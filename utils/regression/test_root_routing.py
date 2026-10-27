#!/usr/bin/env python3
"""Validate and benchmark root-owned routing on Bow at Banff (at most 8 ranks)."""
import argparse
import csv
from datetime import datetime, timedelta
import hashlib
import json
import os
from pathlib import Path
import re
import statistics
import subprocess
import time

import numpy as np
from netCDF4 import Dataset


def prepare_case(data, case, hours, method="3", dt=3600, remap=True, buffered=False):
    case.mkdir(parents=True, exist_ok=False)
    fm = (data / "fileManager.txt").read_text()
    start = datetime.strptime(re.search(r"simStartTime\s+'([^']+)'", fm)[1], "%Y-%m-%d %H:%M")
    end = (start + timedelta(hours=hours - 1)).strftime("%Y-%m-%d %H:%M")
    for key, value in {"simEndTime": end, "settingsPath": str(data) + "/",
                       "forcingPath": str(data / "forcing") + "/", "outputPath": str(case) + "/"}.items():
        fm = re.sub(rf"{key}\s+'[^']+'", f"{key} '{value}'", fm)
    if buffered:
        decisions = (data / "modelDecisions.txt").read_text()
        decisions = re.sub(r"^write_buff\s+.*$", "", decisions, flags=re.M)
        (case / "decisions.txt").write_text(decisions + "\nwrite_buff writeFullSeries\n")
        relative = os.path.relpath(case / "decisions.txt", data)
        fm = re.sub(r"decisionsFile\s+'[^']+'", f"decisionsFile '{relative}'", fm)
    (case / "fileManager.txt").write_text(fm)
    config = (data / "mizu.toml").read_text()
    for key in ("namelist_path", "hfabric_path"):
        config = re.sub(rf"^{key}\s*=.*$", f'{key} = "{data}/"', config, flags=re.M)
    config = re.sub(r'^methods\s*=.*$', f'methods = "{method}"', config, flags=re.M)
    config = re.sub(r'^dt\s*=.*$', f'dt = {dt}', config, flags=re.M)
    config = re.sub(r'^hfabric_newfile\s*=.*\n', "", config, flags=re.M)
    config = re.sub(r'^seg_outlet\s*=.*$', "seg_outlet = -9999", config, flags=re.M)
    if not remap:
        config = config.split("[remapping]")[0]
    (case / "mizu.toml").write_text(config)


def read_results(case, hours, expected_ranks):
    """Compare both timestep and daily fields, concatenated by global GRU range."""
    results = {}
    # The existing full-series writer produces timestep output only.
    buffered = (case / "decisions.txt").exists() and "writeFullSeries" in (case / "decisions.txt").read_text()
    frequencies = ("timestep",) if buffered else ("timestep", "day")
    for frequency in frequencies:
        parts, axes, shared = {}, {}, {}
        paths = sorted(case.glob(f"*_{frequency}.nc"))
        assert len(paths) == expected_ranks, (case, frequency, len(paths), expected_ranks)
        routing_files = 0
        for path in paths:
            with Dataset(path) as nc:
                expected_steps = hours if frequency == "timestep" else hours // 24
                assert len(nc.dimensions["time"]) == expected_steps, path
                if "Q_reach" in nc.variables:
                    routing_files += 1
                    assert len(nc.dimensions["seg"]) == 51
                for key, var in nc.variables.items():
                    if key in {"wallClockTime", "seg_owner"}:
                        continue
                    value = np.asarray(var[:])
                    axis = next((i for i, dim in enumerate(var.dimensions) if dim in {"gru", "hru"}), None)
                    if axis is None:
                        if key in shared:
                            np.testing.assert_array_equal(value, shared[key], err_msg=f"{path}:{key}")
                        else:
                            shared[key] = value
                    else:
                        axes[key] = axis
                        parts.setdefault(key, []).append(value)
        assert routing_files == (1 if frequency == "timestep" else 0), (case, frequency, routing_files)
        for key, values in parts.items():
            shared[key] = np.concatenate(values, axis=axes[key])
        assert shared["gruId"].size == 51
        if frequency == "timestep":
            for key in ("Q_reach", "q_basin"):
                assert np.isfinite(shared[key]).all(), (case, key)
        results.update({f"{frequency}:{key}": value for key, value in shared.items()})
    return results


def compare(reference, candidate):
    assert reference.keys() == candidate.keys(), reference.keys() ^ candidate.keys()
    for key in reference:
        np.testing.assert_array_equal(candidate[key], reference[key], err_msg=key)
        assert candidate[key].dtype == reference[key].dtype, key
        assert candidate[key].tobytes() == reference[key].tobytes(), f"bitwise mismatch: {key}"


def run(binary, ranks, case, launcher, timeout, root_design=True):
    command = ([launcher, "-n", str(ranks)] if ranks else []) + [str(binary), "-m", str(case / "fileManager.txt"),
                                                               "-c", str(case / "mizu.toml"), "-s", "check"]
    started = time.perf_counter()
    with (case / "run.log").open("w") as output:
        result = subprocess.run(command, stdout=output, stderr=subprocess.STDOUT, timeout=timeout)
    wall = time.perf_counter() - started
    log = (case / "run.log").read_text()
    assert result.returncode == 0 and "finished simulation successfully" in log, (case, result.returncode)
    if root_design:
        assert log.count("---- Read river network data ---") == 1
        assert log.count("Root routing: rank 0;") == 1
        if ranks > 1:
            assert "Root SUMMA GRUs: 1\n" in log
    rows = []
    pattern = re.compile(r"^\s*(.*?)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s*$")
    timing_text = log.rsplit("MPI timing summary:", 1)[-1] if "MPI timing summary:" in log else ""
    for line in timing_text.splitlines():
        match = pattern.match(line)
        if match:
            rows.append(dict(component=match[1].strip(), mean_s=float(match[2]), min_s=float(match[3]),
                             max_s=float(match[4]), mean_percent_total=float(match[5])))
    totals = next((r for r in rows if r["component"] == "total elapsed"), None)
    root_compute = re.search(r"Root routing computation \(s\):\s*([\d.]+)", log)
    return {"wall_s": wall, "elapsed_s": totals["max_s"] if totals else wall,
            "root_compute_s": float(root_compute[1]) if root_compute else None}, rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    repo = Path(__file__).resolve().parents[2]
    parser.add_argument("--data-dir", type=Path, default=repo.parent / "bow_at_banff")
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--binary", type=Path, default=repo / "bin/summa_mpi.exe")
    parser.add_argument("--serial-binary", type=Path, default=repo / "bin/summa.exe")
    parser.add_argument("--reference-bin", type=Path, help="Optional previous distributed MPI binary")
    parser.add_argument("--ranks", type=int, nargs="+", default=list(range(1, 9)))
    parser.add_argument("--hours", type=int, default=168)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--mpiexec", default="mpiexec")
    parser.add_argument("--timeout", type=int, default=300)
    args = parser.parse_args()
    if not args.ranks or min(args.ranks) < 1 or max(args.ranks) > 8 or 1 not in args.ranks:
        parser.error("Use rank counts from 1 through 8, including 1 as the scaling reference")
    if args.hours < 1 or args.repeats < 1:
        parser.error("hours and repeats must be positive")
    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    # Pin Open MPI ranks to cores for repeatable small-catchment measurements.
    os.environ.setdefault("OMPI_MCA_hwloc_base_binding_policy", "core")
    os.environ.setdefault("OMPI_MCA_rmaps_base_mapping_policy", "core")
    data, work = args.data_dir.resolve(), args.work_dir.resolve()
    work.mkdir(parents=True, exist_ok=False)
    case = work / "serial_reference"
    prepare_case(data, case, args.hours)
    run(args.serial_binary.resolve(), 0, case, args.mpiexec, args.timeout)
    reference = read_results(case, args.hours, 1)
    print(f"Serial reference: {args.hours} hours, 51 GRUs; {len(reference)} output variables", flush=True)
    versions = {"root-routing": args.binary.resolve()}
    if args.reference_bin:
        versions["distributed-reference"] = args.reference_bin.resolve()
    runs, timings = [], []
    for repeat in range(1, args.repeats + 1):
        # Alternate order between repetitions to reduce systematic drift.
        order = list(versions.items())[::1 if repeat % 2 else -1]
        for ranks in sorted(set(args.ranks)):
            for version, binary in order:
                case = work / f"{version}_n{ranks}_r{repeat}"
                prepare_case(data, case, args.hours)
                measured, rows = run(binary, ranks, case, args.mpiexec, args.timeout, version == "root-routing")
                compare(reference, read_results(case, args.hours, ranks))
                runs.append(dict(version=version, ranks=ranks, repeat=repeat, **measured))
                timings.extend(dict(version=version, ranks=ranks, repeat=repeat, **row) for row in rows)
                print(f"PASS {version}: ranks={ranks}, repeat={repeat}, elapsed={measured['elapsed_s']:.3f}s; bit-for-bit identical", flush=True)
    # Check the other routing kernels, substeps, buffering, and differing ID order.
    for method, remap, buffered in [("3", False, False), ("4", True, False), ("5", True, True)]:
        reference_case = work / f"serial_m{method}"
        prepare_case(data, reference_case, 24, method=method, dt=900, remap=remap, buffered=buffered)
        run(args.serial_binary.resolve(), 0, reference_case, args.mpiexec, args.timeout)
        method_reference = read_results(reference_case, 24, 1)
        for ranks in (2, 8):
            case = work / f"mpi{ranks}_m{method}"
            prepare_case(data, case, 24, method=method, dt=900, remap=remap, buffered=buffered)
            run(args.binary.resolve(), ranks, case, args.mpiexec, args.timeout)
            compare(method_reference, read_results(case, 24, ranks))
        print(f"PASS method={method}, dt=900, remap={remap}, buffered={buffered}: serial/2/8 identical", flush=True)
    summary = []
    for version in versions:
        base = statistics.median(r["elapsed_s"] for r in runs if r["version"] == version and r["ranks"] == 1)
        for ranks in sorted(set(args.ranks)):
            values = [r["elapsed_s"] for r in runs if r["version"] == version and r["ranks"] == ranks]
            median = statistics.median(values)
            summary.append(dict(version=version, ranks=ranks, median_s=median, min_s=min(values), max_s=max(values),
                                speedup=base / median, efficiency=base / median / ranks))
    for name, records in [("runs.csv", runs), ("timings.csv", timings), ("scaling.csv", summary)]:
        with (work / name).open("w", newline="") as output:
            writer = csv.DictWriter(output, fieldnames=list(records[0]))
            writer.writeheader(); writer.writerows(records)
    manifest = {"data_dir": str(data), "hours": args.hours, "repeats": args.repeats,
                "ranks": args.ranks, "versions": {k: str(v) for k, v in versions.items()},
                "hostname": os.uname().nodename, "threads_per_rank": 1,
                "binding_policy": os.environ.get("OMPI_MCA_hwloc_base_binding_policy"),
                "binary_sha256": {k: hashlib.sha256(v.read_bytes()).hexdigest() for k, v in versions.items()},
                "serial_reference_variables": sorted(reference), "runs": runs, "scaling": summary}
    (work / "results.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print("PASS: all consistency and scaling checks complete; results:", work, flush=True)
    for row in summary:
        print(row, flush=True)


if __name__ == "__main__":
    main()
