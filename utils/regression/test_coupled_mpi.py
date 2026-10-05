#!/usr/bin/env python3
"""Integration regression using the Bow at Banff data beside this repository.

Build with USE_MPI=ON and USE_MIZUROUTE=ON, then run this script.
Requires numpy, netCDF4 and an MPI launcher. All runs have a timeout.
"""
import argparse
import os
from datetime import datetime, timedelta
from pathlib import Path
import re
import subprocess

import numpy as np
from netCDF4 import Dataset


def run_case(root, work, name, ranks, methods, dt, remap=True, buffered=False, augmented=False):
    case = work / name
    case.mkdir(parents=True, exist_ok=False)
    fm = (root / 'bow_at_banff/fileManager.txt').read_text()
    start = datetime.strptime(re.search(r"simStartTime\s+'([^']+)'", fm)[1], '%Y-%m-%d %H:%M')
    end = (start + timedelta(hours=23)).strftime('%Y-%m-%d %H:%M')
    fm = re.sub(r"simEndTime\s+'[^']+'", f"simEndTime '{end}'", fm)
    fm = re.sub(r"outputPath\s+'[^']+'", f"outputPath '{case}/'", fm)
    if buffered:
        decisions = (root / 'bow_at_banff/modelDecisions.txt').read_text()
        decisions = re.sub(r'^write_buff\s+.*$', '', decisions, flags=re.M)
        (case / 'decisions.txt').write_text(decisions + '\nwrite_buff writeFullSeries\n')
        relative = os.path.relpath(case / 'decisions.txt', root / 'bow_at_banff')
        fm = re.sub(r"decisionsFile\s+'[^']+'", f"decisionsFile '{relative}'", fm)
    (case / 'fileManager.txt').write_text(fm)
    config = (root / 'bow_at_banff/mizu.toml').read_text()
    config = re.sub(r'methods\s*=\s*"[^"]*"', f'methods = "{methods}"', config)
    config = re.sub(r'^dt\s*=.*$', f'dt = {dt}', config, flags=re.M)
    # Keep topology generation isolated; exercise rank-zero writing separately.
    config = re.sub(r'^hfabric_newfile\s*=.*\n', '', config, flags=re.M)
    if augmented:
        relative = os.path.relpath(case / 'topology_new.nc', root / 'bow_at_banff')
        config = config.replace('[hydrofabric]', f'[hydrofabric]\nhfabric_newfile = "{relative}"')
    else:
        config = re.sub(r'^seg_outlet\s*=.*$', 'seg_outlet = -9999', config, flags=re.M)
    if not remap:
        config = config.split('[remapping]')[0]
    (case / 'mizu.toml').write_text(config)
    binary = root / 'summa-mizuroute/bin' / ('summa_mpi.exe' if ranks else 'summa.exe')
    command = ([ARGS.mpiexec, '-n', str(ranks)] if ranks else []) + [str(binary), '-m', str(case/'fileManager.txt'), '-c', str(case/'mizu.toml'), '-s', 'coupling']
    with (case / 'run.log').open('w') as log:
        result = subprocess.run(command, cwd=root, stdout=log, stderr=subprocess.STDOUT, timeout=180)
    if result.returncode or 'finished simulation successfully' not in (case/'run.log').read_text():
        raise RuntimeError(f'{name} failed; see {case / "run.log"}')
    if augmented:
        with Dataset(case / 'topology_new.nc') as nc:
            assert len(nc.dimensions['seg']) == 51
    routed = []
    runoff = {}
    ownership = {}
    saved_owners = None
    for rank, path in enumerate(sorted(case.glob('*_timestep.nc'))):
        with Dataset(path) as nc:
            assert len(nc.dimensions['time']) == 24, path
            for i, gru in enumerate(nc['gruId'][:]):
                assert int(gru) not in runoff
                ownership[int(gru)] = rank
                runoff[int(gru)] = np.asarray(nc['averageRoutedRunoff'][:, i])
            if 'Q_reach' in nc.variables:
                saved_owners = np.asarray(nc['seg_owner'][:])
                assert nc.routing_mpi_size == (ranks or 1)
                fields = {}
                for key in ['Q_reach', 'q_basin', 'seg', 'mizu_hru', 'method', 'time']:
                    v = nc[key][:]
                    assert not np.ma.getmaskarray(v).any() and np.isfinite(v).all(), (path, key)
                    fields[key] = np.asarray(v)
                assert fields['Q_reach'].shape[1] == 51
                routed.append(fields)
    assert len(routed) == 1, 'Only rank zero should write routing output'
    assert len(runoff) == 51
    out = routed[0]
    np.testing.assert_array_equal(saved_owners, [ownership[int(k)] for k in out['seg']])
    out['runoff'] = np.stack([runoff[k] for k in sorted(runoff)], axis=1)
    expected = np.stack([runoff[int(k)] for k in out['mizu_hru']], axis=1)
    np.testing.assert_array_equal(out['q_basin'], expected)
    print(f'{name}: 24 steps; complete routing output; correct basin mapping', flush=True)
    return out


def compare(reference, candidate):
    for key in reference:
        np.testing.assert_array_equal(candidate[key], reference[key], err_msg=key)
    print('  All routing fields and GRU runoff are bit-for-bit identical.', flush=True)


def check_rank_failure(root, work):
    """One failed writer must abort peers waiting for routing communication."""
    source = work / 'mpi2_kw'
    case = work / 'mpi2_rank_failure'
    case.mkdir()
    fm = (source / 'fileManager.txt').read_text().replace(str(source), str(case))
    (case / 'fileManager.txt').write_text(fm)
    (case / 'mizu.toml').write_text((source / 'mizu.toml').read_text())
    # A directory at rank one's output filename forces a local NetCDF error.
    (case / 'run1_coupling_G27-51_timestep.nc').mkdir()
    command = [ARGS.mpiexec, '-n', '2', str(root / 'summa-mizuroute/bin/summa_mpi.exe'),
               '-m', str(case / 'fileManager.txt'), '-c', str(case / 'mizu.toml'), '-s', 'coupling']
    with (case / 'run.log').open('w') as log:
        result = subprocess.run(command, cwd=root, stdout=log, stderr=subprocess.STDOUT, timeout=30)
    assert result.returncode != 0, 'Expected rank-local output failure'
    assert 'ERROR [rank 1]' in (case / 'run.log').read_text()
    print('mpi2_rank_failure: rank-local error terminated all ranks without hanging.', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project-root', type=Path, default=Path(__file__).resolve().parents[3])
    parser.add_argument('--work-dir', type=Path, required=True, help='New directory for regression artifacts')
    parser.add_argument('--mpiexec', default='mpiexec')
    ARGS = parser.parse_args()
    root, work = ARGS.project_root.resolve(), ARGS.work_dir.resolve()
    work.mkdir(parents=True, exist_ok=False)
    ref = run_case(root, work, 'serial_kw', 0, '3', 3600)
    for n in (1, 2, 4):
        compare(ref, run_case(root, work, f'mpi{n}_kw', n, '3', 3600))
    compare(ref, run_case(root, work, 'mpi2_augmented', 2, '3', 3600, augmented=True))
    compare(ref, run_case(root, work, 'mpi2_buffered', 2, '3', 3600, buffered=True))
    for method in ('3', '4', '5'):
        ref = run_case(root, work, f'serial_m{method}_substeps', 0, method, 900)
        compare(ref, run_case(root, work, f'mpi3_m{method}_substeps', 3, method, 900))
    # No-remap must match by ID even though topology and GRUs have different order.
    compare(ref, run_case(root, work, 'mpi2_no_remap', 2, '5', 900, remap=False))
    check_rank_failure(root, work)
    print('PASS: serial/MPI, uneven partitions, all supported methods, substeps, and ID ordering.', flush=True)
