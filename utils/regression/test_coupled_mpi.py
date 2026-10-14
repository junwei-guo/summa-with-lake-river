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


def run_case(root, work, name, ranks, methods, dt, remap=True, buffered=False, augmented=False, mapping="identity"):
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
    with Dataset(root / 'bow_at_banff/attributes.nc') as nc:
        gru_ids = np.asarray(nc['gruId'][:])
    with Dataset(root / 'bow_at_banff/topology.nc') as nc:
        basin_ids = np.asarray(nc['hruId'][:])
    sources = {int(h): [(int(h), 1.0)] for h in basin_ids}
    if mapping != "identity":
        if mapping == "mixed":
            # A dominant local source plus a source on another SUMMA rank.
            sources = {int(h): [(int(h), 0.75), (int(gru_ids[(np.where(gru_ids == h)[0][0] + 17) % 51]), 0.25)]
                       for h in basin_ids}
        elif mapping == "single_source":
            # All routing belongs to rank zero; other ranks have empty routing arrays.
            sources = {int(h): [(int(gru_ids[0]), 1.0)] for h in basin_ids}
        else:
            raise ValueError(mapping)
        remap_file = case / 'mapping.nc'
        with Dataset(remap_file, 'w', format='NETCDF3_CLASSIC') as nc:
            nc.createDimension('polyid', len(basin_ids))
            nc.createDimension('data', sum(map(len, sources.values())))
            nc.createVariable('polyid', 'i4', ('polyid',))[:] = basin_ids
            nc.createVariable('nOverlaps', 'i4', ('polyid',))[:] = [len(sources[int(h)]) for h in basin_ids]
            nc.createVariable('qhru_id', 'i4', ('data',))[:] = [g for h in basin_ids for g, w in sources[int(h)]]
            nc.createVariable('weight', 'f8', ('data',))[:] = [w for h in basin_ids for g, w in sources[int(h)]]
        relative = os.path.relpath(remap_file, root / 'bow_at_banff')
        config = re.sub(r'^remap_file\s*=.*$', f'remap_file = "{relative}"', config, flags=re.M)
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
    expected_owner = {}
    for h, contributors in sources.items():
        shares = np.zeros(ranks or 1)
        for g, weight in contributors:
            shares[ownership[g]] += weight
        expected_owner[h] = int(np.argmax(shares))
    np.testing.assert_array_equal(saved_owners, [expected_owner[int(k)] for k in out['seg']])
    # Verify allocations, including a rank with zero owned reaches/HRUs.
    storage = re.findall(r'Routing local storage rank (\d+): reaches=(\d+), ghost fluxes=(\d+), HRUs=(\d+)',
                         (case / 'run.log').read_text())
    assert len(storage) == (ranks or 1)
    # Startup must not recreate the expensive legacy topology on every rank.
    log_text = (case / 'run.log').read_text()
    assert log_text.count('---- Read river network data ---') == 1
    assert log_text.count('Routing initialization root:') == 1
    compact_ranks = re.findall(r'Routing compact initialization rank (\d+):', log_text)
    assert sorted(map(int, compact_ranks)) == list(range(ranks or 1))
    with Dataset(root / 'bow_at_banff/topology.nc') as nc:
        edges = list(zip(map(int, nc['segId'][:]), map(int, nc['downSegId'][:])))
    for rank_s, nreach, nghost, nhru in storage:
        rank = int(rank_s)
        owned = sum(o == rank for o in expected_owner.values())
        ghosts = sum(expected_owner[u] != rank and expected_owner.get(d) == rank for u, d in edges)
        assert (int(nreach), int(nghost), int(nhru)) == (owned, ghosts, owned), storage

    out['runoff'] = np.stack([runoff[k] for k in sorted(runoff)], axis=1)
    expected = np.stack([sum(w * runoff[g] for g, w in sources[int(k)]) for k in out['mizu_hru']], axis=1)
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


def check_initialization_failure(root, work):
    """A root-reader error must propagate while workers await connectivity."""
    source = work / 'mpi2_kw'
    case = work / 'mpi3_init_failure'
    case.mkdir()
    (case / 'fileManager.txt').write_text(
        (source / 'fileManager.txt').read_text().replace(str(source), str(case)))
    config = re.sub(r'^hfabric_file\s*=.*$', 'hfabric_file = "missing-topology.nc"',
                    (source / 'mizu.toml').read_text(), flags=re.M)
    (case / 'mizu.toml').write_text(config)
    command = [ARGS.mpiexec, '-n', '3', str(root / 'summa-mizuroute/bin/summa_mpi.exe'),
               '-m', str(case / 'fileManager.txt'), '-c', str(case / 'mizu.toml')]
    with (case / 'run.log').open('w') as log:
        result = subprocess.run(command, cwd=root, stdout=log, stderr=subprocess.STDOUT, timeout=30)
    assert result.returncode != 0
    assert 'missing-topology.nc' in (case / 'run.log').read_text()
    print('mpi3_init_failure: root-reader error terminated all ranks without hanging.', flush=True)


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
    for mapping in ('mixed', 'single_source'):
        ref = run_case(root, work, f'serial_{mapping}', 0, '3', 900, mapping=mapping)
        compare(ref, run_case(root, work, f'mpi4_{mapping}', 4, '3', 900, mapping=mapping))
    check_rank_failure(root, work)
    check_initialization_failure(root, work)
    print('PASS: serial/MPI, uneven partitions, all supported methods, substeps, and ID ordering.', flush=True)
