# MPI with routing on rank 0

The `feature/root-only-routing` experiment starts at
`93e50f76b3afa5be470216283fbc783917b1ba03`. It ports the timing report from
`bfdc8e8f7558d9af6e45f3e90cc7a98d01851c07` without adopting its distributed
river-routing implementation.

For a coupled run with more than one MPI rank, rank 0 owns the complete river
network and computes one SUMMA GRU. The remaining ranks divide the other GRUs
evenly. Keeping one GRU avoids changing SUMMA initialization, forcing, and output
paths that require a nonempty local domain. A single-rank run computes all GRUs.
Uncoupled SUMMA partitioning is unchanged. Coupled runs require at least as many
GRUs as MPI ranks; partial-domain `-g` runs remain unsupported.

During initialization, GRU IDs are gathered to rank 0. Only rank 0 reads and
allocates the routing topology, remapping, parameters, and state. At each land
model timestep, `MPI_Gatherv` collects GRU runoff in the same global order. Rank 0
then runs the existing serial routing kernels, including routing substeps, and
writes routing fields into its SUMMA output file. Workers can advance to their
next physics step after the gather completes. There are no exchanges of river
boundary states or gathers of routing output. Direct runoff mapping matches IDs
rather than assuming identical array order, and treats inactive/missing runoff
as zero using the standalone mizuRoute threshold.

The design reserves almost one core for routing. It therefore needs enough
workers and enough routing work to offset the loss of one physics worker. In
particular, two ranks have only one worker for almost all GRUs. Root also retains
the full routing network in memory, so this experiment does not partition that
memory across ranks.

## Build

On the cluster with Environment Modules:

```bash
bash build/cmake/build.mpi.fir.bash
```

The wrapper loads GCC, Open MPI, NetCDF Fortran, and OpenBLAS, then builds
`bin/summa.exe` and `bin/summa_mpi.exe`. `build/cmake/build.mpi.bash` also supports
an already configured compiler environment.

## Timing interpretation

The report separates GRU physics from `routing (inclusive)`. Inclusive routing
contains the runoff gather, its waiting time, and, on rank 0, river computation.
Workers' routing time is principally time inside the gather. The `runoff MPI
exchange` row is a subset of inclusive routing, not another cost to add to it.
Boundary and output MPI rows are zero in this architecture. The unaveraged
`Root routing computation (s)` footer measures the serial river calculation
after the gather, including remapping and routing substeps. Routing file writes
are accounted for in `output/stats`.

Mean percentages describe average time across ranks, which can hide the root's
cost. Use the maximum `total elapsed` for scaling comparisons. It includes
initialization and the final synchronization, and excludes MPI startup, the
timing report, cleanup sleep, and MPI finalization. `runs.csv` additionally
records subprocess wall time.

## Bow at Banff validation and benchmark

`utils/regression/test_root_routing.py` uses the adjacent `bow_at_banff` dataset:
51 GRUs, 118 HRUs, and 51 river segments. It creates isolated configuration and
output directories without modifying the input dataset. Python needs NumPy and
netCDF4. For example:

```bash
python utils/regression/test_root_routing.py \
  --work-dir /tmp/bow-root-routing \
  --hours 168 --repeats 3 --ranks 1 2 3 4 5 6 7 8
```

Optionally supply `--reference-bin /path/to/previous/summa_mpi.exe` to measure a
previous implementation alongside the new one. Preserve its matching shared
library and library search path. Each version is compared with the new serial
executable. All MPI runs are sequential, use at most eight ranks, and default
to one thread per rank with Open MPI core binding. Use an appropriate scheduler
allocation for those ranks.

The script compares every deterministic timestep and daily output variable,
including runoff and routed discharge, with exact array equality after
concatenating GRU/HRU partitions. It excludes `wallClockTime` and the old
`seg_owner` metadata, whose values intentionally depend on execution. It also
checks methods 3, 4, and 5 with 900-second routing substeps, direct ID mapping,
remapping, and buffered output. CSV and JSON files retain measurements and
median speedup/efficiency for review. These short catchment tests do not establish
performance or memory requirements for continental simulations.

Measured results for the initial experiment are in the
[Bow at Banff benchmark report](../benchmarks/root-routing-bow-at-banff-2026-10-26/README.md).
