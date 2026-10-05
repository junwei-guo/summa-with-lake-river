# MPI SUMMA–mizuRoute coupling

The coupled MPI executable advances SUMMA GRUs and river reaches on the same
communicator. Routing runs on all ranks: it is not a serial routing operation on
rank zero. Existing serial runs remain supported.

## Build

Enable both `USE_MPI=ON` and `USE_MIZUROUTE=ON`. For this Homebrew installation:

```bash
FC=/opt/homebrew/bin/gfortran LIBRARY_LINKS=-llapack \
  bash build/cmake/build.mpi.bash
```

Run this command from the repository, or use the absolute path to the build
script. The script builds `bin/summa.exe` and `bin/summa_mpi.exe` using
`build/cmake_mpi`; it does not delete other build trees. Set `BUILD_TYPE=Debug`
for bounds checks or `BUILD_JOBS` to change build concurrency. With
`LIBRARY_LINKS` unset, the existing CMake OpenBLAS discovery is used.

## Run

From the parent Bow at Banff project directory:

```bash
bash run-bb-mpi.sh 4
```

Equivalent command:

```bash
mpiexec -n 4 summa-mizuroute/bin/summa_mpi.exe \
  -m bow_at_banff/fileManager.txt -c bow_at_banff/mizu.toml -s coupling_mpi
```

The existing file manager controls dates and output location. Use no more ranks
than SUMMA GRUs (51 for Bow at Banff). The `-g` partial-domain option remains
rejected, because omitted upstream GRUs can invalidate routing. This change does
not add routing restart support or enable previously unsupported routing methods.
Select exactly one of the existing supported methods: `"3"` (Eulerian kinematic
wave), `"4"` (Muskingum–Cunge), or `"5"` (diffusive wave).

## Ownership and communication

1. Gather GRU IDs and counts using the communicator and actual rank ownership
   supplied by SUMMA. Gather each rank's runoff vector once per SUMMA timestep.
2. Assign each reach to the rank supplying the largest share of its local basin
   area, using the existing remapping and HRU-to-reach weights. Ties select the
   lowest rank. With the Bow at Banff one-to-one GRU/basin mapping, reach ownership
   exactly matches the corresponding SUMMA GRU ownership. Basinless connectors
   inherit downstream ownership where available, then upstream ownership, with
   rank zero as the final fallback.
3. Build upstream-to-downstream dependency stages. Each rank advances only its
   owned reaches at each stage. Broadcast discharge for cross-rank river edges
   before downstream processing, separately for each routing substep. Internal
   method state stays on the rank advancing that reach.
4. Combine the owned reach discharges for output. Rank zero alone writes the full
   network's routing variables and any augmented topology file.

The no-remapping path also matches basin IDs to GRU IDs explicitly; it no longer
assumes that NetCDF arrays have identical ordering. Remapping must cover all
routing basins and refer to GRUs in the selected full domain.

## Output

Each rank writes its existing local SUMMA output, with a GRU-range suffix. Only
rank zero's timestep file contains the complete routing network. For the default
Bow at Banff prefix and four ranks this is:

```
output/run1_coupling_mpi_G01-13_timestep.nc
```

Routing variables remain `Q_reach(time, seg, method)` in m3/s and
`q_basin(time, mizu_hru)` in m/s. New `seg_owner(seg)` gives the zero-based MPI rank
advancing each reach, and the global attribute `routing_mpi_size` records the
number of ranks. The `seg` and `mizu_hru` coordinates identify the global routing
domain, which is larger than rank zero's local SUMMA `gru` dimension. Do not
concatenate routing variables from the other rank files. Daily files retain
existing SUMMA-only behavior.

## Validation

The integration regression uses the Bow at Banff inputs in the parent project
and writes to a new, isolated directory:

```bash
python3 utils/regression/test_coupled_mpi.py \
  --work-dir ../output/mpi-regression
```

Dependencies: Python `numpy`, `netCDF4`, and `mpiexec`. The directory must not
already exist, so previous validation output cannot be overwritten. The suite
compares complete 24-hour runs bit for bit, checks basin-to-GRU mapping and reach
ownership, checks finite and unmasked routing output, and exercises:

- serial and 1, 2, 4 MPI ranks, including uneven partitions;
- each supported routing method with 900-second substeps;
- buffered output and per-timestep output;
- the configured outlet subset and rank-zero topology generation;
- differing ID order with remapping disabled;
- a rank-local output failure, with a timeout to detect MPI hangs.

Model errors in the MPI driver now abort the MPI job rather than leaving other
ranks blocked in collectives.

## Current scaling limits

This implementation replicates topology, allocated routing state/ghost storage,
and output buffers across ranks. Runoff is all-gathered, boundary discharges are
broadcast, and discharge output is all-reduced. It distributes reach computation,
but does not reduce per-rank routing storage in proportion to rank count. A long
main stem limits concurrency, and preserving SUMMA's partition can leave routing
work unevenly distributed. Correctness tests are not a speedup benchmark;
large-domain scaling needs separate measurement. Sparse ghost storage and
neighbor-only communication would be further work.
