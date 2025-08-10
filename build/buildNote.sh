
# -------------------------------
# On Anvil:
module --force purge
module load gcc  openmpi
module load netcdf-fortran
module load openblas
echo "Loaded modules:"
module list

git clone --branch dev-lake-river --single-branch https://github.com/junwei-guo/summa-with-lake-river.git
cd summa-with-lake-river/build
make mach=anvil

# -------------------------------
# On ComputeCanada Fir:





# -------------------------------
# On UofC ARC:


# -------------------------------
#On Fedora Linux:
git clone --branch dev-lake-river --single-branch https://github.com/junwei-guo/summa-with-lake-river.git
cd summa-with-lake-river/build

make mach=xps
