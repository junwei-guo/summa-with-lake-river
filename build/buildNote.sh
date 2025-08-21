
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

# -------------------------------
#On Anvil with the original serial SUMMA:

git clone --branch dev-lake-river --single-branch https://github.com/junwei-guo/summa-with-lake-river.git
cd summa-with-lake-river/build
git checkout 8321632d84135f70a4d6d31e2caa34208b7b9c00

## Anvil setting:
# module --force purge
# module load gcc  openmpi
# module load netcdf-fortran
# module load openblas


F_MASTER = ../
FC = gfortran

EBROOTGCC=${GCC_HOME}
EBROOTOPENBLAS=${OPENBLAS_HOME}
EBROOTNETCDFMINFORTRAN=${NETCDF_FORTRAN_HOME}

FC_EXE = ${EBROOTGCC}/bin/gfortran
CC = ${EBROOTGCC}/bin/gcc

INCLUDES =  -I${EBROOTOPENBLAS}/include -I${EBROOTNETCDFMINFORTRAN}/include 
LIBRARIES = -L${EBROOTNETCDFMINFORTRAN}/lib -lnetcdff  -L${EBROOTOPENBLAS}/lib -lopenblas  