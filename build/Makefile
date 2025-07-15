
F_MASTER = ../
FC = mpif90
ROUTING = /home/guo/Dropbox/summampi/summaSource/mizuRoute-in-summa/route/build

ifeq ($(mach), xps)

	EBROOTNETCDFMINFORTRAN=/usr/lib64/gfortran/modules				#PATH for netcdf
	EBROOTOPENBLAS=/usr/lib64										#PATH for openblas
	EBROOTOPENMPI=/usr/lib64/openmpi								#PATH for openmpi
	FC_EXE = /usr/lib64/openmpi/bin/mpif90
	CC = /usr/bin/gcc
	INCLUDES = -I/usr/include -I/usr/include/openblas -I$(EBROOTNETCDFMINFORTRAN) -I${EBROOTOPENMPI}/lib
	LIBRARIES = -L$(EBROOTNETCDFMINFORTRAN) -lnetcdff -L${EBROOTOPENBLAS} -lblas -L${EBROOTOPENBLAS} -lopenblas  -p  -g -O0 

else ifeq ($(mach), cc)
	FC_EXE = ${EBROOTOPENMPI}/bin/mpif90
	CC = ${EBROOTGCC}/bin/gcc

	INCLUDES =  -I${EBROOTOPENBLAS}/include -I${EBROOTNETCDFMINFORTRAN}/include -I${EBROOTOPENMPI}/lib
	LIBRARIES = -L${EBROOTNETCDFMINFORTRAN}/lib -lnetcdff  -L${EBROOTOPENBLAS}/lib -lopenblas  


else ifeq ($(mach), anvil)

	EBROOTOPENMPI=${MPI_HOME}
	EBROOTGCC=${GCC_HOME}
	EBROOTOPENBLAS=${OPENBLAS_HOME}
	EBROOTNETCDFMINFORTRAN=${NETCDF_FORTRAN_HOME}

	FC_EXE = ${EBROOTOPENMPI}/bin/mpif90
	CC = ${EBROOTGCC}/bin/gcc

	INCLUDES =  -I${EBROOTOPENBLAS}/include -I${EBROOTNETCDFMINFORTRAN}/include -I${EBROOTOPENMPI}/lib
	LIBRARIES = -L${EBROOTNETCDFMINFORTRAN}/lib -lnetcdff  -L${EBROOTOPENBLAS}/lib -lopenblas  

else ifeq  ($(mach), mac)



	EBROOTOPENMPI=/Users/guo/homebrew/Cellar/open-mpi/5.0.7
	EBROOTOPENBLAS=/Users/guo/homebrew/Cellar/openblas/0.3.29
	EBROOTNETCDFMINFORTRAN=/Users/guo/homebrew/Cellar/netcdf-fortran/4.6.2

	FC_EXE = /Users/guo/homebrew/bin/mpif90
	CC = /Users/guo/homebrew/Cellar/gcc/14.1.0_2/bin/gcc-14

	INCLUDES =  -I${EBROOTOPENBLAS}/include -I${EBROOTNETCDFMINFORTRAN}/include -I${EBROOTOPENMPI}/lib
	LIBRARIES = -L${EBROOTNETCDFMINFORTRAN}/lib -lnetcdff  -L${EBROOTOPENBLAS}/lib -lopenblas  

else 

	# Same as xps
	EBROOTNETCDFMINFORTRAN=/usr/lib64/gfortran/modules				#PATH for netcdf
	EBROOTOPENBLAS=/usr/lib64										#PATH for openblas
	EBROOTOPENMPI=/usr/lib64/openmpi								#PATH for openmpi
	FC_EXE = /usr/lib64/openmpi/bin/mpif90
	CC = /usr/bin/gcc
	INCLUDES = -I/usr/include -I/usr/include/openblas -I$(EBROOTNETCDFMINFORTRAN) -I${EBROOTOPENMPI}/lib
	LIBRARIES = -L$(EBROOTNETCDFMINFORTRAN) -lnetcdff -L${EBROOTOPENBLAS} -lblas -L${EBROOTOPENBLAS} -lopenblas  -p  -g -O0 

endif



# define open MP flags
isOpenMP =
FLAGS_OMP =
LIBOPENMP =


ifeq ($(mode), debug)

# Debug runs
FLAGS_NOAH  = -p  -g -O0 -ffree-form -ffree-line-length-none -fmax-errors=0 -fbacktrace -Wno-unused -Wno-unused-dummy-argument
FLAGS_COMM  = -p  -g -O0 -Wall -ffree-line-length-none -fmax-errors=0 -fbacktrace -fcheck=bounds 
FLAGS_SUMMA = -p  -g -O0 -Wall -ffree-line-length-none -fmax-errors=0 -fbacktrace -fcheck=bounds

else 

FLAGS_NOAH = -O3 -ffree-form -ffree-line-length-none -fmax-errors=0 
FLAGS_COMM = -O3 -ffree-line-length-none -fmax-errors=0 
FLAGS_SUMMA = -O3 -ffree-line-length-none -fmax-errors=0

endif

# #MUST include mizuRoute path
# FLAGS_SUMMA+=-I${ROUTING} 

# #Option 1: include mizuRoute .o and .mod files
# #INCLUDES+=-I${ROUTING} 

# #Option 2: link static libraries
# LIBRARIES+=$(ROUTING)/libmizuroute.a
# #========================================================================
# PART 1: Define directory paths
#========================================================================

# Core directory that contains source code
F_KORE_DIR = $(F_MASTER)/build/source

# Location of the compiled modules
MOD_PATH = $(F_MASTER)/build

# Define the directory for the executables
EXE_PATH = $(F_MASTER)/bin

#========================================================================
# PART 2: Assemble all of the SUMMA sub-routines
#========================================================================

# Define directories
DRIVER_DIR = $(F_KORE_DIR)/driver
HOOKUP_DIR = $(F_KORE_DIR)/hookup
NETCDF_DIR = $(F_KORE_DIR)/netcdf
DSHARE_DIR = $(F_KORE_DIR)/dshare
NUMREC_DIR = $(F_KORE_DIR)/numrec
NOAHMP_DIR = $(F_KORE_DIR)/noah-mp
ENGINE_DIR = $(F_KORE_DIR)/engine

#  utilities
SUMMA_NRUTIL= \
		summa_mpi.f90 \
		nrtype.f90 \
		f2008funcs.f90 \
		nr_utility.f90
NRUTIL = $(patsubst %, $(ENGINE_DIR)/%, $(SUMMA_NRUTIL))

#
# Numerical recipes procedures
# NOTE: all numerical recipes procedures are now replaced with free versions
SUMMA_NRPROC= \
		expIntegral.f90 \
		spline_int.f90
NRPROC = $(patsubst %, $(ENGINE_DIR)/%, $(SUMMA_NRPROC))

# Hook-up modules (set files and directory paths)
SUMMA_HOOKUP= \
		ascii_util.f90 \
		summaFileManager.f90
HOOKUP = $(patsubst %, $(HOOKUP_DIR)/%, $(SUMMA_HOOKUP))

# Data modules
SUMMA_DATAMS= \
		multiconst.f90 \
		var_lookup.f90 \
		data_types.f90 \
		globalData.f90 \
		flxMapping.f90 \
		get_ixname.f90 \
		popMetadat.f90 \
		outpt_stat.f90
DATAMS = $(patsubst %, $(DSHARE_DIR)/%, $(SUMMA_DATAMS))

# utility modules
SUMMA_UTILMS= \
		time_utils.f90 \
		mDecisions.f90 \
		snow_utils.f90 \
		soil_utils.f90 \
		updatState.f90 \
		matrixOper.f90
UTILMS = $(patsubst %, $(ENGINE_DIR)/%, $(SUMMA_UTILMS))

# Model guts
SUMMA_MODGUT= \
		MODGUT = $(patsubst %, $(ENGINE_DIR)/%, $(SUMMA_MODGUT))

# Solver
SUMMA_SOLVER= \
		vegPhenlgy.f90 \
		diagn_evar.f90 \
		stomResist.f90 \
		groundwatr.f90 \
		vegSWavRad.f90 \
		vegNrgFlux.f90 \
		ssdNrgFlux.f90 \
		vegLiqFlux.f90 \
		snowLiqFlx.f90 \
		soilLiqFlx.f90 \
		bigAquifer.f90 \
		computFlux.f90 \
		computResid.f90 \
		computJacob.f90 \
		eval8summa.f90 \
		summaSolve.f90 \
		systemSolv.f90 \
		varSubstep.f90 \
		opSplittin.f90 \
		coupled_em.f90 \
		run_oneHRU.f90 \
		run_oneGRU.f90
SOLVER = $(patsubst %, $(ENGINE_DIR)/%, $(SUMMA_SOLVER))

# Define routines for SUMMA preliminaries
SUMMA_PRELIM= \
		conv_funcs.f90 \
		sunGeomtry.f90 \
		convE2Temp.f90 \
		allocspace.f90 \
		checkStruc.f90 \
		childStruc.f90 \
		ffile_info.f90 \
		read_attrb.f90 \
		read_pinit.f90 \
		pOverwrite.f90 \
		read_param.f90 \
		paramCheck.f90 \
		check_icond.f90
PRELIM = $(patsubst %, $(ENGINE_DIR)/%, $(SUMMA_PRELIM))

SUMMA_NOAHMP= \
		module_model_constants.F \
		module_sf_noahutl.F \
		module_sf_noahlsm.F \
		module_sf_noahmplsm.F

NOAHMP = $(patsubst %, $(NOAHMP_DIR)/%, $(SUMMA_NOAHMP))

# Define routines for the SUMMA model runs
SUMMA_MODRUN = \
		indexState.f90 \
		getVectorz.f90 \
		updateVars.f90 \
		var_derive.f90 \
		read_force.f90 \
		derivforce.f90 \
		snowAlbedo.f90 \
		canopySnow.f90 \
		tempAdjust.f90 \
		snwCompact.f90 \
		layerMerge.f90 \
		layerDivide.f90 \
		volicePack.f90 \
		qTimeDelay.f90
MODRUN = $(patsubst %, $(ENGINE_DIR)/%, $(SUMMA_MODRUN))

# Define routines for the solver
SUMMA_MSOLVE = \

# Define NetCDF routines
SUMMA_NETCDF = \
		netcdf_util.f90 \
		def_output.f90 \
		modelwrite.f90 \
		read_icond.f90
NETCDF = $(patsubst %, $(NETCDF_DIR)/%, $(SUMMA_NETCDF))

# ... stitch together common programs
COMM_ALL = $(NRUTIL) $(NRPROC) $(HOOKUP) $(DATAMS) $(UTILMS)

# ... stitch together SUMMA programs
SUMMA_ALL = $(NETCDF) $(PRELIM) $(MODRUN) $(SOLVER)

# Define the driver routine
SUMMA_DRIVER= \
		summa_type.f90 \
		summa_util.f90 \
		summa_alarms.f90 \
		summa_globalData.f90 \
		summa_defineOutput.f90 \
		summa_init.f90 \
		summa_setup.f90 \
		summa_restart.f90 \
		summa_forcing.f90 \
		summa_modelRun.f90 \
		summa_writeOutput.f90 \
		lake/data_parameters.f90\
		lake/flake_configure.f90\
		lake/flake_derivedtypes.f90\
		lake/flake_parameters.f90\
		lake/flake_albedo_ref.f90\
		lake/flake_paramoptic_ref.f90\
		lake/flake.f90\
		lake/SfcFlx.f90\
		lake/src_flake_interface_1D.f90\
		route/nrtype_mizuroute.f90 \
		route/nr_utility_mizuroute.f90 \
		route/ascii_util_mizuroute.f90 \
		route/ncio_utils.f90 \
		route/gamma_func.f90 \
		route/public_var.f90 \
		route/dataTypes.f90 \
		route/var_lookup_mizuroute.f90 \
		route/time_utils_mizuroute.f90 \
		route/datetime_data.f90 \
		route/csv_data.f90 \
		route/gageMeta_data.f90 \
		route/obs_data.f90 \
		route/globalData_mizuroute.f90 \
		route/popMetadat_mizuroute.f90 \
		route/allocation.f90 \
		route/network_topo.f90 \
		route/process_param.f90 \
		route/process_ntopo.f90 \
		route/pfafstetter.f90 \
		route/domain_decomposition.f90 \
		route/remap.f90 \
		route/read_control.f90 \
		route/read_param_mizuroute.f90 \
		route/read_streamSeg.f90 \
		route/write_streamSeg.f90 \
		route/read_runoff.f90 \
		route/get_basin_runoff.f90 \
		route/read_remap.f90 \
		route/read_restart.f90 \
		route/write_restart.f90 \
		route/write_simoutput.f90 \
		route/model_finalize.f90 \
		route/data_assimilation.f90 \
		route/accum_runoff.f90 \
		route/basinUH.f90 \
		route/irf_route.f90 \
		route/kwt_route.f90 \
		route/dfw_route.f90 \
		route/kw_route.f90 \
		route/mc_route.f90 \
		route/main_route.f90 \
		route/model_setup.f90 \
		route/route_interface.f90 \
		summa_lake.f90 \
		summa_routing.f90 \
		summa_driver.f90
DRIVER = $(patsubst %, $(DRIVER_DIR)/%, $(SUMMA_DRIVER))

# Define the executable
DRIVER__EX = summa.exe

# Define version number
VERSIONFILE = $(DRIVER_DIR)/summaversion.inc
BULTTIM = $(shell date)
GITBRCH = $(shell git describe --long --all --always | sed -e's/heads\///')
GITHASH = $(shell git rev-parse HEAD)
VERSION = $(shell git show-ref --tags | grep $GITHASH | sed 's/.*tags\///' | grep . || echo "undefined")

#========================================================================
# PART 3: Checks
#======================================================================
# make sure that the paths are defined. These are just some high level checks
ifndef F_MASTER
 $(error F_MASTER is undefined)
endif
ifndef FC
 $(error FC is undefined: Specify your compiler)
endif
ifndef FC_EXE
 $(error FC_EXE is undefined: Specify your compiler executable)
endif
ifndef FLAGS_SUMMA
 $(error Specify flags for your compiler: $(FC))
endif
ifndef INCLUDES
 $(error INCLUDES is undefined)
endif
ifndef LIBRARIES
 $(error LIBRARIES is undefined)
endif

#========================================================================
# PART 4: compilation
#======================================================================

# Compile
all: compile_noah compile_comm compile_summa link clean install

check:
	$(info)
	$(info Displaying make variables:)
	$(info   F_MASTER   : $(F_MASTER))
	$(info   EXE_PATH   : $(EXE_PATH))
	$(info   FC         : $(FC))
	$(info   INCLUDES   : $(INCLUDES))
	$(info   LIBRARIES  : $(LIBRARIES))
	$(info   FLAGS_NOAH : $(FLAGS_NOAH))
	$(info   FLAGS_COMM : $(FLAGS_COMM))
	$(info   FLAGS_SUMMA: $(FLAGS_SUMMA))
	$(info)

# update version information
update_version:
	echo "character(len=64), parameter     :: summaVersion = '${VERSION}'" > $(VERSIONFILE)
	echo "character(len=64), parameter     :: buildTime = '${BULTTIM}'" >> $(VERSIONFILE)
	echo "character(len=64), parameter     :: gitBranch = '${GITBRCH}'" >> $(VERSIONFILE)
	echo "character(len=64), parameter     :: gitHash = '${GITHASH}'" >> $(VERSIONFILE)

# compile Noah-MP routines
compile_noah:
	$(FC_EXE) $(FLAGS_NOAH) -c $(NRUTIL) $(NOAHMP)

# compile common routines
compile_comm:
	$(FC_EXE) $(FLAGS_COMM) -c $(COMM_ALL) $(INCLUDES)

# compile SUMMA routines
compile_summa: update_version
	$(FC_EXE) $(FLAGS_SUMMA) -c $(SUMMA_ALL) $(DRIVER) \
	$(INCLUDES)

# link routines
# link:
# 	$(FC_EXE) -g *.o ${ROUTING}/src/*.o $(LIBRARIES) -o $(DRIVER__EX)

link:
	$(FC_EXE) -g *.o $(LIBRARIES) -o $(DRIVER__EX)

# Remove object files
clean:
	rm -f *.o
	rm -f *.mod
	rm -f soil_veg_gen_parm__genmod.f90

# Copy the executable to the bin directory
install:
	@mkdir -p $(EXE_PATH)
	@mv $(DRIVER__EX) $(EXE_PATH)
	$(info $(DRIVER__EX) successfully installed in $(EXE_PATH))
