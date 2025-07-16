!! ======================================================================================================
!! mizuRoute interface to SUMMA
!!
!! ======================================================================================================
MODULE route_interface


  ! ****************************************************
  ! provide access to external data, subroutines
  ! ****************************************************
  ! numeric types
  USE nrtype_mizuroute                                       ! variable types, etc.
  ! subroutines: model set up
  USE model_setup,         ONLY: init_model       ! model setupt - reading control file, populate metadata, read parameter file
  USE model_setup,         ONLY: init_data        ! initialize river reach data
  USE model_setup,         ONLY: update_time      ! Update simulation time information at each time step
  ! subroutines: routing
  USE main_route_module,   ONLY: main_route       ! main routing routine
  ! subroutines: model I/O
  USE get_runoff        ,  ONLY: get_hru_runoff   !
  USE write_simoutput,     ONLY: prep_output      !
  USE write_simoutput,     ONLY: output           !
  USE write_restart,       ONLY: main_restart     ! write netcdf restart file
  USE model_finalize,      ONLY: finalize
  USE model_finalize,      ONLY: handle_error
  USE globalData_mizuRoute,  ONLY: runoff_data             ! data structure to hru runoff data
  USE dataTypes,           ONLY: runoff                 ! runoff data type

  !SUMMA modules
  USE globalData,only:gru_struc                              ! gru->hru mapping structure
  USE globalData,only:index_map                              ! hru->gru mapping structure
  USE globalData,only:time_meta
  USE globalData,only:forc_meta
  USE globalData,only:attr_meta
  USE globalData,only:bvar_meta
  USE data_types,only:gru_doubleVec
  !USE summa_type,only:bvarStruct
  USE summaFileManager,only:OUTPUT_PATH,OUTPUT_PREFIX         ! define output file
  USE globalData,only:output_fileSuffix                       ! suffix for the output & state files (optional summa argument)
  USE summa_type,only:summa1_type_dec                       ! master summa data type
  USE var_lookup,only:iLookFORCE
  USE var_lookup,only:iLookATTR
  USE var_lookup,only:iLookflux
  USE var_lookup,only:iLookBVAR

  implicit none

  public :: route_stand_alone_run
  public :: route_sync_run_initialize
  public :: route_sync_run_stepping

  ! ******
  ! define variables
  ! ************************
  character(len=strLen)         :: cfile_name          ! name of the control file
  integer(i4b)                  :: ierr                ! error code
  character(len=strLen)         :: cmessage            ! error message of downwind routine
  integer(i4b)                  :: iens = 1
  logical(lgt)                  :: finished=.false.
  !Timing
  integer*8                     :: cr, startTime, endTime
  real(dp)                      :: elapsedTime

  integer(i4b),dimension(:),allocatable :: summa_hru_id_mapping

  CONTAINS 

    SUBROUTINE route_stand_alone_run()
      ! ******
      ! system_clock rate
      call system_clock(count_rate=cr)

      ! ******
      ! get command-line argument defining the full path to the control file
      ! ***********************************
      !call getarg(1,cfile_name)
      if(len_trim(cfile_name)==0) call handle_error(50,'need to supply name of the control file as a command-line argument')

      ! *****
      ! *** model setup
      !    - read control files and namelist
      !    - broadcast to all processors
      ! ************************
      call init_model(cfile_name, ierr, cmessage)
      if(ierr/=0) call handle_error(ierr, cmessage)

      ! *****
      ! *** data initialization
      !    - river topology, properties, river network domain decomposition
      !    - runoff data (datetime, domain)
      !    - runoff remapping data
      !    - channel states
      ! ***********************************
      call init_data(ierr, cmessage)
      if(ierr/=0) call handle_error(ierr, cmessage)

      ! ***********************************
      ! start of time-stepping simulation
      ! ***********************************
      do while (.not.finished)

        call prep_output(ierr, cmessage)
        if(ierr/=0) call handle_error(ierr, cmessage)

        call system_clock(startTime)
        call get_hru_runoff(ierr, cmessage)
        if(ierr/=0) call handle_error(ierr, cmessage)
        call system_clock(endTime)
        elapsedTime = real(endTime-startTime, kind(dp))/real(cr)
        write(*,"(A,1PG15.7,A)") '   elapsed-time [read_ro] = ', elapsedTime, ' s'

        call system_clock(startTime)
        call main_route(iens, ierr, cmessage)
        if(ierr/=0) call handle_error(ierr, cmessage)
        call system_clock(endTime)
        elapsedTime = real(endTime-startTime, kind(dp))/real(cr)
        write(*,"(A,1PG15.7,A)") '   elapsed-time [routing] = ', elapsedTime, ' s'

        call system_clock(startTime)
        call output(ierr, cmessage)
        if(ierr/=0) call handle_error(ierr, cmessage)
        call system_clock(endTime)
        elapsedTime = real(endTime-startTime, kind(dp))/real(cr)
        write(*,"(A,1PG15.7,A)") '   elapsed-time [output] = ', elapsedTime, ' s'

        call main_restart(ierr, cmessage)
        if(ierr/=0) call handle_error(ierr, cmessage)

        call update_time(finished, ierr, cmessage)
        if(ierr/=0) call handle_error(ierr, cmessage)

      end do

      call finalize()

    end subroutine route_stand_alone_run


    SUBROUTINE route_sync_run_initialize()
      ! ******
      ! system_clock rate
      call system_clock(count_rate=cr)

      ! ******
      ! get command-line argument defining the full path to the control file
      ! ***********************************
      !call getarg(1,cfile_name)
      if(len_trim(cfile_name)==0) call handle_error(50,'need to supply name of the control file as a command-line argument')

      ! *****
      ! *** model setup
      !    - read control files and namelist
      !    - broadcast to all processors
      ! ************************
      call init_model(cfile_name, ierr, cmessage)
      if(ierr/=0) call handle_error(ierr, cmessage)

      ! *****
      ! *** data initialization
      !    - river topology, properties, river network domain decomposition
      !    - runoff data (datetime, domain)
      !    - runoff remapping data
      !    - channel states
      ! ***********************************
      call init_data(ierr, cmessage)
      if(ierr/=0) call handle_error(ierr, cmessage)

    end SUBROUTINE route_sync_run_initialize


    SUBROUTINE route_sync_run_stepping(summa1_struc,modelTimeStep,numtim,data_step)
      type(summa1_type_dec),intent(inout)   :: summa1_struc       ! master summa data structure
      integer(i4b),intent(in)               :: modelTimeStep      ! current model time step
      integer(i4b),intent(in)               :: numtim             ! number of model time steps
      real(dp),intent(in)                   :: data_step          ! time step size
      type(runoff)                          :: runoff1
      type(runoff)                          :: runoff2


      call prep_output(ierr, cmessage)
      if(ierr/=0) call handle_error(ierr, cmessage)


      runoff1 = runoff_data
      !call system_clock(startTime)
      call get_hru_runoff(ierr, cmessage)
      if(ierr/=0) call handle_error(ierr, cmessage)
      ! call system_clock(endTime)
      ! elapsedTime = real(endTime-startTime, kind(dp))/real(cr)
      ! write(*,"(A,1PG15.7,A)") '   elapsed-time [read_ro] = ', elapsedTime, ' s'
      runoff2 = runoff_data

      call pass_summa_runoff_to_mizuroute(summa1_struc,modelTimeStep,numtim,data_step)

      !call system_clock(startTime)
      call main_route(iens, ierr, cmessage)
      if(ierr/=0) call handle_error(ierr, cmessage)
      ! call system_clock(endTime)
      ! elapsedTime = real(endTime-startTime, kind(dp))/real(cr)
      ! write(*,"(A,1PG15.7,A)") '   elapsed-time [routing] = ', elapsedTime, ' s'

      !call system_clock(startTime)
      call output(ierr, cmessage)
      if(ierr/=0) call handle_error(ierr, cmessage)
      ! call system_clock(endTime)
      ! elapsedTime = real(endTime-startTime, kind(dp))/real(cr)
      ! write(*,"(A,1PG15.7,A)") '   elapsed-time [output] = ', elapsedTime, ' s'

      call main_restart(ierr, cmessage)
      if(ierr/=0) call handle_error(ierr, cmessage)

      call update_time(finished, ierr, cmessage)
      if(ierr/=0) call handle_error(ierr, cmessage)

    END SUBROUTINE route_sync_run_stepping


    SUBROUTINE pass_summa_runoff_to_mizuroute(summa1_struc,modelTimeStep,numtim,data_step)
      type(summa1_type_dec),intent(inout)   :: summa1_struc       ! master summa data structure
      integer(i4b),intent(in)               :: modelTimeStep      ! current model time step
      integer(i4b),intent(in)               :: numtim             ! number of model time steps
      real(dp),intent(in)                   :: data_step          ! time step size
      integer(i4b)                          :: nHRU,nGRU,jGRU,iGRU,gru_id
      real(8),dimension(:),allocatable      :: summa_runoff

      type(runoff)                          :: runoff0
      integer(i4b),dimension(:),allocatable :: hru_id_runoff


      nHRU = sum(gru_struc%hruCount)  
      nGRU = size(gru_struc)

      runoff0 = runoff_data
      hru_id_runoff = runoff_data%hru_id

      ! Create mapping table to avoid searching indices at each time step:
      if (modelTimeStep==1) then
        allocate(summa_hru_id_mapping(nGRU))
        summa_hru_id_mapping = -1
        do iGRU=1,nGRU 
          gru_id = hru_id_runoff(iGRU)
          do jGRU=1,nGRU 
            if (gru_struc(jGRU)%gru_id==gru_id) then
                  summa_hru_id_mapping(iGRU)=jGRU
            end if
          end do 
        end do 
      end if 

      do iGRU=1,nGRU
        jGRU=summa_hru_id_mapping(iGRU)
        summa_runoff=summa1_struc%bvarStruct%gru(jGRU)%var(iLookBVAR%averageRoutedRunoff)%dat(:)
        !summa_runoff=summa1_struc%bvarStat%gru(jGRU)%var(iLookBVAR%averageRoutedRunoff)%dat(:)

        runoff_data%basinRunoff(iGRU) = summa_runoff(1)

      end do 

    END SUBROUTINE
END MODULE route_interface


! **************************************************************************************************
! TODO:
! - avoid read run1_timestep.nc in initialization
! - pass time stepping info from summa to mizuroute instead of reading from mizuroute.control file



! **************************************************************************************************
! GUESSING:
! - Hello team, I have a question about mizuRoute. I am wondering if HRU in mizuRoute is equivalent to GRU in summa? It appears that mizuRoute utilizes the averaged HRU runoff within the same GRU.
