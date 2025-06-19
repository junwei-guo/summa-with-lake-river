

!This module stores single lake information.
!This module contains the interface to update the lake processes.
module SingleLake
      !use ode_jaco
      !use matrix_module
      !use numSolver
      !use rksuite_90
      !use rksuite_90_prec
      !use minpack_module
      !use minpack_capi
      !use numerical_differentiation_module
      !use WaterAirProperties
      use summa_mpi
      use flake_derivedtypes
      use data_parameters
      ! INTEGER, PARAMETER       ::                                         &
      !     ireals    = SELECTED_REAL_KIND (12,200),                       &
      !                   ! number of desired significant digits for
      !                   ! real variables
      !                   ! corresponds to 8 byte real variables

      !     iintegers = KIND  (1)
      !                   ! kind-type parameter of the integer values
      !                   ! corresponds to the default integers

      type :: SGLake
          character(len=50) :: name
          integer ::  hruID       !
          integer ::  gruID
          logical ::  activeLake  ! Define if this lake is active

          integer ::  localLakeIdx
          integer ::  localGRUIdx
          integer ::  localHRUIdx
          integer ::  Hylak_id    ! Unique HydroLake ID

          !Lake information
          real ::     hruArea             ! Fixed at the program initialization
          real ::     area
          real ::     depth
          real ::     volume
          real ::     maxVolume           ! Maximum volume of the lake
          real ::     latitude 
          real ::     longitude
          !real ::     fetch
  
  
          !Meteorological forcing data from the HRU
          real ::     time              ! time since time reference       (s)
          real ::     pptrate           ! precipitation rate              (kg m-2 s-1)
          real ::     airTemp           ! air temperature                 (K)
          real ::     spechum           ! specific humidity               (g/g)
          real ::     windspd           ! windspeed                       (m/s)
          real ::     SWRadAtm          ! downwelling shortwave radiaiton (W m-2)
          real ::     LWRadAtm          ! downwelling longwave radiation  (W m-2)
          real ::     airpres           ! pressure                        (Pa)


          !Derived forcing data
          real ::     rainfall          ! computed rainfall rate (kg m-2 s-1)
          real ::     snowfall          ! computed snowfall rate (kg m-2 s-1)
          
          
          !Coupling to other land unit
          real ::     groundFlowRate      !positive if the flow is coming into the lake, m3/s
          real ::     groundFlowTemp
          real ::     surfaceRunOffRate   !positive if the flow is coming into the lake, m3/s
          real ::     surfaceRunOffTemp
  
  
          !Coupling to other HRU (need MPI)


          !Sediment

          real(kind=ireals) :: &
          depth_bs                          , & ! Depth of the thermally active layer of the bottom sediments [m]
          T_bs                               ! Temperature at the outer edge of 
                                                ! the thermally active layer of the bottom sediments [K]


          !Iterative properties
          real(kind=ireals) :: &
          albedo_water                        , & ! Water surface albedo with respect to the solar radiation
          albedo_ice                          , & ! Ice surface albedo with respect to the solar radiation
          albedo_snow                             ! Snow surface albedo with respect to the solar radiation

          TYPE (opticpar_medium):: & 
          opticpar_water                       , & ! Optical characteristics of water
          opticpar_ice                         , & ! Optical characteristics of ice
          opticpar_snow                            ! Optical characteristics of snow 
  
  
          !Output data
          real(kind=ireals) :: &
          T_snow_out                        , & ! Temperature at the air-snow interface [K] 
          T_ice_out                         , & ! Temperature at the snow-ice or air-ice interface [K]
          T_mnw_out                         , & ! Mean temperature of the water column [K]
          T_wML_out                         , & ! Mixed-layer temperature [K]
          T_bot_out                         , & ! Temperature at the water-bottom sediment interface [K]
          T_B1_out                          , & ! Temperature at the bottom of the upper layer of the sediments [K]
          C_T_out                           , & ! Shape factor (thermocline)
          h_snow_out                        , & ! Snow thickness [m]
          h_ice_out                         , & ! Ice thickness [m]
          h_ML_out                          , & ! Thickness of the mixed-layer [m]
          H_B1_out                          , & ! Thickness of the upper layer of bottom sediments [m]
          T_sfc_n                               ! Updated surface temperature [K]  

  
  
  
      !Parameter constants 
      contains
          procedure :: UpdateInternalProcesses
   
      end type SGLake
  
      contains
  
 
      subroutine UpdateInternalProcesses(this,timeIncrement)
          class(SGLake) :: this
          real(kind = ireals), intent(in) :: timeIncrement
 
          real(kind = ireals) ::&

          !  Input (procedure arguments)
          dMsnowdt_in                       , & ! The rate of snow accumulation [kg m^{-2} s^{-1}]
          I_atm_in                          , & ! Solar radiation flux at the surface [W m^{-2}]
          Q_atm_lw_in                       , & ! Long-wave radiation flux from the atmosphere [W m^{-2}]
          height_u_in                       , & ! Height above the lake surface where the wind speed is measured [m]
          height_tq_in                      , & ! Height where temperature and humidity are measured [m]
          U_a_in                            , & ! Wind speed at z=height_u_in [m s^{-1}]
          T_a_in                            , & ! Air temperature at z=height_tq_in [K]
          q_a_in                            , & ! Air specific humidity at z=height_tq_in
          P_a_in                            , & ! Surface air pressure [N m^{-2} = kg m^{-1} s^{-2}]

          depth_w                           , & ! The lake depth [m]
          fetch                             , & ! Typical wind fetch [m]
          depth_bs                          , & ! Depth of the thermally active layer of the bottom sediments [m]
          T_bs                              , & ! Temperature at the outer edge of 
                                                ! the thermally active layer of the bottom sediments [K]
          par_Coriolis                      , & ! The Coriolis parameter [s^{-1}]
          del_time                          , & ! The model time step [s]

          T_snow_in                        , & ! Temperature at the air-snow interface [K] 
          T_ice_in                         , & ! Temperature at the snow-ice or air-ice interface [K]
          T_mnw_in                         , & ! Mean temperature of the water column [K]
          T_wML_in                         , & ! Mixed-layer temperature [K]
          T_bot_in                         , & ! Temperature at the water-bottom sediment interface [K]
          T_B1_in                          , & ! Temperature at the bottom of the upper layer of the sediments [K]
          C_T_in                           , & ! Shape factor (thermocline)
          h_snow_in                        , & ! Snow thickness [m]
          h_ice_in                         , & ! Ice thickness [m]
          h_ML_in                          , & ! Thickness of the mixed-layer [m]
          H_B1_in                          , & ! Thickness of the upper layer of bottom sediments [m]
          T_sfc_p                          , & ! Surface temperature at the previous time step [K]  

          !  Input/Output (procedure arguments)

          albedo_water                        , & ! Water surface albedo with respect to the solar radiation
          albedo_ice                          , & ! Ice surface albedo with respect to the solar radiation
          albedo_snow                         , & ! Snow surface albedo with respect to the solar radiation




          !  Output (procedure arguments)


          T_snow_out                        , & ! Temperature at the air-snow interface [K] 
          T_ice_out                         , & ! Temperature at the snow-ice or air-ice interface [K]
          T_mnw_out                         , & ! Mean temperature of the water column [K]
          T_wML_out                         , & ! Mixed-layer temperature [K]
          T_bot_out                         , & ! Temperature at the water-bottom sediment interface [K]
          T_B1_out                          , & ! Temperature at the bottom of the upper layer of the sediments [K]
          C_T_out                           , & ! Shape factor (thermocline)
          h_snow_out                        , & ! Snow thickness [m]
          h_ice_out                         , & ! Ice thickness [m]
          h_ML_out                          , & ! Thickness of the mixed-layer [m]
          H_B1_out                          , & ! Thickness of the upper layer of bottom sediments [m]
          T_sfc_n                               ! Updated surface temperature [K]  

          TYPE (opticpar_medium):: & 
          opticpar_water                       , & ! Optical characteristics of water
          opticpar_ice                         , & ! Optical characteristics of ice
          opticpar_snow                         ! Optical characteristics of snow 


          ! Forcing data
          dMsnowdt_in   = this%snowfall
          !dMsnowdt_in   = 0.0             !Turn off the snow accumulation rate. 
          I_atm_in      = this%SWRadAtm 
          Q_atm_lw_in   = this%LWRadAtm
          height_u_in   = 2.0             !Unknown, just put 2.0 m here.
          height_tq_in  = 2.0             !Unknown, just put 2.0 m here.
          U_a_in        = this%windspd
          T_a_in        = this%airTemp
          q_a_in        = this%spechum
          P_a_in        = this%airpres

          depth_w       = this%depth
          fetch         = sqrt(this%area/3.14) * 2.0    !Assuming the lake is a circle, and the fetch is its diameter
          depth_bs      = this%depth_bs                           !Turn off the sediment layer thermal module
          T_bs          = this%T_bs                               !Turn off the sediment layer thermal module
          par_Coriolis  = 2.0*7.2921D-5*sind(this%latitude)        !Read latitude of the lake
          del_time      = timeIncrement



          ! Initail conditions
          T_snow_in    = this%T_snow_out
          T_ice_in     = this%T_ice_out
          T_mnw_in     = this%T_mnw_out
          T_wML_in     = this%T_wML_out
          T_bot_in     = this%T_bot_out
          T_B1_in      = this%T_B1_out
          C_T_in       = this%C_T_out
          h_snow_in    = this%h_snow_out
          h_ice_in     = this%h_ice_out
          h_ML_in      = this%h_ML_out
          H_B1_in      = this%H_B1_out
          T_sfc_p      = this%T_sfc_n


          ! Initial properties
          ! albedo_water = this%albedo_water
          ! albedo_snow  = this%albedo_snow
          ! albedo_ice   = this%albedo_ice

          ! opticpar_water  = this%opticpar_water
          ! opticpar_snow   = this%opticpar_snow
          ! opticpar_ice    = this%opticpar_ice


          call flake_interface ( dMsnowdt_in, I_atm_in, Q_atm_lw_in, height_u_in, height_tq_in,     &
                U_a_in, T_a_in, q_a_in, P_a_in,                                    &
                
                depth_w, fetch, depth_bs, T_bs, par_Coriolis, del_time,            &
                T_snow_in,  T_ice_in,  T_mnw_in,  T_wML_in,  T_bot_in,  T_B1_in,   &
                C_T_in,  h_snow_in,  h_ice_in,  h_ML_in,  H_B1_in, T_sfc_p,        &
                
                albedo_water,   albedo_ice,   albedo_snow,                         &
                opticpar_water, opticpar_ice, opticpar_snow,                       &

                T_snow_out, T_ice_out, T_mnw_out, T_wML_out, T_bot_out, T_B1_out,  & 
                C_T_out, h_snow_out, h_ice_out, h_ML_out, H_B1_out, T_sfc_n )


          ! Updated properties
          ! this%albedo_water = albedo_water
          ! this%albedo_snow  = albedo_snow
          ! this%albedo_ice   = albedo_ice

          ! this%opticpar_water = opticpar_water
          ! this%opticpar_snow  = opticpar_snow
          ! this%opticpar_ice   = opticpar_ice
      

          ! Updated output parameters
          this%T_snow_out   = T_snow_out
          this%T_ice_out    = T_ice_out
          this%T_mnw_out    = T_mnw_out
          this%T_wML_out    = T_wML_out
          this%T_bot_out    = T_bot_out
          this%T_B1_out     = T_B1_out
          this%C_T_out      = C_T_out
          this%h_snow_out   = h_snow_out
          this%h_ice_out    = h_ice_out
          this%h_ML_out     = h_ML_out
          this%H_B1_out     = H_B1_out
          this%T_sfc_n      = T_sfc_n
          

      end subroutine UpdateInternalProcesses
 
          
 
  
 
  

end module SingleLake
  
!This module stores all lake information.
module SubGridLake
      !implicit none
      use NetCDF
      use SingleLake
      USE globalData,only:gru_struc                              ! gru->hru mapping structure
      USE globalData,only:index_map                              ! hru->gru mapping structure
      USE globalData,only:time_meta
      USE globalData,only:forc_meta
      USE globalData,only:attr_meta
      USE summaFileManager,only:OUTPUT_PATH,OUTPUT_PREFIX         ! define output file
      USE globalData,only:output_fileSuffix                       ! suffix for the output & state files (optional summa argument)
      USE summa_type,only:summa1_type_dec                       ! master summa data type
      USE var_lookup,only:iLookFORCE
      USE var_lookup,only:iLookATTR
      USE var_lookup,only:iLookflux
      !Information for all lakes
      type :: LakeSystem
          integer     :: nLake
          type(SGLake),allocatable :: allSGLakes(:)
      contains 
          procedure :: InitAllLake              !Set lake id, gru id, hru id, depth, surface area, etc.
          procedure :: LakeInitialConditions    !Set initial condition of the lakes. 
          procedure :: GetLakeByHRUID
          procedure :: SetForcing
          procedure :: UpdateAllLake
          procedure :: LandCoupling
          procedure :: WriteResult
      end type LakeSystem
    contains
  
  
      subroutine InitAllLake(this)
          class(LakeSystem):: this
          integer :: nHRU, nGRU, nLake, nLakeInFile,localLakeIdx
          integer :: ncid        ! NetCDF file ID
          integer :: varid       ! Variable ID
          integer :: retval      ! Return value for error handling
          real, allocatable :: lakeID(:), Hylak_id(:), lakeGRUID(:),lakeHRUID(:),lakeElev(:),lakeArea(:),lakeDepth(:),lakeVol(:),lakeLat(:),lakeLong(:)


          character(len=*), parameter :: filename = 'lakeinput.nc'
      
          ! Open the NetCDF file in read-only mode
          retval = nf90_open(filename, NF90_NOWRITE, ncid)


          ! Get the variable ID for the variable we want to read

          retval = nf90_inq_dimid(ncid, 'lakeID', varid)
          retval = nf90_inquire_dimension(ncid, varid, len=nLakeInFile)
          allocate(lakeID(nLakeInFile))
          retval = nf90_get_var(ncid, varid, lakeID)

          allocate(Hylak_id(nLakeInFile))
          allocate(lakeGRUID(nLakeInFile))
          allocate(lakeHRUID(nLakeInFile))
          allocate(lakeElev(nLakeInFile))
          allocate(lakeArea(nLakeInFile))          
          allocate(lakeDepth(nLakeInFile))
          allocate(lakeVol(nLakeInFile))
          allocate(lakeLat(nLakeInFile))
          allocate(lakeLong(nLakeInFile))

          retval = nf90_inq_varid(ncid, 'Hylak_id', varid)
          retval = nf90_get_var(ncid, varid, Hylak_id)
          retval = nf90_inq_varid(ncid, 'GRU_ID', varid)
          retval = nf90_get_var(ncid, varid, lakeGRUID)
          retval = nf90_inq_varid(ncid, 'HRU_ID', varid)
          retval = nf90_get_var(ncid, varid, lakeHRUID)
          retval = nf90_inq_varid(ncid, 'Elevation', varid)
          retval = nf90_get_var(ncid, varid, lakeElev)
          retval = nf90_inq_varid(ncid, 'Lake_area', varid)
          retval = nf90_get_var(ncid, varid, lakeArea)
          retval = nf90_inq_varid(ncid, 'Depth_avg', varid)
          retval = nf90_get_var(ncid, varid, lakeDepth)
          retval = nf90_inq_varid(ncid, 'Vot_total', varid)
          retval = nf90_get_var(ncid, varid, lakeVol)
          retval = nf90_inq_varid(ncid, 'Pour_lat', varid)
          retval = nf90_get_var(ncid, varid, lakeLat)
          retval = nf90_inq_varid(ncid, 'Pour_long', varid)
          retval = nf90_get_var(ncid, varid, lakeLong)

          ! Close the NetCDF file
          retval = nf90_close(ncid)


          nHRU = sum(gru_struc%hruCount)  
          nGRU = size(gru_struc)
          nLake = 0

          !Determin nLake in this MPI rank:

          do iGRU=1,nGRU
            do iHRU = 1,  gru_struc(iGRU)%hruCount
                do iLake = 1, nLakeInFile
                  if (gru_struc(iGRU)%gru_id ==lakeGRUID(iLake) .and. gru_struc(iGRU)%hruInfo(iHRU)%hru_id ==lakeHRUID(iLake)) then 
                    nLake = nLake+1
                  endif 
                end do 
            end do 
          end do

          print *,  "[MSG from rank@("//trim(num2str(idx_rank))//", "//trim(num2str(num_rank))//")]: nLake = "//trim(num2str(nLake))



          allocate(this%allSGLakes(nLake))
          this%nLake=nLake
          localLakeIdx = 0
          do iGRU=1,nGRU
            do iHRU = 1,  gru_struc(iGRU)%hruCount
                do iLake = 1, nLakeInFile
                  if (gru_struc(iGRU)%gru_id ==lakeGRUID(iLake) .and. gru_struc(iGRU)%hruInfo(iHRU)%hru_id ==lakeHRUID(iLake)) then 
                      localLakeIdx=localLakeIdx+1

                      this%allSGLakes(localLakeIdx)%gruID            = gru_struc(iGRU)%gru_id
                      this%allSGLakes(localLakeIdx)%hruID            = gru_struc(iGRU)%hruInfo(iHRU)%hru_id 
                      this%allSGLakes(localLakeIdx)%localLakeIdx     = localLakeIdx
                      this%allSGLakes(localLakeIdx)%localGRUIdx      = iGRU
                      this%allSGLakes(localLakeIdx)%localHRUIdx      = iHRU

                      this%allSGLakes(localLakeIdx)%Hylak_id         = Hylak_id(iLake)
                      this%allSGLakes(localLakeIdx)%area             = lakeArea(iLake) * 1.0D6
                      this%allSGLakes(localLakeIdx)%depth            = lakeDepth(iLake)
                      this%allSGLakes(localLakeIdx)%maxVolume        = lakeVol(iLake)  * 1.0D6
                      this%allSGLakes(localLakeIdx)%latitude         = lakeLat(iLake)
                      this%allSGLakes(localLakeIdx)%longitude        = lakeLong(iLake)






                  endif 
                end do 
            end do 
          end do
      end subroutine InitAllLake
  

      subroutine LakeInitialConditions(this)
        class(LakeSystem):: this
        real  :: T0, H0



        ! Read initial conditions from file or set ICs to some values
        do iLake = 1,this%nLake
          T0 = this%allSGLakes(iLake)%airTemp
          H0 = 0.0

          this%allSGLakes(iLake)%T_snow_out   = T0
          this%allSGLakes(iLake)%T_ice_out    = T0
          this%allSGLakes(iLake)%T_mnw_out    = T0
          this%allSGLakes(iLake)%T_wML_out    = T0
          this%allSGLakes(iLake)%T_bot_out    = T0
          this%allSGLakes(iLake)%T_B1_out     = T0


          this%allSGLakes(iLake)%C_T_out     = 0.5
          this%allSGLakes(iLake)%h_snow_out  = H0
          this%allSGLakes(iLake)%h_ice_out   = H0
          this%allSGLakes(iLake)%h_ML_out    = this%allSGLakes(iLake)%depth*0.1
          this%allSGLakes(iLake)%H_B1_out    = H0
          this%allSGLakes(iLake)%T_sfc_n     = T0


          this%allSGLakes(iLake)%T_bs        = T0
          this%allSGLakes(iLake)%depth_bs    = this%allSGLakes(iLake)%depth*0.01


          ! Initialize water, snow and ice properties

          ! this%allSGLakes(iLake)%albedo_water    = 0.1
          ! this%allSGLakes(iLake)%albedo_snow     = 0.6
          ! this%allSGLakes(iLake)%albedo_ice      = 0.6

          ! this%allSGLakes(iLake)%opticpar_water    = 0.1   !Unknown!!!
          ! this%allSGLakes(iLake)%opticpar_snow     = 0.6   !Unknown!!!
          ! this%allSGLakes(iLake)%opticpar_ice      = 0.6   !Unknown!!!


          this%allSGLakes(iLake)%groundFlowRate       = 0.0
          this%allSGLakes(iLake)%groundFlowTemp       = 0.0
          this%allSGLakes(iLake)%surfaceRunOffRate    = 0.0
          this%allSGLakes(iLake)%groundFlowRate       = 0.0


        end do 


      end subroutine 
  
      subroutine SetForcing(this,summa1_struc,modelTimeStep)
        class(LakeSystem) :: this
        type(SGLake) :: lake
        type(summa1_type_dec),intent(inout)   :: summa1_struc       ! master summa data structure
        integer :: localGRUIdx, localHRUIdx,modelTimeStep
        real,dimension(:),allocatable :: forcing,attribute
        
        real    :: lakeGridSize
        real,dimension(:),allocatable :: lakeDepthVec
        real,dimension(:),allocatable :: lakeTempVec
        do iLake = 1,this%nLake
            localGRUIdx = this%allSGLakes(iLake)%localGRUIdx
            localHRUIdx = this%allSGLakes(iLake)%localHRUIdx
            forcing     = summa1_struc%forcStruct%gru(localGRUIdx)%hru(localHRUIdx)%var(:)
            attribute   = summa1_struc%attrStruct%gru(localGRUIdx)%hru(localHRUIdx)%var(:)

            this%allSGLakes(iLake)%time     = forcing(iLookFORCE%time)
            this%allSGLakes(iLake)%pptrate  = forcing(iLookFORCE%pptrate)
            this%allSGLakes(iLake)%airTemp  = forcing(iLookFORCE%airtemp)
            this%allSGLakes(iLake)%spechum  = forcing(iLookFORCE%spechum)
            this%allSGLakes(iLake)%windspd  = forcing(iLookFORCE%windspd)
            this%allSGLakes(iLake)%SWRadAtm = forcing(iLookFORCE%SWRadAtm)
            this%allSGLakes(iLake)%LWRadAtm = forcing(iLookFORCE%LWRadAtm)
            this%allSGLakes(iLake)%airpres  = forcing(iLookFORCE%airpres)
            this%allSGLakes(iLake)%hruArea  = attribute(iLookATTR%HRUarea)


            this%allSGLakes(iLake)%rainfall     =   &
                summa1_struc%fluxStruct%gru(localGRUIdx)%hru(localHRUIdx)%var(iLookFLUX%scalarRainfall)%dat(1)
            this%allSGLakes(iLake)%snowfall     =   &
                summa1_struc%fluxStruct%gru(localGRUIdx)%hru(localHRUIdx)%var(iLookFLUX%scalarSnowfall)%dat(1)

        end do 

      


      end subroutine SetForcing
  
      subroutine GetLakeByHRUID(this, HRUID, lake)
          class(LakeSystem) :: this
          integer, intent(in) :: HRUID
          type(SGLake), intent(out) :: lake
          logical :: found
  
          found = .false.
          do i = 1, size(this%allSGLakes)
              if (this%allSGLakes(i)%hruID == HRUID) then
                  lake = this%allSGLakes(i)
                  found = .true.
                  return
              end if
          end do
  
      end subroutine GetLakeByHRUID
  
  
      subroutine UpdateAllLake(this, timeIncrement)
          class(LakeSystem):: this
          real(8), intent(in) :: timeIncrement
          do iLake = 1, this%nLake
              call this%allSGLakes(iLake)%UpdateInternalProcesses(timeIncrement)
          end do
      end subroutine UpdateAllLake


  
      subroutine LandCoupling(this, timeIncrement)
        class(LakeSystem):: this
        real(ireals), intent(in) :: timeIncrement

    


      end subroutine LandCoupling

      subroutine WriteResult(this,summa1_struc,modelTimeStep,numtim)
        class(LakeSystem):: this

        character(len = 200) :: fn
        type(summa1_type_dec),intent(inout)   :: summa1_struc       ! master summa data structure
        integer :: modelTimeStep,numtim,nlake
        integer :: ncid, time_dim, lake_dim,airtemp_dim

        integer :: varid,Hylak_varid,gru_varid,hru_varid,area_varid,depth_varid
        ! real, allocatable :: lakeID(:), Hylak_id(:), lakeGRUID(:),lakeHRUID(:),lakeArea(:),lakeDepth(:),airTemp(:),&
        !                     T_sfc(:),h_snow(:),h_ice(:),T_mnw(:),T_wML(:)

        real, dimension(this%nLake) :: lakeID, Hylak_id, lakeGRUID,lakeHRUID,lakeArea,lakeDepth,airTemp,&
                            T_sfc,h_snow,h_ice,T_mnw,T_wML,C_T,h_ML,T_bot
        nlake = this%nLake
        ! allocate(lakeID(nlake))
        ! allocate(Hylak_id(nlake))
        ! allocate(lakeGRUID(nlake))
        ! allocate(lakeHRUID(nlake))
        ! allocate(lakeArea(nlake))
        ! allocate(lakeDepth(nlake))
        ! allocate(airTemp(nlake))
        ! allocate(T_sfc(nlake))

        do ilake = 1, nlake
          lakeID(ilake)     = this%allSGLakes(ilake)%localLakeIdx 
          Hylak_id(ilake)   = this%allSGLakes(ilake)%Hylak_id 
          lakeGRUID(ilake)  = this%allSGLakes(ilake)%gruID 
          lakeHRUID(ilake)  = this%allSGLakes(ilake)%hruID 
          lakeArea(ilake)   = this%allSGLakes(ilake)%area 
          lakeDepth(ilake)  = this%allSGLakes(ilake)%depth 

          airTemp(ilake)    = this%allSGLakes(ilake)%airTemp 
          T_sfc(ilake)      = this%allSGLakes(ilake)%T_sfc_n 
          h_snow(ilake)     = this%allSGLakes(ilake)%h_snow_out 
          h_ice(ilake)      = this%allSGLakes(ilake)%h_ice_out 
          T_mnw(ilake)      = this%allSGLakes(ilake)%T_mnw_out 
          T_wML(ilake)      = this%allSGLakes(ilake)%T_wML_out 
          C_T(ilake)        = this%allSGLakes(ilake)%C_T_out 
          h_ML(ilake)       = this%allSGLakes(ilake)%h_ML_out 
          T_bot(ilake)       = this%allSGLakes(ilake)%T_bot_out 

        end do 

 
        write(fn, '(A,A,"_lakeOutput_G"I0,"-"I0,".nc")') trim(OUTPUT_PATH),trim(OUTPUT_PREFIX),sGRU_this_rank, sGRU_this_rank+nGRU_this_rank-1

        ! Initialize the NetCDF file and define dimensions and variable

        if (modelTimeStep == 1) then
          call check(nf90_create(fn, NF90_CLOBBER, ncid))
          !call check(nf90_def_dim(ncid, "time",       numtim,           time_dim))
          call check(nf90_def_dim(ncid, "time",       nf90_unlimited,    time_dim))
          call check(nf90_def_dim(ncid, "nlake",      nlake,             lake_dim))

          call check(nf90_def_var(ncid, "airtemp",    NF90_REAL8,       [lake_dim, time_dim], varid))
          call check(nf90_def_var(ncid, "T_sfc",      NF90_REAL8,       [lake_dim, time_dim], varid))
          call check(nf90_def_var(ncid, "h_snow",     NF90_REAL8,       [lake_dim, time_dim], varid))
          call check(nf90_def_var(ncid, "h_ice",      NF90_REAL8,       [lake_dim, time_dim], varid))
          call check(nf90_def_var(ncid, "T_mnw",      NF90_REAL8,       [lake_dim, time_dim], varid))
          call check(nf90_def_var(ncid, "T_wML",      NF90_REAL8,       [lake_dim, time_dim], varid))
          call check(nf90_def_var(ncid, "C_T",        NF90_REAL8,       [lake_dim, time_dim], varid))
          call check(nf90_def_var(ncid, "h_ML",       NF90_REAL8,       [lake_dim, time_dim], varid))
          call check(nf90_def_var(ncid, "T_bot",      NF90_REAL8,       [lake_dim, time_dim], varid))


          call check(nf90_def_var(ncid, "Hylak_id",   NF90_INT,          lake_dim,            Hylak_varid))
          call check(nf90_def_var(ncid, "gruID",      NF90_INT,          lake_dim,            gru_varid))
          call check(nf90_def_var(ncid, "hruID",      NF90_INT,          lake_dim,            hru_varid))
          call check(nf90_def_var(ncid, "area",       NF90_REAL8,        lake_dim,            area_varid))
          call check(nf90_def_var(ncid, "depth",      NF90_REAL8,        lake_dim,            depth_varid))

          call check(nf90_enddef(ncid))

          !Put basic info:
          call check(nf90_put_var(ncid, Hylak_varid,  Hylak_id))
          call check(nf90_put_var(ncid, gru_varid,    lakeGRUID))
          call check(nf90_put_var(ncid, hru_varid,    lakeHRUID))
          call check(nf90_put_var(ncid, area_varid,   lakeArea))
          call check(nf90_put_var(ncid, depth_varid,  lakeDepth))


          call check(nf90_close(ncid))
        end if


        call check(nf90_open(fn, NF90_WRITE, ncid))


        call check(nf90_inq_varid(ncid, "airtemp", varid))
        call check(nf90_put_var(ncid, varid, airtemp, start=[1, modelTimeStep], count=[nlake, 1]))


        call check(nf90_inq_varid(ncid, "T_sfc", varid))
        call check(nf90_put_var(ncid, varid, T_sfc, start=[1, modelTimeStep], count=[nlake, 1]))

        call check(nf90_inq_varid(ncid, "h_snow", varid))
        call check(nf90_put_var(ncid, varid, h_snow, start=[1, modelTimeStep], count=[nlake, 1]))


        call check(nf90_inq_varid(ncid, "h_ice", varid))
        call check(nf90_put_var(ncid, varid, h_ice, start=[1, modelTimeStep], count=[nlake, 1]))


        call check(nf90_inq_varid(ncid, "T_mnw", varid))
        call check(nf90_put_var(ncid, varid, T_mnw, start=[1, modelTimeStep], count=[nlake, 1]))


        call check(nf90_inq_varid(ncid, "T_wML", varid))
        call check(nf90_put_var(ncid, varid, T_wML, start=[1, modelTimeStep], count=[nlake, 1]))


        call check(nf90_inq_varid(ncid, "C_T", varid))
        call check(nf90_put_var(ncid, varid, C_T, start=[1, modelTimeStep], count=[nlake, 1]))

        call check(nf90_inq_varid(ncid, "h_ML", varid))
        call check(nf90_put_var(ncid, varid, h_ML, start=[1, modelTimeStep], count=[nlake, 1]))

        call check(nf90_inq_varid(ncid, "T_bot", varid))
        call check(nf90_put_var(ncid, varid, T_bot, start=[1, modelTimeStep], count=[nlake, 1]))


        call check(nf90_close(ncid))

        contains

        ! NetCDF error checking subroutine
        subroutine check(status)
            integer, intent(in) :: status
            if (status /= nf90_noerr) then
                print *, "NetCDF error: ", trim(nf90_strerror(status))
                stop "Program terminated due to a NetCDF error."
            end if
        end subroutine check

      end subroutine WriteResult

      

  
  
end module SubGridLake
  
