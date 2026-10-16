! Distributed reach execution using the host SUMMA GRU ownership.
! Runtime storage is local; ghost fluxes hold immediate remote upstream discharge.
module routing_parallel
  use nrtype, only: i4b, i8b, dp
  use mizuroute_types, only: routing_parallel_data, river_network_data, routing_gather
  use dataTypes, only: remap, RCHPRP
  use var_lookup, only: ixNTOPO
  use nr_utils, only: indexx, match_index
#ifdef ROUTING_MPI
  use mpi
#endif
  implicit none
  private
  ! Temporary global connectivity uses flat arrays, never per-reach descriptors.
  ! Only root constructs the legacy, allocation-heavy initialization structures.
  type :: routing_connectivity
    integer, allocatable :: order(:), down(:), down_id(:), lake(:), lake_model(:)
    integer, allocatable :: up_start(:), upstream(:), good(:), hru_start(:), hru(:)
    real(dp), allocatable :: weight(:)
  end type
  public :: routing_connectivity, share_initial_network
  public :: collect_runoff_ids, collect_runoff, plan_routing, localize_routing, exchange_level, gather_output
contains
  function wall_seconds() result(seconds)
    real(dp) :: seconds
    integer(i8b) :: ticks, rate
    call system_clock(ticks,rate)
    seconds=real(ticks,dp)/real(rate,dp)
  end function


  ! Root exports flat connectivity, then discards legacy topology immediately.
  ! Workers never allocate var(:)%dat(:) for the global hydrofabric.
  subroutine share_initial_network(net, c, init_error, message)
    use public_var, only: is_lake_sim
    type(river_network_data), intent(inout) :: net
    type(routing_connectivity), intent(out) :: c
    integer, intent(inout) :: init_error
    character(*), intent(inout) :: message
    integer :: n, s, first, last, mpierr
    associate(p=>net%driver%parallel)
#ifdef ROUTING_MPI
    if(p%nproc>1) then
      call MPI_Bcast(init_error,1,MPI_INTEGER,0,p%comm,mpierr)
      if(mpierr/=0) then
        init_error=mpierr; return
      endif
      if(init_error/=0) call MPI_Bcast(message,len(message),MPI_CHARACTER,0,p%comm,mpierr)
    endif
#endif
    if(init_error/=0) return
    if(p%rank==0) then
      n=size(net%core%ntopo)
      allocate(c%order(n),c%down(n),c%down_id(n),c%lake(n),c%lake_model(n), &
               c%up_start(n+1),c%hru_start(n+1))
      c%up_start(1)=1; c%hru_start(1)=1
      c%lake=0; c%lake_model=0
      do s=1,n
        c%order(s)=net%core%topology%ntopo(s)%var(ixNTOPO%rchOrder)%dat(1)
        c%down(s)=net%core%ntopo(s)%DREACHI
        c%down_id(s)=net%core%ntopo(s)%DREACHK
        c%up_start(s+1)=c%up_start(s)+size(net%core%ntopo(s)%UREACHI)
        c%hru_start(s+1)=c%hru_start(s)+size(net%core%ntopo(s)%HRUIX)
        if(is_lake_sim) then
          c%lake(s)=merge(1,0,net%core%ntopo(s)%LAKINLT) + &
                    merge(2,0,net%core%ntopo(s)%ISLAKE) + &
                    merge(4,0,net%core%ntopo(s)%LAKETARGVOL)
          c%lake_model(s)=net%core%ntopo(s)%LAKEMODELTYPE
        endif
      enddo
      allocate(c%upstream(c%up_start(n+1)-1),c%good(c%up_start(n+1)-1), &
               c%hru(c%hru_start(n+1)-1),c%weight(c%hru_start(n+1)-1))
      do s=1,n
        first=c%up_start(s); last=c%up_start(s+1)-1
        c%upstream(first:last)=net%core%ntopo(s)%UREACHI
        c%good(first:last)=merge(1,0,net%core%ntopo(s)%goodBas)
        first=c%hru_start(s); last=c%hru_start(s+1)-1
        c%hru(first:last)=net%core%ntopo(s)%HRUIX
        c%weight(first:last)=net%core%ntopo(s)%HRUWGT
      enddo
      deallocate(net%core%ntopo)
      associate(t=>net%core%topology)
        deallocate(t%hru,t%seg,t%hru2seg,t%ntopo,t%pfaf)
      end associate
      write(*,'(a,i0,a,i0)') 'Routing initialization root: reaches=',n, &
        ', HRUs=',size(net%driver%hru_id)
    endif
    call bcast_integer(c%order,p,init_error)
    if(init_error/=0) return
    call bcast_integer(c%down,p,init_error)
    if(init_error/=0) return
    call bcast_integer(c%down_id,p,init_error)
    if(init_error/=0) return
    call bcast_integer(c%lake_model,p,init_error)
    if(init_error/=0) return
    call bcast_integer(c%lake,p,init_error)
    if(init_error/=0) return
    call bcast_integer(c%up_start,p,init_error)
    if(init_error/=0) return
    call bcast_integer(c%upstream,p,init_error)
    if(init_error/=0) return
    call bcast_integer(c%good,p,init_error)
    if(init_error/=0) return
    call bcast_integer(c%hru_start,p,init_error)
    if(init_error/=0) return
    call bcast_integer(c%hru,p,init_error)
    if(init_error/=0) return
    call bcast_integer(net%driver%hru_id,p,init_error)
    if(init_error/=0) return
    call bcast_integer(net%driver%seg_id,p,init_error)
    if(init_error/=0) return
    call bcast_real(c%weight,p,init_error)
    if(init_error/=0) return
    call bcast_real(net%driver%totArea,p,init_error)
    if(init_error/=0) return
    net%core%topology%n_seg=size(net%driver%seg_id)
    net%core%topology%n_hru=size(net%driver%hru_id)
    net%core%topology%is_initialized=.true.
    write(*,'(a,i0,a,i0)') 'Routing compact initialization rank ',p%rank, &
      ': reaches=',size(c%down)
    end associate
  end subroutine

  subroutine bcast_integer(values,p,ierr)
    integer, allocatable, intent(inout) :: values(:)
    type(routing_parallel_data), intent(in) :: p
    integer, intent(out) :: ierr
    integer :: n
    ierr=0
#ifdef ROUTING_MPI
    if(p%nproc>1) then
      if(p%rank==0) n=size(values)
      call MPI_Bcast(n,1,MPI_INTEGER,0,p%comm,ierr)
      if(ierr/=0) return
      if(p%rank/=0) allocate(values(n))
      call MPI_Bcast(values,n,MPI_INTEGER,0,p%comm,ierr)
    endif
#endif
  end subroutine

  subroutine bcast_real(values,p,ierr)
    real(dp), allocatable, intent(inout) :: values(:)
    type(routing_parallel_data), intent(in) :: p
    integer, intent(out) :: ierr
    integer :: n
    ierr=0
#ifdef ROUTING_MPI
    if(p%nproc>1) then
      if(p%rank==0) n=size(values)
      call MPI_Bcast(n,1,MPI_INTEGER,0,p%comm,ierr)
      if(ierr/=0) return
      if(p%rank/=0) allocate(values(n))
      call MPI_Bcast(values,n,MPI_DOUBLE_PRECISION,0,p%comm,ierr)
    endif
#endif
  end subroutine

  ! Scatter scalar columns with bounded temporary storage and typed MPI buffers.
  subroutine scatter_parameters(net,param,ierr)
    use public_var, only: is_lake_sim
    type(river_network_data), intent(in) :: net
    type(RCHPRP), intent(out) :: param(:)
    integer, intent(out) :: ierr
    real(dp), allocatable :: real_send(:), real_recv(:)
    integer, allocatable :: int_send(:), int_recv(:)
    integer :: n
    associate(p=>net%driver%parallel)
    n=0
    if(p%rank==0) n=size(p%reaches%index)
    allocate(real_send(n),real_recv(size(param)),int_send(n),int_recv(size(param)))
    ierr=0
    if(p%rank==0) real_send=net%core%param(p%reaches%index)%R_SLOPE
    call scatter_real_column(real_send,real_recv,p,ierr)
    if(ierr/=0) return
    param%R_SLOPE=real_recv
    if(p%rank==0) real_send=net%core%param(p%reaches%index)%R_MAN_N
    call scatter_real_column(real_send,real_recv,p,ierr)
    if(ierr/=0) return
    param%R_MAN_N=real_recv
    if(p%rank==0) real_send=net%core%param(p%reaches%index)%R_WIDTH
    call scatter_real_column(real_send,real_recv,p,ierr)
    if(ierr/=0) return
    param%R_WIDTH=real_recv
    if(p%rank==0) real_send=net%core%param(p%reaches%index)%R_DEPTH
    call scatter_real_column(real_send,real_recv,p,ierr)
    if(ierr/=0) return
    param%R_DEPTH=real_recv
    if(p%rank==0) real_send=net%core%param(p%reaches%index)%RLENGTH
    call scatter_real_column(real_send,real_recv,p,ierr)
    if(ierr/=0) return
    param%RLENGTH=real_recv
    if(p%rank==0) real_send=net%core%param(p%reaches%index)%R_STORAGE
    call scatter_real_column(real_send,real_recv,p,ierr)
    if(ierr/=0) return
    param%R_STORAGE=real_recv
    if(p%rank==0) real_send=net%core%param(p%reaches%index)%SIDE_SLOPE
    call scatter_real_column(real_send,real_recv,p,ierr)
    if(ierr/=0) return
    param%SIDE_SLOPE=real_recv
    if(p%rank==0) real_send=net%core%param(p%reaches%index)%FLDP_SLOPE
    call scatter_real_column(real_send,real_recv,p,ierr)
    if(ierr/=0) return
    param%FLDP_SLOPE=real_recv
    if(p%rank==0) real_send=net%core%param(p%reaches%index)%UPSAREA
    call scatter_real_column(real_send,real_recv,p,ierr)
    if(ierr/=0) return
    param%UPSAREA=real_recv
    if(p%rank==0) real_send=net%core%param(p%reaches%index)%BASAREA
    call scatter_real_column(real_send,real_recv,p,ierr)
    if(ierr/=0) return
    param%BASAREA=real_recv
    if(p%rank==0) real_send=net%core%param(p%reaches%index)%TOTAREA
    call scatter_real_column(real_send,real_recv,p,ierr)
    if(ierr/=0) return
    param%TOTAREA=real_recv
    if(p%rank==0) real_send=net%core%param(p%reaches%index)%MINFLOW
    call scatter_real_column(real_send,real_recv,p,ierr)
    if(ierr/=0) return
    param%MINFLOW=real_recv
    if(is_lake_sim) then
      if(p%rank==0) real_send=net%core%param(p%reaches%index)%D03_MaxStorage
      call scatter_real_column(real_send,real_recv,p,ierr)
      if(ierr/=0) return
      param%D03_MaxStorage=real_recv
      if(p%rank==0) real_send=net%core%param(p%reaches%index)%D03_Coefficient
      call scatter_real_column(real_send,real_recv,p,ierr)
      if(ierr/=0) return
      param%D03_Coefficient=real_recv
      if(p%rank==0) real_send=net%core%param(p%reaches%index)%D03_Power
      call scatter_real_column(real_send,real_recv,p,ierr)
      if(ierr/=0) return
      param%D03_Power=real_recv
      if(p%rank==0) real_send=net%core%param(p%reaches%index)%D03_S0
      call scatter_real_column(real_send,real_recv,p,ierr)
      if(ierr/=0) return
      param%D03_S0=real_recv
      if(p%rank==0) real_send=net%core%param(p%reaches%index)%HYP_E_emr
      call scatter_real_column(real_send,real_recv,p,ierr)
      if(ierr/=0) return
      param%HYP_E_emr=real_recv
      if(p%rank==0) real_send=net%core%param(p%reaches%index)%HYP_E_lim
      call scatter_real_column(real_send,real_recv,p,ierr)
      if(ierr/=0) return
      param%HYP_E_lim=real_recv
      if(p%rank==0) real_send=net%core%param(p%reaches%index)%HYP_E_min
      call scatter_real_column(real_send,real_recv,p,ierr)
      if(ierr/=0) return
      param%HYP_E_min=real_recv
      if(p%rank==0) real_send=net%core%param(p%reaches%index)%HYP_E_zero
      call scatter_real_column(real_send,real_recv,p,ierr)
      if(ierr/=0) return
      param%HYP_E_zero=real_recv
      if(p%rank==0) real_send=net%core%param(p%reaches%index)%HYP_Qrate_emr
      call scatter_real_column(real_send,real_recv,p,ierr)
      if(ierr/=0) return
      param%HYP_Qrate_emr=real_recv
      if(p%rank==0) real_send=net%core%param(p%reaches%index)%HYP_Erate_emr
      call scatter_real_column(real_send,real_recv,p,ierr)
      if(ierr/=0) return
      param%HYP_Erate_emr=real_recv
      if(p%rank==0) real_send=net%core%param(p%reaches%index)%HYP_Qrate_prim
      call scatter_real_column(real_send,real_recv,p,ierr)
      if(ierr/=0) return
      param%HYP_Qrate_prim=real_recv
      if(p%rank==0) real_send=net%core%param(p%reaches%index)%HYP_Qrate_amp
      call scatter_real_column(real_send,real_recv,p,ierr)
      if(ierr/=0) return
      param%HYP_Qrate_amp=real_recv
      if(p%rank==0) int_send=net%core%param(p%reaches%index)%HYP_Qrate_phs
      call scatter_int_column(int_send,int_recv,p,ierr)
      if(ierr/=0) return
      param%HYP_Qrate_phs=int_recv
      if(p%rank==0) int_send=merge(1,0,net%core%param(p%reaches%index)%HYP_prim_F)
      call scatter_int_column(int_send,int_recv,p,ierr)
      if(ierr/=0) return
      param%HYP_prim_F=int_recv/=0
      if(p%rank==0) real_send=net%core%param(p%reaches%index)%HYP_A_avg
      call scatter_real_column(real_send,real_recv,p,ierr)
      if(ierr/=0) return
      param%HYP_A_avg=real_recv
      if(p%rank==0) int_send=merge(1,0,net%core%param(p%reaches%index)%HYP_Qsim_mode)
      call scatter_int_column(int_send,int_recv,p,ierr)
      if(ierr/=0) return
      param%HYP_Qsim_mode=int_recv/=0
      if(p%rank==0) real_send=net%core%param(p%reaches%index)%H06_Smax
      call scatter_real_column(real_send,real_recv,p,ierr)
      if(ierr/=0) return
      param%H06_Smax=real_recv
      if(p%rank==0) real_send=net%core%param(p%reaches%index)%H06_alpha
      call scatter_real_column(real_send,real_recv,p,ierr)
      if(ierr/=0) return
      param%H06_alpha=real_recv
      if(p%rank==0) real_send=net%core%param(p%reaches%index)%H06_envfact
      call scatter_real_column(real_send,real_recv,p,ierr)
      if(ierr/=0) return
      param%H06_envfact=real_recv
      if(p%rank==0) real_send=net%core%param(p%reaches%index)%H06_S_ini
      call scatter_real_column(real_send,real_recv,p,ierr)
      if(ierr/=0) return
      param%H06_S_ini=real_recv
      if(p%rank==0) real_send=net%core%param(p%reaches%index)%H06_c1
      call scatter_real_column(real_send,real_recv,p,ierr)
      if(ierr/=0) return
      param%H06_c1=real_recv
      if(p%rank==0) real_send=net%core%param(p%reaches%index)%H06_c2
      call scatter_real_column(real_send,real_recv,p,ierr)
      if(ierr/=0) return
      param%H06_c2=real_recv
      if(p%rank==0) real_send=net%core%param(p%reaches%index)%H06_exponent
      call scatter_real_column(real_send,real_recv,p,ierr)
      if(ierr/=0) return
      param%H06_exponent=real_recv
      if(p%rank==0) real_send=net%core%param(p%reaches%index)%H06_denominator
      call scatter_real_column(real_send,real_recv,p,ierr)
      if(ierr/=0) return
      param%H06_denominator=real_recv
      if(p%rank==0) real_send=net%core%param(p%reaches%index)%H06_c_compare
      call scatter_real_column(real_send,real_recv,p,ierr)
      if(ierr/=0) return
      param%H06_c_compare=real_recv
      if(p%rank==0) real_send=net%core%param(p%reaches%index)%H06_frac_Sdead
      call scatter_real_column(real_send,real_recv,p,ierr)
      if(ierr/=0) return
      param%H06_frac_Sdead=real_recv
      if(p%rank==0) real_send=net%core%param(p%reaches%index)%H06_E_rel_ini
      call scatter_real_column(real_send,real_recv,p,ierr)
      if(ierr/=0) return
      param%H06_E_rel_ini=real_recv
      if(p%rank==0) real_send=net%core%param(p%reaches%index)%H06_I_Jan
      call scatter_real_column(real_send,real_recv,p,ierr)
      if(ierr/=0) return
      param%H06_I_Jan=real_recv
      if(p%rank==0) real_send=net%core%param(p%reaches%index)%H06_I_Feb
      call scatter_real_column(real_send,real_recv,p,ierr)
      if(ierr/=0) return
      param%H06_I_Feb=real_recv
      if(p%rank==0) real_send=net%core%param(p%reaches%index)%H06_I_Mar
      call scatter_real_column(real_send,real_recv,p,ierr)
      if(ierr/=0) return
      param%H06_I_Mar=real_recv
      if(p%rank==0) real_send=net%core%param(p%reaches%index)%H06_I_Apr
      call scatter_real_column(real_send,real_recv,p,ierr)
      if(ierr/=0) return
      param%H06_I_Apr=real_recv
      if(p%rank==0) real_send=net%core%param(p%reaches%index)%H06_I_May
      call scatter_real_column(real_send,real_recv,p,ierr)
      if(ierr/=0) return
      param%H06_I_May=real_recv
      if(p%rank==0) real_send=net%core%param(p%reaches%index)%H06_I_Jun
      call scatter_real_column(real_send,real_recv,p,ierr)
      if(ierr/=0) return
      param%H06_I_Jun=real_recv
      if(p%rank==0) real_send=net%core%param(p%reaches%index)%H06_I_Jul
      call scatter_real_column(real_send,real_recv,p,ierr)
      if(ierr/=0) return
      param%H06_I_Jul=real_recv
      if(p%rank==0) real_send=net%core%param(p%reaches%index)%H06_I_Aug
      call scatter_real_column(real_send,real_recv,p,ierr)
      if(ierr/=0) return
      param%H06_I_Aug=real_recv
      if(p%rank==0) real_send=net%core%param(p%reaches%index)%H06_I_Sep
      call scatter_real_column(real_send,real_recv,p,ierr)
      if(ierr/=0) return
      param%H06_I_Sep=real_recv
      if(p%rank==0) real_send=net%core%param(p%reaches%index)%H06_I_Oct
      call scatter_real_column(real_send,real_recv,p,ierr)
      if(ierr/=0) return
      param%H06_I_Oct=real_recv
      if(p%rank==0) real_send=net%core%param(p%reaches%index)%H06_I_Nov
      call scatter_real_column(real_send,real_recv,p,ierr)
      if(ierr/=0) return
      param%H06_I_Nov=real_recv
      if(p%rank==0) real_send=net%core%param(p%reaches%index)%H06_I_Dec
      call scatter_real_column(real_send,real_recv,p,ierr)
      if(ierr/=0) return
      param%H06_I_Dec=real_recv
      if(p%rank==0) real_send=net%core%param(p%reaches%index)%H06_D_Jan
      call scatter_real_column(real_send,real_recv,p,ierr)
      if(ierr/=0) return
      param%H06_D_Jan=real_recv
      if(p%rank==0) real_send=net%core%param(p%reaches%index)%H06_D_Feb
      call scatter_real_column(real_send,real_recv,p,ierr)
      if(ierr/=0) return
      param%H06_D_Feb=real_recv
      if(p%rank==0) real_send=net%core%param(p%reaches%index)%H06_D_Mar
      call scatter_real_column(real_send,real_recv,p,ierr)
      if(ierr/=0) return
      param%H06_D_Mar=real_recv
      if(p%rank==0) real_send=net%core%param(p%reaches%index)%H06_D_Apr
      call scatter_real_column(real_send,real_recv,p,ierr)
      if(ierr/=0) return
      param%H06_D_Apr=real_recv
      if(p%rank==0) real_send=net%core%param(p%reaches%index)%H06_D_May
      call scatter_real_column(real_send,real_recv,p,ierr)
      if(ierr/=0) return
      param%H06_D_May=real_recv
      if(p%rank==0) real_send=net%core%param(p%reaches%index)%H06_D_Jun
      call scatter_real_column(real_send,real_recv,p,ierr)
      if(ierr/=0) return
      param%H06_D_Jun=real_recv
      if(p%rank==0) real_send=net%core%param(p%reaches%index)%H06_D_Jul
      call scatter_real_column(real_send,real_recv,p,ierr)
      if(ierr/=0) return
      param%H06_D_Jul=real_recv
      if(p%rank==0) real_send=net%core%param(p%reaches%index)%H06_D_Aug
      call scatter_real_column(real_send,real_recv,p,ierr)
      if(ierr/=0) return
      param%H06_D_Aug=real_recv
      if(p%rank==0) real_send=net%core%param(p%reaches%index)%H06_D_Sep
      call scatter_real_column(real_send,real_recv,p,ierr)
      if(ierr/=0) return
      param%H06_D_Sep=real_recv
      if(p%rank==0) real_send=net%core%param(p%reaches%index)%H06_D_Oct
      call scatter_real_column(real_send,real_recv,p,ierr)
      if(ierr/=0) return
      param%H06_D_Oct=real_recv
      if(p%rank==0) real_send=net%core%param(p%reaches%index)%H06_D_Nov
      call scatter_real_column(real_send,real_recv,p,ierr)
      if(ierr/=0) return
      param%H06_D_Nov=real_recv
      if(p%rank==0) real_send=net%core%param(p%reaches%index)%H06_D_Dec
      call scatter_real_column(real_send,real_recv,p,ierr)
      if(ierr/=0) return
      param%H06_D_Dec=real_recv
      if(p%rank==0) int_send=net%core%param(p%reaches%index)%H06_purpose
      call scatter_int_column(int_send,int_recv,p,ierr)
      if(ierr/=0) return
      param%H06_purpose=int_recv
      if(p%rank==0) int_send=merge(1,0,net%core%param(p%reaches%index)%H06_I_mem_F)
      call scatter_int_column(int_send,int_recv,p,ierr)
      if(ierr/=0) return
      param%H06_I_mem_F=int_recv/=0
      if(p%rank==0) int_send=merge(1,0,net%core%param(p%reaches%index)%H06_D_mem_F)
      call scatter_int_column(int_send,int_recv,p,ierr)
      if(ierr/=0) return
      param%H06_D_mem_F=int_recv/=0
      if(p%rank==0) int_send=net%core%param(p%reaches%index)%H06_I_mem_L
      call scatter_int_column(int_send,int_recv,p,ierr)
      if(ierr/=0) return
      param%H06_I_mem_L=int_recv
      if(p%rank==0) int_send=net%core%param(p%reaches%index)%H06_D_mem_L
      call scatter_int_column(int_send,int_recv,p,ierr)
      if(ierr/=0) return
      param%H06_D_mem_L=int_recv
    endif
    end associate
  end subroutine

  subroutine scatter_real_column(send,recv,p,ierr)
    real(dp), intent(in) :: send(:)
    real(dp), intent(out) :: recv(:)
    type(routing_parallel_data), intent(in) :: p
    integer, intent(out) :: ierr
    ierr=0
#ifdef ROUTING_MPI
    if(p%nproc>1) then
      call MPI_Scatterv(send,p%reaches%counts,p%reaches%offsets,MPI_DOUBLE_PRECISION, &
                        recv,size(recv),MPI_DOUBLE_PRECISION,0,p%comm,ierr)
    else
#endif
      recv=send
#ifdef ROUTING_MPI
    endif
#endif
  end subroutine

  subroutine scatter_int_column(send,recv,p,ierr)
    integer, intent(in) :: send(:)
    integer, intent(out) :: recv(:)
    type(routing_parallel_data), intent(in) :: p
    integer, intent(out) :: ierr
    ierr=0
#ifdef ROUTING_MPI
    if(p%nproc>1) then
      call MPI_Scatterv(send,p%reaches%counts,p%reaches%offsets,MPI_INTEGER, &
                        recv,size(recv),MPI_INTEGER,0,p%comm,ierr)
    else
#endif
      recv=send
#ifdef ROUTING_MPI
    endif
#endif
  end subroutine

  subroutine collect_runoff_ids(p, local_ids, ids, ierr, message)
    type(routing_parallel_data), intent(inout) :: p
    integer(i8b), intent(in) :: local_ids(:)
    integer(i8b), allocatable, intent(out) :: ids(:)
    integer, intent(out) :: ierr
    character(*), intent(out) :: message
    integer :: r, n, i
    integer, allocatable :: sorted_ids(:)
    ierr=0; message='collect_runoff_ids/'
    allocate(p%counts(p%nproc), p%offsets(p%nproc))
    n=size(local_ids)
    p%counts=n
#ifdef ROUTING_MPI
    if(p%nproc>1) then
      call MPI_Allgather(n,1,MPI_INTEGER,p%counts,1,MPI_INTEGER,p%comm,ierr)
      if(ierr/=0) return
    endif
#else
    if(p%nproc>1) then
      ierr=20; message=trim(message)//'routing was built without MPI'; return
    endif
#endif
    p%offsets(1)=0
    do r=2,p%nproc
      p%offsets(r)=p%offsets(r-1)+p%counts(r-1)
    enddo
    allocate(ids(sum(p%counts)), p%host_owner(sum(p%counts)))
    do r=1,p%nproc
      p%host_owner(p%offsets(r)+1:p%offsets(r)+p%counts(r))=r-1
    enddo
#ifdef ROUTING_MPI
    if(p%nproc>1) then
      call MPI_Allgatherv(local_ids,n,MPI_INTEGER8,ids,p%counts,p%offsets,MPI_INTEGER8,p%comm,ierr)
      if(ierr/=0) return
    else
#endif
      ids=local_ids
#ifdef ROUTING_MPI
    endif
#endif
    allocate(sorted_ids(size(ids)))
    call indexx(ids,sorted_ids)
    do i=2,size(ids)
      if(ids(sorted_ids(i-1))==ids(sorted_ids(i))) then
        ierr=20; message=trim(message)//'duplicate host GRU ID'; return
      endif
    enddo
  end subroutine

  subroutine collect_runoff(p, local_runoff, runoff, ierr, message)
    type(routing_parallel_data), intent(inout) :: p
    real(dp), intent(in) :: local_runoff(:)
    real(dp), intent(out) :: runoff(:)
    integer, intent(out) :: ierr
    character(*), intent(out) :: message
    real(dp), allocatable :: send(:)
    real(dp) :: started
    ierr=0; message='collect_runoff/'
    send=local_runoff(p%send_index)
#ifdef ROUTING_MPI
    if(p%nproc>1) then
      started=wall_seconds()
      call MPI_Alltoallv(send,p%send_counts,p%send_offsets,MPI_DOUBLE_PRECISION, &
                         runoff,p%recv_counts,p%recv_offsets,MPI_DOUBLE_PRECISION,p%comm,ierr)
      p%runoff_mpi_seconds=p%runoff_mpi_seconds+wall_seconds()-started
      return
    endif
#endif
    runoff=send
  end subroutine

  subroutine plan_routing(net, c, mapping, do_remapping, ierr, message)
    type(river_network_data), intent(inout) :: net
    type(routing_connectivity), intent(in) :: c
    type(remap), intent(in) :: mapping
    logical, intent(in) :: do_remapping
    integer, intent(out) :: ierr
    character(*), intent(out) :: message
    real(dp), allocatable :: share(:)
    integer, allocatable :: level(:), used(:), map_start(:), map_count(:)
    integer :: n, h, k, j, r, s, d, u, pos, l, nlevels, host
    ierr=0; message='plan_routing/'
    associate(p=>net%driver%parallel)
    n=net%core%topology%n_seg
    allocate(p%owner(n),level(n),p%runoff_index(size(net%driver%hru_id)))
    allocate(share(p%nproc), map_start(size(net%driver%hru_id)), map_count(size(net%driver%hru_id)))
    map_start=0; map_count=0
    ! No-remap is still ID-based, never dependent on file ordering.
    if(.not.do_remapping) then
      p%runoff_index=match_index(net%core%runoff%hru_id,int(net%driver%hru_id,i8b),ierr,message)
      if(ierr/=0) return
      if(any(p%runoff_index<=0)) then
        ierr=20; message='plan_routing/routing basin ID missing from host GRUs'; return
      endif
    else
      ! Keep ragged mapping offsets, not a ranks-by-basins dense matrix.
      pos=1
      do h=1,size(mapping%hru_ix)
        j=mapping%hru_ix(h)
        if(j>0 .and. j<=size(map_start)) then
          if(map_start(j)>0 .or. mapping%num_qhru(h)<=0) then
            ierr=20; message=trim(message)//'duplicate or empty remapping basin'; return
          endif
          map_start(j)=pos; map_count(j)=mapping%num_qhru(h)
          do k=pos,pos+map_count(j)-1
            host=mapping%qhru_ix(k)
            if(host<=0 .or. host>size(p%host_owner)) then
              ierr=20; message=trim(message)//'remapping refers to an absent host GRU'; return
            endif
          enddo
        endif
        pos=pos+max(0,mapping%num_qhru(h))
      enddo
      if(any(map_start==0)) then
        ierr=20; message=trim(message)//'remapping does not cover every routing basin'; return
      endif
    endif
    ! A reach belongs to the rank with the largest local contributing area.
    ! For matching GRUs/basins this exactly preserves SUMMA ownership.
    p%owner=-1
    do s=1,n
      share=0._dp
      do k=c%hru_start(s),c%hru_start(s+1)-1
        h=c%hru(k)
        if(h<=0) cycle
        if(do_remapping) then
          do j=map_start(h),map_start(h)+map_count(h)-1
            r=p%host_owner(mapping%qhru_ix(j))+1
            share(r)=share(r)+c%weight(k)*mapping%weight(j)
          enddo
        else
          r=p%host_owner(p%runoff_index(h))+1
          share(r)=share(r)+c%weight(k)
        endif
      enddo
      if(sum(share)>0._dp) p%owner(s)=maxloc(share,dim=1)-1
    enddo
    ! Basinless connectors inherit a downstream owner, then an upstream owner.
    do k=n,1,-1
      s=c%order(k)
      d=c%down(s)
      if(p%owner(s)<0 .and. d>0 .and. d<=n) p%owner(s)=p%owner(d)
    enddo
    level=0
    do k=1,n
      s=c%order(k)
      level(s)=1
      do j=c%up_start(s),c%up_start(s+1)-1
        u=c%upstream(j)
        if(u<=0) cycle
        if(u>n) then
          ierr=20; message=trim(message)//'invalid upstream index'; return
        endif
        if(level(u)==0) then
          ierr=20; message=trim(message)//'network is not in upstream-to-downstream order'; return
        endif
        level(s)=max(level(s),level(u)+1)
        if(p%owner(s)<0) p%owner(s)=p%owner(u)
      enddo
      if(p%owner(s)<0) p%owner(s)=0
    enddo
    ! Serial execution keeps the original routing order in one stage.
    if(p%nproc==1) level=1
    nlevels=maxval(level)
    allocate(p%level_start(nlevels+1),p%level_order(n),used(nlevels))
    used=0
    do s=1,n
      used(level(s))=used(level(s))+1
    enddo
    p%level_start(1)=1
    do l=1,nlevels
      p%level_start(l+1)=p%level_start(l)+used(l)
    enddo
    used=0
    do k=1,n
      s=c%order(k)
      l=level(s); pos=p%level_start(l)+used(l)
      p%level_order(pos)=s; used(l)=used(l)+1
    enddo
    if(p%nproc>1) write(*,'(a,i0,a,i0,a,i0)') 'Routing rank ',p%rank,': owned reaches=', &
                                           count(p%owner==p%rank),', stages=',nlevels
    end associate
  end subroutine

  ! Discard global initialization structures before solver-state allocation.
  ! Owned reaches use indices 1:n_local; upstream ghost fluxes follow them.
  subroutine localize_routing(net, c, mapping, do_remapping, ierr, message)
    use dataTypes, only: RCHPRP, RCHTOPO
    type(river_network_data), intent(inout) :: net
    type(routing_connectivity), intent(in) :: c
    type(remap), intent(inout) :: mapping
    logical, intent(in) :: do_remapping
    integer, intent(out) :: ierr
    character(*), intent(out) :: message
    type(RCHPRP), allocatable :: param(:)
    type(RCHTOPO), allocatable :: topo(:)
    type(remap) :: local_map
    integer, allocatable :: reach(:), basin(:), g2l(:), h2l(:), host2local(:), hosts(:)
    integer, allocatable :: basin_owner(:), order(:), starts(:)
    logical, allocatable :: ghost(:), needed(:)
    integer :: n, nh, ns, nb, ng, ne, s, d, u, h, k, j, r, l, pos, total, first, last
    ierr=0; message='localize_routing/'
    associate(p=>net%driver%parallel)
    n=size(p%owner); nh=size(net%driver%hru_id)
    ns=count(p%owner==p%rank); p%n_local=ns
    allocate(reach(ns),g2l(n),h2l(nh),ghost(n),basin_owner(nh))
    reach=pack([(s,s=1,n)],p%owner==p%rank)
    g2l=0; ghost=.false.; basin_owner=-1
    do s=1,n
      do k=c%hru_start(s),c%hru_start(s+1)-1
        h=c%hru(k)
        if(h<=0) cycle
        if(basin_owner(h)>=0 .and. basin_owner(h)/=p%owner(s)) then
          ierr=20; message=trim(message)//'routing HRU contributes to multiple reach owners'; return
        endif
        basin_owner(h)=p%owner(s)
      enddo
    enddo
    ! Preserve output for disconnected routing HRUs too, with one owner.
    where(basin_owner<0) basin_owner=0
    nb=count(basin_owner==p%rank)
    basin=pack([(h,h=1,nh)],basin_owner==p%rank)
    h2l=0
    do h=1,nb
      h2l(basin(h))=h
    enddo
    do j=1,ns
      s=reach(j); g2l(s)=j
      do k=c%up_start(s),c%up_start(s+1)-1
        u=c%upstream(k)
        if(u>0) then
          if(p%owner(u)/=p%rank) ghost(u)=.true.
        endif
      enddo
    enddo
    ng=count(ghost); p%n_ghost=ng; pos=ns
    do s=1,n
      if(.not.ghost(s)) cycle
      pos=pos+1; g2l(s)=pos
    enddo

    ! Boundary messages retain global traversal order, but only local edges.
    ne=0
    do s=1,n
      d=c%down(s)
      if(d<=0 .or. d>n) cycle
      if(p%owner(s)==p%owner(d)) cycle
      if(p%owner(s)==p%rank .or. p%owner(d)==p%rank) ne=ne+1
    enddo
    allocate(p%edge_local(ne),p%edge_peer(ne),p%edge_stage(ne),p%edge_send(ne))
    pos=0
    do l=1,size(p%level_start)-1
      do k=p%level_start(l),p%level_start(l+1)-1
        s=p%level_order(k); d=c%down(s)
        if(d<=0 .or. d>n) cycle
        if(p%owner(s)==p%owner(d)) cycle
        if(p%owner(s)/=p%rank .and. p%owner(d)/=p%rank) cycle
        pos=pos+1
        p%edge_local(pos)=g2l(s); p%edge_stage(pos)=l
        p%edge_send(pos)=p%owner(s)==p%rank
        p%edge_peer(pos)=merge(p%owner(d),p%owner(s),p%edge_send(pos))
      enddo
    enddo
    allocate(order(ns),starts(size(p%level_start)))
    pos=1
    do l=1,size(starts)-1
      starts(l)=pos
      do k=p%level_start(l),p%level_start(l+1)-1
        s=p%level_order(k)
        if(p%owner(s)/=p%rank) cycle
        order(pos)=g2l(s); pos=pos+1
      enddo
    enddo
    starts(size(starts))=pos
    call move_alloc(order,p%level_order)
    call move_alloc(starts,p%level_start)

    ! Root stores coordinates, not global simulation states or history buffers.
    if(p%rank==0) then
      net%driver%output_seg_id=net%driver%seg_id
      net%driver%output_hru_id=net%driver%hru_id
      net%driver%output_area=net%driver%totArea
      net%driver%output_owner=p%owner
    endif
    call setup_gather(p%comm,p%nproc,reach,p%reaches,ierr)
    if(ierr/=0) return
    call setup_gather(p%comm,p%nproc,basin,p%basins,ierr)
    if(ierr/=0) return

    allocate(param(ns),topo(ns))
    call scatter_parameters(net, param, ierr)
    if(ierr/=0) return
    do j=1,ns
      s=reach(j)
      topo(j)%REACHID=net%driver%seg_id(s)
      topo(j)%DREACHI=c%down(s)
      topo(j)%DREACHK=c%down_id(s)
      topo(j)%RHORDER=c%order(s)
      topo(j)%LAKINLT=btest(c%lake(s),0)
      topo(j)%ISLAKE=btest(c%lake(s),1)
      topo(j)%LAKETARGVOL=btest(c%lake(s),2)
      topo(j)%LAKEMODELTYPE=c%lake_model(s)
      first=c%up_start(s); last=c%up_start(s+1)-1
      topo(j)%UREACHI=c%upstream(first:last)
      topo(j)%UREACHK=net%driver%seg_id(c%upstream(first:last))
      topo(j)%goodBas=c%good(first:last)/=0
      first=c%hru_start(s); last=c%hru_start(s+1)-1
      topo(j)%HRUIX=c%hru(first:last)
      topo(j)%HRUID=net%driver%hru_id(c%hru(first:last))
      topo(j)%HRUWGT=c%weight(first:last)
      topo(j)%REACHIX=j
      d=topo(j)%DREACHI
      topo(j)%DREACHI=0
      if(d>0 .and. d<=n) topo(j)%DREACHI=g2l(d)
      do k=1,size(topo(j)%UREACHI)
        u=topo(j)%UREACHI(k)
        if(u>0) topo(j)%UREACHI(k)=g2l(u)
      enddo
      do k=1,size(topo(j)%HRUIX)
        h=topo(j)%HRUIX(k)
        if(h>0) topo(j)%HRUIX(k)=h2l(h)
      enddo
      ! These all-upstream lists/UHs are unused by the supported KW/MC/DW methods.
      if(allocated(topo(j)%RCHLIST)) deallocate(topo(j)%RCHLIST)
      if(allocated(topo(j)%UH)) deallocate(topo(j)%UH)
      if(allocated(topo(j)%pfafCode)) deallocate(topo(j)%pfafCode)
    enddo
    call move_alloc(param,net%core%param)
    call move_alloc(topo,net%core%ntopo)
    net%driver%seg_id=net%driver%seg_id(reach)
    net%driver%hru_id=net%driver%hru_id(basin)
    net%driver%totArea=net%driver%totArea(reach)

    ! Retain only local basin remapping rows, preserving floating-point order.
    allocate(needed(size(p%host_owner)),host2local(size(p%host_owner)))
    needed=.false.
    if(do_remapping) then
      total=0
      do h=1,size(mapping%hru_ix)
        j=mapping%hru_ix(h)
        if(j<=0) cycle
        if(h2l(j)>0) total=total+mapping%num_qhru(h)
      enddo
      allocate(local_map%hru_id(nb),local_map%hru_ix(nb),local_map%num_qhru(nb), &
               local_map%qhru_id(total),local_map%qhru_ix(total),local_map%weight(total))
      pos=1; first=1; r=0
      do h=1,size(mapping%hru_ix)
        last=first+max(0,mapping%num_qhru(h))-1
        j=mapping%hru_ix(h)
        if(j>0) then
          if(h2l(j)>0) then
            r=r+1; total=last-first+1
            local_map%hru_id(r)=mapping%hru_id(h)
            local_map%hru_ix(r)=h2l(j)
            local_map%num_qhru(r)=total
            local_map%qhru_id(pos:pos+total-1)=mapping%qhru_id(first:last)
            local_map%qhru_ix(pos:pos+total-1)=mapping%qhru_ix(first:last)
            local_map%weight(pos:pos+total-1)=mapping%weight(first:last)
            needed(mapping%qhru_ix(first:last))=.true.
            pos=pos+total
          endif
        endif
        first=last+1
      enddo
      mapping=local_map
    else
      p%runoff_index=p%runoff_index(basin)
      needed(p%runoff_index)=.true.
    endif
    hosts=pack([(h,h=1,size(needed))],needed)
    host2local=0
    do h=1,size(hosts)
      host2local(hosts(h))=h
    enddo
    if(do_remapping) then
      mapping%qhru_ix=host2local(mapping%qhru_ix)
      deallocate(p%runoff_index)
    else
      p%runoff_index=host2local(p%runoff_index)
    endif
    net%core%runoff%hru_id=net%core%runoff%hru_id(hosts)

    ! Exchange requests once. Subsequent timesteps transfer just needed runoff.
    allocate(p%recv_counts(p%nproc),p%recv_offsets(p%nproc), &
             p%send_counts(p%nproc),p%send_offsets(p%nproc))
    do r=1,p%nproc
      p%recv_counts(r)=count(p%host_owner(hosts)==r-1)
    enddo
    call offsets_for(p%recv_counts,p%recv_offsets)
#ifdef ROUTING_MPI
    if(p%nproc>1) then
      call MPI_Alltoall(p%recv_counts,1,MPI_INTEGER,p%send_counts,1,MPI_INTEGER,p%comm,ierr)
      if(ierr/=0) return
    else
#endif
      p%send_counts=p%recv_counts
#ifdef ROUTING_MPI
    endif
#endif
    call offsets_for(p%send_counts,p%send_offsets)
    allocate(p%send_index(sum(p%send_counts)))
#ifdef ROUTING_MPI
    if(p%nproc>1) then
      call MPI_Alltoallv(hosts,p%recv_counts,p%recv_offsets,MPI_INTEGER, &
                        p%send_index,p%send_counts,p%send_offsets,MPI_INTEGER,p%comm,ierr)
      if(ierr/=0) return
    else
#endif
      p%send_index=hosts
#ifdef ROUTING_MPI
    endif
#endif
    p%send_index=p%send_index-p%offsets(p%rank+1)
    deallocate(p%owner,p%host_owner,p%counts,p%offsets)
    associate(t=>net%core%topology)
      t%n_seg=ns; t%n_hru=nb
      if(allocated(t%hru)) deallocate(t%hru)
      if(allocated(t%seg)) deallocate(t%seg)
      if(allocated(t%hru2seg)) deallocate(t%hru2seg)
      if(allocated(t%ntopo)) deallocate(t%ntopo)
      if(allocated(t%pfaf)) deallocate(t%pfaf)
    end associate
    end associate
  end subroutine

  subroutine offsets_for(counts, offsets)
    integer, intent(in) :: counts(:)
    integer, intent(out) :: offsets(:)
    integer :: r
    offsets(1)=0
    do r=2,size(counts)
      offsets(r)=offsets(r-1)+counts(r-1)
    enddo
  end subroutine

  subroutine setup_gather(comm, nproc, index, layout, ierr)
    integer, intent(in) :: comm, nproc
    integer, intent(in) :: index(:)
    type(routing_gather), intent(out) :: layout
    integer, intent(out) :: ierr
    integer :: n
    ierr=0; n=size(index)
    allocate(layout%counts(nproc),layout%offsets(nproc))
    layout%counts=0
#ifdef ROUTING_MPI
    if(nproc>1) then
      call MPI_Gather(n,1,MPI_INTEGER,layout%counts,1,MPI_INTEGER,0,comm,ierr)
      if(ierr/=0) return
    else
#endif
      layout%counts=n
#ifdef ROUTING_MPI
    endif
#endif
    call offsets_for(layout%counts,layout%offsets)
    allocate(layout%index(sum(layout%counts)))
#ifdef ROUTING_MPI
    if(nproc>1) then
      call MPI_Gatherv(index,n,MPI_INTEGER,layout%index,layout%counts,layout%offsets, &
                       MPI_INTEGER,0,comm,ierr)
    else
#endif
      layout%index=index
#ifdef ROUTING_MPI
    endif
#endif
  end subroutine

  ! One output timestep at a time: O(global domain) scratch only on rank zero.
  subroutine gather_output(p, layout, local, global, ierr, mpi_seconds)
    type(routing_parallel_data), intent(in) :: p
    type(routing_gather), intent(in) :: layout
    real(dp), intent(in) :: local(:)
    real(dp), allocatable, intent(out) :: global(:)
    integer, intent(out) :: ierr
    real(dp), allocatable :: packed(:)
    real(dp), intent(out) :: mpi_seconds
    real(dp) :: started
    ierr=0; mpi_seconds=0._dp
    allocate(packed(size(layout%index)),global(size(layout%index)))
#ifdef ROUTING_MPI
    if(p%nproc>1) then
      started=wall_seconds()
      call MPI_Gatherv(local,size(local),MPI_DOUBLE_PRECISION,packed,layout%counts,layout%offsets, &
                       MPI_DOUBLE_PRECISION,0,p%comm,ierr)
      mpi_seconds=wall_seconds()-started
      if(ierr/=0) return
    else
#endif
      packed=local
#ifdef ROUTING_MPI
    endif
#endif
    if(p%rank==0) global(layout%index)=packed
  end subroutine

  subroutine exchange_level(net, stage, method, ierr, message)
    type(river_network_data), intent(inout) :: net
    integer, intent(in) :: stage, method
    integer, intent(out) :: ierr
    character(*), intent(out) :: message
#ifdef ROUTING_MPI
    integer :: k, j, n
    real(dp) :: started
    integer, allocatable :: request(:), status(:,:)
    real(dp), allocatable, asynchronous :: values(:)
#endif
    ierr=0; message='exchange_level/'
#ifdef ROUTING_MPI
    associate(p=>net%driver%parallel)
    if(p%nproc==1) return
    n=count(p%edge_stage==stage)
    if(n==0) return
    allocate(request(n),status(MPI_STATUS_SIZE,n),values(n))
    j=0
    do k=1,size(p%edge_stage)
      if(p%edge_stage(k)/=stage) cycle
      j=j+1
      if(p%edge_send(k)) then
        values(j)=net%core%flux(p%edge_local(k))%ROUTE(method)%REACH_Q
        call MPI_Isend(values(j),1,MPI_DOUBLE_PRECISION,p%edge_peer(k),731,p%comm,request(j),ierr)
      else
        call MPI_Irecv(values(j),1,MPI_DOUBLE_PRECISION,p%edge_peer(k),731,p%comm,request(j),ierr)
      endif
      if(ierr/=0) return
    enddo
    started=wall_seconds()
    call MPI_Waitall(n,request,status,ierr)
    p%boundary_wait_seconds=p%boundary_wait_seconds+wall_seconds()-started
    if(ierr/=0) return
    j=0
    do k=1,size(p%edge_stage)
      if(p%edge_stage(k)/=stage) cycle
      j=j+1
      if(.not.p%edge_send(k)) net%core%flux(p%edge_local(k))%ROUTE(method)%REACH_Q=values(j)
    enddo
    end associate
#endif
  end subroutine
end module routing_parallel
