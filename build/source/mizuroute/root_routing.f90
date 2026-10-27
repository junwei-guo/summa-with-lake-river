! Root-owned routing: workers send GRU runoff, but never allocate a river model.
module root_routing
  use nrtype, only: dp, i8b
  use mizuroute_types, only: routing_parallel_data
  use nr_utils, only: indexx
#ifdef ROUTING_MPI
  use mpi
#endif
  implicit none
  private
  public :: collect_root_ids, collect_root_runoff, routing_wall_seconds
contains
  subroutine collect_root_ids(p, local_ids, ids, ierr, message)
    type(routing_parallel_data), intent(inout) :: p
    integer(i8b), intent(in) :: local_ids(:)
    integer(i8b), allocatable, intent(out) :: ids(:)
    integer, intent(out) :: ierr
    character(*), intent(out) :: message
    integer, allocatable :: sorted(:)
    integer :: n, r, j, mpierr
    ierr=0; message='collect_root_ids/'
    n=size(local_ids)
    allocate(p%counts(p%nproc),p%offsets(p%nproc))
    p%counts=0; p%offsets=0
#ifdef ROUTING_MPI
    if(p%nproc>1) then
      call MPI_Gather(n,1,MPI_INTEGER,p%counts,1,MPI_INTEGER,0,p%comm,ierr)
      if(ierr/=0) return
    else
#endif
      p%counts(1)=n
#ifdef ROUTING_MPI
    endif
#endif
    do r=2,p%nproc
      p%offsets(r)=p%offsets(r-1)+p%counts(r-1)
    enddo
    allocate(ids(sum(p%counts)))
#ifdef ROUTING_MPI
    if(p%nproc>1) then
      call MPI_Gatherv(local_ids,n,MPI_INTEGER8,ids,p%counts,p%offsets,MPI_INTEGER8,0,p%comm,ierr)
      if(ierr/=0) return
    else
#endif
      ids=local_ids
#ifdef ROUTING_MPI
    endif
#endif
    if(p%rank==0) then
      allocate(sorted(size(ids)))
      call indexx(ids,sorted)
      do j=2,size(ids)
        if(ids(sorted(j-1))==ids(sorted(j))) then
          ierr=20; message=trim(message)//'duplicate SUMMA GRU ID'; exit
        endif
      enddo
    endif
#ifdef ROUTING_MPI
    if(p%nproc>1) then
      call MPI_Bcast(ierr,1,MPI_INTEGER,0,p%comm,mpierr)
      if(mpierr/=0) then
        ierr=mpierr; return
      endif
      if(ierr/=0) call MPI_Bcast(message,len(message),MPI_CHARACTER,0,p%comm,mpierr)
    endif
#endif
  end subroutine

  subroutine collect_root_runoff(p, local_runoff, runoff, ierr)
    type(routing_parallel_data), intent(inout) :: p
    real(dp), intent(in) :: local_runoff(:)
    real(dp), intent(out) :: runoff(:)
    integer, intent(out) :: ierr
    real(dp) :: started
    ierr=0
#ifdef ROUTING_MPI
    if(p%nproc>1) then
      started=routing_wall_seconds()
      call MPI_Gatherv(local_runoff,size(local_runoff),MPI_DOUBLE_PRECISION, &
                       runoff,p%counts,p%offsets,MPI_DOUBLE_PRECISION,0,p%comm,ierr)
      p%runoff_mpi_seconds=p%runoff_mpi_seconds+routing_wall_seconds()-started
      return
    endif
#endif
    runoff=local_runoff
  end subroutine

  function routing_wall_seconds() result(seconds)
    real(dp) :: seconds
    integer(i8b) :: count, rate
    call system_clock(count,rate)
    seconds=real(count,dp)/real(rate,dp)
  end function
end module
