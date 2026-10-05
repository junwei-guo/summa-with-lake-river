! Distributed reach execution using the host SUMMA GRU ownership.
! Topology/ghost storage is replicated; only the owning rank advances a reach.
module routing_parallel
  use nrtype, only: i4b, i8b, dp
  use mizuroute_types, only: routing_parallel_data, river_network_data
  use dataTypes, only: remap
  use var_lookup, only: ixNTOPO
  use nr_utils, only: indexx, match_index
#ifdef ROUTING_MPI
  use mpi
#endif
  implicit none
  private
  public :: collect_runoff_ids, collect_runoff, plan_routing, exchange_level, collect_streamflow
contains
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
    type(routing_parallel_data), intent(in) :: p
    real(dp), intent(in) :: local_runoff(:)
    real(dp), intent(out) :: runoff(:)
    integer, intent(out) :: ierr
    character(*), intent(out) :: message
    ierr=0; message='collect_runoff/'
#ifdef ROUTING_MPI
    if(p%nproc>1) then
      call MPI_Allgatherv(local_runoff,size(local_runoff),MPI_DOUBLE_PRECISION,runoff, &
                         p%counts,p%offsets,MPI_DOUBLE_PRECISION,p%comm,ierr)
      return
    endif
#endif
    runoff=local_runoff
  end subroutine

  subroutine plan_routing(net, mapping, do_remapping, ierr, message)
    type(river_network_data), intent(inout) :: net
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
      do k=1,size(net%core%ntopo(s)%HRUIX)
        h=net%core%ntopo(s)%HRUIX(k)
        if(h<=0) cycle
        if(do_remapping) then
          do j=map_start(h),map_start(h)+map_count(h)-1
            r=p%host_owner(mapping%qhru_ix(j))+1
            share(r)=share(r)+net%core%ntopo(s)%HRUWGT(k)*mapping%weight(j)
          enddo
        else
          r=p%host_owner(p%runoff_index(h))+1
          share(r)=share(r)+net%core%ntopo(s)%HRUWGT(k)
        endif
      enddo
      if(sum(share)>0._dp) p%owner(s)=maxloc(share,dim=1)-1
    enddo
    ! Basinless connectors inherit a downstream owner, then an upstream owner.
    do k=n,1,-1
      s=net%core%topology%ntopo(k)%var(ixNTOPO%rchOrder)%dat(1)
      d=net%core%ntopo(s)%DREACHI
      if(p%owner(s)<0 .and. d>0 .and. d<=n) p%owner(s)=p%owner(d)
    enddo
    level=0
    do k=1,n
      s=net%core%topology%ntopo(k)%var(ixNTOPO%rchOrder)%dat(1)
      level(s)=1
      do j=1,size(net%core%ntopo(s)%UREACHI)
        u=net%core%ntopo(s)%UREACHI(j)
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
      s=net%core%topology%ntopo(k)%var(ixNTOPO%rchOrder)%dat(1)
      l=level(s); pos=p%level_start(l)+used(l)
      p%level_order(pos)=s; used(l)=used(l)+1
    enddo
    if(p%nproc>1) write(*,'(a,i0,a,i0,a,i0)') 'Routing rank ',p%rank,': owned reaches=', &
                                           count(p%owner==p%rank),', stages=',nlevels
    end associate
  end subroutine

  subroutine exchange_level(net, stage, method, ierr, message)
    type(river_network_data), intent(inout) :: net
    integer, intent(in) :: stage, method
    integer, intent(out) :: ierr
    character(*), intent(out) :: message
    integer :: k, s, d
    ierr=0; message='exchange_level/'
#ifdef ROUTING_MPI
    associate(p=>net%driver%parallel)
    if(p%nproc==1) return
    do k=p%level_start(stage),p%level_start(stage+1)-1
      s=p%level_order(k)
      d=net%core%ntopo(s)%DREACHI
      if(d<=0 .or. d>size(p%owner)) cycle
      if(p%owner(s)==p%owner(d)) cycle
      ! Current discharge is the boundary condition for the supported KW/MC/DW methods.
      call MPI_Bcast(net%core%flux(s)%ROUTE(method)%REACH_Q,1,MPI_DOUBLE_PRECISION, &
                     p%owner(s),p%comm,ierr)
      if(ierr/=0) return
    enddo
    end associate
#endif
  end subroutine

  subroutine collect_streamflow(net,method,sub_idx,ierr,message)
    type(river_network_data), intent(inout) :: net
    integer, intent(in) :: method,sub_idx
    integer, intent(out) :: ierr
    character(*), intent(out) :: message
    ierr=0; message='collect_streamflow/'
#ifdef ROUTING_MPI
    if(net%driver%parallel%nproc>1) then
      ! Exactly one rank contributes each reach; other entries are zero.
      call MPI_Allreduce(MPI_IN_PLACE,net%driver%method(method)%streamflow(:,sub_idx), &
                         net%core%topology%n_seg,MPI_DOUBLE_PRECISION,MPI_SUM, &
                         net%driver%parallel%comm,ierr)
    endif
#endif
  end subroutine
end module routing_parallel
