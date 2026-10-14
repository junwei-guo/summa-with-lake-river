! Initialize routing without forcing/physics to measure startup memory separately.
! Usage: mpirun -n N probe_routing_init /absolute/setup/directory/
program probe_routing_init
  use mpi
  use netcdf
  use nrtype, only: i8b, dp
  use mizuroute_types, only: mizuroute_info, mizuroute_domain
  use init_mizuRoute, only: init_mizuroute_domain
  use routing_parallel, only: collect_runoff_ids
  implicit none
  type(mizuroute_info) :: info
  type(mizuroute_domain) :: domain
  integer(i8b), allocatable :: file_ids(:), host_ids(:)
  integer :: rank, nproc, ierr, ncid, vid, did, n, expected_reaches, first, count, u, ios
  integer :: local_counts(2), global_counts(2), rss_kib, peak_kib, sum_peak_kib
  character(4096) :: setup, message, line
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
  call get_command_argument(1,setup)
  if(len_trim(setup)==0) call MPI_Abort(MPI_COMM_WORLD,1,ierr)
  call check(nf90_open(trim(setup)//'/settings/SUMMA/attributes.nc',nf90_nowrite,ncid))
  call check(nf90_inq_dimid(ncid,'gru',did))
  call check(nf90_inquire_dimension(ncid,did,len=n))
  allocate(file_ids(n))
  call check(nf90_inq_varid(ncid,'gruId',vid))
  call check(nf90_get_var(ncid,vid,file_ids))
  call check(nf90_close(ncid))
  call check(nf90_open(trim(setup)//'/settings/mizuRoute/topology.nc',nf90_nowrite,ncid))
  call check(nf90_inq_dimid(ncid,'seg',did))
  call check(nf90_inquire_dimension(ncid,did,len=expected_reaches))
  call check(nf90_close(ncid))
  count=n/nproc+merge(1,0,rank<mod(n,nproc))
  first=rank*(n/nproc)+min(rank,mod(n,nproc))+1
  domain%river_network%driver%parallel%rank=rank
  domain%river_network%driver%parallel%nproc=nproc
  domain%river_network%driver%parallel%comm=MPI_COMM_WORLD
  call collect_runoff_ids(domain%river_network%driver%parallel,file_ids(first:first+count-1), &
                          host_ids,ierr,message)
  call check(ierr)
  deallocate(file_ids)
  info%do_mizuroute=.true.
  info%do_remapping=.false.
  info%is_print=rank==0
  info%is_gridded=.false.
  info%dt_landmodel=3600._dp
  info%mrout%namelist_path=trim(setup)//'/settings/mizuRoute/'
  info%mrout%namelist_file='mizuroute.param'
  info%mrout%methods='3'
  info%mrout%dt=3600._dp
  info%ntopo%hfabric_path=trim(setup)//'/settings/mizuRoute/'
  info%ntopo%hfabric_file='topology.nc'
  info%ntopo%dname_seg='seg'
  info%ntopo%dname_hru='hru'
  info%ntopo%varname_HRUid='hruId'
  info%ntopo%varname_segId='segId'
  info%ntopo%varname_hruSegId='hruToSegId'
  info%ntopo%varname_downSegId='downSegId'
  info%ntopo%varname_area='area'
  info%ntopo%varname_slope='slope'
  info%ntopo%varname_length='length'
  call init_mizuroute_domain(info,domain,[n,-9999],1,host_ids,1._dp,1._dp,ierr,message)
  if(ierr/=0) then
    write(*,*) rank,trim(message)
    call MPI_Abort(MPI_COMM_WORLD,ierr,ios)
  endif
  associate(net=>domain%river_network)
    if(allocated(net%core%topology%seg)) call MPI_Abort(MPI_COMM_WORLD,2,ios)
    local_counts=[size(net%core%state),size(net%core%runoff%basinRunoff)]
    if(size(net%core%param)/=local_counts(1)) call MPI_Abort(MPI_COMM_WORLD,3,ios)
    call MPI_Reduce(local_counts,global_counts,2,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  end associate
  peak_kib=0; rss_kib=0
  open(newunit=u,file='/proc/self/status',status='old')
  do
    read(u,'(a)',iostat=ios) line
    if(ios/=0) exit
    if(index(line,'VmRSS:')==1) read(line(7:),*) rss_kib
    if(index(line,'VmHWM:')==1) read(line(7:),*) peak_kib
  enddo
  close(u)
  write(*,'(a,i0,a,i0,a,i0)') 'Initialization memory rank ',rank,': RSS_KiB=',rss_kib,', peak_KiB=',peak_kib
  call MPI_Reduce(peak_kib,sum_peak_kib,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  if(rank==0) then
    if(any(global_counts/=[expected_reaches,n])) call MPI_Abort(MPI_COMM_WORLD,4,ios)
    write(*,'(a,2(i0,1x),a,f10.3)') 'PASS initialized reaches/HRUs: ',global_counts, &
      '; sum of per-rank peak GiB: ',real(sum_peak_kib,dp)/1048576._dp
  endif
  call MPI_Finalize(ierr)
contains
  subroutine check(status)
    integer, intent(in) :: status
    if(status/=0) then
      write(*,*) rank,' initialization probe error ',status
      call MPI_Abort(MPI_COMM_WORLD,status,ierr)
    endif
  end subroutine
end program
