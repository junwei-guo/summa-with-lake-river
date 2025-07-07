!This module contains subroutines for summa to implement MPI parallelization.
module summa_mpi

USE MPI
public::mpi_print                       !print message to console in mpi enviroment. 
public::mpi_sum_bcast
public::num2str
public::mpi_initial_load_balancer
public::mpi_dynamic_load_balancer
public::mpi_load_balancer_from_table
public::mpi_print_gru_distribution
public::mpi_write_outputs




integer::sGRU_COMM_WORLD, nGRU_COMM_WORLD  !Start and number of GRUs from the command line input (common world of this program excution).
integer::sGRU_this_rank,  nGRU_this_rank
integer::num_rank, idx_rank, mpi_err
!real::mpiSyncFreq = -1                 !MPI barrier synchronization frequency, unit day. if value is missing or negative, no synchronization

logical             :: mpi_table_err            ! if the mpi_table.txt exist.
character(len=256)  :: mpi_table_filename
logical             :: mpi_table_loaded         ! if the mpi_table.txt is correctly loaded.



contains 

subroutine mpi_print(msg,rank_to_prnt)
    ! This subroutine will print message to terminal only when the MPI rank is rank_to_prnt. 
    ! If rank_to_prnt is < 0, print messages from all ranks.
    INTEGER, OPTIONAL,VALUE         :: rank_to_prnt
    !integer                         :: num_rank, idx_rank, mpi_err
    character(len=*), intent(in)    :: msg

    IF(.NOT. PRESENT(rank_to_prnt)) rank_to_prnt = -1 ! default value

    !call MPI_Comm_size(MPI_COMM_WORLD, num_rank, mpi_err)
    !call MPI_Comm_rank(MPI_COMM_WORLD, idx_rank, mpi_err)
    if (num_rank==1) then
        rank_to_prnt = 0
    end if 
    if (rank_to_prnt>=0) then 
        !Print message at rank_to_prnt only.
        if (idx_rank==rank_to_prnt) then 
            print *, msg
            !write(*, '(A)') msg
        end if
    else 
        !print *, "[MPI message from (idx_rank, num_rank) = (",idx_rank,", ",num_rank,")]: ",msg
        print *,  "[MSG from rank@("//trim(num2str(idx_rank))//", "//trim(num2str(num_rank))//")]: "//msg
    end if 
end

subroutine mpi_even_load_balancer()
    !Assgin the GRUs to each MPI rank in a even manner.
    integer:: nGRU_per_rank,rankTran

    nGRU_per_rank  = int(ceiling(real(nGRU_COMM_WORLD)/real(num_rank))) !Round numbers up
    rankTran       = nGRU_COMM_WORLD-(nGRU_per_rank-1)*num_rank

    if (idx_rank<rankTran) then
        sGRU_this_rank  = sGRU_COMM_WORLD+idx_rank*nGRU_per_rank
        nGRU_this_rank  = nGRU_per_rank
    else 
        sGRU_this_rank  = sGRU_COMM_WORLD+ rankTran*nGRU_per_rank +(idx_rank-rankTran)*(nGRU_per_rank-1)
        nGRU_this_rank  = nGRU_per_rank-1
    endif

end 

subroutine mpi_load_balancer_from_table()

    integer    :: mpi_table_values(4)
    integer    :: mpi_table_line_idx
    integer    :: mpi_table_values_idx
    integer    :: mpi_table_nGRU
    

    ! GRU load balancing using 'mpi_table.txt'
    if (idx_rank==0) print *,  "[MPI manager] "//trim(mpi_table_filename)//" exist. Start reading the MPI load balancing using the table values."
    
    ! Open the file and read sGRU and nGRU values from the table.
    open(unit=10, file=mpi_table_filename, status='old', action='read')    
    do mpi_table_line_idx = 0, idx_rank 
        read(10, *) (mpi_table_values(mpi_table_values_idx), mpi_table_values_idx = 1, 4)
    end do

    sGRU_this_rank  = mpi_table_values(2)  ! start GRU at this MPI rank
    nGRU_this_rank  = mpi_table_values(3)  ! number of GRUs at this MPI rank
    close(10)
    mpi_table_loaded = .TRUE.

    !------ Check if this load balancing is correct:

    mpi_table_nGRU = mpi_sum_bcast(nGRU_this_rank) !total number of GRUs read from the mpi_table
    if (mpi_table_nGRU /= nGRU_COMM_WORLD ) then 
      if (idx_rank==0) print *, "[MPI manager] nGRU from the table does not equal to fileGRU or nGRU from command line. Assgin GRUs to each MPI rank evenly."
      if (idx_rank==0) print *, "[MPI manager] mpi_table_nGRU ="//trim(num2str(mpi_table_nGRU))//", nGRU_COMM_WORLD = "//trim(num2str(nGRU_COMM_WORLD))
      mpi_table_loaded = .FALSE.
    end if
    call MPI_Barrier(MPI_COMM_WORLD, mpi_err)

end

subroutine mpi_dynamic_load_balancer()
    !Re-assgin the GRUs to each MPI rank during the time stepping

end

subroutine mpi_print_gru_distribution(sGRU,nGRU)
    !Print distribution of the GRUs on MPI ranks.
    integer:: sGRU, nGRU 

    !call MPI_Barrier(MPI_COMM_WORLD, mpi_err)
    if (idx_rank==0) print *,"[MPI manager] (sGRU_COMM_WORLD, nGRU_COMM_WORLD) = ("//trim(num2str(sGRU_COMM_WORLD))//", "//trim(num2str(nGRU_COMM_WORLD))//"). " 
    call MPI_Barrier(MPI_COMM_WORLD, mpi_err)

    print *, "[MPI manager] (sGRU_this_rank, nGRU_this_rank)   = (", sGRU,", ",nGRU,") @ rank #",idx_rank,"."
   
    if (nGRU<1) then 
     print *, "[MPI manager] no GRU found @ rank #",idx_rank,"."
    end if 
end

subroutine mpi_write_outputs()
    !Unifying the netcdf outputs

end

function mpi_sum_bcast(locVar) result(var_sum)
    integer :: locVar, var_sum
    integer :: ierr 
    call MPI_Reduce(locVar, var_sum, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(var_sum, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
end  function mpi_sum_bcast


function num2str(number) result(temp_str)
    implicit none
    integer, intent(in) :: number
    character(len=256) :: temp_str,str
    ! Write the number to a string

    write(temp_str , *) number
    
end function num2str



function dou2str(number) result(temp_str)
    implicit none
    real, intent(in) :: number
    character(len=256) :: temp_str,str
    ! Write the number to a string

    write(temp_str , *) number
    
end function dou2str

end module summa_mpi 
