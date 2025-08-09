!This module contains the interface to update the river network routing.

module summa_routing 


    ! Structure to control routing coupling
    type :: routing_coupling_control
        logical :: stand_alone_routing      ! One-way coupling
        logical :: sync_routing 
    contains 
        procedure :: set_stand_alone_routing
        procedure :: set_sync_routing
    end type routing_coupling_control


    ! Pass Summa setting to mizuRoute, avoid mizuRoute reads duplicated input from control file or Summa NC files. 
    type :: setting_pass_through
        real :: summa_time_step         ! in sec.
        real :: summa_start_time
        real :: summa_end_time

    contains 

    end type setting_pass_through

    contains 

        subroutine set_stand_alone_routing(this)
          class(routing_coupling_control):: this

            this%stand_alone_routing = .true.
            this%sync_routing = .false. 

        end subroutine set_stand_alone_routing


        subroutine set_sync_routing(this)
          class(routing_coupling_control):: this

            this%stand_alone_routing = .false.
            this%sync_routing = .true. 

        end subroutine set_sync_routing




end module summa_routing