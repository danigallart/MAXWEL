!subroutine read_previous_solution()
!    
!    USE def_io
!    USE def_variables
!    USE def_vectors
!    
!    implicit none
!    
!    integer :: i,io
!    double precision :: x,y,sol_re,sol_im
!    character(len=120) :: line
!    character(len=1) :: comma
!    
!    allocate(prev_sol_array(NP))
!    
!    !open(unit=result_tot_unit,file=result_tot,status='old',action='read')
!    
!    read(read_tot_unit,'(A)') line
!    
!    do i=1,NP
!        read(read_tot_unit,FMT=*,iostat=io) x,y,sol_re,sol_im
!        prev_sol_array(i)%re=sol_re
!        prev_sol_array(i)%im=sol_im
!    enddo
!    
!    
!    !close(result_tot_unit)
!    
!    
!    end subroutine read_previous_solution