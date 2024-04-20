subroutine exit_writer
use def_io
use def_variables
use def_vectors
implicit none
! local
integer :: jj
character(1) :: coma=','

  write(result_scat_unit,*) 'X, Y, Escat_z(V/m) '
  
  do jj=1,NP
      write(result_scat_unit,'(E15.5,a,E15.5,a,E15.5)') coorx(jj),coma,coory(jj),coma,abs(u_scat(jj))
  enddo
  
    write(result_tot_unit,*) 'X, Y, Etot_z(V/m) '
  
  do jj=1,NP
      write(result_tot_unit,'(E15.5,a,E15.5,a,E15.5)') coorx(jj),coma,coory(jj),coma,abs(u_scat(jj)+u_inc(jj))
  enddo
  
  
end subroutine exit_writer