subroutine exit_writer
use def_io
use def_variables
use def_vectors
implicit none
! local
integer :: jj
character(1) :: coma=','


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (pol == 'TM') then

  write(result_scat_unit,*) 'X,Y,Escat_z'
  
else if (pol == 'TE') then
    
    write(result_scat_unit,*) 'X,Y,Hscat_z'
    
endif 
  

do jj=1,NP
    
    write(result_scat_unit,'(E15.5,a,E15.5,a,E15.5)') coorx(jj),coma,coory(jj),coma,abs(u_scat(jj))
    
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (pol == 'TM') then

    write(result_tot_unit,*) 'X,Y,Etot_z'
    
else if (pol == 'TE') then
    
    write(result_tot_unit,*) 'X,Y,Htot_z'

endif

  
do jj=1,NP
    
    write(result_tot_unit,'(E15.5,a,E15.5,a,E15.5)') coorx(jj),coma,coory(jj),coma,abs(u_tot(jj))
    
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (pol == 'TM') then

    write(result_inc_unit,*) 'X,Y,Etot_z'
    
else if (pol == 'TE') then
    
    write(result_inc_unit,*) 'X,Y,Htot_z'

endif

  
do jj=1,NP
    
    write(result_inc_unit,'(E15.5,a,E15.5,a,E15.5)') coorx(jj),coma,coory(jj),coma,real(u_inc(jj))
    
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (pol == 'TM') then

    write(result_plane_unit,*) 'X,Y,Hre_x,Hre_y,Him_x,Him_y'
    
else if (pol == 'TE') then
    
    write(result_plane_unit,*) 'X,Y,Ere_x,Ere_y,Eim_x,Eim_y'

endif

  
do jj=1,NE
    
    write(result_plane_unit,'(E15.5,a,E15.5,a,E15.5,a,E15.5,a,E15.5,a,E15.5)') coorx_mid(jj),coma,coory_mid(jj),coma,real(plane_field_x(jj)),coma,real(plane_field_y(jj)),coma,imag(plane_field_x(jj)),coma,imag(plane_field_y(jj))
    
enddo
  
end subroutine exit_writer