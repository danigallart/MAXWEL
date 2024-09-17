subroutine exit_writer
use def_io
use def_variables
use def_vectors
implicit none
! local
integer :: jj
character(1) :: coma=','
double precision :: x_val,y_val
double precision :: x_real,x_imag
double precision :: y_real,y_imag
double precision :: angle_x,angle_y


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (pol == 'TM') then

  write(result_scat_unit,*) 'X,Y,Escat_real_z,Escat_imag_z'
  
else if (pol == 'TE') then
    
    write(result_scat_unit,*) 'X,Y,Hscat_real_z,Hscat_imag_z'
    
endif 
  

do jj=1,NP
    
    write(result_scat_unit,'(E15.5,a,E15.5,a,E15.5,a,E15.5)') coorx(jj),coma,coory(jj),coma,real(u_scat(jj)),coma,imag(u_scat(jj))
    
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (pol == 'TM') then

    write(result_tot_unit,*) 'X,Y,Etot_real_z,Etot_imag_z'
    
else if (pol == 'TE') then
    
    write(result_tot_unit,*) 'X,Y,Htot_real_z,Htot_imag_z'

endif

  
do jj=1,NP
    
    write(result_tot_unit,'(E15.5,a,E15.5,a,E15.5,a,E15.5)') coorx(jj),coma,coory(jj),coma,real(u_tot(jj)),coma,imag(u_tot(jj))
    
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (pol == 'TM') then

    write(result_inc_unit,*) 'X,Y,Einc_z'
    
else if (pol == 'TE') then
    
    write(result_inc_unit,*) 'X,Y,Hinc_z'

endif

  
do jj=1,NP
    
    write(result_inc_unit,'(E15.5,a,E15.5,a,E15.5)') coorx(jj),coma,coory(jj),coma,real(u_inc(jj))
    
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (pol == 'TM') then

    write(result_plane_unit,*) 'X,Y,Hreal_x,Himag_x,Hreal_y,Himag_y'
    
else if (pol == 'TE') then
    
    write(result_plane_unit,*) 'X,Y,Ereal_x,Eimag_x,Ereal_y,Eimag_y'

endif

  
do jj=1,NE
    
    x_real = real(plane_field_x(jj))
    x_imag = imag(plane_field_x(jj))
    y_real = real(plane_field_y(jj))
    y_imag = imag(plane_field_y(jj))
    
    write(result_plane_unit,'(E15.5,a,E15.5,a,E15.5,a,E15.5,a,E15.5,a,E15.5)') coorx_mid(jj),coma,coory_mid(jj),coma,x_real,coma,x_imag,coma,y_real,coma,y_imag
    
enddo
  
end subroutine exit_writer