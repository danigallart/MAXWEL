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

    write(result_plane_unit,*) 'X,Y,Htot_x,Htot_y'
    
else if (pol == 'TE') then
    
    write(result_plane_unit,*) 'X,Y,Etot_x,Etot_y'

endif

  
do jj=1,NE
    
    x_real = real(plane_field_x(jj))
    x_imag = imag(plane_field_x(jj))
    y_real = real(plane_field_y(jj))
    y_imag = imag(plane_field_y(jj))
    
    angle_x = atan2(x_imag,x_real)
    angle_y = atan2(y_imag,y_real)
    
    
    if ((angle_x >= 0.0) .and. (angle_x < pi)) then
        
        x_val = abs(plane_field_x(jj))
        
    else
        
        x_val = -abs(plane_field_x(jj))
        
    endif
    
    if ((angle_y >= 0.0) .and. (angle_y < pi)) then
        
        y_val = abs(plane_field_y(jj))
        
    else
        
        y_val = -abs(plane_field_y(jj))
        
    endif
    
    
    
    write(result_plane_unit,'(E15.5,a,E15.5,a,E15.5,a,E15.5)') coorx_mid(jj),coma,coory_mid(jj),coma,x_val,coma,y_val
    
enddo
  
end subroutine exit_writer