subroutine initialise
use def_io
implicit none

character*30 :: time_date


  call fdate(time_date) 

! mensaje de bienvenida
   print*, 'Welcome to MAXWEL code.'
   print*, 'Time & date ',time_date
   print*, 'For any kind of problems, bug report or suggestions'
   print*, 'contact hdomingo@bsc.es'
   
    
! apertura de archivo de entrada generico
      
   open(unit=input_unit,file=inputfile,status='old',err=100)
   open(unit=data_unit,file=output_data,status='unknown')
   open(unit=result_scat_unit,file=result_scat,status='unknown')
   open(unit=result_tot_unit,file=result_tot,status='unknown')
   open(unit=result_inc_unit,file=result_inc,status='unknown')
   open(unit=result_plane_unit,file=result_plane,status='unknown')
   open(unit=control_unit,file=control,status='unknown')
   
   open(unit=stiff_matrix_unit,file=stiff_matrix_file,status='unknown')
   open(unit=connectivity_unit,file=connectivity_file,status='unknown')

   open(unit=nodes_unit,file=nodes_file,status='old')
   open(unit=elements_unit, file=elements_file, status='old')
   
   open(unit=mesh_geo_unit, file=mesh_geometry, status='old')
   open(unit=mesh_param_unit, file=mesh_param, status='old')
   open(unit=mesh_element_unit, file=mesh_element, status='old')
   open(unit=mesh_boundary_unit, file=mesh_boundary, status='old')
   
    return

100 write(6,*) 'error in' // inputfile // 'opening'
	print*,"The program will stop"
	stop ' '
   
   

end subroutine initialise
