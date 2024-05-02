subroutine initialise
use def_io
implicit none

character*30 :: time_date


  call fdate(time_date) 

! mensaje de bienvenida
   print*, 'Welcome to Wave Popagation code.'
   print*, 'Time & date ',time_date
   print*, 'This program is protected by me, myself and I,'
   print*, 'for any kind of problems, bug report or suggestions'
   print*, 'contact my dear supervisor Alejandro Soba :)'
   
    
! apertura de archivo de entrada generico
      
   open(unit=input_unit,file=inputfile,status='old',err=100)
   open(unit=data_unit,file=output_data,status='unknown')
   open(unit=grid_unit,file=grid,status='unknown')
   open(unit=result_scat_unit,file=result_scat,status='unknown')
   open(unit=result_tot_unit,file=result_tot,status='unknown')
   open(unit=control_unit,file=control,status='unknown')
   
   open(unit=stiff_matrix_unit,file=stiff_matrix_file,status='unknown')
   open(unit=connectivity_unit,file=connectivity_file,status='unknown')

   
   open(unit=nodes_unit,file=nodes_file,status='old')
   open(unit=elements_unit, file=elements_file, status='old')
   
   open(unit=pml_bin_nodes_unit,file=pml_bin_nodes_file,status='old')
   open(unit=pml_bout_nodes_unit, file=pml_bout_nodes_file, status='old')
   open(unit=pml_nodes_unit,file=pml_nodes_file,status='old')
   
   open(unit=scatb_nodes_unit, file=scatb_nodes_file, status='old')
   open(unit=scatin_nodes_unit,file=scatin_nodes_file,status='old')
   open(unit=scatin_elements_unit, file=scatin_elements_file, status='old')
   open(unit=scatb_elements_unit,file=scatb_elements_file,status='old')
   
   open(unit=huygb_nodes_unit, file=huygb_nodes_file, status='old')
   open(unit=huygb_elements_unit,file=huygb_elements_file,status='old')
   

   ! 


    return

100 write(6,*) 'error in' // inputfile // 'opening'
	print*,"The program will stop"
	stop ' '
   
   

end subroutine initialise
