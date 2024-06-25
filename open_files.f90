   subroutine open_files
    use def_io
    use def_variables
    implicit none
    
   open(unit=data_unit,file=output_data,status='unknown')
   open(unit=result_scat_unit,file=result_scat,status='unknown')
   open(unit=result_tot_unit,file=result_tot,status='unknown')
   open(unit=result_inc_unit,file=result_inc,status='unknown')
   open(unit=result_plane_unit,file=result_plane,status='unknown')
   open(unit=control_unit,file=control,status='unknown')
   
   open(unit=stiff_matrix_unit,file=stiff_matrix_file,status='unknown')
   open(unit=connectivity_unit,file=connectivity_file,status='unknown')

   if (reader_type == 'toka') then
      open(unit=mesh_geo_unit, file=mesh_geometry, status='old')
      open(unit=mesh_param_unit, file=mesh_param, status='old')
      open(unit=mesh_element_unit, file=mesh_element, status='old')
      open(unit=mesh_boundary_unit, file=mesh_boundary, status='old')
   endif
      
    end subroutine open_files