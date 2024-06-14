module def_io
    
    character*(300) :: inputfile = 'EMWAVE.in'
    character*(300) :: output_data = 'EMWAVE.out'
    
    character*(300) :: result_scat = 'scat_field_results.csv'
    character*(300) :: result_tot = 'total_field_results.csv'
    character*(300) :: result_inc = 'inc_field_results.csv'
    
    character*(300) :: result_plane = 'inplane_field_results.csv'

    character*(300) :: control = 'tolerance_control.dat'
    
    character*(300) :: mesh_file = 'mesh_file.dat'
    
    character*(300) ::  stiff_matrix_file = 'stiff_matrix.txt'
    character*(300) ::  connectivity_file = 'connectivity.txt'
    
    character*(300) :: mesh_geometry = 'data/EquilibiumMesh.geo.dat'
    character*(300) :: mesh_param = 'data/EquilibiumMesh.dom.dat'
    character*(300) :: mesh_element = 'data/EquilibiumMesh.set.dat'
    character*(300) :: mesh_boundary = 'data/EquilibiumMesh.fix.dat'

    character*(300) :: logic_file = 'data/sparse_logic.dat'

    
    integer, parameter :: &
        input_unit = 1, &
        data_unit = 2, &
        result_scat_unit = 3, &
        result_tot_unit = 4, &
        result_inc_unit = 5, &
        result_plane_unit = 6, &
        control_unit = 7, &
        mesh_unit = 8, &
        stiff_matrix_unit = 9, &
        connectivity_unit = 10, &
        mesh_geo_unit = 11, &
        mesh_param_unit = 12, &
        mesh_element_unit = 13, &
        mesh_boundary_unit = 14, &
        logic_unit = 15
    
end module def_io