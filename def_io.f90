module def_io
    
    character*(300) :: inputfile = 'EMWAVE.in'
    character*(300) :: output_data = 'EMWAVE.out'
    
    character*(300) :: result_scat = 'scat_field_results.csv'
    character*(300) :: result_tot = 'total_field_results.csv'
    character*(300) :: result_inc = 'inc_field_results.csv'
    
    character*(300) :: result_plane = 'inplane_field_results.csv'

    character*(300) :: control = 'tolerance_control.dat'
    
    character*(300) :: nodes_file = 'circle_mesh_nodes.dat'
    character*(300) :: elements_file = 'circle_mesh_elements.dat'
    
    character*(300) ::  stiff_matrix_file = 'stiff_matrix.txt'
    character*(300) ::  connectivity_file = 'connectivity.txt'


    
    integer, parameter :: &
        input_unit = 1, &
        data_unit = 2, &
        result_scat_unit = 3, &
        result_tot_unit = 4, &
        result_inc_unit = 5, &
        result_plane_unit = 6, &
        control_unit = 7, &
        nodes_unit = 8, &
        elements_unit = 9, &
        stiff_matrix_unit = 10, &
        connectivity_unit = 11
    
end module def_io