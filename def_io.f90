module def_io
    
    character*(300) :: inputfile = 'EMWAVE.inp'
    character*(300) :: output_data = 'EMWAVE.out'
    character*(300) :: grid = 'grid_data.dat'
    character*(300) :: result_scat = 'scat_field_results.csv'
    character*(300) :: result_tot = 'total_field_results.csv'

    character*(300) :: control = 'tolerance_control.dat'
    
    character*(300) :: nodes_file = 'test_circle_nodes.txt'
    character*(300) :: elements_file = 'test_circle_elements.txt'

    character*(300) ::  pml_bin_nodes_file = 'test_circle_pmlbin_nodes.txt'
    character*(300) ::  pml_bout_nodes_file = 'test_circle_pmlbout_nodes.txt'
    character*(300) ::  pml_nodes_file = 'test_circle_pml_nodes.txt'


    character*(300) ::  scatb_nodes_file = 'test_circle_scatb_nodes.txt'
    character*(300) ::  scatin_nodes_file = 'test_circle_scatin_nodes.txt'
    character*(300) ::  scatin_elements_file = 'test_circle_scatin_elements.txt'
    character*(300) ::  scatb_elements_file = 'test_circle_scatb_elements.txt'


    character*(300) ::  huygb_nodes_file = 'test_circle_huygb_nodes.txt'
    character*(300) ::  huygb_elements_file = 'test_circle_huygb_elements.txt'
    
    character*(300) ::  stiff_matrix_file = 'stiff_matrix.txt'
    character*(300) ::  connectivity_file = 'connectivity.txt'


    
    integer, parameter :: &
        input_unit = 1, &
        data_unit = 2, &
        grid_unit = 3, &
        result_scat_unit = 4, &
        result_tot_unit = 5, &
        control_unit = 7, &
        nodes_unit = 8, &
        elements_unit = 9, &
        pml_bin_nodes_unit = 10, &
        pml_bout_nodes_unit = 11, &
        pml_nodes_unit = 12, &
        scatb_nodes_unit = 13, &
        scatin_nodes_unit = 14, &
        scatin_elements_unit = 15, &
        scatb_elements_unit = 16, &
        huygb_nodes_unit = 17, &
        huygb_elements_unit = 18, &
        stiff_matrix_unit = 19, &
        connectivity_unit = 20
    
end module def_io