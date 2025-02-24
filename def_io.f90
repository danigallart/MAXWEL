module def_io
    
    character*(300) :: inputfile = 'EMWAVE.in'
    
    character*(300) :: result_scat = 'scat_field_results.csv'
    character*(300) :: result_tot = 'total_field_results.csv'
    character*(300) :: result_inc = 'inc_field_results.csv'
    
    character*(300) :: result_plane = 'inplane_field_results.csv'

    character*(300) :: control = 'tolerance_control.dat'
    
    character*(300) :: mesh_file = 'mesh_file.dat'
    
    character*(300) ::  stiff_matrix_file = 'stiff_matrix.txt'
    character*(300) ::  connectivity_file = 'connectivity.txt'
    
    character*(300) :: mesh_geometry        ! = 'data/EqMesh_TRIA3_FINE.geo.dat'
    character*(300) :: mesh_param           ! = 'data/EqMesh_TRIA3_FINE.dom.dat'
    character*(300) :: mesh_element         ! = 'data/EqMesh_TRIA3_FINE.set.dat'
    character*(300) :: mesh_boundary        ! = 'data/EqMesh_TRIA3_FINE.fix.dat'
    character*(300) :: mesh_flux            ! = 'data/EqMesh_TRIA3_FINE.psi.dat'
    character*(300) :: mesh_phys            ! = 'data/EqMesh_TRIA3_FINE.ker.dat'

    character*(300) :: logic_file = 'sparse_logic.dat'

    
    integer, parameter :: &
        input_unit = 1, &
        result_scat_unit = 3, &
        result_tot_unit = 4, &
        result_inc_unit = 5, &
        result_plane_unit = 16, &
        control_unit = 7, &
        mesh_unit = 8, &
        stiff_matrix_unit = 9, &
        connectivity_unit = 10, &
        mesh_geo_unit = 11, &
        mesh_param_unit = 12, &
        mesh_element_unit = 13, &
        mesh_boundary_unit = 14, &
        mesh_flux_unit = 15, &
        mesh_phys_unit = 17, &
        logic_unit = 18
    
end module def_io