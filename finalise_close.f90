subroutine finalise
use def_io
implicit none

    close(data_unit)
    close(result_scat_unit)
    close(result_tot_unit)
    close(result_inc_unit)
    close(control_unit)
    close(result_plane_unit)
    
    close(stiff_matrix_unit)
    close(connectivity_unit)
    
    close(mesh_geo_unit)
    close(mesh_param_unit)
    close(mesh_element_unit)
    close(mesh_boundary_unit)    
    close(mesh_flux_unit)   
    close(mesh_phys_unit)

end subroutine finalise