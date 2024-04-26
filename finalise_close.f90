subroutine finalise
use def_io
implicit none

    close(input_unit)
    close(data_unit)
    close(grid_unit)
    close(result_scat_unit)
    close(result_tot_unit)
    close(control_unit)
    close(nodes_unit)
    close(elements_unit)
    close(pml_bin_nodes_unit)
    close(pml_bout_nodes_unit)
    close(pml_nodes_unit)
    close(scatb_nodes_unit)
    close(scatin_nodes_unit)
    close(scatin_elements_unit)
    close(scatb_elements_unit)
    close(huygb_nodes_unit)
    close(huygb_elements_unit)
    
    close(stiff_matrix_unit)
    close(connectivity_unit)
    

end subroutine finalise