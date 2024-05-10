subroutine finalise
use def_io
implicit none

    close(input_unit)
    close(data_unit)
    close(result_scat_unit)
    close(result_tot_unit)
    close(control_unit)
    close(nodes_unit)
    close(elements_unit)
    
    close(stiff_matrix_unit)
    close(connectivity_unit)
    

end subroutine finalise