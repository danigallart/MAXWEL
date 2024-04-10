subroutine finalice()
use def_variables
implicit none

    close(unit_malla)
    close(unit_cc)
    close(unit_cont)
    close(unit_temp)
    close(unit_burn)
!    close(unit_sal)
!    close(unit_2d)
!    close(unit_gra)
!    close(unit_grid)
!    close(unit_camp)
!    close(unit_cel)
!    close(unit_bio)
!    close(unit_bio2)

end subroutine finalice