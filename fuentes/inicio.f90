subroutine inicio()
use def_variables
implicit none

    open(unit=unit_data,file=filedata)
    open(unit=unit_malla,file=filemalla)
    open(unit=unit_cc,file=filecc)
    open(unit=unit_cont,file=file_aux)
    open(unit=unit_temp,file=file_temp)
    open(unit=unit_burn,file=file_burn)
    open(unit=unit_area,file=file_area,Access = 'append')

!    open(unit=unit_sal,file=file_sal)
!    open(unit=unit_2d,file=file_2D)
!    open(unit=unit_gra,file=file_gra)
!    open(unit=unit_grid,file=file_grid)
!    open(unit=unit_camp,file=file_camp)
!    open(unit=unit_cel,file=file_cel)
!    open(unit=unit_bio,file=file_bio)
!    open(unit=unit_bio2,file=file_bio2)
    
end subroutine inicio