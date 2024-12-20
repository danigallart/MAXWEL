subroutine exit_nosolver
use def_io
use def_variables
use def_vectors
implicit none
! local
integer :: ii
character(1) :: tab=char(9)

write(stiff_matrix_unit,*) size(AD)

do ii=1,size(AD)
write(stiff_matrix_unit,*) AD(ii)
enddo

write(stiff_matrix_unit,*) size(AN)

do ii=1,size(AN)
write(stiff_matrix_unit,*) AN(ii)
enddo

write(stiff_matrix_unit,*) size(IA)

do ii=1,size(IA)
write(stiff_matrix_unit,*) IA(ii)
enddo

write(stiff_matrix_unit,*) size(JA)

do ii=1,size(JA)
write(stiff_matrix_unit,*) JA(ii)
enddo

write(stiff_matrix_unit,*) size(indep_vect)

do ii=1,size(indep_vect)
write(stiff_matrix_unit,*) indep_vect(ii)
enddo

do ii=1,NE
    if (nodpel == 3) then
        write(connectivity_unit,*) conn(ii,1),tab,conn(ii,2),tab,conn(ii,3)
    else if (nodpel == 4) then
        write(connectivity_unit,*) conn(ii,1),tab,conn(ii,2),tab,conn(ii,3),tab,conn(ii,4)
    else if (nodpel == 6) then
        write(connectivity_unit,*) conn(ii,1),tab,conn(ii,2),tab,conn(ii,3),tab,conn(ii,4),tab,conn(ii,5),tab,conn(ii,6)
    else if (nodpel == 8) then
        write(connectivity_unit,*) conn(ii,1),tab,conn(ii,2),tab,conn(ii,3),tab,conn(ii,4),tab,conn(ii,5),tab,conn(ii,6),tab,conn(ii,7),tab,conn(ii,8)
    end if
enddo
  
  
end subroutine exit_nosolver