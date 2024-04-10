module def_solver


  double precision, allocatable :: ad(:),an(:),rhs(:),solucion(:),solucion_ant(:),rhs_ant(:)
  integer,allocatable :: cx(:),ia(:),ja(:)
  
  double precision, allocatable :: masa(:)
  
  complex*16, allocatable :: cplx_ad(:),cplx_an(:),cplx_rhs(:),cplx_solucion(:),cplx_solucion_ant(:),cplx_rhs_ant(:)
  
  integer :: NONULL
    
end module def_solver