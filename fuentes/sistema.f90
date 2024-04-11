subroutine sistema()
use def_variables
use def_constantes
use def_solver
implicit none

if(nopcion==1) then     
  write(6,*) '      voy a control!'
  
  if(problema=='HELMOZT') then
     call control_cplx()
  else
     call control()
  endif
  
else
   write(6,*) '      leo el sistema!'
   call lee_sistema(archi_sistema,nnodes)
endif

if(problema=='PAPA') then 

    if(potencial==500) then
       pendiente =  1.4
    elseif(potencial == 800)then
       pendiente =  1.6
    elseif(potencial == 1000)then
       pendiente =  1.7
    elseif(potencial == 1500)then
       pendiente =  1.85
    elseif(potencial == 1700)then
       pendiente =  1.9
    endif   

    ! seteo en uno para volver a empezar con el tema
    !pendiente = 1.0 

endif

    end subroutine sistema

    
    !! Test comment