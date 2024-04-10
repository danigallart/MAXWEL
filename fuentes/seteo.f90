subroutine seteo()
use def_bio
use def_constantes
use def_variables
implicit none


if(problema=='BERGUES') then

  
  Tblod = 37.0
  Tout  = 37.0 
  h_b   =0.010
  W_b = 7.1e-3 
  C_b = 3.8400
  rho_b = 0.0010600
  k_c=0.015
  rho_c = 0.0079
  Cp_c = 0.5
  alfa_b=0.032        

  
  caso_temporal=0
    
  nope=2
  OrdenAlya=0   


elseif(problema=='PAPA') then

  Tblod = 25.0
  Tout  = 25.0
  h_b   =0.010
  pulso = 1.0E-2
  frec = 1
  W_b   =0.0
  C_b   = 3.8400
  rho_b = 0.001060
  k_c   = 0.000562
  rho_c = 0.0011
  Cp_c  = 3.78
  alfa_b=0.135 ! 0.032 ! 0.135 delr referi puto

  ne = 15
  nlog  =0
  nodpel=27
  nope  =3
  DUPLICY=0
  DUPLICZ=0
  OrdenAlya=1   

  
elseif(problema=='PAPA2D') then

  Tblod = 25.0 
  Tout  = 25.0 
  h_b   =0.010
  W_b = 0.0  
  C_b = 3.8400
  rho_b = 0.0010600
  k_c=0.000562
  rho_c = 0.0011
  
  Cp_c = 3.78

  alfa_b=0.032        

  pulso = 1.0E-2
  frec = 1
  
  ne_ber=80 


  caso_temporal=1  ! prendido el tiempo (PAPA!!!!)

  nope=2
  OrdenAlya=0   

elseif(problema=='RATONES') then

  Tblod = 37.0
  Tout  = 37.0 
  h_b   =0.010
  W_b = 7.1e-3 
  C_b = 3.8400
  rho_b = 0.0010600
  k_c=0.015
  rho_c = 0.0079
  Cp_c = 0.5
  alfa_b=0.032        

  
  
  nope=3
  OrdenAlya=1   
  
   
  pulso = 1.0E-2
  frec = 1
  
  caso_temporal=1  

  
endif



end subroutine seteo