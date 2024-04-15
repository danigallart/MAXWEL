module def_constantes


    integer :: ne=10,    &!30 numero de elemtnos por lado
       &       nlog=0,    & ! si queiro malla logaritmica  
       &       nodpel=27, & ! numero de nodos por elemento   PAPA, CECI
!       &       nodpel=8, & ! numero de nodos por elemento   BOCHUM
       &       nope=3,    & ! nodos por direccion   PAPA, CECI
!       &       nope=2,    & ! nodos por direccion    BOCHUM
       &       DUPLICY=0, &  ! si duplica en y
       &       DUPLICZ=0,  & ! si duplica en z
       &       OrdenAlya=1   ! si ordena como en alya!! PAPA, CECI
!       &       OrdenAlya=0   ! si ordena como en alya!! BOCHUM

   integer :: caso_temporal=1  ! 0esta apgado el tiempo. 1 para la papa en 2d
   
      
   double precision :: sig_electrodo= 2200.0
   double precision :: pendiente=1.9   ! 1.7-- 1000   1.9--1500/1700   1.4--500  1.6--800  
   double precision :: toler= 1e-9 ! tolerancia solver
   integer :: itermax=10000
   
   integer :: HayHYALU=0  ! 0 - 1
   
   double precision :: field_limite= 10.0   ! 12.0 c/hyalu   ! 14 s/hyalu
   double precision :: Area_electroporada=0.0
   
   double precision :: radio=0.15, &  ! cm
      &                XX_alto=3.8    ! cm


   integer :: ncada_Xnod=4
   integer :: nmallo_Xnod=4

   integer :: ncerx=20,   &
   &       ncery=20,   &
   &       ncerz=20

   integer :: OPTED=0
   integer :: AGUJAS=5 ! 2 - 4 - 6 - 8

    integer :: NE0
 ! papa  
   integer :: poselex=5
   integer :: poseley=10
   integer :: poselez=10

! Nahuel
!   integer :: poselex=1
!   integer :: poseley=7
!   integer :: poselez=5

! bergues3d
   double precision :: alfazone
   double precision :: pos_cir1=0.328
   double precision :: pos_cir2=0.656

   !integer :: ncada_Xnod=16 
   !integer :: nmallo_Xnod=4

   
   double precision :: escala=0.05   ! escalado log


   double precision :: limite_Basico = 9.0
   double precision :: limite_Acido = 5.5
   double precision :: AreaDOMINIO

   double precision :: temper40=40.0

    end module def_constantes
    
    !New branch mod test.