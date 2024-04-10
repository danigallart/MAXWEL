module def_variables

character(10):: problema
character(10):: nombre
character(20):: archi_malla='    '
character(20):: archi_tumor='    '
character(20):: archi_sistema='    '

double precision,allocatable :: coor_x(:),coor_y(:),coor_z(:),grad_x(:),grad_y(:),grad_z(:),cer(:),gradxel_x(:), &
     &                          gradxel_y(:),gradxel_z(:),tempera(:),tempera_ant(:),jcurrent(:)

integer, allocatable :: conect(:,:), material(:),vec_tierra(:),vec_poten(:),nnodtierra(:),nnodpoten(:),electro_elem(:,:),electro_nod(:,:)
integer, allocatable :: vec_temper(:),nnodtemper(:),vecino(:)

integer, allocatable :: quemado(:)
integer :: cuento_quemado

double precision,allocatable ::  Zgrid(:)
integer, allocatable :: grid(:,:),grid_z(:,:,:),grilla2d(:)
integer :: nzgrid,Ngrid_z,Ngrid_y,Ngrid_x

integer :: nelements,nnodes,nod_tierra,nod_poten,mat_externo,mat_agua,mat_normal,mat_tumor,mat_aislante,mat_electro,  &
   &       elec_elem,nod_centro,nod_temper

integer :: nopcion
integer :: nop_malla


integer :: nnodes_tum
double precision,allocatable :: coor_x_tum(:),coor_y_tum(:),coor_z_tum(:)

integer :: ncerx0,ncery0,ncerz0
integer :: ndimx,ndimy,ndimz

double precision :: tam_zonex,tam_zoney,tam_zonez
double precision :: tamano_elem= 0.1  ! unidades de entrada

integer :: tipo_elec,elec_alt=1,elec_sep=1
double precision :: rad_zone1,rad_zone2,sigma1,sigma2,sigma3,sigma4,rad1,rad2,rad3,rad4,alt1,alt2,alt3
double precision :: potencial, tierra,electropor,electro_no

double precision :: sigma_ave,area_trans,temper_max,sigma_max,E_maximo

double precision :: tiempo_total,pas_ti2,pas_ti1,tiempo_sumado
double precision,allocatable :: historiat(:),histo_gen(:),histo_pot(:)
integer :: pasos_tiempo = 0


! variables apra el probelma Electro
integer :: ne_lado=20 ,    &  ! numero de elementos por lado
   &       nop_el=4           ! nodos por elemento

double precision:: ladox=1.0
double precision:: ladoy=1.0


double precision:: factorx=1.0
double precision:: factory=1.0
double precision:: desplax=0.0
double precision:: desplay=0.0
double precision:: pos_tierra=0.0
double precision:: pos_pot=0.0
double precision:: ancho_inf=0.0
double precision:: ancho_sup=0.0


! variables apra el probelma HELMOZT
complex*16, allocatable  :: cplx_grad_x(:),cplx_grad_y(:),cplx_gradxel_x(:),cplx_gradxel_y(:)
complex*16  :: cplx_ko = 1 ! free space wave number
complex*16  ::cplx_eo=1 ! permittivity free space
complex*16  ::cplx_er=1  !relative permittivity
complex*16  ::cplx_mur=1  !relative permeability

complex*16  ::cplx_sigma=1 ! conductivity
complex*16  ::cplx_w=1   ! ?
complex*16  ::cplx_erc(2,2) ! tensores
complex*16  ::cplx_murc(2,2)



! variables apra el probelma Bergues
integer :: N_electro=2,    &  ! Numero de electrodos en el dominio
   &       ne_ber=20 ,    &  ! numero de elementos grandes
   &       nop_b=16   ,    &  ! nodos por elemento
   &       modo=2             ! modo de arreglo de polaridad electrdos 1=intercalado, 2 enfrentados 

! variables apra el probelma Barbara
integer :: N_electro_tot    ! Numero de electrodos totales
integer :: Base_electro=1    ! Numero de electrodos totales
double precision ::   Alterna        ! numero de repetiiones


   
double precision ::    dis_ele=   1.0,    & ! distancia entre electrodos (mm)
     &                 exentric=  1.0,    & ! 0.0 esfera,  0.6 elipse, 1.0 parabola,  2.0 hyperbola
     &                 zonatot=   40.0,   & ! tamaño de la zona a resolver
     &                 radio_a,           & ! radio y a
     &                 radio_b,           & ! b
     &                 ang30 =0.5236,     &
     &                 ang60=1.0472,      &   
     &                 dprot= 5.26e-3,     &  !1.6E-3,     & ! Difu coef de protones mm2/seg
     &                 dOh= 9.31e-4,      &   !0.8E-3,     & ! Difu coef de oxhidrilos mm2/seg
     &                 font=1.0,          & !valor de la fuentede produccion de H
     &                 sumir=1.0            ! valor de consumo de Ph  

double precision ::  zh = 1.0               !	          Carga de los protones (e) */
double precision ::  zoh = -1.0              !	          Carga de los OH (-e) */
double precision :: Faraday= 96485.34  !C/mol
double precision :: R_cte= 8.314      ! J/K/mol
double precision :: T_cte= 350.0      ! K                 

double precision :: mup=0.00036           ! 36.23E-2  ! Mov ionica proton (mm2 / V s) 
double precision :: muoh= 0.000058  ! Mov ionica OH (mm2 / V s) 


integer,allocatable :: fuentes(:),sumideros(:)
double precision,allocatable :: xfuentes(:),xsumideros(:),yfuentes(:),ysumideros(:),conducta(:)
double precision,allocatable :: Ph_(:)


integer :: unit_data=1,  &! data de entrada
   &       unit_malla=2,   & ! malla de salida
   &       unit_cc=3,   &  ! condiciones de contorno    
   &       unit_cont=10,   & ! para controlas
   &       unit_sal=11,   &  ! resultados
   &       unit_2d=12,   &  
   &       unit_gra=13,   & 
   &       unit_grid=14,   & 
   &       unit_camp=15,  &   ! camp elemental
   &       unit_cel=16,   &  ! para ver celulas tumorales 
   &       unit_bio=17,   &     ! biothemal evolution
   &        unit_bio2=18, &       
   &       unit_ph=19,    &       ! ph
   &       unit_II=20,   &        ! corriente por el dominio 
   &       unit_temp=21,   &     ! temperatura
   &       unit_oh=22,   &    ! OH
   &       unit_burn=23,   &    ! quemado 
   &       unit_area=100           ! areas

character*(30) :: filedata='input.in',  & ! entrada
     &          filemalla='mallado.fem' , &  ! malla
     &          filecc = 'contorno.fem',  &    ! CC
     &          file_aux='control.dat',  & ! controla
     &          file_sal='results.csv',   &   ! salidas de resultados
     &          file_2D='saleplano.csv',  & ! alguna salida especial     
     &          file_gra='gradiente.csv', & ! gradientes     
     &          file_grid='grid2D.dat',  &   ! alguna salida especial     
     &          file_camp='campo.csv',  &     ! campo elemental     
     &          file_cel='cancer.dat',    &
     &          file_bio='biothermal.csv', &
     &          file_bio2='biomini.csv', &
     &          file_ph='Ph.csv',         &
     &          file_II='Current.csv',    &
     &          file_temp='tempera.his',    &
     &          file_area='area.dat',    &
     &          file_oh='OH.csv',               &
     &          file_burn='quemado.csv'
          
     
end module def_variables