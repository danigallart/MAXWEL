subroutine leemalla()
use def_variables
use def_constantes
use def_bio
implicit none
! local
integer :: kk,kk1,kk2,kk3,jj,ii,j,n,ns(8),ntum_old,ntum,nelet
double precision :: cont,c,xmed,ymed,zmed,xliminf,xlimsup,yliminf,ylimsup,zliminf,zlimsup,posx,posy,posz,cer_med,z,xpaso
integer :: busco_arriba,busco_adelante,busco_atras,busco_derecha,busco_izquierda,nelements_new,nnodes_new
integer,allocatable :: nel_tum(:,:),mataux(:),conect_new(:,:),material_new(:),material_build(:)
double precision,allocatable :: aux_x(:),aux_y(:),aux_z(:),cooy_new(:),coox_new(:),cooz_new(:)
integer :: mat_extra,total,nlimite,ncont,nmat4,nxpos,nypos,mat_potencial,mat_tierra,je
integer,allocatable :: nod_electro_pos(:),nod_electro_neg(:)



open(unit=1111,file=archi_malla)
open(unit=1112,file=archi_tumor)


read(1111,*) ncerx0,ncery0,ncerz0

read(1111,*) nnodes

allocate(coor_x(nnodes),coor_y(nnodes),coor_z(nnodes),cer(nnodes))
do kk=1,nnodes
    read(1111,*) n,coor_x(kk),coor_y(kk),coor_z(kk),cer(kk)
enddo



read(1111,*) nelements
nodpel=8
allocate(conect(nelements,nodpel),material(nelements),material_build(nelements))

material_build=0

do kk=1,nelements
       read(1111,*) n,(conect(kk,jj),jj=1,nodpel)
enddo          

close(1111)


read(1112,*) nnodes_tum

allocate(coor_x_tum(nnodes_tum),coor_y_tum(nnodes_tum),coor_z_tum(nnodes_tum))

do kk=1,nnodes_tum
    read(1112,*) n,coor_x_tum(kk),coor_y_tum(kk),coor_z_tum(kk),c
enddo


close(1112)

mat_agua=0
mat_normal=0
mat_externo=0

! determino el material de cada elemento
do kk=1,nelements
   
   cer_med=0
   do jj=1,nodpel
      j=conect(kk,jj)
      cer_med=cer_med+cer(j)
   enddo 
   cer_med=cer_med/real(nodpel)
   
   if(cer_med<=15) then !! exterior
      mat_externo=mat_externo+1
      material(kk) = 0
   elseif(cer_med>15 .and. cer_med<= 75) then ! liquido
      mat_agua=mat_agua+1
      material(kk) = 1
   elseif(cer_med>75 .and. cer_med<= 110) then ! mate gris
      mat_normal=mat_normal+1
      material(kk) = 2
   elseif(cer_med>110 ) then ! mate blanca
      mat_normal=mat_normal+1
      material(kk) = 3
   else
      material(kk)=123456789
   endif

enddo



posx=(ncerx0-1)*real(ncada_Xnod)
posy=(ncery0-1)*real(ncada_Xnod)
posz=(ncerz0-1)*real(ncada_Xnod)

!xliminf = (ncerx0)*real(ncada_Xnod)-ncerx-1
!xlimsup = (ncerx0)*real(ncada_Xnod)+ncerx-1
!yliminf = (ncery0)*real(ncada_Xnod)-ncery-1
!ylimsup = (ncery0)*real(ncada_Xnod)+ncery-1
!zliminf = (ncerz0)*real(ncada_Xnod)-ncerz-1
!zlimsup = (ncerz0)*real(ncada_Xnod)+ncerz-1

xliminf = (ncerx0)-ncerx-mod(((ncerx0)-ncerx),2) 
xlimsup = (ncerx0)+ncerx-1
yliminf = (ncery0)-ncery-1
ylimsup = (ncery0)+ncery-1
zliminf = (ncerz0)-ncerz-1
zlimsup = (ncerz0)+ncerz-1

xpaso = real(ncada_Xnod)/real(nmallo_Xnod)

Ngrid_z = (zlimsup-zliminf)/xpaso + 1
Ngrid_y = (ylimsup-yliminf)/xpaso + 1
Ngrid_x = (xlimsup-xliminf)/xpaso + 1


allocate(grid_z(Ngrid_z,Ngrid_y,Ngrid_x),Zgrid(Ngrid_z))

do kk=1,Ngrid_z
   Zgrid(kk) = zliminf + (kk-1)*xpaso
enddo



ntum=0
! determino ubicacion del mallado fino
do kk=1,nelements
  xmed=0.0
  ymed=0.0
  zmed=0.0
  do jj=1,nodpel
     xmed=xmed + coor_x(conect(kk,jj))
     ymed=ymed + coor_y(conect(kk,jj))
     zmed=zmed + coor_z(conect(kk,jj))
  enddo
  xmed=xmed/real(nodpel) 
  ymed=ymed/real(nodpel) 
  zmed=zmed/real(nodpel) 

  if(xmed>=xliminf*0.99999 .and. xmed<=xlimsup*1.000001 )     then
       
       if(ymed>=yliminf*0.99999 .and. ymed<=ylimsup*1.000001 )     then
          
          if(zmed>=zliminf*0.99999 .and. zmed<=zlimsup*1.000001 )     then
               
               ntum=ntum+1

               material_build(kk)=10

          endif
       endif
  endif

enddo





nlimite=nmallo_Xnod*nmallo_Xnod*nmallo_Xnod*ntum

allocate(nel_tum(nlimite,8))
allocate(aux_x(nlimite*8),aux_y(nlimite*8),aux_z(nlimite*8),mataux(nlimite))

!expando los elementos de la zona fina
ntum_old=ntum

ntum=0
nelet=0
do kk=1,nelements
  if(material_build(kk)==10) then
     do jj=1,nodpel
        ns(jj)=conect(kk,jj)
     enddo
    call expando(ntum,kk,ns,nel_tum,nlimite,nelet,aux_x,aux_y,aux_z,mataux)
  endif
enddo




! rearmo la malla
nelements_new = nelements + nelet-ntum_old
nnodes_new = nnodes + ntum

allocate(conect_new(nelements_new,8),coox_new(nnodes_new),cooy_new(nnodes_new),cooz_new(nnodes_new),material_new(nelements_new))


ncont=0
do kk=1,nelements

   if(material_build(kk)/=10) then

    do jj=1,nodpel
        j=conect(kk,jj)
        conect_new(kk,jj)=j
    enddo
    material_new(kk)=material(kk)   
   else
    ncont=ncont+1
    do jj=1,nodpel
        j=nel_tum(ncont,jj)
        conect_new(kk,jj)=j
    enddo

    material_new(kk)=mataux(ncont)   
     
   
   endif


enddo



do kk=ncont+1,nelet
    do jj=1,nodpel
        j=nel_tum(kk,jj)
        conect_new(nelements+kk-ncont,jj)=j
    enddo
    material_new(nelements+kk-ncont)=mataux(kk)   
enddo


do j=1,nnodes 
   coox_new(j)=coor_x(j)
   cooy_new(j)=coor_y(j)
   cooz_new(j)=coor_z(j)
enddo


do kk=nnodes+1,nnodes_new 
    coox_new(kk)=aux_x(kk-nnodes)
    cooy_new(kk)=aux_y(kk-nnodes)
    cooz_new(kk)=aux_z(kk-nnodes)

enddo


deallocate(conect,material,coor_x,coor_y,coor_z)
nelements=nelements_new
nnodes=nnodes_new

allocate(conect(nelements,nodpel),material(nelements),coor_x(nnodes),coor_y(nnodes),coor_z(nnodes))


do  kk=1,nelements
   do jj=1,nodpel
       conect(kk,jj)=conect_new(kk,jj)
   enddo
   material(kk)=material_new(kk)
enddo



do  kk=1,nnodes

   coor_x(kk)=coox_new(kk)
   coor_y(kk)=cooy_new(kk)
   coor_z(kk)=cooz_new(kk)
   do jj=1,Ngrid_z
      if(coor_z(kk)==Zgrid(jj)) then
         
         if(coor_y(kk)>=yliminf .and.coor_y(kk)<=ylimsup ) then

             nypos = int8( coor_y(kk)-yliminf ) + 1 
                       
             if(coor_x(kk)>=xliminf .and.coor_x(kk)<=xlimsup ) then
                  
                  nxpos = int8( coor_x(kk)-xliminf )   + 1
                  grid_z(jj,nypos,nxpos)=kk
 !                 write(unit_cc,*) jj,nypos,nxpos, grid_z(jj,nypos,nxpos)
             endif
         endif      
      endif
   enddo
   
enddo

if(Ngrid_x/=nxpos) then
  
   write(6,*) '! warning!! troubles in Grid_XX ',Ngrid_x,nxpos  

endif



deallocate(coox_new,cooy_new,cooz_new,conect_new,material_new,nel_tum,aux_x,aux_y,aux_z,mataux,material_build)


            
! ahora limito el tumor
xliminf=coor_x_tum(1)
xlimsup=coor_x_tum(1)
yliminf=coor_y_tum(1)
ylimsup=coor_y_tum(1)
zliminf=coor_z_tum(1)
zlimsup=coor_z_tum(1)

do kk=2,nnodes_tum
    if(coor_x_tum(kk)<xliminf) xliminf=coor_x_tum(kk)
    if(coor_y_tum(kk)<yliminf) yliminf=coor_y_tum(kk)
    if(coor_z_tum(kk)<zliminf) zliminf=coor_z_tum(kk)


    if(coor_x_tum(kk)>xlimsup) xlimsup=coor_x_tum(kk)
    if(coor_y_tum(kk)>ylimsup) ylimsup=coor_y_tum(kk)
    if(coor_z_tum(kk)>zlimsup) zlimsup=coor_z_tum(kk)

enddo

nmat4=0
do kk=1,nelements
  xmed=0.0
  ymed=0.0
  zmed=0.0
  do jj=1,nodpel
     xmed=xmed + coor_x(conect(kk,jj))
     ymed=ymed + coor_y(conect(kk,jj))
     zmed=zmed + coor_z(conect(kk,jj))
  enddo
  xmed=xmed/real(nodpel) 
  ymed=ymed/real(nodpel) 
  zmed=zmed/real(nodpel) 

  if(xmed>=xliminf*0.99999 .and. xmed<=xlimsup*1.000001 )     then
       
       if(ymed>=yliminf*0.99999 .and. ymed<=ylimsup*1.000001 )     then
          
          if(zmed>=zliminf*0.99999 .and. zmed<=zlimsup*1.000001 )     then
               
               nmat4=nmat4+1
               material(kk)=4
          endif
       endif
  endif

enddo

! ahora determino elementos aislante y nodos del electrodo

allocate(nod_electro_pos(AGUJAS/2),nod_electro_neg(AGUJAS/2))


do kk=1,nnodes
   if(coor_z(kk)==ncerz0) then
      if(coor_y(kk)==ncery0) then
          if(coor_x(kk)==ncerx0) then
             nod_centro = kk ! nodo en donde comienza el tumor
          endif
      endif
   endif


   if(AGUJAS==2) then
    
     if(coor_z(kk)==ncerz0+poselez) then
     
       if(coor_y(kk)==ncery0) then
       
         if(coor_x(kk)==ncerx0-1.0) then
           nod_electro_pos(1) = kk ! nodo en donde comienza el tumor
         elseif(coor_x(kk)==ncerx0+2.0) then
           nod_electro_neg(1) = kk ! nodo en donde comienza el tumor
         endif
       
       endif

     endif
  elseif(AGUJAS==4) then

     if(coor_z(kk)==ncerz0) then
       if(coor_x(kk)==ncerx0) then
         if(coor_y(kk)==ncery0-10.0) then
           nod_electro_pos(1) = kk ! nodo en donde comienza el tumor
         elseif(coor_y(kk)==ncery0+10.0) then
           nod_electro_neg(1) = kk ! nodo en donde comienza el tumor
         endif
       endif
    

       if(coor_y(kk)==ncery0) then
         if(coor_x(kk)==ncerx0-10.0) then
           nod_electro_pos(2) = kk 
         elseif(coor_x(kk)==ncerx0+10.0) then
           nod_electro_neg(2) = kk ! nodo en donde comienza el tumor
         endif
       endif

     endif
  elseif(AGUJAS==6) then
     
    
    if(coor_z(kk)==ncerz0+poselez) then
       
       if(coor_y(kk)==ncery0) then
         if(coor_x(kk)==ncerx0-poselex*0.5) then
           nod_electro_pos(1) = kk 
         elseif(coor_x(kk)==ncerx0+poselex*0.5) then
           nod_electro_neg(1) = kk 
         endif
       endif
       
       if(coor_x(kk)==ncerx0-poselex*0.5) then
         if(coor_y(kk)==ncery0-poseley) then
           nod_electro_neg(2) = kk 
         elseif(coor_y(kk)==ncery0+poseley) then
           nod_electro_neg(3) = kk 
         endif
       endif


       if(coor_x(kk)==ncerx0+poselex*0.5) then
         if(coor_y(kk)==ncery0-poseley) then
           nod_electro_pos(2) = kk 
         elseif(coor_y(kk)==ncery0+poseley) then
           nod_electro_pos(3) = kk 
         endif
       endif
     
     
     
     endif

   
  endif



enddo

if(OPTED==1) then

   allocate(electro_elem((elec_sep+2)*elec_alt*3,8))

! busco elemento con ese nodo en el corner izquierdo-infrior-trasero
   elec_elem=0
   do kk=1,nelements
  
    if(conect(kk,1)==nod_centro) then
      elec_elem=elec_elem+1
      do jj=1,nodpel
        electro_elem(elec_elem,jj)=conect(kk,jj)
        material(kk)=5
      enddo
    endif
   enddo

! ahora busco elemento atras, 
  kk = busco_atras(nelements,nodpel,conect,electro_elem(elec_elem,:)) ! otro aislante
  elec_elem=elec_elem+1
  do jj=1,nodpel
    electro_elem(elec_elem,jj)=conect(kk,jj)
    material(kk)=5
  enddo

  kk = busco_atras(nelements,nodpel,conect,electro_elem(elec_elem,:))  ! electroro pos
  elec_elem=elec_elem+1
  do jj=1,nodpel
    electro_elem(elec_elem,jj)=conect(kk,jj)
    material(kk)=10
  enddo

  kk = busco_adelante(nelements,nodpel,conect,electro_elem(elec_elem-2,:))  ! aislante
  elec_elem=elec_elem+1
  do jj=1,nodpel
    electro_elem(elec_elem,jj)=conect(kk,jj)
    material(kk)=5
  enddo
  
  kk = busco_adelante(nelements,nodpel,conect,electro_elem(elec_elem,:))  ! electrodo neg
  elec_elem=elec_elem+1
  do jj=1,nodpel
    electro_elem(elec_elem,jj)=conect(kk,jj)
    material(kk)=11
  enddo

! electrodos
  kk = busco_arriba(nelements,nodpel,conect,electro_elem(5,:))  ! electrodo neg
  elec_elem=elec_elem+1
  do jj=1,nodpel
    electro_elem(elec_elem,jj)=conect(kk,jj)
    material(kk)=11
  enddo
  
  do ii=2,elec_alt-1
    kk = busco_arriba(nelements,nodpel,conect,electro_elem(elec_elem,:))  ! electrodo neg
    elec_elem=elec_elem+1
    do jj=1,nodpel
      electro_elem(elec_elem,jj)=conect(kk,jj)
      material(kk)=11
    enddo
  enddo
  
  kk = busco_arriba(nelements,nodpel,conect,electro_elem(3,:))  ! electrodo pos
  elec_elem=elec_elem+1
  do jj=1,nodpel
    electro_elem(elec_elem,jj)=conect(kk,jj)
    material(kk)=10
  enddo
  do ii=2,elec_alt-1
    kk = busco_arriba(nelements,nodpel,conect,electro_elem(elec_elem,:))  ! electrodo pos
    elec_elem=elec_elem+1
    do jj=1,nodpel
      electro_elem(elec_elem,jj)=conect(kk,jj)
      material(kk)=10
    enddo
  enddo

! aislante
  kk = busco_arriba(nelements,nodpel,conect,electro_elem(1,:))  
  elec_elem=elec_elem+1
  do jj=1,nodpel
    electro_elem(elec_elem,jj)=conect(kk,jj)
    material(kk)=5
  enddo

    do ii=2,elec_alt-1
      kk = busco_arriba(nelements,nodpel,conect,electro_elem(elec_elem,:))  
      elec_elem=elec_elem+1
      do jj=1,nodpel
        electro_elem(elec_elem,jj)=conect(kk,jj)
        material(kk)=5
      enddo
    enddo


    kk = busco_derecha(nelements,nodpel,conect,electro_elem(1,:))  
    elec_elem=elec_elem+1
    do jj=1,nodpel
      electro_elem(elec_elem,jj)=conect(kk,jj)
      material(kk)=5
    enddo

    do ii=2,elec_alt-1
      kk = busco_arriba(nelements,nodpel,conect,electro_elem(elec_elem,:))  
      elec_elem=elec_elem+1
      do jj=1,nodpel
        electro_elem(elec_elem,jj)=conect(kk,jj)
        material(kk)=5
      enddo
    enddo


    kk = busco_izquierda(nelements,nodpel,conect,electro_elem(1,:))  
    elec_elem=elec_elem+1
    do jj=1,nodpel
      electro_elem(elec_elem,jj)=conect(kk,jj)
      material(kk)=5
    enddo

    do ii=2,elec_alt-1
      kk = busco_arriba(nelements,nodpel,conect,electro_elem(elec_elem,:))  
      elec_elem=elec_elem+1
      do jj=1,nodpel
        electro_elem(elec_elem,jj)=conect(kk,jj)
        material(kk)=5
      enddo
    enddo


    kk = busco_arriba(nelements,nodpel,conect,electro_elem(2,:))  
    elec_elem=elec_elem+1
    do jj=1,nodpel
      electro_elem(elec_elem,jj)=conect(kk,jj)
      material(kk)=5
    enddo

    do ii=2,elec_alt-1
      kk = busco_arriba(nelements,nodpel,conect,electro_elem(elec_elem,:))  
      elec_elem=elec_elem+1
      do jj=1,nodpel
        electro_elem(elec_elem,jj)=conect(kk,jj)
        material(kk)=5
      enddo
    enddo

    kk = busco_derecha(nelements,nodpel,conect,electro_elem(2,:))  
    elec_elem=elec_elem+1
    do jj=1,nodpel
      electro_elem(elec_elem,jj)=conect(kk,jj)
      material(kk)=5
    enddo

    do ii=2,elec_alt-1
      kk = busco_arriba(nelements,nodpel,conect,electro_elem(elec_elem,:))  
      elec_elem=elec_elem+1
      do jj=1,nodpel
        electro_elem(elec_elem,jj)=conect(kk,jj)
        material(kk)=5
      enddo
    enddo

    kk = busco_izquierda(nelements,nodpel,conect,electro_elem(2,:))  
    elec_elem=elec_elem+1
    do jj=1,nodpel
      electro_elem(elec_elem,jj)=conect(kk,jj)
      material(kk)=5
    enddo

    do ii=2,elec_alt-1
      kk = busco_arriba(nelements,nodpel,conect,electro_elem(elec_elem,:))  
      elec_elem=elec_elem+1
      do jj=1,nodpel
        electro_elem(elec_elem,jj)=conect(kk,jj)
        material(kk)=5
      enddo
    enddo


    kk = busco_arriba(nelements,nodpel,conect,electro_elem(4,:))  
    elec_elem=elec_elem+1
    do jj=1,nodpel
      electro_elem(elec_elem,jj)=conect(kk,jj)
      material(kk)=5
    enddo

    do ii=2,elec_alt-1
      kk = busco_arriba(nelements,nodpel,conect,electro_elem(elec_elem,:))  
      elec_elem=elec_elem+1
      do jj=1,nodpel
        electro_elem(elec_elem,jj)=conect(kk,jj)
        material(kk)=5
      enddo
    enddo

    kk = busco_derecha(nelements,nodpel,conect,electro_elem(4,:))  
    elec_elem=elec_elem+1
    do jj=1,nodpel
      electro_elem(elec_elem,jj)=conect(kk,jj)
      material(kk)=5
    enddo

    do ii=2,elec_alt-1
      kk = busco_arriba(nelements,nodpel,conect,electro_elem(elec_elem,:))  
      elec_elem=elec_elem+1
      do jj=1,nodpel
        electro_elem(elec_elem,jj)=conect(kk,jj)
        material(kk)=5
      enddo
    enddo

    kk = busco_izquierda(nelements,nodpel,conect,electro_elem(4,:))  
    elec_elem=elec_elem+1
    do jj=1,nodpel
      electro_elem(elec_elem,jj)=conect(kk,jj)
      material(kk)=5
    enddo

    do ii=2,elec_alt-1
      kk = busco_arriba(nelements,nodpel,conect,electro_elem(elec_elem,:))  
      elec_elem=elec_elem+1
      do jj=1,nodpel
        electro_elem(elec_elem,jj)=conect(kk,jj)
        material(kk)=5
      enddo
    enddo

else  ! si son dos electroditos apartados

   
   allocate(electro_elem(elec_alt*AGUJAS,8))



! busco elemento con ese nodo en el corner izquierdo-infrior-trasero
   elec_elem=0
   do kk=1,nelements
  
    do je=1,AGUJAS/2
      if(conect(kk,1)==nod_electro_pos(je)) then
        elec_elem=elec_elem+1
        do jj=1,nodpel
          electro_elem(elec_elem,jj)=conect(kk,jj)
          material(kk)=10
        enddo
      endif
    enddo
  enddo
   
  do kk=1,nelements
  
    do je=1,AGUJAS/2
      if(conect(kk,1)==nod_electro_neg(je)) then
        elec_elem=elec_elem+1
        do jj=1,nodpel
          electro_elem(elec_elem,jj)=conect(kk,jj)
          material(kk)=11
        enddo
      endif
    enddo


  enddo


   do je=1,AGUJAS/2

     kk = busco_arriba(nelements,nodpel,conect,electro_elem(je,:))  ! electrodo pos
     elec_elem=elec_elem+1
     do jj=1,nodpel
       electro_elem(elec_elem,jj)=conect(kk,jj)
       material(kk)=10
     enddo
  
     do ii=2,elec_alt-1
       kk = busco_arriba(nelements,nodpel,conect,electro_elem(elec_elem,:))  ! electrodo pos
       elec_elem=elec_elem+1
       do jj=1,nodpel
         electro_elem(elec_elem,jj)=conect(kk,jj)
         material(kk)=10
     enddo
    enddo
  enddo

   do je=1,AGUJAS/2

     kk = busco_arriba(nelements,nodpel,conect,electro_elem(AGUJAS/2+je,:))  ! electrodo pos
     elec_elem=elec_elem+1
     do jj=1,nodpel
       electro_elem(elec_elem,jj)=conect(kk,jj)
       material(kk)=11
     enddo
  
     do ii=2,elec_alt-1
       kk = busco_arriba(nelements,nodpel,conect,electro_elem(elec_elem,:))  ! electrodo pos
       elec_elem=elec_elem+1
       do  jj=1,nodpel
         electro_elem(elec_elem,jj)=conect(kk,jj)
         material(kk)=11
       enddo
     enddo
   enddo



endif

! condiciones de contorno 
nod_tierra=0
nod_poten=0
nod_temper=0

mat_normal=0
mat_tumor=0
mat_aislante=0
mat_electro=0
mat_agua=0
mat_externo=0
mat_extra=0
mat_tierra=0
mat_potencial=0
allocate(vec_tierra(nnodes),vec_poten(nnodes),vec_temper(nnodes))

vec_tierra=-1
vec_poten=-1
vec_temper=-1

 
do kk=1,nelements
  if(material(kk) == 0 ) then
     mat_externo=mat_externo+1
  elseif(material(kk) == 1 ) then
     mat_agua=mat_agua+1
  elseif(material(kk) == 2 .or. material(kk) == 3 ) then
     mat_normal=mat_normal+1
  elseif(material(kk) == 4 ) then
     mat_tumor=mat_tumor+1
  elseif(material(kk) == 5 ) then
     mat_aislante=mat_aislante+1
  elseif(material(kk) == 10) then
     mat_electro=mat_electro+1
     mat_potencial = mat_potencial+1
     do jj=1,nodpel
        j=conect(kk,jj)
        vec_poten(j)=potencial
     enddo  
  elseif(material(kk) == 11) then
     mat_electro=mat_electro+1
     mat_tierra = mat_tierra+1
     do jj=1,nodpel
        j=conect(kk,jj)
        vec_tierra(j)=tierra
     enddo  
  else
     mat_extra = mat_extra+1
  endif

enddo


do kk=1,nnodes
   if(vec_tierra(kk)/=-1) then
       nod_tierra=nod_tierra+1
   endif
   if(vec_poten(kk)/=-1) then
       nod_poten=nod_poten+1
   endif
enddo

! examino condiciones de contorno de la ecuacion biotermica
xliminf=0.0
xlimsup=0.0
yliminf=0.0
ylimsup=0.0
zliminf=0.0
zlimsup=0.0

do kk=1,nnodes
   if(coor_x(kk).lt.xliminf) xliminf = coor_x(kk)
   if(coor_x(kk).gt.xlimsup) xlimsup = coor_x(kk)
   if(coor_y(kk).lt.yliminf) yliminf = coor_y(kk)
   if(coor_y(kk).gt.ylimsup) ylimsup = coor_y(kk)
   if(coor_z(kk).lt.zliminf) zliminf = coor_z(kk)
   if(coor_z(kk).gt.zlimsup) zlimsup = coor_z(kk)
enddo

do kk=1,nnodes
   if(coor_x(kk)==xliminf) then
     nod_temper=nod_temper + 1
   elseif(coor_x(kk)==xlimsup)then
     nod_temper=nod_temper + 1
   elseif(coor_y(kk)==ylimsup)then
     nod_temper=nod_temper + 1
    elseif(coor_y(kk)==xliminf)then
     nod_temper=nod_temper + 1
   elseif(coor_z(kk)==zlimsup)then
     nod_temper=nod_temper + 1
   elseif(coor_z(kk)==zliminf)then
     nod_temper=nod_temper + 1
   endif
enddo

allocate(nnodtemper(nod_temper))

nod_temper=0

do kk=1,nnodes
   if(coor_x(kk)==xliminf) then
     nod_temper=nod_temper + 1
     nnodtemper(nod_temper)=kk
     vec_temper(kk)=1
   elseif(coor_x(kk)==xlimsup)then
     nod_temper=nod_temper + 1
     nnodtemper(nod_temper)=kk
     vec_temper(kk)=1
   elseif(coor_y(kk)==ylimsup)then
     nod_temper=nod_temper + 1
     nnodtemper(nod_temper)=kk
     vec_temper(kk)=1
    elseif(coor_y(kk)==xliminf)then
     nod_temper=nod_temper + 1
     nnodtemper(nod_temper)=kk
     vec_temper(kk)=1
   elseif(coor_z(kk)==zlimsup)then
     nod_temper=nod_temper + 1
     nnodtemper(nod_temper)=kk
     vec_temper(kk)=1
   elseif(coor_z(kk)==zliminf)then
     nod_temper=nod_temper + 1
     nnodtemper(nod_temper)=kk
     vec_temper(kk)=1
   endif
enddo



write(unit_cc,* )' materiales ', nelements
write(unit_cc,* )'      externo:   ', mat_externo
write(unit_cc,* )'      agua   :   ', mat_agua
write(unit_cc,* )'      tumor  :   ', mat_tumor
write(unit_cc,* )'      normal :   ', mat_normal
write(unit_cc,* )'      aislante:   ', mat_aislante
write(unit_cc,* )'      electrodo:  ', mat_electro,mat_potencial,mat_tierra
write(unit_cc,* )'      extra:  ', mat_extra

write(unit_cc,* )'  ************************************************  '
write(unit_cc,* )'  ************************************************  '



allocate(nnodtierra(nod_tierra),nnodpoten(nod_poten))
nnodtierra=0
nnodpoten=0

nod_tierra=0
nod_poten=0
do kk=1,nnodes
   if(vec_tierra(kk)/=-1) then
       nod_tierra=nod_tierra+1
       nnodtierra(nod_tierra)=kk
   endif
   if(vec_poten(kk)/=-1) then
       nod_poten=nod_poten+1
       nnodpoten(nod_poten)=kk
   endif

enddo

write(unit_cc,* )'nodos tierra ', nod_tierra
      
do kk=1,nod_tierra
      write(unit_cc,'(2i6,3e15.5)') kk,nnodtierra(kk),coor_x(nnodtierra(kk)), coor_y(nnodtierra(kk)),coor_z(nnodtierra(kk))
enddo

write(unit_cc,* )'nodos poten ', nod_poten
do kk=1,nod_poten
      write(unit_cc,'(2i6,3e15.5)') kk,nnodpoten(kk), coor_x(nnodpoten(kk)),coor_y(nnodpoten(kk)),coor_z(nnodpoten(kk))
enddo

write(unit_cc,* )'nodos temper ', nod_temper
do kk=1,nod_temper
      write(unit_cc,'(2i10,3e15.5)') kk,nnodtemper(kk), coor_x(nnodtemper(kk)),coor_y(nnodtemper(kk)),coor_z(nnodtemper(kk))
enddo


if(mat_extra/=0) then
   print*,'cuidado hay materiales sin clasificar',mat_extra
endif

total = mat_externo+mat_agua + mat_aislante + mat_normal + mat_tumor+mat_electro

print*,'total de elementos',total,nelements



call sacomalla(nnodes,coor_x,coor_y,coor_z,conect,nelements,nodpel,material,unit_malla)
  
 
allocate(grad_x(nnodes),grad_y(nnodes),grad_z(nnodes),gradxel_x(nelements),gradxel_y(nelements),gradxel_z(nelements))
allocate(tempera(nnodes),tempera_ant(nnodes),jcurrent(nelements),quemado(nelements))

tempera=Tblod
tempera_ant=Tblod
quemado=0

end subroutine leemalla


integer function busco_arriba(nelements,nodpel,conect,electro)  
implicit none
integer ::nelements,nodpel,conect(nelements,nodpel),electro(nodpel)
! local
integer :: k

do k=1,nelements
   
   if(conect(k,5)==electro(1) .and. conect(k,6)==electro(2) .and. conect(k,7)==electro(3) .and. conect(k,8)==electro(4) ) then
       busco_arriba = k
   endif

enddo


end function  busco_arriba


integer function busco_atras(nelements,nodpel,conect,electro)  
implicit none
integer ::nelements,nodpel,conect(nelements,nodpel),electro(nodpel)
! local
integer :: k

do k=1,nelements
   
   if(conect(k,2)==electro(1) .and. conect(k,3)==electro(4) .and. conect(k,7)==electro(8) .and. conect(k,6)==electro(5) ) then
       busco_atras = k
   endif

enddo


end function  busco_atras


integer function busco_adelante(nelements,nodpel,conect,electro)  
implicit none
integer ::nelements,nodpel,conect(nelements,nodpel),electro(nodpel)
! local
integer :: k

do k=1,nelements
   
   if(conect(k,1)==electro(2) .and. conect(k,4)==electro(3) .and. conect(k,5)==electro(6) .and. conect(k,8)==electro(7) ) then
       busco_adelante = k
   endif

enddo


end function  busco_adelante


integer function busco_izquierda(nelements,nodpel,conect,electro)  
implicit none
integer ::nelements,nodpel,conect(nelements,nodpel),electro(nodpel)
! local
integer :: k

do k=1,nelements
   
   if(conect(k,4)==electro(1) .and. conect(k,3)==electro(2) .and. conect(k,7)==electro(6) .and. conect(k,8)==electro(5) ) then
       busco_izquierda = k
   endif

enddo


end function  busco_izquierda


integer function busco_derecha(nelements,nodpel,conect,electro)  
implicit none
integer ::nelements,nodpel,conect(nelements,nodpel),electro(nodpel)
! local
integer :: k

do k=1,nelements
   
   if(conect(k,1)==electro(4) .and. conect(k,2)==electro(3) .and. conect(k,5)==electro(8) .and. conect(k,6)==electro(7) ) then
       busco_derecha = k
   endif

enddo


end function  busco_derecha


subroutine expando(ntum,kk,ns,nel_tum,nlimite,nelet,aux_x,aux_y,aux_z,mataux)
use def_variables
use def_constantes
implicit none
integer :: ntum,kk,nlimite,nelet,nel_tum(nlimite,8),ns(8),mataux(nlimite)
double precision :: aux_x(*),aux_y(*),aux_z(*)
! local
integer :: jz,jy,jx,mapa(nmallo_Xnod+1,nmallo_Xnod+1,nmallo_Xnod+1),n_new,buscon_nod,nlimit
double precision :: x_x,y_y,z_z,xpaso

nlimit=nmallo_Xnod+1 
xpaso = real(ncada_Xnod)/real(nmallo_Xnod)

do jz=1,nlimit 
  
   do jy=1,nlimit

      do jx=1,nlimit
         
          x_x=coor_x(ns(1))+xpaso*(jx-1)
          y_y=coor_y(ns(1))+xpaso*(jy-1)
          z_z=coor_z(ns(1))+xpaso*(jz-1)


         if(jz==1) then
            if(jy==1) then
               if(jx==1) then

                 mapa(jx,jy,jz) = ns(1)
               elseif(jx==nlimit) then
                 mapa(jx,jy,jz) = ns(2)
               else
                 
                 n_new = buscon_nod(ntum,x_x,y_y,z_z,aux_x,aux_y,aux_z)
                 if(n_new==0) then
                    ntum=ntum+1
                    mapa(jx,jy,jz) = nnodes+ntum
                    aux_x(ntum)=x_x
                    aux_y(ntum)=y_y
                    aux_z(ntum)=z_z

                 else
                    mapa(jx,jy,jz)= nnodes+n_new
                    
                 endif   
               endif
                                   
            elseif(jy==nlimit) then
               if(jx==1) then

                 mapa(jx,jy,jz) = ns(4)
               elseif(jx==nlimit) then
                 mapa(jx,jy,jz) = ns(3)
               else
                 n_new = buscon_nod(ntum,x_x,y_y,z_z,aux_x,aux_y,aux_z)
                 if(n_new==0) then
                    ntum=ntum+1
                    mapa(jx,jy,jz) = nnodes+ntum
                    aux_x(ntum)=x_x
                    aux_y(ntum)=y_y
                    aux_z(ntum)=z_z
                 else
                   mapa(jx,jy,jz) = nnodes+n_new
                    
                 endif   
               endif
            else
                n_new = buscon_nod(ntum,x_x,y_y,z_z,aux_x,aux_y,aux_z)
                 if(n_new==0) then
                    ntum=ntum+1
                    mapa(jx,jy,jz) = nnodes+ntum
                    aux_x(ntum)=x_x
                    aux_y(ntum)=y_y
                    aux_z(ntum)=z_z
                 else
                   mapa(jx,jy,jz)= nnodes+n_new
                    
                 endif   

            endif
       elseif(jz==nlimit) then
            if(jy==1) then
               if(jx==1) then

                 mapa(jx,jy,jz) = ns(5)
               elseif(jx==nlimit) then
                mapa(jx,jy,jz) = ns(6)
               else
                 n_new = buscon_nod(ntum,x_x,y_y,z_z,aux_x,aux_y,aux_z)
                 if(n_new==0) then
                    ntum=ntum+1
                    mapa(jx,jy,jz) = nnodes+ntum
                    aux_x(ntum)=x_x
                    aux_y(ntum)=y_y
                    aux_z(ntum)=z_z
                 else
                   mapa(jx,jy,jz)= nnodes+n_new
                    
                 endif   
               endif
                                   
            elseif(jy==nlimit) then
               if(jx==1) then

                mapa(jx,jy,jz) = ns(8)
               elseif(jx==nlimit) then
                 mapa(jx,jy,jz) = ns(7)
               else
                n_new = buscon_nod(ntum,x_x,y_y,z_z,aux_x,aux_y,aux_z)
                 if(n_new==0) then
                    ntum=ntum+1
                    mapa(jx,jy,jz) = nnodes+ntum
                    aux_x(ntum)=x_x
                    aux_y(ntum)=y_y
                    aux_z(ntum)=z_z
                 else
                    mapa(jx,jy,jz) = nnodes+n_new
                    
                 endif   
               endif
            else
                 n_new = buscon_nod(ntum,x_x,y_y,z_z,aux_x,aux_y,aux_z)
                 if(n_new==0) then
                    ntum=ntum+1
                    mapa(jx,jy,jz) = nnodes+ntum
                    aux_x(ntum)=x_x
                    aux_y(ntum)=y_y
                    aux_z(ntum)=z_z
                 else
                    mapa(jx,jy,jz) = nnodes+n_new
                    
                 endif   


            endif

         else
              n_new = buscon_nod(ntum,x_x,y_y,z_z,aux_x,aux_y,aux_z)
                 if(n_new==0) then
                    ntum=ntum+1
                    mapa(jx,jy,jz) = nnodes+ntum
                    aux_x(ntum)=x_x
                    aux_y(ntum)=y_y
                    aux_z(ntum)=z_z
                 else
                    mapa(jx,jy,jz) = nnodes+n_new
                    
                 endif   
         
         endif 
          
          

      enddo
   enddo
enddo


! creo elementos
! el primero en el lugar anterior
!conect(kk,1)=mapa(1,1,1)
!conect(kk,2)=mapa(1,2,1)
!conect(kk,3)=mapa(2,2,1)
!conect(kk,4)=mapa(2,1,1)
!conect(kk,5)=mapa(1,1,2)
!conect(kk,6)=mapa(1,2,2)
!conect(kk,7)=mapa(2,2,2)
!conect(kk,8)=mapa(2,1,2)


do jz=1,nmallo_Xnod 
  
   do jy=1,nmallo_Xnod

      do jx=1,nmallo_Xnod
         nelet=nelet+1
         nel_tum(nelet,1) = mapa(jx,jy,jz)
         nel_tum(nelet,2) = mapa(jx+1,jy,jz)
         nel_tum(nelet,3) = mapa(jx+1,jy+1,jz)
         nel_tum(nelet,4) = mapa(jx,jy+1,jz)

         nel_tum(nelet,5) = mapa(jx,jy,jz+1)
         nel_tum(nelet,6) = mapa(jx+1,jy,jz+1)
         nel_tum(nelet,7) = mapa(jx+1,jy+1,jz+1)
         nel_tum(nelet,8) = mapa(jx,jy+1,jz+1)
   
         mataux(nelet)=material(kk)

      enddo
   enddo

enddo

end subroutine expando


integer function buscon_nod(ntum,x_x,y_y,z_z,aux_x,aux_y,aux_z)
use def_variables
implicit none
integer :: ntum
double precision::x_x,y_y,z_z,aux_x(*),aux_y(*),aux_z(*)
! local
integer :: jj

buscon_nod=0
do jj=1,ntum
  if(aux_x(jj)>=x_x*0.99999 .and.aux_x(jj)<=x_x*1.0001 .and. aux_y(jj)>=y_y*0.99999 .and. aux_y(jj)<=y_y*1.00001  &
     &                      .and. aux_z(jj)>=z_z*0.99999 .and. aux_z(jj)<=z_z*1.00001) then
      buscon_nod = jj
  endif
enddo

end function buscon_nod
