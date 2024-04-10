subroutine leo_mallado()
use def_variables
use def_constantes
use def_bio
implicit none
!local
integer :: kk,jj,n,nele_aux,n1,n2,n3,materiales
double precision :: ymed,xmed

nodpel=4

!area_trans=2*(3.141592)*radio*XX_alto
area_trans=4*XX_alto

open(unit=1111,file=archi_malla)

read(1111,*) nnodes

allocate(coor_x(nnodes),coor_y(nnodes),coor_z(nnodes))
do kk=1,nnodes
   read(1111,*) n, coor_x(kk),coor_y(kk),coor_z(kk)
enddo

read(1111,*) materiales
read(1111,*) mat_aislante
read(1111,*) mat_electro
read(1111,*) mat_normal

read(1111,*) nelements

allocate(material(nelements),conect(nelements,nodpel),quemado(nelements))
quemado=0

do kk=1,nelements
   read(1111,*) n,material(kk),conect(kk,1),conect(kk,2),conect(kk,3),conect(kk,4)
enddo


allocate(vec_tierra(nnodes),vec_poten(nnodes),vec_temper(nnodes))

vec_tierra=-1
vec_poten=-1
vec_temper=-1

read(1111,*) nele_aux

do kk=1,nele_aux
    read(1111,*) n1,n2,n3
    vec_temper(n1)=1
    vec_temper(n2)=1
    vec_temper(n3)=1
enddo

close(1111)


! asigno 25 grados a lo que esta fuera de la papa

do kk=1,nelements
   xmed=0.0
   do jj=1,nodpel
      xmed = xmed + coor_x(conect(kk,jj))/real(nodpel)
   enddo
   if(xmed>3.5) then
      do jj=1,nodpel
            vec_temper(conect(kk,jj))=1
      enddo
  elseif(xmed<2.0) then 
      do jj=1,nodpel
            vec_temper(conect(kk,jj))=1
      enddo
  endif
enddo

do kk=1,nelements
    if(material(kk)==2) then  ! electrodo
        ymed=0.0
        do jj=1,nodpel
           ymed=ymed + coor_y(conect(kk,jj))/real(nodpel)
        enddo
        if(ymed>0.0) then 
           do jj=1,nodpel
              vec_poten(conect(kk,jj))=1   ! positivo
           enddo
        else
           do jj=1,nodpel
              vec_tierra(conect(kk,jj))=1   ! negativo
           enddo
        endif
   endif
enddo

nod_tierra=0
nod_poten=0
nod_temper=0
do kk=1,nnodes
   if(vec_tierra(kk)==1) nod_tierra=nod_tierra+1
   if(vec_poten(kk)==1) nod_poten=nod_poten+1
   if(vec_temper(kk)==1) nod_temper=nod_temper+1
enddo

allocate(nnodtierra(nod_tierra),nnodpoten(nod_poten))
allocate(nnodtemper(nod_temper))

nod_tierra=0
nod_poten=0
nod_temper=0
do kk=1,nnodes
   if(vec_tierra(kk)==1) then
      nod_tierra=nod_tierra+1
      nnodtierra(nod_tierra)=kk
   endif
   if(vec_poten(kk)==1) then
      nod_poten=nod_poten+1
      nnodpoten(nod_poten)=kk
   endif
   
   if(vec_temper(kk)==1) then
       nod_temper=nod_temper+1
       nnodtemper(nod_temper)=kk
   endif
enddo


allocate(grad_x(nnodes),grad_y(nnodes),grad_z(nnodes),gradxel_x(nelements),gradxel_y(nelements),gradxel_z(nelements))
allocate(tempera(nnodes),tempera_ant(nnodes),jcurrent(nelements))
tempera=Tout
tempera_ant=Tout
grad_x=0
grad_y=0
grad_z=0

jcurrent=0
gradxel_x=0
gradxel_y=0
gradxel_z=0

! saco data
write(unit_malla,*) '********************* malla leida!'
write(unit_malla,*) nnodes
do kk=1,nnodes
   write(unit_malla,*) n, coor_x(kk),coor_y(kk),coor_z(kk)  
enddo

 write(unit_malla,*) materiales
do kk=1,materiales
    write(unit_malla,*) mat_aislante
    write(unit_malla,*) mat_electro
    write(unit_malla,*) mat_normal
enddo

write(unit_malla,*) nelements

do kk=1,nelements
    write(unit_malla,*) n,material(kk),conect(kk,1),conect(kk,2),conect(kk,3),conect(kk,4)
enddo

write(unit_cc,*) '********************* C Contorno!'

write(unit_cc,*) '  tierra '
do kk=1,nod_tierra
    write(unit_cc,*) nnodtierra(kk),vec_tierra(nnodtierra(kk))
enddo

write(unit_cc,*) '  tierra '
do kk=1,nod_tierra
    write(unit_cc,*) nnodtierra(kk)
enddo

write(unit_cc,*) '  potencial '
do kk=1,nod_poten
    write(unit_cc,*) nnodpoten(kk),vec_poten(nnodpoten(kk))
enddo

write(unit_cc,*) '  temperatura '
do kk=1,nod_temper
    write(unit_cc,*) nnodtemper(kk),vec_temper(nnodtemper(kk))
enddo


end subroutine leo_mallado



subroutine leo_thorio()
use def_variables
use def_constantes
use def_bio
implicit none
!local
integer :: kk,jj,n,nele_aux,n1,n2,n3,materiales
double precision :: x1,y1,z1

nodpel=4
materiales=3

open(unit=1111,file=archi_malla)

read(1111,*) nnodes

allocate(coor_x(nnodes),coor_y(nnodes),coor_z(nnodes))
do kk=1,nnodes
   read(1111,*) n, coor_x(kk),coor_y(kk),coor_z(kk)
enddo

mat_aislante=0
mat_electro=0
mat_normal=0

read(1111,*) nelements

allocate(material(nelements),conect(nelements,nodpel))

do kk=1,nelements
   read(1111,*) n,material(kk),conect(kk,1),conect(kk,2),conect(kk,3),conect(kk,4)
   if(material(kk)==1) mat_normal=mat_normal+1
   if(material(kk)==2) mat_electro=mat_electro+1
   if(material(kk)==3) then
       mat_aislante=mat_aislante+1
   !    material(kk)=1
   endif
enddo


allocate(vec_tierra(nnodes))

vec_tierra=-1

read(1111,*) nele_aux

do kk=1,nele_aux
    read(1111,*) n1,x1,y1,z1
    vec_tierra(n1)=1
enddo

close(1111)


! 

! no separo condiciones para los cilindors internos.
!do kk=1,nelements
!   if(material(kk)==3) then
!      do jj=1,nodpel
!            vec_tierra(conect(kk,jj))=2
!      enddo
!  endif
!enddo


nod_tierra=0
do kk=1,nnodes
   if(vec_tierra(kk)>=1) nod_tierra=nod_tierra+1
enddo

allocate(nnodtierra(nod_tierra))

nod_tierra=0
do kk=1,nnodes
   if(vec_tierra(kk)>=1) then
      nod_tierra=nod_tierra+1
      nnodtierra(nod_tierra)=kk
   endif
enddo


allocate(grad_x(nnodes),grad_y(nnodes),grad_z(nnodes),gradxel_x(nelements),gradxel_y(nelements),gradxel_z(nelements))
grad_x=0
grad_y=0
grad_z=0

gradxel_x=0
gradxel_y=0
gradxel_z=0

! saco data
write(unit_malla,*) '********************* malla leida!'
write(unit_malla,*) nnodes
do kk=1,nnodes
   write(unit_malla,*) kk, coor_x(kk),coor_y(kk),coor_z(kk)  
enddo

 write(unit_malla,*) 'materiales: ',materiales
 write(unit_malla,*) 'Agua: ',mat_aislante
 write(unit_malla,*) 'thorio: ',mat_electro
 write(unit_malla,*) 'Grafito: ',mat_normal

write(unit_malla,*) nelements

do kk=1,nelements
    write(unit_malla,*) kk,material(kk),conect(kk,1),conect(kk,2),conect(kk,3),conect(kk,4)
enddo

write(unit_cc,*) '********************* C Contorno!'

write(unit_cc,*) '  tierra '
do kk=1,nod_tierra
    write(unit_cc,*) nnodtierra(kk),vec_tierra(nnodtierra(kk))
enddo

write(unit_cc,*) '  tierra '
do kk=1,nod_tierra
    write(unit_cc,*) nnodtierra(kk)
enddo


end subroutine leo_thorio