SUBROUTINE CONTROL()
use def_solver
use def_constantes
use def_variables
implicit none
double precision :: x(nodpel),y(nodpel),z(nodpel),esm(nodpel,nodpel),ef(nodpel),adiag,sigma_el,qe
integer :: ns(nodpel)
INTEGER  ::INBWT,nbwt,i,ij,j,nb,kk,NCASE,inode,ipoin,jnode,ii,jj,j1,nicio,k,nno,npas
double precision, allocatable :: ud(:),un(:)
integer, allocatable :: cont(:),ip(:),iu(:),iup(:),ju(:),iut(:),consim(:)
integer:: unit_sist
real*4, ALLOCATABLE :: AGENERAL(:)

! ARMO LAS ESTRUCTURAS LOGICAS SPARCE


        allocate(grad_x(nnodes),grad_y(nnodes),gradxel_x(nelements),gradxel_y(nelements))
        ALLOCATE(AD(nnodes),IA(nnodes+1),CONT(nnodes+1),CX(nnodes+1),solucion(nnodes),rhs(nnodes),consim(nnodes),masa(nnodes),  &
         solucion_ant(nnodes),rhs_ant(nnodes))

    !    allocate(cplx_grad_x(nnodes),cplx_grad_y(nnodes),cplx_gradxel_x(nelements),cplx_gradxel_y(nelements))
  

!ALLOCATE(IP(nnodes),IUP(nnodes),IU(nnodes),UD(nnodes),IUT(nnodes))      


!if(nodpel==16) then
!   call masadiag2d()
!else
!   call masadiag()
!endif


RHS=0.0
rhs_ant=0.0
ad=0.0
solucion=0.0
solucion_ant=0.0
!ud=0.0
ia=0
!ip=0
cx=0
!iu=0
!iut=0
cont=0
consim=0


      INBWT=0
      NBWT =0
      DO KK=1,NE
        DO I=1,nodpel
          NS(I)=conect(KK,I)
        ENDDO
        DO  I=1,nodpel-1
		  IJ=I+1
          DO J=IJ,nodpel
            NB=IABS(NS(I)-NS(J))
            IF(NB.EQ.0) THEN 
			   WRITE(unit_cont,*) 'ELEMENTO  ',KK,' TIENE DOS NODOS IDENTICOS'
			   WRITE(6,*) 'ELEMENTO  ',KK,' TIENE DOS NODOS IDENTICOS'
			   STOP
			ENDIF    
            IF(NB.GT.NBWT) THEN
               INBWT=KK
               NBWT =NB
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      NBWT=NBWT+1

      WRITE(unit_cont,*) ' BANDWIDTH: ',NBWT,'  EN ELEMENTO  ',INBWT
   
   
! DETERMINO LA FORMA SIMBOLICA DEL PROBLEMA SPARCE. PARA ESO ENSAMBLO UN SISTEMA DE UNOS Y CEROS


!ALLOCATE( AGENERAL(nnodes*NBWT) )
!AGENERAL = 0

DO KK=1,nnodes+1
  CX(KK)=0.0
  CONT(KK)=0.0
  IA(KK) = 0
ENDDO


NCASE=0
sigma_el=1.0
qe=1.0
IA(1)=1

do nno=1,nnodes
  consim=0
  DO KK=1,nelements
     npas=0
     DO I=1,nodpel
       NS(I)=conect(KK,I)
       if(nno==ns(i)) npas=1
     ENDDO
     if(npas==1) then   
       ESM=1.0
 

       DO I=1,nodpel
          II=NS(I)
            DO J=1,nodpel
              JJ=NS(J)+1-II
              IF(JJ.GT.0) THEN
                 !J1 = (JJ-1)*nnodes + II - (JJ-1)*(JJ-2)/2
	             if(nno==ns(i).and.nno/=ns(j)) consim(ns(j)) = consim(ns(j))+1          
                 !AGENERAL(J1) = AGENERAL(J1)+ ESM(I,J)
              ENDIF  
	        ENDDO
        ENDDO

      endif
   ENDDO

   do i=1,nnodes
      if(consim(i)/=0) then
         cont(nno)=cont(nno)+1
      endif
   enddo
    IA(nno+1) = IA(nno) + CONT(nno)
enddo


!IA(1)=1
!DO I=1,nnodes
!  IA(I+1) = IA(I) + CONT(I)
!ENDDO

NONULL = IA(nnodes+1)-1 !- nnodes


!ia=0
!cont=0

!NICIO = nnodes

!DO I=1,NBWT  !NP-1
!  DO K=1,nnodes-I
!     IF(AGENERAL(NICIO+K).NE.0.0)  CONT(K) = CONT(K)+1
!  ENDDO
!     NICIO = NICIO + (nnodes - I )
!ENDDO

!IA(1)=1
!DO I=1,nnodes
!  IA(I+1) = IA(I) + CONT(I)
!ENDDO

!NONULL = IA(nnodes+1)-1

ALLOCATE (an(NONULL),ja(NONULL))
an=0.0
ja=0
!NICIO = nnodes

!DO I=1,NBWT-1
!  DO K=1,nnodes-I
!   IF(AGENERAL(NICIO+K).NE.0.0) THEN
!      JA(IA(K)+CX(K)) =   K + I
      !AN(IA(K)+CX(K)) = AGENERAL(NICIO + K)
!      CX(K)=CX(K)+1
!    ENDIF
!   ENDDO
!   NICIO = NICIO + (nnodes - I)
!ENDDO

!DEALLOCATE(AGENERAL) 

!cx=0
!ja=0
write(unit_cont,*) 'nonulos del sistema  ',nnodes,nonull,nonull*(1.0/real(nnodes))*(1/real(nnodes))

do nno=1,nnodes
  consim=0
  DO KK=1,nelements
     npas=0
     DO I=1,nodpel
       NS(I)=conect(KK,I)
       if(nno==ns(i)) npas=1
     ENDDO
     if(npas==1) then   
       ESM=1.0
   
     DO I=1,nodpel
       II=NS(I)
       DO J=1,nodpel
          JJ=NS(J)+1-II
          IF(JJ.GT.0) THEN
            if(nno==ns(i).and.nno/=ns(j)) consim(ns(j)) = consim(ns(j))+1          
          ENDIF  
	   ENDDO
     ENDDO

   endif
  ENDDO
  do i=1,nnodes
      
      if(consim(i)/=0) then
          JA(IA(nno)+CX(nno)) =   i
          CX(nno)=CX(nno)+1
      endif
   enddo
enddo


open(unit=unit_sist,file=archi_sistema)
write(unit_sist,*) nnodes
do kk=1,nnodes

   write(unit_sist,*) kk,masa(kk)

enddo

write(unit_sist,*) nnodes+1
do kk=1,nnodes+1

   write(unit_sist,*) kk,ia(kk),cx(kk)

enddo

write(unit_sist,*) nonull
do kk=1,nonull

   write(unit_sist,*) kk,ja(kk)

enddo
close(unit_sist)
RETURN
    END


    
    
SUBROUTINE control_cplx()
use def_solver
use def_constantes
use def_variables
implicit none
double precision :: x(nodpel),y(nodpel),z(nodpel)
complex*16:: esm_cplx(nodpel,nodpel),ef_cplx(nodpel),adiag,sigma_el,qe
integer :: ns(nodpel)
INTEGER  ::INBWT,nbwt,i,ij,j,nb,kk,NCASE,inode,ipoin,jnode,ii,jj,j1,nicio,k,nno,npas
integer, allocatable :: cont(:),consim(:)
integer:: unit_sist
complex, ALLOCATABLE :: AGENERAL_cplx(:)

! ARMO LAS ESTRUCTURAS LOGICAS SPARCE


  ALLOCATE(cplx_AD(nnodes),IA(nnodes+1),CONT(nnodes+1),CX(nnodes+1),cplx_solucion(nnodes),cplx_rhs(nnodes),consim(nnodes),  &
         cplx_solucion_ant(nnodes),cplx_rhs_ant(nnodes))

  allocate(cplx_grad_x(nnodes),cplx_grad_y(nnodes),cplx_gradxel_x(nelements),cplx_gradxel_y(nelements))
  


!if(nodpel==16) then
!   call masadiag2d()
!else
!   call masadiag()
!endif


cplx_RHS=0.0
cplx_rhs_ant=0.0
cplx_ad=0.0
cplx_solucion=0.0
cplx_solucion_ant=0.0
!ud=0.0
ia=0
!ip=0
cx=0
!iu=0
!iut=0
cont=0
consim=0


      INBWT=0
      NBWT =0
      DO KK=1,NE
        DO I=1,nodpel
          NS(I)=conect(KK,I)
        ENDDO
        DO  I=1,nodpel-1
		  IJ=I+1
          DO J=IJ,nodpel
            NB=IABS(NS(I)-NS(J))
            IF(NB.EQ.0) THEN 
			   WRITE(unit_cont,*) 'ELEMENTO  ',KK,' TIENE DOS NODOS IDENTICOS'
			   WRITE(6,*) 'ELEMENTO  ',KK,' TIENE DOS NODOS IDENTICOS'
			   STOP
			ENDIF    
            IF(NB.GT.NBWT) THEN
               INBWT=KK
               NBWT =NB
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      NBWT=NBWT+1

      WRITE(unit_cont,*) ' BANDWIDTH: ',NBWT,'  EN ELEMENTO  ',INBWT
   
   
! DETERMINO LA FORMA SIMBOLICA DEL PROBLEMA SPARCE. PARA ESO ENSAMBLO UN SISTEMA DE UNOS Y CEROS


!ALLOCATE( AGENERAL(nnodes*NBWT) )
!AGENERAL = 0

DO KK=1,nnodes+1
  CX(KK)=0.0
  CONT(KK)=0.0
  IA(KK) = 0
ENDDO


NCASE=0
sigma_el=1.0
qe=1.0
IA(1)=1

do nno=1,nnodes
  consim=0
  DO KK=1,nelements
     npas=0
     DO I=1,nodpel
       NS(I)=conect(KK,I)
       if(nno==ns(i)) npas=1
     ENDDO
     if(npas==1) then   
       esm_cplx=1.0
 

       DO I=1,nodpel
          II=NS(I)
            DO J=1,nodpel
              JJ=NS(J)+1-II
              IF(JJ.GT.0) THEN
                 !J1 = (JJ-1)*nnodes + II - (JJ-1)*(JJ-2)/2
	             if(nno==ns(i).and.nno/=ns(j)) consim(ns(j)) = consim(ns(j))+1          
                 !AGENERAL(J1) = AGENERAL(J1)+ ESM(I,J)
              ENDIF  
	        ENDDO
        ENDDO

      endif
   ENDDO

   do i=1,nnodes
      if(consim(i)/=0) then
         cont(nno)=cont(nno)+1
      endif
   enddo
    IA(nno+1) = IA(nno) + CONT(nno)
enddo


!IA(1)=1
!DO I=1,nnodes
!  IA(I+1) = IA(I) + CONT(I)
!ENDDO

NONULL = IA(nnodes+1)-1 !- nnodes


!ia=0
!cont=0

!NICIO = nnodes

!DO I=1,NBWT  !NP-1
!  DO K=1,nnodes-I
!     IF(AGENERAL(NICIO+K).NE.0.0)  CONT(K) = CONT(K)+1
!  ENDDO
!     NICIO = NICIO + (nnodes - I )
!ENDDO

!IA(1)=1
!DO I=1,nnodes
!  IA(I+1) = IA(I) + CONT(I)
!ENDDO

!NONULL = IA(nnodes+1)-1

ALLOCATE (cplx_an(NONULL),ja(NONULL))
cplx_an=0.0
ja=0
!NICIO = nnodes

!DO I=1,NBWT-1
!  DO K=1,nnodes-I
!   IF(AGENERAL(NICIO+K).NE.0.0) THEN
!      JA(IA(K)+CX(K)) =   K + I
      !AN(IA(K)+CX(K)) = AGENERAL(NICIO + K)
!      CX(K)=CX(K)+1
!    ENDIF
!   ENDDO
!   NICIO = NICIO + (nnodes - I)
!ENDDO

!DEALLOCATE(AGENERAL) 

!cx=0
!ja=0
write(unit_cont,*) 'nonulos del sistema  ',nnodes,nonull,nonull*(1.0/real(nnodes))*(1/real(nnodes))

do nno=1,nnodes
  consim=0
  DO KK=1,nelements
     npas=0
     DO I=1,nodpel
       NS(I)=conect(KK,I)
       if(nno==ns(i)) npas=1
     ENDDO
     if(npas==1) then   
       ESM_cplx=1.0
   
     DO I=1,nodpel
       II=NS(I)
       DO J=1,nodpel
          JJ=NS(J)+1-II
          IF(JJ.GT.0) THEN
            if(nno==ns(i).and.nno/=ns(j)) consim(ns(j)) = consim(ns(j))+1          
          ENDIF  
	   ENDDO
     ENDDO

   endif
  ENDDO
  do i=1,nnodes
      
      if(consim(i)/=0) then
          JA(IA(nno)+CX(nno)) =   i
          CX(nno)=CX(nno)+1
      endif
   enddo
enddo


RETURN
END subroutine control_cplx
