subroutine sparse_logic
use def_io
use def_variables
use def_vectors
! local
integer, allocatable :: consim(:)
integer :: nno,kk,npas



allocate( AD(NP),IA(NP+1),ncount(NP),indep_vect(NP))
allocate( ICX(NP+1),consim(NP) )
allocate(ns(nodpel))

ad=0.0
indep_vect=0.0

ia=0
ncount=0
icx=0

IA(1)=1

do nno=1,np
  consim=0
  DO KK=1,ne
   npas=0
   
   DO I=1,nodpel
      ns(I)=conn(KK,I)
      if(nno==ns(i)) npas=1
   ENDDO

   if(npas==1) then   

      DO I=1,nodpel
        DO J=1,nodpel
             if(nno==ns(i).and.nno/=ns(j)) consim(ns(j)) = consim(ns(j))+1 
        ENDDO
      ENDDO
    endif
  enddo

  do i=1,NP
     if(consim(i)/=0) then
        ncount(nno)=ncount(nno)+1
     endif
  enddo
  IA(nno+1) = IA(nno) + ncount(nno)
enddo

NONULL = IA(NP+1)

allocate (an(NONULL),ja(NONULL))
An=0.0

do nno=1,np
  consim=0
  DO KK=1,ne
     npas=0
     DO I=1,nodpel
        NS(I)=conn(KK,I)
        if(nno==ns(i)) npas=1
     ENDDO
     if(npas==1) then   
        
        DO I=1,nodpel
          DO J=1,nodpel
               if(nno==ns(i).and.nno/=ns(j)) consim(ns(j)) = consim(ns(j))+1 
          ENDDO
        ENDDO

     endif
  ENDDO
  do i=1,np
      if(consim(i)/=0) then
          JA(IA(nno)+ICX(nno)) =   i
          ICX(nno)=ICX(nno)+1
      endif
   enddo
enddo


write(control_unit,*) 'Nonulos del sistema: ',NONULL

deallocate(consim)

end subroutine sparse_logic


!SUBROUTINE armoMasa()
!USE def_variables
!USE def_vectors
!IMPLICIT NONE
!
!local
!INTEGER :: i, j, k, ns(nodpel)
!DOUBLE PRECISION :: x(ndim,nodpel), dphi(ndim, nodpel, nodpel), weight(nodpel), AJACO(2,2), DETER, rho
!DOUBLE PRECISION :: funrho
!    
!    masa = 0.0
!
!    DO i=1, ne
!
!        DO j=1,nodpel
!            ns(j)=conn(i,j)
!
!            x(1,j)=coorx(ns(j))
!            x(2,j)=coory(ns(j))
!        END DO
!
!        CALL shape_nodos(ndim, nodpel, dphi, weight)
!
!        DO j=1,nodpel 
!            AJACO = 0.0
!            DO k = 1,nodpel
!                AJACO(1,1) = AJACO(1,1) + x(1,k) * dphi(1,j,k)
!                AJACO(1,2) = AJACO(1,2) + x(1,k) * dphi(2,j,k)
!                AJACO(2,1) = AJACO(2,1) + x(2,k) * dphi(1,j,k)
!                AJACO(2,2) = AJACO(2,2) + x(2,k) * dphi(2,j,k)
!            END DO
! 
!            DETER = AJACO(1,1) * AJACO(2,2) - AJACO(2,1) * AJACO(1,2)
!         
!            masa(ns(j)) = masa(ns(j)) + weight(j)*DETER
!            masa(ns(j)) = masa(ns(j)) !* 0.0001 !Correccion de unidades
!
!        END DO
!    END DO
!    
!
!    DO i=1,np
!        IF ( masa(i)<1E-8 ) THEN
!            PAUSE 'NODO CON MASA CERO'
!            STOP
!        END IF
!   END DO

!END SUBROUTINE armoMasa