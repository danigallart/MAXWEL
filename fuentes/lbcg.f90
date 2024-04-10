
SUBROUTINE CG(NP,IA,JA,A_SPA,AD,B,X,TOL,ITMAX,ITER,ERR)
      implicit none 
	DOUBLE PRECISION A_SPA(*),AD(*),X(*),B(*)
      INTEGER IA(*),JA(*)
      DOUBLE PRECISION bnrm, Znrm, bknum, bkden, XANT, akden,bk,ak,zm1nrm,dxnrm,xnrm,alfa,b_norma,r_norma,aknum
	DOUBLE PRECISION TOL,EPS,ERR,FUNNORM
	INTEGER ITMAX,ITER,NP,j,kk

      DOUBLE PRECISION, ALLOCATABLE ::Z(:),P(:),PP(:),R(:),ZZ(:),RR(:)

	ALLOCATE (P(NP),PP(NP),ZZ(NP),Z(NP),R(NP),RR(NP)) 

      
	CALL MATXVECSIM(IA,JA,A_SPA,AD,X,R,NP)

      b_norma = 0.0
	DO KK=1,NP
        R(KK)  = B(KK) - R(KK)
        b_norma =b_norma  + b(kk)*b(kk) 
	ENDDO
      b_norma=sqrt(b_norma) 
      

	DO KK=1,NP
	  Z(KK)=R(KK)/AD(KK)
        rr(kk)=z(kk)
	ENDDO

      ZNRM=FUNNORM(NP,Z)
      
      do while(err.gt.tol .and. iter.lt.itmax)

        iter=iter+1

	  bknum=0.d0
  	  do j=1,NP
          bknum=bknum + Z(j)*R(j)
        enddo
	  CALL MATXVECSIM(IA,JA,A_SPA,AD,rr,Z,NP)
        
	  bkden=0.0
	  do j=1,NP
          bkden=bkden + rr(j)*Z(j)
        enddo
	   
        alfa= bknum/bkden    	  

        do kk=1,np
          x(kk)=x(kk)+alfa*rr(kk)          
	    r(kk)= r(kk) - alfa * z(kk)
	  enddo

     	  DO KK=1,NP
	    Z(KK)=R(KK)/AD(KK)
        ENDDO

    	  akden=bknum
	  aknum=0.0
        do  j=1,NP
          aknum = aknum+ z(j)*r(j)
        enddo

        ak=aknum/akden
          
        DO J=1,NP
          rr(j) = z(j)  + ak*rr(j)
        ENDDO

        r_norma = 0.0
	  DO KK=1,NP
          r_norma =r_norma  + r(kk)*r(kk)    
	  ENDDO
        r_norma=sqrt(r_norma) 

        err=r_norma/b_norma
        !write(2,*) iter,err,r_norma,ak,alfa

      ENDDO


	DEALLOCATE (P,PP,ZZ,Z,R,RR) 
	
	end subroutine cg
      
   DOUBLE PRECISION FUNCTION FUNNORM(N,Y)
   implicit none
   DOUBLE PRECISION:: Y(*)
   integer:: N
   !local
   integer :: isamax,i
   DOUBLE PRECISION:: sum2
   
        !ISAMAX = 1
        !DO  I=1,N
        !  IF(ABS(Y(I)).GT.ABS(Y(ISAMAX))) ISAMAX=I
        !ENDDO


	  !FUNNORM=ABS(Y(ISAMAX))

	  sum2=0.0

!$OMP PARALLEL DO 
!$omp& shared(Y), private(i), reduction(+: sum2) 
        DO  I=1,N
          sum2=sum2 + Y(i)*Y(i)
        END DO
!$OMP END PARALLEL DO

	  FUNNORM=sqrt(sum2)


	RETURN
	END

      

SUBROUTINE LIN_BCG(NP,IA,JA,A_SPA,AD,B,X,unit_cont)
use def_constantes
integer :: np,unit_cont
double precision :: a_spa(*),AD(np),X(np),B(np)
integer :: IA(np+1),JA(*)
!local
DOUBLE PRECISION ::bnrm, Znrm, bknum, bkden, XANT, EPS,akden,bk,ak,err,zm1nrm,dxnrm,xnrm,FUNNORM
DOUBLE PRECISION, ALLOCATABLE ::Z(:),P(:),PP(:),R(:),ZZ(:),RR(:)



aLLOCATE (P(NP),PP(NP),ZZ(NP),Z(NP),R(NP),RR(NP)) 

iter=0
EPS= 1.d-14
err = 1.0
      
call MATXVECSIM(IA,JA,A_SPA,AD,X,R,NP)

	DO KK=1,NP
      R(KK)  = B(KK) - R(KK)
	  RR(KK) = R(KK)   
	ENDDO

      CALL MATXVECSIM(IA,JA,A_SPA,AD,R,RR,NP)


 	DO KK=1,NP
	  Z(KK)=B(KK)/AD(KK)
      if(abs(Z(kk))>1.0e-9) then
         write(6,*) kk,zz(kk),b(kk)/ad(kk)
      endif
    ENDDO

      bnrm=FUNNORM(NP,Z)

	DO KK=1,NP
	  Z(KK)=R(KK)/AD(KK)
      ENDDO

      ZNRM=FUNNORM(NP,Z)

      do while(err.gt.toler .and. iter.lt.itmax)

        iter=iter+1
      
 	  DO KK=1,NP
	    ZZ(KK)=RR(KK)/AD(KK)
        ENDDO
   	  
	  bknum=0.d0
  	  do j=1,NP
          bknum=bknum + Z(j)*RR(j)
        enddo

        if(iter.eq.1) then
          do  j=1,NP
             P (j)=Z(j)
             PP(j)=ZZ(j)
          enddo
        else
          bk=bknum/bkden
          do  j=1,NP
            P(j) = bk*P(j)+Z(j)
            PP(j)= bk*PP(j)+ZZ(j)
          enddo
        endif

        bkden=bknum

	  CALL MATXVECSIM(IA,JA,A_SPA,AD,P,Z,NP)
      
    	  akden=0.d0
        do  j=1,NP
          akden = akden+ z(j)*pp(j)
        enddo

        ak=bknum/akden
          
        CALL MATXVECSIM(IA,JA,A_SPA,AD,PP,ZZ,NP)

        DO J=1,NP
          x(j) = x(j)  + ak*p(j)
          r(j) = r(j)  - ak*z(j)
          RR(j)= RR(j) - ak*zz(j)
        ENDDO

     	  DO KK=1,NP
	    Z(KK)=R(KK)/AD(KK)
        ENDDO

        zm1nrm=znrm
        znrm=FUNNORM(NP,Z)
        
	  if(abs(zm1nrm-znrm).gt.EPS*znrm) then
          
		 dxnrm=abs(ak)*FUNNORM(NP,P)
         err=znrm/abs(zm1nrm-znrm)*dxnrm
                  
      	 xnrm=FUNNORM(NP,X)
          
	     if(err.le.0.5d0*xnrm) then
             err=err/xnrm
         else
             err=znrm/bnrm
         endif

	  else
        
	    err=znrm/bnrm
       
	  endif

    ENDDO

    WRITE(unit_cont,*) 'iteraciones internas ', ITER,err  

dEALLOCATE (P,PP,ZZ,Z,R,RR) 
	
end subroutine LIN_BCG

!DOUBLE PRECISION FUNCTION FUNNORM(N,Y) 
!integer :: n
!DOUBLE PRECISION :: Y(N)
! local
!integer :: isamax,i

!  isamax=1
!  do  i=1,N
!      if(abs(Y(i)).gt.abs(Y(isamax))) isamax=i
!  enddo
        
!  FUNNORM=abs(Y(isamax))

!end function FUNNORM
	
SUBROUTINE MATXVECSIM(IA,JA,AN,AD,B,C,Np)
implicit none
integer :: np
INTEGER :: IA(np+1),JA(*)
DOUBLE PRECISION :: B(np),AN(*),AD(np),C(np)
! local
integer :: k,i,iaf,iai,j

DO K=1,np
     C(K)= AD(K)*B(K)
enddo
      
DO I=1,np
     IAI = IA(I)
     IAF = IA(I+1)- 1


     IF(IAF.GE.IAI) THEN
          DO J=IAI,IAF
            C(I) = C(I) + AN(J)*B(JA(J))
            C(JA(J))=   C(JA(J)) + AN(J)*B(I)
          ENDDO
      ENDIF
enddo
     
end subroutine MATXVECSIM
      


SUBROUTINE CGBesta(NP,IA,JA,A_SPA,AD,B,X,TOL,ITMAX,ITER,ERR)
      implicit none 
	DOUBLE PRECISION A_SPA(*),AD(*),X(*),B(*)
      INTEGER IA(*),JA(*)
      DOUBLE PRECISION bnrm, Znrm, bknum, bkden, XANT, akden,bk,ak,zm1nrm,dxnrm,xnrm,alfa,b_norma,r_norma,aknum
	DOUBLE PRECISION TOL,EPS,ERR,FUNNORM
	INTEGER ITMAX,ITER,NP,j,kk

    
    DOUBLE PRECISION, ALLOCATABLE :: F0(:),U0(:),V0(:),Z0(:),q0(:),RR(:),R0(:),T0(:),S0(:)
    DOUBLE PRECISION :: sigma0,rho0,alfa0,omega0,sigma_1,pi0,fhi0,tau0,rho1,sigma1,delta,beta,gamma,nu,tita,kapa

    
    ALLOCATE (F0(NP),U0(NP),V0(NP),Z0(NP),q0(NP),R0(NP),RR(NP),T0(NP),S0(NP))
      

	CALL MATXVECSIM(IA,JA,A_SPA,AD,X,R0,NP)

    DO KK=1,NP
        R0(KK)  = B(KK) - R0(KK)
	ENDDO



	CALL MATXVECSIM(IA,JA,A_SPA,AD,R0,U0,NP)

	CALL MATXVECSIM(IA,JA,A_SPA,AD,R0,F0,NP)

    Z0=0.0
    q0=0.0
    V0=0.0

    ! CALL TRANSPOS(IA,JA,IAT,JAT,NP,NP) en ppio es simetrica

    
    b_norma = 0.0
	sigma0=0.0
    DO KK=1,NP
        sigma0=sigma0 + R0(kk)*U0(kk)
        RR(kk)=R0(kk)
        b_norma =b_norma  + b(kk)*b(kk) 

	ENDDO
    b_norma=sqrt(b_norma) 
      

	!DO KK=1,NP
	!    Z(KK)=R(KK)/AD(KK)
    !    rr(kk)=z(kk)
	!ENDDO

    !ZNRM=FUNNORM(NP,Z)
     
     rho0=1.0
     alfa0=1.0
     omega0=1.0
     sigma_1=0.0
     pi0=0.0
     fhi0=0.0
     tau0=0.0
     
      
    do while(err.gt.tol .and. iter.lt.itmax)

      iter=iter+1

      rho1= fhi0 - omega0*sigma_1 + omega0*alfa0*pi0    !     p = phi - w * sigmam1 + w * old_alfa * pi;

		!		del = (curr_it == 1) ? p : p/old_tau;
		

		

		

!      delta = rho1/(rho0)*alfa0

      delta=rho1
      if(iter>1) delta =rho1/tau0
      rho0=rho1

      beta = delta/omega0 !		B = del / w;

!      tau0= sigma0 + beta*tau0 - delta*pi0 !		tau = sigma + B * (old_tau - w * pi);
      tau0= sigma0 + beta*tau0 - beta*omega0*pi0  !		tau = sigma + B * (old_tau - w * pi);
     
      alfa0 = rho0/tau0  !		alfa = p / tau;

     
      do j=1,NP
          V0(j)=U0(j) + beta * V0(j) - delta*q0(j) 
      enddo

      CALL MATXVECSIM(IA,JA,A_SPA,AD,V0,q0,NP)

      do j=1,NP
          s0(j)=RR(j) -alfa0 * V0(j)  
          t0(j)=U0(j) -alfa0 * q0(j)  
      enddo

      do j=1,NP
          z0(j)=alfa0*RR(j) + beta * Z0(j) - alfa0*delta*V0(j) 
      enddo

      fhi0=0.0
      pi0=0.0
      gamma =0.0
      nu = 0.0
      tita=0.0
      kapa=0.0
      DO KK=1,NP
        fhi0=fhi0 + R0(kk)*S0(kk)
        pi0=pi0 + R0(kk)*q0(kk)
	    gamma = gamma+f0(kk)*s0(kk)
	    nu = nu + f0(kk)*t0(kk)
	    tita = tita + S0(kk)*t0(kk)
	    kapa = kapa + t0(kk)*t0(kk)
      ENDDO

      omega0 = tita/kapa
      sigma_1=sigma0
      sigma0=gamma-omega0*nu

      DO KK=1,NP
         RR(kk) = S0(kk) - omega0*T0(kk)
         X(kk) = X(kk) + Z0(kk)+omega0*S0(kk)
      enddo


	  
      CALL MATXVECSIM(IA,JA,A_SPA,AD,RR,U0,NP)
        
        
	  
        r_norma = 0.0
	  DO KK=1,NP
          r_norma =r_norma  + x(kk)*x(kk)    
	  ENDDO
        r_norma=sqrt(r_norma) 

        err=r_norma/b_norma
        !write(2,*) iter,err,r_norma,ak,alfa

      ENDDO


	
    end subroutine cgBesta
                         
 
    SUBROUTINE LIN_CG_cplx(NP,IA,JA,A_SPA,AD,B,X,TOL,ITMAX,ITER,ERR)
      implicit none
	COMPLEX*16:: A_SPA(*),AD(*),X(*),B(*)
      INTEGER IA(*),JA(*)
	INTEGER  np,kk,j
      complex*16 bnrm,  bknum, bkden, XANT, akden,bk,ak,zm1nrm,dxnrm,xnrm,aknum,alfa,b_norma
	DOUBLE PRECISION TOL,ERR,FUNNORM_cplx,ZNRM,r_norma
	INTEGER ITMAX,ITER

      COMPLEX*16, ALLOCATABLE :: Z(:),P(:),PP(:),R(:),ZZ(:),RR(:)

	ALLOCATE (P(NP),PP(NP),ZZ(NP),Z(NP),R(NP),RR(NP)) 

       

      ITER=  0
	ERR = 1.0
      
	CALL MATXVECSIM_cplx(IA,JA,A_SPA,AD,X,R,NP)

      b_norma = 0.0
	DO KK=1,NP
        R(KK)  = B(KK) - R(KK)
!        b_norma =b_norma  + b(kk)*conjg(b(kk))
        b_norma =b_norma  + b(kk)*(b(kk))
	ENDDO
      b_norma=sqrt(b_norma) 

	DO KK=1,NP
	  Z(KK)=R(KK)/AD(KK)
        rr(kk)=z(kk)
	ENDDO

      ZNRM=FUNNORM_cplx(NP,Z)
      
      do while(err.gt.tol .and. iter.lt.itmax)

        iter=iter+1
        
	  bknum=0.d0
  	  do j=1,NP
!          bknum=bknum + conjg(Z(j))*R(j)
          bknum=bknum + (Z(j))*R(j)
        enddo

	  CALL MATXVECSIM_cplx(IA,JA,A_SPA,AD,rr,Z,NP)
        
	  bkden=0.0
	  do j=1,NP
!          bkden=bkden + conjg(rr(j))*Z(j)
          bkden=bkden + (rr(j))*Z(j)
        enddo
	   
        alfa=bknum/bkden

        do kk=1,np
          x(kk)=x(kk)+alfa*rr(kk)          
	    r(kk)= r(kk) - alfa * z(kk)
	  enddo

     	  DO KK=1,NP
	    Z(KK)=R(KK)/AD(KK)
        ENDDO
        
	  
    	  akden=bknum
	  aknum=0.0
        do  j=1,NP
!          aknum = aknum+ conjg(z(j))*r(j)
          aknum = aknum+ (z(j))*r(j)
        enddo

        ak=aknum/akden
          
        DO J=1,NP
          rr(j) = z(j)  + ak*rr(j)
        ENDDO

        r_norma = 0.0
	  DO KK=1,NP
          r_norma =r_norma  + r(kk)*conjg(r(kk))    
	  ENDDO
        r_norma=sqrt(r_norma) 

        err=abs(r_norma/b_norma)
        write(2,'(i3,7e15.5)') iter,err,r_norma,ak,alfa

      ENDDO


	DEALLOCATE (P,PP,ZZ,Z,R,RR) 
	
	RETURN
      END


      DOUBLE PRECISION FUNCTION FUNNORM_cplx(N,Y) 
      complex*16 ::  Y(*)

        ISAMAX = 1
        DO  I=1,N
          IF(CDABS(Y(I)).GT.CDABS(Y(ISAMAX))) ISAMAX=I
        ENDDO
        
	  FUNNORM_cplx=CDABS(Y(ISAMAX))

	RETURN
	END

      SUBROUTINE MATXVECSIM_cplx(IA,JA,AN,AD,B,C,N)
      INTEGER  IA(*),JA(*)
      COMPLEX*16 :: B(*),AN(*),AD(*),C(*)

      DO K=1,N
        C(K)= AD(K)*B(K)
!c        write(2,*) 'inicio ', k,c(k)
	ENDDO
      
      DO I=1,N
        IAI = IA(I)
        IAF = IA(I+1)- 1

!c        write(2,*) 'i, iai, iaf ', i,iai,iaf

        IF(IAF.GE.IAI) THEN
          DO J=IAI,IAF
            C(I) = C(I) + AN(J)*B(JA(J))
            C(JA(J))=   C(JA(J)) + AN(J)*B(I)
!c            write(2,*) 'i, j, c ', i,ja(j),c(i),C(JA(J))
          ENDDO
        ENDIF
      ENDDO
      
      RETURN
      END




      

   