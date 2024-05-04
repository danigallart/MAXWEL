subroutine solver
use def_io
use def_variables
use def_vectors
implicit none

integer :: ii,i,j
integer :: kk,iter
double precision :: err
complex*16, allocatable :: Au_inc(:)

allocate(u_inc(NP),u_scat(NP),u_tot(NP),Au_inc(NP))

u_inc = exp(ij*k0*(real(complex_coorx)*cos(phii)+real(complex_coory)*sin(phii))) !No need to use complex coordinates. We just need the real value at each node


call MATXVECSIM_cplx(NP,IA,JA,AN,AD,u_inc,Au_inc)


indep_vect=cmplx(0.0,0.0)

do ii=1,n_scatin
    i = scatin_nodes(ii)
    indep_vect(i) = -Au_inc(i)
enddo

u_scat = cmplx(0.0,0.0)

iter=1
err=10.0

!call CG(NP,IA,JA,AN,AD,indep_vect,u_scat,tol_solver,iter_solver,iter,ERR)
call LIN_CG_cplx(NP, IA, JA, AN, AD, indep_vect, u_scat, tol_solver, err, iter_solver, iter)

u_tot = u_scat+u_inc

do ii=1,n_pml
    i = pml_nodes(ii)
    u_tot(i) = cmplx(0.0,0.0)
enddo


write(control_unit,*) 'CG it.   ',iter, '/ ', iter_solver
write(control_unit,*) 'CG err.  ',ERR

end subroutine solver

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     SUBROUTINE LIN_CG_cplx(NP,IA,JA,A_SPA,AD,B,X,TOL,ERR,ITMAX,ITER)
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
      
	CALL MATXVECSIM_cplx(NP,IA,JA,A_SPA,AD,X,R)

      b_norma = 0.0
	DO KK=1,NP
        R(KK)  = B(KK) - R(KK)
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
          bknum=bknum + (Z(j))*R(j)
        enddo

	  CALL MATXVECSIM_cplx(NP,IA,JA,A_SPA,AD,rr,Z)
        
	  bkden=0.0
	  do j=1,NP
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

      SUBROUTINE MATXVECSIM_cplx(N,IA,JA,AN,AD,B,C)
      INTEGER  IA(*),JA(*)
      COMPLEX*16 :: B(*),AN(*),AD(*),C(*)

      DO K=1,N
        C(K)= AD(K)*B(K)
	ENDDO
      
      DO I=1,N
        IAI = IA(I)
        IAF = IA(I+1)- 1


        IF(IAF.GE.IAI) THEN
          DO J=IAI,IAF
            C(I) = C(I) + AN(J)*B(JA(J))
          ENDDO
        ENDIF
      ENDDO
      
      RETURN
      END