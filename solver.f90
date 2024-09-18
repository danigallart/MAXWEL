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

u_inc = exp(ij*(k0*(real(complex_coorx)*cos(phii)+real(complex_coory)*sin(phii))))!/nu0 !No need to use complex coordinates. We just need the real value at each node

!if (system_sym=='Y') then
    
!  call MATXVECSIM_cplx(NP,IA,JA,AN,AD,u_inc,Au_inc)
  
!elseif (system_sym=='N') then
    
!  call MATXVEC_cplex(IA,JA,AN,AD,u_inc,Au_inc,NP,NONULL,0)
  
!endif

!indep_vect=cmplx(0.0,0.0)

!do kk=1,NE
!    if (material(kk) == 1) then
!        do ii = 1, nodpel 
!            i = conn(kk,ii)
!            indep_vect(i) = -Au_inc(i)
!        enddo
!    endif
!enddo

u_scat = cmplx(0.0,0.0)

iter=1
err=10.0

if (system_sym == 'Y') then

    call LIN_CG_cplx(NP, IA, JA, AN, AD, indep_vect, u_scat, tol_solver, err, iter_solver, iter)
    
elseif (system_sym == 'N') then
    
    CALL linbcg(NP,NONULL,indep_vect,u_scat,tol_solver,iter_solver,iter,err,IA,JA,AD,AN)

endif


u_tot = u_scat+u_inc

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
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE linbcg(N,N_NUL,B,X,tol,itmax,iter,err,IA,JA,AD,SA)
      IMPLICIT none
      integer :: n, n_nul, itmax, maxis,iter
      real*8:: tol, err, eps, FUNNORM_cplx, bnrm, xnrm 
      complex*16 :: b(N),x(N),SA(N_NUL),ad(n)
      integer :: IA(n+1),ja(n_nul)
       COMPLEX*16, ALLOCATABLE :: Z(:),P(:),PP(:),R(:),ZZ(:),RR(:)
      integer:: j,i
      complex*16:: znrm, snrm, bknum, bk, bkden, zm1nrm, ak, akden, aknum,  dxnrm

      
      ALLOCATE (P(N),PP(N),ZZ(N),Z(N),R(N),RR(N)) 
      
      iter=0
      EPS=1E-14 
      
      DO I=1,N
        X(I)=0.0
      ENDDO
      
      call MATXVEC_cplex(IA,JA,SA,AD,X,r,n,n_nul,0)
        
      do  j=1,n
        r(j)=b(j)-r(j)
        rr(j)=r(j)
      enddo 
      call MATXVEC_cplex(IA,JA,SA,AD,r,rr,n,n_nul,0)
             
      bnrm= FUNNORM_cplx(N,b)
     
        do  i=1,n
          z(i)=r(i)/AD(i) 
      enddo 


  100 if (iter.le.itmax) then 
    !   Main loop.
      iter=iter+1
       do  i=1,n
          zz(i)=rr(i)/AD(i) 
       enddo 
      
!      Final 1 indicates us,IJA,SAe of transpose matrix AT.
      bknum=0.d0
      do j=1,n 
!      Calculate coeficient bk and direction vectors p and pp 
         bknum=bknum+z(j)*rr(j)
      enddo 
   
      if(iter.eq.1) then
          do  j=1,n
           p(j)=z(j)
           pp(j)=zz(j)
          enddo 
      else 
 
          bk=bknum/bkden
          do  j=1,n
           p(j)=bk*p(j)+z(j)
           pp(j)=bk*pp(j)+zz(j)
          enddo 
      endif
      
      bkden=bknum 
!      Calculate coeficient ak new iterate x and new residuals r and rr 
      call MATXVEC_cplex(IA,JA,SA,AD,p,z,n,n_nul,0)
       
      akden=0.d0
      
          do  j=1,n
           akden=akden+z(j)*pp(j)
          enddo 
          ak=bknum/akden
          call MATXVEC_cplex(IA,JA,SA,AD,pp,zz,n,n_nul,1)
      
          
          do  j=1,n
           x(j)=x(j)+ak*p(j)
           r(j)=r(j)-ak*z(j)
           rr(j)=rr(j)-ak*zz(j)
          enddo 
          
          do  i=1,n
             z(i)=r(i)/AD(i) 
          enddo 

          
!   Solve A ·z = r and check stopping criterion.

         err=  FUNNORM_cplx(N,r)/bnrm

      if(err.gt.tol) goto 100
      
      endif
      
      
      deALLOCATE (P,PP,ZZ,Z,R,RR) 
      
      
      return
    END
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
      DOUBLE PRECISION FUNCTION FUNNORM_cplx(N,Y) 
      complex*16 ::  Y(*)

        ISAMAX = 1
        DO  I=1,N
          IF(CDABS(Y(I)).GT.CDABS(Y(ISAMAX))) ISAMAX=I
        ENDDO
        
	  FUNNORM_cplx=CDABS(Y(ISAMAX))

	RETURN
    END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
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
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
SUBROUTINE MATXVEC_cplex(IA,JA,AN,AD,B,C,N,n_nul,caso)
      implicit none 
      INTEGER  N,k,i,j,iai,iaf,ii,n_nul,caso
	  INTEGER  IA(N+1),JA(n_nul)
      COMPLEX*16 :: B(n),AN(n_nul),AD(n),C(n)
      COMPLEX*16 :: sum
      INTEGER :: IAI_t(N), IAf_t(N),index1

      
      
      DO K=1,N
        C(K)= AD(K)*B(K)
        IAI_t(k) = IA(k)
        IAF_t(k) = IA(k+1) - 1
      ENDDO

      if(caso==0) then
          DO I=1,N
              sum = C(i)
              DO J=IAI_t(i),IAF_t(i)
                  k= JA(J)
                  sum = sum + AN(j)*B(k) 
              ENDDO
              C(i)=sum
          ENDDO
      else
          DO I=1,N
              DO J=IAI_t(i),IAF_t(i)
                 k= JA(J)
                 c(k) = c(k) + AN(J)*B(i) 
              ENDDO
          ENDDO
         
          
      endif
      
          
          
      RETURN
      END
