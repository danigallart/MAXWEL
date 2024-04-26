subroutine solver
use def_io
use def_variables
use def_vectors
implicit none

integer :: ii,i,j
integer :: kk,iter
double precision :: err
complex*16, allocatable :: Au_inc(:)

allocate(u_inc(NP),u_scat(NP),Au_inc(NP))

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

write(control_unit,*) 'CG it.   ',iter, '/ ', iter_solver
write(control_unit,*) 'CG err.  ',ERR

end subroutine solver

!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine cg_solve_complex(n, ia, ja, an, ad, b, x, tol, err, max_iter, iter)
    
    implicit none    

    !Input variables
    integer :: n
    integer :: ia(*),ja(*)
    complex :: ad(*),an(*),b(*)
    integer :: max_iter, iter
    double precision :: tol, err 
    
    !Output variables
    complex :: x(*) ! solution vector
    
    !Local variables
    complex, allocatable :: r(:), p(:), ap(:), z(:)
    integer :: i
    double precision :: norm_r, norm_b, norm_z, FUNNORM
    complex :: alpha, den_alpha
    double precision :: beta, num_alpha, num_beta, den_beta

    
    allocate(r(n), p(n), ap(n), z(n))

    ! initial residual
    call csrmult_complex(n,ia,ja,an,ad,x,r)
    
    do i=1,n
        r(i) = b(i) - r(i)
        z(i) = r(i)/ad(i)
        p(i) = z(i)
    enddo
    
    norm_b = FUNNORM(n,b)
        
    do while (err>tol .and. iter<max_iter)
        
        
        num_alpha = 0.0
        do i=1,n
            num_alpha = num_alpha + z(i)*conjg(r(i))
        enddo
        
        ! matrix-vector multiplication (A * p)
        call csrmult_complex(n,ia,ja,an,ad,p,ap)
        
        den_alpha = cmplx(0.0,0.0)
        
        do i=1,n
            den_alpha=den_alpha+ap(i)*conjg(p(i))
        enddo
        
        alpha = num_alpha/den_alpha
        
        do i=1,n
            x(i) = x(i)+alpha*p(i)
            r(i) = r(i)-alpha*ap(i)
            z(i) = r(i)/ad(i)
        enddo
        
        den_beta = num_alpha
        num_beta = 0.0
        
        do i=1,n
            num_beta = num_beta + z(i)*conjg(r(i))
        enddo
        
        beta = num_beta/den_beta
        
        do i=1,n
            p(i) = z(i) + beta*p(i)
        enddo
        
        norm_r = FUNNORM(n,r)
        
        err = norm_r/norm_b
        
        iter=iter+1

    enddo
    
    deallocate(r,p,ap,z)

end subroutine cg_solve_complex

    
subroutine csrmult_complex(n, ia, ja, an, ad, x, y)
implicit none
    !Input variables
  integer :: n
  integer :: ia(*),ja(*)
  complex :: ad(*),an(*)
  complex :: x(*)
  
  !Output variables
  complex :: y(*)

  !Local variables
  integer :: i, j

  do i=1,n
      !diagonal elements
      y(i) = ad(i)*x(i)
  enddo
  

  do i = 1, n
    ! loop through rows
    do j = ia(i), ia(i+1)-1 ! loop through non-zero elements in row i
        !non-diagonal elements
        y(i) = y(i) + an(j) * x(ja(j))
    enddo
  enddo

end subroutine csrmult_complex
        
DOUBLE PRECISION FUNCTION FUNNORM(N,Y)
implicit none
COMPLEX Y(*)
double precision sum2
integer N,i

	  sum2=0.0

      DO  I=1,N
          sum2=sum2 + Y(i)*conjg(Y(i))
      END DO
	  FUNNORM=sqrt(sum2)


RETURN
    END

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
        b_norma =b_norma  + b(kk)*conjg(b(kk))
!        b_norma =b_norma  + b(kk)*(b(kk))
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
          bknum=bknum + conjg(Z(j))*R(j)
!          bknum=bknum + (Z(j))*R(j)
        enddo

	  CALL MATXVECSIM_cplx(NP,IA,JA,A_SPA,AD,rr,Z)
        
	  bkden=0.0
	  do j=1,NP
          bkden=bkden + conjg(rr(j))*Z(j)
!          bkden=bkden + (rr(j))*Z(j)
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
          aknum = aknum+ conjg(z(j))*r(j)
!          aknum = aknum+ (z(j))*r(j)
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
!c        write(2,*) 'inicio ', k,c(k)
	ENDDO
      
      DO I=1,N
        IAI = IA(I)
        IAF = IA(I+1)- 1

!c        write(2,*) 'i, iai, iaf ', i,iai,iaf

        IF(IAF.GE.IAI) THEN
          DO J=IAI,IAF
            C(I) = C(I) + AN(J)*B(JA(J))
!            C(JA(J))=   C(JA(J)) + AN(J)*B(I)
!c            write(2,*) 'i, j, c ', i,ja(j),c(i),C(JA(J))
          ENDDO
        ENDIF
      ENDDO
      
      RETURN
      END