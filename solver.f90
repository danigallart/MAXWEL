subroutine solver
use def_io
use def_variables
use def_vectors
implicit none

integer :: i,j
integer :: kk,iter
double precision :: err

allocate(u_inc(NP),u_scat(NP))

u_inc = exp(ij*k0*(coorx+coory)) !No need to use complex coordinates. We just need the real value at each node

do i=1,np
    do j=IA(i),IA(i+1)-1
        if (JA(j)==i) then
            indep_vect(i) = indep_vect(i) + AD(i) * u_inc(i)
        else
            indep_vect(JA(j)) = indep_vect(JA(j)) + AN(j) * u_inc(i)
        endif
        
    enddo
enddo

indep_vect = -indep_vect

u_scat = cmplx(0.0,0.0)

iter=1
err=10.0

!call CG(NP,IA,JA,AN,AD,indep_vect,u_scat,tol_solver,iter_solver,iter,ERR)
call cg_solve_complex(NP, IA, JA, AN, AD, indep_vect, u_scat, tol_solver, err, iter_solver, iter) 

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
    complex, allocatable :: r(:), p(:), ap(:)
    integer :: i
    double precision :: norm_r, norm_b, alpha, beta
    
    allocate(r(n),p(n),ap(n))

    ! initial residual
    call MATXVEC_t(ia,ja,an,ad,x,r,n)
    
    do i=1,n
        r(i) = b(i) - r(i)
    enddo

    ! initial direction
    p = r

    do while (err>tol .and. iter<max_iter)
        
        iter=iter+1
        ! matrix-vector multiplication (A * p)
        call MATXVEC_t(ia,ja,an,ad,p,ap,n)

        ! alpha (dot product) with conjg for complex conjugate
        alpha = dot_product(conjg(r), r) / dot_product(p, ap)

        ! update solution
        do i=1,n
            x(i) = x(i) + alpha * p(i)
        enddo
        ! update residual
        r = r - alpha * ap 

        ! check convergence
        norm_r=0.0
        norm_b=0.0
        
        do i=1,n
        norm_r = norm_r + r(i)*conjg(r(i))
        norm_b = norm_b + b(i)*conjg(b(i))
        enddo
        
        err = sqrt(norm_r/norm_b)
        
        if (err < tol) then
        exit
        endif

        ! beta
        beta = dot_product(conjg(r), r) / dot_product(p, ap)

        ! update direction
        p = r + beta * p
    enddo


end subroutine cg_solve_complex

    
subroutine csrmult_complex(n, ia, ja, ad, an, x, y)

    !Input variables
  integer :: n
  integer :: ia(*),ja(*)
  complex :: ad(*),an(*)
  complex :: x(*)
  
  !Output variables
  complex :: y(*)

  !Local variables
  integer :: i, j

  do i = 1, n
    ! loop through rows
    do j = ia(i), ia(i+1)-1 ! loop through non-zero elements in row i
      if (ja(j) == i) then
        ! diagonal element
        y(i) = y(i) + ad(i) * x(i)
      else
        ! non-diagonal element
        y(ja(j)) = y(ja(j)) + an(j) * x(i)
      endif
    enddo
  enddo

end subroutine csrmult_complex
    
    
    
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE CG(NP,IA,JA,A_SPA,AD,B,X,TOL,ITMAX,ITER,ERR)
implicit none 
COMPLEX A_SPA(*),AD(*),X(*),B(*)
INTEGER IA(*),JA(*)
COMPLEX bknum, bkden, XANT, akden,bk, ak,zm1nrm,dxnrm,xnrm,alfa,aknum
DOUBLE PRECISION TOL,EPS,ERR, b_norm, r_norm, z_norm
DOUBLE PRECISION FUNNORM
INTEGER ITMAX,ITER,NP,j,kk

COMPLEX, ALLOCATABLE ::Z(:),P(:),PP(:),R(:),ZZ(:),RR(:)

ALLOCATE (P(NP),PP(NP),ZZ(NP),Z(NP),R(NP),RR(NP)) 
      
CALL MATXVEC_t(IA,JA,A_SPA,AD,X,R,NP)

b_norm = 0.0
DO KK=1,NP
    R(KK)  = B(KK) - R(KK)
    b_norm = b_norm  + b(kk)*conjg(b(kk)) 
ENDDO
b_norm=sqrt(b_norm) 
      

DO KK=1,NP
   Z(KK)=R(KK)/AD(KK)
   rr(kk)=z(kk)
ENDDO

z_norm=FUNNORM(NP,Z)
    
do while(err>tol .and. iter<itmax)

        iter=iter+1

	  bknum=cmplx(0.0,0.0)
  	  do j=1,NP
          bknum=bknum + Z(j)*R(j)
      enddo
	  CALL MATXVEC_t(IA,JA,A_SPA,AD,rr,Z,NP)
        
	  bkden=cmplx(0.0,0.0)
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
	  aknum=cmplx(0.0,0.0)
      do  j=1,NP
          aknum = aknum+ z(j)*r(j)
      enddo

      ak=aknum/akden
          
      DO J=1,NP
          rr(j) = z(j)  + ak*rr(j)
      ENDDO

      r_norm = cmplx(0.0,0.0)
	  DO KK=1,NP
          r_norm =r_norm  + r(kk)*r(kk)    
	  ENDDO
      r_norm=sqrt(r_norm) 

      err=r_norm/b_norm

ENDDO


DEALLOCATE (P,PP,ZZ,Z,R,RR) 
	
    END subroutine cg
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
SUBROUTINE MATXVEC_t(IA,JA,AN,AD,B,C,N)
implicit none 
INTEGER  N,k,i,j,iai,iaf,ii
INTEGER  IA(*),JA(*)
COMPLEX B(*),AN(*),AD(*),C(*)
!DOUBLE PRECISION :: sum,fac
INTEGER :: IAI_t(N), IAf_t(N),index1

 DO K=1,N
        C(K)= AD(K)*B(K)
        IAI_t(k) = IA(k)
        IAF_t(k) = IA(k+1)- 1
 ENDDO


      DO I=1,N
          DO J=IAI_t(i),IAF_t(i)
             C(i) = C(i) + AN(J)*B(JA(J)) 
          ENDDO
      ENDDO
RETURN
END


DOUBLE PRECISION FUNCTION FUNNORM(N,Y)
implicit none
COMPLEX Y(*),sum2
 integer N,isamax,i

	  sum2=CMPLX(0.0,0.0)

      DO  I=1,N
          sum2=sum2 + Y(i)*conjg(Y(i))
      END DO
	  FUNNORM=sqrt(sum2)


RETURN
END
