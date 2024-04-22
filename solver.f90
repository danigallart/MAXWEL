subroutine solver
use def_io
use def_variables
use def_vectors
implicit none

integer :: ii,i,j
integer :: kk,iter
double precision :: err
complex, allocatable :: Au_inc(:)

allocate(u_inc(NP),u_scat(NP),Au_inc(NP))

u_inc = exp(ij*k0*(real(complex_coorx)*cos(phii)+real(complex_coory)*sin(phii))) !No need to use complex coordinates. We just need the real value at each node


call csrmult_complex(NP,IA,JA,AN,AD,u_inc,Au_inc)


indep_vect=cmplx(0.0,0.0)

do ii=1,n_scatin
    i = scatin_nodes(ii)
    indep_vect(i) = -Au_inc(i)
enddo

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
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
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
