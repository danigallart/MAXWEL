SUBROUTINE assembly
use def_io
use def_variables
use def_vectors

implicit none

!DOUBLE PRECISION :: rigidez_local(nodpel,nodpel),indep_local(nodpel),x(ndimension,nodpel),dx(ndimension), QT
!INTEGER :: kk,ii,ipoin,inode,jnode,j,ij,numbergaus,ii,jj,jj2,IAUX,keje,ns(nodpel), vecFlujoElem(nodpel)
!DOUBLE PRECISION :: adiag,fuente, funQT, solucion_local(nodpel), masa_local(nodpel,nodpel), cp, rho
!DOUBLE PRECISION :: funfuente, funconductividad, funcp, funrho

!integer :: kk,ii,i,j,jj,ns(nodpel)
!double precision, dimension(NE) :: rel_permitivity
!double precision :: pmlbin_coorx(n_pml_bin),pmlbin_coory(n_pml_bin), &
!    pmlbout_coorx(n_pml_bout), pmlbout_coory(n_pml_bout)
!double precision :: x_rval, y_rval, x_cval, y_cval
!complex*16 :: coord_matrix(ndim,nodpel)

integer :: kk, ii, i, j, jj, KEJE, IAUX
integer, allocatable :: ns(:)
double precision :: rel_permitivity_xx, rel_permitivity_xy, &
                                 rel_permitivity_yx, rel_permitivity_yy, &
                                 rel_permitivity_zz

double precision :: rel_permeability_xx, rel_permeability_xy, &
                                 rel_permeability_yx, rel_permeability_yy, &
                                 rel_permeability_zz

double precision, allocatable :: cond

double precision, allocatable :: pmlbin_coorx(:),pmlbin_coory(:), pmlbout_coorx(:), pmlbout_coory(:)
complex*16, allocatable :: local_coords(:,:)

double precision :: x_rval, y_rval, x_cval, y_cval

complex*16, allocatable :: JACOB(:,:),INVJACOB(:,:)
double precision, allocatable :: PHI(:,:),DPHI(:,:)
complex*16 :: DETJACOB, determinant

complex*16, allocatable :: AE(:,:)
complex*16 :: pxxe,pxye,pyxe,pyye,qe
complex*16, allocatable :: pe(:,:)


allocate(pmlbin_coorx(n_pml_bin), pmlbin_coory(n_pml_bin))
allocate(pmlbout_coorx(n_pml_bout), pmlbout_coory(n_pml_bout))
allocate(ns(nodpel))
allocate(local_coords(ndim,nodpel))
allocate(AE(nodpel,nodpel))
allocate(pe(ndim,ndim))

allocate(JACOB(ndim,ndim),INVJACOB(ndim,ndim))
allocate(PHI(Ngauss,nodpel), DPHI(ndim, nodpel))

Ngauss = 3

indep_vect=0.0
AD=0.0
AN=0.0


do ii=1,n_pml_bin
    pmlbin_coorx(ii)=coorx(pml_bin_nodes(ii))
    pmlbin_coory(ii)=coory(pml_bin_nodes(ii))
end do


do ii=1,n_pml_bout
    pmlbout_coorx(ii)=coorx(pml_bout_nodes(ii))
    pmlbout_coory(ii)=coory(pml_bout_nodes(ii))
end do

allocate(complex_coorx(NP), complex_coory(NP))

do ii=1,NP
    complex_coorx(ii)=cmplx(coorx(ii),0.0)
    complex_coory(ii)=cmplx(coory(ii),0.0)
end do

do ii=1,n_pml
    jj = pml_nodes(ii)
    call lcpml(coorx(jj),coory(jj),k0,pmlbin_coorx,pmlbin_coory,pmlbout_coorx,pmlbout_coory,complex_coorx(jj),complex_coory(jj))
    !complex_coorx(jj)=cmplx(x_rval,x_cval)
    !complex_coory(jj)=cmplx(y_rval,y_cval)
end do

do kk=1,NE
    
    rel_permeability_xx = 1.0
    rel_permeability_xy = 0.0
    rel_permeability_yx = 0.0
    rel_permeability_yy = 1.0
    rel_permeability_zz = 1.0
    
    
    if (material(kk) == 1) then
            rel_permitivity_xx=9.0
            rel_permitivity_yy=4.0
            rel_permitivity_zz=2.0
            rel_permitivity_xy=0.0
            rel_permitivity_yx=0.0
            
            cond = 0.0
            
    else
        
            rel_permitivity_xx=1.0
            rel_permitivity_yy=1.0
            rel_permitivity_zz=1.0
            rel_permitivity_xy=0.0
            rel_permitivity_yx=0.0
            
            cond = 0.0
    endif
    
    
    do i=1,nodpel
        ns(i) = conn(kk,i)
        local_coords(1,i) = complex_coorx(ns(i))
        local_coords(2,i) = complex_coory(ns(i))
    enddo
    call shape_gauss(local_coords(1,:),local_coords(2,:),PHI,DPHI,JACOB,INVJACOB,DETJACOB,Ngauss,nodpel,ndim)
    
    if (pol == 'TM') then
        determinant = cmplx(rel_permeability_xx,0.0) * cmplx(rel_permeability_yy,0.0) - cmplx(rel_permeability_xy,0.0) * cmplx(rel_permeability_yx,0.0)
        pxxe = cmplx(rel_permeability_xx,0.0)/determinant
        pyye = cmplx(rel_permeability_yy,0.0)/determinant
        pyxe = cmplx(rel_permeability_xy,0.0)/determinant
        pxye = cmplx(rel_permeability_yx,0.0)/determinant
        qe = -(rel_permitivity_zz-ij*cond/(omg*e0))*k0**2
    else if (pol == 'TE') then
        determinant = cmplx(rel_permitivity_xx,-cond/(omg*e0)) * cmplx(rel_permitivity_yy,-cond/(omg*e0)) - cmplx(rel_permitivity_xy,-cond/(omg*e0)) * cmplx(rel_permitivity_yx,-cond/(omg*e0))
        pxxe = cmplx(rel_permitivity_xx,-cond/(omg*e0))/determinant
        pyye = cmplx(rel_permitivity_yy,-cond/(omg*e0))/determinant
        pyxe = cmplx(rel_permitivity_xy,-cond/(omg*e0))/determinant
        pxye = cmplx(rel_permitivity_yx,-cond/(omg*e0))/determinant
        qe = -cmplx(rel_permeability_zz,0.0)*k0**2
    endif
    
    
    call element_matrix(PHI,DPHI,INVJACOB,DETJACOB, &
        pxxe,pyye,pxye,pyxe,qe, &
        AE,Ngauss,nodpel,ndim)
    
    
    do i=1, nodpel
    AD(ns(i)) = AD(ns(i)) + AE(i,i)
    do j=1, nodpel
        do IAUX = 1,ICX(ns(i))
            KEJE = IA(ns(i))+IAUX-1
            if (JA(KEJE) == ns(j)) then
                AN(KEJE) = AN(KEJE) + AE(i,j)
            endif
        enddo
        
    end do
    enddo
    
end do

! Compute elements





END SUBROUTINE assembly


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! LCPML function for Fortran 90
subroutine lcpml(x, y, k, pmlbin_x, pmlbin_y, pmlbout_x, pmlbout_y, xc, yc)
  use def_io
  use def_variables
  use def_vectors

  implicit none

  ! Input arguments
  real(kind=8) :: x, y, k
  real(kind=8), dimension(n_pml_bin) :: pmlbin_x, pmlbin_y
  real(kind=8), dimension(n_pml_bout) :: pmlbout_x, pmlbout_y
  ! Output arguments
  complex*16 :: xc, yc
  !double precision, intent(out) :: xc_r, yc_r,xc_im, yc_im

  ! LC-PML parameters
  !real(kind=8) :: alpha
  real(kind=8),dimension(n_pml_bin) :: dpml1
  integer :: m, ind
  complex*16 :: alphajk, term
  real(kind=8) :: x0,y0
  double precision :: vpx(n_pml_bout),vpy(n_pml_bout)
  double precision :: npx(n_pml_bout),npy(n_pml_bout)
  double precision :: lp(n_pml_bout)
  double precision :: vx, vy, l, nx, ny, x1, y1,dpml2
  double precision :: ksi
  

  ! Set LC-PML parameters
  !alpha = 7.0 * k
  !alphajk = -7.0
  alphajk = cmplx(0,-7.0)
  m = 3  ! PML decay rate (integer 2 or 3)

  ! Find the point on the inner PML boundary nearest to point P
  dpml1 = sqrt((pmlbin_x - x)**2 + (pmlbin_y - y)**2)
  ksi = MINVAL(dpml1,1)
  ind = MINLOC(dpml1,1)  ! Use minloc function for efficiency
  x0 = pmlbin_x(ind)
  y0 = pmlbin_y(ind)

  ! Find the point on the outer PML boundary in the direction of the unit vector
  vpx = pmlbout_x - x0
  vpy = pmlbout_y - y0
  lp = sqrt(vpx**2 + vpy**2)
  npx = vpx / lp  ! Unit vector from r0 to r1 (x-comp)
  npy = vpy / lp  ! Unit vector from r0 to r1 (y-comp)

  vx = x - x0
  vy = y - y0
  l = sqrt(vx**2 + vy**2)
  nx = vx / l  ! Unit vector from r0 to r (x-comp)
  ny = vy / l  ! Unit vector from r0 to r (y-comp)

  if (l < 1.0e-8) then
    xc = cmplx(x,0.0)
    yc = cmplx(y,0.0)
  else
    ! Find the angle between nx and npx using arccosine
    ind = MINLOC(acos(nx * npx + ny * npy),1)

    ! Find the point on the outer PML boundary closest in direction
    !ind = MINLOC(abs(acos(npx * cos(ksi) + npy * sin(ksi))),1)
    x1 = pmlbout_x(ind)
    y1 = pmlbout_y(ind)

    ! Calculate local PML thickness
    dpml2 = sqrt((x1 - x0)**2 + (y1 - y0)**2)

    ! Complex coordinate term
    term = alphajk * ((ksi**m) / (m * (dpml2**(m-1))))

    ! Complex coordinates
    xc = cmplx(x,0.0) + term * nx
    yc = cmplx(y, 0.0) + term * ny
  end if

end subroutine lcpml

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!Function that computes PHI, DPHI, JACOB, INVJACOB at every element in the gaussian integration points
subroutine shape_gauss(xcoord_e,ycoord_e,phi,dphi,jacob,invjacob,detjacob,Ngauss,nodpel,ndim)


implicit none

!Input variables
integer :: Ngauss,nodpel,ndim
complex*16 :: xcoord_e(nodpel),ycoord_e(nodpel)

!Output variables
double precision, intent(out) :: phi(Ngauss,nodpel), dphi(ndim,nodpel)
complex*16, intent(out) :: jacob(ndim,ndim),invjacob(ndim,ndim)
complex*16, intent(out) :: detjacob



!Local variables
integer :: kgauss
double precision :: ksi, eta
double precision, allocatable :: gauss_pt_ksi(:), gauss_pt_eta(:), gauss_wt(:)


!Allocate local variables
allocate(gauss_pt_ksi(Ngauss),gauss_pt_eta(Ngauss), gauss_wt(Ngauss))

gauss_pt_ksi = (/0.5,0.0,0.5/)
gauss_pt_eta = (/0.0,0.5,0.5/)
gauss_wt = 1.0/6.0

do kgauss=1, Ngauss
ksi = gauss_pt_ksi(kgauss)
eta = gauss_pt_eta(kgauss)

phi(kgauss,1) = 1-ksi-eta
phi(kgauss,2) = ksi
phi(kgauss,3) = eta

end do

dphi(1,1) = -1.0
dphi(2,1) = -1.0
           
dphi(1,2) = 1.0
dphi(2,2) = 0.0
           
dphi(1,3) = 0.0
dphi(2,3) = 1.0

jacob(1,1) = xcoord_e(2)-xcoord_e(1)
jacob(1,2) = ycoord_e(2)-ycoord_e(1)
jacob(2,1) = xcoord_e(3)-xcoord_e(1)
jacob(2,2) = ycoord_e(3)-ycoord_e(1)

detjacob = jacob(1,1)*jacob(2,2)-jacob(1,2)*jacob(2,1)

invjacob(1,1) = ycoord_e(3)-ycoord_e(1)
invjacob(1,2) = -(ycoord_e(2)-ycoord_e(1))
invjacob(2,1) = -(xcoord_e(3)-xcoord_e(1))
invjacob(2,2) = xcoord_e(2)-xcoord_e(1)

invjacob = invjacob/detjacob


    end subroutine shape_gauss
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !Function that computes PHI and DPHI in the absolute coordinates
    
subroutine element_matrix(PHI,DPHI,INVJACOB,DETJACOB,pxxe,pyye,pxye,pyxe,qe,AE,Ngauss,nodpel,ndim)
    

!Input variables
double precision :: PHI(Ngauss,nodpel), DPHI(ndim,nodpel)
complex*16 :: INVJACOB(ndim,ndim)
complex*16 :: pxxe,pyye,pxye,pyxe,qe
complex*16 :: DETJACOB

!Output varaibles
complex*16 :: AE(nodpel,nodpel)

!Local variables
complex*16, allocatable :: AE1(:,:), AE2(:,:)
integer :: i,j,k
double precision, allocatable :: gauss_wt(:)
complex*16 :: gauss_sum

allocate(AE1(nodpel,nodpel), AE2(nodpel,nodpel),gauss_wt(Ngauss))

gauss_wt = 1.0/6.0
AE1 = cmplx(0.0,0.0)
AE2 = cmplx(0.0,0.0)


do i=1,nodpel
    do j=1,nodpel
        gauss_sum = cmplx(0.0,0.0)

        AE1(i,j) = (((INVJACOB(1,1) * dphi(1,i) + INVJACOB(1,2) * dphi(2,i))*pxxe + &
                     (INVJACOB(2,1) * dphi(1,i) + INVJACOB(2,2) * dphi(2,i))*pxye) * &
                     (INVJACOB(1,1) * dphi(1,j) + INVJACOB(1,2) * dphi(2,j)) + &
                    ((INVJACOB(1,1) * dphi(1,i) + INVJACOB(1,2) * dphi(2,i))*pyxe + &
                     (INVJACOB(2,1) * dphi(1,i) + INVJACOB(2,2) * dphi(2,i))*pyye * &
                     (INVJACOB(2,1) * dphi(1,j) + INVJACOB(2,2) * dphi(2,j))))*DETJACOB/2
            
!        AE1(i,j)=(pxe*(INVJACOB(1,1)*dphi(1,i)+INVJACOB(1,2)*dphi(2,i)) * (INVJACOB(1,1)*dphi(1,j)+INVJACOB(1,2)*dphi(2,j)) + &
!            pye*(INVJACOB(2,1)*dphi(1,i)+INVJACOB(2,2)*dphi(2,i)) * (INVJACOB(2,1)*dphi(1,j)+INVJACOB(2,2)*dphi(2,j)))*DETJACOB/2
        
        AE1(j,i) = AE1(i,j)
        
        do k=1,Ngauss
            
        gauss_sum = gauss_sum + gauss_wt(k)*qe*PHI(k,i)*PHI(k,j)*DETJACOB
        
        end do
        
        AE2(i,j) = gauss_sum
        AE2(j,i) = AE2(i,j)
        
    end do
end do

AE = AE1 + AE2
    
end subroutine element_matrix

!subroutine gauss_quadrature(Ngauss,gauss_pt_ksi,gauss_pt_eta, gauss_wt)
!implicit none
!double precision :: gauss_pt_ksi(Ngauss), gauss_pt_eta(Ngauss), gauss_wt(Ngauss)

!Ngauss = 3

!gauss_pt_ksi = (/0.5,0.0,0.5/)
!gauss_pt_eta = (/0.0,0.5,0.5/)
!gauss_wt = 1.0/6.0

!end subroutine gauss_quadrature


