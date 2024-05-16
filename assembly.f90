SUBROUTINE assembly
use def_io
use def_variables
use def_vectors

implicit none

integer :: kk, ii, i, j, jj, KEJE, IAUX

complex*16 :: determinant

complex*16, allocatable :: AE(:,:)

allocate(local_coords(ndim,nodpel))
allocate(AE(nodpel,nodpel))

Ngauss = 3


allocate(JACOB(ndim,ndim),INVJACOB(ndim,ndim))
allocate(PHI(Ngauss,nodpel), DPHI(ndim, nodpel))
allocate(DPHIX(nodpel),DPHIY(nodpel))

indep_vect=0.0
AD=0.0
AN=0.0


do kk=1,NE
    
    rel_permeability_xx = cmplx(1.0,0.0)
    rel_permeability_xy = cmplx(0.0,0.0)
    rel_permeability_yx = cmplx(0.0,0.0)
    rel_permeability_yy = cmplx(1.0,0.0)
    rel_permeability_zz = cmplx(1.0,0.0)
    
    
    if (material(kk) == 1) then
        
        if (plasma == 1) then
            
            cond = 0.0
            
            im_rel = -cond/(omg*e0)
        
            rel_permitivity_xx = cmplx(9.0,0.0)
            rel_permitivity_yy = cmplx(4.0,0.0)
            rel_permitivity_zz = cmplx(2.0,0.0)
            rel_permitivity_xy = cmplx(0.0,0.0)
            rel_permitivity_yx = cmplx(0.0,0.0)
            
            
        else
            
            cond = 0.0
            
            im_rel = -cond/(omg*e0)
        
            rel_permitivity_xx = cmplx(9.0,im_rel)
            rel_permitivity_yy = cmplx(4.0,im_rel)
            rel_permitivity_zz = cmplx(2.0,im_rel)
            rel_permitivity_xy = cmplx(0.0,im_rel)
            rel_permitivity_yx = cmplx(0.0,im_rel)
            
        endif
        
              
    else
                
            rel_permitivity_xx = cmplx(1.0,0.0)
            rel_permitivity_yy = cmplx(1.0,0.0)
            rel_permitivity_zz = cmplx(1.0,0.0)
            rel_permitivity_xy = cmplx(0.0,0.0)
            rel_permitivity_yx = cmplx(0.0,0.0)
            
    endif
    
    
    do i=1,nodpel
        ns(i) = conn(kk,i)
        local_coords(1,i) = complex_coorx(ns(i))
        local_coords(2,i) = complex_coory(ns(i))
    enddo
    
    call shape_gauss(local_coords(1,:),local_coords(2,:),PHI,DPHI,JACOB,INVJACOB,DETJACOB,DPHIX,DPHIY,Ngauss,nodpel,ndim)
    
    if (pol == 'TM') then
        determinant = rel_permeability_xx * rel_permeability_yy - rel_permeability_xy * rel_permeability_yx
        pxxe = rel_permeability_xx/determinant
        pyye = rel_permeability_yy/determinant
        pyxe = rel_permeability_xy/determinant
        pxye = rel_permeability_yx/determinant
        qe = -rel_permitivity_zz*k0**2
    else if (pol == 'TE') then
        determinant = rel_permitivity_xx * rel_permitivity_yy - rel_permitivity_xy * rel_permitivity_yx
        pxxe = rel_permitivity_xx/determinant
        pyye = rel_permitivity_yy/determinant
        pyxe = rel_permitivity_xy/determinant
        pxye = rel_permitivity_yx/determinant
        qe = -rel_permeability_zz*k0**2
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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!Function that computes PHI, DPHI, JACOB, INVJACOB at every element in the gaussian integration points
subroutine shape_gauss(xcoord_e,ycoord_e,phi,dphi,jacob,invjacob,detjacob,dphix,dphiy,Ngauss,nodpel,ndim)


implicit none

!Input variables
integer :: Ngauss,nodpel,ndim
complex*16 :: xcoord_e(nodpel),ycoord_e(nodpel)

!Output variables
double precision, intent(out) :: phi(Ngauss,nodpel), dphi(ndim,nodpel)
complex*16, intent(out) :: jacob(ndim,ndim),invjacob(ndim,ndim)
complex*16, intent(out) :: detjacob
complex*16, intent(out) :: dphix(nodpel),dphiy(nodpel)



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

dphix(1) = INVJACOB(1,1) * dphi(1,1) + INVJACOB(1,2) * dphi(2,1)
dphix(2) = INVJACOB(1,1) * dphi(1,2) + INVJACOB(1,2) * dphi(2,2)
dphix(3) = INVJACOB(1,1) * dphi(1,3) + INVJACOB(1,2) * dphi(2,3)

dphiy(1) = INVJACOB(2,1) * dphi(1,1) + INVJACOB(2,2) * dphi(2,1)
dphiy(2) = INVJACOB(2,1) * dphi(1,2) + INVJACOB(2,2) * dphi(2,2)
dphiy(3) = INVJACOB(2,1) * dphi(1,3) + INVJACOB(2,2) * dphi(2,3)


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

