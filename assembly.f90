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

Ngauss = 4
n_species = 3


allocate(JACOB(ndim,ndim,Ngauss),INVJACOB(ndim,ndim,Ngauss))
allocate(DETJACOB(Ngauss))
allocate(PHI(nodpel,Ngauss), DPHI(ndim, nodpel, Ngauss))
allocate(DPHIX(nodpel,Ngauss),DPHIY(nodpel,Ngauss))
allocate(omg_plasma_species(n_species),omg_cyclotron_species(n_species), &
         mass_species(n_species), charge_species(n_species))
allocate(density_species(n_species,NE), mag_field(NE))

indep_vect=0.0
AD=0.0
AN=0.0

charge_species = (/ -1., 1., 1./)

mass1 = 0.511    !MeV
mass2 = 1875.613 !MeV
mass3 = 2808.921 !MeV

mass_species = (/mass1, mass2, mass3 /)
mass_species = mass_species*1.7827E-30 !kg


do kk=1,NE
    
    !Extract local coordinates of nodes
    do i=1,nodpel
        ns(i) = conn(kk,i)
        local_coords(1,i) = complex_coorx(ns(i))
        local_coords(2,i) = complex_coory(ns(i))
    enddo
    
    if (material(kk) == 1) then
        
        rel_permitivity_xx = mu_scat_xx
        rel_permitivity_yy = mu_scat_yy
        rel_permitivity_zz = mu_scat_zz
        rel_permitivity_xy = mu_scat_xy
        rel_permitivity_yx = mu_scat_yx
        
        if (plasma == 1) then
            
            rel_permitivity_xx = cmplx(0.0,0.0)
            rel_permitivity_yy = cmplx(0.0,0.0)
            rel_permitivity_zz = cmplx(0.0,0.0)
            rel_permitivity_xy = cmplx(0.0,0.0)
            rel_permitivity_yx = cmplx(0.0,0.0)
            
            call density_calculation(deu_tri_frac,local_coords(1,:),local_coords(2,:),density_species(:,kk),nodpel,n_species)
            call magnetic_field_calculation(local_coords(1,:),local_coords(2,:),mag_field(kk),mag_field0,nodpel)
            
            do j = 1,n_species
                
                omg_plasma_species(j) = ( density_species(j,kk) * (charge_species(j) * e_charge)**2 ) / ( e0 * mass_species(j) )
                
                omg_cyclotron_species(j) = ( mag_field(kk) * charge_species(j) * e_charge ) / mass_species(j)
                
                rel_permitivity_xx = rel_permitivity_xx + ( omg_plasma_species(j)**2 ) / (omg**2 - omg_cyclotron_species(j)**2)
                rel_permitivity_yy = rel_permitivity_yy + ( omg_plasma_species(j)**2 ) / (omg**2 - omg_cyclotron_species(j)**2)
                rel_permitivity_xy = rel_permitivity_xy + ( omg_cyclotron_species(j) * omg_plasma_species(j)**2)/(omg * (omg**2 - omg_cyclotron_species(j)**2))
                rel_permitivity_yx = -rel_permitivity_xy
                rel_permitivity_zz = rel_permitivity_zz + ( omg_plasma_species(j)**2 ) / (omg**2)
                
            enddo

            rel_permitivity_xx = 1 - rel_permitivity_xx
            rel_permitivity_yy = 1 - rel_permitivity_yy
            rel_permitivity_zz = 1 - rel_permitivity_zz
            rel_permitivity_xy = -ij * rel_permitivity_xy
            rel_permitivity_yx = -ij * rel_permitivity_yx
            
            
        else
                                    
            im_rel = -cond/(omg*e0)
        
            rel_permitivity_xx = epsilon_scat_xx + ij * im_rel
            rel_permitivity_yy = epsilon_scat_yy + ij * im_rel
            rel_permitivity_zz = epsilon_scat_zz + ij * im_rel
            rel_permitivity_xy = epsilon_scat_xy + ij * im_rel
            rel_permitivity_yx = epsilon_scat_yx + ij * im_rel
            
        endif
        
              
    else
                
            rel_permitivity_xx = cmplx(1.0,0.0)
            rel_permitivity_yy = cmplx(1.0,0.0)
            rel_permitivity_zz = cmplx(1.0,0.0)
            rel_permitivity_xy = cmplx(0.0,0.0)
            rel_permitivity_yx = cmplx(0.0,0.0)
            
            rel_permeability_xx = cmplx(1.0,0.0)
            rel_permeability_yy = cmplx(1.0,0.0)
            rel_permeability_zz = cmplx(1.0,0.0)
            rel_permeability_xy = cmplx(0.0,0.0)
            rel_permeability_yx = cmplx(0.0,0.0)
            
    endif
    
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
    
    
    call element_matrix(PHI,DPHIX,DPHIY,DETJACOB, &
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
double precision, intent(out) :: phi(nodpel,Ngauss), dphi(ndim,nodpel,Ngauss)
complex*16, intent(out) :: jacob(ndim,ndim,Ngauss),invjacob(ndim,ndim,Ngauss)
complex*16, intent(out) :: detjacob(Ngauss)
complex*16, intent(out) :: dphix(nodpel,Ngauss),dphiy(nodpel,Ngauss)



!Local variables
integer :: kgauss, ii
double precision :: ksi, eta
double precision, allocatable :: gauss_pt_ksi(:), gauss_pt_eta(:), gauss_wt(:)


!Allocate local variables
allocate(gauss_pt_ksi(Ngauss),gauss_pt_eta(Ngauss), gauss_wt(Ngauss))

if (Ngauss == 3) then
    gauss_pt_ksi = (/0.5,0.0,0.5/)
    gauss_pt_eta = (/0.0,0.5,0.5/)
    gauss_wt = 1.0/6.0

else if (Ngauss == 4) then
    gauss_pt_ksi = (/ (1./3.), 0.6, 0.2, 0.2/)
    gauss_pt_eta = (/ (1./3.), 0.2, 0.6, 0.2/)
    gauss_wt = (/ -(27./96.), (25./96.), (25./96.), (25./96.)/)
endif


if (nodpel == 3) then

    do kgauss=1, Ngauss
        ksi = gauss_pt_ksi(kgauss)
        eta = gauss_pt_eta(kgauss)

        phi(1,kgauss) = 1-ksi-eta
        phi(2,kgauss) = ksi
        phi(3,kgauss) = eta

        dphi(1,1,kgauss) = -1.0
        dphi(2,1,kgauss) = -1.0
           
        dphi(1,2,kgauss) = 1.0
        dphi(2,2,kgauss) = 0.0
           
        dphi(1,3,kgauss) = 0.0
        dphi(2,3,kgauss) = 1.0

        jacob(1,1,kgauss) = xcoord_e(2)-xcoord_e(1)
        jacob(1,2,kgauss) = ycoord_e(2)-ycoord_e(1)
        jacob(2,1,kgauss) = xcoord_e(3)-xcoord_e(1)
        jacob(2,2,kgauss) = ycoord_e(3)-ycoord_e(1)

        detjacob(kgauss) = jacob(1,1,kgauss)*jacob(2,2,kgauss)-jacob(1,2,kgauss)*jacob(2,1,kgauss)

        invjacob(1,1,kgauss) = (ycoord_e(3)-ycoord_e(1))/detjacob(kgauss)
        invjacob(1,2,kgauss) = -(ycoord_e(2)-ycoord_e(1))/detjacob(kgauss)
        invjacob(2,1,kgauss) = -(xcoord_e(3)-xcoord_e(1))/detjacob(kgauss)
        invjacob(2,2,kgauss) = (xcoord_e(2)-xcoord_e(1))/detjacob(kgauss)
        
        do ii=1,nodpel
            dphix(ii,kgauss) = invjacob(1,1,kgauss) * dphi(1,ii,kgauss) + invjacob(1,2,kgauss) * dphi(2,ii,kgauss)
            dphiy(ii,kgauss) = invjacob(2,1,kgauss) * dphi(1,ii,kgauss) + invjacob(2,2,kgauss) * dphi(2,ii,kgauss)
        enddo
        
        if (detjacob(kgauss)%re <= 0.0) then
            print*, 'WARNING: |J| <= 0', detjacob(kgauss)
        endif
    
    end do
    
else if (nodpel == 6) then
    
    do kgauss=1, Ngauss
        ksi = gauss_pt_ksi(kgauss)
        eta = gauss_pt_eta(kgauss)

        phi(1,kgauss) = (1-ksi-eta)*(1-2*ksi-2*eta)
        phi(2,kgauss) = ksi*(2*ksi-1)
        phi(3,kgauss) = eta*(2*eta-1)
        phi(4,kgauss) = (1-ksi-eta)*4*ksi
        phi(5,kgauss) = 4*ksi*eta
        phi(6,kgauss) = 4*eta*(1-ksi-eta)
        
        dphi(1,1,kgauss) = -3 + 4*ksi + 4*eta
        dphi(2,1,kgauss) = -3 + 4*ksi + 4*eta
        
        dphi(1,2,kgauss) = 4*ksi - 1
        dphi(2,2,kgauss) = 0.0
                
        dphi(1,3,kgauss) = 0.0
        dphi(2,3,kgauss) = 4*eta - 1
        
        dphi(1,4,kgauss) = 4*(1-2*ksi-eta)
        dphi(2,4,kgauss) = -4*ksi
        
        dphi(1,5,kgauss) = 4*eta
        dphi(2,5,kgauss) = 4*ksi
        
        dphi(1,6,kgauss) = -4*eta
        dphi(2,6,kgauss) = 4*(1-ksi-2*eta)

        jacob(1,1,kgauss) = cmplx(0.0,0.0)      
        jacob(1,2,kgauss) = cmplx(0.0,0.0)      
        jacob(2,1,kgauss) = cmplx(0.0,0.0)      
        jacob(2,2,kgauss) = cmplx(0.0,0.0)
        
        do ii = 1,nodpel
            jacob(1,1,kgauss) = jacob(1,1,kgauss) + xcoord_e(ii)*dphi(1,ii,kgauss)
            jacob(1,2,kgauss) = jacob(1,2,kgauss) + ycoord_e(ii)*dphi(1,ii,kgauss)
            jacob(2,1,kgauss) = jacob(2,1,kgauss) + xcoord_e(ii)*dphi(2,ii,kgauss)
            jacob(2,2,kgauss) = jacob(2,2,kgauss) + ycoord_e(ii)*dphi(2,ii,kgauss)
        enddo
        
        detjacob(kgauss) = jacob(1,1,kgauss)*jacob(2,2,kgauss)-jacob(1,2,kgauss)*jacob(2,1,kgauss)
        
        invjacob(1,1,kgauss) = jacob(2,2,kgauss)/detjacob(kgauss)      
        invjacob(1,2,kgauss) = -jacob(1,2,kgauss)/detjacob(kgauss)      
        invjacob(2,1,kgauss) = -jacob(2,1,kgauss)/detjacob(kgauss)      
        invjacob(2,2,kgauss) = jacob(1,1,kgauss)/detjacob(kgauss)
        
        do ii = 1,nodpel
            dphix(ii,kgauss) = invjacob(1,1,kgauss) * dphi(1,ii,kgauss) + invjacob(1,2,kgauss) * dphi(2,ii,kgauss)
            dphiy(ii,kgauss) = invjacob(2,1,kgauss) * dphi(1,ii,kgauss) + invjacob(2,2,kgauss) * dphi(2,ii,kgauss)
        enddo
        
        if (detjacob(kgauss)%re <= 0.0) then
            print*, 'WARNING: |J| <= 0', detjacob(kgauss)
        endif

    end do

endif

    end subroutine shape_gauss
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !Function that computes PHI and DPHI in the absolute coordinates
    
subroutine element_matrix(PHI,DPHIX,DPHIY,DETJACOB,pxxe,pyye,pxye,pyxe,qe,AE,Ngauss,nodpel,ndim)
    

!Input variables
double precision :: PHI(nodpel,Ngauss)
complex*16 :: DPHIX(nodpel,Ngauss), DPHIY(nodpel,Ngauss)
complex*16 :: pxxe,pyye,pxye,pyxe,qe
complex*16 :: DETJACOB(Ngauss)

!Output varaibles
complex*16 :: AE(nodpel,nodpel)

!Local variables
complex*16, allocatable :: AE1(:,:), AE2(:,:)
integer :: i,j,k
double precision, allocatable :: gauss_wt(:)
complex*16 :: gauss_sum1, gauss_sum2

allocate(AE1(nodpel,nodpel), AE2(nodpel,nodpel),gauss_wt(Ngauss))

if (Ngauss == 3) then
    gauss_wt = 1.0/6.0
else if (Ngauss == 4) then
    gauss_wt = (/ -(27./96.), (25./96.), (25./96.), (25./96.)/)
endif

AE1 = cmplx(0.0,0.0)
AE2 = cmplx(0.0,0.0)


do i=1,nodpel
    do j=1,nodpel
        gauss_sum1 = cmplx(0.0,0.0)
        gauss_sum2 = cmplx(0.0,0.0)
        
        do k=1,Ngauss

            gauss_sum1 = gauss_sum1 + gauss_wt(k)*((dphix(i,k)*pxxe + dphiy(i,k)*pxye) * dphix(j,k) + (dphix(i,k)*pyxe + dphiy(i,k)*pyye) * dphiy(j,k))*DETJACOB(k)
        
            gauss_sum2 = gauss_sum2 + gauss_wt(k)*qe*PHI(i,k)*PHI(j,k)*DETJACOB(k)
        
        end do
        
        AE1(i,j) = gauss_sum1
        AE1(j,i) = AE1(i,j)

        AE2(i,j) = gauss_sum2
        AE2(j,i) = AE2(i,j)
        
    end do
end do

AE = AE1 + AE2
    
    end subroutine element_matrix
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
subroutine density_calculation(fraction,complx_coorx,complx_coory,density,nodpel,n_species)

implicit none

integer :: nodpel,n_species
complex*16 :: complx_coorx(nodpel), complx_coory(nodpel)
double precision :: density(n_species)
double precision :: radius_element, xmid, ymid
double precision :: fraction

xmid = sum(real(complx_coorx))/size(complx_coorx)
ymid = sum(real(complx_coory))/size(complx_coory)

radius_element = sqrt(xmid**2 + ymid**2)
    
density(n_species-2) = (1-0.01*radius_element**2)**1.5
density(n_species-1) = fraction * density(n_species-2)
density(n_species) = (1-fraction) * density(n_species-2)

density = density * 5e5
    
end subroutine density_calculation

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    subroutine magnetic_field_calculation(complx_coorx,complx_coory,magnetic_field,max_mag_field,nodpel)

implicit none

integer :: nodpel
complex*16 :: complx_coorx(nodpel), complx_coory(nodpel)
double precision :: magnetic_field, max_mag_field
double precision :: radius_element, xmid, ymid
double precision :: displacement

xmid = sum(real(complx_coorx))/size(complx_coorx)
ymid = sum(real(complx_coory))/size(complx_coory)
    
magnetic_field = max_mag_field/xmid
    
end subroutine magnetic_field_calculation

