SUBROUTINE assembly
use def_io
use def_variables
use def_vectors

implicit none

integer :: kk, ii, i, j, jj, KEJE, IAUX, counter, index1, index2
complex*16 :: determinant

complex*16, allocatable :: AE(:,:)
complex*16, allocatable :: BE(:)
complex*16, allocatable :: integ_line(:),integ_line1(:),integ_line2(:),integ_surf(:)
complex*16 :: AB(3,3)
complex*16 :: BB(3)
complex*16 :: alpha,beta,gamma_1_jin,gamma_2_jin

logical :: mask1(NP),mask2(NP)

complex*16, allocatable :: Au_inc(:)

double precision :: distance,radius

allocate(u_inc(NP),Au_inc(NP))

allocate(local_coords(ndim,nodpel))
allocate(AE(nodpel,nodpel))
allocate(BE(nodpel),integ_line(nodpel),integ_line1(nodpel),integ_line2(nodpel),integ_surf(nodpel))
allocate(ls(2),ls1(2),ls2(2))

n_species = 3
dummy_current = cmplx(0.0,0.0)

allocate(JACOB(ndim,ndim,Ngauss),INVJACOB(ndim,ndim,Ngauss))
allocate(DETJACOB(Ngauss))
allocate(PHI(nodpel,Ngauss), DPHI(ndim, nodpel, Ngauss))
allocate(DPHIX(nodpel,Ngauss),DPHIY(nodpel,Ngauss))
allocate(mass_species(n_species), charge_species(n_species))
allocate(density_species(n_species,NE), mag_field(NE))
allocate (JACOB_1D(ndim,Ngauss),PHI_1D(2,Ngauss),DPHI_1D(2,Ngauss))
allocate (JACOB_1D1(ndim,Ngauss),PHI_1D1(2,Ngauss),DPHI_1D1(2,Ngauss))
allocate (JACOB_1D2(ndim,Ngauss),PHI_1D2(2,Ngauss),DPHI_1D2(2,Ngauss))
allocate (coorx_b(2),coory_b(2))
allocate (coorx_b1(2),coory_b1(2),coorx_b2(2),coory_b2(2))



indep_vect1=cmplx(0.0,0.0)
AD=cmplx(0.0,0.0)
AN=cmplx(0.0,0.0)

mass1 = 0.511    !MeV
mass2 = 1875.613 !MeV
mass3 = 2808.921 !MeV
mass4 = 2808.392 !MeV

if (n_species == 2) then
    charge_species = (/ -1., 1./)
    mass_species = (/mass1, mass2 /)
else if (n_species == 3) then
    charge_species = (/ -1., 1., 1./)
    mass_species = (/mass1, mass2, mass3 /)
else if (n_species == 4) then
    charge_species = (/ -1., 1., 1., 2./)
    mass_species = (/mass1, mass2, mass3 , mass4/)
end if
    
mass_species = mass_species*1.78266E-30 !kg

do kk=1,NE
    
    !Extract local coordinates of nodes
    do i=1,nodpel
        ns(i) = conn(kk,i)
        local_coords(1,i) = complex_coorx(ns(i))
        local_coords(2,i) = complex_coory(ns(i))
    enddo
    
    if (material(kk) == 1) then
        
        rel_permeability_xx = mu_scat_xx
        rel_permeability_yy = mu_scat_yy
        rel_permeability_zz = mu_scat_zz
        rel_permeability_xy = mu_scat_xy
        rel_permeability_yx = mu_scat_yx
        
        if (plasma == 1) then
            
            rel_permitivity_xx = cmplx(0.0,0.0)
            rel_permitivity_yy = cmplx(0.0,0.0)
            rel_permitivity_zz = cmplx(0.0,0.0)
            rel_permitivity_xy = cmplx(0.0,0.0)
            rel_permitivity_yx = cmplx(0.0,0.0)
            
            call density_calculation(density_e_0,deuterium_frac,tritium_frac,helium_3_frac,local_coords(1,:),local_coords(2,:),norm_mag_flux_elements(kk),ka,aa,density_species(:,kk),nodpel,n_species,density_flag)
            call magnetic_field_calculation(local_coords(1,:),local_coords(2,:),norm_mag_flux_elements(kk),major_radius,mag_field_0,mag_field(kk),nodpel,magnetic_flag,elem_shape)
            
            do j = 1,n_species
                
                plasma_freq = sqrt(( density_species(j,kk) * (charge_species(j) * e_charge)**2 ) / ( e0 * mass_species(j) ))
                
                cyclo_freq = ( mag_field(kk) * charge_species(j) * e_charge ) / mass_species(j)
                
                rel_permitivity_xx = rel_permitivity_xx + ( plasma_freq**2 ) / (omg**2 - cyclo_freq**2)
                rel_permitivity_yy = rel_permitivity_yy + ( plasma_freq**2 ) / (omg**2 - cyclo_freq**2)
                rel_permitivity_xy = rel_permitivity_xy + ( cyclo_freq * plasma_freq**2)/(omg * (omg**2 - cyclo_freq**2))
                rel_permitivity_yx = -rel_permitivity_xy
                rel_permitivity_zz = rel_permitivity_zz + ( plasma_freq**2 ) / (omg**2)
                
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
    
    
    call element_matrix(PHI,DPHIX,DPHIY,DETJACOB,pxxe,pyye,pxye,pyxe,qe,AE,Ngauss,nodpel)
!    do i=1,nodpel
!        do j=1,nodpel
!            if (AE(i,j)%im /= 0.0) then
!                print*, kk,material(kk),i,j,AE(i,j)
!            endif
!        end do
!    end do
    
    BE = cmplx(0.0,0.0)
    integ_line = cmplx(0.0,0.0)
    integ_surf = cmplx(0.0,0.0)
    integ_line1 = cmplx(0.0,0.0)
    integ_line2 = cmplx(0.0,0.0)
    mask1 = .FALSE.
    mask2= .FALSE.
    if (antenna_source == 'Y') then
        if (material(kk) == 4) then
            if (pol=='TE') then
                num_bnode_e = count(boundary(ns(:)) == 4)
                if ((num_bnode_e == 2) .and. (nodpel == 3)) then
                    counter = 0
                    do i=1,nodpel
                        if (boundary(ns(i)) == 4) then
                            counter = counter + 1
                            ls(counter) = ns(i)
                            coorx_b(counter) = local_coords(1,i)
                            coory_b(counter) = local_coords(2,i)
                        endif
                    enddo
                call line_integ(coorx_b,coory_b,PHI_1D,pyye,pxxe,dummy_current,current_density1,JACOB_1D,integ_line,Ngauss,nodpel,nodpedge,num_bnode_e,ns,ls,ndim)
                else if ((num_bnode_e == 3) .and. (nodpel == 3)) then
                            ls1 = ns(1:2)
                            ls2 = ns(2:3)
                            
                            coorx_b1 = local_coords(1,1:2)
                            coory_b1 = local_coords(2,2:3)
                            
                            coorx_b2= local_coords(1,1:2)
                            coory_b2 = local_coords(2,2:3)
                            
                    call line_integ(coorx_b1,coory_b1,PHI_1D1,pyye,pxxe,dummy_current,current_density1,JACOB_1D1,integ_line1,Ngauss,nodpel,nodpedge,num_bnode_e,ns,ls1,ndim)
                    call line_integ(coorx_b2,coory_b2,PHI_1D2,pyye,pxxe,dummy_current,current_density1,JACOB_1D2,integ_line2,Ngauss,nodpel,nodpedge,num_bnode_e,ns,ls2,ndim)
                    integ_line = integ_line1 + integ_line2
                    
                else if ((num_bnode_e == 3) .and. (nodpel == 6)) then
                    mask1(ns(1:3)) = .TRUE.
                    mask2(ns(4:6)) = .TRUE.
                    ls1(1)=findloc(boundary,5,1,MASK=mask1)
                    ls1(2)=findloc(boundary,5,1,MASK=mask2)
                    ls2(1)=ls1(2)
                    ls2(2)=findloc(boundary,5,1,MASK=mask1,BACK=.TRUE.)
                    
                    coorx_b1 = local_coords(1,1:2)
                    coory_b1 = local_coords(2,2:3)
                            
                    coorx_b2= local_coords(1,1:2)
                    coory_b2 = local_coords(2,2:3)
                            
                    call line_integ(coorx_b1,coory_b1,PHI_1D1,pyye,pxxe,dummy_current,current_density2,JACOB_1D1,integ_line1,Ngauss,nodpel,nodpedge,num_bnode_e,ns,ls1,ndim)
                    call line_integ(coorx_b2,coory_b2,PHI_1D2,pyye,pxxe,dummy_current,current_density2,JACOB_1D2,integ_line2,Ngauss,nodpel,nodpedge,num_bnode_e,ns,ls2,ndim)
                    integ_line = integ_line1 + integ_line2
                else if ((num_bnode_e == 5) .and. (nodpel == 6)) then

                endif
                call surf_integ(local_coords(1,:),local_coords(2,:),PHI,DPHIX,DPHIY,DETJACOB,pyye,pxxe,dummy_current,current_density1,integ_surf,Ngauss,nodpel)
            else if (pol=='TM') then
             call element_indepvec(PHI,DETJACOB,-ij*k0*nu0*current_density1,BE,Ngauss,nodpel)
            endif
    else if (material(kk) == 5) then
            if (pol=='TE') then
                num_bnode_e = count(boundary(ns(:)) == 5)
                if ((num_bnode_e == 2) .and. (nodpel == 3)) then
                    counter = 0
                    do i=1,nodpel
                        if (boundary(ns(i)) == 5) then
                            counter = counter + 1
                            ls(counter) = ns(i)
                            coorx_b(counter) = local_coords(1,i)
                            coory_b(counter) = local_coords(2,i)
                        endif
                    enddo
                call line_integ(coorx_b,coory_b,PHI_1D,pyye,pxxe,dummy_current,current_density2,JACOB_1D,integ_line,Ngauss,nodpel,nodpedge,num_bnode_e,ns,ls,ndim)
                else if ((num_bnode_e == 3) .and. (nodpel == 3)) then
                            ls1 = ns(1:2)
                            ls2 = ns(2:3)
                            
                            coorx_b1 = local_coords(1,1:2)
                            coory_b1 = local_coords(2,2:3)
                            
                            coorx_b2= local_coords(1,1:2)
                            coory_b2 = local_coords(2,2:3)
                            
                    call line_integ(coorx_b1,coory_b1,PHI_1D1,pyye,pxxe,dummy_current,current_density2,JACOB_1D1,integ_line1,Ngauss,nodpel,nodpedge,num_bnode_e,ns,ls1,ndim)
                    call line_integ(coorx_b2,coory_b2,PHI_1D2,pyye,pxxe,dummy_current,current_density2,JACOB_1D2,integ_line2,Ngauss,nodpel,nodpedge,num_bnode_e,ns,ls2,ndim)
                    integ_line = integ_line1 + integ_line2
                    
                else if ((num_bnode_e == 3) .and. (nodpel == 6)) then
                    mask1(ns(1:3)) = .TRUE.
                    mask2(ns(4:6)) = .TRUE.
                    ls1(1)=findloc(boundary,5,1,MASK=mask1)
                    ls1(2)=findloc(boundary,5,1,MASK=mask2)
                    ls2(1)=ls1(2)
                    ls2(2)=findloc(boundary,5,1,MASK=mask1,BACK=.TRUE.)
                    
                    coorx_b1 = local_coords(1,1:2)
                    coory_b1 = local_coords(2,2:3)
                            
                    coorx_b2= local_coords(1,1:2)
                    coory_b2 = local_coords(2,2:3)
                            
                    call line_integ(coorx_b1,coory_b1,PHI_1D1,pyye,pxxe,dummy_current,current_density2,JACOB_1D1,integ_line1,Ngauss,nodpel,nodpedge,num_bnode_e,ns,ls1,ndim)
                    call line_integ(coorx_b2,coory_b2,PHI_1D2,pyye,pxxe,dummy_current,current_density2,JACOB_1D2,integ_line2,Ngauss,nodpel,nodpedge,num_bnode_e,ns,ls2,ndim)
                    integ_line = integ_line1 + integ_line2
                else if ((num_bnode_e == 5) .and. (nodpel == 6)) then
             
                endif
                call surf_integ(local_coords(1,:),local_coords(2,:),PHI,DPHIX,DPHIY,DETJACOB,pyye,pxxe,dummy_current,current_density2,integ_surf,Ngauss,nodpel)
            else if (pol=='TM') then
             call element_indepvec(PHI,DETJACOB,-ij*k0*nu0*current_density2,BE,Ngauss,nodpel)
            endif
    end if
    endif
    
    
    do i=1, nodpel
        AD(ns(i)) = AD(ns(i)) + AE(i,i)
        indep_vect1(ns(i)) = indep_vect1(ns(i)) + BE(i) + integ_line(i) + integ_surf(i)
        do j=1, nodpel
            do IAUX = 1,ICX(ns(i))
                KEJE = IA(ns(i))+IAUX-1
                if (JA(KEJE) == ns(j)) then
                    AN(KEJE) = AN(KEJE) + AE(i,j)
                endif
            enddo
        enddo
    enddo
    
enddo

do ii = 1,NP
    u_inc(ii) = exp(ij*(k0*(real(complex_coorx(ii))*cos(phii)+real(complex_coory(ii))*sin(phii)))) !No need to use complex coordinates. We just need the real value at each node
enddo
!u_inc%re = cos(k0*(real(complex_coorx)*cos(phii)+real(complex_coory)*sin(phii)))
!u_inc%im = 0.0

indep_vect = cmplx(0.0,0.0)
indep_vect2 = cmplx(0.0,0.0)

if (boundary_type == "ABC") then
    do kk = 1,NE
        index1=0
        index2=0
        inner: do ii = 1,nboun
            if (kk == element_boundary(ii,1)) then
                index1 = findloc(conn(kk,:),element_boundary(ii,2),1)
                index2 = findloc(conn(kk,:),element_boundary(ii,3),1)
                distance = sqrt((coorx(conn(kk,index2))-coorx(conn(kk,index1)))**2+(coory(conn(kk,index2))-coory(conn(kk,index1)))**2)
                exit inner
            endif
        enddo inner
        
        do i=1,nodpel
            ns(i) = conn(kk,i)
        enddo
        
        radius = 2.0
        beta = cmplx(0.0,0.0)!u_inc(conn(kk,index1))
        alpha = cmplx(1.0/(2.0*radius),k0*cos(phii))
        gamma_1_jin = (ij*k0+1./(2.0*radius)-((1./radius)**2)*(1./8.)*(1./(1.0/radius+ij*k0)))
        gamma_2_jin = -0.5/(1.0/radius+ij*k0)
        !gamma_1_jin = ij*k0
        !gamma_2_jin = 0.0!-0.5/(ij*k0)
        AB = cmplx(0.0,0.0)
        BB = cmplx(0.0,0.0)
        !call bc_integ((/coorx(conn(kk,index1)),coorx(conn(kk,index2))/),(/coory(conn(kk,index1)),coory(conn(kk,index2))/),phi_1d,pxxe,pyye,k0,jacob_1d,AB,Ngauss,nodpel,2,2,ndim)
        
        if ((index1 == 1).and.(index2==2)) then
            
            AB(1,1) = AB(1,1) + gamma_1_jin*distance/3.0 - gamma_2_jin/distance
            AB(1,2) = AB(1,2) + gamma_1_jin*distance/6.0 + gamma_2_jin/distance
            AB(2,1) = AB(2,1) + gamma_1_jin*distance/6.0 + gamma_2_jin/distance
            AB(2,2) = AB(2,2) + gamma_1_jin*distance/3.0 - gamma_2_jin/distance
            
            BB(1) = BB(1) + beta*distance/2.0
            BB(2) = BB(2) + beta*distance/2.0
            
        else if ((index1 == 2).and.(index2==3)) then
                        
            AB(2,2) = AB(2,2) + gamma_1_jin*distance/3.0 - gamma_2_jin/distance
            AB(2,3) = AB(2,3) + gamma_1_jin*distance/6.0 + gamma_2_jin/distance
            AB(3,2) = AB(3,2) + gamma_1_jin*distance/6.0 + gamma_2_jin/distance
            AB(3,3) = AB(3,3) + gamma_1_jin*distance/3.0 - gamma_2_jin/distance
            
            BB(2) = BB(2) + beta*distance/2.0
            BB(3) = BB(3) + beta*distance/2.0
            
        else if ((index1 == 3).and.(index2==1)) then
                        
            AB(1,1) = AB(1,1) + gamma_1_jin*distance/3.0 - gamma_2_jin/distance
            AB(1,3) = AB(1,3) + gamma_1_jin*distance/6.0 + gamma_2_jin/distance
            AB(3,1) = AB(3,1) + gamma_1_jin*distance/6.0 + gamma_2_jin/distance
            AB(3,3) = AB(3,3) + gamma_1_jin*distance/3.0 - gamma_2_jin/distance
            
            BB(1) = BB(1) + beta*distance/2.0
            BB(3) = BB(3) + beta*distance/2.0
        
        endif
    
        do i=1, nodpel
        AD(ns(i)) = AD(ns(i)) + AB(i,i)
        indep_vect(ns(i)) = indep_vect(ns(i)) + BB(i)
        do j=1, nodpel
            do IAUX = 1,ICX(ns(i))
                KEJE = IA(ns(i))+IAUX-1
                if (JA(KEJE) == ns(j)) then
                    AN(KEJE) = AN(KEJE) + AB(i,j)
                endif
            enddo
        enddo
    enddo
    enddo
endif


if (plane_wave_source=='Y') then
    if (system_sym=='Y') then
    
    call MATXVECSIM_cplx(NP,IA,JA,AN,AD,u_inc,Au_inc)
  
    elseif (system_sym=='N') then
    
    call MATXVEC_cplex(IA,JA,AN,AD,u_inc,Au_inc,NP,NONULL,0)
  
endif

    do kk=1,NE
        if (material(kk) == 1) then
            do ii = 1, nodpel 
                i = conn(kk,ii)
                indep_vect2(i) = -Au_inc(i)
            enddo
        endif
    enddo
endif

indep_vect = indep_vect1 + indep_vect2

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

if ((Ngauss == 3) .and. (nodpel == 3)) then
    gauss_pt_ksi = (/0.5,0.0,0.5/)
    gauss_pt_eta = (/0.0,0.5,0.5/)
    gauss_wt = 1.0/6.0

else if ((Ngauss == 4) .and. (nodpel == 6)) then
    gauss_pt_ksi = (/ (1./3.), 0.6, 0.2, 0.2/)
    gauss_pt_eta = (/ (1./3.), 0.2, 0.6, 0.2/)
    gauss_wt = (/ -(27./96.), (25./96.), (25./96.), (25./96.)/)

else if ((Ngauss == 4) .and. (nodpel == 4)) then
    gauss_pt_ksi = (/(-1./sqrt(3.)),(1./sqrt(3.)),(-1./sqrt(3.)),(1./sqrt(3.))/)
    gauss_pt_eta = (/(-1./sqrt(3.)),(-1./sqrt(3.)),(1./sqrt(3.)),(1./sqrt(3.))/)
    gauss_wt = (/ 1.0, 1.0, 1.0, 1.0/)
    
else if ((Ngauss == 4) .and. (nodpel == 8)) then
    gauss_pt_ksi = (/(-1./sqrt(3.)),(1./sqrt(3.)),(-1./sqrt(3.)),(1./sqrt(3.))/)
    gauss_pt_eta = (/(-1./sqrt(3.)),(-1./sqrt(3.)),(1./sqrt(3.)),(1./sqrt(3.))/)
    gauss_wt = (/ 1.0, 1.0, 1.0, 1.0/)
    
else if ((Ngauss == 9) .and. (nodpel == 8)) then
    gauss_pt_ksi = (/-sqrt(.6),0.0,sqrt(.6),-sqrt(.6),0.0,sqrt(.6),-sqrt(.6),0.0,sqrt(.6)/)
    gauss_pt_eta = (/-sqrt(.6),-sqrt(.6),-sqrt(.6),0.0,0.0,0.0,sqrt(.6),sqrt(.6),sqrt(.6)/)
    gauss_wt = (/ (25./81.), (40./81.), (25./81.), (40./81.),(64./81.),(40./81.),(25./81.),(40./81.),(25./81.)/)

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
        
        !if (detjacob(kgauss)%re <= 0.0) then
        !    print*, 'WARNING: |J| <= 0', detjacob(kgauss)
        !endif
    
    end do

else if (nodpel == 4) then
    
    do kgauss=1, Ngauss
        ksi = gauss_pt_ksi(kgauss)
        eta = gauss_pt_eta(kgauss)

        phi(1,kgauss) = (1./4.)*(1-ksi)*(1-eta)
        phi(2,kgauss) = (1./4.)*(1+ksi)*(1-eta)
        phi(3,kgauss) = (1./4.)*(1+ksi)*(1+eta)
        phi(4,kgauss) = (1./4.)*(1-ksi)*(1+eta)
        
        dphi(1,1,kgauss) = (1./4.)*(-1+eta)
        dphi(2,1,kgauss) = (1./4.)*(-1+ksi)
        
        dphi(1,2,kgauss) = (1./4.)*(1-eta)
        dphi(2,2,kgauss) = (1./4.)*(-1-ksi)
                
        dphi(1,3,kgauss) = (1./4.)*(1+eta)
        dphi(2,3,kgauss) = (1./4.)*(1+ksi)
        
        dphi(1,4,kgauss) = (1./4.)*(-1-eta)
        dphi(2,4,kgauss) = (1./4.)*(1-ksi)
        

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
        
        !if (detjacob(kgauss)%re <= 0.0) then
        !    print*, 'WARNING: |J| <= 0', detjacob(kgauss)
        !endif

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
        
        !if (detjacob(kgauss)%re <= 0.0) then
        !    print*, 'WARNING: |J| <= 0', detjacob(kgauss)
        !endif

    end do

else if (nodpel == 8) then
    
    do kgauss=1, Ngauss
        ksi = gauss_pt_ksi(kgauss)
        eta = gauss_pt_eta(kgauss)

        phi(1,kgauss) = (1./4.)*(1-ksi)*(1-eta)*(-ksi-eta-1)
        phi(2,kgauss) = (1./4.)*(1+ksi)*(1-eta)*(ksi-eta-1)
        phi(3,kgauss) = (1./4.)*(1+ksi)*(1+eta)*(ksi+eta-1)
        phi(4,kgauss) = (1./4.)*(1-ksi)*(1+eta)*(-ksi+eta-1)
        phi(5,kgauss) = (1./2.)*(1-ksi**2)*(1-eta)
        phi(6,kgauss) = (1./2.)*(1+ksi)*(1-eta**2)
        phi(7,kgauss) = (1./2.)*(1-ksi**2)*(1+eta)
        phi(8,kgauss) = (1./2.)*(1-ksi)*(1-eta**2)
        
        dphi(1,1,kgauss) = (1./4.)*(1-eta)*(2*ksi+eta)
        dphi(2,1,kgauss) = (1./4.)*(1-ksi)*(ksi+2*eta)
        
        dphi(1,2,kgauss) = (1./4.)*(1-eta)*(2*ksi-eta)
        dphi(2,2,kgauss) = (1./4.)*(1+ksi)*(-ksi+2*eta)
                
        dphi(1,3,kgauss) = (1./4.)*(1+eta)*(2*ksi+eta)
        dphi(2,3,kgauss) = (1./4.)*(1+ksi)*(ksi+2*eta)
        
        dphi(1,4,kgauss) = (1./4.)*(1+eta)*(2*ksi-eta)
        dphi(2,4,kgauss) = (1./4.)*(1-ksi)*(-ksi+2*eta)
        
        dphi(1,5,kgauss) = -ksi*(1-eta)
        dphi(2,5,kgauss) = -(1./2.)*(1-ksi**2)
        
        dphi(1,6,kgauss) = (1./2.)*(1-eta**2)
        dphi(2,6,kgauss) = (1+ksi)*(-eta)
        
        dphi(1,7,kgauss) = -ksi*(1+eta)
        dphi(2,7,kgauss) = (1./2.)*(1-ksi**2)
        
        dphi(1,8,kgauss) = -(1./2.)*(1-eta**2)
        dphi(2,8,kgauss) = (1-ksi)*(-eta)

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
        
        !if (detjacob(kgauss)%re <= 0.0) then
        !   print*, 'WARNING: |J| <= 0', detjacob(kgauss)
        !endif

    end do


endif

    end subroutine shape_gauss

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine shape_gauss_1D(xcoord_e,ycoord_e,phi,dphi,jacob_1d,Ngauss,nodpel,nodpedge,num_bnode,ndim)
    
    implicit none

    !Input variables
    integer :: Ngauss,nodpel,ndim, nodpedge,num_bnode
    complex*16 :: xcoord_e(2),ycoord_e(2)

    !Output variables
    double precision, intent(out) :: phi(2,Ngauss), dphi(2,Ngauss)
    complex*16, intent(out) :: jacob_1d(ndim,Ngauss)

    !Local variables
    integer :: kgauss, ii
    double precision :: ksi, gauss_point1, gauss_point2, gauss_weight1, gauss_weight2
    double precision :: gauss_pt_ksi(Ngauss), gauss_wt(Ngauss)


    !Allocate local variables
    !allocate(gauss_pt_ksi(Ngauss), gauss_wt(Ngauss))
    
    
    if (Ngauss == 3) then
        
        gauss_point1 = 0.0
        gauss_point2 = sqrt(0.6)
        gauss_weight1 = 8./9.
        gauss_weight2 = 5./9.
        gauss_pt_ksi = (/gauss_point1, gauss_point2, -gauss_point2/)
        gauss_wt = (/gauss_weight1, gauss_weight2, gauss_weight2/)
        
    else if (Ngauss == 4) then
        
        gauss_point1 = sqrt((3.0/7.0)-(2.0/7.0)*sqrt(6.0/5.0))
        gauss_point2 = sqrt((3.0/7.0)+(2.0/7.0)*sqrt(6.0/5.0))
        gauss_weight1 = (18.+sqrt(30.))/36.
        gauss_weight2 = (18.-sqrt(30.))/36.
        gauss_pt_ksi = (/gauss_point1, -gauss_point1, gauss_point2, -gauss_point2/)
        gauss_wt = (/ gauss_weight1, gauss_weight1, gauss_weight2, gauss_weight2/)
    
    else if (Ngauss == 9) then
        
        gauss_pt_ksi = (/0.0,-0.8360311073266358,0.8360311073266358,-0.9681602395076261,0.9681602395076261,-0.3242534234038089,0.3242534234038089,-0.6133714327005904,0.6133714327005904/)
        gauss_wt = (/ 0.3302393550012598, 0.1806481606948574, 0.1806481606948574, 0.0812743883615744,0.0812743883615744,0.3123470770400029,0.3123470770400029,0.2606106964029354,0.2606106964029354/)

    endif
    
    do kgauss=1, Ngauss
        ksi = gauss_pt_ksi(kgauss)

        phi(1,kgauss) = (1-ksi)*0.5
        phi(2,kgauss) = ksi*0.5

        dphi(1,kgauss) = -0.5
           
        dphi(2,kgauss) = 0.5

        jacob_1d(1,kgauss) = jacob_1d(1,kgauss) + xcoord_e(1)*dphi(1,kgauss) + xcoord_e(2)*dphi(2,kgauss)
        jacob_1d(2,kgauss) = jacob_1d(2,kgauss) + ycoord_e(1)*dphi(1,kgauss) + ycoord_e(2)*dphi(2,kgauss)
    
    end do
    
    
    end subroutine shape_gauss_1D
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    subroutine line_integ(coorx,coory,phi_1d,pxxe,pyye,current_x,current_y,jacob_1d,integ,Ngauss,nodpel,nodpedge,num_bnode,n_global,l_global,ndim)
    
    implicit none
    
    !Input variables
    integer :: Ngauss,nodpel, nodpedge, ndim, num_bnode, n_global(nodpel), l_global(2)
    double precision :: phi_1d(2,Ngauss)
    complex*16 :: pxxe,pyye
    complex*16 :: current_x,current_y
    complex*16 :: jacob_1d(ndim,Ngauss)
    complex*16 :: coorx(2),coory(2)
    double precision :: coorx_mid,coory_mid,radius

    !Output variables
    complex*16, intent(out) :: integ(nodpel)

    !Local variables
    integer :: kgauss, ii, i
    complex*16, allocatable :: integ_x(:),integ_y(:)
    double precision :: dummy_phi(nodpel,Ngauss)
    double precision :: ksi, gauss_point1, gauss_point2, gauss_weight1, gauss_weight2
    double precision :: gauss_pt_ksi(Ngauss), gauss_wt(Ngauss)

    !Allocate local variables
    !allocate(gauss_wt(Ngauss),gauss_pt_ksi(Ngauss),dummy_phi(nodpel,Ngauss))
    allocate(integ_x(nodpel),integ_y(nodpel))
    
    integ = cmplx(0.0,0.0)
    integ_x = cmplx(0.0,0.0)
    integ_y = cmplx(0.0,0.0)
    
    coorx_mid = sum(real(coorx))/size(coorx)
    coory_mid = sum(real(coory))/size(coory)
    radius = sqrt(coorx_mid**2+coory_mid**2)

!    current_x = current_x*coory_mid/radius
!    current_y = -current_y*coorx_mid/radius
    
if (Ngauss == 3) then
        gauss_point1 = 0.0
        gauss_point2 = sqrt(0.6)
        gauss_weight1 = 8./9.
        gauss_weight2 = 5./9.
        gauss_pt_ksi = (/gauss_point1, gauss_point2, -gauss_point2/)
        gauss_wt = (/gauss_weight1, gauss_weight2, gauss_weight2/)

    else if (Ngauss == 4) then
        gauss_point1 = sqrt((3.0/7.0)-(2.0/7.0)*sqrt(6.0/5.0))
        gauss_point2 = sqrt((3.0/7.0)+(2.0/7.0)*sqrt(6.0/5.0))
        gauss_weight1 = (18.+sqrt(30.))/36.
        gauss_weight2 = (18.-sqrt(30.))/36.
        gauss_pt_ksi = (/gauss_point1, -gauss_point1, gauss_point2, -gauss_point2/)
        gauss_wt = (/ gauss_weight1, gauss_weight1, gauss_weight2, gauss_weight2/)
    
    else if (Ngauss == 9) then
        gauss_pt_ksi = (/0.0,-0.8360311073266358,0.8360311073266358,-0.9681602395076261,0.9681602395076261,-0.3242534234038089,0.3242534234038089,-0.6133714327005904,0.6133714327005904/)
        gauss_wt = (/ 0.3302393550012598, 0.1806481606948574, 0.1806481606948574, 0.0812743883615744,0.0812743883615744,0.3123470770400029,0.3123470770400029,0.2606106964029354,0.2606106964029354/)

    endif
    
    call shape_gauss_1D(coorx,coory,gauss_pt_ksi,phi_1d,jacob_1d,Ngauss,nodpel,nodpedge,num_bnode,ndim)
    
    dummy_phi = 0.0
    
    do ii=1,nodpel
        do i = 1,nodpedge
            if (n_global(ii) == l_global(i)) then
                dummy_phi(ii,:) = phi_1d(i,:)
            endif
        enddo
        
        do kgauss=1,Ngauss
            integ_x(ii) = integ_x(ii) + gauss_wt(kgauss)*dummy_phi(ii,kgauss)*pxxe*current_x*jacob_1d(1,kgauss)
            integ_y(ii) = integ_y(ii) + gauss_wt(kgauss)*dummy_phi(ii,kgauss)*pyye*current_y*jacob_1d(2,kgauss)
        enddo
        
    enddo
    
    integ = integ_x + integ_y
    
    !print*, 'Line integral'
    !print*, integ
    
    end subroutine line_integ
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    subroutine bc_integ(coorx,coory,phi_1d,pxxe,pyye,k0,jacob_1d,AB,Ngauss,nodpel,nodpedge,num_bnode,ndim)
    
    implicit none
    
    !Input variables
    integer :: Ngauss,nodpel, nodpedge, ndim, num_bnode
    double precision :: phi_1d(2,Ngauss), k0
    complex*16 :: pxxe,pyye
    complex*16 :: jacob_1d(ndim,Ngauss)
    complex*16 :: coorx(2),coory(2),length,nxe,nye

    !Output variables
    complex*16, intent(out) :: AB(nodpel,nodpel)

    !Local variables
    integer :: kgauss, ii, jj, i
    double precision :: ksi, gauss_point1, gauss_point2, gauss_weight1, gauss_weight2
    double precision :: gauss_pt_ksi(Ngauss), gauss_wt(Ngauss)
    complex*16 :: ij

    !Allocate local variables
    !allocate(gauss_wt(Ngauss),gauss_pt_ksi(Ngauss),dummy_phi(nodpel,Ngauss))
    ij = cmplx(0.0,1.0)
    AB = cmplx(0.0,0.0)
    
if (Ngauss == 3) then
        gauss_point1 = 0.0
        gauss_point2 = sqrt(0.6)
        gauss_weight1 = 8./9.
        gauss_weight2 = 5./9.
        gauss_pt_ksi = (/gauss_point1, gauss_point2, -gauss_point2/)
        gauss_wt = (/gauss_weight1, gauss_weight2, gauss_weight2/)

    else if (Ngauss == 4) then
        gauss_point1 = sqrt((3.0/7.0)-(2.0/7.0)*sqrt(6.0/5.0))
        gauss_point2 = sqrt((3.0/7.0)+(2.0/7.0)*sqrt(6.0/5.0))
        gauss_weight1 = (18.+sqrt(30.))/36.
        gauss_weight2 = (18.-sqrt(30.))/36.
        gauss_pt_ksi = (/gauss_point1, -gauss_point1, gauss_point2, -gauss_point2/)
        gauss_wt = (/ gauss_weight1, gauss_weight1, gauss_weight2, gauss_weight2/)
    
    else if (Ngauss == 9) then
        gauss_pt_ksi = (/0.0,-0.8360311073266358,0.8360311073266358,-0.9681602395076261,0.9681602395076261,-0.3242534234038089,0.3242534234038089,-0.6133714327005904,0.6133714327005904/)
        gauss_wt = (/ 0.3302393550012598, 0.1806481606948574, 0.1806481606948574, 0.0812743883615744,0.0812743883615744,0.3123470770400029,0.3123470770400029,0.2606106964029354,0.2606106964029354/)

    endif
    
    call shape_gauss_1D(coorx,coory,gauss_pt_ksi,phi_1d,jacob_1d,Ngauss,nodpel,nodpedge,num_bnode,ndim)
    
    length = sqrt((coorx(2)-coorx(1))**2+(coory(2)-coory(1))**2)
    
    nxe = (coory(2)-coory(1))/length
    nye = (coorx(2)-coorx(1))/length
    
    
    do ii=1,nodpel
        do jj=1,nodpel
            do kgauss=1,Ngauss
                AB(ii,jj) = AB(ii,jj) + gauss_wt(kgauss)*PHI_1D(ii,kgauss)*PHI_1D(jj,kgauss)*(pxxe*nxe+pyye*nye)*ij*k0*length/2
            enddo
        enddo
    enddo
        
    !print*, 'Line integral'
    !print*, integ
    
    end subroutine bc_integ
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
subroutine surf_integ(coorx,coory,PHI,DPHIX,DPHIY,DETJACOB,pxxe,pyye,current_x,current_y,integ_surface,Ngauss,nodpel)

    !Input variables
    double precision :: PHI(nodpel,Ngauss)
    complex*16 :: DETJACOB(Ngauss)
    complex*16 :: current_x, current_y
    complex*16 :: pxxe, pyye
    complex*16 :: DPHIX(nodpel,Ngauss), DPHIY(nodpel,Ngauss)
    complex*16 :: coorx(nodpel),coory(nodpel)
    double precision :: coorx_mid,coory_mid,radius

    
    !Output varaibles
    complex*16 :: integ_surface(nodpel)
    
    !Local variables
    integer :: i,k
    double precision :: gauss_wt(Ngauss)
    complex*16 :: gauss_sum
    
    coorx_mid = sum(real(coorx))/size(coorx)
    coory_mid = sum(real(coory))/size(coory)
    radius = sqrt(coorx_mid**2+coory_mid**2)

    !current_x = current_x*coory_mid/radius
    !current_y = -current_y*coorx_mid/radius
    
if ((Ngauss == 3) .and. (nodpel == 3)) then
    gauss_wt = 1.0/6.0
else if ((Ngauss == 4) .and. (nodpel == 6)) then
    gauss_wt = (/ -(27./96.), (25./96.), (25./96.), (25./96.)/)
else if ((Ngauss == 4) .and. (nodpel == 4)) then
    gauss_wt = (/ 1.0, 1.0, 1.0, 1.0/)
else if ((Ngauss == 4) .and. (nodpel == 8)) then
    gauss_wt = (/ 1.0, 1.0, 1.0, 1.0/)
else if ((Ngauss == 9) .and. (nodpel == 8)) then
    gauss_wt = (/ (25./81.), (40./81.), (25./81.), (40./81.), (64./81.), (40./81.), (25./81.), (40./81.), (25./81.)/)
endif

integ_surface = cmplx(0.0,0.0)

do i=1,nodpel
    do k=1,Ngauss
        
            integ_surface(i) = integ_surface(i) + gauss_wt(k)*(DPHIY(i,k)*pxxe*current_x-DPHIX(i,k)*pyye*current_y)*DETJACOB(k)
    
        enddo
enddo

    !print*, 'Surface integral'
    !print*, integ_surface

    end subroutine surf_integ
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !Function that computes PHI and DPHI in the absolute coordinates
    
subroutine element_matrix(PHI,DPHIX,DPHIY,DETJACOB,pxxe,pyye,pxye,pyxe,qe,AE,Ngauss,nodpel)
    

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

if ((Ngauss == 3) .and. (nodpel == 3)) then
    gauss_wt = 1.0/6.0
else if ((Ngauss == 4) .and. (nodpel == 6)) then
    gauss_wt = (/ -(27./96.), (25./96.), (25./96.), (25./96.)/)
else if ((Ngauss == 4) .and. (nodpel == 4)) then
    gauss_wt = (/ 1.0, 1.0, 1.0, 1.0/)
else if ((Ngauss == 4) .and. (nodpel == 8)) then
    gauss_wt = (/ 1.0, 1.0, 1.0, 1.0/)
else if ((Ngauss == 9) .and. (nodpel == 8)) then
    gauss_wt = (/ (25./81.), (40./81.), (25./81.), (40./81.), (64./81.), (40./81.), (25./81.), (40./81.), (25./81.)/)
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
        !AE1(j,i) = AE1(i,j)

        AE2(i,j) = gauss_sum2
        !AE2(j,i) = AE2(i,j)
                
    end do
end do

AE = AE1 + AE2
    
    end subroutine element_matrix
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
subroutine element_indepvec(PHI,DETJACOB,source_current,BE,Ngauss,nodpel)

    !Input variables
    double precision :: PHI(nodpel,Ngauss)
    complex*16 :: DETJACOB(Ngauss)
    complex*16 :: source_current
    
    !Output varaibles
    complex*16 :: BE(nodpel)
    
    !Local variables
    integer :: i,k
    double precision, allocatable :: gauss_wt(:)
    complex*16 :: gauss_sum
    
    allocate(gauss_wt(Ngauss))
    
if ((Ngauss == 3) .and. (nodpel == 3)) then
    gauss_wt = 1.0/6.0
else if ((Ngauss == 4) .and. (nodpel == 6)) then
    gauss_wt = (/ -(27./96.), (25./96.), (25./96.), (25./96.)/)
else if ((Ngauss == 4) .and. (nodpel == 4)) then
    gauss_wt = (/ 1.0, 1.0, 1.0, 1.0/)
else if ((Ngauss == 4) .and. (nodpel == 8)) then
    gauss_wt = (/ 1.0, 1.0, 1.0, 1.0/)
else if ((Ngauss == 9) .and. (nodpel == 8)) then
    gauss_wt = (/ (25./81.), (40./81.), (25./81.), (40./81.), (64./81.), (40./81.), (25./81.), (40./81.), (25./81.)/)
endif

BE = cmplx(0.0,0.0)

do i=1,nodpel
    do k=1,Ngauss
        
            BE(i) = BE(i) - gauss_wt(k)*PHI(i,k)*source_current*DETJACOB(k)
    
        enddo
enddo


    end subroutine element_indepvec
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
subroutine density_calculation(density_e_0,deuterium_frac,tritium_frac,helium_3_frac,complx_coorx,complx_coory,magnetic_flux,k,a,density,nodpel,n_species,flag)

implicit none

integer :: nodpel,n_species, flag
complex*16 :: complx_coorx(nodpel), complx_coory(nodpel)
double precision :: density(n_species)
double precision :: radius_element, xmid, ymid
double precision :: deuterium_frac, tritium_frac, helium_3_frac, magnetic_flux, k, a
real(kind=8) :: density_e_0

xmid = sum(real(complx_coorx))/size(complx_coorx)
ymid = sum(real(complx_coory))/size(complx_coory)

radius_element = sqrt(xmid**2 + ymid**2)
    
if (flag == 1) then
    
    density(1) = density_e_0 * (1-0.01*radius_element**2)**1.5
    
else if (flag == 2) then

    density(1) = density_e_0 * (k + (1-k)*(1-magnetic_flux**2)**a)
    
end if

if (n_species == 2) then

    density(2) = density(1)
    
else if (n_species == 3) then
    
    density(2) = deuterium_frac * density(1)
    density(3) = tritium_frac * density(1)
    
else if (n_species == 4) then
    
    density(2) = deuterium_frac * density(1)
    density(3) = tritium_frac * density(1)
    density(4) = helium_3_frac * density(1)
    
end if



end subroutine density_calculation

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
subroutine magnetic_field_calculation(complx_coorx,complx_coory,magnetic_flux,major_radius,mag_field_0,magnetic_field,nodpel,flag,element_shape)

implicit none

integer :: nodpel, flag
complex*16 :: complx_coorx(nodpel), complx_coory(nodpel)
double precision :: magnetic_field, magnetic_flux ,mag_field_0
double precision :: radius_element,major_radius, xmid, ymid
double precision :: area1,area2,area
character*4 :: element_shape

xmid = sum(real(complx_coorx))/size(complx_coorx)
ymid = sum(real(complx_coory))/size(complx_coory)

if (flag == 1) then
    
    magnetic_field = mag_field_0*major_radius/xmid
    !magnetic_field = max_mag_field

else if (flag == 2) then
    if (element_shape == 'tria') then
        area = abs(real(complx_coorx(1))*real(complx_coory(2)-complx_coory(3))+real(complx_coorx(2))*real(complx_coory(3)-complx_coory(1))+real(complx_coorx(3))*real(complx_coory(1)-complx_coory(2)))/2
    else if (element_shape == 'squa') then
        area1 = abs(real(complx_coorx(1))*real(complx_coory(2)-complx_coory(3))+real(complx_coorx(2))*real(complx_coory(3)-complx_coory(1))+real(complx_coorx(3))*real(complx_coory(1)-complx_coory(2)))/2
        area2 = abs(real(complx_coorx(1))*real(complx_coory(3)-complx_coory(4))+real(complx_coorx(3))*real(complx_coory(4)-complx_coory(1))+real(complx_coorx(4))*real(complx_coory(1)-complx_coory(3)))/2
        area=area1+area2
    end if
    
    magnetic_field = magnetic_flux/area
    
end if
    
    end subroutine magnetic_field_calculation
    
