subroutine derivatives
    
    use def_io
    use def_variables
    use def_vectors
    
    implicit none
    
    integer :: ii,kk,i
    complex*16 :: DSOLX, DSOLY
    complex*16 :: inv_tensor_xx,inv_tensor_xy, &
                  inv_tensor_yx,inv_tensor_yy
    complex*16 :: determinant_tensor

    
    double precision :: r
    
    double precision :: min_tol, max_tol
    
    allocate(plane_field_x(NE),plane_field_y(NE))
    allocate(u_inc_mid(NE))
        
    min_tol = 0.95
    max_tol = 1.05
    
    do kk = 1,NE
        
        DSOLX = cmplx(0.0,0.0)
        DSOLY = cmplx(0.0,0.0)
        
        do ii = 1, nodpel
            i = conn(kk,ii)
            
            local_coords(1,ii) = complex_coorx(i)!cmplx(coorx(i),0.0)
            local_coords(2,ii) = complex_coory(i)!cmplx(coory(i),0.0)
            
        enddo
            
        call shape_gauss(local_coords(1,:),local_coords(2,:),PHI,DPHI,JACOB,INVJACOB,DETJACOB,DPHIX,DPHIY,Ngauss,nodpel,ndim)
            
        do ii = 1, nodpel
            
            i = conn(kk,ii)                   
            
            DSOLX = DSOLX + u_tot(i)*DPHIX(ii)
            DSOLY = DSOLY + u_tot(i)*DPHIY(ii)
            
            !DSOLX = DSOLX + u_scat(i)*DPHIX(ii)/(nodpel*nodpel)
            !DSOLY = DSOLY + u_scat(i)*DPHIY(ii)/(nodpel*nodpel)
        enddo
        
            rel_permeability_xx = cmplx(1.0,0.0)
            rel_permeability_yy = cmplx(1.0,0.0)
            rel_permeability_zz = cmplx(1.0,0.0)
            rel_permeability_xy = cmplx(0.0,0.0)
            rel_permeability_yx = cmplx(0.0,0.0)
        
        if (material(kk) == 1) then !Falta Plasma
            
            cond = 0.0
            
            im_rel = -cond/(omg*e0)
        
            rel_permitivity_xx = cmplx(9.0,im_rel)
            rel_permitivity_yy = cmplx(4.0,im_rel)
            rel_permitivity_zz = cmplx(2.0,im_rel)
            rel_permitivity_xy = cmplx(0.0,im_rel)
            rel_permitivity_yx = cmplx(0.0,im_rel)
            
        else
            
            rel_permitivity_xx = cmplx(1.0,0.0)
            rel_permitivity_yy = cmplx(1.0,0.0)
            rel_permitivity_zz = cmplx(1.0,0.0)
            rel_permitivity_xy = cmplx(0.0,0.0)
            rel_permitivity_yx = cmplx(0.0,0.0)  
            
        endif
        
        !u_inc_mid = exp(ij*k0*(coorx_mid*cos(phii)+coory_mid*sin(phii)))
        
        if (pol == 'TM') then
            
            determinant_tensor = rel_permeability_xx * rel_permeability_yy - rel_permeability_xy * rel_permeability_yx
            inv_tensor_xx = rel_permeability_yy/determinant_tensor
            inv_tensor_xy = -rel_permeability_xy/determinant_tensor
            inv_tensor_yx = -rel_permeability_yx/determinant_tensor
            inv_tensor_yy = rel_permeability_xx/determinant_tensor
        
            plane_field_x(kk) = -(inv_tensor_xx*DSOLY - inv_tensor_xy*DSOLX)/(ij*omg*mu0)
            plane_field_y(kk) = -(inv_tensor_yx*DSOLY - inv_tensor_yy*DSOLX)/(ij*omg*mu0)
    
        elseif (pol == 'TE') then
            
            determinant_tensor = rel_permitivity_xx * rel_permitivity_yy - rel_permitivity_xy * rel_permitivity_yx
            inv_tensor_xx = rel_permitivity_yy/determinant_tensor
            inv_tensor_xy = -rel_permitivity_xy/determinant_tensor
            inv_tensor_yx = -rel_permitivity_yx/determinant_tensor
            inv_tensor_yy = rel_permitivity_xx/determinant_tensor
        
            plane_field_x(kk) = (inv_tensor_xx*DSOLY - inv_tensor_xy*DSOLX)/(ij*omg*e0)
            plane_field_y(kk) = (inv_tensor_yx*DSOLY - inv_tensor_yy*DSOLX)/(ij*omg*e0)
    
        endif
        
        r = sqrt(coorx_mid(kk)**2 + coory_mid(kk)**2)
         
        if ((r < rpmlin * max_tol) .and. (r > rpmlin * min_tol)) then
        
            plane_field_x(kk) = cmplx(0.0,0.0)
            plane_field_y(kk) = cmplx(0.0,0.0)
        
        endif
        
        !if (material(kk) == 3) then
        !    
        !    plane_field_x(kk) = cmplx(0.0,0.0)
        !    plane_field_y(kk) = cmplx(0.0,0.0)
        !
        !endif
        
    enddo
       
    
    end subroutine derivatives