module def_vectors
    
    integer, allocatable :: conn(:,:)                                                                                       !Connectivity matrix
    
    double precision, allocatable :: coorx(:), coory(:)                                                                     !Array of x and y coordinates for nodes
    double precision, allocatable :: coorx_mid(:), coory_mid(:)
    complex*16, allocatable :: complex_coorx(:), complex_coory(:)                                                           !Complex coordiantes for the PML region
    
    complex*16, allocatable :: u_scat(:), u_inc(:), u_tot(:), &                                                             !Scattered wave and incident wave
                               indep_vect1(:),indep_vect2(:),indep_vect(:), &                                               !Independent vector of the linear equation system
                               grad(:,:), gradxel(:,:), &                                                                   !Gradients of solutions
                               AD(:), AN(:)                                                                                 !Diagonal and non-diagonal elements of stiff matrix
    
    complex*16, allocatable :: u_inc_mid(:)
    
    complex*16, allocatable :: plane_field_x(:), plane_field_y(:)

    integer, allocatable :: IA(:), JA(:), ncount(:), icx(:)                                                                 !IA: rows of non-zero elements in a CSR format, JA: columns of non-zero elements in a CSR format

    integer, allocatable :: material(:)
    integer, allocatable :: boundary(:)
    integer, allocatable :: element_boundary(:,:), boundary_alya(:,:)
    
    complex*16, allocatable :: JACOB(:,:,:),INVJACOB(:,:,:)
    double precision, allocatable :: PHI(:,:),DPHI(:,:,:)
    complex*16, allocatable :: DPHIX(:,:),DPHIY(:,:)
    complex*16, allocatable :: DETJACOB(:)
    
    complex*16, allocatable :: JACOB_1D(:,:), JACOB_1D1(:,:), JACOB_1D2(:,:)
    double precision, allocatable :: PHI_1D(:,:),DPHI_1D(:,:), PHI_1D1(:,:),DPHI_1D1(:,:), PHI_1D2(:,:),DPHI_1D2(:,:)
    
    integer, allocatable :: ns(:), ls(:)
    integer, allocatable :: ls1(:), ls2(:)
    complex*16, allocatable :: local_coords(:,:)
    complex*16, allocatable :: coorx_b(:), coory_b(:)
    complex*16, allocatable :: coorx_b1(:), coory_b1(:)
    complex*16, allocatable :: coorx_b2(:), coory_b2(:)
    
    double precision, allocatable :: mass_species(:)
    double precision, allocatable :: charge_species(:)
    double precision, allocatable :: norm_mag_flux_nodes(:)
    double precision, allocatable :: norm_mag_flux_elements(:)
    double precision, allocatable :: mag_field(:)
    double precision, allocatable :: density_species(:,:)
    
    double precision, allocatable :: source_coorx(:),source_coory(:)
    integer, allocatable :: source_node(:), source_element(:)
    !complex*16, allocatable :: current_density(:)
    
    
end module def_vectors

    