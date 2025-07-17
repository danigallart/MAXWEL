module def_vectors
    
    integer, allocatable :: conn(:,:)                                                                                       !Connectivity matrix
    
    double precision, allocatable :: coorx(:), coory(:)                                                                     !Array of x and y coordinates for nodes
    double precision, allocatable :: coorx_mid(:), coory_mid(:)
    complex*16, allocatable :: complex_coorx(:), complex_coory(:)                                                           !Complex coordiantes for the PML region
    
    complex*16, allocatable :: u_scat(:), u_inc(:), u_tot(:), &                                                             !Scattered wave and incident wave
                               indep_vect1(:),indep_vect2(:),indep_vect(:), &                                               !Independent vector of the linear equation system
                               grad(:,:), gradxel(:,:), &                                                                   !Gradients of solutions
                               AD(:), AN(:)                                                                                 !Diagonal and non-diagonal elements of stiff matrix
    
    complex*16, allocatable :: u_inc_mid(:)                                                                                 !Midpoint of the incident wave, used for postprocessing
    
    complex*16, allocatable :: plane_field_x(:), plane_field_y(:)                                                           !Plane wave field in x and y directions, result of postprocessing

    integer, allocatable :: IA(:), JA(:), ncount(:), icx(:)                                                                 !IA: rows of non-zero elements in a CSR format, JA: columns of non-zero elements in a CSR format

    integer, allocatable :: material(:)                                                                                     !Material of the element, 1: vacuum, 2: plasma, 3: PML, 4: Strap 1, 5: Strap 2
    integer, allocatable :: boundary(:)                                                                                     !Nodal boundary conditions
    integer, allocatable :: element_boundary(:,:), boundary_alya(:,:)                                                       !Boundaries by elements and by nodes, used for Alya mesh
    
    complex*16, allocatable :: JACOB(:,:,:),INVJACOB(:,:,:)                                                                 !Jacobian and inverse Jacobian matrices for the elements, used in shape functions
    double precision, allocatable :: PHI(:,:),DPHI(:,:,:)                                                                   !Shape functions and their derivatives for the elements, used in shape functions
    complex*16, allocatable :: DPHIX(:,:),DPHIY(:,:)                                                                        !Derivatives of shape functions in x and y directions, used in shape functions
    complex*16, allocatable :: DETJACOB(:)                                                                                  !Determinant of the Jacobian matrix for the elements, used in shape functions
    
    complex*16, allocatable :: JACOB_1D(:,:), JACOB_1D1(:,:), JACOB_1D2(:,:)                                                !Jacobian matrices for the 1D elements, used in shape functions
    double precision, allocatable :: PHI_1D(:,:),DPHI_1D(:,:), PHI_1D1(:,:),DPHI_1D1(:,:), PHI_1D2(:,:),DPHI_1D2(:,:)       !Shape functions and their derivatives for the 1D elements, used in shape functions
    
    integer, allocatable :: ns(:), ls(:)
    integer, allocatable :: ls1(:), ls2(:)
    complex*16, allocatable :: local_coords(:,:)
    complex*16, allocatable :: coorx_b(:), coory_b(:)                                                                       !Coordinates of the boundary nodes
    complex*16, allocatable :: coorx_b1(:), coory_b1(:)
    complex*16, allocatable :: coorx_b2(:), coory_b2(:)
    
    double precision, allocatable :: mass_species(:)                                                                        !Mass of the species, used in plasma calculations
    double precision, allocatable :: charge_species(:)                                                                      !Charge of the species, used in plasma calculations
    double precision, allocatable :: norm_mag_flux_nodes(:)                                                                 !Normalised magnetic flux on the nodes, used in plasma calculations
    double precision, allocatable :: norm_mag_flux_elements(:)                                                              !Normalised magnetic flux on the elements, used in plasma calculations
    double precision, allocatable :: mag_field(:)                                                                           !Axial magnetic field on the elements, used in plasma calculations
    double precision, allocatable :: density_species(:,:)                                                                   !Density of the species on the elements, used in plasma calculations
        
    
end module def_vectors

    