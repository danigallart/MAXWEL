module def_vectors
    
    integer, allocatable :: conn(:,:)                                                                                       !Connectivity matrix
    
    double precision, allocatable :: coorx(:), coory(:)                                                                     !Array of x and y coordinates for nodes
    double precision, allocatable :: coorx_mid(:), coory_mid(:)
    complex*16, allocatable :: complex_coorx(:), complex_coory(:)                                                           !Complex coordiantes for the PML region
    
    complex*16, allocatable :: u_scat(:), u_inc(:), u_tot(:), &                                                             !Scattered wave and incident wave
                               indep_vect(:), &                                                                             !Independent vector of the linear equation system
                               grad(:,:), gradxel(:,:), &                                                                   !Gradients of solutions
                               AD(:), AN(:)                                                                                 !Diagonal and non-diagonal elements of stiff matrix
    
    complex*16, allocatable :: u_inc_mid(:)
    
    complex*16, allocatable :: plane_field_x(:), plane_field_y(:)

    integer, allocatable :: IA(:), JA(:), ncount(:), icx(:)                                                                 !IA: rows of non-zero elements in a CSR format, JA: columns of non-zero elements in a CSR format

    integer, allocatable :: material(:)
    integer, allocatable :: boundary(:)
    
    complex*16, allocatable :: JACOB(:,:),INVJACOB(:,:)
    double precision, allocatable :: PHI(:,:),DPHI(:,:)
    complex*16, allocatable :: DPHIX(:),DPHIY(:)
    
    integer, allocatable :: ns(:)
    complex*16, allocatable :: local_coords(:,:)
    
    
end module def_vectors

    