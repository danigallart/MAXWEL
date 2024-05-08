module def_vectors
    
    integer, allocatable :: conn(:,:)                                                                               !Connectivity matrix
    !integer, allocatable :: condiciones_fijas(:,:)
    
    double precision, allocatable :: coorx(:), coory(:)                                                                     !Array of x and y coordinates for nodes
    double precision, allocatable :: coorx_mid(:), coory_mid(:)
    complex*16, allocatable :: complex_coorx(:), complex_coory(:)                                                           !Complex coordiantes for the PML region
    
    complex*16, allocatable :: u_scat(:), u_inc(:), u_tot(:), &                                                             !Scattered wave and incident wave
                               indep_vect(:), &                                                                             !Independent vector of the linear equation system
                               grad(:,:), gradxel(:,:), &                                                                   !Gradients of solutions
                               AD(:), AN(:)                                                                                 !Diagonal and non-diagonal elements of stiff matrix
   
    integer, allocatable :: IA(:), JA(:), ncount(:), icx(:)                                                                 !IA: rows of non-zero elements in a CSR format, JA: columns of non-zero elements in a CSR format

    integer, allocatable :: pml_bin_nodes(:), pml_bout_nodes(:), pml_nodes(:)                                               !Inner boundary, outer boundary and general nodes of PML region
    integer, allocatable :: scatb_nodes(:),scatin_nodes(:),scatin_elements(:),scatb_elements(:)                             !Boundary nodes, inner nodes, boundary elements and inner elements of scaterer element
    integer, allocatable :: huygb_nodes(:),huygb_elements(:)                                                                !Huygens surface, only used to compute RCS and farfield.
    integer, allocatable :: material(:)
    integer, allocatable :: boundary(:)
    
    
end module def_vectors

    