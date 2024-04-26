module def_vectors
    
    integer, allocatable :: conectividad(:,:)
    integer, allocatable :: condiciones_fijas(:,:)
    
    double precision, allocatable :: coorx(:), coory(:)
    complex*16, allocatable :: complex_coorx(:), complex_coory(:)
    
    complex*16, allocatable :: u_scat(:), u_inc(:), &
        indep_vect(:), &
        grad(:,:), gradxel(:,:), &
        AD(:), AN(:)
   
    integer, allocatable :: IA(:), JA(:), ncount(:), icx(:)

    !double precision, allocatable :: valor_condicion_fija(:)
    integer, allocatable :: material(:)
    DOUBLE PRECISION, ALLOCATABLE :: masa(:)
    
    DOUBLE PRECISION, ALLOCATABLE :: indep_local_ant(:), acum_ant(:)

    integer, allocatable :: pml_bin_nodes(:), pml_bout_nodes(:), pml_nodes(:)
    integer, allocatable :: scatb_nodes(:),scatin_nodes(:),scatin_elements(:),scatb_elements(:)
    integer, allocatable :: huygb_nodes(:),huygb_elements(:)

end module def_vectors

    