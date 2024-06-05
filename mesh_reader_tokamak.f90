subroutine mesh_reader_tokamak
    use def_io
    use def_variables
    use def_vectors
    
    implicit none
    
    integer :: istat, ii, jj
        
    double precision :: r
    
    double precision :: min_tol, max_tol
    
    integer :: out_bound1, out_bound2, out_elem
    integer :: ini_1, fini_1, ini_2, fini_2, ini_3, fini_3
    integer :: nod1, nod2, node
    integer :: i, j, n, m
    logical, allocatable :: pml_flag(:)
    
    character*(120) :: text_line
    character*(20) :: label
    
    nodpel = 3 !Nodpel could be read from files (element-wise) but it's not necessary since all elements have the same number of nodes
    
    !Extracting parameters (NE,NP,ndim) from mesh file
    
    do while(text_line /= 'END_DIMENSIONS') 
        read(mesh_param_unit,'(A120)') text_line
        text_line = trim(adjustl(text_line))
        if (index(text_line, 'NODAL_POINTS=') == 1) then
            read(text_line, '(a20,i10)') label,NP
        elseif (index(text_line, 'ELEMENTS=') == 1) then
            read(text_line, '(a20,i10)') label,NE
        elseif(index(text_line, 'SPACE_DIMENSIONS=') == 1) then
            read(text_line, '(a20,i10)') label,ndim
        endif
    enddo
    
    text_line = ''
    
    !Allocate vectors once the parameters are stored
    
    allocate(coorx(NP), coory(NP))
    allocate(conn(NE,nodpel))
    allocate(coorx_mid(NE),coory_mid(NE))
    allocate(material(NE), boundary(NP))
    allocate(pml_flag(NP))
    
    boundary = 0
    pml_flag = .FALSE.

    !Skip irrelevant data
    
    do while(text_line /= 'END_NODES_PER_ELEMENT')
        read(mesh_geo_unit,*) text_line
    enddo
    
    read(mesh_geo_unit, '(a120)') text_line
    
    !Store connectivity matrix
    
    if ( trim(adjustl(text_line)) == 'ELEMENTS') then
        read(mesh_geo_unit,'(A120)') text_line
        do while(text_line /= 'END_ELEMENTS')
            read(text_line,*) ii, conn(ii,nodpel-2), conn(ii,nodpel-1), conn(ii,nodpel)
            read(mesh_geo_unit,'(A120)') text_line
            text_line = trim(adjustl(text_line))
        enddo
    endif
    
    read(mesh_geo_unit, '(a120)') text_line
    
    !Store coordinates array
    
    if ( trim(adjustl(text_line)) == 'COORDINATES') then
    read(mesh_geo_unit,'(A120)') text_line
        do while(text_line /= 'END_COORDINATES')
            read(text_line,*) ii, coorx(ii), coory(ii)
            read(mesh_geo_unit,'(A120)') text_line
            text_line = trim(adjustl(text_line))
        enddo
    endif

    read(mesh_geo_unit, '(a120)') text_line
    
    ! Store material array
    
    read(mesh_element_unit,'(A120)') text_line
    
    if ( trim(adjustl(text_line)) == 'ELEMENTS') then
        read(mesh_element_unit,'(A120)') text_line
        do while(text_line /= 'END_ELEMENTS')
            read(text_line,*) ii, material(ii)
            read(mesh_element_unit,'(A120)') text_line
            text_line = trim(adjustl(text_line))
        enddo
    endif
    
    ! Store boundary array 
    
    read(mesh_boundary_unit,'(A120)') text_line
    if ( trim(adjustl(text_line)) == 'ON_NODES') then
    read(mesh_boundary_unit,'(A120)') text_line
    do while(text_line /= 'END_ON_NODES')
        read(text_line,*) ii, boundary(ii)
        read(mesh_boundary_unit,'(A120)') text_line
        text_line = trim(adjustl(text_line))
    enddo
    endif

    
    do ii=1,NE
        coorx_mid(ii) = sum(coorx(conn(ii,:)))/size(coorx(conn(ii,:)))
        coory_mid(ii) = sum(coory(conn(ii,:)))/size(coory(conn(ii,:)))
        if (material(ii) == 3) then
            do jj = 1, nodpel
                node = conn(ii,jj)
                pml_flag(node) = (boundary(node) /= 2) .and. (boundary(node) /=3)
            enddo
        endif
    enddo
    
    n_scatb = count(boundary == 1, dim=1)
    n_pml_bin = count(boundary == 2, dim=1)
    n_pml_bout = count(boundary == 3, dim=1)

allocate(complex_coorx(NP), complex_coory(NP))

complex_coorx = cmplx(0.0,0.0)
complex_coory = cmplx(0.0,0.0)

call lcpml_tokamak(coorx, coory, k0, boundary, pml_flag, n_pml_bin, n_pml_bout, NP, complex_coorx, complex_coory)


    end subroutine mesh_reader_tokamak
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! LCPML function for Fortran 90
subroutine lcpml_tokamak(x, y, k, boundary_array, flag_array, pml_bin_dim, pml_bout_dim, n_nod, xc, yc)
  use def_io
  use def_variables
  use def_vectors

  implicit none

  ! Input arguments
  integer :: pml_bin_dim, pml_bout_dim, n_nod
  real(kind=8) :: k
  real(kind=8), dimension(n_nod) :: x, y
  integer, dimension(n_nod) :: boundary_array
  logical, dimension(n_nod) :: flag_array
  ! Output arguments
  complex*16, dimension(n_nod) :: xc, yc
  !double precision, intent(out) :: xc_r, yc_r,xc_im, yc_im

  ! LC-PML parameters
  !real(kind=8) :: alpha
  real(kind=8),dimension(pml_bin_dim) :: dpml1
  integer :: m, ind, ii
  complex*16 :: alphajk, term
  real(kind=8) :: x0, y0
  double precision :: vpx(pml_bout_dim),vpy(pml_bout_dim)
  double precision :: npx(pml_bout_dim),npy(pml_bout_dim)
  double precision :: lp(pml_bout_dim)
  double precision :: vx, vy, l, nx, ny, x1, y1,dpml2
  double precision :: ksi
  double precision, dimension(pml_bin_dim) :: pmlbin_x, pmlbin_y
  double precision, dimension(pml_bout_dim) :: pmlbout_x, pmlbout_y
  integer :: count1, count2

count1 = 0
count2 = 0
  
do ii=1,n_nod
    if (boundary_array(ii) == 2) then
        count1 = count1 + 1
        pmlbin_x(count1)=coorx(ii)
        pmlbin_y(count1)=coory(ii)
    else if (boundary_array(ii) == 3) then
        count2 = count2 + 1
        pmlbout_x(count2)=coorx(ii)
        pmlbout_x(count2)=coory(ii)
    endif
enddo
  
do ii=1,n_nod
    if ((boundary_array(ii) == 3) .or. flag_array(ii)) then
        !Set LC-PML parameters
        !alpha = 7.0 * k
        alphajk = cmplx(0.0,-7.0)
        !alphajk = alpha/cmplx(0.0,k)
        m = 3  ! PML decay rate (integer 2 or 3)

        ! Find the point on the inner PML boundary nearest to point P
        dpml1 = sqrt((pmlbin_x - x(ii))**2 + (pmlbin_y - y(ii))**2)
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

        vx = x(ii) - x0
        vy = y(ii) - y0
        l = sqrt(vx**2 + vy**2)
        nx = vx / l  ! Unit vector from r0 to r (x-comp)
        ny = vy / l  ! Unit vector from r0 to r (y-comp)

        if (l < 1.0e-8) then
            xc(ii) = cmplx(x(ii),0.0)
            yc(ii) = cmplx(y(ii),0.0)
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
            xc(ii) = cmplx(x(ii),0.0) + term * nx
            yc(ii) = cmplx(y(ii), 0.0) + term * ny
        endif
        
    else
        xc(ii) = cmplx(x(ii),0.0)
        yc(ii) = cmplx(y(ii),0.0)
    end if
enddo


end subroutine lcpml_tokamak