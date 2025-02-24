subroutine mesh_reader
    use def_io
    use def_variables
    use def_vectors
    
    implicit none
    
    integer :: istat, ii, jj, node
            
    double precision :: min_tol, max_tol
    logical, allocatable :: pml_flag(:)
    double precision, allocatable :: dist(:)
    
    character*(120) :: text_line

    open(unit=mesh_unit,file=mesh_file,status='old')
    
    !Extracting parameters (NE,NP) from mesh file

    do while(text_line /= 'END_DIMENSIONS') 
        read(mesh_unit,'(A120)') text_line
        text_line = trim(adjustl(text_line))
        if (text_line == 'NODAL_POINTS') then
            read(mesh_unit,'(A120)') text_line
            text_line = trim(adjustl(text_line))
            read(text_line, '(i10)') NP
        elseif (text_line == 'NUM_ELEMENTS') then
            read(mesh_unit,'(A120)') text_line
            text_line = trim(adjustl(text_line))
            read(text_line, '(i10)') NE
        endif
    enddo
    
    !Allocate vectors once the parameters are stored
    
    allocate(coorx(NP), coory(NP), dist(NP))
    allocate(conn(NE,nodpel))
    allocate(coorx_mid(NE),coory_mid(NE))
    allocate(material(NE), boundary(NP))
    allocate(pml_flag(NP))
    allocate(complex_coorx(NP), complex_coory(NP))
    
    boundary = 0
    pml_flag = .FALSE.
    
    read(mesh_unit, '(a120)') text_line
    
    !Store connectivity matrix
    if (nodpel == 3) then
        if ( trim(adjustl(text_line)) == 'ELEMENTS') then
            read(mesh_unit,'(A120)') text_line
            do while(text_line /= 'END_ELEMENTS')
                read(text_line,*) ii, conn(ii,nodpel-2), conn(ii,nodpel-1), conn(ii,nodpel)
                read(mesh_unit,'(A120)') text_line
                text_line = trim(adjustl(text_line))
            enddo
        endif

    else if (nodpel == 4) then
        if ( trim(adjustl(text_line)) == 'ELEMENTS') then
            read(mesh_unit,'(A120)') text_line
            do while(text_line /= 'END_ELEMENTS')
                read(text_line,*) ii, conn(ii,nodpel-3), conn(ii,nodpel-2), conn(ii,nodpel-1), conn(ii,nodpel)
                read(mesh_unit,'(A120)') text_line
                text_line = trim(adjustl(text_line))
            enddo
        endif
      
    else if (nodpel == 6) then
        if ( trim(adjustl(text_line)) == 'ELEMENTS') then
            read(mesh_unit,'(A120)') text_line
            do while(text_line /= 'END_ELEMENTS')
                read(text_line,*) ii, conn(ii,nodpel-5), conn(ii,nodpel-4), conn(ii,nodpel-3), conn(ii,nodpel-2), conn(ii,nodpel-1), conn(ii,nodpel)
                read(mesh_unit,'(A120)') text_line
                text_line = trim(adjustl(text_line))
            enddo
        endif

    else if (nodpel == 8) then
        if ( trim(adjustl(text_line)) == 'ELEMENTS') then
            read(mesh_unit,'(A120)') text_line
            do while(text_line /= 'END_ELEMENTS')
                read(text_line,*) ii, conn(ii,nodpel-7), conn(ii,nodpel-6), conn(ii,nodpel-5), conn(ii,nodpel-4), conn(ii,nodpel-3), conn(ii,nodpel-2), conn(ii,nodpel-1), conn(ii,nodpel)
                read(mesh_unit,'(A120)') text_line
                text_line = trim(adjustl(text_line))
            enddo
        endif
    endif
        
    !Store coordinates array
   
    read(mesh_unit, '(a120)') text_line
    
    if ( trim(adjustl(text_line)) == 'COORDINATES') then
    read(mesh_unit,'(A120)') text_line
        do while(text_line /= 'END_COORDINATES')
            read(text_line,*) ii, coorx(ii), coory(ii)
            read(mesh_unit,'(A120)') text_line
            text_line = trim(adjustl(text_line))
        enddo
    endif
    
    ! Store material array
    
    read(mesh_unit,'(A120)') text_line
    
    if ( trim(adjustl(text_line)) == 'MATERIALS') then
        read(mesh_unit,'(A120)') text_line
        do while(text_line /= 'END_MATERIALS')
            read(text_line,*) ii, material(ii)
            read(mesh_unit,'(A120)') text_line
            text_line = trim(adjustl(text_line))
        enddo
    endif
    
    ! Store boundary array
    
    read(mesh_unit,'(A120)') text_line

    if ( trim(adjustl(text_line)) == 'ON_NODES') then
    read(mesh_unit,'(A120)') text_line
    do while(text_line /= 'END_ON_NODES')
        read(text_line,*) ii, boundary(ii)
        read(mesh_unit,'(A120)') text_line
        text_line = trim(adjustl(text_line))
    enddo
    endif
    
    
    close(mesh_unit)
    
    coorx = coorx !+ major_radius
    coory = coory
    coorx_mid = coorx_mid !+ major_radius
    coory_mid = coory_mid
    complex_coorx = cmplx(0.0,0.0)
    complex_coory = cmplx(0.0,0.0)
    
    do ii=1,NE
        coorx_mid(ii) = sum(coorx(conn(ii,:)))/size(coorx(conn(ii,:)))
        coory_mid(ii) = sum(coory(conn(ii,:)))/size(coory(conn(ii,:)))
        if (boundary_type == 'PML') then
            if (material(ii) == 3) then
                do jj = 1, nodpel
                    node = conn(ii,jj)
                    pml_flag(node) = (boundary(node) /= 2) .and. (boundary(node) /=3)
                enddo
            endif
        endif
    enddo
    
    
    if (boundary_type == 'PML') then
        n_pml_bin = count(boundary == 2, dim=1)
        n_pml_bout = count(boundary == 3, dim=1)
        
        call lcpml(coorx, coory, k0, boundary, pml_flag, n_pml_bin, n_pml_bout, NP, complex_coorx, complex_coory)

    else if (boundary_type=='ABC') then
        
        complex_coorx%re = coorx
        complex_coory%re = coory
        
    end if
    
!do jj = 1, source_num
!    dist = sqrt((coorx-source_coorx(jj))**2+(coory-source_coory(jj))**2)
!    source_node(jj) = minloc(dist,1)
!end do
    
    
    end subroutine mesh_reader
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! LCPML function for Fortran 90
subroutine lcpml(x, y, k, boundary_array, flag_array, pml_bin_dim, pml_bout_dim, n, xc, yc)
  use def_io
  use def_variables
  use def_vectors

  implicit none

  ! Input arguments
  integer :: pml_bin_dim, pml_bout_dim, n
  real(kind=8) :: k
  real(kind=8), dimension(n) :: x, y
  integer, dimension(n) :: boundary_array
  logical, dimension(n) :: flag_array
  ! Output arguments
  complex*16, dimension(n) :: xc, yc
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

do ii=1,n
    if (boundary_array(ii) == 2) then
        count1 = count1 + 1
        pmlbin_x(count1) = x(ii)
        pmlbin_y(count1) = y(ii)
    else if (boundary_array(ii) == 3) then
        count2 = count2 + 1
        pmlbout_x(count2) = x(ii)
        pmlbout_y(count2) = y(ii)
    endif
enddo

do ii=1,n
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


end subroutine lcpml