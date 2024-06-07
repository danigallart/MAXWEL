subroutine mesh_reader
    use def_io
    use def_variables
    use def_vectors
    
    implicit none
    
    integer :: istat, ii, jj, node
        
    double precision :: r
    
    double precision :: min_tol, max_tol
    logical, allocatable :: pml_flag(:)
        
    ! Node's coordinates
    
    READ(nodes_unit,*,IOSTAT=istat) delh,NP
    
    allocate(coorx(NP), coory(NP))
    do ii=1,NP
        read(nodes_unit,*,iostat=istat) coorx(ii),coory(ii)
    end do
            
    ! Conectivity matrix
    
    nodpel = 3
    
    read(elements_unit,*,iostat=istat) NE
    
    allocate(conn(NE,nodpel))
    DO ii=1,NE
        read(elements_unit,*,iostat=istat)  conn(ii,nodpel-2), conn(ii,nodpel-1), conn(ii,nodpel)
    END DO
        
    allocate(coorx_mid(NE),coory_mid(NE))
    allocate(material(NE), boundary(NP))
    allocate(pml_flag(NP))
    
    plasma_radius = r_scat
    free_space_dim = plasma_radius/2
    pmldim = plasma_radius/2
    huygdim = delh
    
    rpmlin = plasma_radius + free_space_dim
    rpmlout = rpmlin + pmldim
    rhuyg = plasma_radius + huygdim
    
    min_tol = 1 - boundary_tol
    max_tol = 1 + boundary_tol
    
    n_scatb = 0
    n_pml_bin = 0
    n_pml_bout = 0
    n_huygb = 0
    n_pml = 0
    
    pml_flag = .FALSE.
        
do ii=1,NP
    r = sqrt(coorx(ii)**2 + coory(ii)**2)
    
    if ((r < plasma_radius * max_tol) .and. (r > plasma_radius * min_tol)) then                       ! Boundary of scatterer
        
        boundary(ii) = 1
        n_scatb = n_scatb + 1
        
    elseif ((r < rhuyg * max_tol) .and. (r > rhuyg * min_tol)) then                                   ! Huygens surface
        
        boundary(ii) = 2
        n_huygb = n_huygb + 1
        
    elseif ((r < rpmlin * max_tol) .and. (r > rpmlin * min_tol)) then                                 ! Inner boundary of PML region
        
        boundary(ii) = 3
        n_pml_bin = n_pml_bin + 1
        
    elseif ((r < rpmlout * max_tol) .and. (r > rpmlout * min_tol)) then                               ! Outer boundary of PML region
        
        boundary(ii) = 4
        n_pml_bout = n_pml_bout + 1
        n_pml = n_pml + 1
        
    elseif ((r < rpmlout) .and. (r > rpmlin)) then                                                   ! PML node
        
        boundary(ii) = 5
        n_pml = n_pml + 1
    else                                                                                    ! No boundary node
        boundary(ii) = 0
    endif
enddo
    
do ii=1,NE
    coorx_mid(ii) = sum(coorx(conn(ii,:)))/size(coorx(conn(ii,:)))
    coory_mid(ii) = sum(coory(conn(ii,:)))/size(coory(conn(ii,:)))

    r = sqrt(coorx_mid(ii)**2 + coory_mid(ii)**2)
    
    if (r < plasma_radius) then                                                             ! Scatterer element
        material(ii) = 1
    elseif ((r < rpmlout) .and. (r > rpmlin)) then                                          ! PML element
        material(ii) = 3
        do jj = 1, nodpel
            node = conn(ii,jj)
            pml_flag(node) = (boundary(node) /= 4) .and. (boundary(node) /=3)
        enddo
    else                                                                                    ! Vacuum element
        material(ii) = 2
    endif
enddo

allocate(complex_coorx(NP), complex_coory(NP))

complex_coorx = cmplx(0.0,0.0)
complex_coory = cmplx(0.0,0.0)

coorx = coorx + major_radius
coory = coory
coorx_mid = coorx_mid + major_radius
coory_mid = coory_mid

call lcpml(coorx,coory,k0,boundary,pml_flag,n_pml_bin,n_pml_bout,NP,complex_coorx,complex_coory)



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
  real(kind=8) :: x0,y0
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
    if (boundary_array(ii) == 3) then
        count1 = count1 + 1
        pmlbin_x(count1)=coorx(ii)
        pmlbin_y(count1)=coory(ii)
    else if (boundary_array(ii) == 4) then
        count2 = count2 + 1
        pmlbout_x(count2)=coorx(ii)
        pmlbout_y(count2)=coory(ii)
    endif
enddo

do ii=1,n
    if ((boundary_array(ii) == 4) .or. flag_array(ii)) then
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