subroutine mesh_reader
    use def_io
    use def_variables
    use def_vectors
    
    implicit none
    
    integer :: istat, ii, jj
        
    double precision :: plasma_radius, fsdim, pmldim, huygdim, &
                        rpmlin, rpmlout, rhuyg
    double precision :: r
    
    double precision :: min_tol, max_tol
    
    
    double precision, allocatable :: pmlbin_coorx(:),pmlbin_coory(:), pmlbout_coorx(:), pmlbout_coory(:)
    
    ! Node's coordinates
    
    READ(nodes_unit,*,IOSTAT=istat) delh,NP
    
    allocate(coorx(NP), coory(NP))
    do ii=1,NP
        read(nodes_unit,*,iostat=istat) coorx(ii),coory(ii)
    end do
    
    CLOSE(nodes_unit)
        
    ! Conectivity matrix
    
    nodpel = 3
    
    read(elements_unit,*,iostat=istat) NE
    
    allocate(conn(NE,nodpel))
    DO ii=1,NE
        read(elements_unit,*,iostat=istat)  conn(ii,nodpel-2), conn(ii,nodpel-1), conn(ii,nodpel)
    END DO
    
    CLOSE(elements_unit)
    
    ! PML region and boundary nodes
    
    read(pml_bin_nodes_unit,*,iostat=istat) n_pml_bin
    
    allocate(pml_bin_nodes(n_pml_bin))
    
    do ii=1,n_pml_bin
        read(pml_bin_nodes_unit,*,iostat=istat) pml_bin_nodes(ii)
    end do
    
    close(pml_bin_nodes_unit)
    
    read(pml_bout_nodes_unit,*,iostat=istat) n_pml_bout
    
    allocate(pml_bout_nodes(n_pml_bout))
    
    do ii=1,n_pml_bout
        read(pml_bout_nodes_unit,*,iostat=istat) pml_bout_nodes(ii)
    end do
    
    close(pml_bout_nodes_unit)
    
    read(pml_nodes_unit,*,iostat=istat) n_pml
    
    allocate(pml_nodes(n_pml))
    
    do ii=1,n_pml
        read(pml_nodes_unit,*,iostat=istat) pml_nodes(ii)
    end do
    
    close(pml_nodes_unit)
    
    ! Scattering region and boundaries
    
    read(scatb_nodes_unit,*,iostat=istat) n_scatb
    
    allocate(scatb_nodes(n_scatb))
    
    do ii=1,n_scatb
        read(scatb_nodes_unit,*,iostat=istat) scatb_nodes(ii)
    end do
    
    close(scatb_nodes_unit)
    
    read(scatin_nodes_unit,*,iostat=istat) n_scatin
    
    allocate(scatin_nodes(n_scatin))
    
    do ii=1,n_scatin
        read(scatin_nodes_unit,*,iostat=istat) scatin_nodes(ii)
    end do
    
    close(scatin_nodes_unit)
    
    read(scatin_elements_unit,*,iostat=istat) m_scatin
    
    allocate(scatin_elements(m_scatin))
    
    do ii=1,m_scatin
        read(scatin_elements_unit,*,iostat=istat) scatin_elements(ii)
    end do
    
    close(scatin_elements_unit)
    
    read(scatb_elements_unit,*,iostat=istat) m_scatb
    
    allocate(scatb_elements(m_scatb))
    
    do ii=1,m_scatb
        read(scatb_elements_unit,*,iostat=istat) scatb_elements(ii)
    end do
    
    close(scatb_elements_unit)

    ! Huygens region
    
    read(huygb_nodes_unit,*,iostat=istat) n_huygb
    
    allocate(huygb_nodes(n_huygb))
    
    do ii=1,n_huygb
        read(huygb_nodes_unit,*,iostat=istat) huygb_nodes(ii)
    end do
    
    close(huygb_nodes_unit) 
    
    
    read(huygb_elements_unit,*,iostat=istat) m_huygb
    
    allocate(huygb_elements(m_huygb))
    
    do ii=1,m_huygb
        read(huygb_elements_unit,*,iostat=istat) huygb_elements(ii)
    end do
    
    close(huygb_elements_unit) 
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    allocate(coorx_mid(NE),coory_mid(NE))
    
    do ii=1,NE
    coorx_mid(ii) = sum(coorx(conn(ii,:)))/size(coorx(conn(ii,:)))
    coory_mid(ii) = sum(coory(conn(ii,:)))/size(coory(conn(ii,:)))
    enddo
    
    allocate(material(NE), boundary(NP))
    
    plasma_radius = r_scat * lambda0
    fsdim = lambda0/2
    pmldim = lambda0/2
    huygdim = delh
    
    rpmlin = plasma_radius + fsdim
    rpmlout = rpmlin + pmldim
    rhuyg = plasma_radius + huygdim
    
    min_tol = 1 - boundary_tol
    max_tol = 1 + boundary_tol
    
    n_scatb = 0
    n_pml_bin = 0
    n_pml_bout = 0
    n_huygb = 0
    n_pml = 0
        
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
    r = sqrt(coorx_mid(ii)**2 + coory_mid(ii)**2)
    if (r < plasma_radius) then                                                             ! Scatterer element
        material(ii) = 1
    elseif ((r < rpmlout) .and. (r > rpmlin)) then                                          ! PML element
        material(ii) = 3
    else                                                                                    ! Vacuum element
        material(ii) = 2
    endif
enddo

allocate(pmlbin_coorx(n_pml_bin), pmlbin_coory(n_pml_bin))
allocate(pmlbout_coorx(n_pml_bout), pmlbout_coory(n_pml_bout))


do ii=1,n_pml_bin
    pmlbin_coorx(ii)=coorx(pml_bin_nodes(ii))
    pmlbin_coory(ii)=coory(pml_bin_nodes(ii))
end do


do ii=1,n_pml_bout
    pmlbout_coorx(ii)=coorx(pml_bout_nodes(ii))
    pmlbout_coory(ii)=coory(pml_bout_nodes(ii))
end do

allocate(complex_coorx(NP), complex_coory(NP))

do ii=1,NP
    complex_coorx(ii)=cmplx(coorx(ii),0.0)
    complex_coory(ii)=cmplx(coory(ii),0.0)
end do

do ii=1,n_pml
    jj = pml_nodes(ii)
    call lcpml(coorx(jj),coory(jj),k0,pmlbin_coorx,pmlbin_coory,pmlbout_coorx,pmlbout_coory,complex_coorx(jj),complex_coory(jj))
    !complex_coorx(jj)=cmplx(x_rval,x_cval)
    !complex_coory(jj)=cmplx(y_rval,y_cval)
end do
        


    end subroutine mesh_reader
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! LCPML function for Fortran 90
subroutine lcpml(x, y, k, pmlbin_x, pmlbin_y, pmlbout_x, pmlbout_y, xc, yc)
  use def_io
  use def_variables
  use def_vectors

  implicit none

  ! Input arguments
  real(kind=8) :: x, y, k
  real(kind=8), dimension(n_pml_bin) :: pmlbin_x, pmlbin_y
  real(kind=8), dimension(n_pml_bout) :: pmlbout_x, pmlbout_y
  ! Output arguments
  complex*16 :: xc, yc
  !double precision, intent(out) :: xc_r, yc_r,xc_im, yc_im

  ! LC-PML parameters
  !real(kind=8) :: alpha
  real(kind=8),dimension(n_pml_bin) :: dpml1
  integer :: m, ind
  complex*16 :: alphajk, term
  real(kind=8) :: x0,y0
  double precision :: vpx(n_pml_bout),vpy(n_pml_bout)
  double precision :: npx(n_pml_bout),npy(n_pml_bout)
  double precision :: lp(n_pml_bout)
  double precision :: vx, vy, l, nx, ny, x1, y1,dpml2
  double precision :: ksi
  

  ! Set LC-PML parameters
  !alpha = 7.0 * k
  !alphajk = -7.0
  alphajk = cmplx(0,-7.0)
  m = 3  ! PML decay rate (integer 2 or 3)

  ! Find the point on the inner PML boundary nearest to point P
  dpml1 = sqrt((pmlbin_x - x)**2 + (pmlbin_y - y)**2)
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

  vx = x - x0
  vy = y - y0
  l = sqrt(vx**2 + vy**2)
  nx = vx / l  ! Unit vector from r0 to r (x-comp)
  ny = vy / l  ! Unit vector from r0 to r (y-comp)

  if (l < 1.0e-8) then
    xc = cmplx(x,0.0)
    yc = cmplx(y,0.0)
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
    xc = cmplx(x,0.0) + term * nx
    yc = cmplx(y, 0.0) + term * ny
  end if

end subroutine lcpml