module def_variables
    
	implicit none
    
    character(10) :: problem                                                                        !For future implementations where there are more than one available problem to solve
    character(20) :: file_mesh = ''                                                                 !Name of the file where the mesh is stored
    
    !integer, parameter :: geo_type = 3                                                             !1. plane stress, 2. plane strain, 3. axisymmetric 4. 3D 
    integer, parameter :: ndim = 2                                                                  !Number of physical dimensions
    integer :: NE, NP                                                                               !NE: Number of elements, NP: Number of nodes
    integer :: Nodpel                                                                               !Nodes per element
    integer :: Ngauss                                                                               !Number of gaussian integration points
    
    !INTEGER, parameter :: element_type = 1                                                         !1. triangle 2.quadrilateral
    !INTEGER, parameter :: pol_order = 1                                                            !1. linear 2. quadratic 3. cubic
    
    !integer :: Ndir = 2
    integer :: NONULL                                                                               !Number of non-null matrix elements
        
    double precision :: tol_solver = 1e-9
    integer :: iter_solver = 300000
    double precision :: boundary_tol
    
    !INTEGER, parameter :: problem_type = 2                                                         ! 1. stationary 2. time_implicit
    
    !DOUBLE PRECISION :: dt = 1e-2, total_t = 30.0
    
    !integer :: Nmat = 1
    
    integer :: n_pml_bin                                                                            !PML inside boundary nodes
    integer :: n_pml_bout                                                                           !PML outside boundary nodes
    integer :: n_pml                                                                                !PML nodes
    integer :: n_scatb                                                                              !Scatering boundary nodes
    integer :: n_huygb
    
    double precision :: delh
    
	! Constants
    real(kind=8), parameter :: c0 = 3.0d8                                                           ! m/sec, velocity of light in free space
	real(kind=8), parameter :: pi = 4*ATAN(1.)                                                      ! pi is pi, constant of nature
	complex*16, parameter :: ij = cmplx(0, 1)                                                       ! Imaginary unit, sqrt(-1)


    real(kind=8), parameter :: nu0 = 120.0 * pi                                                     ! ohm, intrinsic impedance of free space
    real(kind=8), parameter :: e0 = (1e-9) / (36.0 * pi)                                            ! F/m, permittivity of free space
    real(kind=8), parameter :: mu0 = 4.0 * pi * 1e-7                                                ! H/m, permeability of free space
    

    ! Input parameters
    real(kind=8) :: freq                                                                            ! MHz, frequency
    real(kind=8) :: freq_hz                                                                         ! Hz, frequency (converted from MHz)
    real(kind=8) :: lambda0                                                                         ! meter, wavelength
    real(kind=8) :: k0                                                                              ! 1/meter, wavenumber
    real(kind=8) :: omg                                                                             ! rad/sec, radial frequency
    real(kind=8) :: phii                                                                            ! rad, angle of incident field
    real(kind=8) :: r_scat                                                                          ! radius of scaterer in units of lambda0
    character(len=2) :: pol                                                                         ! TE: Transversal electric, TM: Transversal magnetic
    

    
    end module def_variables
