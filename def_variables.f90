module def_variables
    
	implicit none

    integer, parameter :: geo_type = 3 !1. plane stress, 2. plane strain, 3. axisymmetric 4. 3D 
    integer, parameter :: ndim = 2
    integer :: NE, NP
    integer :: Nodpel
    integer :: Ngauss
    
    INTEGER, parameter :: element_type = 1 !1. triangle 2.quadrilateral
    INTEGER, parameter :: pol_order = 1 !1. linear 2. quadratic 3. cubic
    
    integer :: Ndir = 2
    integer :: NONULL
        
    double precision :: tol_solver = 1e-9
    integer :: iter_solver = 300000
    
    INTEGER, parameter :: problem_type = 2 ! 1. stationary 2. time_implicit
    
    DOUBLE PRECISION :: dt = 1e-2, total_t = 30.0
    
    !integer :: Nmat = 1
    
    integer :: n_pml_bin
    integer :: n_pml_bout
    integer :: n_pml
    
    integer :: n_scatb
    integer :: n_scatin
    integer :: m_scatin
    integer :: m_scatb
    
    integer :: n_huygb
    integer :: m_huygb
    
	! Constants
    real(kind=8), parameter :: c0 = 3.0d8 ! m/sec, velocity of light in free space
	real(kind=8), parameter :: pi = 4*ATAN(1.) ! m/sec, velocity of light in free space
	complex*16, parameter :: ij = cmplx(0, 1)   ! sqrt(-1)   


    real(kind=8), parameter :: nu0 = 120.0 * pi ! ohm, intrinsic impedance of free space
    real(kind=8), parameter :: e0 = (1e-9) / (36.0 * pi)  ! F/m, permittivity of free space
    real(kind=8), parameter :: mu0 = 4.0 * pi * 1e-7 ! H/m, permeability of free space
    

    ! Input parameters
    real(kind=8) :: freq = 300.0  ! MHz, frequency
    real(kind=8) :: freq_hz ! Hz, frequency (converted from MHz)
    real(kind=8) :: lambda0 ! meter, wavelength
    real(kind=8) :: k0 ! 1/meter, wavenumber
    real(kind=8) :: omg ! rad/sec, radial frequency
    real(kind=8) :: phii = 0 !rad, angle of incident field
    character(len=6) :: pol = 'TE' 
    

    
    end module def_variables
