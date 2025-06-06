module def_variables
    
	implicit none
    
    character(10) :: problem                                                                        !For future implementations where there are more than one available problem to solve
    character(4)  :: reader_type = 'read'
    character*(4) :: elem_type, elem_shape
    character*(3) :: boundary_type
    
    integer :: ndim = 2                                                                             !Number of physical dimensions
    integer :: NE, NP                                                                               !NE: Number of elements, NP: Number of nodes
    integer :: Nodpel                                                                               !Nodes per element
    integer :: nodpedge, num_bnode_e                                                                !Nodes per edge of the element
    integer :: Ngauss                                                                               !Number of gaussian integration points
    integer :: nboun                                                                                !Number of boundaries
    
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
    integer :: n_huygb
    
    double precision :: delh
        
    complex*16 :: rel_permitivity_xx, rel_permitivity_xy, &
                  rel_permitivity_yx, rel_permitivity_yy, &
                  rel_permitivity_zz

    complex*16 :: rel_permeability_xx, rel_permeability_xy, &
                  rel_permeability_yx, rel_permeability_yy, &
                  rel_permeability_zz
    
    complex*16 :: epsilon_scat_xx, epsilon_scat_xy, &
                  epsilon_scat_yx, epsilon_scat_yy, &
                  epsilon_scat_zz

    complex*16 :: mu_scat_xx, mu_scat_xy, &
                  mu_scat_yx, mu_scat_yy, &
                  mu_scat_zz

    double precision :: cond, im_rel
    
    complex*16 :: pxxe,pxye,pyxe,pyye,qe
    
	! Constants
    real(kind=8), parameter :: c0 = 299792458.0                                                     ! m/s, velocity of light in free space
	real(kind=8), parameter :: pi = 4.0*ATAN(1.)                                                    ! pi is pi, constant of nature
	complex*16, parameter :: ij = cmplx(0, 1)                                                       ! Imaginary unit, sqrt(-1)


    real(kind=8), parameter :: nu0 = 376.730313412                                                  ! ohm, intrinsic impedance of free space
    real(kind=8), parameter :: e0 = 8.8542e-12                                                      ! F/m, permittivity of free space
    real(kind=8), parameter :: mu0 = 12.566370e-7                                                   ! H/m, permeability of free space
    

    ! Input parameters
    real(kind=8) :: freq                                                                            ! MHz, frequency
    real(kind=8) :: freq_hz                                                                         ! Hz, frequency (converted from MHz)
    real(kind=8) :: lambda0                                                                         ! meter, wavelength
    real(kind=8) :: k0                                                                              ! 1/meter, wavenumber
    real(kind=8) :: omg                                                                             ! rad/sec, radial frequency
    real(kind=8) :: phii                                                                            ! rad, angle of incident field
    real(kind=8) :: r_scat                                                                          ! radius of scaterer in units of lambda0
    real(kind=8) :: density_e_0                                                                     ! Central plasma electron density
    character(len=2) :: pol                                                                         ! TE: Transversal electric, TM: Transversal magnetic
    character(len=1) :: read_logic                                                                  ! Reads logic (AD,AN), Yes or No
    character(len=1) :: system_sym                                                                  ! Asumes symmetric system, Yes or No
    character(len=1) :: antenna_source                                                              ! Radiation from antenna, Yes(On) or No(Off)
    character(len=1) :: plane_wave_source                                                           ! Radiation from plane wave, Yes(On) or No(Off)
    
    integer :: plasma, density_flag, magnetic_flag                                                  ! Plasma flag
    integer :: n_species
    double precision :: plasma_freq                                                                 ! Plasma frequency
    double precision :: cyclo_freq                                                                  ! Cyclotron frequency for electrons, deuterons and tritons
    double precision :: mass1, mass2, mass3, mass4                                                  ! Mass for electrons, deuterons and tritons
    double precision :: deuterium_frac                                                              ! Fraction of deuterium, n_d/(n_e)
    double precision :: tritium_frac                                                                ! Fraction of tritium, n_t/(n_e)
    double precision :: helium_3_frac                                                               ! Fraction of helium-3, n_he3/(n_e)
    double precision :: ka,aa                                                                       ! Parameters of density function
    
    double precision :: mag_field_0                                                                 ! Axial magnetic field
    double precision :: major_radius                                                                ! Tokamak major radius                                                   
    double precision, parameter :: e_charge = 1.60217662e-19                                        ! Elementary charge in Coulombs

    
    complex*16 :: current_density1,current_density2,dummy_current                                   ! A/m^2, current density
    
    double precision :: plasma_radius, free_space_dim, pmldim, huygdim, &
                        rpmlin, rpmlout, rhuyg
    
    
    end module def_variables
