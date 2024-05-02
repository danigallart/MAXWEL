!  $2D_harmonic_main.f90 
!
!  FUNCTIONS:
!  $2D_harmonic_main - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: $2D_harmonic_main
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program $2D_harmonic_main
    
    USE def_io
    USE def_variables
    USE def_vectors

    implicit none
    
    INTEGER :: itime
    DOUBLE PRECISION :: time

    ! Convert frequency from MHz to Hz
    freq_hz = freq * 1e6

    ! Calculate wavelength and other parameters
    lambda0 = c0 / freq_hz
    k0 = 2.0 * pi / lambda0
    omg = 2.0 * pi * freq_hz
    
        CALL initialise()
        
        CALL mesh_reader()
        
        CALL sparse_logic()
        
        CALL assembly()
        
        CALL solver()
        
        if (.TRUE.) then
            CALL exit_writer()
        else
            CALL exit_nosolver() 
        endif
    
        CALL finalise()

    end program $2D_harmonic_main

