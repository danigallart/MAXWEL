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
    
        CALL initialise()
        
        CALL reader()

    ! Convert frequency from MHz to Hz
    freq_hz = freq * 1e6

    ! Calculate wavelength and other parameters
    lambda0 = c0 / freq_hz
    k0 = 2.0 * pi / lambda0
    omg = 2.0 * pi * freq_hz
        
        print*, "Mesh reader"
        if (.FALSE.) then
            CALL mesh_reader()
        else
            CALL mesh_reader_tokamak()
        endif
        
        print*, "Sparse logic"
        
        CALL sparse_logic()
        
        print*, 'Assembly'
        
        CALL assembly()
        
        print*, 'Solver'
        
        if (.TRUE.) then
            CALL solver()
        else
            CALL read_solution()
        endif
        
        print*, 'Derivatives'
        
        CALL derivatives()
        
        print*, 'Exit'
        
        if (.TRUE.) then
            CALL exit_writer() !Exit data file with solution and coordinates
        else
            CALL exit_nosolver() !Exit data file with linear equation system
        endif
    
        CALL finalise()

    end program $2D_harmonic_main

