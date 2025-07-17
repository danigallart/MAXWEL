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
        
        CALL open_files()

    ! Convert frequency from MHz to Hz
    freq_hz = freq * 1e6

    ! Calculate wavelength and other parameters
    lambda0 = c0 / freq_hz
    k0 = 2.0 * pi / lambda0
    omg = 2.0 * pi * freq_hz
    
        print*, "Mesh reader"
        if (reader_type == 'read') then
            ! Legacy mesh reader, all files should have Alya problem type, hence reader_type is 'toka'
            CALL mesh_reader()
        elseif (reader_type == 'toka') then
            !Alya problem type can be downloaded from Alya repository and installed in GiD
            CALL mesh_reader_tokamak()
        endif
        
        print*, "Sparse logic"
        ! Sparse logic for scalar nodal system
        CALL sparse_logic()
        
        print*, 'Assembly'
        ! The computation of the global system matrix and right-hand side vector and its assembly
        CALL assembly()
        
        print*, 'Solver'
        ! Bi-conjugated gradient solver for linear equation system
        CALL solver()
        
        print*, 'Derivatives'
        ! Calculation of inplane fields, either Ex, Ey or Hx, Hy
            CALL derivatives()
        
        print*, 'Exit'
        
        if (.TRUE.) then
            CALL exit_writer() !Exit data file with solution and coordinates
        else
            CALL exit_nosolver() !Exit data file with linear equation system
        endif
    
        CALL finalise()
        
    end program $2D_harmonic_main

