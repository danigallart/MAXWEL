module def_bio



! temperaturas fijas
double precision :: Tblod = 37.0,   & !37 °C   25 para la papa!!!
       &            Tout  = 37.0,   &    !37 °C 
       &            h_b   =0.010,    &  ! W mm-1 K-1
       &            pulso = 1.0E-2, &  !E-4 seg
       &            frec = 1,      & ! hertz
       &            W_b =0.0,   &    !7.1e-3,       &  ! 1/s 7.1e-3           0.0 para la papa!!
       &            C_b = 3.8400,     &  ! J/gK 
       &            rho_b = 0.0010600,    &   !g/mm3
       &            k_c=0.015,             &
       &            rho_c = 0.0079,       &
       &            Cp_c = 0.5,        &
       &            alfa_b= 0.135 ! 0.032        ! 1/°C
end module def_bio