program cerebro
use def_variables
implicit none
integer :: i
double precision :: time


    write(6,*) ' inicio del programa cerebro'
    call inicio()
    write(6,*) ' leo datos'
    call lectura()
    

    if(nop_malla == 1 ) then
        write(6,*) ' Mallo el dominio generico'
        if(problema=='ELECTRO' .or. problema=='HELMOZT') then
            call mallo_cuadrado()
        else
            call mallagen()
        endif
        
    elseif(nop_malla==2) then
       if(archi_malla/='  ') then
           write(6,*) ' armo mallado desde los archivos! '
           if(problema=='CEREBRO') then
              call leemalla()
           else
              call leo_mallado()
           endif
       else
           write(6,*) 'el archivo con la malla esta vacio!'
           stop ' '
       endif
       
    elseif(nop_malla==3 .or. nop_malla==5) then
        write(6,*) ' Mallo el dominio electrodo'
        call malla3d()
    elseif(nop_malla==4) then
        write(6,*) ' Mallo el dominio bidimensional con electrodos= ',N_electro
        call malla_Bergues()
    elseif(nop_malla==6) then
        write(6,*) ' leo esfera thorio '
        call leo_thorio()
    elseif(nop_malla==7) then
        call malla_raton()
    endif

    write(6,*) 'armo o leo sistema'
    call sistema()
   
   if(problema=='BERGUES') then  ! electroos bidimensionales
      
           call poisson_2d()
          ! write(6,*) ' poisson2d  ',electropor,electro_no

           call bio_therma2d()
           call Ph_evol_2d()
           !call OH_evol_2d()
   
   elseif(problema=='BERGUES3D') then  ! electroos bidimensionales
         tiempo_sumado=0.0
         time=0.0
         pasos_tiempo=1

         do i=1,pasos_tiempo
            call poisson(i,time,histo_gen(i),histo_pot(i))
             write(6,*) ' calculo el potencial ',i,time
         enddo


   elseif(problema=='PAPA2D') then  

          write(6,*) 'Inicio loop temporal'
        
         tiempo_sumado=0.0
         time=0.0
         do i=1,pasos_tiempo

            time=time + histo_gen(i)

            call modifico_fuentes(i,time)

            if(histo_pot(i)/=0.0) then !  hay potencial

                write(6,*) ' calculo el potencial ',i,time

                call poisson_2d_time(i,time,histo_gen(i),histo_pot(i))

            endif
    
            write(6,*) ' resuelvo biothermal! ',i,time
            call bio_therma2d_time(i,time,histo_gen(i),histo_pot(i))
    
         enddo

    elseif(problema=='PAPA_NEW')  then ! papa con nuevo electrodo 3D

        write(6,*) 'Inicio loop temporal'
        
        tiempo_sumado=0.0
        time=0.0
        do i=1,pasos_tiempo

            time=time + histo_gen(i)
        
            if(histo_pot(i)/=0.0) then !  hay potencial

                write(6,*) ' calculo el potencial ',i,time

                call poisson(i,time,histo_gen(i),histo_pot(i))

            endif
    
            write(6,*) ' resuelvo biothermal! ',i,time
            call biotherma(i,time,histo_gen(i),histo_pot(i))
    
        enddo
    elseif(problema=='THORIO')  then ! thorio 3D

        write(6,*) 'resuelvo tempertarua Thorio'
        
        call thorio_temper()
        
    elseif(problema=='ELECTRO')  then ! electrostatic 2D
        
        
        call poisson_2d_electro

    elseif(problema=='HELMOZT')  then ! Helmozt 2D
        
        
        call Helmozt_2d
        
   else  ! CEREBRO o ELECTRODOS o ratones

        write(6,*) 'Inicio loop temporal'
        
        tiempo_sumado=0.0
        time=0.0
        do i=1,pasos_tiempo

            time=time + histo_gen(i)
        
            if(histo_pot(i)/=0.0) then !  hay potencial

                write(6,*) ' calculo el potencial ',i,time

                call poisson(i,time,histo_gen(i),histo_pot(i))

            endif
    
            write(6,*) ' resuelvo biothermal! ',i,time
            call biotherma(i,time,histo_gen(i),histo_pot(i))
    
        enddo
   endif      

    write(6,*) ' cierro todo y me voy'
    call finalice()

    write(6,*) ' chau!'

    
end program cerebro