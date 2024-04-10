subroutine lectura()
use def_variables
implicit none
character*(1)::Z 
character*(120):: textinput,opcion
integer :: leng,last,suma,kk,pasos_aux,j,npas,naux
double precision :: pas

write(Z,'(A1)') Z'09'

do while(textinput /= 'END_DATA')
     read(unit_data,'(A120)') textinput
	 call upcase(textinput)
	 leng=len_trim(textinput) 
     
	 last=0
	 suma=0
	 opcion='  '
	 do while(last<leng)
        last=last+1
		if(textinput(last:last)/=' '.and.textinput(last:last)/=':'.and.textinput(last:last)/='='.and.textinput(last:last)/=Z) then
		  if(textinput(last:last)=='#') then
		    last=leng
		  else
		    suma=suma+1
		    opcion(suma:suma)= textinput(last:last)
		  endif
		endif

	 enddo
     
	 leng=len_trim(opcion)
	   last=0	 
	   do while(last<leng)
			last=last+1
			if(opcion(1:last)=='PROBLEMA') then
				read(opcion(last+1:leng),'(a10)') problema
				write(unit_cont,*) ' PROBLEMA: ', problema
				last=leng+1
			endif
        enddo
	   last=0	 
	   do while(last<leng)
			last=last+1
			if(opcion(1:last)=='NOMBRE') then
				read(opcion(last+1:leng),'(a10)') nombre
				write(unit_cont,*) ' NOMBRE: ', nombre
				last=leng+1
			endif
        enddo
       last=0	 
	   do while(last<leng)
			last=last+1
			if(opcion(1:last)=='MALLA') then
				read(opcion(last+1:leng),'(i3)') nop_malla
				write(unit_cont,*) ' MALLA: ', nop_malla
				last=leng+1
			endif

        enddo
	 last=0	 
	   do while(last<leng)
			last=last+1
			if(opcion(1:last)=='ARCHIVO1') then
				read(opcion(last+1:leng),'(a20)') archi_malla
				write(unit_cont,*) ' ARCHIVO1: ', archi_malla
				last=leng+1
			endif
        enddo
	 last=0	 
	   do while(last<leng)
			last=last+1
			if(opcion(1:last)=='ARCHIVO2') then
				read(opcion(last+1:leng),'(a20)') archi_tumor
				write(unit_cont,*) ' ARCHIVO2: ', archi_tumor
				last=leng+1
			endif
        enddo
       last=0	 
	   do while(last<leng)
			last=last+1
			if(opcion(1:last)=='OPCION') then
				read(opcion(last+1:leng),'(i3)') nopcion
				write(unit_cont,*) ' OPCION: ', nopcion
				last=leng+1
			endif
        enddo
        
         last=0	 
	   do while(last<leng)
			last=last+1
			if(opcion(1:last)=='ARCHIVO3') then
				read(opcion(last+1:leng),'(a20)') archi_sistema
				write(unit_cont,*) ' ARCHIVO3: ', archi_sistema
				last=leng+1
			endif
       enddo
	   last=0	 
	   do while(last<leng)
			last=last+1
			if(opcion(1:last)=='LADOX') then
				read(opcion(last+1:leng),'(e5.0)') tam_zonex
				write(unit_cont,*) ' LADOX: ', tam_zonex
				last=leng+1
			endif
        enddo
	   last=0	 
	   do while(last<leng)
			last=last+1
			if(opcion(1:last)=='LADOY') then
				read(opcion(last+1:leng),'(e5.0)') tam_zoney
				write(unit_cont,*) ' LADOY: ', tam_zoney
				last=leng+1
			endif
        enddo
	   last=0	 

       
       
       last=0	 
	   do while(last<leng)
			last=last+1
			if(opcion(1:last)=='TUMORX') then
				read(opcion(last+1:leng),'(e5.0)') tam_zonex
				write(unit_cont,*) ' TUMORX: ', tam_zonex
				last=leng+1
			endif
        enddo
	   last=0	 
	   do while(last<leng)
			last=last+1
			if(opcion(1:last)=='TUMORY') then
				read(opcion(last+1:leng),'(e5.0)') tam_zoney
				write(unit_cont,*) ' TUMORY: ', tam_zoney
				last=leng+1
			endif
        enddo
	   last=0	 
	   do while(last<leng)
			last=last+1
			if(opcion(1:last)=='TUMORZ') then
				read(opcion(last+1:leng),'(e5.0)') tam_zonez
				write(unit_cont,*) ' TUMORZ: ', tam_zonez
				last=leng+1
			endif
        enddo
         last=0	 
	   do while(last<leng)
			last=last+1
			if(opcion(1:last)=='TAMA_EL') then
				read(opcion(last+1:leng),'(e5.0)') tamano_elem
				write(unit_cont,*) ' TAMA_EL: ', tamano_elem
				last=leng+1
			endif
        enddo

	   last=0	 
	   do while(last<leng)
			last=last+1
			if(opcion(1:last)=='ZONE1') then
				read(opcion(last+1:leng),'(e5.0)') rad_zone1
				write(unit_cont,*) ' ZONE1: ', rad_zone1
				last=leng+1
			endif
        enddo
        last=0	 
	   do while(last<leng)
			last=last+1
			if(opcion(1:last)=='ZONE2') then
				read(opcion(last+1:leng),'(e5.0)') rad_zone2
				write(unit_cont,*) ' ZONE2: ', rad_zone2
				last=leng+1
			endif
        enddo
         last=0	 
	   do while(last<leng)
			last=last+1
			if(opcion(1:last)=='PAPA_SIG') then
				read(opcion(last+1:leng),'(e15.5)') sigma3
				write(unit_cont,*) ' PAPA_SIG: ', sigma3
				last=leng+1
			endif
        enddo
        last=0	 
         last=0	 
	   do while(last<leng)
			last=last+1
			if(opcion(1:last)=='AISLA_SIG') then
				read(opcion(last+1:leng),'(e15.5)') sigma1
				write(unit_cont,*) ' AISLA_SIG: ', sigma1
				last=leng+1
			endif
        enddo
        last=0	 
         last=0	 
	   do while(last<leng)
			last=last+1
			if(opcion(1:last)=='ELEC_SIG') then
				read(opcion(last+1:leng),'(e15.5)') sigma2
				write(unit_cont,*) ' ELEC_SIG: ', sigma2
				last=leng+1
			endif
        enddo
        last=0	 


         last=0	 
	   do while(last<leng)
			last=last+1
			if(opcion(1:last)=='SIGMAT') then
				read(opcion(last+1:leng),'(e15.5)') sigma1
				write(unit_cont,*) ' SIGMAT: ', sigma1
				last=leng+1
			endif
        enddo
        last=0	 
	   do while(last<leng)
			last=last+1
			if(opcion(1:last)=='SIGMAN') then
				read(opcion(last+1:leng),'(e15.5)') sigma2
				write(unit_cont,*) ' SIGMAN: ', sigma2
				last=leng+1
			endif
        enddo
        last=0	 
	   do while(last<leng)
			last=last+1
			if(opcion(1:last)=='SIGMAL') then
				read(opcion(last+1:leng),'(e15.5)') sigma3
				write(unit_cont,*) ' SIGMAL: ', sigma3
				last=leng+1
			endif
        enddo
        last=0	 
	   do while(last<leng)
			last=last+1
			if(opcion(1:last)=='SIGMAE') then
				read(opcion(last+1:leng),'(e15.5)') sigma4
				write(unit_cont,*) ' SIGMAE: ', sigma4
				last=leng+1
			endif
        enddo
        last=0	 
	   do while(last<leng)
			last=last+1
			if(opcion(1:last)=='POTENCIAL') then
				read(opcion(last+1:leng),'(e5.0)') potencial
				write(unit_cont,*) ' POTENCIAL: ', potencial
				last=leng+1
			endif
        enddo
        last=0	 
	   do while(last<leng)
			last=last+1
			if(opcion(1:last)=='TIERRA') then
				read(opcion(last+1:leng),'(e5.0)') tierra
				write(unit_cont,*) ' TIERRA: ', tierra
				last=leng+1
			endif
        enddo
        last=0	 
	   do while(last<leng)
			last=last+1
			if(opcion(1:last)=='TIPOE') then
				read(opcion(last+1:leng),'(i3)') tipo_elec
				write(unit_cont,*) ' TIPOE: ', tipo_elec
				last=leng+1
			endif
        enddo
       last=0	 
	   do while(last<leng)
			last=last+1
			if(opcion(1:last)=='RADIOS') then
				read(opcion(last+1:leng),'(e15.5,e15.5,e15.5,e15.5)') rad1,rad2,rad3,rad4
				write(unit_cont,*) ' radios: ', rad1,rad2,rad3,rad4
				last=leng+1
			endif
        enddo

         last=0	 
	   do while(last<leng)
			last=last+1
			if(opcion(1:last)=='ELEC_ALT') then
				read(opcion(last+1:leng),'(i3)') elec_alt
				write(unit_cont,*) ' ELEC_ALT: ', elec_alt
				last=leng+1
			endif
        enddo

         last=0	 
	   do while(last<leng)
			last=last+1
			if(opcion(1:last)=='ELEC_SEP') then
				read(opcion(last+1:leng),'(i3)') elec_sep
				write(unit_cont,*) ' ELEC_SEP: ', elec_sep
				last=leng+1
			endif
        enddo

       last=0	 
	   do while(last<leng)
			last=last+1
			if(opcion(1:last)=='ALTURAS') then
				read(opcion(last+1:leng),'(e15.5,e15.5,e15.5)') alt1,alt2,alt3
				write(unit_cont,*) ' alturas: ', alt1,alt2,alt3
				last=leng+1
			endif
        enddo
         last=0	 
	   do while(last<leng)
			last=last+1
			if(opcion(1:last)=='TIEMPO') then
				read(opcion(last+1:leng),'(E15.5)') tiempo_total
				write(unit_cont,*) ' TIEMPO: ', tiempo_total
				last=leng+1
			endif
        enddo

         last=0	 
	   do while(last<leng)
			last=last+1
			if(opcion(1:last)=='PASO_TI1') then
				read(opcion(last+1:leng),'(E15.5)') pas_ti1
				write(unit_cont,*) ' PASO_TI1: ', pas_ti1
				last=leng+1
			endif
        enddo
       
         last=0	 
	   do while(last<leng)
			last=last+1
			if(opcion(1:last)=='PASO_TI2') then
				read(opcion(last+1:leng),'(E15.5)') pas_ti2
				write(unit_cont,*) ' PASO_TI2: ', pas_ti2
				last=leng+1
			endif
        enddo
       
         last=0	 
	   do while(last<leng)
			last=last+1
			if(opcion(1:last)=='HISTORIA') then
				read(opcion(last+1:leng),'(I3)') pasos_tiempo
				write(unit_cont,*) ' HISTORIA: ', pasos_tiempo
				last=leng+1
                allocate(historiat(pasos_tiempo))
                naux=pasos_tiempo/tiempo_total
                do kk=1,naux
                   read(unit_data,*) historiat(kk)
                   write(unit_cont,*) kk,historiat(kk)
				enddo
                do kk=naux+1,pasos_tiempo
                   if(mod(kk,2)==0) then
                      historiat(kk) = historiat(2)
                   else
                      historiat(kk) = historiat(1)
                   endif
                enddo

			endif
        enddo
        
        last=0	 
	    do while(last<leng)
			last=last+1
			if(opcion(1:last)=='NELECTRO') then
				read(opcion(last+1:leng),'(i3)') N_electro
				write(unit_cont,*) ' NELECTRO: ', N_electro
				last=leng+1
			endif
        enddo
        
        last=0	 
	    do while(last<leng)
			last=last+1
			if(opcion(1:last)=='DISELE') then
				read(opcion(last+1:leng),'(E15.5)') dis_ele
				write(unit_cont,*) ' DIS_ELE: ', dis_ele
				last=leng+1
			endif
        enddo
        
        last=0	 
	    do while(last<leng)
			last=last+1
			if(opcion(1:last)=='EXENTRIC') then
				read(opcion(last+1:leng),'(E15.5)') exentric
				write(unit_cont,*) ' EXENTRIC: ', exentric
				last=leng+1
			endif
        enddo
        last=0	 
	    do while(last<leng)
			last=last+1
			if(opcion(1:last)=='MODO') then
				read(opcion(last+1:leng),'(I4)') modo
				write(unit_cont,*) ' MODO: ', modo
				last=leng+1
			endif
        enddo
        last=0	 
	    do while(last<leng)
			last=last+1
			if(opcion(1:last)=='ALTERNA') then
				read(opcion(last+1:leng),'(E15.5)') Alterna
				write(unit_cont,*) ' Alterna: ', Alterna
				last=leng+1
			endif
        enddo

   enddo
   

call seteo()

! armo historia

if( pasos_tiempo>0) then
    pasos_aux=0
    
    
    
    do kk=1,pasos_tiempo
      if(mod(kk,2)==0) then ! pasos pares
         npas=  historiat(kk)/pas_ti2

         pasos_aux =  pasos_aux + npas
      else
         npas=  historiat(kk)/pas_ti1
         pasos_aux =  pasos_aux + npas
      endif
    enddo

    allocate(histo_gen(pasos_aux),histo_pot(pasos_aux))
    pasos_aux=0

    do kk=1,pasos_tiempo
     if(mod(kk,2)==0) then ! pasos pares
         npas=  historiat(kk)/pas_ti2
     
         pas = historiat(kk)/npas

         do j=1,npas
            pasos_aux=pasos_aux+1
            histo_gen(pasos_aux)=pas
            histo_pot(pasos_aux)=0
         enddo   
      else
         npas=  historiat(kk)/pas_ti1
         pas = historiat(kk)/npas
         do j=1,npas
            pasos_aux=pasos_aux+1
            histo_gen(pasos_aux)=pas
            histo_pot(pasos_aux)=1
         enddo   
      endif

    enddo


    pasos_tiempo=pasos_aux

    write(unit_cont,*) 'Historia de potencia general ',pasos_tiempo
    do kk=1,pasos_tiempo

       write(unit_cont,*) kk,histo_gen(kk),histo_pot(kk)

    enddo
endif

close(unit_data)


! creo directorio para salidas




end subroutine lectura

subroutine lee_sistema(archi_sistema,nnodes)
use def_solver
implicit none
character(20)::archi_sistema
integer :: nnodes
!local
integer::kk,n,j
integer :: unit_sist=1111

open(unit=unit_sist,file=archi_sistema)


read(unit_sist,*) n
if(n/=nnodes) then
   write(6,*) 'sistema incorrecto!!',n
   stop ' '
endif

allocate(masa(nnodes))

do kk=1,n
   read(unit_sist,*) j,masa(kk)
enddo

read(unit_sist,*) n
allocate(ia(nnodes+1),cx(nnodes+1))
if(n/=nnodes+1) then
   write(6,*) 'sistema incorrecto!! '
   stop ' '
endif

do kk=1,n
   read(unit_sist,*) j,ia(kk),cx(kk)
enddo

read(unit_sist,*) nonull
allocate(ad(nnodes),rhs(nnodes),rhs_ant(nnodes),solucion(nnodes),solucion_ant(nnodes))
rhs_ant=0.0
rhs=0.0
solucion_ant=0.0
solucion=0.0

allocate(ja(nonull),an(nonull))
do kk=1,nonull
   read(unit_sist,*) j,ja(kk)
enddo
close(unit_sist)

end subroutine lee_sistema

subroutine upcase(word)
    !-----------------------------------------------------------------------
    !
    !     This routine converts wopos to upper case 
    !
    !-----------------------------------------------------------------------
    implicit none
    character(*), intent(inout) :: word
    integer                     :: iposi,ioctv

    do iposi=1,60                                   ! process all positions
       ioctv=ichar(word(iposi:iposi))               ! octal value
       if(o'141'<=ioctv.and.o'172'>=ioctv) then ! it is a lower case
          ioctv=ioctv-o'40'                          ! equivalent upper case
          word(iposi:iposi)=char(ioctv)              ! convert it to upcase
       end if
    end do ! iposi=1,5

  end subroutine upcase


