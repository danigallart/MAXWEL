subroutine reader()
use def_io
use def_variables
use def_vectors
implicit none
character*(1)::Z 
character*(120):: textinput,option,complex_val
integer :: leng,last,suma,kk,pasos_aux,j,npas,naux
double precision :: pas



 open(unit=input_unit,file=inputfile,status='old',err=100)

 
do while(textinput /= 'end_data')
     read(input_unit,'(A120)') textinput
	 leng=len_trim(textinput) 
	 last=0
	 suma=0
	 option='  '
	 do while(last<leng)
        last=last+1
		if(textinput(last:last)/=' '.and.textinput(last:last)/=':'.and.textinput(last:last)/='='.and.textinput(last:last)/=Z) then
		  if(textinput(last:last)=='#') then
		    last=leng
		  else
		    suma=suma+1
		    option(suma:suma)= textinput(last:last)
		  endif
		endif
	 enddo
	 leng=len_trim(option)
	   last=0	 
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='problem') then
				read(option(last+1:leng),'(a10)') problem
				last=leng+1
			endif
        enddo
       last=0	  
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='tol_solver') then
				read(option(last+1:leng),'(e5.0)') tol_solver
				last=leng+1
			endif
       enddo
       
        last=0	 
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='iter_solver') then
				read(option(last+1:leng),'(i10)') iter_solver
				last=leng+1
			endif
       enddo
        last=0	 
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='freq') then
				read(option(last+1:leng),'(e10.4)') freq
				last=leng+1
			endif
       enddo
        last=0	 
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='phii') then
				read(option(last+1:leng),'(f5.0)') phii
				last=leng+1
			endif
       enddo
        last=0	 
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='pol') then
				read(option(last+1:leng),'(a2)') pol
				last=leng+1
			endif
       enddo
        last=0	 
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='r_scat') then
				read(option(last+1:leng),'(f1.0)') r_scat
				last=leng+1
			endif
       enddo
       
        last=0	 
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='boundary_tol') then
				read(option(last+1:leng),'(e5.0)') boundary_tol
				last=leng+1
			endif
       enddo
       
        !last=0	 
	   !do while(last<leng)
			!last=last+1
			!if(option(1:last)=='deuterium_frac') then
				!read(option(last+1:leng),'(f5.0)') deu_tri_frac
				!last=leng+1
			!endif
       !enddo
       
        last=0
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='mesh_reader') then
				read(option(last+1:leng),'(a4)') reader_type
				last=leng+1
			endif
        enddo
       
		last=0
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='elem_type') then
				read(option(last+1:leng),'(a4)') elem_type
				last=leng+1
            endif
       enddo
       
       		last=0
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='elem_shape') then
				read(option(last+1:leng),'(a4)') elem_shape
				last=leng+1
            endif              
       enddo
       
		if ((elem_type=='line') .and. elem_shape=='tria') then
                nodpel=3
                nodpedge=2 
                Ngauss=3
        else if ((elem_type=='quad') .and. elem_shape=='tria') then
                nodpel=6
                nodpedge=3 
                Ngauss=4
        else if ((elem_type=='line') .and. elem_shape=='squa') then
                nodpel=4
                nodpedge=2 
                Ngauss=4
        else if ((elem_type=='quad') .and. elem_shape=='squa') then
                nodpel=8
                nodpedge=3 
                Ngauss=4 
        endif

		last=0
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='mesh_file') then
				read(option(last+1:leng),'(A)') mesh_file
                mesh_file = trim(adjustl(mesh_file))
				last=leng+1
			endif
        enddo
       
		last=0	 
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='read_logic') then
				read(option(last+1:leng),'(a1)') read_logic
				last=leng+1
			endif
       enddo
              
		last=0
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='system_sym') then
				read(option(last+1:leng),'(a1)') system_sym
				last=leng+1
			endif
       enddo
       
		last=0
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='logic_file') then
				read(option(last+1:leng),'(A)') logic_file
                logic_file = trim(adjustl(logic_file))
				last=leng+1
			endif
       enddo
       
		last=0
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='antenna_source') then
				read(option(last+1:leng),'(a1)') antenna_source
				last=leng+1
			endif
       enddo
       
		last=0
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='plane_wave_source') then
				read(option(last+1:leng),'(a1)') plane_wave_source
				last=leng+1
			endif
       enddo
       
        last=0	 
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='plasma') then
				read(option(last+1:leng),'(i1)') plasma
				last=leng+1
			endif
       enddo
       
            last=0	 
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='ka') then
				read(option(last+1:leng),'(f10.8)') ka
				last=leng+1
			endif
       enddo
       
            last=0	 
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='aa') then
				read(option(last+1:leng),'(f5.3)') aa
				last=leng+1
			endif
       enddo
       
        last=0	 
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='density_flag') then
				read(option(last+1:leng),'(i1)') density_flag
				last=leng+1
			endif
       enddo

	   last=0	 
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='density_e_0') then
				read(option(last+1:leng),'(e10.4)') density_e_0
				last=leng+1
			endif
       enddo
       
        last=0	 
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='deuterium_frac') then
				read(option(last+1:leng),'(f5.0)') deuterium_frac
				last=leng+1
			endif
       enddo

	   last=0	 
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='tritium_frac') then
				read(option(last+1:leng),'(f5.0)') tritium_frac
				last=leng+1
			endif
       enddo

	   last=0	 
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='helium_3_frac') then
				read(option(last+1:leng),'(f5.0)') helium_3_frac
				last=leng+1
			endif
       enddo

        last=0	 
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='magnetic_flag') then
				read(option(last+1:leng),'(i1)') magnetic_flag
				last=leng+1
			endif
       enddo

	   last=0	 
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='mag_field_0') then
				read(option(last+1:leng),'(f10.8)') mag_field_0
				last=leng+1
			endif
       enddo

	   last=0	 
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='major_radius') then
				read(option(last+1:leng),'(f10.8)') major_radius
				last=leng+1
			endif
       enddo
                  
        last=0	 
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='conductivity') then
				read(option(last+1:leng),'(f3.0)') cond
				last=leng+1
			endif
       enddo
           
        last=0	 
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='epsilon_xx_re') then
				read(option(last+1:leng),'(f10.5)') epsilon_scat_xx%re
				last=leng+1
			endif
       enddo
               last=0	 
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='epsilon_yy_re') then
				read(option(last+1:leng),'(f10.5)') epsilon_scat_yy%re
				last=leng+1
			endif
       enddo
               last=0	 
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='epsilon_zz_re') then
				read(option(last+1:leng),'(f10.5)') epsilon_scat_zz%re
				last=leng+1
			endif
       enddo
               last=0	 
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='epsilon_xy_re') then
				read(option(last+1:leng),'(f5.0)') epsilon_scat_xy%re
				last=leng+1
			endif
       enddo
               last=0	 
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='epsilon_yx_re') then
				read(option(last+1:leng),'(f5.0)') epsilon_scat_yx%re
				last=leng+1
			endif
       enddo
       
               last=0	 
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='epsilon_xx_im') then
				read(option(last+1:leng),'(f5.0)') epsilon_scat_xx%im
				last=leng+1
			endif
       enddo
               last=0	 
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='epsilon_yy_im') then
				read(option(last+1:leng),'(f5.0)') epsilon_scat_yy%im
				last=leng+1
			endif
       enddo
               last=0	 
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='epsilon_zz_im') then
				read(option(last+1:leng),'(f5.0)') epsilon_scat_zz%im
				last=leng+1
			endif
       enddo
               last=0	 
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='epsilon_xy_im') then
				read(option(last+1:leng),'(f10.5)') epsilon_scat_xy%im
				last=leng+1
			endif
       enddo
               last=0	 
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='epsilon_yx_im') then
				read(option(last+1:leng),'(f10.5)') epsilon_scat_yx%im
				last=leng+1
			endif
       enddo
       
               last=0	 
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='mu_xx_re') then
				read(option(last+1:leng),'(f5.0)') mu_scat_xx%re
				last=leng+1
			endif
       enddo
               last=0	 
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='mu_yy_re') then
				read(option(last+1:leng),'(f5.0)') mu_scat_yy%re
				last=leng+1
			endif
       enddo
               last=0	 
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='mu_zz_re') then
				read(option(last+1:leng),'(f5.0)') mu_scat_zz%re
				last=leng+1
			endif
       enddo
               last=0	 
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='mu_xy_re') then
				read(option(last+1:leng),'(f5.0)') mu_scat_xy%re
				last=leng+1
			endif
       enddo
               last=0	 
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='mu_yx_re') then
				read(option(last+1:leng),'(f5.0)') mu_scat_yx%re
				last=leng+1
			endif
       enddo
                    last=0	 
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='mu_xx_im') then
				read(option(last+1:leng),'(f5.0)') mu_scat_xx%im
				last=leng+1
			endif
       enddo
               last=0	 
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='mu_yy_im') then
				read(option(last+1:leng),'(f5.0)') mu_scat_yy%im
				last=leng+1
			endif
       enddo
               last=0	 
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='mu_zz_im') then
				read(option(last+1:leng),'(f5.0)') mu_scat_zz%im
				last=leng+1
			endif
       enddo
               last=0	 
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='mu_xy_im') then
				read(option(last+1:leng),'(f5.0)') mu_scat_xy%im
				last=leng+1
			endif
       enddo
               last=0	 
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='mu_yx_im') then
				read(option(last+1:leng),'(f5.0)') mu_scat_yx%im
				last=leng+1
			endif
       enddo
            last=0	 
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='current_source1_re') then
				read(option(last+1:leng),'(f5.0)') current_density1%re
				last=leng+1
			endif
       enddo
            last=0	 
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='current_source1_im') then
				read(option(last+1:leng),'(f5.0)') current_density1%im
				last=leng+1
			endif
       enddo
       
       last=0	 
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='current_source2_re') then
				read(option(last+1:leng),'(f5.0)') current_density2%re
				last=leng+1
			endif
       enddo
            last=0	 
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='current_source2_im') then
				read(option(last+1:leng),'(f5.0)') current_density2%im
				last=leng+1
			endif
       enddo
       
       last=0	 
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='boundary_type') then
				read(option(last+1:leng),'(a3)') boundary_type
				last=leng+1
			endif
       enddo
enddo


close(input_unit)

return

100 write(6,*) 'error in' // inputfile // 'opening'
	print*,"The program will stop"
	stop ' '
    
end subroutine reader

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


