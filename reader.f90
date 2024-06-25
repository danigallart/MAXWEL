subroutine reader()
use def_io
use def_variables
implicit none
character*(1)::Z 
character*(4):: elem_type
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
			if(option(1:last)=='file mesh') then
				read(option(last+1:leng),'(a20)') file_mesh
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
				read(option(last+1:leng),'(i6)') iter_solver
				last=leng+1
			endif
       enddo
        last=0	 
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='freq') then
				read(option(last+1:leng),'(e6.3)') freq
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
       
        last=0	 
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='deuterium_frac') then
				read(option(last+1:leng),'(f5.0)') deu_tri_frac
				last=leng+1
			endif
       enddo
       
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
				read(option(last+1:leng),'(a3)') elem_type
				last=leng+1
            endif
            if (elem_type=='lin ') then
                nodpel=3
            else if (elem_type=='quad ') then
                nodpel=6
            endif
                
        enddo

		last=0
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='mesh_file') then
				read(option(last+1:leng),*) mesh_file
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
			if(option(1:last)=='logic_file') then
				read(option(last+1:leng),*) logic_file
                logic_file = trim(adjustl(logic_file))
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
			if(option(1:last)=='conductivity') then
				read(option(last+1:leng),'(f2.0)') cond
				last=leng+1
			endif
       enddo
           
        last=0	 
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='epsilon_xx') then
				read(option(last+2:leng-1),'(f2.0,f2.0)') epsilon_scat_xx%re,epsilon_scat_xx%im
				last=leng+1
			endif
       enddo
               last=0	 
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='epsilon_yy') then
				read(option(last+2:leng-1),'(f2.0,f2.0)') epsilon_scat_yy%re,epsilon_scat_yy%im
				last=leng+1
			endif
       enddo
               last=0	 
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='epsilon_zz') then
				read(option(last+2:leng-1),'(f2.0,f2.0)') epsilon_scat_zz%re,epsilon_scat_zz%im
				last=leng+1
			endif
       enddo
               last=0	 
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='epsilon_xy') then
				read(option(last+2:leng-1),'(f2.0,f2.0)') epsilon_scat_xy%re,epsilon_scat_xy%im
				last=leng+1
			endif
       enddo
               last=0	 
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='epsilon_yx') then
				read(option(last+2:leng-1),'(f2.0,f2.0)') epsilon_scat_yx%re,epsilon_scat_yx%im
				last=leng+1
			endif
       enddo
       
               last=0	 
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='mu_xx') then
				read(option(last+2:leng-1),'(f2.0,f2.0)') mu_scat_xx%re,mu_scat_xx%im
				last=leng+1
			endif
       enddo
               last=0	 
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='mu_yy') then
				read(option(last+2:leng-1),'(f2.0,f2.0)') mu_scat_yy%re,mu_scat_yy%im
				last=leng+1
			endif
       enddo
               last=0	 
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='mu_zz') then
				read(option(last+2:leng-1),'(f2.0,f2.0)') mu_scat_zz%re,mu_scat_zz%im
				last=leng+1
			endif
       enddo
               last=0	 
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='mu_xy') then
				read(option(last+2:leng-1),'(f2.0,f2.0)') mu_scat_xy%re,mu_scat_xy%im
				last=leng+1
			endif
       enddo
               last=0	 
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='mu_yx') then
				read(option(last+2:leng-1),'(f2.0,f2.0)') mu_scat_yx%re,mu_scat_yx%im
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


