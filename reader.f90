subroutine reader()
use def_io
use def_variables
implicit none
character*(1)::Z 
character*(120):: textinput,option
integer :: leng,last,suma,kk,pasos_aux,j,npas,naux
double precision :: pas

!write(Z,'(A1)') Z'09'

do while(textinput /= 'end_data')
     read(input_unit,'(A120)') textinput
!	 call upcase(textinput)
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
!	   do while(last<leng)
!			last=last+1
!			if(option(1:last)=='MALLA') then
!				read(option(last+1:leng),'(i3)') nop_malla
!				last=leng+1
!			endif
!        enddo
!	 last=0	 
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
				read(option(last+1:leng),'(f5.0)') freq
				last=leng+1
			endif
       enddo
        last=0	 
	   do while(last<leng)
			last=last+1
			if(option(1:last)=='phi') then
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
    enddo

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


