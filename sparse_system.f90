subroutine sparse_logic
use def_io
use def_variables
use def_vectors
! local
integer, allocatable :: consim(:)
integer :: nno,kk,npas
character*(120) :: text_line


allocate( AD(NP),IA(NP+1),ncount(NP),indep_vect(NP))
allocate( ICX(NP+1),consim(NP) )
allocate(ns(nodpel))

open(unit=logic_unit, file=logic_file, status='unknown')


ad=0.0
indep_vect=0.0

if (read_logic == 'N') then

ia=0
ncount=0
icx=0


IA(1)=1

do nno=1,np
  consim=0
  DO KK=1,ne
   npas=0
   
   DO I=1,nodpel
      ns(I)=conn(KK,I)
      if(nno==ns(i)) npas=1
   ENDDO

   if(npas==1) then   

      DO I=1,nodpel
        DO J=1,nodpel
             if(nno==ns(i).and.nno/=ns(j)) consim(ns(j)) = consim(ns(j))+1 
        ENDDO
      ENDDO
    endif
  enddo

  do i=1,NP
     if(consim(i)/=0) then
        ncount(nno)=ncount(nno)+1
     endif
  enddo
  IA(nno+1) = IA(nno) + ncount(nno)
enddo

NONULL = IA(NP+1)

allocate (an(NONULL),ja(NONULL))
An=0.0

do nno=1,np
  consim=0
  DO KK=1,ne
     npas=0
     DO I=1,nodpel
        NS(I)=conn(KK,I)
        if(nno==ns(i)) npas=1
     ENDDO
     if(npas==1) then   
        
        DO I=1,nodpel
          DO J=1,nodpel
               if(nno==ns(i).and.nno/=ns(j)) consim(ns(j)) = consim(ns(j))+1 
          ENDDO
        ENDDO

     endif
  ENDDO
  do i=1,np
      if(consim(i)/=0) then
          JA(IA(nno)+ICX(nno)) =   i
          ICX(nno)=ICX(nno)+1
      endif
   enddo
enddo

write(logic_unit,*) 'CSR FORMAT'
write(logic_unit,*) 'NONULL = ', NONULL
write(logic_unit,*) 'ROWS (IA)'

do i=1,NP+1
    write(logic_unit,*) i, IA(i)
enddo
write(logic_unit,*) 'END_ROWS'


write(logic_unit,*) 'COLUMNS (JA)'

do i=1,NONULL
    write(logic_unit,*) i, JA(i)
enddo
write(logic_unit,*) 'END_COLUMNS'

write(logic_unit,*) 'COUNTER (ICX)'

do i=1,NP+1
    write(logic_unit,*) i, ICX(i)
enddo
write(logic_unit,*) 'END_COUNTER'


else if (read_logic == 'Y') then
        
    do while(text_line /= 'ROWS (IA)') 
        read(logic_unit,'(A120)') text_line
        text_line = trim(adjustl(text_line))
        if (index(text_line, 'NONULL =') /= 0) then
            read(text_line, '(a12,i10)') label,NONULL
        endif
    enddo
    
    allocate (an(NONULL),ja(NONULL))
    An=0.0

    if (text_line == 'ROWS (IA)') then
        read(logic_unit,'(A120)') text_line
        text_line = trim(adjustl(text_line))
    do while(text_line /= 'END_ROWS')
        read(text_line, '(i10,i10)') i,IA(i)
        read(logic_unit,'(A120)') text_line
        text_line = trim(adjustl(text_line))
    enddo
    endif
    
    read(logic_unit,'(A120)') text_line
    text_line = trim(adjustl(text_line))

    
    if (text_line == 'COLUMNS (JA)') then
        read(logic_unit,'(A120)') text_line
        text_line = trim(adjustl(text_line))
    do while(text_line /= 'END_COLUMNS')
        read(text_line, '(i10,i10)') i,JA(i)
        read(logic_unit,'(A120)') text_line
        text_line = trim(adjustl(text_line))
    enddo
    endif
    
    read(logic_unit,'(A120)') text_line
    text_line = trim(adjustl(text_line))
    
    if (text_line == 'COUNTER (ICX)') then
        read(logic_unit,'(A120)') text_line
        text_line = trim(adjustl(text_line))
    do while(text_line /= 'END_COUNTER') 
        read(text_line, '(i10,i10)') i,icx(i)
        read(logic_unit,'(A120)') text_line
        text_line = trim(adjustl(text_line))
    enddo
    endif

    
end if
    
close(logic_unit) 

write(control_unit,*) 'Nonulos del sistema: ',NONULL

deallocate(consim)

end subroutine sparse_logic
