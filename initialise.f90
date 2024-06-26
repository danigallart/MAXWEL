subroutine initialise
use def_io
implicit none

character*30 :: time_date


  call fdate(time_date) 

! Welcoming message
   print*, 'Welcome to MAXWEL code.'
   print*, 'Time & date ',time_date
   print*, 'For any kind of problems, bug report or suggestions'
   print*, 'contact hdomingo@bsc.es'
   
end subroutine initialise
