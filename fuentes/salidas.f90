subroutine salida_sol_2d(ncase,number,solucion)
use def_variables
use def_constantes
implicit none
integer :: number,ncase
double precision :: solucion(nnodes)
! local 
integer :: kk,jj,ii,ngrid,nzgrid2,casex,ncamp,j,i,nconta
double precision :: camp,xmed,ymed,zmed,col,fhy,area_Acida,area_Basica,max_sol,phmed
character *1:: aa=','
character *3::num
character (30) :: arch1,arch2,arch3,arch4,arch5,arch6




if(number<10) then
    write(num,'(i1)') number
elseif(number>=10 .and. number<100) then
    write(num,'(i2)') number
else
    write(num,'(i3)') number
endif


if(ncase==1) then
    arch1=num//file_sal
    arch3=num//file_gra
 
      open(unit=unit_sal,file=arch1)
      open(unit=unit_gra,file=arch3)


!    ngrid=0

    write(unit_sal,*) 'x,    y,    V'
    max_sol=0.0
    do kk=1,nnodes
   
       write(unit_sal,'(e15.5,a,e15.5,a,e15.5,a,e15.5)') coor_x(kk),aa,coor_y(kk),aa,solucion(kk)
      if(solucion(kk).gt.max_sol) then
           max_sol=solucion(kk)
           nconta=kk
       endif
   
    enddo

    ! salida gradientes
    write(unit_gra,*)  ' x,    y,    E(V/cm), sigma(S/cm)'
   
    write(unit_cont,*) 'Vmax:  ',max_sol

    ! salidas por elemeneto
    !electropor =0.0
    !electro_no=0.0
     max_sol=0.0
    do kk=1,nelements
        xmed=0
        ymed=0
        DO I=1,nodpel
           j=conect(kk,I)
            xmed=xmed + coor_x(j)/real(nodpel)
            ymed=ymed + coor_y(j)/real(nodpel)
        enddo

        camp = dsqrt(gradxel_x(kk)*gradxel_x(kk) + gradxel_y(kk)*gradxel_y(kk))
     
!        if(camp>=field_limite) then
!           electropor = electropor + 1.0
!        else
!           electro_no = electro_no + 1.0
!        endif 

!         write(unit_gra,'(e15.5,a,e15.5,a,e15.5,a,e15.5)')  xmed,aa,ymed,aa,camp*10,aa,conducta(kk)*10
         write(unit_gra,'(e15.5,a,e15.5,a,e15.5)')  xmed,aa,ymed,aa,camp

        if(camp.gt.max_sol) then
           max_sol=camp
           nconta=kk
       endif

    enddo

    close(unit_sal)
    close(unit_gra)
    ! write(unit_cc,'(4e15.5)')  electropor,electro_no,electropor/nelements,electro_no/nelements
     write(unit_cont,*) 'Emax:  ',max_sol,nconta

elseif(ncase==2 .or. ncase==5) then

    arch4=num//file_bio
 
      open(unit=unit_bio,file=arch4)



    write(unit_bio,*) 'x,    y,    T'
    max_sol=0.0
    do kk=1,nnodes
   
       write(unit_bio,'(e15.5,a,e15.5,a,e15.5,a,e15.5)') coor_x(kk),aa,coor_y(kk),aa,solucion(kk)
       if(solucion(kk).gt.max_sol) then
           max_sol=solucion(kk)
           nconta=kk
       endif
   
    enddo

    
     close(unit_bio)

     write(unit_cont,*) 'Tmax:  ',max_sol,coor_x(nconta),coor_y(nconta)

    
elseif(ncase==3) then

    arch4=num//file_ph
 
      open(unit=unit_ph,file=arch4)



    write(unit_ph,*) 'x,    y,    Ph'

    area_Acida=0.0
    area_basica=0.0
    do kk=1,nnodes
   
       !write(unit_ph,'(e15.5,a,e15.5,a,e15.5,a,e15.5)') coor_x(kk),aa,coor_y(kk),aa, -log10((solucion(kk)+1e-12)*1000000.0/6.02E23)
       write(unit_ph,'(e15.5,a,e15.5,a,e15.5,a,e15.5)') coor_x(kk),aa,coor_y(kk),aa, solucion(kk)
       if(-log10((solucion(kk)+1e-12)*1000000.0/6.02E23) > limite_Basico) then
           area_Acida=area_Acida + 1
       endif
       if (-log10((solucion(kk)+1e-12)*1000000.0/6.02E23) < limite_acido) then
           area_Basica=area_Basica + 1
       endif
   
    enddo

    
     close(unit_ph)
    
    write(unit_cont,*) 'PH:  ',area_Acida,area_Acida/real(nnodes),area_Basica,area_Basica/real(nnodes)


    area_Acida=0.0
    area_Basica=0.0
    do kk=1,nelements
        xmed=0
        ymed=0
        phmed=0.0
        DO I=1,nodpel
           j=conect(kk,I)
            xmed=xmed + coor_x(j)/real(nodpel)
            ymed=ymed + coor_y(j)/real(nodpel)
            phmed = phmed + solucion(j)/real(nodpel)
        enddo

       if(-log10(phmed*1000000.0/6.02E23) > limite_Basico) then
           area_Acida=area_Acida + 1
       endif
       if (-log10(phmed*1000000.0/6.02E23) < limite_acido) then
           area_Basica=area_Basica + 1
       endif
     
     enddo

     area_Acida = area_Acida/real(nelements)*AreaDOMINIO
     area_Basica = area_Acida/real(nelements)*AreaDOMINIO

     write(unit_area,*)  'Area_pH: ', area_Acida,area_Basica

elseif(ncase==4) then

    arch4=num//file_oh
 
      open(unit=unit_oh,file=arch4)


    write(unit_oh,*) 'x,    y,    Oh'

    area_Acida=0.0
    area_basica=0.0
    do kk=1,nnodes
   
      ! write(unit_oh,'(e15.5,a,e15.5,a,e15.5,a,e15.5)') coor_x(kk),aa,coor_y(kk),aa, 14.0+ log10((solucion(kk)+1e-12)*1000000.0/6.02E23)
       write(unit_oh,'(e15.5,a,e15.5,a,e15.5,a,e15.5)') coor_x(kk),aa,coor_y(kk),aa, solucion(kk)
      
       if(14.0 + log10(solucion(kk)*1000000.0/6.02E23) > limite_Basico) then
           area_Acida=area_Acida + 1
       endif
       if (14.0 + log10(solucion(kk)*1000000.0/6.02E23) < limite_acido) then
           area_Basica=area_Basica + 1
       endif
   
    enddo

    
     close(unit_oh)
    
    write(unit_cont,*) 'OH:  ',area_Acida,area_Acida/real(nnodes),area_Basica,area_Basica/real(nnodes)

    area_Acida=0.0
    area_Basica=0.0
    do kk=1,nelements
        xmed=0
        ymed=0
        phmed=0.0
        DO I=1,nodpel
           j=conect(kk,I)
            xmed=xmed + coor_x(j)/real(nodpel)
            ymed=ymed + coor_y(j)/real(nodpel)
            phmed = phmed + solucion(j)/real(nodpel)
        enddo

       if(14.0+ log10(phmed*1000000.0/6.02E23) > limite_Basico) then
           area_Acida=area_Acida + 1
       endif
       if (14.0+ log10(phmed*1000000.0/6.02E23) < limite_acido) then
           area_Basica=area_Basica + 1
       endif
     
     enddo

     area_Acida = area_Acida/real(nelements)*AreaDOMINIO
     area_Basica = area_Acida/real(nelements)*AreaDOMINIO

    !write(unit_area,*)  'Area_OH: ', area_Acida,area_Basica


endif

end subroutine salida_sol_2d

subroutine salida_sol_2d_cplx(ncase,number,solucion)
use def_variables
use def_constantes
implicit none
integer :: number,ncase
complex*16 :: solucion(nnodes)
! local 
integer :: kk,jj,ii,ngrid,nzgrid2,casex,ncamp,j,i,nconta
double precision :: camp,xmed,ymed,zmed,col,fhy,area_Acida,area_Basica,max_sol,phmed
character *1:: aa=','
character *3::num
character (30) :: arch1,arch2,arch3,arch4,arch5,arch6




if(number<10) then
    write(num,'(i1)') number
elseif(number>=10 .and. number<100) then
    write(num,'(i2)') number
else
    write(num,'(i3)') number
endif

   arch1=num//file_sal
    arch3=num//file_gra
 
      open(unit=unit_sal,file=arch1)
      open(unit=unit_gra,file=arch3)


!    ngrid=0

    write(unit_sal,*) 'x,    y,    V'
    do kk=1,nnodes
   
       write(unit_sal,'(e15.5,a,e15.5,a,e15.5,a,F15.5,a,F15.5)') coor_x(kk),aa,coor_y(kk),aa,real(solucion(kk)),aa,imag(solucion(kk))
      
    enddo

    ! salida gradientes
    write(unit_gra,*)  ' x,    y,    Rsol,    Isol'
   

    ! salidas por elemeneto
    
    do kk=1,nelements
        xmed=0
        ymed=0
        DO I=1,nodpel
           j=conect(kk,I)
            xmed=xmed + coor_x(j)/real(nodpel)
            ymed=ymed + coor_y(j)/real(nodpel)
        enddo

        camp = sqrt(cplx_gradxel_x(kk)*cplx_gradxel_x(kk) + cplx_gradxel_y(kk)*cplx_gradxel_y(kk))
     
!        if(camp>=field_limite) then
!           electropor = electropor + 1.0
!        else
!           electro_no = electro_no + 1.0
!        endif 

!         write(unit_gra,'(e15.5,a,e15.5,a,e15.5,a,e15.5)')  xmed,aa,ymed,aa,camp*10,aa,conducta(kk)*10
         write(unit_gra,'(e15.5,a,e15.5,a,e15.5)')  xmed,aa,ymed,aa,camp


    enddo

    close(unit_sal)
    close(unit_gra)


end subroutine salida_sol_2d_cplx


subroutine salida_temper(solucion)
use def_variables
use def_constantes
use def_bio
implicit none
double precision :: solucion(nnodes)
! local 
integer :: kk,jj,ii,ngrid,nzgrid2,casex,ncamp,j,i,j_ave,s_ave,jel,mat
character *1:: aa=','


 open(unit=unit_sal,file=file_sal)

write(unit_sal,*) 'x,    y,    z,    V'
do kk=1,nnodes
   write(unit_sal,'(e15.5,a,e15.5,a,e15.5,a,e15.5)') coor_x(kk),aa,coor_y(kk),aa,coor_z(kk),aa,solucion(kk)
enddo

close(unit_sal)
end subroutine salida_temper


subroutine salida_sol(number,time,solucion)
use def_variables
use def_constantes
use def_bio
implicit none
double precision :: time,solucion(nnodes)
integer :: number
! local 
integer :: kk,jj,ii,ngrid,nzgrid2,casex,ncamp,j,i,j_ave,s_ave,jel,mat

double precision :: camp,xmed,ymed,zmed,col,fhy,Iave,jmaxima,sigma_el,Tmed,funsigma1,funsigma4

character *1:: aa=','
character *3::num
character (30) :: arch1,arch2,arch3,arch4,arch5,arch6,arch7

! saco la temperatura en movie

if(number<10) then
    write(num,'(i1)') number
elseif(number>=10 .and. number<100) then
    write(num,'(i2)') number
else
    write(num,'(i3)') number
endif

arch1=num//file_sal
arch2=num//file_2D
arch3=num//file_gra
arch4=num//file_grid
arch5=num//file_camp
arch6=num//file_cel
arch7=num//file_II
 
 
  open(unit=unit_sal,file=arch1)
  open(unit=unit_2d,file=arch2)
  open(unit=unit_gra,file=arch3)
  open(unit=unit_grid,file=arch4)
  open(unit=unit_camp,file=arch5)
  open(unit=unit_cel,file=arch6)
  open(unit=unit_II,file=arch7)


ngrid=0

write(unit_sal,*) 'x,    y,    z,    V'
!write(unit_2D,*) 'x,    y,    z,    V'
do kk=1,nnodes
   
   write(unit_sal,'(e15.5,a,e15.5,a,e15.5,a,e15.5)') coor_x(kk),aa,coor_y(kk),aa,coor_z(kk),aa,solucion(kk)
   
   !if(coor_y(kk)>=ncery0 .and. coor_y(kk)<=ncery0+ncery) then
   !  if(coor_x(kk)>=ncerx0-ncerx .and. coor_x(kk)<=ncerx0+ncerx) then
   !    if(coor_z(kk)>=ncerz0-ncerz .and. coor_z(kk)<=ncerz0+ncerz) then
   !      write(unit_2D,'(e15.5,a,e15.5,a,e15.5,a,e15.5)') coor_x(kk),aa,coor_y(kk),aa,coor_z(kk),aa,solucion(kk)  
   !    endif
   !  endif
   !endif 

   
enddo

! sale grid
if(problema=='ELECTRO' ) then
   nzgrid2=nzgrid
   
   write(unit_grid,*) nzgrid,nzgrid2

  do kk=1,nzgrid
     write(unit_grid,*) (coor_y(grid(kk,jj)),jj=1,nzgrid2)
  enddo

  do kk=1,nzgrid
     write(unit_grid,*) (coor_z(grid(kk,jj)),jj=1,nzgrid2)
  enddo

  do kk=1,nzgrid
     write(unit_grid,*) (solucion(grid(kk,jj)),jj=1,nzgrid2)
  enddo

  do kk=1,nzgrid
     write(unit_grid,*) (grad_x(grid(kk,jj)),jj=1,nzgrid2)
  enddo

  do kk=1,nzgrid
     write(unit_grid,*) (grad_y(grid(kk,jj)),jj=1,nzgrid2)
  enddo

  do kk=1,nzgrid
     write(unit_grid,*) (grad_z(grid(kk,jj)),jj=1,nzgrid2)
  enddo

else ! problema=='CEREBRO'

!write(unit_grid,*) Ngrid_z,Ngrid_y,Ngrid_x
!
!do kk=1,Ngrid_z
!   write(unit_grid,*) kk,Zgrid(kk)
!
!   do jj=1,Ngrid_y
!      write(unit_grid,*) (coor_x(grid_z(kk,jj,ii)),ii=1,Ngrid_x)
!   enddo
!   do jj=1,Ngrid_y
!      write(unit_grid,*) (coor_y(grid_z(kk,jj,ii)),ii=1,Ngrid_x)
!   enddo
!
!   do jj=1,Ngrid_y
!      write(unit_grid,*) (grad_x(grid_z(kk,jj,ii)),ii=1,Ngrid_x)
!   enddo
!   do jj=1,Ngrid_y
!      write(unit_grid,*) (grad_y(grid_z(kk,jj,ii)),ii=1,Ngrid_x)
!   enddo
!   do jj=1,Ngrid_y
!      write(unit_grid,*) (grad_z(grid_z(kk,jj,ii)),ii=1,Ngrid_x)
!   enddo
!
!   do jj=1,Ngrid_y
!      write(unit_grid,*) (solucion(grid_z(kk,jj,ii)),ii=1,Ngrid_x)
!   enddo
!
!enddo   
!
endif

! salida campo por elemento
write(unit_gra,*)  ' x,    y,    z,       E (V/cm), Q'

do kk=1,nelements
        xmed=0
        ymed=0
        zmed=0
        DO I=1,nodpel
           j=conect(kk,I)
            xmed=xmed + coor_x(j)/real(nodpel)
            ymed=ymed + coor_y(j)/real(nodpel)
            zmed=zmed + coor_z(j)/real(nodpel)
        enddo

        camp = dsqrt(gradxel_x(kk)*gradxel_x(kk) + gradxel_y(kk)*gradxel_y(kk)+ gradxel_z(kk)*gradxel_z(kk))
     
         write(unit_gra,'(e15.5,a,e15.5,a,e15.5,a,e15.5,a,i6)')  xmed,aa,ymed,aa,zmed,aa,camp*10,aa,quemado(kk)


           write(unit_2D,'(2e15.5)') xmed,camp  

enddo
    
    
! do kk=1,nnodes
!     camp = dsqrt(grad_x(kk)*grad_x(kk) + grad_y(kk)*grad_y(kk)+grad_z(kk)*grad_z(kk))
     
!     if(problema=='CEREBRO') then
        
      !  if(coor_y(kk)>=ncery0-ncery .and. coor_y(kk)<=ncery0+ncery) then
      !    if(coor_x(kk)>=ncerx0-ncerx .and. coor_x(kk)<=ncerx0+ncerx) then
      !      if(coor_z(kk)>=ncerz0-ncerz .and. coor_z(kk)<=ncerz0+ncerz) then
!              write(unit_gra,'(e15.5,a,e15.5,a,e15.5,a,e15.5)')  coor_x(kk),aa,coor_y(kk),aa,coor_z(kk),aa,camp*10
      !      endif
      !    endif
      !  endif

 !       if(coor_y(kk)>=xzona*0.49 .and. coor_y(kk)<=xzona*0.51 .and. coor_z(kk)>=xzona*0.49.and. coor_z(kk)<=xzona*0.51) then

!           write(unit_2D,'(3e15.5)') coor_x(kk),solucion(kk),camp*10  

 !       endif



  !  else 
  !      write(unit_gra,'(e15.5,a,e15.5,a,e15.5,a,e15.5)')  coor_x(kk),aa,coor_y(kk),aa,coor_z(kk),aa,camp
  !  endif
!enddo

write(unit_camp,*)  ' x,    y,    z,        E (V/cm)'
write(unit_II,*)  ' x,    y,    z,        I (A)'

write(unit_cel,*) 'tiempo ',time


electropor =0.0
electro_no=0.0
Iave=0.0
j_ave=0

jmaxima=0.0
sigma_ave=0.0
s_ave=0
sigma_max=0

do kk=1,nelements
   
   camp = dsqrt(gradxel_x(kk)*gradxel_x(kk) + gradxel_y(kk)*gradxel_y(kk)+gradxel_z(kk)*gradxel_z(kk))
   
   

    mat=material(kk)   

   xmed=0
   ymed=0
   zmed=0
   Tmed=0.0
   do jj=1,nodpel
        j=conect(kk,jj)
        xmed = xmed + coor_x(j)/real(nodpel)
        ymed = ymed + coor_y(j)/real(nodpel)
        zmed = zmed + coor_z(j)/real(nodpel)
        Tmed=Tmed + + tempera(j)/real(nodpel)
   enddo


   if(problema/='PAPA_NEW')  then
        if(problema=='PAPA')  then
           if(vecino(kk)==1) then
              camp=pendiente*camp
           endif
       endif 
       if(mat==4) then
           sigma_el = sigma1
       elseif(mat==2 .or. mat==3) then
           sigma_el = sigma2
       elseif(mat==1 ) then
           sigma_el = sigma3
       elseif(mat==0.or.mat==5 ) then
           sigma_el = sigma4
       elseif(mat==10 .or. mat==11 ) then
           sigma_el = sig_electrodo
       endif
   else
       if(mat==1) then
           sigma_el = sigma1
       elseif(mat==2) then
           sigma_el = sigma2
       elseif(mat==3 ) then
           sigma_el = sigma3
      endif
      
   
   endif  
   
   
   if(problema=='PAPA' .or. problema=='PAPA_NEW') then
       sigma_el = sigma_el * funsigma1(mat,camp,Tmed,alfa_b,time,potencial)
   elseif(problema=='RATONES') then
       sigma_el = sigma_el * funsigma4(mat,camp,HayHYALU)   
   endif

   if(HayHYALU==0) then
      fhy=0.0
   else
      fhy=2.0
   endif  

   if(camp>field_limite-fhy) then
      
      electropor = electropor + 1.0
   
   else
      electro_no = electro_no + 1.0
   
   endif

        write(unit_camp,'(e15.5,a,e15.5,a,e15.5,a,e15.5)')   xmed,aa,ymed,aa,zmed,aa,camp*10

!        if(material(kk)==4) then
!          write(unit_cel,'(4e15.5)')  xmed,ymed,zmed,1.0
!        endif

         write(unit_II,'(e15.5,a,e15.5,a,e15.5,a,e15.5)')  xmed,aa,ymed,aa,zmed,aa,Jcurrent(kk)*area_trans ! (paso a cm2)
       if(problema=='BOCHUM') then

            j_ave=j_ave+1;
            Iave = Iave + Jcurrent(kk)*area_trans ! *elec_sep*tam_zonez

       elseif(problema=='PAPA') then
        
         if(xmed>=tam_zonex*0.5-poselex .and. xmed<=tam_zonex*0.5+poselex .and. ymed>=tam_zoney*0.5-elec_sep*0.5 .and. ymed<=tam_zonex*0.5+elec_sep*0.5 ) then
            
            j_ave=j_ave+1;
            Iave = Iave + Jcurrent(kk)*area_trans ! (paso a cm2) *elec_sep*tam_zonez
            
            s_ave=s_ave+1
            sigma_ave = sigma_ave + sigma_el
         endif          
         if(material(kk)/=11.and.material(kk)/=10) then
       
               if(sigma_el > sigma_max) then
                    sigma_max=sigma_el
               endif
               if(jmaxima<jcurrent(kk)) then
                      jmaxima=jcurrent(kk)
               endif
               if(E_maximo<camp) then
                      E_maximo=camp
               endif

          endif
       elseif(problema=='PAPA_NEW') then      
           
            if(mat==3) then
                j_ave=j_ave+1;
                Iave = Iave + Jcurrent(kk)*area_trans 
            
                s_ave=s_ave+1
                sigma_ave = sigma_ave + sigma_el

                if(sigma_el > sigma_max) then
                    sigma_max=sigma_el
               endif
               if(jmaxima<jcurrent(kk)) then
                      jmaxima=jcurrent(kk)
               endif
               if(E_maximo<camp) then
                      E_maximo=camp
               endif


           endif        
            
       endif
    
  
  

enddo

 !write(unit_cc,*) '% electroporado:  ',electropor,electro_no,electropor/nelements,electro_no/nelements
 write(unit_cc,*) 'Iave   ',time,E_maximo,jmaxima*area_trans,sigma_max,temper_max

 !write(6,'(4e15.5)')  electropor,electro_no,electropor/nelements,electro_no/nelements


 close(unit_sal)
 close(unit_2d)
 close(unit_gra)
 close(unit_grid)
 close(unit_camp)
 close(unit_cel)
 close(unit_II)



 write(unit_burn,*) 'tiempo ',time,cuento_quemado
 do kk=1,nelements
    if(material(kk)==3 .and. quemado(kk)==1) then
         write(unit_burn,*) kk
    endif
 enddo
 

end subroutine salida_sol


subroutine salida_sol_bio(number,solucion,time,dtime)
use def_variables
use def_constantes
implicit none
double precision :: solucion(nnodes),dtime,time
integer :: number

! local 
integer :: kk,jj,ii,n,j
character *1:: aa=','
character *4::num
character (30) :: arch1,arch2
double precision :: xmed,ymed,zmed,tmed

! saco la temperatura en movie

if(number<10) then
    write(num,'(i1)') number
elseif(number>=10 .and. number<100) then
    write(num,'(i2)') number
else
    write(num,'(i3)') number
endif


arch1=num//file_bio
arch2=num//file_bio2

 open(unit=unit_bio,file=arch1)
 open(unit=unit_bio2,file=arch2)
    

    write(unit_bio,*) 'x,    y,    z,    T'

   ! write(unit_bio2,*) 'x,    y,    z,    T'
    

do kk=1,nelements
   

   xmed=0
   ymed=0
   zmed=0
   Tmed=0.0
   do jj=1,nodpel
        j=conect(kk,jj)
        xmed = xmed + coor_x(j)/real(nodpel)
        ymed = ymed + coor_y(j)/real(nodpel)
        zmed = zmed + coor_z(j)/real(nodpel)
        Tmed=Tmed + + solucion(j)/real(nodpel)
   enddo


   if(tmed>temper_max .and. material(kk)  <10) then
       temper_max=tmed
   endif

   
   write(unit_bio,'(e15.5,a,e15.5,a,e15.5,a,e15.5)') xmed,aa,ymed,aa,zmed,aa,Tmed
   
enddo


close(unit_bio)
!close(unit_bio2)

! reviso materiales
!do kk=1,nelements
!  if(material(kk) == 11 .or.  material(kk) == 10 ) then
      
!      write(unit_cont,*) 'material ',kk,material(kk)
!      do jj=1,nodpel
!        n= conect(kk,jj)
!        write(unit_cont,*) 'nodo: ',n,tempera(n)
!        write(unit_cont,*)  '            ',coor_x(n),coor_y(n),coor_z(n) 
!      enddo

!  endif

!enddo

write(unit_cont,*) 'Tmax: ',time,temper_max

      



end subroutine salida_sol_bio