subroutine poisson(number,time,dtime,caso)
use def_variables
use def_constantes
use def_solver
use def_bio
implicit none
double precision :: time,dtime,caso
integer :: number
!local
double precision :: x(nodpel),y(nodpel),z(nodpel),esm(nodpel,nodpel),ef(nodpel),adiag,sigma_el,err,qe
integer :: ns(nodpel)
integer  NPASO,KK,JEL,I,II,ncase,inode,ipoin,jnode,jj,iaux,keje,jj2,mat,j,iter,kmax
double precision :: error,epsil,denom,numer,sol(nodpel),funsigma1,funsigma2,funsigma3,funsigma4,funsigma5,funsigma6,   &
    &               campoxl,Tmed,Pmax,jmaxima,relaja,sigma_max2,jcurr_ave 
integer :: nconta,NCOTA


error=1.0
nconta=0
epsil = 0.1
NCOTA=20
solucion_ant=0.0
relaja=0.5

do while(error>epsil .and. nconta< NCOTA )
    
    nconta=nconta+1
    an=0
    ad=0
    rhs=0
    sigma_ave=0.0
        
    DO JEL=1,nelements

        mat=material(jel)     
        if(problema=='PAPA_NEW') then  
            if(mat==1 ) then
                sigma_el = sigma1
            elseif(mat==2 ) then
                sigma_el = sigma2
            elseif(mat==3 ) then
                sigma_el = sigma3
            endif

        else
            if(mat==4 ) then
                sigma_el = sigma1
            elseif(mat==2 .or. mat==3) then
                sigma_el = sigma2
            elseif(mat==1 ) then
                sigma_el = sigma3
            elseif(mat==0 .or. mat==5 ) then
                sigma_el = sigma4
            elseif(mat==10) then
                sigma_el = sig_electrodo
            elseif(mat==11 ) then
                sigma_el = sig_electrodo
            endif
        endif
        Tmed=0.0
        qe=0.0
        DO I=1,nodpel
            ns(I)=conect(JEL,I)
	        j=NS(I)
            X(i)=coor_x(j)
            y(i)=coor_y(j)
            z(i)=coor_z(j)
            sol(i)=solucion_ant(j)
            Tmed=Tmed + tempera(j)/real(nodpel)
        ENDDO

        campoxl = dsqrt(gradxel_x(jel)*gradxel_x(jel) +gradxel_y(jel)*gradxel_y(jel) +gradxel_z(jel)*gradxel_z(jel))

        if(problema=='CEREBRO'.or. problema=='PAPA' .or. problema=='PAPA_NEW' .or. problema=='BERGUES3D') then
           
           
!           sigma_el = sigma_el * funsigma4(mat,campoxl,HayHYALU) ! Bochum Hyalu
         ! if(mat<10 ) then
   !            sigma_el = sigma_el * funsigma5(mat,campoxl,Tmed,alfa_b,time) !+ sig0 ! para la papa    
           if(mat==4) then
               sigma_el = sigma_el * funsigma1(mat,campoxl,Tmed,alfa_b,time,potencial) !+ sig0 ! para la papa    
           endif 
               !sigma_el = sigma_el * funsigma6(mat,campoxl,HayHYALU,time)  ! para Nahuel
         ! endif
           
           
             ! funsigma1(mat,campoxl,Tmed,alfa_b)  ! cerebro?
!            sigma_el = sigma_el *funsigma3(mat,campoxl,HayHYALU)  !*funsigma2(mat,campoxl,Tmed,alfa_b) cerebro?
        elseif(problema=='RATONES') then
             sigma_el = sigma_el * funsigma4(mat,campoxl,HayHYALU)
        endif
        
        gradxel_x(jel)=0.0
        gradxel_y(jel)=0.0
        gradxel_z(jel)=0.0

        if(nodpel==27) then
            CALL ARMADO(NCASE,JEL,X,Y,Z,ns,nodpel,ESM,EF,sigma_el,qe,OrdenAlya)
        elseif(nodpel==8) then
            CALL ARMADO8(NCASE,JEL,X,Y,Z,ns,nodpel,ESM,EF,sigma_el,qe,sol,gradxel_x(jel),gradxel_y(jel),gradxel_z(jel))
        elseif(nodpel==4) then 
           CALL ARMADO4(NCASE,JEL,X,Y,Z,ns,nodpel,ESM,EF,sigma_el,qe,sol,gradxel_x(jel),gradxel_y(jel),gradxel_z(jel))
        endif

! INTRODUZCO EN LAS MATRICES LAS CONDICIONES DE CONTORNO

        do inode=1,nodpel
            ipoin=NS(inode)
            if(  vec_tierra(ipoin)/=-1 ) then
              adiag=ESM(inode,inode)
              do jnode=1,nodpel
                 ESM(inode,jnode)=0.0
                 EF (jnode)=EF(jnode)-ESM(jnode,inode)*tierra
                 ESM(jnode,inode)=0.0
              end do
              ESM(inode,inode)= adiag
              EF(inode)       = adiag* tierra
            end if
            if(  vec_poten(ipoin)/=-1 ) then
              adiag=ESM(inode,inode)
              do jnode=1,nodpel
                 ESM(inode,jnode)=0.0
                 EF (jnode)=EF(jnode)-ESM(jnode,inode)* potencial
                 ESM(jnode,inode)=0.0
              end do
              ESM(inode,inode)= adiag
              EF(inode)       = adiag* potencial
            end if
        end do


! ENSAMBLO
        DO II=1,nodpel
            RHS(NS(II))=RHS(NS(II)) + EF(II)
	        AD(NS(II)) = AD(NS(II)) + ESM(II,II)
            	     
	        DO JJ=1,nodpel
                JJ2=NS(JJ)+1-NS(II)
	            IF(JJ2.GT.0) THEN
		 
	                DO IAUX = 1,CX(NS(II))
         
		                KEJE = IA(NS(II))+IAUX-1  
	          	  	        
       		            IF( JA(KEJE) .EQ. NS(II) + JJ2 -1) THEN           
 		       
		                    AN( KEJE ) = AN(KEJE ) +  ESM(II,JJ)
                
		                ENDIF
          
		           ENDDO

	            ENDIF
	
	        ENDDO
        ENDDO	  
    

    ENDDO

     
    solucion=RHS
    iter=0
    err=1.0
    call CG(nnodes,IA,JA,AN,AD,RHS,solucion,toler,itermax,ITER,ERR)
    !call CGBesta(nnodes,IA,JA,AN,AD,RHS,solucion,toler,itermax,ITER,ERR)
      
   ! if(problema=='CEREBRO' .or. problema=='BOCHUM') then
       error=0.0
       denom=0
       numer=0.0
       Pmax=0
       do kk=1,nnodes
          numer=numer + (solucion(kk)-solucion_ant(kk))*(solucion(kk)-solucion_ant(kk))
          denom = denom +  solucion(kk)*solucion(kk)
          solucion(kk)=relaja*solucion(kk)+(1-relaja)*solucion_ant(kk)

          solucion_ant(kk)=solucion(kk)

          if(solucion(kk).gt.Pmax) then
             Pmax=solucion(kk)
             kmax=kk
          endif

       enddo
       error=dsqrt(numer/denom)
   ! else
   !    error=epsil*0.5
   ! endif

   if(problema=='PAPA_NEW') then
   
      call campo_tetra(nnodes,nelements,nodpel,solucion,material,conect,coor_x,coor_y,coor_z,sigma1,sigma2,sigma3,sigma4,grad_x,grad_y,grad_z,  &
        &    gradxel_x,gradxel_y,gradxel_z,jcurrent,time,tempera,problema,jmaxima,sigma_ave,tam_zonex,tam_zoney,elec_sep,potencial,sigma_max2,  &
        &    E_maximo,vecino,quemado,cuento_quemado)
 
   else
    
      call campo(nnodes,nelements,nodpel,solucion,material,conect,coor_x,coor_y,coor_z,sigma1,sigma2,sigma3,sigma4,grad_x,grad_y,grad_z,  &
        &    gradxel_x,gradxel_y,gradxel_z,jcurrent,time,tempera,problema,jmaxima,sigma_ave,tam_zonex,tam_zoney,elec_sep,potencial,sigma_max,E_maximo,vecino,jcurr_ave,grilla2d,quemado)
 
   endif

    WRITE(unit_cont,'(A25,i3,e15.5)') 'iteraciones poisson ', nconta, error
    WRITE(6,'(A9,i3,5e15.5)') 'it. pois ', nconta, error,jmaxima*area_trans,sigma_max,E_maximo,Pmax    
    WRITE(unit_cont,'(A25,i3,4e15.5)') 'it. int. de poisson ', nconta, error,jmaxima*area_trans,sigma_ave,jcurr_ave    

enddo ! end while

WRITE(unit_cc,'(A25,i3,5e15.5)') 'Maximos', nconta, E_maximo,jmaxima*area_trans,sigma_max,temper_max,jmaxima
WRITE(unit_cc,'(A25,3e15.5)') 'Medios', sigma_ave,jcurr_ave,Area_electroporada   

if(caso==1) then
   call salida_sol(number,time,solucion)
endif


end subroutine poisson	
      
double precision function funsigma5(mat,E,Tmed,alfa_b,t)
implicit none
integer :: mat
double precision :: E,Tmed,alfa_b,t,sig
! local
double precision :: Tcut=25.0, rampa_p,term_temp  


rampa_p=2.5

! termino temporar
term_temp = rampa_p*t
if(term_temp> 20) term_temp= 20.0  

if(mat==4 .or. mat==2 .or. mat==3) then

   funsigma5 = 1  + 22.0*exp(-exp(-0.01*(E*10-250))) + alfa_b*(Tmed-Tcut) + term_temp 
!   funsigma5 = 1  + 11.66*exp(-exp(-0.01*(E*10-250))) + alfa_b*(Tmed-Tcut) + term_temp 

else

   funsigma5 = 1 

endif
   

   
end function funsigma5




double precision function funsigma3(mat,E,Hyalu)
implicit none
integer :: mat,Hyalu
double precision :: E,Tmed,alfa_b

double precision :: Erev,Eirrev,BB,a_min,factorH


if(Hyalu==0) then
   factorH=1.0
else
   factorH=0.6
   
endif

Erev=24.0*factorH
Eirrev=46.0*factorH
a_min=(Erev+Eirrev)*0.5

BB=10000.0


if(mat==4 .or. mat==2 .or. mat==3) then
   if(E> Eirrev) then
      funsigma3 = 3.5*factorH  
   elseif(E>= Erev .and. E<=Eirrev) then
      funsigma3 = 1.0*factorH  + 2.5*factorH/( 1.0 + exp(-(E-a_min)/BB) )  
   else
      funsigma3 = 1.0*factorH   
   endif
else
   funsigma3 = 1.0 
endif


end function funsigma3


double precision function funsigma4(mat,E,Hyalu)
implicit none
integer :: mat,Hyalu
double precision :: E
double precision :: Erev,Eirrev,factorH


if(Hyalu==0) then
   factorH=1.0
else
   factorH=0.6
endif


Erev=24.0*factorH
Eirrev=46.0*factorH

if(mat==4 .or. mat==2 .or. mat==3) then
   if(E> Eirrev) then
      funsigma4 = 3.5  
   elseif(E>= Erev .and. E<=Eirrev) then
         funsigma4 = 1.0  + 2.5 * ( E-Erev)/(Eirrev-Erev)
   else
      funsigma4 = 1.0 
   endif
else
   funsigma4 = 1 
endif


end function funsigma4


double precision function funsigma1(mat,E,Tmed,alfa_b,time,poten)
implicit none
integer :: mat
double precision :: E,Tmed,alfa_b,time,poten
! local
double precision :: rampa_p,term_temp,Tmed0=25.0,DT
!double precision :: factor = 13.868, Emax=120.0,Emin=10.0  ! mm
double precision :: a_t,b_t,c_t,d_t,e_t,f_t;
double precision :: factor = 13.868, Emax=1200.0,Emin=100.0  ! cm


! seteo temperatura 
DT=Tmed-Tmed0
if(DT<0.0) DT=0.0

a_t = 8.183043616717336E+00
b_t = -1.296793778941757E-02
c_t = 2.862554346859271E-02
d_t = 4.872344433745502E-06
e_t = 1.459048339809790E-08
f_t = 8.413328229027457E-05

!d_t = 4.472344433745502E-06
!e_t = 1.20E-08

!if(poten < 1000.0) then
!   term_temp=0.0
!else
!   term_temp =  a_t + b_t*poten + c_t*time + d_t*poten**2 + e_t*time**2 + f_t*poten*time
!endif

if(E < 1000.0) then
   term_temp = 0.0
else
   term_temp = a_t + b_t*E + c_t*time + d_t*E**2 + e_t*time**2 + f_t*E*time
endif

 
!if(mat==4 .or. mat==2 .or. mat==3) then
if(mat==4) then
   if(E>= Emax) then
      funsigma1 = 1  + factor   + alfa_b*DT  +term_temp 
   elseif(E>= Emin .and. E<Emax) then
      funsigma1 = 1  + factor * ( E-Emin)/(Emax-Emin)    + alfa_b*DT  + term_temp
   else
      funsigma1 = 1  + alfa_b*DT  +term_temp 
   endif

else
   funsigma1 = 1.0 
endif
   

end function funsigma1


double precision function funsigma2(mat,E,Tmed,alfa_b)
implicit none
integer :: mat
double precision :: E,Tmed,alfa_b

if(mat==4 .or. mat==2 .or. mat==3) then
   if(E> 35.0) then
      funsigma2 = 1  + 1.0   + alfa_b*(Tmed-37)  
   elseif(E>= 20.0 .and. E<=35.0) then
      funsigma2 = 1  + 1.0 * ( E-20.0)/(35.0-20.0)    + alfa_b*(Tmed-37)  
   else
      funsigma2 = 1  + alfa_b*(Tmed-37)   
   endif
else
   funsigma2 = 1 
endif


end function funsigma2

double precision function funsigma6(mat,E,Hyalu,time)
implicit none
integer :: mat,Hyalu
double precision :: E,time
!local
integer :: Npulso

Npulso= int4(time) + 1

funsigma6 =  1 + E*(time)**0.25 * (Npulso)**0.4

end function funsigma6

