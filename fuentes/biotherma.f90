subroutine biotherma(number,time_tot,deltat,caso)
use def_variables
use def_bio
use def_constantes
use def_solver
implicit none
double precision :: time_tot,deltat,caso
integer :: number

! local 
double precision :: x(nodpel),y(nodpel),z(nodpel),esm(nodpel,nodpel),ef(nodpel),adiag,sigma_el,err,qe
double precision :: des(nodpel),masa_mini(nodpel),ef_ant(nodpel)
integer :: ns(nodpel)
integer  NPASO,KK,JEL,I,II,ncase,inode,ipoin,jnode,jj,iaux,keje,jj2,mat,j,iter
double precision :: error,epsil,denom,numer,sol(nodpel),funsigma1,campoxl,funsigma5,funsigma4
integer :: nconta,NCOTA,nconta2,ncase2,kmax
double precision ::  coef_pul,dtime,a1,a2,time,error2,epsil2,landa,coef1,demed,tmed,joule,tmax

! para integracion temporal en FEM
double precision :: th2 = 1.0,dt ! coeficiente de integracion en el tiempo  

!double precision,allocatable :: rhs_ant(:)
double precision,allocatable :: tempera0(:), rhs_aux(:)


allocate(tempera0(nnodes),rhs_aux(nnodes))
tempera0=tempera
tempera_ant=tempera



!error=1.0
nconta=0
!NCOTA=10000

!epsil = 1e-4


dt=0.1
if(dt > deltat)  dt=deltat

A2 = (1.0-th2)*dt
A1 = th2*dt

dtime=0.0

do while(dtime < deltat )
    
   dtime=dtime + dt
   nconta=nconta+1
!   error2=1.0
!   epsil2=1e-3
!   nconta2=0
!   ncota=20


!   if(dtime.le.historiat(1)) then
   if(caso==1) then
       joule = 1.0
   else
       joule= 0.0
   endif

   
  ! do while(error2>epsil2 .and. nconta2 < ncota)
    an=0.0
    ad=0.0
    rhs=0.0
    rhs_aux=0.0

    nconta2= nconta2 + 1

        
    DO JEL=1,nelements
       
       mat=material(jel)        

       if(problema=='PAPA_NEW') then
        
        if(mat==3) then ! la papa
            
                sigma_el = sigma3
            
                k_c = 0.00562 ! 0.000565          
                rho_c = 1.1 ! 0.0010390     
                cp_c = 3.78    ! 3.680 
                QE=  0.002161   !0.000010437 
                coef1 =  0.0 ! W_b*C_b*rho_b
                coef_pul =  1.0 ! 0.0001 ! 1.0 
         elseif(mat==1) then
                sigma_el = sigma1
                k_c =  0.002  ! W/cmK          
                rho_c = 1.2  ! g/cm3      
                cp_c = 1.4  ! J/g/K 
                QE= 0.0
                coef1 = 0.0
                coef_pul = 0.0 
         elseif(mat==2) then
                sigma_el = sigma2
                k_c =  0.15      
                rho_c = 7.900     
                cp_c = 0.500 
                QE= 0.0 
                coef1 = 0.0
                coef_pul = 0.0 
         endif

       else
           ! aca todos los parametros!!!
         
            if(mat==4 ) then  ! tumor
                sigma_el = sigma1
                k_c = 0.00055                         
                rho_c = 0.0010300    
                cp_c = 3.750 
                QE= 0.000010437 
                coef1 =  W_b*C_b*rho_b
            
                coef_pul = 1.0 

            elseif(mat==2 .or. mat==3) then ! la papa
            
                sigma_el = sigma2
            
                k_c = 0.000562 ! 0.000565          
                rho_c = 0.0011 ! 0.0010390     
                cp_c = 3.78    ! 3.680 
                QE=  0.000002161   !0.000010437 
                coef1 =  W_b*C_b*rho_b
                coef_pul =  1.0 ! 0.0001 ! 1.0 
       
            elseif(mat==1 ) then
                sigma_el = sigma3
                k_c = 0.000565          
                rho_c = 0.0010390     
                cp_c = 3.6800 
                QE= 0.0
                coef1 = 0.0
                coef_pul = 0.0 

            elseif(mat==0 .or. mat==5 ) then
                sigma_el = sigma4
                k_c =  0.000024          
                rho_c = 0.001000      
                cp_c = 0.003
                QE= 0.0
                coef1 = 0.0
                coef_pul = 0.0 

            elseif(mat==10 .or. mat==11 ) then
                sigma_el = sig_electrodo/1000
                k_c =  0.015      
                rho_c = 0.007900     
                cp_c = 0.500 
                QE= 0.0 
                coef1 = 0.0
                coef_pul = 0.0 

            endif
        endif

        demed=0.0
        tmed=0.0
        DO I=1,nodpel
            ns(I)=conect(JEL,I)
	        j=NS(I)
            X(i)=coor_x(j)
            y(i)=coor_y(j)
            z(i)=coor_z(j)

            !des(I)=tempera0(j)   
            des(I)=tempera_ant(j)   
            
            ef_ant(i)= rhs_ant(j)

            masa_mini(i)=masa(j)
            tmed=tmed + tempera(j)/real(nodpel)
            demed=tmed
        ENDDO


        campoxl = (gradxel_x(jel)*gradxel_x(jel) +gradxel_y(jel)*gradxel_y(jel) +gradxel_z(jel)*gradxel_z(jel))*joule 
       
       ! if(campoxl>2000 .and. mat ==3) then
       !     write(6,*) jel,campoxl,gradxel_x(jel),gradxel_y(jel),gradxel_z(jel)
       ! endif
       if(problema=='PAPA_NEW' .or. problema=='PAPA') then
           sigma_el = sigma_el*funsigma1(mat,dsqrt(campoxl),Tmed,alfa_b,time_tot,potencial) 
       else
           sigma_el = sigma_el * funsigma4(mat,campoxl,HayHYALU)  
       endif

        landa =  rho_c*cp_c
        ! tema de condicion de electrodos como disipadores!!!
      
        if(problema=='PAPA_NEW') then  
            call ARMADObio_tetra(NCASE,JEL,X,Y,Z,des,ns,nodpel,nope,ESM,EF,sigma_el,qe,landa,a1,a2,k_c,coef1,Tblod,coef_pul,campoxl,demed,masa_mini,ef_ant)
        else
            CALL ARMADObio(NCASE,JEL,X,Y,Z,des,ns,nodpel,nope,ESM,EF,sigma_el,qe,landa,a1,a2,k_c,coef1,Tblod,coef_pul,campoxl,demed,masa_mini,ef_ant)
        endif
! INTRODUZCO EN LAS MATRICES LAS CONDICIONES DE CONTORNO

        do inode=1,nodpel
            ipoin=NS(inode)
            if(  vec_temper(ipoin)/=-1 ) then
              adiag=ESM(inode,inode)
              do jnode=1,nodpel
                 ESM(inode,jnode)=0.0
                 EF (jnode)=EF(jnode)-ESM(jnode,inode)*Tout
                 ESM(jnode,inode)=0.0
              end do
              ESM(inode,inode)= adiag
              EF(inode)       = adiag* Tout
            end if
        end do

! ENSAMBLO
        DO II=1,nodpel
            RHS(NS(II))=RHS(NS(II)) + EF(II)
	        AD(NS(II)) = AD(NS(II)) + ESM(II,II)
            RHS_aux(NS(II))=RHS_aux(NS(II)) + EF_ant(II)
            	     
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


    tempera=RHS
    rhs_ant=rhs_aux
    iter=0
    err=1.0
    call CG(nnodes,IA,JA,AN,AD,RHS,tempera,toler/1000,itermax,ITER,ERR)
 
    error2=0.0
    denom=0
    numer=0.0

    temper_max =0.0
    do kk=1,nnodes
        numer=numer + (tempera(kk)-tempera_ant(kk))*(tempera(kk)-tempera_ant(kk))
        denom = denom +  tempera(kk)*tempera(kk)

        tempera_ant(kk)=tempera(kk)

!        if(tempera(kk).gt.temper_max) then
!           temper_max=tempera(kk)
!           kmax=kk
!        endif
    enddo
    error2=dsqrt(numer/denom)


    do kk=1,nelements

      if(material(kk)/=10 .and.material(kk)/=11  ) then
           tmed=0.0
           DO I=1,nodpel
                j=conect(kk,I)
                tmed=tmed + tempera(j)
           enddo
           tmed=tmed/real(nodpel)
            if(tmed.gt.temper_max) then
               temper_max=tmed
            endif
      endif
    enddo   

      
    WRITE(unit_cont,'(A25,i3,3e15.5)') 'iteraciones termicas ',nconta2,error2,time_tot,temper_max  
    WRITE(*,'(A5,i3,4e15.5)') 'it. ',nconta,error2,dtime,time_tot,temper_max
    
    
   !enddo  ! en while no lineal de T  
   
    
     tiempo_sumado=tiempo_sumado+dt
     write(unit_temp,*) tiempo_sumado,temper_max

enddo ! end while de temperatura

call salida_sol_bio(number,tempera,time_tot,deltat)

deallocate(tempera0,rhs_aux)

end subroutine biotherma
