subroutine bio_therma2d()
use def_variables
use def_bio
use def_constantes
use def_solver
implicit none
double precision :: time_tot,deltat,caso,area_quemada
integer :: number

! local 
double precision :: x(nodpel),y(nodpel),esm(nodpel,nodpel),ef(nodpel),adiag,sigma_el,err,qe
double precision :: des(nodpel),masa_mini(nodpel),ef_ant(nodpel)
integer :: ns(nodpel)
integer  NPASO,KK,JEL,I,II,ncase,inode,ipoin,jnode,jj,iaux,keje,jj2,mat,j,iter
double precision :: error,epsil,denom,numer,sol(nodpel),funsig2d,campoxl
integer :: nconta,NCOTA,NCOTA2,nconta2,ncase2,kmax
double precision ::  coef_pul,dtime,a1,a2,time,error2,epsil2,landa,coef1,demed,tmed,joule,tmax,th2,mu

double precision,allocatable :: tempera0(:),rhs_aux(:)

ncase = 2

allocate(tempera0(nnodes),rhs_aux(nnodes))
tempera0=tempera
rhs_ant=0.0

th2=0.8
deltat=1.0

A2 = (1.0-th2)*deltat
A1 = th2*deltat

error=1.0
nconta=0
NCOTA=10000

time=1.0
dtime=0.0

epsil = 1e-4

!do while(dtime < time .and. nconta<NCOTA )
    
   dtime=dtime + deltat
   nconta=nconta+1



    error2=1.0
    epsil2=1e-3
    nconta2=0
    ncota2=10
      
    do while(error2>epsil2 .and. nconta2 < ncota2)
        an=0.0
        ad=0.0
        rhs=0.0
        rhs_aux=0.0

        nconta2= nconta2 + 1

        
        DO JEL=1,nelements
       
           ! aca todos los parametros!!!
            mat=material(jel)        
            if(mat==1 ) then  ! tumor
                sigma_el = sigma1
                k_c = 0.00055                         
                rho_c = 0.0010300    
                cp_c = 3.750 
                QE= 0.000010437 
                coef1 =  W_b*C_b*rho_b
                coef_pul = 1.0 

            elseif(mat==2 .or. mat==3) then
                sigma_el = sigma2
                k_c = 0.000565          
                rho_c = 0.0010390     
                cp_c = 3.680 
                QE=  0.000010437 
                coef1 =  W_b*C_b*rho_b
                coef_pul = 1.0 
       
            endif


            tmed=0.0
            DO I=1,nodpel
                ns(I)=conect(JEL,I)
	            j=NS(I)
                X(i)=coor_x(j)
                y(i)=coor_y(j)

                des(I)=tempera0(j)   
                des(I)=tempera_ant(j)   
            
                ef_ant(i)= rhs_ant(j)
                
               tmed=tmed + tempera(j)/real(nodpel)
            ENDDO

            
            campoxl = (gradxel_x(jel)*gradxel_x(jel) +gradxel_y(jel)*gradxel_y(jel) )
            sigma_el = sigma_el*funsig2d(mat,dsqrt(campoxl),Tmed)

        
            landa =  rho_c*cp_c
            ! tema de condicion de electrodos como disipadores!!!
    
      
            CALL ARMADO16(ncase,JEL,X,Y,ns,nodpel,ESM,EF,sigma_el,qe,sol,gradxel_x(jel),gradxel_y(jel),landa,k_c,mu,coef1,Tblod,coef_pul,campoxl,ef_ant,tmed,des,a1,a2,masa_mini)

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
            	rhs_aux(NS(II)) =  rhs_aux(NS(II)) + EF_ant(ii)
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
        iter=0
        err=1.0
        call CG(nnodes,IA,JA,AN,AD,RHS,tempera,toler,itermax,ITER,ERR)
 
        error2=0.0
        denom=0
        numer=0.0

        Tmax =0.0
        do kk=1,nnodes
            numer=numer + (tempera(kk)-tempera_ant(kk))*(tempera(kk)-tempera_ant(kk))
            denom = denom +  tempera(kk)*tempera(kk)
            tempera_ant(kk)=tempera(kk)

            if(tempera(kk).gt.tmax) then
               tmax=tempera(kk)
               kmax=kk
            endif
        enddo
        error2=dsqrt(numer/denom)
        temper_max=tmax 
      
          WRITE(unit_cont,'(A25,i3,3e15.5,i10)') 'iteraciones termicas No lineales',nconta2,error2,tmax,kmax  
          WRITE(*,'(A25,i3,3e15.5,i10)') 'iteraciones termicas No lineales ',nconta2,error2,tmax,kmax  
    
    
   enddo  ! en while no lineal de T  
   
   rhs_ant=rhs_aux

   error=0.0
   denom=0
   numer=0.0
   area_quemada=0
   do kk=1,nnodes
        numer=numer + (tempera(kk)-tempera0(kk))*(tempera(kk)-tempera0(kk))
        denom = denom +  tempera(kk)*tempera(kk)
        tempera0(kk)=tempera(kk)

        if(tempera(kk)>temper40) then
           area_quemada = area_quemada + 1
        endif
        

   enddo
   area_quemada=area_quemada/real(nnodes)


  error=dsqrt(numer/denom)
    
    
  ! if(mod(nconta,5)==0) then
!      WRITE(unit_cont,*) 'iteraciones termicas ', nconta2, error2,dtime  
!      WRITE(6,*) '****iteraciones temporales ', dtime, nconta, error
  !    WRITE(6,*) 'iteraciones termicas del loop bio ', nconta, error
  ! endif

 !  if(dtime==historiat(1)) then
 !    ncase2=1
 !     call salida_sol_bio(solucion,dtime,ncase2)
  
 !  endif
    

!enddo ! end while de tiempo

number=1
call salida_sol_2d(ncase,number,tempera)

!deallocate(vector_ant)
deallocate(tempera0,rhs_aux)


write(unit_area,*) potencial,modo,exentric,area_quemada


end subroutine bio_therma2d



subroutine bio_therma2d_time(number,time_tot,deltat,caso)
use def_variables
use def_bio
use def_constantes
use def_solver
implicit none
double precision :: time_tot,deltat,caso
integer :: number

! local 
double precision :: x(nodpel),y(nodpel),esm(nodpel,nodpel),ef(nodpel),adiag,sigma_el,err,qe
double precision :: des(nodpel),masa_mini(nodpel),ef_ant(nodpel)
integer :: ns(nodpel)
integer  NPASO,KK,JEL,I,II,ncase,inode,ipoin,jnode,jj,iaux,keje,jj2,mat,j,iter
double precision :: error,epsil,denom,numer,sol(nodpel),funsigma1,campoxl
integer :: nconta,NCOTA,NCOTA2,nconta2,ncase2,kmax
double precision ::  coef_pul,dtime,a1,a2,time,error2,epsil2,landa,coef1,demed,tmed,joule,tmax,th2,dt,mu

double precision,allocatable :: tempera0(:),rhs_aux(:)

ncase = 5

allocate(tempera0(nnodes),rhs_aux(nnodes))
tempera0=tempera
rhs_ant=0.0

th2= 0.5

dt=0.25 !1
if(dt > deltat)  dt=deltat
    

A2 = (1.0-th2)*dt
A1 = th2*dt



error=1.0
nconta=0
NCOTA=10000
dtime=0.0


do while(dtime < deltat )
    
   dtime=dtime + dt
   nconta=nconta+1

   if(caso==1) then
       joule = 1.0
   else
       joule= 0.0
   endif    

   nconta=nconta+1



 an=0.0
 ad=0.0
 rhs=0.0
 rhs_aux=0.0

 DO JEL=1,nelements
       
       ! aca todos los parametros!!!
        mat=material(jel)        
        

        if(mat==4 ) then  ! tumor
            sigma_el = sigma1
            k_c = 0.00055                         
            rho_c = 0.0010300    
            cp_c = 3.750 
            QE= 0.000010437 
            coef1 =  W_b*C_b*rho_b
            
            coef_pul = 17.0 

        elseif(mat==2 .or. mat==3 .or. mat==1) then ! la papa
            
            sigma_el = sigma2
            
            k_c = 0.000562 ! 0.000565          
            rho_c = 0.0011 ! 0.0010390     
            cp_c = 3.78    ! 3.680 
            QE=  0.000002161   !0.000010437 
            coef1 =  W_b*C_b*rho_b
            coef_pul =  1.0 ! 0.0001 ! 1.0 

        elseif(mat==0 .or. mat==5 ) then
            sigma_el = sigma4
            k_c =  0.000024          
            rho_c = 0.001000      
            cp_c = 0.003
            QE= 0.0
            coef1 = 0.0
            coef_pul = 0.0 

        elseif(mat==10 .or. mat==11 ) then
            sigma_el = sig_electrodo*100
            k_c =  0.015      
            rho_c = 0.007900     
            cp_c = 0.500 
            QE= 0.0 
            coef1 = 0.0
            coef_pul = 0.0 

        endif



            tmed=0.0
            DO I=1,nodpel
                ns(I)=conect(JEL,I)
	            j=NS(I)
                X(i)=coor_x(j)
                y(i)=coor_y(j)

                !des(I)=tempera0(j)   
                des(I)=tempera_ant(j)   
            
                ef_ant(i)= rhs_ant(j)
                
                masa_mini(i)=masa(j)

               tmed=tmed + tempera(j)/real(nodpel)
            ENDDO

            
            campoxl = (gradxel_x(jel)*gradxel_x(jel) +gradxel_y(jel)*gradxel_y(jel) )*joule
            
            sigma_el = sigma_el*funsigma1(mat,dsqrt(campoxl),Tmed,alfa_b,time_tot,potencial) 
        
            landa =  rho_c*cp_c
            ! tema de condicion de electrodos como disipadores!!!
    
      
            CALL ARMADO16(ncase,JEL,X,Y,ns,nodpel,ESM,EF,sigma_el,qe,sol,gradxel_x(jel),gradxel_y(jel),landa,k_c,mu,coef1,Tblod,coef_pul,campoxl,ef_ant,tmed,des,a1,a2,masa_mini)

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
        iter=0
        err=1.0
        call CG(nnodes,IA,JA,AN,AD,RHS,tempera,toler,itermax,ITER,ERR)
 
        error2=0.0
        denom=0
        numer=0.0

        Tmax =0.0
        do kk=1,nnodes
            numer=numer + (tempera(kk)-tempera_ant(kk))*(tempera(kk)-tempera_ant(kk))
            denom = denom +  tempera(kk)*tempera(kk)
            tempera_ant(kk)=tempera(kk)

            if(tempera(kk).gt.tmax) then
               tmax=tempera(kk)
               kmax=kk
            endif
        enddo
        error2=dsqrt(numer/denom)
        temper_max=tmax 
      
          WRITE(unit_cont,'(A33,i3,3e15.5,i10)') 'iteraciones termicas No lineales',nconta2,error2,time_tot,tmax
          WRITE(*,'(A33,i3,3e15.5,i10)') 'iteraciones termicas No lineales ',nconta2,error2,time_tot,tmax
    
     
       rhs_ant=rhs_aux
   
   

enddo ! end while de tiempo



call salida_sol_2d(ncase,number,tempera)

deallocate(tempera0,rhs_aux)

end subroutine bio_therma2d_time
