subroutine Ph_evol_2d()
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
double precision ::  coef_pul,dtime,a1,a2,time,error2,epsil2,landa,coef1,demed,tmed,joule,tmax,th2
double precision :: fuentePH,sumideroPH,mu
double precision,allocatable :: Ph0(:),rhs_aux(:)

allocate(Ph_(nnodes),Ph0(nnodes),rhs_aux(nnodes))

ncase = 3
number=1

rhs_ant=0.0
Ph_=6.02e+10  ! at/mmm3  esto equivale a ph=7 inicial!!!  (Mol/dm3)
Ph0= 6.02e+10 ! at/mmm3  esto equivale a ph=7 inicial!!!

fuentePH = 6.0200e+013  ! At/mmm3 equivale a PH=4
sumideroPH= 6.0200e+05 ! 6.0200e+05  ! At/mmm3 equivale a PH=12


time=3600.0

th2=1.0
deltat=60.0 ! 0.01

A2 = (1.0-th2)*deltat
A1 = th2*deltat

error=1.0
nconta=0
NCOTA=10

dtime=0.0

do while(dtime <= time )
    
   dtime=dtime + deltat
   nconta=nconta+1

      
   an=0.0
   ad=0.0
   rhs=0.0
   rhs_aux=0.0

          
   DO JEL=1,nelements
       
           ! aca todos los parametros!!!
            mat=material(jel)        
            if(mat==1 ) then  ! tumor
                k_c =  dprot                         
            elseif(mat==2 .or. mat==3) then
                k_c =  dprot                         
            endif


            DO I=1,nodpel
                ns(I)=conect(JEL,I)
	            j=NS(I)
                X(i)=coor_x(j)
                y(i)=coor_y(j)

                des(I)=Ph0(j)   
            
                ef_ant(i)= rhs_ant(j)
                masa_mini(i) = masa(j)

            ENDDO

         mu = zh*mup        
         landa=1.0
         qe=0.0
         
         CALL ARMADO16(ncase,JEL,X,Y,ns,nodpel,ESM,EF,sigma_el,qe,sol,gradxel_x(jel),gradxel_y(jel),landa,k_c,mu,coef1,Tblod,coef_pul,campoxl,ef_ant,tmed,des,a1,a2,masa_mini)

         do inode=1,nodpel
            ipoin=NS(inode)
            if(  vec_tierra(ipoin)/=-1 ) then
              adiag=ESM(inode,inode)
              do jnode=1,nodpel
                 ESM(inode,jnode)=0.0
                 EF (jnode)=EF(jnode)-ESM(jnode,inode)*sumir*sumideroPH
                 ESM(jnode,inode)=0.0
              end do
              ESM(inode,inode)= adiag
              EF(inode)       = adiag* sumir*sumideroPh
            end if
            if(  vec_poten(ipoin)/=-1 ) then
              adiag=ESM(inode,inode)
              do jnode=1,nodpel
                 ESM(inode,jnode)=0.0
                 EF (jnode)=EF(jnode)-ESM(jnode,inode)* font*fuentePH
                 ESM(jnode,inode)=0.0
              end do
              ESM(inode,inode)= adiag
              EF(inode)       = adiag* font*fuentePh
            end if
        end do

    ! ENSAMBLO
        DO II=1,nodpel
             RHS(NS(II))=RHS(NS(II)) + EF(II)
	         AD(NS(II)) = AD(NS(II)) + ESM(II,II)
             rhs_aux(NS(II))= rhs_aux(NS(II)) + ef_ant(ii)
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


      Ph_=RHS
      iter=0
      err=1.0
      call CG(nnodes,IA,JA,AN,AD,RHS,Ph_,toler/100000,itermax,ITER,ERR)
 
      
   rhs_ant=rhs_aux

   error=0.0
   denom=0
   numer=0.0
   do kk=1,nnodes
        numer=numer + (Ph_(kk)-Ph0(kk))*(Ph_(kk)-Ph0(kk))
        denom = denom +  Ph_(kk)*Ph_(kk)
        Ph0(kk)=Ph_(kk)
   enddo
  error=dsqrt(numer/denom)
    
    
!     WRITE(6,'(A25,i3,2e15.5)') 'Ph: iteraciones de tiempo ',nconta,error,Ph_(5000)  
   
   if(mod(nconta,NCOTA)==0) then
     WRITE(unit_cont,'(A25,i3,2e15.5)') 'Ph: iteraciones de tiempo  ',nconta,error,dtime 
     WRITE(6,'(A25,i3,2e15.5)') 'Ph: iteraciones de tiempo  ',nconta,error,dtime  
     call salida_sol_2d(ncase,number,Ph_)
     number=number+1

   endif

enddo ! end while de tiempo


deallocate(Ph0,Ph_,rhs_aux)

end subroutine Ph_evol_2d
