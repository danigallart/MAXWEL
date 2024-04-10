subroutine thorio_temper()
use def_variables
use def_constantes
use def_solver
use def_bio
implicit none
!local
double precision :: x(nodpel),y(nodpel),z(nodpel),esm(nodpel,nodpel),ef(nodpel),adiag,sigma_el,err,qe
integer :: ns(nodpel)
integer  NPASO,KK,JEL,I,II,ncase,inode,ipoin,jnode,jj,iaux,keje,jj2,mat,j,iter,kmax
double precision :: error,epsil,denom,numer,sol(nodpel),funpoten,funconduc,   &
    &               potencia,Tmed,Pmax,jmaxima,relaja,sigma_max2 
integer :: nconta,NCOTA
integer :: caso_test
double precision :: rprop,poten

error=1.0
nconta=0
epsil = 1e-2
NCOTA=20
solucion_ant=potencial
solucion = potencial

!caso_test=0 ! dado
caso_test=1 ! esfera

poten = 600.0 ! W/cm3

if(caso_test==0) then
   rprop=1.0
else
   rprop=1.1927/2.25
endif

do while(error>epsil .and. nconta< NCOTA )
    
    nconta=nconta+1
    an=0
    ad=0
    rhs=0
    sigma_ave=0.0
        
    DO JEL=1,nelements

        mat=material(jel)     
        
        Tmed=0.0
        
        DO I=1,nodpel
            ns(I)=conect(JEL,I)
	        j=NS(I)
            X(i)=coor_x(j)
            y(i)=coor_y(j)
            z(i)=coor_z(j)
            sol(i)=solucion_ant(j)
            Tmed=Tmed + solucion(j)/real(nodpel)
        ENDDO
        
        sigma_el = funconduc(mat,Tmed,rprop,caso_test)
        potencia = funpoten(mat,rprop,poten,caso_test)

        
        gradxel_x(jel)=0.0
        gradxel_y(jel)=0.0
        gradxel_z(jel)=0.0

        CALL ARMADO4(NCASE,JEL,X,Y,Z,ns,nodpel,ESM,EF,sigma_el,potencia,sol,gradxel_x(jel),gradxel_y(jel),gradxel_z(jel))

! INTRODUZCO EN LAS MATRICES LAS CONDICIONES DE CONTORNO

        do inode=1,nodpel
            ipoin=NS(inode)
            if(  vec_tierra(ipoin)/=-1 ) then
              adiag=ESM(inode,inode)
              do jnode=1,nodpel
                 ESM(inode,jnode)=0.0
                 EF (jnode)=EF(jnode)-ESM(jnode,inode)*potencial
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
    
    error=0.0
    denom=0
    numer=0.0
    Pmax=0
    do kk=1,nnodes
          numer=numer + (solucion(kk)-solucion_ant(kk))*(solucion(kk)-solucion_ant(kk))
          denom = denom +  solucion(kk)*solucion(kk)
          !solucion(kk)=relaja*solucion(kk)+(1-relaja)*solucion_ant(kk)

          solucion_ant(kk)=solucion(kk)

          if(solucion(kk).gt.Pmax) then
             Pmax=solucion(kk)
             kmax=kk
          endif

    enddo
    error=dsqrt(numer/denom)
   
    WRITE(unit_cont,'(A25,i3,e15.5)') 'iteraciones poisson ', nconta, error
    WRITE(6,'(A9,i3,2e15.5)') 'it. pois ', nconta, error, Pmax
    WRITE(unit_cont,'(A25,i3,2e15.5,i8)') 'it. int. de poisson ', nconta, error,Pmax,kmax 

enddo ! end while

WRITE(unit_cc,'(A25,i3,2e15.5,i8)') 'it. int. de poisson ', nconta, Pmax, kmax  

call salida_temper(solucion)

end subroutine thorio_temper	


double precision function funconduc(mat,Tmed,rprop,caso)
implicit none
integer :: mat,caso
double precision :: rprop, Tmed
! local
double precision :: AA,BB,Tk,cond_graf

cond_graf=0.4

if(mat==1) then
   funconduc = cond_graf 
elseif(mat==2) then
   AA=-1.34
   BB=0.02
   Tk=Tmed + 273.0
   funconduc =  1.0/(AA + BB*Tk)
   funconduc = (rprop)*funconduc + cond_graf*(1-rprop)
elseif(mat==3) then  ! Helio  Noagua
   Tk=Tmed + 273.0
   !funconduc = (1.7590E-08*Tmed**3 - 1.2041E-05*Tmed**2 + 2.2318E-03*Tmed + 5.5889E-01)/100.0
    funconduc = 2.639E-5 * Tk**0.7085
endif

!if(caso==1 .and. mat==2) then
!    funconduc = (rprop)*funconduc + cond_graf*(1-rprop)
!endif
   
end function funconduc

double precision function funpoten(mat,rprop,poten,caso)   ! W/cm3
implicit none
integer :: mat,caso
double precision :: rprop,poten

if(mat==1) then
   funpoten = 0.0
elseif(mat==2) then
   if(caso==0) then
      funpoten = poten 
   else
      funpoten = poten * rprop
   
   endif    
elseif(mat==3) then
   funpoten = 0.0 

endif
   
   
end function funpoten

