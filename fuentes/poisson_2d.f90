subroutine poisson_2d_electro()
use def_variables
use def_solver
use def_constantes
implicit none
!local
double precision :: x(nodpel),y(nodpel),esm(nodpel,nodpel),ef(nodpel),adiag,sigma_el,err,qe
integer :: ns(nop_el)
integer  NPASO,KK,JEL,I,II,ncase,inode,ipoin,jnode,jj,iaux,keje,jj2,mat,j,iter,kmax
double precision :: error,epsil,denom,numer,sol(nodpel),funsig2d,campoxl,Tmed,Pmax,mu  
!double precision :: landa,coef1,coef_pul,ef_ant(nodpel),des(nodpel),a1,a2,masa_mini(nodpel)
integer :: nconta,NCOTA,number
character *1:: aa=','
double precision :: xmed,ymed,camp
ncase = 1



error=1.0
nconta=0
epsil = 1.e-3
NCOTA=30
solucion_ant=0.0
do while(error>epsil .and. nconta< NCOTA )
    
    nconta=nconta+1
    an=0
    ad=0
    rhs=0
        
    DO JEL=1,nelements

        mat=material(jel)        
        if(mat==2 ) then
            sigma_el = sigma2   ! dentro de los elecrodos
        elseif(mat==1 ) then
            sigma_el = sigma1   ! fuera de los electrodos
        endif

        DO I=1,nodpel
            ns(I)=conect(JEL,I)
	        j=NS(I)
            X(i)=coor_x(j)
            y(i)=coor_y(j)
            sol(i)=solucion_ant(j)
            qe=0.0  ! si hay carga en el elemento!
        ENDDO

           
        campoxl = dsqrt(gradxel_x(jel)*gradxel_x(jel) +gradxel_y(jel)*gradxel_y(jel))
        

        CALL ARMADOelectro(ncase,JEL,X,Y,ns,nodpel,nope,ESM,EF,sigma_el,qe)
        
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
           endif            
            
            if(  vec_temper(ipoin)/=-1 ) then
             
                  adiag=ESM(inode,inode)
                  do jnode=1,nodpel
                     ESM(inode,jnode)=0.0
                     EF (jnode)=EF(jnode)-ESM(jnode,inode)* 0.0
                     ESM(jnode,inode)=0.0
                  end do
                  ESM(inode,inode)= adiag
                  EF(inode)       = adiag* 0.0
           endif            
            
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
       do kk=1,nnodes
          numer=numer + (solucion(kk)-solucion_ant(kk))*(solucion(kk)-solucion_ant(kk))
          denom = denom +  solucion(kk)*solucion(kk)
          
          solucion(kk)=0.25*solucion(kk)+0.75*solucion_ant(kk)

          solucion_ant(kk)=solucion(kk)

          

       enddo
       error=dsqrt(numer/denom)

    WRITE(unit_cont,'(A25,i3,2e15.5,i6)') 'iteraciones poisson ', nconta, error,Pmax,kmax
    WRITE(6,*) 'iteraciones internas del loop poisson ', nconta, error    
    
    call campo2d_electro(nnodes,nelements,nodpel,solucion,material,conect,coor_x,coor_y,sigma1,sigma2,grad_x,grad_y,gradxel_x,gradxel_y)

enddo ! end while

!call campo2d(nnodes,nelements,nodpel,solucion,material,conect,coor_x,coor_y,sigma1,sigma2,sigma4,grad_x,grad_y,  &
!      &    gradxel_x,gradxel_y)


number=1
call salida_sol_2d(ncase,number,solucion)

 do kk=1,nelements
        xmed=0
        ymed=0
        DO I=1,nodpel
           j=conect(kk,I)
            xmed=xmed + coor_x(j)/real(nodpel)
            ymed=ymed + coor_y(j)/real(nodpel)
        enddo
        camp = dsqrt(gradxel_x(kk)*gradxel_x(kk) + gradxel_y(kk)*gradxel_y(kk))
       
     
       !  write(1111,'(e15.5,a,e15.5,a,e15.5,a,e15.5)')  xmed,aa,ymed,aa,conducta(kk),aa,camp


    enddo


end subroutine poisson_2d_electro	
      

    
    
    subroutine poisson_2d()
use def_variables
use def_constantes
use def_solver
use def_bio
implicit none
double precision :: time,dtime,caso
integer :: number
!local
double precision :: x(nodpel),y(nodpel),esm(nodpel,nodpel),ef(nodpel),adiag,sigma_el,err,qe
integer :: ns(nodpel)
integer  NPASO,KK,JEL,I,II,ncase,inode,ipoin,jnode,jj,iaux,keje,jj2,mat,j,iter,kmax
double precision :: error,epsil,denom,numer,sol(nodpel),funsig2d,campoxl,Tmed,Pmax,mu  
double precision :: landa,coef1,coef_pul,ef_ant(nodpel),des(nodpel),a1,a2,masa_mini(nope)
integer :: nconta,NCOTA
character *1:: aa=','
double precision :: xmed,ymed,camp
ncase = 1

error=1.0
nconta=0
epsil = 1e-1
NCOTA=30
solucion_ant=0.0

do while(error>epsil .and. nconta< NCOTA )
    
    nconta=nconta+1
    an=0
    ad=0
    rhs=0
        
    DO JEL=1,nelements

        mat=material(jel)        
        if(mat==2 ) then
            sigma_el = sigma2
        elseif(mat==1 ) then
            sigma_el = sigma1
        elseif(mat==3 ) then
            sigma_el =0.000001   
        endif

        Tmed=0.0
        qe=0.0
        DO I=1,nodpel
            ns(I)=conect(JEL,I)
	        j=NS(I)
            X(i)=coor_x(j)
            y(i)=coor_y(j)
            sol(i)=solucion_ant(j)
            Tmed=Tmed + tempera(j)/real(nodpel)
        ENDDO

           
        campoxl = dsqrt(gradxel_x(jel)*gradxel_x(jel) +gradxel_y(jel)*gradxel_y(jel))
        
        
        if(mat<3) then   
           sigma_el = sigma_el *funsig2d(mat,campoxl,Tmed)
        endif
        conducta(jel)=sigma_el
        

        CALL ARMADO16(ncase,JEL,X,Y,ns,nodpel,ESM,EF,sigma_el,qe,sol,gradxel_x(jel),gradxel_y(jel),landa,k_c,mu,coef1,Tblod,coef_pul,campoxl,ef_ant,tmed,des,a1,a2,masa_mini)

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
             
            if(vec_poten(ipoin)==1) then 
                  adiag=ESM(inode,inode)
                  do jnode=1,nodpel
                     ESM(inode,jnode)=0.0
                     EF (jnode)=EF(jnode)-ESM(jnode,inode)* potencial
                     ESM(jnode,inode)=0.0
                  end do
                  ESM(inode,inode)= adiag
                  EF(inode)       = adiag* potencial
           elseif(vec_poten(ipoin)==2) then 
                  adiag=ESM(inode,inode)
                  do jnode=1,nodpel
                     ESM(inode,jnode)=0.0
                     EF (jnode)=EF(jnode)-ESM(jnode,inode)* potencial*-1.0
                     ESM(jnode,inode)=0.0
                  end do
                  ESM(inode,inode)= adiag
                  EF(inode)       = adiag* potencial*-1.0
           endif            
            
            
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
       do kk=1,nnodes
          numer=numer + (solucion(kk)-solucion_ant(kk))*(solucion(kk)-solucion_ant(kk))
          denom = denom +  solucion(kk)*solucion(kk)
          
          solucion(kk)=0.25*solucion(kk)+0.75*solucion_ant(kk)

          solucion_ant(kk)=solucion(kk)

          

       enddo
       error=dsqrt(numer/denom)

    WRITE(unit_cont,'(A25,i3,2e15.5,i6)') 'iteraciones poisson ', nconta, error,Pmax,kmax
    WRITE(6,*) 'iteraciones internas del loop poisson ', nconta, error    
    
    call campo2d(nnodes,nelements,nodpel,solucion,material,conect,coor_x,coor_y,sigma1,sigma2,sigma4,grad_x,grad_y,  &
      &    gradxel_x,gradxel_y)

enddo ! end while

!call campo2d(nnodes,nelements,nodpel,solucion,material,conect,coor_x,coor_y,sigma1,sigma2,sigma4,grad_x,grad_y,  &
!      &    gradxel_x,gradxel_y)


number=1
call salida_sol_2d(ncase,number,solucion)

 do kk=1,nelements
        xmed=0
        ymed=0
        DO I=1,nodpel
           j=conect(kk,I)
            xmed=xmed + coor_x(j)/real(nodpel)
            ymed=ymed + coor_y(j)/real(nodpel)
        enddo
        camp = dsqrt(gradxel_x(kk)*gradxel_x(kk) + gradxel_y(kk)*gradxel_y(kk))
       
     
       !  write(1111,'(e15.5,a,e15.5,a,e15.5,a,e15.5)')  xmed,aa,ymed,aa,conducta(kk),aa,camp


    enddo


end subroutine poisson_2d	
      



double precision function funsig2d(mat,E,Tmed)
implicit none
integer :: mat
double precision :: E,Tmed,alfa_b,T0,E0,E1

T0=37.0
alfa_b=0.032
E0 = 20.0
E1 = 30.0

   if(E> E1) then
      funsig2d = 1  + 2.5   + alfa_b*(Tmed-T0)  
   elseif(E>= E0 .and. E<=E1) then
      funsig2d = 1  + 2.5 * ( E-E0)/(E1-E0)    + alfa_b *(Tmed-T0)  
   else
      funsig2d = 1  + alfa_b *(Tmed-T0)   
   endif
   

end function funsig2d




subroutine poisson_2d_time(number,time,dtime,caso)
use def_variables
use def_constantes
use def_solver
use def_bio
implicit none
double precision :: time,dtime,caso
integer :: number
!local
double precision :: x(nodpel),y(nodpel),esm(nodpel,nodpel),ef(nodpel),adiag,sigma_el,err,qe
integer :: ns(nodpel)
integer  NPASO,KK,JEL,I,II,ncase,inode,ipoin,jnode,jj,iaux,keje,jj2,mat,j,iter,kmax
double precision :: error,epsil,denom,numer,sol(nodpel),funsig2d,campoxl,Tmed,Pmax,mu  
double precision :: landa,coef1,coef_pul,ef_ant(nodpel),des(nodpel),a1,a2,masa_mini(nope)
integer :: nconta,NCOTA
character *1:: aa=','
double precision :: xmed,ymed,camp,funsigma1

ncase = 1

error=1.0
nconta=0
epsil = 1e-2
NCOTA=100
!solucion_ant=0.0

do while(error>epsil .and. nconta< NCOTA )
    
    nconta=nconta+1
    an=0
    ad=0
    rhs=0
        
    DO JEL=1,nelements

        mat=material(jel)        
        if(mat==4 ) then
            sigma_el = sigma1
        elseif(mat==2 .or. mat==3) then
            sigma_el = sigma2
        elseif(mat==1 ) then
            sigma_el = sigma2
        elseif(mat==0 .or. mat==5 ) then
            sigma_el = sigma4
        elseif(mat==10) then
            sigma_el = sig_electrodo
        elseif(mat==11 ) then
            sigma_el = sig_electrodo
        endif

        Tmed=0.0
        qe=0.0
        DO I=1,nodpel
            ns(I)=conect(JEL,I)
	        j=NS(I)
            X(i)=coor_x(j)
            y(i)=coor_y(j)
            sol(i)=solucion_ant(j)
            Tmed=Tmed + tempera(j)/real(nodpel)
        ENDDO

           
        campoxl = dsqrt(gradxel_x(jel)*gradxel_x(jel) +gradxel_y(jel)*gradxel_y(jel))

        sigma_el = sigma_el * funsigma1(mat,campoxl,Tmed,alfa_b,time,potencial) 
        
        conducta(jel)=sigma_el

        CALL ARMADO16(ncase,JEL,X,Y,ns,nodpel,ESM,EF,sigma_el,qe,sol,gradxel_x(jel),gradxel_y(jel),landa,k_c,mu,coef1,Tblod,coef_pul,campoxl,ef_ant,tmed,des,a1,a2,masa_mini)

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
             
            if(vec_poten(ipoin)==1) then 
                  adiag=ESM(inode,inode)
                  do jnode=1,nodpel
                     ESM(inode,jnode)=0.0
                     EF (jnode)=EF(jnode)-ESM(jnode,inode)* potencial
                     ESM(jnode,inode)=0.0
                  end do
                  ESM(inode,inode)= adiag
                  EF(inode)       = adiag* potencial
           elseif(vec_poten(ipoin)==2) then 
                  adiag=ESM(inode,inode)
                  do jnode=1,nodpel
                     ESM(inode,jnode)=0.0
                     EF (jnode)=EF(jnode)-ESM(jnode,inode)* tierra ! potencial*-1.0
                     ESM(jnode,inode)=0.0
                  end do
                  ESM(inode,inode)= adiag
                  EF(inode)       = adiag* tierra ! potencial*-1.0
           endif            
            
            
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
       do kk=1,nnodes
          numer=numer + (solucion(kk)-solucion_ant(kk))*(solucion(kk)-solucion_ant(kk))
          denom = denom +  solucion(kk)*solucion(kk)
          
          solucion(kk)=0.25*solucion(kk)+0.75*solucion_ant(kk)

          solucion_ant(kk)=solucion(kk)

          

       enddo
       error=dsqrt(numer/denom)

    WRITE(unit_cont,'(A25,i3,e15.5)') 'iteraciones poisson ', nconta, error
    WRITE(6,*) 'iteraciones internas del loop poisson ', nconta, error    
    
    call campo2d(nnodes,nelements,nodpel,solucion,material,conect,coor_x,coor_y,sigma1,sigma2,sigma4,grad_x,grad_y,  &
      &    gradxel_x,gradxel_y)

enddo ! end while

!call campo2d(nnodes,nelements,nodpel,solucion,material,conect,coor_x,coor_y,sigma1,sigma2,sigma4,grad_x,grad_y,  &
!      &    gradxel_x,gradxel_y)


call salida_sol_2d(ncase,number,solucion)

 
end subroutine poisson_2d_time	
      
