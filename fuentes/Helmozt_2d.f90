subroutine Helmozt_2d()
use def_variables
use def_solver
use def_constantes
implicit none
!local
double precision :: x(nodpel),y(nodpel),error,epsil,denom,numer,xmed,ymed,err
complex*16 :: esm_cplx(nodpel,nodpel),ef_cplx(nodpel),adiag_cplx,sigma_el_cplx,qe_cplx,sol(nodpel),ak_cplx,ij,funF
integer :: ns(nop_el)
integer  NPASO,KK,JEL,I,II,ncase,inode,ipoin,jnode,jj,iaux,keje,jj2,mat,j,iter,kmax
!complex*16 :: error,epsil,denom,numer
integer :: nconta,NCOTA,number
character *1:: aa=','


ij=(0,1) !sqrt(-1.0)

ncase = 1

error=1.0
nconta=0
epsil = 1.e-3
NCOTA=30
cplx_solucion_ant=0.0
do while(error>epsil .and. nconta< NCOTA )
    
    nconta=nconta+1
    cplx_an=0
    cplx_ad=0
    cplx_rhs=0
        
    DO JEL=1,nelements

        mat=material(jel)        
        
        !if(mat==2 ) then
        !    sigma_el = sigma2   ! dentro de los elecrodos
        !elseif(mat==1 ) then
        !    sigma_el = sigma1   ! fuera de los electrodos
        !endif
        ! cplx_ko,cplx_eo,cplx_er,cplx_sigma,cplx_w,cplx_erc(2,2),cplx_murc(2,2)
        
        sigma_el_cplx = 1./cplx_mur
        
        ak_cplx = cplx_ko**2 * (cplx_er- ij* cplx_sigma/(cplx_eo*cplx_w))

        xmed=0.0
        ymed=0.0
        DO I=1,nodpel
            ns(I)=conect(JEL,I)
	        j=NS(I)
            X(i)=coor_x(j)
            y(i)=coor_y(j)
            sol(i)=cplx_solucion_ant(j)
            xmed=xmed + X(i)/nodpel
            ymed=ymed + Y(i)/nodpel
        ENDDO
        qe_cplx=funF(mat,xmed,ymed)  ! si hay carga en el elemento!

           
        !campoxl = dsqrt(cplx_gradxel_x(jel)*cplx_gradxel_x(jel) +cplx_gradxel_y(jel)*cplx_gradxel_y(jel))

        CALL ARMADOhelmozt(ncase,JEL,X,Y,ns,nodpel,nope,ESM_cplx,EF_cplx,sigma_el_cplx,ak_cplx,qe_cplx)
        
        ! INTRODUZCO EN LAS MATRICES LAS CONDICIONES DE CONTORNO
           
        
        do inode=1,nodpel
            ipoin=NS(inode)
            if(  vec_tierra(ipoin)/=-1 ) then
              adiag_cplx=ESM_cplx(inode,inode)
              do jnode=1,nodpel
                 ESM_cplx(inode,jnode)=0.0
                 EF_cplx (jnode)=EF_cplx(jnode)-ESM_cplx(jnode,inode)*tierra
                 ESM_cplx(jnode,inode)=0.0
              end do
              ESM_cplx(inode,inode)= adiag_cplx
              EF_cplx(inode)       = adiag_cplx* tierra
            end if
            if(  vec_poten(ipoin)/=-1 ) then
             
                  adiag_cplx=ESM_cplx(inode,inode)
                  do jnode=1,nodpel
                     ESM_cplx(inode,jnode)=0.0
                     EF_cplx (jnode)=EF_cplx(jnode)-ESM_cplx(jnode,inode)* potencial
                     ESM_cplx(jnode,inode)=0.0
                  end do
                  ESM_cplx(inode,inode)= adiag_cplx
                  EF_cplx(inode)       = adiag_cplx* potencial
           endif            
            
            if(  vec_temper(ipoin)/=-1 ) then
             
                  adiag_cplx=ESM_cplx(inode,inode)
                  do jnode=1,nodpel
                     ESM_cplx(inode,jnode)=0.0
                     EF_cplx (jnode)=EF_cplx(jnode)-ESM_cplx(jnode,inode)* 0.0
                     ESM_cplx(jnode,inode)=0.0
                  end do
                  ESM_cplx(inode,inode)= adiag_cplx
                  EF_cplx(inode)       = adiag_cplx* 0.0
           endif            
            
        end do


! ENSAMBLO
        DO II=1,nodpel
            cplx_RHS(NS(II))=cplx_RHS(NS(II)) + EF_cplx(II)
	        cplx_AD(NS(II)) = cplx_AD(NS(II)) + ESM_cplx(II,II)
            	     
	        DO JJ=1,nodpel
                JJ2=NS(JJ)+1-NS(II)
	            IF(JJ2.GT.0) THEN
		 
	                DO IAUX = 1,CX(NS(II))
         
		                KEJE = IA(NS(II))+IAUX-1  
	          	  	        
       		            IF( JA(KEJE) .EQ. NS(II) + JJ2 -1) THEN           
 		       
		                    cplx_AN( KEJE ) = cplx_AN(KEJE ) +  ESM_cplx(II,JJ)
                
		                ENDIF
          
		           ENDDO

	            ENDIF
	
	        ENDDO
        ENDDO	  
    

    ENDDO

    !do kk=1,Nnodes
    !    Uin(kk)=funinc(coor_x(kk),coor_y(kk))
    !enddo
    
    ! cplx_rhs=A*cplx_Uin
   ! CALL MATXVECSIM_cplx(IA,JA,cplx_AN,cplx_AD,cplx_Uin,cplx_RHS,nnodes)  
    
     
    cplx_solucion=cplx_RHS
    iter=0
    err=1.0

    call LIN_CG_cplx(nnodes,IA,JA,cplx_AN,cplx_AD,cplx_RHS,cplx_solucion,toler,itermax,ITER,ERR)  
       error=0.0
       denom=0
       numer=0.0
       do kk=1,nnodes
          numer=numer + (cplx_solucion(kk)-cplx_solucion_ant(kk))*(cplx_solucion(kk)-cplx_solucion_ant(kk))
          denom = denom +  cplx_solucion(kk)*cplx_solucion(kk)
          
          cplx_solucion(kk)=0.25*cplx_solucion(kk)+0.75*cplx_solucion_ant(kk)

          cplx_solucion_ant(kk)=cplx_solucion(kk)

          

       enddo
       error=dsqrt(numer/denom)

    WRITE(unit_cont,'(A25,i3,2e15.5,i6)') 'iteraciones poisson ', nconta, error,kmax
    WRITE(6,*) 'iteraciones internas del loop poisson ', nconta, error    
    
    !call campo2d_electro(nnodes,nelements,nodpel,solucion,material,conect,coor_x,coor_y,sigma1,sigma2,grad_x,grad_y,gradxel_x,gradxel_y)

enddo ! end while

!call campo2d(nnodes,nelements,nodpel,solucion,material,conect,coor_x,coor_y,sigma1,sigma2,sigma4,grad_x,grad_y,  &
!      &    gradxel_x,gradxel_y)


number=1
call salida_sol_2d_cplx(ncase,number,cplx_solucion)
 

    end subroutine Helmozt_2d	

    
    
    !   CALCULA LAS MATRICES DE CADA ELEMENTO
 SUBROUTINE armadohelmozt(NCASE,NLE,X,Y,ns,nodpel,nope,ESM,EF,sigma_el,ak,qe)
implicit none
INTEGER :: nodpel,nope,NS(nodpel),NCASE
DOUBLE PRECISION :: X(nodpel),Y(nodpel),Z(nodpel)
complex*16:: EF(nodpel),ESM(nodpel,nodpel),sigma_el,qe,ak

integer kk,jj,i,j,ii,K,NLE,pgaus,I2
	 
DOUBLE PREcIsION:: PHI(nodpel),DPHIX(nodpel),DPHIY(nodpel),DPHIZ(nodpel),AJACO(2,2),AJACOI(2,2),DXHI(nodpel), &
     DTHE(nodpel),DPSI(nodpel),GAUSSPT(nope),GAUSSWT(nope),XHI,THE,PSI,DETER,S11,S22,S33,CNST,SUM,demed,factor,efaux


if(nodpel==4) then
   GAUSSPT(1)= -0.57735027
   GAUSSPT(2)=  0.57735027
   GAUSSWT(1)=1.
   GAUSSWT(2)=1.
elseif(nodpel==9) then
   GAUSSPT(1)= -0.774596669241483377035853079956
   GAUSSPT(2)= 0.0
   GAUSSPT(3)= 0.774596669241483377035853079956

   GAUSSWT(1)= 0.5555555556
    GAUSSWT(2)=0.8888888889 
    GAUSSWT(3)=0.5555555556 
endif


ESM=0.0
EF=0.0
  pgaus=0
  DO KK=1,nope
    DO JJ=1,nope
        pgaus=pgaus+1

		XHI = GAUSSPT(kk)
        THE = GAUSSPT(jj)
      
!  CALCULO LA FUNCION PHI, SU DERIVADA EN FUNCCION DE XHI Y DE THE

        if(nodpel==4) then
            call funciones4(XHI,THE,PHI,DPHIX,DPHIY,AJACO,AJACOI,DXHI,DTHE) 
        endif
  
!   JACOBIANO Y SU INVERSA
        AJACO=0.0

          DO K=1,nodpel
              AJACO(1,1)=AJACO(1,1)+DXHI(K)*X(K)
              AJACO(1,2)=AJACO(1,2)+DXHI(K)*Y(K)
              AJACO(2,1)=AJACO(2,1)+DTHE(K)*X(K)
              AJACO(2,2)=AJACO(2,2)+DTHE(K)*Y(K)
          ENDDO

          DO I=1,nodpel
            DPHIX(I)=AJACOI(1,1)*DXHI(I) + AJACOI(1,2)*DTHE(I) 
            DPHIY(I)=AJACOI(2,1)*DXHI(I) + AJACOI(2,2)*DTHE(I) 
          ENDDO
 
!   DERIVADAS DE LAS FUNCIONES DE INTERPOLACION EN FUNCION DE X-Y-Z
            DETER= AJACO(1,1)*AJACO(2,2)-AJACO(1,2)*AJACO(2,1)
               
            IF(DETER.EQ.0.0) THEN
                      write(6,*) 'Error! determinante cero elemento (shapes) ',pgaus
                      STOP ' '
            ENDIF

            AJACOI(1,1)=AJACO(2,2)/DETER
            AJACOI(2,2)=AJACO(1,1)/DETER
            AJACOI(1,2)=-AJACO(1,2)/DETER
            AJACOI(2,1)=-AJACO(2,1)/DETER

           

!       DIFERENCIAL DE INTEGRACION
          CNST = DETER * GAUSSWT(JJ)*GAUSSWT(KK)


!     MATRIZ DE RIGIDEZ LOCAL Y VETOR INDEPENDIENTE
            DO I=1,nodpel
                 DO J=1,nodpel
                   S11= DPHIX(I)*DPHIX(J)
                   S22= DPHIY(I)*DPHIY(J)
            
                  ESM(I,J)=ESM(I,J) - (S11+S22)*CNST*sigma_el
                  
                  ESM(I,J)=ESM(I,J) + ak*PHI(I)*PHI(j)*CNST ! ojo ak!
                ENDDO

                EF(I)=EF(I) + qe*phi(i)*cnst

               ENDDO

      ENDDO
    ENDDO




    end subroutine armadohelmozt
    
    ! si exsite una funcion F diferente a la del libro
    
    complex*16 function funF(mat,x,y)
    double precision :: x,y
    integer :: mat
    
      funF= 1.0
    end function funF
    

