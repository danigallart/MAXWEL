subroutine campo2d_electro(nnodes,nelements,nodpel,solution,material,conect,coor_x,coor_y,sigma1,sigma2,grad_x,grad_y,gradxel_x,gradxel_y)
implicit none
integer :: nnodes,nelements,nodpel,conect(nelements,nodpel),material(nelements)
double precision :: grad_x(nnodes),grad_y(nnodes)
double precision :: gradxel_x(nelements),gradxel_y(nelements)

double precision :: solution(nnodes),coor_x(nnodes),coor_y(nnodes),sigma2 ,sigma1
! local
INTEGER :: NS(nodpel),jel,mat
DOUBLE PRECISION :: X(nodpel),Y(nodpel),sol(nodpel),sigma_el,Ex_el,ey_el,ex(nodpel),ey(nodpel)


DOUBLE PREcIsION:: PHI(nodpel),DPHIX(nodpel),DPHIY(nodpel),AJACO(3,3),AJACOI(3,3),DXHI(nodpel), &
     DTHE(nodpel),XHI,THE,DETER,S11,S22,CNST,denom
DOUBLE PREcIsION,allocatable :: GAUSSPT(:),GAUSSWT(:)

integer kk,jj,i,j,ii,K,NLE,pgaus,I2,ngaus
DOUBLE PRECISION :: t,s,a,c,s1,s2,s3,s4,t1,t2,t3,t4
	 


if(nodpel==4) then
  allocate(gausspt(2),gausswt(2))
  ngaus=2
  GAUSSPT(1)=-0.57735027
  GAUSSPT(2)= 0.57735027
  GAUSSWT(1)= 1.0
  GAUSSWT(2)= 1.0 

elseif(nodpel==9) then
  allocate(gausspt(3),gausswt(3))
  ngaus=3
  GAUSSPT(1)=-0.774596669241483377035853079956
  GAUSSPT(2)= 0.0
  GAUSSPT(3)= 0.774596669241483377035853079956
  GAUSSWT(1)= 0.5555555556
  GAUSSWT(2)= 0.8888888889 
  GAUSSWT(3)=0.5555555556 

elseif(nodpel==16) then
  ngaus=4
  allocate(gausspt(ngaus),gausswt(ngaus))
  GAUSsPT(1)=-0.861136311594053
  GAUsSPT(2)=-0.339981043584856
  GAUsSPT(3)= 0.339981043584856
  GAUsSPT(4)= 0.861136311594053
        
  GAUsSWT(1)=0.347854845137454
  GAUsSWT(2)=0.652145154862546
  GAUsSWT(3)=0.652145154862546
  GAUsSWT(4)=0.347854845137454

endif

grad_x=0.0
grad_y=0.0
DO JEL=1,nelements
   
   DO I=1,nodpel
        ns(I)=conect(JEL,I)
	    j=NS(I)
        X(i)=coor_x(j)
        y(i)=coor_y(j)
        sol(i)=solution(j)
        Ex(i)=0.0
        Ey(i)=0.0
   ENDDO
        
   Ex_el=0
   Ey_el=0

   pgaus=0
   DO KK=1,ngaus
     DO JJ=1,ngaus
         pgaus=pgaus+1

		 t = GAUSSPT(KK)
         s = GAUSSPT(JJ)
      


         call funciones4(t,s,PHI,DPHIX,DPHIY,AJACO,AJACOI,DXHI,DTHE) 
         
         !   JACOBIANO Y SU INVERSA
        AJACO=0.0

          DO K=1,nodpel
              AJACO(1,1)=AJACO(1,1)+DXHI(K)*X(K)
              AJACO(1,2)=AJACO(1,2)+DXHI(K)*Y(K)
              AJACO(2,1)=AJACO(2,1)+DTHE(K)*X(K)
              AJACO(2,2)=AJACO(2,2)+DTHE(K)*Y(K)
          ENDDO

          DETER= AJACO(1,1)*AJACO(2,2)-AJACO(1,2)*AJACO(2,1)
          AJACOI(1,1)=AJACO(2,2)/DETER
          AJACOI(2,2)=AJACO(1,1)/DETER
          AJACOI(1,2)=-AJACO(1,2)/DETER
          AJACOI(2,1)=-AJACO(2,1)/DETER

          DO I=1,nodpel
            DPHIX(I)=AJACOI(1,1)*DXHI(I) + AJACOI(1,2)*DTHE(I) 
            DPHIY(I)=AJACOI(2,1)*DXHI(I) + AJACOI(2,2)*DTHE(I) 
          ENDDO

!          CNST = DETER * GAUSSWT(JJ)*GAUSSWT(KK)*GAUSSWT(II)
         ! gradientes
         DO I=1,nodpel
        
             Ex(i)=Ex(i)+ DPHIX(i)*sol(i)
             Ey(i)=Ey(i)+ DPHIY(i)*sol(i)
             
             Ex_el=Ex_el+ DPHIX(i)*sol(i)
             Ey_el=Ey_el+ DPHIY(i)*sol(i)
         ENDDO


      ENDDO

    ENDDO
      
    do i=1,nodpel
        j=ns(i)
        if(i<5) denom=4.0
        if(i>=5 .and. i<=12) denom=2.0
        if(i>12) denom=1.0

        grad_x(j) =  grad_x(j) - Ex(i)/denom
        grad_y(j) =  grad_y(j) - Ey(i)/denom
    enddo
    
    gradxel_x(jel) = -Ex_el/real(ngaus*ngaus)
    gradxel_y(jel) = -Ey_el/real(ngaus*ngaus)
    
enddo



deallocate(gausspt,gausswt)

end subroutine campo2d_electro

    
    
    subroutine campo2d(nnodes,nelements,nodpel,solution,material,conect,coor_x,coor_y,sigma2,sigma3,sigma4,  &
   &        grad_x,grad_y,gradxel_x,gradxel_y)
implicit none
integer :: nnodes,nelements,nodpel,conect(nelements,nodpel),material(nelements)
double precision :: grad_x(nnodes),grad_y(nnodes)
double precision :: gradxel_x(nelements),gradxel_y(nelements)

double precision :: solution(nnodes),coor_x(nnodes),coor_y(nnodes),sigma2 ,sigma3,sigma4
! local
INTEGER :: NS(nodpel),jel,mat
DOUBLE PRECISION :: X(nodpel),Y(nodpel),sol(nodpel),sigma_el,Ex_el,ey_el,ex(nodpel),ey(nodpel)


DOUBLE PREcIsION:: PHI(nodpel),DPHIX(nodpel),DPHIY(nodpel),AJACO(3,3),AJACOI(3,3),DXHI(nodpel), &
     DTHE(nodpel),XHI,THE,DETER,S11,S22,CNST,denom
DOUBLE PREcIsION,allocatable :: GAUSSPT(:),GAUSSWT(:)

integer kk,jj,i,j,ii,K,NLE,pgaus,I2,ngaus
DOUBLE PRECISION :: t,s,a,c,s1,s2,s3,s4,t1,t2,t3,t4
	 


if(nodpel==8) then
  allocate(gausspt(2),gausswt(2))
  ngaus=2
  GAUSSPT(1)=-0.57735027
  GAUSSPT(2)= 0.57735027
  GAUSSWT(1)= 1.0
  GAUSSWT(2)= 1.0 

elseif(nodpel==27) then
  allocate(gausspt(3),gausswt(3))
  ngaus=3
  GAUSSPT(1)=-0.774596669241483377035853079956
  GAUSSPT(2)= 0.0
  GAUSSPT(3)= 0.774596669241483377035853079956
  GAUSSWT(1)= 0.5555555556
  GAUSSWT(2)= 0.8888888889 
  GAUSSWT(3)=0.5555555556 

elseif(nodpel==16) then
  ngaus=4
  allocate(gausspt(ngaus),gausswt(ngaus))
  GAUSsPT(1)=-0.861136311594053
  GAUsSPT(2)=-0.339981043584856
  GAUsSPT(3)= 0.339981043584856
  GAUsSPT(4)= 0.861136311594053
        
  GAUsSWT(1)=0.347854845137454
  GAUsSWT(2)=0.652145154862546
  GAUsSWT(3)=0.652145154862546
  GAUsSWT(4)=0.347854845137454

endif

grad_x=0.0
grad_y=0.0
DO JEL=1,nelements

   
   
   DO I=1,nodpel
        ns(I)=conect(JEL,I)
	    j=NS(I)
        X(i)=coor_x(j)
        y(i)=coor_y(j)
        sol(i)=solution(j)
        Ex(i)=0.0
        Ey(i)=0.0
   ENDDO
        
   Ex_el=0
   Ey_el=0

   pgaus=0
   DO KK=1,ngaus
     DO JJ=1,ngaus
         pgaus=pgaus+1

		 t = GAUSSPT(KK)
         s = GAUSSPT(JJ)
      

            a =81.0/256.0
            c =1.0/3.0
            s1=1.0+s
            s2=c+s
            s3=c-s
            s4=1.0-s
            t1=1.0+t
            t2=c+t
            t3=c-t
            t4=1.0-t
            phi(1) =   a*s2*s3*s4*t2*t3*t4                   ! 4    10    9    3
            phi( 2) =   a*s1*s2*s3*t2*t3*t4                   ! 
            phi( 3) =   a*s1*s2*s3*t1*t2*t3                   ! 
            phi( 4) =   a*s2*s3*s4*t1*t2*t3                   ! 11   16   15    8
            phi( 5) =-3.0*a*s1*s3*s4*t2*t3*t4              !
            phi( 6) =-3.0*a*s1*s2*s4*t2*t3*t4              !
            phi( 7) =-3.0*a*s1*s2*s3*t1*t3*t4              ! 12   13   14    7
            phi( 8) =-3.0*a*s1*s2*s3*t1*t2*t4              !
            phi( 9) =-3.0*a*s1*s2*s4*t1*t2*t3              !
            phi(10) =-3.0*a*s1*s3*s4*t1*t2*t3              ! 1     5    6    2
            phi(11) =-3.0*a*s2*s3*s4*t1*t2*t4                 
            phi(12) =-3.0*a*s2*s3*s4*t1*t3*t4
            phi(13) = 9.0*a*s1*s3*s4*t1*t3*t4
            phi(14) = 9.0*a*s1*s2*s4*t1*t3*t4
            phi(15) = 9.0*a*s1*s2*s4*t1*t2*t4
            phi(16) = 9.0*a*s1*s3*s4*t1*t2*t4
            DXHI( 1)=  a *t2*t3*t4*(-s2*s3-s2*s4+s3*s4)
            DXHI( 2)=  a *t2*t3*t4*(-s1*s2+s1*s3+s2*s3)
            DXHI( 3)=  a *t1*t2*t3*(-s1*s2+s1*s3+s2*s3)
            DXHI( 4)=  a *t1*t2*t3*(-s2*s3-s2*s4+s3*s4)
            DXHI( 5)=-3.0*a*t2*t3*t4*(-s1*s3-s1*s4+s3*s4)
            DXHI( 6)=-3.0*a*t2*t3*t4*(-s1*s2+s1*s4+s2*s4)
            DXHI( 7)=-3.0*a*t1*t3*t4*(-s1*s2+s1*s3+s2*s3)
            DXHI( 8)=-3.0*a*t1*t2*t4*(-s1*s2+s1*s3+s2*s3)
            DXHI( 9)=-3.0*a*t1*t2*t3*(-s1*s2+s1*s4+s2*s4)
            DXHI(10)=-3.0*a*t1*t2*t3*(-s1*s3-s1*s4+s3*s4)
            DXHI(11)=-3.0*a*t1*t2*t4*(-s2*s3-s2*s4+s3*s4)
            DXHI(12)=-3.0*a*t1*t3*t4*(-s2*s3-s2*s4+s3*s4)
            DXHI(13)= 9.0*a*t1*t3*t4*(-s1*s3-s1*s4+s3*s4)
            DXHI(14)= 9.0*a*t1*t3*t4*(-s1*s2+s1*s4+s2*s4)
            DXHI(15)= 9.0*a*t1*t2*t4*(-s1*s2+s1*s4+s2*s4)
            DXHI(16)= 9.0*a*t1*t2*t4*(-s1*s3-s1*s4+s3*s4)
            DTHE( 1)=  a   *s2*s3*s4*(-t2*t3-t2*t4+t3*t4)
            DTHE( 2)=  a   *s1*s2*s3*(-t2*t3-t2*t4+t3*t4)
            DTHE( 3)=  a   *s1*s2*s3*(-t1*t2+t1*t3+t2*t3)
            DTHE( 4)=  a   *s2*s3*s4*(-t1*t2+t1*t3+t2*t3)
            DTHE( 5)= -3.0*a *s1*s3*s4*(-t2*t3-t2*t4+t3*t4)
            DTHE( 6)= -3.0*a *s1*s2*s4*(-t2*t3-t2*t4+t3*t4)
            DTHE( 7)= -3.0*a *s1*s2*s3*(-t1*t3-t1*t4+t3*t4)
            DTHE( 8)= -3.0*a *s1*s2*s3*(-t1*t2+t1*t4+t2*t4)
            DTHE( 9)= -3.0*a *s1*s2*s4*(-t1*t2+t1*t3+t2*t3)
            DTHE(10)= -3.0*a *s1*s3*s4*(-t1*t2+t1*t3+t2*t3)
            DTHE(11)= -3.0*a *s2*s3*s4*(-t1*t2+t1*t4+t2*t4)
            DTHE(12)= -3.0*a *s2*s3*s4*(-t1*t3-t1*t4+t3*t4)
            DTHE(13)=  9.0*a *s1*s3*s4*(-t1*t3-t1*t4+t3*t4)
            DTHE(14)=  9.0*a *s1*s2*s4*(-t1*t3-t1*t4+t3*t4)
            DTHE(15)=  9.0*a *s1*s2*s4*(-t1*t2+t1*t4+t2*t4)
            DTHE(16)=  9.0*a *s1*s3*s4*(-t1*t2+t1*t4+t2*t4)

!   JACOBIANO Y SU INVERSA
        AJACO=0.0

          DO K=1,nodpel
              AJACO(1,1)=AJACO(1,1)+DXHI(K)*X(K)
              AJACO(1,2)=AJACO(1,2)+DXHI(K)*Y(K)
              AJACO(2,1)=AJACO(2,1)+DTHE(K)*X(K)
              AJACO(2,2)=AJACO(2,2)+DTHE(K)*Y(K)
          ENDDO

          DETER= AJACO(1,1)*AJACO(2,2)-AJACO(1,2)*AJACO(2,1)
          AJACOI(1,1)=AJACO(2,2)/DETER
          AJACOI(2,2)=AJACO(1,1)/DETER
          AJACOI(1,2)=-AJACO(1,2)/DETER
          AJACOI(2,1)=-AJACO(2,1)/DETER

          DO I=1,nodpel
            DPHIX(I)=AJACOI(1,1)*DXHI(I) + AJACOI(1,2)*DTHE(I) 
            DPHIY(I)=AJACOI(2,1)*DXHI(I) + AJACOI(2,2)*DTHE(I) 
          ENDDO

!          CNST = DETER * GAUSSWT(JJ)*GAUSSWT(KK)*GAUSSWT(II)
         ! gradientes
         DO I=1,nodpel
        
             Ex(i)=Ex(i)+ DPHIX(i)*sol(i)
             Ey(i)=Ey(i)+ DPHIY(i)*sol(i)
             
             Ex_el=Ex_el+ DPHIX(i)*sol(i)
             Ey_el=Ey_el+ DPHIY(i)*sol(i)
         ENDDO


      ENDDO

    ENDDO
      
    do i=1,nodpel
        j=ns(i)
        if(i<5) denom=4.0
        if(i>=5 .and. i<=12) denom=2.0
        if(i>12) denom=1.0

        grad_x(j) =  grad_x(j) - Ex(i)/denom
        grad_y(j) =  grad_y(j) - Ey(i)/denom
    enddo
    
    gradxel_x(jel) = -Ex_el/real(ngaus*ngaus)
    gradxel_y(jel) = -Ey_el/real(ngaus*ngaus)
    
enddo


!grad_x=0
!grad_y=0
! do jel=1,nelements
!     do i=1,nodpel
         
!        grad_x(conect(jel,i)) =  grad_x(conect(jel,i)) + gradxel_x(jel)/real(nodpel) 
!        grad_y(conect(jel,i)) =  grad_y(conect(jel,i)) + gradxel_y(jel)/real(nodpel) 
!    enddo
!enddo

deallocate(gausspt,gausswt)

end subroutine campo2d




subroutine campo(nnodes,nelements,nodpe,solution,material,conect,coor_x,coor_y,coor_z,sigma1,sigma2,sigma3,sigma4,  &
   &        grad_x,grad_y,grad_z,gradxel_x,gradxel_y,gradxel_z,jcurrent,time,tempera,problema,jmaxima,sigma_ave,    &
   &        tam_zonex,tam_zoney,elec_sep,potencial,sigma_max,E_maximo,vecino,jcurr_ave,grilla2d,quemado)
use def_bio
use def_constantes
implicit none
character(10):: problema
integer :: nnodes,nelements,nodpe,conect(nelements,nodpe),material(nelements),elec_sep,vecino(nelements),grilla2d(nelements),quemado(nelements)
double precision :: grad_x(nnodes),grad_y(nnodes),grad_z(nnodes),tempera(nnodes)
double precision :: gradxel_x(nelements),gradxel_y(nelements),gradxel_z(nelements),jcurrent(nelements),sigma_ave,potencial,sigma_max,E_maximo

double precision :: solution(nnodes),coor_x(nnodes),coor_y(nnodes),coor_z(nnodes),sigma1,sigma2 ,sigma3,sigma4,   &
      &             campoxl,time,jmaxima,tam_zonex,tam_zoney,jcurr_ave

! local
INTEGER :: NS(nodpe),jel,mat,mate,cuento2d,cuento_quemado
DOUBLE PRECISION :: X(nodpe),Y(nodpe),Z(nodpe),sol(nodpe),sigma_el,Ex_el,ey_el,ez_el,ex(nodpe),ey(nodpe),ez(nodpe)


DOUBLE PREcIsION:: PHI(nodpe),DPHIX(nodpe),DPHIY(nodpe),DPHIZ(nodpe),AJACO(3,3),AJACOI(3,3),DXHI(nodpe), &
     DTHE(nodpe),DPSI(nodpe),XHI,THE,PSI,DETER,S11,S22,S33,CNST,denom,funsigma4,funsigma5,Tmed,funsigma1,xmed,ymed,zmed,funsigma6
DOUBLE PREcIsION,allocatable :: GAUSSPT(:),GAUSSWT(:)


integer kk,jj,i,j,ii,K,NLE,pgaus,I2,ngaus,s_ave


if(nodpe==8) then
  allocate(gausspt(2),gausswt(2))
  ngaus=2
  GAUSSPT(1)=-0.57735027
  GAUSSPT(2)= 0.57735027
  GAUSSWT(1)= 1.0
  GAUSSWT(2)= 1.0 

elseif(nodpe==27) then
  allocate(gausspt(3),gausswt(3))
  ngaus=3
  GAUSSPT(1)=-0.774596669241483377035853079956
  GAUSSPT(2)= 0.0
  GAUSSPT(3)= 0.774596669241483377035853079956
  GAUSSWT(1)= 0.5555555556
  GAUSSWT(2)= 0.8888888889 
  GAUSSWT(3)=0.5555555556 

endif

grad_x=0.0
grad_y=0.0
grad_z=0.0
jmaxima=0.0
sigma_ave=0.0
mate=0
s_ave=0
sigma_max=0.0
E_maximo=0.0
jcurr_ave=.0
cuento2d=0
cuento_quemado=0

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
   endif
   xmed=0.0
   ymed=0.0
   zmed=0.0
  Tmed=0.0
   DO I=1,nodpe
        ns(I)=conect(JEL,I)
	    j=NS(I)
        X(i)=coor_x(j)
        y(i)=coor_y(j)
        z(i)=coor_z(j)
        sol(i)=solution(j)
        Ex(i)=0.0
        Ey(i)=0.0
        Ez(i)=0.0
        Tmed=Tmed + tempera(j)/real(nodpe)
        xmed=xmed + x(i)/real(nodpe)
        ymed=ymed + y(i)/real(nodpe)
        zmed=zmed + z(i)/real(nodpe)

   ENDDO
        
   Ex_el=0
   Ey_el=0
   Ez_el=0

   pgaus=0
   DO KK=1,ngaus
     DO JJ=1,ngaus
       DO II=1,ngaus
         pgaus=pgaus+1

		 XHI = GAUSSPT(KK)
         THE = GAUSSPT(JJ)
         PSI = GAUSSPT(II)
      
      if(ngaus==2) then
        call funciones8(XHI,THE,PSI,PHI,DPHIX,DPHIY,DPHIZ,AJACO,AJACOI,DXHI,DTHE,DPSI ) 
     else         
        if(OrdenAlya==1) then
           call funcORDENALYA(XHI,THE,PSI,PHI,DPHIX,DPHIY,DPHIZ,AJACO,AJACOI,DXHI,DTHE,DPSI )
        else
           call funciones(XHI,THE,PSI,PHI,DPHIX,DPHIY,DPHIZ,AJACO,AJACOI,DXHI,DTHE,DPSI ) !!ojo orden distinto!!
        endif
    endif



!   JACOBIANO Y SU INVERSA
        AJACO=0.0

          DO K=1,nodpe
              AJACO(1,1)=AJACO(1,1)+DXHI(K)*X(K)
              AJACO(1,2)=AJACO(1,2)+DXHI(K)*Y(K)
              AJACO(1,3)=AJACO(1,3)+DXHI(K)*Z(K)
              AJACO(2,1)=AJACO(2,1)+DTHE(K)*X(K)
              AJACO(2,2)=AJACO(2,2)+DTHE(K)*Y(K)
              AJACO(2,3)=AJACO(2,3)+DTHE(K)*Z(K)
              AJACO(3,1)=AJACO(3,1)+DPSI(K)*X(K)
              AJACO(3,2)=AJACO(3,2)+DPSI(K)*Y(K)
              AJACO(3,3)=AJACO(3,3)+DPSI(K)*Z(K)
          ENDDO

          CALL DETERM(AJACO, AJACOI, DETER)

          DO I=1,nodpe
            DPHIX(I)=AJACOI(1,1)*DXHI(I) + AJACOI(1,2)*DTHE(I) + AJACOI(1,3)*DPSI(I)
            DPHIY(I)=AJACOI(2,1)*DXHI(I) + AJACOI(2,2)*DTHE(I) + AJACOI(2,3)*DPSI(I)
            DPHIZ(I)=AJACOI(3,1)*DXHI(I) + AJACOI(3,2)*DTHE(I) + AJACOI(3,3)*DPSI(I)
          ENDDO

          CNST = DETER * GAUSSWT(JJ)*GAUSSWT(KK)*GAUSSWT(II)
         ! gradientes
        
          DO I=1,nodpe
        
             Ex(i)=Ex(i)+ DPHIX(i)*sol(i)
             Ey(i)=Ey(i)+ DPHIY(i)*sol(i)
             Ez(i)=Ez(i)+ DPHIZ(i)*sol(i)

             Ex_el=Ex_el+ DPHIX(i)*sol(i)/real(nodpe)
             Ey_el=Ey_el+ DPHIY(i)*sol(i)/real(nodpe)
             Ez_el=Ez_el+ DPHIZ(i)*sol(i)/real(nodpe)

          ENDDO
      
       ENDDO
     ENDDO
   ENDDO
      
   do i=1,nodpe
        j=ns(i)
        if(i<9) denom=8.0
        if(i>=9 .and. i<=20) denom=4.0
        if(i>=21 .and. i<=26) denom=2.0
        if(i==27) denom=1.0


        grad_x(j) =  grad_x(j) - Ex(i)/denom
        grad_y(j) =  grad_y(j) - Ey(i)/denom
        grad_z(j) =  grad_z(j) - Ez(i)/denom
   enddo
   gradxel_x(jel) = -Ex_el
   gradxel_y(jel) = -Ey_el
   gradxel_z(jel) = -Ez_el

   
        campoxl = dsqrt(gradxel_x(jel)*gradxel_x(jel) +gradxel_y(jel)*gradxel_y(jel) +gradxel_z(jel)*gradxel_z(jel))
        
        if(problema=='PAPA') then

            if(vecino(jel)==1) then
               campoxl=pendiente*campoxl
            endif
            !       sigma_el = sigma_el * funsigma5(mat,campoxl,Tmed,alfa_b,time) !+ sig0 ! para la papa    
           sigma_el = sigma_el * funsigma1(mat,campoxl,Tmed,alfa_b,time,potencial) !+ sig0 ! para la papa    
    !        sigma_el = sigma_el * funsigma6(mat,campoxl,HayHYALU,time)  ! para Nahuel
    
        elseif(problema=='RATONES') then
              sigma_el = sigma_el * funsigma4(mat,campoxl,HayHYALU)   
              jcurrent(jel)  = -sigma_el*gradxel_x(jel)
             
             if(mat/=11.and.mat/=10) then
                 !s_ave=s_ave+1
                 !sigma_ave = sigma_ave + sigma_el
                 
                 if(grilla2d(jel)==1) then
                     s_ave=s_ave+1
                     sigma_ave = sigma_ave + sigma_el
                     cuento2d=cuento2d+1 
                     jcurr_ave=jcurr_ave+jcurrent(jel)
                 endif
              endif  
       
           if(mat==4 .and. abs(campoxl)>=field_limite) then
              quemado(jel)=1
              cuento_quemado=cuento_quemado + 1
           endif         
        elseif(problema=='BERGUES3D') then
        
            if(mat==4) then
               sigma_el = sigma_el * funsigma1(mat,campoxl,Tmed,alfa_b,time,potencial)
            endif
                       
        endif
        
         
       if(problema=='BOCHUM') then

              jcurrent(jel)  = -sigma_el*gradxel_x(jel)
              if(mat/=11.and.mat/=10) then
                 s_ave=s_ave+1
                 sigma_ave = sigma_ave + sigma_el
              endif  
       
       elseif(problema=='CEREBRO' .or. problema=='PAPA') then
           
                if(xmed>=tam_zonex*0.5-poselex .and. xmed<=tam_zonex*0.5+poselex .and. ymed>=tam_zoney*0.5-elec_sep*0.5 .and. ymed<=tam_zoney*0.5+elec_sep*0.5 ) then
                    s_ave=s_ave+1
                    sigma_ave = sigma_ave + sigma_el
                endif              

!           jcurrent(jel)  = sigma_el* gradxel_y(jel)
           jcurrent(jel)  = sigma_el* campoxl

       endif




      if(mat/=11.and.mat/=10) then
            if(sigma_el > sigma_max) then
                sigma_max=sigma_el
            endif

           if(jmaxima<jcurrent(jel)) then
               jmaxima=jcurrent(jel)
           endif

           if(E_maximo<campoxl) then
                E_maximo=campoxl
           endif

      endif
   
   
enddo

sigma_ave=sigma_ave/s_ave
jcurr_ave=jcurr_ave/cuento2d
Area_electroporada = real(cuento_quemado)/real(NE)

grad_x=0
grad_y=0
grad_z=0
 do jel=1,nelements
     do i=1,nodpe
        grad_x(conect(jel,i)) =  grad_x(conect(jel,i)) + gradxel_x(jel)/real(nodpe) 
        grad_y(conect(jel,i)) =  grad_y(conect(jel,i)) + gradxel_y(jel)/real(nodpe) 
        grad_z(conect(jel,i)) =  grad_z(conect(jel,i)) + gradxel_z(jel)/real(nodpe) 
    enddo
enddo

deallocate(gausspt,gausswt)

end subroutine campo



subroutine campo_tetra(nnodes,nelements,nodpe,solution,material,conect,coor_x,coor_y,coor_z,sigma1,sigma2,sigma3,sigma4,  &
   &        grad_x,grad_y,grad_z,gradxel_x,gradxel_y,gradxel_z,jcurrent,time,tempera,problema,jmaxima,sigma_ave,    &
   &        tam_zonex,tam_zoney,elec_sep,potencial,sigma_max,E_maximo,vecino,quemado,cuento_quemado)
use def_bio
use def_constantes
implicit none
character(10):: problema
integer :: nnodes,nelements,nodpe,cuento_quemado,conect(nelements,nodpe),material(nelements),elec_sep,vecino(nelements),quemado(nelements)
double precision :: grad_x(nnodes),grad_y(nnodes),grad_z(nnodes),tempera(nnodes)
double precision :: gradxel_x(nelements),gradxel_y(nelements),gradxel_z(nelements),jcurrent(nelements),sigma_ave,potencial,sigma_max,E_maximo

double precision :: solution(nnodes),coor_x(nnodes),coor_y(nnodes),coor_z(nnodes),sigma1,sigma2 ,sigma3,sigma4,   &
      &             campoxl,time,jmaxima,tam_zonex,tam_zoney

! local
INTEGER :: NS(nodpe),jel,mat,mate
DOUBLE PRECISION :: X(nodpe),Y(nodpe),Z(nodpe),sol(nodpe),sigma_el,Ex_el,ey_el,ez_el,ex(nodpe),ey(nodpe),ez(nodpe)
DOUBLE PREcIsION::  posgp(3,nodpe),weigp(nodpe),sl,tl,zl


DOUBLE PREcIsION:: PHI(nodpe),DPHIX(nodpe),DPHIY(nodpe),DPHIZ(nodpe),AJACO(3,3),AJACOI(3,3),DXHI(nodpe), &
     DTHE(nodpe),DPSI(nodpe),XHI,THE,PSI,DETER,S11,S22,S33,CNST,denom,funsigma4,funsigma5,Tmed,funsigma1,xmed,ymed,zmed,funsigma6
DOUBLE PREcIsION,allocatable :: GAUSSPT(:),GAUSSWT(:)

integer kk,jj,i,j,ii,K,NLE,pgaus,I2,ngaus,s_ave


  ngaus=4
  DXHI=0
DTHE=0
DPSI=0
PHI=0
        posgp(1,1)= 0.0 
        posgp(2,1)= 0.0
        posgp(3,1)= 0.0
        posgp(1,2)= 1.0
        posgp(2,2)= 0.0
        posgp(3,2)= 0.0
        posgp(1,3)= 0.0
        posgp(2,3)= 1.0
        posgp(3,3)= 0.0
        posgp(1,4)= 0.0
        posgp(2,4)= 0.0
        posgp(3,4)= 1.0
        weigp(  1)= 1.0/24.0
        weigp(  2)= 1.0/24.0
        weigp(  3)= 1.0/24.0
        weigp(  4)= 1.0/24.0
  

grad_x=0.0
grad_y=0.0
grad_z=0.0
jmaxima=0.0
sigma_ave=0.0
mate=0
s_ave=0
sigma_max=0.0
E_maximo=0.0
cuento_quemado=0

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
   endif
   xmed=0.0
   ymed=0.0
   zmed=0.0
  Tmed=0.0
   DO I=1,nodpe
        ns(I)=conect(JEL,I)
	    j=NS(I)
        X(i)=coor_x(j)
        y(i)=coor_y(j)
        z(i)=coor_z(j)
        sol(i)=solution(j)
        Ex(i)=0.0
        Ey(i)=0.0
        Ez(i)=0.0
        Tmed=Tmed + tempera(j)/real(nodpe)
        xmed=xmed + x(i)/real(nodpe)
        ymed=ymed + y(i)/real(nodpe)
        zmed=zmed + z(i)/real(nodpe)

   ENDDO
        
   Ex_el=0
   Ey_el=0
   Ez_el=0

 
   DO KK=1,ngaus
 

 	  sl = posgp(1,kk)
      tl = posgp(2,kk)
      zl = posgp(3,kk)

!  CALCULO LA FUNCION PHI, SU DERIVADA EN FUNCCION DE XHI Y DE THE

       ! call funciones8(XHI,THE,PSI,PHI,DPHIX,DPHIY,DPHIZ,AJACO,AJACOI,DXHI,DTHE,DPSI ) !!ojo orden distinto!!
       ! call funcORDENALYA(XHI,THE,PSI,PHI,DPHIX,DPHIY,DPHIZ,AJACO,AJACOI,DXHI,DTHE,DPSI )
      PHI(   1) = 1.0-sl-tl-zl
      PHI(   2) = sl
      PHI(   3) = tl
      PHI(   4) = zl
      DXHI(1) =-1.0
      DTHE(1) =-1.0
      DPSI(1) =-1.0
      DXHI(2) = 1.0
      DTHE(3) = 1.0
      DPSI(4) = 1.0


!   JACOBIANO Y SU INVERSA
        AJACO=0.0

          DO K=1,nodpe
              AJACO(1,1)=AJACO(1,1)+DXHI(K)*X(K)
              AJACO(1,2)=AJACO(1,2)+DXHI(K)*Y(K)
              AJACO(1,3)=AJACO(1,3)+DXHI(K)*Z(K)
              AJACO(2,1)=AJACO(2,1)+DTHE(K)*X(K)
              AJACO(2,2)=AJACO(2,2)+DTHE(K)*Y(K)
              AJACO(2,3)=AJACO(2,3)+DTHE(K)*Z(K)
              AJACO(3,1)=AJACO(3,1)+DPSI(K)*X(K)
              AJACO(3,2)=AJACO(3,2)+DPSI(K)*Y(K)
              AJACO(3,3)=AJACO(3,3)+DPSI(K)*Z(K)
          ENDDO

          CALL DETERM(AJACO, AJACOI, DETER)

          DO I=1,nodpe
            DPHIX(I)=AJACOI(1,1)*DXHI(I) + AJACOI(1,2)*DTHE(I) + AJACOI(1,3)*DPSI(I)
            DPHIY(I)=AJACOI(2,1)*DXHI(I) + AJACOI(2,2)*DTHE(I) + AJACOI(2,3)*DPSI(I)
            DPHIZ(I)=AJACOI(3,1)*DXHI(I) + AJACOI(3,2)*DTHE(I) + AJACOI(3,3)*DPSI(I)
          ENDDO

          CNST = DETER * weigp(kk)
         ! gradientes
        
          DO I=1,nodpe
        
             Ex(i)=Ex(i)+ DPHIX(i)*sol(i)
             Ey(i)=Ey(i)+ DPHIY(i)*sol(i)
             Ez(i)=Ez(i)+ DPHIZ(i)*sol(i)

             Ex_el=Ex_el+ DPHIX(i)*sol(i)/real(nodpe)
             Ey_el=Ey_el+ DPHIY(i)*sol(i)/real(nodpe)
             Ez_el=Ez_el+ DPHIZ(i)*sol(i)/real(nodpe)

          ENDDO
      
   ENDDO
      
   
   gradxel_x(jel) = -Ex_el
   gradxel_y(jel) = -Ey_el
   gradxel_z(jel) = -Ez_el

   
      campoxl = dsqrt(gradxel_x(jel)*gradxel_x(jel) +gradxel_y(jel)*gradxel_y(jel) +gradxel_z(jel)*gradxel_z(jel))

      sigma_el = sigma_el * funsigma1(mat,campoxl,Tmed,alfa_b,time,potencial) !+ sig0 ! para la papa    

      if(mat==3) then
           s_ave=s_ave+1
           sigma_ave = sigma_ave + sigma_el
      endif              

      jcurrent(jel)  = sigma_el* campoxl

      if(mat==3 .and. jcurrent(jel)>0.15) then
          quemado(jel)=1
          cuento_quemado=cuento_quemado + 1
      endif
    
      if(mat==3) then
          
           if(sigma_el > sigma_max) then
               sigma_max=sigma_el
           endif

           if(jmaxima<jcurrent(jel)) then
               jmaxima=jcurrent(jel)
           endif

           if(E_maximo<campoxl) then
               E_maximo=campoxl
           endif

      endif

  
   
enddo

sigma_ave=sigma_ave/s_ave

grad_x=0
grad_y=0
grad_z=0
 do jel=1,nelements
     do i=1,nodpe
        grad_x(conect(jel,i)) =  grad_x(conect(jel,i)) + gradxel_x(jel)/real(nodpe) 
        grad_y(conect(jel,i)) =  grad_y(conect(jel,i)) + gradxel_y(jel)/real(nodpe) 
        grad_z(conect(jel,i)) =  grad_z(conect(jel,i)) + gradxel_z(jel)/real(nodpe) 
    enddo
enddo





end subroutine campo_tetra


