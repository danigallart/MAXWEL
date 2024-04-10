!   CALCULA LAS MATRICES DE CADA ELEMENTO
 SUBROUTINE ARMADOelectro(NCASE,NLE,X,Y,ns,nodpel,nope,ESM,EF,sigma_el,qe)
implicit none
INTEGER :: nodpel,nope,NS(nodpel),NCASE
DOUBLE PRECISION :: X(nodpel),Y(nodpel),Z(nodpel),des(nodpel),EF(nodpel),ESM(nodpel,nodpel),sigma_el,qe

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
            
                  ESM(I,J)=ESM(I,J) + (S11+S22)*CNST*sigma_el
                ENDDO

                EF(I)=EF(I) + qe*phi(i)*cnst

               ENDDO

      ENDDO
    ENDDO




end subroutine armadoelectro



    
    
    
    
    
    SUBROUTINE ARMADO16(ncaso,NLE,X,Y,ns,NOPE,ESM,EF,sigma_el,qe,sol,Ex_el,Ey_el,landa,k_c,mu,coef1,T_b,coef_pul,campoxl,ef_ant,tmed,des,a1,a2,masa_mini)
implicit none
INTEGER :: nope,NS(NOPE),ncaso
DOUBLE PRECISION :: X(NOPE),Y(NOPE),EF(NOPE),ESM(NOPE,NOPE),sigma_el,qe,sol(nope),Ex_el,Ey_el,landa,k_c,mu,coef1,T_b,coef_pul,campoxl,ef_ant(nope),tmed,des(nope),a1,a2
DOUBLE PRECISION ::masa_mini(nope)

DOUBLE PREcIsION:: PHI(16),DPHIX(16),DPHIY(16),AJACO(2,2),AJACOI(2,2),DXHI(16), &
     DTHE(16),GAUSSPT(4),GAUSSWT(4),XHI,THE,DETER,S11,S22,CNST,ESK(NOPE,NOPE)

DOUBLE PRECISION :: t,s,a,c,s1,s2,s3,s4,t1,t2,t3,t4,factor,SUM,efaux
integer kk,jj,i,j,ii,K,NLE,pgaus,I2,Ngaus,kgaus


ESM=0.0
EF=0.0
ESK=0.0

        GAUSsPT(1)=-0.861136311594053
        GAUsSPT(2)=-0.339981043584856
        GAUsSPT(3)= 0.339981043584856
        GAUsSPT(4)= 0.861136311594053
        
        GAUsSWT(1)=0.347854845137454
        GAUsSWT(2)=0.652145154862546
        GAUsSWT(3)=0.652145154862546
        GAUsSWT(4)=0.347854845137454
        
        Ngaus = NOPE/4

        kgaus=0 
        do jj=1,Ngaus
          do ii=1,Ngaus
            kgaus=kgaus+1
            t=GAUsSPT(jj)
            s=GAUsSPT(ii)

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


            AJACO=0.0
            AJACOI=0.0
            DO K=1,nope
                    AJACO(1,1)=AJACO(1,1)+dxhi(K)*X(K)
                    AJACO(1,2)=AJACO(1,2)+dxhi(K)*Y(K)
                    AJACO(2,1)=AJACO(2,1)+DTHE(K)*X(K)
                    AJACO(2,2)=AJACO(2,2)+DTHE(K)*Y(K)
            ENDDO

            DETER= AJACO(1,1)*AJACO(2,2)-AJACO(1,2)*AJACO(2,1)
               
            IF(DETER.EQ.0.0) THEN
                      write(6,*) 'Error! determinante cero elemento (shapes) ',kgaus
                      STOP ' '
            ENDIF

            AJACOI(1,1)=AJACO(2,2)/DETER
            AJACOI(2,2)=AJACO(1,1)/DETER
            AJACOI(1,2)=-AJACO(1,2)/DETER
            AJACOI(2,1)=-AJACO(2,1)/DETER

 
!   DERIVADAS DE LAS FUNCIONES DE INTERPOLACION EN FUNCION DE X-Y-Z

            DO i=1,nope
                    dphix(i)=AJACOI(1,1)*dxhi(i) + AJACOI(1,2)*dthe(i) 
                    dphiy(i)=AJACOI(2,1)*dxhi(i) + AJACOI(2,2)*dthe(i) 
            ENDDO

            CNST= deter * GAUsSWT(ii)* GAUsSWT(jj)

           if(ncaso==1) then ! voltaje
                DO I=1,nope
                 DO J=1,nope
                   S11= DPHIX(I)*DPHIX(J)
                   S22= DPHIY(I)*DPHIY(J)
            
                  ESM(I,J)=ESM(I,J) + (S11+S22)*CNST*sigma_el
                ENDDO
               ! Ex_el=Ex_el- DPHIX(i)*sol(i)/real(nope)
               ! Ey_el=Ey_el- DPHIY(i)*sol(i)/real(nope)

                EF(I)=0.0 ! EF(I) + qe*phi(i)*cnst

               ENDDO

           elseif(ncaso==2) then ! biotherma

           !     MATRIZ DE RIGIDEZ LOCAL Y VETOR INDEPENDIENTE

              DO I=1,nope
                DO J=1,nope
                  S11= DPHIX(I)*DPHIX(J)
                  S22= DPHIY(I)*DPHIY(J)

                    ESM(I,J)=ESM(I,J) + (S11+S22)*CNST*k_c 
              
                    ESM(I,J)=ESM(I,J) + coef1 * PHI(i)*PHI(j) *cnst  
                   
                   !ESK(i,j)= ESK(i,j) + PHI(i)*PHI(j) *cnst*landa

               
                ENDDO
             
                 factor= qe*phi(i)*cnst + coef1*T_b*phi(i)*cnst + coef_pul*campoxl*sigma_el*phi(i)*cnst
            
                 EF(I)=EF(I) +factor

              ENDDO
           elseif(ncaso>=3 .and. ncaso <5) then ! Ph

           !     MATRIZ DE RIGIDEZ LOCAL Y VETOR INDEPENDIENTE

              DO I=1,nope
                DO J=1,nope
                  S11= DPHIX(I)*DPHIX(J)
                  S22= DPHIY(I)*DPHIY(J)

                    ESM(I,J)=ESM(I,J) + (S11+S22)*CNST*k_c 
              
                   
                !   ESK(i,j)= ESK(i,j) + PHI(i)*PHI(j) *cnst

                    ESM(I,J)=ESM(I,J) + mu*CNST* phi(j) *( Ex_el*dphix(I) + Ey_el*dphiy(I) ) 

                ENDDO
                
                  ESK(i,i)= ESK(i,i) + masa_mini(i)*LANDA*cnst
              
            
                 EF(I)=EF(I)+qe*PHI(i)*cnst


              ENDDO
           
           
          elseif(ncaso==5) then ! biotherma


             
              DO I=1,nope
                DO J=1,nope
                  S11= DPHIX(I)*DPHIX(J)
                  S22= DPHIY(I)*DPHIY(J)

                    ESM(I,J)=ESM(I,J) + (S11+S22)*CNST*k_c 
                    
                ENDDO
                 ESK(i,i)= ESK(i,i) + masa_mini(i)*cnst*landa
             
                 factor= qe*phi(i)*cnst + coef1*T_b*phi(i)*cnst + coef_pul*campoxl*sigma_el*phi(i)*cnst
            
                 EF(I)=EF(I) +factor

              ENDDO

           
           
           endif
          
          
          
          
          enddo
       enddo 



if(ncaso>=3) then !
    ! modificado para el tiempo
    DO K=1,NOPE
       SUM=0.0
       DO J=1,nope
           SUM      = SUM + (ESK(K,J) - A2*ESM(K,J)) * des(J)  
           ESM(K,J) = ESK(K,J) + A1*ESM(K,J)
       enddo
       efaux=EF(K)
       EF(K) = SUM + (A1*EF(K)+A2*EF_ant(K))  
!       EF(K) = SUM + (A1*EF(K)+A2*EF(K))  
       EF_ant(K)=efaux
    enddo

endif


end subroutine armado16


SUBROUTINE ARMADO8(NCASE,NLE,X,Y,Z,ns,NOPE,ESM,EF,sigma_el,qe,sol,Ex_el,Ey_el,Ez_el)
implicit none
INTEGER :: nope,NS(NOPE),NCASE
DOUBLE PRECISION :: X(NOPE),Y(NOPE),Z(NOPE),EF(NOPE),ESM(NOPE,NOPE),sigma_el,qe,sol(nope),Ex_el,Ey_el,Ez_el

DOUBLE PREcIsION:: PHI(8),DPHIX(8),DPHIY(8),DPHIZ(8),AJACO(3,3),AJACOI(3,3),DXHI(8), &
     DTHE(8),DPSI(8),GAUSSPT(2),GAUSSWT(2),XHI,THE,PSI,DETER,S11,S22,S33,CNST

integer kk,jj,i,j,ii,K,NLE,pgaus,I2
	 
DATA GAUSSPT/ -0.57735027, 0.57735027/
DATA GAUSSWT/ 1.0, 1.0 /


ESM=0.0
EF=0.0
  pgaus=0
  DO KK=1,2
    DO JJ=1,2
      DO II=1,2
        pgaus=pgaus+1

		XHI = GAUSSPT(KK)
        THE = GAUSSPT(JJ)
        PSI = GAUSSPT(II)
      
!  CALCULO LA FUNCION PHI, SU DERIVADA EN FUNCCION DE XHI Y DE THE

        call funciones8(XHI,THE,PSI,PHI,DPHIX,DPHIY,DPHIZ,AJACO,AJACOI,DXHI,DTHE,DPSI ) !!ojo orden distinto!!
       ! call funcORDENALYA(XHI,THE,PSI,PHI,DPHIX,DPHIY,DPHIZ,AJACO,AJACOI,DXHI,DTHE,DPSI )
  
!   JACOBIANO Y SU INVERSA
        AJACO=0.0

          DO K=1,NOPE
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

          IF(DETER.EQ.0.0) THEN
              WRITE(6,*) ' SUBRUTINE: ARMADO '
              WRITE(6,*) ' ELEMENTO: ',NLE,' PUNTO ',KK,JJ,II,' TIENE DETERMINANTE CERO',DETER
              STOP ' '
          ENDIF

 
!   DERIVADAS DE LAS FUNCIONES DE INTERPOLACION EN FUNCION DE X-Y-Z

          DO I=1,NOPE
            DPHIX(I)=AJACOI(1,1)*DXHI(I) + AJACOI(1,2)*DTHE(I) + AJACOI(1,3)*DPSI(I)
            DPHIY(I)=AJACOI(2,1)*DXHI(I) + AJACOI(2,2)*DTHE(I) + AJACOI(2,3)*DPSI(I)
            DPHIZ(I)=AJACOI(3,1)*DXHI(I) + AJACOI(3,2)*DTHE(I) + AJACOI(3,3)*DPSI(I)
          ENDDO

!       DIFERENCIAL DE INTEGRACION
          CNST = DETER * GAUSSWT(JJ)*GAUSSWT(KK)*GAUSSWT(II)

!     MATRIZ DE RIGIDEZ LOCAL Y VETOR INDEPENDIENTE

          DO I=1,nope
            DO J=1,nope
              S11= DPHIX(I)*DPHIX(J)
              S22= DPHIY(I)*DPHIY(J)
              S33= DPHIZ(I)*DPHIZ(J)

              ESM(I,J)=ESM(I,J) + (S11+S22+S33)*CNST*sigma_el
            ENDDO
            Ex_el=Ex_el- DPHIX(i)*sol(i)/real(nope)
            Ey_el=Ey_el- DPHIY(i)*sol(i)/real(nope)
            Ez_el=Ez_el- DPHIZ(i)*sol(i)/real(nope)

            !EF(I)=EF(I) + qe*phi(i)*cnst

          ENDDO
        ENDDO
      ENDDO
    ENDDO


    

end subroutine armado8


SUBROUTINE ARMADO4(NCASE,NLE,X,Y,Z,ns,NOPE,ESM,EF,sigma_el,qe,sol,Ex_el,Ey_el,Ez_el)
implicit none
INTEGER :: nope,NS(NOPE),NCASE
DOUBLE PRECISION :: X(NOPE),Y(NOPE),Z(NOPE),EF(NOPE),ESM(NOPE,NOPE),sigma_el,qe,sol(nope),Ex_el,Ey_el,Ez_el

DOUBLE PREcIsION:: PHI(4),DPHIX(4),DPHIY(4),DPHIZ(4),AJACO(3,3),AJACOI(3,3),DXHI(4), &
     DTHE(4),DPSI(4),GAUSSPT(4),GAUSSWT(4),XHI,THE,PSI,DETER,S11,S22,S33,CNST
DOUBLE PREcIsION::  posgp(3,nope),weigp(nope),sl,tl,zl

integer kk,jj,i,j,ii,K,NLE,pgaus,I2
	 
!DATA GAUSSPT/ -0.57735027, 0.57735027/
!DATA GAUSSWT/ 1.0, 1.0 /


ESM=0.0
EF=0.0
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

DO KK=1,nope
    
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

          DO K=1,NOPE
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

          IF(DETER.EQ.0.0) THEN
              WRITE(6,*) ' SUBRUTINE: ARMADO '
              WRITE(6,*) ' ELEMENTO: ',NLE,' PUNTO ',KK,JJ,II,' TIENE DETERMINANTE CERO',DETER
              STOP ' '
          ENDIF
 
!   DERIVADAS DE LAS FUNCIONES DE INTERPOLACION EN FUNCION DE X-Y-Z

          DO I=1,NOPE
              DPHIX(I)=AJACOI(1,1)*DXHI(I) + AJACOI(1,2)*DTHE(I) + AJACOI(1,3)*DPSI(I)
              DPHIY(I)=AJACOI(2,1)*DXHI(I) + AJACOI(2,2)*DTHE(I) + AJACOI(2,3)*DPSI(I)
              DPHIZ(I)=AJACOI(3,1)*DXHI(I) + AJACOI(3,2)*DTHE(I) + AJACOI(3,3)*DPSI(I)
          ENDDO

!       DIFERENCIAL DE INTEGRACION
          CNST = DETER * weigp(kk)

!     MATRIZ DE RIGIDEZ LOCAL Y VETOR INDEPENDIENTE

          DO I=1,nope
            DO J=1,nope
              S11= DPHIX(I)*DPHIX(J)
              S22= DPHIY(I)*DPHIY(J)
              S33= DPHIZ(I)*DPHIZ(J)

              ESM(I,J)=ESM(I,J) + (S11+S22+S33)*CNST*sigma_el
            ENDDO
            Ex_el=Ex_el- DPHIX(i)*sol(i)/real(nope)
            Ey_el=Ey_el- DPHIY(i)*sol(i)/real(nope)
            Ez_el=Ez_el- DPHIZ(i)*sol(i)/real(nope)

            EF(I)=EF(I) + qe*phi(i)*cnst

          ENDDO
ENDDO


    

end subroutine armado4




SUBROUTINE ARMADO(NCASE,NLE,X,Y,Z,ns,NOPE,ESM,EF,sigma_el,qe,OrdenAlya)
implicit none
INTEGER :: nope,NS(NOPE),NCASE,OrdenAlya
DOUBLE PRECISION :: X(NOPE),Y(NOPE),Z(NOPE),EF(NOPE),ESM(NOPE,NOPE),sigma_el,qe

DOUBLE PREcIsION:: PHI(27),DPHIX(27),DPHIY(27),DPHIZ(27),AJACO(3,3),AJACOI(3,3),DXHI(27), &
     DTHE(27),DPSI(27),GAUSSPT(3),GAUSSWT(3),XHI,THE,PSI,DETER,S11,S22,S33,CNST

integer kk,jj,i,j,ii,K,NLE,pgaus,I2
	 
DATA GAUSSPT/ -0.774596669241483377035853079956, 0.0, 0.774596669241483377035853079956/
DATA GAUSSWT/ 0.5555555556, 0.8888888889 ,0.5555555556 /

ESM=0.0
EF=0.0
  pgaus=0
  DO KK=1,3
    DO JJ=1,3
      DO II=1,3
        pgaus=pgaus+1

		XHI = GAUSSPT(KK)
        THE = GAUSSPT(JJ)
        PSI = GAUSSPT(II)
      
!  CALCULO LA FUNCION PHI, SU DERIVADA EN FUNCCION DE XHI Y DE THE
        if(OrdenAlya==1) then
           call funcORDENALYA(XHI,THE,PSI,PHI,DPHIX,DPHIY,DPHIZ,AJACO,AJACOI,DXHI,DTHE,DPSI )
        else
           call funciones(XHI,THE,PSI,PHI,DPHIX,DPHIY,DPHIZ,AJACO,AJACOI,DXHI,DTHE,DPSI ) !!ojo orden distinto!!
        endif
  
!   JACOBIANO Y SU INVERSA
        AJACO=0.0

          DO K=1,NOPE
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

          IF(DETER.EQ.0.0) THEN
              WRITE(6,*) ' SUBRUTINE: ARMADO '
              WRITE(6,*) ' ELEMENTO: ',NLE,' PUNTO ',KK,JJ,II,' TIENE DETERMINANTE CERO',DETER
              STOP ' '
          ENDIF

 
!   DERIVADAS DE LAS FUNCIONES DE INTERPOLACION EN FUNCION DE X-Y-Z

          DO I=1,27
            DPHIX(I)=AJACOI(1,1)*DXHI(I) + AJACOI(1,2)*DTHE(I) + AJACOI(1,3)*DPSI(I)
            DPHIY(I)=AJACOI(2,1)*DXHI(I) + AJACOI(2,2)*DTHE(I) + AJACOI(2,3)*DPSI(I)
            DPHIZ(I)=AJACOI(3,1)*DXHI(I) + AJACOI(3,2)*DTHE(I) + AJACOI(3,3)*DPSI(I)
          ENDDO

!       DIFERENCIAL DE INTEGRACION
          CNST = DETER * GAUSSWT(JJ)*GAUSSWT(KK)*GAUSSWT(II)

!     MATRIZ DE RIGIDEZ LOCAL Y VETOR INDEPENDIENTE

          DO I=1,nope
            DO J=1,nope
              S11= DPHIX(I)*DPHIX(J)
              S22= DPHIY(I)*DPHIY(J)
              S33= DPHIZ(I)*DPHIZ(J)

              ESM(I,J)=ESM(I,J) + (S11+S22+S33)*CNST*sigma_el

            ENDDO

            EF(I)=EF(I) + qe*phi(i)*cnst

          ENDDO
        ENDDO
      ENDDO
    ENDDO


end subroutine armado


SUBROUTINE ARMADObio(NCASE,NLE,X,Y,Z,des,ns,nodpel,NOPE,ESM,EF,sigma_el,qe,landa,a1,a2,k_c,coef1,t_b,coef_pul,campoxl,demed,masa_mini,ef_ant)
implicit none
INTEGER :: nodpel,nope,NS(nodpel),NCASE
DOUBLE PRECISION :: X(nodpel),Y(nodpel),Z(nodpel),des(nodpel),EF(nodpel),ESM(nodpel,nodpel),ESK(nodpel,nodpel),sigma_el,qe,landa,a1,a2, &
    &               k_c,coef1,t_b,coef_pul,campoxl,masa_mini(nodpel),ef_ant(nodpel)

integer kk,jj,i,j,ii,K,NLE,pgaus,I2
	 
DOUBLE PREcIsION:: PHI(nodpel),DPHIX(nodpel),DPHIY(nodpel),DPHIZ(nodpel),AJACO(3,3),AJACOI(3,3),DXHI(nodpel), &
     DTHE(nodpel),DPSI(nodpel),GAUSSPT(nope),GAUSSWT(nope),XHI,THE,PSI,DETER,S11,S22,S33,CNST,SUM,demed,factor,efaux




if(nope==2) then
   GAUSSPT(1)= -0.57735027
   GAUSSPT(2)=  0.57735027
   GAUSSWT(1)=0.
   GAUSSWT(2)=0.

elseif(nope==3) then
   GAUSSPT(1)= -0.774596669241483377035853079956
   GAUSSPT(2)= 0.0
   GAUSSPT(3)= 0.774596669241483377035853079956

   GAUSSWT(1)= 0.5555555556
    GAUSSWT(2)=0.8888888889 
    GAUSSWT(3)=0.5555555556 

endif

ESM=0.0
ESK=0.0
EF=0.0
  pgaus=0
  DO KK=1,nope
    DO JJ=1,nope
      DO II=1,nope
        pgaus=pgaus+1

		XHI = GAUSSPT(II)
        THE = GAUSSPT(JJ)
        PSI = GAUSSPT(KK)
      
!  CALCULO LA FUNCION PHI, SU DERIVADA EN FUNCCION DE XHI Y DE THE

        if(nodpel==8) then
            call funciones8(XHI,THE,PSI,PHI,DPHIX,DPHIY,DPHIZ,AJACO,AJACOI,DXHI,DTHE,DPSI ) 
        elseif(nodpel==27) then
            call funcORDENALYA(XHI,THE,PSI,PHI,DPHIX,DPHIY,DPHIZ,AJACO,AJACOI,DXHI,DTHE,DPSI )
        endif
  
!   JACOBIANO Y SU INVERSA
        AJACO=0.0

          DO K=1,nodpel
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

          IF(DETER.EQ.0.0) THEN
              WRITE(6,*) ' SUBRUTINE: ARMADO '
              WRITE(6,*) ' ELEMENTO: ',NLE,' PUNTO ',KK,JJ,II,' TIENE DETERMINANTE CERO',DETER
              STOP ' '
          ENDIF

 
!   DERIVADAS DE LAS FUNCIONES DE INTERPOLACION EN FUNCION DE X-Y-Z

          DO I=1,nodpel
            DPHIX(I)=AJACOI(1,1)*DXHI(I) + AJACOI(1,2)*DTHE(I) + AJACOI(1,3)*DPSI(I)
            DPHIY(I)=AJACOI(2,1)*DXHI(I) + AJACOI(2,2)*DTHE(I) + AJACOI(2,3)*DPSI(I)
            DPHIZ(I)=AJACOI(3,1)*DXHI(I) + AJACOI(3,2)*DTHE(I) + AJACOI(3,3)*DPSI(I)
          ENDDO

!       DIFERENCIAL DE INTEGRACION
          CNST = DETER * GAUSSWT(JJ)*GAUSSWT(KK)*GAUSSWT(II)


!     MATRIZ DE RIGIDEZ LOCAL Y VETOR INDEPENDIENTE

          DO I=1,nodpel
            DO J=1,nodpel
              S11= DPHIX(I)*DPHIX(J)
              S22= DPHIY(I)*DPHIY(J)
              S33= DPHIZ(I)*DPHIZ(J)

              ESM(I,J)=ESM(I,J) + (S11+S22+S33)*CNST*k_c 
              
              ESK(I,J)=ESK(I,J) + landa * PHI(i)*PHI(j) *cnst  ! caso tejido vivo

            ENDDO
           ! matrix de masa
             ESK(i,i)= ESK(i,i) + masa_mini(i)*LANDA*CNST
             
             factor= (qe +coef1*(T_b) + coef_pul*campoxl*sigma_el)*phi(i)*cnst ! caso tejido vivo
             
             !factor= (qe + campoxl*sigma_el)*phi(i)*cnst 
            
             EF(I)=EF(I) + factor

          ENDDO


        ENDDO
      ENDDO
    ENDDO


! modificado para el tiempo
DO K=1,nodpel
   SUM=0.0
   DO J=1,nodpel
       SUM      = SUM + (ESK(K,J) - A2*ESM(K,J)) * des(J)  
       ESM(K,J) = ESK(K,J) + A1*ESM(K,J)
   enddo
!   EF(K) = SUM + (A1*EF(K)+A2*EF(K))  
   efaux=ef(k)
   EF(K) = SUM + (A1*EF(K)+A2*EF_ANT(K))  
   ef_ant(k)=efaux
enddo


end subroutine armadobio



SUBROUTINE ARMADObio_tetra(NCASE,NLE,X,Y,Z,des,ns,nodpel,NOPE,ESM,EF,sigma_el,qe,landa,a1,a2,k_c,coef1,t_b,coef_pul,campoxl, &
   &                       demed,masa_mini,ef_ant)
implicit none
INTEGER :: nodpel,nope,NS(nodpel),NCASE
DOUBLE PRECISION :: X(nodpel),Y(nodpel),Z(nodpel),des(nodpel),EF(nodpel),ESM(nodpel,nodpel),ESK(nodpel,nodpel),sigma_el,qe,landa,a1,a2, &
    &               k_c,coef1,t_b,coef_pul,campoxl,masa_mini(nodpel),ef_ant(nodpel)

integer kk,jj,i,j,ii,K,NLE,pgaus,I2
	 DOUBLE PREcIsION::  posgp(3,nodpel),weigp(nodpel),sl,tl,zl

DOUBLE PREcIsION:: PHI(nodpel),DPHIX(nodpel),DPHIY(nodpel),DPHIZ(nodpel),AJACO(3,3),AJACOI(3,3),DXHI(nodpel), &
     DTHE(nodpel),DPSI(nodpel),GAUSSPT(nodpel),GAUSSWT(nodpel),XHI,THE,PSI,DETER,S11,S22,S33,CNST,SUM,demed,factor,efaux


ESM=0.0
EF=0.0
ESK=0.0
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

DO KK=1,nodpel
    
	 sl = posgp(1,kk)
     tl = posgp(2,kk)
     zl = posgp(3,kk)


!  CALCULO LA FUNCION PHI, SU DERIVADA EN FUNCCION DE XHI Y DE THE
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

          DO K=1,nodpel
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

          IF(DETER.LE.0.0) THEN
              WRITE(6,*) ' SUBRUTINE: ARMADO '
              WRITE(6,*) ' ELEMENTO: ',NLE,' PUNTO ',KK,JJ,II,' TIENE DETERMINANTE CERO',DETER
              STOP ' '
          ENDIF

 
!   DERIVADAS DE LAS FUNCIONES DE INTERPOLACION EN FUNCION DE X-Y-Z

          DO I=1,nodpel
            DPHIX(I)=AJACOI(1,1)*DXHI(I) + AJACOI(1,2)*DTHE(I) + AJACOI(1,3)*DPSI(I)
            DPHIY(I)=AJACOI(2,1)*DXHI(I) + AJACOI(2,2)*DTHE(I) + AJACOI(2,3)*DPSI(I)
            DPHIZ(I)=AJACOI(3,1)*DXHI(I) + AJACOI(3,2)*DTHE(I) + AJACOI(3,3)*DPSI(I)
          ENDDO

!       DIFERENCIAL DE INTEGRACION
          CNST = DETER * weigp(kk)

!     MATRIZ DE RIGIDEZ LOCAL Y VETOR INDEPENDIENTE

          DO I=1,nodpel
            DO J=1,nodpel
              S11= DPHIX(I)*DPHIX(J)
              S22= DPHIY(I)*DPHIY(J)
              S33= DPHIZ(I)*DPHIZ(J)

              ESM(I,J)=ESM(I,J) + (S11+S22+S33)*CNST*k_c 
              
              ESK(I,J)=ESK(I,J) + landa * PHI(i)*PHI(j) *cnst  

            ENDDO
           ! matrix de masa
             !ESK(i,i)= ESK(i,i) + masa_mini(i)*LANDA*CNST
             
             
             factor= (qe + campoxl*sigma_el)*phi(i)*cnst 
            ! factor= qe*phi(i)*cnst 
            
             EF(I)=EF(I) + factor
            

          ENDDO
ENDDO

! modificado para el tiempo
DO K=1,nodpel
   SUM=0.0
   DO J=1,nodpel
       SUM      = SUM + (ESK(K,J) - A2*ESM(K,J)) * des(J)  
       ESM(K,J) = ESK(K,J) + A1*ESM(K,J)
   enddo
   efaux=ef(k)
   EF(K) = SUM + (A1*EF(K)+A2*EF_ANT(K))  
   ef_ant(k)=efaux
enddo


    end subroutine armadobio_tetra


    
    
subroutine funciones4(s,t,PHI,DPHIX,DPHIY,AJACO,AJACOI,DXHI,DTHE)   
IMPLICIT none
double precision PHI(4),DPHIX(4),DPHIY(4),AJACO(2,2),AJACOI(2,2),DXHI(4),DTHE(4),s,t
double precision sm,tm,zm,sq,tp,zp


             
             
                 sm = 0.5*(1.0-s)
                 tm = 0.5*(1.0-t)
                 sq = 0.5*(1.0+s)
                 tp = 0.5*(1.0+t)
                 

                  phi(1)    = sm*tm
                  DXHI(1)  =-0.5*tm
                  DTHE(1)  =-0.5*sm
                 

                 phi(2)     = sq*tm
                 dxhi( 2) = 0.5*tm
                 dthe( 2) =-0.5*sq
                 
                 phi(3)     = sq*tp
                 DXHI( 3) = 0.5*tp
                 dthe( 3) = 0.5*sq
                 
                 phi(4)     = sm*tp
                 DXHI( 4) =-0.5*tp
                 dthe( 4) = 0.5*sm
	 
	 
return
end

    

subroutine funciones8(s,t,z,PHI,DPHIX,DPHIY,DPHIZ,AJACO,AJACOI,DXHI,DTHE,DPSI )   
IMPLICIT none
double precision PHI(8),DPHIX(8),DPHIY(8),DPHIZ(8),AJACO(3,3),AJACOI(3,3),DXHI(8), &
          DTHE(8),DPSI(8),s,t,z
double precision sm,tm,zm,sq,tp,zp


             
             
                 sm = 0.5*(1.0-s)
                 tm = 0.5*(1.0-t)
                 zm = 0.5*(1.0-z)
                 sq = 0.5*(1.0+s)
                 tp = 0.5*(1.0+t)
                 zp = 0.5*(1.0+z)
                 

                  phi(1)    = sm*tm*zm
                  DXHI(1)  =-0.5*tm*zm
                  DTHE(1)  =-0.5*sm*zm
                  DPSI(1)=-0.5*sm*tm


                 phi(2)     = sq*tm*zm
                 dxhi( 2) = 0.5*tm*zm
                 dthe( 2) =-0.5*sq*zm
                 dpsi( 2) =-0.5*sq*tm
                 
                 phi(3)     = sq*tp*zm
                 DXHI( 3) = 0.5*tp*zm
                 dthe( 3) = 0.5*sq*zm
                 dpsi( 3) =-0.5*sq*tp
                 
                 phi(4)     = sm*tp*zm
                 DXHI( 4) =-0.5*tp*zm
                 dthe( 4) = 0.5*sm*zm
                 dpsi( 4) =-0.5*sm*tp
                 
                 phi( 5)    = sm*tm*zp
                 DXHI( 5) =-0.5*tm*zp
                 dthe( 5) =-0.5*sm*zp
                 dpsi( 5) = 0.5*sm*tm
                 
                 phi(6)    = sq*tm*zp 
                 DXHI(6) = 0.5*tm*zp
                 dthe( 6) =-0.5*sq*zp
                 dpsi( 6) = 0.5*sq*tm
                 
                 phi(  7)   = sq*tp*zp
                 DXHI(7) = 0.5*tp*zp
                 dthe( 7) = 0.5*sq*zp
                 dpsi( 7) = 0.5*sq*tp

                 phi(8)     = sm*tp*zp
                 DXHI( 8) =-0.5*tp*zp
                 dthe( 8) = 0.5*sm*zp
                 dpsi( 8) = 0.5*sm*tp

	 
	 
return
end


subroutine funciones(s,t,z,PHI,DPHIX,DPHIY,DPHIZ,AJACO,AJACOI,DXHI,DTHE,DPSI )   
IMPLICIT none
double precision PHI(27),DPHIX(27),DPHIY(27),DPHIZ(27),AJACO(3,3),AJACOI(3,3),DXHI(27), &
          DTHE(27),DPSI(27),s,t,z
double precision s1,z1,t1,sl,tl,zl,sq,tp,zp,s2,t2,z2,s3,t3,z3,s4,t4,z4


	 sl=s*(s-1.0)
     tl=t*(t-1.0)
     zl=z*(z-1.0)
     sq=s*(s+1.0)
     tp=t*(t+1.0)
     zp=z*(z+1.0)
     s1= 2.0*s-1.0
     t1= 2.0*t-1.0
     z1= 2.0*z-1.0
     s2= 1.0-s*s
     t2= 1.0-t*t
     z2= 1.0-z*z
     s3= 1.0+2.0*s
     t3= 1.0+2.0*t
     z3= 1.0+2.0*z
     s4=-2.0*s
     t4=-2.0*t
     z4=-2.0*z
     
	 phi (1) = 0.125*sl*tl*zl
     DXHI(1) = 0.125*s1*tl*zl
     DTHE(1) = 0.125*sl*t1*zl
     DPSI(1) = 0.125*sl*tl*z1
     
	 phi (2) = 0.25*s2*tl*zl
     DXHI(2) = 0.25*s4*tl*zl
     DTHE(2) = 0.25*s2*t1*zl
     DPSI(2) = 0.25*s2*tl*z1
     
     phi (3) = 0.125*sq*tl*zl
     DXHI(3) = 0.125*s3*tl*zl
     DTHE(3) = 0.125*sq*t1*zl
     DPSI(3) = 0.125*sq*tl*z1
     
 	 phi (4) = 0.25*sq*t2*zl
     DXHI(4) = 0.25*s3*t2*zl
     DTHE(4) = 0.25*sq*t4*zl
     DPSI(4) = 0.25*sq*t2*z1
     

	 phi (5) = 0.125*sq*tp*zl
     DXHI(5) = 0.125*s3*tp*zl
     DTHE(5) = 0.125*sq*t3*zl
     DPSI(5) = 0.125*sq*tp*z1
     
     phi (6) = 0.25*s2*tp*zl
     DXHI(6) = 0.25*s4*tp*zl
     DTHE(6) = 0.25*s2*t3*zl
     DPSI(6) = 0.25*s2*tp*z1

	 phi (7) = 0.125*sl*tp*zl
     DXHI(7) = 0.125*s1*tp*zl
     DTHE(7) = 0.125*sl*t3*zl
     DPSI(7) = 0.125*sl*tp*z1
       
     phi (8) = 0.25*sl*t2*zl
     DXHI(8) = 0.25*s1*t2*zl
     DTHE(8) = 0.25*sl*t4*zl
     DPSI(8) = 0.25*sl*t2*z1
     
     phi (9) = 0.5*s2*t2*zl
     DXHI(9) = 0.5*s4*t2*zl
     DTHE(9) = 0.5*s2*t4*zl
     DPSI(9) = 0.5*s2*t2*z1

		 
     phi (10) = 0.25*sl*tl*z2
     DXHI(10) = 0.25*s1*tl*z2
     DTHE(10) = 0.25*sl*t1*z2
     DPSI(10) = 0.25*sl*tl*z4

     phi (11) = 0.5*s2*tl*z2
     DXHI(11) = 0.5*s4*tl*z2
     DTHE(11) = 0.5*s2*t1*z2
     DPSI(11) = 0.5*s2*tl*z4
     
	 	      
	 phi (12) = 0.25*sq*tl*z2
     DXHI(12) = 0.25*s3*tl*z2
     DTHE(12) = 0.25*sq*t1*z2
     DPSI(12) = 0.25*sq*tl*z4
     
	 phi (13) = 0.5*sq*t2*z2
     DXHI(13) = 0.5*s3*t2*z2
     DTHE(13) = 0.5*sq*t4*z2
     DPSI(13) = 0.5*sq*t2*z4
     
	 
	 phi (14) = 0.25*sq*tp*z2
     DXHI(14) = 0.25*s3*tp*z2
     DTHE(14) = 0.25*sq*t3*z2
     DPSI(14) = 0.25*sq*tp*z4
     
	 phi (15) = 0.5*s2*tp*z2
     DXHI(15) = 0.5*s4*tp*z2
     DTHE(15) = 0.5*s2*t3*z2
     DPSI(15) = 0.5*s2*tp*z4
     
	 phi (16) = 0.25*sl*tp*z2
     DXHI(16) = 0.25*s1*tp*z2
     DTHE(16) = 0.25*sl*t3*z2
     DPSI(16) = 0.25*sl*tp*z4
     
	 phi (17) = 0.5*sl*t2*z2
     DXHI(17) = 0.5*s1*t2*z2
     DTHE(17) = 0.5*sl*t4*z2
     DPSI(17) = 0.5*sl*t2*z4
     
	 phi (18) = s2*t2*z2
     DXHI(18) = s4*t2*z2
     DTHE(18) = s2*t4*z2
     DPSI(18) = s2*t2*z4
	 
	 phi (19) = 0.125*sl*tl*zp
     DXHI(19) = 0.125*s1*tl*zp
     DTHE(19) = 0.125*sl*t1*zp
     DPSI(19) = 0.125*sl*tl*z3

     phi (20) = 0.25*s2*tl*zp
     DXHI(20) = 0.25*s4*tl*zp
     DTHE(20) = 0.25*s2*t1*zp
     DPSI(20) = 0.25*s2*tl*z3

	 phi (21) = 0.125*sq*tl*zp
     DXHI(21) = 0.125*s3*tl*zp
     DTHE(21) = 0.125*sq*t1*zp
     DPSI(21) = 0.125*sq*tl*z3
     
     phi (22) = 0.25*sq*t2*zp
     DXHI(22) = 0.25*s3*t2*zp
     DTHE(22) = 0.25*sq*t4*zp
     DPSI(22) = 0.25*sq*t2*z3
     
	 phi (23) = 0.125*sq*tp*zp
     DXHI(23) = 0.125*s3*tp*zp
     DTHE(23) = 0.125*sq*t3*zp
     DPSI(23) = 0.125*sq*tp*z3
     
	 phi (24) = 0.25*s2*tp*zp
     DXHI(24) = 0.25*s4*tp*zp
     DTHE(24) = 0.25*s2*t3*zp
     DPSI(24) = 0.25*s2*tp*z3
     
	 phi (25) = 0.125*sl*tp*zp
     DXHI(25) = 0.125*s1*tp*zp
     DTHE(25) = 0.125*sl*t3*zp
     DPSI(25) = 0.125*sl*tp*z3
     
	 phi (26) = 0.25*sl*t2*zp
     DXHI(26) = 0.25*s1*t2*zp
     DTHE(26) = 0.25*sl*t4*zp
     DPSI(26) = 0.25*sl*t2*z3
     
	 phi (27) = 0.5*s2*t2*zp
     DXHI(27) = 0.5*s4*t2*zp
     DTHE(27) = 0.5*s2*t4*zp
     DPSI(27) = 0.5*s2*t2*z3
	
	
return
end

subroutine funcORDENALYA(s,t,z,PHI,DPHIX,DPHIY,DPHIZ,AJACO,AJACOI,DXHI,DTHE,DPSI )   
IMPLICIT none
double precision PHI(27),DPHIX(27),DPHIY(27),DPHIZ(27),AJACO(3,3),AJACOI(3,3),DXHI(27), &
          DTHE(27),DPSI(27),s,t,z
double precision s1,z1,t1,sl,tl,zl,sq,tp,zp,s2,t2,z2,s3,t3,z3,s4,t4,z4

	 sl=s*(s-1.0)
     tl=t*(t-1.0)
     zl=z*(z-1.0)
     sq=s*(s+1.0)
     tp=t*(t+1.0)
     zp=z*(z+1.0)
     s1= 2.0*s-1.0
     t1= 2.0*t-1.0
     z1= 2.0*z-1.0
     s2= 1.0-s*s
     t2= 1.0-t*t
     z2= 1.0-z*z
     s3= 1.0+2.0*s
     t3= 1.0+2.0*t
     z3= 1.0+2.0*z
     s4=-2.0*s
     t4=-2.0*t
     z4=-2.0*z
     
	 phi (1) = 0.125*sl*tl*zl
     DXHI(1) = 0.125*s1*tl*zl
     DTHE(1) = 0.125*sl*t1*zl
     DPSI(1) = 0.125*sl*tl*z1
     
     phi (2) = 0.125*sq*tl*zl
     DXHI(2) = 0.125*s3*tl*zl
     DTHE(2) = 0.125*sq*t1*zl
     DPSI(2) = 0.125*sq*tl*z1
     
	 phi (3) = 0.125*sq*tp*zl
     DXHI(3) = 0.125*s3*tp*zl
     DTHE( 3) = 0.125*sq*t3*zl
     DPSI( 3) = 0.125*sq*tp*z1
     
	 phi ( 4) = 0.125*sl*tp*zl
     DXHI( 4) = 0.125*s1*tp*zl
     DTHE( 4) = 0.125*sl*t3*zl
     DPSI( 4) = 0.125*sl*tp*z1
     
	 phi ( 5) = 0.125*sl*tl*zp
     DXHI( 5) = 0.125*s1*tl*zp
     DTHE( 5) = 0.125*sl*t1*zp
     DPSI( 5) = 0.125*sl*tl*z3
     
	 phi ( 6) = 0.125*sq*tl*zp
     DXHI( 6) = 0.125*s3*tl*zp
     DTHE( 6) = 0.125*sq*t1*zp
     DPSI( 6) = 0.125*sq*tl*z3
     
	 phi ( 7) = 0.125*sq*tp*zp
     DXHI( 7) = 0.125*s3*tp*zp
     DTHE( 7) = 0.125*sq*t3*zp
     DPSI( 7) = 0.125*sq*tp*z3
     
	 phi ( 8) = 0.125*sl*tp*zp
     DXHI( 8) = 0.125*s1*tp*zp
     DTHE( 8) = 0.125*sl*t3*zp
     DPSI( 8) = 0.125*sl*tp*z3
     
	 phi ( 9) = 0.25*s2*tl*zl
     DXHI( 9) = 0.25*s4*tl*zl
     DTHE( 9) = 0.25*s2*t1*zl
     DPSI( 9) = 0.25*s2*tl*z1
     
	 phi (10) = 0.25*sq*t2*zl
     DXHI(10) = 0.25*s3*t2*zl
     DTHE(10) = 0.25*sq*t4*zl
     DPSI(10) = 0.25*sq*t2*z1
     
	 phi (11) = 0.25*s2*tp*zl
     DXHI(11) = 0.25*s4*tp*zl
     DTHE(11) = 0.25*s2*t3*zl
     DPSI(11) = 0.25*s2*tp*z1
     
	 phi (12) = 0.25*sl*t2*zl
     DXHI(12) = 0.25*s1*t2*zl
     DTHE(12) = 0.25*sl*t4*zl
     DPSI(12) = 0.25*sl*t2*z1
     
	 phi (13) = 0.25*sl*tl*z2
     DXHI(13) = 0.25*s1*tl*z2
     DTHE(13) = 0.25*sl*t1*z2
     DPSI(13) = 0.25*sl*tl*z4
     
	 phi (14) = 0.25*sq*tl*z2
     DXHI(14) = 0.25*s3*tl*z2
     DTHE(14) = 0.25*sq*t1*z2
     DPSI(14) = 0.25*sq*tl*z4
     
	 phi ( 15) = 0.25*sq*tp*z2
     DXHI(15) = 0.25*s3*tp*z2
     DTHE(15) = 0.25*sq*t3*z2
     DPSI(15) = 0.25*sq*tp*z4
     
	 phi ( 16) = 0.25*sl*tp*z2
     DXHI(16) = 0.25*s1*tp*z2
     DTHE(16) = 0.25*sl*t3*z2
     DPSI(16) = 0.25*sl*tp*z4
     
	 phi ( 17) = 0.25*s2*tl*zp
     DXHI(17) = 0.25*s4*tl*zp
     DTHE(17) = 0.25*s2*t1*zp
     DPSI(17) = 0.25*s2*tl*z3
     
	 phi ( 18) = 0.25*sq*t2*zp
     DXHI(18) = 0.25*s3*t2*zp
     DTHE(18) = 0.25*sq*t4*zp
     DPSI(18) = 0.25*sq*t2*z3
     
	 phi ( 19) = 0.25*s2*tp*zp
     DXHI(19) = 0.25*s4*tp*zp
     DTHE(19) = 0.25*s2*t3*zp
     DPSI(19) = 0.25*s2*tp*z3
     
	 phi ( 20) = 0.25*sl*t2*zp
     DXHI(20) = 0.25*s1*t2*zp
     DTHE(20) = 0.25*sl*t4*zp
     DPSI(20) = 0.25*sl*t2*z3
     
	 phi ( 21) = 0.5*s2*t2*zl
     DXHI(21) = 0.5*s4*t2*zl
     DTHE(21) = 0.5*s2*t4*zl
     DPSI(21) = 0.5*s2*t2*z1
     
	 phi ( 22) = 0.5*s2*tl*z2
     DXHI(22) = 0.5*s4*tl*z2
     DTHE(22) = 0.5*s2*t1*z2
     DPSI(22) = 0.5*s2*tl*z4
     
	 phi (23) = 0.5*sq*t2*z2
     DXHI(23) = 0.5*s3*t2*z2
     DTHE(23) = 0.5*sq*t4*z2
     DPSI(23) = 0.5*sq*t2*z4
     
	 phi (  24) = 0.5*s2*tp*z2
     DXHI(24) = 0.5*s4*tp*z2
     DTHE(24) = 0.5*s2*t3*z2
     DPSI(24) = 0.5*s2*tp*z4
     
	 phi ( 25) = 0.5*sl*t2*z2
     DXHI(25) = 0.5*s1*t2*z2
     DTHE(25) = 0.5*sl*t4*z2
     DPSI(25) = 0.5*sl*t2*z4
     
	 phi (  26) = 0.5*s2*t2*zp
     DXHI(26) = 0.5*s4*t2*zp
     DTHE(26) = 0.5*s2*t4*zp
     DPSI(26) = 0.5*s2*t2*z3
     
	 phi (27) = s2*t2*z2
     DXHI(27) = s4*t2*z2
     DTHE(27) = s2*t4*z2
     DPSI(27) = s2*t2*z4
   
return
end


subroutine DETERM(AJACO, AJACOI, DETER)   
IMPLICIT none
double precision :: deter, AJACO(3,3),AJACOI(3,3)
!local
double precision :: t1,t2,t3,denom
 
     t1  = AJACO(2,2)*AJACO(3,3) - AJACO(3,2)*AJACO(2,3)
     t2  =-AJACO(2,1)*AJACO(3,3) + AJACO(3,1)*AJACO(2,3)
     t3  = AJACO(2,1)*AJACO(3,2) - AJACO(3,1)*AJACO(2,2)
     deter = AJACO(1,1)*t1 + AJACO(1,2)*t2 + AJACO(1,3)*t3
     if(deter .eq. 0.0) return
     
	 denom = 1.0/deter
     AJACOI(1,1) = t1*denom
     AJACOI(2,1) = t2*denom
     AJACOI(3,1) = t3*denom
     AJACOI(2,2) = ( AJACO(1,1)*AJACO(3,3) - AJACO(3,1)*AJACO(1,3))*denom
     AJACOI(3,2) = (-AJACO(1,1)*AJACO(3,2) + AJACO(1,2)*AJACO(3,1))*denom
     AJACOI(3,3) = ( AJACO(1,1)*AJACO(2,2) - AJACO(2,1)*AJACO(1,2))*denom
     AJACOI(1,2) = (-AJACO(1,2)*AJACO(3,3) + AJACO(3,2)*AJACO(1,3))*denom
     AJACOI(1,3) = ( AJACO(1,2)*AJACO(2,3) - AJACO(2,2)*AJACO(1,3))*denom
     AJACOI(2,3) = (-AJACO(1,1)*AJACO(2,3) + AJACO(2,1)*AJACO(1,3))*denom


end subroutine determ

