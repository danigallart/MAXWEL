      PROGRAM TDFIELD 
C     PROGRAMA TDFIELD RESUELVE EL PROBLEMA GOBERNADO POR LA ECUACION
C                
C            -d(b.du/dx)/dx - d(b.du/dy)/dy + G.u = Q 
C     CON CONDICIONES DE CONTORNO DE DIRICHLET Y NEUMANN 
      COMMON/ELMATX/ESM(4,4),EST(4,4),EF(4),X(4),Y(4),ESU(4),ESV(4)
      COMMON/HCV/IDBC(50,2),DBC(50,2),NDBC
      COMMON/MATL/DXE,DYE,GE,QE,ALANDA
      COMMON/TLE/TITLE(20)
      COMMON/AV/A(200000),JGF,JGSM,NP,NBW,NN 
      
      DIMENSION NEL(500,4),XC(500),YC(500),MAT(500),TE(500)
      DIMENSION NS(4),MASA(200000),DES(1000)
      
      DIMENSION TELEM(1000),TEANT(1000),TEREAL(1000),TELEANT(1000)

      NTESTI=1
      CALL INICIO
      
      CALL LECTURA(NP,NE,XC,YC,NEL,MAT,NDBC,IDBC,DBC)

      INBW=0
      NBW =0
      DO 20 KK=1,NE
      
        DO I=1,4
          NS(I)=NEL(KK,I)
        ENDDO
      
      DO 21 I=1,3
        IJ=I+1
        DO21 J=IJ,4
        NB=IABS(NS(I)-NS(J))
        IF(NB.EQ.0) WRITE(5,26) KK
   26   FORMAT(/,'ELEMENTO',I3,'TIENE DOS NODOS CON EL MISMO NUMERO')
        IF(NB.LE.NBW) GOTO21
        INBW=KK
        NBW=NB
   21 CONTINUE
   20 CONTINUE
      
      NBW=NBW+1
      
      WRITE(5,27) NBW,INBW
   27 FORMAT(//,'BANDWIDTH ESTA EN',I4,'EN ELEMENTO',I4)
      WRITE(6,*) 'CALCULA EL BANDWITH',NBW,INBW
C
C     INICIALIZACION DEL VECTOR COLUMNA A()
C
C
      JGF=NP
      JGSM=JGF+NP
      JEND=JGSM+NP*NBW
      IF(JEND.GT.200000) THEN
        WRITE(6,23)
   23   FORMAT(//,'DIMENSION DEL VECTOR EXCEDIDA,EJECUCION TERMINADA')
        STOP
      ENDIF
      
      DO K=1,NP
        DES(K)=0.0
        TELEM(K)=0.0
        TEREAL(K)=650
        TEANT(K)=0.0
      ENDDO  
     
C     COMIENZO DEL LOOP TEMPORAL
      CALL INITIME(NTEM,TIME,NT,DT,ALF,EPSLN,NLIMITE,NPASO)
      
 1032 NPASO=NPASO+1
      
      IF( MOD(NPASO,10).EQ.0.0) THEN
        WRITE(6,*) ' PASO TEMPORAL:        ',NPASO
      ENDIF
      
      A2 = (1-ALF)*DT
      A1 = ALF*DT

C       SATEO DE LA NO LINEALIDAD
        NPANL=0
        EPSNL=0.001
        NTOTNL=50
        ERRORNL=1
        
        DO I=1,NE
          TELEANT(I)=TELEM(I)
        ENDDO  
        
        DO WHILE(ERRORNL.GT.EPSNL)
          
          NPANL=NPANL+1
C          WRITE(6,*) '&&&&&   &&&&&PASO NO LINEAL DE TEMPERATURA:',NPANL
          
          IF(NPANL.GT.NTOTNL) THEN
             WRITE(6,*) 'SOLUCION DEL PASO DE TIEMPO:',NPASO
             WRITE(6,*) 'NO ALCANZA LIMITE NO LINEAL:',NPANL
             STOP ' '
          ENDIF
          
          DO I=1,JEND
            MASA(I)=0.0
            A(I)=0.0
          ENDDO


C
C     GENERACION DEL SISTEMA DE ECUACIONES
          DO32 KK=1,NE
      
            DO I=1,4
              NS(I)=NEL(KK,I)
              J=NS(I)
              ESU(I)=DES(J)
              X(I)=XC(J)
              Y(I)=YC(J)
            ENDDO
      
            TEMP=TEREAL(KK)+TELEM(KK)
            DXE=FUNDX(MAT(KK),TEMP)
            DYE=FUNDY(MAT(KK),TEMP)
            GE=FUNGE(MAT(KK))
            QE=FUNQE(MAT(KK))
            ALANDA=FUNLANDA(MAT(KK),TEMP)
      
            CALL ELSTMF(KK,NTESTI)    
      
            CALL TIMER(NTEM,A1,A2,KK)
C
C    PROCEDIMIENTO DE RIGIDEZ DIRECTO
            DO 33 I=1,4
              II=NS(I)
              A(JGF+II)=A(JGF+II)+EF(I)
              DO34J=1,4
                JJ=NS(J)+1-II
                IF(JJ.LE.0) GOTO34
                J1=JGSM+(JJ-1)*NP+II-(JJ-1)*(JJ-2)/2
                A(J1)=A(J1)+ESM(I,J)
                MASA(J1-JGSM)=MASA(J1-JGSM) + EST(I,J)
   34         CONTINUE
   33       CONTINUE
   32     CONTINUE   
C
C     ****************************************************
C     * MODIFICACION Y SOLUCION DEL SISTEMA DE ECUACIONES*
C     * LOS DATOS VIENEN DE LA SUBRUTINA MODIFY          *
C     ****************************************************

C      WRITE(6,*) 'ENSAMBLO '

          CALL MODIFY
      
C      WRITE(6,*) ' VIENE DE MODIFY Y LLAMA A DCMPBD'

          CALL DCMPBD

C      WRITE(6,*) ' VIENE DE DCMPBD Y LLAMA A SLVBD'

          CALL SLVBD

C      WRITE(6,*) ' SALE DE SLVBD'

C     SALIDA DE LAS SOLUCIONES 
C      
          
          DENOM=0
          ANUME=0
          DO I=1,NP
            DENOM=DENOM + A(I)*A(I)
            ANUME= ANUME + (A(I)-TEANT(I))**2
            TEANT(I)=TE(I)
            TE(I)=A(I)    
          ENDDO  

          DO J=1,NE
             TELEM(J)=(TE(NEL(J,1))+TE(NEL(J,2))+TE(NEL(J,3))+
     *              TE(NEL(J,4)))/4
          ENDDO  

          
          ERRORNL=SQRT(ANUME/DENOM)
          
        ENDDO  
        
C  ACA ENTRA EL CALCULO ELASTICO BASADO EN EL TELEM(KK) DEL PASO Y EL DELTA
C  DE TEMPERATURA PARA EL PASO TEMPORAL QUE SE RELACIONA AL CAMBIO DE DIMENSION       
        
        CALL ELASTICO(NP,NE,XC,YC,MAT,NEL,TELEM,TEREAL,TELEANT,NPASO)
         
        IF(MOD(NPASO,10).EQ.0) THEN
           WRITE(1,*)' VALORES DE LAS TEMPERATURAS NODALES'
           DO 113 I=1,NP
             WRITE(8,150) NPASO/10,I,TE(I)
  150        FORMAT(2I5,2X,E15.5)   
  113        WRITE(1,111) NPASO,I,TE(I),TEANT(I)
  111        FORMAT(2I3,3X,2E15.5)
           
           WRITE(1,*) ' VALORES DE LAS TEMPERATURAS ELEMENTALES '
           DO J=1,NE
              WRITE(1,'(I3,3X,E15.5)') J,TEREAL(J)+TELEM(J)
           ENDDO  
        ENDIF
        
        NT   = NT   + 1
        TIME = TIME + DT
C
        DIFF=0.0
        SOLN=0.0
        DO I=1,NP
           SOLN=SOLN + TE(I)*TE(I)
           DIFF=DIFF + (TE(I)-DES(I)) * (TE(I)-DES(I))
           DES(I)=TE(I)
        ENDDO   
      
        PERCNT=1.0
        IF(SOLN.NE.0.0) PERCNT = DSQRT(DIFF/SOLN)
      
        IF(PERCNT.LE.EPSLN) THEN
            WRITE(6,*) ' LA SOLUCION ALCANZO UN ESTADO ESTACIONARIO'
            WRITE(1,300) NPASO,TIME
 300        FORMAT(' ESTADO ESTACIONARIO PASO :',I3,'  TIEMPO:',E12.5)       
            CALL SALIDA(NP,NE,XC,YC,NEL,TE)
            STOP ' '
        ENDIF
      
      
        IF(NPASO.GE.NLIMITE) THEN
            CALL SALIDA(NP,NE,XC,YC,NEL,TE)
            STOP
        ELSE
           IF(MOD(NPASO,20).EQ.0) CALL SALIDA(NP,NE,XC,YC,NEL,TE)
           GOTO 1032
        ENDIF
   
      END
C     ***********************************************************
C     ***********************************************************      
C     *****************  SUBRUTINAS  ****************************
C     *****************  SUBRUTINAS  ****************************
C     ***********************************************************
C     ***********************************************************
C     *********************************************************** 
C     INICIO Y APERTURA DE ARCHIVOS
      SUBROUTINE INICIO
      
c      OPEN(UNIT=2,FILE='PRUT4.FEM',STATUS='OLD',ERR=3001)
c      OPEN (UNIT=3,FILE='PRUT4.TXT',STATUS='OLD',ERR=3002)

      OPEN(1,FILE='temper.txt',STATUS='UNKNOWN')
      
      OPEN(UNIT=2,FILE= 'UMO.FEM',STATUS='OLD',ERR=3001)
      OPEN (UNIT=3,FILE='UMO.TXT',STATUS='OLD',ERR=3002)

      OPEN(UNIT=4, FILE='Temdat4.TXT',STATUS='UNKNOWN',ERR=3004)
      OPEN(UNIT=8, FILE='TEMTIME4.NDP',STATUS='UNKNOWN',ERR=3005)
      
      OPEN(UNIT=9,FILE='SALIELS.NDP')
      OPEN(UNIT=10,FILE='DATELAS.TXT')
      
      OPEN(UNIT=12,FILE='STRESS.TXT')
      OPEN(UNIT=13,FILE='STRAIN.TXT')

      OPEN(UNIT=15,FILE='ZDUDAS.TXT')

      RETURN
 3001 WRITE(6,2001)
      STOP ' '
 3002 WRITE(6,2002)
      STOP ' '
 3004 WRITE(6,2004)
      STOP ' '
 3005 WRITE(6,2005)
      STOP ' '

 2001 FORMAT(/' ERROR ABRIENDO ARCHIVO DE .FEM')
 2002 FORMAT(/' ERROR ABRIENDO ARCHIVO CON LOS DATOS ')   
 2004 FORMAT(/' ERROR ABRIENDO ARCHIVO AUXILIAR')  
 2005 FORMAT(/' ERROR ABRIENDO ARCHIVO TEMTIME.NDP')  

      END
C
C
C     LECTURA DE DATOS
      SUBROUTINE LECTURA(NP,NE,XC,YC,NEL,MAT,NDBC,IDBC,DBC)
      CHARACTER*12 TEXTO
      LOGICAL FLAG 
      DIMENSION XC(500),YC(500),NEL(500,4),IDBC(50,2),DBC(50,2)
      DIMENSION MAT(500),MA(4)
C  
C
      WRITE(6,*) 'LEYENDO DATOS...'
C
      FLAG=.TRUE.
      DOWHILE(FLAG)
      READ(2,'(A12)',ERR=3010,END=3010) TEXTO
      IF(TEXTO.EQ.'*COORDINATES ') FLAG=.FALSE.
      END DO
      READ(2,*,ERR=3020,END=3020) NP
      K=0
      DO1I=1,NP
    1 READ(2,*,ERR=3021,END=3021) K,XC(I),YC(I)
C
C
C
C     
C     %SALIDA DEL TITULO Y DATOS DE ENCABEZAMIENTO
C     
      
      WRITE (1,4) NP
    4 FORMAT(' NUMERO DE NODOS NP =',I5)
      WRITE(1,11)
   11 FORMAT(1X,'NODO',3X,'COORDENADAS DE NODOS'/)
      WRITE(1,12) (I,XC(I),YC(I),I=1,NP)
   12 FORMAT(1X,I4,3X,2E15.5)
C
C
C    
c
C
      FLAG=.TRUE.
      DOWHILE(FLAG)
      READ(2,'(A12)',ERR=3013,END=3013) TEXTO
      IF(TEXTO.EQ.'*ELEMENT_GRO')FLAG=.FALSE.
      END DO
    
      READ(2,*,ERR=3022,END=3022) NMAT
      NE=0
      DO21 I=1,NMAT
      READ(2,*,ERR=3012,END=3012) N,NMA
      MA(I)=NMA
      DO22 J=NE+1,NE+NMA
   22 MAT(J)=I 
   21 NE=NE+NMA
      
      FLAG=.TRUE.
      DOWHILE(FLAG)
      READ(2,'(A10)',ERR=3014,END=3014) TEXTO
      IF(TEXTO.EQ.'*INCIDENCE') FLAG=.FALSE.
      END DO
      DO 100 KK=1,NE
  100 READ(2,*,ERR=3023,END=3023) N,(NEL(KK,I),I=1,4)
      WRITE (1,14) NE
   14 FORMAT(' NUMERO DE ELEMENTOS NE =',I5)
      WRITE(1,110)
  110 FORMAT(1X,'ELEMENTO',3X,'MATERIAL',3X,'NODOS')
      DO101 J=1,NE
  101 WRITE(1,120) J,MAT(J),(NEL(J,I),I=1,4)
  120 FORMAT(I3,3X,I3,3X,4I4)
C
C     AQUI SE LEEN LAS CONDICIONES DE NEUMAN. 
C     
C     NDBC ES EL NUMERO DE CONDICIONES DE CONTORNO DE NEUMAN
C     IDBC(I,1) GUARDA EL NUMERO DE ELEMENTO
C     IDBC(I,2) GUARDA EL LADO SOBRE EL QUE SE DEFINE LA CONDICION
C     DBC(I,1)   GUARDA EL COEFICIENTE M*LONGITUD DEL LADO
C     DBC(I,2)   FUARDA EL COEFICIENTE S*LONGITUD DEL LADO
C     D*dPHI/dx = -M*PHIb + S   CONDICION DE NEUMAN
C
C
C      READ(3,*) NDBC
      NDBC=0
      WRITE(1,300) NDBC
  300 FORMAT(/,'NUMERO DE CONDICIONES DE NEUMAN:',I2)

      IF(NDBC.NE.0) THEN 
      DO 301 I=1,NDBC
      READ(3,*) IDBC(I,1),IDBC(I,2),DBC(I,1),DBC(I,2)
  301 WRITE(1,302) IDBC(I,1),IDBC(I,2),DBC(I,1),DBC(I,2)
  302 FORMAT(/,'ELEMENTO:',I2,2X,'LADO?:',I2,2X,/,'ML:',E15.5,2X,/,
     *'SL:',E15.5)
      ENDIF
      
      RETURN
C
C       
3010  WRITE(6,2010)
      STOP ' '
3011  WRITE(6,2011)
      STOP ' '
3012  WRITE(6,2012)
      STOP ' '
3013  WRITE(6,2013)
      STOP ' '    
3014  WRITE(6,2014)
      STOP ' '
3020  WRITE(6,2020)
      STOP ' '
3021  WRITE(6,2021)
      STOP ' '
3022  WRITE(6,2022)
      STOP ' '
3023  WRITE(6,2023)
      STOP ' '
C    
2010  FORMAT(' ERROR LEYENDO TEXTO coordenadas')
2011  FORMAT(' ERROR LEYENDO TEXTO numero de nodo')
2012  FORMAT(' ERROR LEYENDO  NUMERO DE GRUPO')
2013  FORMAT(' ERROR LEYENDO TEXTO numero de elementos')
2014  FORMAT(' ERROR LEYENDO TEXTO nodo1,nodo2...')
2020  FORMAT(' ERROR LEYENDO numero de nodos')
2021  FORMAT(' ERROR LEYENDO coordenadas')
2022  FORMAT(' ERROR LEYENDO numero de elementos')
2023  FORMAT(' ERROR LEYENDO conectividad')
      END      
C
      SUBROUTINE SALIDA(NP,NE,XC,YC,NEL,TE)
      DIMENSION XC(500),YC(500),NEL(500,3),TE(500)
     
      OPEN(7,FILE='SALIDA4.NDP',STATUS='UNKNOWN')
      WRITE(7,100)
      WRITE(7,110) NP
      WRITE(7,120) (I,XC(I),YC(I),I=1,NP)
      WRITE(7,130) NE
      DO 1000 KK=1,NE
 1000 WRITE(7,140) (NEL(KK,I),I=1,4)
      WRITE(7,150)
      WRITE(7,*)'1'
      WRITE(7,*) '1 UNO'
      WRITE(7,*)'       '
      DO I=1,NP
      WRITE (7,*) I,TE(I)
      END DO
      RETURN
  100 FORMAT('*DIMENSION'/'2')
  110 FORMAT('*COORDINATES'/I6)
  120 FORMAT( I6,1X,E12.6,1X,E12.6)
  130 FORMAT('*ELEMENT GROUPS'/'1'/'1',I6,'  Tri3'/'*INCIDENCES'/)
  140 FORMAT(4I6)
  150 FORMAT(/'*SCALAR FIELDS')
      END
C
C
C    %SUBRUTINA ELSTMF, CALCULA LAS MATRICES DE CADA ELEMENTO
      SUBROUTINE ELSTMF(NLE,NTESTI)
      COMMON/ELMATX/ESM(4,4),EST(4,4),EF(4),X(4),Y(4),ESU(4),ESV(4)
      COMMON/HCV/IDBC(50,2),DBC(50,2),NDBC
      COMMON/MATL/DXE,DYE,GE,QE,ALANDA

      DIMENSION PHI(4),DPHIX(4),DPHIY(4),AJACO(2,2),AJACOI(2,2),DXHI(4)
     *,DTHE(4),GAUSSPT(2),GAUSSWT(2),XNODE(4),YNODE(4),S11T(4,4),
     *S22T(4,4)
      
      DATA GAUSSPT/ -0.57735027, 0.57735027 /
      DATA GAUSSWT/ 1.0, 1.0 /    
      DATA XNODE / -1, 1, 1, -1 /
      DATA YNODE / -1, -1, 1, 1 /


C      WRITE(4,*) NLE,DXE,DYE,GE,QE,ALANDA
       
      DO I=1,4
        DO J=1,4
          ESM(I,J) =0.0
          EST(I,J) =0.0
          S11T(I,J)=0.0
          S22T(I,J)=0.0
        ENDDO
          EF(I)=0.0
      ENDDO

      DO KK=1,2
        DO JJ=1,2
C    XHI-THE PUNTOS DE INTEGRACION DE GAUSS
          XHI = GAUSSPT(KK)
          THE = GAUSSPT(JJ)
          
C    CALCULO LA FUNCION PHI, SU DERIVADA EN FUNCCION DE XHI Y DE THE      

          DO II=1,4
            XHI0 = 1.0 + XHI*XNODE(II)
            THE0 = 1.0 + THE*YNODE(II)
            PHI(II) = 0.25 * XHI0*THE0
            DXHI(II)= 0.25 * XNODE(II)*THE0
            DTHE(II)= 0.25 * YNODE(II)*XHI0
          ENDDO

C   JACOBIANO Y SU INVERSA      
          DO J=1,2  
            DO I=1,2
              AJACO(I,J)=0.0
              AJACO(I,J)=0.0
            ENDDO
          ENDDO    
            
          DO K=1,4
              AJACO(1,1)=AJACO(1,1)+DXHI(K)*X(K) 
              AJACO(1,2)=AJACO(1,2)+DXHI(K)*Y(K)
              AJACO(2,1)=AJACO(2,1)+DTHE(K)*X(K) 
              AJACO(2,2)=AJACO(2,2)+DTHE(K)*Y(K)
          ENDDO
          
          
          DETER= AJACO(1,1)*AJACO(2,2)-AJACO(1,2)*AJACO(2,1)
          

          IF(DETER.EQ.0.0) THEN
              WRITE(6,*) ' ELEMENTO: ',NLE,' PUNTO ',KK,JJ,' TIENE 
     *                         DETERMINANTE CERO'
              STOP ' '
          ENDIF

          AJACOI(1,1)=AJACO(2,2)/DETER
          AJACOI(2,2)=AJACO(1,1)/DETER
          AJACOI(1,2)=-AJACO(1,2)/DETER
          AJACOI(2,1)=-AJACO(2,1)/DETER
          
C   DERIVADAS DE LAS FUNCIONES DE INTERPOLACION EN FUNCION DE X-Y        
          
c          DO I=1,4
c            DPHIX(I)=0.0
c            DPHIY(I)=0.0
c          ENDDO  
          DO I=1,4
            DPHIX(I)=AJACOI(1,1)*DXHI(I) + AJACOI(1,2)*DTHE(I)
            DPHIY(I)=AJACOI(2,1)*DXHI(I) + AJACOI(2,2)*DTHE(I)
          ENDDO
          
C       DIFERENCIAL DE INTEGRACION
          CNST = DETER * GAUSSWT(JJ)*GAUSSWT(KK)

C     CX Y CY SIRVEN CUANDO LOS COEFICIENTES DX Y QE DEPENDEN DE LAS COORDENADAS    
c          CX=0.0
c          CY=0.0
c          
c          DO I=1,4
c            CX=CX+X(I)*PHI(I)
c            CY=CY+Y(I)*PHI(I)
c          ENDDO  
C     MATRIZ DE RIGIDEZ LOCAL VETOR INDEPENDIENTE  MATRIZ DE MASA
   
          DO I=1,4
            DO J=1,4
              S11= DPHIX(I)*DPHIX(J)*CNST
              S22= DPHIY(I)*DPHIY(J)*CNST
              
              S11T(I,J)=S11T(I,J) + S11
              S22T(I,J)=S22T(I,J) + S22
              
              ESM(I,J)=ESM(I,J)+ DXE*S11+DYE*S22
              
              EST(I,J)=EST(I,J)+PHI(I)*PHI(J)*CNST*ALANDA
            ENDDO
            
            EF(I)=EF(I) + CNST*PHI(I)*QE  
            
          ENDDO    
        ENDDO
      ENDDO       


      
C  GUARDA RESULTADOS EN FILE    
      IF(NLE.EQ.NTESTI) THEN
      
C        WRITE(4,10) NLE
C   10   FORMAT(2X,'ELEM:',I4,/,2X,'VECTOR FUERZA',2X,'MATRIZ RIGIDEZ')
C        DO I=1,4
C          WRITE(4,20) EF(I),(ESM(I,J),J=1,4)
C   20     FORMAT(2X,E12.5,2X,4E15.5) 
C        ENDDO   

C        WRITE(4,*) 'MATRIZ DE MASA :',NLE
   
C        DO I=1,4
C          WRITE(4,25) (EST(I,J),J=1,4)
C   25     FORMAT(2X,4E15.5) 
C        ENDDO   
        
        
C        DO I=1,4
C          WRITE(4,30) (S11T(I,J),J=1,4)
C   30     FORMAT(4E15.5) 
C        ENDDO
C        DO I=1,4
C          WRITE(4,40) (S22T(I,J),J=1,4)
C   40     FORMAT(4E15.5) 
C        ENDDO
      
      ENDIF
      
      
     
c
      RETURN
      
    5 WRITE(1,6) KK
    6 FORMAT(//2X,'AREA DEL ELEMENTO',I3,'ES MENOR QUE 0.0001')
      STOP
      END
C
C     SUBRUTINA TIMER
      SUBROUTINE TIMER(NTEM,A1,A2,NLE)
      COMMON/ELMATX/ESM(4,4),EST(4,4),EF(4),X(4),Y(4),ESU(4),ESV(4)
     
      IF(NTEM.NE.0) THEN
               
        DO 30 K=1,4
        SUM=0.0
        DO 20 J=1,4
        SUM   = SUM + (EST(K,J) - A2*ESM(K,J)) * ESU(J)  
  20    ESM(K,J) = EST(K,J) + A1*ESM(K,J)
  30    EF(K) = (A1+A2)*EF(K) + SUM

      
c      WRITE(4,8) NLE
c    8 FORMAT(2X,'DATOS TIMER',/,'ELEMENTO:',I4,/,
c     *2X,'VECTOR FUERZA',2X,'MATRIZ RIGIDEZ')
c      DO 9 I=1,4
c    9 WRITE(4,10) EF(I),(ESM(I,J),J=1,4)
c   10 FORMAT(2X,E12.5,2X,4E15.5)
c      WRITE(4,17) NLE
c   17 FORMAT('MATRIZ DE MASA',I3)  
c      DO19I=1,4
c   19 WRITE(4,21) (EST(I,J),J=1,4)
c  21  FORMAT(4E15.5)

      ENDIF

      RETURN
      END

C
C    SUBRUTINA MODIFY INSERTA LAS CONDICIONES DE CONTORNO
      SUBROUTINE MODIFY
      COMMON/AV/A(200000),JGF,JGSM,NP,NBW,NN
c
C
c      WRITE(1,201)
c  201 FORMAT(//,'FUENTES AND VALORES SINK')
      
      READ(3,*) N
      IF(N.NE.0) THEN

      DO I=1,N
        READ(3,*,ERR=3020,END=3020) NOD,TEM
        A(JGF+NOD)=A(JGF+NOD)+TEM
      ENDDO           
      
      ENDIF
c  202 WRITE(1,501) I,NOD,TEM
c  501 FORMAT(I6,1X,i6,1x,e15.5)   
  
C
C     %ENTRADA DE VALORES DE NODOS PRESCRIPTOS
C         
C          
C         ID - GRADOS DE LIBERTAD DE LOS DESPLAZAMIENTOS CONOCIDOS 
C         BD - VALORES DE LOS DESPLAZAMIENTOS
C
C
C      WRITE(1,208)
C  208 FORMAT(' VALORES DE LAS TEMPERATURAS CONOCIDAS') 
      READ(3,*) NDES
      DO 216 KK =1,NDES
      READ(3,*,ERR=3040,END=3040) ID,BD
C
C     MODIFICACION DE LA MATRIZ GLOBAL DE RIGIDEZ Y DE ELVECTOR FUERZA
C     GLOBAL USANDO EL METODO DE DELETEAR FILAS Y COLUMNAS
      K=ID-1
      DO211 J=2,NBW
      M=ID+J-1
      IF(M.GT.NP) GOTO210
      IJ=JGSM +(J-1)*NP+ID-(J-1)*(J-2)/2
      A(JGF+M)=A(JGF+M)-A(IJ)*BD
      A(IJ)=0.0
  210 IF(K.LE.0) GOTO 211
      KJ=JGSM +(J-1)*NP+K-(J-1)*(J-2)/2
      A(JGF+K)=A(JGF+K)-A(KJ)*BD
      A(KJ)=0.0
      K=K-1
  211 CONTINUE
      A(JGF+ID)=A(JGSM+ID)*BD
  221 CONTINUE   
  216 CONTINUE
C      WRITE(1,209) ID,BD
C  209 FORMAT(I3,E15.5)
C
      REWIND(3)
      RETURN
C
C
C     
 3020 WRITE(6,2020)
      STOP ' '
 3040 WRITE(6,2040)
      STOP' '
C
 2020 FORMAT(' ERROR LEYENDO GRADOS DE LIBERTAD Y VALORES DE LA FUERZA')
 2040 FORMAT(' ERROR LEYENDO GRADOS DE LIBERTAD DE DESPLAZAMIENTOS')
      END
C           
      SUBROUTINE DCMPBD
      COMMON/AV/A(200000),JGF,JGSM,NP,NBW,NN
c     DESCOMPOSICION DE UNA MATRIZ BANDED EN UNA TRIANGULAR SUPERIOR
C     USANDO EL METODO DE ELIMINACION GAUSSIANA
      NP1 = NP-1
      DO226I=1,NP1
      MJ = I+NBW-1
      IF(MJ.GT.NP) MJ=NP
      NJ = I+1
      MK = NBW
      IF((NP-I+1).LT.NBW) MK=NP-I+1
      ND=0
      DO 225 J=NJ,MJ
      MK = MK-1
      ND = ND+1
      NL = ND+1
      DO 225 K=1,MK
      NK = ND+K
      JK = JGSM + (K-1)*NP + J - ((K-1)*(K-2)/2.)
      INL= JGSM + (NL-1)*NP + I - ((NL-1)*(NL-2)/2.)
      INK= JGSM + (NK-1)*NP + I - ((NK-1)*(NK-2)/2.)
      II = JGSM + I
  225 A(JK)=A(JK) - (A(INL)*A(INK)/A(II))
  226 CONTINUE 
C     
      RETURN
      END
C
C
      SUBROUTINE SLVBD
      COMMON/AV/A(200000),JGF,JGSM,NP,NBW,NN
      NP1=NP-1
C
C     %DESCOMPOSICION DEL VECTOR FUERZA GLOBAL
C     
      DO 250 I=1,NP1
      MJ=I+NBW-1
      IF(MJ.GT.NP) MJ=NP
      NJ=I+1
      L=1
      DO 250 J=NJ,MJ
      L=L+1
      IL=JGSM + (L-1)*NP+I-(L-1)*(L-2)/2 
  250 A(JGF+J)=A(JGF+J)-A(IL)*A(JGF+I)/A(JGSM+I)
C
C
C     BACKWARD SUBSTITUTION PARA DETERMINAR LOS VALORES NODALES
C
C

      A(NP)=A(JGF+NP)/A(JGSM+NP)
      DO 252 K=1,NP1
      I=NP-K
      MJ=NBW
      IF((I+NBW-1).GT.NP) MJ=NP-I+1
      SUM=0.0
      DO 251 J=2,MJ
      N=I+J-1
      IJ=JGSM+(J-1)*NP+I-(J-1)*(J-2)/2
  251 SUM=SUM+A(IJ)*A(N)
  252 A(I)=(A(JGF+I)-SUM)/A(JGSM+I)
      
      RETURN
      END   
C
C
      SUBROUTINE INITIME(NTEM,TIME,NT,DT,ALF,EPSLN,NLIMITE,NPASO)

C     NTEM INDICA SI EL PROBLEMA ES TEMPORAL (NTME =1) O NO (NTIM=0)
      NTEM=1
C     TIME TIEMPO INICIAL
      TIME=0.0           
C     NT NO SE QUE ES 
      NT=0  
C     DT ES EL INTERVALO DE TIEMPO          
      DT=0.001                      
C     ALF COEFFICIENTE ALFA
      ALF=0.5                 
      EPSLN=0.00001
C     NLIMITE ES LA CANTIDAD DE PASOS MAXIMA PERMITIDA 
      NLIMITE = 2000                                    
      NPASO=0
      RETURN
      END


      SUBROUTINE ELASTICO(NP,NE,XC,YC,MAT,NEL,TELEM,TEREAL,TELEANT,
     * NPASO)
      
      DIMENSION XC(500),YC(500),MAT(500),NEL(500,4),XSOL(500)
      DIMENSION TELEM(1000),TEREAL(1000),TELEANT(1000)
      
      DIMENSION AEL(1000000),NS(8),NS2(4)
      DIMENSION ESM1(8,8),X1(4),Y1(4),ET(4),FT(8)
      DIMENSION DESPLA(100),NODOD(100),D(500,3,3)

      DIMENSION STRES(500,4,4),STRAIN(500,4,4),ALFA(500),STRE_EL(500,4)
      
      MAXI=1000000   
      INBWE=0
      NBWE =0
      DO KK=1,NE
         DO I=1,4
            NS2(I)=NEL(KK,I)
         ENDDO
         DO I=1,3
            IJ=I+1
            DO J=IJ,4
               NB=IABS(NS2(I)-NS2(J))
               IF(NB.GT.NBWE) THEN
                  INBWE=KK
                  NBWE =NB
               ENDIF   
            ENDDO
         ENDDO
      ENDDO
      
      NBWE=(NBWE+1)*2
      WRITE(10,30) NBWE,INBWE
  30  FORMAT(//,' BANDWIDTH ELASTICO ESTA EN',I4,'EN ELEMENTO',I4)

      NPE=2*NP
      JGFE=NPE
      JGSME=JGFE+NPE
      JENDE=JGSME+NPE*NBWE
      
      IF(JENDE.GT.MAXI) THEN
        WRITE(6,*) 'MAXI EXCEDIDO',JGSME,NBWE,NBWE*NPE,JENDE
        STOP' '
      ENDIF
      
      DO K=1,JENDE
        AEL(K)=0.0
      ENDDO  
      
C COND de CONTORNO 
C      WRITE(10,*) 'CONDICIONES DE CONTORNO '
      NDES=0
      DO K=1,NP
        IF(XC(K).EQ.0.0) THEN
          NDES=NDES+1
          DESPLA(NDES)=0
          NODOD(NDES)=2*K-1
C          write(10,*) NDES,NODOD(NDES),DESPLA(NDES)
        ENDIF
            
        IF(YC(K).EQ.0.0) THEN
          NDES=NDES+1
          DESPLA(NDES)=0
          NODOD(NDES)=2*K
C          write(10,*) NDES,NODOD(NDES),DESPLA(NDES)
        ENDIF
        
      ENDDO                   
      
      NTESTI=1
      
      DO KK=1,NE
        
        TEMP = TEREAL(KK) + TELEM(KK)
        
        ALFA(KK)=FUNAL(MAT(KK),TEMP)
        AL=ALFA(KK)
        T = TELEM(KK)-TELEANT(KK)
        
        EM=FUNEM(MAT(KK),TEMP)
        PR=FUNPOS(MAT(KK),TEMP)
        
C        WRITE(15,*) NPASO,KK,TEMP,T,AL,EM,PR
               
        DO I=1,4
          J=NEL(KK,I)
          NS(2*I-1)=J*2-1
          NS(2*I)=J*2
          X1(I)=XC(J)
          Y1(I)=YC(J)
        ENDDO

        CALL ELSTMX(KK,MAT(KK),D,NTESTI,ESM1,X1,Y1,ET,FT,AL,T,EM,PR)

        DO I=1,8
           II=NS(I)
           AEL(JGFE+II)=AEL(JGFE+II)+FT(I)
           DO J=1,8
              JJ=NS(J)+1-II
              IF(JJ.GT.0) THEN
                 J1=JGSME +(JJ-1)*NPE+II-(JJ-1)*(JJ-2)/2
                 AEL(J1)=AEL(J1)+ESM1(I,J)
              ENDIF
           ENDDO
        ENDDO
      
      ENDDO
      
      CALL MODIFYII(AEL,NDES,DESPLA,NODOD,NPE,JGFE,JGSME,NBWE)

      CALL DCMPBD_EL(AEL,JGFE,JGSME,NPE,NBWE)
      CALL SLVBD_EL(AEL,JGFE,JGSME,NPE,NBWE)
      
      DO J=1,NPE
        XSOL(J)=AEL(J)
      ENDDO  

      CALL CALCULO(NE,XC,YC,ALFA,TELEM,NEL,D,STRES,STRAIN,MAT,XSOL,
     *TELEANT,NPASO,NTESTI,NPE,TEREAL,STRE_EL) 

      IF(MOD(NPASO,10).EQ.0.0) THEN
          CALL SALIDAELASTICA(NPE,NE,XC,YC,NEL,STRE_EL,NPASO)
      ENDIF
      
      RETURN
      END

C     *********************SUBRUTINA ELSTMX************************
      SUBROUTINE ELSTMX(NLE,MAT,D,NTESTI,ESM1,X1,Y1,ET,FT,AL,T,EM,PR)

      DIMENSION ESM1(8,8),X1(4),Y1(4),ET(4),FT(8)
      DIMENSION D(500,3,3)
      
      DIMENSION PHI(4),DPHIX(4),DPHIY(4),AJACO(2,2),AJACOI(2,2),DXHI(4)
     *,DTHE(4),GAUSSPT(2),GAUSSWT(2),XNODE(4),YNODE(4)

      DATA GAUSSPT/ -0.57735027, 0.57735027 /
      DATA GAUSSWT/ 1.0, 1.0 /    
      DATA XNODE / -1, 1, 1, -1 /
      DATA YNODE / -1, -1, 1, 1 /

      
C    PLANE STRAIN
      R=EM/(1+PR)/(1-PR*2)
      D(NLE,1,1)=(1-PR)*R
      D(NLE,2,2)=D(NLE,1,1)
      D(NLE,3,3)=R*(1-2*PR)/2.
      D(NLE,1,2)=R*PR
      D(NLE,2,1)=D(NLE,1,2)
      D(NLE,1,3)=0.0
      D(NLE,3,1)=0.0
      D(NLE,2,3)=0.0
      D(NLE,3,2)=0.0
      
      DO I=1,8
        DO J=1,8
          ESM1(I,J) =0.0
        ENDDO
      ENDDO

C    %VECTOR STRESS TERMINO ET
       ET(1)=AL*T
       ET(2)=AL*T
       ET(3)=0  
           
      DO KK=1,2
        DO JJ=1,2
          XHI = GAUSSPT(KK)
          THE = GAUSSPT(JJ)
          
          DO II=1,4
            XHI0 = 1.0 + XHI*XNODE(II)
            THE0 = 1.0 + THE*YNODE(II)
            PHI(II) = 0.25 * XHI0*THE0
            DXHI(II)= 0.25 * XNODE(II)*THE0
            DTHE(II)= 0.25 * YNODE(II)*XHI0
          ENDDO

          DO J=1,2  
            DO I=1,2
              AJACO(I,J)=0.0
              AJACO(I,J)=0.0
            ENDDO
          ENDDO    
            
          DO K=1,4
              AJACO(1,1)=AJACO(1,1)+DXHI(K)*X1(K) 
              AJACO(1,2)=AJACO(1,2)+DXHI(K)*Y1(K)
              AJACO(2,1)=AJACO(2,1)+DTHE(K)*X1(K) 
              AJACO(2,2)=AJACO(2,2)+DTHE(K)*Y1(K)
          ENDDO
          
          DETER= AJACO(1,1)*AJACO(2,2)-AJACO(1,2)*AJACO(2,1)

          IF(DETER.EQ.0.0) THEN
              WRITE(6,*) ' ELEMENTO: ',NLE,' PUNTO ',KK,JJ,' TIENE 
     *                         DETERMINANTE CERO'
              STOP ' '
          ENDIF

          AJACOI(1,1)= AJACO(2,2)/DETER
          AJACOI(2,2)= AJACO(1,1)/DETER
          AJACOI(1,2)=-AJACO(1,2)/DETER
          AJACOI(2,1)=-AJACO(2,1)/DETER
          
          DO I=1,4
            DPHIX(I)=AJACOI(1,1)*DXHI(I) + AJACOI(1,2)*DTHE(I)
            DPHIY(I)=AJACOI(2,1)*DXHI(I) + AJACOI(2,2)*DTHE(I)
          ENDDO
          
          CNST = DETER * GAUSSWT(JJ)*GAUSSWT(KK)

          IW=1    
          DO I=1,4
            JW=1
            DO J=1,4
              S11= DPHIX(I)*DPHIX(J)*CNST
              S22= DPHIY(I)*DPHIY(J)*CNST
              S12= DPHIX(I)*DPHIY(J)*CNST 
              S21= DPHIY(I)*DPHIX(J)*CNST 
              
              ESM1(IW,JW)    = ESM1(IW,JW)    + D(NLE,1,1)*S11 + 
     *             D(NLE,3,3)*S22
              ESM1(IW,JW+1)  = ESM1(IW,JW+1)  + D(NLE,1,2)*S12 + 
     *             D(NLE,3,3)*S21
              ESM1(IW+1,JW)  = ESM1(IW+1,JW)  + D(NLE,1,2)*S21 + 
     *             D(NLE,3,3)*S12
              ESM1(IW+1,JW+1)= ESM1(IW+1,JW+1)+ D(NLE,3,3)*S11 + 
     *             D(NLE,2,2)*S22
              
              JW=2*J+1
            ENDDO  
            FT(IW)= ((DPHIX(I)*D(NLE,1,1) + DPHIY(I)*D(NLE,3,1))*ET(1)+
     *           (DPHIX(I)*D(NLE,1,2) + DPHIY(I)*D(NLE,3,2))*ET(2))*CNST
                       
            FT(IW+1)=((DPHIY(I)*D(NLE,2,1) + DPHIX(I)*D(NLE,3,1))*ET(1)+
     *           (DPHIY(I)*D(NLE,2,2) + DPHIX(I)*D(NLE,3,2))*ET(2))*CNST
            
            IW=2*I+1
          ENDDO    
        ENDDO
      ENDDO       

C    
C    %SALIDA DE LOS ELEMENTOS DE LA MATRIZ RIGIDEZ
C      WRITE(10,100)  NLE,D(NLE,1,1),D(NLE,1,2),D(NLE,3,3)
C  100 FORMAT(//,'MATRIZ D DEL ELEMENTO',I3,/'D11=D22,D12=D21,D33',//,
C     *3E15.5)    

      IF(NLE.EQ.NTESTI) THEN
  
C        WRITE(10,110) NLE
C  110   FORMAT(' MATRIZ DE RIGIDEZ',I3)
      
C        DO I=1,8
C          WRITE(10,120) (ESM1(I,J),J=1,8)
C  120     FORMAT(8E10.2)
C        ENDDO
      
        WRITE(10,130) NLE
  130   FORMAT(' FUERZAS TERMICA: ',I3)
      
        DO I=1,8
          WRITE(10,140) FT(I)
  140     FORMAT(E15.5)
        ENDDO
      
      ENDIF

      RETURN  
      END
C

C    ******************SUBRUTINA MODIFYII*************************
C      INSERTA LAS CONDICIONES DE CONTORNO DEL PROBLEMA ELASTICO
C
      SUBROUTINE MODIFYII(AEL,NDES,DESPLA,NODOD,NPE,JGFE,JGSME,NBWE)
      DIMENSION AEL(1000000),DESPLA(100),NODOD(100)
      
      DO KK =1,NDES
          ID=NODOD(KK)
          BD=DESPLA(KK)
          K=ID-1
          DO J=2,NBWE
            M=ID+J-1
            IF(M.LE.NPE) THEN
               IJ=JGSME +(J-1)*NPE+ID-(J-1)*(J-2)/2
               AEL(JGFE+M)=AEL(JGFE+M)-AEL(IJ)*BD
               AEL(IJ)=0.0
            ENDIF   
            IF(K.GT.0) THEN
               KJ=JGSME +(J-1)*NPE+K-(J-1)*(J-2)/2
               AEL(JGFE+K)=AEL(JGFE+K)-AEL(KJ)*BD
               AEL(KJ)=0.0
               K=K-1
            ENDIF   
          ENDDO
          AEL(JGFE+ID)=AEL(JGSME+ID)*BD
      ENDDO
      
      RETURN
      END

      SUBROUTINE DCMPBD_EL(AEL,JGF,JGSM,NP,NBW)
      DIMENSION AEL(1000000)

      NP1 = NP-1
      DO226I=1,NP1
      MJ = I+NBW-1
      IF(MJ.GT.NP) MJ=NP
      NJ = I+1
      MK = NBW
      IF((NP-I+1).LT.NBW) MK=NP-I+1
      ND=0
      DO 225 J=NJ,MJ
      MK = MK-1
      ND = ND+1
      NL = ND+1
      DO 225 K=1,MK
      NK = ND+K
      JK = JGSM + (K-1)*NP + J - ((K-1)*(K-2)/2.)
      INL= JGSM + (NL-1)*NP + I - ((NL-1)*(NL-2)/2.)
      INK= JGSM + (NK-1)*NP + I - ((NK-1)*(NK-2)/2.)
      II = JGSM + I
  225 AEL(JK)=AEL(JK) - (AEL(INL)*AEL(INK)/AEL(II))
  226 CONTINUE 
      RETURN
      END


      SUBROUTINE SLVBD_EL(AEL,JGF,JGSM,NP,NBW)
      DIMENSION AEL(1000000)

      NP1=NP-1
      DO 250 I=1,NP1
      MJ=I+NBW-1
      IF(MJ.GT.NP) MJ=NP
      NJ=I+1
      L=1
      DO 250 J=NJ,MJ
      L=L+1
      IL=JGSM + (L-1)*NP+I-(L-1)*(L-2)/2 
  250 AEL(JGF+J)=AEL(JGF+J)-AEL(IL)*AEL(JGF+I)/AEL(JGSM+I)

      AEL(NP)=AEL(JGF+NP)/AEL(JGSM+NP)
      DO 252 K=1,NP1
      I=NP-K
      MJ=NBW
      IF((I+NBW-1).GT.NP) MJ=NP-I+1
      SUM=0.0
      DO 251 J=2,MJ
      N=I+J-1
      IJ=JGSM+(J-1)*NP+I-(J-1)*(J-2)/2
  251 SUM=SUM+AEL(IJ)*AEL(N)
  252 AEL(I)=(AEL(JGF+I)-SUM)/AEL(JGSM+I)
      
      RETURN
      END   


      SUBROUTINE CALCULO(NE,XC,YC,ALFA,TELEM,NEL,D,STRE,STRA,
     *MAT,XSOL,TELEANT,NPASO,NETEST,NPE,TEREAL,STRE_EL) 
      
      DIMENSION  XC(500),YC(500),ALFA(500),NEL(500,4)
      DIMENSION D(500,3,3),U(8),STRA(500,4,4),STRE(500,4,4),
     *STRE_EL(500,4) 
      
      DIMENSION PHI(4),DPHIX(4),DPHIY(4),AJACO(2,2),AJACOI(2,2),DXHI(4)
     *,DTHE(4),GAUSSPT(2),XNODE(4),YNODE(4)

      DIMENSION MAT(500),TELEM(1000),XSOL(500),TELEANT(1000),
     *TEREAL(1000),X1(4),Y1(4),ET(4)
     
      DIMENSION STRESS(4,4),STRAS(4,4),STRQUIV(4),STEQUIV(4)
      
      DATA GAUSSPT/ -0.57735027, 0.57735027 /
      DATA XNODE / -1, 1, 1, -1 /
      DATA YNODE / -1, -1, 1, 1 /
      
      NN=NPE/2
      
      DO I=1,NN
        XC(I)=XC(I)+XSOL(2*I-1)
        YC(I)=YC(I)+XSOL(2*I)
C        WRITE(10,20) I,XSOL(I*2-1),XSOL(I*2),XC(I),YC(I)
C 20     FORMAT(I3,3X,E15.5,3X,E15.5,3X,E15.5,3X,E15.5)
      ENDDO        
      
      
C
C     CALCULO DE LOS STRESS Y STRAIN DE LOS ELEMENTOS Y VALORES 
C     PRINCIPALES DEL STRESS
C
      DO KK=1,NE
        
        DO I=1,4
          DO J=1,4
            STRESS(I,J)=0.0
            STRAS(I,J)=0.0
          ENDDO 
        ENDDO    
        
        AL = ALFA(KK)
        TEMP= TEREAL(KK) + TELEM(KK)
        PR = FUNPOS(MAT(KK),TEMP)  
        
        T = TELEM(KK)-TELEANT(KK)
        
        ET(1)=AL*T*0.25
        ET(2)=AL*T*0.25
        
        DO I=1,4
          J=NEL(KK,I)
          X1(I)=XC(J)
          Y1(I)=YC(J)
          U(2*I-1)= XSOL(2*J-1)
          U(2*I)  = XSOL(2*J)
        ENDDO
      

        JJK=0
        DO K1=1,2
          DO K2=1,2
            JJK=JJK+1
            XHI = GAUSSPT(K1)
            THE = GAUSSPT(K2)
          
            DO II=1,4
               XHI0 = 1.0 + XHI*XNODE(II)
               THE0 = 1.0 + THE*YNODE(II)
               PHI(II) = 0.25 * XHI0*THE0
               DXHI(II)= 0.25 * XNODE(II)*THE0
               DTHE(II)= 0.25 * YNODE(II)*XHI0
            ENDDO

            DO J=1,2  
              DO I=1,2
                AJACO(I,J)=0.0
              ENDDO
            ENDDO    
            
            DO K=1,4
               AJACO(1,1)=AJACO(1,1)+DXHI(K)*X1(K) 
               AJACO(1,2)=AJACO(1,2)+DXHI(K)*Y1(K)
               AJACO(2,1)=AJACO(2,1)+DTHE(K)*X1(K) 
               AJACO(2,2)=AJACO(2,2)+DTHE(K)*Y1(K)
            ENDDO
          
            DETER= AJACO(1,1)*AJACO(2,2)-AJACO(1,2)*AJACO(2,1)

            IF(DETER.EQ.0.0) THEN
                WRITE(6,*) ' ELEMENTO: ',KK,' PUNTO ',K1,K2,' TIENE 
     *                          DETERMINANTE CERO'
                STOP ' '
            ENDIF

            AJACOI(1,1)=AJACO(2,2)/DETER
            AJACOI(2,2)=AJACO(1,1)/DETER
            AJACOI(1,2)=-AJACO(1,2)/DETER
            AJACOI(2,1)=-AJACO(2,1)/DETER
          
            DO I=1,4
              DPHIX(I) = AJACOI(1,1)*DXHI(I) + AJACOI(1,2)*DTHE(I)
              DPHIY(I) = AJACOI(2,1)*DXHI(I) + AJACOI(2,2)*DTHE(I)
            ENDDO
          
            CX=0.0
            CY=0.0
            DO I=1,4
              CX=CX+X1(I)*PHI(I)
              CY=CY+Y1(I)*PHI(I)
            ENDDO  

            UX=0.0
            UY=0.0
            VX=0.0
            VY=0.0
            DO K3=1,4
                UX = UX + U(2*K3-1) * DPHIX(K3)
                UY = UY + U(2*K3-1) * DPHIY(K3)
                VX = VX + U(2*K3)   * DPHIX(K3)
                VY = VY + U(2*K3)   * DPHIY(K3)
            ENDDO
            
            STRAS(JJK,1) =  UX
            STRAS(JJK,2) =  VY
            STRAS(JJK,3) =  UY + VX

            STRESS(JJK,1)=(D(KK,1,1)*(UX-ET(1)) + D(KK,1,2)*(VY-ET(2)))
            STRESS(JJK,2)=(D(KK,1,2)*(UX-ET(1)) + D(KK,2,2)*(VY-ET(2)))
            STRESS(JJK,3)=D(KK,3,3)*(UY+VX)
            STRESS(JJK,4)=PR*(STRESS(JJK,1)+STRESS(JJK,2))

            IF(KK.EQ.NETEST) THEN
               WRITE(10,95) UX,VY,UX+VY
   95          FORMAT(//,'EXX=',E15.5,3X,'EYY=',E15.5,3X,'EXY=',E15.5)

               WRITE(10,*) ' ET:  ', ET(1),ET(2),AL,T
               WRITE(10,*) ' U-ET:  ', UX-ET(1),VY-ET(2)
                  
               WRITE(10,97) STRESS(JJK,1),STRESS(JJK,2),STRESS(JJK,3),
     *         STRESS(JJK,4)
   97          FORMAT(//,'SXX=',E15.5,3X,'SYY=',E15.5,3X,'SXY=',E15.5,
     *         'SZZ=',E15.5)  
            ENDIF
          ENDDO
        ENDDO       

        DO I=1,4
          DO J=1,4
            STRE(KK,I,J)=STRE(KK,I,J)+STRESS(I,J)
            STRA(KK,I,J)=STRA(KK,I,J)+STRAS(I,J)      
          ENDDO
        ENDDO
        
        do j=1,4 
          stre_el(kk,1) = stre_el(kk,1) + stre(kk,j,1)*0.25
          stre_el(kk,2) = stre_el(kk,2) + stre(kk,j,2)*0.25
          stre_el(kk,3) = stre_el(kk,3) + stre(kk,j,3)*0.25
          stre_el(kk,4) = stre_el(kk,4) + stre(kk,j,4)*0.25
        enddo  
        
        IF(KK.EQ.NETEST) THEN    

           WRITE(10,*) KK,NPASO
           DO J=1,4
             STEQUIV(J) = SQRT((STRA(KK,J,1)-STRA(KK,J,2))**2 + 
     *    STRA(KK,J,1)**2 + STRA(KK,J,2)**2 + 6* STRA(KK,J,3)**2 )
             STRQUIV(J) = SQRT((STRE(KK,J,1)-STRE(KK,J,2))**2 + 
     * (STRE(KK,J,2)-STRE(KK,J,4))**2 + (STRE(KK,J,4)-STRE(KK,J,1))**2 )
             WRITE(10,*) STRQUIV(J),STEQUIV(J)
           ENDDO
        
        ENDIF

        IF(MOD(NPASO,10).EQ.0) THEN
         
          WRITE(12,*) KK
          DO J=1,4
            WRITE(12,*) J,(STRE(KK,J,I),I=1,4)
          ENDDO

          WRITE(13,*) KK
          DO J=1,4
            WRITE(13,*) J,(STRA(KK,J,I),I=1,4)
          ENDDO
        
        ENDIF    

           
      ENDDO
        
      RETURN
      END
      

      SUBROUTINE SALIDAELASTICA(NPE,NE,XC,YC,NEL,STRE,NPASO)
      DIMENSION XC(500),YC(500),NEL(500,4),STRE(500,4)
    
      
  100 FORMAT('*DIMENSION'/'2')
  110 FORMAT('*COORDINATES'/I6)
  120 FORMAT( I6,1X,E12.6,1X,E12.6)
  130 FORMAT('*ELEMENT GROUPS'/'1'/'1',I6,'  Tri3'/'*INCIDENCES'/)
  140 FORMAT(3I6)
  150 FORMAT(/'*SCALAR FIELDS')
      
      WRITE(9,*) NPASO
      
      NN=NPE/2
      WRITE(9,*) NN
      WRITE(9,120) (I,XC(I),YC(I),I=1,NN)
      WRITE(9,*) NE
      DO KK=1,NE
        WRITE(9,*) (NEL(KK,I),I=1,4)
      ENDDO  
      DO I=1,NE 
        WRITE(9,*) STRE(I,1),STRE(I,2),STRE(I,3),STRE(I,4)
      END DO

      RETURN
      END




C
C     FUNCIONES LANDA = CALOR ESPECIFICO*DENSIDAD
      FUNCTION FUNLANDA(MAT,TEMP)
      
c      circonio
c      Densidad, g/cm3: 3.8
c      Calor específico (cal/g/°C): 0.066 *4.18 =  0.27588  J/g/k
      
c      U-MO
c      Densidad 17.23 g/cm3 
      
      
c      UO2
c      Densidad   10.97 g/cm3
c      Calor especifico
      A1=0.30227  !J/G/K
      A2=8.463E-6 !J/G/K^2
      A3=8.741E+4 !J/G
      TITA = 548.68 !K
      EA = 18531.7  !K
      
      TE=TEMP + 273

      IF(MAT.EQ.1)  FUNLANDA = (A1*(TITA**2)*EXP(TITA/TE)/(TE**2*
     *(EXP(TITA/TE)-1)) + 2*A2*TE+ A3*EA*EXP(-EA/TE)/TE**2)*17.23

C      0.388*17.23  ! 6.69J/cm3/K (USE EL CP DEL UO2)??
      
      IF(MAT.EQ.2)  FUNLANDA = 1.048344     ! J/cm3/K 
                          
      RETURN
      END

      FUNCTION FUNDX(MAT,TEMP)
      TE=TEMP+273
      A=3.4944
      B=2.243E-2
      C=6.157E+7
      W=1.41
      AK=8.6142E-5          
      IF(MAT.EQ.1) FUNDX=(1/(A+B*TE)) + (C/TE**2)*EXP(-W/(AK*TE))                               

      IF(MAT.EQ.2) FUNDX=0.16487
      RETURN
      END
      
      FUNCTION FUNDY(MAT,TEMP)
      TE=TEMP+273
      A=3.4944
      B=2.243E-2
      C=6.157E+7
      W=1.41
      AK=8.6142E-5          
      IF(MAT.EQ.1) FUNDY=(1/(A+B*TE)) + (C/TE**2)*EXP(-W/(AK*TE))                               
      IF(MAT.EQ.2) FUNDX=0.16487
      RETURN
      END
      
      FUNCTION FUNGE(MAT)
      IF(MAT.EQ.1) FUNGE=0.0
      IF(MAT.EQ.2) FUNGE=0.0
      IF(MAT.EQ.3) FUNGE=0.0
      IF(MAT.EQ.4) FUNGE=0.0
      RETURN
      END
      
      FUNCTION FUNQE(MAT)
      IF(MAT.EQ.1) FUNQE=0.0
      IF(MAT.EQ.2) FUNQE=0.0
      IF(MAT.EQ.3) FUNQE=0.0
      IF(MAT.EQ.4) FUNQE=0.0
      RETURN
      END                                   
      
      
C    COEFICIENTE DE DILATACION TERMICA  
      FUNCTION FUNAL(MAT,T)

C UMO      
      IF(MAT.EQ.1) FUNAL=  7.91E-6 + 1.21E-8*(T+273)
      
C PARA EL UO2
C      (-4.972D-4+((T+273)*7.107D-6)+(((T+273)**2)*2.583D-9))/(T+273)

      IF(MAT.EQ.2) FUNAL= 6.72E-6 - (2.07E-3)/(T+273)    
      RETURN
      END

C    MODULO DE YOUNG  
      FUNCTION FUNEM(MAT,T)
c uo2
c      D=0.95
c      IF(MAT.EQ.1) FUNEM=2.26E+7*(1-(1.131E-4)*T)*(1- 2.62*(1-D))
c uRANIO METALICO
      IF(MAT.EQ.1) FUNEM=17582250
      
      IF(MAT.EQ.2) FUNEM=1.236E+7-6.221E+3*(T+273)
      RETURN
      END      

C    MODULO DE POISSON
      FUNCTION FUNPOS(MAT,T)
c  uo2
c      IF(MAT.EQ.1) FUNPOS=0.316
c  uranio metalico
      IF(MAT.EQ.1) FUNPOS=0.2
      IF(MAT.EQ.2) FUNPOS=0.325
      
      RETURN
      END 
      
      
      FUNCTION FUNYIELD(MAT,T)
      
      IF(MAT.EQ.1) FUNYIELD=22754
      IF(MAT.EQ.2) FUNYIELD=6.578E+4 * (1-1.686E-3*T+7.748E-7*T**2) 
      
      RETURN
      END
      
     
      
      
                                                                                                                        