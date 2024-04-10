      PROGRAM TDFIELD 
C     PROGRAMA TDFIELD RESUELVE EL PROBLEMA GOBERNADO POR LA ECUACION
C                
C            -d(b.du/dx)/dx - d(b.du/dy)/dy + G.u = Q 
C     CON CONDICIONES DE CONTORNO DE DIRICHLET Y NEUMANN 
      COMMON/ELMATX/ESM(3,3),EST(3,3),EF(3),X(3),Y(3),ESU(3),ESV(3)
      COMMON/HCV/IDBC(50,2),DBC(50,2),NDBC
      COMMON/MATL/DXE,DYE,GE,QE,LANDA
      COMMON/TLE/TITLE(20)
      COMMON/AV/A(200000),JGF,JGSM,NP,NBW,NN 
      DIMENSION NEL(500,3),XC(500),YC(500),MAT(500),TE(500),TELEM(500)
      DIMENSION NS(3),PHI(3),MASA(200000),DES(1000),VEL(1000)
      DIMENSION B(3),C(3),XEIG(1000)
C
C   
      CALL INICIO
C      
      CALL LECTURA(NP,NE,XC,YC,NEL,MAT,NDBC,IDBC,DBC)
C
C
C     CREACION E INICIALIZACION DE UN VECTOR
C
C
      INBW=0
      NBW =0
      DO 20 KK=1,NE
      DO 25 I=1,3
   25 NS(I)=NEL(KK,I)
      DO 21 I=1,2
      IJ=I+1
      DO21 J=IJ,3
      NB=IABS(NS(I)-NS(J))
      IF(NB.EQ.0) WRITE(5,26) KK
   26 FORMAT(/,'ELEMENTO',I3,'TIENE DOS NODOS CON EL MISMO NUMERO')
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
   22 WRITE(6,23)
   23 FORMAT(//,'DIMENSION DEL VECTOR EXCEDIDA,EJECUCION TERMINADA')
      STOP
      ENDIF
      
      DO 90 K=1,NP
   90 DES(I)=0.0
     
C     COMIENZO DEL LOOP TEMPORAL
      CALL INITIME(NTEM,TIME,NT,DT,ALF,EPSLN,NLIMITE)
      IF(NTEM.NE.0)  THEN
      NPASO=0
      WRITE(8,1199) NLIMITE,NP,DT
 1199 FORMAT(I3,1X,I3,1X,E15.5)    
      ENDIF
      
      K=1
 1032 NPASO=NPASO+1
      
      IF(NPASO .EQ. K*10) THEN
            
      WRITE(6,*) ' PASO TEMPORAL:        ',NPASO
      K=K+1
      ENDIF
      
      A2 = (1-ALF)*DT
      A1 = ALF*DT

      DO 24 I=1,JEND
      MASA(I)=0.0
   24 A(I)=0.0
C
C     GENERACION DEL SISTEMA DE ECUACIONES
      DO32 KK=1,NE

      DO 31 I=1,3
      NS(I)=NEL(KK,I)
      J=NS(I)
      ESU(I)=DES(J)
C      ESV(I)=VEL(I)
      X(I)=XC(J)
  31  Y(I)=YC(J)
      
C     COEFICIENTES DE LOS ELEMENTOS
C     NMTL ES ALGO ASI COMO NUMERO DE MATERIAL. EN CADA UNO CAMBIA LA 
C     CONDICION DE CONTORNO
      
      DXE=FUNDX(MAT(KK))
      DYE=FUNDY(MAT(KK))
      GE=FUNGE(MAT(KK))
      QE=FUNQE(MAT(KK))
      LANDA=FUNLANDA()
      
c      IF(NPASO.EQ.1) 
      CALL ELSTMF(KK)    
      
      CALL TIMER(NTEM,A1,A2,KK)
C
C    PROCEDIMIENTO DE RIGIDEZ DIRECTO
      DO 33 I=1,3
      II=NS(I)
      A(JGF+II)=A(JGF+II)+EF(I)
      DO34J=1,3
      JJ=NS(J)+1-II
      IF(JJ.LE.0) GOTO34
      J1=JGSM+(JJ-1)*NP+II-(JJ-1)*(JJ-2)/2
      A(J1)=A(J1)+ESM(I,J)
      MASA(J1-JGSM)=MASA(J1-JGSM) + EST(I,J)
   34 CONTINUE
   33 CONTINUE
   32 CONTINUE   
C
C     ****************************************************
C     * MODIFICACION Y SOLUCION DEL SISTEMA DE ECUACIONES*
C     * LOS DATOS VIENEN DE LA SUBRUTINA MODIFY          *
C     ****************************************************
C    
      CALL MODIFY
C                          
C      WRITE(6,*) ' VIENE DE MODIFY Y LLAMA A DCMPBD'
C
      CALL DCMPBD
C     
C      WRITE(6,*) ' VIENE DE DCMPBD Y LLAMA A SLVBD'
C
      CALL SLVBD
C
C      WRITE(6,*) ' SALE DE SLVBD'
C
C     SALIDA DE LAS SOLUCIONES 
C      
      WRITE(1,112)
  112 FORMAT(' VALORES DE LAS TEMPERATURAS NODALES'/)
      DO 113 I=1,NP
      TE(I)=A(I)
      WRITE(8,150) NPASO,I,TE(I)
  150 FORMAT(2I3,2X,2E15.5)   
  113 WRITE(1,111) NPASO,I,TE(I)
  111 FORMAT(2I3,3X,2E15.5)
      
C      WRITE(1,116)
C  116 FORMAT(' VALORES DE LAS TEMPERATURAS ELEMENTALES'/)
C      
C      DO 114J=1,NE
C      TELEM(J)=(TE(NEL(J,1))+TE(NEL(J,2))+TE(NEL(J,3)))/3
C  114 WRITE(1,115) J,TELEM(J)
C  115 FORMAT(I3,3X,E15.5,3X,E15.5)

      NT   = NT   + 1
      TIME = TIME + DT
C
C     EVALUACION DEL VOLUMEN BAJO LA TEMPERATURA DE SUPERFICIE Y EL 
C     GRADIENTE ELEMENTAL
C 
C      ILINE=0
C      DO83 KK=1,NE
C      IF(ILINE.GT.0) GOTO 110
C
C     OUTPUT OF TEH CORRECT GRADIENT HEADING
C      WRITE(1,43) 
C   43 FORMAT(//,'ELEMENT:',4X,'LOCATION','Q(X)',2X,'Q(Y)')
C
C  110 ILINE=ILINE+2
C      IF(ILINE.GT.50) ILINE=0
C      
C     RETRIEVAL OF THE NODAL COORDINATES, THE NODE NUMBER, AND THE TEMP.
C      SP=0.0
C      DO40I=1,3
C      NS(I)=NEL(KK,I)
C      J=NS(I)
C      X(I)=XC(J)
C      Y(I)=YC(J)
C      PHI(I)=A(J)
C   40 SP=SP+PHI(I)
C   
C     DXE=FUNDX(MAT(KK))
C      DYE=FUNDY(MAT(KK))
      
C     EVALUACION DEL GRADIENTE EN LOS ELEMENTOS
C      B(1)=Y(2)-Y(3)
C      B(2)=Y(3)-Y(1)
C      B(3)=Y(1)-Y(2)
C      C(1)=X(3)-X(2)
C      C(2)=X(1)-X(3)
C      C(3)=X(2)-X(1)
C      AR2=X(2)*Y(3)+X(3)*Y(1)+X(1)*Y(2)-X(2)*Y(1)-X(3)*Y(2)-X(1)*Y(3)
C      GRADX=(B(1)*PHI(1)+B(2)*PHI(2)+B(3)*PHI(3))/AR2
C      GRADY=(C(1)*PHI(1)+C(2)*PHI(2)+C(3)*PHI(3))/AR2
C      GRADX=DXE*GRADX*(-1.)
C      GRADY=DYE*GRADY*(-1.)   
C      
C      WRITE(1,52) KK,GRADX,GRADY
C   52 FORMAT(/,I3,'CENTRE',3X,2E15.5)  
      
   
C      VOL=VOL+SP*AR2/6
C   83 CONTINUE           
      
      
      DIFF=0.0
      SOLN=0.0
      DO I=1,NP
         SOLN=SOLN + A(I)*A(I)
         DIFF=DIFF + (A(I)-DES(I)) * (A(I)-DES(I))
         DES(I)=A(I)
      ENDDO   
      
      PERCNT=1.0
      IF(SOLN.NE.0.0) PERCNT = DSQRT(DIFF/SOLN)
      
      IF(PERCNT.LE.EPSLN) THEN
          WRITE(6,*) ' LA SOLUCION ALCANZO UN ESTADO ESTACIONARIO'
          WRITE(1,300) NPASO,TIME
 300      FORMAT(' ESTADO ESTACIONARIO PASO :',I3,'  TIEMPO:',E12.5)       
          CALL SALIDA(NP,NE,XC,YC,NEL,TE)
          STOP ' '
      ENDIF
      
      
      IF(NPASO.GE.NLIMITE) THEN
          CALL SALIDA(NP,NE,XC,YC,NEL,TE)
          STOP
      ELSE      
          GOTO 1032
      ENDIF
   
C      CALL SALIDA(NP,NE,XC,YC,NEL,TE)

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
      OPEN(UNIT=2,FILE='prutem.fem',STATUS='OLD',ERR=3001)
      
      OPEN (UNIT=3,FILE='prutem.txt',STATUS='OLD',ERR=3002)
      OPEN(UNIT=4, FILE='Temdat.TXT',STATUS='UNKNOWN',ERR=3004)
      OPEN(UNIT=8, FILE='TEMTIME.NDP',STATUS='UNKNOWN',ERR=3005)

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
      DIMENSION XC(500),YC(500),NEL(500,3),IDBC(50,2),DBC(50,2)
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
      OPEN(1,FILE='temper.txt',STATUS='UNKNOWN')
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
  100 READ(2,*,ERR=3023,END=3023) N,(NEL(KK,I),I=1,3)
      WRITE (1,14) NE
   14 FORMAT(' NUMERO DE ELEMENTOS NE =',I5)
      WRITE(1,110)
  110 FORMAT(1X,'ELEMENTO',3X,'MATERIAL',3X,'NODOS')
      DO101 J=1,NE
  101 WRITE(1,120) J,MAT(J),(NEL(J,I),I=1,3)
  120 FORMAT(I3,3X,I3,3X,3I4)
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
     
      OPEN(7,FILE='SALIDA.NDP',STATUS='UNKNOWN')
      WRITE(7,100)
      WRITE(7,110) NP
      WRITE(7,120) (I,XC(I),YC(I),I=1,NP)
      WRITE(7,130) NE
      DO 1000 KK=1,NE
 1000 WRITE(7,140) (NEL(KK,I),I=1,3)
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
  140 FORMAT(3I6)
  150 FORMAT(/'*SCALAR FIELDS')
      END
C
C
C    %SUBRUTINA ELSTMF, CALCULA LAS MATRICES DE CADA ELEMENTO
      SUBROUTINE ELSTMF(KK)
      COMMON/ELMATX/ESM(3,3),EST(3,3),EF(3),X(3),Y(3),ESU(3),ESV(3)
      COMMON/HCV/IDBC(50,2),DBC(50,2),NDBC
      COMMON/MATL/DXE,DYE,GE,QE,LANDA
C      DIMENSION ES(4,4),ET(4,4),EG(4,4)
      DIMENSION B(3),C(3)
C      REAL LG
      
C      DATA ES/2.,-2.,-1.,1.,-2.,2.,1.,-1.,-1.,2.,-2.,1.,-1.,-2.,2./
C      DATA ET/2.,1.,-1.,-2.,1.,2.,-2.,-1.,-1.,-2.,2.,1.,-2.,-1.,1.,2./
C      DATA EG/4.,2.,1.,2.,2.,4.,2.,1.,1.,2.,4.,2.,2.,1.,2.,4./
      WRITE(4,100) KK,DXE,DYE,GE,QE
  100 FORMAT(/,I3,4E15.5) 
    
      B(1)=Y(2)-Y(3)
      B(2)=Y(3)-Y(1)
      B(3)=Y(1)-Y(2)
      C(1)=X(3)-X(2)
      C(2)=X(1)-X(3)
      C(3)=X(2)-X(1)
      AR2=X(2)*Y(3)+X(3)*Y(1)+X(1)*Y(2)-X(2)*Y(1)-X(3)*Y(2)-X(1)*Y(3)
      IF(ABS(AR2).LT.0.0001) GOTO5
      DO1 I=1,3
      EF(I)=QE*AR2/6.
      DO1 J=1,3
      A=1.0
      IF(I.EQ.J) A=2.0
   1  ESM(I,J)=((DXE*B(I)*B(J)+DYE*C(I)*C(J))/(AR2*2.))+A*GE*AR2/24.
      
C     S*Lij/2 * [1,1,0]tranapuesto se suma a la fuerza
C     M*Lij/6 * [2,1,0][1,2,0][0,0,0] SE SUMA A MATRIZ RIGIDEZ            
      IF(NDBC.EQ.0) GOTO 7
       
      DO11 I=1,NDBC
                  
      IF(IDBC(I,1).NE.KK) GOTO11
      J=IDBC(I,2)
      K=J+1
      IF(J.EQ.3) K=1
      EF(J)=EF(J)+DBC(I,2)/2.
      EF(K)=EF(K)+DBC(I,2)/2.
      ESM(J,J)=ESM(J,J)+DBC(I,1)/3.
      ESM(J,K)=ESM(J,K)+DBC(I,1)/6.
      ESM(K,J)=ESM(J,K)
      ESM(K,K)=ESM(K,K)+DBC(I,1)/3.
   11 CONTINUE
   
C     MATRIZ DE MASA
   7  DO50 K=1,3    
      DO50 J=1,3
      A=2.0
      IF(K.EQ.J) A=1.0
   50 EST(K,J)=LANDA*AR2/(12*A)
     
c      WRITE(4,8) KK
c    8 FORMAT(2X,'ELEMENTO:',I4,/,2X,'VECTOR FUERZA',2X,'MATRIZ RIGIDEZ')
c      DO 9 I=1,3
c    9 WRITE(4,10) EF(I),(ESM(I,J),J=1,3)
c   10 FORMAT(2X,E12.5,2X,3E15.5)
c      WRITE(4,17) KK
c   17 FORMAT('MATRIZ DE MASA',I3)  
c     DO19I=1,3
c  19 WRITE(4,21) (EST(I,J),J=1,3)
c 21  FORMAT(3E15.5)
c
      RETURN
      
    5 WRITE(1,6) KK
    6 FORMAT(//2X,'AREA DEL ELEMENTO',I3,'ES MENOR QUE 0.0001')
      STOP
      END
C
C     SUBRUTINA TIMER
      SUBROUTINE TIMER(NTEM,A1,A2,NLE)
      COMMON/ELMATX/ESM(3,3),EST(3,3),EF(3),X(3),Y(3),ESU(3),ESV(3)
     
      IF(NTEM.NE.0) THEN
               
        DO 30 K=1,3
        SUM=0.0
        DO 20 J=1,3
        SUM   = SUM + (EST(K,J) - A2*ESM(K,J)) * ESU(J)  
  20    ESM(K,J) = EST(K,J) + A1*ESM(K,J)
  30    EF(K) = (A1+A2)*EF(K) + SUM

      
      WRITE(4,8)  NLE
    8 FORMAT(2X,'DATOS TIMER',/,'ELEMENTO:',I4,/,
     *2X,'VECTOR FUERZA',2X,'MATRIZ RIGIDEZ')
      DO 9 I=1,3
    9 WRITE(4,10) EF(I),(ESM(I,J),J=1,3)
   10 FORMAT(2X,E12.5,2X,3E15.5)
      WRITE(4,17) KK
   17 FORMAT('MATRIZ DE MASA',I3)  
      DO19I=1,3
   19 WRITE(4,21) (EST(I,J),J=1,3)
  21  FORMAT(3E15.5)
      
      ENDIF
      RETURN
      END

C
C    SUBRUTINA MODIFY INSERTA LAS CONDICIONES DE CONTORNO
      SUBROUTINE MODIFY
      COMMON/AV/A(200000),JGF,JGSM,NP,NBW,NN
c
C
      WRITE(1,201)
  201 FORMAT(//,'FUENTES AND VALORES SINK')
      
      READ(3,*) NN
      IF(NN.EQ.0) GOTO 217

      DO 202 I=1,NN
      READ(3,*,ERR=3020,END=3020) NOD,TEM
      A(JGF+NOD)=A(JGF+NOD)+TEM
  202 WRITE(1,501) I,NOD,TEM
  501 FORMAT(I6,1X,i6,1x,e15.5)   
  
C
C     %ENTRADA DE VALORES DE NODOS PRESCRIPTOS
C         
C          
C         ID - GRADOS DE LIBERTAD DE LOS DESPLAZAMIENTOS CONOCIDOS 
C         BD - VALORES DE LOS DESPLAZAMIENTOS
C
C
  217 WRITE(1,208)
  208 FORMAT(' VALORES DE LAS TEMPERATURAS CONOCIDAS') 
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
  216 WRITE(1,209) ID,BD
  209 FORMAT(I3,E15.5)
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
      SUBROUTINE INITIME(NTEM,TIME,NT,DT,ALF,EPSLN,NLIMITE)

C     NTEM INDICA SI EL PROBLEMA ES TEMPORAL (NTME =1) O NO (NTIM=0)
      NTEM=1
C     TIME TIEMPO INICIAL
      TIME=0.0           
C     NT NO SE QUE ES 
      NT=0  
C     DT ES EL INTERVALO DE TIEMPO          
      DT=0.01                      
C     ALF COEFFICIENTE ALFA
      ALF=0.5                 
      EPSLN=0.001
C     NLIMITE ES LA CANTIDAD DE PASOS MAXIMA PERMITIDA 
      NLIMITE = 100                                    
      RETURN
      END
C
C     FUNCIONES
      FUNCTION FUNLANDA()
      FUNLANDA = 1
      RETURN
      END

C
      FUNCTION FUNDX(MAT)
      IF(MAT.EQ.1) FUNDX=1
      IF(MAT.EQ.2) FUNDX=0.00336
      IF(MAT.EQ.3) FUNDX=0.22
      IF(MAT.EQ.4) FUNDX=1
      RETURN
      END
      
      FUNCTION FUNDY(MAT)
      IF(MAT.EQ.1) FUNDY=1
      IF(MAT.EQ.2) FUNDY=0.00336
      IF(MAT.EQ.3) FUNDY=0.22
      IF(MAT.EQ.4) FUNDY=1
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
      IF(MAT.EQ.1) FUNQE=1.0
      IF(MAT.EQ.2) FUNQE=0.0
      IF(MAT.EQ.3) FUNQE=0.0
      IF(MAT.EQ.4) FUNQE=0.0
      RETURN
      END
      
                                                                                                                        