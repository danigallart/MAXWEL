      PROGRAM TDtime 
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
      DIMENSION NS(3),AMASA(200000),DES(1000),VEL(1000),XEIG(1000)
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
C     ************ A U T O V A L O R E S **************
C     PROBLEMA DE AUTOVALORES
      DO 724 I=1,JEND
      AMASA(I)=0.0
  724 A(I)=0.0
C
C     GENERACION DEL SISTEMA DE ECUACIONES

      DO732 KK=1,NE
      DO 731 I=1,3
      NS(I)=NEL(KK,I)
      J=NS(I)
      ESU(I)=DES(J)
C      ESV(I)=VEL(I)
      X(I)=XC(J)
 731  Y(I)=YC(J)
      
C     COEFICIENTES DE LOS ELEMENTOS
C     NMTL ES ALGO ASI COMO NUMERO DE MATERIAL. EN CADA UNO CAMBIA LA 
C     CONDICION DE CONTORNO
      
      DXE=FUNDX(MAT(KK))
      DYE=FUNDY(MAT(KK))
      GE=FUNGE(MAT(KK))
      QE=FUNQE(MAT(KK))
      LANDA=FUNLANDA()
      
      CALL ELSTMF(KK,NPASO)    

      DO 733 I=1,3
      II=NS(I)
      A(JGF+II)=A(JGF+II)+EF(I)
      DO734J=1,3
      JJ=NS(J)+1-II
      IF(JJ.LE.0) GOTO734
      J1=JGSM+(JJ-1)*NP+II-(JJ-1)*(JJ-2)/2
      A(J1)=A(J1)+ESM(I,J)
      AMASA(J1-JGSM)=AMASA(J1-JGSM) + EST(I,J)
  734 CONTINUE
  733 CONTINUE
  732 CONTINUE   


      CALL MODIFY(NPASO)   
      
      CALL MODIFYII(AMASA,NBW,NP,JEND)

      CALL AUTOV(A,AMASA,XEIG,NP,JGSM,JEND)

       VMAX=ABS(XEIG(1))
       DO 222 K=1,NP
       WRITE(1,221) K,XEIG(K)
221    FORMAT( 'VALOR PROPIO: ',I3,2X,E15.5)   
       IF(ABS(XEIG(K)).GT.VMAX) VMAX = ABS(XEIG(K))
222    CONTINUE
       TCRITICO=0
       IF(ALF.NE.0.5) THEN
           TCRITICO = 2/(1-2*ALF)/VMAX
           DTS = TCRITICO / 2
       ENDIF
       WRITE(1,223) VMAX,TCRITICO,DTS
223    FORMAT('MAXIMO AUTOVALOR: ',E15.5,/,'TCRITICO :',E15.5,/,
     *'DT SUGERIDO: ',E15.5)

      WRITE(6,*) ' FIN DEL PROBLEMA DE AUTOVALORES'

C     **************TRANSIENT PROBLEM ******************    
C     COMIENZO DEL LOOP TEMPORAL
      CALL INITIME(NTEM,TIME,NT,DT,ALF,EPSLN,NLIMITE)

      IF(NTEM.NE.0)  THEN

      IF(DT.LT.DTS .OR. DT.GE.(2*DTS) .AND. ALF.LT.0.5 )  DT=DTS
      
      NPASO=0
      WRITE(8,1199) NLIMITE,NP,DT
 1199 FORMAT(I3,1X,I3,1X,E15.5)    
      ENDIF
      
      NKPASO=1
 1032 NPASO=NPASO+1
      
      IF(NPASO .EQ. NKPASO*5) THEN
c            
      WRITE(6,*) ' PASO TEMPORAL:        ',NPASO
      NKPASO=NKPASO+1
      ENDIF
      
      A2 = (1-ALF)*DT
      A1 = ALF*DT

      DO 24 I=1,JEND
      AMASA(I)=0.0
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
      
      CALL ELSTMF(KK,NPASO)    
      
      CALL TIMER(KK,NTEM,A1,A2,NPASO)
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
      AMASA(J1-JGSM)=AMASA(J1-JGSM) + EST(I,J)
   34 CONTINUE
   33 CONTINUE
   32 CONTINUE   

C
C     ****************************************************
C     * MODIFICACION Y SOLUCION DEL SISTEMA DE ECUACIONES*
C     * LOS DATOS VIENEN DE LA SUBRUTINA MODIFY          *
C     ****************************************************
C    
      CALL MODIFY(NPASO)
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
      WRITE(1,112) NPASO
  112 FORMAT(' VALORES DE LAS TEMPERATURAS NODALES DEL PASO:',I3,/)
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
      DO 290 I=1,NP
      IF(NTEM.NE.0) THEN
         SOLN=SOLN + A(I)*A(I)
         DIFF=DIFF + (A(I)-DES(I)) * (A(I)-DES(I))
C      WRITE(6,*) 'DIFF  ',DIFF,' SOLN  ',SOLN,DSQRT(DIFF/SOLN)
      ENDIF
 290  DES(I)=A(I)
     
      IF(NTEM.NE.0 .AND. NT.GT.1) THEN
      PERCNT = DSQRT(DIFF/SOLN)
      IF(PERCNT.LE.EPSLN.AND.NK1.EQ.0) THEN
          WRITE(6,*) ' LA SOLUCION ALCANZO UN ESTADO ESTACIONARIO'
          WRITE(1,300) NPASO,TIME
 300  FORMAT(' ESTADO ESTACIONARIO PASO :',I3,'  TIEMPO:',E12.5)       
      NK1=1
      ENDIF
      ENDIF
      
      
      IF(NTEM.NE.0) THEN
          IF(NPASO.GE.NLIMITE) THEN
                CALL SALIDA(NP,NE,XC,YC,NEL,TE)
                STOP
          ELSE      
               GOTO 1032
          ENDIF
      ENDIF
   
      CALL SALIDA(NP,NE,XC,YC,NEL,TE)
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
      CHARACTER*70 NFEM,NFILE
      
    3 FORMAT(A70)
  100 FORMAT(' ARCHIVO .FEM                   :',$)
  110 FORMAT(' ARCHIVO CON LOS DATOS          :',$)    
C    
      WRITE(6,100)
      READ(5,3) NFEM 
      OPEN(UNIT=2,FILE=NFEM,STATUS='OLD',ERR=3001)
      
      WRITE(6,110)
      READ(5,3) NFILE
      OPEN (UNIT=3,FILE=NFILE,STATUS='OLD',ERR=3002)
C
      OPEN(UNIT=4, FILE='Temdat.TXT',STATUS='UNKNOWN',ERR=3004)
      OPEN(UNIT=8, FILE='TEMTIME.NDP',STATUS='UNKNOWN',ERR=3005)

C              
      RETURN
C      
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
      SUBROUTINE ELSTMF(KK,NPASO)
      COMMON/ELMATX/ESM(3,3),EST(3,3),EF(3),X(3),Y(3),ESU(3),ESV(3)
      COMMON/HCV/IDBC(50,2),DBC(50,2),NDBC
      COMMON/MATL/DXE,DYE,GE,QE,LANDA
C      DIMENSION ES(4,4),ET(4,4),EG(4,4)
      DIMENSION B(3),C(3)
C      REAL LG
      
C      DATA ES/2.,-2.,-1.,1.,-2.,2.,1.,-1.,-1.,2.,-2.,1.,-1.,-2.,2./
C      DATA ET/2.,1.,-1.,-2.,1.,2.,-2.,-1.,-1.,-2.,2.,1.,-2.,-1.,1.,2./
C      DATA EG/4.,2.,1.,2.,2.,4.,2.,1.,1.,2.,4.,2.,2.,1.,2.,4./
C      WRITE(4,100) KK,DXE,DYE,GE,QE
C  100 FORMAT(/,I3,4E15.5) 
    
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
      
      IF(NPASO.EQ.1) THEN
      WRITE(4,8) KK
    8 FORMAT(2X,'ELEMENTO:',I4,/,2X,'VECTOR FUERZA',2X,'MATRIZ RIGIDEZ')
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
      
    5 WRITE(1,6) KK
    6 FORMAT(//2X,'AREA DEL ELEMENTO',I3,'ES MENOR QUE 0.0001')
      STOP
      END
C
C     SUBRUTINA TIMER
      SUBROUTINE TIMER(KK,NTEM,A1,A2,NPASO)
      COMMON/ELMATX/ESM(3,3),EST(3,3),EF(3),X(3),Y(3),ESU(3),ESV(3)
     
      IF(NTEM.NE.0) THEN
               
        DO 30 K=1,3
        SUM=0.0
        DO 20 J=1,3
        SUM   = SUM + (EST(K,J) - A2*ESM(K,J)) * ESU(J)  
  20    ESM(K,J) = EST(K,J) + A1*ESM(K,J)
  30    EF(K) = (A1+A2)*EF(K) + SUM

      IF(NPASO.EQ.1) THEN
      
      WRITE(4,8) KK
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
      
      ENDIF
      RETURN
      END

C
C    SUBRUTINA MODIFY INSERTA LAS CONDICIONES DE CONTORNO
      SUBROUTINE MODIFY(NPASO)
      COMMON/AV/A(200000),JGF,JGSM,NP,NBW,NN
c
C
      IF(NPASO.EQ.1) THEN
      WRITE(1,201)
  201 FORMAT(//,'FUENTES AND VALORES SINK')
      ENDIF
      
      READ(3,*) NN
      
      IF(NN.NE.0) THEN

      DO 202 I=1,NN
      READ(3,*,ERR=3020,END=3020) NOD,TEM
      A(JGF+NOD)=A(JGF+NOD)+TEM

      IF(NPASO.EQ.1) THEN
      WRITE(1,501) I,NOD,TEM
  501 FORMAT(I6,1X,i6,1x,e15.5)   
      ENDIF

  202 CONTINUE
      ENDIF
C
C     %ENTRADA DE VALORES DE NODOS PRESCRIPTOS
C         
C          
C         ID - GRADOS DE LIBERTAD DE LOS DESPLAZAMIENTOS CONOCIDOS 
C         BD - VALORES DE LOS DESPLAZAMIENTOS
C
C

      IF(NPASO.EQ.1) THEN
      WRITE(1,208)
  208 FORMAT(' VALORES DE LAS TEMPERATURAS CONOCIDAS') 
      ENDIF
      
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
      IF(NPASO.EQ.1) THEN
      WRITE(1,209) ID,BD
  209 FORMAT(I3,E15.5)
      ENDIF
  216 CONTINUE
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
C    SUBRUTINA MODIFYII INSERTA LAS CONDICIONES DE CONTORNO EN LA MATRIZ 
C    DE MASA PARA EL PROBLEMA DE AUTOVALORES
      SUBROUTINE MODIFYII(A,NBW,NP,JEND)
      DIMENSION A(200000)
C
C      
      READ(3,*) NN
      DO 202 I=1,NN
  202 READ(3,*,ERR=3020,END=3020) NOD,TEM
      
      READ(3,*) NDES
      DO 216 KK =1,NDES
      READ(3,*,ERR=3040,END=3040) ID,BD
      K=ID-1
      DO211 J=2,NBW
      M=ID+J-1
      IF(M.GT.NP) GOTO210
      IJ=(J-1)*NP+ID-(J-1)*(J-2)/2
      A(IJ)=0.0
  210 IF(K.LE.0) GOTO 211
      KJ=(J-1)*NP+K-(J-1)*(J-2)/2
      A(KJ)=0.0
      K=K-1
  211 CONTINUE
  221 CONTINUE   
  216 CONTINUE
C
C      DO20 K=1,JEND
C  20  WRITE(1,10) A(K)
C  10  FORMAT('MASA MODIFICADA: ', E15.5)  

      REWIND(3)
      RETURN
C     
 3020 WRITE(6,2020)
      STOP ' '
 3040 WRITE(6,2040)
      STOP' '
C
 2020 FORMAT(' ERROR LEYENDO C DE C. DE FUERZA EN MODIFYII')
 2040 FORMAT(' ERROR LEYENDO C. DE C. DE DESPLAZAMIENTOS EN MODIFYII')
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
      NLIMITE = 50                                    
      RETURN
      END
C
C     SUBRUTINA AUTOV
      SUBROUTINE AUTOV(A,AMASA,XEIG,NP,JGSM,JEND)
      DIMENSION AM(50,50),AVN(50,50),AC(50,50),AUX(50,50),XEIG(50)
      DIMENSION A(200000),AMASA(200000)
      
      OPEN(UNIT=11,FILE='EIGEN.TXT')
      DIM=NP
      
C      DO210K=1,JEND
C 210  WRITE(11,211) 'P',A(K),AMASA(K)
C 211  FORMAT(A2,2E15.5)
C      WRITE(11,65)
      
      NT=NP*(NP+1)/2.
      
      DO200 K=1,NP
      DO200 J=K,NP
      NPOS =JGSM + NT - (NP-(J-K))*(NP-(J-K)+1)/2. + K
      AVN(K,J) = A(NPOS)
      AVN(J,K) = AVN(K,J)
      AM(K,J) = AMASA(NPOS-JGSM)
 200  AM(J,K) = AM(K,J)
      
C      DO201 K=1,NP
C      DO201 J=1,NP
C 201  WRITE(11,202) AVN(K,J),AM(K,J)
C 202  FORMAT(2E15.5)   
      
C      WRITE(11,24) NP,JEND
C  24  FORMAT('ESTO ES LO QUE LE MANDO A MATINV',2I3) 
  
C      DO51I=1,NP  
C      DO51J=1,NP
C   51  WRITE(11,61) I,J,AM(I,J)
C   61  FORMAT(2I3,2X,E15.5)
      
C      WRITE(11,65)
C   65 FORMAT(//)
     
      CALL MATINV(AM,NP,NP,IRET,0)   
      
      IF(IRET.EQ.-1) THEN
      WRITE(6,*) ' MATRIZ M SINGULAR'
      STOP
      ENDIF
      
C      DO5I=1,NP
C      DO5J=1,NP
C   5  WRITE(11,6) AVN(I,J)
C   6  FORMAT(E15.5)

C      WRITE(11,65)
      
      DO 66 K=1,NP
      DO 66 J=1,NP
      AUX(K,J)=0
      DO 67 I=1,NP
   67 AUX(K,J) =AUX(K,J) + AM(K,I)*AVN(I,J)
   66 CONTINUE 
   
      DO 77 K=1,NP
      DO 77 I=1,NP
   77 AVN(K,I) =AUX(K,I)
      
C      DO53I=1,NP          
C      DO53J=1,NP
C   53  WRITE(11,63) AVN(I,J)
C  63  FORMAT(E10.4)
C      WRITE(1,65)      

      DO1KI=1,NP
      DO1KJ=1,NP
   1  AC(KI,KJ)=AVN(KI,KJ) 

      BREAK=0.0
      DO100KK=1,10
      
      AMAX=ABS(AVN(1,2))
      I =1
      JJ=2
      
      DO50KI=1,NP-1
      DO50KJ=KI+1,NP
      IF(ABS(AVN(KI,KJ)).GT.AMAX) THEN
            AMAX = ABS(AVN(KI,KJ))
            I =KI
            JJ=KJ
      ENDIF      
  50  CONTINUE    
C      DO10 I=1,DIM-1
C      DO10 JJ=I+1,DIM
      IF(BREAK.EQ.-1.0) GOTO100

      W=AVN(I,JJ)
      E=AVN(I,I)-AVN(JJ,JJ)
      F=SQRT(E**2 + 4*(W**2))
      
      IF( ABS(-E+F).GT.ABS(-E-F) ) THEN
          D=F-E
      ELSE
          D=-F-E
      ENDIF
      
      S=ABS(D)/SQRT( D**2 + 4*(W**2) )
      C=2*S*W/D
      
      GUAR1 = S**2 + C**2
      GUAR2 = S*C*(AVN(I,I)-AVN(JJ,JJ) ) + (S**2) * AVN(JJ,I) - (C**2)*
     * AVN(I,JJ)
      
      WRITE(11,7) W,E,F,D,S,C,GUAR1,GUAR2
   7  FORMAT('W:',E15.5,/,'E:',E15.5,/,'F',E15.5,/,'D',E15.5,/,'S',E15.5
     *,/,'C',E15.5,/,'GUAR1',E15.5,/,'GUAR2',/,E15.5)
      
      
      AC(I,I) = (C**2)*AVN(I,I) + C*S*AVN(I,JJ) +C*S*AVN(JJ,I) + 
     * (S**2)*AVN(JJ,JJ)
      AC(JJ,JJ) = (C**2)*AVN(JJ,JJ) - C*S*AVN(I,JJ) - C*S*AVN(JJ,I) + 
     * (S**2)*AVN(I,I)
      AC(I,JJ) = C*S*(AVN(I,I) - AVN(JJ,JJ) ) + (S**2)*AVN(JJ,I) - 
     * (C**2)*AVN(I,JJ)
      AC(JJ,I) = C*S*(AVN(I,I) - AVN(JJ,JJ) ) + (S**2)*AVN(I,JJ) - 
     * (C**2)*AVN(JJ,I)     
          
      
      
      DO 20 K=1,NP
      IF( I.NE.K .AND. JJ.NE.K ) THEN
          AC(K,I) = C*AVN(K,I) + S*AVN(K,JJ)
          AC(K,JJ)= C*AVN(K,JJ)+ S*AVN(K,I)
          AC(I,K) = C*AVN(I,K) + S*AVN(JJ,K)
          AC(JJ,K)= -C*AVN(JJ,K)+ S*AVN(I,K)
      ENDIF     
   20 CONTINUE

      DO25KI=1,NP
      TRAA=TRAA+AVN(KI,KI)
      TRAC=TRAC+AC(KI,KI)
      DO25KJ=1,NP
      NA=NA+AVN(KI,KJ)**2
      NAC=NAC+AC(KI,KJ)**2
   25 AVN(KI,KJ)=AC(KI,KJ) 

C      DO26 KI=1,NP   
C      DO26 KJ=1,NP
C   26 WRITE(11,27) AVN(KI,KJ)
C   27 FORMAT(E15.5) 
      
      WRITE(11,29) TRAA, TRAC, SQRT(NA), SQRT(NAC) 
   29 FORMAT(//,4E15.5) 
      
      IF(DABS(SQRT(NA)-SQRT(NAC)).LE.0.00001) THEN
      WRITE(11,11) KK
  11  FORMAT (' CONVERGENCIA REQUERIDA ALCANZADA EN EL PASO: ',I3)
      BREAK=-1.0
      ENDIF
      
C   10 CONTINUE
C      DO32I=1,DIM
c   32 WRITE(6,*) A(I,I)
c      STOP
c      ENDIF
  100 CONTINUE 
  
      DO30I=1,NP
      XEIG(I)=AVN(I,I)
   30 WRITE(11,32) XEIG(I)
   32 FORMAT(E15.5)
   
      RETURN
      END
C     
C     SUBRUTINA MATINV

      SUBROUTINE MATINV(A,N,M,IRET,IENT)                                MAT00010
C                                                                       MAT00020
C CALCULA LA MATRIZ INVERSA DE UNA MATRIZ CUADRADA REAL PERMITIENDO REA-MAT00030
C LIZAR LA INVERSION EN DOS ENTRADAS : EN LA PRIMERA SE REALIZA EL PROCEMAT00040
C SO DE ELIMINACION Y EN LA SEGUNDA, LA SUSTITUCION HACIA ATRAS.        MAT00050
C                                                                       MAT00060
C          -DESCRIPCION DE LOS PARAMETROS-                              MAT00070
C                                                                       MAT00080
C  A: MATRIZ DE DIMENSION M*M CONTENIENDO LA MATRIZ A INVERTIR.         MAT00090
C     AL RETORNO, MATRIZ INVERTIDA.                                     MAT00100
C  N: ORDEN DE LA MATRIZ A INVERTIR.(ES UN PARAMETRO REQUERIDO POR      MAT00110
C    LA RUTINA DIFSUB- SI EL USUARIO LO DESEA PUEDE RECOMPILAR          MAT00120
C                      ESTA RUTINA, ELIMINANDOLO)                       MAT00130
C  M: ORDEN DE LA MATRIZ A INVERTIR.                                    MAT00140
C  IRET : CODIGO DE RETORNO.                                            MAT00150
C          0 : LA EJECUCION FUE EXITOSA.                                MAT00160
C         -1 : LA MATRIZ ES SINGULAR.                                   MAT00170
C  IENT : CODIGO DE ENTRADA.                                            MAT00180
C          0 : SI QUIERE REALIZARSE LA INVERSION EN UNA UNICA ENTRADA.  MAT00190
C          1 : PRIMERA ENTRADA. SI QUIERE REALIZARSE LA INVERSION EN    MAT00200
C              DOS ENTRADAS, HARA EL PROCESO DE ELIMINACION.            MAT00210
C          2 : SEGUNDA ENTRADA. SI QUIERE REALIZARSE LA INVERSION EN    MAT00220
C              DOS ENTRADAS, COMPLETARA LA INVERSION REALIZANDO         MAT00230
C              LA SUSTITUCION HACIA ATRAS.                              MAT00240
C                                                                       MAT00250
      DIMENSION A(50,50),IAUXI(1600),IAUXJ(1600)                          MAT00260
      
C      WRITE(11,60)
C  60  FORMAT(' ESTO ES LO QUE LE LLEGA A MATINV')    
C      DO51I=1,N  
C      DO51J=1,N
C   51  WRITE(11,61) I,J,A(I,J)
C   61  FORMAT(2I2,2X,E15.5)
      
      
      IRET=0                                                            MAT00270
      IF(IENT.EQ.2)GO TO 20                                             MAT00280
      DO 10 K=1,M                                                       MAT00290
      BIGA=A(K,K)                                                       MAT00300
      IAUXI(K)=K                                                        MAT00310
      IAUXJ(K)=K                                                        MAT00320
      DO 1 J=K,M                                                        MAT00330
      DO 1 I=K,M                                                        MAT00340
      IF(ABS(BIGA).GE.ABS(A(I,J)))GO TO 1                               MAT00350
      BIGA=A(I,J)                                                       MAT00360
      IAUXI(K)=I                                                        MAT00370
      IAUXJ(K)=J                                                        MAT00380
1     CONTINUE                                                          MAT00390
      I=IAUXI(K)                                                        MAT00400
      IF(I.LE.K)GO TO 2                                                 MAT00410
      DO 3 J=1,M                                                        MAT00420
      HOLD=-A(K,J)                                                      MAT00430
      A(K,J)=A(I,J)                                                     MAT00440
3     A(I,J)=HOLD                                                       MAT00450
2     J=IAUXJ(K)                                                        MAT00460
      IF(J.LE.K)GO TO 4                                                 MAT00470
      DO 5 I=1,M                                                        MAT00480
      HOLD=-A(I,K)                                                      MAT00490
      A(I,K)=A(I,J)                                                     MAT00500
5     A(I,J)=HOLD                                                       MAT00510
4     IF(BIGA.NE.0)GO TO 6                                              MAT00520
      IRET=-1                                                           MAT00530
      RETURN                                                            MAT00540
6     DO 7 I=1,M                                                        MAT00550
      IF(I.EQ.K)GO TO 7                                                 MAT00560
      A(I,K)=A(I,K)/(-BIGA)                                             MAT00570
7     CONTINUE                                                          MAT00580
      DO 8 I=1,M                                                        MAT00590
      HOLD=A(I,K)                                                       MAT00600
      DO 8 J=1,M                                                        MAT00610
      IF(I.EQ.K)GO TO 8                                                 MAT00620
      IF(J.EQ.K)GO TO 8                                                 MAT00630
      A(I,J)=HOLD*A(K,J)+A(I,J)                                         MAT00640
8     CONTINUE                                                          MAT00650
      DO 9 J=1,M                                                        MAT00660
      IF(J.EQ.K)GO TO 9                                                 MAT00670
      A(K,J)=A(K,J)/BIGA                                                MAT00680
9     CONTINUE                                                          MAT00690
10    A(K,K)=1./BIGA                                                    MAT00700
      IF(IENT.EQ.1)RETURN                                               MAT00710
20    K=M                                                               MAT00720
21    K=K-1                                                             MAT00730
      IF(K.LE.0)RETURN                                                  MAT00740
      I=IAUXI(K)                                                        MAT00750
      IF(I.LE.K)GO TO 22                                                MAT00760
      DO 23 J=1,M                                                       MAT00770
      HOLD=A(J,K)                                                       MAT00780
      A(J,K)=-A(J,I)                                                    MAT00790
23    A(J,I)=HOLD                                                       MAT00800
22    J=IAUXJ(K)                                                        MAT00810
      IF(J.LE.K)GO TO 21                                                MAT00820
      DO 24 I=1,M                                                       MAT00830
      HOLD=A(K,I)                                                       MAT00840
      A(K,I)=-A(J,I)                                                    MAT00850
24    A(J,I)=HOLD                                                       MAT00860
      GO TO 21                                                          MAT00870
      END                                                               MAT00880

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
      
