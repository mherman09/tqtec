
      PROGRAM MAIN
C     MAIN PROGRAM TQTec.FOR
C     BUR-ERO-THRUST:UPPER OR LOWER PLATE 
C     CALCULATES THE ONE-DIMENSIONAL TRANSIENT THERMAL FIELD WITHIN
C     AN AREA THAT UNDERGOES EPISODES OF:  BURIAL, EROSION, AND 
C     THRUSTING.  THIS PROGRAM CAN BE USED TO MONITOR POINTS THAT
C     MOVE FROM THE UPPER TO LOWER (or v.v) PLATES.
C     --------------------------------------------------------------
C     VARIABLE LIST:  MAIN
C     RESP = RESPONSE TO:  WANT TO CREATE A NEW INPUT DATA FILE?
C     INFIL = NAME OF INPUT DATA FILE
C     OUTFIL = NAME OF OUTPUT DATA FILE
C     Q1 = TOTAL TIME (LATER, Q1=(TOTAL TIME/TIME STEP LENGTH)=
C        NUMBER OF TIME STEPS
C     M1 = TIME/OUTPUT (LATER, M1=((TIME/OUTPUT)/TIME STEP LENGTH)=
C        NUMBER OF TIME STEPS BETWEEN DISPLAY INTERVALS
C     W(1) = SURFACE TEMPERATURE
C     G1 = SURFACE HEAT FLOW
C     Y(I) = INITIAL DEPTH OF POINTS
C     N = NUMBER OF SPACE STEPS (NODES)
C     H1 = VERTICAL HEIGHT OF NODES
C     K1 = TIME STEP
C     C1 = CONDUCTIVITY
C     A1 = HEAT PRODUCTION
C     B1 = DEPTH OF HEAT PRODUCTION
C     D1 = DIFFUSIVITY
C     W1 = D1*K1/C1 = (A SCALING FACTOR USED TO GET TEMPERATURES)
C     II(1) = H1
C     II(2) = K1
C     II(3) = A1
C     II(4) = B1
C     II(5) = Q1
C     II(6) = D1
C     R1 = D1*K1/(H1*H1) PART OF WHAT IS INSERTED IN THE TRIDIAGONAL
C        MATRIX WHEN FINITE DIFFERENCE IS USED TO SOLVE THE PROBLEM
C     W(3) = TEMPERATURE AT BOTTOM NODE + CHANGE IN TEMPERATURE 
C        WITHOUT HEAT PRODUCTION = TEMP AT NODE N+1
C     V = A COUNTER IN THE P-MATRIX, V'S ARE THE TIME STEPS
C     P(V) = SEE P(I) IN SUBROUTINE HISTI
C     Q(V) = SURFACE HEAT FLOW AT EACH TIME STEP
C     TTI(V,K) = TTI CALCULATED AT EACH TIME STEP
C     --------------------------------------------------------------
      CHARACTER INFIL*20, OUTFIL*20, RESP*1, OUT1*25
      INTEGER P, V, E1, Q1,THTYPE,THNUMB
      REAL II, K1
      COMMON /COM1/ R(50000,2,10), A(5000,3), B(5000), C(5000), D(5000),
     *   E(5000), H(5000), P(50000), Q(50000), T(5000), II(10), Z(5,3),
     *   Y(10), NN(4), W(3),THTYPE(5)
      COMMON /COM2/ H1, A1, B1, C1, W1, E1, R1, N
      COMMON /COM3/ COND(5000),BCOND(50000)
      COMMON /COM4/ INL,TOP(50),THICK(50),ACOND(50)
      WRITE (*,*) 'DO YOU WANT TO CREATE A NEW INPUT DATA_FILE
     *   (Y/N)?'
      READ (*,130) RESP
C      IF (RESP.EQ.'Y') THEN
      if (resp.eq.'Y'.or.resp.eq.'y') then
         WRITE (*,*) 'NAME OF INPUT DATA_FILE: ? '
         READ (*,140) INFIL
         OPEN (UNIT=8, FILE=INFIL)
         CALL INPUT(INFIL)
      ELSE
         WRITE (*,*) 'NAME OF INPUT DATA_FILE: ? '
         READ (*,140) INFIL
         OPEN (UNIT=8, FILE=INFIL)
      ENDIF
      WRITE(*,*) 'NAME OF OUTPUT DATA_FILE: ? '
      READ(*,140) OUTFIL
      OPEN(UNIT=7, FILE=OUTFIL)
      REWIND 8
      READ (8,110) Q1,M1,W(1),G1,C1,A1,B1
      READ (8,150) INL
      IF (INL.GT.0) THEN
         DO 2, I=1,INL
            READ(8,160)TOP(I),THICK(I),ACOND(I)
2        CONTINUE
      ENDIF
      READ (8,120) (Y(I), I=1,10)
      N=1200
      DO 1, I=1,N
         COND(I)=C1
1     CONTINUE
      H1=0.05
      K1=0.005
C      A1=0.0
C      B1=10.0
      D1=32.0
      W1=D1*K1/C1
      II(1)=H1
      II(2)=K1
      II(3)=A1
      II(4)=B1
      II(5)=Q1
      II(6)=D1
      II(7)=W(1)
      Q1=NINT(FLOAT(Q1)/K1)
      M1=M1/K1
      R1=D1*K1/(H1*H1)
      V=0
      CALL INIT(G1)
	  CALL HIST
5     W(3)=B(N)+W(2)
      CALL MAT
	  CALL TRID
	   V=V+1
      WRITE(*,*) 'Did cycle', V
	   IF(E1.NE.0)WRITE(*,*)'*****ERROR EXIT FROM TRID*****'
      IF(P(V).EQ.1)CALL BURIAL(V)
      IF(P(V).EQ.2)CALL EROS
      IF(P(V).GE.3) THEN
          THNUMB = P(V)-2
         IF (THTYPE(THNUMB).EQ.1) CALL THSTUP(V)
         IF (THTYPE(THNUMB).EQ.2) CALL THSTLP(V)
      ENDIF
      Q(V)=(B(10)-B(5))/(5.0*H1)
      SURCON=0.0
      coninv = 0.0
      DO 6, I=1,5
C         coninv = coninv + (1.0/COND(I))
		 		 SURCON=SURCON+COND(I+4)
6     CONTINUE
      SURCON=SURCON/5.0
C      SURCON = (1.0/coninv)/25.0
      Q(V)=Q(V)*SURCON
      DO 10 I=1,10
         IARG=NINT(Y(I))
         IF(IARG.EQ.0) THEN
           R(V,1,I)=W(1)
         ELSEIF(IARG.LT.0) THEN
           R(V,1,I)=0.0
         ELSEIF(IARG.GT.0) THEN
           R(V,1,I)=B(IARG)
         ENDIF
         R(V,2,I)=Y(I)
C         EMP = (R(V,1,I)-105)/10
C         IF (V.EQ.1) THEN
C           TTI(V,I)= II(2)*2**EMP 
C         ELSE
C           TTI(V,I)= TTI(V-1,I) + II(2)*2**EMP
C         ENDIF
10    CONTINUE
      IF(V.GE.Q1)THEN
        CALL OUTPUT(Q1,M1,OUTFIL)
      ELSE
        GOTO 5
      ENDIF
110   FORMAT(/,2I10,5F10.4)
120   FORMAT(10F8.4)
130   FORMAT(A1)
140   FORMAT(A20)
150   FORMAT(I10)
160   FORMAT(3F10.4)
      CLOSE (8)
      CLOSE (7)
      STOP
      END
C     **************************************************************
C     **************************************************************
C     PROGRAM INPUT.FOR
C     INPUT DATA SET FOR LOWPL.FOR AND UPRPL.FOR
      SUBROUTINE INPUT(FILNM)
C     --------------------------------------------------------------
C     VARIABLES:  INPUT
C     IQ1 = TOTAL TIME
C     IM1 = TIME/OUTPUT
C     AW2 = SURFACE TEMPERATURE
C     AG1 = SURFACE HEAT FLOW
C     AY(I) = INITIAL DEPTH OF POINTS 1 TO 10
C     FILNM = NAME OF INPUT FILE
C     INBP = NUMBER OF BURIAL PERIODS
C     AN(1) = BEGINNING OF A BURIAL, EROSION, OR THRUST PERIOD AFTER
C        START OF MODEL
C     AN(2) = DURATION OF BURIAL OR EROSION PERIOD
C     AN(3) = TOTAL BURIAL OR EROSION OF THAT PERIOD
C        (NOTE:  THERE ARE AN'S FOR EACH BURIAL OR EROSION PERIOD)
C     INUEP = NUMBER OF UPLIFT/EROSION PERIODS
C     INTP = NUMBER OF THRUST PERIODS
C     AZ(I,1) = INITIAL BASE OF THRUST FOR EACH THRUST EVENT
C     AZ(I,2) = INITIAL DEPTH OF THRUST
C     AZ(I,3) = INITIAL THICKNESS OF THRUST
C     --------------------------------------------------------------
      CHARACTER FILNM*20
	  COMMON /COM2/ H1, A1, B1, C1, W1, E1, R1, N
      REAL AZ(5,3), AY(10), AN(4), AM1, TOP(5000), THICK(5000), 
     *     ACOND(5000)
      WRITE(*,*)'TOTAL TIME FOR MODEL (MA): ?'
      READ(*,*) IQ1
C      WRITE(*,*) 'TIME INTERVAL FOR OUTPUT TO BE DISPLAYED: ?'
C      READ(*,*) IM1
      IM1=5
      WRITE(*,*)'TEMPERATURE AT UPPER SURFACE BOUNDARY: ?'
      READ(*,*) AW2
      WRITE(*,*) 'SURFACE HEAT FLOW (mW/m2): ?'
      READ(*,*) AG1
      WRITE(*,*) 'INITIAL (BASEMENT) THERMAL CONDUCTIVITY (W/m K): ?'
      READ(*,*) C1
	  WRITE(*,*) 'Surface Heat Production (uW/m3):?'
	  READ(*,*) A1
      WRITE(*,*)'DO YOU WANT TO ACCOUNT FOR VARIATIONS IN'
      WRITE(*,*)'THERMAL CONDUCTIVITY AT THE START OF THE'
      WRITE(*,*)'MODEL? (1=YES)'
      READ(*,*)IREPLY
      IF(IREPLY.EQ.1)THEN
         WRITE(*,*)'NUMBER OF LAYERS TO INPUT CONDUCTIVITY FOR?'
         READ(*,*) INL
         DO 9, I=1,INL
            WRITE(*,*)'DEPTH (KM) OF TOP OF LAYER',I,'?'
            READ(*,*)TOP(I)
            WRITE(*,*)'THICKNESS OF LAYER (KM) ',I,'?'
            READ(*,*)THICK(I)
            WRITE(*,*)'CONDUCTIVITY (W/m K) OF LAYER',I,'?'
            READ(*,*)ACOND(I)
9        CONTINUE
      ENDIF
      DO 10 I=1,10
         WRITE(*,*) 'INITIAL DEPTH (KM) OF POINT',I,' ? '
         READ(*,*) AY(I)
10    CONTINUE
      WRITE(8,100)FILNM
      WRITE(8,110)IQ1,IM1,AW2,AG1,C1,A1
      IF(IREPLY.EQ.1)THEN
      WRITE(8,130)INL
         DO 11, I=1,INL
            WRITE(8,140)TOP(I),THICK(I),ACOND(I)
11       CONTINUE
      ELSE
      INL=0
      WRITE(8,130)INL
      ENDIF
      WRITE(8,120) (AY(I),I=1,10)
      WRITE(*,*)'NUMBER OF BURIAL PERIODS: ?'
      READ(*,*)INBP
      WRITE(8,130)INBP
      DO 20 I=1,INBP
         WRITE(*,*)'BEGINNING (MA AFTER START) OF BURIAL',I,' ? '
         READ(*,*)AN(1)
         WRITE(*,*)'DURATION (MA) OF BURIAL',I,' ? '
         READ(*,*)AN(2)
         WRITE(*,*)'TOTAL BURIAL (KM) EPISODE',I,' ? '
         READ(*,*)AN(3)
         WRITE(*,*)'THERMAL CONDUCTIVITY (W/m K) OF SEDIMENTS'
         WRITE(*,*)'BURIAL EPISODE',I,'?'
         READ(*,*)AN(4)
         WRITE(8,150)(AN(J),J=1,4)
20    CONTINUE
      WRITE(*,*)'NUMBER OF UPLIFT/EROSION PERIODS: ?'
      READ(*,*)INUEP
      WRITE(8,130)INUEP
      DO 30 I=1,INUEP
         WRITE(*,*)'BEGINNING (MA AFTER START) OF UPLIFT',I,' ? '
         READ(*,*)AN(1)
         WRITE(*,*)'DURATION (MA) OF UPLIFT',I,' ? '
         READ(*,*)AN(2)
         WRITE(*,*)'TOTAL UPLIFT (KM) PERIOD',I,' ? '
         READ(*,*)AN(3)
         WRITE(8,140)(AN(J),J=1,3)
30    CONTINUE
      WRITE(*,*)'NUMBER OF THRUST PERIODS: ?'
      READ(*,*)INTP
      WRITE(8,130)INTP
      DO 40 I=1,INTP
         WRITE(*,*)'TIME OF THRUSTING (MA) OF PERIOD',I,' ? '
         READ(*,*)AN(1)
         WRITE(*,*)'POINTS IN UPPER (1) OR LOWER (2) PLATE ?'
         READ(*,*)AN(2)
         WRITE(*,*)'INITIAL BASE (KM) OF THRUST: ?'
         READ(*,*)AZ(I,1)
         WRITE(*,*)'INITIAL DEPTH (KM) OF THRUST FAULT: ?'
         READ(*,*)AZ(I,2)
         WRITE(*,*)'INITIAL  THICKNESS (KM) OF THRUST: ?'
         READ(*,*)AZ(I,3)
         WRITE(8,160)AN(1),AN(2),(AZ(I,J),J=1,3)
40    CONTINUE
100   FORMAT(A20)
110   FORMAT(2I10,4F10.4)
120   FORMAT(10F8.4)
130   FORMAT(I10)
140   FORMAT(3F10.4)
150   FORMAT(4F10.4)
160   FORMAT(5F10.4)
      RETURN
      END
C     **************************************************************
C     **************************************************************
C     PROGRAM INIT.FOR
      SUBROUTINE INIT(G1)
C     --------------------------------------------------------------
C     VARIABLES:  INIT
C     W(2)=(SURFACE HEAT FLOW - (HEAT PRODUCTION*DEPTH OF HEAT PROD)/
C        CONDUCTIVITY = REDUCED HEAT FLOW/CONDUCTIVITY (K/m)=dT/dz
C        GRADIENT WITHOUT HEAT PRODUCTION.  LATER,
C        W(2)=W(2)*NODE SIZE=CHANGE IN TEMPERATURE WITHOUT HEAT PROD.
C     Y(I) = INITIAL DEPTH OF POINTS/NODE SIZE = INITIAL DEPTH IN
C        TERMS OF NODES
C     G1 = SURFACE HEAT FLOW * NODE SIZE / CONDUCTIVITY = CHANGE
C        IN TEMPERATURE WITH HEAT PRODUCTION (K)
C     ARG1 = z/D  FOR HEAT PRODUCTION (z=DEPTH, D=DEPTH OF HEAT PROD.)
C     H(I) = A*exp(-z/D) = EXPONENTIALLY DECAYING HEAT PROD. 
C        IN TERMS OF NODES (I'S) HERE.
C     B(I) = To + qo*z/k + Ao*D*D/k(z/D + exp(z/D) - 1.0)
C        STEADY STATE, ONE-DIMENSIONAL HEAT CONDUCTION EQUATION,
C        IN TERMS OF I'S HERE.
C     Q(1) = k(T(node1) - T(surface)) / dz = SURFACE HEAT FLOW
C     --------------------------------------------------------------
      INTEGER P,V,E1,TN,BN,THTYPE
      REAL II
      COMMON /COM1/ R(50000,2,10), A(5000,3), B(5000), C(5000), D(5000),
     *   E(5000), H(5000), P(50000), Q(50000), T(5000), II(10), Z(5,3),
     *   Y(10), NN(4), W(3),THTYPE(5)
      COMMON /COM2/ H1,A1,B1,C1,W1,E1,R1,N
      COMMON /COM3/ COND(5000),BCOND(50000)
	  COMMON /COM4/ INL,TOP(50),THICK(50),ACOND(50)
      IF(INL.EQ.0)GOTO 7
      DO 5, I=1,INL
         TN=INT(TOP(I)/H1)
         BN=INT((TOP(I)+THICK(I))/H1)
         IF(TN+1.EQ.BN)THEN 
            COND(TN+1)=ACOND(I)
         ELSE
            DO 6, J=TN+1,BN
               COND(J)=ACOND(I)
6           CONTINUE
         ENDIF
5     CONTINUE
7     W(2)=(G1-(A1*B1))/C1
      W(2)=W(2)*H1
      DO 10 I=1,10
         Y(I)=Y(I)/H1
10    CONTINUE
      DO 20 I=1,N
         ARG1=FLOAT(I)*H1/B1
C        H(I) IS THE VOLUMETRIC HEAT PRODUCTION AT EACH NODE, BASED
C           ON EXPONENTIALLY DECAYING HEAT PRODUCTION BUT DISCRETIZED
         H(I)=A1*EXP(-ARG1)
20    CONTINUE
C     COMPUTE TEMPERATURE IN LAYER 1, SURFACE TO NODE 1
      B(1) = W(1) + G1*H1/COND(1) - H(1)*H1**2/(2.*COND(1))
      HTHT = H(1)*H1
      QBASE = G1-HTHT
      DO 30 I=2,N
         B(I)=B(I-1)+QBASE*H1/COND(I)-H(I)*H1**2/(2.*COND(I))
         HTHT = H(I)*H1
         QBASE=QBASE-HTHT
30    CONTINUE
      Q(1)=(COND(1)*(B(1)-W(1)))/H1
      RETURN
      END
C     **************************************************************
C     **************************************************************
C     PROGRAM HIST(HISTORY).FOR
      SUBROUTINE HIST
C     --------------------------------------------------------------
C     VARIABLES:  HIST
C     NBP = # OF BURIAL PERIODS
C     NUEP = # OF UPLIFT AND EROSION PERIODS
C     NTP = # OF THRUST PERIODS
C     XN(1) = BEGINNING OF BURIAL, EROSION, OR THRUST FOR EACH PERIOD
C     XN(2) = DURATION OF BURIAL, EROSION PERIOD
C     XN(3) = TOTAL BURIAL OR EROSION FOR EACH PERIOD
C     NN(1) = TIME STEP # WHERE BURIAL OR EROSION BEGINS
C     NN(2) = DURATION OF BURIAL OR EROSION IN TERMS OF TIME STEPS.
C     NN(3) = TOTAL BURIAL OR EROSION IN TERMS OF NODES.
C     R2 = RATE OF BURIAL IN TERMS OF NODES AND TIME STEPS
C     JBEG = TIME NODE WHERE BURIAL OR EROSION FIRST TAKES PLACE 
C        = NN(1)+1
C     JEND = LAST TIME NODE FOR BURIAL FOR EACH BURIAL OR EROSION
C        PERIOD
C     ARG = INTERMEDIATE VALUE TO DETERMINE WHETHER OR NOT BURIAL OR 
C        EROSION TAKES PLACE.  IF AMOUNT BETWEEN TWO TIME NODES IS
C        GREATER THAN ONE THEN ONE NODE IS ADDED (OR SUBTRACTED).  IF 
C        NOT, IT IS TESTED UNTIL ENOUGH TO MAKE UP ONE NODE IS ADDED 
C        OR SUBTRACTED.  ONLY ONE NODE IS ADDED OR SUBTRACTED AT A 
C        TIME.
C     P(J) = J'S ARE THE TIME STEPS.  WHEN P(J)=1 BURIAL OCCURS.
C        WHEN P(J)=2 EROSION OCCURS.  WHEN P(J)>2 IT IS A THRUST
C        EVENT, WHERE J=3 IS THE FIRST THRUST, J=4 IS THE SECOND ETC.
C     ---------------------------------------------------------------
      INTEGER P,E1,F,THTYPE
      REAL XN(4),II
      COMMON /COM1/ R(50000,2,10), A(5000,3), B(5000), C(5000), D(5000),
     *   E(5000), H(5000), P(50000), Q(50000), T(5000), II(10), Z(5,3),
     *   Y(10), NN(4), W(3),THTYPE(5)
      COMMON /COM2/ H1,A1,B1,C1,W1,E1,R1,N
      COMMON /COM3/ COND(5000),BCOND(50000)
C     HISTORY OF BURIAL PERIODS
      READ(8,100)NBP
      DO 20 I=1,NBP
         READ(8,120) (XN(J),J=1,4)
         NN(1)=INT(XN(1)/II(2))
         NN(2)=INT(XN(2)/II(2))
         NN(3)=INT(XN(3)/II(1))
         F=0
         R2=FLOAT(NN(3))/FLOAT(NN(2))
         JBEG=NN(1)+1
         JEND=NN(1)+NN(2)
         DO 10 J=JBEG,JEND
            ARG=((J-NN(1))*R2)-FLOAT(F)
            IF(ARG.GE.1.0)THEN
               P(J)=1
               BCOND(J)=XN(4)
               F=F+1
            ENDIF
10       CONTINUE
20    CONTINUE
C     HISTORY OF UPLIFT/EROSION PERIODS
      READ(8,100)NUEP
      DO 40 I=1,NUEP
         READ(8,110)(XN(J),J=1,3)
         NN(1)=INT(XN(1)/II(2))
         NN(2)=INT(XN(2)/II(2))
         NN(3)=INT(XN(3)/II(1))
         F=0
         R2=FLOAT(NN(3))/FLOAT(NN(2))
         JBEG=NN(1)+1
         JEND=NN(1)+NN(2)
         DO 30 J=JBEG,JEND
            ARG=((J-NN(1))*R2)-FLOAT(F)
            IF(ARG.GE.1.0)THEN
               P(J)=2
               F=F+1
            ENDIF
30    CONTINUE
40    CONTINUE
C     HISTORY OF THRUST PERIODS
      READ(8,100)NTP
      DO 50 I=1,NTP
         READ(8,130)XN(1),XN(2),(Z(I,J),J=1,3)
         NN(1)=INT(XN(1)/II(2))
         P(NN(1))=2+I
         THTYPE(I)= XN(2)
50    CONTINUE
100   FORMAT(I10)
110   FORMAT(3F10.4)
120   FORMAT(4F10.4)
130   FORMAT(5F10.4)
      RETURN
      END
C     **************************************************************
C     **************************************************************
C     PROGRAM MAT.FOR
      SUBROUTINE MAT
C     --------------------------------------------------------------
C     VARIABLES:  MAT
C     B(I) = INITIAL TEMPERATURES, LATER AT END OF ROUTINE, BECOMES
C        NEW TEMPERATURES.
C     T(I) = NEW TEMPERATURES
C     --------------------------------------------------------------
      INTEGER P,E1,THTYPE
      REAL II
      DIMENSION BET(5000), GAM(5000)
      COMMON /COM1/ R(50000,2,10), A(5000,3), B(5000), C(5000), D(5000),
     *   E(5000), H(5000), P(50000), Q(50000), T(5000), II(10), Z(5,3),
     *   Y(10), NN(4), W(3),THTYPE(5)
      COMMON /COM2/ H1,A1,B1,C1,W1,E1,R1,N
      COMMON /COM3/ COND(5000),BCOND(50000)
      DO 5, I=2,N
         BET(I)=1/(1/COND(I)+1/COND(I-1))
         GAM(I)=1/COND(I)
5     CONTINUE
C     FOR UPPERMOST NODE
      BET(1)=1/(1/COND(1)+1/COND(1))
      GAM(1)=1/COND(1)
      C(1)=-R1*BET(1)*GAM(1)
      E(1)=-R1*BET(2)*GAM(1)
      D(1)=1+R1*(BET(2)*GAM(1)+BET(1)*GAM(1))
      A(1,1)=-C(1)
      A(1,3)=-E(1)
      A(1,2)=2-D(1)
      DO 10 I=2,N-1
         C(I)=-R1*BET(I)*GAM(I)
         E(I)=-R1*BET(I+1)*GAM(I)
         D(I)=1+R1*(BET(I)*GAM(I)+BET(I+1)*GAM(I))
         A(I,1)=-C(I)
         A(I,3)=-E(I)
         A(I,2)=2-D(I)
10    CONTINUE
C     FINAL NODE
      C(N)=-R1*BET(N)*GAM(N)
      E(N)=C(N)
      D(N)=1+R1*2*BET(N)*GAM(N)
      A(N,1)=-C(N)
      A(N,3)=-E(N)
      A(N,2)=2-D(N)
      T(1)=A(1,2)*B(1)+A(1,3)*B(2)+2.0*A(1,1)*W(1)+II(6)*II(2)*H(1)/
     *     COND(1) 
      DO 20 I=2,N-1
         T(I)=A(I,1)*B(I-1)+A(I,2)*B(I)+A(I,3)*B(I+1)+II(6)*II(2)*H(I)/
     *          COND(I)
20    CONTINUE
      T(N)=A(N,1)*B(N-1)+A(N,2)*B(N)+2.0*A(N,3)*W(3)+W1*H(N)
      DO 30 I=1,N
         B(I)=T(I)
30    CONTINUE
      RETURN
      END
C     **************************************************************
C     **************************************************************
C     PROGRAM TRID.FOR
      SUBROUTINE TRID
C     --------------------------------------------------------------
      INTEGER P,E1,THTYPE
      REAL II
      COMMON /COM1/ R(50000,2,10), A(5000,3), B(5000), C(5000), D(5000),
     *   E(5000), H(5000), P(50000), Q(50000), T(5000), II(10), Z(5,3),
     *   Y(10), NN(4), W(3),THTYPE(5)
      COMMON /COM2/ H1,A1,B1,C1,W1,E1,R1,N
      E1=0
      C(1)=D(1)
      IF(N-1.LT.1)GOTO 40
      D(1)=E(1)
      E(1)=0.0
      E(N)=E(1)
      DO 30 K=1,N-1
         IF(ABS(C(K+1)).LT.ABS(C(K)))GOTO 10
         T1=C(K+1)
         C(K+1)=C(K)
         C(K)=T1
         T1=D(K+1)
         D(K+1)=D(K)
         D(K)=T1
         T1=E(K+1)
         E(K+1)=E(K)
         E(K)=T1
         T1=B(K+1)
         B(K+1)=B(K)
         B(K)=T1
10       IF(C(K).NE.0.0)GOTO 20
         E1=K
         RETURN
20       T1=(-C(K+1)/C(K))
         C(K+1)=D(K+1)+(T1*D(K))
         D(K+1)=E(K+1)+(T1*E(K))
         E(K+1)=0.0
         B(K+1)=B(K+1)+(T1*B(K))
30    CONTINUE
40    IF(C(N).NE.0.0) GOTO 50
      E1=N
      RETURN
50    B(N)=B(N)/C(N)
      IF(N.EQ.1) RETURN
      B(N-1)=(B(N-1)-(D(N-1)*B(N)))/C(N-1)
      IF((N-2).LT.1)RETURN
      DO 60 L=1,N-2
         K=N-2-L+1
         B(K)=(B(K)-(D(K)*B(K+1))-(E(K)*B(K+2)))/C(K)
60    CONTINUE
      RETURN
      END
C     **************************************************************
C     **************************************************************
C     PROGRAM BURIAL.FOR
      SUBROUTINE BURIAL(V)
C     --------------------------------------------------------------
      INTEGER P,E1,V,THTYPE
      REAL II
      COMMON /COM1/ R(50000,2,10), A(5000,3), B(5000), C(5000), D(5000),
     *   E(5000), H(5000), P(50000), Q(50000), T(5000), II(10), Z(5,3),
     *   Y(10), NN(4), W(3),THTYPE(5)
      COMMON /COM2/ H1,A1,B1,C1,W1,E1,R1,N
      COMMON /COM3/ COND(5000),BCOND(50000)
      DO 10 I=1,N-1
         J=N-I
         B(J+1)=B(J)
         H(J+1)=H(J)
         COND(J+1)=COND(J)
10    CONTINUE
      B(1)=W(1)
      H(1)=0.0
      COND(1)=BCOND(V)
      DO 20 I=1,10
         Y(I)=Y(I)+1.0
20    CONTINUE
      RETURN
      END
C     **************************************************************
C     **************************************************************
C     PROGRAM EROS(EROSION).FOR
      SUBROUTINE EROS
C     --------------------------------------------------------------
      INTEGER P,E1,THTYPE
      REAL II
      COMMON /COM1/ R(50000,2,10), A(5000,3), B(5000), C(5000), D(5000),
     *   E(5000), H(5000), P(50000), Q(50000), T(5000), II(10), Z(5,3),
     *   Y(10), NN(4), W(3),THTYPE(5)
      COMMON /COM2/ H1,A1,B1,C1,W1,E1,R1,N
      COMMON /COM3/ COND(5000),BCOND(50000)
      DO 10 I=2,N
         B(I-1)=B(I)
         H(I-1)=H(I)
         COND(I-1)=COND(I)
10    CONTINUE
      B(N)=B(N-1)+W(2)
      H(N)=0.0
      COND(N)=C1
      DO 20 I=1,10
         Y(I)=Y(I)-1.0
C        IF(Y(I).LE.0.0) Y(I)=0.0
20    CONTINUE
      RETURN
      END
C     **************************************************************
C     **************************************************************
C     PROGRAM THSTUP.FOR:  THRUST FOR UPPER PLATE BOUNDARY
      SUBROUTINE THSTUP(V)
C     --------------------------------------------------------------
      INTEGER P,E1,V,THTYPE
      REAL II
      DIMENSION UPLC(300), UPLH(300), UPLB(300)
      COMMON /COM1/ R(50000,2,10), A(5000,3), B(5000), C(5000), D(5000),
     *   E(5000), H(5000), P(50000), Q(50000), T(5000), II(10), Z(5,3),
     *   Y(10), NN(4), W(3),THTYPE(5)
      COMMON /COM2/ H1,A1,B1,C1,W1,E1,R1,N
      COMMON /COM3/ COND(5000),BCOND(50000)
      K=P(V)-2
C     NN(1)=THICKNESS PRIOR TO THRUSTING
C     NN(2)=DEPTH OF EMPLACEMENT
C     NN(3)=FINAL THICKNESS OF THRUST SHEET
      NN(1)=INT(Z(K,1)/II(1))
      NN(2)=INT(Z(K,2)/II(1))
      NN(3)=INT(Z(K,3)/II(1))
C     COPY THE PART OF THE ARRAY THAT WILL BE THRUSTED AND ERODE OFF 
C     AMOUNT THAT GETS ERODED DURING THRUST EVENT.
      IERO=NN(1)-NN(3)
      DO 1, I=IERO+1,NN(1)
         UPLC(I-IERO)=COND(I)
         UPLH(I-IERO)=H(I)
         UPLB(I-IERO)=B(I)
1     CONTINUE
C     REMOVE THE OLD ARRAY (LOWER PLATE) DOWN TO THE DEPTH OF 
C     EMPLACEMENT AND MOVE THE REST OF THE ARRAY DOWN TO MAKE ROOM
C     FOR THE UPPER PLATE.
      DO 2, I=NN(2)+1,N
         COND(I-NN(2))=COND(I)
         H(I-NN(2))=H(I)
         B(I-NN(2))=B(I)
2     CONTINUE
      DO 3, I=N,NN(3)+1,-1
         COND(I)=COND(I-NN(3))
         H(I)=H(I-NN(3))
         B(I)=B(I-NN(3))
3     CONTINUE
C     PUT THE TWO ARRAYS TOGETHER
      DO 4, I=1,NN(3)
         COND(I)=UPLC(I)
         H(I)=UPLH(I)
         B(I)=UPLB(I)
4     CONTINUE
C     MOVE POINTS OF INTEREST AROUND
C     FOR UPPER PLATE:
      DO 20 I=1,10
         Y(I)=Y(I)-FLOAT(IERO)
         IF(Y(I).LE.0.0) Y(I)=0.0
20    CONTINUE
      DO 50 J=1,10
         B(1)=(W(1)+B(2))/2.0
         B(2)=(B(1)+B(2)+B(3))/3.0
         DO 40 I=3, 2*NN(1)
            B(I)=(B(I-2)+B(I-1)+B(I)+B(I+1)+B(I+2))/5.0
40       CONTINUE
50    CONTINUE
      RETURN
      END
C     **************************************************************
C     **************************************************************
C     PROGRAM OUTPUT.FOR
      SUBROUTINE OUTPUT(Q1,M1,OUTFIL)
C     --------------------------------------------------------------
      CHARACTER OUTFIL*20
      INTEGER P,E1,Q1,THTYPE
      REAL II
      COMMON /COM1/ R(50000,2,10), A(5000,3), B(5000), C(5000), D(5000),
     *   E(5000), H(5000), P(50000), Q(50000), T(5000), II(10), Z(5,3),
     *   Y(10), NN(4), W(3),THTYPE(5)
      COMMON /COM2/ H1,A1,B1,C1,W1,E1,R1,N
      COMMON /COM3/ COND(5000),BCOND(50000)
      WRITE(7,100)OUTFIL
      DO 5, J=1,10
         WRITE(7,110) II(J)
5     CONTINUE
C      DO 8 I=1,N
C		 WRITE (*,117) H(I), COND(I)
C 8     CONTINUE	
       WRITE (*,117) W(1), W(2)	 
      DO 10, J=2,Q1,2
         WRITE(7,115)Q(J)
10    CONTINUE
      DO 15, K=1,10
         DO 20, J=1,2
            DO 25, I=2,Q1,2
               WRITE(7,120) R(I,J,K)
25          CONTINUE
20       CONTINUE
15    CONTINUE
C      DO 19, K = 1,5
C         DO 18, I = 2, Q1, 2
C            WRITE (7,120) TTI(I,K)
C18       CONTINUE
C19    CONTINUE  
      DO 30, I=1,10
         WRITE(7,130) Y(I)
         WRITE(*,130) Y(I)
C         WRITE(*,120) TTI(Q1,I)
30    CONTINUE
100   FORMAT(A20)
110   FORMAT(F7.3)
115   FORMAT(F6.2)
117   FORMAT(2F6.2)
120   FORMAT(F7.1)
130   FORMAT(F11.4)
      RETURN
      END
C     **************************************************************
C     **************************************************************
C     PROGRAM THSTLP.FOR:  THRUST FOR UPPER PLATE BOUNDARY
      SUBROUTINE THSTLP(V)
C     --------------------------------------------------------------
      INTEGER P,E1,V,THTYPE
      REAL II
      DIMENSION UPLC(300), UPLH(300), UPLB(300)
      COMMON /COM1/ R(50000,2,10), A(5000,3), B(5000), C(5000), D(5000),
     *   E(5000), H(5000), P(50000), Q(50000), T(5000), II(10), Z(5,3),
     *   Y(10), NN(4), W(3),THTYPE(5)
      COMMON /COM2/ H1,A1,B1,C1,W1,E1,R1,N
      COMMON /COM3/ COND(5000),BCOND(50000)
      K=P(V)-2
C     NN(1) = THICKNESS PRIOR TO THRUSTING
C     NN(2) = DEPTH OF EMPLACEMENT
C     NN(3) = FINAL THICKNESS OF THRUST SHEET
      NN(1)=INT(Z(K,1)/II(1))
      NN(2)=INT(Z(K,2)/II(1))
      NN(3)=INT(Z(K,3)/II(1))
      I1=NN(3)-NN(2)
C     COPY THE PAART OF THE ARRAY THAT WILL BE THRUSTED AND ERODE OFF
C     THE AMOUNT THAT GETS ERODED DURING THE THRUST EVENT.
      IERO = NN(1) - NN(3)
      DO 1, I=IERO+1, NN(1)
         UPLC(I-IERO)=COND(I)
         UPLH(I-IERO)=H(I)
         UPLB(I-IERO)=B(I)
1     CONTINUE
C     REMOVE THE OLD ARRAY (LOWER PLATE) DOWN TO THE DEPTH OF
C     EMPLACEMENT AND MOVE THE REST OF THE ARRAY DOWN TO MAKE ROOM 
C     FOR THE UPPER PLATE.
      DO 2, I=NN(2)+1, N
         COND(I-NN(2))=COND(I)
         H(I-NN(2))=H(I)
         B(I-NN(2))=B(I)
2     CONTINUE
      DO 3, I=N,NN(3)+1,-1
         COND(I)=COND(I-NN(3))
         H(I)=H(I-NN(3))
         B(I)=B(I-NN(3))
3     CONTINUE
C     PUT THE TWO ARRAYS TOGETHER
      DO 4, I=1,NN(3)
         COND(I)=UPLC(I)
         H(I)=UPLH(I)
         B(I)=UPLB(I)
4     CONTINUE
C     MOVE POINTS OF INTEREST AROUND
C     FOR LOWER PLATE:
      DO 20, I=1,10
         Y(I)=Y(I)+FLOAT(I1)
20    CONTINUE
      DO 50 J=1,10
         B(1)=(W(1)+B(2))/2.0
         B(2)=(B(1)+B(2)+B(3))/3.0
         DO 40 I=3, 2*NN(1)
            B(I)=(B(I-2)+B(I-1)+B(I)+B(I+1)+B(I+2))/5.0
40       CONTINUE
50    CONTINUE
      RETURN
      END
      
         
      
