      program readTQTec
C     reads the output from program TQTec
      CHARACTER INFILE*20, DUMMY*20, OUTFILE*20
      REAL II
      INTEGER Q1
      DIMENSION S(7),E(7),C(7),A(7),DIFF(5000),PRODN(5000),
     *   S2(9),E2(9),C2(9),A2(9),S3(13),E3(13),C3(13),A3(13)
      COMMON /COM1/ II(10),Q(50000),R(50000,2,10),Y(10)
      WRITE(*,*)'NAME OF INPUT FILE'
      READ(*,100) INFILE
      OPEN (UNIT=8, FILE=INFILE)
      REWIND 8
      R1=1.99
C     XMIN=0.0
      READ (8,100) DUMMY
      DO 5 J=1,10
         READ(8,110) II(J)
5     CONTINUE
      XMIN=II(5)
      Q1=NINT(II(5)/(2*II(2)))
	  write (*,*) Q1
      DO 10 J=1,Q1
         READ(8,115) Q(J)
10    CONTINUE
      DO 25 K=1,10
         DO 25 J=1,2
            DO 25 I=1,Q1
               READ(8,120) R(I,J,K)
25    CONTINUE
C      DO 30 K=1,5
C         DO 30 I=1,Q1
C            READ(8,120) TTI(I,K)
C 30    CONTINUE
      DO 40 I=1,10
         READ(8,130) Y(I)
40    CONTINUE
100   FORMAT(A20)
110   FORMAT(F7.3)
115   FORMAT(F6.2)
120   FORMAT(F7.1)
130   FORMAT(F11.4)
170   FORMAT(I1)
      CLOSE(8)
C
C     CHOOSE PLOT DATA DESIRED
C
500   WRITE(*,*)'WHICH PLOT DATA DO YOU WANT? (ONE DIGIT ONLY)'
      WRITE(*,*)'TEMPERATURE =(1)'
      WRITE(*,*)'DEPTH       =(2)'
C      WRITE(*,*)'TTI         =(3)'
C      WRITE(*,*)'PROD. TYPE 1=(4)'
C      WRITE(*,*)'PROD. TYPE 2=(5)'
C      WRITE(*,*)'PROD. TYPE 3=(6)'
      WRITE(*,*)'SURFACE Q   =(3)'
	  WRITE(*,*) 'TIME@OUTPUT  = (4)'
      WRITE(*,*)'QUIT?       =(9)'
      READ(*,170) IPLT
      IF (IPLT.EQ.1) GOTO 510 
      IF (IPLT.EQ.2) GOTO 520
      IF (IPLT.EQ.3) GOTO 570 
	  IF (IPLT.EQ.4) GOTO 580
C      IF (IPLT.EQ.4) GOTO 540
C      IF (IPLT.EQ.5) GOTO 550
C      IF (IPLT.EQ.6) GOTO 560
C      IF (IPLT.EQ.7) GOTO 570
      IF (IPLT.EQ.9) GOTO 590
C
C     WRITE TEMPERATURE PLOT DATA
C
510   WRITE(*,*) 'TYPE NAME OF OUTPUT FILE FOR TEMP Data'
      READ(*,100) OUTFILE
      OPEN (UNIT=9,FILE=OUTFILE)
      WRITE(9,175) -XMIN,II(2)*2,Q1
      DO 50 L=1,Q1
		WRITE(9,160)(R(L,1,I),I=1,10)
50    CONTINUE
      CLOSE(9)
      GOTO 500
C
C     WRITE DEPTH PLOT DATA
C
520   WRITE(*,*)'TYPE NAME OF OUTPUT FILE FOR DEPTH Data'
      READ(*,100) OUTFILE
      OPEN (UNIT=9, FILE=OUTFILE)
      WRITE(9,175) -XMIN,II(2)*2,Q1
      DO 60 L=1,Q1
          WRITE(9,160)(-1*R(L,2,I)*II(1),I=1,10)
60    CONTINUE
      CLOSE (9)
      GOTO 500
C
C     WRITE TTI PLOT DATA
C
C530   WRITE(*,*)'TYPE NAME OF OUTPUT FILE FOR TTI PLOTS'
C      READ(*,100) OUTFILE
C      OPEN (UNIT=9,FILE=OUTFILE)
C         WRITE(9,175) -XMIN,II(2)*2,Q1,LOG10(1.),LOG10(100000.)
C      DO 70 L=1,5
C         DO 71 M=1,Q1
C            IF (TTI(M,L).LE.1.0) THEN
C               TTI(M,L)=0.0
C            ELSE
C               TTI(M,L)=LOG10(TTI(M,L))
C            ENDIF 
C71       CONTINUE
C         DO 75 K=10,Q1,10
C75          WRITE(9,165)(TTI(I,L),I=K-9,K)
C         IF (K.NE.Q1.AND.K.NE.(Q1+10)) THEN
C            KNEW=ABS(Q1-K)
C            WRITE(9,165)(TTI(I,L),I=Q1-KNEW+1,Q1)
C         ENDIF
C70    CONTINUE
C      CLOSE(9)
C      GOTO 500
C
C     WRITE PRODUCTION TYPE I PLOT DATA
C
540   WRITE(*,*)'TYPE NAME OF OUTPUT FILE FOR PROD TYPE I PLOTS'
      READ(*,100) OUTFILE
      OPEN(UNIT=9,FILE=OUTFILE)
         WRITE(9,175) -XMIN,II(2)*2,Q1,0.0,1000.0
      DATA E/48,50,52,54,56,58,60/,
     *   C/8,10,20,26,810,11,13/,
     *   A/7*3E29/
      DO 80 L=1,5
         DO 81 I=1,7
            S(I)=0.0
81       CONTINUE
         DO 88 ILOOP=1,Q1
            TEMP=0
            DO 89 M=1,7 
               S(M)=S(M)+EXP(-E(M)*1000/(R1*(R(ILOOP,1,L)+273)))
     *            *A(M)*II(2)*2
               TEMP=TEMP+C(M)*(1-EXP(-S(M)))
89          CONTINUE
            PRODN(ILOOP)=TEMP
88       CONTINUE
         DIFF(1)=0
         DO 84 I=2,Q1
           DIFF(I)=(PRODN(I)-PRODN(I-1))/II(2)*2
84       CONTINUE
            DO 85 K=10,Q1,10
85             WRITE(9,160)(PRODN(I),I=K-9,K)
            IF(K.LT.Q1) WRITE(9,160)(PRODN(I),I=K+1,Q1)
         IF (K.NE.Q1.AND.K.NE.(Q1+10)) THEN
            KNEW=ABS(Q1-K)
            WRITE(9,160)(PRODN(I),I=Q1-KNEW+1,Q1)
         ENDIF
80    CONTINUE
      CLOSE(9)
      GOTO 500
C
C     WRITE PRODUCTION TYPE II PLOT DATA
C
550   WRITE(*,*) 'TYPE NAME OF OUTPUT FILE FOR PROD TYPE II PLOTS'
      READ(*,100) OUTFILE
      OPEN(UNIT=9, FILE =OUTFILE)
      WRITE(9,175) -XMIN,II(2)*2,Q1,0.0,700.0
      DATA E2/40,46,48,50,52,54,56,58,60/,
     *C2/6,4,9,32,132,302,104,35,6/,
     *A2/9*3E29/
      DO 280 L=1,5
         DO 281 I=1,9
            S2(I)=0.0
281      CONTINUE
         DO 288 ILOOP=1,Q1
         TEMP=0
         DO 289 M=1,9
            S2(M)=S2(M)+EXP(-E2(M)*1000/(R1*(R(ILOOP,1,L)+273)))
     *           *A2(M)*II(2)*2
            TEMP=TEMP+C2(M)*(1-EXP(-S2(M)))
289      CONTINUE
         PRODN(ILOOP)=TEMP
288   CONTINUE
      DIFF(1)=0
      DO 284 I=2,Q1
         DIFF(I)=(PRODN(I)-PRODN(I-1))/II(2)*2
284   CONTINUE
         DO 285 K=10,Q1,10
285      WRITE(9,160)(PRODN(I),I=K-9,K)
         IF (K.NE.Q1.AND.K.NE.(Q1+10)) THEN
            KNEW=ABS(Q1-K)
            WRITE(9,160)(PRODN(I),I=Q1-KNEW+1,Q1)
         ENDIF
280   CONTINUE
      CLOSE(9)
      GOTO 500
C
C     WRITE PRODUCTION TYPE III PLOT DATA
C     
560   WRITE(*,*) 'TYPE NAME OF OUTPUT FILE FOR PROD TYPE III PLOTS'
      READ(*,100) OUTFILE
      OPEN (UNIT=9,FILE=OUTFILE)
         WRITE(9,175) -XMIN,II(2)*2,Q1,0.0,250.0
      DATA E3/50,52,54,56,58,60,62,64,66,68,70,72,74/,
     *C3/1,5,6,42,82,60,23,12,7,5,3,2,2/,
     *A3/13*3E29/
      DO 380 L=1,5
         DO 381 I=1,13
            S3(I)=0.0
381      CONTINUE
         DO 388 ILOOP=1,Q1
            TEMP=0
            DO 389 M=1,13
               S3(M)=S3(M)+EXP(-E3(M)*1000/(R1*(R(ILOOP,1,L)+273)))
     *               *A3(M)*II(2)*2 
               TEMP=TEMP+C3(M)*(1-EXP(-S3(M)))
389         CONTINUE
            PRODN(ILOOP)=TEMP
388      CONTINUE
         DIFF(1)=0
         DO 384 I=2,Q1
            DIFF(I)=(PRODN(I)-PRODN(I-1))/II(2)*2
384      CONTINUE
         DO 385 K=10,Q1,10
385      WRITE(9,160)(PRODN(I),I=K-9,K)
         IF (K.NE.Q1.AND.K.NE.(Q1+10)) THEN
            KNEW=ABS(Q1-K)
            WRITE(9,160)(PRODN(I),I=Q1-KNEW+1,Q1)
         ENDIF
380   CONTINUE
      CLOSE(9)
      GOTO 500
C
C     WRITE SURFACE HEAT FLOW DATA
C
570   WRITE(*,*) 'TYPE NAME OF OUTPUT FILE FOR SURF. Q Data'
      READ(*,100) OUTFILE
      OPEN (UNIT=9,FILE=OUTFILE)
      WRITE(9,175) -XMIN,II(2)*2,Q1
	  DO 575 I=1,Q1
575      WRITE(9,160) Q(I)
	  CLOSE(9)
      GOTO 500
C      Q(2)=(Q(1)+Q(2)+Q(3))/3.0 
C      Q(Q1-1)=(Q(Q1)+Q(Q1-1)+Q(Q1-2))/3.0
C      DO 571 I=3,Q1-2
C         Q(I)=(Q(I-2)+Q(I-1)+Q(I)+Q(I+1)+Q(I+2))/5
C 571   CONTINUE
C
C     WRITE TIME DATA
C
580   WRITE(*,*)'TYPE NAME OF OUTPUT FILE FOR TIME Data'
      READ(*,100) OUTFILE
      OPEN (UNIT=9, FILE=OUTFILE)
      WRITE(9,175) -XMIN,II(2)*2,Q1
      DO 75 L=1,Q1
          WRITE(9,160)(L*II(2)*2,I=1,10)
75    CONTINUE
      CLOSE (9)
      GOTO 500
C
150   FORMAT(I2,2F7.1,I7)
160   FORMAT(10F8.3)
165   FORMAT(10F7.3)
175   FORMAT(2F10.3,I6)
590   CONTINUE
      STOP
      END 
      
      
