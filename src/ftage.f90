  ! FT DISTRIBUTION MODELLING
  ! MATT LEGG
  !
  ! FILENAME = ftmodel.f90
  ! 
  ! Produces F.T. length distributions for 10 sediments of variable depth;
  ! Also gives final track count and F.T. age;
  !
  ! Designed to accept input temperature files written using TQTec4.f and
  ! readTQTec.f (Furlong)
  !____________________________________________________________________________________

  ! MAIN PROGRAM
   CHARACTER(20) INFILE, OUTFILE, TERMFILE, TTIFILE
   INTEGER Q1
   Dimension Q2(50000), T(50000,10), TERM(50000,10), TERM2(50000,10), RJ(50000,10), &
	 PL(50000,10,20), TEMP(50000,10), TERM3(50000,10), TERM4(50000,10), &
	 GR(50000,10), TIME(50000,10), Y(10,20), Y2(10,20), LLO(20), PPO(20),&
	 TOT0(10), TOT1(10), TOT2(10), AGE0(10), AGE1(10), AGE2(10), AGE3(10), &
	 FS(10,20), TOT3(50000,10), TOT4(10), FTOT(10,20), FTOT2(10), FTOT5(10),&
	 FTOT4(10,20), PIT(10), FTAGE(10), PAGE(50000), LINE(10), LINE2(10), &
	 TTI(50000,10), RN(50000,10), VRO(10)
   INTEGER :: OPENSTATUS, INPUTSTATUS, I, J, L, LINE, LINE2
   REAL :: Q2, T, TERM, TERM2, RJ, PL, TEMP, TERM3, TERM4, GR, TIME, Y, Y2,&
		   LLO, PPO, TOT0, TOT1, TOT2, AGE0, AGE1, AGE2, AGE3, TOT3, TOT4,&
		   FTOT, FTOT2, FTOT4, FTOT5, PIT, FTAGE, PAGE, TTI, RN, VRO
   REAL :: K, R, H, N, Q, A, DTMA 
   !INITIAL F.T. DISTRIBUTION FROM EXPERIMENTAL DATA OF Green et al. (1986):
   REAL :: LO(20) = (/15.,15.,15.,15.,15.,16.,16.,16.,16.,16.,16.,16., &
					  17.,17.,17.,17.,17.,17.,17.,18./)
   REAL :: LOO(20) = (/0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,&
						0.,0.,0./)
   REAL	:: MM(20) = (/1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,&
					  15.,16.,17.,18.,19.,20./)
   REAL :: NUM(10) = (/0.,0.,0.,0.,0.,0.,0.,0.,0.,0./)
					  
  !VARIABLES AND CONSTANTS EXPLAINED:
  !TEMP = TEMPERATURE [DEGREES C]
  !TERM = ARRAY OF TERMS TO SUM FOR RJ EQUATION
  !LO = INITIAL LENGTH
  !PL = FINAL LENGTH
  !RJ = PL/LO

  !DEFINE CONSTANTS:
   K=3.2976E-27		![kcal/K] ; Boltzmann's Constant
   R=0.00199		![kcal/mol/K] ; Universal Gas Constant
   H=1.58E-37		![kcal*s] ; Planck's Constant
   N=0.206			!dimensionless ; Exponent in power-law for 
					!initial radial defect distribution
   Q=40.6		    ![kcal/mol] ; Activation energy to eliminate defect
   A=1.81		    ![um] ; Rate constant for axial shortening
   LOAVE=16.2		!Ave. initial F.T. length
   PST=0.893		!Ratio of present day ave. length for spontaneous tracks
					!in Durango Apatite to LOAVE 
					!(Donelick and Miller, 1991)
   
 
  WRITE (*,*) "Enter Temperature data FileName: "
  READ (*,100) INFILE

  OPEN (UNIT=9, FILE =INFILE)
  REWIND 9

	READ(9,175) XMIN,DTMA,Q1

  	DO 30 I=1,Q1
	      READ(9,160) (T(I,J),J=1,10)
30 	CONTINUE   

	DO 34 J=1,10
		DO 35 I=1,Q1
			TEMP(I,J)=T(I,J)+273.
35		CONTINUE
34	CONTINUE 

	CLOSE(9)

	WRITE(*,175) XMIN,DTMA,Q1

	WRITE (*,*) TEMP(1,1), TEMP(2,2), TEMP(15500,10)

	!CONVERT TIME FROM MA TO SECONDS
        DT = DTMA*1000000*365*24*60*60
	WRITE (*,*) "DT=", DT
	
	!CREATE TERM Q2(I), THE AGE OF EACH TIME STEP
	DO I=1,Q1
		Q2(I)=-XMIN-(REAL(I)-1)*DTMA
	END DO

	!CALCULATE AGE OF EACH TIME STEP FOR PRINTOUT
	!DO I=1,Q1
	!	PAGE(I)=

	!CALCULATE ANNEALING FOR EACH TIME STEP: 
	!CARLSON (1990), EQUATION 4
	DO I=1,Q1
		DO J=1,10
			TERM(I,J)=TEMP(I,J)*DT*EXP(-Q/(R*TEMP(I,J)))
		END DO	 ! J LOOP
	END DO	 ! I LOOP

	!TURN ARRAY UPSIDE-DOWN TO SUM FROM I TO Q1:
	DO I=1,Q1
		DO J=1,10
				TERM2(I,J)=TERM(Q1-I+1,J)
		END DO ! J LOOP
	END DO ! I LOOP

	!SUM ANNEALING TERMS FROM I TO Q1:
	DO I=1,Q1
		DO J=1,10
			TERM3(I,J)=TERM2(I,J)+TERM2(I-1,J)
			TERM2(I,J)=TERM3(I,J)
		END DO
	END DO

	!PUT ARRAY RIGHT-SIDE-UP:
	DO I=1,Q1
		DO J=1,10
				TERM4(I,J)=TERM3(Q1-I+1,J)
		END DO ! J LOOP
	END DO ! I LOOP	
	
	WRITE (*,*) "term4;1,2: ", TERM4(1,1), TERM4(2,1)
	WRITE (*,*) "term4;100,101: ", TERM4(100,1), TERM4(101,1)
	WRITE (*,*) "term4;30000: ", TERM4(30000,1)

	!CALCULATE THE FINAL LENGTH PL:	
	!CARLSON (1990), EQUATION 4
	DO I=1,Q1
        DO J=1,10
			DO L=1,20
				PL(I,J,L)=LO(L)-(A)*((K/H)**N)*(TERM4(I,J)**N)
				IF(TEMP(I,J).LT.278) PL(I,J,L)=LOO(L)-(A)*((K/H)**N)*&
											  (TERM4(I,J)**N)
			END DO	!L LOOP
        END DO	!J LOOP
	END DO  !I LOOP
    print *,PL(6000,1,1)
    print *,PL(14000,1,1)
    print *,PL(28000,1,1)
	
	!PRINT *, "PL(50,1,10)=", PL(50,1,10)
	
	!ELIMINATE NEGATIVE LENGTHS 
	DO I=1,Q1
		DO J=1,10
			DO L=1,20
				IF(PL(I,J,L).LE.0)  PL(I,J,L)=0
				IF(PL(I,J,L).GT.0)	PL(I,J,L)=PL(I,J,L)
			END DO !L LOOP
		END DO	!J LOOP
	END DO  !I LOOP	
	
	!PRINT *, "PL(50,1,10)final=", PL(50,1,10)

	!DEFINE STARTING HISTOGRAM MATRIX
	DO J=1,10
		DO M=1,20
			Y(J,M) = 0
		END DO
	END DO
	
	!CREATE HISTOGRAM OF LENTHS IN BINS 1um WIDE FROM 1 to 20
	DO J=1,10
		DO I=1,Q1
			DO L=1,20
				IF(PL(I,J,L).GT.0.5.AND.PL(I,J,L).LE.1.5)   Y(J,1)=Y(J,1)+1
				IF(PL(I,J,L).GT.1.5.AND.PL(I,J,L).LE.2.5)   Y(J,2)=Y(J,2)+1
				IF(PL(I,J,L).GT.2.5.AND.PL(I,J,L).LE.3.5)   Y(J,3)=Y(J,3)+1
				IF(PL(I,J,L).GT.3.5.AND.PL(I,J,L).LE.4.5)   Y(J,4)=Y(J,4)+1
				IF(PL(I,J,L).GT.4.5.AND.PL(I,J,L).LE.5.5)   Y(J,5)=Y(J,5)+1
				IF(PL(I,J,L).GT.5.5.AND.PL(I,J,L).LE.6.5)   Y(J,6)=Y(J,6)+1
				IF(PL(I,J,L).GT.6.5.AND.PL(I,J,L).LE.7.5)   Y(J,7)=Y(J,7)+1
				IF(PL(I,J,L).GT.7.5.AND.PL(I,J,L).LE.8.5)   Y(J,8)=Y(J,8)+1
				IF(PL(I,J,L).GT.8.5.AND.PL(I,J,L).LE.9.5)   Y(J,9)=Y(J,9)+1
				IF(PL(I,J,L).GT.9.5.AND.PL(I,J,L).LE.10.5)  Y(J,10)=Y(J,10)+1
				IF(PL(I,J,L).GT.10.5.AND.PL(I,J,L).LE.11.5) Y(J,11)=Y(J,11)+1
				IF(PL(I,J,L).GT.11.5.AND.PL(I,J,L).LE.12.5) Y(J,12)=Y(J,12)+1
				IF(PL(I,J,L).GT.12.5.AND.PL(I,J,L).LE.13.5) Y(J,13)=Y(J,13)+1
				IF(PL(I,J,L).GT.13.5.AND.PL(I,J,L).LE.14.5) Y(J,14)=Y(J,14)+1
				IF(PL(I,J,L).GT.14.5.AND.PL(I,J,L).LE.15.5) Y(J,15)=Y(J,15)+1
				IF(PL(I,J,L).GT.15.5.AND.PL(I,J,L).LE.16.5) Y(J,16)=Y(J,16)+1
				IF(PL(I,J,L).GT.16.5.AND.PL(I,J,L).LE.17.5) Y(J,17)=Y(J,17)+1
				IF(PL(I,J,L).GT.17.5.AND.PL(I,J,L).LE.18.5) Y(J,18)=Y(J,18)+1
				IF(PL(I,J,L).GT.18.5.AND.PL(I,J,L).LE.19.5) Y(J,19)=Y(J,19)+1
				IF(PL(I,J,L).GT.19.5.AND.PL(I,J,L).LE.20.5) Y(J,20)=Y(J,20)+1				
			END DO !L LOOP
		END DO  !I LOOP	
	END DO	!J LOOP
	
	!SUM TRACKS IN UNSEGMENTED DISTRIBUTION
	DO J=1,10
		DO M=1,20
			IF(M.EQ.1) TOT0(J)=Y(J,M)
			IF(M.GT.1) TOT0(J)=Y(J,M)+TOT0(J)
		END DO	!M LOOP
	END DO	!J LOOP
	
	!CALCULATE RATIO OF BIN LENGTH TO AVE. INITIAL LENGTH (16.2 um), LLO
	DO M=1,20
		LLO(M) = MM(M)/16.2
	END DO	!M LOOP
	
	!CALCULATE CORRECTION FACTOR,PPO, FOR INTERMEDIATE LENGTHS
	DO M=1,20
		PPO(M) = (2.862*LLO(M))-1.2104
	END DO	!M LOOP
	
	!ACCOUNT FOR SEGMENTED LENGTHS FOR BINS 0-11 
	!A LA CARLSON (1990) EMPERICAL FIT TO FIG. 8B
	DO J=1,10
		DO M=1,20
				IF(M.LT.12) Y(J,M)=Y(J,M)-((-0.1*M)+1.2)*Y(J,M)
				IF(M.GE.12) Y(J,M)=Y(J,M)
				IF(Y(J,M).LE.0.) Y(J,M)=0.
				IF(M.LT.12) Y2(J,M)=Y2(J,M)-((-0.1*M)+1.2)*Y2(J,M)
				IF(M.GE.12) Y2(J,M)=Y2(J,M)
		END DO	!M LOOP
	END DO	!J LOOP

	!CORRECT HISTOGRAM FOR ETCHING/USER BIAS (INCLUDED BIN 11!)
	!WILLETT (1997) EQUATION 4
	DO J=1,10
		DO M=1,20
			IF(M.LE.6) Y2(J,M)=0
			IF(M.GT.6.AND.M.LE.11) Y2(J,M)=Y(J,M)*((2.862*LLO(M))-1.2104)
			IF(M.GT.11) Y2(J,M)=Y(J,M)
			IF(Y2(J,M).LT.0) Y2(J,M)=0
		END DO !M LOOP
	END DO !J LOOP	

	!SUM TRACKS IN SEGMENTED DISTRIBUTIONS 
	!(UNCORRECTED=TOT1, CORRECTED=TOT2)
	DO J=1,10
		DO M=1,20
			IF(M.EQ.1) TOT1(J)=Y(J,M)
			IF(M.GT.1) TOT1(J)=Y(J,M)+TOT1(J)
			IF(M.EQ.1) TOT2(J)=Y2(J,M)
			IF(M.GT.1) TOT2(J)=Y2(J,M)+TOT2(J)
		END DO	!M LOOP
	END DO	!J LOOP
	
	!AND MULTIPLY BY TIME STEP DURATION (DTMA)
	DO I=1,Q1
		DO J=1,10
			DO M=1,20
				IF(M.EQ.1) TOT3(I,J)=PL(I,J,M)
				IF(M.GT.1.AND.M.LT.20) TOT3(I,J)=PL(I,J,M)+TOT3(I,J)
				IF(M.EQ.20) TOT3(I,J)=(((TOT3(I,J)+PL(I,J,M))/20)/16.2)*DTMA
			END DO !M LOOP
		END DO !J LOOP
	END DO !I LOOP
	
	!CALCULATE F.T. AGE
	!KETCHAM ET AL.(2005), EQUATION 14
	DO I=1,Q1
		DO J=1,10
			IF(I.EQ.1) AGE3(J)=TOT3(I,J)
			IF(I.GT.1.AND.I.LT.Q1) AGE3(J)=AGE3(J)+TOT3(I,J)
			IF(I.EQ.Q1) AGE3(J)=(AGE3(J)+TOT3(I,J))*(1/PST)
		END DO !J LOOP
	END DO !I LOOP
	
	!CALCULATE TOTAL CORRECTED TRACK LENGTH IN ORDER TO &
	!NORMALIZE F.T. AGE
	DO J=1,10
		DO M=1,20
			IF(M.EQ.1) FTOT(J,M)=Y2(J,M)*1.
			IF(M.EQ.2) FTOT(J,M)=Y2(J,M)*2.
			IF(M.EQ.3) FTOT(J,M)=Y2(J,M)*3.
			IF(M.EQ.4) FTOT(J,M)=Y2(J,M)*4.
			IF(M.EQ.5) FTOT(J,M)=Y2(J,M)*5.
			IF(M.EQ.6) FTOT(J,M)=Y2(J,M)*6.
			IF(M.EQ.7) FTOT(J,M)=Y2(J,M)*7.
			IF(M.EQ.8) FTOT(J,M)=Y2(J,M)*8.
			IF(M.EQ.9) FTOT(J,M)=Y2(J,M)*9.
			IF(M.EQ.10) FTOT(J,M)=Y2(J,M)*10.
			IF(M.EQ.11) FTOT(J,M)=Y2(J,M)*11.
			IF(M.EQ.12) FTOT(J,M)=Y2(J,M)*12.
			IF(M.EQ.13) FTOT(J,M)=Y2(J,M)*13.
			IF(M.EQ.14) FTOT(J,M)=Y2(J,M)*14.
			IF(M.EQ.15) FTOT(J,M)=Y2(J,M)*15.
			IF(M.EQ.16) FTOT(J,M)=Y2(J,M)*16.
			IF(M.EQ.17) FTOT(J,M)=Y2(J,M)*17.
			IF(M.EQ.18) FTOT(J,M)=Y2(J,M)*18.
			IF(M.EQ.19) FTOT(J,M)=Y2(J,M)*19.
			IF(M.EQ.20) FTOT(J,M)=Y2(J,M)*20.
		END DO !M LOOP
	END DO !J LOOP
	
	DO J=1,10
		DO M=1,20
			IF(M.EQ.1) FTOT2(J)=FTOT(J,M)
			IF(M.GT.1) FTOT2(J)=FTOT2(J)+FTOT(J,M)
		END DO !M LOOP
	END DO !J LOOP	
	
	!CALCULATE TOTAL UNCORRECTED TRACK LENGTH
	DO J=1,10
		DO M=1,20
			IF(M.EQ.1) FTOT4(J,M)=Y(J,M)*1.
			IF(M.EQ.2) FTOT4(J,M)=Y(J,M)*2.
			IF(M.EQ.3) FTOT4(J,M)=Y(J,M)*3.
			IF(M.EQ.4) FTOT4(J,M)=Y(J,M)*4.
			IF(M.EQ.5) FTOT4(J,M)=Y(J,M)*5.
			IF(M.EQ.6) FTOT4(J,M)=Y(J,M)*6.
			IF(M.EQ.7) FTOT4(J,M)=Y(J,M)*7.
			IF(M.EQ.8) FTOT4(J,M)=Y(J,M)*8.
			IF(M.EQ.9) FTOT4(J,M)=Y(J,M)*9.
			IF(M.EQ.10) FTOT4(J,M)=Y(J,M)*10.
			IF(M.EQ.11) FTOT4(J,M)=Y(J,M)*11.
			IF(M.EQ.12) FTOT4(J,M)=Y(J,M)*12.
			IF(M.EQ.13) FTOT4(J,M)=Y(J,M)*13.
			IF(M.EQ.14) FTOT4(J,M)=Y(J,M)*14.
			IF(M.EQ.15) FTOT4(J,M)=Y(J,M)*15.
			IF(M.EQ.16) FTOT4(J,M)=Y(J,M)*16.
			IF(M.EQ.17) FTOT4(J,M)=Y(J,M)*17.
			IF(M.EQ.18) FTOT4(J,M)=Y(J,M)*18.
			IF(M.EQ.19) FTOT4(J,M)=Y(J,M)*19.
			IF(M.EQ.20) FTOT4(J,M)=Y(J,M)*20.
		END DO !M LOOP
	END DO !J LOOP
	
	DO J=1,10
		DO M=1,20
			IF(M.EQ.1) FTOT5(J)=FTOT4(J,M)
			IF(M.GT.1) FTOT5(J)=FTOT5(J)+FTOT4(J,M)
		END DO !M LOOP
	END DO !J LOOP	
	
	DO J=1,10
		DO M=1,20
			IF(M.EQ.1) FTOT2(J)=FTOT(J,M)
			IF(M.GT.1) FTOT2(J)=FTOT2(J)+FTOT(J,M)
		END DO !M LOOP
	END DO !J LOOP				
	
	!CORRECT F.T. AGE BY MULTIPLYING BY THE RATIO OF FINAL & 
	!TRACK LENGTH TO INITIAL AVERAGE (16.2um)
	DO J=1,10
		FTAGE(J)=(FTOT2(J)/FTOT5(J))*AGE3(J)
	END DO !J LOOP
	
	!CALCULATE RETENTION AGES
	DO J=1,10
		AGE0(J)=(TOT0(J)/20)*DTMA
		AGE1(J)=(TOT1(J)/20)*DTMA
		AGE2(J)=(TOT2(J)/20)*DTMA
	END DO	!J LOOP
	
	!PRINT *, "Uncorrected Rentention Age=", AGE1
	!PRINT *, "Corrected Retention Age=", AGE2
	!PRINT *, "Uncorrected F.T. AGE=", AGE3
	!PRINT *, "Corrected F.T. AGE =", FAGE
		
	DO J=1,10
		LINE(J)  = NINT((-XMIN/DTMA)-(FTAGE(J)/DTMA))
		LINE2(J) = NINT((-XMIN/DTMA)-(AGE2(J)/DTMA))
	END DO	!J LOOP
	
	
	PRINT *, "LINE2=", LINE2
	
	!WRITE FINAL LENGTHS TO OUTPUT FILE:
	WRITE(*,*)'TYPE NAME OF F.T. OUT FILE'
	READ *, OUTFILE
	OPEN (UNIT=11, FILE=OUTFILE)
	!	WRITE(11,*) "UNCORRECTED F.T. AGE="
	!	WRITE(11,185) 0,(AGE3(J),J=1,10)	
		WRITE(11,*) "CORRECTED F.T. AGE="
		WRITE(11,185) 0,(FTAGE(J),J=1,10)	
	!	WRITE(11,*) "UNSEGMENTED TOTAL="
	!	WRITE(11,180) 0,(TOT0(J),J=1,10)
	!	WRITE(11,*) "UNSEGMENTED RETENTION AGE="
	!	WRITE(11,185) 0,(AGE0(J),J=1,10)		
	!	WRITE(11,*) "UNCORRECTED HISTOGRAM"
	!	WRITE(11,*) "TOTAL="
	!	WRITE(11,180) 0,(TOT1(J),J=1,10)
	!	WRITE(11,*) "UNCORRECTED RETENTION AGE="
	!	WRITE(11,185) 0,(AGE1(J),J=1,10)
	!	WRITE(11,*) "Histogram:"
	!DO 50	M=1,20
	!	WRITE(11,180) M,(Y(J,M),J=1,10) 
	!50 CONTINUE
	!	WRITE(11,*) "CORRECTED HISTOGRAM"
	!	WRITE(11,*) "TOTAL="
	!	WRITE(11,180) 0,(TOT2(J),J=1,10)
		WRITE(11,*) "CORRECTED RETENTION AGE="
		WRITE(11,185) 0,(AGE2(J),J=1,10)
		WRITE(11,*) "CORRESPONDING TEMP="
		WRITE(11,185) 0,(T(LINE2(J),J),J=1,10)
		WRITE(11,*) "CORRECTED HISTOGRAM:" 			
	DO 60	M=1,20
		WRITE(11,180) M,(Y2(J,M),J=1,10) 
	60 CONTINUE
		WRITE(11,*) "TEMPERATURE DATA:"
	!PRINT EVERY 20th TERM IN TIME-TEMP HISTORY FOR GRAPHING
	DO 70	I=1,Q1,20
		WRITE(11,190) Q2(I),(T(I,J),J=1,10)
	70 CONTINUE
	CLOSE(11)
	!WRITE TTI/VRO FILE:
	DO I=1,Q1
		DO J=1,10
			IF(T(I,J).LT.30.)                       RN(I,J) = -500.
			IF(T(I,J).GE.30..AND.T(I,J).LT.40.)     RN(I,J) = -7.
			IF(T(I,J).GE.40..AND.T(I,J).LT.50.)     RN(I,J) = -6.
			IF(T(I,J).GE.50..AND.T(I,J).LT.60.)     RN(I,J) = -5.
			IF(T(I,J).GE.60..AND.T(I,J).LT.70.)     RN(I,J) = -4.
			IF(T(I,J).GE.70..AND.T(I,J).LT.80.)     RN(I,J) = -3.
			IF(T(I,J).GE.80..AND.T(I,J).LT.90.)     RN(I,J) = -2.
			IF(T(I,J).GE.90..AND.T(I,J).LT.100.)    RN(I,J) = -1.
			IF(T(I,J).GE.100..AND.T(I,J).LT.110.)   RN(I,J) =  0
			IF(T(I,J).GE.110..AND.T(I,J).LT.120.)   RN(I,J) =  1.
			IF(T(I,J).GE.120..AND.T(I,J).LT.130.)   RN(I,J) =  2.
			IF(T(I,J).GE.130..AND.T(I,J).LT.140.)   RN(I,J) =  3.
			IF(T(I,J).GE.140..AND.T(I,J).LT.150.)   RN(I,J) =  4.
			IF(T(I,J).GE.150..AND.T(I,J).LT.160.)   RN(I,J) =  5.
			IF(T(I,J).GE.160..AND.T(I,J).LT.170.)   RN(I,J) =  6.
			IF(T(I,J).GE.170..AND.T(I,J).LT.180.)   RN(I,J) =  7.	
		END DO !J LOOP
	END DO !I LOOP
	
	DO I=1,Q1
		DO J=1,10
			IF(I.EQ.1) TTI(I,J) = DTMA*(2.**RN(I,J))
			IF(I.GT.1) TTI(I,J) = TTI(I-1,J)+(DTMA*(2.**RN(I,J)))
		END DO !J LOOP
	END DO !I LOOP
	
	!Calculate VRO from TTI, based on power law fit to data &
	!In Table 4 (Waples, 1980)
	DO J=1,10
		VRO(J) = (TTI(Q1,J)/40.996)**(0.2177)
	END DO !J LOOP
	
	PRINT *, 'TTI1,6000=', TTI(600,1), 'VRO(1)=', VRO(1)
	
	!WRITE TTI/VRO OUTPUT FILE:
	WRITE(*,*)'TYPE NAME OF TTI OUT FILE'
	READ *, TTIFILE
	OPEN (UNIT=13, FILE=TTIFILE)
	WRITE (13,*) TTIFILE
	WRITE (13,*) "VRo="
		WRITE(13,190) 0.,(VRO(J),J=1,10)
	WRITE (13,*) "TTI="
		WRITE(13,190) 0.,(TTI(Q1,J),J=1,10)
	WRITE (13,*) "TTI list:"
	DO 80	I=1,Q1,20
		WRITE(13,190) Q2(I),(TTI(I,J),J=1,10)
	80 CONTINUE
	CLOSE(13)

100      FORMAT(A20)
160      FORMAT(10F8.3)
162		 FORMAT(10F9.1)
165      FORMAT(10F7.3)
170		 FORMAT(10ES10.3)
175      FORMAT(2F10.3,I6)
180		 FORMAT(I7,10F8.0)
185		 FORMAT(I7,10F8.2)
190		 FORMAT(F7.3,10F8.2)
        STOP
        END
		
		!REFERENCES
		
		! 1) Carlson WD (1990)  Mechanisms and kinetics of apatite
		!	 fission-track annealing; Am. Mineral., 75: 1120-1139
		!
		! 2) Domelick RA, Miller DS (1991)  Enhanced TINT fission track densities
		!    in low density apatites using 252Cf-derived fission fragment tracks:
		!    A model and experimental observations; Nucl. Tracks, 18: 301-307 
		!
		! 2) Green PF, Duddy IR, GLeadow AJW, Tingate PR, Laslett GM
		!    (1986)  Thermal Annealing of Fission Tracks in Apatite
		!	 1. A Qualitative Description; Chem. Geol. (Isotope Geoscience
		!    Section), 59: 237-253
		!
		! 3) Ketcham RA (2005)  Forward and Inverse Modeling of Low-
		!	 Temperature Thermochronometry Data; Rev. Mineral. Geochem.,
		!    58: 275-314
		
		
		
		
		
		
		
		
		
		
