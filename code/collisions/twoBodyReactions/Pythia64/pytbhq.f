C---------------------------------------------------------
C 2)  Q QBAR ->TBH^+
       SUBROUTINE PYTBHQ(Q1,Q2,P1,P2,P3,MT,MB,RMB,MHP,AMP2)
C
C AMP2(OUTPUT) =MATRIX ELEMENT (AMPLITUDE**2) FOR Q QBAR->TB H^+
C (NB SAME STRUCTURE AS FOR PYTBHG ROUTINE ABOVE)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      DOUBLE PRECISION MW2,MT,MB,MHP,MW
      DIMENSION Q1(4),Q2(4),P1(4),P2(4),P3(4)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYMSSM/IMSS(0:99),RMSS(0:99)
      COMMON/PYCTBH/ ALPHA,ALPHAS,SW2,MW2,TANB,VTB,V,A
      SAVE /PYDAT1/,/PYDAT2/,/PYMSSM/,/PYCTBH/
C !THE RELEVANT INPUT PARAMETERS ABOVE ARE NEEDED FOR CALCULATION
C BUT ARE NOT DEFINED HERE SO THAT ONE MAY CHOOSE/VARY THEIR VALUES:
C ACCORDINGLY, WHEN CALLING THESE SUBROUTINES, PLEASE SUPPLY VIA
C THIS COMMON/PARAM/ YOUR PREFERRED ALPHA, ALPHAS,..AND TANB VALUES
C
C THE NORMALIZED V,A COUPLINGS ARE DEFINED BELOW AND USED BOTH
C IN THIS ROUTINE AND IN THE TOP WIDTH CALCULATION PYTBHB(..).
C
      DIMENSION YY(2,2)
 
      PI = 4*DATAN(1.D0)
      MW = DSQRT(MW2)
 
C COLLECTING THE RELEVANT OVERALL FACTORS:
C 3X3 INITIAL QUARK COLOR AVERAGE, 2X2 QUARK SPIN AVERAGE
      PS=1.D0/(3.D0*3.D0 *2.D0*2.D0)
C COUPLING CONSTANT (OVERALL NORMALIZATION)
      FACT=(4.D0*PI*ALPHA)*(4.D0*PI*ALPHAS)**2/SW2/2.D0
C NB ALPHA IS E^2/4/PI, BUT BETTER DEFINED IN TERMS OF G_FERMI:
C ALPHA= DSQRT(2.D0)*GF*SW2*MW**2/PI
C ALPHAS IS ALPHA_STRONG;
C SW2 IS SIN(THETA_W)**2.
C
C      VTB=.998D0
C VTB IS TOP-BOTTOM CKM MATRIX ELEMENT (APPROXIMATE VALUE HERE)
C
      V = ( MT/MW/TANB +RMB/MW*TANB)/2.D0
      A = (-MT/MW/TANB +RMB/MW*TANB)/2.D0
C V AND A ARE (NORMALIZED) VECTOR AND AXIAL TBH^+ COUPLINGS
C
C REDEFINING P2 INGOING FROM OVERALL MOMENTUM CONSERVATION
C (BECAUSE P2 INGOING WAS USED IN OUR GRAPH CALCULATION CONVENTIONS)
      DO 100 KK=1,4
        P2(KK)=P3(KK)-Q1(KK)-Q2(KK)+P1(KK)
  100 CONTINUE
C DEFINING VARIOUS RELEVANT 4-SCALAR PRODUCTS:
      S = 2*PYTBHS(Q1,Q2)
      P1Q1=PYTBHS(Q1,P1)
      P1Q2=PYTBHS(P1,Q2)
      P2Q1=PYTBHS(P2,Q1)
      P2Q2=PYTBHS(P2,Q2)
      P1P2=PYTBHS(P1,P2)
C
C   TOP WIDTH CALCULATION
      CALL PYTBHB(MT,MB,MHP,BR,GAMT)
C   GAMT IS THE TOP WIDTH: T->BH^+ AND/OR T->B W^+
C THEN DEFINE TOP (RESONANT) PROPAGATOR:
      A1INV= S -2*P1Q1 -2*P1Q2
      A1 =A1INV/(A1INV**2+ (GAMT*MT)**2)
C (I.E. INTRODUCE THE TOP WIDTH IN A1 TO REGULARISE THE POLE)
C  NB  A12 = A1*A1 BUT WITH CORRECT WIDTH TREATMENT
      A12 = 1.D0/(A1INV**2+ (GAMT*MT)**2)
      A2 =1.D0/(S +2*P2Q1 +2*P2Q2)
C NOTE A2 IS B PROPAGATOR, DOES NOT NEED A WIDTH
C  NOW COMES THE AMP**2:
C NB COLOR FACTOR (COMING FORM GRAPHS) ALREADY INCLUDED IN
C THE EXPRESSIONS BELOW
      YY(1, 1) = -16*A**2*A2**2*MB*MT+
     &64*A**2*A2**2*P1Q2*P2Q1**2/S**2+
     &128*A**2*A2**2*MB*MT*P2Q1*P2Q2/S**2-
     &128*A**2*A2**2*P1P2*P2Q1*P2Q2/S**2-
     &64*A**2*A2**2*P1Q1*P2Q1*P2Q2/S**2-
     &64*A**2*A2**2*P1Q2*P2Q1*P2Q2/S**2+
     &64*A**2*A2**2*P1Q1*P2Q2**2/S**2-
     &32*A**2*A2**2*MB**3*MT/S+32*A**2*A2**2*MB**2*P1P2/S+
     &32*A**2*A2**2*MB**2*P1Q1/S+32*A**2*A2**2*MB**2*P1Q2/S-
     &32*A**2*A2**2*P1P2*P2Q1/S-32*A**2*A2**2*P1Q1*P2Q1/S-
     &32*A**2*A2**2*P1P2*P2Q2/S-32*A**2*A2**2*P1Q2*P2Q2/S+
     &16*A2**2*MB*MT*V**2+64*A2**2*P1Q2*P2Q1**2*V**2/S**2-
     &128*A2**2*MB*MT*P2Q1*P2Q2*V**2/S**2-
     &128*A2**2*P1P2*P2Q1*P2Q2*V**2/S**2-
     &64*A2**2*P1Q1*P2Q1*P2Q2*V**2/S**2-
     &64*A2**2*P1Q2*P2Q1*P2Q2*V**2/S**2+
     &64*A2**2*P1Q1*P2Q2**2*V**2/S**2
      YY(1, 1)=YY(1, 1)+32*A2**2*MB**3*MT*V**2/S+
     &32*A2**2*MB**2*P1P2*V**2/S+
     &32*A2**2*MB**2*P1Q1*V**2/S+32*A2**2*MB**2*P1Q2*V**2/S-
     &32*A2**2*P1P2*P2Q1*V**2/S-32*A2**2*P1Q1*P2Q1*V**2/S-
     &32*A2**2*P1P2*P2Q2*V**2/S-32*A2**2*P1Q2*P2Q2*V**2/S
      YY(1, 1)=2*YY(1, 1)
 
      YY(1, 2) = -32*A**2*A1*A2*MB*MT+
     &128*A**2*A1*A2*MB*MT*P1Q2*P2Q1/S**2-
     &128*A**2*A1*A2*P1P2*P1Q2*P2Q1/S**2+
     &64*A**2*A1*A2*P1Q1*P1Q2*P2Q1/S**2-
     &64*A**2*A1*A2*P1Q2**2*P2Q1/S**2+
     &64*A**2*A1*A2*P1Q2*P2Q1**2/S**2+
     &128*A**2*A1*A2*MB*MT*P1Q1*P2Q2/S**2-
     &128*A**2*A1*A2*P1P2*P1Q1*P2Q2/S**2-
     &64*A**2*A1*A2*P1Q1**2*P2Q2/S**2+
     &64*A**2*A1*A2*P1Q1*P1Q2*P2Q2/S**2-
     &64*A**2*A1*A2*P1Q1*P2Q1*P2Q2/S**2-
     &64*A**2*A1*A2*P1Q2*P2Q1*P2Q2/S**2+
     &64*A**2*A1*A2*P1Q1*P2Q2**2/S**2-
     &64*A**2*A1*A2*MB*MT*P1P2/S+
     &64*A**2*A1*A2*P1P2**2/S+32*A**2*A1*A2*MB**2*P1Q1/S+
     &32*A**2*A1*A2*P1P2*P1Q1/S+32*A**2*A1*A2*MB**2*P1Q2/S+
     &32*A**2*A1*A2*P1P2*P1Q2/S-32*A**2*A1*A2*MT**2*P2Q1/S
      YY(1, 2)=YY(1, 2)-32*A**2*A1*A2*P1P2*P2Q1/S-
     &64*A**2*A1*A2*P1Q1*P2Q1/S-
     &32*A**2*A1*A2*MT**2*P2Q2/S-32*A**2*A1*A2*P1P2*P2Q2/S-
     &64*A**2*A1*A2*P1Q2*P2Q2/S+32*A1*A2*MB*MT*V**2-
     &128*A1*A2*MB*MT*P1Q2*P2Q1*V**2/S**2 -
     &128*A1*A2*P1P2*P1Q2*P2Q1*V**2/S**2+
     &64*A1*A2*P1Q1*P1Q2*P2Q1*V**2/S**2-
     &64*A1*A2*P1Q2**2*P2Q1*V**2/S**2+
     &64*A1*A2*P1Q2*P2Q1**2*V**2/S**2-
     &128*A1*A2*MB*MT*P1Q1*P2Q2*V**2/S**2-
     &128*A1*A2*P1P2*P1Q1*P2Q2*V**2/S**2-
     &64*A1*A2*P1Q1**2*P2Q2*V**2/S**2+
     &64*A1*A2*P1Q1*P1Q2*P2Q2*V**2/S**2-
     &64*A1*A2*P1Q1*P2Q1*P2Q2*V**2/S**2-
     &64*A1*A2*P1Q2*P2Q1*P2Q2*V**2/S**2+
     &64*A1*A2*P1Q1*P2Q2**2*V**2/S**2+
     &64*A1*A2*MB*MT*P1P2*V**2/S+64*A1*A2*P1P2**2*V**2/S
      YY(1, 2)=YY(1, 2)+32*A1*A2*MB**2*P1Q1*V**2/S+
     &32*A1*A2*P1P2*P1Q1*V**2/S+
     &32*A1*A2*MB**2*P1Q2*V**2/S+32*A1*A2*P1P2*P1Q2*V**2/S-
     &32*A1*A2*MT**2*P2Q1*V**2/S-32*A1*A2*P1P2*P2Q1*V**2/S-
     &64*A1*A2*P1Q1*P2Q1*V**2/S-32*A1*A2*MT**2*P2Q2*V**2/S-
     &32*A1*A2*P1P2*P2Q2*V**2/S-64*A1*A2*P1Q2*P2Q2*V**2/S
 
 
      YY(2, 2) =-16*A**2*A12*MB*MT+
     &128*A**2*A12*MB*MT*P1Q1*P1Q2/S**2-
     &128*A**2*A12*P1P2*P1Q1*P1Q2/S**2+
     &64*A**2*A12*P1Q1*P1Q2*P2Q1/S**2-
     &64*A**2*A12*P1Q2**2*P2Q1/S**2-64*A**2*A12*P1Q1**2*P2Q2/S**2+
     &64*A**2*A12*P1Q1*P1Q2*P2Q2/S**2-32*A**2*A12*MB*MT**3/S+
     &32*A**2*A12*MT**2*P1P2/S+32*A**2*A12*P1P2*P1Q1/S+
     &32*A**2*A12*P1P2*P1Q2/S-32*A**2*A12*MT**2*P2Q1/S-
     &32*A**2*A12*P1Q1*P2Q1/S-32*A**2*A12*MT**2*P2Q2/S-
     &32*A**2*A12*P1Q2*P2Q2/S+16*A12*MB*MT*V**2-
     &128*A12*MB*MT*P1Q1*P1Q2*V**2/S**2-
     &128*A12*P1P2*P1Q1*P1Q2*V**2/S**2+
     &64*A12*P1Q1*P1Q2*P2Q1*V**2/S**2-
     &64*A12*P1Q2**2*P2Q1*V**2/S**2-64*A12*P1Q1**2*P2Q2*V**2/S**2+
     &64*A12*P1Q1*P1Q2*P2Q2*V**2/S**2+32*A12*MB*MT**3*V**2/S+
     &32*A12*MT**2*P1P2*V**2/S+32*A12*P1P2*P1Q1*V**2/S+
     &32*A12*P1P2*P1Q2*V**2/S-32*A12*MT**2*P2Q1*V**2/S
      YY(2, 2)=YY(2, 2)-32*A12*P1Q1*P2Q1*V**2/S-
     &32*A12*MT**2*P2Q2*V**2/S-
     &32*A12*P1Q2*P2Q2*V**2/S
      YY(2, 2)=2*YY(2, 2)
 
      RES=YY(1,1)+2*YY(1,2)+YY(2,2)
      AMP2=  FACT*PS*VTB**2*RES
 
      END