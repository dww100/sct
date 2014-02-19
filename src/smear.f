C--------------------------------------------------------------------------
C
C Copyright 1981-2014 University College London
C
C Licensed under the Apache License, Version 2.0 (the "License");
C you may not use this file except in compliance with the License.
C You may obtain a copy of the License at
C
C    http://www.apache.org/licenses/LICENSE-2.0
C
C Unless required by applicable law or agreed to in writing, software
C distributed under the License is distributed on an "AS IS" BASIS,
C WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
C See the License for the specific language governing permissions and
C limitations under the License.
C
C--------------------------------------------------------------------------

       SUBROUTINE SMEAR(WAS,WAV,ALF,QDIG)
C Robbed without conscience directly from sctpl (ICCC version). Had
C to do a little bit of work, but I DIDN'T TOUCH THE MATHS!!
C Converted to double precision arithmetic: Thus,
C SQRT becomes DSQRT
C ALOG becomes DLOG
C EXP becomes DEXP
C SCATT() becomes I(), and increased the Q and I to 512
C I becomes II
C ************************************************************************
C *      SUBROUTINE SMEAR(SCATT,QDIG,Q,WAS,WAV,ALF)                      *
C *      TAKES THEORET CURVE IN SCATT() AND RETURNS IT SMEARED           *
C *       CORRECTED 24.8.82 : G() IS DIFFERENT AT EACH Q                 *
C *                                                                      *
C *      DIMENSION SCATT(450),Q(450),G(301),QQ(301)                      *
C *      DIMENSION SIG(301),SUM(301),F(301)                              *
C ************************************************************************
       IMPLICIT NONE
       INCLUDE 'SCTPL.COM'
       DOUBLE PRECISION WAS, WAV, ALF
       DOUBLE PRECISION G(301), QQ(301), SIG(301),SUM(301),F(301)
       DOUBLE PRECISION PI, WAB, CON, BON, AON, AA, D, QDIG, VV, SUMMM
       DOUBLE PRECISION TWO, EIGHT
       INTEGER IX, NQ, NQT, IQ, IQ1, IQ2, IQP, KOUNT, II, J, M
C
C       FIX PARAMETERS
      DO 1  IX=1,301
         G(IX)    =0.0 
         QQ(IX)   =0.0 
         SIG(IX)  =0.0 
         SUM(IX)  =0.0 
         F(IX)    =0.0 
 1    CONTINUE
      TWO = 2.0
      EIGHT = 8.0
      NQ=100
      NQT=201
      D=179.0
      PI=3.1415926536
C
C       CALCULATION OF THE SIGMA FOR THE SMEARING GAUSSIAN
C
C       LIKE CHAUVIN BUT 0.25/D**2 REMOVED AND 8LN2 NOT 2LN2
C       WAS AND ALPHA ARE FROM CUSACK JMB 1981 145, 539-541 
C       FOR D=1055,175,140 C=1060,560,280; A IS 0.003,0.010,0.016
C       FOR D11,D17; WAS IS 0.13,0.16 (THEOR 0.08,0.10)
C
      WAB=WAV/(2.0*PI)
      CON=(WAS*WAB)**2
      CON=CON*4.0
      BON= (ALF)**2
C      AON= WAB*DSQRT(8.0*DLOG(2.0))
      AON= WAB*DSQRT(EIGHT*DLOG(TWO))
C
C       CALCULATION OF LOOKUP TABLE FOR SIGMA SIG()
C
      DO 10 IQ=1,150
         IQ1=151-IQ
         IQ2=151+IQ
         QQ(IQ1)=-IQ*QDIG
         QQ(IQ2)= IQ*QDIG
         SIG(IQ2)=(DSQRT(CON*QQ(IQ2)*QQ(IQ2) + BON))/AON
         SIG(IQ1)=SIG(IQ2)
 10   CONTINUE
      QQ(151)=0.0
      SIG(151)=2.0*SIG(152) - SIG(153) 
C
C       APPLY SMEARING IN G() TO THEO F() TO GIVE SMEARED SUM()
C
      DO 25 IQ=1,100
C Replaced by I() 17/6/91 A.S.Nealis
C       F(151+IQ)=SCATT(IQ)
C       F(151-IQ)=SCATT(IQ)
        F(151+IQ)=I(PLOTCNT,IQ)
        F(151-IQ)=I(PLOTCNT,IQ)
 25    CONTINUE
C       SUBSTITUTION AT ZERO Q IS 1000
      F(151)=1000.0
C       NEED TO FILL UP F BEYOND CALC SCATT CURVE (ARBITRARY)
      DO 26 IQ=101,150
C Replaced by I() 17/6/91 A.S.Nealis
C       F(151+IQ)=SCATT(100)
C       F(151-IQ)=SCATT(100)
         F(151+IQ)=I(PLOTCNT,100)
         F(151-IQ)=I(PLOTCNT,100)
 26   CONTINUE
C
C       CONVOLUTION OF G() AND F()
C
      DO 40 IQ=1,301
C
C       CALCULATION OF GAUSSIAN G()
C
         DO 20 IQP=1,301
CCCCCCC NOTE DIFFERENCE FROM CHAUVIN AND CUSACK - IQ NOT IQP FOR AA, VV
            AA=1.0/(SIG(IQ )*DSQRT(2.0*PI))
            G(IQP)=0.0
            VV=(QQ(IQ)-QQ(IQP))/SIG(IQ )
            VV=VV*VV
            IF(VV.GT.76.0) GO TO 20
            G(IQP)=AA*DEXP(-VV/2.0) 
 20     CONTINUE
CCCCCCC SUM(IQ      =F(51)*G(51) + F(251)*G(251)
        DO 30 IQP=1,301
           SUM(IQ ) = SUM(IQ ) + ((F(IQP)*G(IQP)))
 30     CONTINUE
          SUM(IQ)=SUM(IQ)*QDIG
 40   CONTINUE
      SUMMM=SUM(151)
CCCCCCC DO 50 IQ=1,301
CCCCCCC SUM(IQ      =(SUM(IQ)*1000.0)/SUMMM
C50CCCC CONTINUE
C
C
C       OUTPUT FOR THE UNSMEARED CURVE, THEN PART OF THE GAUSSIAN
C
C
C      WRITE(6,505) 
C505   FORMAT(/" UNSMEARED SCATTERING CURVE"/
C    1   " TABULAR OUTPUT IS Q ANGSTROM-1 "/
C    2   "      )DEBYE CURVE  "/
C    3   "      )LOG DEBYE CURVE X 1000")
C      DO 510 M=1,100
C      SCATTL(M)=ALOG10(SCATT(M))
C510   LINE(M)=IFIX(SCATTL(M)*1000.0)
C      WRITE(6,550) 
550    FORMAT(" ")
C      DO 600 N=1,10
C      KOUNT=(N-1)*10
C      WRITE(6,570)(Q(KOUNT+J),J=1,10)
C570   FORMAT(1X,10F7.4)
C      WRITE(6,580)(SCATT(KOUNT+J),J=1,10)
C580   FORMAT(1X,10F7.0)
C      WRITE(6,590)(LINE(KOUNT+J),J=1,10)
C590   FORMAT(1X,10I7,/)
C600   CONTINUE
c       WRITE(6,610)WAV,WAS,ALF
c       WRITE(6,506)QQ(151),QQ(201),QQ(251)
       DO 650 II=1,3 
          IQ=101+(II*50)
          DO 620 IQP=1,301
CCCCCCC NOTE DIFFERENCE FROM CHAUVIN AND CUSACK - IQ NOT IQP FOR AA, VV
             AA=1.0/(SIG(IQ )*DSQRT(2.0*PI))
             G(IQP)=0.0
             VV=(QQ(IQ)-QQ(IQP))/SIG(IQ )
             VV=VV*VV
             IF(VV.GT.76.0) GO TO 620
             G(IQP)=AA*DEXP(-VV/2.0)
 620      CONTINUE
c         IF(II.EQ.1) WRITE(6,609)QQ (151)
c         IF(II.EQ.1) WRITE(6,609)SIG(151)
c          WRITE(6,608)G  (151)
c          IF(II.EQ.1) WRITE(6,550)
          DO 615 M=1,10
             KOUNT=(M-1)*10+ 151
c             IF (II.EQ.1) WRITE(6,607)(QQ (J+KOUNT),J=1,10)
c             IF (II.EQ.1) WRITE(6,607)(SIG(J+KOUNT),J=1,10)
c             WRITE(6,606)(G  (J+KOUNT),J=1,10)
c             IF(II.EQ.1) WRITE(6,550)
 615      CONTINUE
c          WRITE(6,550) 
 650   CONTINUE
c       WRITE(6,700) F(151),SUMMM
       DO 900 IQ=152,251
          II=IQ-151
C Replaced by I() 17/6/91 A.S.Nealis
C       SCATT(II)=SUM(IQ)
          I(PLOTCNT,II)=SUM(IQ)
 900   CONTINUE

       RETURN

 506   FORMAT(/" Q , SIGMA AND GAUSSIAN VALUES AS FOLLOWS",
     1        /" THE THREE RANGES REFER TO Q VALUES OF ",3F7.4)
 606   FORMAT(1X,10F7.3)
 607   FORMAT(1X,10F7.4)
 608   FORMAT(1X,F7.3)
 609   FORMAT(1X,F7.4)
 610      FORMAT(/" SMEARING APPLIED"/
     1            " WAVELENGTH      ",F8.4/
     1            " WAVELENGTH SPREAD ",F8.4/
     1            " BEAM DIVERGENCE      ",F8.4) 
 700   FORMAT(/" INTENSITY AT ZERO Q :  BEFORE GAUSSIAN ="
     1  ,F10.3/ "               :  AFTER GAUSSIAN  ="
     2  ,F10.3/)
       END
