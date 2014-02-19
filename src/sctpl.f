      PROGRAM SCTPL4
C
C Minor update to the SGI version of SCTPL
C Written by A.S. NEALIS 6 JUNE 1991
C Minor refactoring and editing to compile
C using gfortran by David W. Wright 2013
C
C ENHANCEMENTS FROM SCTPL3: 
C INCREASE MAX NUMBER OF POINTS TO 512
C	MORE THAN ONE SET OF PLOTS PER GRAPH
C	ENTER MODEL AND EXPTAL DATA SEPARATELY
C

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

      IMPLICIT NONE
      INCLUDE 'SCTPL.COM'
      COMMON /LINEFIT/ R, W, SI, SIGSI, SL, SIGSL, I0, DI0, RG, DRG
      DOUBLE PRECISION R(6), W(7), SI, SIGSI, SL, SIGSL, I0(10),
     *                 DI0(10), RG(10), DRG(10), WEIGHTED(7680), BTM(25)
      INTEGER MODE, II, NUMPLTS, CURVES, NUMQ, LUN, ERORCD, IRC,
     *        PLSTYL, DASH1, GAP1, DASH2, GAP2, I1, I2, J, RDRMOD,
     *        CHOICE
      COMMON /BOOLRF/ RFMODE
      LOGICAL AXMODE, RFMODE
      CHARACTER * 76 PLTITL

C A.S. Nealis 25 Aug 1992 - Add a facility to log the RG, RXS plot info

      INTEGER LOGFIL
      CHARACTER LOGNAM*25, PLONAM*3

C A.S. Nealis 25 Aug 1992 - Add a facility to log the RG, RXS plot info

      LOGFIL = 3
      LOGNAM = "Sctpl6_Summary"

C Initialise the contents of SCTPL.COM
C SJP 21-Jan-2000 - stripped out all graphics calls so that it only 
C generates text Logged information

 2    CALL INIT
      RFMODE = .FALSE.

      PRINT 100
      LUN = 1
      RDRMOD = 0
      AXMODE = .TRUE.
      CHOICE = -1
C MODE MUST be 0!! else LINFIT will freak out
      PLOTCNT = 1
      MODE = 0
 1    IF(PLOTCNT.GT.MAXPLOT) THEN
        PRINT *, "MORE THAN 10 PLOTS"
	      STOP
      ENDIF

      CALL DATYPE
      IF (CHOICE.EQ.-1) CALL TITLGT(PLTITL)

      IF (PLTYPE(PLOTCNT).EQ.0) THEN
         CALL MPLTGT(CURVES, BTM)
         CALL MDATGT(CURVES, ERORCD, NUMQ, WEIGHTED)
         CALL WGHTDAV(WEIGHTED, BTM, NUMQ, CURVES)
	       CALL ISMEAR(NUMQ)
      ELSE IF (PLTYPE(PLOTCNT).EQ.1) THEN
         CALL GETFILE(1)
      ENDIF

 3    IF ((CHOICE.EQ.-1).OR.(CHOICE.EQ.1).OR.(CHOICE.EQ.7)) THEN
         CALL PTYPGT
      ENDIF

      IF ((CHOICE.EQ.2).AND.(PTYPE(PLOTCNT).EQ.3)) CALL CONSGT

      CALL PLOTQ(I1, I2)

      IF (ISLINFIT(PLOTCNT).EQ.1) CALL FITQ(I1, I2, CHOICE)

      CALL PSTYLE

      IF((PLSTYL.EQ.1).OR.(PLSTYL.EQ.2)) CALL SYMBGT
      IF((PLSTYL.EQ.0).OR.(PLSTYL.EQ.2)) CALL LINEGT(DASH1, GAP1,
     *                                               DASH2, GAP2)
C Convert the data according to the axis types specified
C
      CALL AXTRANS
C If line fitting is required
C
      IF (ISLINFIT(PLOTCNT).EQ.1) THEN
         CALL LINFIT(MODE, SI, SIGSI, SL, SIGSL, W(3))
C Calculate RG/RXS and error on RG/RXS from LINFIT answers
         R(6) = DSQRT(DBLE(1 + PTYPE(PLOTCNT)) * DABS(SL))
         R(5) = 0.5 * SIGSL * DSQRT(DBLE(1 + PTYPE(PLOTCNT)) / DABS(SL))
	       RG(PLOTCNT) = R(6)
	       DRG(PLOTCNT) = R(5)
C Calculate I(0)/I(Q).Q and error on I(0)/I(Q).Q from LINFIT answers
         W(6) = DEXP(SI)
         W(5) = W(6) * SIGSI
	       I0(PLOTCNT) = W(6)
	       DI0(PLOTCNT) = W(5)
C Calculate the I value at the 'start' and 'end' of the fitted line
         IVALF1(PLOTCNT) = SL * QVALF1(PLOTCNT) + SI
         IVALFL(PLOTCNT) = SL * QVALFL(PLOTCNT) + SI
C Calculate Q.RG Range
         R(3) = R(6) * DSQRT( Q( PLOTCNT, FSTFITPT(PLOTCNT) ) )
         R(4) = R(6) * DSQRT( Q( PLOTCNT, LSTFITPT(PLOTCNT) ) )
      ENDIF

C Do the first plot
C-SJP      CALL PAPER(1)
      CALL GRINIT
      IF(CHOICE.EQ.-1) CURPLOT = PTYPE(PLOTCNT)
C If we are starting with a W.A. plot, then use default axes
      IF ((PTYPE(PLOTCNT).EQ.3).OR.(PTYPE(PLOTCNT).EQ.4)) THEN
C A.S. Nealis 26/3/92: changed so axes are preserved between plots
C        CALL MAP(0.0, 0.16, 1., 9.)
C-SJP         CALL MAP(RQMIN1, RQMAX1, RIMIN1, RIMAX1)
C-SJP	 CALL SCALES
	 IF (CHOICE.EQ.2) THEN
	    DO 20 II = 1, PLOTCNT
	       IF (PTYPE(II).EQ.CURPLOT) CALL PLGRAF(II)
 20         CONTINUE
         ELSE
	    CALL PLGRAF(PLOTCNT)
         ENDIF
C-SJP	 CALL PICNOW
         CALL TITLES(RDRMOD, PLTITL)
         CALL WIDANG(RDRMOD)
C-SJP	 CALL PICNOW
	 CALL TOALPH
      ELSE
C-SJP	 CALL BROKEN(DASH1, GAP1, DASH2, GAP2)
C A.S. Nealis: 26/3/92
C CALL SETAXES(RDRMOD, RQMIN1, RQMAX1, RIMIN1, RIMAX1)
	 CALL SETAXES(RDRMOD)
C-SJP         CALL MAP(0.0, RQMAX1, RIMIN1, RIMAX1)
C-SJP	 CALL SCALES
C-SJP	 CALL PICNOW
	 IF (CHOICE.EQ.2) THEN
	    DO 30 II = 1, PLOTCNT
	       IF (PTYPE(II).EQ.CURPLOT) CALL PLGRAF(II)
               IF (ISLINFIT(PLOTCNT).EQ.1) CALL DRWLIN(II)
 30         CONTINUE
         ELSE
	    CALL PLGRAF(PLOTCNT)
            IF (ISLINFIT(PLOTCNT).EQ.1) CALL DRWLIN(PLOTCNT)
         ENDIF
	 CALL GUINLAB(RDRMOD)
         CALL TITLES(RDRMOD, PLTITL)
C-SJP	 CALL PICNOW
	 CALL TOALPH
      ENDIF

C A.S. Nealis 25 Aug 1992 - Add a facility to log the RG, RXS plot info
      OPEN(LOGFIL, FILE=LOGNAM, ACCESS="APPEND", STATUS="UNKNOWN")
      IF(CURPLOT.EQ.0) PLONAM = "RTH"
      IF(CURPLOT.EQ.1) PLONAM = "RXS"
      IF(CURPLOT.EQ.2) PLONAM = "RG "
      IF (CURPLOT.GE.0.OR.CURPLOT.LE.2) WRITE(LOGFIL, 110)
     *     FILENAME(PLOTCNT), PLONAM, PLTITL,
     *     DSQRT(Q(PLOTCNT,FSTFITPT(PLOTCNT))),
     *     DSQRT(Q(PLOTCNT,LSTFITPT(PLOTCNT))),
     *     FSTFITPT(PLOTCNT), LSTFITPT(PLOTCNT),
     *     RG(PLOTCNT), DRG(PLOTCNT), I0(PLOTCNT), DI0(PLOTCNT)
c write to the screen (for debugging only)
c      IF (CURPLOT.GE.0.OR.CURPLOT.LE.2) WRITE(*, 110)
c     *     FILENAME(PLOTCNT), PLONAM, PLTITL,
c     *     DSQRT(Q(PLOTCNT,FSTFITPT(PLOTCNT))),
c     *     DSQRT(Q(PLOTCNT,LSTFITPT(PLOTCNT))),
c     *     FSTFITPT(PLOTCNT), LSTFITPT(PLOTCNT),
c     *     RG(PLOTCNT), DRG(PLOTCNT), I0(PLOTCNT), DI0(PLOTCNT)
      CLOSE(LOGFIL)

C Call the Decision Driverrr... (you'll love it, it's a Way Of Life)
C (decide.f)
      CALL DECIDE(CHOICE, RDRMOD, PLTITL)
      IF (CHOICE.EQ.1) GOTO 1
      IF (CHOICE.EQ.2) GOTO 1
      IF (CHOICE.EQ.3) GOTO 1
      IF (CHOICE.EQ.4) GOTO 2
      IF (CHOICE.EQ.7) GOTO 3

      STOP
C-SJP 100  FORMAT(/,/,/,/,/,/,/,/,/,/,/,/,/,/,/,/,/,/,/,/,/,

 100  FORMAT(/,
     *       "sctpl6g: SG version 4.11 21-Jan-2000",/)
 110  FORMAT(A14,",",A3,",",A76,",",F5.4,"->",F5.4,",",I3,"->",I3,",",
     *       F9.5,",",F9.5,",",E12.6",",
     *       E12.6)

      END

      SUBROUTINE TITLGT(PLTITL)
      IMPLICIT NONE
      INTEGER II
      CHARACTER * 76 PLTITL

      DO 10 II = 1, 76
         PLTITL(II:II) = ' '
 10   CONTINUE
      PRINT 100
      READ(*,'(A)') PLTITL

      RETURN
 100  FORMAT("Enter title of plot:")
      END

      SUBROUTINE DRWLIN(II)
      IMPLICIT NONE
      INCLUDE 'SCTPL.COM'
      REAL Q1, Q2, I1, I2
      INTEGER II

      Q1 = REAL(QVALF1(II))
      Q2 = REAL(QVALFL(II))
      I1 = REAL(IVALF1(II))
      I2 = REAL(IVALFL(II))
C-SJP      CALL BROKEN(1,1,1,1)
C-SJP      CALL POSITN(Q1, I1)
C-SJP      CALL JOIN(Q2, I2)

      RETURN

      END
