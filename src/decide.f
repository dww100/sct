      SUBROUTINE DECIDE(CHOICE, RDRMOD, PLTITL)
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
     *                 DI0(10), RG(10), DRG(10)
      COMMON /RFD/ RFCTR, CONST
      COMMON /RFI/ PLOT1, PLOT2, NUMPT2
      COMMON /BOOL2/ AUTORF, IGNORE
      COMMON /BOOLRF/ RFMODE
      LOGICAL AUTORF, IGNORE(MAXPTS)
      INTEGER PLOT1, PLOT2
      DOUBLE PRECISION RFCTR(11), CONST(11)
C     COMMON /PLTONE/ RQMIN1, RQMAX1, RIMIN1, RIMAX1
      COMMON /BOOL1/ AXMODE
      LOGICAL AXMODE, RFMODE
C     REAL RQMIN1, RQMAX1, RIMIN1, RIMAX1
C A.S Nealis 25/3/92: Added next statement to iron out a bug in errors
C after a call to RFACTOR()
      DOUBLE PRECISION CONS, EUPT2(MAXPTS), ELPT2(MAXPTS)
      DOUBLE PRECISION WEIGHTED(7680), BTM(25), QPT2(MAXPTS)
      DOUBLE PRECISION IPT2(MAXPTS)
      DOUBLE PRECISION VALUE(6)
      INTEGER RDRMOD, DASH1, GAP1, DASH2, GAP2, NVAL, IRC, CHOICE
      INTEGER MODE, II, NUMPLTS, CURVES, NUMQ, LUN, ERORCD,
     *        PLSTYL, I1, I2, J, K, NUMPT2, FPLTPT2, LPLTPT2, L,
C A.S Nealis 25/3/92: Added next INTEGER to iron out a bug in errors
C after a call to RFACTOR
     *        OPLTCNT
      CHARACTER ANSWER * 1, PLTITL *76

      CURPLOT = PTYPE(PLOTCNT)
      ANSWER(1:1) = ' '
 1    CALL TOALPH
      PRINT 100
      CALL GETFVAL(5, VALUE, NVAL, IRC)
      IF (IRC.EQ.2) THEN
         CALL ERRMSG(" Numeric input expected ")
	 GOTO 1
      ENDIF
      IF (NVAL.GT.1) THEN
         CALL ERRMSG(" Only one value please ")
	 GOTO 1
      ENDIF
      CHOICE = IDINT(VALUE(1))
      IF ((CHOICE.LT.1).OR.(CHOICE.GT.12)) THEN
         CALL ERRMSG(" Input out of range ")
	 GOTO 1
      ENDIF
      IF (CHOICE.EQ.1) THEN
         PLOTCNT = PLOTCNT + 1
	 RDRMOD = 0
	 AXMODE = .TRUE.
	 RETURN
      ELSE IF (CHOICE.EQ.2) THEN
         CURPLOT = PTYPE(PLOTCNT)
         PLOTCNT = PLOTCNT + 1
         PTYPE(PLOTCNT) = CURPLOT
	 IF((CURPLOT.GE.0).AND.(CURPLOT.LE.2)) ISLINFIT(PLOTCNT) = 1
	 RDRMOD = 1
	 AXMODE = .FALSE.
	 RETURN
      ELSE IF (CHOICE.EQ.3) THEN
 2       PRINT 110
	 CALL GETFVAL(5, VALUE, NVAL, IRC)
	 IF (IRC.EQ.2) THEN
	    CALL ERRMSG(" Numeric input expected ")
	    GOTO 2
         ENDIF
	 IF (NVAL.NE.4) CALL ERRMSG(" 4 values please ")
	 IF (VALUE(2).LE.VALUE(1)) THEN
	    CALL ERRMSG(" Qmax <= Qmin: try again! ")
	    GOTO 2
         ENDIF
	 IF (VALUE(4).LE.VALUE(3)) THEN
	    CALL ERRMSG(" Imax <= Imin: try again! ")
	    GOTO 2
         ENDIF
	 RQMIN1 = VALUE(1)
	 RQMAX1 = VALUE(2)
	 RIMIN1 = VALUE(3)
	 RIMAX1 = VALUE(4)
C Prompt for Qmin, Qmax, Imin, Imax
         AXMODE = .FALSE.
	 CALL GRINIT
	 IF(RDRMOD.EQ.0) THEN
	    CALL PLOTONE(PLOTCNT, RDRMOD)
	 ELSE IF(RDRMOD.EQ.1) THEN
C A.S. Nealis: 26/3/92
C    CALL PLOTSOME(AXMODE, RDRMOD, RQMIN1, RQMAX1, RIMIN1, RIMAX1)
	    CALL PLOTSOME(AXMODE)
	 ENDIF
	 CALL TITLES(RDRMOD, PLTITL)
	 IF (CURPLOT.LT.3) CALL GUINLAB(RDRMOD)
	 IF (CURPLOT.GT.2) CALL WIDANG(RDRMOD)
C-SJP	 CALL PICNOW
	 CALL TOALPH
	 AXMODE = .TRUE.
	 GOTO 1
      ELSE IF (CHOICE.EQ.4) THEN
         CALL INIT
	 RDRMOD = 0
	 PLOTCNT = 0
	 RETURN
      ELSE IF ((CHOICE.EQ.5).OR.(CHOICE.EQ.6)) THEN
C Get the plotnumbers to be compared, one must be model, the other data
         RDRMOD = 1
	 AXMODE = .FALSE.
	 RFMODE = .TRUE.
 3       PRINT 120
	 CALL GETFVAL(5, VALUE, NVAL, IRC)
	 IF (IRC.EQ.2) THEN
	    CALL ERRMSG(" Numeric input expected ")
	    GOTO 3
         ENDIF
	 IF (NVAL.NE.2) THEN
	    CALL ERRMSG(" Two plot numbers are required ")
	    GOTO 3
         ENDIF
	 PLOT1 = IDINT(VALUE(1))
	 PLOT2 = IDINT(VALUE(2))
C Make sure PLOT1 corresponds to a model, W.A. plot
	 IF ((PLTYPE(PLOT1).NE.0).OR.(PTYPE(PLOT1).NE.3)) THEN
	    CALL ERRMSG(" Wrong plot type for first plot. Try again ")
	    GOTO 3
         ENDIF
C Make sure PLOT2 corresponds to a model, W.A. plot
	 IF ((PLTYPE(PLOT2).NE.1).OR.(PTYPE(PLOT2).NE.3)) THEN
	    CALL ERRMSG(" Wrong plot type for second plot. Try again ")
	    GOTO 3
         ENDIF
         IF (CHOICE.EQ.5) AUTORF = .FALSE.
         IF (CHOICE.EQ.6) AUTORF = .TRUE.
c CALL RFACTOR
	 CALL RFACTOR(CONS)
C Now we know what to fit, let's cut out the 'crap', ie plot only the
C 'not ignored' xptl (PLOT2) data.
C First copy Q(PLOT2,) & I(PLOT2,) into QPT2() & IPT2(). Then cut out
C ignored points, and modify NUMPTS(PLOT2) accordingly.
C
C A.S. Nealis 25/3/92: Fix error bar bug after a call to RFACTOR
         OPLTCNT = PLOTCNT
         PLOTCNT = PLOT2
         CALL UNAXTR(PLOT2)
         CONSTANT(PLOT2) = CONS
         CALL AXTRANS
         PLOTCNT = OPLTCNT

         DO 10 II = 1, MAXPTS
	    QPT2(II) = Q(PLOT2, II)
	    IPT2(II) = I(PLOT2, II)
            EUPT2(II) = ERORU(PLOT2, II)
            ELPT2(II) = ERORL(PLOT2, II)
 10      CONTINUE
         L = 0
         K = NUMPTS(PLOT2)
         DO 20 II = 1, NUMPTS(PLOT2)
	    IF (.NOT.IGNORE(II)) THEN
               L = L + 1
	       Q(PLOT2, L) = Q(PLOT2, II)
	       I(PLOT2, L) = I(PLOT2, II)
               ERORU(PLOT2, l) = ERORU(PLOT2, II)
               ERORL(PLOT2, l) = ERORL(PLOT2, II)
	    ELSE
	        K = K - 1
	    ENDIF
 20      CONTINUE
         NUMPT2 = NUMPTS(PLOT2)
	 NUMPTS(PLOT2) = K
	 FPLTPT2 = FSTPLTPT(PLOT2)
	 FSTPLTPT(PLOT2) = 1
	 LPLTPT2 = LSTPLTPT(PLOT2)
	 LSTPLTPT(PLOT2) = K

	 CALL GRINIT
C A.S. Nealis 26/3/92: Generalise so that can set defaults for axes
C that stay!
C        CALL MAP(0.0, 0.16, 1., 9.)
C-SJP         CALL MAP(RQMIN1, RQMAX1, RIMIN1, RIMAX1)
C-SJP	 CALL SCALES
C-SJP	 CALL MARKER(SYMBLIST(PLOT2))
	 CALL PLGRAF(PLOT1)
	 CALL PLGRAF(PLOT2)
	 CALL TITLES(RDRMOD, PLTITL)
	 CALL RFCTXT
C-SJP	 CALL PICNOW
C Restore the data from QPT2(), IPT2() to Q(), I()
         DO 30 II = 1, MAXPTS
            Q(PLOT2, II) = QPT2(II)
            I(PLOT2, II) = IPT2(II)
            ERORU(PLOT2, II) = EUPT2(II) 
            ERORL(PLOT2, II) = ELPT2(II) 
 30      CONTINUE
         NUMPTS(PLOT2) = NUMPT2
	 FSTPLTPT(PLOT2) = FPLTPT2
	 LSTPLTPT(PLOT2) = LPLTPT2
	 RFMODE = .FALSE.
         GOTO 1
      ELSE IF (CHOICE.EQ.7) THEN
         PRINT 130
         CALL GETFVAL(5, VALUE, NVAL, IRC)
         IF (IRC.EQ.2) THEN
            CALL ERRMSG(" Numeric input expected ")
	    GOTO 1
         ENDIF
         IF (NVAL.GT.1) THEN
            CALL ERRMSG(" Only one value please ")
	    GOTO 1
         ENDIF
         PLOTNUM = IDINT(VALUE(1))
         IF ((PLOTNUM.LT.1).OR.(PLOTNUM.GT.PLOTCNT)) THEN
            CALL ERRMSG(" Input out of range ")
	    GOTO 1
         ENDIF
	 IF (SMEARD(IDINT(VALUE(1))).EQ.0) THEN
    	    CALL UNAXTR(IDINT(VALUE(1)))
         ELSE
	    CALL ERRMSG(" This data has been smeared! Can't unsmear ")
	    PRINT 140
	    READ(*,'(A)') ANSWER(1:1)
	    IF ((ANSWER(1:1).EQ.'Y').OR.(ANSWER(1:1).EQ.'y')) THEN
	       CALL UNAXTR(IDINT(VALUE(1)))
            ELSE
	       GOTO 1
            ENDIF
         ENDIF
	 RETURN
      ELSE IF (CHOICE.EQ.8) THEN
         CALL LOTUS
	 GOTO 1
      ELSE IF (CHOICE.EQ.9) THEN
         CALL PRTSCR
	 GOTO 1
      ELSE IF (CHOICE.EQ.10) THEN
C Get number of plot to redisplay
         PRINT 160
         CALL GETFVAL(5, VALUE, NVAL, IRC)
         IF (IRC.EQ.2) THEN
            CALL ERRMSG(" Numeric input expected ")
	    GOTO 1
         ENDIF
         IF (NVAL.GT.1) THEN
            CALL ERRMSG(" Only one value please ")
	    GOTO 1
         ENDIF
         PLOTNUM = IDINT(VALUE(1))
         IF ((PLOTNUM.LT.1).OR.(PLOTNUM.GT.PLOTCNT)) THEN
            CALL ERRMSG(" Input out of range ")
	    GOTO 1
         ENDIF
         CURPLOT = PTYPE(PLOTNUM)
	 PRINT 170
	 READ(5,'(A)') ANSWER(1:1)
	 CALL GRINIT
	 IF((ANSWER(1:1).EQ.' ').OR.(ANSWER(1:1).EQ.'y').OR.
     *        (ANSWER(1:1).EQ.'Y')) THEN
	    RDRMOD = 1
	    AXMODE = .TRUE.
C A.S. Nealis: 26/3/92
C  	    CALL PLOTSOME(AXMODE, RDRMOD, RQMIN1, RQMAX1,
C    *                   RIMIN1, RIMAX1)
   	    CALL PLOTSOME(AXMODE)
	 ELSE
	    RDRMOD = 0
	    AXMODE = .TRUE.
	    CALL PLOTONE(PLOTNUM, RDRMOD)
	 ENDIF
	 CALL TITLES(RDRMOD, PLTITL)
	 IF(CURPLOT.EQ.3) CALL WIDANG(RDRMOD)
	 IF((CURPLOT.GE.0).AND.(CURPLOT.LE.2)) CALL GUINLAB(RDRMOD)
c-SJP	 CALL PICNOW
	 GOTO 1
      ELSE IF (CHOICE.EQ.11) THEN
C Get number of plot to remove/undisplay
         PRINT 150
         CALL GETFVAL(5, VALUE, NVAL, IRC)
         IF (IRC.EQ.2) THEN
            CALL ERRMSG(" Numeric input expected ")
	    GOTO 1
         ENDIF
         IF (NVAL.GT.1) THEN
            CALL ERRMSG(" Only one value please ")
	    GOTO 1
         ENDIF
         PLOTNUM = IDINT(VALUE(1))
         IF ((PLOTNUM.LT.1).OR.(PLOTNUM.GT.PLOTCNT)) THEN
            CALL ERRMSG(" Input out of range ")
	    GOTO 1
         ENDIF
         ISPLOT(PLOTNUM) = 0
	 RDRMOD = 1
	 AXMODE = .TRUE.
	 IF (CURPLOT.EQ.PTYPE(PLOTNUM)) THEN
            CALL GRINIT
	    IF (RDRMOD.EQ.0) THEN
  	       CALL PLOTONE(PLOTNUM, RDRMOD)
	    ELSE IF (RDRMOD.EQ.1) THEN
C A.S. Nealis: 26/3/92
C  	       CALL PLOTSOME(AXMODE, RDRMOD, RQMIN1, RQMAX1, RIMIN1,
C    *                       RIMAX1)
	       CALL PLOTSOME(AXMODE)
	    ENDIF
C-SJP	    CALL PICNOW
	    CALL TOALPH
	 ENDIF
	 GOTO 1
      ELSE IF (CHOICE.EQ.12) THEN
C-SJP         CALL GREND
	 STOP
      ENDIF
      RETURN
 100  FORMAT("Take your pick:",/,
     *     "1) Get another;",/,
     *     "2) Add another;",/,
     *     "3) Change axes' ranges;",/,
     *     "4) Start again;",/,
     *     "5) Manual R-factor;",/,
     *     "6) Auto R-factor;",/,
     *     "7) Change a given plot's attributes;",/,
     *     "8) Save selected plots as ASCII;",/,
     *     "9) Screen dump;",/,
     *     "10) Redisplay a given plot;",/,
     *     "11) Remove a plot;",/,
     *     "12) Quit: ",$)
 110  FORMAT("Enter new values for Qmin, Qmax, Imin, Imax.")
 120  FORMAT("Enter plot sequence number for model, then xptal file: ",
     *     /,$)
 130  FORMAT("Enter number of plot to modify ",$)
 140  FORMAT("Do you want to continue regardless? ",$)
 150  FORMAT("Enter number of plot to undisplay ",$)
 160  FORMAT("Which plot to redisplay? ",$)
 170  FORMAT("Do you want to see ALL plots of this type? [y] ",$)
      END
