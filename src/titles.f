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
      SUBROUTINE TITLES(RDRMOD, TITLE)
      IMPLICIT NONE
      INCLUDE 'SCTPL.COM'
      REAL FXPOS, FXPOS1, DFXPOS, FYPOS1, FYPOS2, FYPOS3, OLDFXP, XSTART,
     *     YSTART
      COMMON /RFI/ PLOT1, PLOT2, NUMPT2
      INTEGER SYMBL, II, COUNT, LIST(6), XTITLN, YTITLN, RDRMOD, TITLEN,
     *        PLOT1, PLOT2, NUMPT2
      CHARACTER DATE*9, TIME*8, FLNAME(10)*14, TITLE*76, XAXTIT*7,
     *          YAXTIT*12
      COMMON /BOOLRF/ RFMODE
      LOGICAL RFMODE

      FXPOS = 0.1
      FXPOS1 = 0.1 + 20 * .001 * 17
      DFXPOS = 0.22
      FYPOS1 = 0.045
      FYPOS2 = .02
      FYPOS3 = .03
      TITLEN = 1
C Find last non-blank character in TITLE. Set TITLEN to this value
      DO 10 II = 76, 1, -1
         IF (TITLE(II:II).NE.' ') THEN
	    TITLEN = II
	    GOTO 1
         ENDIF
 10   CONTINUE
C Decide on axes' titles, based on plot type CURPLOT
 1    IF ((CURPLOT.GE.0).AND.(CURPLOT.LE.2)) THEN
         XAXTIT(1:7) = "Q*Q A-2"
	 XTITLN = 7
      ELSE
         XAXTIT(1:7) = "Q (A-1)"
	 XTITLN = 7
      ENDIF
      IF (CURPLOT.EQ.0) THEN
         YAXTIT(1:12) = "Ln(I(Q).Q.Q)"
	 YTITLN = 12
      ELSE IF (CURPLOT.EQ.1) THEN
         YAXTIT(1:10) = "Ln(I(Q).Q)"
	 YTITLN = 10
      ELSE IF (CURPLOT.EQ.2) THEN
         YAXTIT(1:8) = "Ln(I(Q))"
	 YTITLN = 8
      ELSE IF (CURPLOT.EQ.3) THEN
         YAXTIT(1:8) = "Ln(I(Q))"
	 YTITLN = 8
      ELSE IF (CURPLOT.EQ.4) THEN
         YAXTIT(1:4) = "I(Q)"
	 YTITLN = 4
      ELSE IF (CURPLOT.EQ.5) THEN
         YAXTIT(1:8) = "I(Q).Q.Q"
	 YTITLN = 8
      ENDIF
C Sort out which plots to plot
C Check we're NOT doing R-factors
      IF (.NOT.RFMODE) THEN
         IF (RDRMOD.EQ.1) THEN
            COUNT = 0
            DO 30 II = 1, MAXPLOT
c A.S. Nealis 17/7/91
c               IF (PTYPE(II).EQ.CURPLOT) THEN
               IF ((ISPLOT(II).EQ.1).AND.(PTYPE(II).EQ.CURPLOT)) THEN
	          COUNT = COUNT + 1
	          LIST(COUNT) = II
	       ENDIF
 30         CONTINUE
         ELSE
            COUNT = 1
	    LIST(1) = PLOTCNT
         ENDIF
      ELSE
C we ARE doing R-factors
         COUNT = 2
	 LIST(1) = PLOT1
	 LIST(2) = PLOT2
      ENDIF
      IF (COUNT.GT.6) THEN
C-SJP         CALL PICNOW
C-SJP	 CALL TOALPH
         CALL ERRMSG(" More than 6 plots ! Only showing first 6 ")
	 COUNT = 6
      ENDIF
C Draw a box around max drawing area
C-SJP      CALL PSPACE(0.0, 1.33, 0.0, 1.0)
C-SJP      CALL BORDER
C-SJP      CALL MAP(0.0, 1.33, 0.0, 1.0)

C Draw plot title
C-SJP      CALL CTRMAG(20)
C-SJP      CALL CTRORI(0.)
C-SJP      CALL PCSCEN(.665, .97,TITLE(1:TITLEN))
C Add axes' labels
C-SJP      CALL CTRORI(0.)
C-SJP      CALL PCSCEN(0.80, 0.04, XAXTIT(1:XTITLN))
C-SJP      CALL CTRORI(90.)
C-SJP      CALL PCSCEN(0.02, 0.5, YAXTIT(1:YTITLN))
C-SJP      CALL CTRORI(0.)

C Draw plot number, name and symbol
C If COUNT < 3 use bigger characters
      IF (COUNT.LT.3) THEN
C-SJP         CALL CTRMAG(20)
C-SJP	 DO 40 II = 1, COUNT
C-SJP	    IF (II.EQ.1) CALL PLOTNI(FXPOS, FYPOS3, LIST(II))
C-SJP	    IF (II.EQ.2) CALL PLOTNI(FXPOS1, FYPOS3, LIST(II))
C-SJP	    CALL TYPECS(' ')
C-SJP	    CALL TYPECS(FILENAME(LIST(II)))
C-SJP	    IF (STYLE(LIST(II)).NE.0) THEN
C-SJP	       CALL TYPECS(' ')
C-SJP	       CALL TYPENC(SYMBLIST(LIST(II)))
C-SJP            ENDIF
C-SJP 40      CONTINUE
C-SJP      ELSE
C Too many plots for a large key, so use a smaller one.
C-SJP         DO 20 II = 1, COUNT, 2
C-SJP            CALL CTRMAG(13)
C-SJP            CALL PLOTNI(FXPOS, FYPOS1, LIST(II))
C-SJP	    CALL TYPECS(' ')
C-SJP	    CALL TYPECS(FILENAME(LIST(II)))
C-SJP	    IF (STYLE(LIST(II)).NE.0) THEN
C-SJP  	       SYMBL = SYMBLIST(LIST(II))
C-SJP	       CALL TYPECS(' ')
C-SJP               CALL CTRMAG(20)
C-SJP	       CALL TYPENC(SYMBL)
C-SJP	    ENDIF
C DRAW PLOT NUM, NAME AND SYMBOL JUST BELOW
C-SJP            IF ((II + 1).LE.COUNT) THEN
C-SJP               CALL CTRMAG(13)
C-SJP               CALL PLOTNI(FXPOS, FYPOS2, LIST(II + 1))
C-SJP	       CALL TYPECS(' ')
C-SJP  	       CALL TYPECS(FILENAME(II + 1))
C-SJP	       IF (STYLE(LIST(II + 1)).NE.0) THEN
C-SJP	          SYMBL = SYMBLIST(LIST(II + 1))
C-SJP	          CALL TYPECS(' ')
C-SJP                  CALL CTRMAG(20)
C-SJP	          CALL TYPENC(SYMBL)
C-SJP	       ENDIF
C-SJP            ENDIF
C-SJP            CALL CTRMAG(13)
            FXPOS = FXPOS + DFXPOS
 20      CONTINUE
      ENDIF
C-SJP      CALL CTRMAG(10)
C-SJP      CALL PICNOW
      RETURN
      END
