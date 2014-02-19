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

      SUBROUTINE RFCTXT
      IMPLICIT NONE
      INCLUDE 'SCTPL.COM'
      COMMON /RFD/ RFCTR, CONST
      COMMON /RFI/ PLOT1, PLOT2, NUMPT2
      DOUBLE PRECISION RFCTR(11), CONST(11)
      COMMON /CHARS/ DATE, TIME, FLNAME, TITLE
      CHARACTER DATE*9, TIME*8, FLNAME(10)*14, TITLE*80
      COMMON /BOOL2/ AUTORF, IGNORE
      COMMON /BOOLRF/ RFMODE
      LOGICAL AUTORF, IGNORE(512), RFMODE
      REAL A, FXPOS, DFXPOS, FYPOS1, FYPOS2, OLDFXP, XSTART, YSTART,
     *     YDIFF, XOFFST, IGNX1, IGNY1, DIGNX1, RMAG, DY, RRFCTR(11),
     *     RCONST(11)
      INTEGER II, SYMBL, NAFTPT, J, PLOT1, PLOT2, NUMPT2

C Get the date and time
C-SJP      CALL ENQDAT(DATE)
C-SJP      CALL ENQTIM(TIME)
C Copy the DOUBLE PRECISION RFCTR, CONST into REAL arrays
      DO 10 II = 1, 11
         RRFCTR(II) = REAL(RFCTR(II))
	 RCONST(II) = REAL(CONST(II))
 10   CONTINUE
C Draw a box around max drawing area
C-SJP      CALL PSPACE(0.0, 1.33, 0.0, 1.0)
C-SJP      CALL BORDER
C-SJP      CALL MAP(0.0, 1.33, 0.0, 1.0)
C add text and posit info
C-SJP      CALL CTRORI(0.0)
C-SJP      CALL CTRMAG(20)
C-SJP      CALL PLOTCS(.965, .92, 'RFACTORS')
C-SJP      CALL CTRMAG(15)
C-SJP      CALL PLOTCS(.965, .89, DATE)
C-SJP      CALL SPACE(2)
C-SJP      CALL TYPECS(TIME)
C-SJP      CALL CSPACE(0.0, 1.33, .0, 1.0)
C-SJP      CALL CTRORI(0.)
C-SJP      CALL PLOTCS(.02, .97,TITLE(1:76))
      XSTART = .98
      XOFFST = .013
      YSTART = .60
C GAP BETWEEN LINES
      DY = 0.035
C TITLES OF MODEL, ETC
C-SJP      CALL CTRMAG(20)
C-SJP      CALL PLOTCS(.965, .83, ' XPTL DATA FILE:')
C-SJP      CALL TYPECS(' ')
C-SJP      CALL PLOTNI(XSTART + XOFFST, .80, PLOT2)
C-SJP      CALL TYPECS(' ')
C-SJP      CALL TYPECS(FILENAME(PLOT2))
      SYMBL = SYMBLIST(PLOT2)
C-SJP      CALL TYPECS(' ')
C-SJP      CALL CTRMAG(20)
C-SJP      CALL TYPENC(SYMBL)
C-SJP      CALL PLOTCS(XSTART + XOFFST, .77, 'Constant = ')
C-SJP      CALL TYPENF(REAL(CONSTANT(PLOT2)), 2)
C-SJP      CALL PLOTCS(.965, .73, ' MODEL DATA FILE:')
C-SJP      CALL CTRMAG(20)
C-SJP      CALL TYPECS(' ')
C-SJP      CALL PLOTNI(XSTART + XOFFST, .70, PLOT1)
C-SJP      CALL TYPECS(' ')
C-SJP      CALL TYPECS(FILENAME(PLOT1))
C-SJP      CALL PLOTCS(XSTART + XOFFST, .67, 'Constant = ')
C-SJP      CALL TYPENF(REAL(CONSTANT(PLOT1)), 2)
C-SJP      CALL PLOTCS(.965, .63, ' CONST   RFACTOR')
      DO 30 II = 1, 11
C-SJP         CALL CTRMAG(20)
C-SJP	 CALL PLOTNF(XSTART + 3.0*XOFFST, YSTART, RCONST(II), 2)
C-SJP	 CALL SPACE(2)
C-SJP	 CALL TYPENF(RRFCTR(II), 4)
	 YSTART = YSTART - DY
 30   CONTINUE
C-SJP      CALL PLOTCS(XSTART + 3.0*XOFFST, YSTART - 0.02, 'PTS DELETED')
      IGNY1 = YSTART - .05
      DIGNX1 = .05
C-SJP      CALL CTRMAG(15)
      J = 0
      IGNX1 = XSTART - 2. * XOFFST
      DO 40 II = 1, NUMPT2
         IF (IGNORE(II)) THEN
	    J = J + 1
	    IGNX1 = IGNX1 + DIGNX1
C-SJP            CALL PLOTNI(IGNX1, IGNY1, II)
	    IF (J.EQ.6) THEN
	       J = 0
               IGNY1 = IGNY1 - .02
               IGNX1 = XSTART - 2. * XOFFST
	    ENDIF
	 ENDIF
 40   CONTINUE
C-SJP      CALL CTRMAG(10)
C-SJP      CALL PICNOW
c      CALL TOALPH
      RETURN
      END
