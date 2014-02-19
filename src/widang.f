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
      SUBROUTINE WIDANG(RDRMOD)
      IMPLICIT NONE
      INCLUDE 'SCTPL.COM'
      REAL XSTART, XOFFST, YSTART, DY, YDIFF
      INTEGER II, COUNT, LIST(6), SYMBL, RDRMOD
      CHARACTER DATE*9, TIME*8

C Get the date and time
C-SJP      CALL ENQDAT(DATE)
C-SJP      CALL ENQTIM(TIME)
c Initialise numbers for the key
      XSTART = .98
      XOFFST = .013
      YSTART = .853846
C GAP BETWEEN LINES
      DY = 0.035
C EXTRA GAP BETWEEN EACH GROUP OF DATA
      YDIFF = 0.02
C Sort out which plots to plot
      IF (RDRMOD.EQ.1) THEN
         COUNT = 0
         DO 30 II = 1, MAXPLOT
           IF ((ISPLOT(II).EQ.1).AND.(PTYPE(II).EQ.CURPLOT)) THEN
	      COUNT = COUNT + 1
	      LIST(COUNT) = II
	   ENDIF
 30     CONTINUE
      ELSE
         COUNT = 1
	 LIST(1) = PLOTCNT
      ENDIF
C Now go ahead and plot them!
      IF (COUNT.GT.6) THEN
C-SJP         CALL PICNOW
	 CALL TOALPH
         CALL ERRMSG(" More than 6 plots ! Only showing first 6 ")
	 COUNT = 6
      ENDIF
C Draw a box around max drawing area
C-SJP      CALL PSPACE(0.0, 1.33, 0.0, 1.0)
C-SJP      CALL BORDER
C-SJP      CALL MAP(0.0, 1.33, 0.0, 1.0)
C-SJP      CALL CTRMAG(20)
C-SJP      CALL PLOTCS(.965, .92, 'WIDE ANGLE')
C-SJP      CALL CTRMAG(15)
C-SJP      CALL PLOTCS(.965, .89, DATE)
C-SJP      CALL SPACE(2)
C-SJP      CALL TYPECS(TIME)
C add text and posit info
C-SJP      DO 35 II = 1, COUNT
C-SJP         CALL CTRMAG(20)
C-SJP	 CALL TYPECS(' ')
C-SJP	 CALL PLOTNI(XSTART + XOFFST, YSTART, LIST(II))
C-SJP	 CALL TYPECS(' ')
C-SJP	 CALL TYPECS(FILENAME(LIST(II)))
C-SJP	 SYMBL = SYMBLIST(LIST(II))
C-SJP	 IF(STYLE(II).NE.0) THEN
C-SJP	    CALL TYPECS(' ')
C-SJP	    CALL CTRMAG(20)
C-SJP	    CALL TYPENC(SYMBL)
C-SJP	 ENDIF
C NEXT LINE
C-SJP         CALL CTRMAG(20)
C-SJP	 YSTART = YSTART - DY
C-SJP	 CALL PLOTCS(XSTART, YSTART, 'PTS')
C-SJP	 CALL TYPENI(FSTPLTPT(LIST(II)))
C-SJP	 CALL TYPENI(LSTPLTPT(LIST(II)))
C NEXT LINE
C-SJP	 YSTART = YSTART - DY
C-SJP	 CALL PLOTCS(XSTART, YSTART, 'CONSTANT')
C-SJP         CALL TYPENF(REAL(CONSTANT(LIST(II))),2)
C EXTRA SPACES BETWEEN EACH SET
C-SJP	 YSTART = YSTART - DY - YDIFF
C-SJP 35   CONTINUE

C-SJP 1    CALL CTRMAG(10)
C-SJP      CALL PICNOW
c      CALL TOALPH

      RETURN
      END

