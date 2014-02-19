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
      SUBROUTINE GETFILE(MODE)
      IMPLICIT NONE
      INCLUDE 'SCTPL.COM'
      DOUBLE PRECISION POINT(1536), VALUE(6)
      REAL RINC, COUNT
      INTEGER NPTS, PCNT, N, LUN, IRC, NVAL, II, MODE, MAXPNT
      LOGICAL HAMBURG

      N = 1
      PCNT = PLOTCNT
      LUN = 1
      COUNT = 0.
      HAMBURG = .FALSE.
      IF (MODE.EQ.0) MAXPNT = MAXCRV * MAXPTS
      IF (MODE.EQ.1) MAXPNT = MAXPTS
      CALL OPNFIL(LUN, IRC, FILENAME, PLOTCNT)
 
 1    IF (INT(COUNT).LE.MAXPNT) THEN
         CALL GETFVAL(LUN,VALUE,NVAL,IRC)
         IF (NVAL.EQ.1) RINC = 1./3.
         IF (NVAL.EQ.2) RINC = 1.
         IF (NVAL.EQ.3) RINC = 1.
         IF (NVAL.EQ.4) RINC = 2.
         IF (NVAL.EQ.6) RINC = 2.
         IF ((INT(COUNT + RINC).GT.MAXPNT).AND.(NVAL.GT.3)) THEN
            IF (INT(COUNT + RINC - 1).LT.MAXPNT) THEN
               IF (NVAL.EQ.4) NVAL = 2
	       IF (NVAL.EQ.6) NVAL = 3
	    ENDIF
         ENDIF
         IF (IRC.EQ.0) THEN
	    N = INT(COUNT + 1)
            IF (NVAL.EQ.1) THEN
               HAMBURG = .TRUE.
	       POINT(N) = VALUE(1)
	       COUNT = COUNT + RINC
               GOTO 1
            ELSE IF(NVAL.EQ.2) THEN
		Q(PCNT,N) = VALUE(1)
		I(PCNT,N) = VALUE(2)
		COUNT = COUNT + RINC
		GOTO 1
	     ELSE IF(NVAL.EQ.3) THEN
		Q(PCNT,N) = VALUE(1)
		I(PCNT,N) = VALUE(2)
		EROR(PCNT,N) = VALUE(3)
		COUNT = COUNT + RINC
		GOTO 1
	     ELSE IF(NVAL.EQ.4) THEN
		Q(PCNT,N) = VALUE(1)
		I(PCNT,N) = VALUE(2)
		Q(PCNT,N+1) = VALUE(4)
		I(PCNT,N+1) = VALUE(5)
		COUNT = COUNT + RINC
		GOTO 1
	     ELSE IF(NVAL.EQ.6) THEN
		Q(PCNT,N) = VALUE(1)
		I(PCNT,N) = VALUE(2)
		EROR(PCNT,N) = VALUE(3)
		Q(PCNT,N+1) = VALUE(4)
		I(PCNT,N+1) = VALUE(5)
		EROR(PCNT,N+1) = VALUE(6)
		COUNT = COUNT + RINC
		GOTO 1
	    ENDIF
         ENDIF
      ENDIF
      NPTS = INT(COUNT)
      NUMPTS(PLOTCNT) = NPTS
      IF(HAMBURG) THEN
         DO 20 II = 1, NPTS, 3
	    Q(PCNT,NPTS) = POINT(II)
	    I(PCNT,NPTS) = POINT(II + 1)
	    EROR(PCNT,NPTS) = POINT(II + 2)
 20      CONTINUE
      ENDIF

      CLOSE(LUN)

      RETURN
      END
