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

      SUBROUTINE LINEGT(DASH1, GAP1, DASH2, GAP2)
      IMPLICIT NONE
      DOUBLE PRECISION VALUE(6)
      INTEGER LINETYP, NVAL, IRC
      INTEGER DASH1, GAP1, DASH2, GAP2

 1    PRINT 100
      CALL GETFVAL(5, VALUE, NVAL, IRC)
      IF (IRC.EQ.2) THEN
         LINETYP = 0
         DASH1 = 0
         GAP1 = 0
         DASH2 = 0
         GAP2 = 0
         RETURN
      ENDIF
      IF(NVAL.GT.1) THEN
         CALL ERRMSG("Only one value expected")
	 GOTO 1
      ENDIF
      IF ((VALUE(1).LT.0).OR.(VALUE(1).GT.3)) THEN
         CALL ERRMSG("Value out of range")
	 GOTO 1
      ELSE
         LINETYP = IDINT(VALUE(1))
      ENDIF
      IF(LINETYP.EQ.1) THEN
	 DASH1 = 5
	 GAP1 = 5
	 DASH2 = 5
	 GAP2 = 5
      ELSE IF(LINETYP.EQ.2) THEN
	 DASH1 = 10
	 GAP1 = 10
	 DASH2 = 10
	 GAP2 = 10
      ELSE IF(LINETYP.EQ.3) THEN
         PRINT 110
	 READ(*,*) DASH1, GAP1, DASH2, GAP2
      ENDIF

      RETURN
 100  FORMAT("Line style:",/,
     *       "0=full, 1=short dash, 2=long dash [0], 3=custom ",
     *       "(see GHOST BROKEN call) [0]: ",$)
 110  FORMAT("Enter DASH1, GAP1, DASH2, GAP2:")
      END
