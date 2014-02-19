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

      SUBROUTINE PSTYLE
      IMPLICIT NONE
      INCLUDE 'SCTPL.COM'
      DOUBLE PRECISION VALUE(6)
      INTEGER DEFSTY, NVAL, IRC

      IF(PLTYPE(PLOTCNT).EQ.0) DEFSTY = 0
      IF(PLTYPE(PLOTCNT).EQ.1) DEFSTY = 1

 1    PRINT 100, DEFSTY
      CALL GETFVAL(5, VALUE, NVAL, IRC)
      IF (IRC.EQ.2) THEN
         STYLE(PLOTCNT) = DEFSTY
      ELSE IF (NVAL.NE.1) THEN
         CALL ERRMSG("Only one value expected")
	 GOTO 1
      ELSE
         STYLE(PLOTCNT) = IDINT(VALUE(1))
      ENDIF

      RETURN
 100  FORMAT("Choose a plot style:",/,
     *       "Join the dots, no symbols [0]: Symbols only [1]:",
     *	     " Dots and symbols [2].",/,"Default is [",I1,"]: ",$)
      END
