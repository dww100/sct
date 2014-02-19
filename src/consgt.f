      SUBROUTINE CONSGT

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
      DOUBLE PRECISION VALUE(6)
      INTEGER IRC, NVAL

 1    PRINT 110
      CALL GETFVAL(5, VALUE, NVAL, IRC)
      IF (IRC.EQ.2) CONSTANT(PLOTCNT) = 0.
      IF (NVAL.GT.1) THEN
         CALL ERRMSG("One value expected")
	 GOTO 1
      ENDIF

      CONSTANT(PLOTCNT) = VALUE(1)

      RETURN
 110  FORMAT("Constant to apply [0]: ",$)
      END
