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
      SUBROUTINE UNAXTR(INDEX)
      IMPLICIT NONE
      INCLUDE 'SCTPL.COM'
      INTEGER J, INDEX
C
C Undo the Q and I from Q**2, Ln(I) etc, of the
C current plot back to Q and I
C
      DO 10 J = 1, NUMPTS(INDEX)
         IF (PTYPE(INDEX).EQ.0) THEN
            I(INDEX, J) = DEXP(I(INDEX, J))/Q(INDEX, J)/Q(INDEX, J)
	    Q(INDEX, J) = DSQRT(Q(INDEX, J))
         ELSE IF (PTYPE(INDEX).EQ.1) THEN
	    Q(INDEX, J) = DSQRT(Q(INDEX, J))
            I(INDEX, J) = DEXP(I(INDEX, J))/Q(INDEX, J)
         ELSE IF (PTYPE(INDEX).EQ.2) THEN
            I(INDEX, J) = DEXP(I(INDEX, J))
	    Q(INDEX, J) = DSQRT(Q(INDEX, J))
         ELSE IF (PTYPE(INDEX).EQ.3) THEN
            I(INDEX, J) = DEXP(I(INDEX, J) - CONSTANT(INDEX))
         ELSE IF (PTYPE(INDEX).EQ.5) THEN
            I(INDEX, J) = I(INDEX, J)/Q(INDEX, J)
	    Q(INDEX, J) = DSQRT(Q(INDEX, J))
	 ENDIF
 10   CONTINUE

      RETURN
      END
