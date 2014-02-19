      SUBROUTINE AXTRANS
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
      INTEGER J
C
C Convert the Q and I of the current plot to the required form
C
      DO 10 J = 1, NUMPTS(PLOTCNT)
	 ERORU(PLOTCNT,J) = 0.0
	 ERORL(PLOTCNT,J) = 0.0
         IF (PTYPE(PLOTCNT).EQ.0) THEN
            I(PLOTCNT,J) = I(PLOTCNT,J)*Q(PLOTCNT,J)*Q(PLOTCNT,J)
	    I(PLOTCNT,J) = DLOG(DABS(I(PLOTCNT,J)))
	    Q(PLOTCNT,J) = Q(PLOTCNT,J)*Q(PLOTCNT,J)
         ELSE IF (PTYPE(PLOTCNT).EQ.1) THEN
	    ERORU(PLOTCNT,J) = DLOG((DABS(I(PLOTCNT,J))+ EROR(PLOTCNT,J))
     *                         * Q(PLOTCNT,J))
	    ERORU(PLOTCNT,J) = ERORU(PLOTCNT,J) + CONSTANT(PLOTCNT)
	    ERORL(PLOTCNT,J) = ERORL(PLOTCNT,J) + CONSTANT(PLOTCNT)
            I(PLOTCNT,J) = I(PLOTCNT,J)*Q(PLOTCNT,J)
	    I(PLOTCNT,J) = DLOG(DABS(I(PLOTCNT,J)))
	    Q(PLOTCNT,J) = Q(PLOTCNT,J)*Q(PLOTCNT,J)
         ELSE IF (PTYPE(PLOTCNT).EQ.2) THEN
	    ERORU(PLOTCNT,J) = DLOG(DABS(I(PLOTCNT,J) + EROR(PLOTCNT,J)))
	    ERORL(PLOTCNT,J) = 2.*DLOG(DABS(I(PLOTCNT,J))) - 
     *                         ERORU(PLOTCNT,J)
	    ERORU(PLOTCNT,J) = ERORU(PLOTCNT,J) + CONSTANT(PLOTCNT)
	    ERORL(PLOTCNT,J) = ERORL(PLOTCNT,J) + CONSTANT(PLOTCNT)
            I(PLOTCNT,J) = DLOG(DABS(I(PLOTCNT,J)))
	    Q(PLOTCNT,J) = Q(PLOTCNT,J)*Q(PLOTCNT,J)
         ELSE IF (PTYPE(PLOTCNT).EQ.3) THEN
	    ERORU(PLOTCNT,J) = DLOG(DABS(I(PLOTCNT,J) + EROR(PLOTCNT,J)))
	    ERORL(PLOTCNT,J) = 2.*DLOG(DABS(I(PLOTCNT,J))) - 
     *                         ERORU(PLOTCNT,J)
c A.S. Nealis 25/3/92
c    ERORU(PLOTCNT,J) = ERORU(PLOTCNT,J)
c    ERORL(PLOTCNT,J) = ERORL(PLOTCNT,J)
	    ERORU(PLOTCNT,J) = ERORU(PLOTCNT,J) + CONSTANT(PLOTCNT)
	    ERORL(PLOTCNT,J) = ERORL(PLOTCNT,J) + CONSTANT(PLOTCNT)
            I(PLOTCNT,J) = DLOG(DABS(I(PLOTCNT,J))) + CONSTANT(PLOTCNT)
         ELSE IF (PTYPE(PLOTCNT).EQ.5) THEN
            I(PLOTCNT,J) = I(PLOTCNT,J)*Q(PLOTCNT,J)*Q(PLOTCNT,J)
	 ENDIF
 10   CONTINUE

      IF (PTYPE(PLOTCNT).EQ.0) THEN
	 QVALF1(PLOTCNT) = QVALF1(PLOTCNT)*QVALF1(PLOTCNT)
	 QVALFL(PLOTCNT) = QVALFL(PLOTCNT)*QVALFL(PLOTCNT)
      ELSE IF (PTYPE(PLOTCNT).EQ.1) THEN
	 QVALF1(PLOTCNT) = QVALF1(PLOTCNT)*QVALF1(PLOTCNT)
	 QVALFL(PLOTCNT) = QVALFL(PLOTCNT)*QVALFL(PLOTCNT)
      ELSE IF (PTYPE(PLOTCNT).EQ.2) THEN
	 QVALF1(PLOTCNT) = QVALF1(PLOTCNT)*QVALF1(PLOTCNT)
	 QVALFL(PLOTCNT) = QVALFL(PLOTCNT)*QVALFL(PLOTCNT)
      ENDIF

      RETURN
      END
