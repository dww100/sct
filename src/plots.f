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

        SUBROUTINE PLGRAF(II)
	IMPLICIT NONE
	INTEGER II, J, K
	REAL EU(512), EL(512), const
	INCLUDE 'SCTPL.COM'

        K = 1
        CONST = CONSTANT(II)
C-SJP	CALL MARKER(SYMBLIST(II))
        DO 10 J = FSTPLTPT(II), LSTPLTPT(II)
	   QTEMP(K) = REAL(Q(II, J))
	   ITEMP(K) = REAL(I(II, J))
	   EU(K) = REAL(ERORU(II,J))
	   EL(K) = REAL(ERORL(II,J))
	   K = K + 1
 10     CONTINUE
C-SJP        IF(PTYPE(II).EQ.3.OR.PTYPE(II).EQ.2.OR.PTYPE(II).EQ.1) THEN
c        IF(PTYPE(II).EQ.3.OR.PTYPE(II).EQ.2) THEN
C-SJP           IF (STYLE(II).EQ.0) THEN
C-SJP	      CALL PTJOIN(QTEMP, ITEMP, 1, K - 1, 0)
C-SJP           ELSE IF (STYLE(II).EQ.1) THEN
C-SJP	      CALL PTPLOT(QTEMP, ITEMP, 1, K - 1, 0)
C-SJP	      CALL ERRBAR(QTEMP, ITEMP, CONST, EL, EU, 1, K - 1, 0,1)
C-SJP           ELSE IF (STYLE(II).EQ.2) THEN
C-SJP	      CALL PTGRAF(QTEMP, ITEMP, 1, K - 1, 0)
C-SJP	      CALL ERRBAR(QTEMP, ITEMP, CONST, EL, EU, 1, K - 1, 0,1)
C-SJP	   ENDIF
C-SJP        ELSE
C-SJP           IF (STYLE(II).EQ.0) CALL PTJOIN(QTEMP, ITEMP, 1, K - 1, 0)
C-SJP           IF (STYLE(II).EQ.1) CALL PTPLOT(QTEMP, ITEMP, 1, K - 1, 0)
C-SJP           IF (STYLE(II).EQ.2) CALL PTGRAF(QTEMP, ITEMP, 1, K - 1, 0)
C-SJP	ENDIF

	RETURN
	END

        SUBROUTINE ERRBAR(Q, I, CONS, EL, EU, FROM, TO, DO)
c       SUBROUTINE ERRBAR(Q, I, EL, EU, FROM, TO, DO)
C A.S. Nealis 29/7/91: GHOST ERRBAR not supplied, have to write my own
C (shame), shouldn't be difficult!
C
	IMPLICIT NONE
        REAL Q(512), I(512), EU(512), EL(512), CONS
	INTEGER FROM, TO, DO, II

C-SJP        DO 10 II = FROM, TO
C-SJP	   CALL POSITN(Q(II), 2 * I(II) - EU(II))
C-SJP	   CALL JOIN(Q(II), EU(II))
 10     CONTINUE

        RETURN
	END
