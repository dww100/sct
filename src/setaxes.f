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

C A.S. Nealis: 26/3/92
      SUBROUTINE SETAXES(RDRMOD)
C The dimensions of QBOUND and IBOUND are to be 2 * MAXPLOT
C RQMIN1, RQMAX1, RIMIN1, RIMAX1: REAL of QMIN1, QMAX1, IMIN1, IMAX1
C                               : for calls to SCALES
      IMPLICIT NONE
      INCLUDE 'SCTPL.COM'
      DOUBLE PRECISION QBOUND(20), IBOUND(20), QMAX1, QMIN1, IMAX1,
     *                  IMIN1, QQ(MAXPTS), II(MAXPTS)
C     REAL RQMIN1, RQMAX1, RIMIN1, RIMAX1
      INTEGER RDRMOD, J, K, L, M

      IF (RDRMOD.EQ.0) THEN
         K = 1
         DO 10 J = FSTPLTPT(PLOTCNT), LSTPLTPT(PLOTCNT)
	    QQ(K) = Q(PLOTCNT, J)
	    II(K) = I(PLOTCNT, J)
	    K = K + 1
 10      CONTINUE
         CALL MINMAX(QQ, QMIN1, QMAX1, K - 1)
         CALL MINMAX(II, IMIN1, IMAX1, K - 1)
      ELSE IF (RDRMOD.EQ.1) THEN
         M = 1
         DO 20 J = 1, PLOTCNT
	    IF (PTYPE(J).EQ.CURPLOT) THEN
	       L = 1
C Copy valid Q and I values of each valid plot into QQ and II, one at a time
               DO 30 K = FSTPLTPT(J), LSTPLTPT(J)
	          QQ(L) = Q(J, K)
	          II(L) = I(J, K)
		  L = L + 1
 30	       CONTINUE
C Each valid plot is sorted into max and min Q and I
	       CALL MINMAX(QQ, QBOUND(M), QBOUND(M+1), L - 1)
	       CALL MINMAX(II, IBOUND(M), IBOUND(M+1), L - 1)
	       M = M + 2
	    ENDIF
 20      CONTINUE
         CALL MINMAX(QBOUND, QMIN1, QMAX1, M-1)
         CALL MINMAX(IBOUND, IMIN1, IMAX1, M-1)
      ENDIF
      RQMIN1 = REAL(QMIN1)
      RQMAX1 = REAL(QMAX1)
      RIMIN1 = REAL(IMIN1)
      RIMAX1 = REAL(IMAX1)
      RETURN
      END
