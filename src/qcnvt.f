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

      SUBROUTINE QCNVT(QVAL1, QVALL, I1ST, ILAST, ERROR)
      IMPLICIT NONE
      INCLUDE 'SCTPL.COM'
      DOUBLE PRECISION QVAL1, QVALL
      INTEGER ERROR, I1ST, ILAST, iii

C If zero was entered for first point, make it's value the first in the file
      IF (QVAL1.EQ.0.) QVAL1 = Q(PLOTCNT, 1)
      ERROR = -1
      IF (Q(PLOTCNT,1).GT.QVAL1) RETURN
      ERROR = -2
      IF (Q(PLOTCNT,1).GT.QVALL) RETURN
      ERROR = -3
      IF (Q(PLOTCNT,NUMPTS(PLOTCNT)).LT.QVAL1) RETURN
      ERROR = -4
      IF (Q(PLOTCNT,NUMPTS(PLOTCNT)).LT.QVALL)THEN
         QVALL = Q(PLOTCNT,NUMPTS(PLOTCNT))
         RETURN
      ENDIF
      I1ST = 1
 1    IF (Q(PLOTCNT,I1ST).LT.QVAL1) THEN
         I1ST = I1ST + 1
	 IF(I1ST.LT.MAXPTS) GOTO 1
      ENDIF
      ILAST = 1
 2    IF (Q(PLOTCNT,ILAST).LE.QVALL) THEN
         ILAST = ILAST + 1
	 IF(ILAST.LE.NUMPTS(PLOTCNT)) GOTO 2
      ENDIF
      ILAST = ILAST - 1
      ERROR = 0

      RETURN
      END
