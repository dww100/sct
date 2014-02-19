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

      SUBROUTINE PLOTQ(I1, I2)
      IMPLICIT NONE
      INCLUDE 'SCTPL.COM'
      DOUBLE PRECISION QVPL1, QVPLL, VALUE(6)
      INTEGER ERCDE, I1, I2, NVAL, IRC
      CHARACTER MESAGE*80

 1    IF ((PTYPE(PLOTCNT).LT.3).OR.(PTYPE(PLOTCNT).GT.4)) THEN
         PRINT 100
         CALL GETFVAL(5, VALUE, NVAL, IRC)
	 IF ((NVAL.NE.2).AND.(IRC.NE.2)) THEN
	    CALL ERRMSG("Two numbers required")
	    GOTO 1
         ENDIF
	 QVPL1 = VALUE(1)
	 QVPLL = VALUE(2)
C      ELSE IF (PTYPE(PLOTCNT).EQ.3) THEN
      ELSE IF ((PTYPE(PLOTCNT).EQ.3).OR.(PTYPE(PLOTCNT).EQ.4)) THEN
C Take default Q range to plot (ie all the data)
 2       PRINT 110
         CALL GETFVAL(5, VALUE, NVAL, IRC)
	 IF (IRC.EQ.2) THEN
  	    FSTPLTPT(PLOTCNT) = 1
	    LSTPLTPT(PLOTCNT) = NUMPTS(PLOTCNT)
C Store also the first and last Q-values to be plotted
	    QVALP1(PLOTCNT) = Q(PLOTCNT, 1)
	    QVALPL(PLOTCNT) = Q(PLOTCNT, NUMPTS(PLOTCNT))
	    RETURN
	 ELSE IF (NVAL.NE.2) THEN
	    CALL ERRMSG("Two numbers required")
	    GOTO 2
         ELSE IF (NVAL.EQ.2) THEN
	    QVPL1 = VALUE(1)
	    QVPLL = VALUE(2)
         ENDIF
      ENDIF
      CALL QCNVT(QVPL1, QVPLL, I1, I2, ERCDE)
      IF (ERCDE.EQ.-1) CALL ERRMSG("Q1 < min Q in file!")
      IF (ERCDE.EQ.-2) CALL ERRMSG("Q2 < min Q in file!")
      IF (ERCDE.EQ.-3) CALL ERRMSG("Q1 > max Q in file!")
      IF (ERCDE.EQ.-4) CALL ERRMSG("Q2 > max Q in file!")
c      IF (ERCDE.EQ.-4) THEN
c         MESAGE(1:45) = "Q2 > max Q in file! Have set last fit point ="
c         MESAGE(46:64) = " last point in file"
c         CALL ERRMSG(MESAGE(1:64))
c         ERCDE = 0
c      ENDIF
      IF (QVPLL.LE.QVPL1) CALL ERRMSG(" Strange choice of values! ")
      IF (ERCDE.LT.0) GOTO 1
      FSTPLTPT(PLOTCNT) = I1
      LSTPLTPT(PLOTCNT) = I2
C Store also the first and last Q-values to be plotted
      QVALP1(PLOTCNT) = Q(PLOTCNT, I1)
      QVALPL(PLOTCNT) = Q(PLOTCNT, I2)

      RETURN
 100  FORMAT("Which Q range to plot (Q1, Q2; Q1 < Q2)? ",$)
 110  FORMAT("Which Q range to plot (Q1, Q2; Q1 < Q2) [Def=All]? ",$)
      END
