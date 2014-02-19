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

cA.S. Nealis 19/11/91: Revised to allow point selection after 1st run through
c ie if CHOICE = 7
      SUBROUTINE FITQ(I1, I2, CHOICE)
      IMPLICIT NONE
      INCLUDE 'SCTPL.COM'
      DOUBLE PRECISION QVPL1, QVPLL, QVFT1, QVFTL
      INTEGER ERCDE, I1, I2, I3, I4, I5, CHOICE
      CHARACTER MESAGE * 80

 2    IF (CHOICE.EQ.7) THEN
         PRINT 110
	 READ(*,*) I3, I4
	 IF (I3.LE.0) I3 = 1
	 IF (I4.LE.0) I4 = 1
	 IF (IABS(I3 - I4).LT.2) THEN
	    CALL ERRMSG("Q1 and Q2 too close together")
	    GOTO 2
         ENDIF
         IF (I4.LT.I3) THEN
	    I5 = I3
	    I3 = I4
	    I4 = I5
         ENDIF
	 IF (I4.GT.NUMPTS(PLOTCNT)) I4 = NUMPTS(PLOTCNT)
      ELSE
 1       PRINT 100
         READ(*,*) QVFT1, QVFTL
         CALL QCNVT(QVFT1, QVFTL, I3, I4, ERCDE)
         IF (ERCDE.EQ.-1) CALL ERRMSG("Q3 < min Q in file!")
         IF (ERCDE.EQ.-2) CALL ERRMSG("Q4 > min Q in file!")
         IF (ERCDE.EQ.-3) CALL ERRMSG("Q3 < max Q in file!")
         IF (ERCDE.EQ.-4) THEN
            MESAGE(1:43) = "Q4 > max Q in file! Have set last fit point"
            MESAGE(43:61) = " = last plot point"
            CALL ERRMSG(MESAGE(1:64))
	    I4 = LSTPLTPT(PLOTCNT)
            ERCDE = 0
         ENDIF
         IF (ERCDE.LT.0) GOTO 1
         IF (I4.LE.I3) THEN
            CALL ERRMSG("Strange choice of values!")
	    GOTO 1
         ELSE IF ((I3.LT.I1).OR.(I4.GT.I2)) THEN
            CALL ERRMSG("Fit range and Plot range overlap")
	    GOTO 1
         ENDIF
      ENDIF
      FSTFITPT(PLOTCNT) = I3
      LSTFITPT(PLOTCNT) = I4
C Store also the first and last Q-values to be plotted
      QVALF1(PLOTCNT) = Q(PLOTCNT, I3)
      QVALFL(PLOTCNT) = Q(PLOTCNT, I4)

      RETURN
 100  FORMAT("Which Q range to fit a line to (Q3, Q4; Q3 < Q4)? ",$)
 110  FORMAT("Which points to fit between (P1, P2; P2 < P1)? ",$)
      END
