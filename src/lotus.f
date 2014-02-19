C***********************************************************************
C A.S.NEALIS - 10 MAY '88
C ROUTINE TO CONVERT DATA AS CREATED BY SCTCOM OR IN EITHER DARESBURY OR
C GRENOBLE FORMAT, AND OUTPUT AS FIVE SEPARATE, ADJACENT COLUMNS
C FOR LOTUS 1-2-3 READIN.
C
C  OUTPUT IS OF FORM
C
C      WRITE(10,100) QEX,IEX,ERR,Q,SCATT
C 100  FORMAT(3E13.6,F10.7,E20.14)
C
C A.S.NEALIS - 12 MAY '88 - ADAPTED FROM LOTUS1, WRITTEN 10 MAY '88
C
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
C
c      SUBROUTINE LOTUS(CONS,NEX,NXGUI,NXPT,NG,X1,X2,Y1,Y2,SIMQ,SIMI,NCR)
      SUBROUTINE LOTUS
      IMPLICIT NONE
      INCLUDE 'SCTPL.COM'
      INTEGER PL2WRI(MAXPLOT), J, K, L, NUMPTMAX, LUN
      CHARACTER * 80 FLNAME

      LUN = 2
      PRINT 100
      READ(*,'(A)') FLNAME
      OPEN(UNIT=LUN,FILE=FLNAME,ACCESS='SEQUENTIAL',STATUS='NEW')

      J = 1
      PL2WRI(1) = 1
 1    PRINT 110
      IF ((J.LE.PLOTCNT).AND.(PL2WRI(J).GE.0)) THEN
         PRINT 120
	 READ(*,*) PL2WRI(J)
	 J = J + 1
	 GOTO 1
      ENDIF

C Find which plot has most points in it
      NUMPTMAX = NUMPTS(PL2WRI(1))
      DO 10 K = 1, J - 1
         IF (NUMPTS(PL2WRI(K)).GT.NUMPTMAX) NUMPTMAX = NUMPTS(PL2WRI(K))
 10   CONTINUE

C Write the xptal data, then the models
C
      DO 20 K = 1, NUMPTMAX
         DO 30 L = 1, J - 1
	    IF (PLTYPE(L).EQ.1) WRITE(LUN,130) Q(PL2WRI(L),K),
     *                          I(PL2WRI(L),K), EROR(PL2WRI(L),K)
	    IF ((PLTYPE(L).LT.0).OR.(PLTYPE(L).GT.1))
     *   CALL ERRMSG("BEEP! BEEP! Something weird's happening in LOTUS")
 30      CONTINUE
         DO 40 L = 1, J - 1
	    IF (PLTYPE(L).EQ.0) WRITE(LUN,140) Q(PL2WRI(L),K),
     *                          I(PL2WRI(L),K)
	    IF ((PLTYPE(L).LT.0).OR.(PLTYPE(L).GT.1))
     *   CALL ERRMSG("BEEP! BEEP! Something weird's happening in LOTUS")
 40      CONTINUE
         WRITE(LUN,150)
 20   CONTINUE
      CLOSE(LUN)

      RETURN
 100  FORMAT("Enter name of file to write the plots to: ",$)
 110  FORMAT("Enter plots to save, by number, terminate with 0 or",/,
     *       " a negative number")
 120  FORMAT("Enter next plot: ",$)
 130  FORMAT(3E13.6,$)
 140  FORMAT(F10.7,E21.14,$)
 150  FORMAT(/,$)
      END
