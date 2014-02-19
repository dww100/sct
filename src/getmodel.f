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
      SUBROUTINE MPLTGT(CURVES, BTM)
      IMPLICIT NONE
      INCLUDE 'SCTPL.COM'
      INTEGER CURVES, II, NVAL, IRC
      DOUBLE PRECISION BTM(25), VALUE(6)

 1    PRINT 110
      CALL GETFVAL(5, VALUE, NVAL, IRC)
      IF (IRC.EQ.2) THEN
         CURVES = 1
      ELSE IF (NVAL.GT.1) THEN
         CALL ERRMSG("One value expected")
	 GOTO 1
      ELSE
         CURVES = IDINT(VALUE(1))
      ENDIF
 2    IF ((CURVES.LT.1).OR.(CURVES.GT.MAXCRV)) THEN
         CALL ERRMSG("Too many curves requested!")
	 GOTO 2
      ENDIF
      IF (CURVES.EQ.1) BTM(1) = 1.0
      IF (CURVES.GT.1) THEN
         PRINT *, "Enter weights for curves:"
         DO 10 II = 1, CURVES
	    PRINT 120, II
	    READ(*,*) BTM(II)
 10      CONTINUE
      ENDIF
      RETURN
 110  FORMAT("Number of curves? [1] ",$)
 120  FORMAT("Weight for curve  ",i2,": ",$)
      END

      SUBROUTINE MDATGT(CURVES, ERORCD, NUMQ, WEIGHTED)
C ERORCD = 0: O.K.
C        = 1: Problem (Hamburg format encountered)
C
      IMPLICIT NONE
      INCLUDE 'SCTPL.COM'
      DOUBLE PRECISION WEIGHTED(7680), VALUE(6)
      INTEGER CURVES, LUN, IRC, MODE, N, PCNT, MAXPNT, NVAL, NUMQ,
     *        ERORCD
      REAL COUNT, RINC

      N = 1
      PCNT = PLOTCNT
      LUN = 1
      COUNT = 0.
      MAXPNT = MAXCRV * MAXPTS
      CALL OPNFIL(LUN, IRC, FILENAME, PLOTCNT)

C Read in up to MAXPNT Intensities. Later, Q values are read in. How many
C Q values depends on the value of CURVES
C
 2    IF (INT(COUNT).LE.MAXPNT) THEN
         CALL GETFVAL(LUN,VALUE,NVAL,IRC)
         IF (NVAL.EQ.1) RINC = 1./3.
         IF (NVAL.EQ.2) RINC = 1.
         IF (NVAL.EQ.3) RINC = 1.
         IF (NVAL.EQ.4) RINC = 2.
         IF (NVAL.EQ.6) RINC = 2.
         IF ((INT(COUNT + RINC).GT.MAXPNT).AND.(NVAL.GT.3)) THEN
            IF (INT(COUNT + RINC - 1).LT.MAXPNT) THEN
               IF (NVAL.EQ.4) NVAL = 2
	       IF (NVAL.EQ.6) NVAL = 3
	    ENDIF
         ENDIF
         IF (IRC.EQ.0) THEN
	    N = INT(COUNT + 1)
            IF (NVAL.EQ.1) THEN
	       CALL ERRMSG("Hamburg format not supported for this op")
	       ERORCD = 1
	       RETURN
            ELSE IF(NVAL.EQ.2) THEN
		WEIGHTED(N) = VALUE(2)
		COUNT = COUNT + RINC
		GOTO 2
	     ELSE IF(NVAL.EQ.3) THEN
		WEIGHTED(N) = VALUE(2)
		COUNT = COUNT + RINC
		GOTO 2
	     ELSE IF(NVAL.EQ.4) THEN
		WEIGHTED(N) = VALUE(2)
		WEIGHTED(N+1) = VALUE(5)
		COUNT = COUNT + RINC
		GOTO 2
	     ELSE IF(NVAL.EQ.6) THEN
		WEIGHTED(N) = VALUE(2)
		WEIGHTED(N+1) = VALUE(5)
		COUNT = COUNT + RINC
		GOTO 2
	    ENDIF
         ENDIF
      ENDIF
         REWIND(LUN)
	 NUMQ = COUNT/CURVES
	 NUMPTS(PLOTCNT) = NUMQ
	 COUNT = 0.
C read in numq q-values
 3    IF (INT(COUNT).LE.NUMQ) THEN
         CALL GETFVAL(LUN,VALUE,NVAL,IRC)
         IF (NVAL.EQ.1) RINC = 1./3.
         IF (NVAL.EQ.2) RINC = 1.
         IF (NVAL.EQ.3) RINC = 1.
         IF (NVAL.EQ.4) RINC = 2.
         IF (NVAL.EQ.6) RINC = 2.
         IF ((INT(COUNT + RINC).GT.MAXPNT).AND.(NVAL.GT.3)) THEN
            IF (INT(COUNT + RINC - 1).LT.MAXPNT) THEN
               IF (NVAL.EQ.4) NVAL = 2
	       IF (NVAL.EQ.6) NVAL = 3
	    ENDIF
         ENDIF
         IF (IRC.EQ.0) THEN
	    N = INT(COUNT + 1)
            IF (NVAL.EQ.1) THEN
	       CALL ERRMSG("Hamburg format not supported for this op")
	       ERORCD = 1
	       RETURN
            ELSE IF(NVAL.EQ.2) THEN
	       Q(PLOTCNT,N) = VALUE(1)
	       COUNT = COUNT + RINC
               GOTO 3
	    ELSE IF(NVAL.EQ.3) THEN
	       Q(PLOTCNT,N) = VALUE(1)
	       COUNT = COUNT + RINC
	       GOTO 3
	    ELSE IF(NVAL.EQ.4) THEN
	       Q(PLOTCNT,N) = VALUE(1)
	       Q(PLOTCNT,N+1) = VALUE(3)
	       COUNT = COUNT + RINC
	       GOTO 3
	    ELSE IF(NVAL.EQ.6) THEN
	       Q(PLOTCNT,N) = VALUE(1)
	       Q(PLOTCNT,N+1) = VALUE(4)
               COUNT = COUNT + RINC
	       GOTO 3
	    ENDIF
          ENDIF
       ENDIF
       REWIND(LUN)
       CLOSE (LUN)
       RETURN
       END

      SUBROUTINE ISMEAR(NUMQ)
      IMPLICIT NONE
      INCLUDE 'SCTPL.COM'
      DOUBLE PRECISION WAS, WAV, ALF, QDIG
      INTEGER NUMQ
      CHARACTER * 1 ISME
C Get smearing parameters if smearing requested
      PRINT 130
      READ(*,'(A1)') ISME
      IF ((ISME.EQ.'Y').OR.(ISME.EQ.'y')) THEN
         QDIG = Q(PLOTCNT, 100)/100.
         IF (NUMQ.NE.100) THEN
	    CALL ERRMSG("Can't smear, because not 100 points")
	    RETURN
         ENDIF
         PRINT 140
	 READ(*,*) WAS, WAV, ALF
C Perform smearing on the averaged 'model' spectrum
	 CALL SMEAR(WAS, WAV, ALF, QDIG)
C Set flag for this plot to tell rest of prog we have smeared data
         SMEARD(PLOTCNT) = 1
      ENDIF

      RETURN
 140  FORMAT("Wavelength Spread (dL/L) :",/,"Wavelength A",/,
     *        "Beam Divergence? ",$)
 130  FORMAT("Smearing (y/n)? ",$)
      END

