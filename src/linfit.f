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

      SUBROUTINE LINFIT(MODE,A,SIGMAA,B,SIGMAB,R)
      IMPLICIT NONE
      INCLUDE 'SCTPL.COM'
      DOUBLE PRECISION A, SIGMAA, B, SIGMAB, R, SUM, SUMX, SUMX2
      DOUBLE PRECISION SUMY, SUMY2, SUMXY, X1, Y1, C, DELTA, VAR
      DOUBLE PRECISION WEIGHT, SIGMAY(MAXPTS)
      INTEGER NPTS, N1, N2, II, MODE, N
C A.S Nealis 18 June 1991: 
C Stripped out some 'deadwood' from the version implemented in sctpl;
C eg MODE ALWAYS ZERO!! else would have problems with 1./SIGMAY(I),
C etc etc..
C
C Converted to Double Precision arithmetic, and some variable names were
C changed to protect the innocent and maintain consistency with the rest
C of the sctpl4 (SG implementation)
C 
C      SUBROUTINE LINFIT(X,Y,SIGMAY,N1,N2,MODE,A,SIGMAA,B,SIGMAB,R)
C      DIMENSION X(450),Y(450),SIGMAY(450)
C
C      LEAST-SQUARES FIT TO A STRAIGHT LINE
C      VERSION FROM BEVINGTON 
C      IMPLEMENTED SC 26/5/1981
C
C
C
C      MODE:  DETERMINES METHOD OF WEIGHTING LEAST-SQUARES FIT
C           +1 WEIGHT(I)=1./SIGMAY(I)**2
C            0 WEIGHT(I)=1.
C           -1 WEIGHT(I)=1./Y(I)
C         A:  Y INTERCEPT OF FITTED STRAIGHT LINE 
C         B:  SLOPE OF FITTED STRAIGHT LINE
C    SIGMAA:  STANDARD DEVIATION OF A
C    SIGMAB:  STANDARD DEVIATION OF B
C
C
C      ACCUMULATE WEIGHTED SUMS
C
c      DO 1 N=1,450
      DO 1 N=1,MAXPTS
      SIGMAY(N)=0.0 
1     CONTINUE
c      NPTS = N2 - N1 + 1
C A.S. Nealis 18/6/91: Force MODE = 0 - it was never used any other way
C in SCTPL on ICCC.
      MODE = 0
      SUM = 0.
      SUMX = 0.
      SUMY = 0.
      SUMX2 = 0.
      SUMY2 = 0.
      SUMXY = 0.
C      DO 50 I = N1, N2
C         X1 = X(I)
C         Y1 = Y(I)
      NPTS = LSTFITPT(PLOTCNT) - FSTFITPT(PLOTCNT) + 1
      DO 50 II = FSTFITPT(PLOTCNT), LSTFITPT(PLOTCNT)
         X1 = Q(PLOTCNT, II)
         Y1 = I(PLOTCNT, II)
C         IF (MODE) 31,36,38
         IF (MODE .LT. 0) THEN
            GOTO 31
         ELSE IF (MODE .EQ. 0) THEN
            GOTO 36
         ELSE
            GOTO 38
         END IF
C 31      IF (Y1) 34,36,32
31       IF (Y1 .LT. 0) THEN
            GOTO 34
         ELSE IF (MODE .EQ. 0) THEN
            GOTO 36
         ELSE
            GOTO 32
         ENDIF
 32      WEIGHT = 1. / Y1
         GO TO 41
 34      WEIGHT = 1. / (-Y1)
         GO TO 41
 36      WEIGHT = 1.
         GO TO 41
c 38      CONTINUE
 38      WEIGHT = 1. / SIGMAY(II)**2
c 38      WEIGHT = 1. / SIGMAY(I)**2
 41      SUM = SUM + WEIGHT
         SUMX = SUMX + WEIGHT * X1
         SUMY = SUMY + WEIGHT * Y1
         SUMX2 = SUMX2 + WEIGHT * X1 * X1
         SUMXY = SUMXY + WEIGHT * X1 * Y1
         SUMY2 = SUMY2 + WEIGHT * Y1 * Y1
 50   CONTINUE
C
C      CALCULATE COEFFICIENTS AND STANDARD DEVIATIONS
C
 51   DELTA = SUM * SUMX2 - SUMX * SUMX
      A = (SUMX2 * SUMY - SUMX * SUMXY) / DELTA
 53   B = (SUMXY * SUM - SUMX * SUMY) / DELTA
C61   IF (MODE) 62,64,62
61    IF (MODE .LT. 0) THEN
          GOTO 62
      ELSE IF (MODE .EQ. 0) THEN
          GOTO 64
      ELSE
          GOTO 62
      ENDIF
 62   VAR = 1.
      GO TO 67
 64   C = DBLE(NPTS-2)
c 64   C = FLOAT(NPTS-2)
      VAR = (SUMY2 + A * A * SUM + B * B * SUMX2
     *    -2. * (A * SUMY + B * SUMXY - A * B * SUMX)) / C 
 67   SIGMAA = DSQRT (VAR * SUMX2 / DELTA)
 68   SIGMAB = DSQRT (VAR * SUM / DELTA)
 71   R = (SUM * SUMXY - SUMX * SUMY) /
     *     DSQRT(DELTA * (SUM * SUMY2 - SUMY * SUMY))
 
      RETURN
      END 
