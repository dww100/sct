C****************************************************
C A.S. NEALIS 23/5/89
C
C Incorporated into sctpl4 (SG version) 18/6/91
C READ(8,) Changed to READ(*,)
C
C IPLT1, IPLT2 are the exponentiated I(PLOT1,) [model] and I(PLOT2,) [xptal]
C data.
C
C Subroutine for calculating crystallographic r-factor for a scatt model vs
C   exptal data.
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
      SUBROUTINE RFACTOR(CONS)
      IMPLICIT NONE
      INCLUDE 'SCTPL.COM'
      COMMON /RFD/ RFCTR, CONST
      COMMON /RFI/ PLOT1, PLOT2, NUMPT2
      COMMON /BOOL2/ AUTORF, IGNORE
      COMMON /BOOLRF/ RFMODE
      CHARACTER ANSWER*1, FLNAME*8, OLDFLNM*8
      DOUBLE PRECISION SIMQ(100), SIMI(100), QEXPT(512), IEXPT(512),
     *                 ICALC(512), CONST(11), RFCTR(11), RINC, CON,
     *                 SUMIOBS, SUMIOIC, CONS, DQMIN, DELTAQ, RFACT,
     *                 OLDCONS, DELTAC, SUMIXPT, SUMICAL, COUNT,
     *                 AVIEXPT, AVICALC, RMIN, IPLT1(MAXPTS),
     *                 IPLT2(MAXPTS),FQMAX
      INTEGER II, J, K, NPTS, IIGNOR1, IIGNOR2, ICNT, PLOT1, PLOT2,
     *        NUMPT2, OPLTCNT
      LOGICAL IGNORE(512), AUTORF, RFMODE
C A.S. Nealis 17/1/90
C Added these declarations for iterative r-fact finding
      DOUBLE PRECISION R(4)
      INTEGER SCORE

C RINC is set to define the increment in the constant used to calculate
C the range of rfactors (RINC = 0.01 USUALLY)
      RINC = 0.01
      CONS = CONSTANT(PLOT2)
      OLDCONS = CONS
      DELTAC =  5.12
      DO 50 II = 1, MAXPTS
         IGNORE(II) = .FALSE.
 50   CONTINUE
  
      IF(.NOT.AUTORF) THEN
C These next lines belonged to the old "manual" version of rfactor
         PRINT *, 'Value of constant to use & increment'
         READ(*,*) CONST(6), RINC
         CONST(1) = CONST(6)
      ENDIF
C******************************************
C Match exptal data with simulated data
C
      J = 0
      DO 1000 II = 1, NUMPTS(PLOT2)
         IF(Q(PLOT2, II).LE.Q(PLOT1, NUMPTS(PLOT1))) J = J + 1
 1000 CONTINUE
C*******************************************
C Sort out which points to ignore
C
 60   PRINT *, 'Enter range of xptal points to ignore,'
      PRINT *, '0 or -ve number to stop.'
      READ(*,*) IIGNOR1,IIGNOR2
      IF(IIGNOR2.LT.IIGNOR1) GOTO 60
      IF((IIGNOR1.LE.0).OR.(IIGNOR2.LE.0)) THEN
         PRINT *, 'Enter y to stop edit, anything else to continue ',
     *               'the edit.'
         READ(*,'(A)') ANSWER
         IF((ANSWER.NE.'Y').AND.(ANSWER.NE.'y')) GOTO 60
      ENDIF
      IF((IIGNOR1.LE.J).AND.(IIGNOR1.GT.0).AND.
     *   (IIGNOR2.LE.J).AND.(IIGNOR2.GT.0)) THEN
         DO 70 II = IIGNOR1, IIGNOR2
            IF(IGNORE(II)) THEN
               PRINT *, 'Channel ',II,' already ignored. type y to'
     *              ,' unignore or anything else to leave it ignored.'
               READ(*,'(A)') ANSWER
               IF(ANSWER.EQ.'Y') IGNORE(II) = .FALSE.
            ELSE
               IGNORE(II) = .TRUE.
            ENDIF
 70      CONTINUE
         GOTO 60
      ENDIF
      IF((IIGNOR1.GT.J).OR.(IIGNOR2.GT.J)) THEN
         CALL ERRMSG("Invalid data range.")
         GOTO 60
       ENDIF

       WRITE(*,*) 'Fit Rfactor Based on what Max Q?'
       READ(*,*) FQMAX
C*******************************************
C For each qexpt, find corresponding qsim. store corresponding simi in
C   ICALC. NB- only done for 0 <= QEXPT <= .16
C
      DQMIN = DABS (Q(PLOT2, 1) - Q(PLOT1, 1))
      DO 10 II = 1, J
         IF(.NOT.IGNORE(II)) THEN
C Set initial minimum to a very big number
            DQMIN = 1.0E+13
            DO 20 K = 1, 100
               DELTAQ = DABS (Q(PLOT2, II) - Q(PLOT1, K))
               IF(DELTAQ.LT.DQMIN) THEN 
                  DQMIN = DELTAQ
                  ICALC(II) = I(PLOT1, K)
                  
                ENDIF
 20         CONTINUE
c            write(3,*) Q(PLOT2,II),ICALC(II),I(PLOT2,II)
          ENDIF
 10   CONTINUE

      IF(AUTORF) THEN
C A.S. Nealis 18/1/90
C Calculate out the average intensities for ICALC and IEXPT
C then use the difference to get a good starting point for the
C iteration
         SUMIXPT = 0.
         SUMICAL = 0.
         COUNT = 0.
         DO 80 II = 1, J
            IF(.NOT.IGNORE(II)) SUMIXPT = SUMIXPT + I(PLOT2, II)
            COUNT = COUNT + 1
 80      CONTINUE
         AVIEXPT = SUMIXPT/COUNT
         DO 85 II = 1, 100
            SUMICAL = SUMICAL + I(PLOT1, II)
            AVICALC = SUMICAL/100.
 85      CONTINUE
         CONST(1) = AVIEXPT - AVICALC
 90      IF(DELTAC.LT.RINC) DELTAC = RINC
         DO 100 II = 1, 3
            CONST(II+1) = CONST(II) + DELTAC
 100     CONTINUE

         IF(DELTAC.EQ.RINC) GOTO 600
         DO 500 II = 1, 4
            R(II) = 0.
            SUMIOBS = 0.0
            SUMIOIC = 0.0
         DO 550 ICNT = 1, J
         IF ((.NOT.IGNORE(ICNT)).AND.(Q(PLOT2,ICNT).LE.FQMAX)) THEN
                  SUMIOIC=SUMIOIC+DABS(EXP(CONST(II)-CONS+
     *                    I(PLOT2,ICNT))-EXP(ICALC(ICNT)))
                  SUMIOBS=SUMIOBS+DABS(EXP(CONST(II)-CONS+
     *                    I(PLOT2,ICNT)))
	       RFACT = SUMIOIC/SUMIOBS
						      
               ENDIF
 550        CONTINUE
 500     R(II) = RFACT

         IF (R(1).EQ.R(2)) THEN
            DELTAC = DELTAC/3.
         ELSE IF (R(2).EQ.R(3)) THEN
            CONST(1) = CONST(2)
            DELTAC = DELTAC/3.
         ELSE IF (R(3).EQ.R(4)) THEN
            CONST(1) = CONST(3)
            DELTAC = DELTAC/3.
         ENDIF
 600     CONTINUE
C Got to find out which points the minimum is between
C
         SCORE = 0
         IF(R(1).GT.R(2)) SCORE = SCORE + 1
         IF(R(1).LT.R(2)) SCORE = SCORE - 1
         IF(R(2).GT.R(3)) SCORE = SCORE + 1
         IF(R(2).LT.R(3)) SCORE = SCORE - 1
         IF(R(3).GT.R(4)) SCORE = SCORE + 1
         IF(R(3).LT.R(4)) SCORE = SCORE - 1

         IF(SCORE.EQ.1) THEN
            CONST(1) = CONST(2)
            CONST(6) = CONST(2)
            DELTAC = DELTAC/2.
         ELSE IF(SCORE.EQ.-1) THEN
            CONST(1) = CONST(1)
            CONST(6) = CONST(1)
C I know it looks pointless, but I think it makes for more readable code.
            DELTAC = DELTAC/2.
         ELSE IF(SCORE.EQ.3) THEN
            CONST(1) = CONST(3)
            CONST(6) = CONST(3)
         ELSE IF(SCORE.EQ.-3) THEN
            CONST(1) = CONST(1) - DELTAC
         ENDIF

         IF(DELTAC.GT.RINC) THEN
            GOTO 90
         ELSE IF(DELTAC.LT.RINC) THEN
            DELTAC = RINC
         ENDIF
      ENDIF
C*******************************************
C Calculate rfact over a range, centred about the chosen const.
C
      DO 2000 II = 1, 11
         CON = CONST(6) + RINC*(II-6)
         CONST(II) = CON
         SUMIOBS = 0.0
         SUMIOIC = 0.0
         DO 1300 ICNT = 1, J
         IF ((.NOT.IGNORE(ICNT)).AND.(Q(PLOT2,ICNT).LE.FQMAX)) THEN
                  SUMIOIC=SUMIOIC+DABS(EXP(CONST(II)-CONS+
     *                    I(PLOT2,ICNT))-EXP(ICALC(ICNT)))
                  SUMIOBS=SUMIOBS+DABS(EXP(CONST(II)-CONS+
     *                    I(PLOT2,ICNT)))
	       RFACT = SUMIOIC/SUMIOBS
            ENDIF
 1300    CONTINUE
         RFCTR(II) = RFACT
         write(3,*) SUMIOBS,SUMIOIC,CON
 2000  CONTINUE

       write(3,*) CON,RINC,CONS
       do 1301 ICNT=1,J
       write(3,*) IGNORE(ICNT),Q(PLOT2,ICNT),I(PLOT2,ICNT),
     *  ICALC(ICNT)
1301   CONTINUE      

c      C A.S. Nealis 18/1/90
C  Find the const corresponding to the smallest RFCTR(II)
C
      IF (AUTORF) THEN
         RMIN = RFCTR(1)
         DO 99 II = 1, 11
            IF(RMIN.GT.RFCTR(II)) CONS = CONST(II)
            IF(RMIN.GT.RFCTR(II)) RMIN = RFCTR(II)
 99      CONTINUE
      ELSE
	 CONS = CONST(6)
      ENDIF
      RETURN
      END 
