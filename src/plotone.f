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

      SUBROUTINE PLOTONE(PNUM, RDRMOD)
      IMPLICIT NONE
      INCLUDE 'SCTPL.COM'
C     COMMON /PLTONE/ RQMIN1, RQMAX1, RIMIN1, RIMAX1
      COMMON /BOOL1/ AXMODE
C     REAL RQMIN1, RQMAX1, RIMIN1, RIMAX1
      INTEGER RDRMOD, PLSTYL, PNUM
      LOGICAL AXMODE

      ISPLOT(PNUM) = 1
C A.S. Nealis: 26/3/92
C     IF (AXMODE) CALL SETAXES(RDRMOD, RQMIN1, RQMAX1, RIMIN1, RIMAX1)
      IF (AXMODE) CALL SETAXES(RDRMOD)
C-SJP      CALL MAP(0.0, RQMAX1, RIMIN1, RIMAX1)
C-SJP      CALL SCALES
C-SJP      CALL MARKER(SYMBLIST(PNUM))
      CALL PLGRAF(PNUM)
      IF (ISLINFIT(PNUM).EQ.1) CALL DRWLIN(PNUM)

      RETURN
      END

      SUBROUTINE MINMAX(ARRAY,MIN,MAX,N)
      IMPLICIT NONE
      INCLUDE 'SCTPL.COM'
      DOUBLE PRECISION ARRAY(MAXPTS), MIN, MAX
      INTEGER N, II
  
      MIN = ARRAY(1)
      MAX = ARRAY(1)
      DO 10 II = 1, N
         IF(MIN.GT.ARRAY(II)) MIN = ARRAY(II)
         IF(MAX.LT.ARRAY(II)) MAX = ARRAY(II)
10    CONTINUE
  
      RETURN
      END 

C A.S. Nealis: 26/3/92
C     SUBROUTINE PLOTSOME(AXMODE, RDRMOD,
C    *                    RQMIN1, RQMAX1, RIMIN1, RIMAX1)
      SUBROUTINE PLOTSOME(AXMODE)
      IMPLICIT NONE
      INCLUDE 'SCTPL.COM'
C     REAL RQMIN1, RQMAX1, RIMIN1, RIMAX1
      INTEGER RDRMOD, II
      LOGICAL AXMODE

      IF (AXMODE) THEN
C A.S. Nealis: 26/3/92
C        CALL SETAXES(RDRMOD, RQMIN1, RQMAX1, RIMIN1, RIMAX1)
         CALL SETAXES(RDRMOD)
C-SJP         CALL MAP(0.0, RQMAX1, RIMIN1, RIMAX1)
C-SJP      ELSE
C-SJP         CALL MAP(RQMIN1, RQMAX1, RIMIN1, RIMAX1)
      ENDIF
C-SJP      CALL SCALES
      DO 10 II = 1, PLOTCNT
C-SJP         IF ((ISPLOT(II).NE.0).AND.(PTYPE(II).EQ.CURPLOT)) THEN
C-SJP            CALL MARKER(SYMBLIST(II))
C-SJP            CALL PLGRAF(II)
            IF (ISLINFIT(II).EQ.1) CALL DRWLIN(II)
C-SJP         ENDIF
 10   CONTINUE

      RETURN
      END
