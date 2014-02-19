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

      SUBROUTINE PTYPGT
      IMPLICIT NONE
      INCLUDE 'SCTPL.COM'
      DOUBLE PRECISION VALUE(6)
      INTEGER IRC, NVAL, DEFTYP

      ISLINFIT(PLOTCNT) = 0
      IF(PLTYPE(PLOTCNT).EQ.0) DEFTYP = 3
      IF(PLTYPE(PLOTCNT).EQ.1) DEFTYP = 2

 1    PRINT 100, DEFTYP
      CALL GETFVAL(5, VALUE, NVAL, IRC)
C Accept default plot type
      IF (IRC.EQ.2) THEN
         PTYPE(PLOTCNT) = DEFTYP
      ELSE IF (NVAL.GT.1) THEN
         CALL ERRMSG("One value expected")
	 GOTO 1
      ELSE
         PTYPE(PLOTCNT) = IDINT(VALUE(1))
      ENDIF
      IF ((PTYPE(PLOTCNT).LT.0).OR.(PTYPE(PLOTCNT).GT.5)) THEN
         CALL ERRMSG("Value out of range")
         GOTO 1
      ENDIF
 2    IF (PTYPE(PLOTCNT).EQ.3) THEN
         CALL CONSGT
      ENDIF
      IF ((PTYPE(PLOTCNT).GE.0).AND.(PTYPE(PLOTCNT).LE.2))
     *             ISLINFIT(PLOTCNT) = 1

      RETURN
 100  FORMAT("Type of plot:",/,
     *       " 0 Ln(IQQ) / Q*Q: Thickness RG",/,
     *       " 1 Ln(IQ) / Q*Q: Cross-sectional RG",/,
     *       " 2 Ln(I) / Q*Q: Overall RG; Guinier plot",/,
     *       " 3 Ln(I) / Q: Wide Angle",/,
     *       " 4 I / Q",/,
     *       " 5 I*Q*Q / Q: Volume plot",/,
     *       "Choose one of the above [",I1,"]: ",$)
      END
