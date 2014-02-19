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

      SUBROUTINE INIT
      IMPLICIT NONE
      INCLUDE 'SCTPL.COM'
      INTEGER II, J
      DATA SYMBLIST /224,227,229,225,228,240,243,245,241,244/

      DO 10 II = 1, MAXPLOT
         FSTPLTPT(II) = 0
	 LSTPLTPT(II) = 0
	 CONSTANT(II) = 0
	 QVALP1(II) = 0.
	 QVALPL(II) = 0.
	 QVALF1(II) = 0.
	 QVALFL(II) = 0.
	 ISLINFIT(II) = 0
	 FSTFITPT(II) = 0
	 PLTYPE(II) = 0
	 SYMBTYP(II) = 0
	 LINETYP(II) = 0
	 NUMPTS(II) = 0
	 PTYPE(II) = 0
	 LSTFITPT(II) = 0
	 SMEARD(II) = 0
	 ISPLOT(II) = 1
         DO 20 J = 1, MAXPTS
            Q(II,J) = 0.
	    I(II,J) = 0.
	    EROR(II,J) = 0.
	    ERORU(II,J) = 0.
	    ERORL(II,J) = 0.
 20      CONTINUE
 10   CONTINUE
      DO 30 II = 1, MAXPTS
         QTEMP(II) = 0.
	 ITEMP(II) = 0.
 30   CONTINUE
      RQMIN1 = 0.0
      RQMAX1 = 0.3
      RIMIN1 = 1.
      RIMAX1 = 9.
      RETURN
      END
