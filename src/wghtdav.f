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
      SUBROUTINE WGHTDAV(ARRAY, WEIGHTS, NUMQ, CURVES)
      IMPLICIT NONE
      INCLUDE 'SCTPL.COM'
      DOUBLE PRECISION ARRAY(7680), WEIGHTS(25)
      INTEGER NUMQ, II, J, CURVES
      
      DO 20 J = 1, CURVES
         DO 10 II = 1, NUMQ
            I(PLOTCNT, II) = I(PLOTCNT, II) + ARRAY((J-1)*NUMQ+II)
     *                       * WEIGHTS(J)
 10      CONTINUE
 20   CONTINUE
      END
      
