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
      SUBROUTINE TOALPH
      IMPLICIT NONE
C 
C Purpose: Switch from graphics mode to alpha mode on a Pericom
C          terminal.
C
C Calls   0:
C 
C-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
C Local variables:
C 
      CHARACTER*1 CTRLX
C
C      DATA CTRLX /24/
       CTRLX = CHAR(24)
C 
C-----------------------------------------------------------------------
C 
C======PRINT ERROR MESSAGE AND RING BELL IN INVERSE VIDEO
C 
      PRINT 1000 ,CTRLX
      RETURN
C 
1000  FORMAT (80A1)
      END
