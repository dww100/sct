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
      SUBROUTINE PRTSCR
      IMPLICIT NONE
C 
C Purpose: Print TEKT graphics screen of a PERICOM MG100 terminal.
C
C Calls   0:
C 
C-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
C Local variables:
C 
      INTEGER     L,LEN
      CHARACTER*1 ESC,SEVEN,ZERO,TWENTY9,TWENTY3,NINETY2,TWENTY4
C 
C ESC    : Escape character
C
C      DATA ESC/27/,SEVEN/55/,ZERO/48/,TWENTY9/29/,TWENTY3/23/
C     *     NINETY2/92/,TWENTY4/24/
      ESC = CHAR(27)
      SEVEN = CHAR(55)
      ZERO = CHAR(48)
      TWENTY9 = CHAR(29)
      TWENTY3 = CHAR(23)
      NINETY2 = CHAR(92)
      TWENTY4 = CHAR(24)
C
C-----------------------------------------------------------------------
C
C======Switch to Graphics mode, dump screen, return to ALPHA mode!
C
      PRINT *, ESC,NINETY2,TWENTY9,ESC,TWENTY3,TWENTY4
C
      END


