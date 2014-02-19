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

      SUBROUTINE ERRMSG (MESAGE)
      IMPLICIT NONE
C 
C Purpose: Print error message in inverse video on a VT100 emulator
C          terminal.
C 
      CHARACTER*(*) MESAGE
C 
C MESAGE : Error message to be displayed
C 
C Calls   0:
C 
C-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
C Local variables:
C 
      INTEGER     L,LEN
      CHARACTER*1 BELL,ESC,SEVEN,ZERO
C 
C BELL   : Bell ascii character
C ESC    : Escape character
C 
C      DATA BELL/CHAR(7)/ , ESC/27/ , SEVEN/55/ , ZERO/48/
      BELL = CHAR(7)
      ESC = CHAR(27)
      SEVEN = CHAR(55)
      ZERO = CHAR(48)
C 
C-----------------------------------------------------------------------
C 
C======PRINT ERROR MESSAGE AND RING BELL IN INVERSE VIDEO
C 
      PRINT 1000 , ESC , '[' , SEVEN , 'm' ,
     &            (MESAGE(L:L),L=1,LEN(MESAGE)),
     &             BELL , ESC , '[' , ZERO , 'm'
      RETURN
C 
1000  FORMAT (' ',80A1)
      END


