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

C A.S. Nealis: 25/3/92 - Modified so can reopen files without a user
C prompt (originally to fix a bug in the handling of errors after
C calls to RFACTOR()
      SUBROUTINE OPNFIL(LUN, IRC, FILENAME, PLOTCNT)
      IMPLICIT NONE
      INTEGER LUN, IRC, PLOTCNT
      CHARACTER MESAGE*26, FILENAME(10)*14
      
      IRC = 1
      MESAGE = "Error opening file to read"
      
 1    CONTINUE
      PRINT 100
      READ(*,'(A)') FILENAME(PLOTCNT)
      IRC = 1
      OPEN(UNIT=LUN, FILE=FILENAME(PLOTCNT), STATUS='OLD', ERR=10)
      REWIND(LUN)

      RETURN
 10   CALL ERRMSG(MESAGE)
      GOTO 1

 100  FORMAT("Enter name of file containing data: ",$)
      END
