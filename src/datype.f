      SUBROUTINE DATYPE
      
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

      IMPLICIT NONE
      INCLUDE 'SCTPL.COM'
      CHARACTER*1 DATATYPE
      
 1    PRINT 100
      READ(*,'(A)') DATATYPE(1:1)
      IF ((DATATYPE(1:1).EQ.'M').OR.(DATATYPE(1:1).EQ.'m')) THEN
         PLTYPE(PLOTCNT) = 0
c         CALL GETMODEL
      ELSEIF ((DATATYPE(1:1).EQ.'X').OR.(DATATYPE(1:1).EQ.'x')) THEN
         PLTYPE(PLOTCNT) = 1
c	 print *, pltype
c         CALL GETXPTL
      ELSE
         CALL ERRMSG("Try again")
         GOTO 1
      ENDIF

      RETURN
 100  FORMAT("Enter data type (model [m] or exptl [x]): ",$)
      END
