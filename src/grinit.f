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
      SUBROUTINE GRINIT
      IMPLICIT NONE
      REAL XMIN, XMAX, DX, YMIN, YMAX, DY

c      XMIN = 0.05
      XMIN = 0.075
      DX = 0.87288961
c      XMAX = 1.0
      XMAX = XMIN + DX
      YMIN = 0.1
c      YMIN = 0.05
      DY = 0.850853242
      YMAX = YMIN + DY
c      CALL PAPER(1)
C-SJP      CALL FRAME
C-SJP      CALL PSPACE(XMIN, XMAX, YMIN, YMAX)
C-SJP      CALL BORDER
      RETURN
      END
