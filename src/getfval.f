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
      SUBROUTINE GETFVAL (ITERM,VALUE,NVAL,IRC)
      IMPLICIT NONE
C
C Purpose : Read record from fort.LUN and decode it into numeric
C           values, indicating if an alpha character or <ctrl-Z>
C           has been input. This is essentially a free-format
C           internal file read.
C 
C Limitations: Number of values should not exceed dimensions of
C              "values" in calling module
C
      DOUBLE PRECISION    VALUE(6)
      INTEGER ITERM,NVAL,IRC
C
C ITERM  : Terminal I/O stream
C VALUE  : Data values entered at terminal
C NVAL   : nos. of values entered
C IRC    : Return code 0 - successful
C                      1 - EOF encountered
C                      2 - alph chars encountered
C
C Updates:
C 26/07/84 SZ  Routine now recognises and reads values in E format
C 09/05/85 GRM Renamed from original routine RDVALU to GETVAL so as to
C              include the terminal input & conversion to Fortran-77.
C
C 07/06/91 ASN Robbed from otoko getval.f and adapted.
C 25/06/91 ASN Missed something! Convert all implicit fncns to 
C              DOUBLE PRECISION, ie SIGN -> DSIGN; REAL -> DBLE
C Calls  0:
C
C-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
C Local Variables:
C
      DOUBLE PRECISION   EXPO,SEXPO,T,W,Z,Y,S
      INTEGER      J,N,IGETS,ICOUNT, iii, numspace
      LOGICAL	   DIGIT,EXPON
      CHARACTER    SPACE,DOT,MINUS,ZERO,ECHAR,DCHAR,PLUS,COMMA,ELCHAR
      CHARACTER*80 LINE
C 
C LINE   : String to hold terminal input
C EXPON  : when true then exponent is being read
C BLANK  :
C DOT    :
C MINUS  :
C ZERO   :
C EXPO   : Value of exponent
C SEXPO  : Sign of exponent
C 
      DATA SPACE/' '/ , ZERO/'0'/  , MINUS/'-'/ , DOT/'.'/
      DATA ECHAR/'E'/ , DCHAR/'D'/ , PLUS/'+'/  , COMMA/','/
      DATA ELCHAR/'e'/
C 
C-----------------------------------------------------------------------
C 
C========INITIALIZE
C
      NUMSPACE = 0
      IRC=1
      NVAL=0
      T=10.0
      W=0.0
      Z=1.0
      Y=1.0
      S=1.0
      EXPO=0.0
      SEXPO=1.0
      EXPON=.FALSE.
      DIGIT=.FALSE.
      DO 20 ICOUNT = 1, 6
 20       VALUE(ICOUNT) = 0.0
C
C========READ RECORD FROM FILE
C
      READ(ITERM,'(A)',END=998, ERR=999) LINE
      DO 30 III = 1, LEN(LINE)
         IF(LINE(III:III).EQ.SPACE) NUMSPACE = NUMSPACE + 1
 30   CONTINUE
      IF(NUMSPACE.GE.80) GOTO 999
C 
C========SEARCH FOR DECIMAL DIGITS UNTIL END OF RECORD
C 
      DO 10 J=1,LEN(LINE)
         N=ICHAR (LINE(J:J))-ICHAR (ZERO)
C 
C========STORE DECIMAL DIGIT
C 
         IF (N.GE.0.AND.N.LE.9) THEN
            DIGIT=.TRUE.
            IF (EXPON) THEN
               EXPO=EXPO*10.0+DBLE(N)
            ELSE
               Z=Z*Y
               W=W*T+Z*DBLE(N)
            ENDIF
C 
C========CHECK FOR DECIMAL POINT
C 
         ELSEIF (LINE(J:J).EQ.DOT.AND.(.NOT.EXPON)) THEN
            Y=0.1
            T=1.0
C 
C========CHECK FOR MINUS SIGN
C 
         ELSEIF (LINE(J:J).EQ.MINUS) THEN
            IF (EXPON) THEN
               SEXPO=-1.0
            ELSE
               S=-1.0
            ENDIF
C 
C========CHECK FOR PLUS
C 
         ELSEIF (LINE(J:J).EQ.PLUS) THEN
            CONTINUE
C 
C========CHECK FOR EXPONENT
C 
         ELSEIF ((LINE(J:J).EQ.ECHAR.OR.LINE(J:J).EQ.ELCHAR.OR.
     1       LINE(J:J).EQ.DCHAR).AND.DIGIT) THEN
            IF (.NOT.EXPON) THEN
               EXPON=.TRUE.
            ELSE
               IRC=1
               GOTO 999
            ENDIF
C 
C========ASSIGN NUMERIC VALUE IF SPACE OR COMMA DELIMITER
C 
         ELSEIF ((LINE(J:J).EQ.COMMA.OR.LINE(J:J).EQ.SPACE).AND.
     1            DIGIT) THEN
            NVAL=NVAL+1
            VALUE(NVAL)=DSIGN(W,S)*10**(DSIGN(EXPO,SEXPO))
            T=10.0
            W=0.0
            Z=1.0
            Y=1.0
            S=1.0
            EXPO=0.0
            SEXPO=1.0
            EXPON=.FALSE.
            DIGIT=.FALSE.
C
C========CHECK FOR SPACES
C
         ELSEIF (LINE(J:J).NE.SPACE) THEN
            CALL ERRMSG ('Error: Numeric input expected')
            IRC=2
            GOTO 999
         ENDIF
10    CONTINUE
      IRC=0
      IF(NUMSPACE.EQ.80) IRC=2
      RETURN
998   IRC=1
999   IF(NUMSPACE.GE.80) IRC=2
      RETURN
      END
