PROGRAM brksph2pdb

USE SJP_UTIL
IMPLICIT NONE
  
REAL, DIMENSION(6000,3) :: COORDS
REAL, DIMENSION(6000) :: R
CHARACTER*50    INFIL,OUTFIL,OUTFIL2,QUERY,T
DATA T/'ATOM      1  C1  SER '/
INTEGER NSPHERE,RESNO,C, IOS 
LOGICAL VALID

! Identify program to user
WRITE (*,'(1X,A)') 'Program BRKSP2PDB version 1.00 '
WRITE (*,*)
! Open input file
CALL GET_FILENAME ('Input File:', INFIL)

CALL READ_SPHERE_FILE (10, INFIL, COORDS, R, NSPHERE)

!Open new output file

VALID = .FALSE.

! get user to specify output file location
DO WHILE (.NOT. VALID)
    WRITE(*,'(1X,A)') 'Output File'
    READ (*,'(A)') OUTFIL 
    IF (OUTFIL .EQ. ' ') THEN
        WRITE (*,'(1X,A)') '**ERROR : Please supply a filename - '
    ELSE
        VALID = .TRUE.
    ENDIF
END DO

! Open output file
OPEN(UNIT=11,FILE=OUTFIL,IOSTAT=IOS,STATUS='NEW')
IF (IOS .NE. 0) THEN
    CLOSE (UNIT=10)
    CLOSE (UNIT=11)
    STOP '**ABORT : Problem with output file'
END IF

! Write out the spheres as PDB records (SER residue, C1 atom)
DO C = 1, NSPHERE
      WRITE(11,'(1A22,1I4,1F12.3,2F8.3,1F6.2)') T,C,COORDS(C,1),COORDS(C,1),COORDS(C,1),R(C)
END DO 

CLOSE (UNIT=10)
CLOSE (UNIT=11)
STOP 'Finished !'

END PROGRAM 
