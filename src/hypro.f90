PROGRAM hypro

! A program adapted from hydrate to give a shell of spheres
! around a brktos sphere file
! Originally produced (Fortran 77) by S. J. Perkins
! Fortran 90 version coded by David W. Wright November 2013

REAL                X,Y,Z,R
CHARACTER*50        INFIL, OUTFIL, OUTFIL2, QUERY, T, E, S
DATA T/'ATOM      1  C1  SER '/
INTEGER             RESNO, NOOFSPHERES, C, IERR

NOOFSPHERES = 1

C = 1

! Identify program to user

WRITE (*,'(1X,A)') 'Program HYPRO version 1.00 awa970218'
WRITE (*,*)

! TODO: some error checking on the input/output here would be nice

! Get input filename from user and open

WRITE(*,'(1X,A)') 'Input file :'
READ (*,'(A)')  INFIL
OPEN(UNIT=10,FILE=INFIL,ERR=5000,STATUS='OLD')

! Open new output file

WRITE(*,'(1X,A)') 'Output File'
READ (*,'(A)') OUTFIL

OPEN(UNIT=11,FILE=OUTFIL,ERR=6000,STATUS='NEW')


! Get the number hydration spheres to add to each existing sphere from user
WRITE(*,'(1X,A)') 'How many spheres [1,7,13,19,27 Recomended]'
READ (*,'(I2)') NOOFSPHERES

IF (NOOFSPHERES .EQ. NULL) THEN
    WRITE (*,'(1X,A)') 'Default set to [1]: '
ENDIF



! Read input file until all entries read

DO
    READ (UNIT=10,FMT='(4F10.2)',IOSTAT=IERR) X,Y,Z,R

    IF (IERR .LT. 0) THEN
        EXIT
    END IF

! Write spheres in up to 27 places:
! 1 = original sphere positions, others are the corners midpounts and
! centres of the sides of a cube centred on the original sphere

    IF(NOOFSPHERES.GE.1) THEN
        write(11,'(1A20,1I6,1F12.3,2F8.3,1F6.2)') T,C,X,Y,Z,R
        C = C + 1
    ENDIF

    IF(NOOFSPHERES.GE.2) THEN
        write(11,'(1A20,1I6,1F12.3,2F8.3,1F6.2)') T,C,X-R*2,Y,Z,R
        C = C + 1
    ENDIF
    IF(NOOFSPHERES.GE.3) THEN
        write(11,'(1A20,1I6,1F12.3,2F8.3,1F6.2)') T,C,X+R*2,Y,Z,R
        C = C + 1
    ENDIF
    IF(NOOFSPHERES.GE.4) THEN
        write(11,'(1A20,1I6,1F12.3,2F8.3,1F6.2)') T,C,X,Y,Z+R*2,R
        C = C + 1
    ENDIF
    IF(NOOFSPHERES.GE.5) THEN
        write(11,'(1A20,1I6,1F12.3,2F8.3,1F6.2)') T,C,X,Y,Z-R*2,R
        C = C + 1
    ENDIF
    IF(NOOFSPHERES.GE.6) THEN
        write(11,'(1A20,1I6,1F12.3,2F8.3,1F6.2)') T,C,X,Y+R*2,Z,R
        C = C + 1
    ENDIF
    IF(NOOFSPHERES.GE.7) THEN
        write(11,'(1A20,1I6,1F12.3,2F8.3,1F6.2)') T,C,X,Y-R*2,Z,R
        C = C + 1
    ENDIF


    IF(NOOFSPHERES.GE.8) THEN
        write(11,'(1A20,1I6,1F12.3,2F8.3,1F6.2)') T,C,X-R*2,Y,Z+R*2,R
        C = C + 1
    ENDIF
    IF(NOOFSPHERES.GE.9) THEN
        write(11,'(1A20,1I6,1F12.3,2F8.3,1F6.2)') T,C,X+R*2,Y,Z-R*2,R
        C = C + 1
    ENDIF
    IF(NOOFSPHERES.GE.10) THEN
        write(11,'(1A20,1I6,1F12.3,2F8.3,1F6.2)') T,C,X,Y+R*2,Z+R*2,R
        C = C + 1
    ENDIF
    IF(NOOFSPHERES.GE.11) THEN
        write(11,'(1A20,1I6,1F12.3,2F8.3,1F6.2)') T,C,X,Y+R*2,Z-R*2,R
        C = C + 1
    ENDIF
    IF(NOOFSPHERES.GE.12) THEN
        write(11,'(1A20,1I6,1F12.3,2F8.3,1F6.2)') T,C,X-R*2,Y-R*2,Z,R
        C = C + 1
    ENDIF
    IF(NOOFSPHERES.GE.13) THEN
        write(11,'(1A20,1I6,1F12.3,2F8.3,1F6.2)') T,C,X+R*2,Y-R*2,Z,R
        C = C + 1
    ENDIF

    IF(NOOFSPHERES.GE.14) THEN
        write(11,'(1A20,1I6,1F12.3,2F8.3,1F6.2)') T,C,X-R*2,Y,Z-R*2,R
        C = C + 1
    ENDIF
    IF(NOOFSPHERES.GE.15) THEN
        write(11,'(1A20,1I6,1F12.3,2F8.3,1F6.2)') T,C,X+R*2,Y,Z+R*2,R
        C = C + 1
    ENDIF
    IF(NOOFSPHERES.GE.16) THEN
        write(11,'(1A20,1I6,1F12.3,2F8.3,1F6.2)') T,C,X,Y-R*2,Z-R*2,R
        C = C + 1
    ENDIF
    IF(NOOFSPHERES.GE.17) THEN
        write(11,'(1A20,1I6,1F12.3,2F8.3,1F6.2)') T,C,X,Y-R*2,Z+R*2,R
        C = C + 1
    ENDIF
    IF(NOOFSPHERES.GE.18) THEN
        write(11,'(1A20,1I6,1F12.3,2F8.3,1F6.2)') T,C,X-R*2,Y+R*2,Z,R
        C = C + 1
    ENDIF
    IF(NOOFSPHERES.GE.19) THEN
        write(11,'(1A20,1I6,1F12.3,2F8.3,1F6.2)') T,C,X+R*2,Y+R*2,Z,R
        C = C + 1
    ENDIF

    IF(NOOFSPHERES.GE.20) THEN
        write(11,'(1A20,1I6,1F12.3,2F8.3,1F6.2)') T,C,X+R*2,Y+R*2,Z+R*2,R
        C = C + 1
    ENDIF
    IF(NOOFSPHERES.GE.21) THEN
        write(11,'(1A20,1I6,1F12.3,2F8.3,1F6.2)') T,C,X+R*2,Y-R*2,Z+R*2,R
        C = C + 1
    ENDIF
    IF(NOOFSPHERES.GE.22) THEN
        write(11,'(1A20,1I6,1F12.3,2F8.3,1F6.2)') T,C,X-R*2,Y+R*2,Z-R*2,R
        C = C + 1
    ENDIF
    IF(NOOFSPHERES.GE.23) THEN
        write(11,'(1A20,1I6,1F12.3,2F8.3,1F6.2)') T,C,X-R*2,Y-R*2,Z-R*2,R
        C = C + 1
    ENDIF
    IF(NOOFSPHERES.GE.24) THEN
        write(11,'(1A20,1I6,1F12.3,2F8.3,1F6.2)') T,C,X+R*2,Y+R*2,Z-R*2,R
        C = C + 1
    ENDIF
    IF(NOOFSPHERES.GE.25) THEN
        write(11,'(1A20,1I6,1F12.3,2F8.3,1F6.2)') T,C,X-R*2,Y-R*2,Z+R*2,R
        C = C + 1
    ENDIF
    IF(NOOFSPHERES.GE.26) THEN
        write(11,'(1A20,1I6,1F12.3,2F8.3,1F6.2)') T,C,X-R*2,Y+R*2,Z+R*2,R
        C = C + 1
    ENDIF
    IF(NOOFSPHERES.GE.27) THEN
        write(11,'(1A20,1I6,1F12.3,2F8.3,1F6.2)') T,C,X+R*2,Y-R*2,Z-R*2,R
        C = C + 1
    ENDIF

    C = 1

END DO

CLOSE (UNIT=10)
CLOSE (UNIT=11)
STOP 'Finished !'

END PROGRAM