PROGRAM BRKTOS

! Program to convert protein coordinates in PDB format into corresponding
! Debye spheres for small angle scattering calculations

! History of the original Fortran 77 code:
! Original version Stephen J. Perkins 1982
! Revised by SJP until 1987
! Further revised by K.F. SMITH 17/7/91
! SGI version 1.00 A. S. Nealis 11/02/1992
! - Gives same volume for CA only and complete files
! SGI version 1.10 M. O. Mayans 22/02/1994 
! - Added carbohydrate support
! SGI version 1.15 M. K. Boehm ??/01/1995
! - Added SIA as a carbohydrate

! Linux Fortran 90 version 1.20
! David W. Wright 30/10/2013

USE SJP_UTIL
IMPLICIT NONE

DOUBLE PRECISION X(3)
INTEGER INCX, INCY, INCZ
! Always make these numbers = dimensions of INC(,,) below!!!
DATA INCX, INCY, INCZ /400, 400, 400/
INTEGER INC(400,400,400), PDBFIL, SPHFIL, RESNUM, LAST_RESN, STDERR
INTEGER ICUT, IPRINT, JUMP, NRES, NATOMM, N, I, J, K, ISUM
INTEGER IBALL, ILOOP, NAT, MX, MY, MZ, NX, NXX, NY, NYY, NZ, NZZ
INTEGER IERR
DOUBLE PRECISION BOX, RADIS, VOLUME, XMA, XMI, YMA, YMI, ZMA, ZMI
DOUBLE PRECISION XB1, XB2, YB1, YB2, ZB1, ZB2, XX, YY, ZZ
DOUBLE PRECISION VOL, V
CHARACTER INFILE*80, OUTFIL*80, ATOM*4, AANAM*3
DOUBLE PRECISION VOLL

! Read format for ATOM cards in PDBS
CHARACTER(LEN=*), PARAMETER  :: ATMFMT = "(A4,13X,A3,2X,I4,4X,3F8.3)"

! Write format for unrecognised amino acid
CHARACTER(LEN=*), PARAMETER  :: ERR1FMT = "(A,A3,A,I4,A,A10)"
! Write format for incorrect box dimension
CHARACTER(LEN=*), PARAMETER  :: ERR2FMT = "(A,I4,A,I4)"

! Name file unit numbers for particular uses
! PDBFIL and SPHFIL are the unit nos for OPEN statements
PDBFIL = 8
SPHFIL = 9
STDERR = 0

CALL GET_FILENAME ('Enter name of Brookhaven format file (PDB): ', INFILE)
OPEN(UNIT=PDBFIL, FILE=INFILE, ACCESS='SEQUENTIAL', STATUS='OLD')

! Prompt user for parameters

! Output types (IPRINT):
! 0: Coords + info (to file)
! 1: Coords only
! 2: info only (to screen)
! 3: Coords (to file) + info (to screen)

WRITE(*,"(A)") "Enter side of box, cutoff, output type and grid type"
READ(5,*) BOX, ICUT, IPRINT, JUMP

IF ( (IPRINT .LT. 0) .OR. (IPRINT .GT. 3) ) THEN
    WRITE(STDERR,'("IPRINT reset to 0")')
    IPRINT = 0
ENDIF

IF ( (JUMP .NE. 0) .OR. (JUMP .NE. 1) ) THEN
    WRITE(STDERR,'("JUMP reset to 0")')
    JUMP = 0
ENDIF

IF(IPRINT .NE. 2) THEN
    PRINT *, 'Enter filename for output coords'
    READ(*,'(A)') OUTFIL
    OPEN(UNIT=SPHFIL, FILE=OUTFIL, STATUS="NEW")
ELSE
    ! Redirect info only to stdout
    SPHFIL = 6
ENDIF

! Initialize grid dimensions

XMA = -99999.0
XMI =  99999.0
YMA = -99999.0
YMI =  99999.0
ZMA = -99999.0
ZMI =  99999.0
XB1 = 0.D0
XB2 = 0.D0
YB1 = 0.D0
YB2 = 0.D0
ZB1 = 0.D0
ZB2 = 0.D0

IF(JUMP .NE. 1) THEN
    PRINT *, "Enter Number of boxes for x, y, z-axes"
    READ(5,*) MX, MY, MZ
    PRINT *, "Enter minimum x, y, z coords for cubic grid"
    READ(5,*) XMI, YMI, ZMI
ENDIF

RADIS = BOX/2.0 

! Two passes are made through the protein file.
! 1. Gets statistics including volume and number of residues
! 2. Counts atoms in the grid used to create sphere model

! Read the protein coordinates from file
REWIND PDBFIL

! Read in PDB file until ATOM found in cols 1-4

DO

    READ (PDBFIL,'(A4)', IOSTAT=IERR) ATOM(1:4)
    IF (IERR .LT. 0) THEN
        WRITE(STDERR, "(A)") "End of file reached with no ATOM recods read."
        STOP
    END IF

    IF (ATOM .NE. "ATOM") THEN
        BACKSPACE PDBFIL
        EXIT
    END IF

END DO


! Initialize statistics
NATOMM = 0
NRES = 0
VOLUME = 0.D0

! Track last residue number
LAST_RESN = 0

! Read in all ATOM records
! Calculate total volume, no. residues and maximum and minimum dimensions
! along the x, y and z axes
! Total volume calculated as the sum of residue volumes
DO

    READ(PDBFIL, ATMFMT, IOSTAT=IERR) ATOM, AANAM, RESNUM, X
    IF (IERR .LT. 0) EXIT
    
    IF (ATOM .EQ. "ATOM") THEN
        NATOMM = NATOMM + 1
        XMA = MAX(XMA, X(1))
        XMI = MIN(XMI, X(1))
        YMA = MAX(YMA, X(2))
        YMI = MIN(YMI, X(2))
        ZMA = MAX(ZMA, X(3))
        ZMI = MIN(ZMI, X(3))
        IF (RESNUM .NE. LAST_RESN) THEN
           NRES = NRES + 1
           V = RES_VOLUME(AANAM(1:3))
           IF (V .EQ. -1) THEN
               WRITE(STDERR, ERR1FMT) "Unknown Amino Acid ", AANAM(1:3), "at residue ", NRES, "in", INFILE
               STOP
           END IF
           VOLUME = VOLUME + V
           LAST_RESN = RESNUM
        END IF
    END IF
    
END DO

MX = INT((XMA-XMI)/BOX) + 1
MY = INT((YMA-YMI)/BOX) + 1
MZ = INT((ZMA-ZMI)/BOX) + 1

IF(MX .GT. INCX) THEN
    WRITE(STDERR, ERR2FMT) "Side of box in X direction = ", MX, ", max allowed = ", INCX
    STOP
ELSE IF(MY .GT. INCY) THEN
    WRITE(STDERR, ERR2FMT) "Side of box in Y direction = ", MY, ", max allowed = ", INCY
    STOP
ELSE IF(MZ .GT. INCZ) THEN
    WRITE(STDERR, ERR2FMT) "Side of box in Z direction = ", MZ, ", max allowed = ", INCZ
    STOP
ENDIF

! Re-read file to populate our grid       

REWIND PDBFIL

! Read in PDB file until ATOM found in cols 1-4
DO 

    READ (PDBFIL,'(A4)', IOSTAT=IERR) ATOM(1:4)

    IF (ATOM .NE. "ATOM") THEN
        BACKSPACE PDBFIL
        EXIT
    END IF

END DO

! Zero the total count and the grid square atom counts
ISUM = 0
DO I = 1, INCX
    DO J = 1, INCY
        DO K = 1, INCZ
            INC(I,J,K) = 0
        END DO
    END DO
END DO

DO

    ! Read the PDB and locate each atom record in a grid cube
    READ(PDBFIL, ATMFMT, IOSTAT=IERR) ATOM, AANAM, RESNUM, X
    IF (IERR .LT. 0) EXIT
    
    IF(ATOM .EQ. "ATOM") THEN
    
        ! Locate atom in a grid square along each axis on at a time
    
        DO NX = 1, MX
            XB1 = XMI + BOX*(NX-1)
            XB2 = XMI + BOX*NX
            IF (X(1) .GT. XB1 .AND. X(1) .LT. XB2) THEN
               NXX = NX
               EXIT
            ENDIF
        END DO
        
        DO NY = 1, MY
            YB1 = YMI + BOX*(NY-1)
            YB2 = YMI + BOX*NY
            IF (X(2) .GT. YB1 .AND. X(2) .LT. YB2) THEN
               NYY = NY
               EXIT
            ENDIF
        END DO
        
        DO NZ = 1, MZ
            ZB1 = ZMI + BOX*(NZ-1)
            ZB2 = ZMI + BOX*NZ
            IF (X(3) .GT. ZB1 .AND. X(3) .LT. ZB2) THEN
               NZZ = NZ
               EXIT
            ENDIF
        END DO       
        
        ! Increment atom count in appropriate grid square and the total
        INC(NXX,NYY,NZZ) = INC(NXX,NYY,NZZ) + 1
        ISUM = ISUM + 1        
        
    END IF

END DO

! Add a sphere to all grid spaces with > ICUT atoms in them

IBALL = 0

DO I = 1, MX
    DO J = 1, MY
        DO K = 1, MZ
            IF (INC(I,J,K) .GE. ICUT) THEN
                IF( IPRINT .EQ. 0 .OR. IPRINT .EQ. 1 .OR. IPRINT.EQ.3) THEN
                     XX = BOX*I + XMI
                     YY = BOX*J + YMI
                     ZZ = BOX*K + ZMI
                     WRITE(SPHFIL,"(4F10.2)") XX, YY, ZZ, RADIS
                ENDIF
                IBALL = IBALL + 1
            ENDIF
        END DO
    END DO
END DO

IF (IPRINT .EQ. 1) STOP

VOL = BOX*BOX*BOX*IBALL

IF (IPRINT .EQ. 3) SPHFIL = 6

WRITE(SPHFIL,"(A)") "brktos: This Version 1.20 from 30/10/2013."
WRITE(SPHFIL,"(A,2I10)") " No. Of Atoms ", NATOMM, ISUM
WRITE(SPHFIL,"(A,I7)") " Total Of Amino Acid Residues ", NRES
WRITE(SPHFIL,"(A,F5.1)") " One Side Of The Box        ", BOX
WRITE(SPHFIL,"(A,I5)") " Cutoff For Creating A Ball ", ICUT
WRITE(SPHFIL,"(A,I5)") " No Of Balls                ", IBALL
WRITE(SPHFIL,"(A,F10.2)") " Volume Of Cubes            ", VOL
WRITE(SPHFIL,"(A,F10.2)") " One Side Of The Box        ", VOLUME
WRITE(SPHFIL,"(A,6F7.2)") " Max And Min X,Y,Z Are ", XMA, XMI, YMA, YMI, ZMA, ZMI
WRITE(SPHFIL,"(A,3I4)") " Increments MX, MY, MZ Are ", MX, MY, MZ

STOP

CONTAINS

DOUBLE PRECISION FUNCTION RES_VOLUME(WORD)
    CHARACTER WORD * 3
    DOUBLE PRECISION AV(28)
    CHARACTER NAMAA(28)*3
    DATA NAMAA/"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY", &
        "HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR", &
        "TRP","TYR","VAL","GAL","MAN","NAG","FUC","SIA","NGA", &
        "GLA","GLC"/
    ! Chothis volumes used here (see sluv)
    DATA AV/91.5, 202.1, 135.2, 124.5, 105.6, 161.1, 155.1, 66.4, &
        167.3, 168.8, 167.9, 171.3, 170.8, 203.4, 129.3,  99.1, 122.1, &
        237.6, 203.6, 141.7, 166.8, 170.8, 222.0, 160.8, 326.3, 232.9, &
        155.1, 171.9/
    ! GLA is added with the same value as GLU 
    ! also added GLC as glucose - volume from sluv awa 950313
    ! NB PDB file may need glucose to be changed to GLC 
    ! NGA added (N-acetyl galactosamine), volume from SLUV - MK BOEHM 6/5/97      
    DO I = 1, 28 
         IF (WORD .EQ. NAMAA(I)) THEN
            RES_VOLUME = AV(I)
            RETURN
         ENDIF
    END DO
    
    ! Handle unrecognized residues
    VOLL = -1.
    RETURN
END FUNCTION

END PROGRAM