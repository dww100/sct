MODULE SJP_UTIL
! Module containing common code for SAS analysis programs

!-------------------------------------------------------------------------------

! Copyright 1981-2014 University College London

! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at

!    http://www.apache.org/licenses/LICENSE-2.0

! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

!-------------------------------------------------------------------------------

! Original author: Stephen J. Perkins
! Rewritten: David W. Wright
! Date: 08 October 2013

IMPLICIT NONE

CONTAINS

SUBROUTINE GET_FILENAME (MSG, FILENAME)

    ! Write a prompt and read a filename from a user

    CHARACTER(LEN=*), INTENT(IN) :: MSG
    CHARACTER(LEN=*), INTENT(OUT) :: FILENAME

    LOGICAL VALID
    INTEGER COUNT

    COUNT = 0

    VALID = .FALSE.
    DO WHILE (.NOT. VALID)

        WRITE(*,'(1X,A)') MSG
        READ (*,'(A)') FILENAME
        IF (FILENAME .NE. ' ') THEN
            VALID = .TRUE.
        ELSE
            WRITE (*,'(1X,A)') 'Error: Please supply a filename'
            COUNT = COUNT +1
        END IF

        IF (COUNT .GT. 3) STOP

    END DO

    RETURN

END SUBROUTINE GET_FILENAME

SUBROUTINE READ_SPHERE_FILE (UNO, FILENAME, COORDS, RADII, NSPHERE)

    ! Read in coordinates and radii from sphere file
    ! UNO specifies unit number to use for file
    ! Assume first three columns = coordinates and 4th = sphere radius
    ! NSPHERE = number of spheres read in (hard limit of 6000)

    INTEGER, INTENT(IN) :: UNO
    CHARACTER(LEN=*), INTENT(IN) :: FILENAME
    REAL, DIMENSION(6000,3), INTENT (OUT) :: COORDS
    REAL, DIMENSION(6000), INTENT (OUT) :: RADII
    INTEGER, INTENT(OUT) :: NSPHERE

    LOGICAL :: FILEEXISTS
    INTEGER :: IERR, J

    INQUIRE(FILE=FILENAME, EXIST=FILEEXISTS)

    IF (.NOT. FILEEXISTS) THEN
        WRITE(*,*) 'Specified file does not exist: ', FILENAME
        STOP
    ELSE
        OPEN(UNIT=UNO, ACCESS='SEQUENTIAL', FILE=FILENAME)

        NSPHERE = 1

        DO
            READ(UNO,FMT="(4F10.2)", IOSTAT=IERR) (COORDS(NSPHERE,J),J=1,3), RADII(NSPHERE)

            IF (IERR .LT. 0) THEN
                EXIT
            ELSE IF (IERR .NE. 0) THEN
                WRITE(*,*) 'Abort: Error reading file: ', FILENAME
                STOP
            END IF

            NSPHERE = NSPHERE + 1
            IF (NSPHERE .EQ. 6000) THEN
                WRITE(*,*) "Error: Number of spheres too great (> 6000)"
                STOP
            END IF
        END DO

    END IF

    CLOSE(UNIT=1)

    NSPHERE = NSPHERE - 1

END SUBROUTINE

SUBROUTINE READ_SCATTER_FILE (FILENAME, QMIN, Q, I, NO)

    ! Read in scattering data
    ! Assume that the file has Q and I as the first two columns
    ! NO = number of values read from file

    CHARACTER(LEN=*), INTENT(IN) :: FILENAME
    DOUBLE PRECISION, INTENT(IN) :: QMIN
    DOUBLE PRECISION, DIMENSION(*), INTENT(INOUT) :: Q , I
    INTEGER, INTENT(OUT) :: NO

    DOUBLE PRECISION :: TMPQ, TMPI
    LOGICAL :: FILEEXISTS
    INTEGER :: IERR, UNIT_NO

    UNIT_NO = 12

    INQUIRE(FILE=FILENAME, EXIST=FILEEXISTS)
    IF (.NOT. FILEEXISTS) THEN
        WRITE(*,*) 'Specified file does not exist: ', FILENAME
        STOP
    ELSE
        OPEN(UNIT=UNIT_NO, FILE=FILENAME, STATUS='OLD', ACTION='READ')
        NO = 1
        DO
            READ(UNIT_NO, *, IOSTAT=IERR) TMPQ, TMPI
            IF (IERR .LT. 0) THEN
                NO = NO - 1
                EXIT
            ELSE IF (IERR .NE. 0) THEN
                WRITE(*,*) 'Abort: Error reading file: ', FILENAME
                STOP
            ELSE IF (TMPQ .GE. QMIN) THEN
                Q(NO) = TMPQ
                I(NO) = TMPI
                NO = NO + 1
            END IF
        END DO

        CLOSE(UNIT_NO)
    END IF

    RETURN

END SUBROUTINE READ_SCATTER_FILE

SUBROUTINE QRANGE_MATCH(QOBS, QCALC, ICALC, XNO, CNO, IMATCH, MNO)

    ! Match modelled I values with experimental Q range

    ! QOBS = Experimental Q values
    ! QCALC, ICALC = Modelled Q and I values
    ! XNO = no. experimental values, CNO = no. modelled values
    ! IMATCH = Modelled I values matched to experimental Q values
    ! MNO = no. matched values
    DOUBLE PRECISION, DIMENSION(*), INTENT(IN) :: QOBS, QCALC, ICALC
    INTEGER, INTENT(IN) :: XNO, CNO
    DOUBLE PRECISION, DIMENSION(*), INTENT(INOUT) :: IMATCH
    INTEGER, INTENT(OUT) :: MNO

    INTEGER :: I, J
    DOUBLE PRECISION :: QMINDIFF, DELTAQ

    ! Find the last experimental Q value that overlaps with the modelled ones
    MNO = 0
    DO I = 1, (XNO)
        IF ( QOBS(I) .LE. QCALC(CNO) ) MNO = MNO + 1
    END DO

    ! Find the corresponding modelled Q values for each experimental one
    ! Store the corresponding I value in IMATCH
    DO I = 1, MNO
        ! Set initial Q minimum difference to a very big number
        QMINDIFF = 1.0E+13
        DO J = 0, (CNO)
            DELTAQ = ABS( QOBS(I) - QCALC(J) )
            IF ( DELTAQ .LT. QMINDIFF ) THEN
                QMINDIFF = DELTAQ
                IMATCH(I) = ICALC(J)
            ENDIF
        END DO
    END DO

    RETURN

END SUBROUTINE QRANGE_MATCH

DOUBLE PRECISION FUNCTION AVERAGE_DATA (DATAIN, N)

    ! Calculate average of values in DATAIN for indices 0 to N
    DOUBLE PRECISION, DIMENSION(*), INTENT(IN) :: DATAIN
    INTEGER, INTENT(IN) :: N

    INTEGER :: I
    DOUBLE PRECISION :: TOTAL

    TOTAL = 0
    DO I = 1, N
        TOTAL = TOTAL + DATAIN(I)
    END DO

    AVERAGE_DATA = TOTAL / N

END FUNCTION AVERAGE_DATA

DOUBLE PRECISION FUNCTION CALC_RFACTOR (QOBS, IOBS, ICALC, N, QMIN, QMAX, CON, VERBOSE)

    ! Calculate the R factor comparing IOBS (experiment) to
    ! ICALC (calculated curve)
    ! R factor is used by analogy with crystallography where:
    ! R = sum (abs(F_obs - F_calc)) / sum (abs(F_obs))

    DOUBLE PRECISION, DIMENSION(*), INTENT(IN) :: QOBS, IOBS, ICALC
    INTEGER, INTENT(IN) :: N
    DOUBLE PRECISION, INTENT(IN) :: QMIN, QMAX
    DOUBLE PRECISION, INTENT(INOUT) :: CON
    LOGICAL, INTENT(IN) :: VERBOSE

    DOUBLE PRECISION :: DELTAC, ORFAC, RFNUM, RFDEN, RFACTOR
    INTEGER :: NDX

    ! Initialize the update and R factor
    DELTAC = CON / 10.
    RFACTOR = 1000000.

    DO WHILE ( ( ABS(DELTAC) ) .GT. ( CON / 10000. ) )

        ORFAC = RFACTOR
        RFACTOR = 0.
        RFNUM = 0.
        RFDEN = 0.

        DO NDX = 1, N
            IF ( ( QOBS(NDX) .LE. QMAX ) .AND. ( QOBS(NDX) .GT. QMIN ) ) THEN
                RFNUM = RFNUM + ABS( IOBS(NDX) - (CON * ICALC(NDX)) )
                RFDEN = RFDEN + ABS( IOBS(NDX) )
            END IF
        END DO

        RFACTOR = RFNUM / RFDEN

        IF (VERBOSE) THEN
            WRITE(*,*) 'BEST YET:', DELTAC, CON, RFACTOR
        END IF

        IF (RFACTOR .LT. ORFAC) THEN
            CON = CON + DELTAC
        ELSE
            DELTAC = DELTAC * (-0.5)
            CON = CON + DELTAC
        ENDIF

    END DO

    IF (VERBOSE) THEN
        WRITE(*,*) 'FINISHED WITH:', DELTAC, CON, RFACTOR
    END IF

    CALC_RFACTOR = RFACTOR * 100

END FUNCTION CALC_RFACTOR

DOUBLE PRECISION FUNCTION CALC_PEARSON (QOBS, IOBS, ICALC, N, QMIN, QMAX, CON, VERBOSE)

    ! Calculate the Pearson Chi squared comparing IOBS (experiment) to
    ! ICALC (calculated curve)
    !
    ! X^2 = sum ((F_calc - F_obs)^2) / F_obs)

    DOUBLE PRECISION, DIMENSION(*), INTENT(IN) :: QOBS, IOBS, ICALC
    INTEGER, INTENT(IN) :: N
    DOUBLE PRECISION, INTENT(IN) :: QMIN, QMAX
    DOUBLE PRECISION, INTENT(INOUT) :: CON
    LOGICAL, INTENT(IN) :: VERBOSE

    DOUBLE PRECISION :: DELTAC, OC2, C2NUM, CHI2
    DOUBLE PRECISION :: SCALEDOBS
    INTEGER :: NDX

    ! Initialize the update and R factor
    DELTAC = CON / 10.
    CHI2 = 1000000.

    DO WHILE ( ( ABS(DELTAC) ) .GT. ( CON / 10000. ) )

        OC2 = CHI2
        CHI2 = 0.
        C2NUM = 0.

        DO NDX = 1, N
            IF ( ( QOBS(NDX) .LE. QMAX ) .AND. ( QOBS(NDX) .GT. QMIN ) ) THEN

                ! Looking to compare I/I(0) to make comparisons between 
                ! different datasets make sense
                SCALEDOBS = IOBS(NDX) / (CON * ICALC(1))
                C2NUM = ( (ICALC(NDX)/ICALC(1)) - SCALEDOBS )**2
                CHI2 = CHI2 + (C2NUM / SCALEDOBS)
            END IF
        END DO

        IF (VERBOSE) THEN
            WRITE(*,*) 'BEST YET:', DELTAC, CON, CHI2
        END IF

        IF (CHI2 .LT. OC2) THEN
            CON = CON + DELTAC
        ELSE
            DELTAC = DELTAC * (-0.5)
            CON = CON + DELTAC
        ENDIF

    END DO

    IF (VERBOSE) THEN
        WRITE(*,*) 'FINISHED WITH:', DELTAC, CON, CHI2
    END IF

    CALC_PEARSON = CHI2

END FUNCTION CALC_PEARSON

DOUBLE PRECISION FUNCTION CALC_CHI2 (QOBS, IOBS, OBSERR, ICALC, N, QMIN, QMAX, CON, VERBOSE)

    ! Calculate the reduced chi-squared statistic comparing IOBS (experiment) to
    ! ICALC (calculated curve)
    !
    ! Chi^2 = 1/(N-1) * sum ((F_obs - F_calc)^2 / OBSERR^2)
    ! Where N = number of observations

    DOUBLE PRECISION, DIMENSION(*), INTENT(IN) :: QOBS, OBSERR, IOBS, ICALC
    INTEGER, INTENT(IN) :: N
    DOUBLE PRECISION, INTENT(IN) :: QMIN, QMAX
    DOUBLE PRECISION, INTENT(INOUT) :: CON
    LOGICAL, INTENT(IN) :: VERBOSE

    DOUBLE PRECISION :: DELTAC, OCHI2, TMPCHI2, CHI2
    DOUBLE PRECISION, DIMENSION(N) :: OBSERR2
    INTEGER :: NDX

    DO NDX = 1, N
        OBSERR2(NDX) = OBSERR(NDX)**2
    END DO

    ! Initialize the update and Chi^2
    DELTAC = CON / 10.
    CHI2 = 1000000.

    DO WHILE ( ( ABS(DELTAC) ) .GT. ( CON / 10000. ) )

        OCHI2 = CHI2
        CHI2 = 0.

        DO NDX = 1, N
            IF ( ( QOBS(NDX) .LE. QMAX ) .AND. ( QOBS(NDX) .GT. QMIN ) ) THEN
                TMPCHI2 = ( (CON * ICALC(NDX)) - IOBS(NDX) )**2
                TMPCHI2 = TMPCHI2 / OBSERR2(NDX)
                CHI2 = CHI2 + TMPCHI2
            END IF
        END DO

        CHI2 = CHI2 / (NDX - 1)

        IF (VERBOSE) THEN
            WRITE(*,*) 'BEST YET:', DELTAC, CON, CHI2
        END IF

        IF (CHI2 .LT. OCHI2) THEN
            CON = CON + DELTAC
        ELSE
            DELTAC = DELTAC * (-0.5)
            CON = CON + DELTAC
        ENDIF

    END DO

    IF (VERBOSE) THEN
        WRITE(*,*) 'FINISHED WITH:', DELTAC, CON, CHI2
    END IF

    CALC_CHI2 = CHI2

END FUNCTION CALC_CHI2

END MODULE SJP_UTIL
