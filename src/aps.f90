PROGRAM APS

! Program to calculate the radius of gyration of sphere models
! Hard limit of 5000 spheres

! Original Fortran 77 code:
! S.J. Perkins July 1981
! Cleaned up by SJP 9 August 1985
! Converted to use with SGI IRIX 3.3.2 by A.S. Nealis 2 Oct 1991

! Rewritten in Fortran 90 by David W. Wright November 2013

! TODO: For reasons currently beyond my comprehension this results in X having a
! start index of 0 and D of 1, this needs to be fixed as it is odd
! and makes the code difficult to understand
DIMENSION X(5000,3),D(5000)

WRITE(6,"(A/)") "APS: VERSION 2.50"
WRITE(6,"(A/)") "CALCULATION OF Rg FOR UP TO 5000 SPHERE MODELS"

!
!  Read in sphere coordinates and radii (also count the number of spheres)
!

REWIND 5

NSPHERE = 0

DO
    READ(5,"(4F10.2)", IOSTAT=IERR) (X(NSPHERE,J),J=1,3), RADIUS
    IF (IERR .LT. 0) THEN
        EXIT
    ELSE IF (IERR .NE. 0) THEN
            WRITE(*,*) 'Abort: Error reading file: ', FILENAME
            STOP
    END IF
    NSPHERE = NSPHERE + 1
END DO

CLOSE(5)


WRITE(6,"(A,I7/)") "COORDINATE READ-IN DONE : TOTAL IS ", NSPHERE

!
! Calculate Rg and volume of model boxes and spheres and total volumes
!

! Sphere radius r: Rg^2 = 3/5 * r^2 -> Rg = (3/5)^1/2 * r

SPHERE_RG = SQRT(3.0 / 5.0) * RADIUS
! Cube side l: Rg^2 = l^2 / 4 -> Rg = l/2
! Here l = 2 * r -> Rg = r -> Rg^2 = r^2
BOX_RG2 = RADIUS**2

PI = 3.1415927
RADIUS3 = RADIUS**3

SPHERE_VOLUME = (4.0 * PI * RADIUS3) / 3.0
VOLUME_MODEL = SPHERE_VOLUME * REAL(NSPHERE)

! Box value = (2 * r)^3 = 8 * r^3
BOX_VOLUME = 8.0 * RADIUS3
TOTAL_BOX_VOLUME = BOX_VOLUME * NSPHERE

WRITE(6,"(A,F5.2)") "SPHERE RADIUS (HALF CUBE SIDE): ", RADIUS
WRITE(6,"(A,F5.2)") "RADIUS OF GYRATION OF SPHERE: ", SPHERE_RG
WRITE(6,"(A,F5.2)") "VOLUME OF SINGLE SPHERE = ", SPHERE_VOLUME
WRITE(6,"(A,F13.2)") "TOTAL VOLUME OF ALL SPHERES   = ", VOLUME_MODEL
WRITE(6,"(A,F5.2)") "RADIUS OF GYRATION OF CUBE: ", RADIUS
WRITE(6,"(A,F7.2)") "VOLUME OF SINGLE CUBE = ", BOX_VOLUME
WRITE(6,"(A,F13.2/)") "TOTAL VOLUME OF ALL BOXES   = ", TOTAL_BOX_VOLUME

!
! Calculation of mean X Y Z coordinates & printout
!

XSUM = 0.0
YSUM = 0.0
ZSUM = 0.0

DO N = 0, NSPHERE - 1
    XSUM = XSUM + X(N,1)
    YSUM = YSUM + X(N,2)
    ZSUM = ZSUM + X(N,3)
END DO

WRITE(6,"(A,I8)") "NO OF COORDINATES = ", NSPHERE
WRITE(6,"(A,3F10.3)") "SUM OF ALL X,Y,Z COORDS = ", XSUM, YSUM, ZSUM

XMEAN = XSUM / NSPHERE
YMEAN = YSUM / NSPHERE
ZMEAN = ZSUM / NSPHERE

WRITE(6,"(A,3F10.3)")"MEAN X,Y,Z COORDS ARE ", XMEAN, YMEAN, ZMEAN

!
! Calculate Rg for the model
!

RADSQ = 0.0
RRAD = 0.0

DO N = 0, NSPHERE - 1

    D(N+1) = ((X(N,1)-XMEAN)**2 + (X(N,2)-YMEAN)**2 + (X(N,3)-ZMEAN)**2)
    RADSQ = RADSQ + D(N+1)
    ARAD = SQRT(D(N+1))
    RRAD=RRAD + ARAD

END DO

UNIT_RADIUS = RRAD / NSPHERE
WRITE(6, "(A, F19.3)") "RADIUS OF UNIT WEIGHT= ", UNIT_RADIUS
RGYR = SQRT( (RADSQ / NSPHERE) + BOX_RG2)

WRITE(6,"(A,F19.3)") "SUM OF RADII SQ      = ", RADSQ
WRITE(6,"(A,F19.3)") "RADIUS OF GYRATION   = ", RGYR

STOP

END PROGRAM
