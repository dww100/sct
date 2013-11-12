PROGRAM FRCC
! Calculate frictional coefficient from sphere models of proteins
IMPLICIT NONE

! Original version of frcc by S. J. PERKINS July 1983
! Following Bloomfield's Biopolymers 1977 series
! and his program supplied as "GENDIA"
! uses the diagonal approximation
! Edited by K.F. Smith 1991-1992
! IRIX 6.5 VERSION by SJP for 2000 spheres

! Fortran 90 version by David W. Wright November 2013

INTEGER MAX_SPH
PARAMETER (MAX_SPH = 2000)
DOUBLE PRECISION R(MAX_SPH,3),G(MAX_SPH),T(MAX_SPH,MAX_SPH,3)
DOUBLE PRECISION E(MAX_SPH)
DOUBLE PRECISION DIS(3),EIG(3)
DOUBLE PRECISION PI, EI, ER, EPS, S1, S2, S3, S4, S5, E2, E3
DOUBLE PRECISION EL2, D2, DISTAN, EE, FACT, A, B, F1, FR1, FR2
DOUBLE PRECISION V, F2, EA, F, FRT
INTEGER IND, NMAX, IP, NEPS, NSTOP, N, J, NSPHERE, K
INTEGER L, I, IT, IERR
LOGICAL TEST

PI = 3.1415926536
IND = 1
NMAX = 100
IP = 3
NEPS = 4

WRITE(6,"(A,/)") " FRCC: THIS VERSION 11.11.2013"

DO N=1,2000
    READ(5, "(4F10.2)", IOSTAT=IERR) (R(N,J),J=1,3),E(N)

    IF (IERR .LT. 0) THEN
        EXIT
    ELSE IF (IERR .NE. 0) THEN
         WRITE(*,*) 'Abort: Error reading input'
         STOP
    END IF

END DO

NSPHERE = N - 1
WRITE(6,"(/,A,I5,/)") " TOTAL SPHERES = ", NSPHERE

EPS = 10.0**(-NEPS)
S1 = 0.0
S2 = 0.0
S3 = 0.0
S4 = 0.0
S5 = 0.0

DO K = 1, NSPHERE

    E2 = E(K)**2
    S1 = S1 + E(K)
    S3 = S3 + E2
    E3 = E2 * E(K)
    S4 = S4 + E3

    DO L = 1, NSPHERE
        EL2 = E(L)**2
        IF (K .LT. L) THEN
            D2 = 0.0
            DO I = 1, 3
                DIS(I) = R(L,I) - R(K,I)
                D2 = D2 + DIS(I)**2
            END DO
            DISTAN = SQRT(D2)
            EE = E(K) * E(L)
            S2 = S2 + EE / DISTAN
            S5 = S5 + EE * EE / DISTAN
            FACT = FLOAT(IND) * (E2 + EL2) / D2
            A = 1.0 + FACT / 3.0
            B = 1.0 - FACT

            ! Calculation of the diagonal part of the Oseen tensor only
            DO I=1,IP
                V=A+(DIS(I)*DIS(I)*B/D2)
                T(K,L,I)=V*6.0/(8.0*DISTAN)
                T(L,K,I)=T(K,L,I)
            END DO
        END IF

    END DO
END DO

S2 = 2.0 * S2
S5 = 2.0 * S5
F1 = S1 / (1.0 + S2 / S1)
F2 = S3**2 / (S4 + S5)
FR1 = 6.0 * PI * 0.01002 * F1 * 1.0E-8
FR2 = 6.0 * PI * 0.01002 * F2 * 1.0E-8

WRITE(6,*)
WRITE(6,*) " KIRKWOOD"
WRITE(6,*)
WRITE(6,"(A,F9.4,E12.4)") " FRICTIONAL COEFFICIENT = ", F1, FR1
WRITE(6,*)
WRITE(6,*) " BLOOMFIELD "
WRITE(6,*)
WRITE(6,"(A,F9.4,E12.4)") " FRICTIONAL COEFFICIENT = ", F2, FR2
WRITE(6,*)
WRITE(6,"(A, I5)") " TOTAL SPHERES = ", NSPHERE

! ITERATIONS START HERE

DO I = 1, NSPHERE
    DO J = 1, NSPHERE
        IF (I .NE. J) THEN
            DO K = 1, IP
                T(I,J,K) = T(I,J,K) * E(J)
            END DO
        END IF
    END DO
END DO

DO K = 1, IP

    WRITE(6,"(A,I2,/)") " COMPONENT ", K

    DO I = 1, NSPHERE
        G(I) = 1.0
    END DO

    EA=0.0

    TEST = .TRUE.
    IT = 1

    DO WHILE (TEST)

        EIG(K)=0.0

        DO I = 1, NSPHERE
            G(I)=1.0
            DO J=1,NSPHERE
                IF (I .NE. J) G(I)=G(I)-T(I,J,K)*G(J)
            END DO

            EIG(K)=EIG(K) + E(I)*G(I)

        END DO

        WRITE(6,"(I5,E11.4)") IT, EIG(K)

        IF (ABS((EIG(K) - EA) / EIG(K)) .GE. EPS) THEN
            EA = EIG(K)
        ELSE
            EXIT
        END IF

        IF (IT .EQ. NMAX) THEN
            WRITE(9,"(A)") " MAXIMUM OF ITERATIONS REACHED"
            EXIT
        END IF

        IT = IT + 1

    END DO


END DO

F = 3.0 / (1.0 / EIG(1) +1.0 / EIG(2) + 1.0 / EIG(3))
FRT = 6.0 * PI * 0.01002 * F * 1.0E-8
WRITE(6,"(A,E10.4)") " FRICTIONAL COEFFICIENT (DIA APPROX)= ", F

WRITE(6,"(A,E10.4)") " PROPER FRICTIONAL COEFFICIENT      = ", FRT


END PROGRAM
