PROGRAM SCT
IMPLICIT NONE

! Calculate the scattering curve from a sphere model

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

! Original version (Fortran 77):
! S.J. PERKINS   JULY 1981
! Adapted from "PCO" for atomic coordinates
! Rewritten 8.7.1981 for spheres, cleaned up in 1982
! Amended for smearing 16.8.1982
! Smearing and plotting improved 24.8.1982
! Amended 27.8.1982 for direct SINQR/QR,GAUSS NORM, etc.
! Adapted to ICCC CDC 10.3.1983
! Adapted 8.4.84 for output of exp data at end of cal data
! Ported to IRIX 24/8/90 (SILICON GRAPHICS)
! Edited to not cycle through the program - A.S. Nealis 18/8/92

! Updated to Fortran 90 by David W. Wright (November 2013)

CHARACTER * 80 MODEL, EXPDATA, SCATTCV, OUTPUT
CHARACTER * 60 TITLE
INTEGER TODAY(3), NOW(8)
DOUBLE PRECISION XX(3),RSQ(3)
DOUBLE PRECISION X(6000,3)
DOUBLE PRECISION Q(102),SCATT(102),SCATTL(102)
DOUBLE PRECISION ARE(102),SINQR(400)
DOUBLE PRECISION RTERM(400)

DOUBLE PRECISION FFAC(102),SUMM(102)
DOUBLE PRECISION QDIG, DIG, RR, QMAX
DOUBLE PRECISION WAS, WAV, ALF, QR, DELTA, SQ
DOUBLE PRECISION AT, ATS, DENOM, R, W, FAC, QXX, GXX, QXYZ
DOUBLE PRECISION AXX, BXX, CXX, DXX, EXX, FXX
INTEGER IDRSQ(400),I1,LSCAT(102)
INTEGER LINE(102), IS, NAMS, NG, NSC, IX, IY, NP, LP, NLOOP
INTEGER IDEFAU, IWHI, IGUI, ISME, NPT, NXGUI, N
INTEGER NATOM, ISTAR, J, M, KOUNT, NITWIT, NITTT, NNN, ILOOP
INTEGER IMAXIM, IDT, ICHECK, NN, IERR

LOGICAL :: FILEEXISTS

COMMON /F/QMAX,R(7),W(7),IS(10),NAMS(20),NG,NSC

!
! Set parameters
!

DO IX = 1,102
    Q(IX) = 0.0
    SCATT(IX) = 0.0
    SCATTL(IX) = 0.0
    ARE(IX) = 0.0
    FFAC(IX) = 0.0
    SUMM(IX) = 0.0
    LSCAT(IX) = 0
    LINE(IX) = 0
END DO

DO IY = 1,400
    SINQR(IY) = 0.0
    RTERM(IY) = 0.0
    IDRSQ(IY) = 0
END DO

NP=100
LP=1
NLOOP=0

WRITE(*,*) "Enter name of file containing model coords"

READ(*,'(A80)') MODEL(1:80)

WRITE(*,*) "Enter name of file for the calc curve"
READ(*,'(A80)') SCATTCV(1:80)
WRITE(*,*) "Enter name of file containing more info"
READ(*,'(A80)') OUTPUT(1:80)

! STATUS = 'UNKNOWN' will overwrite any existing files
INQUIRE(FILE=MODEL, EXIST=FILEEXISTS)

IF (.NOT. FILEEXISTS) THEN
    WRITE(*,*) 'Specified file does not exist: ', MODEL
    STOP
ELSE
    OPEN(UNIT=1,ACCESS='SEQUENTIAL',FILE=MODEL)
END IF

OPEN(UNIT=2,ACCESS='SEQUENTIAL',FILE=OUTPUT,STATUS='UNKNOWN')
OPEN(UNIT=4,ACCESS='SEQUENTIAL',FILE=SCATTCV,STATUS='UNKNOWN')

WRITE(2,"(A)") " CALCULATION OF DEBYE CURVE ORIGINALLY WRITTEN 8.7.1981;"
WRITE(2,"(A)") " FORTRAN 90 VERSION  NOV 2013"


!TODO: reintroduce dating of files
!CALL IDATE(TODAY)
!CALL ITIME(NOW)

!FRMT=FORMAT ('Date ', i2.2, '/', i2.2, '/', i4.4, '; time', i2.2, ':', i2.2, ':', i2.2 )
!WRITE(2,FMT=FRMT) TODAY(1), TODAY(2), TODAY(3), NOW

WRITE(*,"(A)") " TITLE ? "
READ(*,"(A60)") TITLE

WRITE(2,"(A,A60)") " TITLE % - ", TITLE

WRITE(*,"(A)") " DEFAULT PARAMETERS? Y=1 N=0"
READ (*,*) IDEFAU

NPT = 58
IWHI = 1

IF (IDEFAU .NE. 0) THEN

    IGUI = 1
    ISME = 0
    DIG = 0.4
    QMAX = 0.16
    RR = 3.77

ELSE

    WRITE(*,"(A)") " PLOTS?     Y=1 N=0"
    WRITE(*,"(A)") " SMEARING?     Y=1 N=0"
    READ(*,*) IGUI, ISME

    IF (ISME .NE. 0) THEN
        WRITE(*,"(A)") " WAVELENGTH SPREAD? 0.13? "
        WRITE(*,"(A)") " WAVELENGTH? 5.97?  "
        WRITE(*,"(A)") " BEAM DIVERGENCE? 0.016?  NOT 0 "
        READ(*,*) WAS, WAV, ALF
    END IF

    WRITE(*,"(A)") " R TO BE CHECKED FULLY"
    WRITE(*,"(A)") " INTERVALS OF R IN ANGSTROMS?"
    READ(*,*) DIG

    WRITE(*,"(A)") " MAXIMUM Q VALUE?"
    WRITE(*,"(A)") " RADIUS TO BE LOW AND CHECKED RADIUS OF ONE SPHERE? "
    READ(*,*) QMAX, RR

END IF

!
! Set up Q range and increments
!

QDIG = QMAX / DBLE(NP)

DO N = 1, NP
    Q(N) = (DBLE(N)) * QDIG
END DO

!
! Loop to read in spheres
!

NATOM = 1

DO
    READ(1, FMT="(3F10.2)", IOSTAT=IERR) (X(NATOM,J),J=1,3)

    IF (IERR .LT. 0) THEN
        EXIT
    ELSE IF (IERR .NE. 0) THEN
        WRITE(*,*) 'Abort: Error reading file: ', MODEL
        STOP
    END IF

    NATOM = NATOM + 1
    IF (NATOM .EQ. 6000) THEN
        WRITE(*,*) "NUMBER OF SPHERES TOO GREAT"
        STOP
    END IF

END DO

CLOSE(UNIT=1)

NATOM = NATOM - 1

WRITE(2,*) "COORDINATE READ-IN COMPLETED O.K."

!
! Calculate form factor for spherical subunits
!

DO M=1, 100

    QR = Q(M) * RR
    FAC = (3.0*((DSIN(QR)) - (QR*DCOS(QR))))/ (QR**3)
    FFAC (M) = FAC**2

END DO

WRITE(2,"(/A,F10.3)") "RADIUS OF SINGLE SPHERE IS ", RR

WRITE(2,"(/A)") "SQUARED FORM FACTOR OF SPHERE IN INCREMENTS OF Q"

DO N = 1, 20
    KOUNT = (N-1) * 5
    WRITE (2,"(1X,5F14.9)") (FFAC(KOUNT+J),J=1,5)
END DO

! Calculate lookup table for R

DO N = 1, 400
    RTERM (N) = DBLE(N) * DIG
END DO

! Calculates D(R^2) vs R^2 for all spheres
! Placed into range IDRSQ as function of R
! Note: use of R not R^2
!
! This is the time consuming part
!

NITWIT = NATOM - 1

DO N = 1, NITWIT

    DO J=1,3
        XX(J) = X(N,J)
    END DO

    NITTT = N + 1

    DO NNN = NITTT, NATOM

        DO J = 1, 3
            DELTA = X(NNN,J) - XX(J)
            RSQ(J) = DELTA**2
        END DO

        SQ = RSQ(1) + RSQ(2) + RSQ(3)

        IF (SQ > 0) THEN
            SQ = DSQRT(SQ)
        ELSE
           EXIT
        ENDIF

        ILOOP = IFIX(REAL(SQ/DIG) + 0.5)

        IF (ILOOP .GT. 400) THEN
            WRITE(2,"(A,F15.2)") " R VECTOR LOST - VALUE IS ", SQ
        ELSE
            IDRSQ(ILOOP) = IDRSQ(ILOOP) + 1
        ENDIF

    END DO

END DO

! Sum up total of vectors calculated and maximum range

I1=0
IMAXIM=0
DO N = 1, 400
    I1 = I1 + IDRSQ(N)
    IF(IDRSQ(N) .GT. 0) IMAXIM = N
END DO

! For a given R sum up the scattering factors
! Calculate the scattering function for a giveb Q sum product with R
! and WT(?) it with scattering factors

DO IDT = 1, 100
    SCATT(IDT) = 0.0
END DO

AT = DBLE(NATOM)
ATS = AT**2
ATS = 1.0/ATS


DO M = 1,100
    SUMM(M) = 0.0
    DO N = 1, IMAXIM
        QR = Q(M) * RTERM(N)
        SINQR(N) = DSIN(QR)/ QR
        SUMM(M) = SUMM(M) +  (IDRSQ(N)*SINQR(N))
    END DO

END DO

DO M=1,100
    SCATT(M) = FFAC(M) * ATS * (AT + 2.0 * SUMM(M))*1000.0
END DO

! SCATT contains calculated scattering curve
! Apply smearing if ISME == 1

IF (ISME .EQ. 1) CALL SMEAR(SCATT,QDIG,LINE,SCATTL,Q,WAS,WAV,ALF )

! Output scattering curve arrays

DENOM = SCATT(1)

DO M = 1, 100
    IF (SCATT(M) .LE. 0.0) THEN
        WRITE(2,"(A,I5,3X,4F10.2,F10.5)") " NEGATIVE/ZERO VALUE", M, SCATT(M), SUMM(M), FFAC(M), AT, ATS
        SCATT(M)=0.0
    END IF
    SCATTL(M) = DLOG10(SCATT(M))
    LSCAT(M) = IDINT(SCATT(M)*1000000.0)
    LINE(M) = IDINT(SCATTL(M) * 1000.0)
END DO

! Write summary of ranges


WRITE(2,"(A)") " SCATTERING RANGES SUMMED FOR EACH R"
WRITE(2,"(A,F10.2,A)") " R  STEPS  ARE                 ", DIG, " A   "
WRITE(2,"(A, I7)") " MAXIMUM RANGE IN INCREMENTS IS      ", IMAXIM
WRITE(2,"(A)") " R  INCREMENTS ARE IN UNITS OF  (N)*(DIG) A"
WRITE(2,"(A)") " TABULAR OUTPUT IS R ANGSTROM"
WRITE(2,"(A)") "      SUM OCCURENCE"

DO N=1,40
    KOUNT = (N-1)*10
    WRITE(2,"(1X,10F7.2)")(RTERM(KOUNT+J),J=1,10)
    WRITE(2,"(1X,10I7,/)")(IDRSQ( KOUNT+J ),J=1,10)
END DO

WRITE(2,"(A,I10,/)") " TOTAL OF CROSS-TERMS OF MIXED SPHERES IS ", I1
WRITE(*,"(A,I10,/)") " TOTAL OF CROSS-TERMS OF MIXED SPHERES IS ", I1

ICHECK = ((NATOM**2) - NATOM) / 2

WRITE(2,"(A,I10)") " VALUE OF (NATOM**2-NATOM)/2 IS           ", ICHECK
WRITE(*,"(A,I10)") " VALUE OF (NATOM**2-NATOM)/2 IS           ", ICHECK

! Write summary of SUM(DEBYE)

WRITE(2,"(/,/)")
WRITE(2,"(A,F6.4,A,/)") " SCATTERING CURVE OVER RANGE OF Q = ", QMAX, " ANGSTROM-1"
WRITE(2,"(A,/)") " TABULAR OUTPUT IS Q ANGSTROM-1"
WRITE(2,"(A,/)") "                     DEBYE SUM OVER R"
WRITE(2,"(A,/)") "                     DEBYE CURVE"
WRITE(2,"(A,/)") "                     LOG DEBYE CURVE X 1000"

DO N=1,10
    KOUNT = (N-1)*10
    WRITE(2,"(1X,10F7.4)") (Q(KOUNT+J),J=1,10)
    WRITE(2,"(1X,10F7.0)") (SUMM(KOUNT+J),J=1,10)
    WRITE(2,"(1X,10F7.0)") (SCATT(KOUNT+J),J=1,10)
    WRITE(2,"(1X,10I7,/)") (LINE(KOUNT+J),J=1,10)
END DO

! Output of Q and intensities

DO NN=1, 100
    WRITE(4,"(F9.6,1X,E20.14)") Q(NN), SCATT(NN)
END DO

CLOSE (UNIT=2)
CLOSE (UNIT=4)

STOP

CONTAINS

SUBROUTINE SMEAR(SCATT,QDIG,LINE,SCATTL,Q,WAS,WAV,ALF)
IMPLICIT NONE

! Takes theoretical curve in SCATT and returns it smeared

DOUBLE PRECISION SCATT(102),Q(102),G(301),QQ(301)
DOUBLE PRECISION SIG(301),ASUM(301),F(301),SCATTL(102)
DOUBLE PRECISION QDIG, WAS, WAV, ALF, D, PI, WAB
DOUBLE PRECISION AON, BON, CON, VV, AA, SUMMM, ONE
INTEGER LINE(102), IX, NQ, NQT, IQ, IQ1, IQ2, IQP, M, N , KOUNT
INTEGER J, I

! Set parameters

DO IX=1,301
    G(IX) = 0.0
    QQ(IX) = 0.0
    SIG(IX) = 0.0
    ASUM(IX) = 0.0
    F(IX) = 0.0
END DO

NQ = 100
NQT = 201
D = 179.0
ONE = 1.0
PI = 4.0*DATAN(ONE)

! Calculation of the sigma for the smearing gaussian
! Based on CHAUVIN but 0.25/D^2 removed and 8LN2 not 2LN2
! WAS and ALPHA are from Cusack JMB 1981 145, 539-541

WAB = WAV / (2.0*PI)
CON = (WAS * WAB)**2
CON = CON * 4.0
BON = (ALF)**2
AON = WAB * DSQRT(DFLOAT(8) * DLOG(DFLOAT(2)))

! Calculation of lookup table for sigma (SIG)

DO IQ = 1, 150
    IQ1 = 151-IQ
    IQ2 = 151+IQ
    QQ(IQ1) = -IQ * QDIG
    QQ(IQ2) = IQ * QDIG
    SIG(IQ2) = (DSQRT(CON*QQ(IQ2)*QQ(IQ2) + BON))/AON
    SIG(IQ1) = SIG(IQ2)
END DO

QQ(151) = 0.0
SIG(151) = 2.0 * SIG(152) - SIG(153)

! Apply smearing in G to theoretical F to give smeared SUM
DO IQ=1,100
    F(151 + IQ) = SCATT(IQ)
    F(151 - IQ) = SCATT(IQ)
END DO

! substitution at zero Q is 1000
F(151)=1000.0

! Padding of F beyound the calculated scattering curve (arbitrary)

DO IQ=101,150
    F(151+IQ)=SCATT(100)
    F(151-IQ)=SCATT(100)
END DO

! Convolution of G and F

DO IQ = 1, 301

!Calculation of Gaussian G
    DO IQP=1,301
        !NOTE: Difference from Chauvin and Cusack - IQ not IQP for AA, VV

        AA = 1.0 / (SIG(IQ) * DSQRT(2.0 * PI))
        G(IQP) = 0.0
        VV = (QQ(IQ) - QQ(IQP)) / SIG(IQ)
        VV = VV**2
        IF (VV .LE. 76.0) G(IQP) = AA * EXP(-VV/2.0)

    END DO

    DO IQP = 1, 301
       ASUM(IQ ) = ASUM(IQ ) + ((F(IQP)*G(IQP)))
    END DO

    ASUM(IQ)=ASUM(IQ)*QDIG

END DO

SUMMM=ASUM(151)

! Output for the unsmeared curve, then part of the gaussian
WRITE(2,"(/,A)") " UNSMEARED SCATTERING CURVE"
WRITE(2,"(A)") " TABULAR OUTPUT IS Q ANGSTROM-1 "
WRITE(2,"(A)") "      )DEBYE CURVE  "
WRITE(2,"(A)") "      )LOG DEBYE CURVE X 1000"

DO M=1,100
    SCATTL(M)=DLOG10(SCATT(M))
    LINE(M)=IFIX(REAL(SCATTL(M))*1000.0)
END DO

WRITE(2,"(/)")

DO N = 1,10
    KOUNT = (N-1)*10
    WRITE(2,"(1X,10F7.4)") (Q(KOUNT+J),J=1,10)
    WRITE(2,"(1X,10F7.0)") (SCATT(KOUNT+J),J=1,10)
    WRITE(2,"(1X,10I7,/)") (LINE(KOUNT+J),J=1,10)
END DO

WRITE(2,"(/,A,/)") " SMEARING APPLIED"
WRITE(2,"(A,F8.4)") " WAVELENGTH      ", WAV
WRITE(2,"(A,F8.4)") " WAVELENGTH SPREAD ", WAS
WRITE(2,"(A,F8.4)") " BEAM DIVERGENCE      ", ALF

WRITE(2,"(/,A)") " Q , SIGMA AND GAUSSIAN VALUES AS FOLLOWS"
WRITE(2,"(A,3F7.4)") " THE THREE RANGES REFER TO Q VALUES OF ", QQ(151), QQ(201), QQ(251)

DO I=1,3

    IQ = 101 + (I*50)
    DO IQP = 1,301
    ! Note difference from Chauvin and Cusack - IQ not IQP for AA, VV
       AA = 1.0 / (SIG(IQ ) * DSQRT(2.0 * PI))
       G(IQP) = 0.0
       VV=(QQ(IQ)-QQ(IQP))/SIG(IQ)
       VV=VV*VV
       IF(VV .LE. 76.0) G(IQP)=AA*EXP(-VV/2.0)
    END DO

    IF (I .EQ. 1) WRITE(2,"(1X,F7.4)") QQ (151)
    IF (I .EQ. 1) WRITE(2,"(1X,F7.4)") SIG(151)

    WRITE(6,"(1X,F7.3)") G(151)

    IF (I .EQ. 1) WRITE(2,"(/)")

    DO M = 1,10
        KOUNT = (M-1) * 10 + 151
        IF (I .EQ. 1) WRITE(2,"(1X,10F7.4)") (QQ (J+KOUNT),J=1,10)
        IF (I .EQ. 1) WRITE(2,"(1X,10F7.4)") (SIG(J+KOUNT),J=1,10)
        WRITE(6,"(1X,10F7.3)") (G  (J+KOUNT),J=1,10)
        IF(I .EQ. 1) WRITE(2,"(/)")
    END DO

    WRITE(2,"(/)")

END DO

WRITE(2,"(/,A,F10.3)") " INTENSITY AT ZERO Q %  BEFORE GAUSSIAN =", F(151)
WRITE(2,"(/,A,F10.3)") "               %  AFTER GAUSSIAN  =", SUMMM

DO IQ = 152, 251
    I = IQ - 151
    SCATT(I) = ASUM(IQ)
END DO

RETURN

END SUBROUTINE

END PROGRAM
