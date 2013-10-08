Program RfacXN

! RfacXn

! Takes two small angle scattering curves and compares them by calculating the R factor.
! First curve provided is assumed to be a computationally calculated curve, 
! the second experimental in origin.
! Minimimum and maximum Q values to consider must be provided.

! Original author: Stephen J. Perkins
! Rewritten: David W. Wright
! Date: 08 October 2013

USE SJP_UTIL
IMPLICIT NONE

INTEGER XNO, CNO, MATCHNO, IOERR
REAL    ICALC(0:512), IOBS(0:512), IMATCH(0:512)
REAL    QCALC(0:512), QOBS(0:512)
REAL    QMIN, AVIXPT, AVICAL, RFACTOR
REAL    CON, QMAX
CHARACTER*50    INFIL1, INFIL2, OUTFIL
LOGICAL VERBOSE

VERBOSE = .TRUE.

WRITE (*,'(1X,A)') 'Program Rfac version 1.50'
WRITE (*,*)

! Get the filnames for the input calculated and experimental data 
CALL GET_FILENAME ('Provide the name of the file containing the calculated data: ', INFIL1)
CALL GET_FILENAME ('Provide the name of the file containing the experimental data: ', INFIL2)

WRITE(*,*) 'Minimum value of Q to be used in fit?'
READ(*,*)  QMIN
	  
CALL READ_SCATTER_FILE (INFIL1, QMIN, QCALC, ICALC, CNO)
CALL READ_SCATTER_FILE (INFIL2, QMIN, QOBS, IOBS, XNO)

! Read in the output filename and open
! Note: data is appended to the output file
CALL GET_FILENAME ('Provide the output filename: ', OUTFIL)

OPEN(UNIT=11, FILE=OUTFIL, IOSTAT=IOERR, STATUS='UNKNOWN', ACCESS='APPEND')
IF (IOERR .NE. 0) THEN
    WRITE(*,*) 'Abort: Error opening file: ', OUTFIL
    STOP
END IF
      
! Match the observed (experimental) and calculated values to a common Q range
CALL QRANGE_MATCH(QOBS, QCALC, ICALC, XNO, CNO, IMATCH, MATCHNO)

! Calculate averages of the experimental and calculated values over the matched Q range
AVIXPT = AVERAGE_DATA(IOBS, MATCHNO)
AVICAL = AVERAGE_DATA(IMATCH, MATCHNO)

IF (VERBOSE) THEN
    WRITE(*,*) 'Mean Observed=', AVIXPT
    WRITE(*,*) 'Mean Calculated=', AVICAL
    WRITE(*,*) 'Initial Guess=', AVIXPT/AVICAL
    WRITE(*,*) 'Max Xptl Q=', QOBS(XNO - 1)
    WRITE(*,*) 'Max Calc Q=', QCALC(CNO - 1)
END IF
       
WRITE(*,*) 'What is the maximum Q value to be used in the fit?'
READ(*,*)  QMAX

! Use the ratio of the experimental and calculated curves to give a initial guess 
! of the concentration       
CON = AVIXPT / AVICAL
RFACTOR = CALC_RFACTOR (QOBS, IOBS, IMATCH, MATCHNO, QMIN, QMAX, CON, VERBOSE)

WRITE(11,*) INFIL1, QMAX, CON, RFACTOR
CLOSE(11)

END PROGRAM
