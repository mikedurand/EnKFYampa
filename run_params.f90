SUBROUTINE RUN_PARAMS(N_YR,N_T,N_M,N_YP,N_Y,N_YS,N_U,N_R,N_ZR,N_ZP,N_Z,&
  MEASTIMES,CU,CV,LAM_U,ISTART,IEND,NPROCS,N_A,ALPHA_BAR,ZUSE,&
  N_F,FREQ,THETA,MSTART,MEND,C_ALPHA,LAM_ALPHA,PRECIP_SCALE,&
  MEAS_SWITCH,GEN_SWITCH,SPATIAL_PIC,CVEG,LAM_VEG,&
  INVERT_PRECIP,RANK,ATM_SWITCH,RANGE_U,RANGE_A,GEN_SEEDS,ENS_NUM,LAM_B,C_B)

IMPLICIT NONE

INTEGER,INTENT(IN)::NPROCS,RANK
INTEGER,INTENT(IN)::N_YR,N_T,N_M,N_YP,N_Y,N_YS,N_U,N_R,N_ZR,N_Z,N_A,N_F
INTEGER,INTENT(IN) :: N_ZP(N_ZR)
LOGICAL,INTENT(IN)::GEN_SWITCH
INTEGER,INTENT(OUT) :: MEASTIMES(N_M),ZUSE(N_Z),ISTART(NPROCS),IEND(NPROCS),&
  MSTART(NPROCS),MEND(NPROCS),ATM_SWITCH,GEN_SEEDS(5)
INTEGER I,J,P,NEXTRA,NBASE,N_DAYS_STRT,K,ENS_NUM,M,mpi_err
REAL :: CV1,CV2,CV3
REAL :: DIST(N_YP,N_YP),COV_U(N_U),VAR_U(N_U),COV_ALPHA(N_A),VAR_ALPHA(N_A),&
  COV_B,COV_VEG,VAR_VEG,VAR_B,MU_VEG,COV_RADT,VAR_RADT,COV_CHI_P,VAR_CHI_P
REAL :: ABLT_PRCP_SCALE
INTEGER,INTENT(OUT) :: MEAS_SWITCH,SPATIAL_PIC
REAL,INTENT(OUT) :: CU(N_YP,N_YP,N_U),CV(N_Z,N_Z),LAM_U(N_YP,N_U),&
  ALPHA_BAR(N_A),FREQ(N_F),THETA(N_F),C_ALPHA(N_YP,N_YP,N_A),&
  LAM_ALPHA(N_YP,N_A),PRECIP_SCALE,CVEG,LAM_VEG, INVERT_PRECIP,&
  RANGE_U,RANGE_A,LAM_B,C_B
 ! LAM_RADT,C_RADT,PSCALE,RANGE_B,RANGE_RADT,RANGE_CHI_P,C_CHI_P,&
 ! LAM_CHI_P
  
CHARACTER(32) N_UC,N_AC,N_YPC

! 0) READ RUN PAMARETER FILE
!    THIS IS NOW SET UP FOR SPECIFIC FORMATS FOR EACH LINE TO AVOID MISREADS.
!    IN ORDER TO READ A CERTAIN NUMBER OF NUMBERS ON A SINGLE LINE, I HAVE 
!    SET UP A CHARACTER CONCATENATION OF N_UC AND N_AC WITH THE REST OF THE
!    FORMAT STATEMENT FOR SOME OF THE LINES.  N_UC AND N_AC ARE THE CHARACTER
!    EQUIVALENT OF THE INTEGERS N_U AND N_A.  SEE JOURNAL, 6/28/06

OPEN(UNIT=1,FILE='run_params.in',STATUS='OLD')
!FACTOR TO SCALE PRECIPITATION: 0.25 FOR DEPTH ~2M, 0.125 FOR DEPTH ~0.5 M
READ(1,*)
READ(1,'(F4.2)') PRECIP_SCALE

!SWITCH TO SET WHICH MEAS. ARE USED: 0 FOR ALL MEAS, 1 FOR PM ONLY, 2 FOR PM+Alb
READ(1,*) 
READ(1,'(I2)') MEAS_SWITCH

!NUMBER OF DAYS TO WAIT BEFORE BEGINNING UPDATING
READ(1,*) 
READ(1,'(I3)') N_DAYS_STRT

!COVARIANCE VALUES FOR MEASUREMENT ERROR (NUMBER OF ELEMENTS MUST MATCH N_ZR IN SIZES.IN)
READ(1,*) 
READ(1,'(3(F6.4,1X))') CV1,CV2,CV3

!READ RANGE VALUES FOR COMPUTING SPATIAL COVARIANCE MATRIX
READ(1,*) 
READ(1,'(2(F6.2,1X))') RANGE_U,RANGE_A!,RANGE_B,RANGE_RADT,RANGE_CHI_P

!COEFF. OF VAR. VALUES FOR ENSEMBLE FORCING (NUMBER ELEMENTS MUST MATCH N_U)
WRITE(N_UC,*) N_U !FORM A CHARACTER VARIABLE WITH VALUE OF N_U
READ(1,*) 
READ(1,'('// N_UC // '(F4.2,1X))') (COV_U(I),I=1,N_U)

!COEFFICIENT OF VARIATION VALUE FOR INPTERPOLATION / DISAGGREGATION WEIGHTS B
READ(1,*) 
READ(1,'(F4.2)') COV_B

!COEFFICIENT OF VARIATION VALUE FOR ABSORBED RADIATION
READ(1,*) 
READ(1,'(F4.2)') COV_RADT

!COEFFICIENT OF VARIATION VALUE FOR SUBGRID PRECIP. COEFFICIENT OF VARIATION
READ(1,*) 
READ(1,'(F4.2)') COV_CHI_P

!COEFFICIENT OF VARIATION VALUES FOR ENSEMBLE PARAMETERS (MUST MATCH N_A)
WRITE(N_AC,*) N_A !FORM A CHARACTER VARIABLE WITH VALUE OF N_A
READ(1,*) 
READ(1,'('// N_AC // '(F4.2,1X))') (COV_ALPHA(I),I=1,N_A)

!MEAN MODEL PARAMETERS
READ(1,*) 
READ(1,'('// N_AC // '(E8.1,1X))') (ALPHA_BAR(I),I=1,N_A)

!COEFFICIENT OF VARIATION AND BIAS VALUES FOR ENSEMBLE VEGETATION
READ(1,*) 
READ(1,'(2(F4.1,1X))') COV_VEG, MU_VEG

!DISTRIBUTED FORCING (0) OR FORCING PIXEL NUM. FOR SPATIALLY UNIFORM CASE (1-9)
READ(1,*) 
READ(1,'(F4.2)') ABLT_PRCP_SCALE
!FACTOR TO MULTIPLY PRECIPITATION AFTER FIRST MEASUREMENT TIME
READ(1,*) 
READ(1,'(I1)') SPATIAL_PIC
!VARIABLE TO MULTIPLY PRECIPITATION LAPSE RATE: -1 TO INVERT TRUTH, 1 OTHERWISE
READ(1,*) 
READ(1,'(F4.1)') INVERT_PRECIP
!ATMOSPHERIC MODEL SWITCH: 1~ON, 0~OFF
READ(1,*) 
READ(1,'(I1)') ATM_SWITCH
!GENERATOR RANDOM NUMBER SEEDS (PRECIP,GRAIN GROWTH, NEW GRAIN)
READ(1,*) 
READ(1,'(5(I5,1X))') (GEN_SEEDS(I),I=1,5)
!ENSENBLE NUMBER
READ(1,*) 
READ(1,'(I5)') ENS_NUM
CLOSE(1)

! 1) MEASUREMENT TIMES: 

!! OPTION A: ASSUME A SATELLITE OVERPASS TIME OF 1 AM, 1 PM.
!!    THE '7' VALUE BELOW IS DUE TO THE FACT THAT THE DATASET BEGINS AT MIDNIGHT
!!    UTC TIME WHICH CORRESPONDS TO 6PM MOUNTAIN TIME. USED IN THE LSOS WORK, 
!!    AND IN THE BAYESIAN RECONSTRUCTION
!J=0
!DO I=7+N_DAYS_STRT*24,N_T,12 !Wait N_DAYS_STARTS before beginning updates
!  J=J+1
!  IF(J.LE.N_M)THEN
!    MEASTIMES(J)=I
!  END IF
!END DO

! OPTION B: READ MEASTIMES FROM FILE: NOT USED IN THE FILTER RUNS, THUS FAR
OPEN(UNIT=3,FILE='meas_times.in',STATUS='OLD')
DO M=1,N_M
  READ(3,'(I12)') MEASTIMES(M)
END DO
CLOSE(3)

! 2) DEFINE MEASUREMENT ERROR COVARIANCES
ZUSE=1
!CV1   !PM BRIGHTNESS TEMPS
!CV2   !ALBEDO
!CV3   !SURFACE TEMP

!old way... no inter-pixel correlation
DO I=1,N_Z
  DO J=1,N_Z
    IF (I.EQ.J) THEN
      IF (I.LE.N_ZP(1)) THEN
        CV(I,J)=CV1
      ELSEIF (I.LE.(N_ZP(2)+N_ZP(1))) THEN
        CV(I,J)=CV2
      ELSE
        CV(I,J)=CV3
      END IF
    ELSE
      CV(I,J)=0.0
    END IF
  END DO
END DO

! 3) READ DISTANCE MATRIX FROM FILE
OPEN(UNIT=2,FILE='dist.in',STATUS='OLD')
WRITE(N_YPC,*) N_YP !FORM A CHARACTER VARIABLE WITH VALUE OF N_YP
DO I=1,N_YP
  READ(2,'('// N_YPC // '(F10.4))') (DIST(I,J),J=1,N_YP)
END DO
CLOSE(2)

! 4) DEFINE FORCING ERROR PROPERTIES
! 4.1 FORCING
! COMPUTE VARIANCE FOR NORMAL TRANSFORMATION OF LOGNORMAL VARIABLE
VAR_U=LOG(1+COV_U**2)
! COMPUTE MEAN FOR NORMAL TRANSFORMATION OF LOGNORMAL VARIABLE. NOTE THAT
!   LAMBDA MUST HAVE A SPATIAL DIMENSION, EVEN THOUGH WE ARE NOT MODELING
!   A VARIABLE SPATIAL MEAN, IN ORDER FOR DIMENSIONS TO MATCH UP IN MVNRND.
DO I=1,N_U
  DO P=1,N_YP
    LAM_U(P,I)=-0.5*VAR_U(I) !ASSUMES MEAN OF ONE; COULD ADD 'LOG(1)' 
  END DO
END DO
! COMPUTE SPATIAL COVARIANCE FOR NORMAL TRANSFORMATION OF LOGNORMAL VARIABLE.
!   THIS CALCULATION ASSUMES AN EXPONENTIAL VARIOGRAM AND ASSIGNS THE 
!   COVARIANCE BASED ON THE RANGE OF THE VARIOGRAM AND VARIOGRAM-COVARIANCE
!   RELATIONS FROM ISAAKS AND SRIVASTAVA (1986)

DO I=1,N_U
  DO J=1,N_YP
    DO K=1,N_YP
      CU(K,J,I)=VAR_U(I)*EXP(-3.*DIST(K,J)/RANGE_U)
    END DO
  END DO
END DO

! 4.2 PARAMETERS

! GRAIN GROWTH PARAMETERS
!ALPHA_BAR(1)=PARAMETER GOVERNING DRY SNOW GRAIN GROWTH, LAYER 1
!ALPHA_BAR(2)=PARAMETER GOVERNING DRY SNOW GRAIN GROWTH, LAYER 2
!ALPHA_BAR(3)=PARAMETER GOVERNING DRY SNOW GRAIN GROWTH, LAYER 3
!ALPHA_BAR(4)=PARAMETER GOVERNING WET SNOW GRAIN GROWTH
!ALPHA_BAR(5)=ADDITIVE PARAMETER GOVERNING NEW SNOW GRAIN SIZE CALCS
!ALPHA_BAR(6)=EXPONENTIAL PARAM. GOVERNING NEW SNOW GRAIN SIZE CALCS
!ALPHA_BAR(7)=PRECIPITATION LAPSE RATE PARAMETER
!ALPHA_BAR(8)=TEMPERATURE LAPSE RATE PRAMETER DEGREES / METER

! COMPUTE VARIANCE FOR NORMAL TRANSFORMATION OF LOGNORMAL VARIABLE
VAR_ALPHA=LOG(1.0+COV_ALPHA**2)

! COMPUTE MEAN FOR NORMAL TRANSFORMATION OF LOGNORMAL VARIABLE. NOTE THAT
!   LAMBDA MUST HAVE A SPATIAL DIMENSION, EVEN THOUGH WE ARE NOT MODELING
!   A VARIABLE SPATIAL MEAN, IN ORDER FOR DIMENSIONS TO MATCH UP IN MVNRND.
DO I=1,N_A
  DO P=1,N_YP
    LAM_ALPHA(P,I)=-0.5*VAR_ALPHA(I) !ASSUMES MEAN OF ONE; COULD ADD 'LOG(1)' 
  END DO
END DO

! COMPUTE SPATIAL COVARIANCE FOR NORMAL TRANSFORMATION OF LOGNORMAL VARIABLE.
!   THIS CALCULATION ASSUMES AN EXPONENTIAL VARIOGRAM AND ASSIGNS THE 
!   COVARIANCE BASED ON THE RANGE OF THE VARIOGRAM AND VARIOGRAM-COVARIANCE
!   RELATIONS FROM ISAAKS AND SRIVASTAVA (1986)

DO I=1,N_A
  DO K=1,N_YP
    DO J=1,N_YP
      C_ALPHA(K,J,I)=VAR_ALPHA(I)*EXP(-3*DIST(K,J)/RANGE_A)
      !DON'T ALLOW SPATIAL VARIABILITY IN PRECIPITATION LAPSE RATE
      IF (I.EQ.5) C_ALPHA(K,J,I)=VAR_ALPHA(I)
    END DO
  END DO
END DO

! 4.3 VEGETATION 
VAR_VEG=LOG(1+COV_VEG**2)
CVEG=VAR_VEG !THESE ARE IDENTICAL FOR 1-D
LAM_VEG=-0.5*VAR_VEG !ASSUMES MEAN OF ONE; COULD ADD 'LOG(1)' 

! 4.4 INTERPOLATION WEIGHTS FOR FORCING DATA : B
VAR_B=LOG(1+COV_B**2)
C_B=VAR_B !THESE ARE IDENTICAL FOR 1-D
LAM_B=-0.5*VAR_B !ASSUMES MEAN OF ONE; COULD ADD 'LOG(1)' 

!!! 4.5 FACTOR TO PERTURB ABSORBED RADIATION IN SSIB
!VAR_RADT=LOG(1+COV_RADT**2)
!C_RADT=VAR_RADT !THESE ARE IDENTICAL FOR 1-D
!LAM_RADT=-0.5*VAR_RADT !ASSUMES MEAN OF ONE; COULD ADD 'LOG(1)' 

! 4.6 FACTOR TO PERTURB ABSORBED RADIATION IN SSIB
!VAR_CHI_P=LOG(1+COV_CHI_P**2)
!C_CHI_P=VAR_CHI_P !THESE ARE IDENTICAL FOR 1-D
!LAM_CHI_P=-0.5*VAR_CHI_P !ASSUMES MEAN OF ONE; COULD ADD 'LOG(1)' 

! 5) CALCULATE ISTART AND IEND: INDECES TO SUBDIVIDE PIXELS
! FOR ALGORITHM EXPLANATION, SEE 'dIVISION OF pIXELS.SXW' IN rEPORTS

NEXTRA=MOD(N_YP,NPROCS)
NBASE=N_YP/NPROCS
DO I=1,NPROCS
  ! ASSIGN ISTART
  IF(I.EQ.1)THEN
    ISTART(I)=1
  ELSE 
    ISTART(I)=IEND(I-1)+1
  END IF

  ! ASSIGN IEND
  IEND(I)=ISTART(I)+NBASE-1
  IF(I.LE.NEXTRA)THEN
    IEND(I)=IEND(I)+1
  END IF
END DO 


! 6) CALCULATE MSTART AND MEND: INDECES TO SUBDIVIDE MEASUREMENTS
! FOR ALGORITHM EXPLANATION, SEE CALCULATION OF ISTART AND IEND

NEXTRA=MOD(N_M,NPROCS)
NBASE=N_M/NPROCS
DO I=1,NPROCS
  ! ASSIGN ISTART
  IF(I.EQ.1)THEN
    MSTART(I)=1
  ELSE 
    MSTART(I)=MEND(I-1)+1
  END IF

  ! ASSIGN IEND
  MEND(I)=MSTART(I)+NBASE-1
  IF(I.LE.NEXTRA)THEN
    MEND(I)=MEND(I)+1
  END IF
END DO 

! 7) DEFINE FREQUENCIES AND ANGLES FOR TEST
FREQ=(/6.925,10.65,18.7,23.8,36.5,89./)
THETA=55.

END SUBROUTINE RUN_PARAMS
