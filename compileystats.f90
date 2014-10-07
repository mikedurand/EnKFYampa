! ----------------------------------------------------------------------
!
  SUBROUTINE COMPILEYSTATS(Y,YSTATS,N_R,N_TY,N_Y,N_YS,N_RL,RANK,MY_COMM,ERR)
!
! ----------------------------------------------------------------------
!
! SUBROUTINE TO COMPUTE MEAN, MAX, MIN AND VARIANCE OF ARRAYS.  VARIANCE IS 
! COMPUTED BY EXPLICITLY SUMMING THE SQUARE DEVIANCE FROM THE MEAN ACROSS 
! THE ARRAY USING A DO-LOOP AND NORMALIZING.  THE MEAN IS COMPUTED BY 
! FINDING A LOCAL SUM FOR ALL THE REPLICATES ON THIS PROCESS, THEN REDUCING
! THE LOCAL SUMS ACROSS ALL OF THE PROCESSES.  THE MAX AND MIN ARE COMPUTED
! BY FINDING MAX AND MIN, THEN REDUCING.  THE VARIANCE IS COMPUTED BY 
! FINDING THE LOCAL SUM OF SQUARED DEVIANCE FROM THE MEAN, THEN SUMMING THE
! SUM OF SQUARED DEVIANCE ACROSS ALL THE PROCESSES USING REDUCE.
!
! Y - INPUT ARRAY FOR WHICH WE WANT STATISTICS
! YSTATS - OUTPUT ARRAY CONTAINING STATISTICS
! N_R - NUMBER OF REPLICATES
! N_TY - LENGTH OF Y VECTOR IN TIME DIMENSION
! N_Y - TOTAL NUMBER OF STATES
! N_YS - NUMBER OF Y STATISTICS TO BE RECORDED
! MY_COMM - LOCAL COPY OF GLOBAL MPI_COMM_WORLD COMMUNICATOR
! ERR - GLOBAL ERROR VARIABLE 
! N_RL - NUMBER OF REPLICATES ON EACH PROCESS
!
! CALLS: NONE
! CALLED BY: ENKF
! 
! VERSION HISTORY: WRITTEN - MD - JUNE 2005
!             ADAPTED FOR MPI - MD - AUG 2005

IMPLICIT NONE
INCLUDE 'mpif.h'

INTEGER,INTENT(IN)::N_R,N_TY,N_Y,N_YS,N_RL
INTEGER,INTENT(IN)::MY_COMM,ERR,RANK
REAL,INTENT(IN)::Y(N_Y,N_TY,N_RL)
REAL,INTENT(INOUT)::YSTATS(N_Y,N_TY,N_YS)

real,dimension(:,:),allocatable :: localsum,globalsum,localmax,&
  localmin,localsse,globalsse 
real :: mymax(n_y), mymin(n_y)

INTEGER T,I,K,allocate_status,tmax(n_y),tmin(n_y),imax(1),imin(1)

allocate(  LOCALSUM(N_Y,N_TY),GLOBALSUM(N_Y,N_TY),LOCALMAX(N_Y,N_TY),&
  LOCALMIN(N_Y,N_TY),LOCALSSE(N_Y,N_TY),GLOBALSSE(N_Y,N_TY),&
  stat=allocate_status )

! FIRST COMPUTE THE LOCAL SUM, MAX AND MIN
DO T=1,N_TY
  DO I=1,N_Y
    LOCALSUM(I,T)=SUM(Y(I,T,:))
    LOCALMAX(I,T)=MAXVAL(Y(I,T,:))
    LOCALMIN(I,T)=MINVAL(Y(I,T,:))
  END DO
END DO
mymax=maxval(localmax,dim=2)
tmax=maxloc(localmax,dim=2)
imax=maxloc(mymax)
mymin=minval(localmin,dim=2)
tmin=minloc(localmin,dim=2)
imin=minloc(mymin)
!print *, 'after local sums, rank=',rank,'maximum=',maxval(mymax),&
!  'minimum=',minval(mymin),'max time=',tmax(imax(1)),'min time=',tmin(imin(1)),&
!  'time series length =',n_ty,'max state=',imax(1),'min state=',imin(1)

CALL MPI_ALLREDUCE(LOCALSUM,GLOBALSUM,N_TY*N_Y,MPI_REAL,MPI_SUM,MY_COMM,ERR)
!print *, 'after allreduce cmd, rank=',rank
YSTATS(:,:,1)=GLOBALSUM/N_R

!print *, 'before max and min calcs,rank=',rank
CALL MPI_REDUCE(LOCALMAX,YSTATS(:,:,2),N_TY*N_Y,MPI_REAL,MPI_MAX,0,MY_COMM,ERR)
CALL MPI_REDUCE(LOCALMIN,YSTATS(:,:,3),N_TY*N_Y,MPI_REAL,MPI_MIN,0,MY_COMM,ERR)

!print *, 'before variance calc,rank=',rank
! COMPUTE LOCAL VARIANCE, THEN GLOBAL VARIANCE USING REDUCE
IF(N_R>1)THEN
  LOCALSSE=0.0
  DO T=1,N_TY
    DO I=1,N_Y
      DO K=1,N_RL
        LOCALSSE(I,T)=LOCALSSE(I,T)+(Y(I,T,K)-YSTATS(I,T,1))**2
      END DO
    END DO
  END DO

  CALL MPI_REDUCE(LOCALSSE,GLOBALSSE,N_TY*N_Y,MPI_REAL,MPI_SUM,0,MY_COMM,ERR)
  YSTATS(:,:,4)=GLOBALSSE/(N_R-1)
ELSE
  YSTATS(:,:,4)=9999
END IF

deallocate(LOCALSUM,GLOBALSUM,LOCALMAX,LOCALMIN,LOCALSSE,GLOBALSSE)

END SUBROUTINE COMPILEYSTATS
