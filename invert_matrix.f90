! -------------------------------------------------------------------------
!
SUBROUTINE INVERT_MATRIX(MA,INV,N)
!
! -------------------------------------------------------------------------
!     
!     MATRIX INVERSION ALGORITHM OBTAINED FROM GOTOP FORTRAN 90 TEXT,
!     ISBN:957-566-172-9
!
!     COPIED FROM THE COMPANION CD TO THAT TEXT BY MIKE, 1 APRIL 2005

IMPLICIT NONE

INTEGER,INTENT(IN) :: N
REAL,INTENT(IN) :: MA(N,N)
REAL,INTENT(OUT) :: INV(N,N)
REAL,dimension(:,:),allocatable :: temp
INTEGER I,J
 
allocate(temp(n,n))

DO I=1,N
  DO J=1,N
    TEMP(I,J)=MA(I,J)
    INV(I,J)=0.
  END DO
  INV(I,I)=1.
END DO

!print *, 'before upper call'
CALL UPPER(TEMP,INV,N)
!print *, 'before lower call'
CALL LOWER(TEMP,INV,N)

!print *, 'before inv calc'
DO I=1,N
  DO J=1,N
    INV(I,J)=INV(I,J)/TEMP(I,I)
  END DO
END DO

!print *, 'before deallocation'
deallocate(temp)

CONTAINS

SUBROUTINE UPPER(M,S,N)
INTEGER,INTENT(IN):: N
INTEGER I,J,K
REAL E
REAL,INTENT(INOUT):: M(N,N)
REAL,INTENT(INOUT):: S(N,N)

DO I=1,N-1
  DO J=I+1,N            
    E=M(J,I)/M(I,I)
    DO K=1,N
      M(J,K)=M(J,K)-M(I,K)*E
      S(J,K)=S(J,K)-S(I,K)*E
    END DO
  END DO
END DO

END SUBROUTINE UPPER

SUBROUTINE LOWER(M,S,N)
INTEGER,INTENT(IN):: N
REAL,INTENT(INOUT):: M(N,N)
REAL,INTENT(INOUT):: S(N,N)
INTEGER I,J,K
REAL E

DO I=N,2,-1
  DO J=I-1,1,-1         
    E=M(J,I)/M(I,I)
    DO K=1,N
      M(J,K)=M(J,K)-M(I,K)*E
      S(J,K)=S(J,K)-S(I,K)*E 
    END DO
  END DO
END DO

END SUBROUTINE LOWER

END SUBROUTINE INVERT_MATRIX
