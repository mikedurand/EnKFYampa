SUBROUTINE READ_VEG(N_YP,VEG_MAP)

IMPLICIT NONE

INTEGER, INTENT(IN) :: N_YP
INTEGER, INTENT(OUT) :: VEG_MAP(N_YP)
INTEGER :: I

OPEN(UNIT=1,FILE='veg_map.in',STATUS='OLD')
DO I=1,N_YP
  READ(1,'(I4)') VEG_MAP(I)
END DO
CLOSE(1)

END SUBROUTINE READ_VEG
