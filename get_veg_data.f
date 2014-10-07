C This subroutine written by Mike to provide a way to call the vegin.f
C subroutine where all of the vegetation data is stored, and pass that
C data to a F90 (or F77) calling subroutine or program.
C 
C TRAN_VEG_NIR and REF_VEG_NIR are both isolated from the upper-storey,
C   NIR, direct component of the overall Ref0 and Tran0 arrays: index 1,2,1

      SUBROUTINE GET_VEG_DATA(POROSITY,LAI,TRAN_VEG_NIR,REF_VEG_NIR)

      REAL :: POROSITY(13),LAI(13,12,2),TRAN_VEG_NIR(13),REF_VEG_NIR(13)

      INCLUDE 'comsib.in'

      CALL VEGIN

C     Assign porosity and leaf area index
      POROSITY=POROS0
      LAI=ZLT0

C     Assign NIR reflectivity and transmissivity
      TRAN_VEG_NIR=TRAN0(1:13,1,2,1)
      REF_VEG_NIR=REF0(1:13,1,2,1)

      END SUBROUTINE
