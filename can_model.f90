subroutine can_model(ctrl,freq,theta,vegin,tb_can,t_can)

! this canopy model is described in Reports/Models/vegetation and atmospheric/
! combined.sxw.  first the individual leaf dielectric constant is computed,
! then the canopy transmissivity is computed.  the emissivity is obtained 
! from 1 - transmissivity, and is then used to compute canopy brightness
! temperature.  taken from HUT RTM.
!
!   vegin(1)=vegetation water salinity: [ppt]
!   vegin(2)=canopy temperature [K]
!   vegin(3)=vegetation gravimetric water content [frac]
!   vegin(4)=needle or leaf thickness [m]
!   vegin(5)=vegetation biomass [kg/m2]
!
! by mike 9/22/05
!  modified by mike, 06/01/2006, for new inputs

implicit none

integer,intent(in) :: ctrl(9)
real,intent(in) :: freq(ctrl(9)),theta(ctrl(9)),vegin(ctrl(7))
real,intent(inout):: tb_can(2,ctrl(9)),t_can(2,ctrl(9))
integer :: i
complex :: eps_veg

do i=1,ctrl(9)
  !compute leaf dielectric constant
  call epsleaf(freq(i),vegin(1),vegin(2),vegin(3),eps_veg)

  !compute canopy transmissivity
  call can_tran(freq(i),theta(i),eps_veg,vegin(4),vegin(5),t_can(1,i),&
    t_can(2,i))

  !temporary: set all transmissivities to 0.6 a la langlois et al. 2011
  !note this is only applicable for 37 GHz, which is what rhaesung is
  !simulating right now. november 2014. by mtd
  t_can(1,i)=0.6
  t_can(2,i)=0.6

  !compute canopy brightness temperature from emissivity (1-transmissivity),
  !  and the vegetation physical temperature
  tb_can(1,i)=(1-t_can(1,i))*vegin(2)
  tb_can(2,i)=(1-t_can(2,i))*vegin(2)
end do

end subroutine can_model
! -------------------------------------------------------------------------
!
SUBROUTINE EPSLEAF(F,S,T,MG,EVEG)
!
! -------------------------------------------------------------------------
!
!   CODE ORIGINALLY OBTAINED IN MATLAB FROM PULLIAINEN
!
!   CALCULATES THE DIELECTRIC CONSTANT OF FRESH LEAVES
!
!   EVEG = EPSLEAF(F,S,T,MG)
!      EICE:  DIELECTRIC CONSTANT OF FRESH LEAVES
!      F: FREQUENCY IN GHZ
!      S: VEGETATION WATER SALINITY IN PROMILLES (PARTS PER THOUS.) 
!      T: VEGETATION TEMPERATURE IN K
!      MG: GRAVIMETRIC WATER CONTENT (FRACTION)      
!
!   VERSION HISTORY:
!      1.0     ? ?.?.?
!      2.0    MD 15 APR 05 TRANSLATED TO FORTRAN FROM MATLAB
!   
!   USES: EPSWS
!
!   THIS CODE IS DERIVED FROM MATZLER (94), "MICROWAVE (1-100 GHZ) 
!     DIELECTRIC MODEL OF LEAVES." IN IEEE.
!

REAL,INTENT(IN) :: F,S,T,MG
COMPLEX,INTENT(OUT) :: EVEG
REAL TC,FHZ,MD
COMPLEX ESW

!   UNIT CONVERSION
TC=T-273.15
FHZ=F*10**9

MD=1.-MG

CALL EPSWS(FHZ,S,TC,ESW)

EVEG=0.522*(1-1.32*MD)*ESW+0.51+3.84*MD

END SUBROUTINE EPSLEAF

! -------------------------------------------------------------------------
!
SUBROUTINE EPSWS(F,S,T,ESW)
!
! -------------------------------------------------------------------------
!
!   CODE ORIGINALLY OBTAINED IN MATLAB FROM PULLIAINEN
!
!   CALCULATES THE DIELECTRIC CONSTANT OF SALT WATER
!
!   ESW = EPSLEAF(F,S,T)
!      EICE:  DIELECTRIC CONSTANT OF FRESH LEAVES
!      F: FREQUENCY IN HZ
!      S: VEGETATION WATER SALINITY IN PROMILLES (PARTS PER THOUS.) 
!      T: VEGETATION TEMPERATURE IN C
!
!   VERSION HISTORY:
!      1.0     ? ?.?.?
!      2.0    MD 15 APR 05 TRANSLATED TO FORTRAN FROM MATLAB
!   
!   USES: NONE
!
!   THIS CODE IS DERIVED FROM ULABLY ET AL. (1986) MICROWAVE REMOTE SENSING, 
!     ACTIVE AND PASSIVE, VOL III... A SIMILAR MODEL CAN BE FOUND IN KLEIN
!     AND SWIFT (1977) IN IEEE. 

REAL,INTENT(IN) :: F,S,T
COMPLEX,INTENT(OUT) :: ESW
REAL PI,E0,EW_INF,N,A,EW0_T,EW0,B,T0,TW,D,ALFA,SIGMA_25,SIGMA

! CONSTANTS
PI=3.14159
E0= 8.854E-12  ! SEE KLEIN AND SWIFT, P.106
EW_INF=4.9     ! SEE KLEIN AND SWIFT, P.100

! COMPARE NEXT LINE WITH (19) IN KLEIN AND SWIFT:
N = S*(1.707E-2+1.205E-5*S+4.058E-9*S*S)
! COMPARE NEXT LINE WITH (15) IN KLEIN AND SWIFT:
A = 1.00 - 0.2551*N + 5.151E-2 * N*N - 6.889E-3 * N*N*N
! COMPARE NEXT LINE WITH (14) IN KLEIN AND SWIFT:
EW0_T = 87.74 - 0.40008*T + 9.398E-4 * T*T + 1.410E-6 * T*T*T
! COMPARE NEXT LINE WITH (13) IN KLEIN AND SWIFT:
EW0 = EW0_T * A

! NOTE: HERE, I'M INTEGRATING THE OLD SUBROUTINE 'TAU.M' INTO THE NEXT
!   THREE LINES AND COMMENTS OF THIS CODE
! COMPARE NEXT LINE WITH (18) IN KLEIN AND SWIFT:
B=0.1463E-2*N*T+1.00-0.04896*N-0.02967*N*N+5.6441E-3*N*N*N
! COMPARE NEXT LINE WITH (17) IN KLEIN AND SWIFT:
T0=1/(2*PI)*(1.1109E-10-3.824E-12*T+6.938E-14*T*T-5.096E-16*T*T*T)
! COMPARE NEXT LINE WITH (16) IN KLEIN AND SWIFT:
TW = T0* B

! NOTE: HERE, I'M INTEGRATING THE OLD SUBROUTINE 'SIGMA.M' INTO THE NEXT 
!   4 LINES OF THIS CODE
D=25-T
ALFA=2.033E-2+1.266E-4*D+2.464E-6*D*D-S*(1.849E-5-2.551E-7*D+2.551E-8*D*D)
SIGMA_25=S*(0.182521-1.46192E-3*S+2.09324E-5*S*S-1.28205E-7*S*S*S)
SIGMA=SIGMA_25*EXP(-D*ALFA)

! COMPARE NEXT THREE LINES WITH (5) IN KLEIN AND SWIFT:
EW_R=EW_INF+(EW0-EW_INF)/(1+(2*PI*F*TW)**2)
EW_I=(EW0-EW_INF)*2*PI*F*TW/(1+(2*PI*F*TW)**2)+SIGMA/(2*PI*E0*F)

ESW=CMPLX(EW_R,(-1*EW_I))       

END SUBROUTINE EPSWS

! -------------------------------------------------------------------------
!
SUBROUTINE CAN_TRAN(FREQ,THETAD,EPS_VEG,D,B,TCAN_H,TCAN_V)
!
! -------------------------------------------------------------------------
!
!
!   CODE ORIGINALLY OBTAINED IN MATLAB FROM PULLIAINEN
!
!   CALCULATES THE CANOPY TRANSMISSIVITY
!
!   [TCAN_H,TCAN_V] = EPSLEAF(FREQ,THETAD,EPS_VEG,D,N,L,H)
!      TCAN_H,TCAN_V:  HORIZONTAL AND VERTICAL CANOPY TRANSMISSIVITY
!      FREQ: FREQUENCY IN GHZ
!      THETAD: OBSERVATION ANGLE IN DEGREES
!      EPS_VEG: VEGETATION DIELECTRIC CONSTANT (COMPLEX)
!      D: NEEDLE DIAMETER [M]
!      B: VEGETATION BIOMASS [KG/M2]
!
!   VERSION HISTORY:
!      1.0     ? ?.?.?
!      2.0    MD 15 APR 05 TRANSLATED TO FORTRAN FROM MATLAB
!      3.0    MD 01 JUN 06 MODIFIED TO TAKE BIOMASS AND DIAM. INPUTS ONLY
!   
!   USES: NONE
!
!   THIS CODE IS DERIVED FROM WEGMULLER ET AL. (95), SEE REF IN MAIN COMMENTS
!

REAL,INTENT(IN) :: FREQ,THETAD,D,B
COMPLEX,INTENT(IN) :: EPS_VEG
REAL,INTENT(OUT) :: TCAN_H,TCAN_V
REAL PI,MJU0,EPS0,FHZ,THETA,C,K0,KZ0,R_H,R_V,T_H,T_V,A_H,A_V,&
       TAU_H,TAU_V,RHO_VEG
COMPLEX RH,RV,KZ1

! CONSTANTS
PI=3.14159
MJU0=PI*4E-7
EPS0=8.8542E-12
RHO_VEG=950 !KG/M3
      
! UNIT CONVERSIONS
FHZ=FREQ*10**9
THETA=THETAD*PI/180

! WAVE SPEED
C=(MJU0*EPS0)**(-0.5)

K0 = 2*PI*FHZ*SQRT(MJU0*EPS0)
KZ0 = 2*PI*FHZ/C*COS(THETA)
KZ1 = 2*PI*FHZ/C*SQRT(EPS_VEG-SIN(THETA)**2)

RH = (KZ0-KZ1)/(KZ0+KZ1) 
RV = (EPS_VEG*KZ0-KZ1)/(EPS_VEG*KZ0+KZ1)

EPS_VDD=ABS(AIMAG(EPS_VEG))

! HORIZONTAL AND VERTICAL REFLECTIVITY OF A SINGLE LEAF (16)
R_H=ABS(RH*(1-EXP((0.0,1.0)*(-2)*KZ1*D))/(1-RH**2* &
      EXP((0.0,1.0)*(-2)*KZ1*D)))**2
R_V=ABS(RV*(1-EXP((0.0,1.0)*(-2)*KZ1*D))/(1-RV**2* &
      EXP((0.0,1.0)*(-2)*KZ1*D)))**2

! HORIZONTAL (17) AND VERTICAL (18) TRANSMISSIVITY OF A SINGLE LEAF 
T_H = ABS(4*KZ0*KZ1*EXP((0.0,1.0)*(KZ0-KZ1)*D)/((KZ0+KZ1)**2* &
       (1-RH**2*EXP((0.0,1.0)*(-2)*KZ1*D))))**2 

T_V = ABS(4*EPS_VEG*KZ0*KZ1*EXP((0.0,1.0)*(KZ0-KZ1)*D)/ &
       ((EPS_VEG*KZ0+KZ1)**2*(1-RV**2*EXP((0.0,1.0)*(-2)*KZ1*D))))**2

! ABSORPTIVITY OF A SINGLE LEAF (19)
A_H = 1-T_H-R_H
A_V = 1-T_V-R_V

! CANOPY OPACITY
TAU_H = B/RHO_VEG/D*A_H
TAU_V = B/RHO_VEG/D*A_V

! CANOPY TRANSMISSIVITY (25)
TCAN_H=EXP(-TAU_H)
TCAN_V=EXP(-TAU_V)

END SUBROUTINE CAN_TRAN
