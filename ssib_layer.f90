subroutine ssib_layer(snowdepth,dz)

! subroutine to determine layer depths based on ssib3 layering scheme
! copied directly from the ssib3.f code
! by mike 29 sep 05

implicit none
real,intent(in) :: snowdepth
real,intent(out):: dz(3)

IF(snowdepth.gt.0.05.and.snowdepth.le.0.06) THEN     
  dz(1)=0.02
  dz(2)=0.02
  dz(3)=snowdepth- dz(1)- dz(2)
ELSE IF ( snowdepth.gt.0.06.and.snowdepth.le.0.08) then
  dz(3)=0.02
  dz(2)=0.02
  dz(1)=snowdepth- dz(3)- dz(2)
ELSE IF ( snowdepth.gt.0.08.and.snowdepth.le.0.62) then
  dz(3)=0.02
  dz(2)=(snowdepth- dz(3))*0.33333333 
  dz(1)=(snowdepth- dz(3))*0.66666667
ELSE IF ( snowdepth.gt.0.62) then
  dz(3)=0.02
  dz(2)=0.20
  dz(1)=snowdepth- dz(3)- dz(2)        
End IF 

end subroutine ssib_layer 
