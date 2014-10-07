
subroutine rad_xfer_mod(ctrl,freq,theta,aux_in,snow_in,can_in,atm_in,&
  tb_ssca_out,pixel,replicate,rank,meas,v_c,ierr) 

! this 'main' code radiative transfer model implements an atmospheric model,
! canopy model, and snow/soil model using the scheme in 'new memls inputs.pdf'.
! this code performs 'ncalcs' calculations for one set of physical (i.e. snow,
! canopy, etc.) parameters.  the loop over different frequencies is handled 
! in the individual radiative transfer models.  soil, snow, and canopy 
! calculations are performed at two polarizations, while atmospheric 
! calculations are done without polarization distinction; cosmic brightness 
! temperature is specified as 2.7 K without frequency or polarization
! distinction.  the specific versions of the codes used are described in the 
! headers of each individual program codes.  definitions of the inputs other 
! than the ctrl array are given in the individual program codes as well. 
! the snow boundary condition brightness temperature is computed according to 
! a scheme described in Reports/Models/vegetation and atmospheric/
! new memls inputs.sxw. 
!
! modified by mike, 04 sep 2006, to account for fractional vegetation cover
!
! inputs:
!
! ctrl(1) - ncalcs: number of calculations to be done
! ctrl(2) - n_lyrs: number of snow layers
! ctrl(3) - atm_switch: binary atmosphere switch - 1 = on, 0 = off
! ctrl(4) - can_switch: binary canopy switch - 1 = on, 0 = off
! ctrl(5) - n_aux_ins: length of soil and snow auxiliary array (columns)
! ctrl(6) - n_snow_ins: length of snow input array (columns)
! ctrl(7) - n_can_ins: length of canopy input array (columns)
! ctrl(8) - n_atm_ins: length of atmosphere input array (columns)
! ctrl(9) - n_freq: number of frequencies to compute tbs 
! ...
! v_c     - fractional vegetation cover
!
! outputs:
!
! tb_ssca_out(i,1) - t.o.a. pred. horiz. brightness temp for i=1,n_freq
! tb_ssca_out(i,2) - t.o.a. pred. vert. brightness temp for i=1,n_freq 
!
! calls: atm_model, can_model, ss_model
! called from: interfacez.f90


implicit none

integer, intent(in) :: ctrl(9),ierr
real,intent(in) :: freq(ctrl(9)),theta(ctrl(9)),aux_in(ctrl(5)),&
  snow_in(ctrl(2),ctrl(6)),can_in(ctrl(7)),atm_in(ctrl(8)),v_c
real,intent(out) :: tb_ssca_out(2,ctrl(9))
integer :: ctrlc(9)
real :: tb_cosmic
real,dimension(:),allocatable :: tb_atm,tran_atm
real,dimension(:,:),allocatable :: tb_can,tran_can,tb_snow_bc,tb_snow
integer i,j,k,nlyr
integer,intent(in) :: pixel,replicate,rank,meas 

allocate( tb_atm(ctrl(9)),tran_atm(ctrl(9)),tb_can(2,ctrl(9)),&
  tran_can(2,ctrl(9)),tb_snow_bc(2,ctrl(9)),tb_snow(2,ctrl(9)) )

! constants
tb_cosmic=2.7 !value from hut rtms

! atmospheric calculations
if(ctrl(3).eq.1) then
  call atm_model(ctrl(9),ctrl(8),freq,theta,atm_in,&
    tb_atm,tran_atm)
else
  tb_atm=0.
  tran_atm=1.
end if

! canopy calculations
if(ctrl(4).eq.1) then
  call can_model(ctrl,freq,theta,can_in,tb_can,tran_can)
else
  tb_can=0.
  tran_can=1.
end if

! compute memls boundary condition
do j=1,2
  do k=1,ctrl(9)
    !this old way assumes full vegetation cover
    !tb_snow_bc(j,k)=(tb_cosmic*tran_atm(k)+tb_atm(k))*tran_can(j,k)+&
    !  tb_can(j,k)
    !the new way allows for fractional coverage.  i derived it in the written
    !notebook on june 22, 2006, p.75
    tb_snow_bc(j,k)=(1.+(tran_can(j,k)-1)*v_c)*(tb_cosmic*tran_atm(k)+&
      tb_atm(k))+v_c*tb_can(j,k)
  end do
end do

!snow and soil calculations
call ss_model(ctrl,freq,theta,snow_in, tb_snow_bc,aux_in,tb_snow,pixel,&
  replicate,rank,meas)

!if(meas.eq.85) print *, 'radxfer_mod: tb_snow_bc=',tb_snow_bc,'tb_snow=',&
!  tb_snow,'snow_in=',snow_in

!compute satellite-observed brightness temperature
do j=1,2
  do k=1,ctrl(9)
    !this old way assumes full vegetation cover
    !tb_ssca_out(j,k)=tb_atm(k)+tb_can(j,k)*tran_atm(k)+&
    !  tb_snow(j,k)*tran_can(j,k)*tran_atm(k)
    !the new way allows for fractional vegetation coverage.  i derived it in
    !  the written notebook on june 22, 2006, p. 76. 
    tb_ssca_out(j,k)=((1.+(tran_can(j,k)-1.)*v_c)*tb_snow(j,k)+v_c*&
      tb_can(j,k))*tran_atm(k)+tb_atm(k)
  end do
end do

deallocate(tb_atm,tran_atm,tb_can,tran_can,tb_snow_bc,tb_snow)

end subroutine rad_xfer_mod 
