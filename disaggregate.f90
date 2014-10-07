subroutine disaggregate(n_ps,n_pf,n_u,nsteps,nldas_forcing,elevation,zf,zs,u,&
  precip_lapse,temp_lapse,mikesmap,gen_switch,invert_precip)

!need to comment this subroutine
!invert_precip - equal to negative one to invert the lapse rate for the 
!  generator only; set equal to positive one to keep the default lapse rate
!  specified in the precip_lapse variable

implicit none

!inputs
integer,intent(in) :: n_ps, n_pf, n_u, nsteps
real,intent(in) :: nldas_forcing(n_u,nsteps,n_pf),elevation(n_ps),zf(n_pf,2),&
  zs(n_ps,2),precip_lapse,temp_lapse,invert_precip
logical,intent(in) :: gen_switch
!output
real,intent(out) :: u(n_u,nsteps,n_ps)
integer,intent(out) :: mikesmap(625)
!locals
integer :: i,j,k,p,num_i,i_i,min_loc
integer,dimension(:),allocatable :: my_pix
real :: dist(625,9), w(625,9), elev_forc(9),tmp(9),wind(nsteps)

! 1) COMPUTE DISTANCE FROM EACH STATE PIXEL CENTER TO EACH NLDAS PIXEL CENTER
do i=1,n_ps
  do j=1,n_pf
    dist(i,j)=( (zs(i,1)-zf(j,1))**2+(zs(i,2)-zf(j,2))**2 )**0.5
  end do
end do

! 2) ESTIMATE ELEVATION AT EACH NLDAS PIXEL
! Before, we used the closest state elevation as the nldas elevation.
! Now, we find the mean elevation of all the state pixels associated
! with a given nldas pixel.

! 2.1) Figure out which state pixels are closest to each nldas pixel
do i=1,n_ps
  do i_i=1,9
    tmp(i_i)=dist(i,i_i)
  enddo
  call get_min_loc(tmp,min_loc,n_pf)
  mikesmap(i)=min_loc !This maps each state pixel to an NLDAS pixel
end do

! 2.2) Find the average elevation of all state associated w/ each NLDAS pixel
do i=1,n_pf
  !2.2.1) Count the number of state pixels associated w/ NLDAS pixel i
  num_i=0
  do j=1,n_ps
    if(mikesmap(j).eq.i)then
      num_i=num_i+1
    end if
  end do
  !2.2.2) Fill my_pix, which has the indeces of all state pixels associated
  allocate(my_pix(num_i))
  k=0
  do j=1,n_ps
    if(mikesmap(j).eq.i)then
      k=k+1
      my_pix(k)=j
    end if  
  end do
  !2.2.3) Compute mean elevation of all state pixels whose indeces are in my_pix
  elev_forc(i)=0. 
  do j=1,num_i
    elev_forc(i)=elev_forc(i)+elevation(my_pix(j))/num_i
  end do
  deallocate(my_pix)
end do 

! 3) LOCAL PIXEL LOOP: INTERPOLATE FORCING,COMPUTE SNOWFALL DATA

do p=1,n_ps
  !3.1
  ! assign forcing data for pixel p based on nearest neighbor 
  do i_i=1,9
    tmp(i_i)=dist(p,i_i)
  enddo
  call get_min_loc(tmp,min_loc,n_pf)
  do i=1,n_u
    do j=1,nsteps
      u(i,j,p)=nldas_forcing(i,j,min_loc)
    end do
  end do 

  !3.2

  ! perturb precipitation values based on elevation...
  !  the way the elevation gradient is set up gives d(elevation)<0 if the 
  !  state pixel is higher than the forcing pixel and >0 otherwise.  thus
  !  the elevation lapse rate must be >0 here, so that the addition on the
  !  RHS is >0 for state pixels higher than forcing pixels; the precipitation
  !  lapse rate must be <0. 

  do j=1,nsteps
    ! precipitation
    if(u(3,j,p).gt.0.0)then
      if(gen_switch.eq..true.) then
        !if this is the generator, allow lapse rate inversion based on
        !  the invert_precip variable
        u(3,j,p)=u(3,j,p)-invert_precip*precip_lapse*(elev_forc(min_loc)-&
          elevation(p))
      else
        u(3,j,p)=u(3,j,p)-precip_lapse*(elev_forc(min_loc)-elevation(p))
      end if
      if(u(3,j,p).lt.0.0) u(3,j,p)=1.25E-5
    end if

    ! temperature
    if(u(6,j,p).gt.0.0)then
      u(6,j,p)=u(6,j,p)-temp_lapse*(elev_forc(min_loc)-elevation(p))
    end if
  end do

end do

end subroutine disaggregate
