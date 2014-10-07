subroutine interpolate(x,y,xhat,yhat,n)

!should interpolate between any points in (x,y). if values of xhat lie outside
!values of x, then values of yhat are determined by extrpolating either x(1)
!and x(2), or between x(n) and x(n-1). it is assumed that x increases 
!monotonically

!by mike, august 7, 2006
!tested versus a matlab test program in Models/DISORT/interpolation.m

implicit none

!in/out
integer,intent(in) :: n
real,intent(in) :: x(n),y(n),xhat
real,intent(out) :: yhat
!locals
integer :: k,iclose,i
real :: dx(n)

if(xhat.lt.x(1))then
  !extrpolate low
  yhat=y(2)+(y(1)-y(2))*(xhat-x(2))/(x(1)-x(2))
elseif(xhat.gt.x(n))then
  !extrpolate high
  yhat=y(n)+(y(n)-y(n-1))*(xhat-x(n))/(x(n)-x(n-1))
else
  !interpolate
  !1) find closest value of x to xhat
  !1.a) compute distances from x to xhat
  do k=1,n
    dx(k)=abs(x(k)-xhat)
  end do
  !1.b) find smallest distance
  call get_min_loc(dx,iclose,n)
  !2) set i for interpolation based on whether xhat(j) is bigger or 
  !    smaller than x(iclose)
  if(xhat.lt.x(iclose))then
    i=iclose-1
  else
    i=iclose
  end if
  yhat=y(i+1)+(y(i+1)-y(i))*(xhat-x(i+1))/(x(i+1)-x(i))
end if

end subroutine interpolate
