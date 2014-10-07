subroutine get_min_loc(a,min_loc_a,n_a)

implicit none

integer,intent(in)::n_a
real,intent(in)::a(n_a)
integer,intent(out)::min_loc_a

real :: min_a
integer :: i

min_a=a(1)
min_loc_a=1

do i=2,n_a
  if(a(i)<min_a)then
    min_a=a(i)
    min_loc_a=i
  end if
end do

end subroutine get_min_loc
