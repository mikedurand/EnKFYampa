!This subroutine calculates serial dates based on year, month, day and hour
!  based on Matlab's definition obtained by numbering each day.  Jan 1, 0000 
!  is considered day 1, Jan 2, 0000 day 2, etc.  The calculation is made
!  more involved because of the fact that every fourth year is a leap year, 
!  except that every hundredth year is not a leap year, except that every 
!  four hundredth year is a leap year.  Thus, each cycle of one hundred years
!  does not have the same number of days, but each cycle of four hundred years
!  does have the same number of days (146097 days).  The approach followed
!  here is to compute the number of elapsed 400 year cycles, then the number
!  of elapsed 100 year cycles, 4 year cycles and 1 year cycles, computing the
!  proper number of elapsed days for each of the cycles.
!
!See http://scienceworld.wolfram.com/astronomy/LeapYear.html for more info. 
!  Coded by Mike October 25, 2005, based on code written and tested in Matlab. 
!
!Uses: Std2J, LeapYear

Subroutine Serial_Date(Year,Month,Day,Hour,N)

Implicit None

Integer,Intent(In) :: Year,Month,Day,Hour
Double Precision,Intent(Out) :: N
Integer :: Jday,Ndays400,Ndays100,Ndays4,Ndays1,n400c,n100c,n4c,nyr400,nyr100,&
  nyr4,N1,N2,N3,N4,N5,Year00,Year1,NdaysYr00,NdaysYr1,dummy(12)

Call Std2J(Year,Month,Day,Jday)

Ndays400=146097 !Number of days in the 400 year cycle
Ndays100=36524  !Number of days in 100 year cycle if first year is NOT leap yr
Ndays4=1461     !Number of days in 4 year cycle where first year IS leap yr
Ndays1=365      !Number of days in non-leap year

!Compute N1, Number of days in completely elapsed 400 year cycles since 0000
n400c=year/400    !This is the number of elapsed cycles.  Fortran rounds this
                  !  down to the next integer value.
N1=n400c*Ndays400

!Compute N2, Number of days in completely elapsed 100 year cycles since the end
! of the last 400 year cycle.  The key to the below algorithm is that the first
! 100 year cycle ALWAYS has the first year as a leap year, and subsequent 100
! year cycles in the same 400 year period NEVER have the first year as a leap
! year
nyr400=Year-n400c*400 ! This is the number of elapsed years since last 400 year
                      !   cycle ended
n100c=nyr400/100      ! This is the number of elapsed 100 year cycles.  Fortran
                      !   rounds this down to the next integer.

If(n100c.Eq.0) Then
  N2=0
Else
  !Note that the first year in first century of the four hundred year period 
  ! is a leap year, but the first year in the subsequent centuries is NOT a 
  ! leap year
  N2=(Ndays100+1)+Ndays100*(n100c-1)
End If

!Compute N3, the number of days in complete 4 year cycles since the end of the 
!  last 100 year cycle.  The key here is to determine whether or not the first 
!  year of this century is a leap year or not when accounting for all of the
!  days.

!First, determine the first year of the current century and decide whether or
!  or not it is a leap year
nyr100=year-n400c*400-n100c*100 !This is the number of elapsed years since the
                                !  the last 100 year cycle ended
n4c=nyr100/4   !This is the number of elapsed 4 year cycles.  Fortran rounds 
               ! this down to the next integer
Year00=year-nyr100

Call LeapYear(Year00,NdaysYr00,dummy) !NdaysYr00 is the number of days in the
                                      !  first year of the current century
If(n4c.Eq.0) Then
  N3=0
ElseIf(n4c.Gt.0) Then
  !The first term in parentheses below indicates the total number of days
  !  in the first four year cycle since the beginning of the century, while
  !  the second term indicates the number of days in all subsequent four year
  !  cycles. 
  N3=(NdaysYr00+Ndays1*3)+(n4c-1)*Ndays4
End If

!Compute N4, the number of days in complete 1 year cycles since the end of
!  the last 4 year cycle.  The key here is to determine whether or not the 
!  first year of the current 4 year cycle is a leap year or not when 
!  accounting for all of the days.

nyr4=year-n400c*400-n100c*100-n4c*4 ! This is the number of elapsed years since
                                   !  beginning of this 4 year cycle 
Year1=year-nyr4 !This is first year of this cycle
Call LeapYear(Year1,NdaysYr1,dummy) !NdaysYr00 is the number of days in the
                                      !  first year of the current century
If(nyr4.Eq.0) Then
  N4=0
ElseIf(nyr4.Gt.0) Then
  !The first term below indicates the number of days in the first year 
  !  of the four year cycle, while the second term indicates the number 
  !  of days in all subsequent years until the present. 
  N4=NdaysYr1+(nyr4-1)*Ndays1
End IF

!Compute N5, which is simply the current Julian day
N5=Jday

!Compute N as the sum of N1 - N5, plus the hour converted to fraction of a day
N=Dble(N1+N2+N3+N4+N5)+Dble(Hour/24.0)

End Subroutine Serial_Date
