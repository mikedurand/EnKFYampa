!Function to calculate the number of days in a year and number of days in a
!  and number of days in each month as a function of the year

Subroutine LeapYear(Year,Ndays,DaysInMonth)

Implicit None

Integer,Intent(In) :: Year
Integer,Intent(Out):: Ndays,DaysInMonth(12)

DaysInMonth=(/31,-999,31,30,31,30,31,31,30,31,30,31/)

If(Mod(Year,400).Eq.0) Then
  !Year is a Leap Year
  DaysInMonth(2)=29
ElseIf(Mod(Year,100).Eq.0) Then
  !Year is not a Leap Year
  DaysInMonth(2)=28
ElseIf(Mod(Year,4).Eq.0) Then
  !Year is a Leap Year
  DaysInMonth(2)=29
Else
  !Year is not a Leap Year
  DaysInMonth(2)=28
End If

Ndays=Sum(DaysInMonth)

End Subroutine LeapYear
