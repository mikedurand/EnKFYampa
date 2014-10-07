! Subroutine to compute the Julian Day from the month, day of the month and
!  the year

Subroutine Std2J(Year,Month,Mday,Jday)

Implicit None

Integer,Intent(In) :: Month,Mday,Year
Integer,Intent(Out) :: Jday
Integer :: DaysInMonth(12)

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

Jday=Mday+Sum(DaysInMonth(1:Month-1))

End Subroutine Std2J
