!**********************************************************************
module strmod
!**********************************************************************
!  This module contains the generic subroutines for concatenating characters,
!  integers and reals for creating file names and other combined characters.
!
use types, only : rprec

save
private

public :: strcat

interface strcat
  module procedure strcat_aa, strcat_ai, strcat_ar
end interface

character(*), parameter :: int_fmt='(i0)'
character(*), parameter :: real_fmt='(f9.4)'

!**********************************************************************
contains
!**********************************************************************

!**********************************************************************
subroutine strcat_aa(str1, str2)
!**********************************************************************
implicit none

character(*), intent(INOUT) :: str1
character(*), intent(IN) :: str2

str1 = trim(adjustl(str1)) // str2

return
end subroutine strcat_aa

!**********************************************************************
subroutine strcat_ai(str1, i1)
!**********************************************************************
implicit none

character(*), intent(INOUT) :: str1
integer, intent(IN) :: i1
character(120) :: str2

write (str2, int_fmt) i1

call strcat(str1,trim(adjustl(str2)))

end subroutine strcat_ai

!**********************************************************************
subroutine strcat_ar(str1, r1)
!**********************************************************************
implicit none

character(*), intent(INOUT) :: str1
real(rprec), intent(IN) :: r1
character(120) :: str2

write (str2,real_fmt) r1

call strcat(str1,trim(adjustl(str2)))

return
end subroutine strcat_ar

!!**********************************************************************
!subroutine strcat_aai(str1, str2, i1)
!!**********************************************************************
!implicit none

!character(*), intent(INOUT) :: str1
!character(*), intent(IN) :: str2
!integer, intent(IN) :: i1

!call strcat(str1,str2)
!call strcat(str1,i1)

!end subroutine strcat_aai

!!**********************************************************************
!subroutine strcat_aar(str1, str2, r1)
!!**********************************************************************
!implicit none

!character(*), intent(INOUT) :: str1
!character(*), intent(IN) :: str2
!real(rprec), intent(IN) :: r1

!call strcat(str1,str2)
!call strcat(str1,r1)

!return
!end subroutine strcat_aar

!!**********************************************************************
!subroutine strcat_aiai(str1, i1, str2, i2)
!!**********************************************************************
!implicit none

!character(*), intent(INOUT) :: str1
!character(*), intent(IN) :: str2
!real(rprec), intent(IN) :: i1,i2

!call strcat(str1,i1)
!call strcat(str1,str2)
!call strcat(str1,i2)

!return
!end subroutine strcat_aiai

!!**********************************************************************
!subroutine strcat_arar(str1, r1, str2, r2)
!!**********************************************************************
!implicit none

!character(*), intent(INOUT) :: str1
!character(*), intent(IN) :: str2
!real(rprec), intent(IN) :: r1,r2

!call strcat(str1,r1)
!call strcat(str1,str2)
!call strcat(str1,r2)

!return
!end subroutine strcat_arar

!!**********************************************************************
!subroutine strcat_ararar(str1, r1, str2, r2, str3, r3)
!!**********************************************************************
!implicit none

!character(*), intent(INOUT) :: str1
!character(*), intent(IN) :: str2,str3
!real(rprec), intent(IN) :: r1,r2,r3

!call strcat(str1,r1)
!call strcat(str1,str2)
!call strcat(str1,r2)
!call strcat(str1,str3)
!call strcat(str1,r3)

!return
!end subroutine strcat_ararar

!!**********************************************************************
!subroutine strcat_aiaiai(str1, i1, str2, i2, str3, i3)
!!**********************************************************************
!implicit none

!character(*), intent(INOUT) :: str1
!character(*), intent(IN) :: str2, str3
!real(rprec), intent(IN) :: i1,i2, i3

!call strcat(str1,i1)
!call strcat(str1,str2)
!call strcat(str1,i2)
!call strcat(str1,str3)
!call strcat(str1,i3)

!return
!end subroutine strcat_aiaiai



end module strmod
