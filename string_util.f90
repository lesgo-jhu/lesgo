module string_util
implicit none

public

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--count occurences of substring in string
!--does not count occurence to overlap, e.g. eee is 1 occurence of ee,
!  *unless* optional argument overlap = 'yes', then its 2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function countocc (string, substring, overlap)
implicit none

integer :: countocc

character (*), intent (in) :: string, substring
logical, intent (in), optional :: overlap

integer :: i, p, m

!---------------------------------------------------------------------

m = len (substring) - 1
if (present (overlap)) then
  if (overlap) m = 0
end if

countocc = 0
i = index (string, substring)
p = 1
do while (i /= 0)
  countocc = countocc + 1
  p = p + i + m
  i = index (string(p:), substring)
end do

end function countocc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module string_util
