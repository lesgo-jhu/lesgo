module string_util
implicit none

public

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--eats leading and intermediate whitespace, fill trailing space with
!  blanks
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine eat_whtspc (buff, whtspc)
implicit none

character (*), intent (inout) :: buff
character (*), intent (in), optional :: whtspc  !--override default

character (*), parameter :: whtspc_default = achar (9) // achar (32)
                            !--add more characters here if needed
character (1), parameter :: fill_char = ' '

character (1) :: tmp (len (buff))
character (1) :: fill (len (buff))

!---------------------------------------------------------------------

fill = fill_char
tmp = transfer (buff, tmp)

if (present (whtspc)) then
  tmp = pack (tmp, scan (tmp, whtspc) == 0, fill)
else
  tmp = pack (tmp, scan (tmp, whtspc_default) == 0, fill)
end if

buff = transfer (tmp, buff)

end subroutine eat_whtspc

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
