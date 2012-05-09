!**********************************************************************
module string_util
!**********************************************************************
!
!  This module contains the generic subroutines and functions for
!  manipulating strings
!
implicit none

save
private

public :: string_concat, &
     numtostr, &
     eat_whitespace, &
     uppercase, &
     split_string, &
     count_string_occur

interface string_concat
  module procedure strcat_aa, strcat_ai, strcat_ar
end interface

! Explicit interface for overloaded function to convert
! reals and integer to strings
interface numtostr
  module procedure numtostr_r, numtostr_i
end interface

character (*), parameter :: mod_name = 'string_util'

character(*), parameter :: int_fmt='(i0)'
character(*), parameter :: real_fmt='(f9.4)'

contains

!**********************************************************************
subroutine strcat_aa(str1, str2)
!**********************************************************************
use types, only : rprec
implicit none

character(*), intent(INOUT) :: str1
character(*), intent(IN) :: str2

str1 = trim(adjustl(str1)) // str2

return
end subroutine strcat_aa

!**********************************************************************
subroutine strcat_ar(str1, r1)
!**********************************************************************
use types, only : rprec
implicit none

character(*), intent(INOUT) :: str1
real(rprec), intent(IN) :: r1
character(120) :: str2

write (str2,real_fmt) r1

call string_concat(str1,trim(adjustl(str2)))

return
end subroutine strcat_ar

!**********************************************************************
subroutine strcat_ai(str1, i1)
!**********************************************************************
use types, only : rprec
implicit none

character(*), intent(INOUT) :: str1
integer, intent(IN) :: i1
character(120) :: str2

write (str2,int_fmt) i1

call string_concat(str1,trim(adjustl(str2)))

return
end subroutine strcat_ai

!**********************************************************************
function numtostr_r( a, n ) result(c)
!**********************************************************************
!
! This function converts the real variable a to a string b
!
! Inputs
! a : real, scalar value to convert
! n : length of string to return 
!
use types, only : rprec
implicit none

real(rprec), intent(in) :: a
integer, intent(in) :: n
integer, parameter :: l=2*kind(a)+2
character(l) :: b
character(n) :: c
character(25) :: fmt

write(*,*) 'a : ', a
write( fmt, '("(f",i0,".",i0,")")' ) l, l/2
write(*,*) 'fmt : ', fmt

write(b,fmt) a
write(*,*) 'b : ', b

b = trim(adjustl(b))
c = b(:n)
write(*,*) 'c : ', c
return
end function numtostr_r

!**********************************************************************
function numtostr_i( a, n ) result(c)
!**********************************************************************
!
! This function converts the real variable a to a string b
!
! Inputs
! a : real, scalar value to convert
! n : length of string to return 
!
implicit none

integer, intent(in) :: a
integer, intent(in) :: n
integer, parameter :: l = 2*kind(a)+2
character(l) :: b
character(n) :: c

write(b,'(i0)') a
b = adjustl(b)
c = b(:n)

return
end function numtostr_i

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine eat_whitespace (buff, whtspc)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! eats leading and intermediate whitespace, fill trailing space with
! blanks
!

implicit none

character (*), intent (inout) :: buff
character (*), intent (in), optional :: whtspc  !--override default

character (*), parameter :: whtspc_default = achar (9) // achar (32)
                            !--add more characters here if needed
character (1), parameter :: fill_char = ' '

character (1) :: tmp (len (buff))
character (1) :: fill (len (buff))

fill = fill_char
tmp = transfer (buff, tmp)

if (present (whtspc)) then
  tmp = pack (tmp, scan (tmp, whtspc) == 0, fill)
else
  tmp = pack (tmp, scan (tmp, whtspc_default) == 0, fill)
end if

buff = transfer (tmp, buff)

end subroutine eat_whitespace

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function uppercase(str) result(ucstr)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! convert specified string to upper case
!
character (len=*):: str
character (len=len_trim(str)):: ucstr
integer :: i, ilen, iav, ioffset, iqc, iquote

ilen=len_trim(str)
ioffset=iachar('A')-iachar('a')
iquote=0
ucstr=str
do i=1,ilen
  iav=iachar(str(i:i))
  if(iquote==0 .and. (iav==34 .or.iav==39)) then
    iquote=1
    iqc=iav
    cycle
  end if
  if(iquote==1 .and. iav==iqc) then
    iquote=0
    cycle
  end if
  if (iquote==1) cycle
  if(iav >= iachar('a') .and. iav <= iachar('z')) then
    ucstr(i:i)=achar(iav+ioffset)
  else
    ucstr(i:i)=str(i:i)
  end if
end do
return

end function uppercase

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine split_string( string, delim, nseg, sarray )
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! This subroutine splits 'string' based on the specified delimiter
! 'delim'. The number of segments 'nseg' is determined from the string
! based on the delimiter. The output string vector 'sarray' is allocated
! in this subroutine. 
!
use param, only : CHAR_BUFF_LENGTH
use messages

implicit none

character (*), parameter :: sub_name = mod_name // '.split_string'
 
character(*), intent(in) :: string, delim
integer, intent(out) :: nseg

character(CHAR_BUFF_LENGTH), allocatable, dimension(:), intent(inout) :: sarray

! String buffers (assuming length)
character(1024) :: buff_old, buff
integer :: pos, istop, n

! Length of the delimiter
integer :: delim_len

! Get the length of the delimiter (excluding trailing whitespace)
delim_len = len_trim(delim)

! First make sure string is not empty
if( len_trim(string) == 0 ) call error( sub_name, 'specified string is empty')

! Get the number of segments based on the delimiter count
nseg = count_string_occur (string, delim) + 1

! Now allocate string vector
allocate(sarray(nseg))

! Initialize position of delimiter (rewind a bit)
pos=-delim_len+1
! Initialize stop flag
istop=0
! Initialize string buffer
buff = string

n = 0
do 
   n = n + 1
   
   ! If there are more segments than what is expected stop searching
   if( n > nseg ) exit

   ! Save old buffer
   buff_old = buff
   ! New buffer starts after first delimiter position of old buffer
   buff = trim( adjustl( buff_old(pos+delim_len:) ) )
   ! Get the position of the first delmiter in the new buffer
   pos = index( buff, delim )
   
   if( pos > 0 ) then
      sarray(n) = buff(:pos-1)
   else
      ! Assuming this is the last segment
      sarray(n) = trim( adjustl( buff ) )
      exit
   endif

end do 

if( n < nseg ) then
   call error( sub_name, 'number of found segments less than specified number')
elseif( n > nseg ) then
   call error( sub_name, 'number of found segments greater than specified number')
endif 

return
end subroutine split_string

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function count_string_occur (string, substring, overlap) result( countocc )
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! This function counts the number of occurences of 'substring' in
! 'string' --does not count occurence to overlap, e.g. eee is 1 occurence
! of ee, *unless* optional argument overlap = 'yes', then its 2

implicit none

integer :: countocc

character (*), intent (in) :: string, substring
logical, intent (in), optional :: overlap

integer :: i, p, m

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

end function count_string_occur

end module string_util

