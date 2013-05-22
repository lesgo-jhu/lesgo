!!
!!  Copyright 2009,2011,2012 Johns Hopkins University
!!
!!  Licensed under the Apache License, Version 2.0 (the "License"); you may not 
!!  use this file except in compliance with the License. You may obtain a copy of
!!  the License at:
!!
!!    http://www.apache.org/licenses/LICENSE-2.0
!!
!!  Unless required by applicable law or agreed to in writing, software 
!!  distributed under the License is distributed on an "AS IS" BASIS, WITHOUT 
!!  WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the 
!!  License for the specific language governing permissions and limitations under
!!  the License.
!!

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
     string_splice, &
     numtostr, &
     eat_whitespace, &
     uppercase, &
     split_string, &
     count_string_occur

! Concatenates the string. By default eats all trailing whitespace
interface string_concat
  module procedure string_concat_a, string_concat_i, string_concat_r, &
                   string_concat_ai, string_concat_ar, &
                   string_concat_aia, string_concat_ara, &
                   string_concat_aiai, string_concat_arar, &
                   string_concat_aiaia, string_concat_arara, &
                   string_concat_araia, &
                   string_concat_aiaiai, string_concat_ararar, &
                   string_concat_aiaiaia, string_concat_ararara
end interface

interface string_splice
   module procedure string_splice_aa, string_splice_ai, string_splice_ar, &
        string_splice_aia, string_splice_ara, &
        string_splice_aiai, string_splice_arar, &
        string_splice_aiaia, string_splice_arara, &
        string_splice_araia, &
        string_splice_aiaiai, string_splice_ararar, &
        string_splice_aiaiaia, string_splice_ararara
end interface string_splice

! Explicit interface for overloaded function to convert
! reals and integer to strings
interface numtostr
  module procedure numtostr_r, numtostr_i
end interface

character (*), parameter :: mod_name = 'string_util'

character(*), parameter :: iformat='(i0)'
character(*), parameter :: rformat='(f18.5)'

integer, parameter :: BUFF_LENGTH = 64

contains

!///////////////////////////////////////////////////////////////////////////////
!/// STRING_CONCAT
!///////////////////////////////////////////////////////////////////////////////

!**********************************************************************
subroutine string_concat_a(str, str1)
!**********************************************************************
use types, only : rprec
implicit none

character(*), intent(INOUT) :: str
character(*), intent(IN) :: str1

str = trim(adjustl(str)) // trim(adjustl(str1))

return
end subroutine string_concat_a

!**********************************************************************
subroutine string_concat_r(str, r)
!**********************************************************************
use types, only : rprec
implicit none

character(*), intent(INOUT) :: str
real(rprec), intent(IN) :: r
character(BUFF_LENGTH) :: buff

write(buff,rformat) r
call string_concat( str, buff )

return
end subroutine string_concat_r

!**********************************************************************
subroutine string_concat_i(str, i)
!**********************************************************************
use types, only : rprec
implicit none

character(*), intent(INOUT) :: str
integer, intent(IN) :: i
character(BUFF_LENGTH) :: buff

write(buff,iformat) i
call string_concat( str, buff )

return
end subroutine string_concat_i

!**********************************************************************
subroutine string_concat_ai(str, str1, i1)
!**********************************************************************
use types, only : rprec
implicit none

character(*), intent(INOUT) :: str
character(*), intent(IN) :: str1
integer, intent(IN) :: i1

call string_concat(str,str1)
call string_concat(str,i1)

return
end subroutine string_concat_ai

!**********************************************************************
subroutine string_concat_ar(str, str1, r1)
!**********************************************************************
use types, only : rprec
implicit none

character(*), intent(INOUT) :: str
character(*), intent(IN) :: str1
real(rprec), intent(IN) :: r1

call string_concat(str,str1)
call string_concat(str,r1)

return
end subroutine string_concat_ar

!**********************************************************************
subroutine string_concat_aia(str, str1, i1, str2)
!**********************************************************************
use types, only : rprec
implicit none

character(*), intent(INOUT) :: str
character(*), intent(IN) :: str1, str2
integer, intent(IN) :: i1

call string_concat(str,str1)
call string_concat(str,i1)
call string_concat(str,str2)

return
end subroutine string_concat_aia

!**********************************************************************
subroutine string_concat_ara(str, str1, r1, str2)
!**********************************************************************
use types, only : rprec
implicit none

character(*), intent(INOUT) :: str
character(*), intent(IN) :: str1, str2
real(rprec), intent(IN) :: r1

call string_concat(str,str1)
call string_concat(str,r1)
call string_concat(str,str2)

return
end subroutine string_concat_ara

!**********************************************************************
subroutine string_concat_aiaia(str, str1, i1, str2, i2, str3)
!**********************************************************************
use types, only : rprec
implicit none

character(*), intent(INOUT) :: str
character(*), intent(IN) :: str1, str2, str3
integer, intent(IN) :: i1, i2

call string_concat(str,str1)
call string_concat(str,i1)
call string_concat(str,str2)
call string_concat(str,i2)
call string_concat(str,str3)

return
end subroutine string_concat_aiaia

!**********************************************************************
subroutine string_concat_arara(str, str1, r1, str2, r2, str3)
!**********************************************************************
use types, only : rprec
implicit none

character(*), intent(INOUT) :: str
character(*), intent(IN) :: str1, str2, str3
real(rprec), intent(IN) :: r1, r2

call string_concat(str,str1)
call string_concat(str,r1)
call string_concat(str,str2)
call string_concat(str,r2)
call string_concat(str,str3)

return
end subroutine string_concat_arara

!**********************************************************************
subroutine string_concat_aiai(str, str1, i1, str2, i2)
!**********************************************************************
use types, only : rprec
implicit none

character(*), intent(INOUT) :: str
character(*), intent(IN) :: str1, str2
integer, intent(IN) :: i1, i2

call string_concat(str,str1)
call string_concat(str,i1)
call string_concat(str,str2)
call string_concat(str,i2)

return
end subroutine string_concat_aiai

!**********************************************************************
subroutine string_concat_arar(str, str1, r1, str2, r2 )
!**********************************************************************
use types, only : rprec
implicit none

character(*), intent(INOUT) :: str
character(*), intent(IN) :: str1, str2
real(rprec), intent(IN) :: r1, r2

call string_concat(str,str1)
call string_concat(str,r1)
call string_concat(str,str2)
call string_concat(str,r2)

return
end subroutine string_concat_arar

!**********************************************************************
subroutine string_concat_araia(str, str1, r1, str2, i1, str3)
!**********************************************************************
use types, only : rprec
implicit none

character(*), intent(INOUT) :: str
character(*), intent(IN) :: str1, str2, str3
real(rprec) :: r1
integer, intent(IN) :: i1

call string_concat(str,str1)
call string_concat(str,r1)
call string_concat(str,str2)
call string_concat(str,i1)
call string_concat(str,str3)

return
end subroutine string_concat_araia

!**********************************************************************
subroutine string_concat_aiaiai(str, str1, i1, str2, i2, str3, i3)
!**********************************************************************
use types, only : rprec
implicit none

character(*), intent(INOUT) :: str
character(*), intent(IN) :: str1, str2, str3
integer, intent(IN) :: i1, i2, i3

call string_concat(str,str1)
call string_concat(str,i1)
call string_concat(str,str2)
call string_concat(str,i2)
call string_concat(str,str3)
call string_concat(str,i3)

return
end subroutine string_concat_aiaiai

!**********************************************************************
subroutine string_concat_ararar(str, str1, r1, str2, r2, str3, r3 )
!**********************************************************************
use types, only : rprec
implicit none

character(*), intent(INOUT) :: str
character(*), intent(IN) :: str1, str2, str3
real(rprec), intent(IN) :: r1, r2, r3

call string_concat(str,str1)
call string_concat(str,r1)
call string_concat(str,str2)
call string_concat(str,r2)
call string_concat(str,str3)
call string_concat(str,r3)

return
end subroutine string_concat_ararar

!**********************************************************************
subroutine string_concat_aiaiaia(str, str1, i1, str2, i2, str3, i3, str4)
!**********************************************************************
use types, only : rprec
implicit none

character(*), intent(INOUT) :: str
character(*), intent(IN) :: str1, str2, str3, str4
integer, intent(IN) :: i1, i2, i3

call string_concat(str,str1)
call string_concat(str,i1)
call string_concat(str,str2)
call string_concat(str,i2)
call string_concat(str,str3)
call string_concat(str,i3)
call string_concat(str,str4)

return
end subroutine string_concat_aiaiaia

!**********************************************************************
subroutine string_concat_ararara(str, str1, r1, str2, r2, str3, r3, str4)
!**********************************************************************
use types, only : rprec
implicit none

character(*), intent(INOUT) :: str
character(*), intent(IN) :: str1, str2, str3, str4
real(rprec), intent(IN) :: r1, r2, r3

call string_concat(str,str1)
call string_concat(str,r1)
call string_concat(str,str2)
call string_concat(str,r2)
call string_concat(str,str3)
call string_concat(str,r3)
call string_concat(str,str4)

return
end subroutine string_concat_ararara

!///////////////////////////////////////////////////////////////////////////////
!/// STRING_SPLICE
!///////////////////////////////////////////////////////////////////////////////

!*******************************************************************************
subroutine string_splice_aa( s, s1, s2 )
!*******************************************************************************
implicit none

character(*), intent(inout) :: s
character(*), intent(in) :: s1, s2

s = s1 // s2

return
end subroutine string_splice_aa

!**********************************************************************
subroutine string_splice_ar(s, s1, r1)
!**********************************************************************
use types, only : rprec
implicit none

character(*), intent(INOUT) :: s
character(*), intent(in) :: s1
real(rprec), intent(IN) :: r1
character(BUFF_LENGTH) :: b1

write(b1,rformat) r1
!call string_splice(s, s1, trim(adjustl(buff)) )
s = s1 // trim(adjustl(b1))

return
end subroutine string_splice_ar

!**********************************************************************
subroutine string_splice_ai(s, s1, i1)
!**********************************************************************
implicit none

character(*), intent(inout) :: s
character(*), intent(in) :: s1
integer, intent(in) :: i1
character(BUFF_LENGTH) :: b1

write(b1,iformat) i1
!call string_splice(s, s1, trim(adjustl(adjustl(b)))
s = s1 // trim(adjustl(b1))

return
end subroutine string_splice_ai

!**********************************************************************
subroutine string_splice_aia(s, s1, i1, s2 )
!**********************************************************************
implicit none

character(*), intent(inout) :: s
character(*), intent(in) :: s1, s2
integer, intent(in) :: i1
character(BUFF_LENGTH) :: b1

write(b1,iformat) i1

s = s1 // trim(adjustl(b1)) // s2

return
end subroutine string_splice_aia

!**********************************************************************
subroutine string_splice_ara(s, s1, r1, s2 )
!**********************************************************************
use types, only : rprec
implicit none

character(*), intent(inout) :: s
character(*), intent(in) :: s1, s2
real(rprec), intent(in) :: r1
character(BUFF_LENGTH) :: b1

write(b1,rformat) r1

s = s1 // trim(adjustl(b1)) // s2

return
end subroutine string_splice_ara

!**********************************************************************
subroutine string_splice_aiai(s, s1, i1, s2, i2 )
!**********************************************************************
implicit none

character(*), intent(inout) :: s
character(*), intent(in) :: s1, s2
integer, intent(in) :: i1, i2
character(BUFF_LENGTH) :: b1, b2

write(b1,iformat) i1
write(b2,iformat) i2

s = s1 // trim(adjustl(b1)) // s2 // trim(adjustl(b2))

return
end subroutine string_splice_aiai

!**********************************************************************
subroutine string_splice_arar(s, s1, r1, s2, r2 )
!**********************************************************************
use types, only : rprec
implicit none

character(*), intent(inout) :: s
character(*), intent(in) :: s1, s2
real(rprec), intent(in) :: r1, r2
character(BUFF_LENGTH) :: b1, b2

write(b1,rformat) r1
write(b2,rformat) r2

s = s1 // trim(adjustl(b1)) // s2 // trim(adjustl(b2))

return
end subroutine string_splice_arar

!**********************************************************************
subroutine string_splice_aiaia(s, s1, i1, s2, i2, s3 )
!**********************************************************************
implicit none

character(*), intent(inout) :: s
character(*), intent(in) :: s1, s2, s3
integer, intent(in) :: i1, i2
character(BUFF_LENGTH) :: b1, b2

write(b1,iformat) i1
write(b2,iformat) i2

s = s1 // trim(adjustl(b1)) // s2 // trim(adjustl(b2)) // s3

return
end subroutine string_splice_aiaia

!**********************************************************************
subroutine string_splice_arara(s, s1, r1, s2, r2, s3 )
!**********************************************************************
use types, only : rprec
implicit none

character(*), intent(inout) :: s
character(*), intent(in) :: s1, s2, s3
real(rprec), intent(in) :: r1, r2
character(BUFF_LENGTH) :: b1, b2

write(b1,rformat) r1
write(b2,rformat) r2

s = s1 // trim(adjustl(b1)) // s2 // trim(adjustl(b2)) // s3

return
end subroutine string_splice_arara

!**********************************************************************
subroutine string_splice_araia(s, s1, r1, s2, i2, s3 )
!**********************************************************************
use types, only : rprec
implicit none

character(*), intent(inout) :: s
character(*), intent(in) :: s1, s2, s3
real(rprec), intent(in) :: r1
integer, intent(in) :: i2
character(BUFF_LENGTH) :: b1, b2

write(b1,rformat) r1
write(b2,iformat) i2

s = s1 // trim(adjustl(b1)) // s2 // trim(adjustl(b2)) // s3

return
end subroutine string_splice_araia


!**********************************************************************
subroutine string_splice_aiaiai(s, s1, i1, s2, i2, s3, i3)
!**********************************************************************
implicit none

character(*), intent(inout) :: s
character(*), intent(in) :: s1, s2, s3
integer, intent(in) :: i1, i2, i3
character(BUFF_LENGTH) :: b1, b2, b3

write(b1,iformat) i1
write(b2,iformat) i2
write(b3,iformat) i3

s = s1 // trim(adjustl(b1)) // s2 // trim(adjustl(b2)) // s3 // trim(adjustl(b3))

return
end subroutine string_splice_aiaiai

!**********************************************************************
subroutine string_splice_ararar(s, s1, r1, s2, r2, s3, r3)
!**********************************************************************
use types, only : rprec
implicit none

character(*), intent(inout) :: s
character(*), intent(in) :: s1, s2, s3
real(rprec), intent(in) :: r1, r2, r3
character(BUFF_LENGTH) :: b1, b2, b3

write(b1,rformat) r1
write(b2,rformat) r2
write(b3,rformat) r3

s = s1 // trim(adjustl(b1)) // s2 // trim(adjustl(b2)) // s3 // trim(adjustl(b3))

return
end subroutine string_splice_ararar

!**********************************************************************
subroutine string_splice_aiaiaia(s, s1, i1, s2, i2, s3, i3, s4)
!**********************************************************************
implicit none

character(*), intent(inout) :: s
character(*), intent(in) :: s1, s2, s3, s4
integer, intent(in) :: i1, i2, i3
character(BUFF_LENGTH) :: b1, b2, b3

write(b1,iformat) i1
write(b2,iformat) i2
write(b3,iformat) i3

s = s1 // trim(adjustl(b1)) // s2 // trim(adjustl(b2)) // s3 // trim(adjustl(b3)) // s4

return
end subroutine string_splice_aiaiaia

!**********************************************************************
subroutine string_splice_ararara(s, s1, r1, s2, r2, s3, r3, s4)
!**********************************************************************
use types, only : rprec
implicit none

character(*), intent(inout) :: s
character(*), intent(in) :: s1, s2, s3, s4
real(rprec), intent(in) :: r1, r2, r3
character(BUFF_LENGTH) :: b1, b2, b3

write(b1,rformat) r1
write(b2,rformat) r2
write(b3,rformat) r3

s = s1 // trim(adjustl(b1)) // s2 // trim(adjustl(b2)) // s3 // trim(adjustl(b3)) // s4

return
end subroutine string_splice_ararara

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

