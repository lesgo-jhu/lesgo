!!
!!  Copyright (C) 2009-2013  Johns Hopkins University
!!
!!  This file is part of lesgo.
!!
!!  lesgo is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License as published by
!!  the Free Software Foundation, either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  lesgo is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with lesgo.  If not, see <http://www.gnu.org/licenses/>.
!!

module messages
use types, only : rp => rprec
#ifdef PPMPI
  use mpi
#endif
implicit none

save
private

public :: enter_sub, exit_sub, error, warn, mesg
public :: msg, n_l

integer, parameter :: n_blanks = 32
integer, parameter :: n_msg = 1024

character (n_blanks), parameter :: blanks = repeat (' ', n_blanks)
character (2), parameter :: n_l = achar (10) // ' '
                            ! carriage return and a space

character (n_msg) :: msg  ! experiment
character (64) :: fmt

integer, parameter :: lun = 6  ! system dependent

integer :: call_level = 0

interface error
  module procedure error_a, error_ai, error_al, error_aia, error_ai_array,  &
                   error_aiar, error_ar, error_ar_array
end interface

interface mesg
  module procedure message_a, message_ai, message_aiai, message_aiar,  &
                   message_al, message_ar, message_aii, message_air,   &
                   message_ai_array, message_aiai_array,               &
                   message_ar_array, message_aiar_array
end interface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine enter_sub (name)

character (*), intent (in) :: name

integer :: n

!---------------------------------------------------------------------

call_level = call_level + 1

n = min (n_blanks, call_level-1)

write (lun, '(1x,a)') blanks(1:n) // name // ': started'


end subroutine enter_sub

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine exit_sub (name)

character (*), intent (in) :: name

integer :: n

!---------------------------------------------------------------------

n = min (n_blanks, call_level-1)

write (lun, '(1x,a)') blanks(1:n) // name // ': done'

call_level = call_level - 1

end subroutine exit_sub
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine message_a (name, msg)

character (*), intent (in) :: name
character (*), intent (in) :: msg

write (lun, '(1x,a)') name // ': ' // trim (msg)

end subroutine message_a

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine message_ai (name, msg, i)

character (*), intent (in) :: name
character (*), intent (in) :: msg
integer, intent (in) :: i

write (lun, '(1x,a,1x,i0)') name // ': ' // trim (msg), i

end subroutine message_ai

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine message_aiai (name, msg1, i1, msg2, i2)

character (*), intent (in) :: name
character (*), intent (in) :: msg1, msg2
integer, intent (in) :: i1, i2

write (lun, '(2(1x,a,1x,i0))') name // ': ' // trim (msg1), i1,  &
                                               trim (msg2), i2

end subroutine message_aiai

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine message_aiar (name, msg1, i, msg2, r)

character (*), intent (in) :: name
character (*), intent (in) :: msg1, msg2
integer, intent (in) :: i
real (rp), intent (in) :: r

fmt = '(1x,a,1x,i0,1x,a,1x,es11.4)'

write (lun, fmt) name // ': ' // trim (msg1), i, trim (msg2), r

end subroutine message_aiar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine message_al (name, msg, l)

character (*), intent (in) :: name
character (*), intent (in) :: msg
logical, intent (in) :: l

write (lun, '(1x,a,1x,l1)') name // ': ' // trim (msg), l

end subroutine message_al

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine message_aii (name, msg, i, j)

character (*), intent (in) :: name
character (*), intent (in) :: msg
integer, intent (in) :: i, j

fmt = '(1x,a,1x,i0,",",1x,i0)'
                    !--comma to separate output
write (lun, fmt) name // ': ' // trim (msg), i, j

end subroutine message_aii

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine message_air (name, msg, i, r)

character (*), intent (in) :: name
character (*), intent (in) :: msg
integer, intent (in) :: i
real (rp), intent (in) :: r

fmt = '(1x,a,1x,i0,",",1x,es11.4)'
                    !--comma to separate output 
write (lun, fmt) name // ': ' // trim (msg), i, r

end subroutine message_air

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine message_ai_array (name, msg, i_arr)

character (*), intent (in) :: name
character (*), intent (in) :: msg
integer, intent (in) :: i_arr(:)

integer :: n

!---------------------------------------------------------------------

n = size (i_arr)
write (fmt, *) '(1x,a,', n, '(1x,i0))'
write (lun, fmt) name // ': ' // trim (msg), i_arr

end subroutine message_ai_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine message_aiai_array (name, msg1, i, msg2, i_arr)

character (*), intent (in) :: name
character (*), intent (in) :: msg1, msg2
integer, intent (in) :: i
integer, intent (in) :: i_arr(:)

integer :: n

!---------------------------------------------------------------------

n = size (i_arr)
write (fmt, *) '(1x,a,i0,a,', n, '(1x,i0))'
write (lun, fmt) name // ': ' // trim (msg1), i, trim(msg2), i_arr

end subroutine message_aiai_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine message_ar (name, msg, r)

character (*), intent (in) :: name
character (*), intent (in) :: msg
real (rp), intent (in) :: r

write (lun, '(1x,a,1x,es11.4)') name // ': ' // trim (msg), r

end subroutine message_ar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine message_ar_array (name, msg, r_arr)

character (*), intent (in) :: name
character (*), intent (in) :: msg
real (rp), intent (in) :: r_arr(:)

integer :: n

!---------------------------------------------------------------------

n = size (r_arr)
!write (fmt, '(a,i0,a)') '(1x,a,', n, '(1x,es11.4))'
write (lun, *) name // ': ' // trim (msg), r_arr

end subroutine message_ar_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine message_aiar_array (name, msg1, i, msg2, r_arr)

character (*), intent (in) :: name
character (*), intent (in) :: msg1, msg2
integer, intent (in) :: i
real (rp), intent (in) :: r_arr(:)

integer :: n

!---------------------------------------------------------------------

n = size (r_arr)
write (fmt, *) '(1x,a,i0,a,', n, '(1x,es11.4))'
write (lun, fmt) name // ': ' // trim (msg1), i, trim (msg2), r_arr

end subroutine message_aiar_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine warn (name, msg)

character (*), intent (in) :: name
character (*), intent (in) :: msg

write (lun, '(1x,a)') '*****WARNING*****'
write (lun, '(1x,a)') 'In ' // name // ':'
write (lun, '(1x,a)') trim (msg)
write (lun, '(1x,a)') '*****************'

end subroutine warn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine error_a (name, msg)

character (*), intent (in) :: name
character (*), intent (in) :: msg

write (lun, '(1x,a)') '*****ERROR*****'
write (lun, '(1x,a)') 'In ' // name // ':'
write (lun, '(1x,a)') trim (msg)
write (lun, '(1x,a)') '***************'
write (lun, '(1x,a)') 'Program aborted'

stop

end subroutine error_a
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine error_ai (name, msg, i)

character (*), intent (in) :: name
character (*), intent (in) :: msg
integer, intent (in) :: i

write (lun, '(1x,a)') '*****ERROR*****'
write (lun, '(1x,a)') 'In ' // name // ':'
write (lun, '(1x,a,1x,i0)') trim (msg), i
write (lun, '(1x,a)') '***************'
write (lun, '(1x,a)') 'Program aborted'

stop

end subroutine error_ai

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine error_ai_array (name, msg, i_arr)

character (*), intent (in) :: name
character (*), intent (in) :: msg
integer, intent (in) :: i_arr(:)

integer :: n

!---------------------------------------------------------------------

n = size (i_arr)
write (fmt, *) '(1x,a,', n, '(1x,i0))'

write (lun, '(1x,a)') '*****ERROR*****'
write (lun, '(1x,a)') 'In ' // name // ':'
write (lun, fmt) trim (msg), i_arr
write (lun, '(1x,a)') '***************'
write (lun, '(1x,a)') 'Program aborted'

stop

end subroutine error_ai_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine error_aia (name, msg1, i, msg2)

character (*), intent (in) :: name
character (*), intent (in) :: msg1, msg2
integer, intent (in) :: i

write (lun, '(1x,a)') '*****ERROR*****'
write (lun, '(1x,a)') 'In ' // name // ':'
write (lun, '(1x,a,1x,i0,1x,a)') trim (msg1), i, trim (msg2)
write (lun, '(1x,a)') '***************'
write (lun, '(1x,a)') 'Program aborted'

stop

end subroutine error_aia

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine error_aiar (name, msg1, i, msg2, r)

character (*), intent (in) :: name
character (*), intent (in) :: msg1, msg2
integer, intent (in) :: i
real (rp), intent (in) :: r

write (lun, '(1x,a)') '*****ERROR*****'
write (lun, '(1x,a)') 'In ' // name // ':'
write (lun, '(1x,a,1x,i0,1x,a,1x,es11.4)') trim (msg1), i, trim (msg2), r
write (lun, '(1x,a)') '***************'
write (lun, '(1x,a)') 'Program aborted'

stop

end subroutine error_aiar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine error_al (name, msg, l)

character (*), intent (in) :: name
character (*), intent (in) :: msg
logical, intent (in) :: l

write (lun, '(1x,a)') '*****ERROR*****'
write (lun, '(1x,a)') 'In ' // name // ':'
write (lun, '(1x,a,1x,l1)') trim (msg), l
write (lun, '(1x,a)') '***************'
write (lun, '(1x,a)') 'Program aborted'

stop

end subroutine error_al

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine error_ar (name, msg, r)

character (*), intent (in) :: name
character (*), intent (in) :: msg
real (rp), intent (in) :: r

write (lun, '(1x,a)') '*****ERROR*****'
write (lun, '(1x,a)') 'In ' // name // ':'
write (lun, '(1x,a,1x,es11.4)') trim (msg), r
write (lun, '(1x,a)') '***************'
write (lun, '(1x,a)') 'Program aborted'

stop

end subroutine error_ar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine error_ar_array (name, msg, r_arr)

character (*), intent (in) :: name
character (*), intent (in) :: msg
real (rp), intent (in) :: r_arr(:)

integer :: n

!---------------------------------------------------------------------

n = size (r_arr)
write (fmt, *) '(1x,a,', n, '(1x,es11.4))'

write (lun, '(1x,a)') '*****ERROR*****'
write (lun, '(1x,a)') 'In ' // name // ':'
write (lun, fmt) trim (msg), r_arr
write (lun, '(1x,a)') '***************'
write (lun, '(1x,a)') 'Program aborted'

stop

end subroutine error_ar_array
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module messages
