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

module debug_mod
use types, rp=>rprec
implicit none

save
private
public :: DEBUG_write

character (*), parameter :: mod_name = 'debug_mod'
!character (*), parameter :: r_fmt = 'g20.13'
character (*), parameter :: r_fmt = 'g12.5'
character (*), parameter :: c_fmt = '"(",' // r_fmt // ',", ",' //  &
                                    r_fmt // ',")"'

integer, parameter :: fname_len = 64

logical, parameter :: use_tecplot = .false.
logical, parameter :: write_header = .true.

character (fname_len) :: file

logical :: opn

interface DEBUG_write
  module procedure DEBUG_write_carray1, DEBUG_write_carray2,  &
                   DEBUG_write_carray3,                       &
                   DEBUG_write_rarray1, DEBUG_write_rarray2,  &
                   DEBUG_write_rarray3,                       &
                   DEBUG_write_r
end interface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function fname (tag)
use param, only : jt_total, coord
implicit none

character (fname_len) :: fname

character (*), intent (in) :: tag

!---------------------------------------------------------------------

$if ($MPI)
  write (fname, '(3a,i0,a,i0)') 'debug.', trim (tag), '.', jt_total,  &
                                '.MPI.c', coord
$else
  write (fname, '(3a,i0)') 'debug.', trim (tag), '.', jt_total
$endif

if (len (trim (fname)) == fname_len) then
  write (*, *) mod_name // ': warning fname_len is too small'
end if

end function fname

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine DEBUG_write_carray1 (a, tag)
implicit none

complex (rp), intent (in) :: a(:)

character (*), intent (in) :: tag

character (*), parameter :: fmt = '(i0,1x,":",' // c_fmt // ')'

integer :: i

!---------------------------------------------------------------------

inquire (unit=1, opened=opn)
if (opn) write (*, *) mod_name // ' : unit 1 is already open'

file = fname (tag)
open (1, file=file)

if (write_header) then
  write (1, '(a,1x,i0)') '#debug_mod', size (a)
end if

if (use_tecplot) then
  write (1, *) 'zone, f=point, i=', size (a)
end if

do i = 1, size (a)
  write (1, fmt) i, a(i)
end do

close (1)

end subroutine DEBUG_write_carray1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine DEBUG_write_carray2 (a, tag)
implicit none

complex (rp), intent (in) :: a(:, :)

character (*), intent (in) :: tag

character (*), parameter :: fmt = '(2(i0,1x),":",' // c_fmt // ')'

integer :: i, j

!---------------------------------------------------------------------

inquire (unit=1, opened=opn)
if (opn) write (*, *) mod_name // ' : unit 1 is already open'

file = fname (tag)
open (1, file=file)

if (write_header) then
  write (1, '(a,2(1x,i0))') '#debug_mod', size (a, 1), size (a, 2)
end if

if (use_tecplot) then
  write (1, *) 'zone, f=point, i=', size (a, 1), ', j=', size (a, 2)
end if

do j = 1, size (a, 2)
  do i = 1, size (a, 1)
    write (1, fmt) i, j, a(i, j)
  end do
end do

close (1)

end subroutine DEBUG_write_carray2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine DEBUG_write_carray3 (a, tag)
implicit none

complex (rp), intent (in) :: a(:, :, :)

character (*), intent (in) :: tag

character (*), parameter :: fmt = '(3(i0,1x),":",' // c_fmt // ')'

integer :: i, j, k

!---------------------------------------------------------------------

inquire (unit=1, opened=opn)
if (opn) write (*, *) mod_name // ' : unit 1 is already open'

file = fname (tag)
open (1, file=file)

if (write_header) then
  write (1, '(a,3(1x,i0))') '#debug_mod', size (a, 1), size (a, 2), size (a, 3)
end if

if (use_tecplot) then
  write (1, *) 'zone, f=point, i=', size (a, 1), ', j=', size (a,2),  &
               ', k=', size (a, 3)
end if

do k = 1, size (a, 3)
  do j = 1, size (a, 2)
    do i = 1, size (a, 1)
      write (1, fmt) i, j, k, a(i, j, k)
    end do
  end do
end do

close (1)

end subroutine DEBUG_write_carray3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine DEBUG_write_rarray1 (a, tag)
implicit none

real (rp), intent (in) :: a(:)

character (*), intent (in) :: tag

character (*), parameter :: fmt = '(i0,1x,":",' // r_fmt // ')'

integer :: i

!---------------------------------------------------------------------

inquire (unit=1, opened=opn)
if (opn) write (*, *) mod_name // ' : unit 1 is already open'

file = fname (tag)
open (1, file=file)

if (write_header) then
  write (1, '(a,1x,i0)') '#debug_mod', size (a, 1)
end if

if (use_tecplot) then
  write (1, *) 'zone, f=point, i=', size (a)
end if

do i = 1, size (a)
  write (1, fmt) i, a(i)
end do

close (1)

end subroutine DEBUG_write_rarray1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine DEBUG_write_rarray2 (a, tag)
implicit none

real (rp), intent (in) :: a(:, :)

character (*), intent (in) :: tag

character (*), parameter :: fmt = '(2(i0,1x),":",' // r_fmt // ')'

integer :: i, j

!---------------------------------------------------------------------

inquire (unit=1, opened=opn)
if (opn) write (*, *) mod_name // ' : unit 1 is already open'

file = fname (tag)
open (1, file=file)

if (write_header) then
  write (1, '(a,2(1x,i0))') '#debug_mod', size (a, 1), size (a, 2)
end if

if (use_tecplot) then
  write (1, *) 'zone, f=point, i=', size (a, 1), ', j=', size (a, 2)
end if

do j = 1, size (a, 2)
  do i = 1, size (a, 1)
    write (1, fmt) i, j, a(i, j)
  end do
end do

close (1)

end subroutine DEBUG_write_rarray2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine DEBUG_write_rarray3 (a, tag)
implicit none

real (rp), intent (in) :: a(:, :, :)

character (*), intent (in) :: tag

character (*), parameter :: fmt = '(3(i0,1x),":",' // r_fmt // ')'

integer :: i, j, k

!---------------------------------------------------------------------

inquire (unit=1, opened=opn)
if (opn) write (*, *) mod_name // ' : unit 1 is already open'

file = fname (tag)
open (1, file=file)

if (write_header) then
  write (1, '(a,3(1x,i0))') '#debug_mod', size (a, 1), size (a, 2), size (a, 3)
end if

if (use_tecplot) then
  write (1, *) 'zone, f=point, i=', size (a, 1), ', j=', size (a, 2),  &
               ', k=', size (a, 3)
end if

do k = 1, size (a, 3)
  do j = 1, size (a, 2)
    do i = 1, size (a, 1)
      write (1, fmt) i, j, k, a(i, j, k)
    end do
  end do
end do

close (1)

end subroutine DEBUG_write_rarray3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine DEBUG_write_r (a, tag)
implicit none

real (rp), intent (in) :: a

character (*), intent (in) :: tag

character (*), parameter :: fmt = '(' // r_fmt // ')'

!---------------------------------------------------------------------

inquire (unit=1, opened=opn)
if (opn) then
  write (*, *) mod_name // ' : unit 1 is already open'
  stop
end if

file = fname (tag)
open (1, file=file)

if (write_header) then
  write (1, '(a)') '#debug_mod'
endif

write (1, fmt) a

close (1)

end subroutine DEBUG_write_r

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module debug_mod
