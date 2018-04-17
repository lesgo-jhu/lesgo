!!
!!  Copyright (C) 2016-2017  Johns Hopkins University
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

!*******************************************************************************
subroutine init_random_seed()
!*******************************************************************************
!
! This subroutine initializes the intrinsic pseudo random number
! generator. Call init_random_seed before any calls of random_number.
!
! This subroutine is taken from gcc:
! https://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html#RANDOM_005fSEED
!
#ifdef __INTEL_COMPILER
use ifport, only : getpid
#endif
use iso_fortran_env, only : int64
implicit none
integer, allocatable :: seed(:)
integer :: i, n, un, istat, dt(8), pid
integer(int64) :: t

call random_seed(size = n)
allocate(seed(n))
! First try if the OS provides a random number generator
open(newunit=un, file="/dev/urandom", access="stream",                         &
     form="unformatted", action="read", status="old", iostat=istat)
if (istat == 0) then
   read(un) seed
   close(un)
else
   ! Fallback to XOR:ing the current time and pid. The PID is
   ! useful in case one launches multiple instances of the same
   ! program in parallel.
    call system_clock(t)
   if (t == 0) then
      call date_and_time(values=dt)
      t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000                     &
           + dt(2) * 31_int64 * 24 * 60 * 60 * 1000                            &
           + dt(3) * 24_int64 * 60 * 60 * 1000                                 &
           + dt(5) * 60 * 60 * 1000                                            &
           + dt(6) * 60 * 1000 + dt(7) * 1000                                  &
           + dt(8)
   end if
   pid = getpid()
   t = ieor(t, int(pid, kind(t)))
   do i = 1, n
      seed(i) = lcg(t)
   end do
end if

call random_seed(put=seed)

contains

!*******************************************************************************
function lcg(s)
!*******************************************************************************
!
! This is a simple pseudo random number generator. It might not be good
! enough for real work, but is sufficient for seeding a better PRNG.
!
integer :: lcg
integer(int64) :: s
if (s == 0) then
    s = 104729
else
    s = mod(s, 4294967296_int64)
end if
s = mod(s * 279470273_int64, 4294967291_int64)
lcg = int(mod(s, int(huge(0), int64)), kind(0))

end function lcg

end subroutine init_random_seed
