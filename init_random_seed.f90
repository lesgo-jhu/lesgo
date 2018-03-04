!!
!!  Copyright (C) 2016  Johns Hopkins University
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

!*************************************************************
subroutine init_random_seed()
!*************************************************************
!Restarts of queries the state of the pseudorandom number 
!generator used by RANDOM_NUMBER
!
!Call init_random_seed before any calls of random_number.
!
! This subroutine is adapted from gcc:
!https://gcc.gnu.org/onlinedocs/gcc-4.2.3/gfortran/RANDOM_005fSEED.html#RANDOM_005fSEED
!
implicit none
integer :: i, n, clock
integer, dimension(:), allocatable :: seed

call random_seed(size=n)
allocate(seed(n))

call system_clock(count=clock)

seed = clock + 37 * (/ (i -1, i = 1, n) /)
call random_seed(put=seed)

deallocate(seed)

end subroutine init_random_seed
