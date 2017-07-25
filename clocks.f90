!!
!!  Copyright (C) 2011-2016  Johns Hopkins University
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

!*********************************************************************
module clock_m
!*********************************************************************
!
! This module provides the clock data type (object) and the
! subroutines/functions that act on instances of the clock data type.
!
use types, only : rprec
implicit none

save 
private

public clock_t

type clock_t
    real(rprec) :: start_time
    real(rprec) :: stop_time
    real(rprec) :: time
contains
    procedure, public :: start
    procedure, public :: stop
end type clock_t

contains

!*********************************************************************
subroutine start( this )
!*********************************************************************
#ifdef PPMPI
use mpi, only : mpi_wtime
#endif
implicit none

class(clock_t), intent(inout) :: this

#ifdef PPMPI
this % start_time = mpi_wtime()
#else
call cpu_time( this % start_time )
#endif

return
end subroutine start

!*********************************************************************
subroutine stop( this )
!*********************************************************************
#ifdef PPMPI
use mpi, only : mpi_wtime
#endif
implicit none

class(clock_t), intent(inout) :: this

#ifdef PPMPI
this % stop_time = mpi_wtime()
#else
call cpu_time( this % stop_time )
#endif

! Compute the clock time
this % time = this % stop_time - this % start_time

return

end subroutine stop

end module clock_m
