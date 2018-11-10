!!
!!  Copyright (C) 2009-2017  Johns Hopkins University
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
module types
!*********************************************************************
!
! This module provides generic types
!
use iso_c_binding
use mpi
implicit none

public

! Specify precision (make compliant with FFTW)
integer, parameter :: rprec = c_double
integer, parameter :: cprec = c_double_complex
integer, parameter :: MPI_RPREC = MPI_DOUBLE_PRECISION
integer, parameter :: MPI_CPREC = MPI_DOUBLE_COMPLEX

type point3D_t
    real(rprec), dimension(3) :: xyz
end type point3D_t

end module types
