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

!*******************************************************************************
subroutine rmsdiv(rms)
!*******************************************************************************
!
! This subroutine calculates the velocity divergence metric. Currently using the
! L_1 norm ov the velocity divergence.
!
use param, only : rprec
use param
use sim_param, only : dudx, dvdy, dwdz

implicit none
integer :: jx, jy, jz, jz_max
real(rprec) :: rms
#ifdef PPMPI
real(rprec) :: rms_global
#endif

! Initialize variables
rms = 0._rprec
jz_max = nz - 1

! Calculate L1 norm of velocity divergence
do jz = 1, jz_max
do jy = 1, ny
do jx = 1, nx
    rms = rms + abs( dudx(jx,jy,jz) + dvdy(jx,jy,jz) + dwdz(jx,jy,jz) )
end do
end do
end do
rms = rms / (nx*ny*(jz_max))

#ifdef PPMPI
! Transfer between processors
call mpi_reduce(rms, rms_global, 1, MPI_RPREC, MPI_SUM, 0, comm, ierr)
if (rank == 0) then
    rms = rms_global/nproc
end if
#endif

end subroutine rmsdiv
