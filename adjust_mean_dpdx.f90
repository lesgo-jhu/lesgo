!!
!!  Copyright (C) 2010-2016  Johns Hopkins University
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
module adjust_mean_dpdx
!*******************************************************************************

use types, only : rprec
use param
use grid_m
use messages
use string_util
#ifdef PPMPI
use mpi_defs, only : MPI_SYNC_DOWNUP, mpi_sync_real_array 
#endif
use sim_param, only : u 

implicit none

save
private

public :: adjust_dpdx

real(rprec) :: ubulk_old

! Commonly used indices
integer :: i,j,k

contains

!*******************************************************************************
subroutine adjust_dpdx()
!*******************************************************************************

implicit none

real(rprec) :: ubulk_target, ubulk
real(rprec) :: ubulk_global

ubulk_target = 1.0_rprec

ubulk = 0.0_rprec
do k=1, nz-1
   do j=1, ny
      do i= 1, nx
         ubulk = ubulk + u(i,j,k)*dz
      enddo
   enddo
enddo
ubulk = ubulk / (real(nx)*real(ny)*L_z)

#ifdef PPMPI
call mpi_allreduce(ubulk,ubulk_global,1,MPI_RPREC,MPI_SUM, &
     MPI_COMM_WORLD,ierr)
#endif

! initialize if restarting
if (jt == 1 .and. jt_total /= 1) then
   call read_restart
endif

if (jt_total == 1) then
   mean_p_force_x = 0._rprec
   ubulk_old = ubulk_global
endif

mean_p_force_x = mean_p_force_x - ( 1.5_rprec * ( ubulk_global - ubulk_target ) / dt &
     - 0.5_rprec * ( ubulk_old - ubulk_target ) / dt )

ubulk_old = ubulk_global

if (jt_total == nsteps .and. coord == 0) then
   call write_restart
endif

end subroutine adjust_dpdx


!*******************************************************************************
subroutine write_restart()
!*******************************************************************************

open(2,file='constant_mass_restart.txt')
write(2,*) ubulk_old, mean_p_force_x
close(2)

end subroutine write_restart

!*******************************************************************************
subroutine read_restart()
!*******************************************************************************

open(2,file='constant_mass_restart.txt')
read(2,*) ubulk_old, mean_p_force_x
close(2)

end subroutine read_restart


end module adjust_mean_dpdx

