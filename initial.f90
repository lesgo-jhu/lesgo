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

!*******************************************************************************
subroutine initial()
!*******************************************************************************
use types,only:rprec
use param
use sim_param, only : u,v,w,RHSx,RHSy,RHSz,theta,q
use sgs_param, only : Cs_opt2, F_LM, F_MM, F_QN, F_NN
$if ($DYN_TN)
use sgs_param, only:F_ee2,F_deedt2,ee_past
$endif
$if (TURBINES)
use sim_param,only:fxa
$endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Tony ATM
!$if ($ATM)
!use sim_param,only:fxa, fya, fza
!$endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Tony ATM
$if ($LVLSET)
use sim_param,only:fxa,fya,fza
use sim_param,only:fx,fy,fz
$endif
use string_util, only : string_concat
$if ($MPI)
  use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWNUP
$endif

implicit none

logical, parameter :: use_add_random = .false.

character (64) :: fname

$if ($DYN_TN)
logical :: exst
character (64) :: fname_dyn_tn
$endif

integer::jz

$if ($TURBINES)
fxa=0._rprec
$endif

$if ($LVLSET)
fx=0._rprec;fy=0._rprec;fz=0._rprec
fxa=0._rprec; fya=0._rprec; fza=0._rprec
$endif

$if ($DYN_TN)
!Will be over-written if read from dyn_tn.out files
ee_past = 0.1_rprec; F_ee2 = 10.0_rprec; F_deedt2 = 10000.0_rprec
fname_dyn_tn = path // 'dyn_tn.out'
  $if ($MPI)
  call string_concat( fname_dyn_tn, '.c', coord )
  $endif
$endif

fname = checkpoint_file

$if ($MPI)
call string_concat( fname, '.c', coord )
$endif

!TSopen(12,file=path//'vel_sc.out',form='unformatted')

if(initu)then

    $if ($READ_BIG_ENDIAN)
    open(12,file=fname,form='unformatted', convert='big_endian')
    $elseif ($READ_LITTLE_ENDIAN)
    open(12,file=fname,form='unformatted', convert='little_endian')
    $else
    open(12,file=fname,form='unformatted')
    $endif
  
    if(.not. USE_MPI .or. (USE_MPI .and. coord == 0) ) write(*,*) '--> Reading initial velocity field from file'

    read(12) u(:, :, 1:nz), v(:, :, 1:nz), w(:, :, 1:nz),           &
             RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),  &
             Cs_opt2(:,:,1:nz), F_LM(:,:,1:nz), F_MM(:,:,1:nz),     &
             F_QN(:,:,1:nz), F_NN(:,:,1:nz)        

    $if ($DYN_TN)
    ! Read dynamic timescale running averages from file

      if (cumulative_time) then

        inquire (file=fname_dyn_tn, exist=exst)
        if (exst) then

            $if ($READ_BIG_ENDIAN)
            open(13,file=fname_dyn_tn,form='unformatted', convert='big_endian')
            $elseif ($READ_LITTLE_ENDIAN)
            open(13,file=fname_dyn_tn,form='unformatted', convert='little_endian')
            $else
            open(13,file=fname_dyn_tn,form='unformatted')
            $endif

            read(13) F_ee2(:,:,1:nz), F_deedt2(:,:,1:nz), ee_past(:,:,1:nz)

        else

            write(*,*) trim(fname_dyn_tn), ' not found - using default values'

        endif

      endif

    $endif

    !call energy (ke)

    do jz=1,nz
      write(6,7780) jz, sum (u(1:nx, :, jz)) / (nx * ny),  &
                        sum (v(1:nx, :, jz)) / (nx * ny),  &
                        sum (w(1:nx, :, jz)) / (nx * ny)
    end do
7780 format('jz, ubar, vbar, wbar:',(1x,I3,1x,F9.4,1x,F9.4,1x,F9.4))

  close(12) 
  close(13)

else
  if (dns_bc) then
     if (coord == 0) write(*,*) '--> Creating initial velocity field with DNS BCs'
     call ic_dns()
  else
    if (coord == 0) write(*,*) '--> Creating initial fields'
       if (coord == 0) write(*,*) '----> Creating initial velocity field'
       call ic()
  end if
end if

$if ($MPI)

  call mpi_sync_real_array( u, 0, MPI_SYNC_DOWNUP )
  call mpi_sync_real_array( v, 0, MPI_SYNC_DOWNUP ) 
  call mpi_sync_real_array( w, 0, MPI_SYNC_DOWNUP ) 
  
  if (coord == 0) then
    !--set 0-level velocities to BOGUS
    u(:, :, lbz) = BOGUS
    v(:, :, lbz) = BOGUS
    w(:, :, lbz) = BOGUS
  end if
$endif

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine add_random ()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

real (rprec), parameter :: rms = 0.2_rprec

integer :: i, j, k
integer :: seed

real (rprec) :: noise
real (rprec) :: ran3

!---------------------------------------------------------------------

seed = -80

do k = 1, nz
  do j = 1, ny
    do i = 1, nx
      noise=rms/.289_rprec*(ran3(seed)-0.5_rprec)
      u(i, j, k) = u(i, j, k) + noise
      noise=rms/.289_rprec*(ran3(seed)-0.5_rprec)
      v(i, j, k) = v(i, j, k) + noise
      noise=rms/.289_rprec*(ran3(seed)-0.5_rprec)
      w(i, j, k) = w(i, j, k) + noise
    end do
  end do
end do

end subroutine add_random

end subroutine initial
