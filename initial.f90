!*******************************************************************************
subroutine initial()
!*******************************************************************************
use types,only:rprec
use param
use sim_param,only:path,u,v,w,RHSx,RHSy,RHSz,theta,q
!use sgsmodule , only : Cs_opt2, Cs_opt2_avg, F_LM, F_MM, F_QN, F_NN 
use sgsmodule , only : Cs_opt2, F_LM, F_MM, F_QN, F_NN
$if ($DYN_TN)
use sgsmodule, only:F_ee2,F_deedt2,ee_past
$endif

use scalars_module,only:RHS_T ! added by VK
use scalars_module2,only:ic_scal ! added by VK

! VK -label 122 assigned to vel_sc.out for reading input files in case of
! scalars
!!!!XXXXXXXXXX--------Added by Vijayant----XXXXXXX!!!!!

use immersedbc,only:fx,fy,fz
use immersedbc,only:fxa,fya,fza

$if ($MPI)
  use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWNUP
$endif

implicit none

$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif

logical, parameter :: use_add_random = .false.

character (64) :: fname, temp

integer::i,jz

!real (rprec) :: ke

!Cs_opt2_avg=0._rprec
fx=0._rprec;fy=0._rprec;fz=0._rprec
fxa=0._rprec; fya=0._rprec; fza=0._rprec

$if ($DYN_TN)
!Eventually want to read these in from file, but for now, just initialize to zero
ee_past = 0.1_rprec; F_ee2 = 10.0_rprec; F_deedt2 = 10000.0_rprec
$endif

fname = path // 'vel.out'
$if ($MPI)
  write (temp, '(".c",i0)') coord
  fname = trim (fname) // temp
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
  
  if(initsc) then
    if(.not. USE_MPI .or. (USE_MPI .and. coord == 0) ) write(*,*) '--> Reading initial velocity and temperature from file'
    read(12) u,v,w,RHSx,RHSy,RHSz,Cs_opt2,F_LM,F_MM,theta,RHS_T
!TS INITIALIZE THE ZERO CONCENTRATION FIELD IF jt_total=0
    if(sflux_flag)then
       open(1,file=path//'run')
       read(1,*)i
       close(1)
       if(i.eq.0)then
          theta=0._rprec;RHS_T=0._rprec
       endif
    endif
    do jz=1,nz
    write(6,7781) jz, sum (u(1:nx, :, jz)) / (nx * ny),  &
                      sum (v(1:nx, :, jz)) / (nx * ny),  &
                      sum (w(1:nx, :, jz)) / (nx * ny),  &
                      sum (theta(1:nx, :, jz)) / (nx * ny)
    end do
7781 format('jz, ubar, vbar, wbar, Tbar:',(1x,I3,1x,F9.4,1x,F9.4,1x,F9.4,1x,F9.4))
  else

    if(.not. USE_MPI .or. (USE_MPI .and. coord == 0) ) write(*,*) '--> Reading initial velocity field from file'

    select case (model)
      case (1)
        read (12) u(:, :, 1:nz), v(:, :, 1:nz), w(:, :, 1:nz),       &
                  RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz)
      case (2:3)
        read(12) u(:, :, 1:nz),v(:, :, 1:nz),w(:, :, 1:nz),             &
                 RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),  &
                 Cs_opt2(:,:,1:nz)              
      case (4)
        if (inilag) then  !--not sure if Cs_opt2 should be there, just quickie
          read (12) u(:, :, 1:nz),v(:, :, 1:nz),w(:, :, 1:nz),             &
                    RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),  &
                    Cs_opt2(:,:,1:nz)
        else
        $if($DYN_TN)
          read(12) u(:, :, 1:nz),v(:, :, 1:nz),w(:, :, 1:nz),             &
                   RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),  &
                   Cs_opt2(:,:,1:nz), F_LM(:,:,1:nz), F_MM(:,:,1:nz),     &
                   F_ee2(:,:,1:nz), F_deedt2(:,:,1:nz), ee_past(:,:,1:nz)
        $else
          read(12) u(:, :, 1:nz),v(:, :, 1:nz),w(:, :, 1:nz),             &
                   RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),  &
                   Cs_opt2(:,:,1:nz), F_LM(:,:,1:nz), F_MM(:,:,1:nz)   
        $endif
        end if
      case (5)
        if (inilag) then  !--not sure if Cs_opt2 should be there, just quickie
          read (12) u(:, :, 1:nz),v(:, :, 1:nz),w(:, :, 1:nz),             &
                    RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),  &
                    Cs_opt2(:,:,1:nz)
        else
          read(12) u(:, :, 1:nz),v(:, :, 1:nz),w(:, :, 1:nz),             &
                   RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),  &
                   Cs_opt2(:,:,1:nz), F_LM(:,:,1:nz), F_MM(:,:,1:nz),     &
                   F_QN(:,:,1:nz), F_NN(:,:,1:nz)
        end if
      case default
        write (*, *) 'initial: invalid model number'
    end select

    !call energy (ke)

    do jz=1,nz
      write(6,7780) jz, sum (u(1:nx, :, jz)) / (nx * ny),  &
                        sum (v(1:nx, :, jz)) / (nx * ny),  &
                        sum (w(1:nx, :, jz)) / (nx * ny)
    end do
7780 format('jz, ubar, vbar, wbar:',(1x,I3,1x,F9.4,1x,F9.4,1x,F9.4))
  end if

  close(12) 

else
  if (dns_bc) then
     if(.not. USE_MPI .or. (USE_MPI .and. coord == 0) ) write(*,*) '--> Creating initial velocity field with DNS BCs'
     call ic_dns()
  else
    if(.not. USE_MPI .or. (USE_MPI .and. coord == 0) ) write(*,*) '--> Creating initial fields'
    if (S_FLAG) then
       if(.not. USE_MPI .or. (USE_MPI .and. coord == 0) ) write(*,*) '----> Creating initial velocity & scalar fields'
       call ic_scal()
    else
       if(.not. USE_MPI .or. (USE_MPI .and. coord == 0) ) write(*,*) '----> Creating initial velocity field'
       call ic()
    end if
  end if
end if

$if ($MPI)
  !--synchronize the overlapping parts nz-1 (coord) -> 0 (coord + 1) 
  !  and 1 (coord + 1) -> nz (coord)
  !call mpi_sendrecv (u(1, 1, nz-1), ld*ny, MPI_RPREC, up, 1,  &
  !                   u(1, 1, 0), ld*ny, MPI_RPREC, down, 1,   &
  !                   comm, status, ierr)
  !call mpi_sendrecv (v(1, 1, nz-1), ld*ny, MPI_RPREC, up, 2,  &
  !                   v(1, 1, 0), ld*ny, MPI_RPREC, down, 2,   &
  !                   comm, status, ierr)
  !call mpi_sendrecv (w(1, 1, nz-1), ld*ny, MPI_RPREC, up, 3,  &
  !                   w(1, 1, 0), ld*ny, MPI_RPREC, down, 3,   &
  !                   comm, status, ierr)
  !call mpi_sendrecv (u(1, 1, 1), ld*ny, MPI_RPREC, down, 4,  &
  !                   u(1, 1, nz), ld*ny, MPI_RPREC, up, 4,   &
  !                   comm, status, ierr)
  !call mpi_sendrecv (v(1, 1, 1), ld*ny, MPI_RPREC, down, 5,  &
  !                   v(1, 1, nz), ld*ny, MPI_RPREC, up, 5,   &
  !                   comm, status, ierr)
  !call mpi_sendrecv (w(1, 1, 1), ld*ny, MPI_RPREC, down, 6,  &
  !                   w(1, 1, nz), ld*ny, MPI_RPREC, up, 6,   &
  !                   comm, status, ierr)   

  call mpi_sync_real_array( u, MPI_SYNC_DOWNUP )
  call mpi_sync_real_array( v, MPI_SYNC_DOWNUP ) 
  call mpi_sync_real_array( w, MPI_SYNC_DOWNUP ) 
  
$endif

if (USE_MPI .and. coord == 0) then
  !--set 0-level velocities to BOGUS
  u(:, :, $lbz) = BOGUS
  v(:, :, $lbz) = BOGUS
  w(:, :, $lbz) = BOGUS
end if

!  Open vel.out (lun_default in io) for final output
$if ($WRITE_BIG_ENDIAN)
open(11,file=fname,form='unformatted', convert='big_endian')
$elseif ($WRITE_LITTLE_ENDIAN)
open(11,file=fname,form='unformatted', convert='little_endian')
$else
open(11,file=fname,form='unformatted')
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
