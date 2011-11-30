!*******************************************************************************
subroutine initial()
!*******************************************************************************
use types,only:rprec
use param
use sim_param,only:path,u,v,w,RHSx,RHSy,RHSz,theta,q
use sgs_param, only : Cs_opt2, F_LM, F_MM, F_QN, F_NN
$if ($DYN_TN)
use sgs_param, only:F_ee2,F_deedt2,ee_past
$endif

use sim_param,only:fx,fy,fz
use sim_param,only:fxa,fya,fza

$if ($MPI)
  use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWNUP
$endif

implicit none

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
  

    if(.not. USE_MPI .or. (USE_MPI .and. coord == 0) ) write(*,*) '--> Reading initial velocity field from file'

    select case (sgs_model)
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
        write (*, *) 'initial: invalid sgs_model number'
    end select

    !call energy (ke)

    do jz=1,nz
      write(6,7780) jz, sum (u(1:nx, :, jz)) / (nx * ny),  &
                        sum (v(1:nx, :, jz)) / (nx * ny),  &
                        sum (w(1:nx, :, jz)) / (nx * ny)
    end do
7780 format('jz, ubar, vbar, wbar:',(1x,I3,1x,F9.4,1x,F9.4,1x,F9.4))

  close(12) 

else
  if (dns_bc) then
     if(.not. USE_MPI .or. (USE_MPI .and. coord == 0) ) write(*,*) '--> Creating initial velocity field with DNS BCs'
     call ic_dns()
  else
    if(.not. USE_MPI .or. (USE_MPI .and. coord == 0) ) write(*,*) '--> Creating initial fields'
       if(.not. USE_MPI .or. (USE_MPI .and. coord == 0) ) write(*,*) '----> Creating initial velocity field'
       call ic()
  end if
end if

$if ($MPI)

  call mpi_sync_real_array( u, 0, MPI_SYNC_DOWNUP )
  call mpi_sync_real_array( v, 0, MPI_SYNC_DOWNUP ) 
  call mpi_sync_real_array( w, 0, MPI_SYNC_DOWNUP ) 
  
$endif

if (USE_MPI .and. coord == 0) then
  !--set 0-level velocities to BOGUS
  u(:, :, lbz) = BOGUS
  v(:, :, lbz) = BOGUS
  w(:, :, lbz) = BOGUS
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
