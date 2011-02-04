! this is the w-node version
!--MPI: requires u, v 0:nz, except bottom process 1:nz
subroutine interpolag_Ssim ()
! This subroutine takes the arrays F_LM and F_MM from the previous  
!   timestep and essentially moves the values around to follow the 
!   corresponding particles. The (x, y, z) value at the current 
!   timestep will be the (x-u*dt, y-v*dt, z-w*dt) value at the 
!   previous timestep.  Since particle motion does not conform to
!   the grid, an interpolation will be required.

use types,only:rprec
use param
!use param,only:ld,nx,ny,nz,dx,dy,dz,dt,cs_count,cfl_count, jt,  &
!               USE_MPI, coord, nproc, BOGUS
use sgsmodule
use messages
implicit none

character (*), parameter :: sub_name = 'interpolag_Ssim'
$if ($DEBUG)
logical, parameter :: DEBUG = .false.
$endif

integer :: jx, jy, jz
integer :: jjx, jjy, jjz
integer :: jz_min
$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif

real(kind=rprec) :: frac_x,frac_y,frac_z
real(kind=rprec) :: comp_x,comp_y,comp_z
real(kind=rprec), dimension(nx+2,ny+2,nz+2) :: FF_LM,FF_MM
real(kind=rprec) :: lagran_dt, consta
real (rprec) :: lcfl

integer :: addx, addy, addz
integer :: jxaddx, jyaddy, jzaddz

!---------------------------------------------------------------------
$if ($VERBOSE)
call enter_sub (sub_name)
$endif

! To account for particles passing through the boundaries and/or between MPI 
!   coords, "padded" arrays FF_LM and FF_MM will be created & used in place of 
!   F_LM and F_MM during interpolation.  Since particles cannot pass a grid spacing
!   during the Lag timestep this padding only needs "one" extra in each direction.
!   Finally, these padded arrays are shifted to allow for the zero index.

! Create FF_LM and FF_MM (repeat values for periodic BCs and use MPI_sendrecv to 
!   share values among coords)
    !$omp parallel do default(shared) private(jx,jy,jz)	
    do jz = 1, nz
      do jy = 1, ny
        do jx = 1, nx
          FF_LM(jx+1, jy+1, jz+1) = F_LM(jx, jy, jz)
          FF_MM(jx+1, jy+1, jz+1) = F_MM(jx, jy, jz)
        end do
      end do
    end do
    !$omp end parallel do

    !This is a bit like witch craft but the following lines do take care
    !of all the edges including the corners
    FF_LM(1,:,:) = FF_LM(nx+1,:,:)
    FF_LM(nx+2,:,:) = FF_LM(2,:,:)

    FF_LM(:,1,:) = FF_LM(:,ny+1,:) 
    FF_LM(:,ny+2,:) = FF_LM(:,2,:) 

    $if ($MPI)
      !--send F_LM @ jz=nz-1 to F_LM @ jz=0'
      !  i.e. FF_LM @ jz=nz to FF_LM @ jz=1'
      call mpi_sendrecv (FF_LM(1, 1, nz), (nx+2)*(ny+2), MPI_RPREC, up, 1,   &
                         FF_LM(1, 1, 1), (nx+2)*(ny+2), MPI_RPREC, down, 1,  &
                         comm, status, ierr)
      !--F_LM @ jz=nz and F_LM @ jz=1' should already be in sync (test?)
      !  i.e. FF_LM @ jz=nz+1 and FF_LM @ jz=2'
    $endif

    if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
      FF_LM(:,:,1) = FF_LM(:,:,2) 
    end if

    $if ($MPI)
      !--send F_LM @ jz=2 to F_LM @ jz=nz+1'
      !  i.e. FF_LM @ jz=3 to FF_LM @ jz=nz+2'
    call mpi_sendrecv (FF_LM(1, 1, 3), (nx+2)*(ny+2), MPI_RPREC, down, 2,   &
                     FF_LM(1, 1, nz+2), (nx+2)*(ny+2), MPI_RPREC, up, 2,  &
                     comm, status, ierr)
    $endif

    if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
      FF_LM(:,:,nz+2) = FF_LM(:,:,nz+1) 
    end if

    FF_MM(1,:,:) = FF_MM(nx+1,:,:)
    FF_MM(nx+2,:,:) = FF_MM(2,:,:)

    FF_MM(:,1,:) = FF_MM(:,ny+1,:) 
    FF_MM(:,ny+2,:) = FF_MM(:,2,:) 

    $if ($MPI)
      !--send F_MM @ jz=nz-1 to F_MM @ jz=0'
      !  i.e. FF_MM @ jz=nz to FF_MM @ jz=1'
      call mpi_sendrecv (FF_MM(1, 1, nz), (nx+2)*(ny+2), MPI_RPREC, up, 3,   &
                         FF_MM(1, 1, 1), (nx+2)*(ny+2), MPI_RPREC, down, 3,  &
                         comm, status, ierr)
      !--F_MM @ jz=nz and F_MM @ jz=1' should already be in sync (test?)
      !  i.e. FF_MM @ jz=nz+1 and FF_MM @ jz=2'
    $endif

    if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
      FF_MM(:,:,1) = FF_MM(:,:,2) 
    end if

    $if ($MPI)
      !--send F_MM @ jz=2 to F_MM @ jz=nz+1'
      !  i.e. FF_MM @ jz=3 to FF_MM @ jz=nz+2'
      call mpi_sendrecv (FF_MM(1, 1, 3), (nx+2)*(ny+2), MPI_RPREC, down, 4,   &
                         FF_MM(1, 1, nz+2), (nx+2)*(ny+2), MPI_RPREC, up, 4,  &
                         comm, status, ierr)
    $endif

    if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
      FF_MM(:,:,nz+2) = FF_MM(:,:,nz+1) 
    end if
    ! end of witch craft

    $if ($DEBUG)
    if (DEBUG) write (*, *) 'interpolag_Ssim: after FF setup'
    $endif

! Put u_lag and and v_lag on cs nodes
    if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
      jz_min = 2
    else
      jz_min = 1
    end if

    !--must do backwards due to data depenency on plane below current one
    do jz = nz, jz_min, -1
      u_lag(:, :, jz) = 0.5_rprec * (u_lag(:, :, jz) + u_lag(:, :, jz-1))
      v_lag(:, :, jz) = 0.5_rprec * (v_lag(:, :, jz) + v_lag(:, :, jz-1))
    end do

    $if ($MPI)
      !--not sure if 0-level is needed
      w_lag(:, :, $lbz) = BOGUS
    $endif

    if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
      w_lag (:,:,1) = .25_rprec*w_lag (:,:,2)
    end if

! Computes the 3-D inverse displacement arrays that describe 
!   the location where the point was at the previous step
! Previously u_lag represented a velocity  
! Now it represents the number of grid spacings the particle has traveled in the past cs_count timesteps
!   (motion in the positive x-direction will have a negative u_lag)
    u_lag = -u_lag*dt/dx    ! note: Previous u_lag is sum of all u's from past cs_count timesteps
    v_lag = -v_lag*dt/dy    !       so this is equivalent to taking the average then multiplying
    w_lag = -w_lag*dt/dz    !       by cs_count*dt to account for all the time that has elapsed.
                            !       Dividing by dx makes u_lag the number of gridspaces (not a distance)

    $if ($MPI)
      !--perhaps can remove 0-level altogether
      u_lag(:, :, $lbz) = BOGUS
      v_lag(:, :, $lbz) = BOGUS
      w_lag(:, :, $lbz) = BOGUS
    $endif

    if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
      !	because first plane is on u,v,p nodes
      !	this corrects for the fact that the first cell
      !	in the z direction has height dz/2
      !	it doubles the zp fraction if this fraction relates to the cell 
      !     beneath it
      do jy=1,ny
        do jx=1,nx
          w_lag(jx,jy,2) = w_lag(jx,jy,2) + min (w_lag(jx,jy,2), 0._rprec)
          w_lag(jx,jy,1) = w_lag(jx,jy,1) + max (w_lag(jx,jy,1), 0._rprec)
        end do
      end do

      !--we may be missing something here, see interpolag_Sdep at analogous place
    end if

! Compute the Lagrangian CFL number and print to screen
!   Note: this is only in the x-direction... not for complex geometry cases
    if (mod (jt, cfl_count) .eq. 0) then
      lcfl = 0._rprec
      do jz = 1, nz
        lcfl = max ( lcfl,  maxval (abs (u_lag(1:nx, :, jz))) )
      end do
      print*, 'Lagrangian CFL condition= ', lcfl
    end if

! Compute (backwards) interpolation of F_LM, F_MM (represented by "padded" arrays FF_LM, FF_MM)
!   and store the interpolated values to F_LM, F_MM arrays 
! This will be an interpolation between (i,j,k) and i+/-1 etc depending on the sign of u_lag
!   (negative means it moved in the positive x-direction so we want i-1 to go back in time)
!   ditto for v_lag and w_lag
    do jz=1,nz
     jjz = jz+1     ! accounts for shift that was used to accomodate zero index
     do jy=1,ny
      jjy = jy+1
      do jx=1,nx
         jjx = jx+1
         ! addx is the value to add to the index jx (will be either +1 or -1)
         addx = int(sign(1._rprec,u_lag(jx,jy,jz)))
         addy = int(sign(1._rprec,v_lag(jx,jy,jz)))
         addz = int(sign(1._rprec,w_lag(jx,jy,jz)))
         jxaddx = jjx + addx    ! interpolation will be between jjx and jxaddx
         jyaddy = jjy + addy
         jzaddz = jjz + addz

         ! computes the relative weights given to F_** in the cube depending on point location
         comp_x = abs(u_lag(jx,jy,jz))  ! weight for jxaddx
         comp_y = abs(v_lag(jx,jy,jz))
         comp_z = abs(w_lag(jx,jy,jz))

         frac_x = 1._rprec - comp_x     ! weight for jjx
         frac_y = 1._rprec - comp_y
         frac_z = 1._rprec - comp_z

         ! computes interpoated F_LM
         F_LM(jx,jy,jz)=frac_x*frac_y*&
              (FF_LM(jjx,jjy,jjz)*frac_z+FF_LM(jjx,jjy,jzaddz)*comp_z)&
              + frac_x*comp_y*&
              (FF_LM(jjx,jyaddy,jjz)*frac_z+FF_LM(jjx,jyaddy,jzaddz)*comp_z)&
              + comp_x*frac_y*&
              (FF_LM(jxaddx,jjy,jjz)*frac_z+FF_LM(jxaddx,jjy,jzaddz)*comp_z)&
              + comp_x*comp_y*&
              (FF_LM(jxaddx,jyaddy,jjz)*frac_z&
              +FF_LM(jxaddx,jyaddy,jzaddz)*comp_z)

         ! computes interpoated F_MM
         F_MM(jx,jy,jz)=&
              frac_x*frac_y*&
              (FF_MM(jjx,jjy,jjz)*frac_z+FF_MM(jjx,jjy,jzaddz)*comp_z)&
              + frac_x*comp_y*&
              (FF_MM(jjx,jyaddy,jjz)*frac_z+FF_MM(jjx,jyaddy,jzaddz)*comp_z)&
              + comp_x*frac_y*&
              (FF_MM(jxaddx,jjy,jjz)*frac_z+FF_MM(jxaddx,jjy,jzaddz)*comp_z)&
              + comp_x*comp_y*&
              (FF_MM(jxaddx,jyaddy,jjz)*frac_z&
              +FF_MM(jxaddx,jyaddy,jzaddz)*comp_z)
      end do
     end do
    end do

! Reset Lagrangian u/v/w variables for use during next set of cs_counts
    u_lag = 0._rprec
    v_lag = 0._rprec
    w_lag = 0._rprec

$if ($VERBOSE)
call exit_sub (sub_name)
$endif

end subroutine interpolag_Ssim
