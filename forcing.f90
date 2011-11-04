!**********************************************************************
subroutine forcing_applied()
!**********************************************************************
!
!  This subroutine acts as a driver for applying pointwise body forces
!  into the domain. Subroutines contained here should modify f{x,y,z}a
!  which are explicitly applied forces. These forces are applied to RHS 
!  in the evaluation of u* so that mass conservation is preserved.
!
use types, only : rprec
use immersedbc, only : fxa, fya, fza

$if ($LVLSET)
$if ($RNS_LS)
use rns_ls, only : rns_forcing_ls
$endif
$endif

$if ($TURBINES)
use turbines, only:turbines_forcing
$endif

implicit none

! Reset applied force arrays
fxa = 0._rprec
fya = 0._rprec
fza = 0._rprec

$if ($LVLSET)
$if ($RNS_LS)
call rns_forcing_ls()
$endif
$endif

$if ($TURBINES)
call turbines_forcing ()
$endif
   
end subroutine forcing_applied

!**********************************************************************
subroutine forcing_induced()
!**********************************************************************
!  
!  These forces are designated as induced forces such that they are 
!  chosen to obtain a desired velocity at time
!  step m+1. If this is not the case, care should be taken so that the forces
!  here are divergence free in order to preserve mass conservation. For 
!  non-induced forces such as explicitly applied forces they should be 
!  placed in forcing_applied.
!  
use types, only : rprec
use immersedbc, only : fx, fy, fz
$if ($LVLSET)
use level_set, only : level_set_forcing
$if($RNS_LS)
use rns_ls, only : rns_elem_force_ls
$endif
$if ($TREES_LS)
  use trees_ls
$endif
$endif
implicit none


! Initialize
fx = 0._rprec
fy = 0._rprec
fz = 0._rprec

$if($LVLSET)

!  Compute the level set IBM forces
call level_set_forcing ()

$if($TREES_LS)
!--this must come after call to level_set_forcing
!--in /a posteriori/ test, this adds SGS branch force
!--in /a priori/ test, this does not modify force
call trees_ls_calc ()
$endif

$endif

return
end subroutine forcing_induced

!**********************************************************************
subroutine inflow_cond ()
!**********************************************************************
!
!  Enforces prescribed inflow condition. Options are either a uniform
!  inflow velocity or an inlet velocity field generated from a precursor
!  simulation. The inflow condition is enforced using either an IBM type
!  forcing or by modulating directly the velocity in the fringe region.
!
use types, only : rprec
use param, only : face_avg, nx, ny, nz, pi, read_inflow_file,      &
                  buff_end, buff_len, &
                  L_x, dt, dx
use param, only : inflow_sample_velocity, inflow_sample_location, coord
use sim_param, only : u, v, w, theta
use immersedbc, only : fx, fy, fz
use io, only : inflow_read
use messages, only : error
implicit none

character (*), parameter :: sub_name = 'inflow_cond'

integer :: i, i_w
integer :: istart, istart_w
integer :: imid
integer :: iend, iend_w

integer :: icount

real (rprec) :: alpha, beta

type vel_sample_type
   integer :: n
   integer :: istart, iend
   integer :: icount, isample
end type vel_sample_type

type( vel_sample_type ) :: vel_sample_t

!--these may be out of 1, ..., nx
iend = floor (buff_end * nx + 1._rprec)
imid = floor (( buff_end - buff_len / 4 ) * nx + 1._rprec)
istart = floor ((buff_end - buff_len) * nx + 1._rprec)

!--wrapped versions
iend_w = modulo (iend - 1, nx) + 1
istart_w = modulo (istart - 1, nx) + 1

!--read from file
if (read_inflow_file) then  !--read vel inflow @ jx = iend_w from file

   call inflow_read ()  !--this sets u, v, w at (iend_w,:,:)

elseif( inflow_sample_velocity ) then  
   
   ! Zero counter
   vel_sample_t % icount=0

   ! Sample size same as buffer region
   vel_sample_t % n = iend - istart + 1

   vel_sample_t % iend = floor (inflow_sample_location * nx ) + 1
   vel_sample_t % istart = vel_sample_t % iend - vel_sample_t % n + 1

   ! if( vel_sample_t % iend < 1 .or. vel_sample_t % iend > nx .or. &
   !      vel_sample_t % istart < 1 .or. vel_sample_t % istart > nx ) then
   !    call error(sub_name,'Sample region extends beyond domain boundaries')
   ! endif

else     

   ! Use laminar inflow
   u(iend_w, :, :) = face_avg
   v(iend_w, :, :) = 0._rprec
   w(iend_w, :, :) = 0._rprec

end if

!--skip istart since we know vel at istart, iend already
do i = istart + 1, iend - 1

  i_w = modulo (i - 1, nx) + 1

  ! Linear profile
  !factor = real ( i - istart, rprec ) / real ( iend - istart, rprec )
  ! Sine profile
  !factor = 0.5_rprec * ( 1._rprec - cos (pi * real (i - istart, rprec)  &
  !                                       / (iend - istart)) )
  ! Sine profile with plateau
  if ( i > imid ) then 
     beta = 1._rprec
  else
     beta = 0.5_rprec * ( 1._rprec - cos (pi * real (i - istart, rprec)  &
          / (imid - istart)) )
  endif

  alpha = 1.0_rprec - beta

  if( inflow_sample_velocity ) then

     vel_sample_t % icount = vel_sample_t % icount + 1
     vel_sample_t % isample = modulo( vel_sample_t % istart + vel_sample_t % icount - 1, nx ) + 1

     u(i_w, 1:ny, 1:nz) = alpha * u(i_w, 1:ny, 1:nz) + beta * u(vel_sample_t % isample, 1:ny, 1:nz) 
     v(i_w, 1:ny, 1:nz) = alpha * v(i_w, 1:ny, 1:nz) + beta * v(vel_sample_t % isample, 1:ny, 1:nz)
     w(i_w, 1:ny, 1:nz) = alpha * w(i_w, 1:ny, 1:nz) + beta * w(vel_sample_t % isample, 1:ny, 1:nz)

  else

     u(i_w, 1:ny, 1:nz) = alpha * u(istart_w, 1:ny, 1:nz) + beta * u(iend_w, 1:ny, 1:nz)
     v(i_w, 1:ny, 1:nz) = alpha * v(istart_w, 1:ny, 1:nz) + beta * v(iend_w, 1:ny, 1:nz)
     w(i_w, 1:ny, 1:nz) = alpha * w(istart_w, 1:ny, 1:nz) + beta * w(iend_w, 1:ny, 1:nz)

  endif

end do

end subroutine inflow_cond

!**********************************************************************
subroutine project ()
!**********************************************************************
!
! provides u, v, w at 1:nz 
!
use param
use sim_param
use immersedbc
use messages
$if($MPI)
use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWNUP
$endif
implicit none

$if ($DEBUG)
logical, parameter :: DEBUG = .false.
$endif

integer :: jx, jy, jz
integer :: jz_min

real (rprec) :: RHS, tconst

character(*), parameter :: sub_name='project'

! Caching
tconst = tadv1 * dt

do jz = 1, nz - 1
  do jy = 1, ny
    do jx = 1, nx
 
      RHS = -tadv1 * dpdx(jx, jy, jz)
      u(jx, jy, jz) = (u(jx, jy, jz) + dt * (RHS + fx(jx, jy, jz)))
      RHS = -tadv1 * dpdy(jx, jy, jz)
      v(jx, jy, jz) = (v(jx, jy, jz) + dt * (RHS + fy(jx, jy, jz))) 

      !if (DEBUG) then
      !  if ( isnan (u(jx, jy, jz)) ) then
      !    write (*, *) $str($context_doc)
      !    write (*, *) 'nan in u at (jx, jy, jz) = ', jx, jy, jz
      !    stop
      !  end if
      !  if ( isnan (v(jx, jy, jz)) ) then
      !    write (*, *) $str($context_doc)
      !    write (*, *) 'nan in v at (jx, jy, jz) = ', jx, jy, jz
      !    stop
      !  end if
      !end if

    end do
  end do
end do

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
  jz_min = 2
else
  jz_min = 1
end if

do jz = jz_min, nz - 1
  do jy = 1, ny
    do jx = 1, nx

      RHS = -tadv1 * dpdz(jx, jy, jz)
      w(jx, jy, jz) = (w(jx, jy, jz) + dt * (RHS + fz(jx, jy, jz)))

      !if (DEBUG) then
      !  if ( isnan (w(jx, jy, jz)) ) then
      !    write (*, *) 'nan in w at (jx, jy, jz) = ', jx, jy, jz
      !    stop
      !  end if
      !end if

    end do
  end do
end do

if ( inflow ) call inflow_cond ()

!--left this stuff last, so BCs are still enforced, no matter what
!  inflow_cond does

$if ($MPI)

! Exchange ghost node information (since coords overlap)                     
call mpi_sync_real_array( u, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( v, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( w, 0, MPI_SYNC_DOWNUP )  
  
$endif

!--enfore bc at top
if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then

  if (force_top_bot .and. inflow) then
    u(:, :, nz) = face_avg
    v(:, :, nz) = 0._rprec
  else
    ! no-stress top
    u(:,:,nz)=u(:,:,nz-1)
    ! no-stress top
    v(:,:,nz)=v(:,:,nz-1)
  end if

  w(:, :, nz)=0._rprec

end if

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
  ! just a test
  !if (lbc_mom == 'stress free') then
  !  if (force_top_bot) then
  !    u(:, :, 1) = face_avg
  !    v(:, :, 1) = 0._rprec
  !  else
  !    u(:, :, 1) = u(:, :, 2)
  !    v(:, :, 1) = v(:, :, 2)
  !  end if
  !end if

  w(:, :, 1)=0._rprec

end if

end subroutine project
