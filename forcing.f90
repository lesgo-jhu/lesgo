! modified this so that is just calculated the force--it does not do the
! time advancement
subroutine forcing ()
!subroutine forcing(jt)
use types,only:rprec
use param
use sim_param
use immersedbc
$if ($TREES_LS)
  use trees_ls
$endif
$if ($LVLSET)
  use level_set, only : level_set_forcing
$endif
$if ($TURBINES)
  use turbines, only:turbines_forcing
$endif
implicit none

!integer,intent(in)::jt
integer::px,py,lx,ly,lz
integer :: jx,jy,jz,i

real (rprec) :: Rx, Ry, Rz 

! start calculation of body forces (fx,fy,fz)
! 'force' is the mean pressure gradient
if (use_bldg) then
   do i=1,n_bldg
     px=bldg_pts(1,i)
     py=bldg_pts(2,i)
     lx=bldg_pts(3,i)
     ly=bldg_pts(4,i)
     lz=bldg_pts(5,i)
     do jz=1,lz
     do jy=py,py+ly
     do jx=px,px+lx

       ! forces after pressure update
       Rx = -tadv1*dpdx(jx,jy,jz)
       Ry = -tadv1*dpdy(jx,jy,jz)
       Rz = -tadv1*dpdz(jx,jy,jz)

       fx(jx,jy,jz) = ((u_des(jx,jy,jz)-u(jx,jy,jz))/dt - Rx)
       fy(jx,jy,jz) = ((v_des(jx,jy,jz)-v(jx,jy,jz))/dt - Ry)
       fz(jx,jy,jz) = ((w_des(jx,jy,jz)-w(jx,jy,jz))/dt - Rz)

     end do
     end do
     end do
   end do
   ! end calculation of forces
endif

$if ($LVLSET)
  call level_set_forcing ()
$endif

$if ($TREES_LS)
  !--this must come after call to level_set_forcing
  !--in /a posteriori/ test, this adds SGS branch force
  !--in /a priori/ test, this does not modify force
  call trees_ls_calc ()
$endif

$if ($TURBINES)
  call turbines_forcing ()
$endif

if ( inflow .and. use_fringe_forcing ) call inflow_cond ()
    
end subroutine forcing

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function fringe_blend ( x )
use types, only : rp => rprec
implicit none

real (rp) :: fringe_blend

real (rp), intent (in) :: x

real (rp) :: arg

!---------------------------------------------------------------------

if ( x <= 0.0_rp ) then
    fringe_blend = 0.0_rp
else if ( x >= 1.0_rp ) then
    fringe_blend = 1.0_rp
else
    arg = 1.0_rp / ( x - 1 ) + 1.0_rp / x
    fringe_blend = 1.0_rp / ( 1.0_rp + exp ( arg ) )
end if

end function fringe_blend

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine inflow_cond ()
use types, only : rprec
use param, only : face_avg, nx, ny, nz, pi, read_inflow_file,      &
                  sflux_flag, buff_end, buff_len, use_fringe_forcing,  &
                  L_x, dt, dx
use sim_param, only : u, v, w, theta
use immersedbc, only : fx, fy, fz
use io, only : inflow_read
implicit none

integer :: i, i_w
integer :: istart, istart_w
integer :: iend, iend_w

real (rprec) :: factor
real (rprec) :: fringe_blend
real (rprec) :: x1, x2
real (rprec) :: delta_r, delta_f

!---------------------------------------------------------------------

!--these may be out of 1, ..., nx
iend = floor (buff_end * nx + 1._rprec)
istart = floor ((buff_end - buff_len) * nx + 1._rprec)

!--wrapped versions
iend_w = modulo (iend - 1, nx) + 1
istart_w = modulo (istart - 1, nx) + 1

!--read from file
if (read_inflow_file) then  !--read vel inflow @ jx = iend_w from file
  call inflow_read ()  !--this sets u, v, w at (iend_w,:,:)
else
  u(iend_w, :, :) = face_avg
  v(iend_w, :, :) = 0._rprec
  w(iend_w, :, :) = 0._rprec
end if

!--skip istart since we know vel at istart, iend already
do i = istart + 1, iend - 1

  i_w = modulo (i - 1, nx) + 1

  !--fringe function approach: uses forces
  if ( use_fringe_forcing ) then
  
      delta_r = 0.5_rprec * buff_len * L_x
      delta_f = 0.125_rprec * buff_len * L_x
      
      x1 = ( i - istart ) * dx / delta_r
      x2 = ( i - iend ) * dx / delta_f + 1.0_rprec
      
      factor = (1.0_rprec / dt) * ( fringe_blend ( x1 ) - fringe_blend ( x2 ) )

      !--assumes fx is cleared after each time step CHECK THIS
      !--this form is only valid for constant inflow (experimental), need
      !  to alter for general inflows (e.g. turbulent)
      fx(i_w, 1:ny, 1:nz) =  &  !fx(i_w, 1:ny, 1:nz) +             &
                            factor * ( u(iend_w, 1:ny, 1:nz)  &
                                       - u(i_w, 1:ny, 1:nz) )
      fy(i_w, 1:ny, 1:nz) =  &  !fy(i_w, 1:ny, 1:nz) +             &
                            factor * ( v(iend_w, 1:ny, 1:nz)  &
                                       - v(i_w, 1:ny, 1:nz) )
      fz(i_w, 1:ny, 1:nz) =  &  !fz(i_w, 1:ny, 1:nz) +             &
                            factor * ( w(iend_w, 1:ny, 1:nz)  &
                                       - w(i_w, 1:ny, 1:nz) )

  else

      !factor = real ( i - istart, rprec ) / real ( iend - istart, rprec )

      factor = 0.5_rprec * ( 1._rprec - cos (pi * real (i - istart, rprec)  &
                                             / (iend - istart)) )

      !if ( i - istart > (iend - istart) / 2 ) then
      !    factor = 1.0_rprec
      !else
      !    factor = 0.5_rprec *                                        &
      !             ( 1._rprec - cos (2*pi * real (i - istart, rprec)  &
      !                                           / (iend - istart)) )
      !end if

      u(i_w, 1:ny, 1:nz) = u(istart_w, 1:ny, 1:nz) + factor *               &
                            (u(iend_w, 1:ny, 1:nz) - u(istart_w, 1:ny, 1:nz))
      v(i_w, 1:ny, 1:nz) = v(istart_w, 1:ny, 1:nz) + factor *               &
                            (v(iend_w, 1:ny, 1:nz) - v(istart_w, 1:ny, 1:nz))
      w(i_w, 1:ny, 1:nz) = w(istart_w, 1:ny, 1:nz) + factor *               &
                            (w(iend_w, 1:ny, 1:nz) - w(istart_w, 1:ny, 1:nz))

      if (sflux_flag) then
         theta(i_w, 1:ny, 1:nz) = 0._rprec + factor *                   &
                                   (theta(iend_w, 1:ny, 1:nz) - 0._rprec)
      end if

  end if

end do

end subroutine inflow_cond

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--provides u, v, w at 1:nz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine project ()
use param
use sim_param
use immersedbc
implicit none

$if ($DEBUG)
logical, parameter :: DEBUG = .false.
$endif

integer :: jx, jy, jz
integer :: jz_min

real (rprec) :: RHS

!---------------------------------------------------------------------

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

if ( inflow .and. (.not. use_fringe_forcing) ) call inflow_cond ()

!--left this stuff last, so BCs are still enforced, no matter what
!  inflow_cond does
$if ($MPI)
  !--send velocity info down & recv velocity info from above
  call mpi_sendrecv (u(1, 1, 1), ld*ny, MPI_RPREC, down, 1,  &
                     u(1, 1, nz), ld*ny, MPI_RPREC, up, 1,   &
                     comm, status, ierr)
  call mpi_sendrecv (v(1, 1, 1), ld*ny, MPI_RPREC, down, 2,  &
                     v(1, 1, nz), ld*ny, MPI_RPREC, up, 2,   &
                     comm, status, ierr)
  call mpi_sendrecv (w(1, 1, 1), ld*ny, MPI_RPREC, down, 3,  &
                     w(1, 1, nz), ld*ny, MPI_RPREC, up, 3,   &
                     comm, status, ierr)                     
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
