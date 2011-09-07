module level_set
use types, rp => rprec
use param
!use param, only : ld, nx, ny, nz, dx, dy, dz, iBOGUS, BOGUS, VERBOSE,   &
!                  vonK, lbc_mom, USE_MPI, coord, nproc, up, down,       &
!                  comm, ierr, MPI_RPREC, 
use test_filtermodule, only : filter_size
use messages

$if ($DEBUG)
use debug_mod
$endif

use level_set_base
implicit none

$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif

save
private

public :: level_set_forcing, level_set_init, level_set_BC, level_set_Cs
public :: level_set_global_CD
public :: level_set_smooth_vel, level_set_lag_dyn
public :: level_set_Cs_lag_dyn
public :: level_set_vel_err

character (*), parameter :: mod_name = 'level_set'

integer, parameter :: nd = 3

$if ($MPI)
  ! Make sure all values (top and bottom) are less than Nz
  integer, parameter :: nphitop = 3
  integer, parameter :: nphibot = 2
  integer, parameter :: nveltop = 1
  integer, parameter :: nvelbot = 1
  integer, parameter :: ntautop = 3
  integer, parameter :: ntaubot = 2
  integer, parameter :: nFMMtop = 1
  integer, parameter :: nFMMbot = 2
$else
  integer, parameter :: nphitop = 0
  integer, parameter :: nphibot = 0
  integer, parameter :: nveltop = 0
  integer, parameter :: nvelbot = 0
  integer, parameter :: ntautop = 0
  integer, parameter :: ntaubot = 0
  integer, parameter :: nFMMtop = 0
  integer, parameter :: nFMMbot = 0
$endif

!--these are the extra overlap arrays required for BC with MPI
real (rp), dimension (ld, ny, nphitop) :: phitop
real (rp), dimension (ld, ny, nphibot) :: phibot
real (rp), dimension (ld, ny, nveltop) :: utop, vtop, wtop
real (rp), dimension (ld, ny, nvelbot) :: ubot, vbot, wbot
real (rp), dimension (ld, ny, ntautop) :: txxtop, txytop, txztop,  &
                                          tyytop, tyztop, tzztop
real (rp), dimension (ld, ny, ntaubot) :: txxbot, txybot, txzbot,  &
                                          tyybot, tyzbot, tzzbot
!--really only needed for Lagrangian SGS models
real (rp), dimension (ld, ny, nFMMbot) :: FMMbot
real (rp), dimension (ld, ny, nFMMtop) :: FMMtop

$if ($DEBUG)
logical, parameter :: DEBUG = .false.
$endif

real (rp) :: phi_cutoff
real (rp) :: phi_0
!real (rp) :: phi(ld, ny, $lbz:nz)
real (rp) :: norm(nd, ld, ny, $lbz:nz) !--normal vector
                                  !--may want to change so only normals
                                  !  near 0-set are stored
!--experimental: desired velocities for IB method
real (rp) :: udes(ld, ny, $lbz:nz), vdes(ld, ny, $lbz:nz), wdes(ld, ny, $lbz:nz)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains

!**********************************************************************
subroutine level_set_vel_err() 
!**********************************************************************
!
!  This subroutine computes the error in the final velocity field (u^{m+1})
!  with respect to the desired IBM velocity (here zero). The averaged 
!  value is written to file
!
use types, only : rprec
use param, only : nx, ny, nz, total_time
use sim_param, only : u, v, w
$if ($MPI)
use mpi
use param, only : up, down, ierr, MPI_RPREC, status, comm, coord
$endif
implicit none

include 'tecryte.h'

character (*), parameter :: sub_name = mod_name // '.level_set_vel_err'
character(*), parameter :: fname_write = path // 'output/level_set_vel_err.dat'

integer :: i,j,k
integer :: uv_err_navg, w_err_navg
real(rprec) :: u_err, v_err, w_err
$if($MPI)
real(rprec) :: u_err_global, v_err_global, w_err_global
$endif

!  Initialize values
uv_err_navg = 0
w_err_navg = 0 
u_err = 0._rprec
v_err = 0._rprec
w_err = 0._rprec

!  Sum over bottom plane
k=1
do j=1, ny
  do i=1, nx

    if( phi(i,j,k) <= 0._rprec ) then
      uv_err_navg = uv_err_navg + 1
      u_err = u_err + abs( u(i,j,k) )
      v_err = v_err + abs( v(i,j,k) )
    endif

  enddo
enddo

!  Sum over rest of planes
do k=2, nz-1
  do j=1, ny
    do i=1, nx

      if( phi(i,j,k) <= 0._rprec ) then
        uv_err_navg = uv_err_navg + 1
        u_err = u_err + abs( u(i,j,k) )
        v_err = v_err + abs( v(i,j,k) )
      endif


      if( phi(i,j,k) + phi(i,j,k-1) <= 0._rprec ) then
        w_err_navg = w_err_navg + 1
        w_err = w_err + abs( w(i,j,k) )
      endif
      
    enddo
  enddo
enddo

if( uv_err_navg == 0 ) then
  u_err = 0._rprec
  v_err = 0._rprec
else
  u_err = u_err / uv_err_navg
  v_err = v_err / uv_err_navg
endif

if( w_err_navg == 0 ) then
  w_err = 0._rprec
else
  w_err = w_err / w_err_navg
endif

$if( $MPI )

  call mpi_reduce (u_err, u_err_global, 1, MPI_RPREC, MPI_SUM, 0, comm, ierr)
  call mpi_reduce (v_err, v_err_global, 1, MPI_RPREC, MPI_SUM, 0, comm, ierr)
  call mpi_reduce (w_err, w_err_global, 1, MPI_RPREC, MPI_SUM, 0, comm, ierr)

  if( rank == 0 ) then
  
    u_err = u_err_global / nproc
    v_err = v_err_global / nproc
    w_err = w_err_global / nproc

    call write_real_data(fname_write, 'append', 'formatted', 2, &
                         (/ total_time, sqrt( u_err**2 + v_err**2 + w_err**2 ) /))

  endif

$else

call write_real_data(fname_write, 'append', 'formatted', 2, &
                     (/ total_time, sqrt( u_err**2 + v_err**2 + w_err**2 ) /))

$endif





return
end subroutine level_set_vel_err

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--sets Cs2 to epsilon inside solid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine level_set_Cs_lag_dyn ()
use sgsmodule, only : Cs_opt2
implicit none

character (*), parameter :: sub_name = mod_name // '.level_set_Cs_lag_dyn'

real (rp), parameter :: eps = 100._rp * epsilon (0._rp)

integer :: i, j, k
integer :: s

real (rp) :: phix

!---------------------------------------------------------------------
$if ($VERBOSE)
call enter_sub (sub_name)
$endif

do k = 1, nz

  if ( ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) .and.  &
       (k == 1) ) then
    s = 0
  else
    s = 1
  end if
  
  do j = 1, ny
    do i = 1, nx

      phix = 0.5_rp * (phi(i, j, k) + phi(i, j, k - s))

      if (phix < 0._rp) Cs_opt2(i, j, k) = eps

    end do
  end do

end do

$if ($VERBOSE)
call exit_sub (sub_name)
$endif

end subroutine level_set_Cs_lag_dyn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--see trees_BC for reference implementation
!--Sij are on w-nodes
!--F_LM, F_MM on w-nodes
!  * smoothes u, v, w, Sij in solid
!  * zeroes F_LM inside solid
!  * applies neumann condition on F_MM at solid surface
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine level_set_lag_dyn (S11, S12, S13, S22, S23, S33)
use sim_param, only : u, v, w
implicit none

real (rp) :: S11(ld, ny, nz), S12(ld, ny, nz), S13(ld, ny, nz),  &
             S22(ld, ny, nz), S23(ld, ny, nz), S33(ld, ny, nz)

character (*), parameter :: sub_name = mod_name // '.level_set_lag_dyn'

logical, parameter :: lag_dyn_modify_beta = .true.

real (rp) :: phi_c

!---------------------------------------------------------------------
$if ($VERBOSE)
call enter_sub (sub_name)
$endif

!--part 1: smooth variables
phi_c = 0._rp

call smooth (phi_c, $lbz, u)
call smooth (phi_c, $lbz, v)
call smooth (phi_c, $lbz, w, 'w')

call smooth (phi_c, 1, S11, 'w')
call smooth (phi_c, 1, S12, 'w')
call smooth (phi_c, 1, S13, 'w')
call smooth (phi_c, 1, S22, 'w')
call smooth (phi_c, 1, S23, 'w')
call smooth (phi_c, 1, S33, 'w')

!--part 2: zero F_LM
call zero_F_LM ()

!--part 3: neumann condition on F_MM
call neumann_F_MM ()

!--part 4 (optional): modify beta
if (lag_dyn_modify_beta) call modify_beta ()

$if ($VERBOSE)
call exit_sub (sub_name)
$endif

end subroutine level_set_lag_dyn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine modify_beta ()
use sgsmodule, only : beta
implicit none

character (*), parameter :: sub_name = mod_name // '.modify_beta'

real (rp), parameter :: c1 = 0.65_rp
real (rp), parameter :: c2 = 0.7_rp

integer :: i, j, k
integer :: s

real (rp) :: delta
real (rp) :: phix
real (rp) :: dmin
real (rp) :: z

!---------------------------------------------------------------------
$if ($VERBOSE)
call enter_sub (sub_name)
$endif

delta = filter_size * (dx * dy * dz)**(1._rp / 3._rp)

do k = 1, nz

  if ( ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) .and.  &
       (k == 1) ) then
    s = 0
  else
    s = 1
  end if
  
  do j = 1, ny
    do i = 1, nx

      phix = 0.5_rp * (phi(i, j, k) + phi(i, j, k-s))

      if (phix > 0._rp) then

        if (lbc_mom == 'stress free') then
          dmin = phix
        else
          z = (k - 1) * dz
          dmin = min (z, phix)
        end if
        
        beta(i, j, k) = 1._rp - c1 * exp (-c2 * dmin / delta)

      end if
      
    end do
  end do

end do

$if ($VERBOSE)
call exit_sub (sub_name)
$endif

end subroutine modify_beta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--requires F_MM to be valid on 1:nz  (interp_F_MM does)
!--only applies neumann condition to 1:nz-1 though, so additional
!  sync between nz and 1' is needed
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine neumann_F_MM ()
use sgsmodule, only : F_MM
implicit none

character (*), parameter :: sub_name = mod_name // '.neumann_F_MM'

integer, parameter :: tag = 1200

integer :: i, j, k
integer :: s

real (rp) :: phix
real (rp) :: phi1, phi2  !--take action when phi1 < phix < phi2
real (rp) :: dphi
real (rp) :: F_MM_xp
real (rp) :: x(nd), xp(nd)
real (rp) :: n_hat(nd)

!---------------------------------------------------------------------
$if ($VERBOSE)
call enter_sub (sub_name)
$endif

phi1 = -sqrt (dx**2 + dy**2 + dz**2)
phi2 = 0._rp

!--for now dphi does not depend on normal
dphi = sqrt (dx**2 + dy**2 + dz**2)

$if ($MPI)
  call mpi_sync_F_MM ()
$endif

do k = 1, nz - 1

  if ( ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) .and.  &
       (k == 1) ) then
    s = 0  !--u-nodes
  else
    s = 1  !--w-nodes
  end if

  do j = 1, ny
    do i = 1, nx

      phix = 0.5_rp * (phi(i, j, k) + phi(i, j, k - s)) 
   
      if ((phi1 < phix) .and. (phix < phi2)) then

        x = (/ (i - 1) * dx, (j - 1) * dy, (k - 0.5_rp * (1+s)) * dz /)
                                            !--select u/w node

        n_hat = 0.5_rp * (norm(:, i, j, k) + norm(:, i, j, k - s))

        !--dphi depends on normal?  should be > phi2 - phi1
        xp = x + dphi * n_hat

        call interp_scal (1, F_MM, nFMMbot, FMMbot, nFMMtop, FMMtop,  &
                          xp, F_MM_xp, 'w')
        
        F_MM(i, j, k) = F_MM_xp
        
      end if
      
    end do
  end do
end do

$if ($MPI)
  !--make F_MM valid at 1:nz (as in core code) by syncing nz to 1'
  !--not sure this is crucial
  call mpi_sendrecv (F_MM(1, 1, 1), ld*ny, MPI_RPREC, down, tag+1,  &
                     F_MM(1, 1, nz), ld*ny, MPI_RPREC, up, tag+1,   &
                     comm, status, ierr)
$endif

$if ($VERBOSE)
call exit_sub (sub_name)
$endif

$if ($MPI)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  contains
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_sync_F_MM ()
  implicit none

  integer, parameter :: tag = 1100

  integer :: datasize
  integer :: kstart

  !-------------------------------------------------------------------
  
  datasize = ld * ny * nFMMtop
  kstart = 2
  call mpi_sendrecv (F_MM(1, 1, kstart), datasize, MPI_RPREC, down, tag+1,  &
                     FMMtop(1, 1, 1), datasize, MPI_RPREC, up, tag+1,       &
                     comm, status, ierr)

  datasize = ld * ny * nFMMbot
  kstart = nz - nFMMbot  !--F_MM is has lbz = 1
  call mpi_sendrecv (F_MM(1, 1, kstart), datasize, MPI_RPREC, up, tag+4,  &
                     FMMbot(1, 1, 1), datasize, MPI_RPREC, down, tag+4,   &
                     comm, status, ierr)

  end subroutine mpi_sync_F_MM
$endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine neumann_F_MM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--zeros F_LM for phi < phi_F_LM
!--only operates on 1:nz-1, MAY NEED 1:nz for complete consistency!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine zero_F_LM ()
use sgsmodule, only : F_LM
implicit none

character (*), parameter :: sub_name = mod_name // '.zero_F_LM'

integer, parameter :: tag = 1300

integer :: i, j, k
integer :: s

real (rp) :: phix
real (rp) :: phi_F_LM

!---------------------------------------------------------------------
$if ($VERBOSE)
call enter_sub (sub_name)
$endif

!phi_F_LM = 0._rp
phi_F_LM = filter_size * dx  !--experimental

do k = 1, nz - 1

  if ( ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) .and.  &
       (k == 1) ) then
    s = 0  !--use phi at u-node
  else
    s = 1
  end if
  
  do j = 1, ny
    do i = 1, nx

      
      phix = 0.5_rp * (phi(i, j, k) + phi(i, j, k-s))
      
      if (phix < phi_F_LM) F_LM(i, j, k) = 0._rp

    end do
  end do

end do

$if ($MPI)
  !--make F_LM valid at 1:nz (as in core code) by syncing nz to 1'
  !--not sure this is crucial
  !--probably can change above loop to 1:nz, since phi is valid there
  call mpi_sendrecv (F_LM(1, 1, 1), ld*ny, MPI_RPREC, down, tag+1,  &
                     F_LM(1, 1, nz), ld*ny, MPI_RPREC, up, tag+1,   &
                     comm, status, ierr)
$endif

$if ($VERBOSE)
call exit_sub (sub_name)
$endif

end subroutine zero_F_LM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine modify_dutdn (i, j, k, tau, phix, x_hat, y_hat, z_hat, node)
use sim_param, only : dudx, dudy, dudz,  &
                      dvdx, dvdy, dvdz,  &
                      dwdx, dwdy, dwdz
implicit none

integer, intent (in) :: i, j, k

real (rp), intent (in) :: tau
real (rp), intent (in) :: phix
real (rp), intent (in) :: x_hat(nd), y_hat(nd), z_hat(nd)

character (*), intent (in), optional :: node  !--'u' or 'w'

character (*), parameter :: sub_name = mod_name // '.modify_dutdn'

logical :: unode

real (rp) :: A(nd, nd)
real (rp) :: G(nd, nd), Gp(nd, nd)
real (rp) :: grad
real (rp) :: phi_min

!---------------------------------------------------------------------
$if ($VERBOSE)
call enter_sub (sub_name)
$endif

A(1, :) = x_hat
A(2, :) = y_hat
A(3, :) = z_hat

if (present (node)) then

  select case (node)
    case ('u'); unode = .true.
    case ('w'); unode = .false.
    case default; call error (sub_name, 'invalid node = ' // node)
  end select

else  !--assume u-node

  unode = .true.
  
end if

if (unode) then
  G(:, 1) = (/ dudx(i, j, k), dvdx(i, j, k),                 &
               0.5_rp * (dwdx(i, j, k) + dwdx(i, j, k + 1)) /)
  G(:, 2) = (/ dudy(i, j, k), dvdy(i, j, k),                 &
               0.5_rp * (dwdy(i, j, k) + dwdy(i, j, k + 1)) /)
  G(:, 3) = (/ 0.5_rp * (dudz(i, j, k) + dudz(i, j, k + 1)),  &
               0.5_rp * (dvdz(i, j, k) + dvdz(i, j, k + 1)),  &
               dwdz(i, j, k) /)
else  !--w-node
  G(:, 1) = (/ 0.5_rp * (dudx(i, j, k) + dudx(i, j, k - 1)),  &
               0.5_rp * (dvdx(i, j, k) + dvdx(i, j, k - 1)),  &
               dwdx(i, j, k) /)
  G(:, 2) = (/ 0.5_rp * (dudy(i, j, k) + dudy(i, j, k - 1)),  &
               0.5_rp * (dvdy(i, j, k) + dvdy(i, j, k - 1)),  &
               dwdy(i, j, k) /)
  G(:, 3) = (/ dudz(i, j, k), dvdz(i, j, k),                 &
               0.5_rp * (dwdz(i, j, k) + dwdz(i, j, k - 1)) /)
end if

!--calculate Gp: gradient in primed (rotated) coordinates
Gp = matmul ( A, matmul (G, transpose (A)) )

!--set the value of dut/dn to be consistent w/log-law
if (phi_cutoff_is_set) then
  phi_min = (0.1_rp * phi_cutoff)  !--experimental
else
  call error (sub_name, 'trying to use unset phi_cutoff')
end if

grad = sqrt (abs (tau)) / (vonK * (max (phi_min, phix) + z0))

Gp(1, 3) = grad

!--rotate back to x, y, z frame
G = matmul ( transpose (A), matmul (Gp, A))

!--still not sure if this is right: only re-assign u-node components
!  on u-nodes, same with w
if (unode) then
  dudx(i, j, k) = G(1, 1)
  dvdx(i, j, k) = G(2, 1)
  dudy(i, j, k) = G(1, 2)
  dvdy(i, j, k) = G(2, 2)
  dwdz(i, j, k) = G(3, 3)
else
  dudz(i, j, k) = G(1, 3)
  dvdz(i, j, k) = G(2, 3)
  dwdx(i, j, k) = G(3, 1)
  dwdy(i, j, k) = G(3, 2)
end if

$if ($VERBOSE)
call exit_sub (sub_name)
$endif

end subroutine modify_dutdn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--this is supposed to replace extrap_tau
!--this avoids using fit3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine extrap_tau_simple ()
use param, only : jt
use sim_param, only : u, v, w, txx, txy, txz, tyy, tyz, tzz
implicit none

character (*), parameter :: sub_name = mod_name // '.extrap_tau_simple'
character (*), parameter :: fprefix = 'output/extrap_tau_simple.'
character (*), parameter :: fmta3r = '(a,3(es12.5,1x))'
character (*), parameter :: fmta3i = '(a,3(i0,1x))'

integer, parameter :: noutput = 200
integer, parameter :: lun = 1

!logical, parameter :: DEBUG = .false.
logical, parameter :: use_output = .false.

real (rp), parameter :: eps = 100._rp * epsilon (0._rp)

character (128) :: fname

integer :: i, j, k
!--experiment to reduce time spent on looping over fluid pts
integer, save :: imn = 1, imx = nx, jmn = 1, jmx = ny
integer :: imn_used, imx_used, jmn_used, jmx_used
integer :: nbad
integer :: kmn

logical :: output
logical :: exst, opn

real (rp) :: phi_c, phi_x
real (rp) :: dphi0, dphi
real (rp) :: phi1, phi2
real (rp) :: txx1, txy1, txz1, tyy1, tyz1, tzz1
real (rp) :: txx2, txy2, txz2, tyy2, tyz2, tzz2
real (rp) :: wgt
real (rp) :: x(nd), x1(nd), x2(nd)
real (rp) :: n_hat(nd)

!---------------------------------------------------------------------
$if ($VERBOSE)
call enter_sub (sub_name)
$endif

if (use_output) then

  if (modulo (jt, noutput) == 0) then

    output = .true.

    write (fname, '(a,i6.6,a)') trim (fprefix), jt, '.dat'

    inquire (unit=lun, exist=exst, opened=opn)

    if (exst .and. (.not. opn)) then
      open (lun, file=fname)
    else
      call error (sub_name, 'problem opening file')
    end if

  else
    output = .false.
  end if

else

  output = .false.

end if

if (phi_cutoff_is_set) then
  phi_c = phi_cutoff
else
  call error (sub_name, 'trying to use uninitialized phi_cutoff')
end if

if (.not. phi_0_is_set) then
  call error (sub_name, 'trying to use uninitialized phi_0')
end if

dphi0 = filter_size * dx  !--experiment

if (output) then
  write (lun, *) 'phi_c = ', phi_c
  write (lun, *) 'phi_0 = ', phi_0
  write (lun, *) 'dphi = ', dphi
end if

!  initial values for _used variables
!  setting to non-previously used values
imn_used = nx
imx_used = 1
jmn_used = ny
jmx_used = 1

nbad = 0

if (output) write (lun, *) 'u-node pass'

!--u-node pass
do k = 1, nz - 1
  do j = jmn, jmx
    do i = imn, imx

      phi_x = phi(i, j, k)

      if ((-phi_c <= phi_x) .and. (phi_x < phi_0)) then

        ! Find bounding box for all "banded" points"
        imn_used = min (imn_used, i)
        imx_used = max (imx_used, i)
        jmn_used = min (jmn_used, j)
        jmx_used = max (jmx_used, j)

        x = (/ (i - 1) * dx, (j - 1) * dy, (k - 0.5_rp) * dz /)

        $if ($DEBUG)
        if (DEBUG) then
          call mesg (sub_name, '(i, j, k) =', (/ i, j, k /))
          call mesg (sub_name, 'x =', x)
        end if
        $endif
        
        n_hat = norm(:, i, j, k)

        dphi = dphi0
        x1 = x + dphi * n_hat  !--this pt is supposed to be in fluid
                               !--this means dphi >= phi_c

        !--check phi(x1) >= 0
        call interp_phi (x1, phi1)

        if (phi1 < phi_0) then
          !--try a larger dphi--experimental
          dphi = 1.5_rp * dphi
          x1 = x + dphi * n_hat
          call interp_phi (x1, phi1)

          if (phi1 < phi_0) then
            !--probably here because of the way the normal (phi)
            !  inside the trees are calculated: its not a true
            !  signed distance function, due to the complex internal
            !  geometry
            !--not sure how to set tau here: for now set to 0._rp
            nbad = nbad + 1

            !call mesg (sub_name, 'bad point (u): ', (/ i, j, k /))

            txx(i, j, k) = 0._rp
            txy(i, j, k) = 0._rp
            txz(i, j, k) = 0._rp
            tyy(i, j, k) = 0._rp
            tyz(i, j, k) = 0._rp
            tzz(i, j, k) = 0._rp

            cycle  !--onto next point

            !call mesg (sub_name, 'at u-node (i, j, k) =', (/ i, j, k /))
            !call mesg (sub_name, 'phi_x =', phi_x)
            !call mesg (sub_name, 'x =', x)
            !call mesg (sub_name, 'n_hat =', n_hat)
            !call mesg (sub_name, 'dphi =', dphi)
            !call mesg (sub_name, 'x1 =', x1)
            !call mesg (sub_name, 'phi_0 =', phi_0)
            !call mesg (sub_name, 'phi1 =', phi1)
            !call error (sub_name, 'phi1 < phi_0')

          end if
        end if

        !if (phi1 < 0._rp) call error (sub_name, 'phi1 =', phi1)

        x2 = x1 + dphi * n_hat  !--this is supposed to be further into fluid
                                !  than x1
                                !--problem for highly wrinkled surface

        !--check phi(x2) >= 0
        call interp_phi (x2, phi2)
        
        !--calculate tij1, for use in extrapolation
        call interp_tij_u (x1, txx1, txy1, tyy1, tzz1)
        !call interp_scal (txx, x1, txx1)
        !call interp_scal (txy, x1, txy1)
        !call interp_scal (tyy, x1, tyy1)
        !call interp_scal (tzz, x1, tzz1)

        !--calculate tij2, for use in extrapolation
        call interp_tij_u (x2, txx2, txy2, tyy2, tzz2)
        !call interp_scal (txx, x2, txx2)
        !call interp_scal (txy, x2, txy2)
        !call interp_scal (tyy, x2, tyy2)
        !call interp_scal (tzz, x2, tzz2)

        !--now extrapolate
        !--simple linear extrap, for now
        !wgt = (phi2 - phi1) / (phi2 - phi_x)
        !--since the two dphi steps are equal, the wgt is 0.5
        wgt = 0.5_rp

        !--the way we define the phis right now, expect wgt ~= 0.5
        !--its a problem is wgt is smaller than say 0.25
        !if (wgt < 0.25_rp) then
        !  write (msg, *) 'extrapolation weight too small (u-node)', n_l,  &
        !                 'coord = ', coord, n_l,                          &
        !                 'wgt = ', wgt, n_l,                              &
        !                 'i, j, k = ', i, j, k, n_l,                      &
        !                 'x = ', x, n_l,                                  &
        !                 'x1 = ', x1, n_l,                                &
        !                 'x2 = ', x2, n_l,                                &
        !                 'phi = ', phi_x, n_l,                            &
        !                 'phi1 = ', phi1, n_l,                            &
        !                 'phi2 = ', phi2, n_l,                            &
        !                 'n_hat = ', n_hat
        !  call error (sub_name, msg)
        !end if

        txx(i, j, k) = (txx1 - (1._rp - wgt) * txx2) / wgt
        txy(i, j, k) = (txy1 - (1._rp - wgt) * txy2) / wgt
        tyy(i, j, k) = (tyy1 - (1._rp - wgt) * tyy2) / wgt
        tzz(i, j, k) = (tzz1 - (1._rp - wgt) * tzz2) / wgt

        if (output) then
          write (lun, fmta3i) 'i, j, k = ', i, j, k
          write (lun, *) 'phi_x = ', phi_x
          write (lun, fmta3r) 'x = ', x
          write (lun, fmta3r) 'n_hat = ', n_hat
          write (lun, fmta3r) 'x1 = ', x1
          write (lun, *) 'phi1 = ', phi1
          write (lun, fmta3r) 'x2 = ', x2
          write (lun, *) 'phi2 = ', phi2
          write (lun, *) 'txx1 = ', txx1
          write (lun, *) 'txy1 = ', txy1
          write (lun, *) 'tyy1 = ', tyy1
          write (lun, *) 'tzz1 = ', tzz1
          write (lun, *) 'txx2 = ', txx2
          write (lun, *) 'txy2 = ', txy2
          write (lun, *) 'tyy2 = ', tyy2
          write (lun, *) 'tzz2 = ', tzz2
          write (lun, *) 'wgt = ', wgt
          write (lun, *) 'txx = ', txx(i, j, k)
          write (lun, *) 'txy = ', txy(i, j, k)
          write (lun, *) 'tyy = ', tyy(i, j, k)
          write (lun, *) 'tzz = ', tzz(i, j, k)
        end if

      end if

    end do
  end do
end do

$if ($DEBUG)
if (DEBUG) call mesg (sub_name, '# u bad points =', nbad)
$endif
nbad = 0

if (output) write (lun, *) 'w-node pass'

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
  kmn = 2
else
  kmn = 1
end if

!--w-node pass
do k = kmn, nz - 1
  do j = jmn, jmx
    do i = imn, imx

      phi_x = 0.5_rp * (phi(i, j, k) + phi(i, j, k - 1))
                                       !--MPI: requires phi(k=0)

      if ((-phi_c <= phi_x) .and. (phi_x < phi_0)) then

        imn_used = min (imn_used, i)
        imx_used = max (imx_used, i)
        jmn_used = min (jmn_used, j)
        jmx_used = max (jmx_used, j)

        x = (/ (i - 1) * dx, (j - 1) * dy, (k - 1) * dz /)

        n_hat = 0.5_rp * (norm(:, i, j, k) + norm(:, i, j, k - 1))
                                             !--MPI: requires phi(k=0)

        if (mag (n_hat) > eps) then
          n_hat = n_hat / mag (n_hat)
        else  !--normal cancelled?
          n_hat = norm(:, i, j, k)
        end if

        dphi = dphi0
        x1 = x + dphi * n_hat

        !--check phi(x1) >= 0
        call interp_phi (x1, phi1)

        if (phi1 < phi_0) then
          !--try a larger dphi--experimental
          dphi = 1.5_rp * dphi
          x1 = x + dphi * n_hat
          call interp_phi (x1, phi1)

          if (phi1 < phi_0) then
          
            nbad = nbad + 1

            !call mesg (sub_name, 'bad point (w): ', (/ i, j, k /))

            txz(i, j, k) = 0._rp
            tyz(i, j, k) = 0._rp

            cycle

            !call mesg (sub_name, 'at w-node (i, j, k) =', (/ i, j, k /))
            !call mesg (sub_name, 'phi_x =', phi_x)
            !call mesg (sub_name, 'x =', x)
            !call mesg (sub_name, 'n_hat =', n_hat)
            !call mesg (sub_name, 'dphi =', dphi)
            !call mesg (sub_name, 'x1 =', x1)
            !call mesg (sub_name, 'phi_0 =', phi_0)
            !call mesg (sub_name, 'phi1 =', phi1)
            !call error (sub_name, 'phi1 < phi_0')

          end if
        end if

        !if (phi1 < 0._rp) call error (sub_name, 'phi1 =', phi1)

        x2 = x1 + dphi * n_hat  !--this is supposed to be further into fluid
                                !  than x1

        !--check phi(x2) >= 0
        call interp_phi (x2, phi2)

        !--calculate tij1, for use in extrapolation
        call interp_tij_w (x1, txz1, tyz1)
        !call interp_scal (txz, x1, txz1, 'w')
        !call interp_scal (tyz, x1, tyz1, 'w')

        !--calculate tij2, for use in extrapolation
        call interp_tij_w (x2, txz2, tyz2)
        !call interp_scal (txz, x2, txz2, 'w')
        !call interp_scal (tyz, x2, tyz2, 'w')

        !--now extrapolate
        !--simple linear extrap, for now
        !wgt = (phi2 - phi1) / (phi2 - phi_x)
        wgt = 0.5_rp

        !--the way we define the phis right now, expect wgt ~= 0.5
        !--its a problem is wgt is smaller than say 0.25
        !if (wgt < 0.25_rp) then
        !  write (msg, *) 'extrapolation weight too small (w-node)', n_l,  &
        !                 'wgt = ', wgt, n_l,                              &
        !                 'i, j, k = ', i, j, k, n_l,                      &
        !                 'x = ', x, n_l,                                  &
        !                 'x1 = ', x1, n_l,                                &
        !                 'x2 = ', x2, n_l,                                &
        !                 'phi = ', phi_x, n_l,                            &
        !                 'phi1 = ', phi1, n_l,                            &
        !                 'phi2 = ', phi2, n_l,                            &
        !                 'n_hat = ', n_hat
        !  call error (sub_name, msg)
        !end if
          
        txz(i, j, k) = (txz1 - (1._rp - wgt) * txz2) / wgt
        tyz(i, j, k) = (tyz1 - (1._rp - wgt) * tyz2) / wgt

        if (output) then
          write (lun, fmta3i) 'i, j, k = ', i, j, k
          write (lun, *) 'phi_x = ', phi_x
          write (lun, fmta3r) 'x = ', x
          write (lun, fmta3r) 'n_hat = ', n_hat
          write (lun, fmta3r) 'x1 = ', x1
          write (lun, *) 'phi1 = ', phi1
          write (lun, fmta3r) 'x2 = ', x2
          write (lun, *) 'phi2 = ', phi2
          write (lun, *) 'txz1 = ', txz1
          write (lun, *) 'tyz1 = ', tyz1
          write (lun, *) 'txz2 = ', txz2
          write (lun, *) 'tyz2 = ', tyz2
          write (lun, *) 'wgt = ', wgt
          write (lun, *) 'txz = ', txz(i, j, k)
          write (lun, *) 'tyz = ', tyz(i, j, k)
        end if

      end if

    end do
  end do
end do

$if ($DEBUG)
if (DEBUG) call mesg (sub_name, '# w bad points =', nbad)
$endif

imn = imn_used
imx = imx_used
jmn = jmn_used
jmx = jmx_used

if (output) then
  write (lun, *) ' '
  close (lun)
end if

$if ($VERBOSE)
call exit_sub (sub_name)
$endif

end subroutine extrap_tau_simple

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--extrapolate fluid tau through log-law at wall into solid
!--this is NOT a drop-in replacement for extrap_tau, instead it
!  replaces interp_tau AND extrap_tau
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine extrap_tau_log ()
use sim_param, only : u, v, w, txx, txy, txz, tyy, tyz, tzz
implicit none

character (*), parameter :: sub_name = mod_name // '.extrap_tau_log'

real (rp), parameter :: eps = 100._rp * epsilon (1._rp)

integer :: i, j, k

real (rp) :: phi_a, phi_b, phi_x
real (rp) :: dphi
real (rp) :: phi1, phi2
real (rp) :: tau_w
real (rp) :: txx1, txy1, txz1, tyy1, tyz1, tzz1
real (rp) :: txx2, txy2, txz2, tyy2, tyz2, tzz2
real (rp) :: txx_im, txy_im, txz_im, tyy_im, tyz_im, tzz_im
real (rp) :: txx_w, txy_w, txz_w, tyy_w, tyz_w, tzz_w
real (rp) :: wgt, wgt_im
real (rp) :: v1n
real (rp) :: x(nd), x1(nd), x2(nd)
real (rp) :: n_hat(nd)
real (rp) :: v1(nd), v1t(nd)
real (rp) :: x_hat(nd), y_hat(nd), z_hat(nd)

!---------------------------------------------------------------------
$if ($VERBOSE)
call enter_sub (sub_name)
$endif

!--assumes phi_cutoff is set. CAREFUL
!--set phi_a
phi_a = -phi_cutoff

!--set phi_b
phi_b = 0._rp

dphi = phi_cutoff

!--u-node pass
do k = 1, nz - 1
  do j = 1, ny
    do i = 1, nx

      phi_x = phi(i, j, k)

      if ((phi_a < phi_x) .and. (phi_x < phi_b)) then

        x = (/ (i - 1) * dx, (j - 1) * dy, (k - 0.5_rp) * dz /)

        n_hat = norm(:, i, j, k)

        x1 = x + (dphi - phi_a) * n_hat

        call interp_vel (x1, v1)

        v1n = dot_product (v1, n_hat)  !--scalar
        v1t = v1 - v1n * n_hat         !--vector

        if (mag (v1t) < eps) then

          txx(i, j, k) = 0._rp
          txy(i, j, k) = 0._rp
          tyy(i, j, k) = 0._rp
          tzz(i, j, k) = 0._rp
       
        else
        
          !--local coordinate system
          x_hat = v1t / mag (v1t)
          y_hat = cross_product (n_hat, x_hat)
          z_hat = n_hat

          call interp_phi (x1, phi1)

          tau_w = -(mag (v1t) * vonk / log (phi1 / z0))**2

          !--next step is to rotate
          !--special case: only nonzero t_{i'j'} is t_{1'3'} & t_{3'1'}
          txx_w = (x_hat(1) * z_hat(1) + z_hat(1) * x_hat(1)) * tau_w
          txy_w = (x_hat(1) * z_hat(2) + z_hat(1) * x_hat(2)) * tau_w
          tyy_w = (x_hat(2) * z_hat(2) + z_hat(2) * x_hat(2)) * tau_w
          tzz_w = (x_hat(3) * z_hat(3) + z_hat(3) * x_hat(3)) * tau_w

          !--calculate tij1, for use in extrapolation
          call interp_tij_u (x1, txx1, txy1, tyy1, tzz1)
          !call interp_scal (txx, x1, txx1)
          !call interp_scal (txy, x1, txy1)
          !call interp_scal (tyy, x1, tyy1)
          !call interp_scal (tzz, x1, tzz1)

          !--now extrapolate
          !--simple linear extrap is OK for now
          wgt = abs (phi1) / (abs (phi_x) + abs (phi1))

          if (wgt >= 0.5_rp) then
          
            txx(i, j, k) = (txx_w - (1._rp - wgt) * txx1) / wgt
            txy(i, j, k) = (txy_w - (1._rp - wgt) * txy1) / wgt
            tyy(i, j, k) = (tyy_w - (1._rp - wgt) * tyy1) / wgt
            tzz(i, j, k) = (tzz_w - (1._rp - wgt) * tzz1) / wgt

          else  !--image pt method

            x2 = x1 + dphi * n_hat
            call interp_phi (x2, phi2)

            call interp_tij_u (x2, txx2, txy2, tyy2, tzz2)
            !call interp_scal (txx, x2, txx2)
            !call interp_scal (txy, x2, txy2)
            !call interp_scal (tyy, x2, tyy2)
            !call interp_scal (tzz, x2, tzz2)

            wgt_im = (abs (phi2) - abs (phi_x)) / (abs (phi2) - abs (phi1))

            txx_im = wgt_im * txx1 + (1._rp - wgt_im) * txx2
            txy_im = wgt_im * txy1 + (1._rp - wgt_im) * txy2
            tyy_im = wgt_im * tyy1 + (1._rp - wgt_im) * tyy2
            tzz_im = wgt_im * tzz1 + (1._rp - wgt_im) * tzz2

            txx(i, j, k) = 2._rp * txx_w - txx_im
            txy(i, j, k) = 2._rp * txy_w - txy_im
            tyy(i, j, k) = 2._rp * tyy_w - tyy_im
            tzz(i, j, k) = 2._rp * tzz_w - tzz_im

          end if

        end if

        $if ($IFORT || $IFC)

          if (isnan (txx(i, j, k))) then
            call mesg (sub_name, 'phi_x =', phi_x)
            call mesg (sub_name, 'x =', x)
            call mesg (sub_name, 'n_hat =', n_hat)
            call mesg (sub_name, 'x1 =', x1)
            call mesg (sub_name, 'v1 =', v1)
            call mesg (sub_name, 'x_hat =', x_hat)
            call mesg (sub_name, 'y_hat =', y_hat)
            call mesg (sub_name, 'z_hat =', z_hat)
            call mesg (sub_name, 'phi1 =', phi1)
            call mesg (sub_name, 'v1n =', v1n)
            call mesg (sub_name, 'v1t =', v1t)
            call mesg (sub_name, 'tau_w =', tau_w)
            call mesg (sub_name, 'wgt =', wgt)
            call mesg (sub_name, 'x2 =', x2)
            call mesg (sub_name, 'phi2 =', phi2)
            call mesg (sub_name, 'txx2 =', txx2)
            call mesg (sub_name, 'wgt_im =', wgt_im)
            call mesg (sub_name, 'txx_im =', txx_im)
            call error (sub_name, 'Nan in txx at (i, j, k) =', (/ i, j, k /))
          end if

          if (isnan (txy(i, j, k))) then
            call error (sub_name, 'NaN in txy at (i, j, k) =', (/ i, j, k /))
          end if

          if (isnan (tyy(i, j, k))) then
            call error (sub_name, 'NaN in tyy at (i, j, k) =', (/ i, j, k /))
          end if

          if (isnan (tzz(i, j, k))) then
            call error (sub_name, 'NaN in tzz at (i, j, k) =', (/ i, j, k /))
          end if

        $endif

      end if

    end do
  end do
end do

!--w-node pass
do k = 2, nz
  do j = 1, ny
    do i = 1, nx

      phi_x = 0.5_rp * (phi(i, j, k) + phi(i, j, k - 1))

      if ((phi_a < phi_x) .and. (phi_x < phi_b)) then

        x = (/ (i - 1) * dx, (j - 1) * dy, (k - 1) * dz /)

        n_hat = 0.5_rp * (norm(:, i, j, k) + norm(:, i, j, k - 1))
        n_hat = n_hat / mag (n_hat)

        x1 = x + (dphi - phi_a) * n_hat

        call interp_vel (x1, v1)

        v1n = dot_product (v1, n_hat)  !--scalar
        v1t = v1 - v1n * n_hat         !--vector

        if (mag (v1t) < eps) then

          txz(i, j, k) = 0._rp
          tyz(i, j, k) = 0._rp
       
        else
        
          !--local coordinate system
          x_hat = v1t / mag (v1t)
          y_hat = cross_product (n_hat, x_hat)
          z_hat = n_hat

          call interp_phi (x1, phi1)

          tau_w = -(mag (v1t) * vonk / log (phi1 / z0))**2

          !--next step is to rotate
          !--special case: only nonzero t_{i'j'} is t_{1'3'} & t_{3'1'}
          txz_w = (x_hat(1) * z_hat(3) + z_hat(1) * x_hat(3)) * tau_w
          tyz_w = (x_hat(2) * z_hat(3) + z_hat(2) * x_hat(3)) * tau_w

          !--calculate tij1, for use in extrapolation
          call interp_tij_w (x1, txz1, tyz1)
          !call interp_scal (txz, x1, txz1, 'w')
          !call interp_scal (tyz, x1, tyz1, 'w')

          !--now extrapolate
          !--simple linear extrap is OK for now
          wgt = abs (phi1) / (abs (phi_x) + abs (phi1))

          if (wgt >= 0.5_rp) then
          
            txz(i, j, k) = (txz_w - (1._rp - wgt) * txz1) / wgt
            tyz(i, j, k) = (tyz_w - (1._rp - wgt) * tyz1) / wgt

          else  !--image pt method

            x2 = x1 + dphi * n_hat
            call interp_phi (x2, phi2)

            call interp_tij_w (x2, txz2, tyz2)
            !call interp_scal (txz, x2, txz2, 'w')
            !call interp_scal (tyz, x2, tyz2, 'w')

            wgt_im = (abs (phi2) - abs (phi_x)) / (abs (phi2) - abs (phi1))

            txz_im = wgt_im * txz1 + (1._rp - wgt_im) * txz2
            tyz_im = wgt_im * tyz1 + (1._rp - wgt_im) * tyz2

            txz(i, j, k) = 2._rp * txz_w - txz_im
            tyz(i, j, k) = 2._rp * tyz_w - tyz_im

          end if

        end if

      end if

    end do
  end do
end do

$if ($VERBOSE)
call exit_sub (sub_name)
$endif

end subroutine extrap_tau_log

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine enforce_un ()
use param, only : jt  !--plus stuff above
use sim_param, only : u, v, w
implicit none

character (*), parameter :: sub_name = mod_name // '.enforce_un'

!logical, parameter :: DEBUG = .false.

integer :: i, j, k

real (rp) :: dphi
real (rp) :: phi0, phi1, phi2, phic
real (rp) :: phi_x2
real (rp) :: x1(nd), x2(nd)
real (rp) :: n_hat(nd)
real (rp) :: vel1(nd), vel2(nd)
real (rp) :: vel1_n, vel2_n

!---------------------------------------------------------------------
$if ($VERBOSE)
call enter_sub (sub_name)
$endif

!--need to experiment with these values
if (phi_0_is_set) then
  phi0 = phi_0
else
  call error (sub_name, 'trying to use uninitialized phi_0')
end if

if (phi_cutoff_is_set) then
  phic = phi_cutoff
else
  call error (sub_name, 'trying to use uninitialized phi_cutoff')
end if

!phic = filter_size * dx
!phic = filter_size * sqrt (dx**2 + dy**2 + dz**2)

!dphi = filter_size * dx
dphi = filter_size * sqrt (dx**2 + dy**2 + dz**2)

udes = huge (1._rp)
vdes = huge (1._rp)
wdes = huge (1._rp)

!--u-node pass
do k = 1, nz - 1
  do j = 1, ny
    do i = 1, nx

      phi1 = phi(i, j, k)

      if ((phi0 < phi1) .and. (phi1 < phic)) then
 
        x1 = (/ (i-1) * dx, (j-1) * dy, (k-0.5_rp) * dz /)

        n_hat = norm(:, i, j, k)

        x2 = x1 + dphi * n_hat

        phi2 = phi1 + dphi  !--assume this is approx. correct

        call interp_phi (x2, phi_x2)

        !if (phi_x2 > phi0) then

          call interp_vel (x1, vel1)
          call interp_vel (x2, vel2)

          vel2_n = dot_product (vel2, n_hat)
        
          vel1_n = vel2_n * (phi1 / phi2)**2

          vel1 = vel1 + (vel2_n * (phi1 / phi2)**2 - vel1_n)* n_hat 

          !--only take u,v components
          udes(i, j, k) = vel1(1)
          vdes(i, j, k) = vel1(2)

          $if ($IFORT || $IFC)

            if (isnan (udes(i, j, k))) then
              call error (sub_name, 'Nan at (i, j, k) =', (/ i, j, k /))
            end if

            if (isnan (vdes(i, j, k))) then
              call error (sub_name, 'NaN at (i, j, k) =', (/ i, j, k /))
            end if

          $endif

        !else
          !--do not know what to put here
        !end if

      end if

    end do
  end do
end do

!--w-node pass
do k = 2, nz - 1  !--(-1) here due to BOGUS
  do j = 1, ny
    do i = 1, nx

      phi1 = 0.5_rp * (phi(i, j, k) + phi(i, j, k-1))

      if ((phi0 < phi1) .and. (phi1 < phic)) then
      
        x1 = (/ (i-1) * dx, (j-1) * dy, (k-1) * dz /)

        n_hat = 0.5_rp * (norm(:, i, j, k) + norm(:, i, j, k-1))
        n_hat = n_hat / mag (n_hat)

        x2 = x1 + dphi * n_hat

        phi2 = phi1 + dphi  !--assume this is approx. correct

        call interp_phi (x2, phi_x2)

        !if (phi_x2 > phi0) then

          call interp_vel (x1, vel1)
          call interp_vel (x2, vel2)

          if (k == nz) then
            call mesg (sub_name, '(i, j) =', (/ i, j /))
            call mesg (sub_name, 'vel2 =', vel2)
            call mesg (sub_name, 'u =', u(i, j, k))
            call mesg (sub_name, 'v =', v(i, j, k))
            call mesg (sub_name, 'w =', w(i, j, k))
          end if

          vel2_n = dot_product (vel2, n_hat)
        
          vel1_n = vel2_n * (phi1 / phi2)**2

          vel1 = vel1 + (vel2_n * (phi1 / phi2)**2 - vel1_n)* n_hat 

          !--only take w-component
          wdes(i, j, k) = vel1(3)

          $if ($IFORT || $IFC)

            if (isnan (wdes(i, j, k))) then
              call error (sub_name, 'NaN at (i, j, k) =', (/ i, j, k /))
            end if

          $endif

        !else
          !--do not know what to put here
        !end if

      end if

    end do
  end do
end do

$if ($VERBOSE)
call exit_sub (sub_name)
$endif

end subroutine enforce_un

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine enforce_log_profile ()
use param, only : jt  !--plus stuff above
use sim_param, only : u, v, w
implicit none

character (*), parameter :: sub_name = mod_name // '.enforce_log_profile'

!logical, parameter :: DEBUG = .false.

integer :: i, j, k

real (rp) :: dphi
real (rp) :: phi0, phi1, phi2, phic
real (rp) :: phi_x2
real (rp) :: x1(nd), x2(nd)
real (rp) :: n_hat(nd)
real (rp) :: vel1(nd), vel2(nd)
real (rp) :: vel1_n, vel2_n
real (rp) :: vel1_p(nd), vel2_p(nd)

!---------------------------------------------------------------------

!--need to experiment with these values
if (phi_0_is_set) then
  phi0 = phi_0
else
  call error (sub_name, 'trying to use uninitialized phi_0')
end if

if (phi_cutoff_is_set) then
  phic = phi_cutoff
else
  call error (sub_name, 'trying to use uninitialized phi_cutoff')
end if

!phic = filter_size * dx
!phic = filter_size * sqrt (dx**2 + dy**2 + dz**2)

!dphi = filter_size * dx
dphi = filter_size * sqrt (dx**2 + dy**2 + dz**2)

udes = huge (1._rp)
vdes = huge (1._rp)
wdes = huge (1._rp)

!--u-node pass
do k = 1, nz - 1
  do j = 1, ny
    do i = 1, nx

      phi1 = phi(i, j, k)

      if ((phi0 < phi1) .and. (phi1 < phic)) then
      
        x1 = (/ (i-1) * dx, (j-1) * dy, (k-0.5_rp) * dz /)

        n_hat = norm(:, i, j, k)

        x2 = x1 + dphi * n_hat

        phi2 = phi1 + dphi  !--assume this is approx. correct

        call interp_phi (x2, phi_x2)

        !if (phi_x2 > phi0) then

          call interp_vel (x2, vel2)

          vel2_n = dot_product (vel2, n_hat)
          vel2_p = vel2 - vel2_n * n_hat
        
          vel1_n = vel2_n * (phi1 / phi2)**2
          vel1_p = vel2_p * (log (1._rp + phi1 / z0) / log (1._rp + phi2 / z0))

          vel1 = vel1_p + vel1_n * n_hat

          !--only take u,v components
          udes(i, j, k) = vel1(1)
          vdes(i, j, k) = vel1(2)

          $if ($IFORT || $IFC)

            if (isnan (udes(i, j, k))) then
              call error (sub_name, 'Nan at (i, j, k) =', (/ i, j, k /))
            end if

            if (isnan (vdes(i, j, k))) then
              call error (sub_name, 'NaN at (i, j, k) =', (/ i, j, k /))
            end if

          $endif

        !else
          !--do not know what to put here
        !end if

      end if

    end do
  end do
end do

!--w-node pass
do k = 2, nz - 1  !--(-1) here due to BOGUS
  do j = 1, ny
    do i = 1, nx

      phi1 = 0.5_rp * (phi(i, j, k) + phi(i, j, k-1))

      if ((phi0 < phi1) .and. (phi1 < phic)) then
      
        x1 = (/ (i-1) * dx, (j-1) * dy, (k-1) * dz /)

        n_hat = 0.5_rp * (norm(:, i, j, k) + norm(:, i, j, k-1))
        n_hat = n_hat / mag (n_hat)

        x2 = x1 + dphi * n_hat

        phi2 = phi1 + dphi  !--assume this is approx. correct

        call interp_phi (x2, phi_x2)

        !if (phi_x2 > phi0) then

          call interp_vel (x2, vel2)

          if (k == nz) then
            call mesg (sub_name, '(i, j) =', (/ i, j /))
            call mesg (sub_name, 'vel2 =', vel2)
            call mesg (sub_name, 'u =', u(i, j, k))
            call mesg (sub_name, 'v =', v(i, j, k))
            call mesg (sub_name, 'w =', w(i, j, k))
          end if

          vel2_n = dot_product (vel2, n_hat)
          vel2_p = vel2 - vel2_n * n_hat
        
          vel1_n = vel2_n * (phi1 / phi2)**2
          vel1_p = vel2_p * (log (1._rp + phi1 / z0) / log (1._rp + phi2 / z0))

          vel1 = vel1_p + vel1_n * n_hat

          !--only take w-component
          wdes(i, j, k) = vel1(3)

          $if ($IFORT || $IFC)

            if (isnan (wdes(i, j, k))) then
              call error (sub_name, 'NaN at (i, j, k) =', (/ i, j, k /))
            end if

          $endif

        !else
          !--do not know what to put here
        !end if

      end if

    end do
  end do
end do

end subroutine enforce_log_profile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--performs tri-linear interpolation to obtain a at point x
!--assumes a is on u-nodes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine interp_scal (albz, a, nbot, abot, ntop, atop, x, a_x, node)
use grid_defs, only : grid_t !autowrap_i, autowrap_j
use functions, only : cell_indx
use messages
implicit none

integer, intent (in) :: albz
real (rp), intent (in) :: a(ld, ny, albz:nz)  !--albz usually 0 or 1
integer, intent (in) :: nbot, ntop
real (rp), intent (in) :: abot(ld, ny, nbot)
real (rp), intent (in) :: atop(ld, ny, ntop)
real (rp), intent (in) :: x(nd)

real (rp), intent (out) :: a_x

character (*), intent (in), optional :: node  !--'u' or 'w'

character (*), parameter :: sub_name = mod_name // '.interp_scal'

integer :: i, j, k
integer :: i1, j1, k1
integer :: ks, ks1

integer, pointer, dimension(:) :: autowrap_i, autowrap_j

real (rp) :: s
real (rp) :: x1, x2, x3
real (rp) :: f1, f2, f3, f4, f5, f6, f7, f8
real (rp) :: w1, w2, w3, w4, w5, w6, w7, w8

real (rp) :: xmod(nd) ! Spatial location of autowrapped point

$if($VERBOSE)
call enter_sub( sub_name )
$endif

nullify(autowrap_i, autowrap_j)

autowrap_i => grid_t % autowrap_i
autowrap_j => grid_t % autowrap_j

!---------------------------------------------------------------------
xmod=x ! Initialize
xmod(1)=modulo(x(1),L_x) ! Ensures i is located in the domain
xmod(2)=modulo(x(2),L_y) ! Ensures j is located in the domain

!--calculate indices
!  Could replace floor() with cell_indx (JSG)
!i = autowrap_i( floor (xmod(1) / dx + 1._rp) )
!j = autowrap_j( floor (xmod(2) / dy + 1._rp) )

i = autowrap_i( cell_indx('i', dx, xmod(1)) )
j = autowrap_j( cell_indx('j', dy, xmod(2)) )

if (present (node)) then

  select case (node)
    case ('u'); s = 0.5_rp
    case ('w'); s = 1._rp
    case default; call error (sub_name, 'invalid node = ' // node)
  end select

else  !-default to u-nodes

  s = 0.5_rp

end if

!  This needs special treatment for u and w
ks = floor (xmod(3) / dz + s)

if ((i < 1) .or. (i > nx)) then
  write (msg, *) 'i out of range, i = ', i, n_l,  &
                 'x = ', x, n_l,                  &
                 '(i, j, ks) = ', i, j, ks
  call error (sub_name, msg)
end if

if ((j < 1) .or. (j > ny)) then
  call error (sub_name, 'j out of range, j =', j)
end if


!--need to bounds check ks
if ( (.not. USE_MPI) .or. (USE_MPI .and. (coord == 0)) ) then
  if (ks < 1) call error (sub_name, 'ks out of range, ks =', ks)
else
  if (ks < -nbot + albz) call error (sub_name,                       &
                                     'out of range, (ks, ksmin) =',  &
                                     (/ ks, -nbot + albz /))
end if

if ( (.not. USE_MPI) .or. (USE_MPI .and. (coord == nproc - 1)) ) then
  if (ks > nz - 1) call error (sub_name, 'ks out of range, ks =', ks)
else
  if (ks > nz + ntop - 1) call error (sub_name,                      &
                                      'out of range (ks, ksmax) =',  &
                                      (/ ks, nz + ntop - 1 /))
end if

!--try to handle boundaries nicely for the +1 indices
!i1 = modulo (i, nx) + 1
!j1 = modulo (j, ny) + 1
i1 = autowrap_i( i + 1 )
j1 = autowrap_j( j + 1 )

ks1 = ks + 1

!--calculate interpolation weights
x1 = modulo (xmod(1), dx) / dx
x2 = modulo (xmod(2), dy) / dy
x3 = xmod(3) / dz - (floor (xmod(3) / dz + s) - s)

w1 = (1._rp - x1) * (1._rp - x2) * (1._rp - x3)
w2 = (    x1    ) * (1._rp - x2) * (1._rp - x3)
w3 = (1._rp - x1) * (    x2    ) * (1._rp - x3)
w4 = (    x1    ) * (    x2    ) * (1._rp - x3)
w5 = (1._rp - x1) * (1._rp - x2) * (    x3    )
w6 = (    x1    ) * (1._rp - x2) * (    x3    )
w7 = (1._rp - x1) * (    x2    ) * (    x3    )
w8 = (    x1    ) * (    x2    ) * (    x3    )

$if ($MPI)

  if (ks < albz) then
    k = nbot + ks + 1 - albz
    f1 = abot(i , j , k)
    f2 = abot(i1, j , k)
    f3 = abot(i , j1, k)
    f4 = abot(i1, j1, k)
  else if (ks > nz) then
    k = ks - nz
    f1 = atop(i , j , k)
    f2 = atop(i1, j , k)
    f3 = atop(i , j1, k)
    f4 = atop(i1, j1, k)
  else
    k = ks
    f1 = a(i , j , k)
    f2 = a(i1, j , k)
    f3 = a(i , j1, k)
    f4 = a(i1, j1, k)
  end if

  if (ks1 < albz) then
    k1 = nbot + ks1 + 1 - albz
    f5 = abot(i , j , k1)
    f6 = abot(i1, j , k1)
    f7 = abot(i , j1, k1)
    f8 = abot(i1, j1, k1)
  else if (ks1 > nz) then
    k1 = ks1 - nz
    f5 = atop(i , j , k1)
    f6 = atop(i1, j , k1)
    f7 = atop(i , j1, k1)
    f8 = atop(i1, j1, k1)
  else
    k1 = ks1
    f5 = a(i , j , k1)
    f6 = a(i1, j , k1)
    f7 = a(i , j1, k1)
    f8 = a(i1, j1, k1)
  end if

$else

  f1 = a(i , j , ks)
  f2 = a(i1, j , ks)
  f3 = a(i , j1, ks)
  f4 = a(i1, j1, ks)
  f5 = a(i , j , ks1)
  f6 = a(i1, j , ks1)
  f7 = a(i , j1, ks1)
  f8 = a(i1, j1, ks1)

$endif

a_x = w1 * f1 + w2 * f2 + w3 * f3 + w4 * f4 +  &
      w5 * f5 + w6 * f6 + w7 * f7 + w8 * f8

nullify(autowrap_i, autowrap_j)

$if($VERBOSE)
call exit_sub( sub_name )
$endif

end subroutine interp_scal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine interp_tij_u (x, txx_x, txy_x, tyy_x, tzz_x)
use sim_param, only : txx, txy, tyy, tzz
use functions, only : cell_indx
use grid_defs, only : grid_t !autowrap_i, autowrap_j
use messages

implicit none

real (rp), intent (in) :: x(nd)

real (rp), intent (out) :: txx_x, txy_x, tyy_x, tzz_x

character (*), parameter :: sub_name = mod_name // '.interp_tij_u'

integer :: i, j, k
integer :: i1, j1, k1
integer :: ku, ku1

integer, pointer, dimension(:) :: autowrap_i, autowrap_j

real (rp) :: s
real (rp) :: x1, x2, x3
real (rp) :: w(8)
real (rp) :: f(8)

real (rp) :: xmod(nd) ! Spatial location of autowrapped point

nullify(autowrap_i, autowrap_j)

$if($VERBOSE)
call enter_sub( sub_name )
$endif

autowrap_i => grid_t % autowrap_i
autowrap_j => grid_t % autowrap_j

!---------------------------------------------------------------------
xmod=x ! Initialize
xmod(1)=modulo(x(1),L_x) ! Ensures i is located in the domain
xmod(2)=modulo(x(2),L_y) ! Ensures j is located in the domain

!--calculate indices
!i = floor (xmod(1) / dx + 1._rp)
!j = floor (xmod(2) / dy + 1._rp)
i = autowrap_i( cell_indx('i', dx, xmod(1)) )
j = autowrap_j( cell_indx('j', dy, xmod(2)) )

s = 0.5_rp
ku = floor (xmod(3) / dz + s)

if ((i < 1) .or. (i > nx)) then
  call error (sub_name, 'i out of range, i =', i)
end if

if ((j < 1) .or. (j > ny)) then
  call error (sub_name, 'j out of range, j =', j)
end if

!--need to bounds check ku
if ( (.not. USE_MPI) .or. (USE_MPI .and. (coord == 0)) ) then
  if (ku < 1) call error (sub_name, 'ku out of range, ku =', ku)
else
  if (ku < -ntaubot) call error (sub_name,                       &
                                 'out of range, (ku, kumin) =',  &
                                 (/ ku, -ntaubot /))
end if

if ( (.not. USE_MPI) .or. (USE_MPI .and. (coord == nproc - 1)) ) then
  if (ku > nz - 1) call error (sub_name, 'ku out of range, ku =', ku)
else
  if (ku > nz + ntautop - 1) call error (sub_name,                     &
                                         'out of range (ku, kumax) =', &
                                         (/ ku, nz + ntautop - 1 /))
end if

!--try to handle boundaries nicely for the +1 indices

!i1 = modulo (i, nx) + 1
!j1 = modulo (j, ny) + 1
i1 = autowrap_i( i + 1 )
j1 = autowrap_j( j + 1 )

ku1 = ku + 1

!--calculate interpolation weights
x1 = modulo (xmod(1), dx) / dx
x2 = modulo (xmod(2), dy) / dy
x3 = xmod(3) / dz - (floor (xmod(3) / dz + s) - s)

w(1) = (1._rp - x1) * (1._rp - x2) * (1._rp - x3)
w(2) = (    x1    ) * (1._rp - x2) * (1._rp - x3)
w(3) = (1._rp - x1) * (    x2    ) * (1._rp - x3)
w(4) = (    x1    ) * (    x2    ) * (1._rp - x3)
w(5) = (1._rp - x1) * (1._rp - x2) * (    x3    )
w(6) = (    x1    ) * (1._rp - x2) * (    x3    )
w(7) = (1._rp - x1) * (    x2    ) * (    x3    )
w(8) = (    x1    ) * (    x2    ) * (    x3    )

call fill_f (txxbot, txxtop, txx)

if ((abs (x(1) - 0.385684_rp) < 0.001_rp) .and.  &
    (abs (x(2) - 0.000000_rp) < 0.001_rp) .and.  &
    (abs (x(3) - 0.507011_rp) < 0.001_rp)) then
  call mesg (sub_name, 'i=', i)
  call mesg (sub_name, 'j=', j)
  call mesg (sub_name, 'ku=', ku)
  call mesg (sub_name, 'ku1=', ku1)
  call mesg (sub_name, 'nz=', nz)
  call mesg (sub_name, 'w=', w)
  call mesg (sub_name, 'f=', f)
  call mesg (sub_name, 'txx(i,j,ku1)=', txx(i, j, ku1))
end if

txx_x = w(1) * f(1) + w(2) * f(2) + w(3) * f(3) + w(4) * f(4) +  &
        w(5) * f(5) + w(6) * f(6) + w(7) * f(7) + w(8) * f(8)

call fill_f (txybot, txytop, txy)

txy_x = w(1) * f(1) + w(2) * f(2) + w(3) * f(3) + w(4) * f(4) +  &
        w(5) * f(5) + w(6) * f(6) + w(7) * f(7) + w(8) * f(8)

call fill_f (tyybot, tyytop, tyy)

tyy_x = w(1) * f(1) + w(2) * f(2) + w(3) * f(3) + w(4) * f(4) +  &
        w(5) * f(5) + w(6) * f(6) + w(7) * f(7) + w(8) * f(8)

call fill_f (tzzbot, tzztop, tzz)

tzz_x = w(1) * f(1) + w(2) * f(2) + w(3) * f(3) + w(4) * f(4) +  &
        w(5) * f(5) + w(6) * f(6) + w(7) * f(7) + w(8) * f(8)

nullify(autowrap_i, autowrap_j)        

$if($VERBOSE)
call exit_sub( sub_name )
$endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fill_f (abot, atop, a)
implicit none

real (rp), intent (in) :: abot(ld, ny, ntaubot)
real (rp), intent (in) :: atop(ld, ny, ntautop)
real (rp), intent (in) :: a(ld, ny, $lbz:nz)  !--since tij are $lbz:nz

!---------------------------------------------------------------------

$if ($MPI)

  if (ku < 0) then
    k = ntaubot + ku + 1
    f(1) = abot(i , j , k)
    f(2) = abot(i1, j , k)
    f(3) = abot(i , j1, k)
    f(4) = abot(i1, j1, k)
  else if (ku > nz) then
    k = ku - nz
    f(1) = atop(i , j , k)
    f(2) = atop(i1, j , k)
    f(3) = atop(i , j1, k)
    f(4) = atop(i1, j1, k)
  else
    k = ku
    f(1) = a(i , j , k)
    f(2) = a(i1, j , k)
    f(3) = a(i , j1, k)
    f(4) = a(i1, j1, k)
  end if

  if (ku1 < 0) then
    k1 = nvelbot + ku1 + 1
    f(5) = abot(i , j , k1)
    f(6) = abot(i1, j , k1)
    f(7) = abot(i , j1, k1)
    f(8) = abot(i1, j1, k1)
  else if (ku1 > nz) then
    k1 = ku1 - nz
    f(5) = atop(i , j , k1)
    f(6) = atop(i1, j , k1)
    f(7) = atop(i , j1, k1)
    f(8) = atop(i1, j1, k1)
  else
    k1 = ku1
    f(5) = a(i , j , k1)
    f(6) = a(i1, j , k1)
    f(7) = a(i , j1, k1)
    f(8) = a(i1, j1, k1)
  end if

$else

  f(1) = a(i , j , ku )
  f(2) = a(i1, j , ku )
  f(3) = a(i , j1, ku )
  f(4) = a(i1, j1, ku )
  f(5) = a(i , j , ku1)
  f(6) = a(i1, j , ku1)
  f(7) = a(i , j1, ku1)
  f(8) = a(i1, j1, ku1)

$endif

end subroutine fill_f

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine interp_tij_u

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine interp_tij_w (x, txz_x, tyz_x)
use sim_param, only : txz, tyz
use functions, only : cell_indx
use grid_defs, only : grid_t !autowrap_i, autowrap_j
use messages

implicit none

real (rp), intent (in) :: x(nd)

real (rp), intent (out) :: txz_x, tyz_x

character (*), parameter :: sub_name = mod_name // '.interp_tij_w'

integer :: i, j, k
integer :: i1, j1, k1
integer :: kw, kw1

integer, pointer, dimension(:) :: autowrap_i, autowrap_j

real (rp) :: s
real (rp) :: x1, x2, x3
real (rp) :: w(8)
real (rp) :: f(8)

real (rp) :: xmod(nd) ! Spatial location of autowrapped point

nullify(autowrap_i, autowrap_j)

$if($VERBOSE)
call enter_sub( sub_name )
$endif

autowrap_i => grid_t % autowrap_i
autowrap_j => grid_t % autowrap_j

!---------------------------------------------------------------------
xmod=x ! Initialize
xmod(1)=modulo(x(1),L_x) ! Ensures i is located in the domain
xmod(2)=modulo(x(2),L_y) ! Ensures j is located in the domain

!--calculate indices
!i = floor (xmod(1) / dx + 1._rp)
!j = floor (xmod(2) / dy + 1._rp)
i = autowrap_i( cell_indx('i', dx, xmod(1)) )
j = autowrap_j( cell_indx('j', dy, xmod(2)) )

s = 1._rp
kw = floor (xmod(3) / dz + s)

if ((i < 1) .or. (i > nx)) then
  call error (sub_name, 'i out of range, i =', i)
end if

if ((j < 1) .or. (j > ny)) then
  call error (sub_name, 'j out of range, j =', j)
end if

!--need to bounds check kw
if ( (.not. USE_MPI) .or. (USE_MPI .and. (coord == 0)) ) then
  if (kw < 1) call error (sub_name, 'kw out of range, kw =', kw)
else
  if (kw < -ntaubot) call error (sub_name,                       &
                                 'out of range, (kw, kwmin) =',  &
                                 (/ kw, -ntaubot /))
end if

if ( (.not. USE_MPI) .or. (USE_MPI .and. (coord == nproc - 1)) ) then
  if (kw > nz - 1) call error (sub_name, 'kw out of range, kw =', kw)
else
  if (kw > nz + ntautop - 1) call error (sub_name,                     &
                                         'out of range (kw, kwmax) =', &
                                         (/ kw, nz + ntautop - 1 /)) 
end if

!--try to handle boundaries nicely for the +1 indices

!i1 = modulo (i, nx) + 1
!j1 = modulo (j, ny) + 1
i1 = autowrap_i( i + 1 )
j1 = autowrap_j( j + 1 )

kw1 = kw + 1

!--calculate interpolation weights
x1 = modulo (xmod(1), dx) / dx
x2 = modulo (xmod(2), dy) / dy
x3 = xmod(3) / dz - (floor (xmod(3) / dz + s) - s)

w(1) = (1._rp - x1) * (1._rp - x2) * (1._rp - x3)
w(2) = (    x1    ) * (1._rp - x2) * (1._rp - x3)
w(3) = (1._rp - x1) * (    x2    ) * (1._rp - x3)
w(4) = (    x1    ) * (    x2    ) * (1._rp - x3)
w(5) = (1._rp - x1) * (1._rp - x2) * (    x3    )
w(6) = (    x1    ) * (1._rp - x2) * (    x3    )
w(7) = (1._rp - x1) * (    x2    ) * (    x3    )
w(8) = (    x1    ) * (    x2    ) * (    x3    )

call fill_f (txzbot, txztop, txz)

txz_x = w(1) * f(1) + w(2) * f(2) + w(3) * f(3) + w(4) * f(4) +  &
        w(5) * f(5) + w(6) * f(6) + w(7) * f(7) + w(8) * f(8)

call fill_f (tyzbot, tyztop, tyz)

tyz_x = w(1) * f(1) + w(2) * f(2) + w(3) * f(3) + w(4) * f(4) +  &
        w(5) * f(5) + w(6) * f(6) + w(7) * f(7) + w(8) * f(8)

nullify(autowrap_i, autowrap_j)

$if($VERBOSE)
call exit_sub( sub_name )
$endif

!**********************************************************************
contains
!**********************************************************************

!**********************************************************************
subroutine fill_f (abot, atop, a)
!**********************************************************************
implicit none

real (rp), intent (in) :: abot(ld, ny, ntaubot)
real (rp), intent (in) :: atop(ld, ny, ntautop)
real (rp), intent (in) :: a(ld, ny, $lbz:nz)  !--since tij are $lbz:nz

!---------------------------------------------------------------------

$if ($MPI)

  if (kw < 0) then
    k = ntaubot + kw + 1
    f(1) = abot(i , j , k)
    f(2) = abot(i1, j , k)
    f(3) = abot(i , j1, k)
    f(4) = abot(i1, j1, k)
  else if (kw > nz) then
    k = kw - nz
    f(1) = atop(i , j , k)
    f(2) = atop(i1, j , k)
    f(3) = atop(i , j1, k)
    f(4) = atop(i1, j1, k)
  else
    k = kw
    f(1) = a(i , j , k)
    f(2) = a(i1, j , k)
    f(3) = a(i , j1, k)
    f(4) = a(i1, j1, k)
  end if

  if (kw1 < 0) then
    k1 = nvelbot + kw1 + 1
    f(5) = abot(i , j , k1)
    f(6) = abot(i1, j , k1)
    f(7) = abot(i , j1, k1)
    f(8) = abot(i1, j1, k1)
  else if (kw1 > nz) then
    k1 = kw1 - nz
    f(5) = atop(i , j , k1)
    f(6) = atop(i1, j , k1)
    f(7) = atop(i , j1, k1)
    f(8) = atop(i1, j1, k1)
  else
    k1 = kw1
    f(5) = a(i , j , k1)
    f(6) = a(i1, j , k1)
    f(7) = a(i , j1, k1)
    f(8) = a(i1, j1, k1)
  end if

$else

  f(1) = a(i , j , kw )
  f(2) = a(i1, j , kw )
  f(3) = a(i , j1, kw )
  f(4) = a(i1, j1, kw )
  f(5) = a(i , j , kw1)
  f(6) = a(i1, j , kw1)
  f(7) = a(i , j1, kw1)
  f(8) = a(i1, j1, kw1)

$endif

end subroutine fill_f

end subroutine interp_tij_w


!**********************************************************************
subroutine interp_phi (x, phi_x)
!**********************************************************************
!
!--performs tri-linear interpolation to obtain phi at point x
!--assumes phi is on u-nodes
!
use functions, only : cell_indx
use grid_defs, only : grid_t !autowrap_i, autowrap_j
use messages
implicit none

real (rp), intent (in) :: x(nd)
real (rp), intent (out) :: phi_x

character (*), parameter :: sub_name = mod_name // '.interp_phi'

integer :: i, j, k, ku
integer :: i1, j1, k1, ku1

integer, pointer, dimension(:) :: autowrap_i, autowrap_j

real (rp) :: x1, x2, x3
real (rp) :: w(8), f(8)

real (rp) :: xmod(nd) ! Spatial location of autowrapped point

nullify(autowrap_i, autowrap_j)

$if($VERBOSE)
call enter_sub( sub_name )
$endif

autowrap_i => grid_t % autowrap_i
autowrap_j => grid_t % autowrap_j

!---------------------------------------------------------------------
xmod=x ! Initialize
xmod(1)=modulo(x(1),L_x) ! Ensures i is located in the domain
xmod(2)=modulo(x(2),L_y) ! Ensures j is located in the domain

!--calculate indices
!i = floor (xmod(1) / dx + 1._rp)
!j = floor (xmod(2) / dy + 1._rp)
i = autowrap_i( cell_indx( 'i', dx, xmod(1) ) )
j = autowrap_j( cell_indx( 'j', dy, xmod(2) ) )

ku = floor (xmod(3) / dz + 0.5_rp) !--assumes phi on u-nodes

!--need to bounds check i, j, ku, kw
if ((i < 1) .or. (i > nx)) then
  write (msg, *) 'i out of range', n_l,            &
                 '(i, j, ku) = ', i, j, ku, n_l,   &
                 'x = ', x
  call error (sub_name, msg)
end if

if ((j < 1) .or. (j > ny)) then
  write (msg, *) 'j out of range', n_l,            &
                 '(i, j, ku) = ', i, j, ku, n_l,   &
                 'x = ', x
  call error (sub_name, msg)
end if

!--need to bounds check ku
!--non-MPI has 0-sized nphibot, nphitop
if ( (.not. USE_MPI) .or. (USE_MPI .and. (coord == 0)) ) then
  if (ku < 1) call error (sub_name, 'ku out of range, ku =', ku)
else
  if (ku < -nphibot) call error (sub_name,                     &
                                 'out of range, (ku, kumin) =',  &
                                 (/ ku, -nphibot /))
end if

if ( (.not. USE_MPI) .or. (USE_MPI .and. (coord == nproc - 1)) ) then
  if (ku > nz - 1) call error (sub_name, 'ku out of range, ku =', ku)
else
  !--(-1) needed since we will use ku + 1 below!!!!
  if (ku > nz + nphitop - 1) call error (sub_name,                      &
                                         'out of range (ku, kumax) =',  &
                                         (/ ku, nz + nphitop - 1 /))
end if

!--try to handle boundaries nicely for the +1 indices
!i1 = modulo (i, nx) + 1
!j1 = modulo (j, ny) + 1
i1 = autowrap_i( i + 1 )
j1 = autowrap_j( j + 1 )

ku1 = ku + 1

!--calculate interpolation weights
!  Computes fraction of dx,dy,dz that point exists from 
!  starting i,j,k of cell
x1 = modulo (xmod(1), dx) / dx
x2 = modulo (xmod(2), dy) / dy
x3 = xmod(3) / dz - (floor (xmod(3) / dz + 0.5_rp) - 0.5_rp)

w(1) = (1._rp - x1) * (1._rp - x2) * (1._rp - x3)
w(2) = (    x1    ) * (1._rp - x2) * (1._rp - x3)
w(3) = (1._rp - x1) * (    x2    ) * (1._rp - x3)
w(4) = (    x1    ) * (    x2    ) * (1._rp - x3)
w(5) = (1._rp - x1) * (1._rp - x2) * (    x3    )
w(6) = (    x1    ) * (1._rp - x2) * (    x3    )
w(7) = (1._rp - x1) * (    x2    ) * (    x3    )
w(8) = (    x1    ) * (    x2    ) * (    x3    )

!--u-nodes
$if ($MPI)

  if (ku < 0) then
    k = nphibot + ku + 1
    f(1) = phibot(i , j , k)
    f(2) = phibot(i1, j , k)
    f(3) = phibot(i , j1, k)
    f(4) = phibot(i1, j1, k)
  else if (ku > nz) then
    k = ku - nz
    f(1) = phitop(i , j , k)
    f(2) = phitop(i1, j , k)
    f(3) = phitop(i , j1, k)
    f(4) = phitop(i1, j1, k)
  else
    k = ku
    f(1) = phi(i , j , k)
    f(2) = phi(i1, j , k)
    f(3) = phi(i , j1, k)
    f(4) = phi(i1, j1, k)
  end if

  if (ku1 < 0) then
    k1 = nphibot + ku1 + 1
    f(5) = phibot(i , j , k1)
    f(6) = phibot(i1, j , k1)
    f(7) = phibot(i , j1, k1)
    f(8) = phibot(i1, j1, k1)
  else if (ku1 > nz) then
    k1 = ku1 - nz
    f(5) = phitop(i , j , k1)
    f(6) = phitop(i1, j , k1)
    f(7) = phitop(i , j1, k1)
    f(8) = phitop(i1, j1, k1)
  else
    k1 = ku1
    f(5) = phi(i , j , k1)
    f(6) = phi(i1, j , k1)
    f(7) = phi(i , j1, k1)
    f(8) = phi(i1, j1, k1)
  end if

$else

  f(1) = phi(i , j , ku)
  f(2) = phi(i1, j , ku)
  f(3) = phi(i , j1, ku)
  f(4) = phi(i1, j1, ku)
  f(5) = phi(i , j , ku1)
  f(6) = phi(i1, j , ku1)
  f(7) = phi(i , j1, ku1)
  f(8) = phi(i1, j1, ku1)

$endif

phi_x = w(1) * f(1) + w(2) * f(2) + w(3) * f(3) + w(4) * f(4) +  &
        w(5) * f(5) + w(6) * f(6) + w(7) * f(7) + w(8) * f(8)

nullify(autowrap_i, autowrap_j)

$if($VERBOSE)
call exit_sub( sub_name )
$endif

end subroutine interp_phi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--performs tri-linear interpolation to obtain velocity field at point x
!--special care is required near bdrys: watch out for BOGUS!
!--this may be out-of-date with use_log_profile
!--this assumes u, v, w, at k=0, nz are not BOGUS, and are useable
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine interp_vel (x, vel)
use sim_param, only : u, v, w
use functions, only : cell_indx
use grid_defs, only : grid_t !autowrap_i, autowrap_j
use messages

implicit none

real (rp), intent (in) :: x(nd)
real (rp), intent (out) :: vel(nd)

character (*), parameter :: sub_name = mod_name // '.interp_vel'

integer :: i, j, k, ku, kw
integer :: i1, j1, k1, ku1, kw1

integer, pointer, dimension(:) :: autowrap_i, autowrap_j

real (rp) :: x1, x2, x3u, x3w
real (rp) :: w1, w2, w3, w4, w5, w6, w7, w8
real (rp) :: f1, f2, f3, f4, f5, f6, f7, f8

real (rp) :: xmod(nd) ! Spatial location of autowrapped point

nullify(autowrap_i, autowrap_j)

$if($VERBOSE)
call enter_sub( sub_name )
$endif

autowrap_i => grid_t % autowrap_i
autowrap_j => grid_t % autowrap_j

!---------------------------------------------------------------------
xmod=x ! Initialize
xmod(1)=modulo(x(1),L_x) ! Ensures i is located in the domain
xmod(2)=modulo(x(2),L_y) ! Ensures j is located in the domain

!--calculate indices
!i = floor (xmod(1) / dx + 1._rp)
!j = floor (xmod(2) / dy + 1._rp)
i = autowrap_i( cell_indx( 'i', dx, xmod(1) ) )
j = autowrap_j( cell_indx( 'j', dy, xmod(2) ) )
ku = floor (xmod(3) / dz + 0.5_rp)  !--assumes phi on u-nodes
kw = floor (xmod(3) / dz + 1._rp)

if ((i < 1) .or. (i > nx)) then
  call error (sub_name, 'i out of range, i =', i)
end if

if ((j < 1) .or. (j > ny)) then
  call error (sub_name, 'j out of range, j =', j)
end if

!--need to bounds check ku, kw
if ( (.not. USE_MPI) .or. (USE_MPI .and. (coord == 0)) ) then
  if (ku < 1) call error (sub_name, 'ku out of range, ku =', ku)
  if (kw < 1) call error (sub_name, 'kw out of range, kw =', kw)
else
  if (ku < -nvelbot) call error (sub_name,                       &
                                 'out of range, (ku, kumin) =',  &
                                 (/ ku, -nvelbot /))
  if (kw < -nvelbot) call error (sub_name,                       &
                                 'out of range, (kw, kwmin) =',  &
                                 (/ kw, -nvelbot /))
end if

if ( (.not. USE_MPI) .or. (USE_MPI .and. (coord == nproc - 1)) ) then
  if (ku > nz - 1) call error (sub_name, 'ku out of range, ku =', ku)
  if (kw > nz - 1) call error (sub_name, 'kw out of range, kw =', kw)
else
  if (ku > nz + nveltop - 1) call error (sub_name,                     &
                                         'out of range (ku, kumax) =', &
                                         (/ ku, nz + nveltop - 1 /))
  if (kw > nz + nveltop - 1) call error (sub_name,                      &
                                         'out of range (kw, kwmax) =',  &
                                         (/ kw, nz + nveltop - 1 /))
end if

!--try to handle boundaries nicely for the +1 indices
!i1 = modulo (i, nx) + 1
!j1 = modulo (j, ny) + 1
i1 = autowrap_i( i + 1 )
j1 = autowrap_j( j + 1 )

ku1 = ku + 1
kw1 = kw + 1

!--calculate interpolation weights
x1 = modulo (xmod(1), dx) / dx
x2 = modulo (xmod(2), dy) / dy
x3u = xmod(3) / dz - (floor (xmod(3) / dz + 0.5_rp) - 0.5_rp)
x3w = modulo (xmod(3), dz) / dz

w1 = (1._rp - x1) * (1._rp - x2) * (1._rp - x3u)
w2 = (    x1    ) * (1._rp - x2) * (1._rp - x3u)
w3 = (1._rp - x1) * (    x2    ) * (1._rp - x3u)
w4 = (    x1    ) * (    x2    ) * (1._rp - x3u)
w5 = (1._rp - x1) * (1._rp - x2) * (   x3u     )
w6 = (    x1    ) * (1._rp - x2) * (   x3u     )
w7 = (1._rp - x1) * (    x2    ) * (   x3u     )
w8 = (    x1    ) * (    x2    ) * (   x3u     )

!--u-nodes
$if ($MPI)

  if (ku < 0) then
    k = nvelbot + ku + 1
    f1 = ubot(i , j , k)
    f2 = ubot(i1, j , k)
    f3 = ubot(i , j1, k)
    f4 = ubot(i1, j1, k)
  else if (ku > nz) then
    k = ku - nz
    f1 = utop(i , j , k)
    f2 = utop(i1, j , k)
    f3 = utop(i , j1, k)
    f4 = utop(i1, j1, k)
  else
    k = ku
    f1 = u(i , j , k)
    f2 = u(i1, j , k)
    f3 = u(i , j1, k)
    f4 = u(i1, j1, k)
  end if

  if (ku1 < 0) then
    k1 = nvelbot + ku1 + 1
    f5 = ubot(i , j , k1)
    f6 = ubot(i1, j , k1)
    f7 = ubot(i , j1, k1)
    f8 = ubot(i1, j1, k1)
  else if (ku1 > nz) then
    k1 = ku1 - nz
    f5 = utop(i , j , k1)
    f6 = utop(i1, j , k1)
    f7 = utop(i , j1, k1)
    f8 = utop(i1, j1, k1)
  else
    k1 = ku1
    f5 = u(i , j , k1)
    f6 = u(i1, j , k1)
    f7 = u(i , j1, k1)
    f8 = u(i1, j1, k1)
  end if

$else

  f1 = u(i , j , ku)
  f2 = u(i1, j , ku)
  f3 = u(i , j1, ku)
  f4 = u(i1, j1, ku)
  f5 = u(i , j , ku1)
  f6 = u(i1, j , ku1)
  f7 = u(i , j1, ku1)
  f8 = u(i1, j1, ku1)

$endif

vel(1) = w1 * f1 + w2 * f2 + w3 * f3 + w4 * f4 +  &
         w5 * f5 + w6 * f6 + w7 * f7 + w8 * f8

!--u-nodes
$if ($MPI)

  if (ku < 0) then
    k = nvelbot + ku + 1
    f1 = vbot(i , j , k)
    f2 = vbot(i1, j , k)
    f3 = vbot(i , j1, k)
    f4 = vbot(i1, j1, k)
  else if (ku > nz) then
    k = ku - nz
    f1 = vtop(i , j , k)
    f2 = vtop(i1, j , k)
    f3 = vtop(i , j1, k)
    f4 = vtop(i1, j1, k)
  else
    k = ku
    f1 = v(i , j , k)
    f2 = v(i1, j , k)
    f3 = v(i , j1, k)
    f4 = v(i1, j1, k)
  end if

  if (ku1 < 0) then
    k1 = nvelbot + ku1 + 1
    f5 = vbot(i , j , k1)
    f6 = vbot(i1, j , k1)
    f7 = vbot(i , j1, k1)
    f8 = vbot(i1, j1, k1)
  else if (ku1 > nz) then
    k1 = ku1 - nz
    f5 = vtop(i , j , k1)
    f6 = vtop(i1, j , k1)
    f7 = vtop(i , j1, k1)
    f8 = vtop(i1, j1, k1)
  else
    k1 = ku1
    f5 = v(i , j , k1)
    f6 = v(i1, j , k1)
    f7 = v(i , j1, k1)
    f8 = v(i1, j1, k1)
  end if

$else

  f1 = v(i , j , ku )
  f2 = v(i1, j , ku )
  f3 = v(i , j1, ku )
  f4 = v(i1, j1, ku )
  f5 = v(i , j , ku1)
  f6 = v(i1, j , ku1)
  f7 = v(i , j1, ku1)
  f8 = v(i1, j1, ku1)

$endif

vel(2) = w1 * f1 + w2 * f2 + w3 * f3 + w4 * f4 +  &
         w5 * f5 + w6 * f6 + w7 * f7 + w8 * f8

!--weights for w-nodes
w1 = (1._rp - x1) * (1._rp - x2) * (1._rp - x3w)
w2 = (    x1    ) * (1._rp - x2) * (1._rp - x3w)
w3 = (1._rp - x1) * (    x2    ) * (1._rp - x3w)
w4 = (    x1    ) * (    x2    ) * (1._rp - x3w)
w5 = (1._rp - x1) * (1._rp - x2) * (   x3w     )
w6 = (    x1    ) * (1._rp - x2) * (   x3w     )
w7 = (1._rp - x1) * (    x2    ) * (   x3w     )
w8 = (    x1    ) * (    x2    ) * (   x3w     )

!--w-nodes
$if ($MPI)
  if (kw < 0) then
    k = nvelbot + kw + 1
    f1 = wbot(i , j , k)
    f2 = wbot(i1, j , k)
    f3 = wbot(i , j1, k)
    f4 = wbot(i1, j1, k)
  else if (kw > nz) then
    k = kw - nz
    f1 = wtop(i , j , k)
    f2 = wtop(i1, j , k)
    f3 = wtop(i , j1, k)
    f4 = wtop(i1, j1, k)
  else
    k = kw
    f1 = w(i , j , k)
    f2 = w(i1, j , k)
    f3 = w(i , j1, k)
    f4 = w(i1, j1, k)
  end if

  if (kw1 < 0) then
    k1 = nvelbot + kw1 + 1
    f5 = wbot(i , j , k1)
    f6 = wbot(i1, j , k1)
    f7 = wbot(i , j1, k1)
    f8 = wbot(i1, j1, k1)
  else if (kw1 > nz) then
    k1 = kw1 - nz
    f5 = wtop(i , j , k1)
    f6 = wtop(i1, j , k1)
    f7 = wtop(i , j1, k1)
    f8 = wtop(i1, j1, k1)
  else
    k1 = kw1
    f5 = w(i , j , k1)
    f6 = w(i1, j , k1)
    f7 = w(i , j1, k1)
    f8 = w(i1, j1, k1)
  end if

$else

  f1 = w(i , j , kw )
  f2 = w(i1, j , kw )
  f3 = w(i , j1, kw )
  f4 = w(i1, j1, kw )
  f5 = w(i , j , kw1)
  f6 = w(i1, j , kw1)
  f7 = w(i , j1, kw1)
  f8 = w(i1, j1, kw1)

$endif

vel(3) = w1 * f1 + w2 * f2 + w3 * f3 + w4 * f4 +  &
         w5 * f5 + w6 * f6 + w7 * f7 + w8 * f8

nullify(autowrap_i, autowrap_j)

$if($VERBOSE)
call exit_sub( sub_name )
$endif

end subroutine interp_vel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--a driver routine for calling smooth routine with stresses
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine smooth_tau (phi_c, txx, txy, txz, tyy, tyz, tzz)
implicit none

real (rp), intent (in) :: phi_c

real (rp), intent (in out), dimension (ld, ny, $lbz:nz) ::  &
                                 txx, txy, txz, tyy, tyz, tzz

character (*), parameter :: sub_name = mod_name // '.level_set_smooth_tau'

!---------------------------------------------------------------------
$if ($VERBOSE)
call enter_sub (sub_name)
$endif

!if (phi_cutoff_is_set) then
!  phi_c = -phi_cutoff - 10._rp * epsilon (0._rp)
!else
!  call error (sub_name, 'trying to use uninitialized phi_cutoff')
!end if

!phi_c = 0._rp  !--anything with phi < 0 is smoothed
!--experimental
!phi_c = -(filter_size * sqrt (dx**2 + dy**2 + dz**2) +  &
!          10._rp * epsilon (0._rp))

call smooth (phi_c, $lbz, txx)
call smooth (phi_c, $lbz, txy)
call smooth (phi_c, $lbz, txz, node='w')
call smooth (phi_c, $lbz, tyy)
call smooth (phi_c, $lbz, tyz, node='w')
call smooth (phi_c, $lbz, tzz)

$if ($VERBOSE)
call exit_sub (sub_name)
$endif

end subroutine smooth_tau

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--basically a driver routine for calling the smooth routine with
!  velocities
!--make sure phi_c is consistent with other bc routines
!--modifies u, v, w.  Be careful.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine level_set_smooth_vel (u, v, w)
implicit none

real (rp), intent (in out), dimension (ld, ny, $lbz:nz) :: u, v, w

character (*), parameter :: sub_name = mod_name // '.level_set_smooth_vel'

real (rp), parameter :: phi_c = 0._rp !--any pt with phi < 0 is smoothed

!---------------------------------------------------------------------
$if ($VERBOSE)
call enter_sub (sub_name)
$endif

call smooth (phi_c, lbound (u, 3), u)
call smooth (phi_c, lbound (v, 3), v)
call smooth (phi_c, lbound (w, 3), w, node='w')

$if ($VERBOSE)
call exit_sub (sub_name)
$endif

end subroutine level_set_smooth_vel
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--smoothes in-place, so be careful
!--uses SOR for laplace equation to acheive smoothing
!--only smoothes region phi < phi0
!--Used to NOT work near boundaries (of grid), but does now as
!--autowrapping of points has been added
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine smooth (phi0, albz, a, node)
use grid_defs, only : grid_t !autowrap_i, autowrap_j
implicit none

real (rp), intent (in) :: phi0
integer, intent (in) :: albz
real (rp), intent (in out) :: a(ld, ny, albz:nz)
character (*), intent (in), optional :: node  !--'u' or 'w'

character (*), parameter :: sub_name = mod_name // '.smooth'

! Moved to level_set_base (renamed to smooth_mode )
!character (*), parameter :: mode = 'xy'  !--'xy', '3d'

integer, parameter :: niter = 5  !--number of SOR iterations

real (rp), parameter :: omega = 1.5_rp  !--SOR parameter

integer :: i, j, k
integer :: s  !--shift for u/w node selection
integer :: nnbr  !--number of neighbors
integer :: iter
integer :: kmin, kmax

integer, pointer, dimension(:) :: autowrap_i, autowrap_j

real (rp) :: phi1
real (rp) :: update

! For autowrapping points: im1 = i-1, ip1 = i+1, etc.
real(rp) :: im1, ip1, jm1, jp1

!---------------------------------------------------------------------
nullify(autowrap_i, autowrap_j)

$if ($VERBOSE)
call enter_sub (sub_name)
$endif

autowrap_i => grid_t % autowrap_i
autowrap_j => grid_t % autowrap_j

if (present (node)) then

  select case (node)
    case ('u'); s = 0
    case ('w'); s = 1
    case default; call error (sub_name, 'invalid node =' // node)
  end select

else

  s = 0  !--u-node by default (assumes phi is on u-node)

end if

select case (smooth_mode)
  case ('xy')
    kmin = 1
    kmax = nz
  case ('3d')
    kmin = 2  !--this is not MPI-enabled
    kmax = nz - 1
    $if($MPI)
    call error (sub_name, 'smooth_mode 3d not MPI compliant')
    $endif
  case default
    call error (sub_name, 'invalid smooth_mode =' // smooth_mode)
end select

do iter = 1, niter
  
  do k = kmin, kmax
    do j = 1, ny
      do i = 1, nx

        if ( ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) .and.  &
             (k == s) ) then
          !--avoid phi(0)
          phi1 = phi(i, j, k)
        else
          phi1 = 0.5_rp * (phi(i, j, k) + phi(i, j, k - s))
        end if

        if (phi1 < phi0) then  !--note its less than

          !  Autowrap boundary points
          im1 = autowrap_i(i-1)
          ip1 = autowrap_i(i+1)
          jm1 = autowrap_j(j-1)
          jp1 = autowrap_j(j+1)

          select case (smooth_mode)
            case ('xy')
              nnbr = 4
              update = (a(im1, j, k) + a(ip1, j, k) +     &
                        a(i, jm1, k) + a(i, jp1, k)) / nnbr
            case ('3d')
              nnbr = 6
              update = (a(im1, j, k) + a(ip1, j, k) +     &
                        a(i, jm1, k) + a(i, jp1, k) +     &
                        a(i, j, k-1) + a(i, j, k+1)) / nnbr
            case default
              call error (sub_name, 'invalid smooth_mode =' // smooth_mode)
          end select

          a(i, j, k) = (1._rp - omega) * a(i, j, k) + omega * update 

        end if

      end do
    end do
  end do

end do

nullify(autowrap_i, autowrap_j)

$if ($VERBOSE)
call exit_sub (sub_name)
$endif

end subroutine smooth

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--this routine relies on user-supplied data about the projected area
!  of the level set object, because it is hard to get exact values from the
!  level set alone
!--may want to put option to append to previous output, if it exists
!--also calculates CL (coeff. of lift)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine level_set_global_CD ()
use param, only : jt, jt_total, dt, L_y, L_z
use immersedbc, only : fx, fy, fz
use sim_param, only : u
implicit none

include 'tecryte.h'

character (*), parameter :: sub_name = mod_name // '.level_set_global_CD'
character (*), parameter :: fCD_out = 'output/global_CD.dat'

integer, parameter :: lun = 99  !--keep open between calls
integer, parameter :: n_calc_CD = 10  !--# t-steps between updates
!integer, parameter :: Ldir = 2
!                      !--lift direction:
!                      !  2 when cyl-axis is z
!                      !  3 when cyl axis is y

!logical, parameter :: DEBUG = .false.

real (rp), parameter :: Ap = 1._rp * 1._rp  !--projected area

logical, save :: file_init = .false.
logical :: opn, exst

real (rp) :: CD, CL
real (rp) :: Uinf   !--velocity scale used in calculation of CD
real (rp) :: fD, fL      !--drag, lift force
real (rp) :: fD_global, fL_global, Uinf_global

!---------------------------------------------------------------------

if (modulo (jt, n_calc_CD) /= 0) return  !--do nothing

fD = -sum (fx(1:nx, :, 1:nz-1)) * dx * dy * dz
     !--(-) since want force ON cylinder
     !--dx*dy*dz is since force is per cell (unit volume)
     !--may want to restrict this sum to points with phi < 0.

select case (Ldir)
  case (2)
    fL = -sum (fy(1:nx, :, 1:nz-1)) * dx * dy * dz
  case (3)
    fL = -sum (fz(1:nx, :, 1:nz-1)) * dx * dy * dz
  case default
    call error (sub_name, 'invalid Ldir =', Ldir)
end select

Uinf = sum (u(1, :, 1:nz-1)) / (ny * (nz - 1))  !--measure at inflow plane

$if ($MPI)

  !--accumulate at coord 0
  call mpi_reduce (fD, fD_global, 1, MPI_RPREC, MPI_SUM,  &
                   rank_of_coord(0), comm, ierr)
  call mpi_reduce (fL, fL_global, 1, MPI_RPREC, MPI_SUM,  &
                   rank_of_coord(0), comm, ierr)
  call mpi_reduce (Uinf, Uinf_global, 1, MPI_RPREC, MPI_SUM,  &
                   rank_of_coord(0), comm, ierr)

  if (coord == 0) Uinf_global = Uinf_global / nproc

$else

  fD_global = fD
  fL_global = fL
  Uinf_global = Uinf

$endif

$if($MPI)
if( coord == 0 ) then
$endif 

  CD = fD_global / (0.5_rp * Ap * Uinf_global**2)
  CL = fL_global / (0.5_rp * Ap * Uinf_global**2)

  if( .not. file_init ) then

    inquire (file=fCD_out, exist=exst, opened=opn)

    !  Check that output is not already opened
    if (opn) call error (sub_name, 'unit', lun, ' is already open')

    if( .not. exst ) then

      call write_tecplot_header_xyline(fCD_out, 'rewind', '"t", "CD", "fD", "CL", "fL", "Uinf"')

    endif

    file_init = .true.

  endif

  call write_real_data(fCD_out, 'append', 'formatted', 6, &
    (/ total_time, CD, fD_global, CL, fL_global, Uinf_global /))

  $if ($DEBUG)
  if (DEBUG) call mesg (sub_name, 'jt_total =', jt_total)
  $endif

$if($MPI)
end if
$endif

end subroutine level_set_global_CD



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine level_set_Cs (delta)
use param, only : vonK, Co, dx, dy, dz, nx, ny, nz, n => nnn 
use sgsmodule, only : Cs_opt2
implicit none

real (rp), intent (in) :: delta

character (*), parameter :: sub_name = mod_name // '.level_set_Cs'

integer :: jx, jy, jz
integer :: jz_min
integer :: node

logical, save :: initialized = .false.

real (rp) :: dmin, z
real (rp) :: phi_tmp

!----------------------------------------------------------------------
$if ($VERBOSE)
call enter_sub (sub_name)
$endif

if (.not. initialized) then

  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then

    jz = 1

    !--u-node at jz = 1
    do jy = 1, ny
      do jx = 1, nx

        if (lbc_mom == 'stress free') then

          dmin = phi(jx, jy, jz)

        else

          z = (jz - 0.5_rp) * dz
          dmin = min (phi(jx, jy, jz), z)

        end if

        !--sorry, this splitting is ugly
        Cs_opt2(jx, jy, jz) = ( Co**(-n) + (delta / vonK / (dmin + z0))**n  &
                              )**(-2._rp / n)

      end do
    end do
        
    jz_min = 2
    
  else

    jz_min = 1

  end if

  !--w-nodes
  do jz = jz_min, nz
    do jy = 1, ny
      do jx = 1, nx

        if (lbc_mom == 'stress free') then

          dmin = 0.5_rp * (phi(jx, jy, jz) + phi(jx, jy, jz - 1))
                           !--MPI: requires phi(k=0)

        else  !--also take wall into account

          if (USE_MPI) then
            z = (coord * (nz - 1) + jz - 1) * dz
          else
            z = (jz - 1) * dz
          end if
          
          phi_tmp = 0.5_rp * (phi(jx, jy, jz) + phi(jx, jy, jz - 1))
                              !--MPI: requires phi(k=0)
          
          dmin = min (phi_tmp, z)  !--min distance to surface

        end if
      
        !--sorry, this splitting is ugly
        Cs_opt2(jx, jy, jz) = ( Co**(-n) + (delta / vonK / (dmin + z0))**n  &
                              )**(-2._rp / n)

      end do
    end do
  end do

  initialized = .true.

end if

$if ($VERBOSE)
call exit_sub (sub_name)
$endif

end subroutine level_set_Cs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--this routine communicates extra the boundary info required for
!  the boundary conditions
!--the thickness of the extra layers is variable, b/c the BCs may
!  change
!--this routine only needed when using MPI
!--this routine is now ridiculously expensive due to all the MPI calls
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
$if ($MPI)
subroutine mpi_sync ()
use sim_param, only : u, v, w, txx, txy, txz, tyy, tyz, tzz
implicit none

character (*), parameter :: sub_name = mod_name // '.mpi_sync'

integer, parameter :: tag = 1000

integer :: datasize
integer :: kstart

logical, save :: phi_synced = .false.

!---------------------------------------------------------------------
$if ($VERBOSE)
call enter_sub (sub_name)
$endif

!--this logic MUST match that in level_set_BC
if (.not. use_log_profile) then

  if (use_extrap_tau_log) then

    !<extrap_tau_log is used>
    call error (sub_name, 'not implemented for use_extrap_tau_log')
    
  else
  
    !--interp_tau: needs sync of u, v, w

    !--this assumes u, v, w already hold 0:nz (bottom 1:nz)
    !--send data down (goes into "top" of process below)
    !--start of send data is k=2, since the corresponds to nz+1 at the
    !  process below the current one
    datasize = ld * ny * nveltop
    kstart = 2
    call mpi_sendrecv (u(1, 1, kstart), datasize, MPI_RPREC, down, tag+1,  &
                       utop(1, 1, 1), datasize, MPI_RPREC, up, tag+1,      &
                       comm, status, ierr)
    call mpi_sendrecv (v(1, 1, kstart), datasize, MPI_RPREC, down, tag+2,  &
                       vtop(1, 1, 1), datasize, MPI_RPREC, up, tag+2,      &
                       comm, status, ierr)
    call mpi_sendrecv (w(1, 1, kstart), datasize, MPI_RPREC, down, tag+3,  &
                       wtop(1, 1, 1), datasize, MPI_RPREC, up, tag+3,      &
                       comm, status, ierr)

    !--send data up (goes into "bottom" of process above)
    !--start of send data is k = nz - 1 - nvelbot, which corresponds to
    !  -nvelbot at the process above the current one
    datasize = ld * ny * nvelbot
    kstart = nz - 1 - nvelbot
    call mpi_sendrecv (u(1, 1, kstart), datasize, MPI_RPREC, up, tag+4,  &
                       ubot(1, 1, 1), datasize, MPI_RPREC, down, tag+4,  &
                       comm, status, ierr)
    call mpi_sendrecv (v(1, 1, kstart), datasize, MPI_RPREC, up, tag+5,  &
                       vbot(1, 1, 1), datasize, MPI_RPREC, down, tag+5,  &
                       comm, status, ierr)
    call mpi_sendrecv (w(1, 1, kstart), datasize, MPI_RPREC, up, tag+6,  &
                       wbot(1, 1, 1), datasize, MPI_RPREC, down, tag+6,  &
                       comm, status, ierr)

    if (use_extrap_tau_simple) then
      !--extrap_tau_simple: needs sync of tij, and phi too (I think)
      !--phi only needs to be "sync"ed once at first call

      if (.not. phi_synced) then

        !--this assumes phi is valid 0:nz
        !--make sure this is consistent with level_set_init

        datasize = ld * ny * nphitop
        kstart = 2
        call mpi_sendrecv (phi(1, 1, kstart), datasize, MPI_RPREC, down,  &
                           tag+7,                                         &
                           phitop(1, 1, 1), datasize, MPI_RPREC, up,      &
                           tag+7, comm, status, ierr)

        datasize = ld * ny * nphibot
        kstart = nz - 1 - nphibot
        call mpi_sendrecv (phi(1, 1, kstart), datasize, MPI_RPREC, up,  &
                           tag+8,                                       &
                           phibot(1, 1, 1), datasize, MPI_RPREC, down,  &
                           tag+8, comm, status, ierr)

      end if

      call mpi_sync_tau ()

    else
      !<extrap_tau>
      call error (sub_name, 'not implemented for use with extrap_tau')
    end if

  end if
  
end if

$if ($VERBOSE)
call exit_sub (sub_name)
$endif

end subroutine mpi_sync
$endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
$if ($MPI)
subroutine mpi_sync_tau ()
use sim_param, only : txx, txy, txz, tyy, tyz, tzz
use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWNUP
implicit none

character (*), parameter :: sub_name = mod_name // '.mpi_sync_tau'

integer, parameter :: tag = 1400

integer :: datasize
integer :: kstart

!---------------------------------------------------------------------

!--tij only valid 1:nz-1
!--first step synchronize tij to make is valid at 0:nz

!--send 1 level down to nz
call mpi_sendrecv (txx(1, 1, 1), ld*ny, MPI_RPREC, down, tag+51,  &
                   txx(1, 1, nz), ld*ny, MPI_RPREC, up, tag+51,   &
                   comm, status, ierr)
call mpi_sendrecv (txy(1, 1, 1), ld*ny, MPI_RPREC, down, tag+52,  &
                   txy(1, 1, nz), ld*ny, MPI_RPREC, up, tag+52,   &
                   comm, status, ierr)
call mpi_sendrecv (txz(1, 1, 1), ld*ny, MPI_RPREC, down, tag+53,  &
                   txz(1, 1, nz), ld*ny, MPI_RPREC, up, tag+53,   &
                   comm, status, ierr)
call mpi_sendrecv (tyy(1, 1, 1), ld*ny, MPI_RPREC, down, tag+54,  &
                   tyy(1, 1, nz), ld*ny, MPI_RPREC, up, tag+54,   &
                   comm, status, ierr)
call mpi_sendrecv (tyz(1, 1, 1), ld*ny, MPI_RPREC, down, tag+55,  &
                   tyz(1, 1, nz), ld*ny, MPI_RPREC, up, tag+55,   &
                   comm, status, ierr)
call mpi_sendrecv (tzz(1, 1, 1), ld*ny, MPI_RPREC, down, tag+56,  &
                         tzz(1, 1, nz), ld*ny, MPI_RPREC, up, tag+56,   &
                         comm, status, ierr)

!--send nz-1 level up to 0 level
call mpi_sendrecv (txx(1, 1, nz-1), ld*ny, MPI_RPREC, up, tag+61,  &
                   txx(1, 1, 0), ld*ny, MPI_RPREC, down, tag+61,   &
                   comm, status, ierr)
call mpi_sendrecv (txy(1, 1, nz-1), ld*ny, MPI_RPREC, up, tag+62,  &
                   txy(1, 1, 0), ld*ny, MPI_RPREC, down, tag+62,   &
                   comm, status, ierr)
call mpi_sendrecv (txz(1, 1, nz-1), ld*ny, MPI_RPREC, up, tag+63,  &
                   txz(1, 1, 0), ld*ny, MPI_RPREC, down, tag+63,   &
                   comm, status, ierr)
call mpi_sendrecv (tyy(1, 1, nz-1), ld*ny, MPI_RPREC, up, tag+64,  &
                   tyy(1, 1, 0), ld*ny, MPI_RPREC, down, tag+64,   &
                   comm, status, ierr)
call mpi_sendrecv (tyz(1, 1, nz-1), ld*ny, MPI_RPREC, up, tag+65,  &
                   tyz(1, 1, 0), ld*ny, MPI_RPREC, down, tag+65,   &
                   comm, status, ierr)
call mpi_sendrecv (tzz(1, 1, nz-1), ld*ny, MPI_RPREC, up, tag+66,  &
                   tzz(1, 1, 0), ld*ny, MPI_RPREC, down, tag+66,   &
                   comm, status, ierr)

!call mpi_sync_real_array( txx, MPI_SYNC_DOWNUP )
!call mpi_sync_real_array( txy, MPI_SYNC_DOWNUP )
!call mpi_sync_real_array( txz, MPI_SYNC_DOWNUP )
!call mpi_sync_real_array( tyy, MPI_SYNC_DOWNUP )
!call mpi_sync_real_array( tyz, MPI_SYNC_DOWNUP )
!call mpi_sync_real_array( tzz, MPI_SYNC_DOWNUP )

!--at this point, tij 0:nz are valid (1:nz at bottom, 0:nz-1 at top)

!--this assumes tij already hold 0:nz
datasize = ld * ny * ntautop
kstart = 2
call mpi_sendrecv (txx(1, 1, kstart), datasize, MPI_RPREC, down, tag+9,  &
                   txxtop(1, 1, 1), datasize, MPI_RPREC, up, tag+9,      &
                   comm, status, ierr)
call mpi_sendrecv (txy(1, 1, kstart), datasize, MPI_RPREC, down,  &
                   tag+10,                                        &
                   txytop(1, 1, 1), datasize, MPI_RPREC, up,      &
                   tag+10, comm, status, ierr)
call mpi_sendrecv (txz(1, 1, kstart), datasize, MPI_RPREC, down,  &
                   tag+11,                                        &
                   txztop(1, 1, 1), datasize, MPI_RPREC, up,      &
                   tag+11, comm, status, ierr)
call mpi_sendrecv (tyy(1, 1, kstart), datasize, MPI_RPREC, down,  &
                   tag+12,                                        &
                   tyytop(1, 1, 1), datasize, MPI_RPREC, up,      &
                   tag+12, comm, status, ierr)
call mpi_sendrecv (tyz(1, 1, kstart), datasize, MPI_RPREC, down,  &
                   tag+13,                                        &
                   tyztop(1, 1, 1), datasize, MPI_RPREC, up,      &
                   tag+13, comm, status, ierr)
call mpi_sendrecv (tzz(1, 1, kstart), datasize, MPI_RPREC, down,  &
                   tag+14,                                        &
                   tzztop(1, 1, 1), datasize, MPI_RPREC, up,      &
                   tag+14, comm, status, ierr)

datasize = ld * ny * ntaubot
kstart = nz - 1 - ntaubot
call mpi_sendrecv (txx(1, 1, kstart), datasize, MPI_RPREC, up,  &
                   tag+15,                                      &
                   txxbot(1, 1, 1), datasize, MPI_RPREC, down,  &
                   tag+15, comm, status, ierr)
call mpi_sendrecv (txy(1, 1, kstart), datasize, MPI_RPREC, up,  &
                   tag+16,                                      &
                   txybot(1, 1, 1), datasize, MPI_RPREC, down,  &
                   tag+16, comm, status, ierr)
call mpi_sendrecv (txz(1, 1, kstart), datasize, MPI_RPREC, up,  &
                   tag+17,                                      &
                   txzbot(1, 1, 1), datasize, MPI_RPREC, down,  &
                   tag+17, comm, status, ierr)
call mpi_sendrecv (tyy(1, 1, kstart), datasize, MPI_RPREC, up,  &
                   tag+18,                                      &
                   tyybot(1, 1, 1), datasize, MPI_RPREC, down,  &
                   tag+18, comm, status, ierr)
call mpi_sendrecv (tyz(1, 1, kstart), datasize, MPI_RPREC, up,  &
                   tag+19,                                      &
                   tyzbot(1, 1, 1), datasize, MPI_RPREC, down,  &
                   tag+19, comm, status, ierr)
call mpi_sendrecv (tzz(1, 1, kstart), datasize, MPI_RPREC, up,  &
                   tag+20,                                      &
                   tzzbot(1, 1, 1), datasize, MPI_RPREC, down,  &
                   tag+20, comm, status, ierr)

end subroutine mpi_sync_tau
$endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--when this is called from sgs_stag:
!  * tij are valid at 1:nz-1 only
!  * u, v, w, are valid at 0:nz (bottom has only 1:nz)
!--upon exit: tij are valid 1:nz-1 only
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine level_set_BC ()
use sim_param, only : txx, txy, txz, tyy, tyz, tzz
implicit none

character (*), parameter :: sub_name = mod_name // '.level_set_BC'

$if ($DEBUG)
logical, parameter :: DEBUG = .false.
$endif

!---------------------------------------------------------------------
$if ($VERBOSE)
call enter_sub (sub_name)
$endif

$if ($MPI)
  call mpi_sync ()
$endif

!--set the phi_cutoff:  all the BC routines should rely on this
if (.not. phi_cutoff_is_set) then

  phi_cutoff = filter_size * 1.1_rp * dx
  !phi_cutoff = filter_size * sqrt (dx**2 + dy**2 + dz**2)
  phi_cutoff_is_set = .true.

end if

if (.not. phi_0_is_set) then

  phi_0 = -100._rp * epsilon (1._rp)
  !phi_0 = 0._rp
  phi_0_is_set = .true.

end if

if (.not. use_log_profile) then  !--skip this if enforce log profile directly

  if (use_extrap_tau_log) then

    call extrap_tau_log ()

  else
  
    $if ($DEBUG)
    if (DEBUG) then
      call DEBUG_write (txx(:, :, 1:nz), 'level_set_BC.a.txx')
      call DEBUG_write (txy(:, :, 1:nz), 'level_set_BC.a.txy')
      call DEBUG_write (txz(:, :, 1:nz), 'level_set_BC.a.txz')
      call DEBUG_write (tyy(:, :, 1:nz), 'level_set_BC.a.tyy')
      call DEBUG_write (tyz(:, :, 1:nz), 'level_set_BC.a.tyz')
      call DEBUG_write (tzz(:, :, 1:nz), 'level_set_BC.a.tzz')
    end if
    $endif

    call interp_tau ()

    $if ($DEBUG)
    if (DEBUG) then
      call DEBUG_write (txx(:, :, 1:nz), 'level_set_BC.b.txx')
      call DEBUG_write (txy(:, :, 1:nz), 'level_set_BC.b.txy')
      call DEBUG_write (txz(:, :, 1:nz), 'level_set_BC.b.txz')
      call DEBUG_write (tyy(:, :, 1:nz), 'level_set_BC.b.tyy')
      call DEBUG_write (tyz(:, :, 1:nz), 'level_set_BC.b.tyz')
      call DEBUG_write (tzz(:, :, 1:nz), 'level_set_BC.b.tzz')
    end if
    $endif
    if (use_extrap_tau_simple) then
      call extrap_tau_simple ()
    else
      call extrap_tau ()
    end if

    $if ($DEBUG)
    if (DEBUG) then
      call DEBUG_write (txx(:, :, 1:nz), 'level_set_BC.c.txx')
      call DEBUG_write (txy(:, :, 1:nz), 'level_set_BC.c.txy')
      call DEBUG_write (txz(:, :, 1:nz), 'level_set_BC.c.txz')
      call DEBUG_write (tyy(:, :, 1:nz), 'level_set_BC.c.tyy')
      call DEBUG_write (tyz(:, :, 1:nz), 'level_set_BC.c.tyz')
      call DEBUG_write (tzz(:, :, 1:nz), 'level_set_BC.c.tzz')
    end if 
    $endif

  end if

end if

if (use_smooth_tau) then
  !--phi_cutoff is set at start of this routine, so no need to check is_set
  call smooth_tau (-phi_cutoff, txx, txy, txz, tyy, tyz, tzz)

  $if ($DEBUG)
  if (DEBUG) then
    call DEBUG_write (txx(:, :, 1:nz), 'level_set_BC.d.txx')
    call DEBUG_write (txy(:, :, 1:nz), 'level_set_BC.d.txy')
    call DEBUG_write (txz(:, :, 1:nz), 'level_set_BC.d.txz')
    call DEBUG_write (tyy(:, :, 1:nz), 'level_set_BC.d.tyy')
    call DEBUG_write (tyz(:, :, 1:nz), 'level_set_BC.d.tyz')
    call DEBUG_write (tzz(:, :, 1:nz), 'level_set_BC.d.tzz')
  end if
  $endif

end if

$if ($VERBOSE)
call exit_sub (sub_name)
$endif

end subroutine level_set_BC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--now, we also need to do extrapolation passes for dead/solid nodes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine extrap_tau ()
use param, only : jt  !--just for debug
use sim_param, only : txx, txy, txz, tyy, tyz, tzz
implicit none

character (*), parameter :: sub_name = mod_name // '.extrap_tau'

character (*), parameter :: fmt1 = '(4(i0,1x))'
integer, parameter :: lun1 = 1

$if ($DEBUG)
logical, parameter :: DEBUG = .false.
$endif
!real (rp), parameter :: phi_0 = 0._rp * z0  !--should be consistent with interp

integer :: i, j, k, m, id
integer :: counter, nlist
integer :: d_2(3)
integer :: pt(nd), s(nd)
integer :: indx(nd, 7), list(nd, 7)
integer :: nxi(nd)

$if ($DEBUG)
logical, save :: first_call = .true.
$endif

real (rp) :: phi_c, phiw, phiw_m
real (rp) :: sphi
real (rp) :: normw(nd)

!---------------------------------------------------------------------
$if ($VERBOSE)
call enter_sub (sub_name)
$endif

if (phi_cutoff_is_set) then
  phi_c = phi_cutoff
else
  call error (sub_name, 'trying to use uninitialized phi_cutoff')
end if

if (.not. phi_0_is_set) then
  call error (sub_name, 'trying to use uninitialized phi_0')
end if

nxi = (/ nx, ny, nz /)

$if ($DEBUG)
if (DEBUG) then
  if (first_call) then
    open (lun1, file='level_set_extrap-u.dat')
    write (lun1, *) 'variables = "i" "j" "k" "nlist"'
    write (lun1, *) 'zone, f=point, i=', nx-2, ', j=', ny-2, ', k=', nz-2
  end if
end if
$endif
!--u-nodes
do k = 2, nz-1
  do j = 2, ny-1
    do i = 2, nx-1

      $if ($DEBUG)
      if (DEBUG .and. first_call) nlist = 0
      $endif

      if ((phi(i, j, k) < phi_0) .and. (phi(i, j, k) >= -phi_c)) then

        pt = (/ i, j, k /)

        !--try to locate three closest fluid nodes and fit plane
        s = nint (sign (1._rp, norm(:, i, j, k)))
        counter = 0
        list = iBOGUS  !--set to bogus value

        do m = 1, 7  !--7 other corners apart from this one

          d_2(1) = modulo (m, 2)  !--this is easier than bit routines
          d_2(2) = modulo (m/2, 2)
          d_2(3) = modulo (m/4, 2)

          indx(:, m) = pt + d_2 * s

          $if ($DEBUG)
          if (DEBUG .and. jt >= 615) then
            call mesg (sub_name, '(i, j, k)=', (/ i, j, k /))
            call mesg (sub_name, 'pt=', pt)
            call mesg (sub_name, 'd_2=', d_2)
            call mesg (sub_name, 's=', s)
            call mesg (sub_name, 'indx(:,', m, ') =', indx(:, m))
          end if
          $endif

          do id = 1, nd
            if ( (indx(id, m) < 1) .or. (indx(id, m) > nxi(id)) ) then
              write (msg, '(3(a,i0),a)') 'indx(', id, ',', m, ')=',     &
                                         indx(id, m), ' is out of bounds'
              call error (sub_name, msg)
            end if
          end do

          sphi = sign (1._rp, phi(indx(1, m), indx(2, m), indx(3, m)))

          if (sphi > 0._rp) then
            counter = counter + 1
            list(:, counter) = indx(:, m)  !--store positive-phi points
          end if
          
        end do

        nlist = counter

        !if ((counter > 0) .and. (counter < 3)) then
        !  !--try to find more points
        !  
        !end do
       
        !--stress extrapolation
        select case (nlist)
          case (0)

            !--do nothing

          case (1)

            !--copy stress
            txx(i, j, k) = txx(list(1, 1), list(2, 1), list(3, 1))
            txy(i, j, k) = txy(list(1, 1), list(2, 1), list(3, 1))
            tyy(i, j, k) = tyy(list(1, 1), list(2, 1), list(3, 1))
            
          case (2)

            !--should try to find more fluid points, this is only "for now"
            txx(i, j, k) = (txx(list(1, 1), list(2, 1), list(3, 1)) +      &
                            txx(list(1, 2), list(2, 2), list(3, 2))) / 2._rp
            txy(i, j, k) = (txy(list(1, 1), list(2, 1), list(3, 1)) +      &
                            txy(list(1, 2), list(2, 2), list(3, 2))) / 2._rp
            tyy(i, j, k) = (tyy(list(1, 1), list(2, 1), list(3, 1)) +      &
                            tyy(list(1, 2), list(2, 2), list(3, 2))) / 2._rp
            tzz(i, j, k) = (tzz(list(1, 1), list(2, 1), list(3, 1)) +      &
                            tzz(list(1, 2), list(2, 2), list(3, 2))) / 2._rp

          case (3:7)

            call fit3 (pt, nlist, list, txx) 
            call fit3 (pt, nlist, list, txy)
            call fit3 (pt, nlist, list, tyy)
            call fit3 (pt, nlist, list, tzz)
            
          case default
            call error (sub_name, 'invalid counter value')
        end select

      end if
 
      $if ($DEBUG)
      if (DEBUG) then
        if (jt >= 615) then
          call mesg (sub_name, 'pt =', pt)
          call mesg (sub_name, 'txx =', txx(pt(1), pt(2), pt(3)))
          call mesg (sub_name, 'txy =', txy(pt(1), pt(2), pt(3)))
          call mesg (sub_name, 'tyy =', tyy(pt(1), pt(2), pt(3)))
          call mesg (sub_name, 'tzz =', tzz(pt(1), pt(2), pt(3)))
        end if
      end if
      $endif
 
      $if ($DEBUG)
      if (DEBUG .and. first_call) then
        write (lun1, fmt1) i, j, k, nlist
      end if
      $endif

    end do
  end do
end do

$if ($DEBUG)
if (DEBUG .and. first_call) then
  close (lun1)
  open (lun1, file='level_set_extrap-w.dat')
  write (lun1, *) 'variables = "i" "j" "k" "nlist"'
  write (lun1, *) 'zone, f=point, i=', nx-2, ', j=', ny-2, ', k=', nz-2
end if
$endif

$if ($DEBUG)
if (DEBUG) call mesg (sub_name, 'done extrapolation (u)')
$endif

!--w-nodes
do k = 2, nz-1
  do j = 2, ny-1
    do i = 2, nx-1

      $if ($DEBUG)
      if (DEBUG .and. first_call) nlist = 0
      $endif
      phiw = (phi(i, j, k) + phi(i, j, k-1)) / 2._rp
      
      if ((phiw < phi_0) .and. (phiw >= -phi_c)) then

        pt = (/ i, j, k /)

        !--try to locate three closest fluid nodes and fit plane
        normw = (norm(:, i, j, k) + norm(:, i, j, k-1)) / 2._rp
        normw = normw / mag (normw)
        
        s = nint (sign (1._rp, normw))
        counter = 0
        list = -1  !--set to bogus value

        do m = 1, 7  !--7 other corners apart from this one

          d_2(1) = modulo (m, 2)  !--this is easier than bit routines
          d_2(2) = modulo (m/2, 2)
          d_2(3) = modulo (m/4, 2)

          indx(:, m) = pt + d_2 * s

          phiw_m = (phi(indx(1, m), indx(2, m), indx(3, m)) +       &
                    phi(indx(1, m), indx(2, m), indx(3, m)-1)) / 2._rp
          sphi = sign (1._rp, phiw_m)

          if (sphi > 0._rp) then
            counter = counter + 1
            list(:, counter) = indx(:, m)  !--store positive-phi points
          end if
          
        end do

        nlist = counter

        !if ((counter > 0) .and. (counter < 3)) then
        !  !--try to find more points
        !  
        !end do
       
        !--stress extrapolation
        select case (nlist)
          case (0)

            !--do nothing

          case (1)

            !--copy stress
            txz(i, j, k) = txz(list(1, 1), list(2, 1), list(3, 1))
            tyz(i, j, k) = tyz(list(1, 1), list(2, 1), list(3, 1))
            
          case (2)

            !--should try to find more fluid points, this is only "for now"
            txz(i, j, k) = (txz(list(1, 1), list(2, 1), list(3, 1)) +      &
                            txz(list(1, 2), list(2, 2), list(3, 2))) / 2._rp
            tyz(i, j, k) = (tyz(list(1, 1), list(2, 1), list(3, 1)) +      &
                            tyz(list(1, 2), list(2, 2), list(3, 2))) / 2._rp

          case (3:7)

            call fit3 (pt, nlist, list, txz) 
            call fit3 (pt, nlist, list, tyz)
            
          case default
            call error (sub_name, 'invalid counter value')
        end select
       
      end if
      
      $if ($DEBUG)
      if (DEBUG) then
        if (jt >= 615) then
          call mesg (sub_name, 'pt =', pt)
          call mesg (sub_name, 'txz =', txz(pt(1), pt(2), pt(3)))
          call mesg (sub_name, 'tyz =', tyz(pt(1), pt(2), pt(3)))
        end if
      end if
      $endif

      $if ($DEBUG)
      if (DEBUG .and. first_call) then
        write (lun1, fmt1) i, j, k, nlist
      end if
      $endif
    end do
  end do
end do

$if ($DEBUG)
if (DEBUG) call mesg (sub_name, 'done extrapolation (w)')
$endif

$if ($DEBUG)
if (DEBUG .and. first_call) then
  close (lun1)
  first_call = .false.
end if
$endif

$if ($VERBOSE)
call exit_sub (sub_name)
$endif

end subroutine extrap_tau

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine interp_tau ()
use param, only : jt  !--in addition to stuff above
use sim_param, only : u, v, w, txx, txy, txz, tyy, tyz, tzz,             &
                      dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
implicit none

character (*), parameter :: sub_name = mod_name // '.interp_tau'
character (*), parameter :: fprefix = 'output/interp_tau.'
character (*), parameter :: fmta3r = '(a,3(es12.5,1x))'
character (*), parameter :: fmta3i = '(a,3(i0,1x))'

integer, parameter :: noutput = 200
integer, parameter :: lun = 1

!logical, parameter :: DEBUG = .false.
logical, parameter :: use_output = .false.

real (rp), parameter :: eps = 100._rp * epsilon (0._rp)
!real (rp), parameter :: phi_0 = 0._rp * z0  !--this is adjustable

character (128) :: fname

integer :: i, j, k
integer :: kmin, kmax
!--experiment for skipping unsed points
integer, save :: imn = 1, imx = nx, jmn = 1, jmx = ny
integer :: imn_used, imx_used, jmn_used, jmx_used

logical :: output
logical :: exst, opn

real (rp) :: kappa
real (rp) :: phi_c, phix
real (rp) :: tau
real (rp) :: n_hat(nd)
real (rp) :: vel(nd), vel_t(nd)
real (rp) :: x_hat(nd), y_hat(nd), z_hat(nd)
real (rp) :: x(nd), xv(nd)

!---------------------------------------------------------------------
$if ($VERBOSE)
call enter_sub (sub_name)
$endif

if (use_output) then

  if (modulo (jt, noutput) == 0) then

    output = .true.

    write (fname, '(a,i6.6,a)') trim (fprefix), jt, '.dat'

    inquire (unit=lun, exist=exst, opened=opn)

    if (exst .and. (.not. opn)) then
      open (lun, file=fname)
    else
      call error (sub_name, 'problem opening file')
    end if

  else
    output = .false.
  end if

else

  output = .false.

end if

if (phi_cutoff_is_set) then
  phi_c = phi_cutoff
else
  call error (sub_name, 'trying to use uninitialized phi_cutoff')
end if

if (.not. phi_0_is_set) then
  call error (sub_name, 'trying to use uninitialized phi_0')
end if

if (output) then
  write (lun, *) 'phi_c = ', phi_c
  write (lun, *) 'phi_0 = ', phi_0
end if

kappa = vonk

!--initial values for _used variables
imn_used = nx
imx_used = 1
jmn_used = ny
jmx_used = 1

if (output) write (lun, *) 'u-node pass'

!--u-node pass
do k = 1, nz-1
  do j = jmn, jmx
    do i = imn, imx

      phix = phi(i, j, k)

      if ((phix >= phi_0) .and. (phix <= phi_c)) then

        imn_used = min (imn_used, i)
        imx_used = max (imx_used, i)
        jmn_used = min (jmn_used, j)
        jmx_used = max (jmx_used, j)

        x = (/ (i - 1) * dx, (j - 1) * dy, (k - 0.5_rp) * dz /)

        n_hat = norm(:, i, j, k)

        !--make sure norm is well-defined
        if (mag (n_hat) < eps) then
          write (msg, *) 'n_hat has 0 magnitude (u) at ',              &
                         (/ i, j, k /), ' [coord =', coord, ']', n_l,  &
                         'n_hat = ', n_hat
          call error (sub_name, msg)
        end if

        if (physBC) then

          !--determine velocity vector at point with phi ~ phi_c
          xv = x + n_hat * (phi_c - phix)

          call interp_vel (xv, vel)

        else
        
          !--determine velocity vector here: beware w may include solid pt
          vel = (/ u(i, j, k), v(i, j, k),               &
                   0.5_rp * (w(i, j, k) + w(i, j, k+1)) /)

        end if
        
        vel_t = vel - dot_product (vel, n_hat) * n_hat  !--tangential part

        if (mag (vel_t) < eps) then
          
          txx(i, j, k) = 0._rp
          txy(i, j, k) = 0._rp
          tyy(i, j, k) = 0._rp
          tzz(i, j, k) = 0._rp
            
        else
         
          !--local coordinate system
          x_hat = vel_t / mag (vel_t)
          y_hat = cross_product (n_hat, x_hat)
          z_hat = n_hat

          if (physBC) then
            tau = -(kappa * mag (vel_t) / log (1._rp + phi_c / z0))**2
          else
            tau = -(kappa * mag (vel_t) / log (1._rp + phix / z0))**2
          end if

          !--special case: only nonzero t_{i'j'} is t_{1'3'} & t_{3'1'}
          txx(i, j, k) = (x_hat(1) * z_hat(1) + z_hat(1) * x_hat(1)) * tau
          txy(i, j, k) = (x_hat(1) * z_hat(2) + z_hat(1) * x_hat(2)) * tau
          tyy(i, j, k) = (x_hat(2) * z_hat(2) + z_hat(2) * x_hat(2)) * tau
          tzz(i, j, k) = (x_hat(3) * z_hat(3) + z_hat(3) * x_hat(3)) * tau

          if (use_modify_dutdn) then
            $if ($MPI)
              call error (sub_name, 'modify_dutdn not MPI-enabled')
            $else
              call modify_dutdn (i, j, k, tau, phix, x_hat, y_hat, z_hat)
            $endif
          end if

        end if

        if (output) then
          write (lun, fmta3i) 'i, j, k = ', i, j, k
          write (lun, *) 'phix = ', phix
          write (lun, fmta3r) 'x = ', x
          write (lun, fmta3r) 'n_hat = ', n_hat
          if (physBC) write (lun, fmta3r) 'xv = ', xv
          write (lun, fmta3r) 'u,v,w @ i  , j  , k : ', u(i, j, k),  &
                              v(i, j, k), w(i, j, k)
          write (lun, fmta3r) 'u,v,w @ i+1, j  , k : ', u(i+1, j, k),  &
                              v(i+1, j, k), w(i+1, j, k)
          write (lun, fmta3r) 'u,v,w @ i-1, j  , k : ', u(i-1, j, k),  &
                              v(i-1, j, k), w(i-1, j, k)
          write (lun, fmta3r) 'u,v,w @ i  , j+1, k : ', u(i, j+1, k),  &
                              v(i, j+1, k), w(i, j+1, k)
          write (lun, fmta3r) 'u,v,w @ i  , j-1, k : ', u(i, j-1, k),  &
                              v(i, j-1, k), w(i, j-1, k)
          write (lun, fmta3r) 'vel = ', vel
          write (lun, fmta3r) 'vel_t = ', vel_t
          if (mag (vel_t) < eps) then
            write (lun, '(a)') 'mag (vel_t) < eps branch taken'
          else
            write (lun, fmta3r) 'x_hat = ', x_hat
            write (lun, fmta3r) 'y_hat = ', y_hat
            write (lun, fmta3r) 'z_hat = ', z_hat
            write (lun, *) 'tau = ', tau
            write (lun, *) 'txx = ', txx(i, j, k)
            write (lun, *) 'txy = ', txy(i, j, k)
            write (lun, *) 'tyy = ', tyy(i, j, k)
            write (lun, *) 'tzz = ', tzz(i, j, k)
            if (use_modify_dutdn) then
              write (lun, *) 'dudx = ', dudx(i, j, k)
              write (lun, *) 'dudy = ', dudy(i, j, k)
              write (lun, *) 'dudz = ', dudz(i, j, k)
              write (lun, *) 'dvdx = ', dvdx(i, j, k)
              write (lun, *) 'dvdy = ', dvdy(i, j, k)
              write (lun, *) 'dvdz = ', dvdz(i, j, k)
              write (lun, *) 'dwdx = ', dwdx(i, j, k)
              write (lun, *) 'dwdy = ', dwdy(i, j, k)
              write (lun, *) 'dwdz = ', dwdz(i, j, k)
            end if
          end if

        end if

      end if

    end do
  end do
end do

$if ($DEBUG)
if (DEBUG) call mesg (sub_name, 'done interpolation (u)')
$endif

if (output) write (lun, *) 'w-node pass'

!--w-node pass

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
  kmin = 2
else
  kmin = 1
end if

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc - 1)) then
  kmax = nz - 1 !--not certain this is what we want here, maybe nz
else
  kmax = nz - 1
end if

do k = kmin, kmax
  do j = jmn, jmx
    do i = imn, imx

      phix = 0.5_rp * (phi(i, j, k) + phi(i, j, k-1)) !--put phi on w-node

      if ((phix >= phi_0) .and. (phix <= phi_c)) then

        imn_used = min (imn_used, i)
        imx_used = max (imx_used, i)
        jmn_used = min (jmn_used, j)
        jmx_used = max (jmx_used, j)

        x = (/ (i - 1) * dx, (j - 1) * dy, (k - 1) * dz /)
        
        !--make sure norms are well defined
        if ((mag (norm(:, i, j, k)) < eps) .or.  &
            (mag (norm(:, i, j, k-1)) < eps)) then
          write (msg, *) 'norm has 0 magnitude (w)', n_l,  &
                         '(i,j,k) =', i, j, k,       n_l,  &
                         'coord =', coord,           n_l,  &
                         'x =', x
          call error (sub_name, msg)
        end if

        !--interpolate norm to w-nodes
        n_hat = 0.5_rp * (norm(:, i, j, k) + norm(:, i, j, k-1))
        n_hat = n_hat / mag (n_hat)

        if (physBC) then
        
          !--determine velocity vector at point with phi ~ phi_c
          xv = x + n_hat * (phi_c - phix)

          call interp_vel (xv, vel)

        else

          !--get velocity at w-node
          vel = (/ 0.5_rp * (u(i, j, k) + u(i, j, k-1)),             &
                   0.5_rp * (v(i, j, k) + v(i, j, k-1)), w(i, j, k) /)

        end if

        vel_t = vel - dot_product (vel, n_hat) * n_hat

        if (mag (vel_t) < eps) then

          txz(i, j, k) = 0._rp
          tyz(i, j, k) = 0._rp
          
        else

          x_hat = vel_t / mag (vel_t)
          y_hat = cross_product (n_hat, x_hat)
          z_hat = n_hat

          if (physBC) then
            tau = -(kappa * mag (vel_t) / log (1._rp + phi_c / z0))**2
          else
            tau = -(kappa * mag (vel_t) / log (1._rp + phix / z0))**2
          end if
          
          txz(i, j, k) = (x_hat(1) * z_hat(3) + z_hat(1) * x_hat(3)) * tau
          tyz(i, j, k) = (x_hat(2) * z_hat(3) + z_hat(2) * x_hat(3)) * tau
          
          if (use_modify_dutdn) then
            $if ($MPI)
              call error (sub_name, 'modify_dutdn not MPI enabled')
            $else
              call modify_dutdn (i, j, k, tau, phix, x_hat, y_hat, z_hat, 'w')
            $endif
          end if

        end if

        if (output) then
          write (lun, fmta3i) 'i, j, k = ', i, j, k
          write (lun, *) 'phix = ', phix
          write (lun, fmta3r) 'x = ', x
          write (lun, fmta3r) 'n_hat = ', n_hat
          if (physBC) write (lun, fmta3r) 'xv = ', xv
          write (lun, fmta3r) 'u,v,w @ i  , j  , k : ', u(i, j, k),  &
                              v(i, j, k), w(i, j, k)
          write (lun, fmta3r) 'u,v,w @ i+1, j  , k : ', u(i+1, j, k),  &
                              v(i+1, j, k), w(i+1, j, k)
          write (lun, fmta3r) 'u,v,w @ i-1, j  , k : ', u(i-1, j, k),  &
                              v(i-1, j, k), w(i-1, j, k)
          write (lun, fmta3r) 'u,v,w @ i  , j+1, k : ', u(i, j+1, k),  &
                              v(i, j+1, k), w(i, j+1, k)
          write (lun, fmta3r) 'u,v,w @ i  , j-1, k : ', u(i, j-1, k),  &
                              v(i, j-1, k), w(i, j-1, k)
          write (lun, fmta3r) 'vel = ', vel
          write (lun, fmta3r) 'vel_t = ', vel_t
          if (mag (vel_t) < eps) then
            write (lun, '(a)') 'mag (vel_t) < eps branch taken'
          else
            write (lun, fmta3r) 'x_hat = ', x_hat
            write (lun, fmta3r) 'y_hat = ', y_hat
            write (lun, fmta3r) 'z_hat = ', z_hat
            write (lun, *) 'tau = ', tau
            write (lun, *) 'txz = ', txz(i, j, k)
            write (lun, *) 'tyz = ', tyz(i, j, k)
            if (use_modify_dutdn) then
              write (lun, *) 'dudx = ', dudx(i, j, k)
              write (lun, *) 'dudy = ', dudy(i, j, k)
              write (lun, *) 'dudz = ', dudz(i, j, k)
              write (lun, *) 'dvdx = ', dvdx(i, j, k)
              write (lun, *) 'dvdy = ', dvdy(i, j, k)
              write (lun, *) 'dvdz = ', dvdz(i, j, k)
              write (lun, *) 'dwdx = ', dwdx(i, j, k)
              write (lun, *) 'dwdy = ', dwdy(i, j, k)
              write (lun, *) 'dwdz = ', dwdz(i, j, k)
            end if
          end if
        end if

      end if

    end do
  end do
end do

!--to be used at next time step
imn = imn_used
imx = imx_used
jmn = jmn_used
jmx = jmx_used

if (output) then
  write (lun, *) ' '
  close (lun)
end if

$if ($MPI)
  !--resync tij, since we have altered it (this is annoying)
  call mpi_sync_tau ()
$endif

$if ($VERBOSE)
call exit_sub (sub_name)
$endif

end subroutine interp_tau

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--nl must be 3 or larger here
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fit3 (p, nl, l, t)
use linear_simple, only : solve_linear
implicit none

integer, intent (in) :: p(nd)
integer, intent (in) :: nl
integer, intent (in) :: l(nd, nl)
real (rp), intent (in out) :: t(ld, ny, nz)

character (*), parameter :: sub_name = mod_name // '.fit3'

$if ($DEBUG)
logical, parameter :: DEBUG = .false.
$endif

integer :: counter
integer :: m, q
integer :: l_plane (nd, 3)  !--list that contains our 3 points for plane
integer :: indx(nd-1, 3)
integer :: dir, d

logical :: failed

real (rp) :: tmp
real (rp) :: A(3, 3), coeff(3), t_known(3)
real (rp) :: dxi(nd-1)

!---------------------------------------------------------------------
$if ($VERBOSE)
call enter_sub (sub_name)
$endif
if (nl < 3) call error (sub_name, 'nl should be >= 3')

counter = 0
dir = -1

do d = 1, nd

  counter = 0
  
  do m = 1, nl
    if (l(d, m) == p(d)) counter = counter + 1
  end do
  
  if (counter >= 3) then
    dir = d
    exit  !--d-loop
  end if
  
end do

!--experiment
!failed = .true.

if (dir == -1) then
  failed = .true.
else
  failed = .false.
end if

if (.not. failed) then

  !--fill l_plane
  l_plane = iBOGUS
  q = 0
  do m = 1, nl
  
    if (l(dir, m) == p(dir)) then

      q = q + 1
      l_plane(:, q) = l(:, m)

      if (q == 3) exit

     end if
  end do

  !--need to exclude dir from the calculation now
  dxi = pack ((/ dx, dy, dz /), mask=(/ 1, 2, 3 /) /= dir)

  do q = 1, 3
    indx(:, q) = pack (l_plane(:, q) - p(:), mask=(/ 1, 2, 3 /) /= dir)
    t_known(q) = t(l_plane(1, q), l_plane(2, q), l_plane(3, q))
  end do

  !--setup 3x3 matrix problem of the form
  !  / 1 x1 y1 \ / a \   / t_known(1) \
  !  | 1 x2 y2 | | b | = | t_known(2) |
  !  \ 1 x3 y3 / \ c /   \ t_known(3) /
  A(:, 1) = 1._rp
  A(:, 2) = indx(1, :) * dxi(1)  !--coords relative to p
  A(:, 3) = indx(2, :) * dxi(2)  !--coords relative to p

  call solve_linear (A, t_known, coeff)

  $if ($DEBUG)
  if (DEBUG) then
    if (coeff(1) > 1.e7_rp) then
      call mesg (sub_name, 'A(1, :) =', A(1, :))
      call mesg (sub_name, 'A(2, :) =', A(2, :))
      call mesg (sub_name, 'A(3, :) =', A(3, :))
      call mesg (sub_name, 't_known =', t_known)
      call mesg (sub_name, 'coeff =', coeff)
      call error (sub_name, 'coeff(1) too big')
    end if
  end if
  $endif

  !--we only need coeff(1)
  t(p(1), p(2), p(3)) = coeff(1)
  
else  !--just average all the points, for now

  tmp = 0._rp

  do m = 1, nl
    tmp = tmp + t(l(1, m), l(2, m), l(3, m))
  end do

  tmp = tmp / m
  
  t(p(1), p(2), p(3)) = tmp

end if

$if ($DEBUG)
if (DEBUG) call mesg (sub_name, 't(p) =', t(p(1), p(2), p(3)))
$endif

$if ($VERBOSE)
call exit_sub (sub_name)
$endif

end subroutine fit3

!**********************************************************************
subroutine level_set_forcing ()
!**********************************************************************
! 
! Set fx, fy, fz at 1:nz-1
!
use param, only : tadv1, dt, BOGUS, dx  !--in addition to param vars above
use sim_param
use immersedbc, only : fx, fy, fz

implicit none

character (*), parameter :: sub_name = mod_name // '.level_set_forcing'

$if ($DEBUG)
logical, parameter :: DEBUG = .false.
$endif

integer :: i, j, k
integer :: k_min

real (rp) :: Rx, Ry, Rz

!---------------------------------------------------------------------
$if ($VERBOSE)
call enter_sub (sub_name)
$endif

!--this is experimental
if (vel_BC) then
  if (use_enforce_un) call enforce_un ()
  if (use_log_profile) call enforce_log_profile ()
end if

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then

  k = 1
  do j = 1, ny
    do i = 1, nx

      call level_set_force_xy()

      !--no w-node part here: the velocity should be zero b/c of BCs

    end do
  end do

  k_min = 2

else

  k_min = 1

end if

do k = k_min, nz - 1
  do j = 1, ny
    do i = 1, nx

      call level_set_force_xy()
      call level_set_force_z()

    end do
  end do
end do

fx(:, :, nz) = BOGUS
fy(:, :, nz) = BOGUS

!Setting 0 at physical top boundary
$if($MPI)
  if( coord == nproc - 1 ) then
    fz(:, :, nz) = 0._rprec
  else
    fz(:, :, nz) = BOGUS
  endif
$else
  fz(:,:,nz) = 0._rprec
$endif

$if ($DEBUG)
if (DEBUG) then
  call DEBUG_write (u(:, :, 1:nz), 'level_set_forcing.u')
  call DEBUG_write (v(:, :, 1:nz), 'level_set_forcing.v')
  call DEBUG_write (w(:, :, 1:nz), 'level_set_forcing.w')
  call DEBUG_write (fx(:, :, 1:nz), 'level_set_forcing.fx')
  call DEBUG_write (fy(:, :, 1:nz), 'level_set_forcing.fy')
  call DEBUG_write (fz(:, :, 1:nz), 'level_set_forcing.fz')
end if
$endif

$if ($VERBOSE)
call exit_sub (sub_name)
$endif

return

contains

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine level_set_force_xy()
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

if (phi(i, j, k) <= 0._rp) then  !--uv-nodes
  
  $if($PC_SCHEME_0)    
  ! Original PC
  Rx = -tadv1 * dpdx(i, j, k)
  Ry = -tadv1 * dpdy(i, j, k)        
  fx(i, j, k) = (-u(i, j, k)/dt - Rx) 
  fy(i, j, k) = (-v(i, j, k)/dt - Ry)

  $elseif($PC_SCHEME_1 or $PC_SCHEME_3)
  ! Updated PC
  fx(i,j,k) = -(u(i,j,k)/dt + tadv1 * RHSx(i, j, k) + tadv2 * RHSx_f(i,j,k) - dpdx_f(i,j,k))
  fy(i,j,k) = -(v(i,j,k)/dt + tadv1 * RHSy(i, j, k) + tadv2 * RHSy_f(i,j,k) - dpdy_f(i,j,k))

  $elseif($PC_SCHEME_2)
  ! Updated PC-2
  fx(i,j,k) = -(u(i,j,k)/dt + tadv1 * RHSx(i, j, k) + tadv2 * RHSx_f(i,j,k))
  fy(i,j,k) = -(v(i,j,k)/dt + tadv1 * RHSy(i, j, k) + tadv2 * RHSy_f(i,j,k))

  $else

  call error(sub_name,'Makefile pressure correction scheme not specified properly')

  $endif

  !  Commented since not compliant with correct pressure scheme (JSG)
!else if (vel_BC) then
       
  !call error(sub_name,'Disabled this feature until forcing updated to be based on u*')

  ! forces after pressure update
  !Rx = -tadv1 * dpdx(i, j, k)
  !Ry = -tadv1 * dpdy(i, j, k)

  !if (udes(i, j, k) < huge (1._rp) / 2) then
    !fx(i, j, k) = ((udes(i, j, k) - u(i, j, k)) / dt - Rx)
  !end if

  !if (vdes(i, j, k) < huge (1._rp) / 2) then
    !fy(i, j, k) = ((vdes(i, j, k) - v(i, j, k)) / dt - Ry)
  !end if

end if

return
end subroutine level_set_force_xy

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine level_set_force_z()
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none


! Original FV
if (phi(i,j,k) + phi(i,j,k-1) <= 0._rp) then  !--w-nodes

  $if($PC_SCHEME_0)
  ! Original PC
  Rz = -tadv1 * dpdz(i, j, k)
  fz(i, j, k) = (-w(i, j, k)/dt - Rz)

  $elseif($PC_SCHEME_1 or $PC_SCHEME_3)
  ! Updated PC
  fz(i,j,k) = -(w(i,j,k)/dt + tadv1 * RHSz(i, j, k) + tadv2 * RHSz_f(i,j,k) - dpdz_f(i,j,k))

  $elseif($PC_SCHEME_2)
  ! Update PC-2
  fz(i,j,k) = -(w(i,j,k)/dt + tadv1 * RHSz(i, j, k) + tadv2 * RHSz_f(i,j,k))

  $else

  call error(sub_name,'Makefile pressure correction scheme not specified properly')

  $endif

!  Commented since not compliant with correct pressure scheme (JSG)
!else if (vel_BC) then
  !call error(sub_name,'Disabled this feature until forcing updated to be based on u*')

  !Rz = -tadv1 * dpdz(i, j, k)

  !if (wdes(i, j, k) < huge (1._rp) / 2._rp) then
    !fz(i, j, k) = ((wdes(i, j, k) - w(i, j, k)) / dt - Rz)
  !end if

end if

return
end subroutine level_set_force_z

end subroutine level_set_forcing

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--calculated centered differences, excetp 1-sided differences at top and
!--bottom edges
!--this assumes f(:, :, nz) is valid
!--will insert BOGUS at nz-level, unless its the top process
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real(rp) function safe_cd (i, j, k, d, f)
use param, only : dx, dy, dz  !--in addition to those above
use grid_defs, only : grid_t ! autowrap_i, autowrap_j
implicit none

integer, intent (in) :: i, j, k
integer, intent (in) :: d  !--d is dimension to difference along
real (rp), intent (in) :: f(ld, ny, $lbz:nz)

character (*), parameter :: sub_name = mod_name // '.safe_cd'

integer :: n

integer, pointer, dimension(:) :: autowrap_i, autowrap_j

real (rp) :: delta

nullify( autowrap_i, autowrap_j )

autowrap_i => grid_t % autowrap_i  
autowrap_j => grid_t % autowrap_j

!---------------------------------------------------------------------

select case (d)
  case (1)
    delta = dx

    !  Commented (JSG)
    !n = nx
    !if (i == 1) then

    !  safe_cd = ( f(i + 1, j, k) - f(i, j, k) ) / delta

    !else if (i == n) then

    !  safe_cd = ( f(i, j, k) - f(i - 1, j, k) ) / delta

    !else

    !  safe_cd = ( f(i + 1, j, k) - f(i - 1, j, k) ) / (2._rp * delta)

    !end if

    !  Using autowrap to take care of edges
    safe_cd = ( f(autowrap_i(i + 1), j, k) - f(autowrap_i(i - 1), j, k) ) / (2._rp * delta)
    
  case (2)
 
    delta = dy
    !  Commented (JSG)
    ! n = ny
    !if (j == 1) then

    !  safe_cd = ( f(i, j + 1, k) - f(i, j, k) ) / delta

    !else if (j == n) then

    !  safe_cd = ( f(i, j, k) - f(i, j - 1, k) ) / delta

    !else

    !  safe_cd = ( f(i, j + 1, k) - f(i, j - 1, k) ) /  (2._rp * delta)

    !end if

    !  Using autowrap to take care of edges
    safe_cd = ( f(i, autowrap_j(j + 1), k) - f(i, autowrap_j(j - 1), k) ) /  (2._rp * delta)

  case (3)
    n = nz
    delta = dz

    if ( ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0))  &
         .and. (k == 1) ) then

      safe_cd = ( f(i, j, k + 1) - f(i, j, k) ) / delta

    else if (k == n) then  !--this should not happen anymore
                           !--this routine only called with 1:nz-1
      call error (sub_name, 'called with k == nz')
      !if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc - 1)) then
      !  safe_cd = ( f(i, j, k) - f(i, j, k - 1) ) / delta
      !else
      !  safe_cd = BOGUS
      !end if

    else

      safe_cd = ( f(i, j, k + 1) - f(i, j, k - 1) ) / (2._rp * delta)

    end if

  case default
    call error (sub_name, 'invalid d =', d)
end select

nullify( autowrap_i, autowrap_j )

end function safe_cd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--try 1-sided derivatives to estimate normal
!--to be used when safe_cd gives indeterminate results
!--ierr < 0 upon failure to calc normal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fix_norm (i, j, k, eps, ntmp, ierr)
implicit none

integer :: i, j, k
real (rp), intent (in) :: eps
real (rp), intent (out) :: ntmp(nd)
integer, intent (out), optional :: ierr !--error code: < 0 upon failure

character (*), parameter :: sub_name = mod_name // '.fix_norm'

integer :: i1, j1, k1
integer :: ii, jj, kk

!---------------------------------------------------------------------

if (present (ierr)) ierr = 0

do k1 = 1, -1, -2

  !if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
  !  if ((k1 == -1) .and. (k == 1)) exit  !--same as cycle here
  !end if

  do j1 = 1, -1, -2
    do i1 = 1, -1, -2

      ii = i + i1
      jj = j + j1
      kk = k + k1

      if ((ii < 1) .or. (ii > nx) .or. (jj < 1) .or. (jj > ny) .or.  &
          (kk < 1) .or. (kk > nz)) cycle
        
      !--hack: factor of i1, j1, k1 in front control sign
      ntmp(1) = i1 * (phi(ii, j , k ) - phi(i, j, k)) / dx
      ntmp(2) = j1 * (phi(i , jj, k ) - phi(i, j, k)) / dy
      ntmp(3) = k1 * (phi(i , j , kk) - phi(i, j, k)) / dz

      if (mag (ntmp) > eps) then
        ntmp = ntmp / mag (ntmp)
        return
      end if

    end do
  end do
end do

!--if get here, then we have failed to get a normal

if (present (ierr)) then
  ierr = -1  !--set error code
else
  call mesg (sub_name, 'failed to get normal at', (/ i, j, k /))
end if

ntmp = (/ 1._rp, 0._rp, 0._rp /)  !--what else to try??

end subroutine fix_norm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--uses centered finite differences to calculate unit normal from 
!  the signed distance function
!--uses 1-sided finite differences near top and bottom boundaries
!--not caring about speed too much here
!--always returns a unit vector
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fill_norm ()
implicit none

character (*), parameter :: sub_name = mod_name // '.fill_norm'

integer, parameter :: bad_size = 100  !--space for storage bad norm locations

real (rp), parameter :: eps = 1000._rp * epsilon (0._rp)

integer :: i, j, k
$if ($MPI)
  integer, parameter :: tag = 50
  integer :: datasize
$endif
integer :: nbad
integer :: bad(3, bad_size)
integer :: ierr

real (rp) :: ntmp(nd)

!---------------------------------------------------------------------
$if ($VERBOSE)
call enter_sub (sub_name)
$endif

nbad = 0
bad = iBOGUS

!--convention: phi(nz) is valid, but norm(nz) is not
do k = 1, nz - 1
  do j = 1, ny
    do i = 1, nx

      ntmp(1) = safe_cd (i, j, k, 1, phi)
      ntmp(2) = safe_cd (i, j, k, 2, phi)
      ntmp(3) = safe_cd (i, j, k, 3, phi)

      if (mag (ntmp) > eps) then
 
        norm(:, i, j, k) = ntmp / mag (ntmp)

        !--double-check
        if (mag (norm(:, i, j, k)) <= eps) then
          call error (sub_name, 'norm failed double-check')
        end if

      else  !--special treatment required

        call fix_norm (i, j, k, eps, ntmp, ierr)
             !--always returns a unit vector
        norm(:, i, j, k) = ntmp

        !--check error code for "bad" value
        if (ierr < 0) then
        
          nbad = nbad + 1
          
          if (nbad <= bad_size) then
            bad(:, nbad) = (/ i, j, k /)
          else
            call error (sub_name, 'bad_size is too small')
          end if
          
        end if
        
      end if

    end do
  end do
end do

if (nbad > 0) then  !--attempt to fix "bad" values

  do i = 1, nbad
    call avg_norm (bad(1, i), bad(2, i), bad(3, i), eps)
         !--assumes bad points are isolated
  end do

end if

!--fill top level with BOGUS
norm(:, :, :, nz) = BOGUS

$if ($MPI)

  norm(:, :, :, 0) = BOGUS

  !--sync 0-level with level below it
  !--this is required to get norm at k=1 for w-nodes
  datasize = size (norm(:, :, :, 1))  !--size of a z-slice
  call mpi_sendrecv (norm(1, 1, 1, nz-1), datasize, MPI_RPREC, up, tag,  &
                     norm(1, 1, 1, 0), datasize, MPI_RPREC, down, tag,   &
                     comm, status, ierr)

$endif

$if ($VERBOSE)
call exit_sub (sub_name)
$endif

end subroutine fill_norm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--calculates norm at i, j, k by averaging neighboring values
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine avg_norm (i, j, k, eps)
implicit none

integer, intent (in) :: i, j, k
real (rp), intent (in) :: eps

character (*), parameter :: sub_name = mod_name // '.avg_norm'

real (rp) :: avg(nd)

!---------------------------------------------------------------------

avg = 0._rp

!--check we are in bounds
if ( ((1 > k) .or. (k > nz-1)) .or. ((1 > j) .or. (j > ny)) .or.  &
     ((1 > i) .or. (i > nx)) ) then
  call error (sub_name, 'out of bounds (i, j, k)=', (/ i, j, k /))
end if

if (1 < k) avg = avg + norm(:, i, j, k-1)
if (k < nz-1) avg = avg + norm(:, i, j, k+1)

if (1 < j) avg = avg + norm(:, i, j-1, k)
if (j < ny) avg = avg + norm(:, i, j+1, k)

if (1 < i) avg = avg + norm(:, i-1, j, k)
if (i < nx) avg = avg + norm(:, i+1, j, k)

!--normalize if successful, complain otherwise
if (mag (avg) > eps) then
  avg = avg / mag (avg)
  norm(:, i, j, k) = avg
else
  call mesg (sub_name, 'failed to get normal at', (/ i, j, k /))
end if

end subroutine avg_norm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function mag (v)
implicit none

real (rp) :: mag

real (rp), intent (in) :: v(nd)

!---------------------------------------------------------------------

mag = sqrt (dot_product (v, v))

end function mag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine level_set_init ()
use param, only : dx, dy, dz  !--in addition to those above
implicit none

character (*), parameter :: sub_name = mod_name // '.level_set_init'

character (*), parameter :: fphi_in_base = 'phi.out'
character (*), parameter :: fnorm_out_base = 'norm.dat'
character (*), parameter :: MPI_suffix = '.c'

integer, parameter :: lun = 1

logical, parameter :: do_write_norm = .true.

character (128) :: fphi_in, fnorm_out

integer :: i, j, k

logical :: exst, opn

real (rp) :: x, y, z

!---------------------------------------------------------------------
$if ($VERBOSE)
call enter_sub (sub_name)
$endif

$if($MPI)
!  Check that the buffer arrays DO NOT extent beyond neighboring processors
if( nphitop >= Nz .or. nphibot >= Nz .or. &
    nveltop >= Nz .or. nvelbot >= Nz .or. &
    ntautop >= Nz .or. ntaubot >= Nz .or. &
    nFMMtop >= Nz .or. nFMMbot >= Nz )  &

  call error( sub_name, 'Buffer array extents beyond neighboring processor')
$endif

inquire (unit=lun, exist=exst, opened=opn)

if (.not. exst) call error (sub_name, 'lun =', lun, 'does not exist')
if (opn) call error (sub_name, 'lun =', lun, 'is already open')

$if ($MPI)
  write (fphi_in, '(a,a,i0)') trim (fphi_in_base), MPI_suffix, coord
$else
  fphi_in = trim (fphi_in_base)
$endif

inquire (file=fphi_in, exist=exst, opened=opn)

if (.not. exst) call error (sub_name, 'file ' // fphi_in // ' does not exist')
if (opn) call error (sub_name, 'file ' // fphi_in // ' is aleady open')

$if ($READ_BIG_ENDIAN)
open (lun, file=fphi_in, form='unformatted', action='read', position='rewind', convert='big_endian')
$elseif ($READ_LITTLE_ENDIAN)
open (lun, file=fphi_in, form='unformatted', action='read', position='rewind', convert='little_endian')
$else
open (lun, file=fphi_in, form='unformatted', action='read', position='rewind')
$endif

read (lun) phi(:, :, $lbz:nz)
           !--phi(:, :, 0) will be BOGUS at coord == 0
           !--for now, phi(:, :, nz) will be valid at coord = nproc - 1
close (lun)

call mesg (sub_name, 'level set function initialized')

!--calculate the normal
!--provides 0:nz-1, except at coord = 0 it provides 1:nz-1
call fill_norm ()

if (do_write_norm) then

  $if ($MPI)
    write (fnorm_out, '(a,a,i0)') trim (fnorm_out_base), MPI_suffix, coord
  $else
    fnorm_out = trim (fnorm_out_base)
  $endif

  !--output normal in ascii, for checking purposes
  !--recall nz-level is BOGUS
  open (lun, file=fnorm_out, action='write')

  write (lun, '(a)') 'variables = "x" "y" "z" "n1" "n2" "n3"'
  write (lun, '(3(a,i0))') 'zone, f=point, i=', nx, ', j=', ny, ', k=', nz-1

  do k = 1, nz-1

    z = (k - 0.5_rp) * dz

    do j = 1, ny

      y = (j - 1) * dy
    
      do i = 1, nx

        x = (i - 1) * dx
      
        write (lun, '(6(1x,es12.5))') x, y, z, norm(:, i, j, k)

      end do
    
    end do

  end do

  close (lun)

end if

$if ($VERBOSE)
call exit_sub (sub_name)
$endif

end subroutine level_set_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function cross_product (a, b)
implicit none

real (rp) :: cross_product(nd)

real (rp), intent (in) :: a(nd), b(nd)

!----------------------------------------------------------------------

cross_product(1) = a(2) * b(3) - a(3) * b(2)
cross_product(2) = a(3) * b(1) - a(1) * b(3)
cross_product(3) = a(1) * b(2) - a(2) * b(1)

end function cross_product

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module level_set
