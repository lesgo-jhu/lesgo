module press_util

implicit none

private
public press_stag_array

contains

!**********************************************************************
subroutine press_stag_array(p_hat,dfdx,dfdy)   
!**********************************************************************
! p_hat contains the physical space pressure on exit
!--provides p_hat, dfdx, dfdy 0:nz-1
!-------------------    
! Boundary Layer version with 4th order derivs in vertical.
!  04 December 1995
!	Mods.
!	12/6: added upper and lower boundary conditions.	
!    12/8: Corrected sign on imag. x and y deriv of RHS.
!	12/17 Forcing pressure equal zero at wall 
!			Removed forcing of <P>=0 for all z.
!		    Will need to change BC when hetero. surface stress.
!	12/17 Added ficticious node below wall (Note solver on Nz+1 x Nz+1
!          prev. version of this subroutine saved as p_bl4.old.12.17
!    12/18 Revised 2st deriv stencil at wall (avg of deriv at -1 and at 1)
!    12/21 Redid FDD to 2nd order accurate.
!    12/22 Now value of P(wall) diagnosed from prev P, to match gradient BC
!....1/13: major changes.
!		Broke out mean pressure for separate solution
!    1/20 back to tridag for natrix solution (same sol'n as LUDCMP...A Keeper!)
!....1/23 Staggered solution
!.........Using Nz+1 levels for computing P, but tossing out level below ground
!....4/1 Changed sign on Div T_iz at wall and lid (five places)
!-------------------          
use types,only:rprec
use param
use messages
use sim_param,only:u,v,w,divtz
use fft
use emul_complex, only : OPERATOR(.MULI.)
$if ($DEBUG)
use debug_mod
$endif

implicit none      

! z indexing must start at 0 for p_hat 
real(rprec), dimension(:, :, 0:), intent(inout) :: p_hat
real(rprec), dimension(:, :, :), intent(inout) :: dfdx,dfdy

real(rprec) , dimension(:, :, :), allocatable :: rH_x,rH_y,rH_z
real(rprec) , dimension(:, :), allocatable :: rtopw, rbottomw

real(rprec) ::const

integer::jx,jy,jz,k

integer :: ir, ii ! Used for complex emulation of real array

character (64) :: fname

$if ($DEBUG)
logical, parameter :: DEBUG = .false.
logical, parameter :: TRI_DEBUG = .false.
$endif

integer :: jz_min

real(rprec), dimension(:, :, :), allocatable :: RHS_col
real(rprec), dimension(:, :, :), allocatable :: a,b,c

real(rprec), dimension(2) :: aH_x, aH_y ! Used to emulate complex scalar

!---------------------------------------------------------------------
const = 1._rprec/(nx*ny)

! Allocate
allocate ( rH_x(ld,ny,lbz:nz), rH_y(ld,ny,lbz:nz), rH_z(ld,ny,lbz:nz) )
allocate ( rtopw(ld,ny), rbottomw(ld,ny) )
allocate ( RHS_col(ld,ny,nz+1) )
allocate ( a(lh,ny,nz+1), b(lh,ny,nz+1), c(lh,ny,nz+1) )

$if ($VERBOSE)
write (*, *) 'started press_stag_array'
$endif

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
  p_hat(:, :, 0) = (0._rprec, 0._rprec)
else
  p_hat(:, :, 0) = BOGUS
end if

!==========================================================================
! Get the right hand side ready 
! Loop over levels    
do jz=1,nz-1  !--experiment: was nz here (see below experiments)
! temp storage for sum of RHS terms.  normalized for fft
! sc: recall that the old timestep guys already contain the pressure
!   term

   rH_x(:, :, jz) = const / tadv1 * (u(:, :, jz) / dt)
   rH_y(:, :, jz) = const / tadv1 * (v(:, :, jz) / dt)
   rH_z(:, :, jz) = const / tadv1 * (w(:, :, jz) / dt)

   call rfftwnd_f77_one_real_to_complex(forw,rH_x(:,:,jz),fftwNull_p)
   call rfftwnd_f77_one_real_to_complex(forw,rH_y(:,:,jz),fftwNull_p)
   call rfftwnd_f77_one_real_to_complex(forw,rH_z(:,:,jz),fftwNull_p)   
end do



$if ($MPI)
!  H_x(:, :, 0) = BOGUS
!  H_y(:, :, 0) = BOGUS
!  H_z(:, :, 0) = BOGUS
  !Careful - only update real values (odd indicies)
  rH_x(1:ld:2,:,0) = BOGUS
  rH_y(1:ld:2,:,0) = BOGUS
  rH_z(1:ld:2,:,0) = BOGUS
$endif

!--experiment
!--this causes blow-up
!H_x(:, :, nz) = BOGUS
!H_y(:, :, nz) = BOGUS
!Careful - only update real values (odd indicies)
rH_x(1:ld:2,:,nz) = BOGUS
rH_y(1:ld:2,:,nz) = BOGUS

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
  !H_z(:, :, nz) = (0._rprec, 0._rprec)
  rH_z(:,:,nz) = 0._rprec
else
  !H_z(:, :, nz) = BOGUS  !--perhaps this should be 0 on top process?
  !Careful - only update real values (odd indicies)
  rH_z(1:ld:2,:,nz) = BOGUS !--perhaps this should be 0 on top process?
end if

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
  rbottomw(:, :) = const * divtz(:, :, 1)
  call rfftwnd_f77_one_real_to_complex (forw, rbottomw(:, :), fftwNull_p)
end if

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
  rtopw(:, :) = const * divtz(:, :, nz)
  call rfftwnd_f77_one_real_to_complex (forw, rtopw(:, :), fftwNull_p)
end if

! set oddballs to 0
! probably can get rid of this if we're more careful below
!H_x(lh, :, 1:nz-1)=0._rprec
!H_y(lh, :, 1:nz-1)=0._rprec
!H_z(lh, :, 1:nz-1)=0._rprec
!H_x(:, ny/2+1, 1:nz-1)=0._rprec
!H_y(:, ny/2+1, 1:nz-1)=0._rprec
!H_z(:, ny/2+1, 1:nz-1)=0._rprec
rH_x(ld-1:ld, :, 1:nz-1)=0._rprec
rH_y(ld-1:ld, :, 1:nz-1)=0._rprec
rH_z(ld-1:ld, :, 1:nz-1)=0._rprec
rH_x(:, ny/2+1, 1:nz-1)=0._rprec
rH_y(:, ny/2+1, 1:nz-1)=0._rprec
rH_z(:, ny/2+1, 1:nz-1)=0._rprec

!--with MPI; topw and bottomw are only on top & bottom processes
!topw(lh, :)=0._rprec
!topw(:, ny/2+1)=0._rprec
!bottomw(lh, :)=0._rprec
!bottomw(:, ny/2+1)=0._rprec
rtopw(ld-1:ld, :)=0._rprec
rtopw(:, ny/2+1)=0._rprec
rbottomw(ld-1:ld, :)=0._rprec
rbottomw(:, ny/2+1)=0._rprec

!==========================================================================
! Loop over (Kx,Ky) to solve for Pressure amplitudes

$if ($DEBUG)
if (TRI_DEBUG) then
  a = BOGUS
  b = BOGUS
  c = BOGUS
  !RHS_col = BOGUS
  !Careful - only update real values (odd indicies)
  RHS_col(1:ld:2,:,:) = BOGUS
end if
$endif

!--switch order of inner/outer loops here
if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then

  !  a,b,c are treated as the real part of a complex array
  a(:, :, 1) = BOGUS  !--was 0._rprec
  b(:, :, 1) = -1._rprec
  c(:, :, 1) = 1._rprec
  !RHS_col(:, :, 1) = -dz * bottomw(:, :)
  RHS_col(:,:,1) = -dz * rbottomw(:,:)

  $if ($DEBUG)
  if (TRI_DEBUG) then
    a(:, :, 1) = BOGUS  !--was 0._rprec
    b(:, :, 1) = 2._rprec
    c(:, :, 1) = 1._rprec
    !RHS_col(:, :, 1) = 1._rprec
    !Careful - only update real values (odd indicies)
    RHS_col(1:ld:2,:,1) = 1._rprec
  end if
  $endif

  jz_min = 2

else

  jz_min = 1

end if

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then

  !--top nodes
  a(:, :, nz+1) = -1._rprec
  b(:, :, nz+1) = 1._rprec
  c(:, :, nz+1) = BOGUS  !--was 0._rprec
  !RHS_col(:, :, nz+1) = -topw(:, :) * dz
  RHS_col(:,:,nz+1) = -dz * rtopw(:,:)

  $if ($DEBUG)
  if (TRI_DEBUG) then
    a(:, :, nz+1) = 1._rprec
    b(:, :, nz+1) = 2._rprec
    c(:, :, nz+1) = BOGUS  !--was 0._rprec
    $if ($MPI)
      !RHS_col(:, :, nz+1) = real (nz+1 + coord * (nz-1), rprec)
      !Careful - only update real values (odd indicies)
      RHS_col(1:ld:2,:,nz+1) = real (nz+1 + coord * (nz-1), rprec)
    $else
      !RHS_col(:, :, nz+1) = real (nz+1, rprec)
      !Careful - only update real values (odd indicies)
      RHS_col(1:ld:2,:,nz+1) =  real (nz+1, rprec)
    $endif
  end if
  $endif

end if

$if ($DEBUG)
if (DEBUG) write (*, *) coord, ' before H send/recv'
$endif

$if ($MPI)
  !--could maybe combine some of these to less communication is needed
  !--fill H_x, H_y, H_z at jz=0 (from nz-1)
  !--cant just change lbz above, since u,v,w (jz=0) are not in sync yet
  !call mpi_sendrecv (H_x(1, 1, nz-1), lh*ny, MPI_CPREC, up, 1,  &
  !                   H_x(1, 1, 0), lh*ny, MPI_CPREC, down, 1,   &
  !                   comm, status, ierr)
  !call mpi_sendrecv (H_y(1, 1, nz-1), lh*ny, MPI_CPREC, up, 2,  &
  !                   H_y(1, 1, 0), lh*ny, MPI_CPREC, down, 2,   &
  !                   comm, status, ierr)
  !call mpi_sendrecv (H_z(1, 1, nz-1), lh*ny, MPI_CPREC, up, 3,  &
  !                   H_z(1, 1, 0), lh*ny, MPI_CPREC, down, 3,   &
  !                   comm, status, ierr)
  call mpi_sendrecv (rH_x(1, 1, nz-1), ld*ny, MPI_RPREC, up, 1,  &
                     rH_x(1, 1, 0), ld*ny, MPI_RPREC, down, 1,   &
                     comm, status, ierr)
  call mpi_sendrecv (rH_y(1, 1, nz-1), ld*ny, MPI_RPREC, up, 2,  &
                     rH_y(1, 1, 0), ld*ny, MPI_RPREC, down, 2,   &
                     comm, status, ierr)
  call mpi_sendrecv (rH_z(1, 1, nz-1), ld*ny, MPI_RPREC, up, 3,  &
                     rH_z(1, 1, 0), ld*ny, MPI_RPREC, down, 3,   &
                     comm, status, ierr)

  !--fill H_x, H_y, H_z at jz=nz (from 1)
  !call mpi_sendrecv (H_x(1, 1, 1), lh*ny, MPI_CPREC, down, 4,  &
  !                   H_x(1, 1, nz), lh*ny, MPI_CPREC, up, 4,   &
  !                   comm, status, ierr)
  !call mpi_sendrecv (H_y(1, 1, 1), lh*ny, MPI_CPREC, down, 5,  &
  !                   H_y(1, 1, nz), lh*ny, MPI_CPREC, up, 5,   &
  !                   comm, status, ierr)
  !call mpi_sendrecv (H_z(1, 1, 1), lh*ny, MPI_CPREC, down, 6,  &
  !                   H_z(1, 1, nz), lh*ny, MPI_CPREC, up, 6,   &
  !                   comm, status, ierr)
  call mpi_sendrecv (rH_z(1, 1, 1), ld*ny, MPI_RPREC, down, 6,  &
                     rH_z(1, 1, nz), ld*ny, MPI_RPREC, up, 6,   &
                     comm, status, ierr)                     
$endif

$if ($DEBUG)
if (DEBUG) then
  write (*, *) coord, ' after H send/recv'
  !call DEBUG_write (H_x(:, :, 1:nz), 'w.H_x')
  !call DEBUG_write (H_y(:, :, 1:nz), 'w.H_y')
  !call DEBUG_write (H_z(:, :, 1:nz), 'w.H_z')
  !call DEBUG_write (topw, 'w.topw')
  !call DEBUG_write (bottomw, 'w.bottomw')
  call DEBUG_write (rH_x(:, :, 1:nz), 'w.H_x')
  call DEBUG_write (rH_y(:, :, 1:nz), 'w.H_y')
  call DEBUG_write (rH_z(:, :, 1:nz), 'w.H_z')
  call DEBUG_write (rtopw, 'w.topw')
  call DEBUG_write (rbottomw, 'w.bottomw')  

end if
$endif

do jz = jz_min, nz
  do jy = 1, ny

    if (jy == ny/2 + 1) cycle

    do jx = 1, lh-1

      if (jx*jy == 1) cycle

      ii = 2*jx   ! imaginary index 
      ir = ii - 1 ! real index

      ! JDA dissertation, eqn(2.85) a,b,c=coefficients and RHS_col=r_m
      !a(jx, jy, jz) = 1._rprec/(dz**2)
      !b(jx, jy, jz) = -(kx(jx, jy)**2 + ky(jx, jy)**2 + 2._rprec/(dz**2))
      !c(jx, jy, jz) = 1._rprec/(dz**2)
      !RHS_col(jx, jy, jz) = eye * (kx(jx, jy) * H_x(jx, jy, jz-1) +   &
      !                             ky(jx, jy) * H_y(jx, jy, jz-1)) +  &
      !                      (H_z(jx, jy, jz) - H_z(jx, jy, jz-1)) / dz
      a(jx, jy, jz) = 1._rprec/(dz**2)
      b(jx, jy, jz) = -(kx(jx, jy)**2 + ky(jx, jy)**2 + 2._rprec/(dz**2))
      c(jx, jy, jz) = 1._rprec/(dz**2)   
      !  Compute eye * kx * H_x 
!      call mult_real_complex_imag( rH_x(ir:ii, jy, jz-1), kx(jx, jy), aH_x )
      aH_x = rH_x(ir:ii, jy, jz-1) .MULI. kx(jx,jy)

      !  Compute eye * ky * H_y
!      call mult_real_complex_imag( rH_y(ir:ii, jy, jz-1), ky(jx, jy), aH_y )           
      aH_y = rH_y(ir:ii, jy, jz-1) .MULI. ky(jx,jy) 

      RHS_col(ir:ii,jy,jz) =  aH_x + aH_y + (rH_z(ir:ii, jy, jz) - rH_z(ir:ii, jy, jz-1)) / dz

      $if ($DEBUG)
      if (TRI_DEBUG) then
        a(jx, jy, jz) = 1._rprec
        b(jx, jy, jz) = 2._rprec
        c(jx, jy, jz) = 1._rprec
        $if ($MPI)
          !RHS_col(jx, jy, jz) = jz + coord * (nz-1)
          !Careful - only update real value
          RHS_col(ir,jy,jz) = jz + coord * (nz-1)
        $else
          !RHS_col(jx, jy, jz) = jz
          !Careful - only update real value
          RHS_col(ir,jy,jz) = jz
        $endif
      end if
      $endif
     
    end do
  end do
end do

!a = 1._rprec
!c = 1._rprec
!b = 2._rprec
!do jz=1,nz+1
!  $if ($MPI)
!    RHS_col(:, :, jz) = jz + coord * (nz-1)
!  $else
!    RHS_col(:, :, jz) = jz
!  $endif
!end do

$if ($DEBUG)
if (DEBUG) then
  write (*, *) coord, ' before tridag_array'
  call DEBUG_write (a, 'v.a')
  call DEBUG_write (b, 'v.b')
  call DEBUG_write (c, 'v.c')
  !call DEBUG_write (RHS_col, 'v.RHS_col')
  call DEBUG_write( RHS_col, 'v.RHS_col')
end if
$endif

!if (DEBUG) then
!    write (*, *) 'jt, pre tridag p_hat(1,1,0) = ', jt, p_hat(1, 1, 0)
!end if

!--this skips zero wavenumber solution, nyquist freqs
!call tridag_array (a, b, c, RHS_col, p_hat)
$if ($MPI)
  !call tridag_array_pipelined (0, a, b, c, RHS_col, p_hat)
  call tridag_array_pipelined( 0, a, b, c, RHS_col, p_hat )
$else
  !call tridag_array (a, b, c, RHS_col, p_hat)
  call tridag_array (a, b, c, RHS_col, p_hat)
$endif

$if ($DEBUG)
if (DEBUG) then
  write (*, *) coord, ' after tridag_array'
  call DEBUG_write (p_hat, 'press_stag_array.c.p_hat')
endif
$endif

!--zero-wavenumber solution
$if ($MPI)
  !--wait for p_hat(1, 1, 1) from "down"
  !call mpi_recv (p_hat(1, 1, 1), 1, MPI_CPREC, down, 8, comm, status, ierr)
  call mpi_recv (p_hat(1:2, 1, 1), 2, MPI_RPREC, down, 8, comm, status, ierr)
$endif

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then

  !p_hat(1, 1, 0) = 0._rprec
  !p_hat(1, 1, 1) = p_hat(1, 1, 0) - dz * bottomw(1, 1)
  p_hat(1:2, 1, 0) = 0._rprec
  p_hat(1:2, 1, 1) = p_hat(1:2,1,0) - dz * rbottomw(1:2,1)

end if

do jz = 2, nz
  ! JDA dissertation, eqn(2.88)
  !p_hat(1, 1, jz) = p_hat(1, 1, jz-1) + H_z(1, 1, jz)*dz
  p_hat(1:2, 1, jz) = p_hat(1:2, 1, jz-1) + rH_z(1:2, 1, jz) * dz
end do

$if ($MPI)
  !--send p_hat(1, 1, nz) to "up"
  !call mpi_send (p_hat(1, 1, nz), 1, MPI_CPREC, up, 8, comm, ierr)
  call mpi_send (p_hat(1:2, 1, nz), 2, MPI_RPREC, up, 8, comm, ierr)
$endif

$if ($MPI)
  !--make sure 0 <-> nz-1 are syncronized
  !-- 1 <-> nz should be in sync already
  !call mpi_sendrecv (p_hat(1, 1, nz-1), lh*ny, MPI_CPREC, up, 2,  &
  !                   p_hat(1, 1, 0), lh*ny, MPI_CPREC, down, 2,   &
  !                   comm, status, ierr)
  call mpi_sendrecv (p_hat(1, 1, nz-1), ld*ny, MPI_RPREC, up, 2,  &
                     p_hat(1, 1, 0), ld*ny, MPI_RPREC, down, 2,   &
                     comm, status, ierr)  

$endif

!--zero the nyquist freqs
!p_hat(lh, :, :) = 0._rprec
!p_hat(:, ny/2+1, :) = 0._rprec
p_hat(ld-1:ld, :, :) = 0._rprec
p_hat(:, ny/2+1, :) = 0._rprec

$if ($DEBUG)
if (DEBUG) call DEBUG_write (p_hat, 'press_stag_array.d.p_hat')
$endif

!=========================================================================== 
!...Now need to get p_hat(wave,level) to physical p(jx,jy,jz)   
!.....Loop over height levels     

$if ($DEBUG)
if (DEBUG) write (*, *) 'press_stag_array: before inverse FFT'
$endif

call rfftwnd_f77_one_complex_to_real(back,p_hat(:,:,0),fftwNull_p)
do jz=1,nz-1  !--used to be nz
do jy=1,ny
do jx=1,lh
  ii = 2*jx
  ir = ii - 1
! complex
   !dfdx(jx,jy,jz)=eye*kx(jx,jy)*p_hat(jx,jy,jz)
   !dfdy(jx,jy,jz)=eye*ky(jx,jy)*p_hat(jx,jy,jz)
   !call mult_real_complex_imag( p_hat(ir:ii,jy,jz), kx(jx,jy), dfdx(ir:ii,jy,jz) )
   !call mult_real_complex_imag( p_hat(ir:ii,jy,jz), ky(jx,jy), dfdy(ir:ii,jy,jz) )
   dfdx(ir:ii,jy,jz) = p_hat(ir:ii,jy,jz) .MULI. kx(jx,jy) 
   dfdy(ir:ii,jy,jz) = p_hat(ir:ii,jy,jz) .MULI. ky(jx,jy) 
   

! note the oddballs of p_hat are already 0, so we should be OK here
end do
end do
call rfftwnd_f77_one_complex_to_real(back,dfdx(:,:,jz),fftwNull_p)
call rfftwnd_f77_one_complex_to_real(back,dfdy(:,:,jz),fftwNull_p)
call rfftwnd_f77_one_complex_to_real(back,p_hat(:,:,jz),fftwNull_p)
end do

!--nz level is not needed elsewhere (although its valid)
dfdx(:, :, nz) = BOGUS
dfdy(:, :, nz) = BOGUS
p_hat(:, :, nz) = BOGUS

$if ($VERBOSE)
write (*, *) 'finished press_stag_array'
$endif

end subroutine press_stag_array

end module press_util
