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
subroutine press_stag_array()
!*******************************************************************************
!
! Calculate the pressure and its derivatives on exit. Everything is in physical
! space on exit.
!
use types, only : rprec
use param
use messages
use sim_param, only : u, v, w, divtz, p, dpdx, dpdy, dpdz
use fft
use emul_complex, only : OPERATOR(.MULI.)

implicit none

real(rprec) :: const, const2, const3, const4
integer :: jx, jy, jz
integer :: ir, ii
integer :: jz_min

real(rprec), save, dimension(:,:,:), allocatable :: rH_x, rH_y, rH_z
real(rprec), save, dimension(:,:), allocatable :: rtopw, rbottomw
real(rprec), save, dimension(:,:,:), allocatable :: RHS_col
real(rprec), save, dimension(:,:,:), allocatable :: a, b, c

logical, save :: arrays_allocated = .false.

real(rprec), dimension(2) :: aH_x, aH_y

! Specifiy cached constants
const = 1._rprec/(nx*ny)
const2 = const/tadv1/dt
const3 = 1._rprec/(dz**2)
const4 = 1._rprec/(dz)

! Allocate arrays
if( .not. arrays_allocated ) then
    allocate ( rH_x(ld,ny,lbz:nz), rH_y(ld,ny,lbz:nz), rH_z(ld,ny,lbz:nz) )
    allocate ( rtopw(ld,ny), rbottomw(ld,ny) )
    allocate ( RHS_col(ld,ny,nz+1) )
    allocate ( a(lh,ny,nz+1), b(lh,ny,nz+1), c(lh,ny,nz+1) )

    arrays_allocated = .true.
endif

if (coord == 0) then
    p(:,:,0) = 0._rprec
#ifdef PPSAFETYMODE
else
    p(:,:,0) = BOGUS
#endif
end if

! Get the right hand side ready
! Loop over levels
! Recall that the old timestep guys already contain the pressure
do jz = 1, nz-1
    rH_x(:,:,jz) = const2 * u(:,:,jz)
    rH_y(:,:,jz) = const2 * v(:,:,jz)
    rH_z(:,:,jz) = const2 * w(:,:,jz)

    call dfftw_execute_dft_r2c(forw, rH_x(:,:,jz), rH_x(:,:,jz))
    call dfftw_execute_dft_r2c(forw, rH_y(:,:,jz), rH_y(:,:,jz))
    call dfftw_execute_dft_r2c(forw, rH_z(:,:,jz), rH_z(:,:,jz))
end do

#if defined(PPMPI) && defined(PPSAFETYMODE)
  !Careful - only update real values (odd indicies)
  rH_x(1:ld:2,:,0) = BOGUS
  rH_y(1:ld:2,:,0) = BOGUS
  rH_z(1:ld:2,:,0) = BOGUS
#endif

#ifdef PPSAFETYMODE
!Careful - only update real values (odd indicies)
rH_x(1:ld:2,:,nz) = BOGUS
rH_y(1:ld:2,:,nz) = BOGUS
#endif

#ifdef PPMPI
if (coord == nproc-1) then
    rH_z(:,:,nz) = const2 * w(:,:,nz)
    call dfftw_execute_dft_r2c(forw, rH_z(:,:,nz), rH_z(:,:,jz))
#ifdef PPSAFETYMODE
else
    rH_z(1:ld:2,:,nz) = BOGUS !--perhaps this should be 0 on top process?
#endif
endif
#else
rH_z(:,:,nz) = const2 * w(:,:,nz)
call dfftw_execute_dft_r2c(forw, rH_z(:,:,nz), rH_z(:,:,jz))
#endif

if (coord == 0) then
    rbottomw(:,:) = const * divtz(:,:,1)
    call dfftw_execute_dft_r2c(forw, rbottomw, rbottomw )
end if

#ifdef PPMPI
if (coord == nproc-1) then
#endif
    rtopw(:,:) = const * divtz(:,:,nz)
    call dfftw_execute_dft_r2c(forw, rtopw, rtopw)
#ifdef PPMPI
endif
#endif


! set oddballs to 0
rH_x(ld-1:ld,:,1:nz-1) = 0._rprec
rH_y(ld-1:ld,:,1:nz-1) = 0._rprec
rH_z(ld-1:ld,:,1:nz-1) = 0._rprec
rH_x(:,ny/2+1,1:nz-1) = 0._rprec
rH_y(:,ny/2+1,1:nz-1) = 0._rprec
rH_z(:,ny/2+1,1:nz-1) = 0._rprec
! should also set to zero for rH_z (nz) on coord == nproc-1
if (coord == nproc-1) then
    rH_z(ld-1:ld,:,nz) = 0._rprec
    rH_z(:,ny/2+1,nz) = 0._rprec
end if

! with MPI; topw and bottomw are only on top & bottom processes
rtopw(ld-1:ld, :) = 0._rprec
rtopw(:, ny/2+1) = 0._rprec
rbottomw(ld-1:ld, :) = 0._rprec
rbottomw(:, ny/2+1) = 0._rprec

! Loop over (Kx,Ky) to solve for Pressure amplitudes
if (coord == 0) then
    !  a, b, and c are treated as the real part of a complex array
#ifdef PPSAFETYMODE
    a(:,:,1) = BOGUS
#endif
    b(:,:,1) = -1._rprec
    c(:,:,1) = 1._rprec
    RHS_col(:,:,1) = -dz * rbottomw(:,:)

    jz_min = 2
else
  jz_min = 1
end if

#ifdef PPMPI
if (coord == nproc-1) then
#endif
    !--top nodes
    a(:,:,nz+1) = -1._rprec
    b(:,:,nz+1) = 1._rprec
#ifdef PPSAFETYMODE
    c(:,:,nz+1) = BOGUS
#endif
    RHS_col(:,:,nz+1) = -dz * rtopw(:,:)
#ifdef PPMPI
endif
#endif

#ifdef PPMPI
    call mpi_sendrecv (rH_x(1, 1, nz-1), ld*ny, MPI_RPREC, up, 1,              &
        rH_x(1, 1, 0), ld*ny, MPI_RPREC, down, 1, comm, status, ierr)
    call mpi_sendrecv (rH_y(1, 1, nz-1), ld*ny, MPI_RPREC, up, 2,              &
        rH_y(1, 1, 0), ld*ny, MPI_RPREC, down, 2, comm, status, ierr)
    call mpi_sendrecv (rH_z(1, 1, nz-1), ld*ny, MPI_RPREC, up, 3,              &
        rH_z(1, 1, 0), ld*ny, MPI_RPREC, down, 3, comm, status, ierr)
    call mpi_sendrecv (rH_z(1, 1, 1), ld*ny, MPI_RPREC, down, 6,               &
        rH_z(1, 1, nz), ld*ny, MPI_RPREC, up, 6, comm, status, ierr)
#endif

do jz = jz_min, nz
do jy = 1, ny
    if (jy == ny/2 + 1) cycle

    do jx = 1, lh-1

        if (jx*jy == 1) cycle

        ii = 2*jx   ! imaginary index
        ir = ii - 1 ! real index

        ! JDA dissertation, eqn(2.85) a,b,c=coefficients and RHS_col=r_m
        a(jx, jy, jz) = const3
        b(jx, jy, jz) = -(kx(jx, jy)**2 + ky(jx, jy)**2 + 2._rprec*const3)
        c(jx, jy, jz) = const3

        !  Compute eye * kx * H_x
        aH_x(1) = -rH_x(ii,jy,jz-1) * kx(jx,jy)
        aH_x(2) =  rH_x(ir,jy,jz-1) * kx(jx,jy)
        aH_y(1) = -rH_y(ii,jy,jz-1) * ky(jx,jy)
        aH_y(2) =  rH_y(ir,jy,jz-1) * ky(jx,jy)

        RHS_col(ir:ii,jy,jz) =  aH_x + aH_y + (rH_z(ir:ii, jy, jz) -           &
            rH_z(ir:ii, jy, jz-1)) *const4

    end do
end do
end do

! this skips zero wavenumber solution, nyquist freqs
call tridag_array (a, b, c, RHS_col, p)

! zero-wavenumber solution
#ifdef PPMPI
! wait for p(1, 1, 1) from "down"
call mpi_recv (p(1:2, 1, 1), 2, MPI_RPREC, down, 8, comm, status, ierr)
#endif

if (coord == 0) then
    p(1:2, 1, 0) = 0._rprec
    p(1:2, 1, 1) = p(1:2,1,0) - dz * rbottomw(1:2,1)
end if

do jz = 2, nz
    ! JDA dissertation, eqn(2.88)
    p(1:2, 1, jz) = p(1:2, 1, jz-1) + rH_z(1:2, 1, jz) * dz
end do

#ifdef PPMPI
! send p(1, 1, nz) to "up"
call mpi_send (p(1:2, 1, nz), 2, MPI_RPREC, up, 8, comm, ierr)
#endif

#ifdef PPMPI
! make sure 0 <-> nz-1 are syncronized
! 1 <-> nz should be in sync already
call mpi_sendrecv (p(1, 1, nz-1), ld*ny, MPI_RPREC, up, 2,                     &
    p(1, 1, 0), ld*ny, MPI_RPREC, down, 2, comm, status, ierr)
#endif

! zero the nyquist freqs
p(ld-1:ld,:,:) = 0._rprec
p(:,ny/2+1,:) = 0._rprec

! Now need to get p(wave,level) to physical p(jx,jy,jz)
! Loop over height levels
call dfftw_execute_dft_c2r(back,p(:,:,0), p(:,:,0))
do jz = 1, nz-1
    do jy = 1, ny
    do jx = 1,lh
        ii = 2*jx
        ir = ii - 1
        dpdx(ir,jy,jz) = -p(ii,jy,jz) * kx(jx,jy)
        dpdx(ii,jy,jz) =  p(ir,jy,jz) * kx(jx,jy)
        dpdy(ir,jy,jz) = -p(ii,jy,jz) * ky(jx,jy)
        dpdy(ii,jy,jz) =  p(ir,jy,jz) * ky(jx,jy)
    end do
    end do

    ! note the oddballs of p are already 0, so we should be OK here
    call dfftw_execute_dft_c2r(back,dpdx(:,:,jz), dpdx(:,:,jz))
    call dfftw_execute_dft_c2r(back,dpdy(:,:,jz), dpdy(:,:,jz))
    call dfftw_execute_dft_c2r(back,p(:,:,jz), p(:,:,jz))
end do

if(coord==nproc-1) call dfftw_execute_dft_c2r(back,p(:,:,nz),p(:,:,nz))

! nz level is not needed elsewhere (although its valid)
#ifdef PPSAFETYMODE
dpdx(:,:,nz) = BOGUS
dpdy(:,:,nz) = BOGUS
if(coord<nproc-1) p(:,:,nz) = BOGUS
#endif

! Final step compute the z-derivative of p
! note: p has additional level at z=-dz/2 for this derivative
dpdz(1:nx, 1:ny, 1:nz-1) = (p(1:nx, 1:ny, 1:nz-1) - p(1:nx, 1:ny, 0:nz-2)) / dz
#ifdef PPSAFETYMODE
if(coord<nproc-1)  dpdz(:,:,nz) = BOGUS
#endif
if(coord==nproc-1) dpdz(1:nx,1:ny,nz) = (p(1:nx,1:ny,nz)-p(1:nx,1:ny,nz-1))/ dz

end subroutine press_stag_array
