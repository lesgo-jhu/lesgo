!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Written by: 
!!
!!   Luis 'Tony' Martinez <tony.mtos@gmail.com> (Johns Hopkins University)
!!
!!   Copyright (C) 2012-2013, Johns Hopkins University
!!
!!   This file is part of The Actuator Turbine Model Library.
!!
!!   LESGO is free software: you can redistribute it 
!!   and/or modify it under the terms of the GNU General Public License as 
!!   published by the Free Software Foundation, either version 3 of the 
!!   License, or (at your option) any later version.
!!
!!   LESGO is distributed in the hope that it will be 
!!   useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
!!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!   GNU General Public License for more details.
!!
!!   You should have received a copy of the GNU General Public License
!!   along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
module hit_inflow
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This module provides basic functionalities to the actuator turbine model
! Real precision variable and dynamic allocation types are stored here

! Real precision from LESGO
use types, only : rprec

! Data from LESGO
use param, only : ny, nz, dt

! Grid definition (LESGO)
use grid_m

implicit none

public :: hit, inflow_HIT, initialize_HIT

type hit_t
    ! This type includes all the information for a HIT case
    real(rprec) :: Lx, Ly, Lz
    integer :: Nx, Ny, Nz

    ! Location of the plane
    real(rprec) :: xloc=0.

    ! The sweeping velocity
    real(rprec) :: U_sweep

    ! The input and output turbulence intensity
    ! By definition up_in/U_sweep = TI_out
    real(rprec) :: up_in, TI_out

    ! The size of the HIT data-set
    real(rprec), allocatable, dimension(:) :: x
    real(rprec), allocatable, dimension(:) :: y
    real(rprec), allocatable, dimension(:) :: z

    ! This is the HIT field (x, y, z)
    real(rprec), allocatable, dimension(:,:,:) :: u
    real(rprec), allocatable, dimension(:,:,:) :: v
    real(rprec), allocatable, dimension(:,:,:) :: w

    ! The plane of data for the inflow (y, z)
    real(rprec), allocatable, dimension(:,:) :: u_plane
    real(rprec), allocatable, dimension(:,:) :: v_plane
    real(rprec), allocatable, dimension(:,:) :: w_plane

    ! Name of the input files
    character(128) :: u_file ! file of u field
    character(128) :: v_file ! file of v field
    character(128) :: w_file ! file of w field

    ! Name of the restart file
    character(128) :: restartFile='restartHIT.dat'

end type hit_t

! Declare turbine array variable
type(hit_t), target :: hit

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine initialize_HIT ()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  This initializes the HIT case by reading input and allocating arrays
!
implicit none

! Grid size in HIT data
real(rprec) :: dx, dy, dz

! Size of the input
integer :: nx_hit, ny_hit, nz_hit

! Index for loop
integer :: i

nx_hit = hit % Nx
ny_hit = hit % Ny
nz_hit = hit % Nz

! Allocate coordinates
allocate( hit % x(nx_hit))
allocate( hit % y(ny_hit))
allocate( hit % z(nz_hit))

! Grid spacing in HIT data
dx = hit % Lx / (nx_hit - 1)
dy = hit % Ly / (ny_hit - 1)
dz = hit % Lz / (nz_hit - 1)

! Create the coordinate arrays x, y, and z

! x
do i = 1, nx_hit
    hit % x(i) = dx * (i - 1)
enddo

! y
do i = 1, ny_hit
    hit % y(i) = dy * (i - 1)
enddo

! z
do i = 1, nz_hit
    hit % z(i) = dz * (i - 1)
enddo

! Allocate velocity input
allocate( hit % u(nx_hit, ny_hit, nz_hit))
allocate( hit % v(nx_hit, ny_hit, nz_hit))
allocate( hit % w(nx_hit, ny_hit, nz_hit))

! Allocate the plane data
allocate(hit % u_plane(ny, nz))
allocate(hit % v_plane(ny, nz))
allocate(hit % w_plane(ny, nz))

! Read the input velocity field
call extract_HIT_data()

! Read the restart file if present
call hit_read_restart()

end subroutine initialize_HIT

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine extract_HIT_data ()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  This extracts the data from the input files
!

implicit none

integer :: readFile=19  ! File number to read
integer :: i, j, k

! Open the velocity field and extract the data
write(*,*) 'Reading u HIT from ', hit % u_file
open( unit=readFile, file=trim(hit % u_file), action='read' )
do i=1, hit % Nx
    do j=1, hit % Ny
        do k=1, hit % Nz
            read(readFile, *) hit % u(i, j, k)
        enddo
    enddo
enddo
close(readFile)

write(*,*) 'Reading v HIT from ', hit % v_file
open( unit=readFile, file=trim(hit % v_file), action='read')
do i=1, hit % Nx
    do j=1, hit % Ny
        do k=1, hit % Nz
            read(readFile, *) hit % v(i, j, k)
        enddo
    enddo
enddo
close(readFile)

write(*,*) 'Reading w HIT from ', hit % w_file
open( unit=readFile, file=trim(hit % w_file), action='read')
do i=1, hit % Nx
    do j=1, hit % Ny
        do k=1, hit % Nz
            read(readFile, *) hit % w(i, j, k)
        enddo
    enddo
enddo
close(readFile)

end subroutine extract_HIT_data

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine compute_HIT_plane_data ()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  This interploates the HIT data for each plane
!

! Inflow velocity
use param, only : inflow_velocity

! Indices for looping in y and z
integer :: j, k

! Compute the location in x
! Sweeping velocity is U_sw = u'_DB / TI
hit % U_sweep = hit % up_in / hit % TI_out

! Update the location of where to sample
hit % xloc = hit % xloc + 1. * dt

! Periodic condition
if (hit % xloc > hit % Lx) then
    hit % xloc = hit % xloc - hit % Lx
endif

! Interpolate data onto plane
do j = 1, Ny
    do k = 1, Nz
        ! Scale the inflow by the inflow velocity
        hit % u_plane(j,k) = inflow_velocity * (1. +  1. / hit % U_sweep *     &
            interpolate3D(hit % xloc, grid % y(j), grid %  z(k),               &
                            hit % x, hit % y, hit % z, hit % u))
        hit % v_plane(j,k) = inflow_velocity * (1. / hit % U_sweep  *          &
            interpolate3D(hit % xloc, grid % y(j), grid %  z(k),               &
                            hit % x, hit % y, hit % z, hit % v))

        hit % w_plane(j,k) = inflow_velocity * (1. / hit % U_sweep *           &
            interpolate3D(hit % xloc, grid % y(j), grid %  zw(k),              &
                            hit % x, hit % y, hit % z, hit % w))
    enddo
enddo

end subroutine compute_HIT_plane_data

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine inflow_HIT ()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Enforces prescribed inflow condition based on an uniform inflow
!  velocity with Homogeneous Isotropic Inflow.
!
use param, only : nx, ny, nz
use sim_param, only : u, v, w
use messages, only : error
use fringe_util
implicit none

integer :: i, i_w
integer :: istart, istart_w
integer :: iplateau
integer :: iend, iend_w

real (rprec) :: alpha, beta

! Compute the velocity at a plane
call compute_HIT_plane_data ()

!--these may be out of 1, ..., nx
call fringe_init( istart, iplateau, iend )

!--wrapped versions
iend_w = modulo (iend - 1, nx) + 1
istart_w = modulo (istart - 1, nx) + 1

! Set end of domain (uniform inflow + turbulence)
u(iend_w, :, :) = hit % u_plane(:,:)
v(iend_w, :, :) = hit % v_plane(:,:)
w(iend_w, :, :) = hit % w_plane(:,:)

!--skip istart since we know vel at istart, iend already
do i = istart + 1, iend - 1

  i_w = modulo (i - 1, nx) + 1

  beta = fringe_weighting( i, istart, iplateau )
  alpha = 1.0_rprec - beta

  u(i_w, 1:ny, 1:nz) = alpha * u(i_w, 1:ny, 1:nz) + beta * hit % u_plane(:,:)
  v(i_w, 1:ny, 1:nz) = alpha * v(i_w, 1:ny, 1:nz) + beta * hit % v_plane(:,:)
  w(i_w, 1:ny, 1:nz) = alpha * w(i_w, 1:ny, 1:nz) + beta * hit % w_plane(:,:)

end do

end subroutine inflow_HIT

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine hit_write_restart()
! This subroutine writes the hit restart information
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

integer :: restartFile=21 ! File to write restart data

! Open the file 
open( unit=restartFile, file=trim(hit % restartFile), status="replace")

write(restartFile,*) 'xloc'

! Store the location x in the file
write(restartFile,*) hit % xloc

close(restartFile)


end subroutine hit_write_restart

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine hit_read_restart()
! This subroutine reads the hit restart information
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

integer :: restartFile=27 ! File to write restart data
logical :: file_exists  ! Flag to check if restart file exists

! Check if restart file exists
inquire(file = trim(hit % restartFile), exist=file_exists)

! If restart file exists read it
if (file_exists) then
    ! Open the file
    open( unit=restartFile, file=trim(hit % restartFile), action="read")

    ! Read past the first line
    read(restartFile,*)

    ! Read the location x in the file
    read(restartFile,*) hit % xloc
    
    close(restartFile)
endif

end subroutine hit_read_restart

!-------------------------------------------------------------------------------
function interpolate3D(xp, yp, zp, x, y, z, u)
! This function does a trilinear intepolation
! xp, yp, zp - 3D point to intepolate
! x, y, z - vectors for doing interpolation
!-------------------------------------------------------------------------------
real(rprec) :: interpolate3D
real(rprec), dimension(:), intent(in) :: x, y, z
real(rprec), dimension(:,:,:), intent(in) :: u
real(rprec), intent(in) ::  xp, yp, zp

real(rprec) :: xd, x0, x1
real(rprec) :: yd, y0, y1
real(rprec) :: zd, z0, z1
integer :: i, i0, i1, nx
integer :: j, j0, j1, ny
integer :: k, k0, k1, nz
real(rprec) :: c00, c01, c10, c11, c0, c1

! Size of the x, y, z vector
nx=size(x)
ny=size(y)
nz=size(z)

!!!!! x
! Point outside (lower bound)
if (xp < x(1)) then 
    ! Pick the first point in the interpolation
    i0 = 1
    i1 = 1
    xd = 1.    
! Point outside (upper bound)
else if (xp > x(nx)) then
    ! Pick the last point in the interpolation
    i0 = nx
    i1 = nx
    xd = 1.
! Point inside
else
    ! Pick the points that are inside the vector
    do i = 2, nx
        if ( ( xp >= x(i-1) ) .and. ( xp <= x(i) ) ) then
            ! Assign the indices between the points
            i0 = i-1
            i1 = i
        endif
    enddo
    ! Points 0 and 1
    x0 = x(i0)
    x1 = x(i1)
    ! The difference
    xd = (xp - x0) / (x1 - x0)
endif

!!!!! y
! Point outside (lower bound)
if (yp < y(1)) then 
    ! Pick the first point in the interpolation
    j0 = 1
    j1 = 1
    yd = 1.    
! Point outside (upper bound)
else if (yp > y(ny)) then
    ! Pick the last point in the interpolation
    j0 = ny
    j1 = ny
    yd = 1.
! Point inside
else
    ! Pick the points that are inside the vector
    do j = 2, ny
        if ( ( yp >= y(j-1) ) .and. ( yp <= y(j) ) ) then
            ! Assign the indices between the points
            j0 = j-1
            j1 = j
        endif
    enddo
    ! Points 0 and 1
    y0 = y(j0)
    y1 = y(j1)
    ! The difference
    yd = (yp - y0) / (y1 - y0)
endif


!!!!! z
! Point outside (lower bound)
if (zp < z(1)) then 
    ! Pick the first point in the interpolation
    k0 = 1
    k1 = 1
    zd = 1.    
! Point outside (upper bound)
else if (zp > z(nz)) then
    ! Pick the last point in the interpolation
    k0 = nz
    k1 = nz
    zd = 1.
! Point inside
else
    ! Pick the points that are inside the vector
    do k = 2, nz
        if ( ( zp >= z(k-1) ) .and. ( zp <= z(k) ) ) then
            ! Assign the indices between the points
            k0 = k-1
            k1 = k
        endif
    enddo
    ! Points 0 and 1
    z0 = z(k0)
    z1 = z(k1)
    ! The difference
    zd = (zp - z0) / (z1 - z0)
endif


! Interpolate along x
c00 = u(i0, j0, k0) * (1.-xd) + u(i1, j0, k0) * xd
c01 = u(i0, j0, k1) * (1.-xd) + u(i1, j0, k1) * xd
c10 = u(i0, j1, k0) * (1.-xd) + u(i1, j1, k0) * xd
c11 = u(i0, j1, k1) * (1.-xd) + u(i1, j1, k1) * xd

! Interpolate along y
c0 = c00 * (1.-yd) + c10 * yd
c1 = c01 * (1.-yd) + c11 * yd

! Interpolate along z
interpolate3D = c0 * (1.-zd) + c1 * zd

return
end function interpolate3D

end module hit_inflow

