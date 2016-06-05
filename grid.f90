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


!**********************************************************************
module grid_defs
!**********************************************************************
use types, only : rprec
implicit none
save
private
!public x, y, z, zw, grid_build, grid_built
!public autowrap_i, autowrap_j
public grid, grid_build, grid_build_fourier

type grid_t
  logical :: built
  real(rprec), pointer, dimension(:) :: x, y, z, zw
  integer, pointer, dimension(:) :: autowrap_i, autowrap_j
end type grid_t

type(grid_t) :: grid
!real(rprec), allocatable, dimension(:) :: x, y, z, zw
! These need to be used in conjunction with modulo


contains

!**********************************************************************
subroutine grid_build()
!**********************************************************************
!
!  This subroutine creates the uv grid for the domain. It uses the x,y,z
!  variables decalared in grid_defs. This subroutine should only be called
!  once. To use x,y,z elsewhere in the code make sure
!  
!  use grid_defs, only : x,y,z 
!  
!  is placed in the routine
!  
use param, only : nx,ny,nz,jzmin,jzmax,dx,dy,dz,coord,lbz,nproc
implicit none

integer :: i,j,k
real(rprec), pointer, dimension(:) :: x,y,z,zw
integer, pointer, dimension(:) :: autowrap_i, autowrap_j

nullify(x,y,z,zw)
nullify(autowrap_i,autowrap_j)

!  x and y go to nx+1, ny+1 respectively for adding
!  the buffered points for periodicity
allocate(grid % x(nx+1),grid % y(ny+1))
allocate(grid % z(lbz:nz), grid % zw(lbz:nz))
allocate(grid % autowrap_i(0:nx+1), grid % autowrap_j(0:ny+1))

! Initialize built
grid % built = .false. 

! Set pointers
x => grid % x
y => grid % y
z => grid % z
zw => grid %zw

autowrap_i => grid % autowrap_i
autowrap_j => grid % autowrap_j

do k=lbz,nz
  $if ($MPI)
  z(k) = (coord*(nz-1) + k - 0.5_rprec) * dz
  $else
  z(k) = (k - 0.5_rprec) * dz
  $endif
enddo
do j=1,ny+1
  y(j) = (j-1)*dy
enddo
do i=1,nx+1
  x(i) = (i - 1)*dx
enddo
zw = z - dz/2._rprec

! Set index autowrapping arrays
autowrap_i(0) = nx
autowrap_j(0) = ny
autowrap_i(nx+1) = 1
autowrap_j(ny+1) = 1
do i=1,nx; autowrap_i(i) = i; enddo
do j=1,ny; autowrap_j(j) = j; enddo
     
grid % built = .true. 

nullify(x,y,z,zw)
nullify(autowrap_i,autowrap_j)

! Set jzmin and jzmax - the levels that this processor "owns"
$if($MPI)
  if (coord == 0) then
    jzmin = 0
    jzmax = nz-1
  elseif (coord == nproc-1) then
    jzmin = 1
    jzmax = nz
  else
    jzmin = 1
    jzmax = nz-1
  endif
$else
  jzmin = 1
  jzmax = nz
$endif

return
end subroutine grid_build 

!**********************************************************************
subroutine grid_build_fourier()
!**********************************************************************
!
!  This subroutine creates the uv grid for the domain. It uses the x,y,z
!  variables decalared in grid_defs. This subroutine should only be called
!  once. To use x,y,z elsewhere in the code make sure
!  
!  use grid_defs, only : x,y,z 
!  
!  is placed in the routine
!  
use param, only : nxp,ny,nz,jzmin,jzmax,dx,dy,dz,coord,lbz,nproc
implicit none

integer :: i,j,k
real(rprec), pointer, dimension(:) :: x,y,z,zw
integer, pointer, dimension(:) :: autowrap_i, autowrap_j

nullify(x,y,z,zw)
nullify(autowrap_i,autowrap_j)

!  x and y go to nx+1, ny+1 respectively for adding
!  the buffered points for periodicity
allocate(grid % x(nxp+1),grid % y(ny+1))
allocate(grid % z(lbz:nz), grid % zw(lbz:nz))
allocate(grid % autowrap_i(0:nxp+1), grid % autowrap_j(0:ny+1))

! Initialize built
grid % built = .false. 

! Set pointers
x => grid % x
y => grid % y
z => grid % z
zw => grid %zw

autowrap_i => grid % autowrap_i
autowrap_j => grid % autowrap_j

do k=lbz,nz
  $if ($MPI)
  z(k) = (coord*(nz-1) + k - 0.5_rprec) * dz
  $else
  z(k) = (k - 0.5_rprec) * dz
  $endif
enddo
do j=1,ny+1
  y(j) = (j-1)*dy
enddo
do i=1,nxp+1
  x(i) = (i - 1)*dx
enddo
zw = z - dz/2._rprec

! Set index autowrapping arrays
autowrap_i(0) = nxp
autowrap_j(0) = ny
autowrap_i(nxp+1) = 1
autowrap_j(ny+1) = 1
do i=1,nxp; autowrap_i(i) = i; enddo
do j=1,ny; autowrap_j(j) = j; enddo
     
grid % built = .true. 

nullify(x,y,z,zw)
nullify(autowrap_i,autowrap_j)

! Set jzmin and jzmax - the levels that this processor "owns"
$if($MPI)
  if (coord == 0) then
    jzmin = 0
    jzmax = nz-1
  elseif (coord == nproc-1) then
    jzmin = 1
    jzmax = nz
  else
    jzmin = 1
    jzmax = nz-1
  endif
$else
  jzmin = 1
  jzmax = nz
$endif

return
end subroutine grid_build_fourier


end module grid_defs
