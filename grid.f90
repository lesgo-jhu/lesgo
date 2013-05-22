!!
!!  Copyright 2009,2010,2011,2012 Johns Hopkins University
!!
!!  Licensed under the Apache License, Version 2.0 (the "License"); you may not 
!!  use this file except in compliance with the License. You may obtain a copy of
!!  the License at:
!!
!!    http://www.apache.org/licenses/LICENSE-2.0
!!
!!  Unless required by applicable law or agreed to in writing, software 
!!  distributed under the License is distributed on an "AS IS" BASIS, WITHOUT 
!!  WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the 
!!  License for the specific language governing permissions and limitations under
!!  the License.
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
public grid, grid_build

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

end module grid_defs
