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
module functions
!*******************************************************************************
use messages
use types, only : rprec
implicit none

save
private

public interp_to_uv_grid, trilinear_interp, trilinear_interp_w, binary_search, &
    bilinear_interp, linear_interp, cell_indx, buff_indx, interp_to_w_grid,    &
    get_tau_wall_bot, get_tau_wall_top, count_lines

character (*), parameter :: mod_name = 'functions'

interface linear_interp
    module procedure :: linear_interp_ss
    module procedure :: linear_interp_sa
    module procedure :: linear_interp_aa
end interface linear_interp

interface bilinear_interp
    module procedure :: bilinear_interp_ss
    module procedure :: bilinear_interp_sa
    module procedure :: bilinear_interp_aa
end interface bilinear_interp

contains

!*******************************************************************************
function interp_to_uv_grid(var, lbz) result(var_uv)
!*******************************************************************************
!  This function interpolates the array var, which resides on the w-grid,
!  onto the uv-grid variable var_uv using linear interpolation. It is
!  important to note that message passing is required for MPI cases and
!  all processors must call this routine. If this subroutine is call from a
!  only a subset of the total processors, the code will hang due to the usage
!  of the syncronous send/recv functions and certain processors waiting
!  to recv data but it never gets there.
!
!  NOTE: It is assumed that the size of var and var_uv are the same as the
!  coord (processor) domain and that k=nz-1 (k=0) and k=1 (k=nz) are overlap
!  nodes - no interpolation is performed for k=0 and k=nz
!
!  It is assumed that the array var has been synced if using MPI.

use types, only : rprec
use param, only : nz
use messages
#ifdef PPMPI
use param, only : coord, nproc
use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWN, MPI_SYNC_DOWNUP
#endif
implicit none

real(rprec), dimension(:,:,lbz:), intent(in) :: var
integer, intent(in) :: lbz
real(rprec), allocatable, dimension(:,:,:) :: var_uv

integer :: sx,sy,ubz

character (*), parameter :: sub_name = mod_name // '.interp_to_uv_grid'

sx = size(var,1)
sy = size(var,2)
ubz = ubound(var,3)

if( ubz .ne. nz ) call error(sub_name, 'Input array must lbz:nz z dimensions.')

allocate(var_uv(sx,sy,lbz:ubz))

!do k=1,ubz-1
!  do j=1,uby
!    do i=1,ubx
!      var_uv(i,j,k) = 0.5_rprec * (var(i,j,k+1) + var(i,j,k))
!    enddo
!  enddo
!enddo
! Perform the interpolation
var_uv(:,:,1:ubz-1) = 0.5_rprec * (var(:,:,2:ubz) + var(:,:,1:ubz-1))

#ifdef PPMPI
!  Take care of top "physical" boundary
if(coord == nproc - 1) var_uv(:,:,ubz) = var_uv(:,:,ubz-1)

!  Sync all overlapping data
if( lbz == 0 ) then
    call mpi_sync_real_array( var_uv, lbz, MPI_SYNC_DOWNUP )
elseif( lbz == 1 ) then
!     call mpi_sendrecv(var_uv(:,:,1), ubx*uby, MPI_RPREC, down, 1,  &
!                       var_uv(:,:,nz), ubx*uby, mpi_rprec, up, 1,   &
!                       comm, status, ierr)
    call mpi_sync_real_array( var_uv, lbz, MPI_SYNC_DOWN )
endif
#else
!  Take care of top "physical" boundary
var_uv(:,:,ubz) = var_uv(:,:,ubz-1)
#endif

!deallocate(var_uv)

!#ifdef PPMPI
!deallocate(buf)
!#endif

end function interp_to_uv_grid

!*******************************************************************************
function interp_to_w_grid(var, lbz) result(var_w)
!*******************************************************************************
!  This function interpolates the array var, which resides on the uv-grid,
!  onto the w-grid variable var_w using linear interpolation. It is
!  important to note that message passing is required for MPI cases and
!  all processors must call this routine. If this subroutine is call from a
!  only a subset of the total processors, the code will hang due to the usage
!  of the syncronous send/recv functions and certain processors waiting
!  to recv data but it never gets there.
!
!  NOTE: It is assumed that the size of var and var_w are the same as the
!  coord (processor) domain and that k=nz-1 (k=0) and k=1 (k=nz) are overlap
!  nodes - no interpolation is performed for k=0 and k=nz
!
!  ALSO: Bogus values at coord==0 and k==0 might lead to problems with the
!  k==1 interpolated value.  Velocities and sgs stresses (txx,txy,tyy,tzz)
!  are exactly zero at this point (k==1 on the w-grid), though, and can
!  therefore be set manually after this interpolation.

use types, only : rprec
use param, only : nz
use messages
#ifdef PPMPI
use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWN, MPI_SYNC_DOWNUP
#endif
implicit none

real(rprec), dimension(:,:,lbz:), intent(in) :: var
integer, intent(in) :: lbz
real(rprec), allocatable, dimension(:,:,:) :: var_w
integer :: sx,sy,ubz
!integer :: i,j,k

character (*), parameter :: sub_name = mod_name // '.interp_to_w_grid'

sx = size(var,1)
sy = size(var,2)
ubz = ubound(var,3)

if( ubz .ne. nz ) call error(sub_name, 'Input array must lbz:nz z dimensions.')

allocate(var_w(sx,sy,lbz:ubz))

! Perform the interpolation - does not work for lbz level
var_w(:,:,lbz+1:ubz) = 0.5_rprec * (var(:,:,lbz:ubz-1) + var(:,:,lbz+1:ubz))


#ifdef PPMPI
!  Sync all overlapping data
if( lbz == 0 ) then
    call mpi_sync_real_array( var_w, lbz, MPI_SYNC_DOWNUP )
elseif( lbz == 1 ) then
    call mpi_sync_real_array( var_w, lbz, MPI_SYNC_DOWN )
endif
#endif

end function interp_to_w_grid

!*******************************************************************************
integer function cell_indx_w(indx,dx,px)
!*******************************************************************************
! This routine takes index=['i' or 'j' or 'k'] and the magnitude of the
!   spacing=[dx or dy or dz] and the [x or y or z] location and returns
!   the value of the lower index (cell index). Also include is implicit
!   wrapping of the spatial location px. For z, the w-grid is used. To
!   use uv-grid for cell index, please see cell_indx() function.
!
!  cell_indx should always be:
!  1<= cell_indx <= Nx
!  or
!  1<= cell_indx <= Ny
!  or
!  lbz <= cell_indx < Nz
!
use types, only : rprec
use grid_m
use messages
use param, only : nx, ny, nz, L_x, L_y, L_z, lbz
implicit none

character (*), intent (in) :: indx
real(rprec), intent(in) :: dx
real(rprec) :: px ! Global value
character (*), parameter :: func_name = mod_name // '.cell_indx'
real(rprec), parameter :: thresh = 1.e-9_rprec
real(rprec), pointer, dimension(:) :: zw

! Nullify pointers
nullify(zw)
! Intialize result
cell_indx_w = -1

if(.not. grid % built) call grid%build()

zw => grid % zw

select case (indx)
    case ('i')
        ! Autowrap spatial point
        px = modulo(px,L_x)
        ! Check lower boundary
        if( abs(px) / L_x < thresh ) then
            cell_indx_w = 1
            ! Check upper boundary
        elseif( abs( px - L_x ) / L_x < thresh ) then
            cell_indx_w = Nx
        else
            ! Returned values 1 < cell_indx < Nx
            cell_indx_w = floor (px / dx) + 1
        endif
    case ('j')
        ! Autowrap spatial point
        px = modulo(px, L_y)

        ! Check lower boundary
        if( abs(px) / L_y < thresh ) then
            cell_indx_w = 1
        ! Check upper boundary
        elseif( abs( px - L_y ) / L_y < thresh ) then
            cell_indx_w = Ny
        else
            ! Returned values 1 < cell_indx < Ny
            cell_indx_w = floor (px / dx) + 1
        endif
    !  Need to compute local distance to get local k index
    case ('k')
        ! Check upper boundary
        if( abs( px - zw(Nz) ) / L_z < thresh ) then
            cell_indx_w = Nz-1
        else
            cell_indx_w = floor ((px - zw(1)) / dx) + 1
        endif
end select

nullify(zw)

end function cell_indx_w


!*******************************************************************************
integer function cell_indx(indx,dx,px)
!*******************************************************************************
! This routine takes index=['i' or 'j' or 'k'] and the magnitude of the
!   spacing=[dx or dy or dz] and the [x or y or z] location and returns
!   the value of the lower index (cell index). Also include is implicit
!   wrapping of the spatial location px
!
!  cell_indx should always be:
!  1<= cell_indx <= Nx
!  or
!  1<= cell_indx <= Ny
!  or
!  lbz <= cell_indx < Nz
!
use types, only : rprec
use grid_m
use messages
use param, only : nx, ny, nz, L_x, L_y, L_z, lbz
implicit none

character (*), intent (in) :: indx
real(rprec), intent(in) :: dx
real(rprec) :: px ! Global value
character (*), parameter :: func_name = mod_name // '.cell_indx'
real(rprec), parameter :: thresh = 1.e-9_rprec
real(rprec), pointer, dimension(:) :: z

! Nullify pointers
nullify(z)
! Intialize result
cell_indx = -1

if(.not. grid % built) call grid%build()

z => grid % z

select case (indx)
    case ('i')
        ! Autowrap spatial point
        px = modulo(px,L_x)

        ! Check lower boundary
        if( abs(px) / L_x < thresh ) then
            cell_indx = 1
        ! Check upper boundary
        elseif( abs( px - L_x ) / L_x < thresh ) then
            cell_indx = Nx
        else
            ! Returned values 1 < cell_indx < Nx
            cell_indx = floor (px / dx) + 1
        endif
    case ('j')
        ! Autowrap spatial point
        px = modulo(px, L_y)

        ! Check lower boundary
        if( abs(px) / L_y < thresh ) then
            cell_indx = 1
        ! Check upper boundary
        elseif( abs( px - L_y ) / L_y < thresh ) then
            cell_indx = Ny
        else
            ! Returned values 1 < cell_indx < Ny
            cell_indx = floor (px / dx) + 1
        endif
    !  Need to compute local distance to get local k index
    case ('k')
        ! Check upper boundary
        if( abs( px - z(Nz) ) / L_z < thresh ) then
            cell_indx = Nz-1
        else
            cell_indx = floor ((px - z(1)) / dx) + 1
        endif
end select

nullify(z)

end function cell_indx

!*******************************************************************************
real(rprec) function trilinear_interp_w(var,lbz,xyz)
!*******************************************************************************
!
!  This subroutine perform trilinear interpolation for a point that
!  exists in the cell with lower dimension (cell index) : istart,jstart,kstart
!  for the point xyz
!
!  istart, jstart, kstart are used to determine the cell location on the
!  w-grid (with exceptions, see below); these are defined in output_init
!
!  This assumes that var is on the w-grid, except (jz=1 .and. coord==0)
!  .or. (jz=nz .and. coord==nproc-1) when these are near wall BC, in which case
!  the var is assumed to be on the uv-grid. This oddity is because this
!  function is designed for interpolating variables for the lagrangian subgrid
!  models and for wall BC lesgo stores nu_t not at the wall but 0.5*dz off it.
!
!  The variable sent to this subroutine should have a lower-bound-on-z
!  (lbz) set as an input so the k-index will match the k-index of z.
!  Before calling this function, make sure the point exists on the coord
!  [ test using: z(1) \leq z_p < z(nz-1) ]
!
use grid_m
use types, only : rprec
use param, only : dx, dy, dz, coord, nproc, lbc_mom, ubc_mom, nz
implicit none

integer, intent(in) :: lbz
real(rprec), dimension(:,:,lbz:), intent(in) :: var
integer :: istart, jstart, kstart, istart1, jstart1, kstart1
real(rprec), intent(in), dimension(3) :: xyz
!integer, parameter :: nvar = 3
real(rprec) :: u1, u2, u3, u4, u5, u6
real(rprec) :: xdiff, ydiff, zdiff
real(rprec), pointer, dimension(:) :: x, y, z, zw
integer, pointer, dimension(:) :: autowrap_i, autowrap_j
nullify(x,y,z,zw)
nullify(autowrap_i, autowrap_j)

x  => grid % x
y  => grid % y
z  => grid % z
zw => grid % zw
autowrap_i => grid % autowrap_i
autowrap_j => grid % autowrap_j

!  Initialize stuff
u1 = 0._rprec; u2 = 0._rprec; u3 = 0._rprec;
u4 = 0._rprec; u5 = 0._rprec; u6 = 0._rprec

! X and Y index not affected by z-grid details
! Determine istart, jstart by calling cell_indx
istart = cell_indx_w('i',dx,xyz(1)) ! 1<= istart <= Nx
jstart = cell_indx_w('j',dy,xyz(2)) ! 1<= jstart <= Ny

! Set +1 values
istart1 = autowrap_i(istart+1) ! Autowrap index
jstart1 = autowrap_j(jstart+1) ! Autowrap index

!  Compute xdiff
xdiff = xyz(1) - x(istart)
!  Compute ydiff
ydiff = xyz(2) - y(jstart)


if(coord==0 .and. lbc_mom>0 .and. xyz(3) < zw(2)) then
    ! special case
    if (xyz(3) < z(1)) then
        ! below first point, zero grad BC
        kstart  = 1
        kstart1 = 1 ! use same z-index for both sides of interp (zero grad)
        zdiff = 0._rprec  ! kstart==kstart1, so u1==u3 .and. u2==u4, so u5==u6
        ! therefore, zdiff is irrelevant in this case
    else
        ! in-between jz=1 and jz=2, which are sep by only 0.5*dz
        kstart  = 1
        kstart1 = 2
        zdiff   = 2._rprec*(xyz(3) - z(kstart)) ! use uvp-grid for z
        ! multiply by 2 to achieve effective dz/2 in last step of interp
    end if
else if (coord==nproc-1 .and. ubc_mom>0 .and. xyz(3) > zw(nz-1)) then
    ! special case
    if (xyz(3) > z(nz-1)) then
        ! above last point, zero grad BC
        kstart  = nz
        kstart1 = nz ! use same z-index for both sides of interp (zero grad)
        zdiff = 0._rprec  ! kstart==kstart1, so u1==u3 .and. u2==u4, so u5==u6
        ! therefore, zdiff is irrelevant in this case
    else
        ! in-between jz=nz-1 and jz=nz, which are sep by only 0.5*dz
        kstart  = nz-1
        kstart1 = nz
        zdiff   = 2._rprec*(xyz(3) - zw(kstart)) ! use w-grid for z
        ! multiply by 2 to achieve effective dz/2 in last step of interp
    end if
else
    ! Determine kstart by calling cell_indx_w
    kstart = cell_indx_w('k',dz,xyz(3)) ! lbz <= kstart < Nz

    ! Set +1 values
    kstart1 = kstart + 1

    !  Compute zdiff
    zdiff = xyz(3) - zw(kstart)
end if

!  Perform the 7 linear interpolations
!  Perform interpolations in x-direction
u1 = var(istart,  jstart,  kstart)  + (xdiff) * (var(istart1, jstart,  kstart) &
    - var(istart,  jstart,  kstart)) / dx
u2 = var(istart,  jstart1, kstart)  + (xdiff) * (var(istart1, jstart1, kstart) &
    - var(istart,  jstart1, kstart)) / dx
u3 = var(istart,  jstart,  kstart1) + (xdiff) * (var(istart1, jstart,  kstart1)&
    - var(istart,  jstart,  kstart1)) / dx
u4 = var(istart,  jstart1, kstart1) + (xdiff) * (var(istart1, jstart1, kstart1)&
    - var(istart,  jstart1, kstart1)) / dx

!  Perform interpolations in y-direction
u5 = u1 + (ydiff) * (u2 - u1) / dy
u6 = u3 + (ydiff) * (u4 - u3) / dy
!  Perform interpolation in z-direction
trilinear_interp_w = u5 + (zdiff) * (u6 - u5) / dz

nullify(x,y,z,zw)
nullify(autowrap_i, autowrap_j)

end function trilinear_interp_w

!*******************************************************************************
real(rprec) function trilinear_interp(var,lbz,xyz)
!*******************************************************************************
!
!  This subroutine perform trilinear interpolation for a point that
!  exists in the cell with lower dimension (cell index) : istart,jstart,kstart
!  for the point xyz
!
!  istart, jstart, kstart are used to determine the cell location on the
!  uv-grid; these are defined in output_init
!
!  Takes care of putting w-grid variables onto the uv-grid; this assumes
!  that var is on the uv-grid
!
!  The variable sent to this subroutine should have a lower-bound-on-z
!  (lbz) set as an input so the k-index will match the k-index of z.
!  Before calling this function, make sure the point exists on the coord
!  [ test using: z(1) \leq z_p < z(nz-1) ]
!
use grid_m
use types, only : rprec
use param, only : dx, dy, dz
implicit none

integer, intent(in) :: lbz
real(rprec), dimension(:,:,lbz:), intent(in) :: var
integer :: istart, jstart, kstart, istart1, jstart1, kstart1
real(rprec), intent(in), dimension(3) :: xyz
!integer, parameter :: nvar = 3
real(rprec) :: u1, u2, u3, u4, u5, u6
real(rprec) :: xdiff, ydiff, zdiff
real(rprec), pointer, dimension(:) :: x, y, z
integer, pointer, dimension(:) :: autowrap_i, autowrap_j

nullify(x,y,z)
nullify(autowrap_i, autowrap_j)

x => grid % x
y => grid % y
z => grid % z
autowrap_i => grid % autowrap_i
autowrap_j => grid % autowrap_j

!  Initialize stuff
u1 = 0._rprec; u2 = 0._rprec; u3 = 0._rprec; u4 = 0._rprec; u5 = 0._rprec; u6=0.

! Determine istart, jstart, kstart by calling cell_indx
istart = cell_indx('i',dx,xyz(1)) ! 1<= istart <= Nx
jstart = cell_indx('j',dy,xyz(2)) ! 1<= jstart <= Ny
kstart = cell_indx('k',dz,xyz(3)) ! lbz <= kstart < Nz

! Extra term with kstart accounts for shift in var k-index if lbz.ne.1
! Set +1 values
istart1 = autowrap_i(istart+1) ! Autowrap index
jstart1 = autowrap_j(jstart+1) ! Autowrap index
kstart1 = kstart + 1

!  Compute xdiff
xdiff = xyz(1) - x(istart)
!  Compute ydiff
ydiff = xyz(2) - y(jstart)
!  Compute zdiff
zdiff = xyz(3) - z(kstart)

!  Perform the 7 linear interpolations
!  Perform interpolations in x-direction
u1 = var(istart,  jstart,  kstart)  + (xdiff) * (var(istart1, jstart,  kstart) &
    - var(istart,  jstart,  kstart)) / dx
u2 = var(istart,  jstart1, kstart)  + (xdiff) * (var(istart1, jstart1, kstart) &
    - var(istart,  jstart1, kstart)) / dx
u3 = var(istart,  jstart,  kstart1) + (xdiff) * (var(istart1, jstart,  kstart1)&
    - var(istart,  jstart,  kstart1)) / dx
u4 = var(istart,  jstart1, kstart1) + (xdiff) * (var(istart1, jstart1, kstart1)&
    - var(istart,  jstart1, kstart1)) / dx
!  Perform interpolations in y-direction
u5 = u1 + (ydiff) * (u2 - u1) / dy
u6 = u3 + (ydiff) * (u4 - u3) / dy
!  Perform interpolation in z-direction
trilinear_interp = u5 + (zdiff) * (u6 - u5) / dz

nullify(x,y,z)
nullify(autowrap_i, autowrap_j)

end function trilinear_interp

!*******************************************************************************
real(rprec) function bilinear_interp_ss(u11,u21,u12,u22,dx,dy,xdiff,ydiff)
!*******************************************************************************
!
!  This function performs linear interpolation
!
!  Inputs:
!  u11          - lower bound value in x direction for lower y
!  u21          - upper bound value in x direction for lower y
!  u12          - lower bound value in x direction for upper y
!  u22          - upper bound value in x direction for upper y
!  dx           - length delta for the grid in x direction
!  dy           - length delta for the grid in y direction
!  xdiff        - distance from the point of interest to the u11 node in x direction
!  xdiff        - distance from the point of interest to the u11 node in y direction
!
use types, only : rprec
implicit none

real(rprec), intent(in) :: u11, u12, u21, u22, dx, dy, xdiff, ydiff
real(rprec) :: v1, v2

v1 = linear_interp(u11, u21, dx, xdiff)
v2 = linear_interp(u12, u22, dx, xdiff)

bilinear_interp_ss = linear_interp(v1,v2,dy,ydiff)

end function bilinear_interp_ss

!*******************************************************************************
function bilinear_interp_sa_nocheck(x, y, v, xq, yq) result(vq)
!*******************************************************************************
!
!  This function performs the linear interpolation fo bilinear_interp_sa
!  and bilinear_interp_aa without checking bounds of the input arrays.
!  It should not be called directly.
!
implicit none
real(rprec), dimension(:), intent(in) :: x, y
real(rprec), dimension(:,:), intent(in) :: v
real(rprec), intent(in) :: xq, yq
real(rprec) :: vq
integer     :: i, j, Nx, Ny

Nx = size(x)
Ny = size(y)
i = binary_search(x, xq)
if (i == 0) then
    vq = linear_interp(y, v(1,:), yq)
else if (i == Nx) then
    vq = linear_interp(y, v(Nx,:), yq)
else
    j = binary_search(y, yq)
    if (j == 0) then
        vq = linear_interp(v(i,1), v(i+1,1), x(i+1)-x(i), xq - x(i))
    else if (j == Ny) then
        vq = linear_interp(v(i,Ny), v(i+1,Ny), x(i+1)-x(i), xq - x(i))
    else
        vq = bilinear_interp_ss( v(i,j), v(i+1,j), v(i,j+1), v(i+1,j+1), &
             x(i+1) - x(i), y(j+1) - y(j), xq - x(i), yq - y(j) )
    end if
end if

end function bilinear_interp_sa_nocheck

!*******************************************************************************
function bilinear_interp_sa(x, y, v, xq, yq) result(vq)
!*******************************************************************************
!
!  This function performs linear interpolation from a set of points
!  defined on a grid (x,y) with values v to a query point (xq, yq)
!
!  Inputs:
!  x            - array of sample points
!  v            - array of values at sample points
!  xq           - query point
!
implicit none
real(rprec), dimension(:), intent(in) :: x, y
real(rprec), dimension(:,:), intent(in) :: v
real(rprec), intent(in) :: xq, yq
real(rprec) :: vq
character(*), parameter :: func_name = mod_name // '.bilinear_interp_sa'

if ( size(v,1) /= size(x) .or. size(v,2) /= size(y)) then
    call error(func_name, 'Array v must have a size of [size(x), size(y)]')
end if

vq = bilinear_interp_sa_nocheck(x, y, v, xq, yq)

end function bilinear_interp_sa

!*******************************************************************************
function bilinear_interp_aa(x, y, v, xq, yq) result(vq)
!*******************************************************************************
!
!  This function performs linear interpolation from a set of points
!  defined on a grid (x,y) with values v to an array of query points
!  (xq, yq)
!
!  Inputs:
!  x            - array of sample points
!  v            - array of values at sample points
!  xq           - array of query points
!
implicit none
real(rprec), dimension(:), intent(in) :: x, y
real(rprec), dimension(:,:), intent(in) :: v
real(rprec), dimension(:), intent(in) :: xq, yq
real(rprec), dimension(:), allocatable :: vq
integer :: i, N
character(*), parameter :: func_name = mod_name // '.bilinear_interp_sa'

if ( size(v,1) /= size(x) .or. size(v,2) /= size(y)) then
    call error(func_name, 'Array v must have a size of [size(x), size(y)]')
end if
if ( size(xq) /= size(yq) ) then
    call error(func_name, 'Arrays xq and yq must be the same size')
end if

N = size(xq)
allocate(vq(N))

do i = 1, N
    vq(i) = bilinear_interp_sa_nocheck(x, y, v, xq(i), yq(i))
end do

end function bilinear_interp_aa

!*******************************************************************************
real(rprec) function linear_interp_ss(u1,u2,dx,xdiff)
!*******************************************************************************
!
!  This function performs linear interpolation
!
!  Inputs:
!  u1           - lower bound value in the increasing index direction
!  u2           - upper bound value in the increasing index direction
!  dx           - length delta for the grid in the correct direction
!  xdiff        - distance from the point of interest to the u1 node
!
use types, only : rprec
implicit none

real(rprec), intent(in) :: u1, u2, dx, xdiff

linear_interp_ss = u1 + (xdiff) * (u2 - u1) / dx

end function linear_interp_ss

!*******************************************************************************
function linear_interp_sa_nocheck(x, v, xq) result(vq)
!*******************************************************************************
!
!  This function performs the linear interpolation fo linear_interp_sa
!  and linear_interp_aa without checking bounds of the input arrays.
!  It should not be called directly.
!
implicit none
real(rprec), dimension(:), intent(in) :: x, v
real(rprec), intent(in) :: xq
real(rprec) :: vq
integer :: i, N

N = size(v)
i = binary_search(x, xq)
if (i == 0) then
    vq = v(1)
else if (i == N) then
    vq = v(N)
else
    vq = linear_interp_ss(v(i), v(i+1), x(i+1)-x(i), (xq - x(i)))
end if

end function linear_interp_sa_nocheck

!*******************************************************************************
function linear_interp_sa(x, v, xq) result(vq)
!*******************************************************************************
!
!  This function performs linear interpolation from a set of points x
!  with values v to a query point xq
!
!  Inputs:
!  x            - array of sample points
!  v            - array of values at sample points
!  xq           - query point
!
implicit none
real(rprec), dimension(:), intent(in) :: x, v
real(rprec), intent(in) :: xq
real(rprec) :: vq
character(*), parameter :: func_name = mod_name // '.linear_interp_sa'

if ( size(v) /= size(x) ) then
    call error(func_name, 'Arrays x and v must be the same size')
end if

vq = linear_interp_sa_nocheck(x, v, xq)

end function linear_interp_sa

!*******************************************************************************
function linear_interp_aa(x, v, xq) result(vq)
!*******************************************************************************
!
!  This function performs linear interpolation from a set of points x
!  with values v to an array of query points xq
!
!  Inputs:
!  x            - array of sample points
!  v            - array of values at sample points
!  xq           - array of query points
!
implicit none
real(rprec), dimension(:), intent(in) :: x, v
real(rprec), dimension(:), intent(in) :: xq
real(rprec), dimension(:), allocatable :: vq
integer :: i, N
character(*), parameter :: func_name = mod_name // '.linear_interp_aa'

! Check array sizes
if ( size(v) /= size(x) ) then
    call error(func_name, 'Arrays x and v must be the same size')
end if

! Allocate output
N = size(xq)
allocate(vq(N))

! For each element of the array perform interpolation
do i = 1, N
    vq(i) = linear_interp_sa_nocheck(x, v, xq(i))
end do

end function linear_interp_aa

!*******************************************************************************
function binary_search(arr,val) result(low)
!*******************************************************************************
!
!  This function performs a binary search on a sorted array. Given the
!  provided value, adjacent low and high indices of the array are found
!  such that the provided value is bracketed. Guaranteed log2(N) search.
!
!  Inputs:
!  arr          - sorted array of values to search
!  val          - value to be bracketed
!
!  Output:
!  low          - lower index of the array bracket
!                 0 if val < arr(1), N if val < arr(N))
!
implicit none

real(rprec), dimension(:) :: arr
real(rprec) :: val
integer :: low, mid, high, N

! Size of array
N = size(arr)

! Check if value is outside bounds
if ( val < arr(1) ) then
    low = 0
    return
end if
if ( val > arr(N) ) then
    low = N
    return
end if

! Otherwise perform bisection
low = 1
high = N
do while (high - low > 1)
    mid = (low + high) / 2
    if ( arr(mid) > val ) then
        high = mid
    elseif ( arr(mid) < val ) then
        low = mid
    else
        low = mid
        return
    endif
end do

end function binary_search

!*******************************************************************************
integer function buff_indx(i,imax)
!*******************************************************************************
!  This function returns the physical index associated with the buffer
!  region for the specified i and imax.
!  For i = imax + 1 -> 1 is returned otherwise i is returned
implicit none

integer, intent(in) :: i,imax

if (i == imax + 1) then
    buff_indx = 1
else
    buff_indx = i
endif

return
end function buff_indx

!*******************************************************************************
function get_tau_wall_bot() result(twall)
!*******************************************************************************
!
! This function provides plane-averaged value of wall stress magnitude
use types, only: rprec
use param, only : nx, ny
use sim_param, only : txz, tyz

implicit none
real(rprec) :: twall, txsum, tysum
integer :: jx, jy

txsum = 0._rprec
tysum = 0._rprec
do jx = 1, nx
do jy = 1, ny
    txsum = txsum + txz(jx,jy,1)
    tysum = tysum + tyz(jx,jy,1)
enddo
enddo

twall = sqrt( (txsum/(nx*ny))**2 + (tysum/(nx*ny))**2  )

end function get_tau_wall_bot


!*******************************************************************************
function get_tau_wall_top() result(twall)
!*******************************************************************************
!
! This function provides plane-averaged value of wall stress magnitude
use types, only: rprec
use param, only : nx, ny, nz
use sim_param, only : txz, tyz

implicit none
real(rprec) :: twall, txsum, tysum
integer :: jx, jy

txsum = 0._rprec
tysum = 0._rprec
do jx = 1, nx
do jy = 1, ny
    txsum = txsum + txz(jx,jy,nz)
    tysum = tysum + tyz(jx,jy,nz)
enddo
enddo

twall = sqrt( (txsum/(nx*ny))**2 + (tysum/(nx*ny))**2  )

end function get_tau_wall_top

!*******************************************************************************
function count_lines(fname) result(N)
!*******************************************************************************
!
! This function counts the number of lines in a file
!
use messages
use param, only : CHAR_BUFF_LENGTH
implicit none
character(*), intent(in) :: fname
logical :: exst
integer :: fid, ios
integer :: N

character(*), parameter :: sub_name = mod_name // '.count_lines'

! Check if file exists and open
inquire (file = trim(fname), exist = exst)
if (.not. exst) then
    call error (sub_name, 'file ' // trim(fname) // 'does not exist')
end if
open(newunit=fid, file=trim(fname), status='unknown', form='formatted',        &
    position='rewind')

! count number of lines and close
ios = 0
N = 0
do
    read(fid, *, IOstat = ios)
    if (ios /= 0) exit
    N = N + 1
end do

! Close file
close(fid)

end function count_lines

end module functions
