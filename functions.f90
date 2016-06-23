!!
!!  Copyright (C) 2009-2016  Johns Hopkins University
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
module functions
!**********************************************************************
use messages
use types, only : rprec
implicit none
save
private
public interp_to_uv_grid,   &
    trilinear_interp,       &
    bilinear_interp,        &
    linear_interp,          &
    cell_indx,              &
    buff_indx,              &
    points_avg_3d,          & 
    plane_avg_3d,           &     
    interp_to_w_grid,       &
    get_tau_wall

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

!**********************************************************************
function interp_to_uv_grid(var, lbz) result(var_uv)
!**********************************************************************
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
use param,only : nz
use messages
#ifdef PPMPI
use param, only : coord, nproc, MPI_RPREC, down, up,  comm, status, ierr
use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWN, MPI_SYNC_DOWNUP
#endif

implicit none

real(rprec), dimension(:,:,lbz:), intent(in) :: var
integer, intent(in) :: lbz
real(rprec), allocatable, dimension(:,:,:) :: var_uv

integer :: sx,sy,ubz

character (*), parameter :: sub_name = mod_name // '.interp_to_uv_grid'

sx=size(var,1)
sy=size(var,2)
ubz=ubound(var,3)

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
  ! call mpi_sendrecv(var_uv(:,:,1), ubx*uby, MPI_RPREC, down, 1,  &
  !                   var_uv(:,:,nz), ubx*uby, mpi_rprec, up, 1,   &
  !                   comm, status, ierr)
  call mpi_sync_real_array( var_uv, lbz, MPI_SYNC_DOWN )
endif                    

#else

!  Take care of top "physical" boundary
var_uv(:,:,ubz) = var_uv(:,:,ubz-1)

#endif
  
return 

!deallocate(var_uv)

!#ifdef PPMPI
!deallocate(buf)
!#endif

end function interp_to_uv_grid

!**********************************************************************
function interp_to_w_grid(var, lbz) result(var_w)
!**********************************************************************
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
use param,only : nz
use messages
#ifdef PPMPI
use param, only : coord, nproc, MPI_RPREC, down, up,  comm, status, ierr
use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWN, MPI_SYNC_DOWNUP
#endif

implicit none

real(rprec), dimension(:,:,lbz:), intent(in) :: var
integer, intent(in) :: lbz
real(rprec), allocatable, dimension(:,:,:) :: var_w
integer :: sx,sy,ubz
!integer :: i,j,k

character (*), parameter :: sub_name = mod_name // '.interp_to_w_grid'

sx=size(var,1)
sy=size(var,2)
ubz=ubound(var,3)

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
  
return 


end function interp_to_w_grid

!**********************************************************************
integer function cell_indx(indx,dx,px)
!**********************************************************************
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

return
end function cell_indx

!**********************************************************************
real(rprec) function trilinear_interp(var,lbz,xyz)
!**********************************************************************
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
use sim_param, only : u,v
use param, only : nx, ny, nz, dx, dy, dz, coord, L_x, L_y
implicit none

integer, intent(in) :: lbz
real(rprec), dimension(:,:,lbz:), intent(in) :: var
integer :: istart, jstart, kstart, istart1, jstart1, kstart1
real(rprec), intent(in), dimension(3) :: xyz

!integer, parameter :: nvar = 3
real(rprec) :: u1,u2,u3,u4,u5,u6
real(rprec) :: xdiff, ydiff, zdiff

real(rprec), pointer, dimension(:) :: x,y,z
integer, pointer, dimension(:) :: autowrap_i, autowrap_j

nullify(x,y,z)
nullify(autowrap_i, autowrap_j)

x => grid % x
y => grid % y
z => grid % z
autowrap_i => grid % autowrap_i
autowrap_j => grid % autowrap_j

!  Initialize stuff
u1=0.; u2=0.; u3=0.; u4=0.; u5=0.; u6=0.

! Determine istart, jstart, kstart by calling cell_indx
istart = cell_indx('i',dx,xyz(1)) ! 1<= istart <= Nx
jstart = cell_indx('j',dy,xyz(2)) ! 1<= jstart <= Ny
kstart = cell_indx('k',dz,xyz(3)) ! lbz <= kstart < Nz

! Extra term with kstart accounts for shift in var k-index if lbz.ne.1
! Set +1 values
istart1 = autowrap_i(istart+1) ! Autowrap index
jstart1 = autowrap_j(jstart+1) ! Autowrap index
kstart1 = kstart + 1


! Can probably bypass storing in uvar and put directly in linear_interp
!  Contains the 6 points that make of the cube
!uvar = 0.
!uvar(1,1,1) = var(istart,  jstart,  kstart)
!uvar(2,1,1) = var(istart1, jstart,  kstart)
!uvar(1,2,1) = var(istart,  jstart1, kstart)
!uvar(2,2,1) = var(istart1, jstart1, kstart)
!uvar(1,1,2) = var(istart,  jstart,  kstart1)
!uvar(2,1,2) = var(istart1, jstart,  kstart1)
!uvar(1,2,2) = var(istart,  jstart1, kstart1)
!uvar(2,2,2) = var(istart1, jstart1, kstart1)

!  Compute xdiff
xdiff = xyz(1) - x(istart)
!  Compute ydiff
ydiff = xyz(2) - y(jstart)
!  Compute zdiff
zdiff = xyz(3) - z(kstart)

!  Perform the 7 linear interpolations
!  Perform interpolations in x-direction 
!u1 = linear_interp(uvar(1,1,1),uvar(2,1,1),dx,xdiff)
!u2 = linear_interp(uvar(1,2,1),uvar(2,2,1),dx,xdiff)
!u3 = linear_interp(uvar(1,1,2),uvar(2,1,2),dx,xdiff)
!u4 = linear_interp(uvar(1,2,2),uvar(2,2,2),dx,xdiff)
!u1 = linear_interp(var(istart,  jstart,  kstart), &
!                   var(istart1, jstart,  kstart), &
!                   dx, xdiff)
!u2 = linear_interp(var(istart,  jstart1, kstart), &
!                   var(istart1, jstart1, kstart), &
!                   dx, xdiff)
!u3 = linear_interp(var(istart,  jstart,  kstart1), &
!                   var(istart1, jstart,  kstart1), &
!                   dx, xdiff)
!u4 = linear_interp(var(istart,  jstart1, kstart1), &
!                   var(istart1, jstart1, kstart1), &
!                   dx, xdiff)
u1=var(istart,  jstart,  kstart)  + (xdiff) * (var(istart1, jstart,  kstart)  - var(istart,  jstart,  kstart)) / dx
u2=var(istart,  jstart1, kstart)  + (xdiff) * (var(istart1, jstart1, kstart)  - var(istart,  jstart1, kstart)) / dx
u3=var(istart,  jstart,  kstart1) + (xdiff) * (var(istart1, jstart,  kstart1) - var(istart,  jstart,  kstart1)) / dx
u4=var(istart,  jstart1, kstart1) + (xdiff) * (var(istart1, jstart1, kstart1) - var(istart,  jstart1, kstart1)) / dx

!  Perform interpolations in y-direction
!u5 = linear_interp(u1,u2,dy,ydiff)
!u6 = linear_interp(u3,u4,dy,ydiff)
u5=u1 + (ydiff) * (u2 - u1) / dy
u6=u3 + (ydiff) * (u4 - u3) / dy
!  Perform interpolation in z-direction
!trilinear_interp = linear_interp(u5,u6,dz,zdiff)
trilinear_interp = u5 + (zdiff) * (u6 - u5) / dz

nullify(x,y,z)
nullify(autowrap_i, autowrap_j)

return
end function trilinear_interp

!~ !**********************************************************************
!~ real(rprec) function trilinear_interp(var,lbz,xyz)
!~ !**********************************************************************
!~ !  This subroutine perform trilinear interpolation for a point that
!~ !  exists in the cell with lower dimension (cell index) : istart,jstart,kstart
!~ !  for the point xyz
!~ !  
!~ !  istart, jstart, kstart are used to determine the cell location on the
!~ !  uv-grid; these are defined in output_init
!~ !
!~ !  Takes care of putting w-grid variables onto the uv-grid; this assumes
!~ !  that var is on the uv-grid
!~ !
!~ !  The variable sent to this subroutine should have a lower-bound-on-z 
!~ !  (lbz) set as an input so the k-index will match the k-index of z.  
!~ !  Before calling this function, make sure the point exists on the coord
!~ !  [ test using: z(1) \leq z_p < z(nz-1) ]
!~ use grid_m
!~ use param, only : dx,dy,dz,L_x,L_y,L_z, nz
!~ implicit none
!~ real(rprec), dimension(:,:,lbz:), intent(in) :: var
!~ integer    , intent(in) :: lbz
!~ real(rprec), intent(in), dimension(3) :: xyz
!~ real(rprec), pointer   , dimension(:) :: x,y,z
!~ integer, pointer, dimension(:) :: autowrap_i, autowrap_j
!~ real(rprec) :: u1,u2,u3,u4,u5,u6,xdiff,ydiff,zdiff,px,py
!~ integer :: istart,jstart,kstart,istart1,jstart1,kstart1
!~ 
!~ nullify(x,y,z,autowrap_i, autowrap_j)
!~ x => grid % x
!~ y => grid % y
!~ z => grid % z
!~ autowrap_i => grid % autowrap_i
!~ autowrap_j => grid % autowrap_j
!~ 
!~ ! Determine istart, jstart, kstart by calling cell_indx
!~ px = modulo(xyz(1),L_x)
!~ istart = floor (px / dx) + 1
!~ py = modulo(xyz(2),L_y)
!~ jstart = floor (py / dy) + 1
!~ if( abs( xyz(3) - z(nz) ) / L_z < 1.e-9 ) then
!~ kstart = nz-1
!~ else
!~ kstart = floor ((xyz(3) - z(1)) / dz) + 1
!~ endif
!~ 
!~ ! Extra term with kstart accounts for shift in var k-index if lbz.ne.1
!~ ! Set +1 values
!~ istart1 = autowrap_i(istart+1) ! Autowrap index
!~ jstart1 = autowrap_j(jstart+1) ! Autowrap index
!~ kstart1 = kstart + 1
!~ 
!~ !  Compute xdiff
!~ xdiff = px - x(istart)
!~ !  Compute ydiff
!~ ydiff = py - y(jstart)
!~ !  Compute zdiff
!~ zdiff = xyz(3) - z(kstart)
!~ 
!~ !  Perform the 7 linear interpolations
!~ !  Perform interpolations in x-direction 
!~ u1=var(istart,  jstart,  kstart)  + (xdiff) * (var(istart1, jstart,  kstart)  - var(istart,  jstart,  kstart)) / dx
!~ u2=var(istart,  jstart1, kstart)  + (xdiff) * (var(istart1, jstart1, kstart)  - var(istart,  jstart1, kstart)) / dx
!~ u3=var(istart,  jstart,  kstart1) + (xdiff) * (var(istart1, jstart,  kstart1) - var(istart,  jstart,  kstart1)) / dx
!~ u4=var(istart,  jstart1, kstart1) + (xdiff) * (var(istart1, jstart1, kstart1) - var(istart,  jstart1, kstart1)) / dx
!~ 
!~ !  Perform interpolations in y-direction
!~ u5=u1 + (ydiff) * (u2 - u1) / dy
!~ u6=u3 + (ydiff) * (u4 - u3) / dy
!~ !  Perform interpolation in z-direction
!~ trilinear_interp = u5 + (zdiff) * (u6 - u5) / dz
!~ 
!~ nullify(x,y,z,autowrap_i, autowrap_j)
!~ 
!~ return
!~ end function trilinear_interp

!**********************************************************************
real(rprec) function bilinear_interp_ss(u11,u21,u12,u22,dx,dy,xdiff,ydiff)
!**********************************************************************
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

return
end function bilinear_interp_ss

!**********************************************************************
function bilinear_interp_sa_nocheck(x, y, v, xq, yq) result(vq)
!**********************************************************************
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
Ny = size(x)
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

!**********************************************************************
function bilinear_interp_sa(x, y, v, xq, yq) result(vq)
!**********************************************************************
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

!**********************************************************************
function bilinear_interp_aa(x, y, v, xq, yq) result(vq)
!**********************************************************************
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

!**********************************************************************
real(rprec) function linear_interp_ss(u1,u2,dx,xdiff)
!**********************************************************************
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

return
end function linear_interp_ss

!**********************************************************************
function linear_interp_sa_nocheck(x, v, xq) result(vq)
!**********************************************************************
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

!**********************************************************************
function linear_interp_sa(x, v, xq) result(vq)
!**********************************************************************
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

!**********************************************************************
function linear_interp_aa(x, v, xq) result(vq)
!**********************************************************************
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

!**********************************************************************
function binary_search(arr,val) result(low)
!**********************************************************************
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

return
end function binary_search

!**********************************************************************
real(rprec) function plane_avg_3d(var, lbz, bp1, bp2, bp3, nzeta, neta)
!**********************************************************************
!
!  This subroutine computes the average of a specified quantity on an arbitrary
!  plane in 3d space. The bounding points, bp{1,2,3} are used to define the plane
!  such that the zeta direction 2 -> 1 and the eta direction 2 -> 3.
!
!  When sending the variable to this subroutine, it is important that the
!  ranges (1:nx,1:ny,1:nz) be stated explicitly to avoid incorrect matching
!  of indices between this variable and the x,y,z arrays.
!

use types, only : rprec
use param, only : Nx, Ny, Nz, dx, dy, dz, L_x, L_y
#ifdef PPMPI
use mpi
use param, only : up, down, ierr, MPI_RPREC, status, comm, coord
#endif
use grid_m
use messages
implicit none

real(rprec), intent(in), dimension(:,:,lbz:) :: var
integer, intent(in) :: lbz   !lower bound on z (lbz) for variable sent
real(rprec), intent(in), dimension(:) :: bp1, bp2, bp3

integer, intent(in) :: nzeta, neta

character (*), parameter :: func_name = mod_name // '.plane_avg_3d'

integer :: i, j, nsum

#ifdef PPMPI
integer :: nsum_global
real(rprec) :: var_sum_global
#endif

real(rprec) :: dzeta, deta, vec_mag, zmin, zmax
real(rprec) :: var_sum

real(rprec), dimension(3) :: zeta_vec, eta_vec, eta, cell_center
real(rprec), dimension(3) :: bp4

real(rprec), pointer, dimension(:) :: z

nullify(z)

!  Build computational mesh if needed
if(.not. grid % built) call grid%build()

z => grid % z

nsum = 0
var_sum=0.

!  vector in zeta direction
zeta_vec = bp1 - bp2
!  vector in eta direction
eta_vec   = bp3 - bp2

!  Compute fourth point of plane
bp4 = bp2 + zeta_vec + eta_vec

zmin = min(bp1(3), bp2(3), bp3(3), bp4(3))
zmax = max(bp1(3), bp2(3), bp3(3), bp4(3))

!  Normalize to create unit vector
vec_mag = sqrt(zeta_vec(1)*zeta_vec(1) + zeta_vec(2)*zeta_vec(2) + zeta_vec(3)*zeta_vec(3))
dzeta = vec_mag/nzeta
zeta_vec = zeta_vec / vec_mag

vec_mag = sqrt(eta_vec(1)*eta_vec(1) + eta_vec(2)*eta_vec(2) + eta_vec(3)*eta_vec(3))
deta = vec_mag/neta
eta_vec = eta_vec / vec_mag

!  Compute cell centers
do j=1,neta
  !  Attempt for cache friendliness
  eta = (j - 0.5)*deta*eta_vec
  do i=1,nzeta
  ! Simple vector addition
    cell_center = bp2 + (i - 0.5)*dzeta*zeta_vec + eta

    if(cell_center(3) >= z(1) .and. cell_center(3) < z(nz)) then
      
      !  Include autowrapping for x and y directions
      !cell_center(1) = modulo(cell_center(1), L_x)
      !cell_center(2) = modulo(cell_center(2), L_y)
        
      !  Perform trilinear interpolation       
      var_sum = var_sum + trilinear_interp(var, lbz, cell_center)
      nsum = nsum + 1

    endif

  enddo
enddo

#ifdef PPMPI
!  Perform averaging; all procs have this info
 call mpi_allreduce(var_sum, var_sum_global, 1, MPI_RPREC, MPI_SUM, comm, ierr)
 call mpi_allreduce(nsum, nsum_global, 1, MPI_INTEGER, MPI_SUM, comm, ierr)

 if(nsum_global == 0) then
 
  write(*,'(1a,i4,f9.4)') 'coord, bp1 : ', coord, bp1
  write(*,'(1a,i4,f9.4)') 'coord, bp2 : ', coord, bp2
  write(*,'(1a,i4,f9.4)') 'coord, bp3 : ', coord, bp3
  
  call error(func_name, 'nsum_global = 0')
  
 endif
 
  !  Average over all procs; assuming distribution is even
 plane_avg_3d = var_sum_global / nsum_global
  
  !write(*,*) 'var_sum_global : ', var_sum_global
  
#else
  
  plane_avg_3d = var_sum / nsum
  
#endif
   
nullify(z)

return

end function plane_avg_3d

!**********************************************************************
real(rprec) function points_avg_3d(var, lbz, npoints, points)
!**********************************************************************
!
!  This subroutine computes the arithmetic average of a specified 
!  quantity defined on a set of arbitrary points
!

use types, only : rprec
use param, only : dx, dy, dz, L_x, L_y, nz
#ifdef PPMPI
use mpi
use param, only : up, down, ierr, MPI_RPREC, status, comm, coord
#endif
use grid_m
use messages
implicit none

real(rprec), intent(in), dimension(:,:,lbz:) :: var
integer, intent(in) :: lbz      !lower bound on z (lbz) for variable sent
integer, intent(in) :: npoints
real(rprec), intent(in), dimension(3,npoints) :: points

character (*), parameter :: func_name = mod_name // '.points_avg_3d'

!integer :: istart, jstart, kstart, nsum
integer :: nsum
integer :: n

#ifdef PPMPI
integer :: nsum_global
real(rprec) :: var_sum_global
#endif

real(rprec) :: var_sum
real(rprec) :: xp, yp, zp

real(rprec), pointer, dimension(:) :: z

nullify(z)

!  Check that points is a column major ordered array of dim-3
!if( size(points,1) .ne. 3 ) call error(func_name, 'points not specified correctly.')

!  Build computational mesh if needed
if(.not. grid % built) call grid%build()

z => grid % z

nsum = 0
var_sum=0.

! Get the number of specified points
!npoint = size(points,2)

do n=1, npoints
  
  zp = points(3,n)
  
  if(zp >= z(1) .and. zp < z(nz)) then
  
    xp = points(1,n)
    yp = points(2,n)
    
    !  Include autowrapping for x and y directions
    !xp = modulo(xp, L_x)
    !yp = modulo(yp, L_y)
    
    !  Perform trilinear interpolation
    !istart = autowrap_i( cell_indx('i', dx, xp) )
    !jstart = autowrap_j( cell_indx('j', dy, yp) )
    !kstart = cell_indx('k', dz, zp)
        
    var_sum = var_sum + trilinear_interp(var, lbz, (/ xp, yp, zp /))
    nsum = nsum + 1
    
  endif
  
enddo

#ifdef PPMPI

!  Perform averaging; all procs have this info
call mpi_allreduce(var_sum, var_sum_global, 1, MPI_RPREC, MPI_SUM, comm, ierr)
call mpi_allreduce(nsum, nsum_global, 1, MPI_INTEGER, MPI_SUM, comm, ierr)

if(nsum_global == 0) then
  
  call error(func_name, 'nsum_global = 0')
  
endif
 
!  Average over all procs; assuming distribution is even
points_avg_3d = var_sum_global / nsum_global
  
#else
  
points_avg_3d = var_sum / nsum
  
#endif
   
nullify(z)

return

end function points_avg_3d

!**********************************************************************
integer function buff_indx(i,imax)
!**********************************************************************
!  This function returns the physical index associated with the buffer 
!  region for the specified i and imax. 
!  For i = imax + 1 -> 1 is returned otherwise i is returned
implicit none

integer, intent(in) :: i,imax

if(i == imax + 1) then
  buff_indx = 1
else
  buff_indx = i
endif
  
return
end function buff_indx

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function get_tau_wall() result(twall)       !!jb
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
do jx=1,nx
   do jy=1,ny
      txsum = txsum + txz(jx,jy,1)
      tysum = tysum + tyz(jx,jy,1)
   enddo
enddo
twall = sqrt( (txsum/(nx*ny))**2 + (tysum/(nx*ny))**2  )

return
end function get_tau_wall

end module functions
