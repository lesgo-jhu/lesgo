!**********************************************************************
module functions
!**********************************************************************
use messages
use types, only : rprec
implicit none
save
private
public trilinear_interp, linear_interp, cell_indx, plane_avg_3D
public buff_indx, points_avg_3D, det2D

character (*), parameter :: mod_name = 'functions'

contains

!**********************************************************************
integer function cell_indx(indx,dx,px)
!**********************************************************************
! This routine takes index=['i' or 'j' or 'k'] and the magnitude of the 
!   spacing=[dx or dy or dz] and the [x or y or z] location and returns
!   the value of the lower index (cell index)
!
use types, only : rprec
use grid_defs, only : z, grid_built, grid_build
use messages 
use param, only : nx, ny, nz
implicit none

character (*), intent (in) :: indx
real(rprec), intent(IN) :: dx
real(rprec), intent(IN) :: px ! Global value

character (*), parameter :: func_name = mod_name // '.cell_indx'

if(.not. grid_built) call grid_build()

select case (indx)
  case ('i')
    cell_indx = floor (px / dx) + 1
    if( cell_indx > Nx .or. cell_indx < 1) call error(func_name, 'Specified point is not in spatial domain - wrap with modulo')
  case ('j')
    cell_indx = floor (px / dx) + 1
    if( cell_indx > Ny .or. cell_indx < 1) call error(func_name, 'Specified point is not in spatial domain - wrap with modulo')
  !  Need to compute local distance to get local k index
  case ('k')
    cell_indx = floor ((px - z(1)) / dx) + 1
    if( cell_indx > Nz .or. cell_indx < 0) call error(func_name, 'Specified point is not in spatial domain')    
  case default
    call error (func_name, 'invalid indx =' // indx)
end select

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
!  uv-grid; these are defined in stats_init
!
!  Takes care of putting w-grid variables onto the uv-grid; this assumes
!  that var is on the uv-grid
!
!  The variable sent to this subroutine should have a lower-bound-on-z 
!  (lbz) set as an input so the k-index will match the k-index of z.  
!  Before calling this function, make sure the point exists on the coord
!  [ test using: z(1) \leq z_p < z(nz-1) ]
!
use grid_defs, only : x,y,z, autowrap_i, autowrap_j
use types, only : rprec
use sim_param, only : u,v
use param, only : nz, dx, dy, dz, coord, L_x, L_y
implicit none

real(rprec), dimension(:,:,:), intent(IN) :: var
integer, intent(IN) :: lbz
integer :: istart, jstart, kstart
real(rprec), intent(IN), dimension(3) :: xyz

real(rprec), dimension(2,2,2) :: uvar
integer, parameter :: nvar = 3
integer :: i,j,k
real(rprec) :: u1,u2,u3,u4,u5,u6
real(rprec) :: xdiff, ydiff, zdiff

!  Initialize stuff
u1=0.; u2=0.; u3=0.; u4=0.; u5=0.; u6=0.

! Wrap x,y if necessary (periodic BCs)
! Determine istart, jstart, kstart by calling cell_indx
istart = cell_indx('i',dx,mod(xyz(1),L_x))
jstart = cell_indx('j',dy,mod(xyz(2),L_y))
kstart = cell_indx('k',dz,xyz(3))
!write(*,*) 'coord,is,js,ks',coord,istart,jstart,kstart

!  Contains the 6 points that make of the cube
uvar = 0.

! Extra term with kstart accounts for shift in var k-index is lbz.ne.1
do k=0,1
  do j=0,1
    do i=0,1
      uvar(i+1,j+1,k+1) = var(istart+i, jstart+j, kstart+(1-lbz)+k)
    enddo
  enddo
enddo

!  Compute xdiff
xdiff = xyz(1) - x(istart)
!  Compute ydiff
ydiff = xyz(2) - y(jstart)
!  Compute zdiff
zdiff = xyz(3) - z(kstart)

!  Perform the 7 linear interpolations
!  Perform interpolations in x-direction 
u1 = linear_interp(uvar(1,1,1),uvar(2,1,1),dx,xdiff)
u2 = linear_interp(uvar(1,2,1),uvar(2,2,1),dx,xdiff)
u3 = linear_interp(uvar(1,1,2),uvar(2,1,2),dx,xdiff)
u4 = linear_interp(uvar(1,2,2),uvar(2,2,2),dx,xdiff)
!  Perform interpolations in y-direction
u5 = linear_interp(u1,u2,dy,ydiff)
u6 = linear_interp(u3,u4,dy,ydiff)
!  Perform interpolation in z-direction
trilinear_interp = linear_interp(u5,u6,dz,zdiff)

return
end function trilinear_interp

!**********************************************************************
real(rprec) function linear_interp(u1,u2,dx,xdiff)
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

real(rprec), intent(IN) :: u1, u2, dx, xdiff

linear_interp = u1 + (xdiff) * (u2 - u1) / dx

return
end function linear_interp

!**********************************************************************
real(rprec) function plane_avg_3D(var, lbz, bp1, bp2, bp3, nzeta, neta)
!**********************************************************************
!
!  This subroutine computes the average of a specified quantity on an arbitrary
!  plane in 3D space. The bounding points, bp{1,2,3} are used to define the plane
!  such that the zeta direction 2 -> 1 and the eta direction 2 -> 3.
!
!  When sending the variable to this subroutine, it is important that the
!  ranges (1:nx,1:ny,1:nz) be stated explicitly to avoid incorrect matching
!  of indices between this variable and the x,y,z arrays.
!

use types, only : rprec
use param, only : Nx, Ny, Nz, dx, dy, dz, L_x, L_y
$if ($MPI)
use mpi
use param, only : up, down, ierr, MPI_RPREC, status, comm, coord
$endif
use grid_defs
use messages
implicit none

real(rprec), intent(IN), dimension(:,:,:) :: var
integer, intent(IN) :: lbz   !lower bound on z ($lbz) for variable sent
real(RPREC), intent(IN), dimension(:) :: bp1, bp2, bp3

INTEGER, INTENT(IN) :: nzeta, neta

character (*), parameter :: func_name = mod_name // '.plane_avg_3D'

integer :: i, j, nsum

$if ($MPI)
integer :: nsum_global
REAL(RPREC) :: var_sum_global
$endif

REAL(RPREC) :: dzeta, deta, vec_mag, zmin, zmax
REAL(RPREC) :: var_sum

real(RPREC), dimension(3) :: zeta_vec, eta_vec, eta, cell_center
real(RPREC), dimension(3) :: bp4

!  Build computational mesh if needed
if(.not. grid_built) call grid_build()

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
      cell_center(1) = modulo(cell_center(1), L_x)
      cell_center(2) = modulo(cell_center(2), L_y)
        
      !  Perform trilinear interpolation       
      var_sum = var_sum + trilinear_interp(var, lbz, cell_center)
      nsum = nsum + 1

    endif

  enddo
enddo

$if ($MPI)
!  Perform averaging; all procs have this info
 call mpi_allreduce(var_sum, var_sum_global, 1, MPI_RPREC, MPI_SUM, comm, ierr)
 call mpi_allreduce(nsum, nsum_global, 1, MPI_INTEGER, MPI_SUM, comm, ierr)

 if(nsum_global == 0) then
 
  write(*,'(1a,1i,3f9.4)') 'coord, bp1 : ', coord, bp1
  write(*,'(1a,1i,3f9.4)') 'coord, bp2 : ', coord, bp2
  write(*,'(1a,1i,3f9.4)') 'coord, bp3 : ', coord, bp3
  
  call error(func_name, 'nsum_global = 0')
  
 endif
 
  !  Average over all procs; assuming distribution is even
 plane_avg_3D = var_sum_global / nsum_global
  
  !write(*,*) 'var_sum_global : ', var_sum_global
  
 $else
  
  plane_avg_3D = var_sum / nsum
  
 $endif
   
return

end function plane_avg_3D

!**********************************************************************
real(rprec) function points_avg_3D(var, lbz, npoints, points)
!**********************************************************************
!
!  This subroutine computes the arithmetic average of a specified 
!  quantity defined on a set of arbitrary points
!

use types, only : rprec
use param, only : dx, dy, dz, L_x, L_y, nz
$if ($MPI)
use mpi
use param, only : up, down, ierr, MPI_RPREC, status, comm, coord
$endif
use grid_defs
use messages
implicit none

real(rprec), intent(IN), dimension(:,:,:) :: var
integer, intent(IN) :: lbz      !lower bound on z ($lbz) for variable sent
integer, intent(IN) :: npoints
real(rprec), intent(IN), dimension(3,npoints) :: points

character (*), parameter :: func_name = mod_name // '.points_avg_3D'

integer :: istart, jstart, kstart, nsum
integer :: n

$if ($MPI)
integer :: nsum_global
real(rprec) :: var_sum_global
$endif

real(rprec) :: var_sum
real(rprec) :: xp, yp, zp

!  Check that points is a column major ordered array of dim-3
!if( size(points,1) .ne. 3 ) call error(func_name, 'points not specified correctly.')

!  Build computational mesh if needed
if(.not. grid_built) call grid_build()

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
    xp = modulo(xp, L_x)
    yp = modulo(yp, L_y)
    
    !  Perform trilinear interpolation
    istart = autowrap_i( cell_indx('i', dx, xp) )
    jstart = autowrap_j( cell_indx('j', dy, yp) )
    kstart = cell_indx('k', dz, zp)
        
    var_sum = var_sum + trilinear_interp(var, lbz, (/ xp, yp, zp /))
    nsum = nsum + 1
    
  endif
  
enddo

$if ($MPI)

!  Perform averaging; all procs have this info
call mpi_allreduce(var_sum, var_sum_global, 1, MPI_RPREC, MPI_SUM, comm, ierr)
call mpi_allreduce(nsum, nsum_global, 1, MPI_INTEGER, MPI_SUM, comm, ierr)

if(nsum_global == 0) then
  
  call error(func_name, 'nsum_global = 0')
  
endif
 
!  Average over all procs; assuming distribution is even
points_avg_3D = var_sum_global / nsum_global
  
$else
  
points_avg_3D = var_sum / nsum
  
$endif
   
return

end function points_avg_3D

!**********************************************************************
real(rprec) function fractal_scale(var, scale_fact, ng)
!**********************************************************************
use types, only : rprec
implicit none

real(rprec), intent(in) :: var, scale_fact
integer, intent(in) :: ng

fractal_scale = var * scale_fact ** ( ng - 1 )

return
end function fractal_scale

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

!**********************************************************************
real(rprec) function det2D(A)
!**********************************************************************
!  Computes the determinant of the 2D matrix A
!
implicit none
real(rprec), intent(IN), dimension(2,2) :: A

det2D = A(1,1) * A(2,2) - A(1,2) * A(2,1)

return
end function det2D

end module functions
