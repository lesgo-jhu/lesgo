!**********************************************************************
module functions
!**********************************************************************
use messages
use types, only : rprec
implicit none
save
private
public interp_to_uv_grid, trilinear_interp, linear_interp, cell_indx, plane_avg_3D
public buff_indx, points_avg_3D, det2D

character (*), parameter :: mod_name = 'functions'

$if ($MPI)
  !--this dimensioning adds a ghost layer for finite differences
  !--its simpler to have all arrays dimensioned the same, even though
  !  some components do not need ghost layer
  $define $lbz 0
$else
  $define $lbz 1
$endif

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

use types, only : rprec
use param,only : nz
use messages
$if ($MPI)
use param, only : coord, nproc, MPI_RPREC, down, up,  comm, status, ierr
use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWNUP
$endif

implicit none

real(rprec), dimension(:,:,lbz:), intent(IN) :: var
integer, intent(in) :: lbz
real(rprec), allocatable, dimension(:,:,:) :: var_uv

integer :: ubx,uby,ubz
integer :: i,j,k

character (*), parameter :: sub_name = mod_name // '.interp_to_uv_grid'

ubx=ubound(var,1)
uby=ubound(var,2)
ubz=ubound(var,3)

if( ubz .ne. nz ) call error( 'interp_to_uv_grid', 'Input array must lbz:nz z dimensions.')

allocate(var_uv(ubx,uby,lbz:ubz))

do k=1,ubz-1
  do j=1,uby
    do i=1,ubx
      var_uv(i,j,k) = 0.5 * (var(i,j,k+1) + var(i,j,k))
    enddo
  enddo
enddo

$if ($MPI)

!  Take care of top "physical" boundary
if(coord == nproc - 1) var_uv(:,:,ubz) = var_uv(:,:,ubz-1)

!  Sync all overlapping data
if( lbz == 0 ) then
  call mpi_sync_real_array( var_uv, MPI_SYNC_DOWNUP )
elseif( lbz == 1 ) then
  call mpi_sendrecv(var_uv(:,:,1), ubx*uby, MPI_RPREC, down, 1,  &
                    var_uv(:,:,nz), ubx*uby, mpi_rprec, up, 1,   &
                    comm, status, ierr)
endif                    

$else

!  Take care of top "physical" boundary
var_uv(:,:,ubz) = var_uv(:,:,ubz-1)

$endif
  
return 

!deallocate(var_uv)

!$if($MPI)
!deallocate(buf)
!$endif

end function interp_to_uv_grid

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
!  $lbz <= cell_indx < Nz
!
use types, only : rprec
use grid_defs
use messages 
use param, only : nx, ny, nz, L_x, L_y, L_z
implicit none

character (*), intent (in) :: indx
real(rprec), intent(IN) :: dx

real(rprec) :: px ! Global value

character (*), parameter :: func_name = mod_name // '.cell_indx'
real(rprec), parameter :: thresh = 1.e-9_rprec

real(rprec), pointer, dimension(:) :: z

! Nullify pointers
nullify(z)
! Intialize result
cell_indx = -1

if(.not. grid_t % built) call grid_build()

z => grid_t % z

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

   if( cell_indx > Nx .or. cell_indx < 1)  then
     write(*,*) 'px, dx, L_x, cell_indx : ', px, dx, L_x, cell_indx
     call error(func_name, 'Specified point is not in spatial domain (x-direction)')
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

   if( cell_indx > Ny .or. cell_indx < 1)  then
      write(*,*) 'px, dx, L_y, cell_indx : ', px, dx, L_y, cell_indx
      call error(func_name, 'Specified point is not in spatial domain (y-direction)')
   endif

  !  Need to compute local distance to get local k index
  case ('k')

      ! Check upper boundary 
    if( abs( px - z(Nz) ) / L_z < thresh ) then

      cell_indx = Nz-1

    else

      cell_indx = floor ((px - z(1)) / dx) + 1

    endif

    if( cell_indx >= Nz .or. cell_indx < $lbz) call error(func_name, 'Specified point is not in spatial domain (z-direction)')    

  case default

    call error (func_name, 'invalid indx =' // indx)

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
use grid_defs, only : grid_t
use types, only : rprec
use sim_param, only : u,v
use param, only : nx, ny, nz, dx, dy, dz, coord, L_x, L_y
implicit none

integer, intent(IN) :: lbz
real(rprec), dimension(:,:,lbz:), intent(IN) :: var
integer :: istart, jstart, kstart, istart1, jstart1, kstart1
real(rprec), intent(IN), dimension(3) :: xyz

!real(rprec), dimension(2,2,2) :: uvar
integer, parameter :: nvar = 3
integer :: i,j,k
real(rprec) :: u1,u2,u3,u4,u5,u6
real(rprec) :: xdiff, ydiff, zdiff

real(rprec), pointer, dimension(:) :: x,y,z
integer, pointer, dimension(:) :: autowrap_i, autowrap_j

nullify(x,y,z)
nullify(autowrap_i, autowrap_j)

x => grid_t % x
y => grid_t % y
z => grid_t % z
autowrap_i => grid_t % autowrap_i
autowrap_j => grid_t % autowrap_j

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
u1 = linear_interp(var(istart,  jstart,  kstart), &
                   var(istart1, jstart,  kstart), &
                   dx, xdiff)
u2 = linear_interp(var(istart,  jstart1, kstart), &
                   var(istart1, jstart1, kstart), &
                   dx, xdiff)
u3 = linear_interp(var(istart,  jstart,  kstart1), &
                   var(istart1, jstart,  kstart1), &
                   dx, xdiff)
u4 = linear_interp(var(istart,  jstart1, kstart1), &
                   var(istart1, jstart1, kstart1), &
                   dx, xdiff)
!  Perform interpolations in y-direction
u5 = linear_interp(u1,u2,dy,ydiff)
u6 = linear_interp(u3,u4,dy,ydiff)
!  Perform interpolation in z-direction
trilinear_interp = linear_interp(u5,u6,dz,zdiff)

nullify(x,y,z)
nullify(autowrap_i, autowrap_j)

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

real(rprec), intent(IN), dimension(:,:,lbz:) :: var
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

real(rprec), pointer, dimension(:) :: z

nullify(z)

!  Build computational mesh if needed
if(.not. grid_t % built) call grid_build()

z => grid_t % z

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

$if ($MPI)
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
 plane_avg_3D = var_sum_global / nsum_global
  
  !write(*,*) 'var_sum_global : ', var_sum_global
  
 $else
  
  plane_avg_3D = var_sum / nsum
  
 $endif
   
nullify(z)

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

real(rprec), intent(IN), dimension(:,:,lbz:) :: var
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

real(rprec), pointer, dimension(:) :: z

nullify(z)

!  Check that points is a column major ordered array of dim-3
!if( size(points,1) .ne. 3 ) call error(func_name, 'points not specified correctly.')

!  Build computational mesh if needed
if(.not. grid_t % built) call grid_build()

z => grid_t % z

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
   
nullify(z)

return

end function points_avg_3D

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
