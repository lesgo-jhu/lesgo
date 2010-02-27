!**********************************************************************
module functions
!**********************************************************************
use messages
implicit none
save
private
public trilinear_interp, linear_interp, interp_to_uv_grid, &
       find_istart, plane_avg_3D

character (*), parameter :: mod_name = 'functions'

contains

!**********************************************************************
subroutine find_istart(x,nx,px,istart,xdiff)
!**********************************************************************
! This routine should be setup to directly compute istart, xdiff from
! modulo function
!
implicit none

integer, intent(IN) :: nx
double precision, dimension(nx), intent(IN) :: x
double precision, intent(IN) :: px
integer, intent(OUT) :: istart
double precision, intent(OUT) :: xdiff

integer :: i

isearch: do i=2,nx
  if(x(i) >= px) then
    istart = i-1
    xdiff = px - x(istart)
    exit isearch
  endif
enddo isearch

return
end subroutine find_istart

!**********************************************************************
double precision function trilinear_interp(cvar,istart,jstart,kstart,xyz)
!**********************************************************************
!
!  This subroutine perform trilinear interpolation for a point that
!  exists in the cell with lower dimension: istart,jstart,kstart
!  for the point xyz
!  
!  istart, jstart, kstart are used to determine the cell location;
!  these are defined in stats_init
!
!  Takes care of putting w on uv grid
!
use grid_defs, only : x,y,z
use types, only : rprec
use sim_param, only : u,v,w
use param, only : nz, dx, dy, dz, coord
implicit none

character(*), intent(IN) :: cvar
double precision, dimension(2,2,2) :: uvar
integer, intent(IN) :: istart, jstart, kstart
double precision, intent(IN), dimension(3) :: xyz
integer, parameter :: nvar = 3
integer :: i,j,k
double precision :: u1,u2,u3,u4,u5,u6
double precision :: xdiff, ydiff, zdiff

!  Initialize stuff
u1=0.
u2=0.
u3=0.
u4=0.
u5=0.
u6=0.

!  Contains the 6 points that make of the cube
uvar = 0.

if(cvar == 'u') then
  uvar(:,:,:) = u(istart:istart+1,jstart:jstart+1,kstart:kstart+1)
elseif(cvar == 'v') then
  uvar(:,:,:) = v(istart:istart+1,jstart:jstart+1,kstart:kstart+1)
elseif(cvar == 'w') then
!  Put w node values on uv grid
  do k=1,2; do j=1,2; do i=1,2
   uvar(i,j,k) = interp_to_uv_grid('w', i + istart - 1, j + jstart - 1, k + kstart - 1)
  enddo; enddo; enddo
else
  write(*,*) 'Error: variable specification not specified properly!'
  stop
endif

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
double precision function linear_interp(u1,u2,dx,xdiff)
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

implicit none

double precision, intent(IN) :: u1, u2, dx, xdiff

linear_interp = u1 + (xdiff) * (u2 - u1) / dx

return
end function linear_interp

!**********************************************************************
double precision function interp_to_uv_grid(cvar,i,j,k)
!**********************************************************************
!  This function computes any values the read in value u(k) and
!  u(k+1) to the uv grid location k
use param,only : nz,ld
use sim_param, only : w, dudz

character(*), intent(IN) :: cvar
integer,intent(IN) :: i,j,k


if(trim(adjustl(cvar)) == 'w') then
  if(k==nz) then
    interp_to_uv_grid = 3./2.*w(i,j,k) - 0.5*w(i,j,k-1)
  else
    interp_to_uv_grid = 0.5*(w(i,j,k)+w(i,j,k+1))
  endif
elseif(trim(adjustl(cvar)) == 'dudz') then 
  if(k==nz) then
    interp_to_uv_grid = 3./2.*dudz(i,j,k) - 0.5*dudz(i,j,k-1)
  else
    interp_to_uv_grid = 0.5*(dudz(i,j,k)+dudz(i,j,k+1))
  endif
else
  write(*,*) 'Error: variable specification not specified properly!'
  stop
endif
return
end function interp_to_uv_grid

$if ($DEVEL)
!**********************************************************************
double precision function plane_avg_3D(cvar, bound_points, nzeta, neta)
!**********************************************************************
!
!  This subroutine computes the average of a specified quantity on an arbitrary
!  plane in 3D space. 
!

use types, only : rprec
use param, only : Nx, Ny, Nz
$if ($MPI)
use mpi
use param, only : up, down, ierr, MPI_RPREC, status, rank_of_coord, coord, comm
$endif
use grid_defs

implicit none

character(*), intent(IN) :: cvar

real(RPREC), intent(IN), dimension(:,:) :: bound_points
INTEGER, INTENT(IN) :: nzeta, neta

integer :: i, j, istart, jstart, kstart, nsum

$if ($MPI)
integer :: isum_send, isum_recieve, iavg_send, iavg_recieve
integer :: nsum_global
REAL(RPREC) :: var_sum_global
$endif

REAL(RPREC) :: dzeta, deta, Lzeta, Leta, vec_mag, zmin, zmax
REAL(RPREC) :: xdiff, ydiff, zdiff, var_sum, var_interp

real(RPREC), dimension(3) :: zeta_vec, eta_vec, eta, cell_center
real(RPREC), dimension(3) :: bp1, bp2, bp3, bp4

real(RPREC), dimension(3,nzeta,neta) :: cell_centers ! nodes of plane

!  Build computational mesh if needed
if(.not. grid_built) call grid_build()

$if ($MPI)
isum_send    = .false.
isum_recieve = .false.
iavg_send    = .false.
iavg_recieve = .false.
$endif

nsum = 0
var_sum=0.

!  Attempt for cache friendliness
bp1 = bound_points(:,1)
bp2 = bound_points(:,2) !  Serves as local origin of (zeta,eta) plane
bp3 = bound_points(:,3)
if(coord == 0) then
  write(*,'(1a,3f12.6)') 'bp1 : ', bp1
  write(*,'(1a,3f12.6)') 'bp2 : ', bp2
  write(*,'(1a,3f12.6)') 'bp3 : ', bp3
  !stop
endif

!  vector in zeta direction
zeta_vec = bp1 - bp2
!  vector in eta direction
eta_vec   = bp3 - bp2

if(coord == 0) then
  write(*,'(1a,3f12.6)') 'zeta_vec : ', zeta_vec
  write(*,'(1a,3f12.6)') 'eta_vec : ', eta_vec
  !stop  stop
endif

!  Compute fourth point of plane
bp4 = bp2 + zeta_vec + eta_vec

zmin = min(bp1(3), bp2(3), bp3(3), bp4(3))
zmax = max(bp1(3), bp2(3), bp3(3), bp4(3))

!  Normalize to create unit vector
vec_mag = sqrt(zeta_vec(1)*zeta_vec(1) + zeta_vec(2)*zeta_vec(2) + zeta_vec(3)*zeta_vec(3))
dzeta = vec_mag/(nzeta-1)
zeta_vec = zeta_vec / vec_mag

vec_mag = sqrt(eta_vec(1)*eta_vec(1) + eta_vec(2)*eta_vec(2) + eta_vec(3)*eta_vec(3))
deta = vec_mag/(neta-1)
eta_vec = eta_vec / vec_mag

if(coord == 0) then
  write(*,'(1a,3f12.6)') 'zeta_vec : ', zeta_vec
  write(*,'(1a,3f12.6)') 'eta_vec : ', eta_vec
  !stop
endif

!  Compute cell centers
do j=1,neta
  !  Attempt for cache friendliness
  eta = (j - 1)*deta*eta_vec
  do i=1,nzeta
    ! Simple vector addition
    cell_centers(:,i,j) = bp2 + (i - 1)*dzeta*zeta_vec + eta
	
!	if(coord == 0) write(*,'(1a,3i,3f12.6)') 'i, j, x, y, z : ', i, j, cell_centers(1,i,j), cell_centers(2,i,j), cell_centers(3,i,j)
  enddo
enddo

!  Check if plane is associated with processor
if(z(nz) <= zmax .or. z(1) >= zmin) then

!  $if ($MPI)
!  !  Check if points lie below current proc
!  if(zmin < z(1)) then
!    isum_send = .true.
!    iavg_recieve = .true.
!  endif

  !  Check if points lie above current proc
!  if(zmax > z(nz)) then
!    isum_recieve = .true.
!    iavg_send = .true.
!  endif
  
!  $endif

  do j=1,neta
    do i=1,nzeta

      if(cell_centers(3,i,j) > z(1) .and. cell_centers(3,i,j) < z(nz)) then
        !  Perform trilinear interpolation
        call find_istart(x,Nx,cell_centers(1,i,j),istart,xdiff)
        call find_istart(y,Ny,cell_centers(2,i,j),jstart,ydiff)
        call find_istart(z,Nz,cell_centers(3,i,j),kstart,zdiff)

        var_interp = trilinear_interp(cvar,istart,jstart,kstart,cell_centers(:,i,j))
        !if(cvar == 'u') write(*,*) 'var_interp : ', cvar, var_interp
        var_sum = var_sum + var_interp
        nsum = nsum + 1
 
      endif
    enddo
  enddo

  $if ($MPI)
!  if(isum_recieve) then
!    CALL MPI_Recv(var_sum_up, 1, MPI_RPREC, up, 1, MPI_COMM_WORLD, status, ierr)
!    CALL MPI_Recv(nsum_up, 1, MPI_RPREC, up, 2, MPI_COMM_WORLD, status, ierr)
	
!     var_sum = var_sum + var_sum_up
!     nsum = nsum + nsum_up
	
!  endif
  
!  if(isum_send) then
!    CALL MPI_Send(var_sum, 1, MPI_RPREC, down, 1, MPI_COMM_WORLD, ierr)
!    CALL MPI_Send(nsum, 1, MPI_RPREC, down, 2, MPI_COMM_WORLD, ierr)
	
!  else
!    ! Should be the bottom most proc 
!	plane_avg_3D = var_sum / nsum
!     call mpi_scatter(plane_avg_3D, 1, MPI_RPREC, plane_avg_3D, 1, MPI_RPREC, rank_of_coord(coord), MPI_COMM_WORLD, ierr)
!  endif
  
!  if(iavg_recieve) CALL MPI_Recv(plane_avg_3D, 1, MPI_RPREC, down, 3, MPI_COMM_WORLD, status, ierr)
!  if(iavg_send) CALL MPI_Send(plane_avg_3D, 1, MPI_RPREC, up, 3, MPI_COMM_WORLD, ierr)
  
  
 $else
  
  plane_avg_3D = var_sum / nsum
  
 $endif
   
else
  write(*,*) 'need to put message here'  
endif

$if ($MPI)
!  Perform averaging and scatter to all procs
 call mpi_allreduce(var_sum, var_sum_global, 1, MPI_RPREC, MPI_SUM, comm, ierr)
 call mpi_allreduce(nsum, nsum_global, 1, MPI_INT, MPI_SUM, comm, ierr)

 !write(*,*) 'coord, var_sum : ', coord, var_sum
  !if(coord==0) write(*,*) 'coord, nsum_global : ', coord, nsum_global
  !  Average over all procs; assuming distribution is even
  plane_avg_3D = var_sum_global / nsum_global
$endif

return

end function plane_avg_3D

$endif

end module functions
