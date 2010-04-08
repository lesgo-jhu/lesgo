!**********************************************************************
module functions
!**********************************************************************
use messages
implicit none
save
private
public trilinear_interp, linear_interp, index_start, plane_avg_3D
	   !interp_to_uv_grid, interp_to_uv_grid_buf, interp_to_w_grid

character (*), parameter :: mod_name = 'functions'

contains

!**********************************************************************
integer function index_start(indx,dx,px)
!**********************************************************************
! This routine does ...
!
use types, only : rprec
use grid_defs, only : z, grid_built, grid_build
implicit none

character (*), intent (in) :: indx
real(rprec), intent(IN) :: dx
real(rprec), intent(IN) :: px ! Global value

character (*), parameter :: func_name = mod_name // '.index_start'

if(.not. grid_built) call grid_build()

select case (indx)
  case ('i'); index_start = floor (px / dx) + 1
  case ('j'); index_start = floor (px / dx) + 1
  !  Need to compute local distance to get local k index
  case ('k'); index_start = floor ((px - z(1)) / dx) + 1
  case default; call error (func_name, 'invalid indx =' // indx)
end select

return
end function index_start

!**********************************************************************
real(rprec) function trilinear_interp(var,istart,jstart,kstart,xyz)
!**********************************************************************
!
!  This subroutine perform trilinear interpolation for a point that
!  exists in the cell with lower dimension: istart,jstart,kstart
!  for the point xyz
!  
!  istart, jstart, kstart are used to determine the cell location on the
!  uv-grid; these are defined in stats_init
!
!  Takes care of putting w-grid variables onto the uv-grid; this assumes
!  that var is on the uv-grid
!
use grid_defs, only : x,y,z
use types, only : rprec
use sim_param, only : u,v
use param, only : nz, dx, dy, dz, coord
implicit none

real(rprec), dimension(:,:,:), intent(IN) :: var
integer, intent(IN) :: istart, jstart, kstart
real(rprec), intent(IN), dimension(3) :: xyz

real(rprec), dimension(2,2,2) :: uvar
integer, parameter :: nvar = 3
integer :: i,j,k
real(rprec) :: u1,u2,u3,u4,u5,u6
real(rprec) :: xdiff, ydiff, zdiff

!  Initialize stuff
u1=0.
u2=0.
u3=0.
u4=0.
u5=0.
u6=0.

!  Contains the 6 points that make of the cube
uvar = 0.

uvar(:,:,:) = var(istart:istart+1,jstart:jstart+1,kstart:kstart+1)

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

double precision, intent(IN) :: u1, u2, dx, xdiff

linear_interp = u1 + (xdiff) * (u2 - u1) / dx

return
end function linear_interp

!!!**********************************************************************
!!real(rprec) function interp_to_uv_grid(cvar,i,j,k)
!!!**********************************************************************
!!!  This function computes any values the read in value u(k) and
!!!  u(k+1) to the uv grid location k
!!use types, only : rprec
!!use param,only : nz,ld
!!use sim_param, only : w, dudz

!!character(*), intent(IN) :: cvar
!!integer,intent(IN) :: i,j,k


!!if(trim(adjustl(cvar)) == 'w') then
!!  if(k==nz) then
!!    interp_to_uv_grid = 3./2.*w(i,j,k) - 0.5*w(i,j,k-1)
!!  else
!!    interp_to_uv_grid = 0.5*(w(i,j,k)+w(i,j,k+1))
!!  endif
!!elseif(trim(adjustl(cvar)) == 'dudz') then 
!!  if(k==nz) then
!!    interp_to_uv_grid = 3./2.*dudz(i,j,k) - 0.5*dudz(i,j,k-1)
!!  else
!!    interp_to_uv_grid = 0.5*(dudz(i,j,k)+dudz(i,j,k+1))
!!  endif
!!else
!!  write(*,*) 'Error: variable specification not specified properly!'
!!  stop
!!endif
!!return
!!end function interp_to_uv_grid

!!!**********************************************************************
!!real(rprec) function interp_to_uv_grid_buf(cvar,i,j,k)
!!!**********************************************************************
!!!  This function computes any values the read in value u(k) and
!!!  u(k+1) to the uv grid location k
!!use types, only : rprec
!!use param,only : nz,ld
!!use sim_param, only : w, dudz
!!$if ($MPI)
!!use mpi
!!use param, only : MPI_RPREC, down, up, comm, status, ierr, nproc, &
!!  coord, coord_of_rank, rank
!!$endif

!!implicit none

!!character(*), intent(IN) :: cvar
!!integer,intent(IN) :: i,j,k
!!real(rprec), dimension(2) :: var

!!$if ($MPI)
!!real(rprec) :: buf
!!$endif

!!if(trim(adjustl(cvar)) == 'w') then
!!  
!!  $if ($MPI)
!!  if(coord_of_rank(rank) < nproc - 1) then
!!    call mpi_recv (buf, 1, MPI_RPREC, up, 100, comm, status, ierr)
!!  endif
!!  if(coord_of_rank(rank) > 0) then
!!    call mpi_send (w(1, 1, 2), 1, MPI_RPREC, down, 100, comm, ierr)
!!  endif
!!  $endif

!!  if(k==nz) then
!!  
!!    $if ($MPI)
!!    if (coord == nproc - 1) then
!!      var = (/ w(i,j,k), w(i,j,k) /)
!!    else
!!      var = (/buf, w(i,j,k) /)
!!	endif
!!    $else

!!    var = (/ w(i,j,k), w(i,j,k) /)    

!!    $endif
!!  else
!!    var = (/ w(i,j,k+1), w(i,j,k)/)
!!  endif

!!elseif(trim(adjustl(cvar)) == 'dudz') then 

!!  $if ($MPI)
!!  
!!  if(coord_of_rank(rank) < nproc - 1) then
!!    call mpi_recv (buf, 1, MPI_RPREC, up, 200, comm, status, ierr)
!!  endif
!!  
!!  if(coord_of_rank(rank) > 0) then
!!    call mpi_send (dudz(1, 1, 2), 1, MPI_RPREC, down, 200, comm, ierr)
!!  endif
!!  $endif

!!  if(k==nz) then
!!  
!!    $if ($MPI) 
!!    if (coord == nproc - 1) then
!!      var = (/ dudz(i,j,k), dudz(i,j,k) /)
!!    else
!!      var = (/buf, dudz(i,j,k) /)
!!	endif
!!    $else

!!    var = (/ dudz(i,j,k), dudz(i,j,k) /)    

!!    $endif
!!  else
!!    var = (/ dudz(i,j,k+1), dudz(i,j,k)/)
!!  endif
!!else
!!  write(*,*) 'Error: variable specification not specified properly!'
!!  stop
!!endif

!!interp_to_uv_grid_buf = 0.5*(var(1) + var(2))

!!return
!!end function interp_to_uv_grid_buf


!!!**********************************************************************
!!real(rprec) function interp_to_w_grid(cvar,i,j,k)
!!!**********************************************************************
!!!  This function computes any values the read in value u(k) and
!!!  u(k-1) to the w grid location k
!!use types, only : rprec
!!use param,only : nz,ld,USE_MPI,coord
!!use sim_param, only : u, v
!!$if ($LVLSET)
!!use level_set, only : phi
!!$endif

!!character(*), intent(IN) :: cvar
!!integer,intent(IN) :: i,j,k
!!real(rprec), dimension(2) :: var

!!if(trim(adjustl(cvar)) == 'u') then
!!  if(k == 1) then
!!    if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
!!	  var = (/ 0., 0./)
!!	else
!!	  var = (/u(i,j,k), u(i,j,k-1)/);
!!	endif
!!  else
!!    var = (/u(i,j,k), u(i,j,k-1)/);
!!  endif

!!elseif(trim(adjustl(cvar)) == 'v') then
!!  if(k == 1) then
!!    if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
!!	  var = (/ 0., 0./)
!!	else
!!	  var = (/v(i,j,k), v(i,j,k-1)/);
!!	endif
!!  else
!!    var = (/v(i,j,k), v(i,j,k-1)/);
!!  endif
!!$if ($LVLSET)
!!elseif(trim(adjustl(cvar)) == 'phi') then
!!  var = (/phi(i,j,k), phi(i,j,k-1)/);
!!$endif
!!else
!!  write(*,*) 'Error: variable specification not specified properly!'
!!  stop
!!endif  

!!interp_to_w_grid = 0.5*(var(1) + var(2))

!!return
!!end function interp_to_w_grid

!**********************************************************************
real(rprec) function plane_avg_3D(var, bp1, bp2, bp3, nzeta, neta)
!**********************************************************************
!
!  This subroutine computes the average of a specified quantity on an arbitrary
!  plane in 3D space. 
!

use types, only : rprec
use param, only : Nx, Ny, Nz, dx, dy, dz, L_x, L_y
$if ($MPI)
use mpi
use param, only : up, down, ierr, MPI_RPREC, status, comm, coord
$endif
use grid_defs

implicit none

real(rprec), intent(IN), dimension(:,:,:) :: var
real(RPREC), intent(IN), dimension(:) :: bp1, bp2, bp3

INTEGER, INTENT(IN) :: nzeta, neta

character (*), parameter :: func_name = mod_name // '.plane_avg_3D'

integer :: i, j, istart, jstart, kstart, nsum

$if ($MPI)
integer :: isum_send, isum_recieve, iavg_send, iavg_recieve
integer :: nsum_global
REAL(RPREC) :: var_sum_global
$endif

REAL(RPREC) :: dzeta, deta, Lzeta, Leta, vec_mag, zmin, zmax
REAL(RPREC) :: xdiff, ydiff, zdiff, var_sum, var_interp

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

!if(coord == 0) then
!  write(*,'(1a,3f12.6)') 'zeta_vec : ', zeta_vec
!  write(*,'(1a,3f12.6)') 'eta_vec  : ', eta_vec
!endif

!!  Check if plane is associated with processor
!if(z(nz) <= zmax .or. z(1) >= zmin) then

!  Compute cell centers
  do j=1,neta
  !  Attempt for cache friendliness
    eta = (j - 0.5)*deta*eta_vec
    do i=1,nzeta
    ! Simple vector addition
      cell_center = bp2 + (i - 0.5)*dzeta*zeta_vec + eta

      if(cell_center(3) >= z(1) .and. cell_center(3) < z(nz)) then
        !  Perform trilinear interpolation
        !  Include autowrapping for x and y directions
        istart = index_start('i', dx, modulo(cell_center(1),L_x))
        jstart = index_start('j', dy, modulo(cell_center(2),L_y))
        kstart = index_start('k', dz, cell_center(3))
        
        var_sum = var_sum + trilinear_interp(var, istart, jstart, kstart, cell_center)
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

!  Perform averaging; all procs have this info
 call mpi_allreduce(var_sum, var_sum_global, 1, MPI_RPREC, MPI_SUM, comm, ierr)
 call mpi_allreduce(nsum, nsum_global, 1, MPI_INTEGER, MPI_SUM, comm, ierr)

  !  Average over all procs; assuming distribution is even
  plane_avg_3D = var_sum_global / nsum_global
  
  !write(*,*) 'var_sum_global : ', var_sum_global
  
 $else
  
  plane_avg_3D = var_sum / nsum
  
 $endif
   
!else
!  write(*,*) 'need to put message here'
!  stop
!endif


return

end function plane_avg_3D

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

end module functions
