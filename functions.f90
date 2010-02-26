!**********************************************************************
module functions
!**********************************************************************
use messages
implicit none
save
private
public trilinear_interp, linear_interp, interp_to_uv_grid

character (*), parameter :: mod_name = 'functions'

contains

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
double precision, dimension(istart:istart+1,jstart:jstart+1,kstart:kstart+1) :: uvar
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
  uvar = u(istart:istart+1,jstart:jstart+1,kstart:kstart+1)
elseif(cvar == 'v') then
  uvar = v(istart:istart+1,jstart:jstart+1,kstart:kstart+1)
elseif(cvar == 'w') then
!  Put w node values on uv grid
  do k=kstart,kstart+1; do j=jstart,jstart+1; do i=istart,istart+1
   uvar(i,j,k) = interp_to_uv_grid('w',i,j,k)
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
u1 = linear_interp(uvar(istart,jstart,kstart),uvar(istart+1,jstart,kstart),dx,xdiff)
u2 = linear_interp(uvar(istart,jstart+1,kstart),uvar(istart+1,jstart+1,kstart),dx,xdiff)
u3 = linear_interp(uvar(istart,jstart,kstart+1),uvar(istart+1,jstart,kstart+1),dx,xdiff)
u4 = linear_interp(uvar(istart,jstart+1,kstart+1),uvar(istart+1,jstart+1,kstart+1),dx,xdiff)
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
subroutine plane_avg_3D(bound_points, nzeta, neta, var, var_avg)
!**********************************************************************
!
!  This subroutine computes the average of a specified quantity on an arbitrary
!  plane in 3D space. 
!

use types, only : rprec
implicit none

real(RPREC), intent(IN), dimension(:,:) :: bound_points
INTEGER, INTENT(IN) :: nzeta, neta
REAL(RPREC), intent(IN), DIMENSION(:,:,:) :: var

REAL(RPREC), intent(OUT) :: var_avg

REAL(RPREC) :: dzeta, deta, Lzeta, Leta, vec_mag

real(RPREC), dimension(3) :: zeta_vec, eta_vec
real(RPREC), allocatable, dimension(3,nzeta,neta) :: cell_centers

!  Attempt for cache friendliness
bp1 = bound_points(:,1)
bp2 = bound_points(:,2) !  Serves as local origin of (zeta,eta) plane
bp3 = bound_points(:,3)

!  vector in zeta direction
zeta_vec = bp1 - bp2
!  vector in eta direction
eta_vec   = bp3 - bp2

dzeta = zeta_vec/nzeta
deta   = eta_vec/neta

!  Normalize to create unit vector
vec_mag = sqrt(zeta_vec(1)*zeta_vec(1) + zeta_vec(2)*zeta_vec(2) + zeta_vec(3)*zeta_vec(3))
zeta_vec = zeta_vec / vec_mag
vec_mag = sqrt(eta_vec(1)*eta_vec(1) + eta_vec(2)*eta_vec(2) + eta_vec(3)*eta_vec(3))
eta_vec = zeta_vec / vec_mag

do j=1,neta
  !  Attempt for cache friendliness
  eta = (j - 0.5)*deta*eta_vec
  do i=1,nzeta
    cell_centers(:,i,j) = bp2 + (i - 0.5)*dzeta*zeta_vec + eta
    zmax = max(cell_centers(3,i,j), zmax)
    zmin = min(cell_centers(3,i,j), zmin)
  enddo
enddo

if(z(nz) <= zmax .or. z(1) >= zmin) then

  if(zmax > z(nz)) then
    isum_recieve = 1
    iavg_send = 1
  endif

  if(zmin < z(1)) then
    isum_send = 1
    iavg_recieve = 1
  endif

  do j=1,neta
    do i=1,nzeta
      if(cell_centers(3,i,j) > z(1) .and. cell_centers(3,i,j) < z(nz)) then
        !  Perform trilinear interpolation
	istart = 
	jstart = 
	kstart = 
	var_sum = 
      endif
      

    enddo
  enddo
    
  



return

end subroutine plane_avg_3D

$endif

end module functions
