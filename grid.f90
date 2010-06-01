
!**********************************************************************
module grid_defs
!**********************************************************************
use types, only : rprec
implicit none
save
private
public x, y, z, zw, grid_build, grid_built

logical :: grid_built

real(rprec), allocatable, dimension(:) :: x, y, z, zw

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
use param, only : nx,ny,nz,dx,dy,dz,coord
implicit none

integer :: i,j,k

grid_built = .false.

$if ($MPI)
  !--this dimensioning adds a ghost layer for finite differences
  !--its simpler to have all arrays dimensioned the same, even though
  !  some components do not need ghost layer
  $define $lbz 0
$else
  $define $lbz 1
$endif

allocate(x(nx+1),y(ny+1),z($lbz:nz),zw($lbz:nz))

do k=$lbz,nz
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
     
grid_built = .true. 

return
end subroutine grid_build 

end module grid_defs
