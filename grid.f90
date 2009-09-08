!**********************************************************************
module grid_defs
!**********************************************************************
use types, only : rprec
implicit none
save
private
public x, y, z, zw, grid_build
real(rprec), allocatable, dimension(:) :: x,y,z,zw

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

allocate(x(nx),y(ny),z(nz),zw(nz))

do k=1,nz
  $if ($MPI)
  z(k) = (coord*(nz-1) + k - 0.5_rprec) * dz
  $else
  z(k) = (k - 0.5_rprec) * dz
  $endif
  do j=1,ny
    y(j) = (j-1)*dy
    do i=1,nx
      x(i) = (i - 1)*dx
    enddo
  enddo
enddo
zw = z - dz/2._rprec
     
return
end subroutine grid_build 

end module grid_defs