!**********************************************************************
subroutine cfl_max(cfl)
!**********************************************************************
use types, only : rprec
use param, only : dt, dx, dy, dz, nx, ny, nz
use sim_param, only : u,v,w

$if($MPI)
use mpi
use param, only : up, down, ierr, MPI_RPREC, status, comm, coord
$endif

implicit none
real(rprec), intent(OUT) :: cfl
$if($PGI)
integer :: i, j, k, icount
real(rprec), allocatable, dimension(:) :: u_inter
$endif

$if($MPI)
real(rprec) :: cfl_buf
$endif

$if($PGI)
allocate(u_inter(nx*ny*(nz-1)*3))
icount=0
do k=1,nz-1
  do j=1,ny
    do i=1,nx
      icount = icount + 1
      u_inter(icount) = abs(u(i,j,k)) / dx
    enddo
  enddo
enddo
do k=1,nz-1
  do j=1,ny
    do i=1,nx
      icount = icount + 1
      u_inter(icount) = abs(v(i,j,k)) / dy
    enddo
  enddo
enddo
do k=1,nz-1
  do j=1,ny
    do i=1,nx
      icount = icount + 1
      u_inter(icount) = abs(w(i,j,k)) / dz
    enddo
  enddo
enddo

$endif

$if($PGI)
cfl = maxval( u_inter )
deallocate( u_inter )
$else
!cfl = maxval( (/ dabs(u(1:nx,1:ny,1:nz-1)) / dx, &
!                 dabs(v(1:nx,1:ny,1:nz-1)) / dy, &
!                 dabs(w(1:nx,1:ny,1:nz-1)) / dz /) )
cfl = maxval( abs(u(1:nx,1:ny,1:nz-1)) / dx )
cfl = maxval( (/ cfl, abs(v(1:nx,1:ny,1:nz-1)) / dy /) )
cfl = maxval( (/ cfl, abs(w(1:nx,1:ny,1:nz-1)) / dz /) )
$endif

cfl = dt * cfl

$if($MPI)
call mpi_allreduce(cfl, cfl_buf, 1, MPI_RPREC, MPI_MAX, comm, ierr)
cfl = cfl_buf
$endif

return
end subroutine cfl_max

$if($CFL_DT)
!**********************************************************************
subroutine cfl_set_dt(dt)
!**********************************************************************
use types, only : rprec
use param, only : cfl, dx, dy, dz, nx, ny, nz
use sim_param, only : u,v,w

$if($MPI)
use mpi
use param, only : up, down, ierr, MPI_RPREC, status, comm, coord
$endif

implicit none

real(rprec), intent(OUT) :: dt

$if($PGI)
integer :: i, j, k, icount
real(rprec), allocatable, dimension(:) :: u_inter
$endif

$if($MPI)
real(rprec) :: dt_buf
$endif

$if($PGI)
allocate(u_inter(nx*ny*(nz-1)*3))
icount=0
do k=1,nz-1
  do j=1,ny
    do i=1,nx
      icount = icount + 1
      u_inter(icount) = dx / abs(u(i,j,k))
    enddo
  enddo
enddo
do k=1,nz-1
  do j=1,ny
    do i=1,nx
      icount = icount + 1
      u_inter(icount) = dy / abs(v(i,j,k))
    enddo
  enddo
enddo
do k=1,nz-1
  do j=1,ny
    do i=1,nx
      icount = icount + 1
      u_inter(icount) = dz / abs(w(i,j,k))
    enddo
  enddo
enddo
$endif

$if($PGI)
dt = minval( u_inter )
deallocate(u_inter)
$else
!dt = minval( (/ dx / abs(u(1:nx,1:ny,1:nz-1)), &
!                dy / abs(v(1:nx,1:ny,1:nz-1)), &
!                dz / abs(w(1:nx,1:ny,1:nz-1)) /) )
dt = minval( dx / abs(u(1:nx,1:ny,1:nz-1)) )
dt = minval((/ dt, dy / abs(v(1:nx,1:ny,1:nz-1)) /) )
dt = minval((/ dt, dz / abs(w(1:nx,1:ny,1:nz-1)) /) )
$endif

dt = cfl * dt

$if($MPI)
call mpi_allreduce(dt, dt_buf, 1, MPI_RPREC, MPI_MIN, comm, ierr)
dt = dt_buf
$endif

return
end subroutine cfl_set_dt
$endif
