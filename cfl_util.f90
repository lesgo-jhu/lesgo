!**********************************************************************
subroutine maxcfl(cfl)
!**********************************************************************
use param, only : dt, dx, dy, dz
use sim_param, only : u,v,w

$if($MPI)
use types, only : rprec
use mpi
use param, only : up, down, ierr, MPI_RPREC, status, comm, coord
$endif

real(rprec), intent(OUT) :: cfl

$if($MPI)
real(rprec) :: cfl_buf
$endif

cfl = maxval( (/ dt * u / dx, dt * v / dy, dt * w / dz /) )

$if($MPI)
call mpi_allreduce(cfl, cfl_buf, 1, MPI_RPREC, MPI_MAX, comm, ierr)
cfl = cfl_buf
$endif

return
end subroutine maxcfl
