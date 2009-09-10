!**********************************************************************
module mpi_defs
!**********************************************************************
use mpi
implicit none

save
public 

integer :: mpierror, mpisize, mpirank, mpicount

contains

!**********************************************************************
subroutine initialize_mpi()
!**********************************************************************
use param, only : nproc
implicit none

!  Initialize mpi communication
call MPI_Init(mpierror)
call MPI_Comm_size(MPI_COMM_WORLD, mpisize, mpierror)
call MPI_Comm_rank(MPI_COMM_WORLD, mpirank, mpierror)

if(mpisize .ne. nproc) then
  write(*,*) 'Error: requested number of processors not equal to specified number!'
  stop
endif

return
end subroutine initialize_mpi

!**********************************************************************
subroutine finalize_mpi()
!**********************************************************************
implicit none
!  Finalize mpi communication
call MPI_FINALIZE(mpierror)
return
end subroutine finalize_mpi

end module mpi_defs

