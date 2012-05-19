!*******************************************************************************
subroutine finalize()
!*******************************************************************************
! 
! This subroutine is called by the main program. It is a driver subroutine for
! calling all the finalize routines of the various lesgo modules.
!
use param, only : coord, sgs_hist_calc
use io, only : closefiles
use sgs_hist
$if($MPI)
use param, only : MPI_COMM_WORLD, ierr
$endif
$if ($LVLSET and $RNS_LS )
use rns_ls, only : rns_finalize_ls
$endif
$if ($TURBINES)
use turbines, only : turbines_finalize
$endif

implicit none

! Close all files opened by calling 'openfiles'
call closefiles()

! Level set:
$if ($LVLSET)

  $if ($RNS_LS)
  call rns_finalize_ls ()
  $endif
$endif

! Turbines:
$if ($TURBINES)
call turbines_finalize ()   ! must come before MPI finalize
$endif   

! SGS variable histograms
if (sgs_hist_calc) then
  call sgs_hist_finalize()
endif

! MPI:
$if ($MPI)
! First make sure everyone in has finished
call mpi_barrier( MPI_COMM_WORLD, ierr )
call mpi_finalize (ierr)
$endif

return
end subroutine finalize
