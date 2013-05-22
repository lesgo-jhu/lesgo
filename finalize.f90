!!
!!  Copyright 2012 Johns Hopkins University
!!
!!  Licensed under the Apache License, Version 2.0 (the "License"); you may not 
!!  use this file except in compliance with the License. You may obtain a copy of
!!  the License at:
!!
!!    http://www.apache.org/licenses/LICENSE-2.0
!!
!!  Unless required by applicable law or agreed to in writing, software 
!!  distributed under the License is distributed on an "AS IS" BASIS, WITHOUT 
!!  WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the 
!!  License for the specific language governing permissions and limitations under
!!  the License.
!!

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
