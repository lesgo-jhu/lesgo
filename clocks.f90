!!
!!  Copyright 2011,2012 Johns Hopkins University
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

!*********************************************************************
module clocks
!*********************************************************************
!
! This module provides the clock data type (object) and the
! subroutines/functions that act on instances of the clock data type.
!
use types, only : rprec
implicit none

save 
private

public clock_t, &
     clock_start, &
     clock_stop

type clock_t
   real(rprec) :: start
   real(rprec) :: stop
   real(rprec) :: time
end type clock_t

!---------------------------------------------------------------------
contains
!---------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine clock_start( this )
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

type(clock_t), intent(inout) :: this

call cpu_time( this % start )

return

end subroutine clock_start

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine clock_stop( this )
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

type(clock_t), intent(inout) :: this

call cpu_time( this % stop )

! Compute the clock time
this % time = this % stop - this % start

return

end subroutine clock_stop

end module clocks