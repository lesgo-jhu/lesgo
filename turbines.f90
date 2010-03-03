module turbines
use types,only:rprec
use param
!use test_filtermodule --subroutines are not actually part of this module (yet)

implicit none

save
private

public :: turbines_forcing

!Filter type: 1->cut off 2->Gaussian 3->Top-hat
integer,parameter::turbines_filter_type=2
real (rprec) :: turbines_filter_size=1.5


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine turbines_forcing()
!locate applicable nodes (another subroutine, probably)
!calculate forcing
!filter forcing
!apply forcing
end subroutine

end module turbines
