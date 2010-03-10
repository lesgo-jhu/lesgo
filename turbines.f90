module turbines
use types,only:rprec
use param
use stat_defs, only:turbine_t
!use test_filtermodule --subroutines are not actually part of this module (yet)

implicit none

save
private

public :: turbines_init, turbines_forcing


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine turbines_init(ifilter,alpha)
!initialize filter
!locate applicable nodes (turbine_nodes)

!Filter type: 1->cut off 2->Gaussian 3->Top-hat
integer,parameter::ifilter
real (rprec) :: alpha

call turbines_nodes()

end subroutine turbines_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine turbines_forcing()
!calculate forcing
!filter forcing
!apply forcing
end subroutine turbines_forcing

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine turbines_nodes()
! use the following to give all nodes 0/1 - where 1 means there is a turbine there
	!turbine_t%dia    
	!turbine_t%height  
	!turbine_t%nloc       
	!turbine_t%xloc(1)  !and more
	!turbine_t%yloc(1)  !and more
end subroutine turbines_nodes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module turbines
