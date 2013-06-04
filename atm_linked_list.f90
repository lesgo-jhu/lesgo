!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
module atm_linked_list
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This module is used to create dynamic array

! Imported modules
use atm_input_util, only : rprec! Precision of real numbers

implicit none

! Define the type used to create a dynamic list
! An element pointer of the same type is created recusrisvely 
! This allows dynamic allocation
type :: DynamicList_
    integer :: id                   ! This identifies the element in the 
                                    ! list 1, 2, 3, 4, ...
    ! Variables that may be used
    character(128) :: name_id       ! Name identifier
    real(rprec) :: vector(3)        ! A vector (a,b,c)
    real(rprec) :: real_number      ! A real number
    integer :: integer_number       ! An integer
    
    type(DynamicList_), pointer :: prev => null()
    type(DynamicList_), pointer :: next => null()
end type DynamicList_

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine initializeDynamicListVector(DynamicList,vector)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialize a dynamic list vector
    implicit none

    type(DynamicList_), pointer, intent(inout) :: DynamicList
    real(rprec), intent(in) :: vector(3)

    allocate( DynamicList )
    DynamicList % id = 1
    DynamicList % vector(1) = vector(1)        ! Modify vector element 1
    DynamicList % vector(2) = vector(2)        ! Modify vector element 2
    DynamicList % vector(3) = vector(3)        ! Modify vector element 3
end subroutine initializeDynamicListVector

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine initializeDynamicListCharacter(DynamicList,name_id)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    implicit none

    type(DynamicList_), pointer, intent(inout) :: DynamicList
    character(128) :: name_id

    allocate( DynamicList )
    DynamicList % id = 1
    DynamicList % name_id = name_id        ! Modify character
end subroutine initializeDynamicListCharacter

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine appendVector(DynamicList,vector)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    implicit none

    real(rprec), intent(in) :: vector(3)
    type(DynamicList_), pointer, intent(inout) :: DynamicList

    allocate(DynamicList % next)
    nullify(DynamicList % next % next)
    DynamicList % next % id = DynamicList % id +1     ! Adds counter
    DynamicList % next % vector(1) = vector(1)        ! Modify vector element 1
    DynamicList % next % vector(2) = vector(2)        ! Modify vector element 2
    DynamicList % next % vector(3) = vector(3)        ! Modify vector element 3

end subroutine appendVector

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine appendCharacter(DynamicList,name_id)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Will append a character to a dynamic list

    implicit none

    character(128) :: name_id
    type(DynamicList_), pointer, intent(inout) :: DynamicList

    allocate(DynamicList % next)
    nullify(DynamicList % next % next)
    DynamicList % next % id = DynamicList % id +1     ! Adds counter
    DynamicList % next % name_id = name_id        ! Modify vector element 1

end subroutine appendCharacter

end module atm_linked_list
