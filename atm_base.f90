!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Written by: 
!!
!!   Luis 'Tony' Martinez <tony.mtos@gmail.com> (Johns Hopkins University)
!!
!!   Copyright (C) 2012-2013, Johns Hopkins University
!!
!!   This file is part of The Actuator Turbine Model Library.
!!
!!   The Actuator Turbine Model is free software: you can redistribute it 
!!   and/or modify it under the terms of the GNU General Public License as 
!!   published by the Free Software Foundation, either version 3 of the 
!!   License, or (at your option) any later version.
!!
!!   The Actuator Turbine Model is distributed in the hope that it will be 
!!   useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
!!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!   GNU General Public License for more details.
!!
!!   You should have received a copy of the GNU General Public License
!!   along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
module atm_base
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This module provides basic functionalities to the actuator turbine model
! Real precision variable and dynamic allocation types are stored here

implicit none

! Precision of real numbers
integer, parameter :: rprec = kind (1.d0)

!private
!public  DynamicList_, initializeDynamicListVector, 
!   & initializeDynamicListCharacter, appendVector, appendCharacter

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
subroutine DynamicListGoTo(DynamicList,i)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Goes to element i on the list
    implicit none
    integer, intent(in) :: i 
    type(DynamicList_), pointer, intent(inout) :: DynamicList

    do while (i .lt. DynamicList % id )
        DynamicList => DynamicList % next
    enddo
    do while (i .gt. DynamicList % id )
        DynamicList => DynamicList % prev
    enddo

end subroutine DynamicListGoTo

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
subroutine appendVector(DynamicList,vector)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    implicit none

    real(rprec), intent(in) :: vector(3)
    type(DynamicList_), pointer, intent(inout) :: DynamicList

    allocate(DynamicList % prev)
    DynamicList % prev => DynamicList

    allocate(DynamicList % next)

    DynamicList  % next% id = DynamicList % id +1     ! Adds counter
    DynamicList  % next% vector(1) = vector(1)        ! Modify vector element 1
    DynamicList  % next% vector(2) = vector(2)        ! Modify vector element 2
    DynamicList  % next% vector(3) = vector(3)        ! Modify vector element 3

    DynamicList => DynamicList % next
end subroutine appendVector

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
subroutine appendCharacter(DynamicList,name_id)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Will append a character to a dynamic list

    implicit none

    character(128) :: name_id
    type(DynamicList_), pointer, intent(inout) :: DynamicList

    allocate(DynamicList % next)
    DynamicList % next % prev => DynamicList

    DynamicList % next % id = DynamicList % id +1     ! Adds counter
    DynamicList % next % name_id = name_id        ! Modify character element 1
    DynamicList => DynamicList % next

end subroutine appendCharacter

end module atm_base
