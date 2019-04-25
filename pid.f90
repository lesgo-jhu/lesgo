!!
!!  Copyright (C) 2019  Johns Hopkins University
!!
!!  This file is part of lesgo.
!!
!!  lesgo is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License as published by
!!  the Free Software Foundation, either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  lesgo is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with lesgo.  If not, see <http://www.gnu.org/licenses/>.
!!

!*******************************************************************************
module pid_m
!*******************************************************************************
! PID-controller
use types, only : rprec
use messages
implicit none

private
public pid_t

type pid_t
    real(rprec) :: e = 0                    ! error (y - y_set)
    real(rprec) :: e_int = 0                ! integral of error
    real(rprec) :: e_prev = 0               ! previous error
    real(rprec) :: y_set = 0                ! output setoint
    real(rprec) :: Kp = 0, Ki = 0, Kd = 0   ! controller gains
contains
    procedure, private :: advance_set, advance_noset
    generic, public :: advance => advance_set, advance_noset
end type pid_t

interface pid_t
    module procedure :: constructor
end interface pid_t

contains

!*******************************************************************************
function constructor(Kp, Ki, Kd, y_set) result(this)
!*******************************************************************************
real(rprec), intent(in) :: Kp, Ki, Kd, y_set
type(pid_t) :: this

if (Kp >= 0) then
    this%Kp = Kp
else
    call error ('pid_t/constructor', 'Kp must be non-negative')
endif

if (Kp >= 0) then
    this%Ki = Ki
else
    call error ('pid_t/constructor', 'Ki must be non-negative')
end if

if (Kd >= 0) then
    this%Kd = Kd
else
    call error ('pid_t/constructor', 'Kd must be non-negative')
end if

this%y_set = y_set

end function constructor

!*******************************************************************************
function advance_noset(this, y, dt) result(u)
!*******************************************************************************
class(pid_t) :: this
real(rprec), intent(in) :: y, dt
real(rprec) :: u

! Save previous error
this%e_prev = this%e

! Calculate new error
this%e = this%y_set - y

! Integrate error
this%e_int = this%e_int + this%e*dt

! Calculate output
u = this%Kp*this%e + this%Ki*this%e_int + this%Kd*(this%e - this%e_prev)/dt

end function advance_noset

!*******************************************************************************
function advance_set(this, y, dt, y_set) result(u)
!*******************************************************************************
class(pid_t) :: this
real(rprec), intent(in) :: y, dt, y_set
real(rprec) :: u

this%y_set = y_set
u = this%advance_noset(y, dt)

end function advance_set

end module pid_m
