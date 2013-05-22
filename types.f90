!!
!!  Copyright 2009,2010,2011,2012 Johns Hopkins University
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
module types
!*********************************************************************
!
! This module provides generic types
!
implicit none

public
! rprec is used to specify precision
$if(DBLPREC)
integer, parameter :: rprec = kind (1.d0)
$else
integer, parameter :: rprec = kind (1.0)
$endif
 
!integer, parameter :: rprec = kind (1.e0)
!integer, parameter :: rprec = selected_real_kind (6)
!integer, parameter :: rprec = selected_real_kind (15)

type vec3d_t
  real(rprec) :: mag
  real(rprec), dimension(3) :: xyz
end type vec3d_t

type vec2d_t
  real(rprec) :: mag
  real(rprec), dimension(2) :: xy
end type vec2d_t

type point3D_t
  real(rprec), dimension(3) :: xyz
end type point3D_t



end module types
