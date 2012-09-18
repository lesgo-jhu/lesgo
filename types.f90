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
