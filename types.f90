module types
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

type vec3d
  real(rprec) :: mag
  real(rprec), dimension(3) :: xyz
end type vec3d

type vec2d
  real(rprec) :: mag
  real(rprec), dimension(2) :: xy
end type vec2d

end module types
