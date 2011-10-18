!*********************************************************************
module types
!*********************************************************************
!
! This module provides generic types and routines that act on them.
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

type vec3d
  real(rprec) :: mag
  real(rprec), dimension(3) :: xyz
end type vec3d

type vec2d
  real(rprec) :: mag
  real(rprec), dimension(2) :: xy
end type vec2d

type point3D
  real(rprec), dimension(3) :: xyz
end type point3D

type clock_type
   real(rprec) :: start
   real(rprec) :: end
end type clock_type

!---------------------------------------------------------------------
contains
!---------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function clock_time( clock_t ) result( time )
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

type(clock_type), intent(in) :: clock_t
real(rprec) :: time

time = clock_t % end - clock_t % start

return
end function clock_time

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine clock_start( clock_t )
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

type(clock_type), intent(inout) :: clock_t

call cpu_time( clock_t % start )

return
end subroutine clock_start

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine clock_end( clock_t )
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

type(clock_type), intent(inout) :: clock_t

call cpu_time( clock_t % end )

return
end subroutine clock_end

end module types
