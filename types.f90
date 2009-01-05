module types
implicit none

public
! rprec is used to specify precision 
!integer, parameter :: rprec = kind (1.e0)
integer, parameter :: rprec = kind (1.d0) 
!integer, parameter :: rprec = selected_real_kind (6)
!integer, parameter :: rprec = selected_real_kind (15)
end module types
