!**********************************************************************
module cylinder_skew_defs
!**********************************************************************
type cylinder_skew
  integer :: ngen
  integer, dimension(:), allocatable :: igen, kbottom, kbottom_inside, ktop, ktop_inside, lun
  double precision, dimension(:), allocatable :: dz_bottom, dz_top
end type cylinder_skew

type(cylinder_skew) :: cylinder_skew_t

end module cylinder_skew_defs