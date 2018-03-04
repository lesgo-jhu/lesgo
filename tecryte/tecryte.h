!///////////////////////////////////////////////////////////////////////////////
interface write_real_data
!///////////////////////////////////////////////////////////////////////////////

!*************************************************************
subroutine write_real_data_single_(fname, write_position, write_format, nvars, vars)
!*************************************************************
implicit none
character(*), intent(in) :: fname, write_position, write_format
integer, intent(in) :: nvars
real(4), intent(in), dimension(:) :: vars
end subroutine write_real_data_single_

!*************************************************************
subroutine write_real_data_double_(fname, write_position, write_format, nvars, vars)
!*************************************************************
implicit none
character(*), intent(in) :: fname, write_position, write_format
integer, intent(in) :: nvars
real(8), intent(in), dimension(:) :: vars
end subroutine write_real_data_double_

!*************************************************************
subroutine write_real_data_fid_single_(fid, write_format, nvars, vars)
!*************************************************************
implicit none
integer, intent(in) :: fid
character(*), intent(in) :: write_format
integer, intent(in) :: nvars
real(4), intent(in), dimension(:) :: vars
end subroutine write_real_data_fid_single_

!*************************************************************
subroutine write_real_data_fid_double_(fid, write_format, nvars, vars)
!*************************************************************
implicit none
integer, intent(in) :: fid
character(*), intent(in) :: write_format
integer, intent(in) :: nvars
real(8), intent(in), dimension(:) :: vars
end subroutine write_real_data_fid_double_

end interface write_real_data

!///////////////////////////////////////////////////////////////////////////////
interface write_real_data_1D
!///////////////////////////////////////////////////////////////////////////////

!*************************************************************
subroutine write_real_data_1D_single_(fname, write_position, write_format, nvars, &
  imax, vars, ibuff, x)
!*************************************************************
implicit none
character(*), intent(in) :: fname, write_position, write_format
integer, intent(in) :: nvars, imax
real(4), intent(in), dimension(nvars*imax) :: vars
integer, intent(in) :: ibuff
real(4), intent(in), dimension(:), optional :: x
end subroutine write_real_data_1D_single_

!*************************************************************
subroutine write_real_data_1D_double_(fname, write_position, write_format, nvars, &
  imax, vars, ibuff, x)
!*************************************************************
implicit none
character(*), intent(in) :: fname, write_position, write_format
integer, intent(in) :: nvars, imax
real(8), intent(in), dimension(nvars*imax) :: vars
integer, intent(in) :: ibuff
real(8), intent(in), dimension(:), optional :: x
end subroutine write_real_data_1D_double_

end interface write_real_data_1D

!///////////////////////////////////////////////////////////////////////////////
interface write_real_data_2D
!///////////////////////////////////////////////////////////////////////////////

!*************************************************************
subroutine write_real_data_2D_single_(fname, write_position, write_format, nvars, &
  imax, jmax, vars, ibuff, x, y)
!*************************************************************
implicit none
character(*), intent(in) :: fname, write_position, write_format
integer, intent(in) :: nvars, imax, jmax
real(4), intent(in), dimension(nvars*imax*jmax) :: vars
integer, intent(in) :: ibuff
real(4), intent(in), dimension(:), optional :: x, y
end subroutine write_real_data_2D_single_

!*************************************************************
subroutine write_real_data_2D_double_(fname, write_position, write_format, nvars, &
  imax, jmax, vars, ibuff, x, y)
!*************************************************************
implicit none
character(*), intent(in) :: fname, write_position, write_format
integer, intent(in) :: nvars, imax, jmax
real(8), intent(in), dimension(nvars*imax*jmax) :: vars
integer, intent(in) :: ibuff
real(8), intent(in), dimension(:), optional :: x, y
end subroutine write_real_data_2D_double_

end interface write_real_data_2D

!///////////////////////////////////////////////////////////////////////////////
interface write_real_data_3D
!//////////////////////////////////////////////////////////////////////////////

!*************************************************************
subroutine write_real_data_3D_single_(fname, write_position, write_format, nvars, &
  imax, jmax, kmax, vars, ibuff, x,y,z)
!*************************************************************
implicit none
character(*), intent(in) :: fname, write_position, write_format
integer, intent(in) :: nvars, imax, jmax, kmax
real(4), intent(in), dimension(nvars*imax*jmax*kmax) :: vars
integer, intent(in) :: ibuff
real(4), intent(in), dimension(:), optional :: x,y,z
end subroutine write_real_data_3D_single_

!*************************************************************
subroutine write_real_data_3D_double_(fname, write_position, write_format, nvars, &
  imax, jmax, kmax, vars, ibuff, x,y,z)
!*************************************************************
implicit none
character(*), intent(in) :: fname, write_position, write_format
integer, intent(in) :: nvars, imax, jmax, kmax
real(8), intent(in), dimension(nvars*imax*jmax*kmax) :: vars
integer, intent(in) :: ibuff
real(8), intent(in), dimension(:), optional :: x,y,z
end subroutine write_real_data_3D_double_

end interface write_real_data_3D

!///////////////////////////////////////////////////////////////////////////////
interface write_tecplot_header_xyline
!///////////////////////////////////////////////////////////////////////////////

!*************************************************************
subroutine write_tecplot_header_xyline_fname_(fname, write_position, var_list)
!*************************************************************
implicit none
character(*), intent(in) :: fname, write_position, var_list
end subroutine write_tecplot_header_xyline_fname_

!*************************************************************
subroutine write_tecplot_header_xyline_fid_(fid, var_list)
!*************************************************************
implicit none
integer, intent(in) :: fid
character(*), intent(in) :: var_list
end subroutine write_tecplot_header_xyline_fid_

end interface write_tecplot_header_xyline

!///////////////////////////////////////////////////////////////////////////////
interface 
!///////////////////////////////////////////////////////////////////////////////

!*************************************************************
function open_file_(fname, open_position, open_format) result(fid)
!*************************************************************
implicit none
character(*), intent(in) :: fname, open_position, open_format
integer :: fid
end function open_file_

end interface

!///////////////////////////////////////////////////////////////////////////////
interface write_tecplot_header_ND
!///////////////////////////////////////////////////////////////////////////////

!*************************************************************
subroutine write_tecplot_header_ND_(fname, write_position, nvars, &
  domain_size, var_list, zone, data_type, soln_time)
!*************************************************************
implicit none
character(*), intent(in) :: fname, write_position
integer, intent(in) :: nvars
integer, dimension(:), intent(in) :: domain_size
character(*), intent(in) :: var_list, zone
integer, intent(in) :: data_type
real(4), optional, intent(in) :: soln_time
end subroutine write_tecplot_header_ND_

end interface write_tecplot_header_ND
