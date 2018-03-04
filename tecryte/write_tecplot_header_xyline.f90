!*************************************************************
subroutine write_tecplot_header_xyline_fname_(fname, write_position, var_list)
!*************************************************************
!  The purpose of this routine is to write Tecplot header information
!  for xy-line formatted files
! 
!  Inputs:
!  fname (char)     - file name to write to
!  write_position (char) - position in file to write data: append or rewind
!  var_list	(char)  - string containing variable names: Ex. "x", "u"
!
use tecryte
implicit none

character(*), intent(in) :: fname, write_position, var_list

character(*), parameter :: sub_name = mod_name // '.write_tecplot_header_xyline'

!  Check if write position has been specified correctly
call check_write_position(write_position, sub_name)

open (unit = 2,file = fname, status='unknown',form='formatted', &
  action='write',position=write_position)

write(2,'(1a)') 'variables = ' // var_list

close(2)

return
end subroutine write_tecplot_header_xyline_fname_

!*************************************************************
subroutine write_tecplot_header_xyline_fid_(fid, var_list)
!*************************************************************
!  The purpose of this routine is to write Tecplot header information
!  for xy-line formatted files
! 
!  Inputs:
!  fid (int)     - file identifier to write to
!  var_list	(char)  - string containing variable names: Ex. "x", "u"
!
use tecryte
implicit none

integer, intent(in) :: fid
character(*), intent(in) ::  var_list

character(*), parameter :: sub_name = mod_name // '.write_tecplot_header_xyline_fid'

write(fid,'(1a)') 'variables = ' // var_list

return
end subroutine write_tecplot_header_xyline_fid_

