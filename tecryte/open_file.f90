!*************************************************************
function open_file_(fname, open_position, open_format) result(fid)
!*************************************************************
! 
!  This function opens the file 'fname' and returns the file identifier
!  'fid'. It provides a method for openfiles while managing the file identifiers
!  in the backend. There will be no need to explicitly manage these values.
!  
!  Inputs:
!  fname (char)         - file to open to
!  open_position (char) - position to open in file : 'append' or 'rewind'
!  open_format (char)   - Fortran format flag : 'formatted' or 'unformatted'
!
!  Returns: fid (int)  - file identifier
!
use tecryte
implicit none

character(*), intent(in) :: fname, open_position, open_format
character(*), parameter :: sub_name = mod_name // '.open_file'

integer :: fid

call check_write_position(open_position, sub_name)
call check_write_format(open_format, sub_name)

! Update global fid counter
fid_open_file = fid_open_file + 1
! Set current fid to global fid value
fid = fid_open_file
  
open (unit = fid, file = fname, status='unknown',form=open_format, &
  position=open_position)  

return
end function open_file_
