!*************************************************************
subroutine write_real_data_single_(fname, write_position, write_format, nvars, vars)
!*************************************************************
! 
!  This subroutine writes the variables given by vars to the
!  specified file, fname. The number of variables can be arbitrary
!  but must be specified by nvars. The output will be of the form
!  
!  vars(1) vars(2) ... vars(nvars-1) vars(nvars)
!
!  The primary purpose of this routine is to write instantaneous data
!  or something similar to file
!
!  Inputs:
!  fname (char) - file to write to
!  write_position (char) - position to write in file : 'append' or 'rewind'
!  write_format (char) - Fotran format flag : 'formatted' or 'unformatted'
!  nvars (int) - number of variables contained in vars
!  vars (real,vector) - vector contaning variables to write
!
use tecryte
implicit none

character(*), intent(in) :: fname, write_position, write_format
integer, intent(in) :: nvars
real(4), intent(in), dimension(:) :: vars

character(*), parameter :: sub_name = mod_name // '.write_real_data'

character(64) :: frmt

!  Check if file exists
!inquire ( file=fname, exist=exst)
!if (.not. exst) call mesg(sub_name, 'Creating : ' // fname)

call check_write_position(write_position, sub_name)
call check_write_format(write_format, sub_name)
  
open (unit = 2,file = fname, status='unknown',form=write_format, &
  action='write',position=write_position)  
  
!  Write the data
select case(write_format)

case('formatted')

   ! Specify output format; may want to use a global setting
   ! Setting data to be written as single precision
   write (frmt, '("(",i0,"e15.7)")') nvars

   !  Write the data
   write(2,frmt) vars  

case('unformatted')

   write(2) vars

end select

close(2)
  
return
end subroutine write_real_data_single_

!*************************************************************
subroutine write_real_data_fid_single_(fid, write_format, nvars, vars)
!*************************************************************
! 
!  This subroutine writes the variables given by vars to the
!  specified file, fname. The number of variables can be arbitrary
!  but must be specified by nvars. The output will be of the form
!  
!  vars(1) vars(2) ... vars(nvars-1) vars(nvars)
!
!  The primary purpose of this routine is to write instantaneous data
!  or something similar to file
!
!  Inputs:
!  fid (integer) - file id to write to (assumed to alread be opened)
!  write_format (char) - Fotran format flag : 'formatted' or 'unformatted'
!  nvars (int) - number of variables contained in vars
!  vars (real,vector) - vector contaning variables to write
!
use tecryte
implicit none

integer, intent(in) :: fid
character(*), intent(in) :: write_format
integer, intent(in) :: nvars
real(4), intent(in), dimension(:) :: vars

character(*), parameter :: sub_name = mod_name // '.write_real_data'

character(64) :: frmt

!  Check if file exists
!inquire ( file=fname, exist=exst)
!if (.not. exst) call mesg(sub_name, 'Creating : ' // fname)

call check_write_format(write_format, sub_name)
  
!  Write the data
select case(write_format)

case('formatted')

   ! Specify output format; may want to use a global setting
   ! Setting data to be written as single precision
   write (frmt, '("(",i0,"e15.7)")') nvars

   !  Write the data
   write(fid,frmt) vars  

case('unformatted')

   write(fid) vars

end select

return
end subroutine write_real_data_fid_single_


!*************************************************************
subroutine write_real_data_double_(fname, write_position, write_format, nvars, vars)
!*************************************************************
! 
!  This subroutine writes the variables given by vars to the
!  specified file, fname. The number of variables can be arbitrary
!  but must be specified by nvars. The output will be of the form
!  
!  vars(1) vars(2) ... vars(nvars-1) vars(nvars)
!
!  The primary purpose of this routine is to write instantaneous data
!  or something similar to file
!
!  Inputs:
!  fname (char) - file to write to
!  write_position (char) - position to write in file : 'append' or 'rewind'
!  write_format (char) - Fotran format flag : 'formatted' or 'unformatted'
!  nvars (int) - number of variables contained in vars
!  vars (real,vector) - vector contaning variables to write
!
use tecryte
implicit none

character(*), intent(in) :: fname, write_position, write_format
integer, intent(in) :: nvars
real(8), intent(in), dimension(:) :: vars

character(*), parameter :: sub_name = mod_name // '.write_real_data'

character(64) :: frmt

!  Check if file exists
!inquire ( file=fname, exist=exst)
!if (.not. exst) call mesg(sub_name, 'Creating : ' // fname)

call check_write_position(write_position, sub_name)
call check_write_format(write_format, sub_name)
  
open (unit = 2,file = fname, status='unknown',form=write_format, &
  action='write',position=write_position)  
  
!  Write the data
select case(write_format)

case('formatted')

   ! Specify output format; may want to use a global setting
   ! Setting data to be written as single precision
   write (frmt, '("(",i0,"e15.7)")') nvars

   !  Write the data
   write(2,frmt) vars  

case('unformatted')

   write(2) vars

end select

close(2)
  
return
end subroutine write_real_data_double_

!*************************************************************
subroutine write_real_data_fid_double_(fid, write_format, nvars, vars)
!*************************************************************
! 
!  This subroutine writes the variables given by vars to the
!  specified file, fname. The number of variables can be arbitrary
!  but must be specified by nvars. The output will be of the form
!  
!  vars(1) vars(2) ... vars(nvars-1) vars(nvars)
!
!  The primary purpose of this routine is to write instantaneous data
!  or something similar to file
!
!  Inputs:
!  fid (integer) - file id to write to (assumed to alread be opened)
!  write_format (char) - Fotran format flag : 'formatted' or 'unformatted'
!  nvars (int) - number of variables contained in vars
!  vars (real,vector) - vector contaning variables to write
!
use tecryte
implicit none

integer, intent(in) :: fid
character(*), intent(in) :: write_format
integer, intent(in) :: nvars
real(8), intent(in), dimension(:) :: vars

character(*), parameter :: sub_name = mod_name // '.write_real_data'

character(64) :: frmt

!  Check if file exists
!inquire ( file=fname, exist=exst)
!if (.not. exst) call mesg(sub_name, 'Creating : ' // fname)

call check_write_format(write_format, sub_name)
  
!  Write the data
select case(write_format)

case('formatted')

   ! Specify output format; may want to use a global setting
   ! Setting data to be written as single precision
   write (frmt, '("(",i0,"e15.7)")') nvars

   !  Write the data
   write(fid,frmt) vars  

case('unformatted')

   write(fid) vars

end select

return
end subroutine write_real_data_fid_double_


