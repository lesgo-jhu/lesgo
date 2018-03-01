!!
!!  Copyright (C) 2010-2017  Johns Hopkins University
!!
!!  This file is part of lesgo.
!!
!!  lesgo is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License as published by
!!  the Free Software Foundation, either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  lesgo is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with lesgo.  If not, see <http://www.gnu.org/licenses/>.
!!

!*******************************************************************************
module open_file_fid_mod
!*******************************************************************************
implicit none

save
public

! Begin file identifier at 1000
integer :: global_fid_count = 1000

contains

!*******************************************************************************
function open_file_fid(fname, open_position, open_format) result(fid)
!*******************************************************************************
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
implicit none

character(*), intent(in) :: fname, open_position, open_format
character(*), parameter :: sub_name = 'open_file_fid'

integer :: fid

! Check write position
select case(open_position)
    case('append')
        !  do nothing
    case('rewind')
        !  do nothing
    case default
        write(*,*) 'Called from ', sub_name
        write(*,*) 'Incorrect write position : ' // open_position
        stop
end select

! Check write format
select case(open_format)
    case('formatted')
        !  do nothing
    case('unformatted')
        !  do nothing
    case default
        write(*,*) 'Called from ', sub_name
        write(*,*) 'Incorrect write format : ' // open_format
        stop
end select

! Update global fid counter
global_fid_count = global_fid_count + 1
! Set current fid to global fid value
fid = global_fid_count

open (unit = fid, file = fname, status='unknown',form=open_format,             &
    position=open_position)

end function open_file_fid

end module open_file_fid_mod
