module tecryte
implicit none

save
public
!public write_real_data, write_real_data_1D, write_real_data_2D, &
!       write_real_data_3D, write_tecplot_header_xyline, &
!       write_tecplot_header_ND

!private rprec

!integer, parameter :: rprec = kind(1.0d0)
character(*), parameter :: mod_name = 'tecryte'

interface strcat
  module procedure strcat_aa, strcat_ai, strcat_ar
end interface

! Begin file identifier at 1000
integer :: fid_open_file = 1000

character(*), parameter :: int_fmt='(i0)'
character(*), parameter :: real_fmt='(f9.4)'

! Format specifier used for writing data elements. The structure of
! the specifier is 
!
!    <N>e<w>.<d>
!
! where <N> is the number of entries per line, <w> the width of the
! line, and <d> the depth. The larger the value of <N> the better in
! most cases, which wll allow for larger blocks of data to be written
! at a given time. Just be sure that <N>*<w> < 4000 as Tecplot only
! supports up to 4000 characters per line.
character(*), parameter :: data_format_single='(128e15.7)'

! double precison data written to file as single precision in order to
! 1) decrease overall write time and 2) save space since the extra
! bits probably aren't necessary.
character(*), parameter :: data_format_double=data_format_single
!character(*), parameter :: data_format_double='(128e23.15)'
contains

!*************************************************************
subroutine tecplot_data_type_string(nvars, data_type, sub_name, tec_dt_str)
!*************************************************************
implicit none

integer, intent(in) ::  nvars, data_type
character, intent(IN) :: sub_name
character(120), intent(OUT) :: tec_dt_str
character(7) :: tec_dt

!character(*), parameter :: sub_name = mod_name // '.tecplot_data_type_string'

integer :: n

!  Specify single or double precision
if(data_type == 1) then
  tec_dt = ' SINGLE'
elseif(data_type == 2) then
  tec_dt = ' DOUBLE'
else
  
  write(*,*) 'Called from ', sub_name
  write(*,*) 'Incorrect data type : ', data_type
  stop

endif

!  Create DT string
tec_dt_str = 'DT=('
do n=1, nvars
  call strcat(tec_dt_str,tec_dt)
enddo
call strcat(tec_dt_str,')')

return
end subroutine tecplot_data_type_string

!*************************************************************
subroutine check_write_position(write_pos, sub_name)
!*************************************************************
!  This subroutine checks whether the write position has been
!  specified properly
!
implicit none

character(*), intent(IN) :: write_pos, sub_name

!  Check if write position has been specified correctly
select case(write_pos)
  case('append')
  !  do nothing
  case('rewind')
  !  do nothing
  case default
    write(*,*) 'Called from ', sub_name
    write(*,*) 'Incorrect write position : ' // write_pos
    stop
end select

return
end subroutine check_write_position

!*************************************************************
subroutine check_write_format(write_fmt, sub_name)
!*************************************************************
!  This subroutine checks whether the write position has been
!  specified properly
!
implicit none

character(*), intent(IN) :: write_fmt, sub_name

!  Check if write position has been specified correctly
select case(write_fmt)
  case('formatted')
  !  do nothing
  case('unformatted')
  !  do nothing
  case default
    write(*,*) 'Called from ', sub_name
    write(*,*) 'Incorrect write format : ' // write_fmt
    stop
end select

return
end subroutine check_write_format

!**********************************************************************
integer function buff_indx(i,imax)
!**********************************************************************
!  This function returns the physical index associated with the buffer 
!  region for the specified i and imax. 
!  For i = imax + 1 -> 1 is returned otherwise i is returned
implicit none

integer, intent(in) :: i,imax

if(i == imax + 1) then
  buff_indx = 1
else
  buff_indx = i
endif

return
end function buff_indx

!**********************************************************************
subroutine strcat_aa(str1, str2)
!**********************************************************************
implicit none

character(*), intent(INOUT) :: str1
character(*), intent(IN) :: str2

str1 = trim(adjustl(str1)) // str2

return
end subroutine strcat_aa

!**********************************************************************
subroutine strcat_ai(str1, i1)
!**********************************************************************
implicit none

character(*), intent(INOUT) :: str1
integer, intent(IN) :: i1
character(120) :: str2

write (str2, int_fmt) i1

call strcat(str1,trim(adjustl(str2)))

end subroutine strcat_ai

!**********************************************************************
subroutine strcat_ar(str1, r1)
!**********************************************************************
implicit none

character(*), intent(INOUT) :: str1
real(8), intent(IN) :: r1
character(120) :: str2

write (str2,real_fmt) r1

call strcat(str1,trim(adjustl(str2)))

return
end subroutine strcat_ar

end module tecryte

