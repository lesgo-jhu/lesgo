!*************************************************************
subroutine write_tecplot_header_nd_(fname, write_position, nvars, &
  domain_size, var_list, zone, data_type, soln_time)
!*************************************************************
!  The purpose of this routine is to write Tecplot header information
!  for 1D, 2D, and 3D data files.
!
!  NOTE: domain_size needs to be specified as a vector even for 1D data.
!  Ex: 
!    1D : (/ Nx /)
!    2D : (/ Nx, Ny /)
!    3D : (/ Nx, Ny, Nz /)
! 
!  Inputs:
!  fname (char)                  - file name to write to
!  write_position (char)         - position in file to write data: append or rewind
!  nvars (int)                   - number of variables
!  domain_size (int, vector)     - vector containing the diminsions of the data.
!  var_list	(char)           - string containing variable names: Ex. '"x", "u"'
!  zone (char)                   - zone name/id
!  date_type (int) 	         - specify Tecplot data type (precision): 1 - single, 2 - double
!  soln_time (real(4), optional) - time stamp
!
use tecryte
implicit none

character(*), intent(in) :: fname, write_position
integer, intent(in) :: nvars
integer, dimension(:), intent(in) :: domain_size
character(*), intent(in) :: var_list, zone
integer, intent(in) :: data_type
real(4), optional, intent(in) :: soln_time

character(*), parameter :: sub_name = mod_name // '.write_tecplot_header_nd'
character(120) :: tec_dt_str, tec_dat_str
integer :: ndims

!  Check if write position has been specified correctly
call check_write_position(write_position, sub_name)

!  Get number of dimensions for data
ndims = size(domain_size,1)

if(ndims == 1) then
  write(tec_dat_str,"(1a,1a,1a,i9,1a,i9)") 'ZONE T="', &
    zone,'", DATAPACKING=BLOCK, i=', domain_size(1)
elseif(ndims == 2) then
  write(tec_dat_str,"(1a,1a,1a,i9,1a,i9,1a,i9)") 'ZONE T="', &
    zone,'", DATAPACKING=BLOCK, i=', domain_size(1),', j=', domain_size(2)
elseif(ndims == 3) then
  write(tec_dat_str,"(1a,1a,1a,i9,1a,i9,1a,i9,1a,i9)") 'ZONE T="', &
    zone,'", DATAPACKING=BLOCK, i=', domain_size(1),', j=', domain_size(2),', k=', domain_size(3)
else
  
  write(*,*) 'Called from ', sub_name
  write(*,*) 'Incorrect number of dimensions : ', ndims
  stop

endif

!  Create Tecplot DT string
call tecplot_data_type_string(nvars, data_type, sub_name, tec_dt_str)

open (unit = 2,file = fname, status='unknown',form='formatted', &
  action='write',position=write_position)

!  Write variable list
write(2,'(1a)') 'variables = ' // var_list
!  Write data layout size information
write(2,'(1a)') tec_dat_str
!  Write Tecplot data type for each variable
write(2,'(1a)') tec_dt_str

if (present (soln_time)) then
write(2,'(1a,e15.6)') 'solutiontime=', soln_time
endif

close(2)

  
return
end subroutine write_tecplot_header_nd_
