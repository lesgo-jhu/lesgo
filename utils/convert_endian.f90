module convert_endian_base
implicit none
integer, parameter :: rp = kind (1.d0)

character(128) :: fbase

logical :: backup=.false.
logical :: integer_data=.false.

! Not implemented yet
logical :: streamio = .false.

! 1 - little to big endian; 2 - big to little endian
integer :: iendian = 1

! Exact size of the variables and number of variables
integer :: nx=1, ny=1 , nz=1, nvar=1

end module convert_endian_base

!************************************************************
module messages
!************************************************************
!
!  This module contains functions for printing spectrum 
!  messages to screen.
!
!  Authors:
!    Jason Graham <jgraha8@gmail.com>
!
implicit none

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine print_error_message( message )
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

character(*), intent(in) :: message

write(*,*) 'convert-endian: ' // message
stop
end subroutine print_error_message

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine print_message( message )
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

character(*), intent(in) :: message

write(*,*) 'convert-endian: ' // message

end subroutine print_message

!************************************************************
subroutine print_help_message()
!************************************************************
!
! Subroutine used to print help message of convert-endian
!
! Authors:
!     Jason Graham <jgraha8@gmail.com>
!
write(*,*) 'Usage:'
write(*,*) '  convert-endian [options] -e {1|2} -n <elements> file'
write(*,*) ' '
write(*,*) 'Options:'
write(*,*) '  -B                       Backup input file (off by default)'
write(*,*) '  -i                       Specifies that data is of type integer (real default)'
write(*,*) '  -nx <nx>                 Set the number of dimensions in x-direction (default 1)'
write(*,*) '  -ny <ny>                 Set the number of dimensions in y-direction (default 1)'
write(*,*) '  -nz <nz>                 Set the number of dimensions in z-direction (default 1)'
write(*,*) ' '
write(*,*) 'Arguments:'
write(*,*) '  -e                       Endian conversion flag: 1 - little to big; 2 - big to little'
write(*,*) '  -n                       Set the number of data elements (or number of variables if '
write(*,*) '                           nx, ny, and/or nz is specified)'
write(*,*) ' '
write(*,*) 'NOTE: the order of arguments are arbitrary'
end subroutine print_help_message

end module messages

!***********************************************************************
program convert_endian
!***********************************************************************
use convert_endian_base
implicit none

character (128) :: fname
character (64) :: read_endian, write_endian, data_access

! Array which holds the data
real(rp), target, allocatable, dimension(:,:,:,:) :: vars


integer, target, allocatable, dimension(:,:,:,:) :: ivars
!---------------------------------------------------------------------


!//////////////////////////////////////
!/// READ DATA                      ///
!//////////////////////////////////////
call read_input_arguments()
!write(*,*) 'nx, ny, nz, nvar : ', nx, ny, nz, nvar
! Allocate space
if( integer_data ) then
   allocate(ivars(nx,ny,nz,nvar))
else
   allocate(vars(nx,ny,nz,nvar))
   vars=huge(1.0_rp)
endif

if(iendian == 1) then
  read_endian = 'little_endian'
  write_endian = 'big_endian'
elseif(iendian == 2) then
  read_endian = 'big_endian'
  write_endian = 'little_endian' 
else
  write(*,*) 'Error: incorrect endian specification.'
  stop
endif

if(streamio) then
   data_access='stream'
else
   data_access='sequential'
endif

! Load the data
write(*,*) 'Converting file : ', trim(adjustl(fbase))
open (1, file=fbase, form='unformatted', convert=trim(adjustl(read_endian)), &
     position='rewind', access=data_access)
if( integer_data ) then
   read(1) ivars
else
   read(1) vars
endif
close (1)

if( backup ) then
   ! Save a copy
   fname = trim(adjustl(fbase)) // '.save'
   write(*,*) 'Backing up input file to : ', trim(adjustl(fname))
   call system('cp ' // trim(adjustl(fbase)) // ' ' // trim(adjustl(fname))) 
endif

! Now over write the original file
write(*,*) 'Converting from ' // trim(adjustl(read_endian)) // ' to ' // trim(adjustl(write_endian))
open (1, file=fbase, form='unformatted', convert=trim(adjustl(write_endian)), &
     position='rewind', access=data_access)

if( integer_data ) then
   write(1) ivars
else
   write(1) vars
endif

close (1)

write(*,*) 'Conversion complete'

stop
end program convert_endian

!************************************************************
subroutine read_input_arguments()
!************************************************************
!
!  Subroutine used to read input arguments from command line
!
!  Written by: 
!    Adrien Thormann
!    Jason Graham <jgraha8@gmail.com>
!
  
use convert_endian_base
use messages
implicit none
integer :: n, narg
character(64), allocatable, dimension(:) :: arg
logical :: input_filename_read = .false.
logical :: exst = .false.

narg = command_argument_count()  !  Get total number of input arguments
if( narg == 0 ) call print_error_message( 'no input file' )

!  Read input arguments
allocate(arg(1:narg))
do n = 1, narg
  call get_command_argument(n, arg(n))
  arg(n) = trim(adjustl(arg(n)))
enddo

!
! March through command line argument array
!
n=0
input: do while (n < narg) 

  n=n+1

  if(arg(n) == '-B') then

     backup = .true.

  elseif( arg(n) == '-i') then
     
     integer_data = .true.

  elseif(arg(n) == '-e') then

     n=n+1
     ! Check that argument specified
     if( n > narg ) call print_error_message('endian conversion flag not specified')
     read(arg(n),*) iendian

  elseif( arg(n) == '-nx' ) then 

     n=n+1
     ! Check that argument specified
     if( n > narg ) call print_error_message('nx not specified')
     read(arg(n),*) nx

  elseif( arg(n) == '-ny' ) then 

     n=n+1
     ! Check that argument specified
     if( n > narg ) call print_error_message('ny not specified')
     read(arg(n),*) ny

  elseif( arg(n) == '-nz' ) then 

     n=n+1
     ! Check that argument specified
     if( n > narg ) call print_error_message('nz not specified')
     read(arg(n),*) nz

  elseif( arg(n) == '-n' ) then 

     n=n+1
     ! Check that argument specified
     if( n > narg ) call print_error_message('nvar not specified')
     read(arg(n),*) nvar


  elseif(arg(n) == '-h' .or. arg(n) == '--help') then
    call print_help_message()
    stop
  else 
     read(arg(n),*) fbase
     input_filename_read = .true.
  endif
enddo input

!  Check that an input file has been read
if( .not. input_filename_read ) call print_error_message('no input file')

! Check that input file exists
inquire(file=fbase,exist=exst)
if( .not. exst ) call print_error_message('input file does not exist')

return
end subroutine read_input_arguments



