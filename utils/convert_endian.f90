module convert_endian_base
implicit none
integer, parameter :: rp = kind (1.d0)

character(128) :: fbase

logical :: backup=.false.
logical :: integer_flag=.false.

! Not implemented yet
logical :: streamio = .false.

! 1 - little to big endian; 2 - big to little endian
integer :: iendian = 1

! Exact size of the variables and number of variables
integer :: nx=1, ny=1 , nz=1, nelem=1

! Number of records in the file
character(128) :: frecords
logical :: multiple_records = .false.
integer :: nrecords=1

type records
   logical :: integer_flag
   integer :: nelem
   integer, pointer, dimension(:) :: ivars
   real(rp), pointer, dimension(:) :: vars   
end type records

type( records ), allocatable, dimension(:) :: records_t

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
write(*,*) '  -f <recordsfile>         Specify input records file (required when multiple records'
write(*,*) '                           are in the input file)'
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

integer :: n
!---------------------------------------------------------------------


!//////////////////////////////////////
!/// READ DATA                      ///
!//////////////////////////////////////
! Load global definitions
call read_input_arguments()

! Initialize the data structures
if( multiple_records) then

   call read_records_file()

else

   allocate( records_t( 1 ) )

   records_t(1) % integer_flag = integer_flag
   records_t(1) % nelem = nx*ny*nz*nelem
   
   if( integer_flag ) then
      allocate(records_t(1) % ivars( records_t(1) % nelem ))
      records_t(1) % ivars = huge(1)
   else
      allocate(records_t(1) % vars( records_t(1) % nelem ))
      records_t(1) % vars=huge(1.0_rp)
   endif
 
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


if( backup ) then
   ! Save a copy
   fname = trim(adjustl(fbase)) // '.save'
   write(*,*) 'Backing up input file to : ', trim(adjustl(fname))
   call system('cp ' // trim(adjustl(fbase)) // ' ' // trim(adjustl(fname))) 
endif

! Load the data
write(*,*) 'Converting file : ', trim(adjustl(fbase))
open (1, file=fbase, form='unformatted', convert=trim(adjustl(read_endian)), &
     position='rewind', access=data_access)

do n = 1, nrecords

   if( records_t(n) % integer_flag ) then
      read(1) records_t(n) % ivars
   else
      read(1) records_t(n) % vars
   endif

enddo

! Rewind file to beginning
REWIND(1)

! Now over write the original file
write(*,*) 'Converting from ' // trim(adjustl(read_endian)) // ' to ' // trim(adjustl(write_endian))
!open (1, file=fbase, form='unformatted', convert=trim(adjustl(write_endian)), &
!     position='rewind', access=data_access)

do n = 1, nrecords

   if( records_t(n) % integer_flag ) then
      write(1) records_t(n) % ivars
   else
      write(1) records_t(n) % vars
   endif

enddo
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
     
     integer_flag = .true.

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

  elseif(arg(n) == '-f') then

     multiple_records =  .true.

     n=n+1
     ! Check that argument specified
     if( n > narg ) call print_error_message('records file not specified')
     read(arg(n),*) frecords

  elseif(arg(n) == '-e') then

     n=n+1
     ! Check that argument specified
     if( n > narg ) call print_error_message('endian conversion flag not specified')
     read(arg(n),*) iendian



  elseif( arg(n) == '-n' ) then 

     n=n+1
     ! Check that argument specified
     if( n > narg ) call print_error_message('number of elements (or variables)  not specified')
     read(arg(n),*) nelem


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

!***********************************************************************
subroutine read_records_file()
!***********************************************************************
use convert_endian_base, only : rp, records_t, nrecords, frecords
implicit none

integer :: n, nelem
logical :: integer_flag

integer :: istat

open(101,file=frecords,position='rewind')

nrecords=0
istat=0
do while (istat == 0)

  read(101,*,IOSTAT=istat) integer_flag, nelem
  if( istat == 0 ) nrecords = nrecords + 1
  
enddo

! Start back at beginning
rewind(101)

! Now allocate records
allocate( records_t( nrecords ) ) 
do n=1, nrecords

   read(101,*) integer_flag, nelem

   records_t( n ) % integer_flag = integer_flag
   records_t( n ) % nelem = nelem

   if( integer_flag ) then
      allocate( records_t( n ) % ivars( nelem ) )
      records_t(n) % ivars = huge(1)
   else
      allocate( records_t( n ) % vars( nelem ) )
      records_t( n ) % vars = huge(1.0_rp) 
   endif
enddo

return
end subroutine read_records_file

