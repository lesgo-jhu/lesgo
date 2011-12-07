module convert_endian_base
implicit none
integer, parameter :: rp = kind (1.d0)

character(128) :: fbase
logical :: backup=.false.

! 1 - little to big endian; 2 - big to little endian
integer :: iendian = 1

! Exact size of the variables and number of variables
integer :: nx=32, ny=32 , nz=32, nvar=3

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
write(*,*) '  convert-endian [options] -e {1|2} -nx <nx> -ny <ny> -nz <nz> -n <n> file'
write(*,*) ' '
write(*,*) 'Options:'
write(*,*) '  -b                       Backup input file'
write(*,*) ' '
write(*,*) 'Arguments:'
write(*,*) '  -e                       Endian conversion flag: 1 - little to endian; 2 - endian to little'
write(*,*) '  -nx                      Set the number of dimensions in x-direction'
write(*,*) '  -ny                      Set the number of dimensions in y-direction'
write(*,*) '  -nz                      Set the number of dimensions in z-direction'
write(*,*) '  -n                       Set the number of variables'
write(*,*) ' '
write(*,*) 'NOTE: the order of arguments are arbitrary'
end subroutine print_help_message

end module messages

program convert_endian
use convert_endian_base
implicit none

!character (*), parameter :: path='./output'
!character (*), parameter :: fbase = path // 'uvw.$FITER.out'
!character (*), parameter :: MPI_suffix = '.c'  !--must be proc number

!logical, parameter :: MPI = .true.


character (64) :: fmt
character (128) :: fname
character (64) :: read_endian, write_endian

integer :: ip
integer :: lbz, ubz
integer :: i, j, k

logical :: fexist

! Array which holds the data
real(rp), allocatable, dimension(:,:,:,:) :: vars

!---------------------------------------------------------------------

!//////////////////////////////////////
!/// READ DATA                      ///
!//////////////////////////////////////
call read_input_arguments()
!write(*,*) 'nx, ny, nz, nvar : ', nx, ny, nz, nvar
! Allocate space
allocate(vars(nx,ny,nz,nvar))
vars=0.0

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

! Load the data
write(*,*) 'Converting file : ', trim(adjustl(fbase))
open (1, file=fbase, form='unformatted', convert=trim(adjustl(read_endian)))
read(1) vars
close (1)

if( backup ) then
   ! Save a copy
   fname = trim(adjustl(fbase)) // '.save'
   write(*,*) 'Backing up input file to : ', trim(adjustl(fname))
   call system('cp ' // trim(adjustl(fbase)) // ' ' // trim(adjustl(fname))) 
endif

! Now over write the original file
write(*,*) 'Converting from ' // trim(adjustl(read_endian)) // ' to ' // trim(adjustl(write_endian))
open (1, file=fbase, form='unformatted', convert=trim(adjustl(write_endian)))
write (1) vars
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

  if(arg(n) == '-b') then

     backup = .true.

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



