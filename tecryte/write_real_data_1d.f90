!*************************************************************
subroutine write_real_data_1D_single_(fname, write_position, write_format, nvars, &
  imax, vars, ibuff, x)
!*************************************************************
! 
!  This subroutine writes the variables given by vars to the
!  specified file, fname. The number of variables can be arbitrary
!  but must be specified by nvars. An example in which vars 
!  is to be specified as:
!    (/ u, v, w /)
!  
!  Inputs:
!  fname (char) - file to write to
!  write_position (char) - postition to write in file : 'append' or 'rewind'
!  write_format (char) - Fotran format flag : 'formatted' or 'unformatted'
!  nvars (int) - number of variables contained in vars
!  imax (int) - size of 1st dimension of variables in vars
!  vars (real, vector) - vector contaning variables to write
!  ibuff (int) - flag for adding buffer region due to periodicity
!     0 - no buffer region
!     1 - buffer on i direction (i = 1 = imax + 1)
!  x (real,vector,optional) - vector containing x coordinates 
!
use tecryte
implicit none

character(*), intent(in) :: fname, write_position, write_format
integer, intent(in) :: nvars, imax
real(4), intent(in), dimension(nvars*imax) :: vars
integer, intent(in) :: ibuff
real(4), intent(in), dimension(:), optional :: x

logical :: coord_pres

character(*), parameter :: sub_name = mod_name // '.write_real_data_1D'

integer :: i,n
integer :: i0, imax_buff

real(4), allocatable, dimension(:,:) :: vars_1d

call check_write_position(write_position, sub_name)
call check_write_format(write_format, sub_name)

!  Check if spatial coordinates specified
if(present(x)) coord_pres = .true.
 
if( ibuff == 0 ) then

  imax_buff = imax
  
elseif( ibuff == 1 ) then

  imax_buff = imax + 1
 
else 

  write(*,*) 'ibuff not specified correctly'
  stop
  
endif

allocate(vars_1d(imax_buff,nvars))  

do n=1,nvars
 
  do i = 1, imax_buff

    i0 = buff_indx(i,imax)
      
    vars_1d(i,n) = vars( (n-1)*imax + i0 )
    
  enddo
  
enddo

open (unit = 2,file = fname, status='unknown',form=write_format, &
  action='write',position=write_position)
   
!  Write the data
select case(write_format)
  case('formatted')
  
    if (coord_pres) then
  
      write(2,data_format_single) x, vars_1d

   else

      write(2,data_format_single) vars_1d
  
    endif
    
  case('unformatted')
  
    if (coord_pres) then

       write(2) x, vars_1d

    else

       write(2) vars_1d

    endif
    
end select

close(2)

deallocate ( vars_1d )

return
end subroutine write_real_data_1D_single_
!*************************************************************
subroutine write_real_data_1D_double_(fname, write_position, write_format, nvars, &
  imax, vars, ibuff, x)
!*************************************************************
! 
!  This subroutine writes the variables given by vars to the
!  specified file, fname. The number of variables can be arbitrary
!  but must be specified by nvars. An example in which vars 
!  is to be specified as:
!    (/ u, v, w /)
!  
!  Inputs:
!  fname (char) - file to write to
!  write_position (char) - postition to write in file : 'append' or 'rewind'
!  write_format (char) - Fotran format flag : 'formatted' or 'unformatted'
!  nvars (int) - number of variables contained in vars
!  imax (int) - size of 1st dimension of variables in vars
!  vars (real, vector) - vector contaning variables to write
!  ibuff (int) - flag for adding buffer region due to periodicity
!     0 - no buffer region
!     1 - buffer on i direction (i = 1 = imax + 1)
!  x (real,vector,optional) - vector containing x coordinates 
!
use tecryte
implicit none

character(*), intent(in) :: fname, write_position, write_format
integer, intent(in) :: nvars, imax
real(8), intent(in), dimension(nvars*imax) :: vars
integer, intent(in) :: ibuff
real(8), intent(in), dimension(:), optional :: x

logical :: coord_pres

character(*), parameter :: sub_name = mod_name // '.write_real_data_1D'

integer :: i,n
integer :: i0, imax_buff

real(8), allocatable, dimension(:,:) :: vars_1d

call check_write_position(write_position, sub_name)
call check_write_format(write_format, sub_name)

!  Check if spatial coordinates specified
if(present(x)) coord_pres = .true.
 
if( ibuff == 0 ) then

  imax_buff = imax
  
elseif( ibuff == 1 ) then

  imax_buff = imax + 1
 
else 

  write(*,*) 'ibuff not specified correctly'
  stop
  
endif

allocate(vars_1d(imax_buff,nvars))  

do n=1,nvars
 
  do i = 1, imax_buff

    i0 = buff_indx(i,imax)
      
    vars_1d(i,n) = vars( (n-1)*imax + i0 )
    
  enddo
  
enddo

open (unit = 2,file = fname, status='unknown',form=write_format, &
  action='write',position=write_position)
   
!  Write the data
select case(write_format)
  case('formatted')
  
    if (coord_pres) then
  
      write(2,data_format_double) x, vars_1d

   else

      write(2,data_format_double) vars_1d
  
    endif
    
  case('unformatted')
  
    if (coord_pres) then

       write(2) x, vars_1d

    else

       write(2) vars_1d

    endif
    
end select

close(2)

deallocate ( vars_1d )

return
end subroutine write_real_data_1D_double_
