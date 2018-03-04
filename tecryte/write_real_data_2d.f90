!*************************************************************
subroutine write_real_data_2D_single_(fname, write_position, write_format, nvars, &
  imax, jmax, vars, ibuff, x, y)
!*************************************************************
! 
!  This subroutine variables the variables given by vars to the
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
!  jmax (int) - size of 2nd dimension of variables in vars 
!  vars (real, vector) - vector contaning variables to write
!  ibuff (int) - flag for adding buffer region due to periodicity
!     0 - no buffer region
!     1 - buffer on i direction
!     2 - buffer on j direction
!     3 - buffer on i and j directions
!  x,y (real, vector, optional) - vectors containing x,y coordinates 
!
use tecryte
implicit none

character(*), intent(in) :: fname, write_position, write_format
integer, intent(in) :: nvars, imax, jmax
real(4), intent(in), dimension(nvars*imax*jmax) :: vars
integer, intent(in) :: ibuff
real(4), intent(in), dimension(:), optional :: x, y

real(4), allocatable, dimension(:,:) :: x_2d, y_2d
real(4), allocatable, dimension(:,:,:) :: vars_2d

logical :: coord_pres

character(*), parameter :: sub_name = mod_name // '.write_real_data_2D'

integer :: i,j,n
integer :: i0, j0, imax_buff, jmax_buff

!  Check that position and format were specified sanely
call check_write_position(write_position, sub_name)
call check_write_format(write_format, sub_name)

!  Check if spatial coordinates specified
coord_pres=.false.
if(present(x) .and. present(y)) coord_pres = .true.

if( ibuff == 0 ) then

  imax_buff = imax
  jmax_buff = jmax
  
elseif( ibuff == 1 ) then

  imax_buff = imax + 1
  jmax_buff = jmax

elseif( ibuff == 2 ) then

  imax_buff = imax
  jmax_buff = jmax + 1
  
elseif( ibuff == 3 ) then

  imax_buff = imax + 1
  jmax_buff = jmax + 1
  
else 

  write(*,*) 'ibuff not specified correctly'
  stop
  
endif

allocate(vars_2d(imax_buff,jmax_buff,nvars))

do n=1,nvars

  do j = 1, jmax_buff
  
    j0 = buff_indx(j,jmax)
    
    do i = 1, imax_buff

      i0 = buff_indx(i,imax)
      
      vars_2d(i,j,n) = vars( (n-1)*imax*jmax + (j0-1)*imax + i0 )
    
    enddo
    
  enddo
  
enddo

if( coord_pres ) then

  allocate(x_2d(imax_buff, jmax_buff))
  allocate(y_2d(imax_buff, jmax_buff))

  do j=1, jmax_buff
    do i=1, imax_buff
      x_2d(i,j) = x(i)
      y_2d(i,j) = y(j)
    enddo
  enddo

endif

open (unit = 2,file = fname, status='unknown',form=write_format, &
  action='write',position=write_position)

!  Write the data
select case(write_format)

  case('formatted')
  
    !  Specify output format; may want to use a global setting
    if (coord_pres) then

      write(2,data_format_single) x_2d, y_2d, vars_2d

    else

      write(2,data_format_single) vars_2d

    endif

  case('unformatted')

    if (coord_pres) then

      write(2) x_2d, y_2d, vars_2d

    else

      write(2) vars_2d

    endif

  
end select

close(2)


deallocate ( vars_2d )

if ( coord_pres ) then
    deallocate ( x_2d )
    deallocate ( y_2d )
endif


return
end subroutine write_real_data_2D_single_
!*************************************************************
subroutine write_real_data_2D_double_(fname, write_position, write_format, nvars, &
  imax, jmax, vars, ibuff, x, y)
!*************************************************************
! 
!  This subroutine variables the variables given by vars to the
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
!  jmax (int) - size of 2nd dimension of variables in vars 
!  vars (real, vector) - vector contaning variables to write
!  ibuff (int) - flag for adding buffer region due to periodicity
!     0 - no buffer region
!     1 - buffer on i direction
!     2 - buffer on j direction
!     3 - buffer on i and j directions
!  x,y (real, vector, optional) - vectors containing x,y coordinates 
!
use tecryte
implicit none

character(*), intent(in) :: fname, write_position, write_format
integer, intent(in) :: nvars, imax, jmax
real(8), intent(in), dimension(nvars*imax*jmax) :: vars
integer, intent(in) :: ibuff
real(8), intent(in), dimension(:), optional :: x, y

real(8), allocatable, dimension(:,:) :: x_2d, y_2d
real(8), allocatable, dimension(:,:,:) :: vars_2d

logical :: coord_pres

character(*), parameter :: sub_name = mod_name // '.write_real_data_2D'

integer :: i,j,n
integer :: i0, j0, imax_buff, jmax_buff

!  Check that position and format were specified sanely
call check_write_position(write_position, sub_name)
call check_write_format(write_format, sub_name)

!  Check if spatial coordinates specified
coord_pres=.false.
if(present(x) .and. present(y)) coord_pres = .true.

if( ibuff == 0 ) then

  imax_buff = imax
  jmax_buff = jmax
  
elseif( ibuff == 1 ) then

  imax_buff = imax + 1
  jmax_buff = jmax

elseif( ibuff == 2 ) then

  imax_buff = imax
  jmax_buff = jmax + 1
  
elseif( ibuff == 3 ) then

  imax_buff = imax + 1
  jmax_buff = jmax + 1
  
else 

  write(*,*) 'ibuff not specified correctly'
  stop
  
endif

allocate(vars_2d(imax_buff,jmax_buff,nvars))

do n=1,nvars

  do j = 1, jmax_buff
  
    j0 = buff_indx(j,jmax)
    
    do i = 1, imax_buff

      i0 = buff_indx(i,imax)
      
      vars_2d(i,j,n) = vars( (n-1)*imax*jmax + (j0-1)*imax + i0 )
    
    enddo
    
  enddo
  
enddo

if( coord_pres ) then

  allocate(x_2d(imax_buff, jmax_buff))
  allocate(y_2d(imax_buff, jmax_buff))

  do j=1, jmax_buff
    do i=1, imax_buff
      x_2d(i,j) = x(i)
      y_2d(i,j) = y(j)
    enddo
  enddo

endif

open (unit = 2,file = fname, status='unknown',form=write_format, &
  action='write',position=write_position)

!  Write the data
select case(write_format)

  case('formatted')
  
    !  Specify output format; may want to use a global setting
    if (coord_pres) then

      write(2,data_format_double) x_2d, y_2d, vars_2d

    else

      write(2,data_format_double) vars_2d

    endif

  case('unformatted')

    if (coord_pres) then

      write(2) x_2d, y_2d, vars_2d

    else

      write(2) vars_2d

    endif

  
end select

close(2)


deallocate ( vars_2d )

if ( coord_pres ) then
    deallocate ( x_2d )
    deallocate ( y_2d )
endif


return
end subroutine write_real_data_2D_double_
