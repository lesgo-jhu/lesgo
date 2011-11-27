!**********************************************************************
module input_mod
!**********************************************************************
implicit none

save 
private

public read_input_conf

character (*), parameter :: input_conf = 'input.conf'
character (*), parameter :: comment = '!'
!character (*), parameter :: ldelim = '('  !--no whitespace allowed
!character (*), parameter :: rdelim = ')'  !--no whitespace allowed
character (*), parameter :: block_entry = '{'
character (*), parameter :: block_exit = '}'
character (*), parameter :: equal = '='
character (*), parameter :: esyntax = 'syntax error at line'

contains


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine read_input_conf ()
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use settings
use messages
use strmod, only : eat_whtspc, uppercase
implicit none

character (*), parameter :: sub = 'read_conf'

integer, parameter :: lun = 1
integer, parameter :: BUFF_LEN = 256

character (BUFF_LEN) :: buff

integer :: block_entry_pos, block_exit_pos, equal_pos

integer :: ios
integer :: line

logical :: exst

! Check that the configuration file exists
inquire (file=input_conf, exist=exst)

if (exst) then
  open (lun, file=input_conf, action='read')
else
  call error (sub, 'file ' // input_conf // ' does not exist')
end if

line = 0
do

  call readline( lun, line, buff, block_entry_pos, block_exit_pos, &
                 equal_pos, ios )
  if (ios /= 0) exit

  if (block_entry_pos == 0) then  !--for now, invalid format if no block entry found
    call error (sub, 'block entry not found on line', line) 
  end if

  ! Find block
  select case (uppercase(buff(1:block_entry_pos-1)))

    case ('GRID')

      call grid_block()

    case ('TIME')

      call time_block()

    case ('AVERAGING')

      call averaging_block()      
      
    case ('OUTPUT')

      call output_block()

    case ('SOLVER')

      call solver_block()      

    case ('FLOW')

      call flow_block()      

    case ('BC')

      call bc_block()
      
    case default
      call error (sub, 'invalid block label ' // buff(1:equal_pos-1) //  &
                  'at line', line)
  end select
  
end do

close (lun)

write(*,*) 'GRID SETTINGS : '
write(*,'(a16,i9)')    'Nx : ', Nx
write(*,'(a16,i9)')    'Ny : ', Ny
write(*,'(a16,f12.6)')  'Lx : ', Lx
write(*,'(a16,f12.6)') 'Ly : ', Ly
write(*,'(a16,l)') 'non_uniform : ', non_uniform
write(*,'(a16,f12.6)') 'alpha : ', alpha

write(*,*) ''
write(*,*) 'TIME SETTINGS : '
write(*,'(a16,f12.6)') 'dt : ', dt
write(*,'(a16,i9)') 'Nmax : ', Nmax

write(*,*) ''
write(*,*) 'AVERAGING SETTINGS : '
write(*,'(a16,l)') 'avg_compute : ', avg_compute
write(*,'(a16,i9)') 'avg_start : ', avg_start

write(*,*) ''
write(*,*) 'OUTPUT SETTINGS : '
write(*,'(a16,i9)') 'output_skip : ', output_skip
write(*,'(a16,a)') 'output_path : ', output_path
write(*,'(a22,l)') 'output_stream_func : ', output_stream_func
write(*,*) ''
write(*,*) 'SOLVER SETTINGS : '
write(*,'(a16,e12.6)') 'eps : ', eps

write(*,*) ''
write(*,*) 'FLOW SETTINGS : '
write(*,'(a16,f12.6)') 'Re : ', Re

write(*,*) ''
write(*,*) 'BC SETTINGS : '
write(*,'(a16,2f12.6)') 'Ue : ', Ue
write(*,'(a16,2f12.6)') 'Uw : ', Uw
write(*,'(a16,2f12.6)') 'Un : ', Un
write(*,'(a16,2f12.6)') 'Us : ', Us

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine grid_block()
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use settings
implicit none

do 

  call readline( lun, line, buff, block_entry_pos, block_exit_pos, &
                 equal_pos, ios )
  if (ios /= 0) call error( sub, 'Bad read in block')

  ! Check if we have found a block exit
  if( block_exit_pos == 1 ) return 

  ! Check that the data entry conforms to correct format
  call checkentry()  

  select case (uppercase(buff(1:equal_pos-1)))

  case ('NX')
    read (buff(equal_pos+1:), *) Nx
  case ('NY')
    read (buff(equal_pos+1:), *) Ny
  case ('LX')
    read (buff(equal_pos+1:), *) Lx
  case ('LY')
    read (buff(equal_pos+1:), *) Ly
  case ('NON_UNIFORM')
    read (buff(equal_pos+1:), *) non_uniform
  case ('ALPHA')
    read (buff(equal_pos+1:), *) alpha
  case default
     
    call error (sub, 'invalid grid block data value ' // buff(1:equal_pos-1) //  &
                'at line', line)
  end select

enddo

return
end subroutine grid_block

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine time_block()
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use settings
implicit none

do 

  call readline( lun, line, buff, block_entry_pos, block_exit_pos, &
                 equal_pos, ios )
  if (ios /= 0) call error( sub, 'Bad read in block')


  ! Check if we have found a block exit
  if( block_exit_pos == 1 ) return 

  ! Check that the data entry conforms to correct format
  call checkentry()

  select case (uppercase(buff(1:equal_pos-1)))

  case ('DT')
    read (buff(equal_pos+1:), *) dt
  case ('NMAX')
    read (buff(equal_pos+1:), *) Nmax
  case default
     
    call error (sub, 'invalid block data value ' // buff(1:equal_pos-1) //  &
                'at line', line)
  end select

enddo

return
end subroutine  time_block

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine averaging_block()
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use settings
implicit none

do 

  call readline( lun, line, buff, block_entry_pos, block_exit_pos, &
                 equal_pos, ios )
  if (ios /= 0) call error( sub, 'Bad read in block')


  ! Check if we have found a block exit
  if( block_exit_pos == 1 ) return 

  ! Check that the data entry conforms to correct format
  call checkentry()

  select case (uppercase(buff(1:equal_pos-1)))

  case ('AVG_COMPUTE')
    read (buff(equal_pos+1:), *) avg_compute
  case ('AVG_START')
    read (buff(equal_pos+1:), *) avg_start
  case default
     
    call error (sub, 'invalid block data value ' // buff(1:equal_pos-1) //  &
                'at line', line)
  end select

enddo

return
end subroutine  averaging_block

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine output_block()
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use settings
implicit none

do 

  call readline( lun, line, buff, block_entry_pos, block_exit_pos, &
                 equal_pos, ios )
  if (ios /= 0) call error( sub, 'Bad read in block')


  if( block_exit_pos == 1 ) return ! Exit block '}' found

    ! Check that the data entry conforms to correct format
  call checkentry()

  select case (uppercase(buff(1:equal_pos-1)))

  case ('OUTPUT_SKIP')
    read (buff(equal_pos+1:), *) output_skip
  case ('OUTPUT_PATH')
    read (buff(equal_pos+1:), *) output_path
  case ('OUTPUT_STREAM_FUNC')
    read (buff(equal_pos+1:), *) output_stream_func
  case default
     
    call error (sub, 'invalid block data value ' // buff(1:equal_pos-1) //  &
                ' at line', line)
  end select

enddo

return
end subroutine  output_block

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine solver_block()
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use settings
implicit none

do 

  call readline( lun, line, buff, block_entry_pos, block_exit_pos, &
                 equal_pos, ios )
  if (ios /= 0) call error( sub, 'Bad read in block')


  ! Check if we have found a block exit
  if( block_exit_pos == 1 ) return

  ! Check that the data entry conforms to correct format
  call checkentry()

  select case (uppercase(buff(1:equal_pos-1)))

  case ('EPS')
    
    read (buff(equal_pos+1:), *) eps

  case default
     
    call error (sub, 'invalid block data value ' // buff(1:equal_pos-1) //  &
                ' at line', line)
  end select

enddo

return
end subroutine  solver_block


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine flow_block()
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use settings
implicit none

do 

  call readline( lun, line, buff, block_entry_pos, block_exit_pos, &
                 equal_pos, ios )
  if (ios /= 0) call error( sub, 'Bad read in block')


  ! Check if we have found a block exit
  if( block_exit_pos == 1 ) return

  ! Check that the data entry conforms to correct format
  call checkentry()

  select case (uppercase(buff(1:equal_pos-1)))

  case ('RE')
    
    read (buff(equal_pos+1:), *) Re

  case default
     
    call error (sub, 'invalid block data value ' // buff(1:equal_pos-1) //  &
                ' at line', line)
  end select

enddo

return
end subroutine  flow_block

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine bc_block()
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use settings
implicit none

do 

  call readline( lun, line, buff, block_entry_pos, block_exit_pos, &
                 equal_pos, ios )
  if (ios /= 0) call error( sub, 'Bad read in block')


  ! Check if we have found a block exit
  if( block_exit_pos == 1 ) return

  ! Check that the data entry conforms to correct format
  call checkentry()

  select case (uppercase(buff(1:equal_pos-1)))

  case ('UE')
    
    read (buff(equal_pos+1:), *) Ue

  case ('UW')
    
    read (buff(equal_pos+1:), *) Uw

  case ('UN')
    
    read (buff(equal_pos+1:), *) Un

  case ('US')
    
    read (buff(equal_pos+1:), *) Us    

  case default
     
    call error (sub, 'invalid block data value ' // buff(1:equal_pos-1) //  &
                ' at line', line)
  end select

enddo

return
end subroutine  bc_block

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine checkentry()
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

if( equal_pos == 0 ) call error( sub, 'Bad read in block at line', line)
!--invalid if nothing after equals
if (len_trim (buff) == equal_pos) call error (sub, 'nothing after equals sign in line', line) 

return
end subroutine checkentry  

end subroutine read_input_conf

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine readline(lun, line, buff, block_entry_pos, &
                    block_exit_pos, equal_pos, ios )
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use strmod
implicit none

integer, intent(in) :: lun
integer, intent(inout) :: line

character(*), intent(inout) :: buff
integer, intent(out) :: block_entry_pos, block_exit_pos, &
                        equal_pos, ios

block_entry_pos = 0
block_exit_pos = 0
equal_pos = 0
ios = -1

do     

  line = line + 1
  read (lun, '(a)', iostat=ios) buff
  if (ios /= 0) exit

  call eat_whtspc (buff)

  if (verify (buff, ' ') == 0) cycle  !--drop blank lines
  
  if (buff (1:len (comment)) == comment) cycle  !--drop comment lines

  block_entry_pos = index( buff, block_entry )
  block_exit_pos  = index( buff, block_exit )
  equal_pos       = index( buff, equal )

  exit

enddo 
return
end subroutine readline

end module input_mod
