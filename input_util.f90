!**********************************************************************
module input_mod
!**********************************************************************
implicit none

save 
private

public read_input_conf

character (*), parameter :: input_conf = 'lesgo.conf'
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
use param
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

  case ('DOMAIN')

     call domain_block()

  case ('MODEL')
         
     call model_block()

  case ('TIME')
     
     call time_block()

  case ('FLOW_COND')

     call flow_cond_block()      
      
  case ('OUTPUT')

     call output_block()

  case default

     call error (sub, 'invalid block label ' // buff(1:equal_pos-1) //  &
          'at line', line)

  end select
  
end do

close (lun)

! write(*,*) 'GRID SETTINGS : '
! write(*,'(a16,i9)')    'Nx : ', Nx
! write(*,'(a16,i9)')    'Ny : ', Ny
! write(*,'(a16,f12.6)')  'Lx : ', Lx
! write(*,'(a16,f12.6)') 'Ly : ', Ly
! write(*,'(a16,l)') 'non_uniform : ', non_uniform
! write(*,'(a16,f12.6)') 'alpha : ', alpha

! write(*,*) ''
! write(*,*) 'TIME SETTINGS : '
! write(*,'(a16,f12.6)') 'dt : ', dt
! write(*,'(a16,i9)') 'Nmax : ', Nmax

! write(*,*) ''
! write(*,*) 'AVERAGING SETTINGS : '
! write(*,'(a16,l)') 'avg_compute : ', avg_compute
! write(*,'(a16,i9)') 'avg_start : ', avg_start

! write(*,*) ''
! write(*,*) 'OUTPUT SETTINGS : '
! write(*,'(a16,i9)') 'output_skip : ', output_skip
! write(*,'(a16,a)') 'output_path : ', output_path
! write(*,'(a22,l)') 'output_stream_func : ', output_stream_func
! write(*,*) ''
! write(*,*) 'SOLVER SETTINGS : '
! write(*,'(a16,e12.6)') 'eps : ', eps

! write(*,*) ''
! write(*,*) 'FLOW SETTINGS : '
! write(*,'(a16,f12.6)') 'Re : ', Re

! write(*,*) ''
! write(*,*) 'BC SETTINGS : '
! write(*,'(a16,2f12.6)') 'Ue : ', Ue
! write(*,'(a16,2f12.6)') 'Uw : ', Uw
! write(*,'(a16,2f12.6)') 'Un : ', Un
! write(*,'(a16,2f12.6)') 'Us : ', Us

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine domain_block()
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
  case ('NZ') 
     read (buff(equal_pos+1:), *) Nz
  case ('Z_I')
    read (buff(equal_pos+1:), *) z_i
  case ('L_X')
    read (buff(equal_pos+1:), *) L_x
  case ('L_Y')
    read (buff(equal_pos+1:), *) L_y
  case ('L_Z')
    read (buff(equal_pos+1:), *) L_z
  case default
     
    call error (sub, 'invalid block data value ' // buff(1:equal_pos-1) //  &
                'at line', line)
  end select

enddo

return
end subroutine domain_block

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine model_block()
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

  case ('MODEL')
     read (buff(equal_pos+1:), *) model
  case ('MODELS')
     read (buff(equal_pos+1:), *) models
  case ('NNN') 
     read (buff(equal_pos+1:), *) nnn
  case ('CS_COUNT')
    read (buff(equal_pos+1:), *) cs_count
  case ('DYN_INIT')
    read (buff(equal_pos+1:), *) DYN_init
  case ('CO')
    read (buff(equal_pos+1:), *) Co
  case ('IFILTER')
    read (buff(equal_pos+1:), *) ifilter
  case ('U_STAR')
    read (buff(equal_pos+1:), *) u_star
  case ('PR')
    read (buff(equal_pos+1:), *) Pr
  case ('VONK')
    read (buff(equal_pos+1:), *) vonk
  case ('CORIOLIS_FORCING')
    read (buff(equal_pos+1:), *) coriolis_forcing
  case ('UG')
    read (buff(equal_pos+1:), *) ug
  case ('VG')
    read (buff(equal_pos+1:), *) vg
  case ('NU_MOLEC')
    read (buff(equal_pos+1:), *) nu_molec
  case ('MOLEC')
    read (buff(equal_pos+1:), *) molec
  case ('SGS')
    read (buff(equal_pos+1:), *) sgs
  case ('DNS_BC')
    read (buff(equal_pos+1:), *) dns_bc

  case default
     
    call error (sub, 'invalid block data value ' // buff(1:equal_pos-1) //  &
                'at line', line)
  end select

enddo

return
end subroutine model_block

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

  case ('NSTEPS')
    read (buff(equal_pos+1:), *) nsteps

  $if($CFL_DT)
  case ('CFL')
    read (buff(equal_pos+1:), *) cfl
  $else
  case('DT')
     read (buff(equal_pos+1:), *) dt
  $endif
  case default
     
    call error (sub, 'invalid block data value ' // buff(1:equal_pos-1) //  &
                'at line', line)
  end select

enddo

return
end subroutine  time_block

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine flow_cond_block()
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

  case ('INITU')
    read (buff(equal_pos+1:), *) initu
  case ('INILAG')
    read (buff(equal_pos+1:), *) inilag
  case ('UBC')
    read (buff(equal_pos+1:), *) ubc
  case ('LBC_MOM')
    Read (buff(equal_pos+1:), *) lbc_mom
  case ('ZO')
    read (buff(equal_pos+1:), *) zo
  case ('INFLOW')
    read (buff(equal_pos+1:), *) inflow
  case ('FRINGE_REGION_END')
    read (buff(equal_pos+1:), *) fringe_region_end
  case ('FRINGE_REGION_LEN')
    read (buff(equal_pos+1:), *) fringe_region_len
  case ('UNIFORM_INFLOW')
    read (buff(equal_pos+1:), *) uniform_inflow
  case ('INFLOW_VELOCITY')
    read (buff(equal_pos+1:), *) inflow_velocity
  case ('FORCE_TOP_BOT')
    read (buff(equal_pos+1:), *) force_top_bot
  case ('USE_MEAN_P_FORCE')
    read (buff(equal_pos+1:), *) use_mean_p_force
  case ('MEAN_P_FORCE')
    read (buff(equal_pos+1:), *) mean_p_force

  case default
     
    call error (sub, 'invalid block data value ' // buff(1:equal_pos-1) //  &
                'at line', line)
  end select

enddo

return
end subroutine  flow_cond_block

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

  case ('WBASE')
    read (buff(equal_pos+1:), *) wbase

  case ('NENERGY')
    read (buff(equal_pos+1:), *) nenergy

  case ('LAG_CFL_COUNT')
    read (buff(equal_pos+1:), *) cfl_count

  case ('TAVG_CALC')
    read (buff(equal_pos+1:), *) tavg_calc
  case ('TAVG_NSTART')
    read (buff(equal_pos+1:), *) tavg_nstart

  case ('POINT_CALC')
    read (buff(equal_pos+1:), *) point_calc
  case ('POINT_NSTART')
    read (buff(equal_pos+1:), *) point_nstart
  case ('POINT_NLOC')
    read (buff(equal_pos+1:), *) point_nloc
  case ('POINT_LOC')
     allocate( point_loc( point_nloc ) )
     read (buff(equal_pos+1:), *) point_loc

  case ('DOMAIN_CALC')
    read (buff(equal_pos+1:), *) domain_calc
  case ('DOMAIN_NSTART')
    read (buff(equal_pos+1:), *) domain_nstart

  case ('XPLANE_CALC')
    read (buff(equal_pos+1:), *) xplane_calc
  case ('XPLANE_NSTART')
    read (buff(equal_pos+1:), *) xplane_nstart
  case ('XPLANE_NLOC')
    read (buff(equal_pos+1:), *) xplane_nloc
  case ('XPLANE_LOC')
     allocate( xplane_loc( xplane_nloc ) )
     read (buff(equal_pos+1:), *) xplane_loc

  case ('YPLANE_CALC')
    read (buff(equal_pos+1:), *) yplane_calc
  case ('YPLANE_NSTART')
    read (buff(equal_pos+1:), *) yplane_nstart
  case ('YPLANE_NLOC')
    read (buff(equal_pos+1:), *) yplane_nloc
  case ('YPLANE_LOC')
     allocate( yplane_loc( yplane_nloc ) )
     read (buff(equal_pos+1:), *) yplane_loc

  case ('ZPLANE_CALC')
    read (buff(equal_pos+1:), *) zplane_calc
  case ('ZPLANE_NSTART')
    read (buff(equal_pos+1:), *) zplane_nstart
  case ('ZPLANE_NLOC')
    read (buff(equal_pos+1:), *) zplane_nloc
  case ('ZPLANE_LOC')
     allocate( zplane_loc( zplane_nloc ) )
     read (buff(equal_pos+1:), *) zplane_loc

  case ('SPECTRA_CALC')
    read (buff(equal_pos+1:), *) spectra_calc
  case ('SPECTRA_NSTART')
    read (buff(equal_pos+1:), *) spectra_nstart
  case ('SPECTRA_NEND')
    read (buff(equal_pos+1:), *) spectra_nend
  case ('SPECTRA_NLOC')
    read (buff(equal_pos+1:), *) spectra_nloc
  case ('SPECTRA_LOC')
     allocate( spectra_loc( spectra_nloc ) )
     read (buff(equal_pos+1:), *) spectra_loc

  case default
     
    call error (sub, 'invalid block data value ' // buff(1:equal_pos-1) //  &
                ' at line', line)
  end select

enddo

return
end subroutine  output_block

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

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine readline(lun, line, buff, block_entry_pos, &
                    block_exit_pos, equal_pos, ios )
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
! This subroutine reads the specified line and determines the attributes
! of the contents of the line.
!
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

  call eat_white_space (buff)

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
