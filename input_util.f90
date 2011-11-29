!**********************************************************************
module input_util
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

! Delimiters used for reading vectors and points
character(*), parameter :: delim_minor=','
character(*), parameter :: delim_major='//'

! Default buffer length for characters of unknown length
integer, parameter :: BUFF_LEN = 256

interface parse_vector
  module procedure parse_vector_real, parse_vector_point3D
end interface

contains

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine read_input_conf ()
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use param
use messages
use string_util, only : eat_whitespace, uppercase
implicit none

character (*), parameter :: sub = 'read_input_conf'

integer, parameter :: lun = 1


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

  write(*,*) buff

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

     call mesg( sub, 'Found unused input block :' // buff(1:equal_pos-1) )

  end select
  
end do

close (lun)

contains

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine domain_block()
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use param
implicit none

write(*,*) 'Im in domain_block'

do 

  call readline( lun, line, buff, block_entry_pos, block_exit_pos, &
                 equal_pos, ios )
  if (ios /= 0) call error( sub, 'Bad read in block')

  if( block_exit_pos == 0 ) then

     ! Check that the data entry conforms to correct format
     call checkentry() 
     
     select case (uppercase(buff(1:equal_pos-1)))

     case ('NX')
        read (buff(equal_pos+1:), *) Nx
     case ('NY')
        read (buff(equal_pos+1:), *) Ny
     case ('NZ') 
        read (buff(equal_pos+1:), *) Nz_tot
     case ('Z_I')
        read (buff(equal_pos+1:), *) z_i
     case ('LX')
        read (buff(equal_pos+1:), *) L_x
     case ('LY')
        read (buff(equal_pos+1:), *) L_y
     case ('LZ')
        read (buff(equal_pos+1:), *) L_z
     case default
        call mesg( sub, 'Found unused data value in GRID block :' // buff(1:equal_pos-1) )
     end select

  elseif ( block_exit_pos == 1 ) then

     ! === Set dependant variables ===
     
     ! Set the processor owned vertical grid spacing
     nz = ( nz_tot - 1 ) / nproc + 1 
     ! Grid size for dealiasing
     nx2 = 3 * nx / 2
     ny2 = 3 * ny / 2
     ! Grid size for FFT's
     lh = nx / 2 + 1
     ld = 2 * lh
     lh_big = nx2 / 2 + 1
     ld_big = 2 * lh_big
     
     ! Grid spacing
     dx = L_x / nx
     dy = L_y / ny
     dz = L_z / ( nz_tot - 1 )

     return

  else

     call error( sub, 'GRID data block not formatted correctly:' // buff(1:equal_pos-1) )

  endif
     
enddo

return
end subroutine domain_block

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine model_block()
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use param
implicit none

do 

  call readline( lun, line, buff, block_entry_pos, block_exit_pos, &
                 equal_pos, ios )
  if (ios /= 0) call error( sub, 'Bad read in block') 

  if( block_exit_pos == 0 ) then

     ! Check that the data entry conforms to correct format
     call checkentry()  

     select case (uppercase(buff(1:equal_pos-1)))

     case ('MODEL')
        read (buff(equal_pos+1:), *) sgs_model
     case ('WALL_DAMP_EXP') 
        read (buff(equal_pos+1:), *) wall_damp_exp
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
     case ('VONK')
        read (buff(equal_pos+1:), *) vonk
     case ('CORIOLIS_FORCING')
        read (buff(equal_pos+1:), *) coriolis_forcing
     case ('CORIOL')
        read (buff(equal_pos+1:), *) coriol
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
        call mesg( sub, 'Found unused data value in MODEL block :' // buff(1:equal_pos-1) )
     end select

  elseif( block_exit_pos == 1 ) then

     return
     
  else
     
     call error( sub, 'MODEL data block not formatted correctly:' // buff(1:equal_pos-1) )

  endif

enddo

return
end subroutine model_block

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine time_block()
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use param
implicit none

do 

  call readline( lun, line, buff, block_entry_pos, block_exit_pos, &
                 equal_pos, ios )
  if (ios /= 0) call error( sub, 'Bad read in block')


  if( block_exit_pos == 0 ) then

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

     case('CUMULATIVE_TIME')
        read (buff(equal_pos+1:), *) cumulative_time
     case default
        call mesg( sub, 'Found unused data value in TIME block :' // buff(1:equal_pos-1) )
     end select

  elseif( block_exit_pos == 1 ) then

     ! Set dependent data
     $if(not $CFL_DT)
     ! Set dimensional time step
     dt_dim = dt * z_i / u_star
     ! Set AB2 integration coefficients
     tadv1 = 1.5_rprec
     tadv2 = 1.0_rprec - tadv1
     $endif
     return

  else

     call error( sub, 'TIME data block not formatted correctly:' // buff(1:equal_pos-1) )

  endif

enddo

return
end subroutine  time_block

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine flow_cond_block()
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use param
implicit none

do 

  call readline( lun, line, buff, block_entry_pos, block_exit_pos, &
                 equal_pos, ios )
  if (ios /= 0) call error( sub, 'Bad read in block')

  if( block_exit_pos == 0 ) then

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

        call mesg( sub, 'Found unused data value in FLOW_COND block :' // buff(1:equal_pos-1) )

     end select

  elseif( block_exit_pos == 1 ) then

     return

  else

     call error( sub, 'FLOW_COND data block not formatted correctly:' // buff(1:equal_pos-1) )

  endif

enddo

return
end subroutine  flow_cond_block

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine output_block()
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use param
implicit none

do 

  call readline( lun, line, buff, block_entry_pos, block_exit_pos, &
                 equal_pos, ios )
  if (ios /= 0) call error( sub, 'Bad read in block')

  if( block_exit_pos == 0 ) then

     ! Check that the data entry conforms to correct format
     call checkentry()

     select case (uppercase(buff(1:equal_pos-1)))

     case ('WBASE')
        read (buff(equal_pos+1:), *) wbase

     case ('NENERGY')
        read (buff(equal_pos+1:), *) nenergy

     case ('LAG_CFL_COUNT')
        read (buff(equal_pos+1:), *) lag_cfl_count

     case ('TAVG_CALC')
        read (buff(equal_pos+1:), *) tavg_calc
     case ('TAVG_NSTART')
        read (buff(equal_pos+1:), *) tavg_nstart
     case ('TAVG_NEND')
        read (buff(equal_pos+1:), *) tavg_nend

     case ('POINT_CALC')
        read (buff(equal_pos+1:), *) point_calc
     case ('POINT_NSTART')
        read (buff(equal_pos+1:), *) point_nstart
     case ('POINT_NEND')
        read (buff(equal_pos+1:), *) point_nend
     case ('POINT_NSKIP')
        read (buff(equal_pos+1:), *) point_nskip
     case ('POINT_NLOC')
        read (buff(equal_pos+1:), *) point_nloc
     case ('POINT_LOC')
        allocate( point_loc( point_nloc ) )
        call parse_vector( buff(equal_pos+1:), point_loc )

     case ('DOMAIN_CALC')
        read (buff(equal_pos+1:), *) domain_calc
     case ('DOMAIN_NSTART')
        read (buff(equal_pos+1:), *) domain_nstart
     case ('DOMAIN_NEND')
        read (buff(equal_pos+1:), *) domain_nend
     case ('DOMAIN_NSKIP')
        read (buff(equal_pos+1:), *) domain_nskip

     case ('XPLANE_CALC')
        read (buff(equal_pos+1:), *) xplane_calc
     case ('XPLANE_NSTART')
        read (buff(equal_pos+1:), *) xplane_nstart
     case ('XPLANE_NEND')
        read (buff(equal_pos+1:), *) xplane_nend
     case ('XPLANE_NSKIP')
        read (buff(equal_pos+1:), *) xplane_nskip
     case ('XPLANE_NLOC')
        read (buff(equal_pos+1:), *) xplane_nloc
     case ('XPLANE_LOC')
        allocate( xplane_loc( xplane_nloc ) )
        call parse_vector( buff(equal_pos+1:), xplane_loc )

     case ('YPLANE_CALC')
        read (buff(equal_pos+1:), *) yplane_calc
     case ('YPLANE_NSTART')
        read (buff(equal_pos+1:), *) yplane_nstart
     case ('YPLANE_NEND')
        read (buff(equal_pos+1:), *) yplane_nend
     case ('YPLANE_NSKIP')
        read (buff(equal_pos+1:), *) yplane_nskip
     case ('YPLANE_NLOC')
        read (buff(equal_pos+1:), *) yplane_nloc
     case ('YPLANE_LOC')
        allocate( yplane_loc( yplane_nloc ) )
        call parse_vector( buff(equal_pos+1:), yplane_loc )

     case ('ZPLANE_CALC')
        read (buff(equal_pos+1:), *) zplane_calc
     case ('ZPLANE_NSTART')
        read (buff(equal_pos+1:), *) zplane_nstart
     case ('ZPLANE_NEND')
        read (buff(equal_pos+1:), *) zplane_nend
     case ('ZPLANE_NSKIP')
        read (buff(equal_pos+1:), *) zplane_nskip
     case ('ZPLANE_NLOC')
        read (buff(equal_pos+1:), *) zplane_nloc
     case ('ZPLANE_LOC')
        allocate( zplane_loc( zplane_nloc ) )
        call parse_vector( buff(equal_pos+1:), zplane_loc )

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
        call parse_vector( buff(equal_pos+1:), spectra_loc )

     case default
        call mesg( sub, 'Found unused data value in OUTPUT block :' // buff(1:equal_pos-1) )
     end select

  elseif( block_exit_pos == 1 ) then

     return

  else

     call error( sub, 'OUTPUT data block not formatted correctly:' // buff(1:equal_pos-1) )

  endif

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
use string_util, only : eat_whitespace
implicit none

integer, intent(in) :: lun
integer, intent(inout) :: line

character(*), intent(inout) :: buff
integer, intent(out) :: block_entry_pos, &
                        block_exit_pos, &
                        equal_pos, ios

block_entry_pos = 0
block_exit_pos = 0
equal_pos = 0
ios = -1

do     

  line = line + 1
  read (lun, '(a)', iostat=ios) buff
  if (ios /= 0) exit

  call eat_whitespace (buff)

  if (verify (buff, ' ') == 0) cycle  !--drop blank lines
  
  if (buff (1:len (comment)) == comment) cycle  !--drop comment lines

  block_entry_pos = index( buff, block_entry )
  block_exit_pos  = index( buff, block_exit )
  equal_pos       = index( buff, equal )

  exit

enddo 
return
end subroutine readline

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine parse_vector_real( string, vector )
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use types, only : rprec
use string_util, only : split_string

implicit none

character(*), intent(in) :: string
real(rprec), dimension(:), intent(inout) :: vector
character(BUFF_LEN), dimension(:), allocatable :: svector

integer :: nelem

! Get the number of elements in the vector
nelem = size(vector,1)
allocate( svector( nelem ) )

call split_string( string, delim_minor, nelem, svector )
read( svector, * ) vector

deallocate(svector)

return
end subroutine parse_vector_real

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine parse_vector_point3D( string, vector )
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use types, only : rprec, point3D
use string_util, only : split_string

implicit none

character(*), intent(in) :: string
type(point3D), dimension(:), intent(inout) :: vector
character(BUFF_LEN), allocatable, dimension(:) :: svector

integer :: n, nelem
real(rprec), dimension(3) :: vector_minor

! Get the number of elements in the vector
nelem = size(vector,1)

allocate( svector( nelem ) )

! Split based on major delimiter
call split_string( string, delim_major, nelem, svector )
! Now parse result string 
do n=1, nelem
   call parse_vector_real( svector(n), vector_minor )
   vector(n) = point3D( (/ vector_minor(1), vector_minor(2), vector_minor(3) /) )
enddo

deallocate(svector)

return
end subroutine parse_vector_point3D

end module input_util
