!!
!!  Copyright 2011,2012,2013 Johns Hopkins University
!!
!!  Licensed under the Apache License, Version 2.0 (the "License"); you may not 
!!  use this file except in compliance with the License. You may obtain a copy of
!!  the License at:
!!
!!    http://www.apache.org/licenses/LICENSE-2.0
!!
!!  Unless required by applicable law or agreed to in writing, software 
!!  distributed under the License is distributed on an "AS IS" BASIS, WITHOUT 
!!  WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the 
!!  License for the specific language governing permissions and limitations under
!!  the License.
!!

!**********************************************************************
module input_util
!**********************************************************************
use types, only : rprec
use param, only : path, CHAR_BUFF_LENGTH
implicit none

save 
private

public :: read_input_conf

character (*), parameter :: mod_name = 'string_util'
 
character (*), parameter :: input_conf = path // 'lesgo.conf'
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

! Thresh hold for evaluating differences in floating point values.
real(rprec), parameter :: thresh = 1.0e-6_rprec

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

character (*), parameter :: sub_name = mod_name // '.read_input_conf'

integer, parameter :: lun = 1


character (CHAR_BUFF_LENGTH) :: buff

integer :: block_entry_pos, block_exit_pos, equal_pos

integer :: ios
integer :: line

logical :: exst

! Check that the configuration file exists
inquire (file=input_conf, exist=exst)

if (exst) then
  open (lun, file=input_conf, action='read')
else
  call error (sub_name, 'file ' // input_conf // ' does not exist')
end if

line = 0
do

  call readline( lun, line, buff, block_entry_pos, block_exit_pos, &
                 equal_pos, ios )

  if (ios /= 0) exit

  if (block_entry_pos == 0) then  !--for now, invalid format if no block entry found
    call error (sub_name, 'block entry not found on line', line) 
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

  $if($LVLSET)
  case ('LEVEL_SET')
     
     call level_set_block()

    $if($RNS_LS)
    case ('RNS_LS')
       call rns_block()
    $endif

    $if($CYL_SKEW_LS)
    case ('CYL_SKEW_LS')
       call cyl_skew_block()
    $endif

  $endif

  case ('SGS_HIST')

     call sgs_hist_block()

  $if ($TURBINES)
  case ('TURBINES')
  
     call turbines_block()
  $endif

  case default

     if(coord == 0) call mesg( sub_name, 'Found unused input block: ' // buff(1:block_entry_pos-1) )
     ! Now need to 'fast-forward' untile we reach the end of the block
     do while ( block_exit_pos == 0 )
        call readline( lun, line, buff, block_entry_pos, block_exit_pos, &
             equal_pos, ios )
        if (ios /= 0) exit ! exit if end of file is reached
     enddo

  end select
  
end do

close (lun)

contains

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine domain_block()
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use types, only : rprec
use param
implicit none

character(*), parameter :: block_name = 'DOMAIN'

integer :: ival_read
real(rprec) :: val_read

do 

  call readline( lun, line, buff, block_entry_pos, block_exit_pos, &
                 equal_pos, ios )
  if (ios /= 0) call error( sub_name, 'Bad read in block')

  if( block_exit_pos == 0 ) then

     ! Check that the data entry conforms to correct format
     call checkentry() 
     
     select case (uppercase(buff(1:equal_pos-1)))

     case ('NPROC')
        read (buff(equal_pos+1:), *) nproc
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
     case ('UNIFORM_SPACING')
        read (buff(equal_pos+1:), *) uniform_spacing
     case default
        if(coord == 0) call mesg( sub_name, 'Found unused data value in ' // block_name // ' block: ' // buff(1:equal_pos-1) )
     end select

  elseif ( block_exit_pos == 1 ) then

     ! === Set dependant variables ===
     $if( not $MPI )
     if( nproc /= 1 ) then
        ival_read = nproc
        ! Reset to 1
        nproc=1
        if( coord == 0 ) &
             call mesg( sub_name, 'Reseting nproc to: ', nproc )          
     endif
     $endif
     
     ! Set the processor owned vertical grid spacing
     nz = floor ( real( nz_tot, rprec ) / nproc ) + 1

     ! Recompute nz_tot to be compliant with computed nz
     ival_read = nz_tot
     nz_tot = ( nz - 1 ) * nproc + 1 
     if( coord == 0 .AND. ival_read /= nz_tot ) &
          call mesg( sub_name, 'Reseting Nz (total) to: ', nz_tot )          
     ! Grid size for dealiasing
     nx2 = 3 * nx / 2
     ny2 = 3 * ny / 2
     ! Grid size for FFT's
     lh = nx / 2 + 1
     ld = 2 * lh
     lh_big = nx2 / 2 + 1
     ld_big = 2 * lh_big

     ! Grid spacing (x direction)
     dx = L_x / nx

     ! Check if we are to enforce uniform grid spacing
     if( uniform_spacing ) then

        ! Adjust L_y
        val_read = L_y
        L_y = ny * dx
        if( coord == 0 .AND. abs( val_read - L_y ) >= thresh ) &
             call mesg( sub_name, 'Reseting Ly to: ', L_y )

        ! Adjust L_z
        val_read = L_z
        L_z = (nz_tot - 1 ) * dx
        if( coord == 0 .AND. abs( val_read - L_z ) >= thresh ) &
             call mesg( sub_name, 'Reseting Lz to: ', L_z )

     endif
     
     ! Grid spacing (y and z directions)
     dy = L_y / ny
     dz = L_z / ( nz_tot - 1 )

     return

  else

     call error( sub_name, block_name // ' data block not formatted correctly: ' // buff(1:equal_pos-1) )

  endif
     
enddo

return
end subroutine domain_block

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine model_block()
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use param
implicit none

character(*), parameter :: block_name = 'MODEL'

do 

  call readline( lun, line, buff, block_entry_pos, block_exit_pos, &
                 equal_pos, ios )
  if (ios /= 0) call error( sub_name, 'Bad read in block') 

  if( block_exit_pos == 0 ) then

     ! Check that the data entry conforms to correct format
     call checkentry()  

     select case (uppercase(buff(1:equal_pos-1)))

     case ('SGS_MODEL')
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

        if(coord == 0) call mesg( sub_name, 'Found unused data value in ' // block_name // ' block: ' // buff(1:equal_pos-1) )

     end select

  elseif( block_exit_pos == 1 ) then

     return
     
  else

     call error( sub_name, block_name // ' data block not formatted correctly: ' // buff(1:equal_pos-1) )

  endif

enddo

return
end subroutine model_block

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine time_block()
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use types, only : rprec
use param
implicit none

character(*), parameter :: block_name = 'TIME'

do 

  call readline( lun, line, buff, block_entry_pos, block_exit_pos, &
                 equal_pos, ios )
  if (ios /= 0) call error( sub_name, 'Bad read in block')


  if( block_exit_pos == 0 ) then

     ! Check that the data entry conforms to correct format
     call checkentry()

     select case (uppercase(buff(1:equal_pos-1)))

     case ('NSTEPS')
        read (buff(equal_pos+1:), *) nsteps

     case ('RUNTIME')
        read (buff(equal_pos+1:), *) runtime

     case ('USE_CFL_DT') 
        read (buff(equal_pos+1:), *) use_cfl_dt
        
     case ('CFL')
        read (buff(equal_pos+1:), *) cfl

     case('DT')
        read (buff(equal_pos+1:), *) dt

     case('CUMULATIVE_TIME')
        read (buff(equal_pos+1:), *) cumulative_time
     case default

        if(coord == 0) call mesg( sub_name, 'Found unused data value in ' // block_name // ' block: ' // buff(1:equal_pos-1) )
     end select

  elseif( block_exit_pos == 1 ) then

     ! Set dependent data
     if( .not. use_cfl_dt ) then
       ! Set dimensional time step
       dt_dim = dt * z_i / u_star
       ! Set AB2 integration coefficients
       tadv1 = 1.5_rprec
       tadv2 = 1.0_rprec - tadv1
     endif

     return

  else

     call error( sub_name, block_name // ' data block not formatted correctly: ' // buff(1:equal_pos-1) )

  endif

enddo

return
end subroutine  time_block

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine flow_cond_block()
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use param
implicit none

character(*), parameter :: block_name = 'FLOW_COND'

real(rprec) :: val_read 

do 

  call readline( lun, line, buff, block_entry_pos, block_exit_pos, &
                 equal_pos, ios )
  if (ios /= 0) call error( sub_name, 'Bad read in block')

  if( block_exit_pos == 0 ) then

     ! Check that the data entry conforms to correct format
     call checkentry()

     select case (uppercase(buff(1:equal_pos-1)))

     case ('INITU')
        read (buff(equal_pos+1:), *) initu
     case ('INILAG')
        read (buff(equal_pos+1:), *) inilag
     case ('UBC')
        Read (buff(equal_pos+1:), *) lbc_mom
     case ('ZO')
        read (buff(equal_pos+1:), *) zo
     case ('INFLOW')
        read (buff(equal_pos+1:), *) inflow
     case ('FRINGE_REGION_END')
        read (buff(equal_pos+1:), *) fringe_region_end
     case ('FRINGE_REGION_LEN')
        read (buff(equal_pos+1:), *) fringe_region_len
     case ('INFLOW_VELOCITY')
        read (buff(equal_pos+1:), *) inflow_velocity
     case ('FORCE_TOP_BOT')
        read (buff(equal_pos+1:), *) force_top_bot
     case ('USE_MEAN_P_FORCE')
        read (buff(equal_pos+1:), *) use_mean_p_force
     case ('EVAL_MEAN_P_FORCE')
           read (buff(equal_pos+1:), *) eval_mean_p_force
     case ('MEAN_P_FORCE')
        read (buff(equal_pos+1:), *) mean_p_force

     case default      

        if(coord == 0) call mesg( sub_name, 'Found unused data value in ' // block_name // ' block: ' // buff(1:equal_pos-1) )

     end select

  elseif( block_exit_pos == 1 ) then

     if( use_mean_p_force .AND. eval_mean_p_force ) then

        val_read = mean_p_force
        ! Evaluate the mean pressure force
        mean_p_force = 1.0_rprec / L_z
        if( coord == 0 .AND. abs( val_read - mean_p_force ) >= thresh )  &
             call mesg( sub_name, 'Reseting mean_p_force to: ', mean_p_force ) 

     endif


     return

  else

     call error( sub_name, block_name // ' data block not formatted correctly: ' // buff(1:equal_pos-1) )

  endif

enddo

return
end subroutine  flow_cond_block

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine output_block()
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use param
implicit none

character(*), parameter :: block_name = 'OUTPUT'

do 

  call readline( lun, line, buff, block_entry_pos, block_exit_pos, &
                 equal_pos, ios )
  if (ios /= 0) call error( sub_name, 'Bad read in block')

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

     case ('CHECKPOINT_DATA')
        read (buff(equal_pos+1:), *) checkpoint_data
     case ('CHECKPOINT_NSKIP')        
        read (buff(equal_pos+1:), *) checkpoint_nskip

     case ('TAVG_CALC')
        read (buff(equal_pos+1:), *) tavg_calc
     case ('TAVG_NSTART')
        read (buff(equal_pos+1:), *) tavg_nstart
     case ('TAVG_NEND')
        read (buff(equal_pos+1:), *) tavg_nend
     case ('TAVG_NSKIP')
        read (buff(equal_pos+1:), *) tavg_nskip

     case ('POINT_CALC')
        read (buff(equal_pos+1:), *) point_calc
     ! Only read if point data is to be recorded. It is important that the
     ! following only be used elsewhere in the code if point_calc=.true.; also
     ! it is required that point_calc be listed before the rest in lesgo.conf
     case ('POINT_NSTART')
        if( point_calc ) read (buff(equal_pos+1:), *) point_nstart
     case ('POINT_NEND')
        if( point_calc ) read (buff(equal_pos+1:), *) point_nend
     case ('POINT_NSKIP')
        if( point_calc ) read (buff(equal_pos+1:), *) point_nskip
     case ('POINT_LOC')
           call parse_vector( buff(equal_pos+1:), point_nloc, point_loc )
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
     case ('XPLANE_LOC')
        call parse_vector( buff(equal_pos+1:), xplane_nloc, xplane_loc )

     case ('YPLANE_CALC')
        read (buff(equal_pos+1:), *) yplane_calc
     case ('YPLANE_NSTART')
        read (buff(equal_pos+1:), *) yplane_nstart
     case ('YPLANE_NEND')
        read (buff(equal_pos+1:), *) yplane_nend
     case ('YPLANE_NSKIP')
        read (buff(equal_pos+1:), *) yplane_nskip
     case ('YPLANE_LOC')
        call parse_vector( buff(equal_pos+1:), yplane_nloc, yplane_loc )

     case ('ZPLANE_CALC')
        read (buff(equal_pos+1:), *) zplane_calc
     case ('ZPLANE_NSTART')
        read (buff(equal_pos+1:), *) zplane_nstart
     case ('ZPLANE_NEND')
        read (buff(equal_pos+1:), *) zplane_nend
     case ('ZPLANE_NSKIP')
        read (buff(equal_pos+1:), *) zplane_nskip
     case ('ZPLANE_LOC')
        call parse_vector( buff(equal_pos+1:), zplane_nloc, zplane_loc )

     case ('SPECTRA_CALC')
        read (buff(equal_pos+1:), *) spectra_calc
     case ('SPECTRA_NSTART')
        read (buff(equal_pos+1:), *) spectra_nstart
     case ('SPECTRA_NEND')
        read (buff(equal_pos+1:), *) spectra_nend
     case ('SPECTRA_NSKIP')
        read (buff(equal_pos+1:), *) spectra_nskip
     case ('SPECTRA_LOC')
        call parse_vector( buff(equal_pos+1:), spectra_nloc, spectra_loc )

     case default

        if(coord == 0) call mesg( sub_name, 'Found unused data value in ' // block_name // ' block: ' // buff(1:equal_pos-1) )
     end select

  elseif( block_exit_pos == 1 ) then

     return

  else

     call error( sub_name, block_name // ' data block not formatted correctly: ' // buff(1:equal_pos-1) )

  endif

enddo

return
end subroutine  output_block

$if($LVLSET)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine level_set_block()
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use level_set_base
implicit none

character(*), parameter :: block_name = 'LEVEL_SET'
do 

  call readline( lun, line, buff, block_entry_pos, block_exit_pos, &
                 equal_pos, ios )
  if (ios /= 0) call error( sub_name, 'Bad read in block')

  if( block_exit_pos == 0 ) then

     ! Check that the data entry conforms to correct format
     call checkentry()

     select case (uppercase(buff(1:equal_pos-1)))

     case ('GLOBAL_CD_CALC') 
        read (buff(equal_pos+1:), *) global_CD_calc
     case ('LDIR')
        read (buff(equal_pos+1:), *) Ldir
     case ('VEL_BC')
        read (buff(equal_pos+1:), *) vel_bc
     case ('USE_LOG_PROFILE')
        Read (buff(equal_pos+1:), *) use_log_profile
     case ('USE_ENFORCE_UN')
        read (buff(equal_pos+1:), *) use_enforce_un
     case ('PHYSBC')
        read (buff(equal_pos+1:), *) physBC
     case ('USE_SMOOTH_TAU')
        read (buff(equal_pos+1:), *) use_smooth_tau
     case ('USE_EXTRAP_TAU_LOG')
        read (buff(equal_pos+1:), *) use_extrap_tau_log
     case ('USE_EXTRAP_TAU_SIMPLE')
        read (buff(equal_pos+1:), *) use_extrap_tau_simple
     case ('USE_MODIFY_DUTDN')
        read (buff(equal_pos+1:), *) use_modify_dutdn
     case ('LAG_DYN_MODIFY_BETA')
        read (buff(equal_pos+1:), *) lag_dyn_modify_beta
     case ('SMOOTH_MODE')
        read (buff(equal_pos+1:), *) smooth_mode
     case ('ZO_LEVEL_SET')
        read (buff(equal_pos+1:), *) zo_level_set
     
     $if($MPI)
     case ('NPHITOP')
        read (buff(equal_pos+1:), *) nphitop
     case ('NPHIBOT')
        read (buff(equal_pos+1:), *) nphibot
     case ('NVELTOP')
        read (buff(equal_pos+1:), *) nveltop
     case ('NVELBOT')
        read (buff(equal_pos+1:), *) nvelbot
     case ('NTAUTOP')
        read (buff(equal_pos+1:), *) ntautop
     case ('NTAUBOT')
        read (buff(equal_pos+1:), *) ntaubot
     case ('NFMMTOP')
        read (buff(equal_pos+1:), *) nFMMtop
     case ('NFMMBOT')
        read (buff(equal_pos+1:), *) nFMMbot
     $endif

     case default

        if(coord == 0) call mesg( sub_name, 'Found unused data value in ' // block_name // ' block: ' // buff(1:equal_pos-1) )

     end select

  elseif( block_exit_pos == 1 ) then

     return

  else

     call error( sub_name, block_name // ' data block not formatted correctly: ' // buff(1:equal_pos-1) )

  endif

enddo

return
end subroutine  level_set_block

$if($RNS_LS)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine rns_block()
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use rns_base_ls
implicit none

character(*), parameter :: block_name = 'RNS'

do 

  call readline( lun, line, buff, block_entry_pos, block_exit_pos, &
                 equal_pos, ios )
  if (ios /= 0) call error( sub_name, 'Bad read in block')

  if( block_exit_pos == 0 ) then

     ! Check that the data entry conforms to correct format
     call checkentry()

     select case (uppercase(buff(1:equal_pos-1)))

     case ('RNS_NTREE') 
        read (buff(equal_pos+1:), *) rns_ntree
     case ('RNS_TREE_LAYOUT')
        read (buff(equal_pos+1:), *) rns_tree_layout
     case ('TEMPORAL_WEIGHT')
        read (buff(equal_pos+1:), *) temporal_weight
     case ('TCONST')
        Read (buff(equal_pos+1:), *) Tconst
     case ('WEIGHT_NSTART')
        read (buff(equal_pos+1:), *) weight_nstart
     case ('TEMPORAL_MODEL')
        read (buff(equal_pos+1:), *) temporal_model
     case ('SPATIAL_MODEL')
        read (buff(equal_pos+1:), *) spatial_model
     case ('OUTPUT_NSKIP')
        read (buff(equal_pos+1:), *) output_nskip
     case ('CD_RAMP_NSTEP')
        read (buff(equal_pos+1:), *) CD_ramp_nstep
     case ('ALPHA_WIDTH')
        read (buff(equal_pos+1:), *) alpha_width
     case ('ALPHA_DIST')
        read (buff(equal_pos+1:), *) alpha_dist
     case ('CHI_CUTOFF')
        read (buff(equal_pos+1:), *) chi_cutoff

     case default

        if(coord == 0) call mesg( sub_name, 'Found unused data value in ' // block_name // ' block: ' // buff(1:equal_pos-1) )

     end select

  elseif( block_exit_pos == 1 ) then

     return

  else

     call error( sub_name, block_name // ' data block not formatted correctly: ' // buff(1:equal_pos-1) )

  endif

enddo

return
end subroutine  rns_block
$endif

$if($CYL_SKEW_LS)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine cyl_skew_block()
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use types, only : rprec
use param, only : pi
use cyl_skew_base_ls, only : zrot_angle, skew_angle, &
                             use_bottom_surf, z_bottom_surf, &
                             use_top_surf, z_top_surf, &
                             use_left_surf, y_left_surf, &
                             use_right_surf, y_right_surf, &
                             ntree, tree_location, &
                             ngen, ngen_reslv, nbranch, d, l, offset, &
                             scale_fact, filter_chi, filt_width
implicit none

character(*), parameter :: block_name = 'CYL_SKEW'

do 

  call readline( lun, line, buff, block_entry_pos, block_exit_pos, &
                 equal_pos, ios )

  if (ios /= 0) call error( sub_name, 'Bad read in block')

  if( block_exit_pos == 0 ) then

     ! Check that the data entry conforms to correct format
     call checkentry()

     select case (uppercase(buff(1:equal_pos-1)))

     case ('ZROT_ANGLE') 
        read (buff(equal_pos+1:), *) zrot_angle
        ! Convert to radians
        zrot_angle = pi * zrot_angle / 180.0_rprec
     case ('SKEW_ANGLE')
        read (buff(equal_pos+1:), *) skew_angle
        ! Convert to radians
        skew_angle = pi * skew_angle / 180.0_rprec
     case ('USE_BOTTOM_SURF')
        read (buff(equal_pos+1:), *) use_bottom_surf
     case ('Z_BOTTOM_SURF')
        Read (buff(equal_pos+1:), *) z_bottom_surf
     case ('USE_TOP_SURF')
        read (buff(equal_pos+1:), *) use_top_surf
     case ('Z_TOP_SURF')
        Read (buff(equal_pos+1:), *) z_top_surf
     case ('USE_RIGHT_SURF')
        read (buff(equal_pos+1:), *) use_right_surf
     case ('Y_RIGHT_SURF')
        Read (buff(equal_pos+1:), *) y_right_surf
     case ('USE_LEFT_SURF')
        read (buff(equal_pos+1:), *) use_left_surf
     case ('Y_LEFT_SURF')
        Read (buff(equal_pos+1:), *) y_left_surf
     case ('TREE_LOCATION')
        call parse_vector(buff(equal_pos+1:), ntree, tree_location )
     case ('NGEN')
        read (buff(equal_pos+1:), *) ngen
     case ('NGEN_RESLV')
        read (buff(equal_pos+1:), *) ngen_reslv
     case ('NBRANCH')
        read (buff(equal_pos+1:), *) nbranch
     case ('D')
        read (buff(equal_pos+1:), *) d
     case ('L')
        read (buff(equal_pos+1:), *) l
     case ('OFFSET')
        read (buff(equal_pos+1:), *) offset
     case ('SCALE_FACT')
        read (buff(equal_pos+1:), *) scale_fact
     case ('FILTER_CHI')
        read (buff(equal_pos+1:), *) filter_chi
     case ('FILT_WIDTH')
        read (buff(equal_pos+1:), *) filt_width

     case default

        if(coord == 0) call mesg( sub_name, 'Found unused data value in ' // block_name // ' block: ' // buff(1:equal_pos-1) )

     end select

  elseif( block_exit_pos == 1 ) then

     return

  else

     call error( sub_name, block_name // ' data block not formatted correctly: ' // buff(1:equal_pos-1) )

  endif

enddo

return
end subroutine  cyl_skew_block
$endif

$endif

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine sgs_hist_block()
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use param
implicit none

character(*), parameter :: block_name = 'SGS_HIST'

do 

  call readline( lun, line, buff, block_entry_pos, block_exit_pos, &
                 equal_pos, ios )
  if (ios /= 0) call error( sub_name, 'Bad read in block')

  if( block_exit_pos == 0 ) then

     ! Check that the data entry conforms to correct format
     call checkentry()

     select case (uppercase(buff(1:equal_pos-1)))

     case ('SGS_HIST_CALC')
        read (buff(equal_pos+1:), *) sgs_hist_calc
     case ('SGS_HIST_CUMULATIVE')
        read (buff(equal_pos+1:), *) sgs_hist_cumulative

     case ('SGS_HIST_NSTART')
        read (buff(equal_pos+1:), *) sgs_hist_nstart
     case ('SGS_HIST_NSKIP')
        read (buff(equal_pos+1:), *) sgs_hist_nskip
     case ('SGS_HIST_LOC')
        call parse_vector( buff(equal_pos+1:), sgs_hist_nloc, sgs_hist_loc )

     case ('CS2_BMIN')
        read (buff(equal_pos+1:), *) cs2_bmin
     case ('CS2_BMAX')
        read (buff(equal_pos+1:), *) cs2_bmax
     case ('CS2_NBINS')
        read (buff(equal_pos+1:), *) cs2_nbins

     case ('TN_BMIN')
        read (buff(equal_pos+1:), *) tn_bmin
     case ('TN_BMAX')
        read (buff(equal_pos+1:), *) tn_bmax
     case ('TN_NBINS')
        read (buff(equal_pos+1:), *) tn_nbins

     case ('NU_BMIN')
        read (buff(equal_pos+1:), *) nu_bmin
     case ('NU_BMAX')
        read (buff(equal_pos+1:), *) nu_bmax
     case ('NU_NBINS')
        read (buff(equal_pos+1:), *) nu_nbins

     case ('EE_BMIN')
        read (buff(equal_pos+1:), *) ee_bmin
     case ('EE_BMAX')
        read (buff(equal_pos+1:), *) ee_bmax
     case ('EE_NBINS')
        read (buff(equal_pos+1:), *) ee_nbins

     case default

        if(coord == 0) call mesg( sub_name, 'Found unused data value in ' // block_name // ' block: ' // buff(1:equal_pos-1) )
     end select

  elseif( block_exit_pos == 1 ) then

     return

  else

     call error( sub_name, block_name // ' data block not formatted correctly: ' // buff(1:equal_pos-1) )

  endif

enddo

return
end subroutine  sgs_hist_block

$if ($TURBINES)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine turbines_block()
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use turbines_base
implicit none

character(*), parameter :: block_name = 'TURBINES'

do 

  call readline( lun, line, buff, block_entry_pos, block_exit_pos, &
                 equal_pos, ios )
  if (ios /= 0) call error( sub_name, 'Bad read in block')

  if( block_exit_pos == 0 ) then

     ! Check that the data entry conforms to correct format
     call checkentry()

     select case (uppercase(buff(1:equal_pos-1)))

     case ('NUM_X')
        read (buff(equal_pos+1:), *) num_x
     case ('NUM_Y')
        read (buff(equal_pos+1:), *) num_y

     case ('DIA_ALL')
        read (buff(equal_pos+1:), *) dia_all
     case ('HEIGHT_ALL')
        read (buff(equal_pos+1:), *) height_all
     case ('THK_ALL')
        read (buff(equal_pos+1:), *) thk_all

     case ('ORIENTATION')
        read (buff(equal_pos+1:), *) orientation
     case ('STAG_PERC')
        read (buff(equal_pos+1:), *) stag_perc

     case ('THETA1_ALL')
        read (buff(equal_pos+1:), *) theta1_all
     case ('THETA2_ALL')
        read (buff(equal_pos+1:), *) theta2_all

     case ('CT_PRIME')
        read (buff(equal_pos+1:), *) Ct_prime
     case ('CT_NOPRIME')
        read (buff(equal_pos+1:), *) Ct_noprime

     case ('T_AVG_DIM')
        read (buff(equal_pos+1:), *) T_avg_dim

     case ('ALPHA')
        read (buff(equal_pos+1:), *) alpha
     case ('TRUNC')
        read (buff(equal_pos+1:), *) trunc
     case ('FILTER_CUTOFF')
        read (buff(equal_pos+1:), *) filter_cutoff
     case ('TURBINE_CUMULATIVE_TIME')
        read (buff(equal_pos+1:), *) turbine_cumulative_time
     case ('TBASE')
        read (buff(equal_pos+1:), *) tbase
     case default

        if(coord == 0) call mesg( sub_name, 'Found unused data value in ' // block_name // ' block: ' // buff(1:equal_pos-1) )
     end select

  elseif( block_exit_pos == 1 ) then

     return

  else

     call error( sub_name, block_name // ' data block not formatted correctly: ' // buff(1:equal_pos-1) )

  endif

enddo

return
end subroutine  turbines_block
$endif !($TURBINES)

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine checkentry()
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

if( equal_pos == 0 ) call error( sub_name, 'Bad read in block at line', line, ': ' // trim(adjustl(buff)))
!--invalid if nothing after equals
if (len_trim (buff) == equal_pos) call error (sub_name, 'nothing after equals sign in line', line) 

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
subroutine parse_vector_real( string, nelem, vector )
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use types, only : rprec
use messages
use string_util, only : split_string

implicit none

character (*), parameter :: sub_name = mod_name // '.parse_vector_real'

character(*), intent(in) :: string
integer, intent(out) :: nelem
real(rprec), allocatable, dimension(:), intent(inout) :: vector
character(CHAR_BUFF_LENGTH), dimension(:), allocatable :: svector

call split_string( string, delim_minor, nelem, svector )

if( allocated( vector ) ) then
   ! Check that things are consistent
   if( nelem /= size( vector ) ) call error( sub_name, 'mismatch in element number and vector size')
else
   ! Now allocate the output vector if not allocated outside of parse_vector_real
   allocate( vector( nelem ) )
endif

! Read the string vector into the vector
read( svector(1:nelem), * ) vector(1:nelem)

deallocate(svector)

return
end subroutine parse_vector_real

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine parse_vector_point3D( string, nelem, vector )
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use types, only : rprec, point3D_t
use messages
use string_util, only : split_string

implicit none

character (*), parameter :: sub_name = mod_name // '.parse_vector_point3d'

character(*), intent(in) :: string
integer, intent(out) :: nelem
type(point3D_t), allocatable, dimension(:), intent(inout) :: vector

character(CHAR_BUFF_LENGTH), allocatable, dimension(:) :: svector

integer :: n, nelem_minor
real(rprec), allocatable, dimension(:) :: vector_minor

allocate( vector_minor(3) )

! Split based on major delimiter (allocates svector)
call split_string( string, delim_major, nelem, svector )

if( allocated( vector ) ) then
   ! Check that things are consistent
   if( nelem /= size( vector ) ) call error( sub_name, 'mismatch in element number and vector size')
else
   ! Now allocate the output vector if not allocated outside of parse_vector_real
   allocate( vector( nelem ) )
endif

! Now parse result string 
do n=1, nelem

   ! Dimension of the minor vector
   nelem_minor = 3

   call parse_vector_real( svector(n), nelem_minor, vector_minor )

   ! Check that the number of elements has not been reset
   if( nelem_minor /= 3 ) call error( sub_name, 'vector not specified correctly')

   vector(n) = point3D_t( (/ vector_minor(1), vector_minor(2), vector_minor(3) /) )

enddo

deallocate(svector)

return
end subroutine parse_vector_point3D

end module input_util
