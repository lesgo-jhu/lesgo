!!
!!  Copyright (C) 2011-2017  Johns Hopkins University
!!
!!  This file is part of lesgo.
!!
!!  lesgo is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License as published by
!!  the Free Software Foundation, either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  lesgo is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with lesgo.  If not, see <http://www.gnu.org/licenses/>.
!!

!*******************************************************************************
module input_util
!*******************************************************************************
use types, only : rprec
use param, only : path, CHAR_BUFF_LENGTH
implicit none

save
private

public :: read_input_conf

character(*), parameter :: mod_name = 'input_util'

character(*), parameter :: input_conf = path // 'lesgo.conf'
character(*), parameter :: comment = '!'
character(*), parameter :: block_entry = '{'
character(*), parameter :: block_exit = '}'
character(*), parameter :: equal = '='
character(*), parameter :: esyntax = 'syntax error at line'

! Delimiters used for reading vectors and points
character(*), parameter :: delim_minor=','
character(*), parameter :: delim_major='//'

! Thresh hold for evaluating differences in floating point values.
real(rprec), parameter :: thresh = 1.0e-6_rprec

interface parse_vector
    module procedure parse_vector_real, parse_vector_point3D
end interface

contains

!*******************************************************************************
subroutine read_input_conf ()
!*******************************************************************************
use param
use messages
use string_util, only : eat_whitespace, uppercase
implicit none

integer, parameter :: lun = 1
character (CHAR_BUFF_LENGTH) :: buff
integer :: block_entry_pos, block_exit_pos, equal_pos
integer :: ios
integer :: line
logical :: exst

character(*), parameter :: sub_name = mod_name // '.read_input_conf'

! Check that the configuration file exists
inquire (file=input_conf, exist=exst)

if (exst) then
    open (lun, file=input_conf, action='read')
else
    call error (sub_name, 'file ' // input_conf // ' does not exist')
end if

line = 0
do
    call readline( lun, line, buff, block_entry_pos, block_exit_pos,           &
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
#ifdef PPLVLSET
        case ('LEVEL_SET')
            call level_set_block()
#endif
#ifdef PPTURBINES
        case ('TURBINES')
            call turbines_block()
#endif
        case default
            if (coord == 0) write(*,*) 'Found unused input block: '//          &
                buff(1:block_entry_pos-1)
            ! Now need to 'fast-forward' until we reach the end of the block
            do while ( block_exit_pos == 0 )
                call readline( lun, line, buff, block_entry_pos,               &
                    block_exit_pos, equal_pos, ios )
                ! exit if end of file is reached
                if (ios /= 0) exit
            enddo
    end select
end do

close (lun)

contains

!*******************************************************************************
subroutine domain_block()
!*******************************************************************************
use types, only : rprec
use param
implicit none

character(*), parameter :: block_name = 'DOMAIN'

integer :: ival_read
integer :: np
real(rprec) :: val_read

do
    call readline( lun, line, buff, block_entry_pos, block_exit_pos,           &
        equal_pos, ios )

    if (ios /= 0) call error( sub_name, 'Bad read in block')

    if (block_exit_pos == 0) then
        ! Check that the data entry conforms to correct format
        call checkentry()
        select case (uppercase(buff(1:equal_pos-1)))
            case ('NPROC')
                read (buff(equal_pos+1:), *) np
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
                if (coord == 0) write(*,*) 'Found unused data value in '       &
                    // block_name // ' block: ' // buff(1:equal_pos-1)
        end select
    elseif (block_exit_pos == 1) then
        ! Check or reset nproc based on MPI setup
#ifndef PPMPI
        ! Reset to nproc to 1 when not using MPI
        if( nproc /= 1 ) then
            nproc = 1
            call mesg(sub_name, 'Reseting nproc to: ', nproc)
        end if
#else
! check if run-time number of processes agrees with nproc parameter
        if (np /= nproc) then
            call error(sub_name, 'Runtime number of procs = ', nproc,          &
                ' is not equal to nproc = ', np)
        endif
#endif

        ! Set the processor owned vertical grid spacing
        nz = floor ( real( nz_tot, rprec ) / nproc ) + 1

        ! Recompute nz_tot to be compliant with computed nz
        ival_read = nz_tot
        nz_tot = ( nz - 1 ) * nproc + 1
        if (coord == 0 .AND. ival_read /= nz_tot )                             &
           write(*,*) 'Reseting Nz (total) to: ', nz_tot

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
        if (uniform_spacing) then
            ! Adjust L_y
            val_read = L_y
            L_y = ny * dx
            if (coord == 0 .AND. abs( val_read - L_y ) >= thresh)              &
                call mesg( sub_name, 'Reseting Ly to: ', L_y )

            ! Adjust L_z
            val_read = L_z
            L_z = (nz_tot - 1 ) * dx
            if (coord == 0 .AND. abs( val_read - L_z ) >= thresh)              &
                call mesg( sub_name, 'Reseting Lz to: ', L_z )
        endif

        ! Grid spacing (y and z directions)
        dy = L_y / ny
        dz = L_z / ( nz_tot - 1 )

        return
    else
        call error( sub_name, block_name //                                    &
            'data block not formatted correctly: ' // buff(1:equal_pos-1) )
    endif
enddo

end subroutine domain_block

!*******************************************************************************
subroutine model_block()
!*******************************************************************************
use param
implicit none

character(*), parameter :: block_name = 'MODEL'

do
    call readline( lun, line, buff, block_entry_pos, block_exit_pos,           &
        equal_pos, ios )

    if (ios /= 0) call error( sub_name, 'Bad read in block')

    if (block_exit_pos == 0) then

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
            case default
                if (coord == 0) write(*,*) 'Found unused data value in '       &
                    // block_name // ' block: ' // buff(1:equal_pos-1)
        end select
    elseif( block_exit_pos == 1 ) then
        return
    else
        call error( sub_name, block_name //                                    &
            ' data block not formatted correctly: ' // buff(1:equal_pos-1) )
    endif
enddo

end subroutine model_block

!*******************************************************************************
subroutine time_block()
!*******************************************************************************
use types, only : rprec
use param
implicit none

character(*), parameter :: block_name = 'TIME'

do
    call readline( lun, line, buff, block_entry_pos, block_exit_pos,           &
        equal_pos, ios )

    if (ios /= 0) call error( sub_name, 'Bad read in block')

    if (block_exit_pos == 0) then

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
                if (coord == 0) write(*,*) 'Found unused data value in '       &
                    // block_name // ' block: ' // buff(1:equal_pos-1)
        end select

    elseif (block_exit_pos == 1) then
        ! Set dependent data
        if (.not. use_cfl_dt) then
            ! Set dimensional time step
            dt_dim = dt * z_i / u_star
            ! Set AB2 integration coefficients
            tadv1 = 1.5_rprec
            tadv2 = 1.0_rprec - tadv1
        endif

        return
    else
        call error( sub_name, block_name //                                    &
            ' data block not formatted correctly: ' // buff(1:equal_pos-1) )
    endif
enddo

end subroutine  time_block

!*******************************************************************************
subroutine flow_cond_block()
!*******************************************************************************
use param

#ifdef PPHIT
! Type hit has all the information inside
use hit_inflow, only : hit
#endif

implicit none

character(*), parameter :: block_name = 'FLOW_COND'

real(rprec) :: val_read

do
    call readline( lun, line, buff, block_entry_pos, block_exit_pos,           &
        equal_pos, ios )

    if (ios /= 0) call error( sub_name, 'Bad read in block')

    if (block_exit_pos == 0) then

        ! Check that the data entry conforms to correct format
        call checkentry()

        select case (uppercase(buff(1:equal_pos-1)))

            case ('INITU')
                read (buff(equal_pos+1:), *) initu
            case ('INILAG')
                read (buff(equal_pos+1:), *) inilag
            case ('LBC_MOM')
                Read (buff(equal_pos+1:), *) lbc_mom
            case ('UBC_MOM')
                Read (buff(equal_pos+1:), *) ubc_mom
            case ('UBOT')
                Read (buff(equal_pos+1:), *) ubot
            case ('UTOP')
                Read (buff(equal_pos+1:), *) utop
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
            case ('USE_MEAN_P_FORCE')
                read (buff(equal_pos+1:), *) use_mean_p_force
            case ('EVAL_MEAN_P_FORCE')
                read (buff(equal_pos+1:), *) eval_mean_p_force
            case ('MEAN_P_FORCE')
                read (buff(equal_pos+1:), *) mean_p_force
            case ('USE_RANDOM_FORCE')
                read (buff(equal_pos+1:), *) use_random_force
            case ('STOP_RANDOM_FORCE')
                read (buff(equal_pos+1:), *) stop_random_force
            case ('RMS_RANDOM_FORCE')
                read (buff(equal_pos+1:), *) rms_random_force

#ifdef PPHIT
            ! Read the input for HIT case

            ! Turbulence intensity input and output
            case('UP_IN')
                read (buff(equal_pos+1:), *) hit % up_in
            case('TI_OUT')
                read (buff(equal_pos+1:), *) hit % TI_out

            ! Domain length of HIT input
            case ('LX_HIT')
                read (buff(equal_pos+1:), *) hit % Lx
            case ('LY_HIT')
                read (buff(equal_pos+1:), *) hit % Ly
            case ('LZ_HIT')
                read (buff(equal_pos+1:), *) hit % Lz

            ! Number of grid points
            case ('NX_HIT')
                read (buff(equal_pos+1:), *) hit % Nx
            case ('NY_HIT')
                read (buff(equal_pos+1:), *) hit % Ny
            case ('NZ_HIT')
                read (buff(equal_pos+1:), *) hit % Nz

            ! Names of the input files
            case ('U_FILE')
                read (buff(equal_pos+1:), *) hit % u_file
            case ('V_FILE')
                read (buff(equal_pos+1:), *) hit % v_file
            case ('W_FILE')
                read (buff(equal_pos+1:), *) hit % w_file
#endif

            case default
                if (coord == 0) write(*,*) 'Found unused data value in '       &
                    // block_name // ' block: ' // buff(1:equal_pos-1)
        end select
    elseif (block_exit_pos == 1) then
        if( use_mean_p_force .AND. eval_mean_p_force ) then
            val_read = mean_p_force
            ! Evaluate the mean pressure force
            mean_p_force = 1.0_rprec / L_z
            if (coord == 0 .AND. abs( val_read - mean_p_force ) >= thresh)     &
                call mesg(sub_name, 'Reseting mean_p_force to: ', mean_p_force)
        endif
        return
    else
        call error( sub_name, block_name //                                    &
        ' data block not formatted correctly: ' // buff(1:equal_pos-1) )
    endif
enddo

end subroutine  flow_cond_block

!*******************************************************************************
subroutine output_block()
!*******************************************************************************
use param
implicit none

character(*), parameter :: block_name = 'OUTPUT'

do
    call readline(lun, line, buff, block_entry_pos, block_exit_pos,            &
        equal_pos, ios )

    if (ios /= 0) call error( sub_name, 'Bad read in block')

    if (block_exit_pos == 0) then

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
            ! Only read if point data is to be recorded. It is important that
            ! the following only be used elsewhere in the code if
            ! point_calc=.true.; also it is required that point_calc be listed
            ! before the rest in lesgo.conf
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
            case default
                if (coord == 0) write(*,*) 'Found unused data value in '       &
                    // block_name // ' block: ' // buff(1:equal_pos-1)
        end select
    elseif (block_exit_pos == 1 ) then
        return
    else
        call error( sub_name, block_name //                                    &
            ' data block not formatted correctly: ' // buff(1:equal_pos-1) )
    endif
enddo

end subroutine  output_block

#ifdef PPLVLSET
!*******************************************************************************
subroutine level_set_block()
!*******************************************************************************
use level_set_base
implicit none

character(*), parameter :: block_name = 'LEVEL_SET'

do
    call readline( lun, line, buff, block_entry_pos, block_exit_pos,           &
        equal_pos, ios )

    if (ios /= 0) call error( sub_name, 'Bad read in block')

    if (block_exit_pos == 0) then

        ! Check that the data entry conforms to correct format
        call checkentry()

        select case (uppercase(buff(1:equal_pos-1)))
            case ('GLOBAL_CA_CALC')
                read (buff(equal_pos+1:), *) global_CA_calc
            case ('GLOBAL_CA_NSKIP')
                read (buff(equal_pos+1:), *) global_CA_nskip
            case ('VEL_BC')
                read (buff(equal_pos+1:), *) vel_bc
            case ('USE_LOG_PROFILE')
                read (buff(equal_pos+1:), *) use_log_profile
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
            case ('USE_TREES')
                read (buff(equal_pos+1:), *) use_trees
#ifdef PPMPI
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
#endif
            case default
                if (coord == 0) write(*,*) 'Found unused data value in '       &
                    // block_name // ' block: ' // buff(1:equal_pos-1)
        end select
    elseif( block_exit_pos == 1 ) then
        return
    else
        call error( sub_name, block_name //                                    &
            ' data block not formatted correctly: ' // buff(1:equal_pos-1) )

    endif
enddo

end subroutine  level_set_block
#endif

#ifdef PPTURBINES
!*******************************************************************************
subroutine turbines_block()
!*******************************************************************************
use turbines
implicit none

character(*), parameter :: block_name = 'TURBINES'

do
    call readline( lun, line, buff, block_entry_pos, block_exit_pos,           &
        equal_pos, ios )

    if (ios /= 0) call error( sub_name, 'Bad read in block')

    if (block_exit_pos == 0) then

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
            case ('CP_PRIME')
                read (buff(equal_pos+1:), *) Cp_prime
            case ('READ_PARAM')
                read (buff(equal_pos+1:), *) read_param
            case ('DYN_THETA1')
                read (buff(equal_pos+1:), *) dyn_theta1
            case ('DYN_THETA2')
                read (buff(equal_pos+1:), *) dyn_theta2
            case ('DYN_BETA')
                read (buff(equal_pos+1:), *) dyn_beta
            case ('DYN_TORQUE_GAIN')
                read (buff(equal_pos+1:), *) dyn_torque_gain
            case ('T_AVG_DIM')
                read (buff(equal_pos+1:), *) T_avg_dim
            case ('ALPHA')
                read (buff(equal_pos+1:), *) alpha
            case ('FILTER_CUTOFF')
                read (buff(equal_pos+1:), *) filter_cutoff
            case ('TBASE')
                read (buff(equal_pos+1:), *) tbase
            case ('RHO')
                read (buff(equal_pos+1:),*) rho
            case ('INERTIA_ALL')
                read (buff(equal_pos+1:),*) inertia_all
            case ('USE_WAKE_MODEL')
                read (buff(equal_pos+1:),*) use_wake_model
            case ('TORQUE_GAIN')
                read (buff(equal_pos+1:),*) torque_gain
            case ('TAU_U_INFTY')
                read (buff(equal_pos+1:),*) tau_U_infty
            case ('SIGMA_DU')
                read (buff(equal_pos+1:),*) sigma_du
            case ('SIGMA_K')
                read (buff(equal_pos+1:),*) sigma_k
            case ('SIGMA_UHAT')
                read (buff(equal_pos+1:),*) sigma_uhat
            case ('NUM_ENSEMBLE')
                read (buff(equal_pos+1:),*) num_ensemble
            case ('USE_RECEDING_HORIZON')
                read (buff(equal_pos+1:),*) use_receding_horizon
            case ('SOLVER')
                read (buff(equal_pos+1:),*) solver
            case ('ADVANCEMENT_BASE')
                read (buff(equal_pos+1:),*) advancement_base
            case ('HORIZON_TIME')
                read (buff(equal_pos+1:),*) horizon_time
            case ('MAX_ITER')
                read (buff(equal_pos+1:),*) max_iter
            case ('SPEED_PENALTY')
                read (buff(equal_pos+1:),*) speed_penalty
            case ('OMEGA_MIN')
                read (buff(equal_pos+1:),*) omega_min
            case ('OMEGA_MAX')
                read (buff(equal_pos+1:),*) omega_max
            case ('BETA_PENALTY')
                read (buff(equal_pos+1:),*) beta_penalty
            case ('BETA_STAR')
                read (buff(equal_pos+1:),*) beta_star
            case ('TSR_PENALTY')
                read (buff(equal_pos+1:),*) tsr_penalty
            case ('LAMBDA_PRIME_STAR')
                read (buff(equal_pos+1:),*) lambda_prime_star
            case default
                if (coord == 0) write(*,*) 'Found unused data value in '       &
                    // block_name // ' block: ' // buff(1:equal_pos-1)
        end select
    elseif (block_exit_pos == 1) then
        return
    else
        call error( sub_name, block_name //                                    &
            ' data block not formatted correctly: ' // buff(1:equal_pos-1) )
    endif
enddo

end subroutine turbines_block
#endif

!*******************************************************************************
subroutine checkentry()
!*******************************************************************************
implicit none

if( equal_pos == 0 ) call error( sub_name,                                     &
    'Bad read in block at line', line, ': ' // trim(adjustl(buff)))
! invalid if nothing after equals
if (len_trim (buff) == equal_pos) call error (sub_name,                        &
    'nothing after equals sign in line', line)

end subroutine checkentry

end subroutine read_input_conf

!*******************************************************************************
subroutine readline(lun, line, buff, block_entry_pos,                          &
    block_exit_pos, equal_pos, ios )
!*******************************************************************************
!
! This subroutine reads the specified line and determines the attributes
! of the contents of the line.
!
use string_util, only : eat_whitespace
implicit none

integer, intent(in) :: lun
integer, intent(inout) :: line
character(*), intent(inout) :: buff
integer, intent(out) :: block_entry_pos, block_exit_pos, equal_pos, ios

block_entry_pos = 0
block_exit_pos = 0
equal_pos = 0
ios = -1

do
    line = line + 1
    read (lun, '(a)', iostat=ios) buff
    if (ios /= 0) exit

    call eat_whitespace (buff)

    if (verify (buff, ' ') == 0) cycle              ! drop blank lines

    if (buff (1:len (comment)) == comment) cycle    ! drop comment lines

    block_entry_pos = index( buff, block_entry )
    block_exit_pos = index( buff, block_exit )
    equal_pos = index( buff, equal )

    exit
enddo

end subroutine readline

!*******************************************************************************
subroutine parse_vector_real( string, nelem, vector )
!*******************************************************************************
use types, only : rprec
use messages
use string_util, only : split_string

implicit none

character(*), parameter :: sub_name = mod_name // '.parse_vector_real'

character(*), intent(in) :: string
integer, intent(out) :: nelem
real(rprec), allocatable, dimension(:), intent(inout) :: vector
character(CHAR_BUFF_LENGTH), dimension(:), allocatable :: svector

call split_string( string, delim_minor, nelem, svector )

if (allocated( vector )) then
    ! Check that things are consistent
    if (nelem /= size( vector ) ) call error(sub_name,                         &
        'mismatch in element number and vector size')
else
    ! Now allocate the output vector if not allocated
    ! outside of parse_vector_real
    allocate( vector( nelem ) )
endif

! Read the string vector into the vector
read( svector(1:nelem), *) vector(1:nelem)

deallocate(svector)

end subroutine parse_vector_real

!*******************************************************************************
subroutine parse_vector_point3D( string, nelem, vector )
!*******************************************************************************
use types, only : rprec, point3D_t
use messages
use string_util, only : split_string

implicit none

character(*), parameter :: sub_name = mod_name // '.parse_vector_point3d'

character(*), intent(in) :: string
integer, intent(out) :: nelem
type(point3D_t), allocatable, dimension(:), intent(inout) :: vector
character(CHAR_BUFF_LENGTH), allocatable, dimension(:) :: svector
integer :: n, nelem_minor
real(rprec), allocatable, dimension(:) :: vector_minor

allocate( vector_minor(3) )

! Split based on major delimiter (allocates svector)
call split_string( string, delim_major, nelem, svector )

if (allocated( vector )) then
    ! Check that things are consistent
    if( nelem /= size( vector ) ) call error( sub_name,                        &
        'mismatch in element number and vector size')
    else
    ! Now allocate the output vector if not allocated
    ! outside of parse_vector_real
    allocate( vector( nelem ) )
endif

! Now parse result string
do n = 1, nelem
   ! Dimension of the minor vector
   nelem_minor = 3

   call parse_vector_real( svector(n), nelem_minor, vector_minor )

   ! Check that the number of elements has not been reset
   if( nelem_minor /= 3 ) call error(sub_name, 'vector not specified correctly')

   vector(n) = point3D_t((/vector_minor(1), vector_minor(2), vector_minor(3)/))
enddo

deallocate(svector)

end subroutine parse_vector_point3D

end module input_util
