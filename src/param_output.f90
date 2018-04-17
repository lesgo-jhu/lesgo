!!
!!  Copyright (C) 2010-2013  Johns Hopkins University
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
subroutine param_output()
!*******************************************************************************
use param
#ifdef PPLVLSET
use level_set_base
#endif

implicit none

integer :: n

character(*), parameter :: fname = path // 'output/lesgo_param.out'

character(*), parameter :: c_fmt = '(a)'
character(*), parameter :: x2c_fmt = '(2a)'

character(*), parameter :: l_fmt = '(a,l)'
character(*), parameter :: x3l_fmt = '(a,3l)'

character(*), parameter :: i_fmt = '(a,i7)'
character(*), parameter :: x2i_fmt = '(a,2i7)'
character(*), parameter :: x3i_fmt = '(a,3i7)'
character(*), parameter :: x4i_fmt = '(a,4i7)'

character(*), parameter :: f_fmt = '(a,e15.7)'
character(*), parameter :: x2f_fmt = '(a,2e15.7)'
character(*), parameter :: x3f_fmt = '(a,3e15.7)'

character(*), parameter :: if_fmt='(a,i7,e15.7)'
character(*), parameter :: ix3f_fmt='(a,i7,3e15.7)'
character(13) :: ch

open (unit = 2,file = fname, status='unknown',form='formatted',                &
  action='write',position='rewind')

write(2,c_fmt) '**********************************************************************'
write(2,c_fmt) 'PARAM'
write(2,c_fmt) '**********************************************************************'
write(2,c_fmt) ''
write(2,c_fmt) '---------------------------------------------------'
write(2,c_fmt) 'COMPUTATIONAL DOMAIN PARAMETERS'
write(2,c_fmt) '---------------------------------------------------'
write(2,i_fmt) 'nproc : ', nproc
write(2,x4i_fmt) 'nx, ny, nz, nz_tot : ', nx, ny, nz, nz_tot
write(2,f_fmt) 'z_i : ', z_i
write(2,x3f_fmt) 'L_x, L_y, L_z : ', L_x, L_y, L_z
write(2,x3f_fmt) 'dx, dy, dz : ', dx, dy, dz
write(2,c_fmt) ''
write(2,c_fmt) '---------------------------------------------------'
write(2,c_fmt) 'MODEL PARAMETERS'
write(2,c_fmt) '---------------------------------------------------'
write(2,x3i_fmt) 'sgs_model, wall_damp_exp : ', sgs_model, wall_damp_exp
write(2,f_fmt) 'Co : ', Co
write(2,i_fmt) 'cs_count : ', cs_count
write(2,i_fmt) 'ifilter : ', ifilter
write(2,x2f_fmt) 'u_star : ', u_star
write(2,f_fmt) 'vonk : ', vonk
write(2,l_fmt) 'coriolis_forcing : ', coriolis_forcing
write(2,x3f_fmt) 'coriol : ', coriol, ug, vg
write(2,f_fmt) 'nu_molec : ', nu_molec
write(2,x3l_fmt) 'molec, sgs : ', molec, sgs
write(2,c_fmt) ''
write(2,c_fmt) '---------------------------------------------------'
write(2,c_fmt) 'TIMESTEP PARAMETERS'
write(2,c_fmt) '---------------------------------------------------'
write(2,i_fmt) 'nsteps : ', nsteps
write(2,i_fmt) 'runtime : ', runtime
#ifdef PPCFL_DT
write(2,f_fmt) 'cfl : ', cfl
#else
write(2,f_fmt) 'dt : ', dt
#endif
write(2,l_fmt) 'cumulative_time : ', cumulative_time
write(2, x2c_fmt) 'fcumulative_time : ', fcumulative_time
write(2,c_fmt) ''
write(2,c_fmt) '---------------------------------------------------'
write(2,c_fmt) 'BOUNDARY/INITIAL CONDITION PARAMETERS'
write(2,c_fmt) '---------------------------------------------------'
write(2,l_fmt) 'initu : ', initu
write(2,l_fmt) 'inilag : ', inilag
write(2,i_fmt) 'lbc_mom : ', lbc_mom
write(2,i_fmt) 'ubc_mom : ', ubc_mom
write(2,f_fmt) 'ubot : ', ubot
write(2,f_fmt) 'utop : ', utop
write(2,f_fmt) 'zo : ', zo
write(2,l_fmt) 'inflow : ', inflow
write(2,f_fmt) 'fringe_region_end : ', fringe_region_end
write(2,f_fmt) 'fringe_region_len : ', fringe_region_len
write(2,f_fmt) 'inflow_velocity : ', inflow_velocity
write(2,l_fmt) 'use_mean_p_force : ', use_mean_p_force
write(2,f_fmt) 'mean_p_force : ', mean_p_force
write(2,l_fmt) 'use_random_force : ', use_random_force
write(2,i_fmt) 'stop_random_force : ', stop_random_force
write(2,f_fmt) 'rms_random_force : ', rms_random_force
write(2,c_fmt) ''
write(2,c_fmt) '---------------------------------------------------'
write(2,c_fmt) 'DATA OUTPUT PARAMETERS'
write(2,c_fmt) '---------------------------------------------------'
#if ( defined(PPWRITE_BIG_ENDIAN) || defined(PPWRITE_LITTLE_ENDIAN) )
write(2,x2c_fmt) 'write_endian : ', write_endian
#elif (defined(__INTEL_COMPILER))
! According to https://software.intel.com/en-us/node/524834:
! Intel Fortran expects numeric data to be in native little endian order.
write(2,x2c_fmt) 'write_endian : ', 'LITTLE_ENDIAN'
#else
inquire(2,convert=ch)
write(2,x2c_fmt) 'write_endian : ', ch
#endif
#if ( defined(PPREAD_BIG_ENDIAN) || defined(PPREAD_LITTLE_ENDIAN) )
write(2,x2c_fmt) 'read_endian : ', write_endian
#elif (defined(__INTEL_COMPILER))
! According to https://software.intel.com/en-us/node/524834:
! Intel Fortran expects numeric data to be in native little endian order.
write(2,x2c_fmt) 'read_endian : ', 'LITTLE_ENDIAN'
#else
inquire(2,convert=ch)
write(2,x2c_fmt) 'read_endian : ', ch
#endif
write(2,i_fmt) 'wbase : ', wbase
write(2,i_fmt) 'nenergy : ', nenergy
write(2,i_fmt) 'lag_cfl_count : ', lag_cfl_count
write(2,l_fmt) 'checkpoint_data : ', checkpoint_data
write(2,i_fmt) 'checkpoint_nskip : ', checkpoint_nskip
write(2,l_fmt) 'tavg_calc : ', tavg_calc
write(2,x3i_fmt) 'tavg_nstart, tavg_nend, tavg_nskip : ', tavg_nstart, tavg_nend, tavg_nskip
write(2,l_fmt) 'point_calc : ', point_calc
write(2,x3i_fmt) 'point_nstart, point_nend, point_nskip : ', point_nstart, point_nend, point_nskip
write(2,i_fmt) 'point_nloc : ', point_nloc
do n=1,point_nloc
  write(2,ix3f_fmt) 'n, point_loc(n)%xyz : ', n, point_loc(n)%xyz
enddo
write(2,l_fmt) 'domain_calc : ', domain_calc
write(2,x3i_fmt) 'domain_nstart, domain_nend, domain_nskip : ', domain_nstart, domain_nend, domain_nskip

write(2,l_fmt) 'xplane_calc : ', xplane_calc
write(2,x3i_fmt) 'xplane_nstart, xplane_nend, xplane_nskip : ', xplane_nstart, xplane_nend, xplane_nskip
write(2,i_fmt) 'xplane_nloc : ', xplane_nloc
do n=1,xplane_nloc
  write(2,if_fmt) 'n, xplane_loc(n) : ', n, xplane_loc(n)
enddo

write(2,l_fmt) 'yplane_calc : ', yplane_calc
write(2,x3i_fmt) 'yplane_nstart, yplane_nend, yplane_nskip : ', yplane_nstart, yplane_nend, yplane_nskip
write(2,i_fmt) 'yplane_nloc : ', yplane_nloc
do n=1,yplane_nloc
  write(2,if_fmt) 'n, yplane_loc(n) : ', n, yplane_loc(n)
enddo

write(2,l_fmt) 'zplane_calc : ', zplane_calc
write(2,x3i_fmt) 'zplane_nstart, zplane_nend, zplane_nskip : ', zplane_nstart, zplane_nend, zplane_nskip
write(2,i_fmt) 'zplane_nloc : ', zplane_nloc
do n=1,zplane_nloc
  write(2,if_fmt) 'n, zplane_loc(n) : ', n, zplane_loc(n)
enddo

#ifdef PPLVLSET
write(2,c_fmt) ''
write(2,c_fmt) '**********************************************************************'
write(2,c_fmt) 'LEVEL_SET_BASE'
write(2,c_fmt) '**********************************************************************'
write(2,c_fmt) ''
write(2,c_fmt) '---------------------------------------------------'
write(2,c_fmt) 'FORCE PARAMETERS'
write(2,c_fmt) '---------------------------------------------------'
write(2,l_fmt) 'global_CA_calc : ', global_CA_calc
write(2,i_fmt) 'global_CA_nskip : ', global_CA_nskip
write(2,c_fmt) ''
write(2,c_fmt) '---------------------------------------------------'
write(2,c_fmt) 'BC PARAMETERS'
write(2,c_fmt) '---------------------------------------------------'
write(2,l_fmt) 'vel_BC : ', vel_BC
write(2,l_fmt) 'use_log_profile : ', use_log_profile
write(2,l_fmt) 'use_enforce_un : ', use_enforce_un
write(2,l_fmt) 'physBC : ', physBC
write(2,f_fmt) 'zo_level_set : ', zo_level_set
write(2,c_fmt) ''
write(2,c_fmt) '---------------------------------------------------'
write(2,c_fmt) 'SMOOTHING PARAMETERS'
write(2,c_fmt) '---------------------------------------------------'
write(2,l_fmt) 'use_smooth_tau : ', use_smooth_tau
write(2,l_fmt) 'use_extrap_tau_log : ', use_extrap_tau_log
write(2,l_fmt) 'use_extrap_tau_simple : ', use_extrap_tau_simple
write(2,l_fmt) 'use_modify_dutdn : ', use_modify_dutdn
write(2,x2c_fmt) 'smooth_mode : ', smooth_mode
write(2,c_fmt) ''
write(2,c_fmt) '---------------------------------------------------'
write(2,c_fmt) 'SGS PARAMETERS'
write(2,c_fmt) '---------------------------------------------------'
write(2,l_fmt) 'lag_dyn_modify_beta : ', lag_dyn_modify_beta
#endif

close(2)

return
end subroutine param_output
