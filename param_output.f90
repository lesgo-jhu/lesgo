!**********************************************************************
subroutine param_output()
!**********************************************************************
use param
$if($RNS_LS)
use rns_base_ls
$endif
$if($CYL_SKEW_LS)
use cyl_skew_base_ls
$endif
$if($LVLSET)
use level_set_base
$endif

implicit none

integer :: n

character(*), parameter :: fname = 'lesgo_param.out'

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

open (unit = 2,file = fname, status='unknown',form='formatted', &
  action='write',position='rewind')

write(2,c_fmt) '**********************************************************************'
write(2,c_fmt) 'PARAM'
write(2,c_fmt) '**********************************************************************'
write(2,c_fmt) ''
write(2,c_fmt) '---------------------------------------------------'
write(2,c_fmt) 'MPI PARAMETERS'
write(2,c_fmt) '---------------------------------------------------'
write(2,i_fmt) 'nproc : ', nproc
write(2,c_fmt) ''
write(2,c_fmt) '---------------------------------------------------'
write(2,c_fmt) 'COMPUTATIONAL DOMAIN PARAMETERS'
write(2,c_fmt) '---------------------------------------------------'
write(2,x4i_fmt) 'nx, ny, nz, nz_tot : ', nx, ny, nz, nz_tot
write(2,f_fmt) 'z_i : ', z_i
write(2,x3f_fmt) 'L_x, L_y, L_z : ', L_x, L_y, L_z
write(2,x3f_fmt) 'dx, dy, dz : ', dx, dy, dz
write(2,c_fmt) ''
write(2,c_fmt) '---------------------------------------------------'
write(2,c_fmt) 'MODEL PARAMETERS'
write(2,c_fmt) '---------------------------------------------------'
write(2,x3i_fmt) 'model, models, nnn : ', model, models, nnn
write(2,f_fmt) 'Co : ', Co
write(2,i_fmt) 'cs_count : ', cs_count
write(2,i_fmt) 'ifilter : ', ifilter
write(2,x2f_fmt) 'u_star, Pr : ', u_star, Pr
write(2,f_fmt) 'vonk : ', vonk
write(2,l_fmt) 'coriolis_forcing : ', coriolis_forcing
write(2,x3f_fmt) 'coriol : ', coriol, ug, vg
write(2,f_fmt) 'nu_molec : ', nu_molec
write(2,x3l_fmt) 'molec, sgs, dns_bc : ', molec, sgs, dns_bc
write(2,c_fmt) ''
write(2,c_fmt) '---------------------------------------------------'
write(2,c_fmt) 'TIMESTEP PARAMETERS'
write(2,c_fmt) '---------------------------------------------------'
write(2,i_fmt) 'nsteps : ', nsteps
$if($CFL_DT)
write(2,f_fmt) 'cfl : ', cfl
$else
write(2,f_fmt) 'dt : ', dt
$endif
write(2,l_fmt) 'cumulative_time : ', cumulative_time
write(2, x2c_fmt) 'fcumulative_time : ', fcumulative_time
write(2,c_fmt) ''
write(2,c_fmt) '---------------------------------------------------'
write(2,c_fmt) 'BOUNDARY/INITIAL CONDITION PARAMETERS'
write(2,c_fmt) '---------------------------------------------------'
write(2,l_fmt) 'initu : ', initu
write(2,l_fmt) 'inilag : ', inilag
write(2,i_fmt) 'ubc : ', ubc
write(2,x2c_fmt) 'lbc_mom : ', lbc_mom
write(2,f_fmt) 'zo : ', zo
write(2,l_fmt) 'inflow : ', inflow
write(2,l_fmt) 'use_fringe_forcing : ', use_fringe_forcing
write(2,f_fmt) 'buff_end : ', buff_end
write(2,f_fmt) 'buff_len : ', buff_len
write(2,f_fmt) 'face_avg : ', face_avg
write(2,l_fmt) 'read_inflow_file : ', read_inflow_file
write(2,l_fmt) 'write_inflow_file : ', write_inflow_file
write(2,i_fmt) 'jt_start_write : ', jt_start_write
write(2,l_fmt) 'force_top_bot : ', force_top_bot
write(2,l_fmt) 'inflow_sample_velocity : ', inflow_sample_velocity
write(2,f_fmt) 'inflow_sample_location : ', inflow_sample_location
write(2,l_fmt) 'use_mean_p_force : ', use_mean_p_force
write(2,f_fmt) 'mean_p_force : ', mean_p_force
write(2,c_fmt) ''
write(2,c_fmt) '---------------------------------------------------'
write(2,c_fmt) 'DATA OUTPUT PARAMETERS'
write(2,c_fmt) '---------------------------------------------------'
write(2,i_fmt) 'wbase : ', wbase
write(2,i_fmt) 'nenergy : ', nenergy
write(2,i_fmt) 'cfl_count : ', cfl_count
write(2,l_fmt) 'tavg_calc : ', tavg_calc
write(2,x2i_fmt) 'tavg_nstart, tavg_nend : ', tavg_nstart, tavg_nend
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

write(2,l_fmt) 'spectra_calc : ', spectra_calc
write(2,x2i_fmt) 'spectra_nstart, spectra_nend : ', spectra_nstart, spectra_nend
write(2,i_fmt) 'spectra_nloc : ', spectra_nloc
do n=1,spectra_nloc
  write(2,if_fmt) 'n, spectra_loc(n) : ', n, spectra_loc(n)
enddo

$if($LVLSET)
write(2,c_fmt) ''
write(2,c_fmt) '**********************************************************************'
write(2,c_fmt) 'LEVEL_SET_BASE'
write(2,c_fmt) '**********************************************************************'
write(2,c_fmt) ''
write(2,c_fmt) '---------------------------------------------------'
write(2,c_fmt) 'FORCE PARAMETERS'
write(2,c_fmt) '---------------------------------------------------'
write(2,l_fmt) 'global_CD_calc : ', global_CD_calc
write(2,i_fmt) 'Ldir : ', Ldir
write(2,c_fmt) ''
write(2,c_fmt) '---------------------------------------------------'
write(2,c_fmt) 'BC PARAMETERS'
write(2,c_fmt) '---------------------------------------------------'
write(2,l_fmt) 'vel_BC : ', vel_BC
write(2,l_fmt) 'use_log_profile : ', use_log_profile
write(2,l_fmt) 'use_enforce_un : ', use_enforce_un
write(2,l_fmt) 'physBC : ', physBC
write(2,f_fmt) 'z0 : ', z0
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
$endif

$if($RNS_LS)
write(2,c_fmt) ''
write(2,c_fmt) '**********************************************************************'
write(2,c_fmt) 'RNS_BASE_LS'
write(2,c_fmt) '**********************************************************************'
write(2,c_fmt) ''
write(2,c_fmt) '---------------------------------------------------'
write(2,c_fmt) 'RNS PARAMETERS'
write(2,c_fmt) '---------------------------------------------------'
write(2,i_fmt) 'rns_ntree : ', rns_ntree
write(2,x2c_fmt) 'rns_tree_layout : ', rns_tree_layout
write(2,i_fmt) 'temporal_weight : ', temporal_weight
write(2,f_fmt) 'Tconst : ', Tconst
write(2,i_fmt) 'weight_nstart : ', weight_nstart
write(2,i_fmt) 'temporal_model : ', temporal_model
write(2,i_fmt) 'spatial_model : ', spatial_model
write(2,i_fmt) 'output_nskip : ', output_nskip
write(2,i_fmt) 'CD_ramp_nstep : ', CD_ramp_nstep
write(2,f_fmt) 'alpha_width : ', alpha_width
write(2,f_fmt) 'alpha_dist : ', alpha_dist
write(2,f_fmt) 'chi_cutoff : ', chi_cutoff
write(2,i_fmt) 'ndim : ', ndim
$endif

$if($CYL_SKEW_LS)
write(2,c_fmt) ''
write(2,c_fmt) '**********************************************************************'
write(2,c_fmt) 'CYL_SKEW_BASE_LS'
write(2,c_fmt) '**********************************************************************'
write(2,c_fmt) ''
write(2,c_fmt) '---------------------------------------------------'
write(2,c_fmt) 'CYL_SKEW TREE PARAMETERS'
write(2,c_fmt) '---------------------------------------------------'
write(2,f_fmt) 'zrot_angle : ', zrot_angle
write(2,f_fmt) 'skew_angle : ', skew_angle
write(2,i_fmt) 'ntree : ', ntree
write(2,i_fmt) 'ngen : ', ngen
write(2,i_fmt) 'ngen_reslv : ', ngen_reslv
write(2,i_fmt) 'nbranch : ', nbranch
write(2,f_fmt) 'd : ', d
write(2,f_fmt) 'l : ', l
write(2,f_fmt) 'offset : ', offset
write(2,f_fmt) 'scale_fact : ', scale_fact
write(2,l_fmt) 'use_bottom_surf : ', use_bottom_surf
write(2,f_fmt) 'z_bottom_surf : ', z_bottom_surf
write(2,l_fmt) 'filter_chi : ', filter_chi
write(2,f_fmt) 'filt_width : ', filt_width
$endif

close(2)

return
end subroutine param_output
