!
!!  Copyright (C) 2009-2013  Johns Hopkins University
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

!///////////////////////////////////////////////////////////////////////////////
module io
!///////////////////////////////////////////////////////////////////////////////
use types,only:rprec
use param, only : ld, nx, ny, nz, nz_tot, path,  &
                  coord, rank, nproc, jt_total, total_time, &
                  total_time_dim, lbz, jzmin, jzmax
use param, only : cumulative_time, fcumulative_time
use sim_param, only : w, dudz, dvdz
use sgs_param,only:Cs_opt2
use string_util
use messages

implicit none

save
private

!!$public openfiles,output_loop,output_final,                   &
!!$     inflow_write, avg_stats
public jt_total, openfiles, closefiles, energy, output_loop, output_final

public output_init

character (*), parameter :: mod_name = 'io'

! Output file id's (see README for assigned values)
integer :: ke_fid

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine openfiles()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use param, only : use_cfl_dt, dt, cfl_f, path
implicit none
include 'tecryte.h'

logical :: exst

! Temporary values used to read time step and CFL from file
real(rprec) :: dt_r, cfl_r

! Open output files using tecryte library
! Kinetic energy (check_ke.dat)
ke_fid = open_file( path // 'output/check_ke.dat', 'append', 'formatted' )

if (cumulative_time) then

  inquire (file=fcumulative_time, exist=exst)
  if (exst) then

    open (1, file=fcumulative_time)   
    read(1, *) jt_total, total_time, total_time_dim, dt_r, cfl_r
    close (1)
    
  else  !--assume this is the first run on cumulative time
    if( coord == 0 ) then
      write (*, *) '--> Assuming jt_total = 0, total_time = 0.0'
    endif
    jt_total = 0
    total_time = 0.0_rprec
    total_time_dim = 0.0_rprec
  end if

end if

! Update dynamic time stepping info if required; otherwise discard.
if( use_cfl_dt ) then
   dt = dt_r
   cfl_f = cfl_r
endif

return
end subroutine openfiles

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine closefiles()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! This subroutine is used to close all open files used by the main
! program. These files are opened by calling 'openfiles'.
!
implicit none

! Close kinetic energy file
close( ke_fid )

return
end subroutine closefiles

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine energy (ke)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use types,only:rprec
use param
use sim_param,only:u,v,w
use messages
$if ($XLF)
  use ieee_arithmetic  !--for NAN checking
$endif
implicit none

include 'tecryte.h'

character (*), parameter :: sub_name = 'energy'

integer, parameter :: NAN_MAX = 10
                      !--write this many NAN's before calling error (to aid
                      !  diagnosis of problem)

$if ($DEBUG)
logical, parameter :: DEBUG = .true.
$endif

!logical, parameter :: flush = .true.

integer :: jx, jy, jz
integer :: nan_count

$if($DEBUG)
logical :: nan
$endif

real(kind=rprec)::KE,temp_w
$if ($MPI)
  real (rprec) :: ke_global
$endif

! Initialize variables
nan_count = 0
ke=0._rprec

z_loop: do jz=1,nz-1
    do jy=1,ny
        do jx=1,nx
            
            temp_w = 0.5_rprec*(w(jx,jy,jz)+w(jx,jy,jz+1))
            ke = ke + (u(jx,jy,jz)**2+v(jx,jy,jz)**2+temp_w**2)
            
            $if ($DEBUG)
            if (DEBUG) then
                $if ($IFORT || $IFC)
                    nan = isnan (ke)
                $elsif ($XLF)
                    !--this is a bit verbose, should make into sub-program
                    if (ieee_support_datatype (ke)) then
                        if (ieee_support_nan (ke)) nan = ieee_is_nan (ke)
                    end if
                $endif
 
                if (nan) then
                    nan_count = nan_count + 1
                    write (*, *) 'NaN in ke at (jx, jy, jz) =', jx, jy, jz
                    write (*, *) 'jt_total = ', jt_total
                    write (*, *) 'u = ', u(jx, jy, jz)
                    write (*, *) 'v = ', v(jx, jy, jz)
                    write (*, *) 'w = ', w(jx, jy, jz)
                    if ( nan_count >= NAN_MAX ) exit z_loop
                end if
            end if
            $endif 
            
        end do
    end do
end do z_loop

! Perform spatial averaging
ke = ke*0.5_rprec/(nx*ny*(nz-1))

! Check if NaN's where found
if ( nan_count > 0 ) call error (sub_name, 'NaN found')

$if ($MPI)

  call mpi_reduce (ke, ke_global, 1, MPI_RPREC, MPI_SUM, 0, comm, ierr)
  if (rank == 0) then  !--note its rank here, not coord
    ke = ke_global/nproc
    !ke = ke_global
    !write (13, *) total_time, ke

  call write_real_data( ke_fid, 'formatted', 2, (/ total_time, ke /))

  end if
  !if (rank == 0) ke = ke_global/nproc  !--its rank here, not coord

$else

!  write (13, *) total_time, ke
  call write_real_data( ke_fid, 'formatted', 2, (/ total_time, ke /))

$endif

!if ( flush ) then
!  if (rank == 0) then
!    close (13)
!    open ( 13, file=path//'output/check_ke.out', position='append' )
!  end if
!end if

end subroutine energy

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine output_loop()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  This subroutine is called every time step and acts as a driver for 
!  computing statistics and outputing instantaneous data. No actual
!  calculations are performed here.
!
use param, only : nsteps, jt_total, dt
use param, only : checkpoint_data, checkpoint_nskip
use param, only : tavg_calc, tavg_nstart, tavg_nend, tavg_nskip
use param, only : spectra_calc, spectra_nstart, spectra_nend, spectra_nskip
use param, only : point_calc, point_nstart, point_nend, point_nskip
use param, only : domain_calc, domain_nstart, domain_nend, domain_nskip
use param, only : xplane_calc, xplane_nstart, xplane_nend, xplane_nskip
use param, only : yplane_calc, yplane_nstart, yplane_nend, yplane_nskip
use param, only : zplane_calc, zplane_nstart, zplane_nend, zplane_nskip
use stat_defs, only : tavg_initialized, spectra_initialized
use stat_defs, only: tavg_dt, spectra_dt
implicit none

! Determine if we are to checkpoint intermediate times
if( checkpoint_data ) then
   ! Now check if data should be checkpointed this time step
   ! Don't checkpoint for last time step since this is done in output_final
!   if ( jt_total < nsteps .and. modulo (jt_total, checkpoint_nskip) == 0) call checkpoint()
   if ( modulo (jt_total, checkpoint_nskip) == 0) call checkpoint()
endif 

!  Determine if time summations are to be calculated
if (tavg_calc) then

  ! Are we between the start and stop timesteps?
  if ((jt_total >= tavg_nstart).and.(jt_total <= tavg_nend)) then

    ! Every timestep (between nstart and nend), add to tavg_dt
    tavg_dt = tavg_dt + dt

    ! Are we at the beginning or a multiple of nstart?
    if ( mod(jt_total-tavg_nstart,tavg_nskip)==0 ) then

      ! Check if we have initialized tavg
      if (.not.tavg_initialized) then
        if (coord == 0) then
          write(*,*) '-------------------------------'   
          write(*,"(1a,i9,1a,i9)") 'Starting running time summation from ', tavg_nstart, ' to ', tavg_nend
          write(*,*) '-------------------------------'   
        endif  ! coord==0
        call tavg_init()
      else
        call tavg_compute ()
      endif  ! init
     
    endif  ! mod of nskip

  endif  ! between nstart and nend

endif  ! tavg_calc
     

!  Determine if spectra are to be calculated
if (spectra_calc) then

  ! Are we between the start and stop timesteps?
  if ((jt_total >= spectra_nstart).and.(jt_total <= spectra_nend)) then

    ! Every timestep (between nstart and nend), add to spectra_dt
    spectra_dt = spectra_dt + dt

    ! Are we at the beginning or a multiple of nstart?
    if ( mod(jt_total-spectra_nstart,spectra_nskip)==0 ) then

      ! Check if we have initialized spectra
      if (.not.spectra_initialized) then
        if (coord == 0) then
         write(*,*) '-------------------------------'
         write(*,"(1a,i9,1a,i9)") 'Starting running spectra calculations from ', spectra_nstart, ' to ', spectra_nend
         write(*,*) '-------------------------------'
        endif  ! coord==0
      
        call spectra_init()
      else
        call spectra_compute ()
      endif  ! init
     
    endif  ! mod of nskip

  endif  ! between nstart and nend

endif  ! spectra_calc


!  Determine if instantaneous point velocities are to be recorded
if(point_calc) then
  if(jt_total >= point_nstart .and. jt_total <= point_nend .and. ( mod(jt_total-point_nstart,point_nskip)==0) ) then
    if(jt_total == point_nstart) then
      if (coord == 0) then   
        write(*,*) '-------------------------------'   
        write(*,"(1a,i9,1a,i9)") 'Writing instantaneous point velocities from ', point_nstart, ' to ', point_nend
        write(*,"(1a,i9)") 'Iteration skip:', point_nskip
        write(*,*) '-------------------------------'
      endif
    endif
    call inst_write(1)
  endif
endif
  
!  Determine if instantaneous domain velocities are to be recorded
if(domain_calc) then
  if(jt_total >= domain_nstart .and. jt_total <= domain_nend .and. ( mod(jt_total-domain_nstart,domain_nskip)==0) ) then
    if(jt_total == domain_nstart) then
      if (coord == 0) then        
        write(*,*) '-------------------------------'
        write(*,"(1a,i9,1a,i9)") 'Writing instantaneous domain velocities from ', domain_nstart, ' to ', domain_nend
        write(*,"(1a,i9)") 'Iteration skip:', domain_nskip
        write(*,*) '-------------------------------'
      endif

    endif
    call inst_write(2)
  endif
endif 

!  Determine if instantaneous x-plane velocities are to be recorded
if(xplane_calc) then
  if(jt_total >= xplane_nstart .and. jt_total <= xplane_nend .and. ( mod(jt_total-xplane_nstart,xplane_nskip)==0) ) then
    if(jt_total == xplane_nstart) then
      if (coord == 0) then        
        write(*,*) '-------------------------------'
        write(*,"(1a,i9,1a,i9)") 'Writing instantaneous x-plane velocities from ', xplane_nstart, ' to ', xplane_nend
        write(*,"(1a,i9)") 'Iteration skip:', xplane_nskip
        write(*,*) '-------------------------------'
      endif

    endif
    call inst_write(3)
  endif
endif  

!  Determine if instantaneous y-plane velocities are to be recorded
if(yplane_calc) then
  if(jt_total >= yplane_nstart .and. jt_total <= yplane_nend .and. ( mod(jt_total-yplane_nstart,yplane_nskip)==0) ) then
    if(jt_total == yplane_nstart) then
      if (coord == 0) then        
        write(*,*) '-------------------------------'
        write(*,"(1a,i9,1a,i9)") 'Writing instantaneous y-plane velocities from ', yplane_nstart, ' to ', yplane_nend
        write(*,"(1a,i9)") 'Iteration skip:', yplane_nskip
        write(*,*) '-------------------------------'
      endif
    endif
    call inst_write(4)
  endif
endif    

!  Determine if instantaneous z-plane velocities are to be recorded
if(zplane_calc) then
  if(jt_total >= zplane_nstart .and. jt_total <= zplane_nend .and. ( mod(jt_total-zplane_nstart,zplane_nskip)==0) ) then
    if(jt_total == zplane_nstart) then
      if (coord == 0) then        
        write(*,*) '-------------------------------'
        write(*,"(1a,i9,1a,i9)") 'Writing instantaneous z-plane velocities from ', zplane_nstart, ' to ', zplane_nend
        write(*,"(1a,i9)") 'Iteration skip:', zplane_nskip
        write(*,*) '-------------------------------'
      endif
    endif
    call inst_write(5)
  endif
endif


return
end subroutine output_loop

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine inst_write(itype)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! This subroutine is used to write all of the instantaneous data from
! lesgo to file. The types of data written are:
! 
!   points   : itype=1
!   domain   : itype=2
!   x-planes : itype=3
!   y-planes : itype=4
!   z-planes : itype=5
!
! For the points and planar data, this subroutine writes using the
! locations specfied from the param module. All of the data is written
! using the tecryte library. If additional instantenous values are
! desired to be written, they should be done so using this subroutine.
!
! REMARK: It may be desired to convert this subroutine to a driver and
! have it call appropriate subroutines to perform the writing, just to
! clean things up a bit
!
use functions, only : linear_interp, trilinear_interp, interp_to_uv_grid, interp_to_w_grid
use param, only : path, output_path
use param, only : point_nloc, point_loc
use param, only : xplane_nloc, xplane_loc
use param, only : yplane_nloc, yplane_loc
use param, only : zplane_nloc, zplane_loc
use grid_defs, only : grid
use sim_param, only : u,v,w,dudx,dvdy,dwdz,theta, &
                      sal, dTdx, dTdy, dTdz, dSdx, dSdy, dSdz, buoyancy, &
                      divtx, divty, divtz   !Eshwan
use sim_param, only: sigma_theta, sigma_sal
$if($DEBUG)
use sim_param, only : p, dpdx, dpdy, dpdz
use sim_param, only : RHSx, RHSy, RHSz
$endif
use stat_defs, only : xplane, yplane, zplane, point
$if($MPI)
use mpi
use param, only : ld, ny, nz, MPI_RPREC, down, up, comm, status, ierr
$endif

$if($LVLSET)
use level_set_base, only : phi
use sim_param, only : fx, fy, fz, fxa, fya, fza
$endif
use param, only : dx,dy,dz
use param, only : sgs_model
use sgs_param, only : F_LM,F_MM,F_QN,F_NN,Tn_all,Beta,Cs_opt2,Nu_t,cs2_clips !Eshwan
use sgs_param, only : kappa_tt,Ds_opt2_t !Eshwan
use sgs_param, only : s_Beta_t, s_Tn_all_t, I_LM_t, I_MM_t, I_QN_t, I_NN_t !Eshwan
use sgs_param, only : ds2_clips_t !Eshwan
use sgs_param, only : kappa_ts,Ds_opt2_s !Eshwan
use sgs_param, only : s_Beta_s, s_Tn_all_s, I_LM_s, I_MM_s, I_QN_s, I_NN_s !Eshwan
use sgs_param, only : ds2_clips_s !Eshwan
implicit none

include 'tecryte.h'      

integer, intent(IN) :: itype

character (*), parameter :: sub_name = mod_name // '.inst_write'

character (64) :: fname
integer :: n, i, j, k
$if(not $OUTPUT_BINARY)
character (64) :: var_list
integer :: nvars
$endif

real(rprec), allocatable, dimension(:,:,:) :: ui, vi, wi
real(rprec), allocatable, dimension(:,:,:) :: w_uv
real(rprec), allocatable, dimension(:,:,:) :: u_w, v_w
real(rprec), allocatable, dimension(:,:,:) :: divtx_w, divty_w

$if($LVLSET)
real(rprec), allocatable, dimension(:,:,:) :: fx_tot, fy_tot, fz_tot
$endif

$if($DEBUG)
real(rprec), allocatable, dimension(:,:,:) :: divvel
$endif

real(rprec), pointer, dimension(:) :: x,y,z,zw

!$if($OUTPUT_EXTRA)
! Arrays used for outputing slices of LDSM variables
real(rprec), allocatable, dimension(:,:) :: F_LM_s,F_MM_s,F_QN_s,F_NN_s,beta_s,Cs_opt2_s,Nu_t_s
real(rprec), allocatable, dimension(:,:,:) :: F_LM_uv,F_MM_uv,F_QN_uv,F_NN_uv,beta_uv,Cs_opt2_uv,Nu_t_uv_interp
!$endif

! Nullify pointers
nullify(x,y,z,zw)

! Set grid pointers
x => grid % x
y => grid % y
z => grid % z
zw => grid % zw

!  Allocate space for the interpolated w values
allocate(w_uv(nx,ny,lbz:nz))

!  Make sure w has been interpolated to uv-grid
w_uv = interp_to_uv_grid(w(1:nx,1:ny,lbz:nz), lbz)

!Interpolate u and v on w-grid to be consistent with values output
allocate(u_w(nx,ny,lbz:nz))
allocate(v_w(nx,ny,lbz:nz))

u_w = interp_to_w_grid(u(1:nx,1:ny,lbz:nz), lbz)
v_w = interp_to_w_grid(v(1:nx,1:ny,lbz:nz), lbz)
u_w(:,:,0:1) = 0.0_rprec
v_w(:,:,0:1) = 0.0_rprec

!Interpolate divtx and divty on w-grid
allocate(divtx_w(nx,ny,lbz:nz))
allocate(divty_w(nx,ny,lbz:nz))

divtx_w = interp_to_w_grid(divtx(1:nx,1:ny,lbz:nz), lbz)
divty_w = interp_to_w_grid(divty(1:nx,1:ny,lbz:nz), lbz)
divtx_w(:,:,0:1) = 0.0_rprec
divty_w(:,:,0:1) = 0.0_rprec

!$if($OUTPUT_EXTRA)
!  Allocate arrays and interpolate to uv grid for LDSM output
if( sgs_model == 4 .or. sgs_model == 5 ) then

  if( itype == 3 .or. itype == 4 .or. itype == 5 ) then

    allocate( F_LM_uv(nx,ny,nz), F_MM_uv(nx,ny,nz) )
    allocate( beta_uv(nx,ny,nz), Cs_opt2_uv(nx,ny,nz) )
    allocate( Nu_t_uv_interp(nx,ny,nz) )

    F_LM_uv = interp_to_uv_grid( F_LM(1:nx,1:ny,1:nz), 1 )
    F_MM_uv = interp_to_uv_grid( F_MM(1:nx,1:ny,1:nz), 1 )
    beta_uv = interp_to_uv_grid( beta(1:nx,1:ny,1:nz), 1 )
    Cs_opt2_uv = interp_to_uv_grid( Cs_opt2(1:nx,1:ny,1:nz), 1 )
    Nu_t_uv_interp = interp_to_uv_grid( Nu_t(1:nx,1:ny,1:nz), 1 )

    if( sgs_model == 5) then

      allocate( F_QN_uv(nx, ny, nz), F_NN_uv(nx,ny,nz) )

      F_QN_uv = interp_to_uv_grid( F_QN(1:nx,1:ny,1:nz), 1 )
      F_NN_uv = interp_to_uv_grid( F_NN(1:nx,1:ny,1:nz), 1 )

    endif

  endif

endif       
!$endif

if(itype==1) then

  do n=1,point_nloc

    !  For parallel runs check if data is on correct proc
    $if ($MPI)
    if(point(n) % coord == coord) then
    $endif

    ! Want to replace with write based on fid
    call write_real_data(point(n) % fid, 'formatted', 4, (/ total_time, &
         trilinear_interp(u(1:nx,1:ny,lbz:nz), lbz, point_loc(n)%xyz), &
         trilinear_interp(v(1:nx,1:ny,lbz:nz), lbz, point_loc(n)%xyz), &
         trilinear_interp(w_uv(1:nx,1:ny,lbz:nz), lbz, point_loc(n)%xyz) /))
    

    $if ($MPI)
    endif
    $endif

  enddo
!  Instantaneous write for entire domain
elseif(itype==2) then

  !Eshwan: write temperature
  !////////////////////////////////////////////
  !/// WRITE TEMPERATURE                    ///
  !////////////////////////////////////////////
  $if ($BINARY)
  call string_splice (fname, output_path // 'output/binary_T.', jt_total,'.dat')  !Eshwan
  $else
  call string_splice(fname, output_path // 'output/T.',jt_total,'.dat')  !Eshwan
  $endif

  $if ($MPI)
  call string_concat (fname, '.c', coord)
  $endif 

  var_list = ' "x", "y", "z", "theta", "dTdx", "dTdy", "dTdz" '
  nvars = 7

  call write_tecplot_header_ND(fname,'rewind',nvars, (/Nx+1,Ny+1,Nz/), &
       trim(adjustl(var_list)), numtostr(coord, 6), 2, real(total_time,4))
  call write_real_data_3D(fname, 'append', 'formatted', 4, nx, ny, nz, &
       (/theta(1:nx,1:ny,1:nz), &
         dTdx(1:nx,1:ny,1:nz), &
         dTdy(1:nx,1:ny,1:nz), &
         dTdz(1:nx,1:ny,1:nz)/), 4, x, y, z(1:nz))

  $if($MPI)
  !Ensure that all processes finish before attempting to write additional files.
  !Otherwise it may flood the system with too many I/O request and crash the process.
  call mpi_barrier (comm, ierr)
  $endif


  !Eshwan: write salinity
  !////////////////////////////////////////////
  !/// WRITE SALINITY                       ///
  !////////////////////////////////////////////
  $if ($BINARY)
  call string_splice (fname, output_path // 'output/binary_sal.', jt_total,'.dat')  !Eshwan
  $else
  call string_splice(fname, output_path // 'output/sal.',jt_total,'.dat')  !Eshwan
  $endif

  $if ($MPI)
  call string_concat (fname, '.c', coord)
  $endif 

  var_list = ' "x", "y", "z", "sal", "dSdx", "dSdy", "dSdz" '
  nvars = 7

  call write_tecplot_header_ND(fname,'rewind',nvars, (/Nx+1,Ny+1,Nz/), &
       trim(adjustl(var_list)), numtostr(coord, 6), 2, real(total_time,4))
  call write_real_data_3D(fname, 'append', 'formatted', 7, nx, ny, nz, &
       (/sal(1:nx,1:ny,1:nz), &
         dSdx(1:nx,1:ny,1:nz), &
         dSdy(1:nx,1:ny,1:nz), &
         dSdz(1:nx,1:ny,1:nz)/), 4, x, y, z(1:nz))

  $if($MPI)
  !Ensure that all processes finish before attempting to write additional files.
  !Otherwise it may flood the system with too many I/O request and crash the process.
  call mpi_barrier (comm, ierr)
  $endif



  !////////////////////////////////////////////
  !/// WRITE STANDARD DEVIATION OF SCALARS  ///
  !////////////////////////////////////////////
  $if ($BINARY)
  call string_splice (fname, output_path // 'output/binary_sigma_zplane.', jt_total,'.dat')  !Eshwan
  $else
  call string_splice(fname, output_path // 'output/sigma_zplane.',jt_total,'.dat')  !Eshwan
  $endif

  $if ($MPI)
  call string_concat (fname, '.c', coord)
  $endif 

  var_list = ' "z", "sigma_theta", "sigma_sal" '
  nvars = 3

  call write_tecplot_header_ND(fname,'rewind',nvars, (/ Nz /), &
       trim(adjustl(var_list)), numtostr(coord, 6), 2, real(total_time,4))
  call write_real_data_1D(fname, 'append', 'formatted', 2, nz, &
       (/sigma_theta(1:nz), &
         sigma_sal(1:nz)/), & 
         0, z(1:nz))

  $if($MPI)
  !Ensure that all processes finish before attempting to write additional files.
  !Otherwise it may flood the system with too many I/O request and crash the process.
  call mpi_barrier (comm, ierr)
  $endif



  !Eshwan: write divergence of shear stress
  !////////////////////////////////////////////
  !/// WRITE DIVERGENCE OF SHEAR STRESS     ///
  !////////////////////////////////////////////
  $if ($BINARY)
  call string_splice (fname, output_path // 'output/binary_divt.', jt_total,'.dat')  !Eshwan
  $else
  call string_splice(fname, output_path // 'output/divt.',jt_total,'.dat')  !Eshwan
  $endif

  $if ($MPI)
  call string_concat (fname, '.c', coord)
  $endif 

  var_list = ' "x", "y", "z", "divtx", "divty", "divtz" '
  nvars = 6

  call write_tecplot_header_ND(fname,'rewind',nvars, (/Nx+1,Ny+1,Nz/), &
       trim(adjustl(var_list)), numtostr(coord, 6), 2, real(total_time,4))
  call write_real_data_3D(fname, 'append', 'formatted', 3, nx, ny, nz, &
       (/divtx_w(1:nx,1:ny,1:nz), &
         divty_w(1:nx,1:ny,1:nz), &
         divtz(1:nx,1:ny,1:nz) /), 4, x, y, zw(1:nz))

  $if($MPI)
  !Ensure that all processes finish before attempting to write additional files.
  !Otherwise it may flood the system with too many I/O request and crash the process.
  call mpi_barrier (comm, ierr)
  $endif


  !Eshwan: write buoyancy
  !////////////////////////////////////////////
  !/// WRITE BUOYANCY                       ///
  !////////////////////////////////////////////
  $if ($BINARY)
  call string_splice (fname, output_path // 'output/binary_buoy.', jt_total,'.dat')  !Eshwan
  $else
  call string_splice(fname, output_path // 'output/buoy.',jt_total,'.dat')  !Eshwan
  $endif

  $if ($MPI)
  call string_concat (fname, '.c', coord)
  $endif 

  var_list = ' "x", "y", "z", "buoyancy" '
  nvars = 4

  call write_tecplot_header_ND(fname,'rewind',nvars, (/Nx+1,Ny+1,Nz/), &
       trim(adjustl(var_list)), numtostr(coord, 6), 2, real(total_time,4))
  call write_real_data_3D(fname, 'append', 'formatted', 1, nx, ny, nz, &
       (/buoyancy(1:nx,1:ny,1:nz)  /), 4, x, y, zw(1:nz))

  $if($MPI)
  !Ensure that all processes finish before attempting to write additional files.
  !Otherwise it may flood the system with too many I/O request and crash the process.
  call mpi_barrier (comm, ierr)
  $endif



  !Eshwan: subgrid-scale parameters
  !////////////////////////////////////////////
  !/// WRITE SGS PARAMETERS                 ///
  !////////////////////////////////////////////
  $if ($BINARY)
  call string_splice (fname, output_path // 'output/binary_sgs.', jt_total,'.dat')  !Eshwan
  $else
  call string_splice(fname, output_path // 'output/sgs.',jt_total,'.dat')  !Eshwan
  $endif

  $if ($MPI)
  call string_concat (fname, '.c', coord)
  $endif 

  var_list = ' "x", "y", "z", "cs2", "Nu_t"' 
  nvars = 5

  call write_tecplot_header_ND(fname,'rewind',nvars, (/Nx+1,Ny+1,Nz/), &
       trim(adjustl(var_list)), numtostr(coord, 6), 2, real(total_time,4))
  call write_real_data_3D(fname, 'append', 'formatted', 2, nx, ny, nz, &
       (/Cs_opt2(1:nx,1:ny,1:nz), &
         Nu_t(1:nx,1:ny,1:nz)/), &
          4, x, y, z(1:nz))

  $if($MPI)
  !Ensure that all processes finish before attempting to write additional files.
  !Otherwise it may flood the system with too many I/O request and crash the process.
  call mpi_barrier (comm, ierr)
  $endif
  


  !Eshwan: subgrid-scale parameters
  !////////////////////////////////////////////
  !/// WRITE SGS PARAMETERS                 ///
  !////////////////////////////////////////////
  !Write Fsub parameters
  $if ($BINARY)
  call string_splice (fname, output_path // 'output/binary_F.', jt_total,'.dat')  !Eshwan
  $else
  call string_splice(fname, output_path // 'output/F.',jt_total,'.dat')  !Eshwan
  $endif

  $if ($MPI)
  call string_concat (fname, '.c', coord)
  $endif 

  var_list = ' "x", "y", "z", "Tn", "Beta", "F_LM", "F_MM", "F_QN", "F_NN", "cs2_clips"' 
  nvars = 10

  call write_tecplot_header_ND(fname,'rewind',nvars, (/Nx+1,Ny+1,Nz/), &
       '"x", "y", "z", "Tn", "Beta", "F_LM", "F_MM", "F_QN", "F_NN", "cs2_clips"', &
        numtostr(coord, 6), 2, real(total_time,4))
  call write_real_data_3D(fname, 'append', 'formatted', 7, nx, ny, nz, &
       (/Tn_all(1:nx,1:ny,1:nz), &
         Beta(1:nx,1:ny,1:nz), &
         F_LM(1:nx,1:ny,1:nz), &
         F_MM(1:nx,1:ny,1:nz), &
         F_QN(1:nx,1:ny,1:nz), &
         F_NN(1:nx,1:ny,1:nz), &
         cs2_clips(1:nx,1:ny,1:nz)/), &
          4, x, y, z(1:nz))
  
  $if($MPI)
  !Ensure that all processes finish before attempting to write additional files.
  !Otherwise it may flood the system with too many I/O request and crash the process.
  call mpi_barrier (comm, ierr)
  $endif


  ! Write scalar sgs parameters
  $if ($BINARY)
  call string_splice (fname, output_path // 'output/binary_ssgs.', jt_total,'.dat')  !Eshwan
  $else
  call string_splice(fname, output_path // 'output/ssgs.',jt_total,'.dat')  !Eshwan
  $endif

  $if ($MPI)
  call string_concat (fname, '.c', coord)
  $endif 

  var_list = ' "x", "y", "z", "ds2_t", "kappa_tt", "ds2_s", "kappa_ts"' 
  nvars = 7

  call write_tecplot_header_ND(fname,'rewind',nvars, (/Nx+1,Ny+1,Nz/), &
       '"x", "y", "z", "ds2_t", "kappa_tt", "ds2_s", "kappa_ts"', &
        numtostr(coord, 6), 2, real(total_time,4))
  call write_real_data_3D(fname, 'append', 'formatted', 4, nx, ny, nz, &
       (/Ds_opt2_t(1:nx,1:ny,1:nz), &
         kappa_tt(1:nx,1:ny,1:nz), &
         Ds_opt2_s(1:nx,1:ny,1:nz), &
         kappa_ts(1:nx,1:ny,1:nz)/), &
          4, x, y, z(1:nz))

  $if($MPI)
  !Ensure that all processes finish before attempting to write additional files.
  !Otherwise it may flood the system with too many I/O request and crash the process.
  call mpi_barrier (comm, ierr)
  $endif



  !Write Isub_t parameters
  $if ($BINARY)
  call string_splice (fname, output_path // 'output/binary_I_t.', jt_total,'.dat')  !Eshwan
  $else
  call string_splice(fname, output_path // 'output/I_t.',jt_total,'.dat')  !Eshwan
  $endif

  $if ($MPI)
  call string_concat (fname, '.c', coord)
  $endif 

  var_list = ' "x", "y", "z", "s_Tn_t", "s_Beta_t", "I_LM_t", "I_MM_t", "I_QN_t", "I_NN_t" "ds2_clips_t"' 
  nvars = 10

  call write_tecplot_header_ND(fname,'rewind',nvars, (/Nx+1,Ny+1,Nz/), &
       '"x", "y", "z", "s_Tn_t", "s_Beta_t", "I_LM_t", "I_MM_t", "I_QN_t", "I_NN_t", "ds2_clips_t"', &
        numtostr(coord, 6), 2, real(total_time,4))
  call write_real_data_3D(fname, 'append', 'formatted', 7, nx, ny, nz, &
       (/s_Tn_all_t(1:nx,1:ny,1:nz), &
         s_Beta_t(1:nx,1:ny,1:nz), &
         I_LM_t(1:nx,1:ny,1:nz), &
         I_MM_t(1:nx,1:ny,1:nz), &
         I_QN_t(1:nx,1:ny,1:nz), &
         I_NN_t(1:nx,1:ny,1:nz), &
         ds2_clips_t(1:nx,1:ny,1:nz)/), &
          4, x, y, z(1:nz))

  $if($MPI)
  !Ensure that all processes finish before attempting to write additional files.
  !Otherwise it may flood the system with too many I/O request and crash the process.
  call mpi_barrier (comm, ierr)
  $endif



  !Write Isub_s parameters
  $if ($BINARY)
  call string_splice (fname, output_path // 'output/binary_I_s.', jt_total,'.dat')  !Eshwan
  $else
  call string_splice(fname, output_path // 'output/I_s.',jt_total,'.dat')  !Eshwan
  $endif

  $if ($MPI)
  call string_concat (fname, '.c', coord)
  $endif 

  var_list = ' "x", "y", "z", "s_Tn_s", "s_Beta_s", "I_LM_s", "I_MM_s", "I_QN_s", "I_NN_s" "ds2_clips_s"' 
  nvars = 10

  call write_tecplot_header_ND(fname,'rewind',nvars, (/Nx+1,Ny+1,Nz/), &
       '"x", "y", "z", "s_Tn_s", "s_Beta_s", "I_LM_s", "I_MM_s", "I_QN_s", "I_NN_s", "ds2_clips_s"', &
        numtostr(coord, 6), 2, real(total_time,4))
  call write_real_data_3D(fname, 'append', 'formatted', 7, nx, ny, nz, &
       (/s_Tn_all_s(1:nx,1:ny,1:nz), &
         s_Beta_s(1:nx,1:ny,1:nz), &
         I_LM_s(1:nx,1:ny,1:nz), &
         I_MM_s(1:nx,1:ny,1:nz), &
         I_QN_s(1:nx,1:ny,1:nz), &
         I_NN_s(1:nx,1:ny,1:nz), &
         ds2_clips_s(1:nx,1:ny,1:nz)/), &
          4, x, y, z(1:nz))

  $if($MPI)
  !Ensure that all processes finish before attempting to write additional files.
  !Otherwise it may flood the system with too many I/O request and crash the process.
  call mpi_barrier (comm, ierr)
  $endif



  !////////////////////////////////////////////
  !/// WRITE VELOCITY                       ///
  !////////////////////////////////////////////

  $if( $BINARY )
  call string_splice( fname, output_path // 'output/binary_vel.', jt_total,'.dat')  !Eshwan
  $else
  call string_splice( fname, output_path // 'output/vel.', jt_total, '.dat')  !Eshwan
  $endif

  $if ($MPI)
  call string_concat( fname, '.c', coord, '.tec' )
  $endif
 
  $if($BINARY)

  ! RICHARD
  open(unit=13,file=fname,form='`unformatted',convert='big_endian', access='direct',recl=nx*ny*nz*rprec)
  write(13,rec=1) u(:nx,:ny,1:nz)
  write(13,rec=2) v(:nx,:ny,1:nz)
  write(13,rec=3) w_uv(:nx,:ny,1:nz)
  close(13)

  $else

    $if($LVLSET)
    var_list = '"x", "y", "z", "u", "v", "w", "phi"'
    nvars = 7
    call write_tecplot_header_ND(fname, 'rewind', nvars, (/ Nx+1, Ny+1, Nz/), &
         trim(adjustl(var_list)), numtostr(coord, 6), 2, real(total_time,4))
    call write_real_data_3D(fname, 'append', 'formatted', 4, nx, ny, nz, &
         (/ u(1:nx,1:ny,1:nz), &
         v(1:nx,1:ny,1:nz), &
         w_uv(1:nx,1:ny,1:nz), &
         phi(1:nx,1:ny,1:nz)/), & 
         4, x, y, z(1:nz))
    $else   
    var_list = '"x", "y", "z", "u", "v", "w"'
    nvars = 6
    call write_tecplot_header_ND(fname, 'rewind', nvars, (/ Nx+1, Ny+1, Nz/), &
         trim(adjustl(var_list)), numtostr(coord, 6), 2, real(total_time,4))
    call write_real_data_3D(fname, 'append', 'formatted', 3, nx,ny,nz, &
         (/ u_w(1:nx,1:ny,1:nz), &
         v_w(1:nx,1:ny,1:nz), &
         w(1:nx,1:ny,1:nz) /), &
         4, x, y, zw(1:nz))    !Eshwan
    $endif
  
  $endif
  

  $if($MPI)
    ! Ensure that all processes finish before attempting to write 
    ! additional files. Otherwise it may flood the system with 
    ! too many I/O requests and crash the process 
    call mpi_barrier( comm, ierr )
  $endif

  !  Output instantaneous force field 
  $if( not BINARY )
  $if($LVLSET)
    !////////////////////////////////////////////
    !/// WRITE FORCES                         ///
    !////////////////////////////////////////////

    ! Compute the total forces 
    call force_tot()

    !  Open file which to write global data
    call string_splice( fname, output_path // 'output/force.', jt_total, '.dat')  !Eshwan


    $if ($MPI)
    call string_concat( fname, '.c', coord )
    $endif

    var_list = '"x", "y", "z", "f<sub>x</sub>", "f<sub>y</sub>", "f<sub>z</sub>", "phi"'
    nvars = 7

    call write_tecplot_header_ND(fname, 'rewind', nvars, (/ Nx+1, Ny+1, Nz/), &
         trim(adjustl(var_list)), numtostr(coord, 6), 2, real(total_time,4))
    call write_real_data_3D(fname, 'append', 'formatted', 4, nx, ny,nz, &
                            (/ fx_tot, fy_tot, fz_tot, &
                            phi(1:nx,1:ny,1:nz) /), 4, x, y, z(1:nz))

    deallocate(fx_tot, fy_tot, fz_tot)

    $if($MPI)
      ! Ensure that all processes finish before attempting to write 
      ! additional files. Otherwise it may flood the system with 
      ! too many I/O requests and crash the process 
      call mpi_barrier( comm, ierr )
    $endif

  $endif
  $endif

  $if($DEBUG)

    !////////////////////////////////////////////
    !/// WRITE VELOCITY DIVERGENCE            ///
    !////////////////////////////////////////////


    !  Output divergence of velocity field
    allocate(divvel(nx,ny,nz))
    divvel = dudx(1:nx,1:ny,1:nz) + dvdy(1:nx,1:ny,1:nz) + dwdz(1:nx,1:ny,1:nz)

    !  Open file which to write global data
    $if( $BINARY )
    call string_splice( fname, output_path // 'output/binary_divvel.', jt_total,'.dat')  !Eshwan
    $else
    call string_splice( fname, output_path // 'output/divvel.', jt_total, '.dat')  !Eshwan
    $endif

    $if ($MPI)
      call string_concat( fname, '.c', coord )
    $endif

    $if($BINARY)
    ! RICHARD
    open(unit=13,file=fname2,form='unformatted',convert='big_endian', access='direct',recl=nx*ny*nz*rprec)
    write(13,rec=1) divvel(:nx,:ny,1:nz)
    close(13)
     
    $else
  
      $if($LVLSET)
      var_list = '"x", "y", "z", "divvel", "phi"'
      nvars = 5
      call write_tecplot_header_ND(fname, 'rewind', nvars, (/ Nx+1, Ny+1, Nz/), &
           trim(adjustl(var_list)), numtostr(coord, 6), 2, real(total_time,4))
      call write_real_data_3D(fname, 'append', 'formatted', 2, nx, ny,nz, &
           (/ divvel, phi(1:nx,1:ny,1:nz) /), 4, x, y, z(1:nz))
      $else   
      var_list = '"x", "y", "z", "divvel"'
      nvars = 4
      call write_tecplot_header_ND(fname, 'rewind', nvars, (/ Nx+1, Ny+1, Nz/), &
           trim(adjustl(var_list)), numtostr(coord, 6), 2, real(total_time,4))
      call write_real_data_3D(fname, 'append', 'formatted', 1, nx, ny,nz, &
           (/ divvel /), 4, x, y, z(1:nz))
      $endif

    $endif

    deallocate(divvel)

    $if($MPI)
      ! Ensure that all processes finish before attempting to write 
      ! additional files. Otherwise it may flood the system with 
      ! too many I/O requests and crash the process 
      call mpi_barrier( comm, ierr )
    $endif

    !////////////////////////////////////////////
    !/// WRITE PRESSURE                       ///
    !////////////////////////////////////////////

    !  Open file which to write global data

    $if($BINARY)
    call string_splice( fname, output_path // 'output/binary_pressure.', jt_total,'.dat')  !Eshwan
    $else
    call string_splice( fname, output_path // 'output/pressure.', jt_total, '.dat')  !Eshwan
    $endif

    $if ($MPI)
    call string_concat( fname, '.c', coord )
    $endif

    call pressure_sync()

    $if($BINARY)

    ! RICHARD
    open(unit=13,file=fname2,form='unformatted',convert='big_endian', access='direct',recl=nx*ny*nz*rprec)
    write(13,rec=1) p(:nx,:ny,1:nz)
    write(13,rec=2) dpdx(:nx,:ny,1:nz)
    write(13,rec=3) dpdy(:nx,:ny,1:nz)
    write(13,rec=4) interp_to_uv_grid(dpdz(:nx,:ny,1:nz),1)
    close(13)

    $else
    
      $if($LVLSET)

      var_list = '"x", "y", "z", "p", "dpdx", "dpdy", "dpdz", "phi"'
      nvars = 8
      call write_tecplot_header_ND(fname, 'rewind', nvars, (/ Nx+1, Ny+1, Nz/), &
           trim(adjustl(var_list)), numtostr(coord,6), 2, real(total_time,4))
      call write_real_data_3D(fname, 'append', 'formatted', 5, nx, ny,nz, &
           (/ p(1:nx,1:ny,1:nz), &
           dpdx(1:nx,1:ny,1:nz), &
           dpdy(1:nx,1:ny,1:nz), &
           interp_to_uv_grid(dpdz(1:nx,1:ny,1:nz),1), &
           phi(1:nx,1:ny,1:nz) /), &
           4, x, y, z(1:nz))

      $else

      var_list = '"x", "y", "z", "p", "dpdx", "dpdy", "dpdz"'
      nvars = 7
      call write_tecplot_header_ND(fname, 'rewind', nvars, (/ Nx+1, Ny+1, Nz/), &
           trim(adjustl(var_list)), numtostr(coord,6), 2, real(total_time,4))
      call write_real_data_3D(fname, 'append', 'formatted', 4, nx, ny,nz, &
        (/ p(1:nx,1:ny,1:nz), &
        dpdx(1:nx,1:ny,1:nz), &
        dpdy(1:nx,1:ny,1:nz), &
        interp_to_uv_grid(dpdz(1:nx,1:ny,1:nz),1) /), &
        4, x, y, z(1:nz))

      $endif

    $endif

    $if($MPI)
      ! Ensure that all processes finish before attempting to write 
      ! additional files. Otherwise it may flood the system with 
      ! too many I/O requests and crash the process 
      call mpi_barrier( comm, ierr )
    $endif
  
    !////////////////////////////////////////////
    !/// WRITE RHS                            ///
    !////////////////////////////////////////////

    !  Open file which to write global data
    $if($BINARY)
    call string_splice( fname, output_path // 'output/binary_RHS.', jt_total,'.dat')  !Eshwan
    $else
    call string_splice( fname, output_path // 'output/RHS.', jt_total, '.dat')  !Eshwan
    $endif

    $if ($MPI)
    call string_concat( fname, '.c', coord )
    $endif
  
    call RHS_sync()

    $if($BINARY)
    ! RICHARD
    open(unit=13,file=fname,form='unformatted',convert='big_endian', access='direct',recl=nx*ny*nz*rprec)
    write(13,rec=1) RHSx(:nx,:ny,1:nz)
    write(13,rec=2) RHSy(:nx,:ny,1:nz)
    write(13,rec=3) dpdy(:nx,:ny,1:nz)
    write(13,rec=4) interp_to_uv_grid(RHSz(:nx,:ny,1:nz),1)
    close(13)
    
    $else

      $if($LVLSET)
      
      var_list = '"x", "y", "z", "RHSx", "RHSy", "RHSz", "phi"'
      nvars = 7
      call write_tecplot_header_ND(fname, 'rewind', nvars, (/ Nx+1, Ny+1, Nz/), &
           trim(adjustl(var_list)), numtostr(coord,6), 2, real(total_time,4))
      call write_real_data_3D(fname, 'append', 'formatted', 4, nx, ny, nz, &
           (/ RHSx(1:nx,1:ny,1:nz), &
           RHSy(1:nx,1:ny,1:nz), &
           interp_to_uv_grid(RHSz(1:nx,1:ny,1:nz),1), &
           phi(1:nx,1:ny,1:nz) /), & 
           4, x, y, z(1:nz))

      $else

      var_list = '"x", "y", "z", "RHSx", "RHSy", "RHSz"'
      nvars = 6
      call write_tecplot_header_ND(fname, 'rewind', nvars, (/ Nx+1, Ny+1, Nz/), &
           trim(adjustl(var_list)), numtostr(coord,6), 2, real(total_time,4))
      call write_real_data_3D(fname, 'append', 'formatted', 3, nx, ny,nz, &
           (/ RHSx(1:nx,1:ny,1:nz), &
           RHSy(1:nx,1:ny,1:nz), &
           interp_to_uv_grid(RHSz(1:nx,1:ny,1:nz),1) /), &
           4, x, y, z(1:nz))

      $endif

    $endif

    $if($MPI)
      ! Ensure that all processes finish before attempting to write 
      ! additional files. Otherwise it may flood the system with 
      ! too many I/O requests and crash the process 
      call mpi_barrier( comm, ierr )
    $endif

  $endif

!  Write instantaneous x-plane values
elseif(itype==3) then

  allocate(ui(1,ny,nz), vi(1,ny,nz), wi(1,ny,nz))

  $if($LVLSET)
  call force_tot()
  $endif

!  Loop over all xplane locations
  do i=1,xplane_nloc

     call string_splice( fname, path // 'output/vel.x-', xplane_loc(i), '.', jt_total, '.dat')
  
     $if ($MPI)
       call string_concat( fname, '.c', coord )
     $endif

    call write_tecplot_header_ND(fname, 'rewind', 6, (/ 1, Ny+1, Nz /), &
      '"x", "y", "z", "u", "v", "w"', numtostr(coord,6), 2, real(total_time,4))  
  
    do k=1,nz
      do j=1,ny

        ui(1,j,k) = linear_interp(u(xplane(i) % istart,j,k), &
             u(xplane(i) % istart+1,j,k), dx, xplane(i) % ldiff)
        vi(1,j,k) = linear_interp(v(xplane(i) % istart,j,k), &
             v(xplane(i) % istart+1,j,k), dx, xplane(i) % ldiff)
        wi(1,j,k) = linear_interp(w_uv(xplane(i) % istart,j,k), &
             w_uv(xplane(i) % istart+1,j,k), dx, &
             xplane(i) % ldiff)
      enddo
    enddo

    call write_real_data_3D(fname, 'append', 'formatted', 3, 1, ny, nz, &
      (/ ui, vi, wi /), 2, (/ xplane_loc(i) /), y, z(1:nz))     

    $if($LVLSET)

    call string_splice( fname, path // 'output/force.x-', xplane_loc(i), '.', jt_total, '.dat')

    $if ($MPI)
    call string_concat( fname, '.c', coord )
    $endif

    call write_tecplot_header_ND(fname, 'rewind', 6, (/ 1, Ny+1, Nz/), &
      '"x", "y", "z", "f<sub>x</sub>", "f<sub>y</sub>", "f<sub>z</sub>"', &
      numtostr(coord,6), 2, real(total_time,4))

    !  Sum both induced forces, f{x,y,z}, and applied forces, f{x,y,z}a
    do k=1,nz
      do j=1,ny

        ui(1,j,k) = linear_interp(fx_tot(xplane(i) % istart,j,k), &
             fx_tot(xplane(i) % istart+1,j,k), dx, xplane(i) % ldiff)

        vi(1,j,k) = linear_interp(fy_tot(xplane(i) % istart,j,k), &
             fy_tot(xplane(i) % istart+1,j,k), dx, xplane(i) % ldiff)

        wi(1,j,k) = linear_interp(fz_tot(xplane(i) % istart,j,k), &
             fz_tot(xplane(i) % istart+1,j,k), dx, xplane(i) % ldiff)

      enddo
    enddo


    call write_real_data_3D(fname, 'append', 'formatted', 3, 1, ny, nz, &
      (/ ui, vi, wi /), 2, (/ xplane_loc(i) /), y, z(1:nz))




    $endif

    !$if($OUTPUT_EXTRA)

    !////////////////////////////////////////////
    !/// WRITE LDSM                           ///
    !////////////////////////////////////////////

    if( sgs_model == 4 ) then

      allocate(F_LM_s(nx,nz),F_MM_s(nx,nz))
      allocate(beta_s(nx,nz),Cs_opt2_s(nx,nz))
      allocate(Nu_t_s(nx,nz))

      do k=1,Nz
        do j=1,Ny
!
          F_LM_s(j,k) = linear_interp(F_LM_uv(xplane(i) % istart,j,k), &
               F_LM_uv(xplane(i) % istart+1,j,k), &
               dx, xplane(i) % ldiff)
          F_MM_s(j,k) = linear_interp(F_MM_uv(xplane(i) % istart,j,k), &
               F_MM_uv(xplane(i) % istart+1,j,k), &
               dx, xplane(i) % ldiff)
          beta_s(j,k) = linear_interp(beta_uv(xplane(i) % istart,j,k), &
               beta_uv(xplane(i) % istart+1,j,k), &
               dx, xplane(i) % ldiff)
          Cs_opt2_s(j,k) = linear_interp(Cs_opt2_uv(xplane(i) % istart,j,k), &
               Cs_opt2_uv(xplane(i) % istart+1,j,k), &
               dx, xplane(i) % ldiff)
          Nu_t_s(j,k) = linear_interp(Nu_t_uv_interp(xplane(i) % istart,j,k), &
               Nu_t_uv_interp(xplane(i) % istart+1,j,k), &
               dx, xplane(i) % ldiff)

        enddo
      enddo

      call string_splice( fname, path // 'output/ldsm.x-', xplane_loc(i), '.', jt_total, '.dat')

      $if ($MPI)
      call string_concat( fname, '.c', coord )
      $endif

      var_list = '"x", "y", "z", "F<sub>LM</sub>", "F<sub>MM</sub>"'
      var_list = trim(adjustl(var_list)) // ', "<greek>b</greek>", "Cs<sup>2</sup>"'
      var_list = trim(adjustl(var_list)) // ', "<greek>n</greek><sub>T</sub>"'

      call write_tecplot_header_ND(fname, 'rewind', 8, (/ 1, Ny+1, Nz/), &
           trim(adjustl(var_list)), numtostr(coord,6), 2, real(total_time,4)) 

      call write_real_data_3D(fname, 'append', 'formatted', 5, 1,ny,nz, &
                              (/ F_LM_s, F_MM_s, beta_s, Cs_opt2_s, Nu_t_s /), &
                              2, (/ xplane_loc(i) /), y, z(1:nz)) 

      deallocate(F_LM_s,F_MM_s,beta_s,Cs_opt2_s,Nu_t_s)

    elseif( sgs_model == 5 ) then

      allocate(F_LM_s(nx,nz),F_MM_s(nx,nz))
      allocate(F_QN_s(nx,nz),F_NN_s(nx,nz))
      allocate(beta_s(nx,nz),Cs_opt2_s(nx,nz))
      allocate(Nu_t_s(nx,nz))

      do k=1,Nz
        do j=1,Ny
!
          F_LM_s(j,k) = linear_interp(F_LM_uv(xplane(i) % istart,j,k), &
               F_LM_uv(xplane(i) % istart+1,j,k), &
               dx, xplane(i) % ldiff)
          F_MM_s(j,k) = linear_interp(F_MM_uv(xplane(i) % istart,j,k), &
               F_MM_uv(xplane(i) % istart+1,j,k), &
               dx, xplane(i) % ldiff)
          F_QN_s(j,k) = linear_interp(F_QN_uv(xplane(i) % istart,j,k), &
               F_QN_uv(xplane(i) % istart+1,j,k), &
               dx, xplane(i) % ldiff)  
          F_NN_s(j,k) = linear_interp(F_NN_uv(xplane(i) % istart,j,k), &
               F_NN_uv(xplane(i) % istart+1,j,k), &
               dx, xplane(i) % ldiff)                                         
          beta_s(j,k) = linear_interp(beta_uv(xplane(i) % istart,j,k), &
               beta_uv(xplane(i) % istart+1,j,k), &
               dx, xplane(i) % ldiff)
          Cs_opt2_s(j,k) = linear_interp(Cs_opt2_uv(xplane(i) % istart,j,k), &
               Cs_opt2_uv(xplane(i) % istart+1,j,k), &
               dx, xplane(i) % ldiff)
          Nu_t_s(j,k) = linear_interp(Nu_t_uv_interp(xplane(i) % istart,j,k), &
               Nu_t_uv_interp(xplane(i) % istart+1,j,k), &
               dx, xplane(i) % ldiff)

        enddo
      enddo

      call string_splice( fname , path // 'output/ldsm.x-', xplane_loc(i), '.', jt_total, '.dat')

      $if ($MPI)
      call string_concat( fname, '.c', coord )
      $endif

      var_list = '"x", "y", "z", "F<sub>LM</sub>", "F<sub>MM</sub>"'
      var_list = trim(adjustl(var_list)) // ', "F<sub>QN</sub>", "F<sub>NN</sub>"'
      var_list = trim(adjustl(var_list)) // ', "<greek>b</greek>", "Cs<sup>2</sup>"'
      var_list = trim(adjustl(var_list)) // ', "<greek>n</greek><sub>T</sub>"'

      call write_tecplot_header_ND(fname, 'rewind', 10, (/ 1, Ny+1, Nz/), &
        trim(adjustl(var_list)), numtostr(coord,6), 2, real(total_time,4)) 

      call write_real_data_3D(fname, 'append', 'formatted', 7, 1,ny,nz, &
                             (/ F_LM_s, F_MM_s, F_QN_s, F_NN_s, beta_s, Cs_opt2_s, Nu_t_s /), &
                             2, (/ xplane_loc(i) /), y, z(1:nz)) 

      deallocate(F_LM_s,F_MM_s,F_QN_s,F_NN_s,beta_s,Cs_opt2_s,Nu_t_s)

    endif   

    !$endif    
    
  enddo   
  
  deallocate(ui,vi,wi)

  $if($LVLSET)
  deallocate( fx_tot, fy_tot, fz_tot )
  $endif

!  Write instantaneous y-plane values
elseif(itype==4) then
  
  allocate(ui(nx,1,nz), vi(nx,1,nz), wi(nx,1,nz))

  $if($LVLSET)
  call force_tot()
  $endif
  
!  Loop over all yplane locations
  do j=1,yplane_nloc

    call string_splice( fname, path // 'output/vel.y-', yplane_loc(j), '.', jt_total, '.dat')
  
    $if ($MPI)
    call string_concat( fname, '.c', coord )
    $endif

    call write_tecplot_header_ND(fname, 'rewind', 6, (/ Nx+1, 1, Nz/), &
      '"x", "y", "z", "u", "v", "w"', numtostr(coord,6), 2, real(total_time,4)) 

    do k=1,nz
      do i=1,nx

        ui(i,1,k) = linear_interp(u(i,yplane(j) % istart,k), &
             u(i,yplane(j) % istart+1,k), dy, yplane(j) % ldiff)
        vi(i,1,k) = linear_interp(v(i,yplane(j) % istart,k), &
             v(i,yplane(j) % istart+1,k), dy, yplane(j) % ldiff)
        wi(i,1,k) = linear_interp(w_uv(i,yplane(j) % istart,k), &
             w_uv(i,yplane(j) % istart+1,k), dy, &
             yplane(j) % ldiff)
        
      enddo
    enddo
    
    call write_real_data_3D(fname, 'append', 'formatted', 3, nx,1,nz, &
      (/ ui, vi, wi /), 1, x, (/ yplane_loc(j) /), z(1:nz))    
  
    $if($LVLSET)

    call string_splice( fname, path // 'output/force.y-', yplane_loc(j), '.', jt_total, '.dat')

    $if ($MPI)
    call string_concat( fname, '.c', coord )
    $endif

    call write_tecplot_header_ND(fname, 'rewind', 6, (/ Nx+1, 1, Nz/), &
      '"x", "y", "z", "fx", "fy", "fz"', numtostr(coord,6), 2, real(total_time,4))  
  
    do k=1,nz
      do i=1,nx

        ui(i,1,k) = linear_interp(fx_tot(i,yplane(j) % istart,k), &
             fx_tot(i,yplane(j) % istart+1,k), dy, yplane(j) % ldiff)
        vi(i,1,k) = linear_interp(fy_tot(i,yplane(j) % istart,k), &
             fy_tot(i,yplane(j) % istart+1,k), dy, yplane(j) % ldiff)
        wi(i,1,k) = linear_interp(fz_tot(i,yplane(j) % istart,k), &
             fz_tot(i,yplane(j) % istart+1,k), dy, yplane(j) % ldiff)

      enddo
    enddo
    
    call write_real_data_3D(fname, 'append', 'formatted', 3, nx,1,nz, &
      (/ ui, vi, wi /), 1, x, (/ yplane_loc(j) /), z(1:nz))       
    
    $endif

    !$if($OUTPUT_EXTRA)

    !////////////////////////////////////////////
    !/// WRITE LDSM                           ///
    !////////////////////////////////////////////

    if( sgs_model == 4 ) then

      allocate(F_LM_s(nx,nz),F_MM_s(nx,nz))
      allocate(beta_s(nx,nz),Cs_opt2_s(nx,nz))
      allocate(Nu_t_s(nx,nz))

      do k=1,Nz
        do i=1,Nx
!
          F_LM_s(i,k) = linear_interp(F_LM_uv(i,yplane(j) % istart,k), &
               F_LM_uv(i,yplane(j) % istart+1,k), &
               dy, yplane(j) % ldiff)
          F_MM_s(i,k) = linear_interp(F_MM_uv(i,yplane(j) % istart,k), &
               F_MM_uv(i,yplane(j) % istart+1,k), &
               dy, yplane(j) % ldiff)
          beta_s(i,k) = linear_interp(beta_uv(i,yplane(j) % istart,k), &
               beta_uv(i,yplane(j) % istart+1,k), &
               dy, yplane(j) % ldiff)
          Cs_opt2_s(i,k) = linear_interp(Cs_opt2_uv(i,yplane(j) % istart,k), &
               Cs_opt2_uv(i,yplane(j) % istart+1,k), &
               dy, yplane(j) % ldiff)
          Nu_t_s(i,k) = linear_interp(Nu_t_uv_interp(i,yplane(j) % istart,k), &
               Nu_t_uv_interp(i,yplane(j) % istart+1,k), &
               dy, yplane(j) % ldiff)

        enddo
      enddo

      call string_splice( fname, path // 'output/ldsm.y-', yplane_loc(j), '.', jt_total, '.dat')

      $if ($MPI)
      call string_concat( fname, '.c', coord )
      $endif

      var_list = '"x", "y", "z", "F<sub>LM</sub>", "F<sub>MM</sub>"'
      var_list = trim(adjustl(var_list)) // ', "<greek>b</greek>", "Cs<sup>2</sup>"'
      var_list = trim(adjustl(var_list)) // ', "<greek>n</greek><sub>T</sub>"'

      call write_tecplot_header_ND(fname, 'rewind', 8, (/ Nx+1, 1, Nz/), &
        trim(adjustl(var_list)), numtostr(coord,6), 2, real(total_time,4)) 

      call write_real_data_3D(fname, 'append', 'formatted', 5, nx,1,nz, &
                              (/ F_LM_s, F_MM_s, beta_s, Cs_opt2_s, Nu_t_s /), &
                              1, x, (/ yplane_loc(j) /), z(1:nz)) 

      deallocate(F_LM_s,F_MM_s,beta_s,Cs_opt2_s,Nu_t_s)

    elseif( sgs_model == 5 ) then

      allocate(F_LM_s(nx,nz),F_MM_s(nx,nz))
      allocate(F_QN_s(nx,nz),F_NN_s(nx,nz))
      allocate(beta_s(nx,nz),Cs_opt2_s(nx,nz))
      allocate(Nu_t_s(nx,nz))

      do k=1,Nz
        do i=1,Nx

          F_LM_s(i,k) = linear_interp(F_LM_uv(i,yplane(j) % istart,k), &
               F_LM_uv(i,yplane(j) % istart+1,k), &
               dy, yplane(j) % ldiff)
          F_MM_s(i,k) = linear_interp(F_MM_uv(i,yplane(j) % istart,k), &
               F_MM_uv(i,yplane(j) % istart+1,k), &
               dy, yplane(j) % ldiff)
          F_QN_s(i,k) = linear_interp(F_QN_uv(i,yplane(j) % istart,k), &
               F_QN_uv(i,yplane(j) % istart+1,k), &
               dy, yplane(j) % ldiff)
          F_NN_s(i,k) = linear_interp(F_NN_uv(i,yplane(j) % istart,k), &
               F_NN_uv(i,yplane(j) % istart+1,k), &
               dy, yplane(j) % ldiff)                                      
          beta_s(i,k) = linear_interp(beta_uv(i,yplane(j) % istart,k), &
               beta_uv(i,yplane(j) % istart+1,k), &
               dy, yplane(j) % ldiff)
          Cs_opt2_s(i,k) = linear_interp(Cs_opt2_uv(i,yplane(j) % istart,k), &
               Cs_opt2_uv(i,yplane(j) % istart+1,k), &
               dy, yplane(j) % ldiff)
          Nu_t_s(i,k) = linear_interp(Nu_t_uv_interp(i,yplane(j) % istart,k), &
               Nu_t_uv_interp(i,yplane(j) % istart+1,k), &
               dy, yplane(j) % ldiff)

        enddo
      enddo

      call string_splice( fname, path // 'output/ldsm.y-', yplane_loc(j), '.', jt_total,'.dat')

      $if ($MPI)
      call string_concat( fname, '.c', coord )
      $endif

      var_list = '"x", "y", "z", "F<sub>LM</sub>", "F<sub>MM</sub>"'
      var_list = trim(adjustl(var_list)) // ',  "F<sub>QN</sub>", "F<sub>NN</sub>"'
      var_list = trim(adjustl(var_list)) // ', "<greek>b</greek>", "Cs<sup>2</sup>"'
      var_list = trim(adjustl(var_list)) // ', "<greek>n</greek><sub>T</sub>"'

      call write_tecplot_header_ND(fname, 'rewind', 10, (/ Nx+1, 1, Nz/), &
        trim(adjustl(var_list)), numtostr(coord,6), 2, real(total_time,4)) 

      call write_real_data_3D(fname, 'append', 'formatted', 7, nx,1,nz, &
                              (/ F_LM_s, F_MM_s, F_QN_s, F_NN_s, beta_s, Cs_opt2_s, Nu_t_s /), &
                              1, x, (/ yplane_loc(j) /), z(1:nz)) 

      deallocate(F_LM_s,F_MM_s,F_QN_s,F_NN_s,beta_s,Cs_opt2_s,Nu_t_s)

    endif   

    !$endif

  enddo  

  deallocate(ui,vi,wi)

  $if($LVLSET)
  deallocate(fx_tot, fy_tot, fz_tot)
  $endif
  
!  Write instantaneous z-plane values
elseif(itype==5) then

  allocate(ui(nx,ny,1), vi(nx,ny,1), wi(nx,ny,1))

  $if($LVLSET)
  call force_tot()
  $endif

!  Loop over all zplane locations
  do k=1,zplane_nloc
  
    $if ($MPI)
    if(zplane(k) % coord == coord) then
    $endif

    $if ($BINARY)
    call string_splice( fname, path // 'output/binary_vel.z-', zplane_loc(k), '.', jt_total, '.dat')
    $else
    call string_splice( fname, path // 'output/vel.z-', zplane_loc(k), '.', jt_total, '.dat')    
    $endif
    
    do j=1,Ny
      do i=1,Nx

        ui(i,j,1) = linear_interp(u(i,j,zplane(k) % istart), &
             u(i,j,zplane(k) % istart+1), &
             dz, zplane(k) % ldiff)
        vi(i,j,1) = linear_interp(v(i,j,zplane(k) % istart), &
             v(i,j,zplane(k) % istart+1), &
             dz, zplane(k) % ldiff)
        wi(i,j,1) = linear_interp(w_uv(i,j,zplane(k) % istart), &
             w_uv(i,j,zplane(k) % istart+1), &
             dz, zplane(k) % ldiff)

      enddo
    enddo
    
    $if ($BINARY)
    open(unit=13,file=fname,form='unformatted',convert='big_endian', access='direct',recl=nx*ny*1*rprec)
    write(13,rec=1) ui(1:nx,1:ny,1)
    write(13,rec=2) vi(1:nx,1:ny,1)
    write(13,rec=3) wi(1:nx,1:ny,1)
    close(13)
    
    $else
    call write_tecplot_header_ND(fname, 'rewind', 6, (/ Nx+1, Ny+1, 1/), &
      '"x", "y", "z", "u", "v", "w"', numtostr(coord,6), 2, real(total_time,4)) 
    call write_real_data_3D(fname, 'append', 'formatted', 3, nx,ny,1, &
    (/ ui, vi, wi /), 4, x, y, (/ zplane_loc(k) /))   
    $endif
    
    $if($LVLSET)

    call string_splice( fname, path // 'output/force.z-', zplane_loc(k), '.', jt_total, '.dat')

    call write_tecplot_header_ND(fname, 'rewind', 6, (/ Nx+1, Ny+1, 1/), &
      '"x", "y", "z", "f<sub>x</sub>", "f<sub>y</sub>", "f<sub>z</sub>"', &
      numtostr(coord,6), 2, real(total_time,4))

    do j=1,Ny
      do i=1,Nx

        ui(i,j,1) = linear_interp(fx_tot(i,j,zplane(k) % istart), &
             fx_tot(i,j,zplane(k) % istart+1), &
             dz, zplane(k) % ldiff) 
        vi(i,j,1) = linear_interp(fy_tot(i,j,zplane(k) % istart), &
             fy_tot(i,j,zplane(k) % istart+1), &
             dz, zplane(k) % ldiff) 
        wi(i,j,1) = linear_interp(fz_tot(i,j,zplane(k) % istart), &
             fz_tot(i,j,zplane(k) % istart+1), &
             dz, zplane(k) % ldiff)

      enddo
    enddo
 
    call write_real_data_3D(fname, 'append', 'formatted', 3, nx,ny,1, &
         (/ ui, vi, wi /), 4, x, y, (/ zplane_loc(k) /) )      
    
    $endif

    !$if($OUTPUT_EXTRA)

    !////////////////////////////////////////////
    !/// WRITE LDSM                           ///
    !////////////////////////////////////////////

    if( sgs_model == 4 ) then

      allocate(F_LM_s(nx,ny),F_MM_s(nx,ny))
      allocate(beta_s(nx,ny),Cs_opt2_s(nx,ny))
      allocate(Nu_t_s(nx,ny))

      do j=1,Ny
        do i=1,Nx
!
          F_LM_s(i,j) = linear_interp(F_LM_uv(i,j,zplane(k) % istart), &
               F_LM_uv(i,j,zplane(k) % istart+1), &
               dz, zplane(k) % ldiff)
          F_MM_s(i,j) = linear_interp(F_MM_uv(i,j,zplane(k) % istart), &
               F_MM_uv(i,j,zplane(k) % istart+1), &
               dz, zplane(k) % ldiff)
          beta_s(i,j) = linear_interp(beta_uv(i,j,zplane(k) % istart), &
               beta_uv(i,j,zplane(k) % istart+1), &
               dz, zplane(k) % ldiff)
          Cs_opt2_s(i,j) = linear_interp(Cs_opt2_uv(i,j,zplane(k) % istart), &
               Cs_opt2_uv(i,j,zplane(k) % istart+1), &
               dz, zplane(k) % ldiff)
          Nu_t_s(i,j) = linear_interp(Nu_t_uv_interp(i,j,zplane(k) % istart), &
               Nu_t_uv_interp(i,j,zplane(k) % istart+1), &
               dz, zplane(k) % ldiff)

        enddo
      enddo

      call string_splice( fname, path // 'output/ldsm.z-', zplane_loc(k), '.', jt_total, '.dat')

      var_list = '"x", "y", "z", "F<sub>LM</sub>", "F<sub>MM</sub>"'
      var_list = trim(adjustl(var_list)) // ', "<greek>b</greek>", "Cs<sup>2</sup>"'
      var_list = trim(adjustl(var_list)) // ', "<greek>n</greek><sub>T</sub>"'

      call write_tecplot_header_ND(fname, 'rewind', 8, (/ Nx+1, Ny+1, 1/), &
        trim(adjustl(var_list)), numtostr(coord,6), 2, real(total_time,4))

      call write_real_data_3D(fname, 'append', 'formatted', 5, nx,ny,1, &
                              (/ F_LM_s, F_MM_s, beta_s, Cs_opt2_s, Nu_t_s /), &
                              4, x, y, (/ zplane_loc(k) /) )

      deallocate(F_LM_s,F_MM_s,beta_s,Cs_opt2_s,Nu_t_s)

     elseif( sgs_model == 5 ) then
 
      allocate(F_LM_s(nx,ny),F_MM_s(nx,ny))
      allocate(F_QN_s(nx,ny),F_NN_s(nx,ny))
      allocate(beta_s(nx,ny),Cs_opt2_s(nx,ny))
      allocate(Nu_t_s(nx,ny))

      do j=1,Ny
        do i=1,Nx
!
          F_LM_s(i,j) = linear_interp(F_LM_uv(i,j,zplane(k) % istart), &
               F_LM_uv(i,j,zplane(k) % istart+1), &
               dz, zplane(k) % ldiff) 
          F_MM_s(i,j) = linear_interp(F_MM_uv(i,j,zplane(k) % istart), &
               F_MM_uv(i,j,zplane(k) % istart+1), &
               dz, zplane(k) % ldiff) 
          F_QN_s(i,j) = linear_interp(F_QN_uv(i,j,zplane(k) % istart), &
               F_QN_uv(i,j,zplane(k) % istart+1), &
               dz, zplane(k) % ldiff)
          F_NN_s(i,j) = linear_interp(F_NN_uv(i,j,zplane(k) % istart), &
               F_NN_uv(i,j,zplane(k) % istart+1), &
               dz, zplane(k) % ldiff)
          beta_s(i,j) = linear_interp(beta_uv(i,j,zplane(k) % istart), &
               beta_uv(i,j,zplane(k) % istart+1), &
               dz, zplane(k) % ldiff)            
          Cs_opt2_s(i,j) = linear_interp(Cs_opt2_uv(i,j,zplane(k) % istart), &
               Cs_opt2_uv(i,j,zplane(k) % istart+1), &
               dz, zplane(k) % ldiff)                                         
          Nu_t_s(i,j) = linear_interp(Nu_t_uv_interp(i,j,zplane(k) % istart), &
               Nu_t_uv_interp(i,j,zplane(k) % istart+1), &
               dz, zplane(k) % ldiff)

        enddo
      enddo      

      call string_splice( fname, path // 'output/ldsm.z-', zplane_loc(k), '.', jt_total, '.dat')

      var_list = '"x", "y", "z", "F<sub>LM</sub>", "F<sub>MM</sub>"'
      var_list = trim(adjustl(var_list)) // ', "F<sub>QN</sub>", "F<sub>NN</sub>"'
      var_list = trim(adjustl(var_list)) // ', "<greek>b</greek>", "Cs<sup>2</sup>"'
      var_list = trim(adjustl(var_list)) // ', "<greek>n</greek><sub>T</sub>"'

      call write_tecplot_header_ND(fname, 'rewind', 10, (/ Nx+1, Ny+1, 1/), &
        trim(adjustl(var_list)), numtostr(coord,6), 2, real(total_time,4))

      call write_real_data_3D(fname, 'append', 'formatted', 7, nx,ny,1, &
                              (/ F_LM_s, F_MM_s, F_QN_s, F_NN_s, beta_s, Cs_opt2_s, Nu_t_s /), &
                              4, x, y, (/ zplane_loc(k) /) )         

      deallocate(F_LM_s,F_MM_s,F_QN_s,F_NN_s,beta_s,Cs_opt2_s,Nu_t_s)

    endif
    !$endif

    $if ($MPI)
    endif
    $endif

  enddo  
  
  deallocate(ui,vi,wi)
  
  $if($LVLSET)
  deallocate(fx_tot, fy_tot, fz_tot)
  $endif

else
  write(*,*) 'Error: itype not specified properly to inst_write!'
  stop
endif

deallocate(w_uv)
nullify(x,y,z,zw)

$if($LVLSET or $DEBUG)
contains
$endif

$if($LVLSET)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine force_tot()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
$if($MPI)
use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWN
$endif
implicit none

! Zero bogus values
fx(:,:,nz) = 0._rprec
fy(:,:,nz) = 0._rprec
fz(:,:,nz) = 0._rprec

!  Sum both the induced and applied forces
allocate(fx_tot(nx,ny,nz), fy_tot(nx,ny,nz), fz_tot(nx,ny,nz))

! Richard: Might not be necessary to do this as the function only seems to be called when LVLSET is activated
$if($TURBINES)
fx_tot = fxa(1:nx,1:ny,1:nz)
fy_tot = 0._rprec
fz_tot = 0._rprec
$elseif($LVLSET)
fx_tot = fx(1:nx,1:ny,1:nz)+fxa(1:nx,1:ny,1:nz)
fy_tot = fy(1:nx,1:ny,1:nz)+fya(1:nx,1:ny,1:nz)
fz_tot = fz(1:nx,1:ny,1:nz)+fza(1:nx,1:ny,1:nz)
$else
fx_tot = 0._rprec
fy_tot = 0._rprec
fz_tot = 0._rprec
$endif

$if($MPI)
!  Sync forces
call mpi_sync_real_array( fx_tot, 1, MPI_SYNC_DOWN )
call mpi_sync_real_array( fy_tot, 1, MPI_SYNC_DOWN )
call mpi_sync_real_array( fz_tot, 1, MPI_SYNC_DOWN )
$endif

! Put fz_tot on uv-grid
fz_tot(1:nx,1:ny,1:nz) = interp_to_uv_grid( fz_tot(1:nx,1:ny,1:nz), 1 )

return
end subroutine force_tot
$endif

$if($DEBUG)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine pressure_sync()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWN
use param, only : ld
implicit none

! Reset bogus values
p(:,:,nz) = p(:,:,nz-1)
dpdx(:,:,nz) = dpdx(:,:,nz-1)
dpdy(:,:,nz) = dpdy(:,:,nz-1)
dpdz(:,:,nz) = dpdz(:,:,nz-1)

$if($MPI)
!  Sync pressure
call mpi_sync_real_array( p, 0 , MPI_SYNC_DOWN )
call mpi_sync_real_array( dpdx, 1 , MPI_SYNC_DOWN )
call mpi_sync_real_array( dpdy, 1 , MPI_SYNC_DOWN )
call mpi_sync_real_array( dpdz, 1 , MPI_SYNC_DOWN )
$endif

return
end subroutine pressure_sync

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine RHS_sync()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use param, only : ld
use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWN
implicit none

! Reset bogus values
RHSx(:,:,nz) = RHSx(:,:,nz-1)
RHSy(:,:,nz) = RHSy(:,:,nz-1)
RHSz(:,:,nz) = RHSz(:,:,nz-1)

$if($MPI)
!  Sync RHS
call mpi_sync_real_array( RHSx, 0 , MPI_SYNC_DOWN )
call mpi_sync_real_array( RHSy, 0 , MPI_SYNC_DOWN )
call mpi_sync_real_array( RHSz, 0 , MPI_SYNC_DOWN )
$endif

return
end subroutine RHS_sync
$endif

end subroutine inst_write

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine checkpoint ()
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use param, only : nz, checkpoint_file, tavg_calc, spectra_calc
$if($MPI)
use param, only : coord
use param, only : comm, ierr
$endif
use sim_param, only : u, v, w, RHSx, RHSy, RHSz
use sim_param, only:  theta, sal, RHS_T, RHS_S  !Eshwan
use sgs_param, only : Cs_opt2, F_LM, F_MM, F_QN, F_NN
use sgs_param, only : Ds_opt2_t, I_LM_t, I_MM_t, I_QN_t, I_NN_t
use sgs_param, only : Ds_opt2_s, I_LM_s, I_MM_s, I_QN_s, I_NN_s
$if($DYN_TN)
use sgs_param, only: F_ee2, F_deedt2, ee_past
$endif
use param, only : jt_total, total_time, total_time_dim, dt
use param, only : use_cfl_dt, cfl
use cfl_util, only : get_max_cfl
use stat_defs, only : tavg_initialized, spectra_initialized
use string_util, only : string_concat
implicit none

character(64) :: fname
real(rprec) :: cfl_w

fname = checkpoint_file
$if ($MPI)
call string_concat( fname, '.c', coord )
$endif

!  Open vel.out (lun_default in io) for final output
$if ($WRITE_BIG_ENDIAN)
open(11,file=fname,form='unformatted', convert='big_endian', status='unknown', position='rewind')
$elseif ($WRITE_LITTLE_ENDIAN)
open(11,file=fname,form='unformatted', convert='little_endian', status='unknown', position='rewind')
$else
open(11,file=fname,form='unformatted', status='unknown', position='rewind')
$endif

write (11) u(:, :, 1:nz), v(:, :, 1:nz), w(:, :, 1:nz),      &
     RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),   &
     Cs_opt2(:,:,1:nz), F_LM(:,:,1:nz), F_MM(:,:,1:nz),      &
     F_QN(:,:,1:nz), F_NN(:,:,1:nz), theta(:,:,1:nz),        &
     RHS_T(:,:,1:nz), Ds_opt2_t(:,:,1:nz), I_LM_t(:,:,1:nz), &
     I_MM_t(:,:,1:nz), I_QN_t(:,:,1:nz), I_NN_t(:,:,1:nz),   &
     sal(:,:,1:nz), RHS_S(:,:,1:nz), Ds_opt2_s(:,:,1:nz),    &
     I_LM_s(:,:,1:nz), I_MM_s(:,:,1:nz), I_QN_s(:,:,1:nz),   &
     I_NN_s(:,:,1:nz)

! Close the file to ensure that the data is flushed and written to file
close(11)

$if($MPI)
call mpi_barrier( comm, ierr )
$endif

$if ($DYN_TN) 
! Write running average variables to file
  fname = path // 'dyn_tn.out'
  $if ($MPI)
  call string_concat( fname, '.c', coord)
  $endif
  
  $if ($WRITE_BIG_ENDIAN)
  open(unit=13,file=fname,form='unformatted',position='rewind',convert='big_endian')
  $elseif ($WRITE_LITTLE_ENDIAN)
  open(unit=13,file=fname,form='unformatted',position='rewind',convert='little_endian')
  $else
  open(unit=13,file=fname,form='unformatted',position='rewind')
  $endif

  write(13) F_ee2(:,:,1:nz), F_deedt2(:,:,1:nz), ee_past(:,:,1:nz)

  close(13)

  $if($MPI)
  call mpi_barrier( comm, ierr )
  $endif
  
$endif

! Checkpoint time averaging restart data
if( tavg_calc .and. tavg_initialized ) call tavg_checkpoint()
! Checkpoint spectra restart data
if( spectra_calc .and. spectra_initialized ) call spectra_checkpoint()

! Write time and current simulation state
! Set the current cfl to a temporary (write) value based whether CFL is
! specified or must be computed
if( use_cfl_dt ) then
   cfl_w = cfl
else
   cfl_w = get_max_cfl()
endif

!  Update total_time.dat after simulation
!if ((cumulative_time) .and. (lun == lun_default)) then
if (coord == 0) then
   !--only do this for true final output, not intermediate recording
   open (1, file=fcumulative_time)

   write(1, *) jt_total, total_time, total_time_dim, dt, cfl_w

   close(1)

end if
!end if

return
end subroutine checkpoint


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine output_final()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use stat_defs, only : tavg, point, tavg_initialized, spectra_initialized
use param, only : tavg_calc, point_calc, point_nloc, spectra_calc
implicit none

! Perform final checkpoing
call checkpoint()

!  Check if average quantities are to be recorded
if(tavg_calc .and. tavg_initialized ) call tavg_finalize()

!  Check if spectra is to be computed
if(spectra_calc .and. spectra_initialized ) call spectra_finalize()

return
end subroutine output_final

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine len_da_file(fname, lenrec, length)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!--finds number of records on existing direct-access unformatted file
!--taken from Clive Page's comp.lang.fortran posting (12/16/2003), 
!  under the topic counting number of records in a Fortran direct file
!--minor changes/renaming
!
use messages
implicit none
character (*), intent(in) :: fname  ! name of existing direct-access file
integer, intent(in)       :: lenrec ! record length (O/S dependent units)
integer, intent(out) :: length      ! number of records.
!
character(*), parameter :: sub_name = mod_name // '.len_dat_file'
character (1) :: cdummy
integer :: lunit, nlo, nhi, mid, kode
logical :: exists, open
!
! find a free unit on which to open the file
!
do lunit = 99, 1, -1
  !--units to skip (compiler dependent)
  select case (lunit)
    case (5:6)
      !--do nothing
    case default
      inquire(unit=lunit, exist=exists, opened=open)
      if(exists .and. .not. open) exit
  end select
end do
open(unit=lunit, file=fname, access="direct", recl=lenrec, iostat=kode)
if(kode /= 0) then
   call mesg( sub_name, 'error in len_da_file: ' // trim(fname) // ' does not exist' )
   return
end if
!
! expansion phase
!
mid = 1
do
  read(lunit, rec=mid, iostat=kode) cdummy
  if(kode /= 0) exit
  mid = 2 * mid
end do
!
! length is between mid/2 and mid, do binary search to refine
!
nlo = mid/2
nhi = mid
do while(nhi - nlo > 1)
  mid = (nlo + nhi) / 2
  read(lunit, rec=mid, iostat=kode) cdummy
  if(kode == 0) then
     nlo = mid
  else
     nhi = mid
  end if
end do
length = nlo
close(unit=lunit)
return
end subroutine len_da_file

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine output_init ()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  This subroutine allocates the memory for arrays used for statistical
!  calculations
use param, only : path
use param, only : L_x,L_y,L_z,dx,dy,dz,nx,ny,nz,nsteps,coord,nproc,lbz,lh
use param, only : point_calc, point_nloc, point_loc
use param, only : xplane_calc, xplane_nloc, xplane_loc
use param, only : yplane_calc, yplane_nloc, yplane_loc
use param, only : zplane_calc, zplane_nloc, zplane_loc
use param, only : spectra_calc, spectra_nloc, spectra_loc
use param, only : tavg_calc
use grid_defs, only : grid
use functions, only : cell_indx
use stat_defs, only : point, xplane, yplane, zplane
use stat_defs, only : tavg, tavg_zplane, spectra, theta_avg, theta_avg_zplane, & 
                      sal_avg, sal_avg_zplane, ustar_avg, t_flux_avg, tstar_avg, & 
                      sal_flux_avg, sstar_avg, sigma_avg_zplane
!$if($OUTPUT_EXTRA)
use stat_defs, only : tavg_sgs,tavg_scalar_sgs
!$endif
use stat_defs, only : type_set
implicit none

include 'tecryte.h'

!character(120) :: cx,cy,cz
character(120) :: var_list, fname
integer :: i,j,k

logical :: exst

real(rprec), pointer, dimension(:) :: x,y,z

nullify(x,y,z)

x => grid % x
y => grid % y
z => grid % z

if( tavg_calc ) then

  allocate(tavg(nx,ny,lbz:nz))
  allocate(tavg_zplane(lbz:nz))
  allocate(theta_avg(nx,ny,lbz:nz))    !Eshwan
  allocate(theta_avg_zplane(lbz:nz))   !Eshwan
  allocate(sal_avg(nx,ny,lbz:nz))      !Eshwan
  allocate(sal_avg_zplane(lbz:nz))     !Eshwan
  allocate(ustar_avg(nx,ny))
  allocate(t_flux_avg(nx,ny))
  allocate(tstar_avg(nx,ny))
  allocate(sal_flux_avg(nx,ny))
  allocate(sstar_avg(nx,ny))  
  allocate(sigma_avg_zplane(lbz:nz))   !Eshwan

  !$if($OUTPUT_EXTRA)
  allocate(tavg_sgs(nx,ny,lbz:nz))
  allocate(tavg_scalar_sgs(nx,ny,lbz:nz))
  !$endif

  ! Initialize the derived types tavg and tavg_zplane  
  do k=1,Nz
    do j=1, Ny
      do i=1, Nx
        call type_set ( tavg(i,j,k), 0._rprec )
        call type_set ( theta_avg(i,j,k), 0._rprec )           !Eshwan
        call type_set ( sal_avg(i,j,k), 0._rprec )           !Eshwan
        !$if($OUTPUT_EXTRA)
        call type_set( tavg_sgs(i,j,k), 0._rprec )
        call type_set( tavg_scalar_sgs(i,j,k), 0._rprec )
        !$endif        
      enddo
    enddo
    
    call type_set( tavg_zplane(k), 0._rprec )
    call type_set( theta_avg_zplane(k), 0._rprec )              !Eshwan
    call type_set( sal_avg_zplane(k), 0._rprec )              !Eshwan
    call type_set( sigma_avg_zplane(k), 0._rprec )            !Eshwan
  enddo
 
  do j=1,Ny
     do i=1,Nx
        ustar_avg(i,j) = 0._rprec
        t_flux_avg(i,j) = 0._rprec
        tstar_avg(i,j) = 0._rprec
        sal_flux_avg(i,j) = 0._rprec
        sstar_avg(i,j) = 0._rprec
     end do
  end do
   
endif

if(spectra_calc) then

  allocate(spectra(spectra_nloc))
!  Initialize 
  spectra(:) % istart = -1
  spectra(:) % ldiff = 0._rprec
  spectra(:) % coord = -1

!  Compute istart and ldiff
  do k=1,spectra_nloc

    !  Initialize sub-arrays
    allocate( spectra(k)%power(lh) )
    spectra(k) % power = 0._rprec

    $if ($MPI)
    if(spectra_loc(k) >= z(1) .and. spectra_loc(k) < z(nz)) then
      spectra(k) % coord = coord
      spectra(k) % istart = cell_indx('k',dz,spectra_loc(k))
      spectra(k) % ldiff = spectra_loc(k) - z(spectra(k) % istart)
    endif
    $else
    spectra(k) % coord = 0
    spectra(k) % istart = cell_indx('k',dz,spectra_loc(k))
    spectra(k) % ldiff = spectra_loc(k) - z(spectra(k) % istart)
    $endif

  enddo

endif

! Initialize information for x-planar stats/data
if(xplane_calc) then

  allocate(xplane(xplane_nloc))

  xplane(:) % istart = -1
  xplane(:) % ldiff = 0.
  
!  Compute istart and ldiff
  do i=1,xplane_nloc
    xplane(i) % istart = cell_indx('i', dx, xplane_loc(i))
    xplane(i) % ldiff = xplane_loc(i) - x(xplane(i) % istart)   
  enddo
    
endif

! Initialize information for y-planar stats/data
if(yplane_calc) then

  allocate(yplane(yplane_nloc))

  yplane(:) % istart = -1
  yplane(:) % ldiff = 0.
  
!  Compute istart and ldiff
  do j=1,yplane_nloc
    yplane(j) % istart = cell_indx('j', dy, yplane_loc(j))
    yplane(j) % ldiff = yplane_loc(j) - y(yplane(j) % istart)
  enddo
    
endif

! Initialize information for z-planar stats/data
if(zplane_calc) then

  allocate(zplane(zplane_nloc))

!  Initialize 
  zplane(:) % istart = -1
  zplane(:) % ldiff = 0. 
  zplane(:) % coord=-1 
  
!  Compute istart and ldiff
  do k=1,zplane_nloc

    $if ($MPI)
    if(zplane_loc(k) >= z(1) .and. zplane_loc(k) < z(nz)) then
      zplane(k) % coord = coord
      zplane(k) % istart = cell_indx('k',dz,zplane_loc(k))
      zplane(k) % ldiff = zplane_loc(k) - z(zplane(k) % istart)
    endif
    $else
    zplane(k) % coord = 0
    zplane(k) % istart = cell_indx('k',dz,zplane_loc(k))
    zplane(k) % ldiff = zplane_loc(k) - z(zplane(k) % istart)
    $endif

  enddo  
  
endif


!  Open files for instantaneous writing
if(point_calc) then

  allocate(point(point_nloc))

  !  Intialize the coord values (-1 shouldn't be used as coord so initialize to this)
  point % coord=-1
  point % fid = -1

  do i=1,point_nloc
    !  Find the processor in which this point lives
    $if ($MPI)
    if(point_loc(i)%xyz(3) >= z(1) .and. point_loc(i)%xyz(3) < z(nz)) then
    $endif

    point(i) % coord = coord
  
    point(i) % istart = cell_indx('i',dx,point_loc(i)%xyz(1))
    point(i) % jstart = cell_indx('j',dy,point_loc(i)%xyz(2))
    point(i) % kstart = cell_indx('k',dz,point_loc(i)%xyz(3))
    
    point(i) % xdiff = point_loc(i)%xyz(1) - x(point(i) % istart)
    point(i) % ydiff = point_loc(i)%xyz(2) - y(point(i) % jstart)
    point(i) % zdiff = point_loc(i)%xyz(3) - z(point(i) % kstart)
    
    call string_splice(fname, path // 'output/vel.x-', point_loc(i)%xyz(1), &
         '.y-', point_loc(i)%xyz(2), &
         '.z-', point_loc(i)%xyz(3),'.dat')
    
    !  Add tecplot header if file does not exist
    inquire (file=fname, exist=exst)
    if (exst) then
       point(i) % fid = open_file( fname, 'append', 'formatted' )
    else

       point(i) % fid = open_file( fname, 'rewind', 'formatted' )
       var_list = '"t", "u", "v", "w"'
       ! Compilation error
       call write_tecplot_header_xyline( point(i) % fid, var_list )

    endif
    
    $if($MPI)
    endif
    $endif

  enddo
endif

nullify(x,y,z)

return
end subroutine output_init

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine tavg_init()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  Load tavg.out files
use param, only : path
use param, only : coord, dt, Nx, Ny, Nz
use messages
use stat_defs, only : tavg, tavg_total_time, tavg_dt, tavg_initialized
use stat_defs, only : operator(.MUL.)
!$if($OUTPUT_EXTRA)
use stat_defs, only : tavg_sgs, tavg_total_time_sgs
!$endif
use param, only : tavg_nstart, tavg_nend
implicit none

character (*), parameter :: sub_name = mod_name // '.tavg_init'
character (*), parameter :: ftavg_in = path // 'tavg.out'
!$if($OUTPUT_EXTRA)
character (*), parameter :: ftavg_sgs_in = path // 'tavg_sgs.out'
!$endif
$if ($MPI)
character (*), parameter :: MPI_suffix = '.c'
$endif
character (128) :: fname

logical :: opn, exst
!integer :: i,j,k

inquire (unit=1, opened=opn)
if (opn) call error (sub_name, 'unit 1 already open')

fname = ftavg_in
$if ($MPI)
call string_concat( fname, MPI_suffix, coord )
$endif

inquire (file=fname, exist=exst)
if (.not. exst) then
  !  Nothing to read in
  if (coord == 0) then
    write(*,*) ' '
    write(*,*)'No previous time averaged data - starting from scratch.'
  endif

    tavg_total_time = 0.0_rprec
    ! note: tavg was already initialized to zero in output_init routine
 
else

    $if ($READ_BIG_ENDIAN)
    open (1, file=fname, action='read', position='rewind',  &
      form='unformatted', convert='big_endian')
    $elseif ($READ_LITTLE_ENDIAN)
    open (1, file=fname, action='read', position='rewind',  &
      form='unformatted', convert='little_endian')  
    $else
    open (1, file=fname, action='read', position='rewind',  &
      form='unformatted')
    $endif

    read (1) tavg_total_time
    read (1) tavg

    close(1)

endif

!------
!$if($OUTPUT_EXTRA)

    fname = ftavg_sgs_in
    $if ($MPI)
    call string_concat( fname, MPI_suffix, coord )
    $endif

    inquire (file=fname, exist=exst)
    if (.not. exst) then
        !  Nothing to read in
        if(coord == 0) then
            write(*,*) ' '
            write(*,*)'No previous time averaged data (sgs) - starting from scratch.'
        endif

        tavg_total_time_sgs = 0._rprec  
    else
        $if ($READ_BIG_ENDIAN)
        open (1, file=fname, action='read', position='rewind',  &
          form='unformatted', convert='big_endian')
        $elseif ($READ_LITTLE_ENDIAN)
        open (1, file=fname, action='read', position='rewind',  &
          form='unformatted', convert='little_endian')  
        $else
        open (1, file=fname, action='read', position='rewind',  &
          form='unformatted')
        $endif

        read (1) tavg_total_time_sgs
        read (1) tavg_sgs

        close(1)    
    endif
    
!$endif

! Initialize tavg_dt
tavg_dt = 0.0_rprec

! Set global switch that tavg as been initialized
tavg_initialized = .true. 

return
end subroutine tavg_init

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine tavg_compute()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  This subroutine collects the stats for each flow 
!  variable quantity
use types, only : rprec
use param, only : tavg_nskip
use stat_defs, only : tavg, tavg_zplane, tavg_total_time, tavg_dt, theta_avg, & 
                      theta_avg_zplane, t_flux_avg, tstar_avg, &
                      ustar_avg, sal_avg, sal_avg_zplane, sal_flux_avg, &
                      sstar_avg, sigma_avg_zplane   !Eshwan
!$if($OUTPUT_EXTRA)
use param, only : sgs_model
use stat_defs, only : tavg_sgs, tavg_total_time_sgs, tavg_scalar_sgs
use sgs_param
!$endif
use param, only : nx,ny,nz,dt,lbz,jzmin,jzmax
use sim_param, only : u,v,w, dudz, dvdz, txx, txy, tyy, txz, tyz, tzz,u_avg, &
                      theta, t_flux, tstar, sal_flux, sstar, ustar, &
                      sal, dTdx, dTdy, dTdz, dSdx, dSdy, dSdz, &
                      sigma_theta, sigma_sal, &
                      divtx, divty, divtz, buoyancy   !Eshwan
$if($TURBINES)
use sim_param, only : fxa
$elseif($LVLSET)
use sim_param, only : fx, fy, fz, fxa, fya, fza
$endif

use functions, only : interp_to_w_grid, interp_to_uv_grid, &
                      interp_eddy_to_uv_grid

implicit none

integer :: i,j,k

real(rprec) :: u_p, v_p, w_p
real(rprec), allocatable, dimension(:,:,:) :: u_w, v_w ! Eshwan
real(rprec), allocatable, dimension(:,:,:) :: w_uv !Eshwan
real(rprec), allocatable, dimension(:,:,:) :: theta_w
real(rprec), allocatable, dimension(:,:,:) :: sal_w
real(rprec), allocatable, dimension(:,:,:) :: divtx_w, divty_w

allocate(u_w(nx,ny,lbz:nz),v_w(nx,ny,lbz:nz))
allocate(w_uv(nx,ny,lbz:nz))
allocate(theta_w(nx,ny,lbz:nz))
allocate(sal_w(nx,ny,lbz:nz))
allocate(divtx_w(nx,ny,lbz:nz),divty_w(nx,ny,lbz:nz))

!  Interpolate velocities to w-grid
u_w(1:nx,1:ny,lbz:nz) = interp_to_w_grid( u(1:nx,1:ny,lbz:nz), lbz )
v_w(1:nx,1:ny,lbz:nz) = interp_to_w_grid( v(1:nx,1:ny,lbz:nz), lbz )

!Interpolate velocity w to uv-grid
w_uv(1:nx,1:ny,lbz:nz) = interp_to_uv_grid( w(1:nx,1:ny,lbz:nz), lbz)

!Interpolate temperature to w-grid
theta_w(1:nx,1:ny,lbz:nz) = interp_to_w_grid( theta(1:nx,1:ny,lbz:nz), lbz )

!Interpolate salinity to w-grid
sal_w(1:nx,1:ny,lbz:nz) = interp_to_w_grid( sal(1:nx,1:ny,lbz:nz), lbz )

!Interpolate divtx and divty to w-grid
divtx_w(1:nx,1:ny,lbz:nz) = interp_to_w_grid( divtx(1:nx,1:ny,lbz:nz), lbz)
divty_w(1:nx,1:ny,lbz:nz) = interp_to_w_grid( divty(1:nx,1:ny,lbz:nz), lbz)

$if($MPI)
k=0
! === w-grid variables === 
sigma_avg_zplane(k)%sigma_theta = sigma_avg_zplane(k)%sigma_theta + sigma_theta(k) * tavg_dt
sigma_avg_zplane(k)%sigma_sal = sigma_avg_zplane(k)%sigma_sal + sigma_sal(k) * tavg_dt 

do j=1,ny
   do i=1,nx

      u_p = u_w(i,j,k)
      v_p = v_w(i,j,k) 
      w_p = w(i,j,k)

      ! === w-grid variables === 
      tavg(i,j,k)%u = tavg(i,j,k)%u + u_p * tavg_dt                    
      tavg(i,j,k)%v = tavg(i,j,k)%v + v_p * tavg_dt                         
      tavg(i,j,k)%w = tavg(i,j,k)%w + w_p * tavg_dt

      tavg(i,j,k)%u2 = tavg(i,j,k)%u2 + u_p * u_p * tavg_dt
      tavg(i,j,k)%v2 = tavg(i,j,k)%v2 + v_p * v_p * tavg_dt
      tavg(i,j,k)%w2 = tavg(i,j,k)%w2 + w_p * w_p * tavg_dt
      tavg(i,j,k)%uv = tavg(i,j,k)%uv + u_p * v_p * tavg_dt
      tavg(i,j,k)%uw = tavg(i,j,k)%uw + u_p * w_p * tavg_dt
      tavg(i,j,k)%vw = tavg(i,j,k)%vw + v_p * w_p * tavg_dt

      tavg(i,j,k)%dudz = tavg(i,j,k)%dudz + dudz(i,j,k) * tavg_dt
      tavg(i,j,k)%dvdz = tavg(i,j,k)%dvdz + dvdz(i,j,k) * tavg_dt
      tavg(i,j,k)%dudz2 = tavg(i,j,k)%dudz2 + dudz(i,j,k) * dudz(i,j,k) * tavg_dt !Eshwan

      tavg(i,j,k)%divtx = tavg(i,j,k)%divtx + divtx_w(i,j,k) *tavg_dt
      tavg(i,j,k)%divty = tavg(i,j,k)%divty + divty_w(i,j,k) *tavg_dt
      tavg(i,j,k)%divtz = tavg(i,j,k)%divtz + divtz(i,j,k) *tavg_dt
      tavg(i,j,k)%buoyancy = tavg(i,j,k)%buoyancy + buoyancy(i,j,k) *tavg_dt

      theta_avg(i,j,k)%theta_w = theta_avg(i,j,k)%theta_w + theta_w(i,j,k) * tavg_dt !Eshwan
      theta_avg(i,j,k)%wtheta = theta_avg(i,j,k)%wtheta + w(i,j,k) * theta_w(i,j,k) * tavg_dt !Eshwan
      theta_avg(i,j,k)%theta2 = theta_avg(i,j,k)%theta2 + theta_w(i,j,k) * theta_w(i,j,k) * tavg_dt !Eshwan
      theta_avg(i,j,k)%dTdz = theta_avg(i,j,k)%dTdz + dTdz(i,j,k) * tavg_dt !Eshwan
      theta_avg(i,j,k)%sgs_tflux = theta_avg(i,j,k)%sgs_tflux + kappa_tt(i,j,k) * dTdz(i,j,k) *tavg_dt !Eshwan
     
      sal_avg(i,j,k)%sal_w = sal_avg(i,j,k)%sal_w + sal_w(i,j,k) * tavg_dt !Eshwan
      sal_avg(i,j,k)%wsal = sal_avg(i,j,k)%wsal + w(i,j,k) * sal_w(i,j,k) * tavg_dt !Eshwan
      sal_avg(i,j,k)%sal2 = sal_avg(i,j,k)%sal2 + sal_w(i,j,k) * sal_w(i,j,k) * tavg_dt !Eshwan
      sal_avg(i,j,k)%dSdz = sal_avg(i,j,k)%dSdz + dSdz(i,j,k) * tavg_dt !Eshwan
      sal_avg(i,j,k)%sgs_sflux = sal_avg(i,j,k)%sgs_sflux + kappa_ts(i,j,k) * dSdz(i,j,k) * tavg_dt !Eshwan  
      
      ! === uv-grid variables ===
      tavg(i,j,k)%u_uv = tavg(i,j,k)%u_uv + u(i,j,k) * tavg_dt !Eshwan
      tavg(i,j,k)%v_uv = tavg(i,j,k)%v_uv + v(i,j,k) * tavg_dt !Eshwan
      tavg(i,j,k)%w_uv = tavg(i,j,k)%w_uv + w_uv(i,j,k) * tavg_dt !Eshwan
      
      tavg(i,j,k)%txx = tavg(i,j,k)%txx + txx(i,j,k) * tavg_dt
      tavg(i,j,k)%tyy = tavg(i,j,k)%tyy + tyy(i,j,k) * tavg_dt
      tavg(i,j,k)%tzz = tavg(i,j,k)%tzz + tzz(i,j,k) * tavg_dt
      tavg(i,j,k)%txy = tavg(i,j,k)%txy + txy(i,j,k) * tavg_dt

      theta_avg(i,j,k)%theta = theta_avg(i,j,k)%theta + theta(i,j,k) * tavg_dt !Eshwan
      theta_avg(i,j,k)%dTdx = theta_avg(i,j,k)%dTdx + dTdx(i,j,k) * tavg_dt !Eshwan
      theta_avg(i,j,k)%dTdy = theta_avg(i,j,k)%dTdy + dTdy(i,j,k) * tavg_dt !Eshwan

      sal_avg(i,j,k)%sal = sal_avg(i,j,k)%sal + sal(i,j,k) * tavg_dt !Eshwan
      sal_avg(i,j,k)%dSdx = sal_avg(i,j,k)%dSdx + dSdx(i,j,k) * tavg_dt !Eshwan
      sal_avg(i,j,k)%dSdy = sal_avg(i,j,k)%dSdy + dSdy(i,j,k) * tavg_dt !Eshwan

      ! === w-grid variables === 
      tavg(i,j,k)%txz = tavg(i,j,k)%txz + txz(i,j,k) * tavg_dt
      tavg(i,j,k)%tyz = tavg(i,j,k)%tyz + tyz(i,j,k) * tavg_dt

      tavg(i,j,k)%cs_opt2 = tavg(i,j,k)%cs_opt2 + Cs_opt2(i,j,k) * tavg_dt
      tavg(i,j,k)%Nu_t = tavg(i,j,k)%Nu_t + Nu_t(i,j,k) * tavg_dt !Eshwan
   enddo
enddo
$endif

do k=1,jzmax  
  ! === w-grid variables === 
  sigma_avg_zplane(k)%sigma_theta = sigma_avg_zplane(k)%sigma_theta + sigma_theta(k) * tavg_dt
  sigma_avg_zplane(k)%sigma_sal = sigma_avg_zplane(k)%sigma_sal + sigma_sal(k) * tavg_dt 
  do j=1,ny
    do i=1,nx
   
      u_p = u_w(i,j,k)
      v_p = v_w(i,j,k) 
      w_p = w(i,j,k)
    
      ! === w-grid variables === 
      tavg(i,j,k)%u = tavg(i,j,k)%u + u_p * tavg_dt                    
      tavg(i,j,k)%v = tavg(i,j,k)%v + v_p * tavg_dt                         
      tavg(i,j,k)%w = tavg(i,j,k)%w + w_p * tavg_dt

      tavg(i,j,k)%u2 = tavg(i,j,k)%u2 + u_p * u_p * tavg_dt
      tavg(i,j,k)%v2 = tavg(i,j,k)%v2 + v_p * v_p * tavg_dt
      tavg(i,j,k)%w2 = tavg(i,j,k)%w2 + w_p * w_p * tavg_dt
      tavg(i,j,k)%uv = tavg(i,j,k)%uv + u_p * v_p * tavg_dt
      tavg(i,j,k)%uw = tavg(i,j,k)%uw + u_p * w_p * tavg_dt
      tavg(i,j,k)%vw = tavg(i,j,k)%vw + v_p * w_p * tavg_dt

      tavg(i,j,k)%dudz = tavg(i,j,k)%dudz + dudz(i,j,k) * tavg_dt
      tavg(i,j,k)%dvdz = tavg(i,j,k)%dvdz + dvdz(i,j,k) * tavg_dt
      tavg(i,j,k)%dudz2 = tavg(i,j,k)%dudz2 + dudz(i,j,k) * dudz(i,j,k) * tavg_dt              !Eshwan

      tavg(i,j,k)%divtx = tavg(i,j,k)%divtx + divtx_w(i,j,k) *tavg_dt
      tavg(i,j,k)%divty = tavg(i,j,k)%divty + divty_w(i,j,k) *tavg_dt
      tavg(i,j,k)%divtz = tavg(i,j,k)%divtz + divtz(i,j,k) *tavg_dt
      tavg(i,j,k)%buoyancy = tavg(i,j,k)%buoyancy + buoyancy(i,j,k) *tavg_dt

      theta_avg(i,j,k)%theta_w = theta_avg(i,j,k)%theta_w + theta_w(i,j,k) * tavg_dt !Eshwan
      theta_avg(i,j,k)%wtheta = theta_avg(i,j,k)%wtheta + w(i,j,k) * theta_w(i,j,k) * tavg_dt !Eshwan
      theta_avg(i,j,k)%theta2 = theta_avg(i,j,k)%theta2 + theta_w(i,j,k) * theta_w(i,j,k) * tavg_dt !Eshwan
      theta_avg(i,j,k)%dTdz = theta_avg(i,j,k)%dTdz + dTdz(i,j,k) * tavg_dt !Eshwan
      theta_avg(i,j,k)%sgs_tflux = theta_avg(i,j,k)%sgs_tflux + kappa_tt(i,j,k) * dTdz(i,j,k) *tavg_dt !Eshwan

      sal_avg(i,j,k)%sal_w = sal_avg(i,j,k)%sal_w + sal_w(i,j,k) * tavg_dt !Eshwan
      sal_avg(i,j,k)%wsal = sal_avg(i,j,k)%wsal + w(i,j,k) * sal_w(i,j,k) * tavg_dt !Eshwan
      sal_avg(i,j,k)%sal2 = sal_avg(i,j,k)%sal2 + sal_w(i,j,k) * sal_w(i,j,k) * tavg_dt !Eshwan
      sal_avg(i,j,k)%dSdz = sal_avg(i,j,k)%dSdz + dSdz(i,j,k) * tavg_dt !Eshwan
      sal_avg(i,j,k)%sgs_sflux = sal_avg(i,j,k)%sgs_sflux + kappa_ts(i,j,k) * dSdz(i,j,k) * tavg_dt !Eshwan  

      ! === uv-grid variables ===
      tavg(i,j,k)%u_uv = tavg(i,j,k)%u_uv + u(i,j,k) * tavg_dt !Eshwan
      tavg(i,j,k)%v_uv = tavg(i,j,k)%v_uv + v(i,j,k) * tavg_dt !Eshwan
      tavg(i,j,k)%w_uv = tavg(i,j,k)%w_uv + w_uv(i,j,k) * tavg_dt !Eshwan
      
      tavg(i,j,k)%txx = tavg(i,j,k)%txx + txx(i,j,k) * tavg_dt
      tavg(i,j,k)%tyy = tavg(i,j,k)%tyy + tyy(i,j,k) * tavg_dt
      tavg(i,j,k)%tzz = tavg(i,j,k)%tzz + tzz(i,j,k) * tavg_dt
      tavg(i,j,k)%txy = tavg(i,j,k)%txy + txy(i,j,k) * tavg_dt

      theta_avg(i,j,k)%theta = theta_avg(i,j,k)%theta + theta(i,j,k) * tavg_dt !Eshwan
      theta_avg(i,j,k)%dTdx = theta_avg(i,j,k)%dTdx + dTdx(i,j,k) * tavg_dt !Eshwan
      theta_avg(i,j,k)%dTdy = theta_avg(i,j,k)%dTdy + dTdy(i,j,k) * tavg_dt !Eshwan
      
      sal_avg(i,j,k)%sal = sal_avg(i,j,k)%sal + sal(i,j,k) * tavg_dt !Eshwan
      sal_avg(i,j,k)%dSdx = sal_avg(i,j,k)%dSdx + dSdx(i,j,k) * tavg_dt !Eshwan
      sal_avg(i,j,k)%dSdy = sal_avg(i,j,k)%dSdy + dSdy(i,j,k) * tavg_dt !Eshwan
 
      ! === w-grid variables === 
      tavg(i,j,k)%txz = tavg(i,j,k)%txz + txz(i,j,k) * tavg_dt
      tavg(i,j,k)%tyz = tavg(i,j,k)%tyz + tyz(i,j,k) * tavg_dt

      tavg(i,j,k)%cs_opt2 = tavg(i,j,k)%cs_opt2 + Cs_opt2(i,j,k) * tavg_dt
      tavg(i,j,k)%Nu_t = tavg(i,j,k)%Nu_t + Nu_t(i,j,k) * tavg_dt !Eshwan
      
    enddo
  enddo
enddo

!Average friction velocity, friction temperature, friction, salinity, and 
!surface scalar fluxes
do j=1,ny
   do i=1,nx
      ustar_avg(i,j) = ustar_avg(i,j) + ustar(i,j)*tavg_dt
      t_flux_avg(i,j) = t_flux_avg(i,j) + t_flux(i,j)*tavg_dt
      tstar_avg(i,j) = tstar_avg(i,j) + tstar(i,j)*tavg_dt
      sal_flux_avg(i,j) = sal_flux_avg(i,j) + sal_flux(i,j)*tavg_dt
      sstar_avg(i,j) = sstar_avg(i,j) + sstar(i,j)*tavg_dt
   enddo
enddo

!$if ($OUTPUT_EXTRA)
if (sgs_model.gt.1) then
   do k=1,jzmax       
      do j=1,ny
         do i=1,nx
            ! === w-grid variables ===
            tavg_sgs(i,j,k)%ee_now = tavg_sgs(i,j,k)%ee_now + ee_now(i,j,k) * tavg_dt
         enddo
      enddo
   enddo
endif

if (sgs_model.gt.3) then
   do k=1,jzmax       
      do j=1,ny
         do i=1,nx
            ! === w-grid variables ===
            tavg_sgs(i,j,k)%Tn = tavg_sgs(i,j,k)%Tn + Tn_all(i,j,k) * tavg_dt !Eshwan
            tavg_sgs(i,j,k)%Beta = tavg_sgs(i,j,k)%Beta + Beta(i,j,k) !Eshwan
            tavg_sgs(i,j,k)%F_LM = tavg_sgs(i,j,k)%F_LM + F_LM(i,j,k) * tavg_dt
            tavg_sgs(i,j,k)%F_MM = tavg_sgs(i,j,k)%F_MM + F_MM(i,j,k) * tavg_dt
            tavg_sgs(i,j,k)%F_QN = tavg_sgs(i,j,k)%F_QN + F_QN(i,j,k) * tavg_dt
            tavg_sgs(i,j,k)%F_NN = tavg_sgs(i,j,k)%F_NN + F_NN(i,j,k) * tavg_dt
            tavg_sgs(i,j,k)%cs2_clips = tavg_sgs(i,j,k)%cs2_clips + cs2_clips(i,j,k)
            
            $if ($DYN_TN)
            tavg_sgs(i,j,k)%F_ee2 = tavg_sgs(i,j,k)%F_ee2 + F_ee2(i,j,k) * tavg_dt
            tavg_sgs(i,j,k)%F_deedt2 = tavg_sgs(i,j,k)%F_deedt2 + F_deedt2(i,j,k) * tavg_dt
            $endif         
         enddo
      enddo
   enddo
endif
!$endif


do k=1,jzmax       
  do j=1,ny
    do i=1,nx
       ! === w-grid variables ===
       tavg_scalar_sgs(i,j,k)%Ds_opt2_t = tavg_scalar_sgs(i,j,k)%Ds_opt2_t + Ds_opt2_t(i,j,k) * tavg_dt
       tavg_scalar_sgs(i,j,k)%kappa_tt = tavg_scalar_sgs(i,j,k)%kappa_tt + kappa_tt(i,j,k) * tavg_dt
       tavg_scalar_sgs(i,j,k)%s_Tn_t = tavg_scalar_sgs(i,j,k)%s_Tn_t + s_Tn_all_t(i,j,k) * tavg_dt
       tavg_scalar_sgs(i,j,k)%s_Beta_t = tavg_scalar_sgs(i,j,k)%s_Beta_t + s_Beta_t(i,j,k) * tavg_dt
       tavg_scalar_sgs(i,j,k)%I_LM_t = tavg_scalar_sgs(i,j,k)%I_LM_t + I_LM_t(i,j,k) * tavg_dt
       tavg_scalar_sgs(i,j,k)%I_MM_t = tavg_scalar_sgs(i,j,k)%I_MM_t + I_MM_t(i,j,k) * tavg_dt
       tavg_scalar_sgs(i,j,k)%I_QN_t = tavg_scalar_sgs(i,j,k)%I_QN_t + I_QN_t(i,j,k) * tavg_dt
       tavg_scalar_sgs(i,j,k)%I_NN_t = tavg_scalar_sgs(i,j,k)%I_NN_t + I_NN_t(i,j,k) * tavg_dt
       tavg_scalar_sgs(i,j,k)%ds2_clips_t = tavg_scalar_sgs(i,j,k)%ds2_clips_t + ds2_clips_t(i,j,k)

       tavg_scalar_sgs(i,j,k)%Ds_opt2_s = tavg_scalar_sgs(i,j,k)%Ds_opt2_s + Ds_opt2_s(i,j,k) * tavg_dt
       tavg_scalar_sgs(i,j,k)%kappa_ts = tavg_scalar_sgs(i,j,k)%kappa_ts + kappa_ts(i,j,k) * tavg_dt
       tavg_scalar_sgs(i,j,k)%s_Tn_s = tavg_scalar_sgs(i,j,k)%s_Tn_s + s_Tn_all_s(i,j,k) * tavg_dt
       tavg_scalar_sgs(i,j,k)%s_Beta_s = tavg_scalar_sgs(i,j,k)%s_Beta_s + s_Beta_s(i,j,k) * tavg_dt
       tavg_scalar_sgs(i,j,k)%I_LM_s = tavg_scalar_sgs(i,j,k)%I_LM_s + I_LM_s(i,j,k) * tavg_dt
       tavg_scalar_sgs(i,j,k)%I_MM_s = tavg_scalar_sgs(i,j,k)%I_MM_s + I_MM_s(i,j,k) * tavg_dt
       tavg_scalar_sgs(i,j,k)%I_QN_s = tavg_scalar_sgs(i,j,k)%I_QN_s + I_QN_s(i,j,k) * tavg_dt
       tavg_scalar_sgs(i,j,k)%I_NN_s = tavg_scalar_sgs(i,j,k)%I_NN_s + I_NN_s(i,j,k) * tavg_dt
       tavg_scalar_sgs(i,j,k)%ds2_clips_s = tavg_scalar_sgs(i,j,k)%ds2_clips_s + ds2_clips_s(i,j,k)
    enddo
 enddo
enddo


deallocate( u_w, v_w )
deallocate( theta_w )
deallocate ( sal_w ) 

! Added to compute the average of u on the grid points and prevent the 
! linear interpolation, which is not exactly valid due to log profile in boundary layer
do k=lbz,jzmax  
  do j=1,ny
    do i=1,nx
    u_p = u(i,j,k)
    u_avg(i,j,k)=u_avg(i,j,k) + u_p*tavg_dt                   
    enddo
  enddo
enddo

! Update tavg_total_time for variable time stepping
tavg_total_time = tavg_total_time + tavg_dt
!$if($OUTPUT_EXTRA)
tavg_total_time_sgs = tavg_total_time_sgs + tavg_dt
!$endif

! Set tavg_dt back to zero for next increment
tavg_dt = 0.0_rprec

return

end subroutine tavg_compute

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine tavg_finalize()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use grid_defs, only : grid !x,y,z
use stat_defs, only : tavg_t, theta_avg_t, tavg_zplane, tavg_total_time, tavg, &
                      theta_avg, theta_avg_zplane, t_flux_avg, &
                      t_flux_avg_zplane, tstar_avg, tstar_avg_zplane, ustar_avg, &
                      ustar_avg_zplane, sal_flux_avg, sal_flux_avg_zplane, &
                      sstar_avg, sstar_avg_zplane, sal_avg_t, sal_avg, & 
                      sal_avg_zplane    !Eshwan
use stat_defs, only : rs_t, rs, rs_zplane, cnpy_zplane, theta_stat_t, theta_stat, &
                      theta_stat_zplane, sal_stat_t, sal_stat, sal_stat_zplane, & 
                      sigma_t, sigma_avg_zplane !Eshwan 
use stat_defs, only : operator(.DIV.), operator(.MUL.)
use stat_defs, only :  operator(.ADD.), operator(.SUB.)
use stat_defs, only : type_set, type_zero_bogus
use stat_defs, only : tavg_interp_to_uv_grid, tavg_interp_to_w_grid
use stat_defs, only : rs_compute, cnpy_tavg_mul, theta_stat_compute, &
                      sal_stat_compute
!$if($OUTPUT_EXTRA)i
use stat_defs, only : tavg_sgs_t, tavg_scalar_sgs_t
use stat_defs, only : tavg_sgs, tavg_scalar_sgs, tavg_total_time_sgs, & !Eshwan
                      tavg_sgs_zplane, tavg_scalar_sgs_zplane
!$endif
use param, only : nx,ny,nz,dx,dy,dz,L_x,L_y,L_z, nz_tot
use param, only : output_path
use sim_param, only : u_avg
$if($MPI)
use mpi
use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWNUP
use param, only : MPI_RPREC, coord_of_rank, &
     rank_of_coord, comm, &
     ierr, down, up, status
use stat_defs, only : rs_t, tavg_t
$endif
$if($LVLSET)
use level_set_base, only : phi
$endif

implicit none

include 'tecryte.h'

character (*), parameter :: sub_name = mod_name // '.tavg_finalize'
character(64) :: fname_vel, &
     fname_vel2, fname_ddz, &
     fname_tau, fname_f, &
     fname_rs, fname_cs, fname_u_vel_grid, &
     fname_theta, fname_theta_stat, & !Eshwan
     fname_sal, fname_sal_stat, &
     fname_divt, fname_buoy !Eshwan
character(64) :: fname_velb, &
     fname_vel2b, fname_ddzb, &
     fname_taub, fname_fb, &
     fname_rsb, fname_csb, fname_u_vel_gridb, &
     fname_thetab, fname_salb, &
     fname_divtb, fname_buoyb !Eshwan
     
character(64) :: fname_vel_zplane, fname_vel2_zplane, &
  fname_ddz_zplane, fname_tau_zplane, fname_f_zplane, &
  fname_rs_zplane, fname_cs_zplane, fname_cnpy_zplane,&
  fname_theta_zplane, fname_theta_stat_zplane, &
  fname_sal_zplane, fname_sal_stat_zplane, &
  fname_divt_zplane, fname_buoy_zplane, &
  fname_sigma_zplane !Eshwan
!$if($OUTPUT_EXTRA)  
character(64) :: fname_sgs_Nu, fname_sgs_Fsub
character(64) :: fname_scalar_sgs
character(64) :: fname_sgs_Isub_t, fname_sgs_Isub_s
  !$if($DYN_TN)
  character(64) :: fname_sgs_ee
  !$endif
character(64) :: fname_sgs_Nu_zplane, fname_scalar_sgs_zplane
character(64) :: fname_sgs_Fsub_zplane
character(64) :: fname_sgs_Isub_t_zplane, fname_sgs_Isub_s_zplane
!$endif  

integer :: i,j,k

$if($MPI)
!character(64) :: temp

integer :: MPI_RS, MPI_CNPY, MPI_TAVG, MPI_THETA_AVG,  &
           MPI_THETA_STAT, MPI_TSTAT_AVG, MPI_SAL_AVG, MPI_SAL_STAT, &
           MPI_SGS, MPI_SCALAR_SGS, MPI_SIGMA !Eshwan
integer :: rs_type(1), rs_block(1), rs_disp(1)
integer :: theta_stat_type(1), theta_stat_block(1), theta_stat_disp(1) !Eshwan
integer :: sal_stat_type(1), sal_stat_block(1), sal_stat_disp(1) !Eshwan
integer :: cnpy_type(1), cnpy_block(1), cnpy_disp(1)
integer :: tavg_type(1), tavg_block(1), tavg_disp(1)
integer :: theta_avg_type(1), theta_avg_block(1), theta_avg_disp(1) !Eshwan
integer :: sal_avg_type(1), sal_avg_block(1), sal_avg_disp(1) !Eshwan
integer :: tavg_sgs_type(1), tavg_sgs_block(1), tavg_sgs_disp(1) !Eshwan
integer :: tavg_scalar_sgs_type(1), tavg_scalar_sgs_block(1), tavg_scalar_sgs_disp(1) !Eshwan
integer :: sigma_avg_type(1), sigma_avg_block(1), sigma_avg_disp(1) !Eshwan

! Definitions for reconstructing z-planar averaged data
integer :: sendsize
integer, allocatable, dimension(:) :: recvsize, recvstride
real(rprec), allocatable, dimension(:) :: z_tot, zw_tot
type(rs_t), allocatable, dimension(:) :: rs_zplane_tot
type(theta_stat_t), allocatable, dimension (:) :: theta_stat_zplane_tot    !Eshwan
type(sal_stat_t), allocatable, dimension (:) :: sal_stat_zplane_tot    !Eshwan
type(rs_t), allocatable, dimension(:) :: cnpy_zplane_tot
type(tavg_t), allocatable, dimension(:) :: tavg_zplane_tot
type(theta_avg_t), allocatable, dimension(:) :: theta_avg_zplane_tot   !Eshwan
type(sal_avg_t), allocatable, dimension(:) :: sal_avg_zplane_tot   !Eshwan
type(tavg_sgs_t), allocatable, dimension(:) :: tavg_sgs_zplane_tot
type(tavg_scalar_sgs_t), allocatable, dimension(:) :: tavg_scalar_sgs_zplane_tot
type(sigma_t), allocatable, dimension(:) :: sigma_avg_zplane_tot !Eshwan
$endif

type(rs_t) :: cnpy_avg
type(tavg_t) :: tavg_avg

real(rprec) :: favg

$if($LVLSET)
real(rprec) :: fx_global, fy_global, fz_global

$endif

real(rprec), pointer, dimension(:) :: x,y,z,zw

favg = real(nx*ny,kind=rprec)

nullify(x,y,z,zw)

allocate(rs(nx,ny,lbz:nz), rs_zplane(nz))
allocate(theta_stat(nx,ny,lbz:nz), theta_stat_zplane(nz))
allocate(sal_stat(nx,ny,lbz:nz), sal_stat_zplane(nz))
allocate(cnpy_zplane(nz))
allocate(tavg_sgs_zplane(lbz:nz))
allocate( tavg_scalar_sgs_zplane(lbz:nz))
x => grid % x
y => grid % y
z => grid % z
zw => grid % zw

! All processors need not do this, but that is ok
!  Set file names
fname_vel = output_path // 'output/vel_avg.dat'
fname_vel2 = output_path // 'output/vel2_avg.dat'
fname_ddz = output_path // 'output/ddz_avg.dat'
fname_tau = output_path // 'output/tau_avg.dat'
fname_f = output_path // 'output/force_avg.dat'
fname_rs = output_path // 'output/rs.dat'
fname_theta_stat = output_path // 'output/T_stat_avg.dat'     !Eshwan
fname_sal_stat = output_path // 'output/S_stat_avg.dat'     !Eshwan
fname_cs = output_path // 'output/cs_opt2.dat'
fname_u_vel_grid = output_path // 'output/u_grid_vel.dat'
fname_theta = output_path // 'output/T_avg.dat'       !Eshwan
fname_sal = output_path // 'output/S_avg.dat'       !Eshwan
fname_divt = output_path // 'output/divt_avg.dat'
fname_buoy = output_path // 'output/buoy_avg.dat'

fname_velb = output_path // 'output/binary_vel_avg.dat'
fname_vel2b = output_path // 'output/binary_vel2_avg.dat'
fname_ddzb = output_path // 'output/binary_ddz_avg.dat'
fname_taub = output_path // 'output/binary_tau_avg.dat'
fname_fb = output_path // 'output/binary_force_avg.dat'
fname_rsb = output_path // 'output/binary_rs.dat'
fname_csb = output_path // 'output/binary_cs_opt2.dat'
fname_u_vel_gridb = output_path // 'output/binary_u_grid_vel.dat'
fname_thetab = output_path // 'output/binary_T.dat'   !Eshwan
fname_salb = output_path // 'output/binary_S.dat'   !Eshwan
fname_divtb = output_path // 'output/binary_divt.dat'
fname_buoyb = output_path // 'output/binary_buoy.dat'

fname_vel_zplane = output_path // 'output/vel_zplane_avg.dat'
fname_vel2_zplane = output_path // 'output/vel2_zplane.dat'
fname_ddz_zplane = output_path // 'output/ddz_zplane_avg.dat'
fname_tau_zplane = output_path // 'output/tau_zplane_avg.dat'
fname_f_zplane = output_path // 'output/force_zplane_avg.dat'
fname_rs_zplane = output_path // 'output/rs_zplane.dat'
fname_theta_stat_zplane = output_path // 'output/T_stat_zplane.dat'
fname_sal_stat_zplane = output_path // 'output/S_stat_zplane.dat'
fname_cnpy_zplane = output_path // 'output/cnpy_zplane.dat'
fname_cs_zplane = output_path // 'output/cs_opt2_zplane.dat'
fname_theta_zplane = output_path // 'output/T_zplane_avg.dat'     !Eshwan
fname_sal_zplane = output_path // 'output/S_zplane_avg.dat'       !Eshwan
fname_divt_zplane = output_path // 'output/divt_zplane.dat'
fname_buoy_zplane = output_path // 'output/buoy_zplane.dat'
fname_sigma_zplane = output_path // 'output/sigma_avg_zplane.dat' !Eshwan 

!$if($OUTPUT_EXTRA)  
fname_sgs_Nu = output_path // 'output/Nu_avg.dat'
fname_sgs_Fsub = output_path // 'output/F_avg.dat'
fname_scalar_sgs = output_path // 'output/ssgs_avg.dat'
fname_sgs_Isub_t = output_path // 'output/I_t_avg.dat'
fname_sgs_Isub_s = output_path // 'output/I_s_avg.dat'
!$if($DYN_TN)
fname_sgs_ee = output_path // 'output/ee_avg.dat'
!$endif
fname_sgs_Nu_zplane = output_path // 'output/Nu_zplane_avg.dat'
fname_scalar_sgs_zplane = output_path // 'output/ssgs_zplane.dat'
fname_sgs_Fsub_zplane = output_path // 'output/F_zplane_avg.dat'
fname_sgs_Isub_t_zplane = output_path // 'output/I_t_zplane_avg.dat'
fname_sgs_Isub_s_zplane = output_path // 'output/I_s_zplane_avg.dat'
!$endif  
 
$if ($MPI)
  $if($BINARY)
  !  For MPI implementation     
  call string_concat( fname_velb, '.c', coord)
  call string_concat( fname_vel2b, '.c', coord)
  call string_concat( fname_ddzb, '.c', coord)
  call string_concat( fname_taub, '.c', coord)
  call string_concat( fname_fb, '.c', coord)
  call string_concat( fname_rsb, '.c', coord)
  call string_concat( fname_csb, '.c', coord)
  call string_concat( fname_u_vel_gridb, '.c',coord)
  call string_concat( fname_thetab, '.c',coord)                  !Eshwan
  call string_concat( fname_salb, '.c',coord)                  !Eshwan
  $else
  call string_concat( fname_vel, '.c', coord)
  call string_concat( fname_vel2, '.c', coord)
  call string_concat( fname_ddz, '.c', coord)
  call string_concat( fname_tau, '.c', coord)
  call string_concat( fname_f, '.c', coord)
  call string_concat( fname_rs, '.c', coord)
  call string_concat( fname_theta_stat, '.c', coord)
  call string_concat( fname_sal_stat, '.c', coord)
  call string_concat( fname_cs, '.c', coord)
  call string_concat( fname_u_vel_grid, '.c', coord)
  call string_concat( fname_theta, '.c', coord)                  !Eshwan
  call string_concat( fname_sal, '.c', coord)                  !Eshwan
  $endif
  
  !$if($OUTPUT_EXTRA)  
  call string_concat( fname_sgs_Nu, '.c', coord)
  call string_concat( fname_sgs_Fsub, '.c', coord)
  call string_concat( fname_scalar_sgs, '.c', coord)
  call string_concat( fname_sgs_Isub_t, '.c', coord)
  call string_concat( fname_sgs_Isub_s, '.c', coord)
 !$if($DYN_TN)
  call string_concat( fname_sgs_ee, '.c', coord)
 !$endif
  !$endif    
$endif

! Final checkpoint all restart data
call tavg_checkpoint()

$if($MPI)
call mpi_barrier( comm, ierr )
$endif

! Zero bogus values
call type_zero_bogus( tavg(:,:,nz) )

$if($MPI)
! Build MPI types for derived types
rs_type = MPI_RPREC
rs_block = 6 ! Number of rs subtypes
rs_disp = 0

! Eshwan: theta flux and variance
theta_stat_type = MPI_RPREC
theta_stat_block = 3
theta_stat_disp = 0

! Eshwan: sal flux and variance
sal_stat_type = MPI_RPREC
sal_stat_block = 3
sal_stat_disp = 0

cnpy_type = MPI_RPREC
cnpy_block = 6 ! Number of cnpy subtypes
cnpy_disp = 0

tavg_type = MPI_RPREC
tavg_block = 30 ! Number of tavg subtypes
tavg_disp = 0

! Eshwan: average theta
theta_avg_type = MPI_RPREC ! Number of theta_avg subtypes
theta_avg_block = 8
theta_avg_disp = 0

! Eshwan: average sal
sal_avg_type = MPI_RPREC ! Number of sal_avg subtypes
sal_avg_block = 8
sal_avg_disp = 0

!Eshwan: sgs 
tavg_sgs_type = MPI_RPREC
tavg_sgs_block = 8
tavg_sgs_disp = 0 

!Eshwan: scalar sgs
tavg_scalar_sgs_type = MPI_RPREC
tavg_scalar_sgs_block = 18
tavg_scalar_sgs_disp = 0

!Eshwan: standard deviation of temperature and salinity
sigma_avg_type = MPI_RPREC
sigma_avg_block = 2
sigma_avg_disp = 0

!  Create MPI type structures consistent with the derived types
call MPI_TYPE_STRUCT(1, rs_block, rs_disp, rs_type, MPI_RS, ierr)
if(ierr /= 0) call error(sub_name,'Error in setting MPI_RS:', ierr)
Call MPI_Type_commit(MPI_RS,ierr)
if(ierr /= 0) call error(sub_name,'Error in committing MPI_RS:', ierr)

call MPI_TYPE_STRUCT(1, theta_stat_block, theta_stat_disp, theta_stat_type, MPI_THETA_STAT, ierr)
if(ierr /= 0) call error(sub_name,'Error in setting MPI_THETA_STAT:', ierr)
Call MPI_Type_commit(MPI_THETA_STAT,ierr)
if(ierr /= 0) call error(sub_name,'Error in committing MPI_THETA_STAT:', ierr)

call MPI_TYPE_STRUCT(1, sal_stat_block, sal_stat_disp, sal_stat_type, MPI_SAL_STAT, ierr)
if(ierr /= 0) call error(sub_name,'Error in setting MPI_SAL_STAT:', ierr)
Call MPI_Type_commit(MPI_SAL_STAT,ierr)
if(ierr /= 0) call error(sub_name,'Error in committing MPI_SAL_STAT:', ierr)

call MPI_TYPE_STRUCT(1, cnpy_block, cnpy_disp, cnpy_type, MPI_CNPY, ierr)
if(ierr /= 0) call error(sub_name,'Error in setting MPI_CNPY:', ierr)
Call MPI_Type_commit(MPI_CNPY,ierr)  
if(ierr /= 0) call error(sub_name,'Error in committing MPI_CNPY:', ierr)
   
call MPI_TYPE_STRUCT(1, tavg_block, tavg_disp, tavg_type, MPI_TAVG, ierr)
if(ierr /= 0) call error(sub_name,'Error in setting MPI_TAVG:', ierr)
Call MPI_Type_commit(MPI_TAVG,ierr)
if(ierr /= 0) call error(sub_name,'Error in committing MPI_TAVG:', ierr)

call MPI_TYPE_STRUCT(1, theta_avg_block, theta_avg_disp, theta_avg_type, MPI_THETA_AVG, ierr)
if(ierr /= 0) call error(sub_name,'Error in setting MPI_THETA_AVG:', ierr)
Call MPI_Type_commit(MPI_THETA_AVG,ierr)
if(ierr /= 0) call error(sub_name,'Error in committing MPI_THETA_AVG:', ierr)

call MPI_TYPE_STRUCT(1, sal_avg_block, sal_avg_disp, sal_avg_type, MPI_SAL_AVG, ierr)
if(ierr /= 0) call error(sub_name,'Error in setting MPI_SAL_AVG:', ierr)
Call MPI_Type_commit(MPI_SAL_AVG,ierr)
if(ierr /= 0) call error(sub_name,'Error in committing MPI_SAL_AVG:', ierr)

call MPI_TYPE_STRUCT(1, tavg_sgs_block, tavg_sgs_disp, tavg_sgs_type, MPI_SGS, ierr)
if(ierr /= 0) call error(sub_name,'Error in setting MPI_SGS:', ierr)
Call MPI_Type_commit(MPI_SGS,ierr)
if(ierr /= 0) call error(sub_name,'Error in committing MPI_SGS:', ierr)

call MPI_TYPE_STRUCT(1, tavg_scalar_sgs_block, tavg_scalar_sgs_disp, tavg_scalar_sgs_type, MPI_SCALAR_SGS, ierr)
if(ierr /= 0) call error(sub_name,'Error in setting MPI_SCALAR_SGS:', ierr)
Call MPI_Type_commit(MPI_SCALAR_SGS,ierr)
if(ierr /= 0) call error(sub_name,'Error in committing MPI_SCALAR_SGS:', ierr)

call MPI_TYPE_STRUCT(1, sigma_avg_block, sigma_avg_disp, sigma_avg_type, MPI_SIGMA, ierr)
if(ierr /= 0) call error(sub_name,'Error in setting MPI_SIGMA:', ierr)
Call MPI_Type_commit(MPI_SIGMA,ierr)
if(ierr /= 0) call error(sub_name,'Error in committing MPI_SIGMA:', ierr)


!  Allocate space only on base processor for assembled z-plane data
! *_tot_t is the assembled data without the overlap nodes (the final stuff that is outputted)
! All processes need to allocate the space even though it is not directly used
allocate(z_tot(nz_tot))
allocate(zw_tot(nz_tot))
allocate(rs_zplane_tot(nz_tot))
allocate(theta_stat_zplane_tot(nz_tot)) !Eshwan
allocate(sal_stat_zplane_tot(nz_tot)) !Eshwan
allocate(cnpy_zplane_tot(nz_tot))
allocate(tavg_zplane_tot(nz_tot))
allocate(theta_avg_zplane_tot(nz_tot)) !Eshwan
allocate(sal_avg_zplane_tot(nz_tot)) !Eshwan
allocate(tavg_sgs_zplane_tot(nz_tot)) !Eshwan
allocate(tavg_scalar_sgs_zplane_tot(nz_tot)) !Eshwan
allocate(sigma_avg_zplane_tot(nz_tot)) !Eshwan
$endif

!  Perform time averaging operation
!  tavg = tavg / tavg_total_time

do k=jzmin,jzmax
  sigma_avg_zplane(k) = sigma_avg_zplane(k) .DIV. tavg_total_time !Eshwan
  do j=1, Ny
    do i=1, Nx
      tavg(i,j,k) = tavg(i,j,k) .DIV. tavg_total_time
      u_avg(i,j,k) =u_avg(i,j,k) / tavg_total_time

      theta_avg(i,j,k) = theta_avg(i,j,k) .DIV. tavg_total_time !Eshwan 
      sal_avg(i,j,k) = sal_avg(i,j,k) .DIV. tavg_total_time !Eshwan 
    enddo
  enddo
enddo

!Average friction quantities and surface fluxes
do j=1,Ny
   do i=1,Nx
      ustar_avg(i,j) = ustar_avg(i,j) / tavg_total_time
      t_flux_avg(i,j) = t_flux_avg(i,j) / tavg_total_time
      tstar_avg(i,j) = tstar_avg(i,j) / tavg_total_time
      sal_flux_avg(i,j) = sal_flux_avg(i,j) / tavg_total_time
      sstar_avg(i,j) = sstar_avg(i,j) / tavg_total_time
   end do
end do

!$if ($OUTPUT_EXTRA)
do k=1,jzmax
  do j=1, Ny
    do i=1, Nx
      tavg_sgs(i,j,k) = tavg_sgs(i,j,k) .DIV. tavg_total_time_sgs
    enddo
  enddo
enddo

do k=1,jzmax
  do j=1, Ny
    do i=1, Nx
      tavg_scalar_sgs(i,j,k) = tavg_scalar_sgs(i,j,k) .DIV. tavg_total_time
    enddo
  enddo
enddo
!$endif

!  Sync entire tavg structure
$if($MPI)
call mpi_sendrecv (tavg(:,:,1), nx*ny, MPI_TAVG, down, 1,  &
                   tavg(:,:,nz), nx*ny, MPI_TAVG, up, 1,   &
                   comm, status, ierr)
call mpi_sendrecv (tavg(:,:,nz-1), nx*ny, MPI_TAVG, up, 2,  &
                   tavg(:,:,0), nx*ny, MPI_TAVG, down, 2,   &
                   comm, status, ierr)

call mpi_sendrecv (theta_avg(:,:,1), nx*ny, MPI_THETA_AVG, down, 3,  &
                   theta_avg(:,:,nz), nx*ny, MPI_THETA_AVG, up, 3,   &
                   comm, status, ierr)                                   !Eshwan
call mpi_sendrecv (theta_avg(:,:,nz-1), nx*ny, MPI_THETA_AVG, up, 4,  &
                   theta_avg(:,:,0), nx*ny, MPI_THETA_AVG, down, 4,   &
                   comm, status, ierr)                                   !Eshwan

call mpi_sendrecv (sal_avg(:,:,1), nx*ny, MPI_SAL_AVG, down, 5,  &
                   sal_avg(:,:,nz), nx*ny, MPI_SAL_AVG, up, 5,   &
                   comm, status, ierr)                                   !Eshwan
call mpi_sendrecv (sal_avg(:,:,nz-1), nx*ny, MPI_SAL_AVG, up, 6,  &
                   sal_avg(:,:,0), nx*ny, MPI_SAL_AVG, down, 6,   &
                   comm, status, ierr)                                   !Eshwan

call mpi_sendrecv (tavg_scalar_sgs(:,:,1), nx*ny, MPI_SCALAR_SGS, down, 7,  &
                   tavg_scalar_sgs(:,:,nz), nx*ny, MPI_SCALAR_SGS, up, 7,   &
                   comm, status, ierr)                                   !Eshwan
call mpi_sendrecv (tavg_scalar_sgs(:,:,nz-1), nx*ny, MPI_SCALAR_SGS, up, 8,  &
                   tavg_scalar_sgs(:,:,0), nx*ny, MPI_SCALAR_SGS, down, 8,   &
                   comm, status, ierr)                                   !Eshwan 

call mpi_sendrecv (tavg_sgs(:,:,1), nx*ny, MPI_SGS, down, 9,  &
                   tavg_sgs(:,:,nz), nx*ny, MPI_SGS, up, 9,   &
                   comm, status, ierr)                                   !Eshwan
call mpi_sendrecv (tavg_sgs(:,:,nz-1), nx*ny, MPI_SGS, up, 10,  &
                   tavg_sgs(:,:,0), nx*ny, MPI_SGS, down, 10,   &
                   comm, status, ierr)                                   !Eshwan 

call mpi_sendrecv (sigma_avg_zplane(1), 1, MPI_SIGMA, down, 11,  &
                   sigma_avg_zplane(nz), 1, MPI_SIGMA, up, 11,   &
                   comm, status, ierr)                                   !Eshwan
call mpi_sendrecv (sigma_avg_zplane(nz-1), 1, MPI_SIGMA, up, 12,  &
                   sigma_avg_zplane(0), 1, MPI_SIGMA, down, 12,   &
                   comm, status, ierr)                                   !Eshwan 


$endif

! Anything with velocity is on the w-grid so these
!   values for coord==0 at k=1 should be zeroed (otherwise bogus)
if (jzmin==0) then  !coord==0 only
  tavg(:,:,0:1)%u = 0.0_rprec
  tavg(:,:,0:1)%v = 0.0_rprec
  tavg(:,:,0:1)%w = 0.0_rprec
  tavg(:,:,0:1)%u2 = 0.0_rprec
  tavg(:,:,0:1)%v2 = 0.0_rprec
  tavg(:,:,0:1)%w2 = 0.0_rprec
  tavg(:,:,0:1)%uw = 0.0_rprec
  tavg(:,:,0:1)%vw = 0.0_rprec
  tavg(:,:,0:1)%uv = 0.0_rprec
  tavg(:,:,0:1)%divtx = 0.0_rprec
  tavg(:,:,0:1)%divty = 0.0_rprec
endif

! Interpolate between grids where necessary
!tavg = tavg_interp_to_uv_grid( tavg )
tavg = tavg_interp_to_w_grid( tavg )

! Anything with velocity is on the w-grid so these
!   values for coord==0 at k=1 should be zeroed (otherwise bogus)
if (jzmin==0) then  !coord==0 only
  tavg(:,:,0:1)%fx = 0.0_rprec 
  tavg(:,:,0:1)%fy = 0.0_rprec   

  tavg(:,:,0:1)%txx = 0.0_rprec    
  tavg(:,:,0:1)%txy = 0.0_rprec    
  tavg(:,:,0:1)%tyy = 0.0_rprec    
  tavg(:,:,0:1)%tzz = 0.0_rprec
endif

!Values for coord==0 at k=1 should be zeroed
if (jzmin==0) then !coord==0 only
    theta_avg(:,:,0:1)%theta_w = 0.0_rprec
    theta_avg(:,:,0:1)%wtheta = 0.0_rprec
    theta_avg(:,:,0:1)%theta2 = 0.0_rprec
end if

!Values for coord==0 at k=1 should be zeroed
if (jzmin==0) then !coord==0 only
    sal_avg(:,:,0:1)%sal_w = 0.0_rprec
    sal_avg(:,:,0:1)%wsal = 0.0_rprec
    sal_avg(:,:,0:1)%sal2 = 0.0_rprec
end if

!Values for coord==0 at k=1 should be zeroed
if (jzmin==0) then !coord==0 only
   sigma_avg_zplane(0:1)%sigma_theta = 0.0_rprec
   sigma_avg_zplane(0:1)%sigma_sal = 0.0_rprec
end if

!  Average over z-planes
do k=1,nz
  
  !  Initialize to 0 for summations
  call type_set( tavg_zplane(k), 0._rprec )
  call type_set( theta_avg_zplane(k), 0._rprec )
  call type_set( sal_avg_zplane(k), 0._rprec )
  call type_set( tavg_sgs_zplane(k), 0._rprec )
  call type_set( tavg_scalar_sgs_zplane(k), 0._rprec ) 
  do j=1, Ny
    do i=1, Nx

      tavg_zplane(k) = tavg_zplane(k) .ADD. tavg(i,j,k)
      theta_avg_zplane(k) = theta_avg_zplane(k) .ADD. theta_avg(i,j,k) !Eshwan
      sal_avg_zplane(k) = sal_avg_zplane(k) .ADD. sal_avg(i,j,k) !Eshwan 
      tavg_sgs_zplane(k) = tavg_sgs_zplane(k) .ADD. tavg_sgs(i,j,k) !Eshwan
      tavg_scalar_sgs_zplane(k) = tavg_scalar_sgs_zplane(k) .ADD. tavg_scalar_sgs(i,j,k) !Eshwan
    enddo
  enddo

  !  Divide by number of summation points 
  tavg_zplane(k) = tavg_zplane(k) .DIV. favg
  theta_avg_zplane(k) = theta_avg_zplane(k) .DIV. favg !Eshwan
  sal_avg_zplane(k) = sal_avg_zplane(k) .DIV. favg !Eshwan
  tavg_sgs_zplane(k) = tavg_sgs_zplane(k) .DIV. favg !Eshwan
  tavg_scalar_sgs_zplane(k) = tavg_scalar_sgs_zplane(k) .DIV. favg !Eshwan
end do  

!Average surface fluxes and friction quantities over horizontal plane
ustar_avg_zplane = 0._rprec
t_flux_avg_zplane = 0._rprec
tstar_avg_zplane = 0._rprec
sal_flux_avg_zplane = 0._rprec
sstar_avg_zplane = 0._rprec
do j=1,Ny
   do i=1,Nx
      ustar_avg_zplane = ustar_avg_zplane + ustar_avg(i,j)
      t_flux_avg_zplane = t_flux_avg_zplane + t_flux_avg(i,j)
      tstar_avg_zplane = tstar_avg_zplane + tstar_avg(i,j)
      sal_flux_avg_zplane = sal_flux_avg_zplane + sal_flux_avg(i,j)
      sstar_avg_zplane = sstar_avg_zplane + sstar_avg(i,j)
   end do
end do
ustar_avg_zplane = ustar_avg_zplane / favg
t_flux_avg_zplane = t_flux_avg_zplane / favg
tstar_avg_zplane = tstar_avg_zplane / favg
sal_flux_avg_zplane = sal_flux_avg_zplane / favg
sstar_avg_zplane = sstar_avg_zplane / favg

if (coord==0) write(*,*) 'ustar_avg_zplane', ustar_avg_zplane
if (coord==0) write(*,*) 't_flux_avg_zplane', t_flux_avg_zplane
if (coord==0) write(*,*) 'tstar_avg_zplane', tstar_avg_zplane
if (coord==0) write(*,*) 'sal_flux_avg_zplane', sal_flux_avg_zplane

! Compute the Reynolds stresses: bar(u_i * u_j) - bar(u_i) * bar(u_j)
rs = rs_compute( tavg , lbz)

! Compute the temperature flux and variance:
! bar(w *theta) - bar(w)*bar(theta), bar(theta*theta) - bar(theta)*bar(theta) -- Eshwan
theta_stat = theta_stat_compute ( tavg, theta_avg, lbz)

! Compute the salinity flux and variance: 
!bar(w *sal) - bar(w)*bar(sal), bar(sal*sal) - bar(sal)*bar(sal) -- Eshwan
sal_stat = sal_stat_compute ( tavg, sal_avg, lbz)

! Compute planar averaged Reynolds stress
do k = 1, nz

  !  Initialize to 0
  call type_set( rs_zplane(k), 0._rprec )
  call type_set( theta_stat_zplane(k), 0._rprec ) !Eshwan
  call type_set( sal_stat_zplane(k), 0._rprec ) !Eshwan
  do j = 1, ny
    do i = 1, nx
      rs_zplane(k) = rs_zplane(k) .ADD. rs(i,j,k) 
      theta_stat_zplane(k) = theta_stat_zplane(k) .ADD. theta_stat(i,j,k) !Eshwan
      sal_stat_zplane(k) = sal_stat_zplane(k) .ADD. sal_stat(i,j,k) !Eshwan
    enddo    
  enddo

  rs_zplane(k) = rs_zplane(k) .DIV. favg
  theta_stat_zplane(k) = theta_stat_zplane(k) .DIV. favg
  sal_stat_zplane(k) = sal_stat_zplane(k) .DIV. favg
enddo


$if($MPI)
call mpi_barrier( comm, ierr )
$endif

!Write output

call write_tecplot_header_ND(fname_vel, 'rewind', 9, (/ Nx+1, Ny+1, Nz/), &
     '"x", "y", "z", "<u>","<v>","<w>", "<u_uv>", "<v_uv>", "<w_uv>"', numtostr(coord, 6), 2)
call write_real_data_3D(fname_vel, 'append', 'formatted', 6, nx, ny, nz, &
     (/ tavg(:,:,1:nz) % u, &
     tavg(:,:,1:nz) % v, &
     tavg(:,:,1:nz) % w, &
     tavg(:,:,1:nz) % u_uv, &
     tavg(:,:,1:nz) % v_uv, &
     tavg(:,:,1:nz) % w_uv /), &
     4, x, y, zw(1:nz))

call write_tecplot_header_ND(fname_vel2, 'rewind', 9, (/ Nx+1, Ny+1, Nz/), &
     '"x", "y", "z", "<u<sup>2</sup>>","<v<sup>2</sup>>","<w<sup>2</sup>>", "<uw>", "<vw>", "<uv>"', &
     numtostr(coord,6), 2)
call write_real_data_3D(fname_vel2, 'append', 'formatted', 6, nx, ny, nz, &
     (/ tavg(:,:,1:nz) % u2, &
     tavg(:,:,1:nz) % v2, &
     tavg(:,:,1:nz) % w2, &
     tavg(:,:,1:nz) % uw, &
     tavg(:,:,1:nz) % vw, &
     tavg(:,:,1:nz) % uv /), &
     4, x, y, zw(1:nz))

call write_tecplot_header_ND(fname_ddz, 'rewind', 6, (/ Nx+1, Ny+1, Nz/), &
     '"x", "y", "z", "<dudz>","<dvdz>", "<dudz2>" ', numtostr(coord,6), 2)
call write_real_data_3D(fname_ddz, 'append', 'formatted', 3, nx, ny, nz, &
     (/ tavg(:,:,1:nz) % dudz, &
     tavg(:,:,1:nz) % dvdz, &
     tavg(:,:,1:nz) % dudz2 /), &
     4, x, y, zw(1:nz))

call write_tecplot_header_ND(fname_tau, 'rewind', 9, (/ Nx+1, Ny+1, Nz/), &
     '"x", "y", "z", "<t<sub>xx</sub>>","<t<sub>xy</sub>>","<t<sub>yy</sub>>", "<t<sub>xz</sub>>", "<t<sub>yz</sub>>", "<t<sub>zz</sub>>"', &
     numtostr(coord,6), 2)
call write_real_data_3D(fname_tau, 'append', 'formatted', 6, nx, ny, nz, &
     (/ tavg(:,:,1:nz) % txx, &
     tavg(:,:,1:nz) % txy, &
     tavg(:,:,1:nz) % tyy, &
     tavg(:,:,1:nz) % txz, &
     tavg(:,:,1:nz) % tyz, &
     tavg(:,:,1:nz) % tzz /), &
     4, x, y, zw(1:nz))

call write_tecplot_header_ND(fname_divt, 'rewind', 6, (/ Nx+1, Ny+1, Nz/), &
     '"x", "y", "z", "<divtx>","<divty>","<divtz>"', numtostr(coord, 6), 2)
call write_real_data_3D(fname_divt, 'append', 'formatted', 3, nx, ny, nz, &
     (/ tavg(:,:,1:nz) % divtx, &
     tavg(:,:,1:nz) % divty, &
     tavg(:,:,1:nz) % divtz /), &
     4, x, y, zw(1:nz))

call write_tecplot_header_ND(fname_buoy, 'rewind', 4, (/ Nx+1, Ny+1, Nz/), &
     '"x", "y", "z", "<buoyancy>"', numtostr(coord, 6), 2)
call write_real_data_3D(fname_buoy, 'append', 'formatted', 1, nx, ny, nz, &
     (/ tavg(:,:,1:nz) % buoyancy /), &
     4, x, y, zw(1:nz))

call write_tecplot_header_ND(fname_f, 'rewind', 6, (/ Nx+1, Ny+1, Nz/), &
     '"x", "y", "z", "<f<sub>x</sub>>","<f<sub>y</sub>>","<f<sub>z</sub>>"', &
     numtostr(coord,6), 2)
call write_real_data_3D(fname_f, 'append', 'formatted', 3, nx, ny, nz, &
     (/ tavg(:,:,1:nz) % fx, &
     tavg(:,:,1:nz) % fy, &
     tavg(:,:,1:nz) % fz /), &
     4, x, y, zw(1:nz))

call write_tecplot_header_ND(fname_rs, 'rewind', 9, (/ Nx+1, Ny+1, Nz/), &
     '"x", "y", "z", "<upup>","<vpvp>","<wpwp>", "<upwp>", "<vpwp>", "<upvp>"', &
     numtostr(coord,6), 2)
call write_real_data_3D(fname_rs, 'append', 'formatted', 6, nx, ny, nz, &
     (/ rs(:,:,1:nz)%up2, &
     rs(:,:,1:nz)%vp2, &
     rs(:,:,1:nz)%wp2, &
     rs(:,:,1:nz)%upwp, &
     rs(:,:,1:nz)%vpwp, &
     rs(:,:,1:nz)%upvp /), &
     4, x, y, zw(1:nz))

!Eshwan:: temperature statistics
call write_tecplot_header_ND(fname_theta_stat, 'rewind', 6, (/ Nx+1, Ny+1, Nz/), &
     '"x", "y", "z",  "<wpthetap>","<thetap2>","<total_tflux>"', &
    numtostr(coord,6), 2) !Eshwan 
call write_real_data_3D(fname_theta_stat, 'append', 'formatted', 3, nx, ny, nz, &
     (/ theta_stat(:,:,1:nz)%wpthetap, &
     theta_stat(:,:,1:nz)%thetap2, &
     theta_stat(:,:,1:nz)%total_tflux /), &
     4, x, y, zw(1:nz)) !Eshwan

!Eshwan:: salinity statistics
call write_tecplot_header_ND(fname_sal_stat, 'rewind', 6, (/ Nx+1, Ny+1, Nz/), &
     '"x", "y", "z",  "<wpsalp>","<salp2>","<total_sflux>"', &
    numtostr(coord,6), 2) !Eshwan 
call write_real_data_3D(fname_sal_stat, 'append', 'formatted', 3, nx, ny, nz, &
     (/ sal_stat(:,:,1:nz)%wpsalp, &
     sal_stat(:,:,1:nz)%salp2, &
     sal_stat(:,:,1:nz)%total_sflux /), &
     4, x, y, zw(1:nz)) !Eshwan

call write_tecplot_header_ND(fname_cs, 'rewind', 4, (/ Nx+1, Ny+1, Nz/), &
     '"x", "y", "z", "<cs2>"', numtostr(coord,6), 2)
call write_real_data_3D(fname_cs, 'append', 'formatted', 1, nx, ny, nz, &
     (/ tavg(:,:,1:nz)% cs_opt2 /), &
     4, x, y, zw(1:nz))

call write_tecplot_header_ND(fname_theta,'rewind', 11, (/Nx+1, Ny+1, Nz/), &        
     '"x", "y", "z", "<theta>", "<theta_w>", "<wtheta>", "<theta2>", "<dTdx>", "<dTdy>", "<dTdz>", "<sgs_tflux>" ', &
     numtostr(coord,6), 2) !Eshwan
call write_real_data_3D(fname_theta, 'append', 'formatted', 8, nx, ny, nz, &       
     (/ theta_avg(:,:,1:nz)%theta, & 
     theta_avg(:,:,1:nz)%theta_w, &
     theta_avg(:,:,1:nz)%wtheta, &
     theta_avg(:,:,1:nz)%theta2, &
     theta_avg(:,:,1:nz)%dTdx, &
     theta_avg(:,:,1:nz)%dTdy, &
     theta_avg(:,:,1:nz)%dTdz, &
     theta_avg(:,:,1:nz)%sgs_tflux /), &
     4 , x, y, zw(1:nz)) !Eshwan/

call write_tecplot_header_ND(fname_sal,'rewind', 11, (/Nx+1, Ny+1, Nz/), &        
     '"x", "y", "z", "<sal>", "<sal_w>", "<wsal>", "<sal2>", "<dSdx>", "<dSdy>", "<dSdz>", "<sgs_sflux>" ', &
     numtostr(coord,6), 2) !Eshwan
call write_real_data_3D(fname_sal, 'append', 'formatted', 8, nx, ny, nz, &       
     (/ sal_avg(:,:,1:nz)%sal, & 
     sal_avg(:,:,1:nz)%sal_w, &
     sal_avg(:,:,1:nz)%wsal, &
     sal_avg(:,:,1:nz)%sal2, &
     sal_avg(:,:,1:nz)%dSdx, &
     sal_avg(:,:,1:nz)%dSdy, &
     sal_avg(:,:,1:nz)%dSdz, &
     sal_avg(:,:,1:nz)%sgs_sflux /), &
     4 , x, y, zw(1:nz)) !Eshwan

!----
!$if($OUTPUT_EXTRA)

  call write_tecplot_header_ND(fname_sgs_Nu, 'rewind', 4, (/ Nx+1, Ny+1, Nz/), &
       '"x", "y", "z", "<Nu_t>"', numtostr(coord,6), 2)
  call write_real_data_3D(fname_sgs_Nu, 'append', 'formatted', 1, nx, ny, nz, &
       (/ tavg(:,:,1:nz) % Nu_t /), &
       4, x, y, zw(1:nz))

  call write_tecplot_header_ND(fname_sgs_Fsub, 'rewind', 10, (/ Nx+1, Ny+1, Nz/), &
       '"x", "y", "z", "<Tn>", "<Beta>", "<F_LM>", "<F_MM>", "<F_QN>", "<F_NN>", "<cs2_clips>"', numtostr(coord,6), 2)
  call write_real_data_3D(fname_sgs_Fsub, 'append', 'formatted', 7, nx, ny, nz, &
       (/ tavg_sgs(:,:,1:nz) % Tn, &
          tavg_sgs(:,:,1:nz) % Beta, &
          tavg_sgs(:,:,1:nz) % F_LM, &
          tavg_sgs(:,:,1:nz) % F_MM, &
          tavg_sgs(:,:,1:nz) % F_QN, &
          tavg_sgs(:,:,1:nz) % F_NN, &
          tavg_sgs(:,:,1:nz) % cs2_clips /), &
       4, x, y, zw(1:nz))

  call write_tecplot_header_ND(fname_sgs_ee, 'rewind', 4, (/ Nx+1, Ny+1, Nz/), &
         '"x", "y", "z", "<ee_now>"', numtostr(coord,6), 2)
  call write_real_data_3D(fname_sgs_ee, 'append', 'formatted', 1, nx, ny, nz, &
         (/ tavg_sgs(:,:,1:nz) % ee_now /), &
         4, x, y, zw(1:nz))
    
  call write_tecplot_header_ND(fname_scalar_sgs, 'rewind', 7, (/ Nx+1, Ny+1, Nz/), &
       '"x", "y", "z", "<ds2_t>", "<ds2_s>",  "<kappa_tt>", "<kappa_ts>"', numtostr(coord,6), 2)
  call write_real_data_3D(fname_scalar_sgs, 'append', 'formatted', 4, nx, ny, nz, &
       (/ tavg_scalar_sgs(:,:,1:nz) % Ds_opt2_t, &
          tavg_scalar_sgs(:,:,1:nz) % Ds_opt2_s, &
          tavg_scalar_sgs(:,:,1:nz) % kappa_tt, &
          tavg_scalar_sgs(:,:,1:nz) % kappa_ts /), &
       4, x, y, zw(1:nz))

  call write_tecplot_header_ND(fname_sgs_Isub_t, 'rewind', 10, (/ Nx+1, Ny+1, Nz/), &
       '"x", "y", "z", "<s_Tn_t>", "<s_Beta_t>", "<I_LM_t>", "<I_MM_t>", "<I_QN_t>", "<I_NN_t>", "<ds2_clips_t>"', numtostr(coord,6), 2)
  call write_real_data_3D(fname_sgs_Isub_t, 'append', 'formatted', 7, nx, ny, nz, &
       (/ tavg_scalar_sgs(:,:,1:nz) % s_Tn_t, &
          tavg_scalar_sgs(:,:,1:nz) % s_Beta_t, &
          tavg_scalar_sgs(:,:,1:nz) % I_LM_t, &
          tavg_scalar_sgs(:,:,1:nz) % I_MM_t, &
          tavg_scalar_sgs(:,:,1:nz) % I_QN_t, &
          tavg_scalar_sgs(:,:,1:nz) % I_NN_t, &
          tavg_scalar_sgs(:,:,1:nz) % ds2_clips_t /), &
       4, x, y, zw(1:nz))

  call write_tecplot_header_ND(fname_sgs_Isub_s, 'rewind', 10, (/ Nx+1, Ny+1, Nz/), &
       '"x", "y", "z", "<s_Tn_s>", "<s_Beta_s>", "<I_LM_s>", "<I_MM_s>", "<I_QN_s>", "<I_NN_s>", "<ds2_clips_s>"', numtostr(coord,6), 2)
  call write_real_data_3D(fname_sgs_Isub_s, 'append', 'formatted', 7, nx, ny, nz, &
       (/ tavg_scalar_sgs(:,:,1:nz) % s_Tn_s, &
          tavg_scalar_sgs(:,:,1:nz) % s_Beta_s, &
          tavg_scalar_sgs(:,:,1:nz) % I_LM_s, &
          tavg_scalar_sgs(:,:,1:nz) % I_MM_s, &
          tavg_scalar_sgs(:,:,1:nz) % I_QN_s, &
          tavg_scalar_sgs(:,:,1:nz) % I_NN_s, &
          tavg_scalar_sgs(:,:,1:nz) % ds2_clips_s /), &
       4, x, y, zw(1:nz))
!$endif
!----

! Construct zplane data 
$if($MPI)

  allocate(recvsize(nproc),recvstride(nproc))

  sendsize=nz-1
  recvsize=nz-1
  recvstride=coord_of_rank*(nz-1)

  ! Set very bottom values
  if( coord == 0 ) then

     z_tot(1) = z(1)
     zw_tot(1) = zw(1)
     rs_zplane_tot(1) = rs_zplane(1)
     theta_stat_zplane_tot(1) = theta_stat_zplane(1) !Eshwan
     sal_stat_zplane_tot(1) = sal_stat_zplane(1) !Eshwan
     cnpy_zplane_tot(1) = cnpy_zplane(1)
     tavg_zplane_tot(1) = tavg_zplane(1)
     theta_avg_zplane_tot(1) = theta_avg_zplane(1) !Eshwan
     sal_avg_zplane_tot(1) = sal_avg_zplane(1) !Eshwan
     tavg_sgs_zplane_tot(1) = tavg_sgs_zplane(1) !Eshwan
     tavg_scalar_sgs_zplane_tot(1) = tavg_scalar_sgs_zplane(1) !Eshwan
     sigma_avg_zplane_tot(1) = sigma_avg_zplane(1) !Eshwan
  endif

  call mpi_gatherv( z(2), sendsize, MPI_RPREC, &
       z_tot(2), recvsize, recvstride, MPI_RPREC, &
       rank_of_coord(0), comm, ierr)

  call mpi_gatherv( zw(2), sendsize, MPI_RPREC, &
       zw_tot(2), recvsize, recvstride, MPI_RPREC, &
       rank_of_coord(0), comm, ierr)
  
  call mpi_gatherv( rs_zplane(2), sendsize, MPI_RS, &
       rs_zplane_tot(2), recvsize, recvstride, MPI_RS, &
       rank_of_coord(0), comm, ierr)

  call mpi_gatherv( theta_stat_zplane(2), sendsize, MPI_THETA_STAT, &
       theta_stat_zplane_tot(2), recvsize, recvstride, MPI_THETA_STAT, &
       rank_of_coord(0), comm, ierr) !Eshwan

  call mpi_gatherv( sal_stat_zplane(2), sendsize, MPI_SAL_STAT, &
       sal_stat_zplane_tot(2), recvsize, recvstride, MPI_SAL_STAT, &
       rank_of_coord(0), comm, ierr) !Eshwan

  call mpi_gatherv( cnpy_zplane(2), sendsize, MPI_CNPY, &
       cnpy_zplane_tot(2), recvsize, recvstride, MPI_CNPY, &
       rank_of_coord(0), comm, ierr)

  call mpi_gatherv( tavg_zplane(2), sendsize, MPI_TAVG, &
       tavg_zplane_tot(2), recvsize, recvstride, MPI_TAVG, &
       rank_of_coord(0), comm, ierr)

  call mpi_gatherv( theta_avg_zplane(2), sendsize, MPI_THETA_AVG, &
       theta_avg_zplane_tot(2), recvsize, recvstride, MPI_THETA_AVG, &
       rank_of_coord(0), comm, ierr) !Eshwan

  call mpi_gatherv( sal_avg_zplane(2), sendsize, MPI_SAL_AVG, &
       sal_avg_zplane_tot(2), recvsize, recvstride, MPI_SAL_AVG, &
       rank_of_coord(0), comm, ierr) !Eshwan

  call mpi_gatherv( tavg_sgs_zplane(2), sendsize, MPI_SGS, &
       tavg_sgs_zplane_tot(2), recvsize, recvstride, MPI_SGS, &
       rank_of_coord(0), comm, ierr)

  call mpi_gatherv( tavg_scalar_sgs_zplane(2), sendsize, MPI_SCALAR_SGS, &
       tavg_scalar_sgs_zplane_tot(2), recvsize, recvstride, MPI_SCALAR_SGS, &
       rank_of_coord(0), comm, ierr)

  call mpi_gatherv( sigma_avg_zplane(2), sendsize, MPI_SIGMA, &
       sigma_avg_zplane_tot(2), recvsize, recvstride, MPI_SIGMA, &
       rank_of_coord(0), comm, ierr)

  deallocate(recvsize, recvstride)

  call MPI_Type_free (MPI_RS, ierr)
  call MPI_Type_free (MPI_THETA_STAT, ierr)    !Eshwan
  call MPI_Type_free (MPI_SAL_STAT, ierr)    !Eshwan
  call MPI_Type_free (MPI_CNPY, ierr)
  call mpi_type_free (MPI_TAVG, ierr)    
  call MPI_Type_free (MPI_THETA_AVG, ierr)   !Eshwan
  call MPI_Type_free (MPI_SAL_AVG, ierr)   !Eshwan
  call mpi_type_free (MPI_SGS, ierr)    
  call mpi_type_free (MPI_SCALAR_SGS, ierr)    
  call mpi_type_free (MPI_SIGMA, ierr) !Eshwan
  ! Write reconstructed data only if bottom processor
  if( coord == 0 ) then
  
    call write_tecplot_header_ND(fname_vel_zplane, 'rewind', 7, (/ Nz_tot /), &
      '"z", "<u>","<v>","<w>", "<u_uv>","<v_uv>","<w_uv>"', numtostr(coord,6), 2) !Eshwan
    call write_real_data_1D(fname_vel_zplane, 'append', 'formatted', 6, Nz_tot, &
      (/ tavg_zplane_tot % u, tavg_zplane_tot % v, tavg_zplane_tot % w, & 
      tavg_zplane_tot % u_uv, tavg_zplane_tot % v_uv, tavg_zplane_tot % w_uv &
      /), 0, zw_tot) !Eshwan

    call write_tecplot_header_ND(fname_vel2_zplane, 'rewind', 7, (/ Nz_tot/), &
      '"z", "<u<sup>2</sup>>","<v<sup>2</sup>>","<w<sup>2</sup>>", "<uw>", "<vw>", "<uv>"', &
      numtostr(coord,6), 2)
    call write_real_data_1D(fname_vel2_zplane, 'append', 'formatted', 6, Nz_tot, &
      (/ tavg_zplane_tot % u2, tavg_zplane_tot % v2, tavg_zplane_tot % w2, &
      tavg_zplane_tot % uw, tavg_zplane_tot % vw, tavg_zplane_tot % uv /), &
      0, zw_tot) 
  
    call write_tecplot_header_ND(fname_ddz_zplane, 'rewind', 4, (/ Nz_tot/), &
      '"z", "<dudz>","<dvdz>", "<dudz2>"', numtostr(coord,6), 2)
    call write_real_data_1D(fname_ddz_zplane, 'append', 'formatted', 3, Nz_tot, &
      (/ tavg_zplane_tot % dudz, tavg_zplane_tot % dvdz, tavg_zplane_tot % dudz2 /), 0, zw_tot)
  
    call write_tecplot_header_ND(fname_tau_zplane, 'rewind', 7, (/Nz_tot/), &
      '"z", "<t<sub>xx</sub>>","<t<sub>xy</sub>>","<t<sub>yy</sub>>", "<t<sub>xz</sub>>", "<t<sub>yz</sub>>", "<t<sub>zz</sub>>"', numtostr(coord,6), 2)  
    call write_real_data_1D(fname_tau_zplane, 'append', 'formatted', 6, Nz_tot, &
      (/ tavg_zplane_tot % txx, tavg_zplane_tot % txy, tavg_zplane_tot % tyy, &
      tavg_zplane_tot % txz, tavg_zplane_tot % tyz, tavg_zplane_tot % tzz /), 0, zw_tot) 
 
    call write_tecplot_header_ND(fname_divt_zplane, 'rewind', 4, (/Nz_tot/), &
      '"z", "<divtx>", "<divty>", "<divtz>"', numtostr(coord,6), 2)  
    call write_real_data_1D(fname_divt_zplane, 'append', 'formatted', 3, Nz_tot, &
      (/ tavg_zplane_tot % divtx, tavg_zplane_tot % divty, tavg_zplane_tot % divtz &
       /), 0, zw_tot) 
 
    call write_tecplot_header_ND(fname_buoy_zplane, 'rewind', 2, (/Nz_tot/), &
      '"z", "<buoyancy>"', numtostr(coord,6), 2)  
    call write_real_data_1D(fname_buoy_zplane, 'append', 'formatted', 1, Nz_tot, &
      (/ tavg_zplane_tot % buoyancy /), 0, zw_tot) 
 
    call write_tecplot_header_ND(fname_f_zplane, 'rewind', 4, (/Nz_tot/), &
      '"z", "<f<sub>x</sub>>","<f<sub>y</sub>>","<f<sub>z</sub>>"', numtostr(coord,6), 2)
    call write_real_data_1D(fname_f_zplane, 'append', 'formatted', 3, Nz_tot, &
      (/ tavg_zplane_tot % fx, tavg_zplane_tot % fy, tavg_zplane_tot % fz /), 0, zw_tot)  
  
    call write_tecplot_header_ND(fname_rs_zplane, 'rewind', 7, (/Nz_tot/), &
      '"z", "<upup>","<vpvp>","<wpwp>", "<upwp>", "<vpwp>", "<upvp>"', numtostr(coord,6), 2)  
    call write_real_data_1D(fname_rs_zplane, 'append', 'formatted', 6, Nz_tot, &
      (/ rs_zplane_tot % up2, rs_zplane_tot%vp2, rs_zplane_tot%wp2, &
      rs_zplane_tot%upwp, rs_zplane_tot%vpwp, rs_zplane_tot%upvp /), 0, zw_tot)    

    call write_tecplot_header_ND(fname_theta_stat_zplane, 'rewind', 4, (/Nz_tot/), &
      '"z", "<wpthetap>","<thetap2>","<total_tflux>"', numtostr(coord,6), 2) !Eshwan
    call write_real_data_1D(fname_theta_stat_zplane, 'append', 'formatted', 3, Nz_tot, &
      (/ theta_stat_zplane_tot % wpthetap, theta_stat_zplane_tot%thetap2, &
      theta_stat_zplane_tot % total_tflux /), 0, zw_tot) !Eshwan

    call write_tecplot_header_ND(fname_sal_stat_zplane, 'rewind', 4, (/Nz_tot/), &
      '"z", "<wpsalp>","<salp2>","<total_sflux>"', numtostr(coord,6), 2) !Eshwan
    call write_real_data_1D(fname_sal_stat_zplane, 'append', 'formatted', 3, Nz_tot, &
      (/ sal_stat_zplane_tot % wpsalp, sal_stat_zplane_tot % salp2, & 
      sal_stat_zplane_tot % total_sflux /), 0, zw_tot) !Eshwan

    call write_tecplot_header_ND(fname_sigma_zplane, 'rewind', 3, (/Nz_tot/), &
      '"z", "<sigma_theta>","<sigma_sal>"', numtostr(coord,6), 2) !Eshwan
    call write_real_data_1D(fname_sigma_zplane, 'append', 'formatted', 2, Nz_tot, &
      (/ sigma_avg_zplane_tot % sigma_theta, sigma_avg_zplane_tot % sigma_sal /), 0, zw_tot) !Eshwan
 
    call write_tecplot_header_ND(fname_cnpy_zplane, 'rewind', 7, (/Nz_tot/), &
      '"z", "<upup>","<vpvp>","<wpwp>", "<upwp>", "<vpwp>", "<upvp>"', numtostr(coord,6), 2)  
    call write_real_data_1D(fname_cnpy_zplane, 'append', 'formatted', 6, Nz_tot, &
      (/ cnpy_zplane_tot % up2, cnpy_zplane_tot%vp2, cnpy_zplane_tot%wp2, &
      cnpy_zplane_tot%upwp, cnpy_zplane_tot%vpwp, cnpy_zplane_tot%upvp /), 0, zw_tot)         
      
    call write_tecplot_header_ND(fname_cs_zplane, 'rewind', 2, (/Nz_tot/), &
      '"z", "<cs2>"', numtostr(coord,6), 2)
    call write_real_data_1D(fname_cs_zplane, 'append', 'formatted', 1, Nz_tot, &
      (/ tavg_zplane_tot % cs_opt2 /), 0, zw_tot)       

    call write_tecplot_header_ND(fname_theta_zplane, 'rewind', 9, (/ Nz_tot /), &
      ' "z", "<theta>", "<theta_w>", "<wtheta>", "<theta2>", "<dTdx>", "<dTdy>", "<dTdz>", "<sgs_tflux>"', numtostr(coord,6), 2) !Eshwan
    call write_real_data_1D(fname_theta_zplane, 'append', 'formatted', 8, Nz_tot, &
      (/ theta_avg_zplane_tot % theta, theta_avg_zplane_tot % theta_w, & 
      theta_avg_zplane_tot % wtheta, theta_avg_zplane_tot % theta2, &
      theta_avg_zplane_tot % dTdx, theta_avg_zplane_tot % dTdy, &
      theta_avg_zplane_tot % dTdz, theta_avg_zplane_tot % sgs_tflux /), 0, zw_tot) !Eshwan                   
  
    call write_tecplot_header_ND(fname_sal_zplane, 'rewind', 9, (/ Nz_tot /), &
      ' "z", "<sal>", "<sal_w>", "<wsal>", "<sal2>", "<dSdx>", "<dSdy>", "<dSdz>", "<sgs_sflux>"', numtostr(coord,6), 2) !Eshwan
    call write_real_data_1D(fname_sal_zplane, 'append', 'formatted', 8, Nz_tot, &
      (/ sal_avg_zplane_tot % sal, sal_avg_zplane_tot % sal_w, & 
      sal_avg_zplane_tot % wsal, sal_avg_zplane_tot % sal2, &
      sal_avg_zplane_tot % dSdx, sal_avg_zplane_tot % dSdy, &
      sal_avg_zplane_tot % dSdz, sal_avg_zplane_tot % sgs_sflux /), 0, zw_tot) !Eshwan                   
 
    call write_tecplot_header_ND(fname_sgs_Nu_zplane, 'rewind', 2, (/ Nz_tot /), &
      ' "z", "<Nu_t>" ', numtostr(coord,6), 2) !Eshwan
    call write_real_data_1D(fname_sgs_Nu_zplane, 'append', 'formatted', 1, Nz_tot, &
      (/ tavg_zplane_tot % Nu_t /), 0, zw_tot) !Eshwan                  

    call write_tecplot_header_ND(fname_sgs_Fsub_zplane, 'rewind', 8, (/ Nz_tot /), &
       '"z", "<Tn>", "<Beta>", "<F_LM>",  "<F_MM>", "<F_QN>", "<F_NN>", "<cs2_clips>"', numtostr(coord,6), 2) !Eshwan
    call write_real_data_1D(fname_sgs_Fsub_zplane, 'append', 'formatted', 7, Nz_tot, &
      (/ tavg_sgs_zplane_tot % Tn, tavg_sgs_zplane_tot % Beta, &
         tavg_sgs_zplane_tot % F_LM, tavg_sgs_zplane_tot % F_MM, &
         tavg_sgs_zplane_tot % F_QN, tavg_sgs_zplane_tot % F_NN, & 
         tavg_sgs_zplane_tot % cs2_clips /), 0, zw_tot) !Eshwan                  

    call write_tecplot_header_ND(fname_scalar_sgs_zplane, 'rewind', 5, (/ Nz_tot /), &
       '"z", "<ds2_t>", "<ds2_s>",  "<kappa_tt>", "<kappa_ts>"', numtostr(coord,6), 2) !Eshwan
    call write_real_data_1D(fname_scalar_sgs_zplane, 'append', 'formatted', 4, Nz_tot, &
      (/ tavg_scalar_sgs_zplane_tot % Ds_opt2_t, tavg_scalar_sgs_zplane_tot % Ds_opt2_s, &
         tavg_scalar_sgs_zplane_tot % kappa_tt, tavg_scalar_sgs_zplane_tot % kappa_ts /), 0, zw_tot) !Eshwan                  

    call write_tecplot_header_ND(fname_sgs_Isub_t_zplane, 'rewind', 8, (/ Nz_tot /), &
       '"z",  "<s_Tn_t>", "<s_Beta_t>", "<I_LM_t>", "<I_MM_t>", "<I_QN_t>", "<I_NN_t>", "<ds2_clips_t>"', numtostr(coord,6), 2) !Eshwan
    call write_real_data_1D(fname_sgs_Isub_t_zplane, 'append', 'formatted', 7, Nz_tot, &
      (/ tavg_scalar_sgs_zplane_tot % s_Tn_t, tavg_scalar_sgs_zplane_tot % s_Beta_t, &
         tavg_scalar_sgs_zplane_tot % I_LM_t, tavg_scalar_sgs_zplane_tot % I_MM_t, &
         tavg_scalar_sgs_zplane_tot % I_QN_t, tavg_scalar_sgs_zplane_tot % I_NN_t, &
         tavg_scalar_sgs_zplane_tot % ds2_clips_t /), 0, zw_tot) !Eshwan                  
 
    call write_tecplot_header_ND(fname_sgs_Isub_s_zplane, 'rewind', 8, (/ Nz_tot /), &
       '"z",  "<s_Tn_s>", "<s_Beta_s>", "<I_LM_s>", "<I_MM_s>", "<I_QN_s>", "<I_NN_s>", "<ds2_clips_s>"', numtostr(coord,6), 2) !Eshwan
    call write_real_data_1D(fname_sgs_Isub_s_zplane, 'append', 'formatted', 7, Nz_tot, &
      (/ tavg_scalar_sgs_zplane_tot % s_Tn_s, tavg_scalar_sgs_zplane_tot % s_Beta_s, &
         tavg_scalar_sgs_zplane_tot % I_LM_s, tavg_scalar_sgs_zplane_tot % I_MM_s, &
         tavg_scalar_sgs_zplane_tot % I_QN_s, tavg_scalar_sgs_zplane_tot % I_NN_s, &
         tavg_scalar_sgs_zplane_tot % ds2_clips_s /), 0, zw_tot) !Eshwan                  

    deallocate(z_tot, zw_tot)      
    deallocate(tavg_zplane_tot)
    deallocate(rs_zplane_tot)
    deallocate(cnpy_zplane_tot)
    deallocate(theta_avg_zplane_tot) !Eshwan
    deallocate(theta_stat_zplane_tot) !Eshwan
    deallocate(sal_avg_zplane_tot) !Eshwan
    deallocate(sal_stat_zplane_tot) !Eshwan
    deallocate(tavg_sgs_zplane_tot) !Eshwan
    deallocate(tavg_scalar_sgs_zplane_tot) !Eshwan
    deallocate(sigma_avg_zplane_tot) !Eshwan
  endif

$else
call write_tecplot_header_ND(fname_vel_zplane, 'rewind', 7, (/ Nz /), &
   '"z", "<u>","<v>","<w>", "<u_uv>","<v_uv>","<w_uv>"', numtostr(coord,6), 2) !Eshwan
call write_real_data_1D(fname_vel_zplane, 'append', 'formatted', 6, nz, &
  (/ tavg_zplane % u, tavg_zplane % v, tavg_zplane % w, & 
  tavg_zplane % u_uv, tavg_zplane % v_uv, tavg_zplane % w_uv &
  /), 0, zw(1:nz)) !Eshwan

call write_tecplot_header_ND(fname_vel2_zplane, 'rewind', 7, (/ Nz/), &
   '"z", "<u<sup>2</sup>>","<v<sup>2</sup>>","<w<sup>2</sup>>", "<uw>", "<vw>", "<uv>"', &
   numtostr(coord,6), 2)
call write_real_data_1D(fname_vel2_zplane, 'append', 'formatted', 6, nz, &
  (/ tavg_zplane % u2, tavg_zplane % v2, tavg_zplane % w2, &
  tavg_zplane % uw, tavg_zplane % vw, tavg_zplane % uv /), 0, zw(1:nz)) 
  
call write_tecplot_header_ND(fname_ddz_zplane, 'rewind', 4, (/ Nz/), &
   '"z", "<dudz>","<dvdz>", "<dudz2>"', numtostr(coord,6), 2)
call write_real_data_1D(fname_ddz_zplane, 'append', 'formatted', 3, nz, &
  (/ tavg_zplane % dudz, tavg_zplane % dvdz, tavg_zplane % dudz2 /), 0, zw(1:nz))
  
call write_tecplot_header_ND(fname_tau_zplane, 'rewind', 7, (/Nz/), &
   '"z", "<t<sub>xx</sub>>","<t<sub>xy</sub>>","<t<sub>yy</sub>>", "<t<sub>xz</sub>>", "<t<sub>yz</sub>>", "<t<sub>zz</sub>>"', &
   numtostr(coord,6), 2)  
call write_real_data_1D(fname_tau_zplane, 'append', 'formatted', 6, nz, &
  (/ tavg_zplane % txx, tavg_zplane % txy, tavg_zplane % tyy, &
  tavg_zplane % txz, tavg_zplane % tyz, tavg_zplane % tzz /), 0, zw(1:nz)) 

call write_tecplot_header_ND(fname_divt_zplane, 'rewind', 4, (/Nz/), &
   '"z", "<divtx>", "<divty>", "<divtz>"', numtostr(coord,6), 2)  
call write_real_data_1D(fname_divt_zplane, 'append', 'formatted', 3, nz, &
  (/ tavg_zplane % divtx, tavg_zplane % divty, tavg_zplane % divtz /), 0, zw(1:nz)) 
 
call write_tecplot_header_ND(fname_buoy_zplane, 'rewind', 2, (/Nz/), &
  '"z", "<buoyancy>"', numtostr(coord,6), 2)  
call write_real_data_1D(fname_buoy_zplane, 'append', 'formatted', 1, nz, &
  (/ tavg_zplane % buoyancy /), 0, zw(1:nz)) 
  
call write_tecplot_header_ND(fname_f_zplane, 'rewind', 4, (/Nz/), &
   '"z", "<f<sub>x</sub>>","<f<sub>y</sub>>","<f<sub>z</sub>>"', numtostr(coord,6), 2)
call write_real_data_1D(fname_f_zplane, 'append', 'formatted', 3, nz, &
  (/ tavg_zplane % fx, tavg_zplane % fy, tavg_zplane % fz /), 0, zw(1:nz))  
  
call write_tecplot_header_ND(fname_rs_zplane, 'rewind', 7, (/Nz/), &
   '"z", "<upup>","<vpvp>","<wpwp>", "<upwp>", "<vpwp>", "<upvp>"', numtostr(coord,6), 2)  
call write_real_data_1D(fname_rs_zplane, 'append', 'formatted', 6, nz, &
  (/ rs_zplane % up2, rs_zplane%vp2, rs_zplane%wp2, &
  rs_zplane%upwp, rs_zplane%vpwp, rs_zplane%upvp /), 0, zw(1:nz))

call write_tecplot_header_ND(fname_theta_stat_zplane, 'rewind', 4, (/Nz/), &
   '"z", "<wpthetap>","<thetap2>","<total_tflux>"', numtostr(coord,6), 2) !Eshwan
call write_real_data_1D(fname_theta_stat_zplane, 'append', 'formatted', 3, nz, &
  (/ theta_stat_zplane % wpthetap, theta_stat_zplane%thetap2, &
  theta_stat_zplane % total_tflux /), 0, zw(1:nz)) !Eshwan

call write_tecplot_header_ND(fname_sal_stat_zplane, 'rewind', 4, (/Nz/), &
   '"z", "<wpsalp>","<salp2>","<total_sflux>"', numtostr(coord,6), 2) !Eshwan
call write_real_data_1D(fname_sal_stat_zplane, 'append', 'formatted', 3, nz, &
  (/ sal_stat_zplane % wpsalp, sal_stat_zplane%salp2, &
  sal_stat_zplane % total_sflux /), 0, zw(1:nz)) !Eshwan

call write_tecplot_header_ND(fname_sigma_zplane, 'rewind', 3, (/Nz/), &
   '"z", "<sigma_theta>","<sigma_sal>"', numtostr(coord,6), 2) !Eshwan
call write_real_data_1D(fname_sigma_zplane, 'append', 'formatted', 2, nz, &
  (/ sigma_avg_zplane_tot % sigma_theta, sigma_avg_zplane_tot % sigma_sal /), 0, zw(1:nz)) !Eshwan

call write_tecplot_header_ND(fname_cnpy_zplane, 'rewind', 7, (/Nz/), &
   '"z", "<upup>","<vpvp>","<wpwp>", "<upwp>", "<vpwp>", "<upvp>"', numtostr(coord,6), 2)  
call write_real_data_1D(fname_cnpy_zplane, 'append', 'formatted', 6, nz, &
  (/ cnpy_zplane % up2, cnpy_zplane%vp2, cnpy_zplane%wp2, &
  cnpy_zplane%upwp, cnpy_zplane%vpwp, cnpy_zplane%upvp /), 0, zw(1:nz))  
  
call write_tecplot_header_ND(fname_cs_zplane, 'rewind', 2, (/Nz/), &
   '"z", "<cs2>"', numtostr(coord,6), 2)
call write_real_data_1D(fname_cs_zplane, 'append', 'formatted', 1, nz, &
  (/ tavg_zplane % cs_opt2 /), 0, zw(1:nz))    

call write_tecplot_header_ND(fname_theta_zplane, 'rewind', 9, (/ Nz /), &
   ' "z", "<theta>", "<theta_w>", "<wtheta>", "<theta2>", "<dTdx>", "<dTdy>", "<dTdz>", "<sgs_tflux>"', numtostr(coord,6), 2) !Eshwan
call write_real_data_1D(fname_theta_zplane, 'append', 'formatted', 8, nz, &
  (/ theta_avg_zplane % theta, theta_avg_zplane % theta_w, &
     theta_avg_zplane % wtheta, theta_avg_zplane % theta2, &
     theta_avg_zplane % dTdx, theta_avg_zplane % dTdy, &
     theta_avg_zplane % dTdz, theta_avg_zplane % sgs_tflux /), 0, zw(1:nz)) !Eshwan                   

call write_tecplot_header_ND(fname_sal_zplane, 'rewind', 9, (/ Nz /), &
   ' "z", "<sal>", "<sal_w>", "<wsal>", "<sal2>", "<dSdx>", "<dSdy>", "<dSdz>", "<sgs_sflux>"', numtostr(coord,6), 2) !Eshwan
call write_real_data_1D(fname_sal_zplane, 'append', 'formatted', 8, nz, &
  (/ sal_avg_zplane % sal, sal_avg_zplane % sal_w, &
     sal_avg_zplane % wsal, sal_avg_zplane % sal2, &
     sal_avg_zplane % dSdx, sal_avg_zplane % dSdy, &
     sal_avg_zplane % dSdz, sal_avg_zplane % sgs_slux /), 0, zw(1:nz)) !Eshwan                   

call write_tecplot_header_ND(fname_sgs_Nu_zplane, 'rewind', 2, (/ Nz /), &
   ' "z", "<Nu_t>"', numtostr(coord,6), 2) !Eshwan
call write_real_data_1D(fname_sgs_Nu_zplane, 'append', 'formatted', 1, Nz, &
  (/ tavg_zplane % Nu_t /), 0, zw(1:nz)) !Eshwan                  

call write_tecplot_header_ND(fname_sgs_Fsub_zplane, 'rewind', 8, (/ Nz /), &
   ' "z", "<Tn>", "<Beta>", "<F_LM>",  "<F_MM>", "<F_QN>", "<F_NN>", "<cs2_clips>"', numtostr(coord,6), 2) !Eshwan
call write_real_data_1D(fname_sgs_Fsub_zplane, 'append', 'formatted', 7, nz, &
   (/ tavg_sgs_zplane % Tn, tavg_sgs_zplane % Beta, &
      tavg_sgs_zplane % F_LM, tavg_sgs_zplane % F_MM, &
      tavg_sgs_zplane % F_QN, tavg_sgs_zplane % F_NN, &
      tavg_sgs_zplane % cs2_clips /), 0, zw(1:nz)) !Eshwan                  

call write_tecplot_header_ND(fname_scalar_sgs_zplane, 'rewind', 5, (/ Nz /), &
   ' "z", "<ds2_t>", "<ds2_s>",  "<kappa_tt>", "<kappa_ts>"', numtostr(coord,6), 2) !Eshwan
call write_real_data_1D(fname_scalar_sgs_zplane, 'append', 'formatted', 4, nz, &
   (/ tavg_scalar_sgs_zplane % Ds_opt2_t, tavg_scalar_sgs_zplane % Ds_opt2_s, &
      tavg_scalar_sgs_zplane % kappa_tt, tavg_scalar_sgs_zplane % kappa_ts /), 0, zw(1:nz)) !Eshwan                  

call write_tecplot_header_ND(fname_sgs_Isub_t_zplane, 'rewind', 8, (/ Nz /), &
   ' "z",  "<s_Tn_t>", "<s_Beta_t>", "<I_LM_t>", "<I_MM_t>", "<I_QN_t>", "<I_NN_t>", "<ds2_clips_t>"', numtostr(coord,6), 2) !Eshwan
call write_real_data_1D(fname_sgs_Isub_t_zplane, 'append', 'formatted', 7, nz, &
   (/ tavg_scalar_sgs_zplane % s_Tn_t, tavg_scalar_sgs_zplane % s_Beta_t, &
      tavg_scalar_sgs_zplane % I_LM_t, tavg_scalar_sgs_zplane % I_MM_t, &
      tavg_scalar_sgs_zplane % I_QN_t, tavg_scalar_sgs_zplane % I_NN_t, &
      tavg_scalar_sgs_zplane % ds2_clips_t /), 0, zw(1:nz)) !Eshwan                  

call write_tecplot_header_ND(fname_sgs_Isub_s_zplane, 'rewind', 8, (/ Nz /), &
   ' "z",  "<s_Tn_s>", "<s_Beta_s>", "<I_LM_s>", "<I_MM_s>", "<I_QN_s>", "<I_NN_s>", "<ds2_clips_s>"', numtostr(coord,6), 2) !Eshwan
call write_real_data_1D(fname_sgs_Isub_s_zplane, 'append', 'formatted', 7, nz, &
   (/ tavg_scalar_sgs_zplane % s_Tn_s, tavg_scalar_sgs_zplane % s_Beta_s, &
      tavg_scalar_sgs_zplane % I_LM_s, tavg_scalar_sgs_zplane % I_MM_s, &
      tavg_scalar_sgs_zplane % I_QN_s, tavg_scalar_sgs_zplane % I_NN_s, &
      tavg_scalar_sgs_zplane % ds2_clips_s /), 0, zw(1:nz)) !Eshwan                  

$endif

deallocate(tavg, tavg_zplane, rs, rs_zplane, cnpy_zplane)
deallocate(theta_stat_zplane)  !Eshwan
deallocate(theta_avg)
deallocate(theta_avg_zplane) !Eshwan
deallocate(sal_stat_zplane)  !Eshwan
deallocate(sal_avg)
deallocate(sal_avg_zplane) !Eshwan
deallocate(sigma_avg_zplane)

!$if($OUTPUT_EXTRA)
deallocate(tavg_sgs)
deallocate(tavg_sgs_zplane, tavg_scalar_sgs_zplane) !Eshwan
!$endif

nullify(x,y,z,zw)

$if($MPI)
! Ensure all writes complete before preceeding
call mpi_barrier( comm, ierr )
$endif

return
end subroutine tavg_finalize

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine tavg_checkpoint()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! This subroutine writes the restart data and is to be called by 'checkpoint'
! for intermediate checkpoints and by 'tavg_finalize' at the end of the
! simulation.
!
use param, only : checkpoint_tavg_file
use stat_defs, only : tavg_total_time, tavg
!$if($OUTPUT_EXTRA)
use param, only : checkpoint_tavg_sgs_file
use stat_defs, only : tavg_total_time_sgs, tavg_sgs
!$endif
implicit none

character(*), parameter :: sub_name = mod_name // '.tavg_checkpoint'

character(64) :: fname
logical :: opn

fname = checkpoint_tavg_file
$if($MPI)
call string_concat( fname, '.c', coord)
$endif

!  Write data to tavg.out
inquire (unit=1, opened=opn)
if (opn) call error (sub_name, 'unit 1 already open')

$if ($WRITE_BIG_ENDIAN)
open (1, file=fname, action='write', position='rewind', &
  form='unformatted', convert='big_endian')
$elseif ($WRITE_LITTLE_ENDIAN)
open (1, file=fname, action='write', position='rewind', &
  form='unformatted', convert='little_endian')
$else
open (1, file=fname, action='write', position='rewind', form='unformatted')
$endif

! write the entire structures
write (1) tavg_total_time
write (1) tavg
close(1)

!----
!$if($OUTPUT_EXTRA)
fname = checkpoint_tavg_sgs_file
  $if($MPI)
  call string_concat( fname, '.c', coord)
  $endif

  !  Write data to tavg_sgs.out
  $if ($WRITE_BIG_ENDIAN)
  open (1, file=fname, action='write', position='rewind', &
       form='unformatted', convert='big_endian')
  $elseif ($WRITE_LITTLE_ENDIAN)
  open (1, file=fname, action='write', position='rewind', &
       form='unformatted', convert='little_endian')
  $else
  open (1, file=fname, action='write', position='rewind', form='unformatted')
  $endif

  ! write the entire structures
  write (1) tavg_total_time_sgs
  write (1) tavg_sgs
  close(1)
!$endif

return
end subroutine tavg_checkpoint

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine spectra_init()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use types, only : rprec
use param, only : path, checkpoint_spectra_file
use messages
use param, only : coord, dt, spectra_nloc, lh, nx
use stat_defs, only : spectra, spectra_total_time, spectra_initialized,spectra_dt
implicit none

character (*), parameter :: sub_name = mod_name // '.spectra_init'
$if ($MPI)
character (*), parameter :: MPI_suffix = '.c'
$endif
character (128) :: fname

logical :: opn, exst
integer :: k

inquire (unit=1, opened=opn)
if (opn) call error (sub_name, 'unit 1 already open')

fname = checkpoint_spectra_file
$if ($MPI)
call string_concat( fname, MPI_suffix, coord )
$endif

inquire (file=fname, exist=exst)
if (.not. exst) then
   !  Nothing to read in
   if (coord == 0) then
      write(*,*) ' '
      write(*,*)'No previous spectra data - starting from scratch.'
   endif

   spectra_total_time = 0._rprec
   ! note: spectra already initialized during output_init

else


   $if ($READ_BIG_ENDIAN)
   open (1, file=fname, action='read', position='rewind',  &
        form='unformatted', convert='big_endian')
   $elseif ($READ_LITTLE_ENDIAN)
   open (1, file=fname, action='read', position='rewind',  &
        form='unformatted', convert='little_endian')
   $else
   open (1, file=fname, action='read', position='rewind',  &
        form='unformatted')
   $endif

   read (1) spectra_total_time
   do k=1, spectra_nloc
      read (1) spectra(k) % power
   enddo

   close(1)

endif

! Initialize spectra_dt
spectra_dt = 0.0_rprec

! Set global switch that spectra as been initialized
spectra_initialized = .true.

return
end subroutine spectra_init

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine spectra_compute()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use types, only : rprec
use sim_param, only : u
use param, only : Nx, Ny, dt, dz, lh
use param, only : spectra_nloc, spectra_nskip
use stat_defs, only : spectra, spectra_total_time, spectra_dt
use functions, only : linear_interp
use fft, only : forw_spectra
implicit none

integer :: i, j, k
real(rprec), allocatable, dimension(:) :: ui, uhat, power

real(rprec) :: nxny

! Interpolation variable
allocate(ui(nx), uhat(nx), power(lh))


nxny = real(nx*ny,rprec)

!  Loop over all spectra locations
do k=1, spectra_nloc
  
  $if ($MPI)
  if(spectra(k) % coord == coord) then
  $endif

  do j=1,Ny

    do i=1,Nx
        ! Interpolate to the specified plane
        ui(i) = linear_interp(u(i,j,spectra(k) % istart),u(i,j,spectra(k) % istart+1), &
                              dz, spectra(k) % ldiff)

    enddo

    ! 1) Compute uhat for the given j
    ui = ui - sum(ui) / Nx ! Remove the mean
    ! Compute FFT
    call rfftw_f77_one(forw_spectra, ui, uhat)
    !  Normalize
    uhat = uhat / Nx

    ! 2) Compute power spectra for given j
    power(1) = uhat(1)**2
    do i=2,lh-1
      power(i) = uhat(i)**2 + uhat(Nx-i+2)**2
    enddo
    power(lh) = uhat(lh)**2 ! Nyquist

    ! Sum jth component, normalize, and include averaging over y 
    spectra(k) % power = spectra(k) % power + spectra_dt * power / nxny

  enddo

  $if ($MPI)
  endif
  $endif

enddo  

!  Update spectra total time
spectra_total_time = spectra_total_time + spectra_dt
  
deallocate(ui, uhat, power)

! Set spectra_dt back to zero for next use
spectra_dt = 0.0_rprec

return
end subroutine spectra_compute

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine spectra_finalize()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use types, only : rprec
use param, only : path
use param, only : lh, spectra_nloc, spectra_loc
$if($MPI)
use param, only : comm, ierr
$endif
use fft, only : kx
use stat_defs, only : spectra, spectra_total_time
implicit none

include 'tecryte.h'

character (*), parameter :: sub_name = mod_name // '.spectra_finalize'
!character(25) :: cl
character (64) :: fname

integer ::  k

! Perform final checkpoint 
call spectra_checkpoint()

!  Loop over all zplane locations
do k=1,spectra_nloc
  
  $if ($MPI)
  if(spectra(k) % coord == coord) then
  $endif

  ! Finalize averaging for power spectra (averaging over time)
  spectra(k) % power = spectra(k) % power / spectra_total_time

  !  Create unique file name
  call string_splice( fname, path // 'output/spectra.z-', spectra_loc(k),'.dat' )

  !  Omitting Nyquist from output
  call write_tecplot_header_ND(fname, 'rewind', 2, (/ lh-1/), &
    '"k", "E(k)"', numtostr(k, 6), 2 ) 
  call write_real_data_1D(fname, 'append', 'formatted', 1, lh-1, &
    (/ spectra(k) % power(1:lh-1) /), 0, (/ kx(1:lh-1,1) /))
    
  $if ($MPI)
  endif
  $endif

enddo  
  
deallocate(spectra)

return

end subroutine spectra_finalize

!*******************************************************************************
subroutine spectra_checkpoint()
!*******************************************************************************
!
! This subroutine writes the restart data and is to be called by 'checkpoint'
! for intermediate checkpoints and by 'spectra_finalize' at the end of the 
! simulation.    
! 
use param, only : checkpoint_spectra_file, spectra_nloc
$if($MPI)
use param, only : comm, ierr
$endif
use stat_defs, only : spectra, spectra_total_time

implicit none

character (*), parameter :: sub_name = mod_name // '.spectra_checkpoint'
character (64) :: fname

logical :: opn
integer :: k

fname = checkpoint_spectra_file
$if($MPI)
call string_concat( fname, '.c', coord)
$endif

!  Write data to spectra.out
inquire (unit=1, opened=opn)
if (opn) call error (sub_name, 'unit 1 already open')

$if ($WRITE_BIG_ENDIAN)
open (1, file=fname, action='write', position='rewind', &
  form='unformatted', convert='big_endian')
$elseif ($WRITE_LITTLE_ENDIAN)
open (1, file=fname, action='write', position='rewind', &
  form='unformatted', convert='little_endian')
$else
open (1, file=fname, action='write', position='rewind', form='unformatted')
$endif

! write the accumulated time and power
write (1) spectra_total_time
do k=1, spectra_nloc
  write (1) spectra(k) % power
enddo
close(1)

$if($MPI)
! Ensure all writes complete before preceeding
call mpi_barrier( comm, ierr )
$endif

return
end subroutine spectra_checkpoint

end module io
