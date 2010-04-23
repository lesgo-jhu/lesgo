module io
use types,only:rprec
use param, only : ld, nx, ny, nz, nz_tot, write_inflow_file, path,  &
                  USE_MPI, coord, rank, nproc, jt_total, total_time, total_time_dim
use sim_param, only : w, dudz
use messages
use strmod

implicit none

save
private

!!$public openfiles,output_loop,output_final,                   &
!!$     inflow_write, avg_stats
public jt_total, openfiles, inflow_read, inflow_write, output_loop, output_final
public mean_u,mean_u2,mean_v,mean_v2,mean_w,mean_w2
public w_uv, dudz_uv, w_uv_tag, dudz_uv_tag, interp_to_uv_grid, stats_init
public write_tecplot_header_xyline, write_tecplot_header_ND
public write_real_data, write_real_data_1D, write_real_data_2D, write_real_data_3D

!interface write_real_data_ND
!  module procedure write_real_data_1D, write_real_data_2D, write_real_data_3D
!end interface

!!$ Region commented by JSG 
!!$integer,parameter::base=2000,nwrite=base
!!$
logical, parameter :: cumulative_time = .true.
character (*), parameter :: fcumulative_time = path // 'total_time.dat'

character (*), parameter :: mod_name = 'io'
!!$
!!$logical, parameter :: use_avg_stats = .false.
!!$integer, parameter :: n_avg_stats = 5000
!!$                      !--interval for updates in avg_stats
!!$character (*), parameter :: end_hdr_avg = '# end header'
!!$
!!$!--write velocity slices (entire planes, instantaneous)
!!$!--this is not related to avgslice
!!$!--this output to be averaged by postprocessing program
!!$logical, parameter :: use_vel_slice = .false.
!!$integer, parameter :: n_vel_slice_write = 20000
!!$
!!$!--write immersed bdry force field
!!$logical, parameter :: write_f = .false.
!!$
!!$!--write whole txz field or not (for checking momentum balance)
!!$!--does not affect avgslice writing of txz
!!$logical, parameter :: write_txz = .false.
!!$
!!$!--write 1d-slices (w/o any processing) at specified stations
!!$logical, parameter :: use_write_1dx = .true.  !--x-slices
!!$logical, parameter :: use_write_1dy = .true.  !--y-slices
!!$integer, parameter :: n_1dx = 1  !--number of x-slices
!!$integer, parameter :: n_1dy = 1  !--number of y-slices
!!$real (rprec) :: y_1dx(n_1dx) = (/ 0.5_rprec /)
!!$                !--y-value for x-slice
!!$real (rprec) :: x_1dy(n_1dy) = (/ 0.5_rprec /)
!!$                !--x-value for y-slice
!!$
!!$!!!!  io_spec=.true. output plan-averaged spectrum
!!$!!!!  time_spec>0 output time series spectrum (need additional calcu.)
!!$logical,parameter::io_spec=.false.
!!$integer,parameter::time_spec=0
!!$integer::n_obs
!!$integer,allocatable::obs_pt(:,:)
!!$
!!$!!!!  io_mean=.true. output small domain time-averaged velocity
!!$logical,parameter::io_mean=.false.
!integer,parameter::jx_pls=1,jx_ple=nx,width=ny/2-1
integer,parameter::jx_pls=1,jx_ple=1,width=1
integer,parameter::jy_pls=ny/2-width,jy_ple=ny/2+width+1
real(kind=rprec),dimension(jx_pls:jx_ple,jy_pls:jy_ple,nz):: &
     mean_u,mean_v,mean_w,mean_u2,mean_v2,mean_w2
!!$
!!$!!!!  io_lambda2
!!$!logical,parameter::io_lambda2=.false.
!!$!real(kind=rprec),dimension(nx,ny,nz)::lam2  !--commented to save mem.

real(rprec), dimension(lbound(w,1):ubound(w,1),lbound(w,2):ubound(w,2),lbound(w,3):ubound(w,3)) :: w_uv
real(rprec), dimension(lbound(dudz,1):ubound(dudz,1),lbound(dudz,2):ubound(dudz,2),lbound(dudz,3):ubound(dudz,3)) :: dudz_uv ! on the uv grid
integer :: w_uv_tag, dudz_uv_tag  
                   
!**********************************************************************
contains
!**********************************************************************

!**********************************************************************
subroutine interp_to_uv_grid(var,var_uv,tag)
!**********************************************************************
!  This function interpolates the array var, which resides on the w-grid,
!  onto the uv-grid variable var_uv using linear interpolation. It is 
!  important to note that message passing is required for MPI cases and 
!  all processors must call this routine. If this subroutine is call from a 
!  only a subset of the total processors, the code will hang due to the usage
!  of the syncronous send/recv functions and certain processors waiting
!  to recv data but it never gets there.
!
use types, only : rprec
use param,only : nz,ld,jt
use sim_param, only : w, dudz
$if ($MPI)
use mpi_defs, only : mpi_sync_real_array
$endif

implicit none

real(rprec), dimension(:,:,:), intent(IN) :: var
real(rprec), dimension(:,:,:), intent(OUT) :: var_uv
integer, intent(INOUT) :: tag
!real(rprec), dimension(2) :: var
integer :: lbx,ubx,lby,uby,lbz,ubz
integer :: i,j,k

character (*), parameter :: sub_name = mod_name // '.interp_to_uv_grid'

if(tag == jt) then
$if ($VERBOSE)
  call mesg(sub_name, 'Interpolation already performed for current time step')
$endif
  return
endif
  
lbx=lbound(var,1); ubx=ubound(var,1)
lby=lbound(var,2); uby=ubound(var,2)
lbz=lbound(var,3); ubz=ubound(var,3)

do k=lbz,ubz-1
  do j=lby,uby
    do i=lbx,ubx
      var_uv(i,j,k) = 0.5 * (var(i,j,k+1) + var(i,j,k))
    enddo
  enddo
enddo

$if ($MPI)

!  Take care of top "physical" boundary
if(coord == nproc - 1) var_uv(:,:,ubz) = var_uv(:,:,ubz-1)

!  Sync all overlapping data
call mpi_sync_real_array(var_uv)

$else

!  Take care of top "physical" boundary
var_uv(:,:,ubz) = var_uv(:,:,ubz-1)

$endif
  
!  set k=nz
!k=ubz
!do j=lby,uby; do i=lbx,ubx
!  $if ($MPI)
!!  Check for top "physical" boundary
!  if (coord == nproc - 1) then
!    var_uv(i,j,ubz) = var(i,j,ubz) 
!  else
!    var_uv(i,j,k) = 0.5*(buf(i,j) + var(i,j,k))
!  endif
!  $else
!  var_uv(i,j,k) = var(i,j,k)
!  $endif
!enddo; enddo

tag = jt ! Set identifying tag 

return 

!!$if($MPI)
!deallocate(buf)
!!$endif

end subroutine interp_to_uv_grid

!!$ Commented by JSG
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$!--write out velocity slices for later postprocessing (avging)
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$subroutine vel_slice ()
!!$use sim_param, only: u, v, w
!!$
!!$implicit none
!!$
!!$integer, parameter :: nx_skip = nx / 4  !--space between x-slices
!!$integer, parameter :: ny_skip = ny / 4  !--space between y-slices
!!$integer, parameter :: nz_skip = (nz_tot-1) / 4  !--space between z-slices
!!$
!!$character (128) :: fname
!!$character (32) :: temp
!!$
!!$integer :: kmin, kk
!!$
!!$!---------------------------------------------------------------------
!!$
!!$!--x-direction slices
!!$write (fname, '(a,i0,a,i6.6,a)') path // 'output/vel.xskip', nx_skip,  &
!!$                                   '.', jt_total, '.out'
!!$
!!$$if ($MPI)
!!$  write (temp, '(".c",i0)') coord
!!$  fname = trim (fname) // temp
!!$$endif
!!$
!!$open(1,file=fname,form='unformatted')
!!$
!!$write (1) u(1:nx:nx_skip, 1:ny, 1:nz), v(1:nx:nx_skip, 1:ny, 1:nz),  &
!!$          w(1:nx:nx_skip, 1:ny, 1:nz)
!!$
!!$close(1)
!!$
!!$!--y-direction slices
!!$write (fname, '(a,i0,a,i6.6,a)') path // 'output/vel.yskip', ny_skip,  &
!!$                                   '.', jt_total, '.out'
!!$
!!$$if ($MPI)
!!$  write (temp, '(".c",i0)') coord
!!$  fname = trim (fname) // temp
!!$$endif
!!$
!!$open(1,file=fname,form='unformatted')
!!$
!!$write (1) u(1:nx, 1:ny:ny_skip, 1:nz), v(1:nx, 1:ny:ny_skip, 1:nz),  &
!!$          w(1:nx, 1:ny:ny_skip, 1:nz)
!!$
!!$close(1)
!!$
!!$!--z-direction slices
!!$write (fname, '(a,i0,a,i6.6,a)') path // 'output/vel.zskip', nz_skip,  &
!!$                                   '.', jt_total, '.out'
!!$
!!$$if ($MPI)
!!$
!!$  !--first global k for corresponding to this process
!!$  !  kg = 1 + kk * nz_skip
!!$  kk = ceiling (real (coord * (nz - 1), rprec) / nz_skip)
!!$  !--this corresponds to local k of kmin:
!!$  kmin = 1 + kk * nz_skip - coord * (nz - 1)
!!$
!!$  write (temp, '(".c",i0)') coord
!!$  fname = trim (fname) // temp
!!$  
!!$$else
!!$
!!$  kmin = 1
!!$
!!$$endif
!!$
!!$if ((1 <= kmin) .and. (kmin <= nz-1)) then  !--no write if out of range
!!$
!!$  open (1, file=fname, form='unformatted')
!!$
!!$  write (1) u(1:nx, 1:ny, kmin:nz-1:nz_skip),  &
!!$            v(1:nx, 1:ny, kmin:nz-1:nz_skip),  &
!!$            w(1:nx, 1:ny, kmin:nz-1:nz_skip)
!!$
!!$  close (1)
!!$
!!$end if
!!$
!!$end subroutine vel_slice

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! file number 1-10 are used for temporary use
! 11-19 are basic output
! 20-40 are avgslice
! use >50 for debugging
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine openfiles()
!!$use param,only:sflux_flag
use sim_param,only:path
implicit none

!--to hold file names
character (64) :: temp
!!$character (64) :: fCS1plan, fCS2plan, fCS4plan, fVISCplan,  &
!!$                  fDISSplan, fCS1Vplan, fCS2Vplan, fCS4Vplan

integer::i

logical :: exst

if (cumulative_time) then

  inquire (file=fcumulative_time, exist=exst)
  if (exst) then
    open (1, file=fcumulative_time)
    read (1, *) jt_total, total_time, total_time_dim
    close (1)
  else  !--assume this is the first run on cumulative time
    write (*, *) 'file ', fcumulative_time, ' not found'
    write (*, *) 'assuming jt_total = 0, total_time = 0., total_time_dim = 0.'
    jt_total = 0
    total_time = 0.
    total_time_dim = 0._rprec
  end if

end if

!--see also energy () subr. for flushing this file
if ((.not. USE_MPI) .or. (USE_MPI .and. rank == 0)) then
  open(13, file=path//'output/check_ke.out', position='append')
end if

!!$if(time_spec.gt.0)then
!!$  open(15,file=path//'output/velspec.out',form='unformatted',position='append')
!!$  if(jt_total.eq.0)rewind(15)
!!$endif
!!$
!!$if(io_mean)then
!!$  open(51,file=path//'output/mean_u.out',form='unformatted',position='append')
!!$  if(jt_total.eq.0)then
!!$    rewind(51)
!!$    write(51)jx_pls,jx_ple,jy_pls,jy_ple
!!$  endif
!!$endif

!!$fCS1plan = path // 'output/CS1plan.out'
!!$fCS2plan = path // 'output/CS2plan.out'
!!$fCS4plan = path // 'output/CS4plan.out'
!!$fVISCplan = path // 'output/VISCplan.out'
!!$fDISSplan = path // 'output/DISSplan.out'
!!$fCS1Vplan = path // 'output/CS1Vplan.out'
!!$fCS2Vplan = path // 'output/CS2Vplan.out'
!!$fCS4Vplan = path // 'output/CS4Vplan.out'
!!$
!!$$if ($MPI)
!!$  !--append coordinate identifiers
!!$  write (temp, '(".c",i0)') coord
!!$  fCS1plan = trim (fCS1plan) // temp
!!$  fCS2plan = trim (fCS2plan) // temp
!!$  fCS4plan = trim (fCS4plan) // temp
!!$  fVISCplan = trim (fVISCplan) // temp
!!$  fDISSplan = trim (fDISSplan) // temp
!!$  fCS1Vplan = trim (fCS1Vplan) // temp
!!$  fCS2Vplan = trim (fCS2Vplan) // temp
!!$  fCS4Vplan = trim (fCS4Vplan) // temp
!!$$endif

!!$open (90, file=fCS1plan, form='unformatted')
!!$open (91, file=fCS2plan, form='unformatted')
!!$open (92, file=fCS4plan, form='unformatted')
!!$open (93, file=fVISCplan, form='unformatted')
!!$open (94, file=fDISSplan, form='unformatted')
!!$open (95, file=fCS1Vplan, form='unformatted')
!!$open (96, file=fCS2Vplan, form='unformatted')
!!$open (97, file=fCS4Vplan, form='unformatted')
!!$!TSif(sflux_flag)open(98,file='/home/bluesky/ytseng/RESEARCH/JHU_LES/CHANNEL/output/series_data.dat',form='unformatted')

!!$if(time_spec.gt.0)then
!!$open(1,file=path//'obs.pt')
!!$read(1,*)n_obs
!!$allocate(obs_pt(1:2,n_obs))
!!$do i=1,n_obs
!!$read(1,*)obs_pt(1:2,i)
!!$enddo
!!$close(1)
!!$endif

end subroutine openfiles

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output_loop(jt)
use param, only : dx, dy, dz, jt_total, dt_dim
use stat_defs, only: tstats_t, point_t, domain_t, yplane_t, zplane_t, zplane_avg_t
use sim_param, only : u, v, w, txx, txy, tyy, txz, tyz, tzz, dudz, dvdz
use grid_defs, only: z

!!$use param,only:output,dt,c_count,S_FLAG,SCAL_init
!!$use sim_param,only:path,u,v,w,dudz,dudx,p,&
!!$     RHSx,RHSy,RHSz,theta, txx, txy, txz, tyy, tyz, tzz
!!$use sgsmodule,only:Cs_opt2!,Cs_opt2_avg  !--not essential, save mem.
!!$use scalars_module2,only:scalar_slice
!!$use immersedbc, only : fx, fy, fz
implicit none
character(len=20)::req
integer,intent(in)::jt
!real(kind=rprec),dimension(ld,ny,nz)::usp,vsp  !--not used
real(kind=rprec),dimension(nz)::u_ndim

character (64) :: fname, temp

integer :: i,j,k
integer::jx,jy,jz

!  Determine if time summations are to be calculated
if(tstats_t%calc) then
!  Check if we are in the time interval for running summations
  if(jt >= tstats_t%nstart .and. jt <= tstats_t%nend) then
    if(.not. tstats_t%started) then
      if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
        write(*,*) '-------------------------------'   
        write(*,"(1a,i9,1a,i9)") 'Starting running time summation from ', tstats_t%nstart, ' to ', tstats_t%nend
        write(*,*) '-------------------------------'   
      endif

      call tstats_init()

    endif
!  Compute running summations
    call tstats_compute ()
  endif 
  
endif

!  Determine if instantaneous point velocities are to be recorded
if(point_t%calc) then
  if(jt >= point_t%nstart .and. jt <= point_t%nend .and. ( jt == point_t%nstart .or. mod(jt,point_t%nskip)==0) ) then
    if(.not. point_t%started) then
      if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then   
        write(*,*) '-------------------------------'   
        write(*,"(1a,i9,1a,i9)") 'Writing instantaneous point velocities from ', point_t%nstart, ' to ', point_t%nend
        write(*,"(1a,i9)") 'Iteration skip:', point_t%nskip
        write(*,*) '-------------------------------'
      endif
      point_t%started=.true.
    endif
    call inst_write(1)
  endif
endif
  
!  Determine if instantaneous domain velocities are to be recorded
if(domain_t%calc) then
  if(jt >= domain_t%nstart .and. jt <= domain_t%nend .and. ( jt == domain_t%nstart .or. mod(jt,domain_t%nskip)==0) ) then
    if(.not. domain_t%started) then
      if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then        
        write(*,*) '-------------------------------'
        write(*,"(1a,i9,1a,i9)") 'Writing instantaneous domain velocities from ', domain_t%nstart, ' to ', domain_t%nend
        write(*,"(1a,i9)") 'Iteration skip:', domain_t%nskip
        write(*,*) '-------------------------------'
      endif
      domain_t%started=.true.
    endif
    call inst_write(2)
  endif
endif 

!  Determine if instantaneous y-plane velocities are to be recorded
if(yplane_t%calc) then
  if(jt >= yplane_t%nstart .and. jt <= yplane_t%nend .and. ( jt == yplane_t%nstart .or. mod(jt,yplane_t%nskip)==0) ) then
    if(.not. yplane_t%started) then
      if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then        
        write(*,*) '-------------------------------'
        write(*,"(1a,i9,1a,i9)") 'Writing instantaneous y-plane velocities from ', yplane_t%nstart, ' to ', yplane_t%nend
        write(*,"(1a,i9)") 'Iteration skip:', yplane_t%nskip
        write(*,*) '-------------------------------'
      endif
      yplane_t%started=.true.
    endif
    call inst_write(3)
  endif
endif    

!  Determine if instantaneous z-plane velocities are to be recorded
if(zplane_t%calc) then
  if(jt >= zplane_t%nstart .and. jt <= zplane_t%nend .and. ( jt == zplane_t%nstart .or. mod(jt,zplane_t%nskip)==0) ) then
    if(.not. zplane_t%started) then
      if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then        
        write(*,*) '-------------------------------'
        write(*,"(1a,i9,1a,i9)") 'Writing instantaneous z-plane velocities from ', zplane_t%nstart, ' to ', zplane_t%nend
        write(*,"(1a,i9)") 'Iteration skip:', zplane_t%nskip
        write(*,*) '-------------------------------'
      endif
      zplane_t%started=.true.
    endif
    call inst_write(4)
  endif
endif

!  Determine if time-averaged z-plane data are to be recorded
if(zplane_avg_t%calc) then
  if(jt >= zplane_avg_t%nstart .and. jt <= zplane_avg_t%nend) then
    if(.not. zplane_avg_t%started) then
      if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then        
        write(*,*) '-------------------------------'
        write(*,"(1a,i9,1a,i9)") 'Writing time-averaged z-plane data from ', zplane_avg_t%nstart, ' to ', zplane_avg_t%nend
        write(*,*) '-------------------------------'
      endif
      zplane_avg_t%started=.true.	
	  !Write Tecplot header files
	  call zplane_avg_compute(1)
	endif

	!Calculate the averaged data for each time step
	  call zplane_avg_compute(2)

    !if(jt > zplane_avg_t%nstart .and. mod(jt,zplane_avg_t%nskip)==0) then
	if(jt == zplane_avg_t%nend) then
	  !Write the averaged data to file and reset the variables	
	  call zplane_avg_compute(3)	  
	endif
  endif
endif

!if(yplane_t%calc .or. zplane_t%calc) call plane_avg_compute(jt)

return
end subroutine output_loop

!**********************************************************************
subroutine tstats_init()
!**********************************************************************

!  Load tstats.out files

use param, only : coord, dt
use messages
use stat_defs, only : tstats_t
implicit none

character (*), parameter :: sub_name = mod_name // '.tstats_init'
character (*), parameter :: ftstats_in = 'tstats.out'
$if ($MPI)
character (*), parameter :: MPI_suffix = '.c'
character (128) :: fname
$endif

logical :: opn, exst

!---------------------------------------------------------------------

inquire (unit=1, opened=opn)
if (opn) call error (sub_name, 'unit 1 already open')

$if ($MPI)
write (fname, '(a,a,i0)') ftstats_in, MPI_suffix, coord
$else
fname = trim(adjustl(ftstats_in))
$endif

inquire (file=fname, exist=exst)
if (.not. exst) then
  !  Nothing to read in
  call mesg(sub_name, 'No previous time statistical data - starting from scratch.')

  tstats_t % total_time = dble(tstats_t % nend - tstats_t % nstart + 1) * dt

  tstats_t % started = .true.
  return
endif
  
open (1, file=fname, action='read', position='rewind', form='unformatted')

read (1) tstats_t % total_time
read (1) tstats_t % u 
read (1) tstats_t % v
read (1) tstats_t % w
read (1) tstats_t % u2
read (1) tstats_t % v2
read (1) tstats_t % w2
read (1) tstats_t % uw
read (1) tstats_t % vw
read (1) tstats_t % uv
read (1) tstats_t % dudz

$if($LVLSET)
$if(RNS_LS)
read (1) tstats_t % fx
read (1) tstats_t % fy
read (1) tstats_t % fz
$endif
$endif

close(1)

! Now initialize all quantities for summation
tstats_t % u = tstats_t % u * tstats_t % total_time
tstats_t % v = tstats_t % v * tstats_t % total_time
tstats_t % w = tstats_t % w * tstats_t % total_time
tstats_t % u2 = tstats_t % u2 * tstats_t % total_time
tstats_t % v2 = tstats_t % v2 * tstats_t % total_time
tstats_t % w2 = tstats_t % w2 * tstats_t % total_time
tstats_t % uw = tstats_t % uw * tstats_t % total_time
tstats_t % vw = tstats_t % vw * tstats_t % total_time
tstats_t % uv = tstats_t % uv * tstats_t % total_time
tstats_t % dudz = tstats_t % dudz * tstats_t % total_time

$if($LVLSET)
$if(RNS_LS)
tstats_t % fx = tstats_t % fx * tstats_t % total_time
tstats_t % fy = tstats_t % fy * tstats_t % total_time
tstats_t % fz = tstats_t % fz * tstats_t % total_time
$endif
$endif

tstats_t % total_time = tstats_t%total_time + dble(tstats_t % nend - tstats_t % nstart + 1) * dt


tstats_t%started=.true.

return
end subroutine tstats_init

!**********************************************************************
subroutine inst_write(itype)
!**********************************************************************
!  This subroutine writes the instantaneous values
!  at specified i,j,k locations
use functions, only : linear_interp, trilinear_interp
use stat_defs, only : point_t, domain_t, yplane_t, zplane_t
use grid_defs, only : x,y,z,zw
use sim_param, only : u,v,w

$if($LVLSET)
use level_set, only : phi
$endif

use param, only : jt_total, dt_dim, nx, ny, nz,dx,dy,dz,z_i,L_x,L_y,L_z,coord
implicit none

integer, intent(IN) :: itype

character(25) :: cl, ct
character (64) :: fname, temp, var_list
integer :: n, fid, i, j, k, nvars

real(rprec) :: ui, vi, wi

! real(rprec) :: dnx, dny, dnz

!  Write point data; assumes files have been opened properly
!  in stats_init

!  Make sure w has been interpolated to uv-grid
call interp_to_uv_grid(w, w_uv, w_uv_tag)

if(itype==1) then

  do n=1,point_t%nloc
!  Files have been opened in stats_init

!  For parallel runs check if data is on correct proc
    $if ($MPI)
    if(point_t%coord(n) == coord) then
    $endif

    call write_real_data(point_t%fname(n), 'append', 4, (/ total_time_dim, &
      trilinear_interp(u,point_t%istart(n),point_t%jstart(n), point_t%kstart(n), point_t%xyz(:,n)), &
      trilinear_interp(v,point_t%istart(n),point_t%jstart(n), point_t%kstart(n), point_t%xyz(:,n)), &
      trilinear_interp(w,point_t%istart(n),point_t%jstart(n), point_t%kstart(n), point_t%xyz(:,n)) /))
	  


    $if ($MPI)
    endif
    $endif

  enddo
!  Instantaneous write for entire domain
elseif(itype==2) then

!  Convert total iteration time to string
  write(ct,*) jt_total
!  Open file which to write global data
  write (fname,*) 'output/uvw.', trim(adjustl(ct)),'.dat'
  fname = trim(adjustl(fname))

  $if ($MPI)
    write (temp, '(".c",i0)') coord
    fname = trim (fname) // temp
  $endif
  


  $if($LVLSET)
  !write(7,*) 'variables = "x", "y", "z", "u", "v", "w", "phi"';
  var_list = '"x", "y", "z", "u", "v", "w", "phi"'
  nvars = 7
  $else
  !write(7,*) 'variables = "x", "y", "z", "u", "v", "w"';
  var_list = '"x", "y", "z", "u", "v", "w"'
  nvars = 6
  $endif
  
  
  call write_tecplot_header_ND(fname, 'rewind', nvars, (/ Nx, Ny, Nz/), var_list, coord, 2, total_time_dim)
  !write_tecplot_header_ND(fname, write_posn, var_list, nvars, zone, data_prec, domain_size, soln_time)
  

  !write(7,"(1a,i9,1a,i3,1a,i3,1a,i3,1a,i3)") 'ZONE T="', &
  !  j,'", DATAPACKING=POINT, i=', Nx,', j=',Ny,', k=', Nz

  !!$if($LVLSET)
  !!write(7,"(1a)") ''//adjustl('DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)')//''
  !!$else
  !!write(7,"(1a)") ''//adjustl('DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)')//''
  !!$endif

  !!write(7,"(1a,f18.6)") 'solutiontime=', total_time_dim
  !open(unit = 7,file = fname, status='old',form='formatted', &
  !  action='write',position='append')
	
  $if($LVLSET)
  call write_real_data_3D(fname, 'append', 'formatted', 4, nx, ny,nz, &
    (/ u(1:nx,1:ny,1:nz), v(1:nx,1:ny,1:nz), w_uv(1:nx,1:ny,1:nz), phi(1:nx,1:ny,1:nz) /), x, y, z(1:nz))
  $else
  call write_real_data_3D(fname, 'append', 'formatted', 3, nx,ny,nz, &
    (/ u(1:nx,1:ny,1:nz), v(1:nx,1:ny,1:nz), w_uv(1:nx,1:ny,1:nz) /), x, y, z(1:nz))
  $endif
	
  !do k=1,nz
  !  do j=1,ny
  !    do i=1,nx
  !      
  !      write(7,*) x(i), y(j), z(k), u(i,j,k), v(i,j,k), w_uv(i,j,k), phi(i,j,k)
  !      $else
  !      write(7,*) x(i), y(j), z(k), u(i,j,k), v(i,j,k), w_uv(i,j,k)
  !      $endif

  !    enddo
  !  enddo
  !enddo
  
  !close(7)
!  Write instantaneous y-plane values
elseif(itype==3) then

!  Loop over all yplane locations
  do j=1,yplane_t%nloc

    write(cl,'(F9.4)') yplane_t%loc(j)
    !  Convert total iteration time to string
    write(ct,*) jt_total
    write(fname,*) 'output/uvw.y-',trim(adjustl(cl)),'.',trim(adjustl(ct)),'.dat'
    fname=trim(adjustl(fname))

    $if ($MPI)
!  For MPI implementation
      write (temp, '(".c",i0)') coord
      fname = trim (fname) // temp
    $endif


    !write(2,*) 'variables = "x", "y", "z", "u", "v", "w"';
    !write(2,"(1a,i9,1a,i3,1a,i3,1a,i3,1a,i3)") 'ZONE T="', &
    !  j,'", DATAPACKING=POINT, i=', Nx,', j=',1,', k=', Nz
    !write(2,"(1a)") ''//adjustl('DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)')//''
    !write(2,"(1a,f18.6)") 'solutiontime=', total_time_dim
	
    call write_tecplot_header_ND(fname, 'rewind', 6, (/ Nx, 1, Nz/), &
	  '"x", "y", "z", "u", "v", "w"', coord, 2, total_time_dim)  
	  
	open (unit = 2,file = fname, status='unknown',form='formatted', &
      action='write',position='append')
	  
	!call write_real_data_2D(fname, 'rewind', 'formatted', 3, nx, nz, &
	!  (/ linear_interp(u(:,yplane_t%istart(j),:), &
    !      u(:,yplane_t%istart(j)+1,:), dy, yplane_t%ldiff(j)), &
	!	  linear_interp(v(:,yplane_t%istart(j),:), &
    !      v(:,yplane_t%istart(j)+1,:), dy, yplane_t%ldiff(j)), &
	!	  linear_interp(w_uv(:,yplane_t%istart(j),:), &
    !      w_uv(:,yplane_t%istart(j)+1,:), dy, &
    !      yplane_t%ldiff(j)) /), x, z
	  
    do k=1,nz
      do i=1,nx

        ui = linear_interp(u(i,yplane_t%istart(j),k), &
          u(i,yplane_t%istart(j)+1,k), dy, yplane_t%ldiff(j))
        vi = linear_interp(v(i,yplane_t%istart(j),k), &
          v(i,yplane_t%istart(j)+1,k), dy, yplane_t%ldiff(j))
        wi = linear_interp(w_uv(i,yplane_t%istart(j),k), &
          w_uv(i,yplane_t%istart(j)+1,k), dy, &
          yplane_t%ldiff(j))

        write(2,*) x(i), yplane_t%loc(j), z(k), ui, vi, wi
      enddo
    enddo
    close(2)
  enddo  

!  Write instantaneous z-plane values
elseif(itype==4) then

!  Loop over all zplane locations
  do k=1,zplane_t%nloc
  
    $if ($MPI)
    if(zplane_t%coord(k) == coord) then
    $endif

    write(cl,'(F9.4)') zplane_t%loc(k)
    !  Convert total iteration time to string
    write(ct,*) jt_total
    write(fname,*) 'output/uvw.z-',trim(adjustl(cl)),'.',trim(adjustl(ct)),'.dat'
    fname=trim(adjustl(fname))

!     $if ($MPI)
! !  For MPI implementation
!       write (temp, '(".c",i0)') coord
!       fname = trim (fname) // temp
!     $endif


    !write(2,*) 'variables = "x", "y", "z", "u", "v", "w"';
    !write(2,"(1a,i9,1a,i3,1a,i3,1a,i3,1a,i3)") 'ZONE T="', &
    !  j,'", DATAPACKING=POINT, i=', Nx,', j=',Ny,', k=', 1
    !write(2,"(1a)") ''//adjustl('DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)')//''
    !write(2,"(1a,f18.6)") 'solutiontime=', total_time_dim
	
    call write_tecplot_header_ND(fname, 'rewind', 6, (/ Nx, Ny, 1/), &
	 '"x", "y", "z", "u", "v", "w"', coord, 2, total_time_dim)
	  	  
    open (unit = 2,file = fname, status='unknown',form='formatted', &
      action='write',position='append')	  

    do j=1,Ny
      do i=1,Nx

        ui = linear_interp(u(i,j,zplane_t%istart(k)),u(i,j,zplane_t%istart(k)+1), &
          dz, zplane_t%ldiff(k))
        vi = linear_interp(v(i,j,zplane_t%istart(k)),v(i,j,zplane_t%istart(k)+1), &
          dz, zplane_t%ldiff(k))
        wi = linear_interp(w_uv(i,j,zplane_t%istart(k)), &
          w_uv(i,j,zplane_t%istart(k)+1), &
          dz, zplane_t%ldiff(k))

        write(2,*) x(i), y(j), zplane_t%loc(k), ui, vi, wi

      enddo
    enddo
    close(2)

    $if ($MPI)
    endif
    $endif

  enddo  


else
  write(*,*) 'Error: itype not specified properly to inst_write!'
  stop
endif
return
end subroutine inst_write

!*************************************************************
subroutine write_real_data(fname, write_pos, nvars, vars)
!*************************************************************
! 
!  This subroutine writes the variables given by vars to the
!  specified file, fname. The number of variables can be arbitrary
!  but must be specified by nvars. The output will be of the form
!  
!  vars(1) vars(2) ... vars(nvars-1) vars(nvars)
!
!  The primary purpose of this routine is to write instantaneous data
!  or something similar to file
!
!  Inputs:
!  fname (char) - file to write to
!  write_pos (char) - position to write in file : 'append' or 'rewind'
!  nvars (int) - number of variables contained in vars
!  vars (real,vector) - vector contaning variables to write
!
implicit none

character(*), intent(in) :: fname, write_pos
integer, intent(in) :: nvars
real(rprec), intent(in), dimension(:) :: vars


character(64) :: frmt
logical :: exst

character(*), parameter :: sub_name = mod_name // '.write_real_data'

!  Check if file exists
!inquire ( file=fname, exist=exst)
!if (.not. exst) call mesg(sub_name, 'Creating : ' // fname)

call check_write_pos(write_pos, sub_name)

open (unit = 2,file = fname, status='unknown',form='formatted', &
  action='write',position='append')
  
!  Specify output format; may want to use a global setting
write (frmt, '("(",i0,"e)")') nvars
  
!  Write the data
write(2,frmt) vars

close(2)
  
return
end subroutine write_real_data

!*************************************************************
subroutine write_real_data_3D(fname, write_pos, write_fmt, nvars, &
  imax, jmax, kmax, vars, x,y,z)
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
!  write_pos (char) - postition to write in file : 'append' or 'rewind'
!  write_fmt (char) - Fotran format flag : 'formatted' or 'unformatted'
!  nvars (int) - number of variables contained in vars
!  imax (int) - size of 1st dimension of variables in vars
!  jmax (int) - size of 2nd dimension of variables in vars  
!  kmax (int)- size of 3rd dimension of variables in vars
!  vars (real, vector) - vector contaning variables to write
!  x,y,z (real, vector, optional) - vectors containing x,y,z coordinates 
!
implicit none

character(*), intent(in) :: fname, write_pos, write_fmt
integer, intent(in) :: nvars, imax, jmax, kmax
real(rprec), intent(in), dimension(nvars*imax*jmax*kmax) :: vars
real(rprec), intent(in), dimension(:), optional :: x,y,z

character(64) :: frmt
logical :: exst, xpres

character(*), parameter :: sub_name = mod_name // '.write_real_data_3D'

integer :: i,j,k,n
real(rprec), allocatable, dimension(:,:,:,:) :: vars_dim

!  Check if file exists
!inquire ( file=fname, exist=exst)
!if (.not. exst) call mesg(sub_name, 'Creating : ' // fname)
call check_write_pos(write_pos, sub_name)
call check_write_fmt(write_fmt, sub_name)

!  Check if spatial coordinates specified
xpres=.false.
if(present(x) .and. present(y) .and. present(z)) xpres = .true.

!  Put data in correct shape for sequential reads/writes
allocate(vars_dim(nvars,imax,jmax,kmax)) 
do n=1,nvars
  do k=1,kmax
    do j=1,jmax
	  do i=1,imax
	    vars_dim(n,i,j,k) = vars((n-1)*imax*jmax*kmax + (k-1)*imax*jmax + (j-1)*imax + i)
	  enddo
	enddo
  enddo
enddo

open (unit = 2,file = fname, status='unknown',form=write_fmt, &
  action='write',position=write_pos)
  
!  Write the data
select case(write_fmt)
  case('formatted')
  
    !  Specify output format; may want to use a global setting
    	
    if (xpres) then
	  
	  write (frmt, '("(3e,",i0,"e)")') nvars
	  
	  do k=1, kmax
	    do j=1, jmax
  	      do i=1, imax
            write(2,frmt) x(i), y(j), z(k), vars_dim(:,i,j,k)
	      enddo 
	    enddo
	  enddo
	  
	else
	  
	  write (frmt, '("(",i0,"e)")') nvars
	 
	  do k=1, kmax
	    do j=1, jmax
  	      do i=1, imax
            write(2,frmt) vars_dim(:,i,j,k)
	      enddo 
	    enddo
	  enddo
	  
	endif

  case('unformatted')
  
    if (xpres) then
	  
	  do k=1, kmax
	    do j=1, jmax
  	      do i=1, imax
            write(2) x(i), y(j), z(k), vars_dim(:,i,j,k)
	      enddo 
	    enddo
	  enddo
	  
	else
	
	  do k=1, kmax
	    do j=1, jmax
  	      do i=1, imax
            write(2) vars_dim(:,i,j,k)
	      enddo 
	    enddo
	  enddo
	  
	endif
	
  !case default
  
  !  call error(sub_name, 'Incorrect write format : ' // write_pos)
	
end select



close(2)
  
deallocate(vars_dim)

return
end subroutine write_real_data_3D

!*************************************************************
subroutine write_real_data_2D(fname, write_pos, write_fmt, nvars, &
  imax, jmax, vars, x, y)
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
!  write_pos (char) - postition to write in file : 'append' or 'rewind'
!  write_fmt (char) - Fotran format flag : 'formatted' or 'unformatted'
!  nvars (int) - number of variables contained in vars
!  imax (int) - size of 1st dimension of variables in vars
!  jmax (int) - size of 2nd dimension of variables in vars 
!  vars (real, vector) - vector contaning variables to write
!  x,y (real, vector, optional) - vectors containing x,y coordinates 
!

implicit none

character(*), intent(in) :: fname, write_pos, write_fmt
integer, intent(in) :: nvars, imax, jmax
real(rprec), intent(in), dimension(nvars*imax*jmax) :: vars
real(rprec), intent(in), dimension(:), optional :: x, y

character(64) :: frmt
logical :: exst, xpres

character(*), parameter :: sub_name = mod_name // '.write_real_data_2D'

integer :: i,j,n
  
real(rprec), allocatable, dimension(:,:,:) :: vars_dim

!  Check if file exists
!inquire ( file=fname, exist=exst)
!if (.not. exst) call mesg(sub_name, 'Creating : ' // fname)
call check_write_pos(write_pos, sub_name)
call check_write_fmt(write_fmt, sub_name)

!  Check if spatial coordinates specified
xpres=.false.
if(present(x) .and. present(y)) xpres = .true.

open (unit = 2,file = fname, status='unknown',form=write_fmt, &
  action='write',position=write_pos)

!  Put data in correct shape for sequential reads/writes
allocate(vars_dim(nvars,imax,jmax)) 
do n=1,nvars
  do j=1,jmax
    do i=1,imax
	  vars_dim(n,i,j) = vars((n-1)*jmax*imax + (j-1)*imax + i) 
	enddo
  enddo
enddo
  
  
!  Write the data
select case(write_fmt)
  case('formatted')
  
    !  Specify output format; may want to use a global setting
    	
    if (xpres) then
	  
	  write (frmt, '("(2e,",i0,"e)")') nvars
	  
	  do j=1, jmax
	    do i=1,imax
          write(2,frmt) x(i), y(j), vars_dim(:,i,j)
	    enddo 
	  enddo
	  
	else
	  
	  write (frmt, '("(",i0,"e)")') nvars
	 
	  do j=1, jmax
	    do i=1,imax
          write(2,frmt) vars_dim(:,i,j)
	    enddo 
	  enddo
	  
	endif

  case('unformatted')
  
    if (xpres) then
	  
	  do j=1, jmax
	    do i=1,imax
          write(2) x(i), y(j), vars_dim(:,i,j)
	    enddo 
	  enddo
	  
	else
	
	  do j=1,jmax
	    do i=1,imax
          write(2) vars_dim(:,i,j)
	    enddo 
      enddo
	  
	endif
	
  !case default
  !  call error(sub_name, 'Incorrect write format : ' // write_pos)
end select

close(2)

deallocate(vars_dim)

return
end subroutine write_real_data_2D

!*************************************************************
subroutine write_real_data_1D(fname, write_pos, write_fmt, nvars, &
  imax, vars, x)
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
!  write_pos (char) - postition to write in file : 'append' or 'rewind'
!  write_fmt (char) - Fotran format flag : 'formatted' or 'unformatted'
!  nvars (int) - number of variables contained in vars
!  imax (int) - size of 1st dimension of variables in vars
!  vars (real, vector) - vector contaning variables to write
!  x (real,vector,optional) - vectors containing x coordinates 
!
implicit none

character(*), intent(in) :: fname, write_pos, write_fmt
integer, intent(in) :: nvars, imax
real(rprec), intent(in), dimension(nvars*imax) :: vars
real(rprec), intent(in), dimension(:), optional :: x

character(64) :: frmt
logical :: exst, xpres

character(*), parameter :: sub_name = mod_name // '.write_real_data_1D'

integer :: i,n

real(rprec), allocatable, dimension(:,:) :: vars_dim

!  Check if file exists
!inquire ( file=fname, exist=exst)
!if (.not. exst) call mesg(sub_name, 'Creating : ' // fname)
call check_write_pos(write_pos, sub_name)
call check_write_fmt(write_fmt, sub_name)

!  Check if spatial coordinates specified
if(present(x)) xpres = .true.

!  Put data in correct shape for sequential reads/writes
allocate(vars_dim(nvars,imax)) 
do n=1,nvars
  do i=1,imax
     vars_dim(n,i) = vars((n-1)*imax + i)
  enddo
enddo

open (unit = 2,file = fname, status='unknown',form=write_fmt, &
  action='write',position=write_pos)
   
!  Write the data
select case(write_fmt)
  case('formatted')
  
    if (xpres) then
	
      write (frmt, '("(1e,",i0,"e)")') nvars
	  
	  do i=1,imax
        write(2,frmt) x(i), vars_dim(:,i)
	  enddo 
	  
	else
	
	  write (frmt, '("(",i0,"e)")') nvars
	  
	  do i=1,imax
        write(2,frmt) vars_dim(:,i)
	  enddo 
	  
	endif

  case('unformatted')
  
    if (xpres) then
	  
	  do i=1,imax
        write(2) x(i), vars_dim(:,i)
	  enddo 
	  
	else
	
	  do i=1,imax
        write(2) x(i), vars_dim(:,i)
	  enddo 
	  
	endif
	
  !case default
  !  call error(sub_name, 'Incorrect write format : ' // write_pos)
end select

close(2)

deallocate(vars_dim)

return
end subroutine write_real_data_1D

!*************************************************************
subroutine write_tecplot_header_xyline(fname, write_pos, var_list)
!*************************************************************
!  The purpose of this routine is to write Tecplot header information
!  for xy-line formatted files
! 
!  Inputs:
!  fname (char)     - file name to write to
!  write_pos (char) - position in file to write data: append or rewind
!  var_list	(char)  - string containing variable names: Ex. "x", "u"
!

implicit none

character(*), intent(in) :: fname, write_pos, var_list
character(120) :: tec_dt_str

character(*), parameter :: sub_name = mod_name // '.write_tecplot_header_xyline'

!  Check if write position has been specified correctly
call check_write_pos(write_pos, sub_name)	

open (unit = 2,file = fname, status='unknown',form='formatted', &
  action='write',position=write_pos)

write(2,'(1a)') 'variables = ' // var_list

close(2)

return
end subroutine write_tecplot_header_xyline


!*************************************************************
subroutine write_tecplot_header_ND(fname, write_pos, nvars, &
  domain_size, var_list, zone, data_type, soln_time)
!*************************************************************
!  The purpose of this routine is to write Tecplot header information
!  for 1D, 2D, and 3D data files.
!
!  NOTE: domain_size needs to be specified as a vector even for 1D data.
!  Ex: 
!    1D : (/ Nx /)
!    2D : (/ Nx, Ny /)
!    3D : (/ Nx, Ny, Nz /)
! 
!  Inputs:
!  fname (char)     - file name to write to
!  write_pos (char) - position in file to write data: append or rewind
!  nvars (int)      - number of variables
!  domain_size (int, vector) - vector containing the diminsions of the data.
!  var_list	(char)  - string containing variable names: Ex. '"x", "u"'
!  zone (int)       - zone number
!  date_type (int) 	- specify Tecplot data type (precision): 1 - single, 2 - double
!  soln_time (real, optional) - time stamp
!
implicit none

character(*), intent(in) :: fname, write_pos
integer, intent(in) :: nvars
integer, dimension(:), intent(in) :: domain_size
character(*), intent(in) :: var_list
integer, intent(in) :: zone, data_type
real(rprec), optional, intent(in) :: soln_time

character(*), parameter :: sub_name = mod_name // '.write_tecplot_header_ND'
character(64) :: tec_dt, posn
character(120) :: tec_dt_str, tec_dat_str
integer :: ndims, n

!  Check if write position has been specified correctly
call check_write_pos(write_pos, sub_name)

!  Get number of dimensions for data
ndims = size(domain_size,1)

if(ndims == 1) then
  write(tec_dat_str,"(1a,i9,1a,i3,1a,i3)") 'ZONE T="', &
    zone,'", DATAPACKING=POINT, i=', domain_size(1)
elseif(ndims == 2) then
  write(tec_dat_str,"(1a,i9,1a,i3,1a,i3,1a,i3)") 'ZONE T="', &
    zone,'", DATAPACKING=POINT, i=', domain_size(1),', j=', domain_size(2)
elseif(ndims == 3) then
  write(tec_dat_str,"(1a,i9,1a,i3,1a,i3,1a,i3,1a,i3)") 'ZONE T="', &
    zone,'", DATAPACKING=POINT, i=', domain_size(1),', j=', domain_size(2),', k=', domain_size(3)
else
  call error(sub_name, 'Incorrect number of dimensions : ', ndims)
endif

!  Create Tecplot DT string
call tecplot_data_type_str(nvars, data_type, sub_name, tec_dt_str)

open (unit = 2,file = fname, status='unknown',form='formatted', &
  action='write',position=write_pos)

!  Write variable list
write(2,'(1a)') 'variables = ' // var_list
!  Write data layout size information
write(2,'(1a)') tec_dat_str
!  Write Tecplot data type for each variable
write(2,'(1a)') tec_dt_str

if (present (soln_time)) then
write(2,'(1a,e)') 'solutiontime=', soln_time
endif

close(2)

  
return
end subroutine write_tecplot_header_ND

!*************************************************************
subroutine tecplot_data_type_str(nvars, data_type, sub_name, tec_dt_str)
!*************************************************************
implicit none

integer, intent(in) ::  nvars, data_type
character, intent(IN) :: sub_name
character(120), intent(OUT) :: tec_dt_str
character(7) :: tec_dt

!character(*), parameter :: sub_name = mod_name // '.tecplot_data_type_str'

integer :: n

!  Specify single or double precision
if(data_type == 1) then
  tec_dt = ' SINGLE'
elseif(data_type == 2) then
  tec_dt = ' DOUBLE'
else
  call error(sub_name, 'Incorrect data type : ', data_type)
endif

!  Create DT string
tec_dt_str = 'DT=('
do n=1, nvars
  call strcat(tec_dt_str,tec_dt)
enddo
call strcat(tec_dt_str,')')

return
end subroutine tecplot_data_type_str

!*************************************************************
subroutine check_write_pos(write_pos, sub_name)
!*************************************************************
!  This subroutine checks whether the write position has been
!  specified properly
!
implicit none

character(*), intent(IN) :: write_pos, sub_name

!  Check if write position has been specified correctly
select case(write_pos)
  case('append')
  !  do nothing
  case('rewind')
  !  do nothing
  case default
    call error(sub_name, 'Incorrect write position : ' // write_pos)
end select

return
end subroutine check_write_pos

!*************************************************************
subroutine check_write_fmt(write_fmt, sub_name)
!*************************************************************
!  This subroutine checks whether the write position has been
!  specified properly
!
implicit none

character(*), intent(IN) :: write_fmt, sub_name

!  Check if write position has been specified correctly
select case(write_fmt)
  case('formatted')
  !  do nothing
  case('unformatted')
  !  do nothing
  case default
    call error(sub_name, 'Incorrect write position : ' // write_fmt)
end select

return
end subroutine check_write_fmt


! !**********************************************************************
! subroutine rs_write()
! !**********************************************************************
! use grid_defs, only : x,y,z
! use stat_defs, only : rs_t
! use param, only : nx,ny,nz,dx,dy,dz,L_x,L_y,L_z
! implicit none
! 
! character (64) :: fname, temp
! integer i,j,k
! 
! !  Set default output file name
! fname = 'output/rs.dat'
! 
! $if ($MPI)
! !  Update fname for MPI implementation     
!   write (temp, '(".c",i0)') coord
!   fname = trim (fname) // temp
! $endif
! 
! !  Add temporary output information
! open(unit = 7,file = fname)
! write(7,*) 'variables= "x", "y", "z", "up2", "vp2", "wp2", "upwp", "vpwp", "upvp"'
! write(7,"(1a,i9,1a,i3,1a,i3,1a,i3,1a,i3)") 'ZONE T="', &
!   1,'", DATAPACKING=POINT, i=', Nx,', j=',Ny, ', k=', Nz
! write(7,"(1a)") ''//adjustl('DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)')//''  
! do k=1,nz
!   do j=1,ny
!     do i=1,nx
! !  Write spatially averaged, temporally averaged quantities   
!      write(7,*) x(i), y(j), z(k), rs_t%up2(i,j,k), rs_t%vp2(i,j,k), rs_t%wp2(i,j,k), &
!        rs_t%upwp(i,j,k), rs_t%vpwp(i,j,k), rs_t%upvp(i,j,k)
!     enddo
!   enddo
! enddo
! close(7)
!   
! return
! end subroutine rs_write

!**********************************************************************
subroutine tstats_finalize()
!**********************************************************************
use grid_defs, only : x,y,z
use stat_defs, only : tstats_t
use param, only : nx,ny,nz,dx,dy,dz,L_x,L_y,L_z
implicit none

character (*), parameter :: sub_name = mod_name // '.tstats_finalize'
character (64) :: fname, temp, fname_dat
integer :: i,j,k

logical :: opn

fname = 'tstats.out'
fname_dat = 'output/tstats.dat'

$if ($MPI)
!  For MPI implementation     
  write (temp, '(".c",i0)') coord
  fname = trim (fname) // temp
  fname_dat = trim (fname_dat) // temp
$endif

!  Perform Averaging operation
tstats_t % u = tstats_t % u / tstats_t % total_time
tstats_t % v = tstats_t % v / tstats_t % total_time
tstats_t % w = tstats_t % w / tstats_t % total_time
tstats_t % u2 = tstats_t % u2 / tstats_t % total_time
tstats_t % v2 = tstats_t % v2 / tstats_t % total_time
tstats_t % w2 = tstats_t % w2 / tstats_t % total_time
tstats_t % uw = tstats_t % uw / tstats_t % total_time
tstats_t % vw = tstats_t % vw / tstats_t % total_time
tstats_t % uv = tstats_t % uv / tstats_t % total_time
tstats_t % dudz = tstats_t % dudz / tstats_t % total_time

$if($LVLSET)
$if(RNS_LS)
tstats_t % fx = tstats_t % fx / tstats_t % total_time
tstats_t % fy = tstats_t % fy / tstats_t % total_time
tstats_t % fz = tstats_t % fz / tstats_t % total_time
$endif
$endif

inquire (unit=1, opened=opn)
if (opn) call error (sub_name, 'unit 1 already open')

open (1, file=fname, action='write', position='rewind', form='unformatted')

write (1) tstats_t % total_time
write (1) tstats_t % u 
write (1) tstats_t % v
write (1) tstats_t % w
write (1) tstats_t % u2
write (1) tstats_t % v2
write (1) tstats_t % w2
write (1) tstats_t % uw
write (1) tstats_t % vw
write (1) tstats_t % uv
write (1) tstats_t % dudz

$if($LVLSET)
$if(RNS_LS)
write (1) tstats_t % fx
write (1) tstats_t % fy
write (1) tstats_t % fz
$endif
$endif

close(1)

 ! call write_real_data_3D(fname,'rewind', 'formatted', 10, nx, ny, nz, (/ tstats_t%u(1:nx,1:ny,1:nz), &
 ! tstats_t%v(1:nx,1:ny,1:nz), tstats_t%w(1:nx,1:ny,1:nz), tstats_t%u2(1:nx,1:ny,1:nz), &
 ! tstats_t%v2(1:nx,1:ny,1:nz), tstats_t%w2(1:nx,1:ny,1:nz), tstats_t%uw(1:nx,1:ny,1:nz), &
 ! tstats_t%vw(1:nx,1:ny,1:nz), tstats_t%uv(1:nx,1:ny,1:nz), tstats_t%dudz(1:nx,1:ny,1:nz) /), &
 ! x, y, z(1:nz))
 
$if($LVLSET)
$if(RNS_LS)

call write_tecplot_header_ND(fname_dat, 'rewind', 16, (/ Nx, Ny, Nz/), &
   '"x", "y", "z", "u","v","w","u2","v2","w2","uw","vw","uv","dudz", "fx", "fy", "fz"', coord, 2)

call write_real_data_3D(fname_dat,'append', 'formatted', 13, nx, ny, nz, (/ tstats_t%u(1:nx,1:ny,1:nz), &
tstats_t%v(1:nx,1:ny,1:nz), tstats_t%w(1:nx,1:ny,1:nz), tstats_t%u2(1:nx,1:ny,1:nz), &
tstats_t%v2(1:nx,1:ny,1:nz), tstats_t%w2(1:nx,1:ny,1:nz), tstats_t%uw(1:nx,1:ny,1:nz), &
tstats_t%vw(1:nx,1:ny,1:nz), tstats_t%uv(1:nx,1:ny,1:nz), tstats_t%dudz(1:nx,1:ny,1:nz), &
tstats_t%fx(1:nx,1:ny,1:nz), tstats_t%fy(1:nx,1:ny,1:nz), tstats_t%fz(1:nx,1:ny,1:nz) /), &
x(1:nx), y(1:ny), z(1:nz))

$endif

$else

call write_tecplot_header_ND(fname_dat, 'rewind', 10, (/ Nx, Ny, Nz/), &
   '"u","v","w","u2","v2","w2","uw","vw","uv","dudz"', coord, 2)

call write_real_data_3D(fname_dat,'append', 'formatted', 10, nx, ny, nz, (/ tstats_t%u(1:nx,1:ny,1:nz), &
tstats_t%v(1:nx,1:ny,1:nz), tstats_t%w(1:nx,1:ny,1:nz), tstats_t%u2(1:nx,1:ny,1:nz), &
tstats_t%v2(1:nx,1:ny,1:nz), tstats_t%w2(1:nx,1:ny,1:nz), tstats_t%uw(1:nx,1:ny,1:nz), &
tstats_t%vw(1:nx,1:ny,1:nz), tstats_t%uv(1:nx,1:ny,1:nz), tstats_t%dudz(1:nx,1:ny,1:nz) /), &
x(1:nx), y(1:ny), z(1:nz))


$endif
 
 !call write_real_data_3D(fname,'rewind', 'unformatted', 1, nx, ny, nz, tstats_t%u(1:nx,1:ny,1:nz))
 !call write_real_data_3D(fname,'append', 'unformatted', 1, nx, ny, nz, tstats_t%v(1:nx,1:ny,1:nz))
 !call write_real_data_3D(fname,'append', 'unformatted', 1, nx, ny, nz, tstats_t%w(1:nx,1:ny,1:nz))
 !call write_real_data_3D(fname,'append', 'unformatted', 1, nx, ny, nz, tstats_t%u2(1:nx,1:ny,1:nz))
 !call write_real_data_3D(fname,'append', 'unformatted', 1, nx, ny, nz, tstats_t%v2(1:nx,1:ny,1:nz))
 !call write_real_data_3D(fname,'append', 'unformatted', 1, nx, ny, nz, tstats_t%w2(1:nx,1:ny,1:nz))
 !call write_real_data_3D(fname,'append', 'unformatted', 1, nx, ny, nz, tstats_t%uw(1:nx,1:ny,1:nz))
 !call write_real_data_3D(fname,'append', 'unformatted', 1, nx, ny, nz, tstats_t%vw(1:nx,1:ny,1:nz))
 !call write_real_data_3D(fname,'append', 'unformatted', 1, nx, ny, nz, tstats_t%uv(1:nx,1:ny,1:nz))
 !call write_real_data_3D(fname,'append', 'unformatted', 1, nx, ny, nz, tstats_t%dudz(1:nx,1:ny,1:nz))

return
end subroutine tstats_finalize

! !**********************************************************************
! subroutine tavg_write()
! !**********************************************************************
! use grid_defs, only : x,y,z
! use stat_defs, only : tavg_t
! use param, only : nx,ny,nz,dx,dy,dz,L_x,L_y,L_z
! implicit none
! 
! character (64) :: fname, temp
! integer :: i,j,k
! 
! fname = 'output/uvw_avg.dat'
! 
! $if ($MPI)
! !  For MPI implementation     
!   write (temp, '(".c",i0)') coord
!   fname = trim (fname) // temp
! $endif
! 
! open(unit = 7,file = fname)
! !  open(unit = 8,file = 'output/avg_dudz.out')
! write(7,*) 'variables= "x", "y", "z", "<u>", "<v>", "<w>"'
! write(7,"(1a,i9,1a,i3,1a,i3,1a,i3,1a,i3)") 'ZONE T="', &
!   1,'", DATAPACKING=POINT, i=', Nx,', j=',Ny, ', k=', Nz
! write(7,"(1a)") ''//adjustl('DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)')//''	
! !  write(8,*) 'variables= "z", "<dudz>/u*"'
! do k=1,nz
!   do j=1,ny
!     do i=1,nx
!        write(7,*) x(i), y(j), z(k), tavg_t%u(i,j,k), tavg_t%v(i,j,k), tavg_t%w(i,j,k)
! !	write(8,*) z(k), sum(tavg_t%dudz(:,:,k))/(dnx*dny)
!     enddo
!   enddo
! enddo
! close(7)
! !  close(8)
! 
! return
! end subroutine tavg_write


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--assumes lun is open and positioned correctly
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine checkpoint (lun)
use param, only : nz, S_FLAG
use sim_param, only : u, v, w, RHSx, RHSy, RHSz, theta
use sgsmodule, only : Cs_opt2, F_LM, F_MM, F_QN, F_NN
use scalars_module, only : RHS_T
implicit none

integer, intent (in) :: lun

!---------------------------------------------------------------------

if (S_FLAG) then !WITH SCALARS
  write (lun) u(:, :, 1:nz), v(:, :, 1:nz), w(:, :, 1:nz),           &
              RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),  &
              Cs_opt2, F_LM, F_MM, F_QN, F_NN, theta, RHS_T
else ! No SCALARS
  write (lun) u(:, :, 1:nz), v(:, :, 1:nz), w(:, :, 1:nz),           &
              RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),  &
              Cs_opt2, F_LM, F_MM, F_QN, F_NN
end if

end subroutine checkpoint

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--lun_opt gives option for writing to different unit, and is used by
!  inflow_write
!--assumes the unit to write to (lun_default or lun_opt is already
!  open for sequential unformatted writing
!--this routine also closes the unit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output_final(jt, lun_opt)
use stat_defs, only : rs_t, tstats_t, point_t, yplane_t, zplane_t

implicit none

integer,intent(in)::jt
integer, intent (in), optional :: lun_opt  !--if present, write to unit lun
integer, parameter :: lun_default = 11
integer::i,fid,jx,jy,jz
integer :: lun

logical :: opn

!---------------------------------------------------------------------

if (present (lun_opt)) then
  lun = lun_opt
else
  lun = lun_default
end if

inquire (unit=lun, opened=opn)

if (.not. opn) then
  write (*, *) 'output_final: lun=', lun, ' is not open'
  stop
end if

rewind (lun)

call checkpoint (lun)

close (lun)


!  Update total_time.dat after simulation
if ((cumulative_time) .and. (lun == lun_default)) then
  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
    !--only do this for true final output, not intermediate recording
    open (1, file=fcumulative_time)
    write (1, *) jt_total, total_time, total_time_dim
    close (1)
  end if
end if

! !  Check if writing statistics
! if(rs_t%calc) then
!   call rs_compute()
!   call rs_write()
! endif 
! 
! !  Check if average quantities are to be recorded
! if(tavg_t%calc) call tavg_write()

!  Check if average quantities are to be recorded
if(tstats_t%calc) call tstats_finalize()
 
!  Close instantaneous velocity files
if(point_t%calc) then

  do i=1,point_t%nloc
    $if ($MPI)
    if(zplane_t%coord(i) == coord) then
    $endif

    fid=3000*i
    close(fid)

    $if ($MPI)
    endif
    $endif

  enddo

endif

! !  Write averaged yplane data
! if(yplane_t%calc .or. zplane_t%calc) call plane_write()

return
end subroutine output_final

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine plane_write()
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !  This subroutine is used to write yplane data to a formatted Fortran
! !  file in Tecplot format
! use grid_defs, only : x,y,z
! use stat_defs, only : yplane_t, zplane_t
! use param, only : Nx, Ny, Nz, dx, dy, dz
! implicit none
! character(50) :: ct
! character (64) :: fname, temp
! integer :: i,j,k
! 
! if(yplane_t%calc) then
!   do j=1,yplane_t%nloc
!     write(ct,'(F9.4)') yplane_t%loc(j)
!     write(fname,*) 'output/uvw_avg.y-',trim(adjustl(ct)),'.dat'
!     fname=trim(adjustl(fname))
! 
!     $if ($MPI)
! !  For MPI implementation
!       write (temp, '(".c",i0)') coord
!       fname = trim (fname) // temp
!     $endif
! 
!     open (unit = 2,file = fname, status='unknown',form='formatted', &
!       action='write',position='rewind')
!     write(2,*) 'variables = "x", "z", "u", "v", "w"';
!     write(2,"(1a,i9,1a,i3,1a,i3,1a,i3)") 'ZONE T="', &
!       j,'", DATAPACKING=POINT, i=', Nx,', j=',Nz
!     write(2,"(1a)") ''//adjustl('DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)')//''
!     do k=1,nz
!       do i=1,nx
!         write(2,*) x(i), z(k), yplane_t%ua(i,j,k), yplane_t%va(i,j,k), yplane_t%wa(i,j,k)
!       enddo
!     enddo
!     close(2)
!   enddo
! endif
! 
! if(zplane_t%calc) then
! 
!   do k=1,zplane_t%nloc
! 
!   $if ($MPI)
!   if(zplane_t%coord(k) == coord) then
!   $endif
! 
!     write(ct,'(F9.4)') zplane_t%loc(k)
!     write(fname,*) 'output/uvw_avg.z-',trim(adjustl(ct)),'.dat'
!     fname=trim(adjustl(fname))
! 
!     open (unit = 2,file = fname, status='unknown',form='formatted', &
!       action='write',position='rewind')
!     write(2,*) 'variables = "x", "y", "u", "v", "w"';
!     write(2,"(1a,i9,1a,i3,1a,i3,1a,i3)") 'ZONE T="', &
!       j,'", DATAPACKING=POINT, i=', Nx,', j=',Ny
!     write(2,"(1a)") ''//adjustl('DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)')//''
! !  Write uvw for each jyp plane
!     do j=1,ny
!       do i=1,nx
!         write(2,*) x(i), y(j), zplane_t%ua(i,j,k), zplane_t%va(i,j,k), zplane_t%wa(i,j,k)
!       enddo
!     enddo
!     close(2)
! 
!     $if ($MPI)
!     endif
!     $endif
! 
!   enddo
! 
! endif
! 
! 
! return
! end subroutine plane_write

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine inflow_write ()
use param, only : jt_total, jt_start_write, buff_end,  &
                  read_inflow_file, write_inflow_file
use sim_param, only : u, v, w
implicit none

character (*), parameter :: sub = 'inflow_write'
character (*), parameter :: inflow_file = 'output/inflow_BC.out'
character (*), parameter :: field_file = 'output/inflow.vel.out'
character (*), parameter :: MPI_suffix = '.c'

integer, parameter :: lun = 80
integer, parameter :: field_lun = 81

$if ($DEBUG)
logical, parameter :: DEBUG = .false.
$endif
character (64) :: fname

integer, save :: rec = 0
integer :: nrec
integer :: iolen
integer :: iend, iend_w

logical, save :: initialized = .false.
logical :: opn, exst

!---------------------------------------------------------------------

!--option check
if ( read_inflow_file .and. write_inflow_file ) then
  write (*, *) sub // ': cannot have read_inflow_file and write_inflow_file'
  stop
end if

!--check consistency with inflow_cond
iend = floor (buff_end * nx + 1._rprec)
iend_w = modulo (iend - 1, nx) + 1

if (.not. initialized) then

  inquire ( unit=lun, exist=exst, opened=opn )
  if ( .not. exst ) then
    write (*, *) sub // ': lun = ', lun, ' does not exist'
    stop
  end if
  if (opn) then
    write (*, *) sub // ': lun = ', lun, ' is already open'
    stop
  end if

  if ( USE_MPI ) then
      write ( fname, '(a,a,i0)' ) trim (inflow_file), MPI_suffix, coord
  else
      write ( fname, '(a)' ) inflow_file
  end if
  
  inquire ( file=fname, exist=exst, opened=opn )
  if (exst .and. opn) then
    write (*, *) sub // ': file = ', trim (fname), ' is already open'
    stop
  end if
  
  !--figure out the record length
  inquire (iolength=iolen) u(iend_w, :, :), v(iend_w, :, :), w(iend_w, :, :)

  !--always add to existing inflow_file
  !--inflow_file records always start at 1
  if ( exst ) then

      !--figure out the number of records already in file
      call len_da_file (fname, iolen, nrec)

      write (*, *) sub // ': #records in ' // trim (fname) // '= ', nrec

      rec = nrec

  else

      rec = 0

  end if
  
  !--using direct-access file to allow implementation of 'inflow recycling'
  !  more easily
  !--may want to put in some simple checks on ny, nz
  open (unit=lun, file=fname, access='direct', action='write',  &
        recl=iolen)

  initialized = .true.

end if

if (jt_total == jt_start_write) then  !--write entire flow field out
  inquire (unit=field_lun, exist=exst, opened=opn)
  if (exst .and. .not. opn) then
    open (unit=field_lun, file=field_file, form='unformatted')
    call output_final (jt_total, field_lun)
  else
    write (*, *) sub // ': problem opening ' // field_file
    stop
  end if
end if

if (jt_total >= jt_start_write) then
  rec = rec + 1
  write (unit=lun, rec=rec) u(iend_w, :, :), v(iend_w, :, :), w(iend_w, :, :)
  $if ($DEBUG)
  if ( DEBUG ) write (*, *) sub // ': wrote record ', rec
  $endif
end if

end subroutine inflow_write

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine inflow_read ()
use param, only : ny, nz, pi, nsteps, jt_total, buff_end
use sim_param, only : u, v, w
implicit none

character (*), parameter :: sub = 'inflow_read'
character (*), parameter :: inflow_file = 'output/inflow_BC.out'
character (*), parameter :: debug_file = 'inflow_read_debug.dat'
character (*), parameter :: MPI_suffix = '.c'

integer, parameter :: lun = 80  !--inflow_write for now
integer, parameter :: lun_DEBUG = 88

integer, parameter :: l_blend = 300  !--length of blending zone (recycling)
                                     !--should correspond to integral scale
                                     !--this is number of t-steps
logical, parameter :: recycle = .false.

$if ($DEBUG)
logical, parameter :: DEBUG = .false.
$endif
character (32) :: fmt
character (64) :: fname

!--check for consistency with sim_param here
!--could define a fortran integer lbz in sim_param, and make it visible
!  here, however, this may have complications elsewhere where the name lbz
!  is used.
$if ( $MPI )
    $define $lbz 0
$else
    $define $lbz 1
$endif

integer :: jy, jz
integer :: iend, iend_w
integer :: i
integer :: iolen
integer, save :: rec
integer, save :: nrec
integer :: recp

logical, save :: init_DEBUG = .false.
logical, save :: initialized = .false.
logical :: exst, opn

real (rprec) :: wgt

real (rprec) :: u_tmp(ny, $lbz:nz), v_tmp(ny, $lbz:nz), w_tmp(ny, $lbz:nz)

!---------------------------------------------------------------------

iend = floor ( buff_end * nx + 1.0_rprec )
iend_w = modulo ( iend - 1, nx ) + 1

if ( .not. initialized ) then

    inquire ( unit=lun, exist=exst, opened=opn )
    if ( .not. exst ) then
        write (*, *) sub // ': lun = ', lun, ' does not exist'
        stop
    end if
    if ( opn ) then
        write (*, *) sub // ': lun = ', lun, ' is already open'
        stop
    end if

    if ( USE_MPI ) then
        write ( fname, '(a,a,i0)' ) trim (inflow_file), MPI_suffix, coord
    else
        write ( fname, '(a)' ) inflow_file
    end if
    
    inquire ( file=fname, exist=exst, opened=opn )
    if ( exst ) then
        if ( opn ) then
            write (*, *) sub // ': file = ', fname, ' is already open'
            stop
        end if
    else
        write (*, *) sub // ': file = ', fname, ' does not exist'
        stop
    end if

    !--can only reach this point if exst and .not. opn
  
    !--figure out the record length
    inquire ( iolength=iolen ) u(iend_w, :, :), v(iend_w, :, :), w(iend_w, :, :)
  
    !--figure out the number of records
    call len_da_file ( fname, iolen, nrec )

    write (*, *) sub // ': number of records = ', nrec

    if ( recycle ) then
        !--check minimum length requirement
        !  checks that there are some points that will be non-blended
        
        if ( 2 * (l_blend - 1) > nrec ) then
            write (*, *) sub // ': ', fname, 'is too short to recycle'
            stop
        end if
    end if

    open ( unit=lun, file=fname, access='direct', action='read',  &
           recl=iolen )

    !--file always starts a record 1, but in continued runs, we may need to
    !  access a record that is not 1 to begin with
    !--actually, with wrap-around of records, it means the reading can start
    !  at any point in the file and be OK
    !--intended use: jt_total = 1 here at start of set of runs reading
    !  from the inflow_file, so the first record read will be record 1
    rec = jt_total - 1

    initialized = .true.

end if

rec = rec + 1

if ( recycle ) then
    rec = modulo ( rec - 1, nrec - l_blend + 1 ) + 1
else
    rec = modulo ( rec - 1, nrec ) + 1
end if

read ( unit=lun, rec=rec ) u(iend_w, :, :), v(iend_w, :, :), w(iend_w, :, :)

$if ($DEBUG)
if ( DEBUG ) write (*, *) sub // ' : read record ', rec
$endif
    
if ( recycle ) then

    if ( rec < l_blend ) then

        recp = nrec - l_blend + 1 + rec

        wgt = 0.5_rprec * ( 1.0_rprec -                              &
                            cos ( pi * real (rec, rprec) / l_blend ) )
            !--wgt = 0+ when rec = 1
            !  wgt = 1- when rec = l_blend

        read ( unit=lun, rec=recp ) u_tmp, v_tmp, w_tmp

        u(iend_w, :, :) = wgt * u(iend_w, :, :) + (1.0_rprec - wgt) * u_tmp
        v(iend_w, :, :) = wgt * v(iend_w, :, :) + (1.0_rprec - wgt) * v_tmp
        w(iend_w, :, :) = wgt * w(iend_w, :, :) + (1.0_rprec - wgt) * w_tmp

    end if

end if

$if ($DEBUG)
if ( DEBUG ) then  !--write out slices as an ascii time series

    if ( .not. init_DEBUG ) then

        inquire ( unit=lun_DEBUG, exist=exst, opened=opn )

        if ( exst .and. (.not. opn) ) then
        
            if ( USE_MPI ) then
                open ( unit=lun_DEBUG, file=debug_file // MPI_suffix )
            else
                open ( unit=lun_DEBUG, file=debug_file )
            end if
        
            write ( lun_DEBUG, '(a)' ) 'variables = "y" "z" "t" "u" "v" "w"'
            write ( lun_DEBUG, '(3(a,i0))' ) 'zone, f=point, i= ', ny,  &
                                             ', j= ', nz,               &
                                             ', k= ', nsteps

        else
            write (*, *) sub // ': problem opening debug file'
            stop
        end if

        init_DEBUG = .true.

    end if

    fmt = '(3(1x,i0),3(1x,es12.5))'
    do jz = 1, nz
        do jy = 1, ny
            write ( lun_DEBUG, fmt ) jy, jz, jt_total, u(iend_w, jy, jz),  &
                                     v(iend_w, jy, jz), w(iend_w, jy, jz)
        end do
    end do

end if
$endif

end subroutine inflow_read

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--finds number of records on existing direct-access unformatted file
!--taken from Clive Page's comp.lang.fortran posting (12/16/2003), 
!  under the topic counting number of records in a Fortran direct file
!--minor changes/renaming
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine len_da_file(fname, lenrec, length)
implicit none
character (*), intent(in) :: fname  ! name of existing direct-access file
integer, intent(in)       :: lenrec ! record length (O/S dependent units)
integer, intent(out) :: length      ! number of records.
!
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
  print *, 'error in len_da_file: ', trim(fname), ' does not exist'
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!**********************************************************************
subroutine stats_init ()
!***************************************************************
!  This subroutine allocates the memory for arrays
!  used for statistical calculations 

use param, only : L_x,L_y,L_z,dx,dy,dz,nx,ny,nz,nsteps,coord,nproc
use stat_defs
use grid_defs
use functions, only : index_start
implicit none

!character(120) :: cx,cy,cz
character(120) :: fname, var_list
integer :: fid, i,j,k

!  All nstart and nend values are based
!  on jt and not jt_total
tstats_t%calc = .true.
tstats_t%nstart = 1
tstats_t%nend = nsteps

!  Turns instantaneous velocity recording on or off
point_t%calc = .false.
point_t%nstart = 1
point_t%nend   = nsteps
point_t%nskip = 1
point_t%nloc = 2
point_t%xyz(:,1) = (/L_x/2., L_y/2., 1.5_rprec/)
point_t%xyz(:,2) = (/L_x/2., L_y/2., 2.5_rprec/)

domain_t%calc = .false.
domain_t%nstart = 1000
domain_t%nend   = nsteps
domain_t%nskip = 1000

!  y-plane stats/data
yplane_t%calc   = .false.
yplane_t%nstart = 1
yplane_t%nend   = nsteps
yplane_t%nskip  = 100
yplane_t%nloc   = 2
yplane_t%loc(1) = 1.0
yplane_t%loc(2) = 3.0

!  z-plane stats/data
zplane_t%calc   = .false.
zplane_t%nstart = 100
zplane_t%nend   = nsteps
zplane_t%nskip  = 100
zplane_t%nloc   = 3
zplane_t%loc(1) = 1.95925
zplane_t%loc(2) = 2.1636
zplane_t%loc(3) = 2.26515

!  z-plane TIME-AVERAGED stats/data
zplane_avg_t%calc   = .false.
zplane_avg_t%nstart = 1
zplane_avg_t%nend   = nsteps


$if ($MPI)
  !--this dimensioning adds a ghost layer for finite differences
  !--its simpler to have all arrays dimensioned the same, even though
  !  some components do not need ghost layer
  $define $lbz 0
$else
  $define $lbz 1
$endif

!  Allocate arrays for variable summation for Reynolds
!  stress calculations
if(tstats_t%calc) then 

  allocate(tstats_t%u(nx, ny, nz))
  allocate(tstats_t%v(nx, ny, nz))
  allocate(tstats_t%w(nx, ny, nz))
  allocate(tstats_t%u2(nx, ny, nz))
  allocate(tstats_t%v2(nx, ny, nz))
  allocate(tstats_t%w2(nx, ny, nz))
  allocate(tstats_t%uw(nx, ny, nz))
  allocate(tstats_t%vw(nx, ny, nz))
  allocate(tstats_t%uv(nx, ny, nz))
  allocate(tstats_t%dudz(nx, ny, nz))
  
  !  Initialize arrays
  tstats_t%u=0.
  tstats_t%v=0.
  tstats_t%w=0.
  tstats_t%u2=0.
  tstats_t%v2=0.
  tstats_t%w2=0.
  tstats_t%uw=0.
  tstats_t%vw=0.
  tstats_t%uv=0.
  tstats_t%dudz=0.
  
  $if($LVLSET)
  $if($RNS_LS)
  
  allocate(tstats_t%fx(nx,ny,nz))
  allocate(tstats_t%fy(nx,ny,nz))
  allocate(tstats_t%fz(nx,ny,nz))
  
  tstats_t%fx = 0._rprec
  tstats_t%fx = 0._rprec
  tstats_t%fx = 0._rprec
  
  $endif
  $endif
  
endif

! if(rs_t%calc) then
!   allocate(rs_t%up2(nx, ny, $lbz:nz))
!   allocate(rs_t%vp2(nx, ny, $lbz:nz))
!   allocate(rs_t%wp2(nx, ny, $lbz:nz))
!   allocate(rs_t%upwp(nx, ny, $lbz:nz))
!   allocate(rs_t%vpwp(nx, ny, $lbz:nz))
!   allocate(rs_t%upvp(nx, ny, $lbz:nz))
!   rs_t%up2=0.
!   rs_t%vp2=0.
!   rs_t%wp2=0.
!   rs_t%upwp=0.
!   rs_t%vpwp=0.
!   rs_t%upvp=0.
! endif

! Initialize information for y-planar stats/data
if(yplane_t%calc) then
!   allocate(yplane_t%ua(Nx,yplane_t%nloc,Nz))
!   allocate(yplane_t%va(Nx,yplane_t%nloc,Nz))
!   allocate(yplane_t%wa(Nx,yplane_t%nloc,Nz))
  
!   yplane_t%fa = 1./(dble(yplane_t%nend - yplane_t%nstart + 1))
!   yplane_t%ua = 0.
!   yplane_t%va = 0.
!   yplane_t%wa = 0.
  yplane_t%istart = -1
  yplane_t%ldiff = 0.
!  Not really needed
  yplane_t%coord = -1
  
!  Compute istart and ldiff
  do j=1,yplane_t%nloc
	  yplane_t%istart(j) = index_start('j', dy, yplane_t%loc(j))
	  yplane_t%ldiff(j) = y(yplane_t%istart(j)) - yplane_t%loc(j)
  enddo
    
endif

! Initialize information for y-planar stats/data
if(zplane_t%calc) then
!   allocate(zplane_t%ua(Nx,Ny,zplane_t%nloc))
!   allocate(zplane_t%va(Nx,Ny,zplane_t%nloc))
!   allocate(zplane_t%wa(Nx,Ny,zplane_t%nloc))
!   zplane_t%fa = 1./(dble(zplane_t%nend - zplane_t%nstart + 1))
!   zplane_t%ua = 0.
!   zplane_t%va = 0.
!   zplane_t%wa = 0.

!  Initialize 
  zplane_t%istart = -1
  zplane_t%ldiff = 0. 
  zplane_t%coord=-1 
  
!  Compute istart and ldiff
  do k=1,zplane_t%nloc

    $if ($MPI)
    if(zplane_t%loc(k) >= z(1) .and. zplane_t%loc(k) < z(nz)) then
      zplane_t%coord(k) = coord

	  zplane_t%istart(k) = index_start('k',dz,zplane_t%loc(k))
	  zplane_t%ldiff(k) = z(zplane_t%istart(k)) - zplane_t%loc(k)
    endif
    $else
    zplane_t%coord(k) = 0
    zplane_t%istart(k) = index_start('k',dz,zplane_t%loc(k))
    zplane_t%ldiff(k) = z(zplane_t%istart(k)) - zplane_t%loc(k)
    $endif

  enddo  
  
endif

!  Intialize the coord values (-1 shouldn't be used as coord so initialize to this)
point_t%coord=-1

!  Open files for instantaneous writing
if(point_t%calc) then

  do i=1,point_t%nloc
!  Find the processor in which this point lives
  $if ($MPI)
    if(point_t%xyz(3,i) >= z(1) .and. point_t%xyz(3,i) < z(nz)) then
    point_t%coord(i) = coord
	  
	  point_t%istart(i) = index_start('i',dx,point_t%xyz(1,i))
	  point_t%jstart(i) = index_start('j',dy,point_t%xyz(2,i))
	  point_t%kstart(i) = index_start('k',dz,point_t%xyz(3,i))
	  
	  point_t%xdiff(i) = x(point_t%istart(i)) - point_t%xyz(1,i)
	  point_t%ydiff(i) = y(point_t%jstart(i)) - point_t%xyz(2,i)
	  point_t%zdiff(i) = z(point_t%kstart(i)) - point_t%xyz(3,i)

	  
    fid=3000*i
	  
	  !  Can't concatenate an empty string
    point_t%fname(i)=''
	  call strcat(point_t%fname(i),'output/uvw_inst-')
	  call strcat(point_t%fname(i), point_t%xyz(1,i))
	  call strcat(point_t%fname(i),'-')
	  call strcat(point_t%fname(i),point_t%xyz(2,i))
	  call strcat(point_t%fname(i),'-')
	  call strcat(point_t%fname(i),point_t%xyz(3,i))
	  call strcat(point_t%fname(i),'.dat')
	
	  
	  !call strcat(point_t%fname(i),point_t%xyz(1,i))
      !write (point_t%fname(i),*) ,trim(adjustl(cx)),'-',  &
      !trim(adjustl(cy)),'-', trim(adjustl(cz)),'.dat'
	   
	  var_list = '"t (s)", "u", "v", "w"'
	  call write_tecplot_header_xyline(point_t%fname(i), 'rewind', var_list)
	    
    endif
  $else
    point_t%coord(i) = 0
	point_t%istart(i) = index_start('i',dx,point_t%xyz(1,i))
	point_t%jstart(i) = index_start('j',dy,point_t%xyz(2,i))
	point_t%kstart(i) = index_start('k',dz,point_t%xyz(3,i))
	
	point_t%xdiff(i) = x(point_t%istart(i)) - point_t%xyz(1,i)
	point_t%ydiff(i) = y(point_t%jstart(i)) - point_t%xyz(2,i)
	point_t%zdiff(i) = z(point_t%kstart(i)) - point_t%xyz(3,i)

    !write(cx,'(F9.4)') point_t%xyz(1,i)
    !write(cy,'(F9.4)') point_t%xyz(2,i)
    !write(cz,'(F9.4)') point_t%xyz(3,i)
	
	fid=3000*i
	
	call strcat(point_t%fname(i),'output/uvw_inst-')
	call strcat(point_t%fname(i), point_t%xyz(1,i))
	call strcat(point_t%fname(i),'-')
	call strcat(point_t%fname(i),point_t%xyz(2,i))
	call strcat(point_t%fname(i),'-')
	call strcat(point_t%fname(i),point_t%xyz(3,i))
	call strcat(point_t%fname(i),'.dat')
	
    !write (point_t%fname(i),*) 'output/uvw_inst-',trim(adjustl(cx)),'-',  &
    !  trim(adjustl(cy)),'-', trim(adjustl(cz)),'.dat'
	var_list = '"t (s)", "u", "v", "w"'
	call write_tecplot_header_xyline(point_t%fname(i), 'rewind', var_list)
	
  $endif
  
  enddo
endif

! Initialize for zplane_avg calculations
if(zplane_avg_t%calc) then 
  allocate(zplane_avg_t%u_avg(nz-1))
  allocate(zplane_avg_t%v_avg(nz-1))
  allocate(zplane_avg_t%w_avg(nz-1))
  allocate(zplane_avg_t%u2_avg(nz-1))
  allocate(zplane_avg_t%v2_avg(nz-1))
  allocate(zplane_avg_t%w2_avg(nz-1))  
  allocate(zplane_avg_t%txx_avg(nz-1))
  allocate(zplane_avg_t%txy_avg(nz-1))
  allocate(zplane_avg_t%tyy_avg(nz-1))
  allocate(zplane_avg_t%txz_avg(nz-1))
  allocate(zplane_avg_t%tyz_avg(nz-1))
  allocate(zplane_avg_t%tzz_avg(nz-1))  
  allocate(zplane_avg_t%dudz_avg(nz-1))
  allocate(zplane_avg_t%dvdz_avg(nz-1))
  allocate(zplane_avg_t%uv_avg(nz-1))
  allocate(zplane_avg_t%uw_avg(nz-1))
  allocate(zplane_avg_t%vw_avg(nz-1))  

  !  Initialize arrays
  zplane_avg_t%u_avg(:)=0.
  zplane_avg_t%v_avg(:)=0.
  zplane_avg_t%w_avg(:)=0.
  zplane_avg_t%u2_avg(:)=0.
  zplane_avg_t%v2_avg(:)=0.
  zplane_avg_t%w2_avg(:)=0.  
  zplane_avg_t%txx_avg(:)=0.
  zplane_avg_t%txy_avg(:)=0.
  zplane_avg_t%tyy_avg(:)=0.
  zplane_avg_t%txz_avg(:)=0.
  zplane_avg_t%tyz_avg(:)=0.
  zplane_avg_t%tzz_avg(:)=0.  
  zplane_avg_t%dudz_avg(:)=0.
  zplane_avg_t%dvdz_avg(:)=0.  
  zplane_avg_t%uv_avg(:)=0.
  zplane_avg_t%uw_avg(:)=0.
  zplane_avg_t%vw_avg(:)=0.    
endif

return
end subroutine stats_init

!***************************************************************
subroutine tstats_compute()

!***************************************************************
!  This subroutine collects the stats for each flow 
!  variable quantity
use types, only : rprec
use stat_defs, only : tstats_t
use param, only : nx,ny,nz, dt
use sim_param, only : u,v,w, dudz
$if($LVLSET)
$if($RNS_LS)
use immersedbc, only : fx, fy, fz
$endif
$endif

implicit none

!use io, only : w_uv, w_uv_tag, dudz_uv, dudz_uv_tag, interp_to_uv_grid
integer :: i,j,k
real(rprec) :: u_p, v_p, w_p, dudz_p

!  Make sure w has been interpolated to uv-grid
call interp_to_uv_grid(w, w_uv, w_uv_tag)
call interp_to_uv_grid(dudz, dudz_uv, dudz_uv_tag)

do k=1,nz
  do j=1,ny
    do i=1,nx
!  Being cache friendly
      u_p = u(i,j,k) * dt 
      v_p = v(i,j,k) * dt 
!  Interpolate each w and dudz to uv grid
      w_p = w_uv(i,j,k) * dt
      dudz_p = dudz_uv(i,j,k) * dt 
          
      tstats_t%u(i,j,k)=tstats_t%u(i,j,k) + u_p 
      tstats_t%v(i,j,k)=tstats_t%v(i,j,k) + v_p 
      tstats_t%w(i,j,k)=tstats_t%w(i,j,k) + w_p 
      tstats_t%u2(i,j,k)=tstats_t%u2(i,j,k) + u_p*u_p 
      tstats_t%v2(i,j,k)=tstats_t%v2(i,j,k) + v_p*v_p 
      tstats_t%w2(i,j,k)=tstats_t%w2(i,j,k) + w_p*w_p 
      tstats_t%uw(i,j,k)=tstats_t%uw(i,j,k) + u_p*w_p 
      tstats_t%vw(i,j,k)=tstats_t%vw(i,j,k) + v_p*w_p 
      tstats_t%uv(i,j,k)=tstats_t%uv(i,j,k) + u_p*v_p 
      tstats_t%dudz(i,j,k)=tstats_t%dudz(i,j,k) + dudz_p 
      
      $if($LVLSET)
      $if($RNS_LS)
      tstats_t % fx(i,j,k) = tstats_t % fx(i,j,k) + fx(i,j,k) * dt 
      tstats_t % fy(i,j,k) = tstats_t % fy(i,j,k) + fy(i,j,k) * dt 
      tstats_t % fz(i,j,k) = tstats_t % fz(i,j,k) + fz(i,j,k) * dt 
      $endif
      $endif
      
      
    enddo
  enddo
enddo

return

end subroutine tstats_compute

!!***************************************************************
!subroutine tavg_compute()
!!***************************************************************
!!  This subroutine collects the stats for each flow 
!!  variable quantity
!use stat_defs, only : tavg_t
!use param, only : nx,ny,nz
!use sim_param, only : u,v,w, dudz
!!use io, only : w_uv, w_uv_tag, dudz_uv, dudz_uv_tag, interp_to_uv_grid
!implicit none
!integer :: i,j,k, navg
!double precision :: w_interp, dudz_interp, fa

!!  Make sure w has been interpolated to uv-grid
!call interp_to_uv_grid(w, w_uv, w_uv_tag)
!call interp_to_uv_grid(dudz, dudz_uv, dudz_uv_tag)

!!  Initialize w_interp and dudz_interp
!w_interp = 0.
!dudz_interp = 0. 

!!  Compute number of times to average over
!navg = tavg_t%nend - tavg_t%nstart + 1
!!  Averaging factor
!fa=1./dble(navg)

!do k=1,nz
!  do j=1,ny
!    do i=1,nx
!!  Interpolate each w and dudz to uv grid
!      w_interp = w_uv(i,j,k)  
!      dudz_interp = dudz_uv(i,j,k)
!      
!      tavg_t%u(i,j,k)=tavg_t%u(i,j,k) + fa*u(i,j,k)
!      tavg_t%v(i,j,k)=tavg_t%v(i,j,k) + fa*v(i,j,k)
!      tavg_t%w(i,j,k)=tavg_t%w(i,j,k) + fa*w_interp
!      tavg_t%u2(i,j,k)=tavg_t%u2(i,j,k) + fa*u(i,j,k)*u(i,j,k)
!      tavg_t%v2(i,j,k)=tavg_t%v2(i,j,k) + fa*v(i,j,k)*v(i,j,k)
!      tavg_t%w2(i,j,k)=tavg_t%w2(i,j,k) + fa*w_interp*w_interp
!      tavg_t%uw(i,j,k)=tavg_t%uw(i,j,k) + fa*u(i,j,k)*w_interp
!      tavg_t%vw(i,j,k)=tavg_t%vw(i,j,k) + fa*v(i,j,k)*w_interp
!      tavg_t%uv(i,j,k)=tavg_t%uv(i,j,k) + fa*u(i,j,k)*v(i,j,k)
!      tavg_t%dudz(i,j,k)=tavg_t%dudz(i,j,k) + fa*dudz_interp
!    enddo
!  enddo
!enddo

!return

!end subroutine tavg_compute

!!***************************************************************
!subroutine rs_compute()
!!***************************************************************
!!  This subroutine computes Reynolds stresses from tavg_t u,v,w
!!  values
!use stat_defs, only : rs_t, tavg_t
!use param, only : nx,ny,nz,dx,dy,dz,L_x,L_y,L_z
!use sim_param, only : u,v,w

!implicit none

!integer :: i,j,k

!do k=1,nz
!  do j=1,ny
!    do i=1,nx
!      rs_t%up2(i,j,k)=tavg_t%u2(i,j,k) - tavg_t%u(i,j,k)*tavg_t%u(i,j,k)
!      rs_t%vp2(i,j,k)=tavg_t%v2(i,j,k) - tavg_t%v(i,j,k)*tavg_t%v(i,j,k)
!	  rs_t%wp2(i,j,k)=tavg_t%w2(i,j,k) - tavg_t%w(i,j,k)*tavg_t%w(i,j,k)
!	  rs_t%upwp(i,j,k)=tavg_t%uw(i,j,k) - tavg_t%u(i,j,k)*tavg_t%w(i,j,k)
!	  rs_t%vpwp(i,j,k)=tavg_t%vw(i,j,k) - tavg_t%v(i,j,k)*tavg_t%w(i,j,k)
!	  rs_t%upvp(i,j,k)=tavg_t%uv(i,j,k) - tavg_t%u(i,j,k)*tavg_t%v(i,j,k)
!	enddo
!  enddo
!enddo
!  
!return
!end subroutine rs_compute

! !**********************************************************************
! subroutine plane_avg_compute(jt)
! !**********************************************************************
! use functions, only : linear_interp, interp_to_uv_grid
! use grid_defs, only : z
! use param, only : dx, dy, dz, Nx, Ny, Nz,coord
! use stat_defs,only:yplane_t, zplane_t
! use sim_param, only : u, v, w
! 
! implicit none
! 
! integer,intent(in)::jt
! 
! integer :: i,j,k
! double precision :: ui, vi ,wi
! 
! !  Determine if y-plane averaging should be performed
! if(yplane_t%calc .and. jt >= yplane_t%nstart .and. jt <= yplane_t%nend) then
!   
! !  Compute average for each y-plane (jya)
!   do j=1,yplane_t%nloc
!     do k=1,Nz
!       do i=1,Nx
!         ui = linear_interp(u(i,yplane_t%istart(j),k), &
!           u(i,yplane_t%istart(j)+1,k), dy, yplane_t%ldiff(j))
!         vi = linear_interp(v(i,yplane_t%istart(j),k), &
!           v(i,yplane_t%istart(j)+1,k), dy, yplane_t%ldiff(j))
!         wi = linear_interp(interp_to_uv_grid('w',i,yplane_t%istart(j),k), &
!           interp_to_uv_grid('w',i,yplane_t%istart(j)+1,k), dy, &
!           yplane_t%ldiff(j))
!         yplane_t%ua(i,j,k) = yplane_t%ua(i,j,k) + yplane_t%fa*ui 
!         yplane_t%va(i,j,k) = yplane_t%va(i,j,k) + yplane_t%fa*vi
!         yplane_t%wa(i,j,k) = yplane_t%wa(i,j,k) + yplane_t%fa*wi
! 
!       enddo
!     enddo
!   enddo
! endif 
! 
! !  Determine if y-plane averaging should be performed
! if(zplane_t%calc .and. jt >= zplane_t%nstart .and. jt <= zplane_t%nend) then
! !  Compute average for each y-plane (jya)
! 
!   do k=1,zplane_t%nloc
! 
!   $if ($MPI)
!   if(zplane_t%coord(k) == coord) then
!   $endif
! 
!     do j=1,Ny
!       do i=1,Nx
! 
!         ui = linear_interp(u(i,j,zplane_t%istart(k)),u(i,j,zplane_t%istart(k)+1), &
!           dz, zplane_t%ldiff(k))
!         vi = linear_interp(v(i,j,zplane_t%istart(k)),v(i,j,zplane_t%istart(k)+1), &
!           dz, zplane_t%ldiff(k))
!         wi = linear_interp(interp_to_uv_grid('w',i,j,zplane_t%istart(k)), &
!           interp_to_uv_grid('w',i,j,zplane_t%istart(k)+1), &
!           dz, zplane_t%ldiff(k))
! 
!         zplane_t%ua(i,j,k) = zplane_t%ua(i,j,k) + zplane_t%fa*ui
!         zplane_t%va(i,j,k) = zplane_t%va(i,j,k) + zplane_t%fa*vi
!         zplane_t%wa(i,j,k) = zplane_t%wa(i,j,k) + zplane_t%fa*wi
! 
!       enddo
!     enddo
! 
!     $if ($MPI)
!     endif
!     $endif
! 
!   enddo
! 
! endif
! 
! return
! end subroutine plane_avg_compute

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!subroutine strcat(str1, str2)
!implicit none

!character(*), intent(INOUT) :: str1
!character(*), intent(IN) :: str2

!str1 = trim(adjustl(str1)) // str2

!return
!end subroutine strcat


!***************************************************************
subroutine zplane_avg_compute(ipart)
!***************************************************************
!  This subroutine averages the following parameters over z-planes:
!  u, v, w, u^2, v^2, w^2, txx, txy, tyy, txz, tyz, tzz, dudz, dvdz
!  and writes them to files avg_*.dat

! ipart=1 to create files and write Tecplot headers
! ipart=2 to compute averages
! ipart=3 to write averages and reset


use param, only : dx, dy, dz, jt_total, dt_dim
use stat_defs, only: zplane_avg_t
use sim_param, only : u, v, w, txx, txy, tyy, txz, tyz, tzz, dudz, dvdz
use grid_defs, only: z

implicit none

integer :: k, ipart, nt_z
character(64) :: frmt
character (64) :: fname, temp, fname2, temp2, fname3, temp3, fname4, temp4, fname5, temp5, fname6, temp6

if (ipart==1) then
!Create files and write Tecplot headers         
		!parameter file
		  if(coord==0) then 
		    fname='output/avg.dat'
			open (unit = 2,file = fname, status='unknown',form='formatted', action='write',position='rewind')
			write (frmt, '("(",i0,"i9)")') 3
			write(2,*) 'nz, nstart, nend'
			write(2,frmt) nz,zplane_avg_t%nstart, zplane_avg_t%nend
			close(2)
			!call write_real_data(fname, 'append', 3, (/ zplane_avg_t%nstart, zplane_avg_t%nend, zplane_avg_t%nskip /))
		  endif
		!u-v-w file  
		  fname2 = 'output/avg_uvw.dat'
			$if ($MPI)
			  write (temp2, '(".c",i0)') coord
			  fname2 = trim (fname2) // temp2
			$endif	 	  
		  call write_tecplot_header_ND(fname2, 'rewind', 4, (/ Nz/),'"z", "u_avg", "v_avg", "w_avg"', coord, 2, total_time_dim)  
		!txx-txy,etc file
		  fname3 = 'output/avg_tau.dat'
			$if ($MPI)
			  write (temp3, '(".c",i0)') coord
			  fname3 = trim (fname3) // temp3
			$endif	 	  
		  call write_tecplot_header_ND(fname3, 'rewind', 7, (/ Nz/), &
		  '"z", "txx_avg", "txy_avg", "tyy_avg", "txz_avg", "tyz_avg", "tzz_avg"', coord, 2, total_time_dim)  
		 !ddz file  
		  fname4 = 'output/avg_ddz.dat'
			$if ($MPI)
			  write (temp4, '(".c",i0)') coord
			  fname4 = trim (fname4) // temp4
			$endif	 	  
		  call write_tecplot_header_ND(fname4, 'rewind', 3, (/ Nz/),'"z", "dudz_avg", "dvdz_avg"', coord, 2, total_time_dim) 
		!u2-v2-w2 file  
		  fname5 = 'output/avg_uvw2.dat'
			$if ($MPI)
			  write (temp5, '(".c",i0)') coord
			  fname5 = trim (fname5) // temp5
			$endif	 	  
		  call write_tecplot_header_ND(fname5, 'rewind', 4, (/ Nz/),'"z", "u2_avg", "v2_avg", "w2_avg"', coord, 2, total_time_dim) 		  
		!uv-uw-vw file  
		  fname6 = 'output/avg_uv_uw_vw.dat'
			$if ($MPI)
			  write (temp6, '(".c",i0)') coord
			  fname6 = trim (fname6) // temp6
			$endif	 	  
		  call write_tecplot_header_ND(fname6, 'rewind', 4, (/ Nz/),'"z", "uv_avg", "uw_avg", "vw_avg"', coord, 2, total_time_dim) 		  

elseif (ipart==2) then
!Compute averages
	nt_z = zplane_avg_t%nend - zplane_avg_t%nstart + 1
    do k=1,nz-1
	    zplane_avg_t%u_avg(k) = zplane_avg_t%u_avg(k) + sum(u(:,:,k))/(nx*ny*nt_z)
	    zplane_avg_t%v_avg(k) = zplane_avg_t%v_avg(k) + sum(v(:,:,k))/(nx*ny*nt_z)
	    zplane_avg_t%w_avg(k) = zplane_avg_t%w_avg(k) + sum(w(:,:,k))/(nx*ny*nt_z)
		
		zplane_avg_t%txx_avg(k) = zplane_avg_t%txx_avg(k) + sum(txx(:,:,k))/(nx*ny*nt_z)
		zplane_avg_t%txy_avg(k) = zplane_avg_t%txy_avg(k) + sum(txy(:,:,k))/(nx*ny*nt_z)
		zplane_avg_t%tyy_avg(k) = zplane_avg_t%tyy_avg(k) + sum(tyy(:,:,k))/(nx*ny*nt_z)
		zplane_avg_t%txz_avg(k) = zplane_avg_t%txz_avg(k) + sum(txz(:,:,k))/(nx*ny*nt_z)
		zplane_avg_t%tyz_avg(k) = zplane_avg_t%tyz_avg(k) + sum(tyz(:,:,k))/(nx*ny*nt_z)
		zplane_avg_t%tzz_avg(k) = zplane_avg_t%tzz_avg(k) + sum(tzz(:,:,k))/(nx*ny*nt_z)
		
		zplane_avg_t%dudz_avg(k) = zplane_avg_t%dudz_avg(k) + sum(dudz(:,:,k))/(nx*ny*nt_z)
	    zplane_avg_t%dvdz_avg(k) = zplane_avg_t%dvdz_avg(k) + sum(dvdz(:,:,k))/(nx*ny*nt_z)
		
		zplane_avg_t%u2_avg(k) = zplane_avg_t%u2_avg(k) + sum(u(:,:,k)*u(:,:,k))/(nx*ny*nt_z)
	    zplane_avg_t%v2_avg(k) = zplane_avg_t%v2_avg(k) + sum(v(:,:,k)*v(:,:,k))/(nx*ny*nt_z)
	    zplane_avg_t%w2_avg(k) = zplane_avg_t%w2_avg(k) + sum(w(:,:,k)*w(:,:,k))/(nx*ny*nt_z)
		
		zplane_avg_t%uv_avg(k) = zplane_avg_t%uv_avg(k) + sum(u(:,:,k)*v(:,:,k))/(nx*ny*nt_z)
	    zplane_avg_t%uw_avg(k) = zplane_avg_t%uw_avg(k) + sum(u(:,:,k)*w(:,:,k))/(nx*ny*nt_z)
	    zplane_avg_t%vw_avg(k) = zplane_avg_t%vw_avg(k) + sum(v(:,:,k)*w(:,:,k))/(nx*ny*nt_z)		
	enddo
elseif (ipart==3) then
!Write averages and reset
	  !u-v-w file
		  call write_real_data_1D(fname2,'append','formatted',3,(nz-1),(/ zplane_avg_t%u_avg,  zplane_avg_t%v_avg,  zplane_avg_t%w_avg /),z(1:nz-1))
		  zplane_avg_t%u_avg(:) = 0.
		  zplane_avg_t%v_avg(:) = 0.
		  zplane_avg_t%w_avg(:) = 0.
	  !tau file
		  call write_real_data_1D(fname3,'append','formatted',6,(nz-1),(/ zplane_avg_t%txx_avg,  &
		  zplane_avg_t%txy_avg, zplane_avg_t%tyy_avg,  zplane_avg_t%txz_avg,  zplane_avg_t%tyz_avg,  zplane_avg_t%tzz_avg /),z(1:nz-1))
		  zplane_avg_t%txx_avg(:) = 0.
		  zplane_avg_t%txy_avg(:) = 0.
		  zplane_avg_t%tyy_avg(:) = 0.
		  zplane_avg_t%txz_avg(:) = 0.
		  zplane_avg_t%tyz_avg(:) = 0.
		  zplane_avg_t%tzz_avg(:) = 0.	 
	  !ddz file
		  call write_real_data_1D(fname4,'append','formatted',2,(nz-1),(/ zplane_avg_t%dudz_avg,  zplane_avg_t%dvdz_avg /),z(1:nz-1))
		  zplane_avg_t%dudz_avg(:) = 0.
		  zplane_avg_t%dvdz_avg(:) = 0.	
	  !u2-v2-w2 file
		  call write_real_data_1D(fname5,'append','formatted',3,(nz-1),(/ zplane_avg_t%u2_avg,  zplane_avg_t%v2_avg,  zplane_avg_t%w2_avg /),z(1:nz-1))
		  zplane_avg_t%u2_avg(:) = 0.
		  zplane_avg_t%v2_avg(:) = 0.
		  zplane_avg_t%w2_avg(:) = 0.	
	  !uv-uw-vw file
		  call write_real_data_1D(fname6,'append','formatted',3,(nz-1),(/ zplane_avg_t%uv_avg,  zplane_avg_t%uw_avg,  zplane_avg_t%vw_avg /),z(1:nz-1))
		  zplane_avg_t%uv_avg(:) = 0.
		  zplane_avg_t%uw_avg(:) = 0.
		  zplane_avg_t%vw_avg(:) = 0.			  
else 
!error
	write(*,*) 'ERROR: Incorrect call zplane_avg_compute(ipart)'
endif


return
end subroutine zplane_avg_compute


end module io
