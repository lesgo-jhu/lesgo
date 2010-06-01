module io
use types,only:rprec
use param, only : ld, nx, ny, nz, nz_tot, write_inflow_file, path,  &
                  USE_MPI, coord, rank, nproc, jt_total, total_time, total_time_dim
use sim_param, only : w, dudz, dvdz
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
real(rprec), dimension(lbound(dvdz,1):ubound(dvdz,1),lbound(dvdz,2):ubound(dvdz,2),lbound(dvdz,3):ubound(dvdz,3)) :: dvdz_uv ! on the uv grid
integer :: w_uv_tag, dudz_uv_tag, dvdz_uv_tag
                   
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
use param, only : tavg_calc, tavg_nstart, tavg_nend
use stat_defs, only: tavg_t, point_t, domain_t, xplane_t, yplane_t, zplane_t
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
if(tavg_calc) then
!  Check if we are in the time interval for running summations
  if(jt >= tavg_nstart .and. jt <= tavg_nend) then
    if(jt == tavg_nstart) then
      if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
        write(*,*) '-------------------------------'   
        write(*,"(1a,i9,1a,i9)") 'Starting running time summation from ', tavg_nstart, ' to ', tavg_nend
        write(*,*) '-------------------------------'   
      endif

      call tavg_init()

    endif
!  Compute running summations
    call tavg_compute ()
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

!  Determine if instantaneous x-plane velocities are to be recorded
if(xplane_t%calc) then
  if(jt >= xplane_t%nstart .and. jt <= xplane_t%nend .and. ( jt == xplane_t%nstart .or. mod(jt,xplane_t%nskip)==0) ) then
    if(.not. xplane_t%started) then
      if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then        
        write(*,*) '-------------------------------'
        write(*,"(1a,i9,1a,i9)") 'Writing instantaneous x-plane velocities from ', xplane_t%nstart, ' to ', xplane_t%nend
        write(*,"(1a,i9)") 'Iteration skip:', xplane_t%nskip
        write(*,*) '-------------------------------'
      endif
      xplane_t%started=.true.
    endif
    call inst_write(5)
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

!!  Determine if time-averaged z-plane data are to be recorded
!if(zplane_avg_t%calc) then
!  if(jt >= zplane_avg_t%nstart .and. jt <= zplane_avg_t%nend) then
!    if(.not. zplane_avg_t%started) then
!      if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then        
!        write(*,*) '-------------------------------'
!        write(*,"(1a,i9,1a,i9)") 'Writing time-averaged z-plane data from ', zplane_avg_t%nstart, ' to ', zplane_avg_t%nend
!        write(*,*) '-------------------------------'
!      endif
!      zplane_avg_t%started=.true.	
!	  !Write Tecplot header files
!	  call zplane_avg_compute(1)
!	endif

!	!Calculate the averaged data for each time step
!	  call zplane_avg_compute(2)

!    !if(jt > zplane_avg_t%nstart .and. mod(jt,zplane_avg_t%nskip)==0) then
!	if(jt == zplane_avg_t%nend) then
!	  !Write the averaged data to file and reset the variables	
!	  call zplane_avg_compute(3)	  
!	endif
!  endif
!endif

!if(yplane_t%calc .or. zplane_t%calc) call plane_avg_compute(jt)

return
end subroutine output_loop

!**********************************************************************
subroutine tavg_init()
!**********************************************************************

!  Load tavg.out files

use param, only : coord, dt, USE_MPI
use messages
use stat_defs, only : tavg_t, tavg_zplane_t, tavg_total_time
use param, only : tavg_nstart, tavg_nend
implicit none

character (*), parameter :: sub_name = mod_name // '.tavg_init'
character (*), parameter :: ftavg_in = 'tavg.out'
$if ($MPI)
character (*), parameter :: MPI_suffix = '.c'
$endif
character (128) :: fname

logical :: opn, exst

!---------------------------------------------------------------------

inquire (unit=1, opened=opn)
if (opn) call error (sub_name, 'unit 1 already open')

$if ($MPI)
write (fname, '(a,a,i0)') ftavg_in, MPI_suffix, coord
$else
fname = trim(adjustl(ftavg_in))
$endif

inquire (file=fname, exist=exst)
if (.not. exst) then
  !  Nothing to read in
  if(.not. USE_MPI .or. (USE_MPI .and. coord == 0)) then
    write(*,*) ' '
    write(*,*)'No previous time averaged data - starting from scratch.'
  endif
  tavg_total_time = dble(tavg_nend - tavg_nstart + 1) * dt

  !tavg_started = .true.
  
  return
endif
  
open (1, file=fname, action='read', position='rewind', form='unformatted')

read (1) tavg_total_time
read (1) tavg_t
read (1) tavg_zplane_t
!!read (1) tavg_t % u 
!!read (1) tavg_t % v
!!read (1) tavg_t % w
!!read (1) tavg_t % u2
!!read (1) tavg_t % v2
!!read (1) tavg_t % w2
!!read (1) tavg_t % uw
!!read (1) tavg_t % vw
!!read (1) tavg_t % uv
!!read (1) tavg_t % dudz

!!$if($LVLSET)
!!$if(RNS_LS)
!!read (1) tavg_t % fx
!!read (1) tavg_t % fy
!!read (1) tavg_t % fz
!!$endif
!!$endif

close(1)

! Now initialize all quantities for summation
tavg_t % u = tavg_t % u * tavg_total_time
tavg_t % v = tavg_t % v * tavg_total_time
tavg_t % w = tavg_t % w * tavg_total_time
tavg_t % u2 = tavg_t % u2 * tavg_total_time
tavg_t % v2 = tavg_t % v2 * tavg_total_time
tavg_t % w2 = tavg_t % w2 * tavg_total_time
tavg_t % uw = tavg_t % uw * tavg_total_time
tavg_t % vw = tavg_t % vw * tavg_total_time
tavg_t % uv = tavg_t % uv * tavg_total_time
tavg_t % dudz = tavg_t % dudz * tavg_total_time
tavg_t % dvdz = tavg_t % dvdz * tavg_total_time

tavg_t % txx = tavg_t % txx * tavg_total_time
tavg_t % txy = tavg_t % txy * tavg_total_time
tavg_t % tyy = tavg_t % tyy * tavg_total_time
tavg_t % txz = tavg_t % txz * tavg_total_time
tavg_t % tyz = tavg_t % tyz * tavg_total_time
tavg_t % tzz = tavg_t % tzz * tavg_total_time

tavg_t % fx = tavg_t % fx * tavg_total_time
tavg_t % fy = tavg_t % fy * tavg_total_time
tavg_t % fz = tavg_t % fz * tavg_total_time

tavg_zplane_t % u = tavg_zplane_t % u * tavg_total_time
tavg_zplane_t % v = tavg_zplane_t % v * tavg_total_time
tavg_zplane_t % w = tavg_zplane_t % w * tavg_total_time
tavg_zplane_t % u2 = tavg_zplane_t % u2 * tavg_total_time
tavg_zplane_t % v2 = tavg_zplane_t % v2 * tavg_total_time
tavg_zplane_t % w2 = tavg_zplane_t % w2 * tavg_total_time
tavg_zplane_t % uw = tavg_zplane_t % uw * tavg_total_time
tavg_zplane_t % vw = tavg_zplane_t % vw * tavg_total_time
tavg_zplane_t % uv = tavg_zplane_t % uv * tavg_total_time
tavg_zplane_t % dudz = tavg_zplane_t % dudz * tavg_total_time
tavg_zplane_t % dvdz = tavg_zplane_t % dvdz * tavg_total_time

tavg_zplane_t % txx = tavg_zplane_t % txx * tavg_total_time
tavg_zplane_t % txy = tavg_zplane_t % txy * tavg_total_time
tavg_zplane_t % tyy = tavg_zplane_t % tyy * tavg_total_time
tavg_zplane_t % txz = tavg_zplane_t % txz * tavg_total_time
tavg_zplane_t % tyz = tavg_zplane_t % tyz * tavg_total_time
tavg_zplane_t % tzz = tavg_zplane_t % tzz * tavg_total_time

tavg_zplane_t % fx = tavg_zplane_t % fx * tavg_total_time
tavg_zplane_t % fy = tavg_zplane_t % fy * tavg_total_time
tavg_zplane_t % fz = tavg_zplane_t % fz * tavg_total_time

tavg_total_time = tavg_total_time + dble(tavg_nend - tavg_nstart + 1) * dt

return
end subroutine tavg_init

!**********************************************************************
subroutine inst_write(itype)
!**********************************************************************
!  This subroutine writes the instantaneous values
!  at specified i,j,k locations
use functions, only : linear_interp, trilinear_interp
use stat_defs, only : point_t, domain_t, xplane_t, yplane_t, zplane_t
use grid_defs, only : x,y,z,zw
use sim_param, only : u,v,w

$if($LVLSET)
use level_set, only : phi
$if($RNS_LS)
use immersedbc, only : fx, fy, fz
$endif
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

    call write_real_data(point_t%fname(n), 'append', 4, (/ total_time, &
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
  write (fname,*) 'output/vel.', trim(adjustl(ct)),'.dat'
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
  
  
  call write_tecplot_header_ND(fname, 'rewind', nvars, (/ Nx+1, Ny+1, Nz/), var_list, coord, 2, total_time)
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
  call write_real_data_3D(fname, 'append', 'formatted', 4, nx, ny, nz, &
    (/ u(1:nx,1:ny,1:nz), v(1:nx,1:ny,1:nz), w_uv(1:nx,1:ny,1:nz), phi(1:nx,1:ny,1:nz) /), & 
    4, x, y, z(1:nz))
  $else
  call write_real_data_3D(fname, 'append', 'formatted', 3, nx,ny,nz, &
    (/ u(1:nx,1:ny,1:nz), v(1:nx,1:ny,1:nz), w_uv(1:nx,1:ny,1:nz) /), &
    4, x, y, z(1:nz))
  $endif
  
  !  Output Instantaneous Force Field for RNS Simulations
  !  Still need to put fz on uv grid may need a better way
  $if($LVLSET)
  $if($RNS_LS)
  !  Open file which to write global data
  write (fname,*) 'output/force.', trim(adjustl(ct)),'.dat'
  fname = trim(adjustl(fname))

  $if ($MPI)
    write (temp, '(".c",i0)') coord
    fname = trim (fname) // temp
  $endif
  
  
    
  !write(7,*) 'variables = "x", "y", "z", "u", "v", "w", "phi"';
  var_list = '"x", "y", "z", "f<sub>x</sub>", "f<sub>y</sub>", "f<sub>z</sub>", "phi"'
  nvars = 7
  call write_tecplot_header_ND(fname, 'rewind', nvars, (/ Nx+1, Ny+1, Nz/), var_list, coord, 2, total_time)
  call write_real_data_3D(fname, 'append', 'formatted', 4, nx, ny,nz, &
  (/ fx(1:nx,1:ny,1:nz), fy(1:nx,1:ny,1:nz), fz(1:nx,1:ny,1:nz), phi(1:nx,1:ny,1:nz) /), &
  4, x, y, z(1:nz))
  
  $endif
  $endif

!  Write instantaneous y-plane values
elseif(itype==3) then

!  Loop over all yplane locations
  do j=1,yplane_t%nloc

    write(cl,'(F9.4)') yplane_t%loc(j)
    !  Convert total iteration time to string
    write(ct,*) jt_total
    write(fname,*) 'output/vel.y-',trim(adjustl(cl)),'.',trim(adjustl(ct)),'.dat'
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
	  '"x", "y", "z", "u", "v", "w"', coord, 2, total_time)  
	  
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
  
  $if($LVLSET)
  $if($RNS_LS)
  
    write(fname,*) 'output/force.y-',trim(adjustl(cl)),'.',trim(adjustl(ct)),'.dat'
    fname=trim(adjustl(fname))

    $if ($MPI)
!  For MPI implementation
      write (temp, '(".c",i0)') coord
      fname = trim (fname) // temp
    $endif
  
    call write_tecplot_header_ND(fname, 'rewind', 6, (/ Nx, 1, Nz/), &
	  '"x", "y", "z", "fx", "fy", "fz"', coord, 2, total_time)  
	  
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

        ui = linear_interp(fx(i,yplane_t%istart(j),k), &
          fx(i,yplane_t%istart(j)+1,k), dy, yplane_t%ldiff(j))
        vi = linear_interp(fy(i,yplane_t%istart(j),k), &
          fy(i,yplane_t%istart(j)+1,k), dy, yplane_t%ldiff(j))
        wi = linear_interp(fz(i,yplane_t%istart(j),k), &
          fz(i,yplane_t%istart(j)+1,k), dy, &
          yplane_t%ldiff(j))

        write(2,*) x(i), yplane_t%loc(j), z(k), ui, vi, wi
      enddo
    enddo
    close(2)
    
    $endif
    $endif
    
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
    write(fname,*) 'output/vel.z-',trim(adjustl(cl)),'.',trim(adjustl(ct)),'.dat'
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
	 '"x", "y", "z", "u", "v", "w"', coord, 2, total_time)
	  	  
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
    
    $if($LVLSET)
    $if($RNS_LS)
    
    write(fname,*) 'output/force.z-',trim(adjustl(cl)),'.',trim(adjustl(ct)),'.dat'
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
	 '"x", "y", "z", "fx", "fy", "fz"', coord, 2, total_time)
	  	  
    open (unit = 2,file = fname, status='unknown',form='formatted', &
      action='write',position='append')	  

    do j=1,Ny
      do i=1,Nx

        ui = linear_interp(fx(i,j,zplane_t%istart(k)),fx(i,j,zplane_t%istart(k)+1), &
          dz, zplane_t%ldiff(k))
        vi = linear_interp(fy(i,j,zplane_t%istart(k)),fy(i,j,zplane_t%istart(k)+1), &
          dz, zplane_t%ldiff(k))
        wi = linear_interp(fz(i,j,zplane_t%istart(k)), &
          fz(i,j,zplane_t%istart(k)+1), &
          dz, zplane_t%ldiff(k))

        write(2,*) x(i), y(j), zplane_t%loc(k), ui, vi, wi

      enddo
    enddo
    close(2)
    
    $endif
    $endif

    $if ($MPI)
    endif
    $endif

  enddo  

!  Write instantaneous x-plane values
elseif(itype==5) then

!  Loop over all xplane locations
  do j=1,xplane_t%nloc

    write(cl,'(F9.4)') xplane_t%loc(j)
    !  Convert total iteration time to string
    write(ct,*) jt_total
    write(fname,*) 'output/vel.x-',trim(adjustl(cl)),'.',trim(adjustl(ct)),'.dat'
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
	
    call write_tecplot_header_ND(fname, 'rewind', 6, (/ 1, Ny, Nz/), &
	  '"x", "y", "z", "u", "v", "w"', coord, 2, total_time)  
	  
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
      do i=1,ny

        ui = linear_interp(u(xplane_t%istart(j),i,k), &
          u(xplane_t%istart(j)+1,i,k), dx, xplane_t%ldiff(j))
        vi = linear_interp(v(xplane_t%istart(j),i,k), &
          v(xplane_t%istart(j)+1,i,k), dx, xplane_t%ldiff(j))
        wi = linear_interp(w_uv(xplane_t%istart(j),i,k), &
          w_uv(xplane_t%istart(j)+1,i,k), dx, &
          xplane_t%ldiff(j))

        write(2,*) xplane_t%loc(j), y(i), z(k), ui, vi, wi
      enddo
    enddo
    close(2)
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
  imax, jmax, kmax, vars, ibuff, x,y,z)
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
!  ibuff (int) - flag for adding buffer region due to periodicity
!     0 - no buffer region
!     1 - buffer on i direction
!     2 - buffer on j direction
!     3 - buffer on k direction
!     4 - buffer on i,j direction
!     5 - buffer on i,k direction
!     6 - buffer on j,k direction
!     7 - buffer on i,j,k directions
!  x,y,z (real, vector, optional) - vectors containing x,y,z coordinates 
!
use functions, only : buff_indx
implicit none

character(*), intent(in) :: fname, write_pos, write_fmt
integer, intent(in) :: nvars, imax, jmax, kmax
real(rprec), intent(in), dimension(nvars*imax*jmax*kmax) :: vars
integer, intent(in) :: ibuff
real(rprec), intent(in), dimension(:), optional :: x,y,z

character(64) :: frmt
logical :: exst, coord_pres

character(*), parameter :: sub_name = mod_name // '.write_real_data_3D'

integer :: i,j,k,n
integer :: i0, j0, k0, imax_buff, jmax_buff, kmax_buff

!integer, allocatable, dimension(:) :: ikey_x, ikey_y, ikey_z
integer, allocatable, dimension(:,:,:,:) :: ikey_vars
!real(rprec), allocatable, dimension(:,:,:,:) :: vars_dim

!  Check if file exists
!inquire ( file=fname, exist=exst)
!if (.not. exst) call mesg(sub_name, 'Creating : ' // fname)
call check_write_pos(write_pos, sub_name)
call check_write_fmt(write_fmt, sub_name)

!  Check if spatial coordinates are specified
coord_pres=.false.
if(present(x) .and. present(y) .and. present(z)) coord_pres = .true.

if( ibuff == 0 ) then

  imax_buff = imax
  jmax_buff = jmax
  kmax_buff = kmax
  
elseif( ibuff == 1 ) then

  imax_buff = imax + 1
  jmax_buff = jmax
  kmax_buff = kmax

elseif( ibuff == 2 ) then

  imax_buff = imax
  jmax_buff = jmax + 1
  kmax_buff = kmax
  
elseif( ibuff == 3 ) then

  imax_buff = imax
  jmax_buff = jmax
  kmax_buff = kmax + 1
  
elseif( ibuff == 4 ) then

  imax_buff = imax + 1
  jmax_buff = jmax + 1
  kmax_buff = kmax  
  
elseif( ibuff == 5 ) then

  imax_buff = imax + 1
  jmax_buff = jmax
  kmax_buff = kmax + 1  
  
elseif( ibuff == 6 ) then

  imax_buff = imax
  jmax_buff = jmax + 1
  kmax_buff = kmax + 1  
  
elseif( ibuff == 7 ) then

  imax_buff = imax + 1
  jmax_buff = jmax + 1
  kmax_buff = kmax + 1  
  
else 

  call error(sub_name, 'ibuff not specified correctly')
  
endif

allocate(ikey_vars(nvars,imax_buff,jmax_buff,kmax_buff)) 


do n=1,nvars

  do k=1,kmax_buff
  
    k0 = buff_indx(k,kmax)
    
    do j = 1, jmax_buff
  
      j0 = buff_indx(j,jmax)
    
      do i = 1, imax_buff

        i0 = buff_indx(i,imax)
      
        ikey_vars(n,i,j,k) = (n-1)*imax*jmax*kmax + (k0-1)*imax*jmax + (j0-1)*imax + i0
        
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
    	
    if (coord_pres) then
	    
	    do k=1, kmax_buff
  	    do j=1, jmax_buff
    	    do i=1, imax_buff
            write(2,'(1e)') x(i)
	        enddo 
	      enddo
	    enddo
      
      do k=1, kmax_buff
  	    do j=1, jmax_buff
    	    do i=1, imax_buff
            write(2,'(1e)') y(j)
	        enddo 
	      enddo
	    enddo
      
	    do k=1, kmax_buff
  	    do j=1, jmax_buff
    	    do i=1, imax_buff
            write(2,'(1e)') z(k)
	        enddo 
	      enddo
	    enddo      
      
    endif      
    
    do n=1, nvars
	    do k=1, kmax_buff
        do j=1, jmax_buff
          do i=1, imax_buff
            write(2,'(1e)') vars(ikey_vars(n,i,j,k))
	        enddo 
	      enddo
	    enddo 
    enddo

  case('unformatted')
  
    if (coord_pres) then
	  
	    do k=1, kmax_buff
  	    do j=1, jmax_buff
    	    do i=1, imax_buff
            write(2) x(i)
	        enddo 
	      enddo
	    enddo
      
      do k=1, kmax_buff
  	    do j=1, jmax_buff
    	    do i=1, imax_buff
            write(2) y(j)
	        enddo 
	      enddo
	    enddo
      
	    do k=1, kmax_buff
  	    do j=1, jmax_buff
    	    do i=1, imax_buff
            write(2) z(k)
	        enddo 
	      enddo
	    enddo  
	  
    endif
    
    do n=1, nvars
	    do k=1, kmax_buff
        do j=1, jmax_buff
          do i=1, imax_buff
            write(2) vars(ikey_vars(n,i,j,k))
	        enddo 
	      enddo
	    enddo 
    enddo
	
end select

close(2)
  
deallocate(ikey_vars)

return
end subroutine write_real_data_3D

!*************************************************************
subroutine write_real_data_2D(fname, write_pos, write_fmt, nvars, &
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
!  write_pos (char) - postition to write in file : 'append' or 'rewind'
!  write_fmt (char) - Fotran format flag : 'formatted' or 'unformatted'
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
use functions, only : buff_indx
implicit none

character(*), intent(in) :: fname, write_pos, write_fmt
integer, intent(in) :: nvars, imax, jmax
real(rprec), intent(in), dimension(nvars*imax*jmax) :: vars
integer, intent(in) :: ibuff
real(rprec), intent(in), dimension(:), optional :: x, y

character(64) :: frmt
logical :: exst, coord_pres

character(*), parameter :: sub_name = mod_name // '.write_real_data_2D'

integer :: i,j,n
integer :: i0, j0, imax_buff, jmax_buff
!integer, allocatable, dimension(:) :: ikey_x, ikey_y
integer, allocatable, dimension(:,:,:) :: ikey_vars

!  Check if file exists
!inquire ( file=fname, exist=exst)
!if (.not. exst) call mesg(sub_name, 'Creating : ' // fname)
call check_write_pos(write_pos, sub_name)
call check_write_fmt(write_fmt, sub_name)

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

  call error(sub_name, 'ibuff not specified correctly')
  
endif

allocate(ikey_vars(nvars,imax_buff,jmax_buff)) 

!if(coord_pres) then

!  allocate(ikey_x(imax_buff), ikey_y(jmax_buff))
!  
!  do i=1, imax_buff
!  
!    i0 = buff_indx(i,imax)
!    
!    ikey_x(i) = i0
!    
!  enddo
!  
!  do j=1, jmax_buff

!    j0 = buff_indx(j,jmax)
!  
!    ikey_y(j) = j0
!    
!  enddo

!endif
  
do n=1,nvars

  do j = 1, jmax_buff
  
    j0 = buff_indx(j,jmax)
    
    do i = 1, imax_buff

      i0 = buff_indx(i,imax)
      
      ikey_vars(n,i,j) = (n-1)*imax*jmax + (j0-1)*imax + i0
    
    enddo
    
  enddo
  
enddo

open (unit = 2,file = fname, status='unknown',form=write_fmt, &
  action='write',position=write_pos)

!  Write the data
select case(write_fmt)
  case('formatted')
  
    !  Specify output format; may want to use a global setting
    	
    if (coord_pres) then

	    do j=1, jmax_buff
  	    do i=1,imax_buff
          write(2,'(1e)') x(i)
	      enddo 
	    enddo
    	  
      do j=1, jmax_buff
  	    do i=1,imax_buff
          write(2,'(1e)') y(j)
	      enddo 
	    enddo
 
    endif
	  
    do n=1, nvars
      do j=1, jmax_buff
	      do i=1,imax_buff
          write(2,'(1e)') vars(ikey_vars(n,i,j))
	      enddo 
	    enddo
    enddo

  case('unformatted')
  
    if (coord_pres) then
	  
	  
	    do j=1, jmax_buff
  	    do i=1,imax_buff
          write(2) x(i)
	      enddo 
	    enddo
    	  
      do j=1, jmax_buff
  	    do i=1,imax_buff
          write(2) y(j)
	      enddo 
	    enddo
 
    endif

    do n=1, nvars
      do j=1, jmax_buff
	      do i=1,imax_buff
          write(2) vars(ikey_vars(n,i,j))
	      enddo 
	    enddo
    enddo    
	
  !case default
  !  call error(sub_name, 'Incorrect write format : ' // write_pos)
end select

close(2)

deallocate(ikey_vars)

return
end subroutine write_real_data_2D

!*************************************************************
subroutine write_real_data_1D(fname, write_pos, write_fmt, nvars, &
  imax, vars, ibuff, x)
!*************************************************************
! 
!  This subroutine writes the variables given by vars to the
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
!  ibuff (int) - flag for adding buffer region due to periodicity
!     0 - no buffer region
!     1 - buffer on i direction (i = 1 = imax + 1)
!  x (real,vector,optional) - vector containing x coordinates 
!
use functions, only : buff_indx
implicit none

character(*), intent(in) :: fname, write_pos, write_fmt
integer, intent(in) :: nvars, imax
real(rprec), intent(in), dimension(nvars*imax) :: vars
integer, intent(in) :: ibuff
real(rprec), intent(in), dimension(:), optional :: x

character(64) :: frmt
logical :: exst, coord_pres

character(*), parameter :: sub_name = mod_name // '.write_real_data_1D'

integer :: i,n
integer :: i0, imax_buff
!integer, allocatable, dimension(:) :: ikey_x
integer, allocatable, dimension(:,:) :: ikey_vars

!  Check if file exists
!inquire ( file=fname, exist=exst)
!if (.not. exst) call mesg(sub_name, 'Creating : ' // fname)
call check_write_pos(write_pos, sub_name)
call check_write_fmt(write_fmt, sub_name)

!  Check if spatial coordinates specified
if(present(x)) coord_pres = .true.

if( ibuff == 0 ) then

  imax_buff = imax
  
elseif( ibuff == 1 ) then

  imax_buff = imax + 1
 
else 

  call error(sub_name, 'ibuff not specified correctly')
  
endif

allocate(ikey_vars(nvars,imax_buff)) 
  
do n=1,nvars
 
  do i = 1, imax_buff

    i0 = buff_indx(i,imax)
      
    ikey_vars(n,i) = (n-1)*imax + i0
    
  enddo
  
enddo

open (unit = 2,file = fname, status='unknown',form=write_fmt, &
  action='write',position=write_pos)
   
!  Write the data
select case(write_fmt)
  case('formatted')
  
    if (coord_pres) then
	  
      do i=1,imax_buff
        write(2,'(1e)') x(i)
      enddo 
	  
    endif
    
    do n=1, nvars
      do i=1,imax_buff
        write(2,'(1e)') vars(ikey_vars(n,i))
      enddo
	  enddo 

  case('unformatted')
  
    if (coord_pres) then
	  
      do i=1,imax_buff
        write(2) x(i)
      enddo 
	  
    endif
    
    do n=1, nvars
      do i=1,imax_buff
        write(2) vars(ikey_vars(n,i))
      enddo
	  enddo 

end select

close(2)

deallocate(ikey_vars)

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
    zone,'", DATAPACKING=BLOCK, i=', domain_size(1)
elseif(ndims == 2) then
  write(tec_dat_str,"(1a,i9,1a,i3,1a,i3,1a,i3)") 'ZONE T="', &
    zone,'", DATAPACKING=BLOCK, i=', domain_size(1),', j=', domain_size(2)
elseif(ndims == 3) then
  write(tec_dat_str,"(1a,i9,1a,i3,1a,i3,1a,i3,1a,i3)") 'ZONE T="', &
    zone,'", DATAPACKING=BLOCK, i=', domain_size(1),', j=', domain_size(2),', k=', domain_size(3)
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

!**********************************************************************
subroutine tavg_finalize()
!**********************************************************************
use grid_defs, only : x,y,z
use stat_defs, only : tavg_t, tavg_zplane_t, tavg, tavg_total_time
use stat_defs, only : rs, rs_t, rs_zplane_t
use param, only : nx,ny,nz,dx,dy,dz,L_x,L_y,L_z, nz_tot
$if($MPI)
use mpi
use mpi_defs, only : mpi_sync_real_array
use param, only : MPI_RPREC, rank_of_coord, comm, ierr
use stat_defs, only : rs_zplane_buf_t, rs_zplane_tot_t
use stat_defs, only : tavg_zplane_buf_t, tavg_zplane_tot_t
$endif
implicit none

character (*), parameter :: sub_name = mod_name // '.tavg_finalize'
character (64) :: temp
character(64) :: fname_out, fname_vel, fname_vel2, fname_ddz, &
  fname_tau, fname_f, fname_rs
character(64) :: fname_vel_zplane, fname_vel2_zplane, &
  fname_ddz_zplane, fname_tau_zplane, fname_f_zplane, fname_rs_zplane

integer :: i,j,k

$if($MPI)
real(rprec), allocatable, dimension(:) :: z_tot
integer :: MPI_RS, MPI_TSTATS
integer :: rs_type(1), rs_block(1), rs_disp(1)
integer :: tavg_type(1), tavg_block(1), tavg_disp(1)

integer :: ip, kbuf_start, kbuf_end, ktot_start, ktot_end
integer, allocatable, dimension(:) :: gather_coord
$endif

logical :: opn


!type rs
!real(rprec) :: up2, vp2, wp2, upwp, vpwp, upvp
!end type rs

!type(rs), pointer, dimension(:,:,:) :: rs_t
!type(rs), pointer, dimension(:) :: rs_zplane_t

!!$if($MPI)
!type(rs), pointer, dimension(:) :: rs_zplane_tot_t
!type(rs), pointer, dimension(:) :: rs_zplane_buf_t

!type(tavg), pointer, dimension(:) :: tavg_zplane_tot_t
!type(tavg), pointer, dimension(:) :: tavg_zplane_buf_t
!!$endif

allocate(rs_t(nx,ny,nz), rs_zplane_t(nz))

$if($MPI)

rs_type = MPI_RPREC
rs_block = 6 ! Number of rs subtypes
rs_disp = 0

tavg_type = MPI_RPREC
tavg_block = 20 ! Number of tavg subtypes
tavg_disp = 0

if(coord == 0) then

  !  Allocate space only on base processor for assembled z-plane data
  ! *_tot_t is the assembled data without the overlap nodes (the final stuff that is outputted)
  ! *_buf_t contains the overlap data and is used to recieve the z-plane data from all other nodes
  allocate(rs_zplane_tot_t(nz_tot), rs_zplane_buf_t(nz*nproc))
  allocate(tavg_zplane_tot_t(nz_tot), tavg_zplane_buf_t(nz*nproc))
  
  allocate(z_tot(nz_tot))
  ! In order to ensure that *_tot_t is assembled correctly we make sure that the processor number 
  ! is consistent with the spatial location
  allocate(gather_coord(nproc))
    
  do k=1, nz_tot
  
    z_tot(k) = (dble(k) - 0.5_rprec) * dz
    
  enddo
  
elseif(coord == nproc - 1) then

  !  Set fx,fy,fz to 0 as these are bogus
  tavg_t(:,:,nz) % fx = 0._rprec
  tavg_t(:,:,nz) % fy = 0._rprec
  tavg_t(:,:,nz) % fz = 0._rprec
  
  tavg_t(:,:,nz) % txx = 0._rprec
  tavg_t(:,:,nz) % txy = 0._rprec
  tavg_t(:,:,nz) % tyy = 0._rprec
  tavg_t(:,:,nz) % txz = 0._rprec
  tavg_t(:,:,nz) % tyz = 0._rprec
  tavg_t(:,:,nz) % tzz = 0._rprec
  
  tavg_zplane_t(nz) % fx = 0._rprec
  tavg_zplane_t(nz) % fy = 0._rprec
  tavg_zplane_t(nz) % fz = 0._rprec
  
  tavg_zplane_t(nz) % txx = 0._rprec
  tavg_zplane_t(nz) % txy = 0._rprec
  tavg_zplane_t(nz) % tyy = 0._rprec  
  tavg_zplane_t(nz) % txz = 0._rprec
  tavg_zplane_t(nz) % tyz = 0._rprec
  tavg_zplane_t(nz) % tzz = 0._rprec
  
endif

$else

  !  Set fx,fy,fz to 0 as these are bogus
  tavg_t(:,:,nz) % fx = 0._rprec
  tavg_t(:,:,nz) % fy = 0._rprec
  tavg_t(:,:,nz) % fz = 0._rprec
  
  tavg_t(:,:,nz) % txx = 0._rprec
  tavg_t(:,:,nz) % txy = 0._rprec
  tavg_t(:,:,nz) % tyy = 0._rprec
  tavg_t(:,:,nz) % txz = 0._rprec
  tavg_t(:,:,nz) % tyz = 0._rprec
  tavg_t(:,:,nz) % tzz = 0._rprec
  
  tavg_zplane_t(nz) % fx = 0._rprec
  tavg_zplane_t(nz) % fy = 0._rprec
  tavg_zplane_t(nz) % fz = 0._rprec
  
  tavg_zplane_t(nz) % txx = 0._rprec
  tavg_zplane_t(nz) % txy = 0._rprec
  tavg_zplane_t(nz) % tyy = 0._rprec  
  tavg_zplane_t(nz) % txz = 0._rprec
  tavg_zplane_t(nz) % tyz = 0._rprec
  tavg_zplane_t(nz) % tzz = 0._rprec

$endif
! All processors need not do this
!  Set file names
fname_out = 'tavg.out'

fname_vel = 'output/vel_avg.dat'
fname_vel2 = 'output/vel2_avg.dat'
fname_ddz = 'output/ddz_avg.dat'
fname_tau = 'output/tau_avg.dat'
fname_f = 'output/force_avg.dat'
fname_rs = 'output/rs.dat'

fname_vel_zplane = 'output/vel_zplane_avg.dat'
fname_vel2_zplane = 'output/vel2_zplane_avg.dat'
fname_ddz_zplane = 'output/ddz_zplane_avg.dat'
fname_tau_zplane = 'output/tau_zplane_avg.dat'
fname_f_zplane = 'output/force_zplane_avg.dat'
fname_rs_zplane = 'output/rs_zplane.dat'

$if ($MPI)
!  For MPI implementation     
  write (temp, '(".c",i0)') coord
  fname_out = trim (fname_out) // temp
  
  fname_vel = trim (fname_vel) // temp
  fname_vel2 = trim (fname_vel2) // temp
  fname_ddz = trim (fname_ddz) // temp
  fname_tau = trim (fname_tau) // temp
  fname_f = trim (fname_f) // temp
  fname_rs = trim (fname_rs) // temp
   
$endif

!  Perform Averaging operation
tavg_t % u = tavg_t % u / tavg_total_time
tavg_t % v = tavg_t % v / tavg_total_time
tavg_t % w = tavg_t % w / tavg_total_time
tavg_t % u2 = tavg_t % u2 / tavg_total_time
tavg_t % v2 = tavg_t % v2 / tavg_total_time
tavg_t % w2 = tavg_t % w2 / tavg_total_time
tavg_t % uw = tavg_t % uw / tavg_total_time
tavg_t % vw = tavg_t % vw / tavg_total_time
tavg_t % uv = tavg_t % uv / tavg_total_time
tavg_t % dudz = tavg_t % dudz / tavg_total_time
tavg_t % dvdz = tavg_t % dvdz / tavg_total_time

tavg_t % txx = tavg_t % txx / tavg_total_time
tavg_t % txy = tavg_t % txy / tavg_total_time
tavg_t % tyy = tavg_t % tyy / tavg_total_time
tavg_t % txz = tavg_t % txz / tavg_total_time
tavg_t % tyz = tavg_t % tyz / tavg_total_time
tavg_t % tzz = tavg_t % tzz / tavg_total_time

tavg_t % fx = tavg_t % fx / tavg_total_time
tavg_t % fy = tavg_t % fy / tavg_total_time
tavg_t % fz = tavg_t % fz / tavg_total_time

!  Average zplane values
tavg_zplane_t % u = tavg_zplane_t % u / tavg_total_time
tavg_zplane_t % v = tavg_zplane_t % v / tavg_total_time
tavg_zplane_t % w = tavg_zplane_t % w / tavg_total_time
tavg_zplane_t % u2 = tavg_zplane_t % u2 / tavg_total_time
tavg_zplane_t % v2 = tavg_zplane_t % v2 / tavg_total_time
tavg_zplane_t % w2 = tavg_zplane_t % w2 / tavg_total_time
tavg_zplane_t % uw = tavg_zplane_t % uw / tavg_total_time
tavg_zplane_t % vw = tavg_zplane_t % vw / tavg_total_time
tavg_zplane_t % uv = tavg_zplane_t % uv / tavg_total_time
tavg_zplane_t % dudz = tavg_zplane_t % dudz / tavg_total_time
tavg_zplane_t % dvdz = tavg_zplane_t % dvdz / tavg_total_time

tavg_zplane_t % txx = tavg_zplane_t % txx / tavg_total_time
tavg_zplane_t % txy = tavg_zplane_t % txy / tavg_total_time
tavg_zplane_t % tyy = tavg_zplane_t % tyy / tavg_total_time
tavg_zplane_t % txz = tavg_zplane_t % txz / tavg_total_time
tavg_zplane_t % tyz = tavg_zplane_t % tyz / tavg_total_time
tavg_zplane_t % tzz = tavg_zplane_t % tzz / tavg_total_time

tavg_zplane_t % fx = tavg_zplane_t % fx / tavg_total_time
tavg_zplane_t % fy = tavg_zplane_t % fy / tavg_total_time
tavg_zplane_t % fz = tavg_zplane_t % fz / tavg_total_time

!  Compute the Reynolds Stresses
do k = 1, nz

  rs_zplane_t(k) % up2 = tavg_zplane_t(k) % u2 - tavg_zplane_t(k) % u * tavg_zplane_t(k) % u
  rs_zplane_t(k) % vp2 = tavg_zplane_t(k) % v2 - tavg_zplane_t(k) % v * tavg_zplane_t(k) % v
  rs_zplane_t(k) % wp2 = tavg_zplane_t(k) % w2 - tavg_zplane_t(k) % w * tavg_zplane_t(k) % w
  rs_zplane_t(k) % upwp = tavg_zplane_t(k) % uw - tavg_zplane_t(k) % u * tavg_zplane_t(k) % w
  rs_zplane_t(k) % vpwp = tavg_zplane_t(k) % vw - tavg_zplane_t(k) % v * tavg_zplane_t(k) % w
  rs_zplane_t(k) % upvp = tavg_zplane_t(k) % uv - tavg_zplane_t(k) % u * tavg_zplane_t(k) % v
  
  do j = 1, ny

    do i = 1, nx
    
      rs_t(i,j,k) % up2 = tavg_t(i,j,k) % u2 - tavg_t(i,j,k) % u * tavg_t(i,j,k) % u
      rs_t(i,j,k) % vp2 = tavg_t(i,j,k) % v2 - tavg_t(i,j,k) % v * tavg_t(i,j,k) % v
      rs_t(i,j,k) % wp2 = tavg_t(i,j,k) % w2 - tavg_t(i,j,k) % w * tavg_t(i,j,k) % w
      rs_t(i,j,k) % upwp = tavg_t(i,j,k) % uw - tavg_t(i,j,k) % u * tavg_t(i,j,k) % w
      rs_t(i,j,k) % vpwp = tavg_t(i,j,k) % vw - tavg_t(i,j,k) % v * tavg_t(i,j,k) % w
      rs_t(i,j,k) % upvp = tavg_t(i,j,k) % uv - tavg_t(i,j,k) % u * tavg_t(i,j,k) % v
      
    enddo
    
  enddo
  
enddo

!  Write data to tavg.out
inquire (unit=1, opened=opn)
if (opn) call error (sub_name, 'unit 1 already open')

open (1, file=fname_out, action='write', position='rewind', form='unformatted')

! write the entire structures
write (1) tavg_total_time
write (1) tavg_t          
write (1) tavg_zplane_t

close(1)

! Construct zplane data 
$if($MPI)
 
  !  Create MPI type structures consistent with the derived types
  call MPI_TYPE_STRUCT(1, rs_block, rs_disp, rs_type, MPI_RS, ierr)
  Call MPI_Type_commit(MPI_RS,ierr)

  call MPI_TYPE_STRUCT(1, tavg_block, tavg_disp, tavg_type, MPI_TSTATS, ierr)
  Call MPI_Type_commit(MPI_TSTATS,ierr)

  call mpi_gather( rs_zplane_t, nz, MPI_RS, rs_zplane_buf_t, nz, &
    MPI_RS, rank_of_coord(0), comm, ierr)
    
  call mpi_gather( tavg_zplane_t, nz, MPI_TSTATS, tavg_zplane_buf_t, nz, &
    MPI_TSTATS, rank_of_coord(0), comm, ierr)   
    
  !  Get the rank that was used for mpi_gather (ensure that assembly of {rs,tavg}_zplane_tot_t is
  !  done in increasing coord
  call mpi_gather( coord, 1, MPI_INTEGER, gather_coord, 1, &
  MPI_INTEGER, rank_of_coord(0), comm, ierr)
  
  call MPI_Type_free (MPI_RS, ierr)
  call mpi_type_free (MPI_TSTATS, ierr)
  
  if(coord == 0) then
  
  !  Need to remove overlapping nodes
  !! Set initial block of data  
  !  rs_zplane_tot_t(1:nz) = rs_zplane_buf_t(1:nz)
  !  tavg_zplane_tot_t(1:nz) = tavg_zplane_buf_t(1:nz)
    
    do ip=1, nproc
      
      kbuf_start = gather_coord(ip)*nz + 1
      kbuf_end   = kbuf_start + nz - 1
      
      ktot_start = kbuf_start - gather_coord(ip)
      ktot_end = ktot_start + nz - 1
                
      rs_zplane_tot_t(ktot_start:ktot_end) = rs_zplane_buf_t(kbuf_start:kbuf_end)
      tavg_zplane_tot_t(ktot_start:ktot_end) = tavg_zplane_buf_t(kbuf_start:kbuf_end)
      
    enddo   
    
    call write_tecplot_header_ND(fname_vel_zplane, 'rewind', 4, (/ Nz_tot /), &
      '"z", "<u>","<v>","<w>"', coord, 2)
    call write_real_data_1D(fname_vel_zplane, 'append', 'formatted', 3, Nz_tot, &
      (/ tavg_zplane_tot_t % u, tavg_zplane_tot_t % v, tavg_zplane_tot_t % w /), 0, z_tot)

    call write_tecplot_header_ND(fname_vel2_zplane, 'rewind', 7, (/ Nz_tot/), &
      '"z", "<u<sup>2</sup>>","<v<sup>2</sup>>","<w<sup>2</sup>>", "<uw>", "<vw>", "<uv>"', coord, 2)
    call write_real_data_1D(fname_vel2_zplane, 'append', 'formatted', 6, Nz_tot, &
      (/ tavg_zplane_tot_t % u2, tavg_zplane_tot_t % v2, tavg_zplane_tot_t % w2, &
      tavg_zplane_tot_t % uw, tavg_zplane_tot_t % vw, tavg_zplane_tot_t % uv /), &
      0, z_tot) 
  
    call write_tecplot_header_ND(fname_ddz_zplane, 'rewind', 3, (/ Nz_tot/), &
      '"z", "<dudz>","<dvdz>"', coord, 2)
    call write_real_data_1D(fname_ddz_zplane, 'append', 'formatted', 2, Nz_tot, &
      (/ tavg_zplane_tot_t % dudz, tavg_zplane_tot_t % dvdz /), 0, z_tot)
  
    call write_tecplot_header_ND(fname_tau_zplane, 'rewind', 7, (/Nz_tot/), &
      '"z", "<t<sub>xx</sub>>","<t<sub>xy</sub>>","<t<sub>yy</sub>>", "<t<sub>xz</sub>>", "<t<sub>yx</sub>>", "<t<sub>zz</sub>>"', coord, 2)  
    call write_real_data_1D(fname_tau_zplane, 'append', 'formatted', 6, Nz_tot, &
      (/ tavg_zplane_tot_t % txx, tavg_zplane_tot_t % txy, tavg_zplane_tot_t % tyy, &
      tavg_zplane_tot_t % txz, tavg_zplane_tot_t % tyz, tavg_zplane_tot_t % tzz /), 0, z_tot) 
  
    call write_tecplot_header_ND(fname_f_zplane, 'rewind', 4, (/Nz_tot/), &
      '"z", "<fx>","<fy>","<fz>"', coord, 2)
    call write_real_data_1D(fname_f_zplane, 'append', 'formatted', 3, Nz_tot, &
      (/ tavg_zplane_tot_t % fx, tavg_zplane_tot_t % fy, tavg_zplane_tot_t % fz /), 0, z_tot)  
  
    call write_tecplot_header_ND(fname_rs_zplane, 'rewind', 7, (/Nz_tot/), &
      '"z", "<upup>","<vpvp>","<wpwp>", "<upwp>", "<vpwp>", "<upvp>"', coord, 2)  
    call write_real_data_1D(fname_rs_zplane, 'append', 'formatted', 6, Nz_tot, &
      (/ rs_zplane_tot_t % up2, rs_zplane_tot_t%vp2, rs_zplane_tot_t%wp2, &
      rs_zplane_tot_t%upwp, rs_zplane_tot_t%vpwp, rs_zplane_tot_t%upvp /), 0, z_tot)    
 
    deallocate(tavg_zplane_tot_t, tavg_zplane_buf_t)
    deallocate(rs_zplane_tot_t, rs_zplane_buf_t)
    deallocate(z_tot, gather_coord)
  
  endif

$else

call write_tecplot_header_ND(fname_vel_zplane, 'rewind', 4, (/ Nz /), &
   '"x", "y", "z", "<u>","<v>","<w>"', coord, 2)
call write_real_data_1D(fname_vel_zplane, 'append', 'formatted', 3, nz, &
  (/ tavg_zplane_t % u, tavg_zplane_t % v, tavg_zplane_t % w /), 0, z(1:nz))

call write_tecplot_header_ND(fname_vel2_zplane, 'rewind', 7, (/ Nz/), &
   '"z", "<u<sup>2</sup>>","<v<sup>2</sup>>","<w<sup>2</sup>>", "<uw>", "<vw>", "<uv>"', coord, 2)
call write_real_data_1D(fname_vel2_zplane, 'append', 'formatted', 6, nz, &
  (/ tavg_zplane_t % u2, tavg_zplane_t % v2, tavg_zplane_t % w2, &
  tavg_zplane_t % uw, tavg_zplane_t % vw, tavg_zplane_t % uv /), 0, z(1:nz)) 
  
call write_tecplot_header_ND(fname_ddz_zplane, 'rewind', 3, (/ Nz/), &
   '"z", "<dudz>","<dvdz>"', coord, 2)
call write_real_data_1D(fname_ddz_zplane, 'append', 'formatted', 2, nz, &
  (/ tavg_zplane_t % dudz, tavg_zplane_t % dvdz /), 0, z(1:nz))
  
call write_tecplot_header_ND(fname_tau_zplane, 'rewind', 7, (/Nz/), &
   '"z", "<t<sub>xx</sub>>","<t<sub>xy</sub>>","<t<sub>yy</sub>>", "<t<sub>xz</sub>>", "<t<sub>yx</sub>>", "<t<sub>zz</sub>>"', coord, 2)  
call write_real_data_1D(fname_tau_zplane, 'append', 'formatted', 6, nz, &
  (/ tavg_zplane_t % txx, tavg_zplane_t % txy, tavg_zplane_t % tyy, &
  tavg_zplane_t % txz, tavg_zplane_t % tyz, tavg_zplane_t % tzz /), 0, z(1:nz)) 
  
call write_tecplot_header_ND(fname_f_zplane, 'rewind', 4, (/Nz/), &
   '"z", "<fx>","<fy>","<fz>"', coord, 2)
call write_real_data_1D(fname_f_zplane, 'append', 'formatted', 3, nz, &
  (/ tavg_zplane_t % fx, tavg_zplane_t % fy, tavg_zplane_t % fz /), 0, z(1:nz))  
  
call write_tecplot_header_ND(fname_rs_zplane, 'rewind', 7, (/Nz/), &
   '"z", "<upup>","<vpvp>","<wpwp>", "<upwp>", "<vpwp>", "<upvp>"', coord, 2)  
call write_real_data_1D(fname_rs_zplane, 'append', 'formatted', 6, nz, &
  (/ rs_zplane_t % up2, rs_zplane_t%vp2, rs_zplane_t%wp2, &
  rs_zplane_t%upwp, rs_zplane_t%vpwp, rs_zplane_t%upvp /), 0, z(1:nz))

$endif

$if ($MPI)
!  Sync data across all nodes for a subset of varibles which contain bogus initialization
!call mpi_sync_real_array(  tavg_t % txx )
!call mpi_sync_real_array(  tavg_t % txy )
!call mpi_sync_real_array(  tavg_t % tyy )
!call mpi_sync_real_array(  tavg_t % txz )
!call mpi_sync_real_array(  tavg_t % tyz )
!call mpi_sync_real_array(  tavg_t % tzz )
$endif

! ----- Write all the 3D data -----

call write_tecplot_header_ND(fname_vel, 'rewind', 6, (/ Nx+1, Ny+1, Nz/), &
   '"x", "y", "z", "<u>","<v>","<w>"', coord, 2)
call write_real_data_3D(fname_vel, 'append', 'formatted', 3, nx, ny, nz, &
  (/ tavg_t % u, tavg_t % v, tavg_t % w /), &
  4, x, y, z(1:nz))

call write_tecplot_header_ND(fname_vel2, 'rewind', 9, (/ Nx+1, Ny+1, Nz/), &
   '"x", "y", "z", "<u<sup>2</sup>>","<v<sup>2</sup>>","<w<sup>2</sup>>", "<uw>", "<vw>", "<uv>"', coord, 2)
call write_real_data_3D(fname_vel2, 'append', 'formatted', 6, nx, ny, nz, &
  (/ tavg_t % u2, tavg_t % v2, tavg_t % w2, tavg_t % uw, tavg_t % vw, tavg_t % uv /), &
  4, x, y, z(1:nz)) 
  
call write_tecplot_header_ND(fname_ddz, 'rewind', 5, (/ Nx+1, Ny+1, Nz/), &
   '"x", "y", "z", "<dudz>","<dvdz>"', coord, 2)
call write_real_data_3D(fname_ddz, 'append', 'formatted', 2, nx, ny, nz, &
  (/ tavg_t % dudz, tavg_t % dvdz /), 4, x, y, z(1:nz))

call write_tecplot_header_ND(fname_tau, 'rewind', 9, (/ Nx+1, Ny+1, Nz/), &
   '"x", "y", "z", "<t<sub>xx</sub>>","<t<sub>xy</sub>>","<t<sub>yy</sub>>", "<t<sub>xz</sub>>", "<t<sub>yx</sub>>", "<t<sub>zz</sub>>"', coord, 2)  
call write_real_data_3D(fname_tau, 'append', 'formatted', 6, nx, ny, nz, &
  (/ tavg_t % txx, tavg_t % txy, tavg_t % tyy, tavg_t % txz, tavg_t % tyz, tavg_t % tzz /), &
  4, x, y, z(1:nz)) 
  
call write_tecplot_header_ND(fname_f, 'rewind', 6, (/ Nx+1, Ny+1, Nz/), &
   '"x", "y", "z", "<fx>","<fy>","<fz>"', coord, 2)
call write_real_data_3D(fname_f, 'append', 'formatted', 3, nx, ny, nz, &
  (/ tavg_t % fx, tavg_t % fy, tavg_t % fz /), &
  4, x, y, z(1:nz))  
  
call write_tecplot_header_ND(fname_rs, 'rewind', 9, (/ Nx+1, Ny+1, Nz/), &
   '"x", "y", "z", "<upup>","<vpvp>","<wpwp>", "<upwp>", "<vpwp>", "<upvp>"', coord, 2)  
call write_real_data_3D(fname_rs, 'append', 'formatted', 6, nx, ny, nz, &
  (/ rs_t % up2, rs_t%vp2, rs_t%wp2, rs_t%upwp, rs_t%vpwp, rs_t%upvp /), &
  4, x, y, z(1:nz))   

deallocate(tavg_t, tavg_zplane_t, rs_t, rs_zplane_t)

return
end subroutine tavg_finalize

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
use stat_defs, only : tavg_t, point_t, yplane_t, zplane_t
use param, only : tavg_calc

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
if(tavg_calc) call tavg_finalize()
 
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
use param, only : tavg_calc
use stat_defs
use grid_defs
use functions, only : cell_indx
implicit none

!character(120) :: cx,cy,cz
character(120) :: fname, var_list
integer :: fid, i,j,k

! --- ALL TIME AVERAGING CONTROLLED IN PARAM ---
!!  All nstart and nend values are based
!!  on jt and not jt_total
!tavg_t%calc = .true.
!tavg_nstart = 1
!tavg_nend = nsteps

!  Turns instantaneous velocity recording on or off
point_t%calc = .true.
point_t%nstart = 1
point_t%nend   = nsteps
point_t%nskip = 10
point_t%nloc = 1
point_t%xyz(:,1) = (/L_x/2. + 2.162162_rprec, L_y/2., 4.461996_rprec/)
!point_t%xyz(:,2) = (/L_x/2., L_y/2., 2.5_rprec/)

domain_t%calc = .true.
domain_t%nstart = 100
domain_t%nend   = nsteps
domain_t%nskip = 100

!  x-plane stats/data
xplane_t%calc   = .true.
xplane_t%nstart = 100
xplane_t%nend   = nsteps
xplane_t%nskip  = 100
xplane_t%nloc   = 2
xplane_t%loc(1) = 1.0
xplane_t%loc(2) = 3.0

!  y-plane stats/data
yplane_t%calc   = .true.
yplane_t%nstart = 100
yplane_t%nend   = nsteps
yplane_t%nskip  = 100
yplane_t%nloc   = 2
yplane_t%loc(1) = 1.0
yplane_t%loc(2) = 3.0

!  z-plane stats/data
zplane_t%calc   = .true.
zplane_t%nstart = 100
zplane_t%nend   = nsteps
zplane_t%nskip  = 100
zplane_t%nloc   = 7
zplane_t%loc(1) = 0.733347
zplane_t%loc(2) = 1.550644
zplane_t%loc(3) = 1.959293
zplane_t%loc(4) = 2.163617
zplane_t%loc(5) = 2.265780
zplane_t%loc(6) = 2.316861
zplane_t%loc(7) = 2.342401

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
if(tavg_calc) then 

  !allocate(tavg_t%u(nx, ny, nz))
  !allocate(tavg_t%v(nx, ny, nz))
  !allocate(tavg_t%w(nx, ny, nz))
  !allocate(tavg_t%u2(nx, ny, nz))
  !allocate(tavg_t%v2(nx, ny, nz))
  !allocate(tavg_t%w2(nx, ny, nz))
  !allocate(tavg_t%uw(nx, ny, nz))
  !allocate(tavg_t%vw(nx, ny, nz))
  !allocate(tavg_t%uv(nx, ny, nz))
  !allocate(tavg_t%dudz(nx, ny, nz))
  allocate(tavg_t(nx,ny,nz))
  allocate(tavg_zplane_t(nz))
  
  !  Initialize arrays
  tavg_t(:,:,:) = tavg(0._rprec, 0._rprec, 0._rprec, 0._rprec, &
    0._rprec, 0._rprec, 0._rprec, 0._rprec, &
    0._rprec, 0._rprec, 0._rprec, 0._rprec, &
    0._rprec, 0._rprec, 0._rprec, 0._rprec, &
    0._rprec, 0._rprec, 0._rprec, 0._rprec)
    
  tavg_zplane_t(:) = tavg(0._rprec, 0._rprec, 0._rprec, 0._rprec, &
    0._rprec, 0._rprec, 0._rprec, 0._rprec, &
    0._rprec, 0._rprec, 0._rprec, 0._rprec, &
    0._rprec, 0._rprec, 0._rprec, 0._rprec, &
    0._rprec, 0._rprec, 0._rprec, 0._rprec)
  
  !!$if($LVLSET)
  !!$if($RNS_LS)
  
  !!allocate(tavg_t%fx(nx,ny,nz))
  !!allocate(tavg_t%fy(nx,ny,nz))
  !!allocate(tavg_t%fz(nx,ny,nz))
  
  !!tavg_t%fx = 0._rprec
  !!tavg_t%fx = 0._rprec
  !!tavg_t%fx = 0._rprec
  
  !!$endif
  !!$endif
  
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

! Initialize information for x-planar stats/data
if(xplane_t%calc) then
  xplane_t%istart = -1
  xplane_t%ldiff = 0.
!  Not really needed
  xplane_t%coord = -1
  
!  Compute istart and ldiff
  do j=1,xplane_t%nloc
	  xplane_t%istart(j) = cell_indx('j', dx, xplane_t%loc(j))
	  xplane_t%ldiff(j) = x(xplane_t%istart(j)) - xplane_t%loc(j)
  enddo
    
endif

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
	  yplane_t%istart(j) = cell_indx('j', dy, yplane_t%loc(j))
	  yplane_t%ldiff(j) = y(yplane_t%istart(j)) - yplane_t%loc(j)
  enddo
    
endif

! Initialize information for z-planar stats/data
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

	  zplane_t%istart(k) = cell_indx('k',dz,zplane_t%loc(k))
	  zplane_t%ldiff(k) = z(zplane_t%istart(k)) - zplane_t%loc(k)
    endif
    $else
    zplane_t%coord(k) = 0
    zplane_t%istart(k) = cell_indx('k',dz,zplane_t%loc(k))
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
	  
	  point_t%istart(i) = cell_indx('i',dx,point_t%xyz(1,i))
	  point_t%jstart(i) = cell_indx('j',dy,point_t%xyz(2,i))
	  point_t%kstart(i) = cell_indx('k',dz,point_t%xyz(3,i))
	  
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
	point_t%istart(i) = cell_indx('i',dx,point_t%xyz(1,i))
	point_t%jstart(i) = cell_indx('j',dy,point_t%xyz(2,i))
	point_t%kstart(i) = cell_indx('k',dz,point_t%xyz(3,i))
	
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

!! Initialize for zplane_avg calculations
!if(zplane_avg_t%calc) then 
!  allocate(zplane_avg_t%u_avg(nz-1))
!  allocate(zplane_avg_t%v_avg(nz-1))
!  allocate(zplane_avg_t%w_avg(nz-1))
!  allocate(zplane_avg_t%u2_avg(nz-1))
!  allocate(zplane_avg_t%v2_avg(nz-1))
!  allocate(zplane_avg_t%w2_avg(nz-1))  
!  allocate(zplane_avg_t%txx_avg(nz-1))
!  allocate(zplane_avg_t%txy_avg(nz-1))
!  allocate(zplane_avg_t%tyy_avg(nz-1))
!  allocate(zplane_avg_t%txz_avg(nz-1))
!  allocate(zplane_avg_t%tyz_avg(nz-1))
!  allocate(zplane_avg_t%tzz_avg(nz-1))  
!  allocate(zplane_avg_t%dudz_avg(nz-1))
!  allocate(zplane_avg_t%dvdz_avg(nz-1))
!  allocate(zplane_avg_t%uv_avg(nz-1))
!  allocate(zplane_avg_t%uw_avg(nz-1))
!  allocate(zplane_avg_t%vw_avg(nz-1))  

!  !  Initialize arrays
!  zplane_avg_t%u_avg(:)=0.
!  zplane_avg_t%v_avg(:)=0.
!  zplane_avg_t%w_avg(:)=0.
!  zplane_avg_t%u2_avg(:)=0.
!  zplane_avg_t%v2_avg(:)=0.
!  zplane_avg_t%w2_avg(:)=0.  
!  zplane_avg_t%txx_avg(:)=0.
!  zplane_avg_t%txy_avg(:)=0.
!  zplane_avg_t%tyy_avg(:)=0.
!  zplane_avg_t%txz_avg(:)=0.
!  zplane_avg_t%tyz_avg(:)=0.
!  zplane_avg_t%tzz_avg(:)=0.  
!  zplane_avg_t%dudz_avg(:)=0.
!  zplane_avg_t%dvdz_avg(:)=0.  
!  zplane_avg_t%uv_avg(:)=0.
!  zplane_avg_t%uw_avg(:)=0.
!  zplane_avg_t%vw_avg(:)=0.    
!endif

return
end subroutine stats_init

!***************************************************************
subroutine tavg_compute()

!***************************************************************
!  This subroutine collects the stats for each flow 
!  variable quantity
use types, only : rprec
use stat_defs, only : tavg_t, tavg_zplane_t
use param, only : nx,ny,nz, dt
use sim_param, only : u,v,w, dudz, dvdz, txx, txy, tyy, txz, tyz, tzz
use immersedbc, only : fx, fy, fz

implicit none

!use io, only : w_uv, w_uv_tag, dudz_uv, dudz_uv_tag, interp_to_uv_grid
integer :: i,j,k
real(rprec) :: u_p, v_p, w_p, dudz_p, dvdz_p, fa

!  Make sure w stuff has been interpolated to uv-grid
call interp_to_uv_grid(w, w_uv, w_uv_tag)
call interp_to_uv_grid(dudz, dudz_uv, dudz_uv_tag)
call interp_to_uv_grid(dvdz, dvdz_uv, dvdz_uv_tag)

do k=1,nz

  !  Include dt for current weighting
  fa = dt / (nx * ny)
  
  tavg_zplane_t(k)%u = tavg_zplane_t(k)%u + fa * sum(u(1:nx,1:ny,k))
  tavg_zplane_t(k)%v = tavg_zplane_t(k)%v + fa * sum(v(1:nx,1:ny,k))
  tavg_zplane_t(k)%w = tavg_zplane_t(k)%w + fa * sum(w(1:nx,1:ny,k))
  tavg_zplane_t(k)%u2 = tavg_zplane_t(k)%u2 + fa * sum(u(1:nx,1:ny,k)*u(1:nx,1:ny,k))
  tavg_zplane_t(k)%v2 = tavg_zplane_t(k)%v2 + fa * sum(v(1:nx,1:ny,k)*v(1:nx,1:ny,k)) 
  tavg_zplane_t(k)%w2 = tavg_zplane_t(k)%w2 + fa * sum(w_uv(1:nx,1:ny,k)*w_uv(1:nx,1:ny,k)) 
  tavg_zplane_t(k)%uw = tavg_zplane_t(k)%uw+ fa * sum(u(1:nx,1:ny,k)*w_uv(1:nx,1:ny,k))
  tavg_zplane_t(k)%vw = tavg_zplane_t(k)%vw + fa * sum(v(1:nx,1:ny,k)*w_uv(1:nx,1:ny,k)) 
  tavg_zplane_t(k)%uv = tavg_zplane_t(k)%uv + fa * sum(u(1:nx,1:ny,k)*v(1:nx,1:ny,k))
  tavg_zplane_t(k)%dudz = tavg_zplane_t(k)%dudz + fa * sum(dudz_uv(1:nx,1:ny,k)) 
  tavg_zplane_t(k)%dvdz = tavg_zplane_t(k)%dvdz + fa * sum(dvdz_uv(1:nx,1:ny,k))
      
  tavg_zplane_t(k)%txx = tavg_zplane_t(k)%txx + fa * sum(txx(1:nx,1:ny,k))
  tavg_zplane_t(k)%txy = tavg_zplane_t(k)%txy + fa * sum(txy(1:nx,1:ny,k))
  tavg_zplane_t(k)%tyy = tavg_zplane_t(k)%tyy + fa * sum(tyy(1:nx,1:ny,k))
  tavg_zplane_t(k)%txz = tavg_zplane_t(k)%txz + fa * sum(txz(1:nx,1:ny,k))
  tavg_zplane_t(k)%tyz = tavg_zplane_t(k)%tyz + fa * sum(tyz(1:nx,1:ny,k))
  tavg_zplane_t(k)%tzz = tavg_zplane_t(k)%tzz + fa * sum(tzz(1:nx,1:ny,k))
      
  tavg_zplane_t(k)%fx = tavg_zplane_t(k)%fx + fa * sum(fx(1:nx,1:ny,k)) 
  tavg_zplane_t(k)%fy = tavg_zplane_t(k)%fy + fa * sum(fy(1:nx,1:ny,k))
  tavg_zplane_t(k)%fz = tavg_zplane_t(k)%fz + fa * sum(fz(1:nx,1:ny,k)) 
  
  do j=1,ny
    do i=1,nx
      
      !  Being cache friendly
      u_p = u(i,j,k)
      v_p = v(i,j,k) 
      w_p = w_uv(i,j,k)
      dudz_p = dudz_uv(i,j,k) 
      dvdz_p = dvdz_uv(i,j,k) 
          
      tavg_t(i,j,k)%u = tavg_t(i,j,k)%u + u_p * dt
      tavg_t(i,j,k)%v = tavg_t(i,j,k)%v + v_p * dt
      tavg_t(i,j,k)%w = tavg_t(i,j,k)%w + w_p * dt
      tavg_t(i,j,k)%u2 = tavg_t(i,j,k)%u2 + u_p*u_p * dt
      tavg_t(i,j,k)%v2 = tavg_t(i,j,k)%v2 + v_p*v_p * dt
      tavg_t(i,j,k)%w2 = tavg_t(i,j,k)%w2 + w_p*w_p * dt
      tavg_t(i,j,k)%uw = tavg_t(i,j,k)%uw+ u_p*w_p * dt
      tavg_t(i,j,k)%vw = tavg_t(i,j,k)%vw + v_p*w_p * dt
      tavg_t(i,j,k)%uv = tavg_t(i,j,k)%uv + u_p*v_p * dt
      tavg_t(i,j,k)%dudz = tavg_t(i,j,k)%dudz + dudz_p * dt
      tavg_t(i,j,k)%dvdz = tavg_t(i,j,k)%dvdz + dvdz_p * dt
      
      tavg_t(i,j,k)%txx = tavg_t(i,j,k)%txx + txx(i,j,k) * dt
      tavg_t(i,j,k)%txy = tavg_t(i,j,k)%txy + txy(i,j,k) * dt
      tavg_t(i,j,k)%tyy = tavg_t(i,j,k)%tyy + tyy(i,j,k) * dt
      tavg_t(i,j,k)%txz = tavg_t(i,j,k)%txz + txz(i,j,k) * dt
      tavg_t(i,j,k)%tyz = tavg_t(i,j,k)%tyz + tyz(i,j,k) * dt
      tavg_t(i,j,k)%tzz = tavg_t(i,j,k)%tzz + tzz(i,j,k) * dt
      
      tavg_t(i,j,k)%fx = tavg_t(i,j,k)%fx + fx(i,j,k) * dt 
      tavg_t(i,j,k)%fy = tavg_t(i,j,k)%fy + fy(i,j,k) * dt 
      tavg_t(i,j,k)%fz = tavg_t(i,j,k)%fz + fz(i,j,k) * dt 
 
      
    enddo
  enddo
enddo

return

end subroutine tavg_compute

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


!!***************************************************************
!subroutine zplane_avg_compute(ipart)
!!***************************************************************
!!  This subroutine averages the following parameters over z-planes:
!!  u, v, w, u^2, v^2, w^2, txx, txy, tyy, txz, tyz, tzz, dudz, dvdz
!!  and writes them to files avg_*.dat

!! ipart=1 to create files and write Tecplot headers
!! ipart=2 to compute averages
!! ipart=3 to write averages and reset


!use param, only : dx, dy, dz, jt_total, dt_dim
!use stat_defs, only: zplane_avg_t
!use sim_param, only : u, v, w, txx, txy, tyy, txz, tyz, tzz, dudz, dvdz
!use grid_defs, only: z

!implicit none

!integer :: k, ipart, nt_z
!character(64) :: frmt
!character (64) :: fname, temp, fname2, temp2, fname3, temp3, fname4, temp4, fname5, temp5, fname6, temp6

!if (ipart==1) then
!!Create files and write Tecplot headers         
!		!parameter file
!		  if(coord==0) then 
!		    fname='output/avg.dat'
!			open (unit = 2,file = fname, status='unknown',form='formatted', action='write',position='rewind')
!			write (frmt, '("(",i0,"i9)")') 3
!			write(2,*) 'nz, nstart, nend'
!			write(2,frmt) nz,zplane_avg_t%nstart, zplane_avg_t%nend
!			close(2)
!			!call write_real_data(fname, 'append', 3, (/ zplane_avg_t%nstart, zplane_avg_t%nend, zplane_avg_t%nskip /))
!		  endif
!		!u-v-w file  
!		  fname2 = 'output/avg_uvw.dat'
!			$if ($MPI)
!			  write (temp2, '(".c",i0)') coord
!			  fname2 = trim (fname2) // temp2
!			$endif	 	  
!		  call write_tecplot_header_ND(fname2, 'rewind', 4, (/ Nz/),'"z", "u_avg", "v_avg", "w_avg"', coord, 2, total_time_dim)  
!		!txx-txy,etc file
!		  fname3 = 'output/avg_tau.dat'
!			$if ($MPI)
!			  write (temp3, '(".c",i0)') coord
!			  fname3 = trim (fname3) // temp3
!			$endif	 	  
!		  call write_tecplot_header_ND(fname3, 'rewind', 7, (/ Nz/), &
!		  '"z", "txx_avg", "txy_avg", "tyy_avg", "txz_avg", "tyz_avg", "tzz_avg"', coord, 2, total_time_dim)  
!		 !ddz file  
!		  fname4 = 'output/avg_ddz.dat'
!			$if ($MPI)
!			  write (temp4, '(".c",i0)') coord
!			  fname4 = trim (fname4) // temp4
!			$endif	 	  
!		  call write_tecplot_header_ND(fname4, 'rewind', 3, (/ Nz/),'"z", "dudz_avg", "dvdz_avg"', coord, 2, total_time_dim) 
!		!u2-v2-w2 file  
!		  fname5 = 'output/avg_uvw2.dat'
!			$if ($MPI)
!			  write (temp5, '(".c",i0)') coord
!			  fname5 = trim (fname5) // temp5
!			$endif	 	  
!		  call write_tecplot_header_ND(fname5, 'rewind', 4, (/ Nz/),'"z", "u2_avg", "v2_avg", "w2_avg"', coord, 2, total_time_dim) 		  
!		!uv-uw-vw file  
!		  fname6 = 'output/avg_uv_uw_vw.dat'
!			$if ($MPI)
!			  write (temp6, '(".c",i0)') coord
!			  fname6 = trim (fname6) // temp6
!			$endif	 	  
!		  call write_tecplot_header_ND(fname6, 'rewind', 4, (/ Nz/),'"z", "uv_avg", "uw_avg", "vw_avg"', coord, 2, total_time_dim) 		  

!elseif (ipart==2) then
!!Compute averages
!	nt_z = zplane_avg_t%nend - zplane_avg_t%nstart + 1
!    do k=1,nz-1
!	    zplane_avg_t%u_avg(k) = zplane_avg_t%u_avg(k) + sum(u(:,:,k))/(nx*ny*nt_z)
!	    zplane_avg_t%v_avg(k) = zplane_avg_t%v_avg(k) + sum(v(:,:,k))/(nx*ny*nt_z)
!	    zplane_avg_t%w_avg(k) = zplane_avg_t%w_avg(k) + sum(w(:,:,k))/(nx*ny*nt_z)
!		
!		zplane_avg_t%txx_avg(k) = zplane_avg_t%txx_avg(k) + sum(txx(:,:,k))/(nx*ny*nt_z)
!		zplane_avg_t%txy_avg(k) = zplane_avg_t%txy_avg(k) + sum(txy(:,:,k))/(nx*ny*nt_z)
!		zplane_avg_t%tyy_avg(k) = zplane_avg_t%tyy_avg(k) + sum(tyy(:,:,k))/(nx*ny*nt_z)
!		zplane_avg_t%txz_avg(k) = zplane_avg_t%txz_avg(k) + sum(txz(:,:,k))/(nx*ny*nt_z)
!		zplane_avg_t%tyz_avg(k) = zplane_avg_t%tyz_avg(k) + sum(tyz(:,:,k))/(nx*ny*nt_z)
!		zplane_avg_t%tzz_avg(k) = zplane_avg_t%tzz_avg(k) + sum(tzz(:,:,k))/(nx*ny*nt_z)
!		
!		zplane_avg_t%dudz_avg(k) = zplane_avg_t%dudz_avg(k) + sum(dudz(:,:,k))/(nx*ny*nt_z)
!	    zplane_avg_t%dvdz_avg(k) = zplane_avg_t%dvdz_avg(k) + sum(dvdz(:,:,k))/(nx*ny*nt_z)
!		
!		zplane_avg_t%u2_avg(k) = zplane_avg_t%u2_avg(k) + sum(u(:,:,k)*u(:,:,k))/(nx*ny*nt_z)
!	    zplane_avg_t%v2_avg(k) = zplane_avg_t%v2_avg(k) + sum(v(:,:,k)*v(:,:,k))/(nx*ny*nt_z)
!	    zplane_avg_t%w2_avg(k) = zplane_avg_t%w2_avg(k) + sum(w(:,:,k)*w(:,:,k))/(nx*ny*nt_z)
!		
!		zplane_avg_t%uv_avg(k) = zplane_avg_t%uv_avg(k) + sum(u(:,:,k)*v(:,:,k))/(nx*ny*nt_z)
!	    zplane_avg_t%uw_avg(k) = zplane_avg_t%uw_avg(k) + sum(u(:,:,k)*w(:,:,k))/(nx*ny*nt_z)
!	    zplane_avg_t%vw_avg(k) = zplane_avg_t%vw_avg(k) + sum(v(:,:,k)*w(:,:,k))/(nx*ny*nt_z)		
!	enddo
!elseif (ipart==3) then
!!Write averages and reset
!	  !u-v-w file
!		  call write_real_data_1D(fname2,'append','formatted',3,(nz-1),(/ zplane_avg_t%u_avg,  zplane_avg_t%v_avg,  zplane_avg_t%w_avg /),z(1:nz-1))
!		  zplane_avg_t%u_avg(:) = 0.
!		  zplane_avg_t%v_avg(:) = 0.
!		  zplane_avg_t%w_avg(:) = 0.
!	  !tau file
!		  call write_real_data_1D(fname3,'append','formatted',6,(nz-1),(/ zplane_avg_t%txx_avg,  &
!		  zplane_avg_t%txy_avg, zplane_avg_t%tyy_avg,  zplane_avg_t%txz_avg,  zplane_avg_t%tyz_avg,  zplane_avg_t%tzz_avg /),z(1:nz-1))
!		  zplane_avg_t%txx_avg(:) = 0.
!		  zplane_avg_t%txy_avg(:) = 0.
!		  zplane_avg_t%tyy_avg(:) = 0.
!		  zplane_avg_t%txz_avg(:) = 0.
!		  zplane_avg_t%tyz_avg(:) = 0.
!		  zplane_avg_t%tzz_avg(:) = 0.	 
!	  !ddz file
!		  call write_real_data_1D(fname4,'append','formatted',2,(nz-1),(/ zplane_avg_t%dudz_avg,  zplane_avg_t%dvdz_avg /),z(1:nz-1))
!		  zplane_avg_t%dudz_avg(:) = 0.
!		  zplane_avg_t%dvdz_avg(:) = 0.	
!	  !u2-v2-w2 file
!		  call write_real_data_1D(fname5,'append','formatted',3,(nz-1),(/ zplane_avg_t%u2_avg,  zplane_avg_t%v2_avg,  zplane_avg_t%w2_avg /),z(1:nz-1))
!		  zplane_avg_t%u2_avg(:) = 0.
!		  zplane_avg_t%v2_avg(:) = 0.
!		  zplane_avg_t%w2_avg(:) = 0.	
!	  !uv-uw-vw file
!		  call write_real_data_1D(fname6,'append','formatted',3,(nz-1),(/ zplane_avg_t%uv_avg,  zplane_avg_t%uw_avg,  zplane_avg_t%vw_avg /),z(1:nz-1))
!		  zplane_avg_t%uv_avg(:) = 0.
!		  zplane_avg_t%uw_avg(:) = 0.
!		  zplane_avg_t%vw_avg(:) = 0.			  
!else 
!!error
!	write(*,*) 'ERROR: Incorrect call zplane_avg_compute(ipart)'
!endif


!return
!end subroutine zplane_avg_compute


end module io
