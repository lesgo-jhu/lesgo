module io
use types,only:rprec
use param, only : ld, nx, ny, nz, nz_tot, write_inflow_file, path,  &
                  USE_MPI, coord, rank, nproc, jt_total
implicit none
private

!!$public openfiles,output_loop,output_final,                   &
!!$     inflow_write, avg_stats
public jt_total, openfiles, inflow_read, inflow_write, output_loop, output_final
public mean_u,mean_u2,mean_v,mean_v2,mean_w,mean_w2

!!$ Region commented by JSG 
!!$integer,parameter::base=2000,nwrite=base
!!$
logical, parameter :: cumulative_time = .true.
character (*), parameter :: fcumulative_time = path // 'total_time.dat'
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
                   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains

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
    read (1, *) jt_total
    close (1)
  else  !--assume this is the first run on cumulative time
    write (*, *) 'file ', fcumulative_time, ' not found'
    write (*, *) 'assuming jt_total = 0'
    jt_total = 0
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
use param, only : dx, dy, dz
use stat_defs, only: tsum_t, point_t, domain_t, yplane_t, zplane_t
use sim_param, only : u, v, w
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
if(tsum_t%calc) then
!  Check if we are in the time interval for running summations
  if(jt >= tsum_t%nstart .and. jt <= tsum_t%nend) then
    if(.not. tsum_t%started) then
      if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
        write(*,*) '-------------------------------'   
        write(*,"(1a,i9,1a,i9)") 'Starting running time summation from ', tsum_t%nstart, ' to ', tsum_t%nend
        write(*,*) '-------------------------------'   
      endif
      tsum_t%started=.true.
    endif
!  Compute running summations
    call tsum_compute ()
  endif 
endif

!  Determine if instantaneous point velocities are to be recorded
if(point_t%calc) then
  if(jt >= point_t%nstart .and. jt <= point_t%nend .and. mod(jt,point_t%nskip)==0) then
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
  if(jt >= domain_t%nstart .and. jt <= domain_t%nend .and. mod(jt,domain_t%nskip)==0) then
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
  if(jt >= yplane_t%nstart .and. jt <= yplane_t%nend .and. mod(jt,yplane_t%nskip)==0) then
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
  if(jt >= zplane_t%nstart .and. jt <= zplane_t%nend .and. mod(jt,zplane_t%nskip)==0) then
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

!if(yplane_t%calc .or. zplane_t%calc) call plane_avg_compute(jt)

return
end subroutine output_loop

!**********************************************************************
subroutine inst_write(itype)
!**********************************************************************
!  This subroutine writes the instantaneous values
!  at specified i,j,k locations
use functions, only : linear_interp, trilinear_interp, interp_to_uv_grid
use stat_defs, only : point_t, domain_t, yplane_t, zplane_t
use grid_defs, only : x,y,z
use sim_param, only : u,v,w
use level_set, only : phi
use param, only : jt_total, dt_dim, nx, ny, nz,dx,dy,dz,z_i,L_x,L_y,L_z,coord
implicit none

integer, intent(IN) :: itype

character(25) :: cl, ct
character (64) :: fname, temp
integer :: n, fid, i, j, k

real(rprec) :: ui, vi, wi

! double precision :: dnx, dny, dnz

!  Write point data; assumes files have been opened properly
!  in stats_init

if(itype==1) then

  do n=1,point_t%nloc
!  Files have been opened in stats_init

!  For parallel runs check if data is on correct proc
    $if ($MPI)
    if(point_t%coord(n) == coord) then
    $endif

      fid=3000*n

      write(fid,*) jt_total*dt_dim, &
      trilinear_interp('u',point_t%istart(n),point_t%jstart(n), point_t%kstart(n), point_t%xyz(:,n)), &
      trilinear_interp('v',point_t%istart(n),point_t%jstart(n), point_t%kstart(n), point_t%xyz(:,n)), &
      trilinear_interp('w',point_t%istart(n),point_t%jstart(n), point_t%kstart(n), point_t%xyz(:,n))


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
  
  open(unit = 7,file = fname, status='unknown',form='formatted', &
        action='write',position='rewind')

  $if($LVLSET)
  write(7,*) 'variables = "x", "y", "z", "u", "v", "w", "phi"';
  $else
  write(7,*) 'variables = "x", "y", "z", "u", "v", "w"';
  $endif

  write(7,"(1a,i9,1a,i3,1a,i3,1a,i3,1a,i3)") 'ZONE T="', &
    j,'", DATAPACKING=POINT, i=', Nx,', j=',Ny,', k=', Nz

  $if($LVLSET)
  write(7,"(1a)") ''//adjustl('DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)')//''
  $else
  write(7,"(1a)") ''//adjustl('DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)')//''
  $endif

  write(7,"(1a,f18.6)") 'solutiontime=', jt_total*dt_dim
  
  do k=1,nz
    do j=1,ny
      do i=1,nx
        $if($LVLSET)
        write(7,*) x(i), y(j), z(k), u(i,j,k), v(i,j,k), interp_to_uv_grid('w',i,j,k), phi(i,j,k)
        $else
        write(7,*) x(i), y(j), z(k), u(i,j,k), v(i,j,k), interp_to_uv_grid('w',i,j,k)
        $endif
      enddo
    enddo
  enddo
  
  close(7)
!  Write instantaneous y-plane values
elseif(itype==3) then

!  Loop over all yplane locations
  do j=1,yplane_t%nloc

    write(cl,'(F9.4)') yplane_t%loc(j)
    !  Convert total iteration time to string
    write(ct,*) jt_total
    write(fname,*) 'output/uvw.y-',trim(adjustl(cl)),'.',trim(adjustl(ct)),'.out'
    fname=trim(adjustl(fname))

    $if ($MPI)
!  For MPI implementation
      write (temp, '(".c",i0)') coord
      fname = trim (fname) // temp
    $endif

    open (unit = 2,file = fname, status='unknown',form='formatted', &
      action='write',position='rewind')
    write(2,*) 'variables = "x", "y", "z", "u", "v", "w"';
    write(2,"(1a,i9,1a,i3,1a,i3,1a,i3,1a,i3)") 'ZONE T="', &
      j,'", DATAPACKING=POINT, i=', Nx,', j=',1,', k=', Nz
    write(2,"(1a)") ''//adjustl('DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)')//''
    write(2,"(1a,f18.6)") 'solutiontime=', jt_total*dt_dim
    do k=1,nz
      do i=1,nx

        ui = linear_interp(u(i,yplane_t%istart(j),k), &
          u(i,yplane_t%istart(j)+1,k), dy, yplane_t%ldiff(j))
        vi = linear_interp(v(i,yplane_t%istart(j),k), &
          v(i,yplane_t%istart(j)+1,k), dy, yplane_t%ldiff(j))
        wi = linear_interp(interp_to_uv_grid('w',i,yplane_t%istart(j),k), &
          interp_to_uv_grid('w',i,yplane_t%istart(j)+1,k), dy, &
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
    write(fname,*) 'output/uvw.z-',trim(adjustl(cl)),'.',trim(adjustl(ct)),'.out'
    fname=trim(adjustl(fname))

!     $if ($MPI)
! !  For MPI implementation
!       write (temp, '(".c",i0)') coord
!       fname = trim (fname) // temp
!     $endif

    open (unit = 2,file = fname, status='unknown',form='formatted', &
      action='write',position='rewind')
    write(2,*) 'variables = "x", "y", "z", "u", "v", "w"';
    write(2,"(1a,i9,1a,i3,1a,i3,1a,i3,1a,i3)") 'ZONE T="', &
      j,'", DATAPACKING=POINT, i=', Nx,', j=',Ny,', k=', 1
    write(2,"(1a)") ''//adjustl('DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)')//''
    write(2,"(1a,f18.6)") 'solutiontime=', jt_total*dt_dim

    do j=1,Ny
      do i=1,Nx

        ui = linear_interp(u(i,j,zplane_t%istart(k)),u(i,j,zplane_t%istart(k)+1), &
          dz, zplane_t%ldiff(k))
        vi = linear_interp(v(i,j,zplane_t%istart(k)),v(i,j,zplane_t%istart(k)+1), &
          dz, zplane_t%ldiff(k))
        wi = linear_interp(interp_to_uv_grid('w',i,j,zplane_t%istart(k)), &
          interp_to_uv_grid('w',i,j,zplane_t%istart(k)+1), &
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
subroutine tsum_write()
!**********************************************************************
use grid_defs, only : x,y,z
use stat_defs, only : tsum_t
use param, only : nx,ny,nz,dx,dy,dz,L_x,L_y,L_z
implicit none

character (64) :: fname, temp
integer :: i,j,k

fname = 'output/tsum.out'

$if ($MPI)
!  For MPI implementation     
  write (temp, '(".c",i0)') coord
  fname = trim (fname) // temp
$endif

open (unit = 7,file = fname, status='unknown',form='unformatted', &
  action='write',position='rewind')
! write(7,*) 'variables= "x", "y", "z", "<u>", "<v>", "<w>"'
! write(7,"(1a,i9,1a,i3,1a,i3,1a,i3,1a,i3)") 'ZONE T="', &
!   1,'", DATAPACKING=POINT, i=', Nx,', j=',Ny, ', k=', Nz
! write(7,"(1a)") ''//adjustl('DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)')//''	
!  write(8,*) 'variables= "z", "<dudz>/u*"'
do k=1,nz
  do j=1,ny
    do i=1,nx
       write(7) x(i), y(j), z(k), tsum_t%u(i,j,k), tsum_t%v(i,j,k), tsum_t%w(i,j,k), &
         tsum_t%u(i,j,k), tsum_t%v(i,j,k), tsum_t%w(i,j,k), tsum_t%u2(i,j,k), &
         tsum_t%v2(i,j,k), tsum_t%w2(i,j,k), tsum_t%uw(i,j,k), &
         tsum_t%vw(i,j,k), tsum_t%uv(i,j,k), tsum_t%dudz(i,j,k)
    enddo
  enddo
enddo
close(7)
!  close(8)

return
end subroutine tsum_write

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
use stat_defs, only : rs_t, tsum_t, point_t, yplane_t, zplane_t

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
    write (1, *) jt_total
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
if(tsum_t%calc) call tsum_write()
 
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
end module io
