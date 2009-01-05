module io
use types,only:rprec
use param, only : ld, nx, ny, nz, nz_tot, write_inflow_file, path,  &
                  use_avgslice, USE_MPI, coord, rank, nproc, jt_total
implicit none
private
public openfiles,output_loop,output_final,                   &
     mean_u,mean_u2,mean_v,mean_v2,mean_w,mean_w2,jt_total,  &
     inflow_read, inflow_write, avg_stats

!integer::jt_total  !--moved into param

integer,parameter::base=2000,nwrite=base

logical, parameter :: cumulative_time = .true.
character (*), parameter :: fcumulative_time = path // 'total_time.dat'

logical, parameter :: use_avg_stats = .false.
integer, parameter :: n_avg_stats = 500
                      !--interval for updates in avg_stats
character (*), parameter :: end_hdr_avg = '# end header'

!--write velocity slices (entire planes, instantaneous)
!--this is not related to avgslice
!--this output to be averaged by postprocessing program
logical, parameter :: use_vel_slice = .true.
integer, parameter :: n_vel_slice_write = 1000

!--write immersed bdry force field
logical, parameter :: write_f = .false.

!--write whole txz field or not (for checking momentum balance)
!--does not affect avgslice writing of txz
logical, parameter :: write_txz = .false.

!--write 1d-slices (w/o any processing) at specified stations
logical, parameter :: use_write_1dx = .true.  !--x-slices
logical, parameter :: use_write_1dy = .true.  !--y-slices
integer, parameter :: n_1dx = 1  !--number of x-slices
integer, parameter :: n_1dy = 1  !--number of y-slices
real (rprec) :: y_1dx(n_1dx) = (/ 0.5_rprec /)
                !--y-value for x-slice
real (rprec) :: x_1dy(n_1dy) = (/ 0.5_rprec /)
                !--x-value for y-slice

!!!!  io_spec=.true. output plan-averaged spectrum
!!!!  time_spec>0 output time series spectrum (need additional calcu.)
logical,parameter::io_spec=.false.
integer,parameter::time_spec=0
integer::n_obs
integer,allocatable::obs_pt(:,:)

!!!!  io_mean=.true. output small domain time-averaged velocity
logical,parameter::io_mean=.false.
!integer,parameter::jx_pls=1,jx_ple=nx,width=ny/2-1
integer,parameter::jx_pls=1,jx_ple=1,width=1
integer,parameter::jy_pls=ny/2-width,jy_ple=ny/2+width+1
real(kind=rprec),dimension(jx_pls:jx_ple,jy_pls:jy_ple,nz)::&
     mean_u,mean_v,mean_w,mean_u2,mean_v2,mean_w2

!!!!  io_lambda2
!logical,parameter::io_lambda2=.false.
!real(kind=rprec),dimension(nx,ny,nz)::lam2  !--commented to save mem.
                   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--write out velocity slices for later postprocessing (avging)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine vel_slice ()
use sim_param, only: u, v, w

implicit none

integer, parameter :: nx_skip = nx / 4  !--space between x-slices
integer, parameter :: ny_skip = ny / 4  !--space between y-slices
integer, parameter :: nz_skip = (nz_tot-1) / 4  !--space between z-slices

character (128) :: fname
character (32) :: temp

integer :: kmin, kk

!---------------------------------------------------------------------

!--x-direction slices
write (fname, '(a,i0,a,i6.6,a)') path // 'output/vel.xskip', nx_skip,  &
                                   '.', jt_total, '.out'

$if ($MPI)
  write (temp, '(".c",i0)') coord
  fname = trim (fname) // temp
$endif

open(1,file=fname,form='unformatted')

write (1) u(1:nx:nx_skip, 1:ny, 1:nz), v(1:nx:nx_skip, 1:ny, 1:nz),  &
          w(1:nx:nx_skip, 1:ny, 1:nz)

close(1)

!--y-direction slices
write (fname, '(a,i0,a,i6.6,a)') path // 'output/vel.yskip', ny_skip,  &
                                   '.', jt_total, '.out'

$if ($MPI)
  write (temp, '(".c",i0)') coord
  fname = trim (fname) // temp
$endif

open(1,file=fname,form='unformatted')

write (1) u(1:nx, 1:ny:ny_skip, 1:nz), v(1:nx, 1:ny:ny_skip, 1:nz),  &
          w(1:nx, 1:ny:ny_skip, 1:nz)

close(1)

!--z-direction slices
write (fname, '(a,i0,a,i6.6,a)') path // 'output/vel.zskip', nz_skip,  &
                                   '.', jt_total, '.out'

$if ($MPI)

  !--first global k for corresponding to this process
  !  kg = 1 + kk * nz_skip
  kk = ceiling (real (coord * (nz - 1), rprec) / nz_skip)
  !--this corresponds to local k of kmin:
  kmin = 1 + kk * nz_skip - coord * (nz - 1)

  write (temp, '(".c",i0)') coord
  fname = trim (fname) // temp
  
$else

  kmin = 1

$endif

if ((1 <= kmin) .and. (kmin <= nz-1)) then  !--no write if out of range

  open (1, file=fname, form='unformatted')

  write (1) u(1:nx, 1:ny, kmin:nz-1:nz_skip),  &
            v(1:nx, 1:ny, kmin:nz-1:nz_skip),  &
            w(1:nx, 1:ny, kmin:nz-1:nz_skip)

  close (1)

end if

end subroutine vel_slice

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! file number 1-10 are used for temporary use
! 11-19 are basic output
! 20-40 are avgslice
! use >50 for debugging
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine openfiles()
use param,only:sflux_flag
use sim_param,only:path
implicit none

!--to hold file names
character (64) :: temp
character (64) :: fCS1plan, fCS2plan, fCS4plan, fVISCplan,  &
                  fDISSplan, fCS1Vplan, fCS2Vplan, fCS4Vplan

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

if(time_spec.gt.0)then
  open(15,file=path//'output/velspec.out',form='unformatted',position='append')
  if(jt_total.eq.0)rewind(15)
endif

if(io_mean)then
  open(51,file=path//'output/mean_u.out',form='unformatted',position='append')
  if(jt_total.eq.0)then
    rewind(51)
    write(51)jx_pls,jx_ple,jy_pls,jy_ple
  endif
endif

fCS1plan = path // 'output/CS1plan.out'
fCS2plan = path // 'output/CS2plan.out'
fCS4plan = path // 'output/CS4plan.out'
fVISCplan = path // 'output/VISCplan.out'
fDISSplan = path // 'output/DISSplan.out'
fCS1Vplan = path // 'output/CS1Vplan.out'
fCS2Vplan = path // 'output/CS2Vplan.out'
fCS4Vplan = path // 'output/CS4Vplan.out'

$if ($MPI)
  !--append coordinate identifiers
  write (temp, '(".c",i0)') coord
  fCS1plan = trim (fCS1plan) // temp
  fCS2plan = trim (fCS2plan) // temp
  fCS4plan = trim (fCS4plan) // temp
  fVISCplan = trim (fVISCplan) // temp
  fDISSplan = trim (fDISSplan) // temp
  fCS1Vplan = trim (fCS1Vplan) // temp
  fCS2Vplan = trim (fCS2Vplan) // temp
  fCS4Vplan = trim (fCS4Vplan) // temp
$endif

open (90, file=fCS1plan, form='unformatted')
open (91, file=fCS2plan, form='unformatted')
open (92, file=fCS4plan, form='unformatted')
open (93, file=fVISCplan, form='unformatted')
open (94, file=fDISSplan, form='unformatted')
open (95, file=fCS1Vplan, form='unformatted')
open (96, file=fCS2Vplan, form='unformatted')
open (97, file=fCS4Vplan, form='unformatted')
!TSif(sflux_flag)open(98,file='/home/bluesky/ytseng/RESEARCH/JHU_LES/CHANNEL/output/series_data.dat',form='unformatted')

if(time_spec.gt.0)then
open(1,file=path//'obs.pt')
read(1,*)n_obs
allocate(obs_pt(1:2,n_obs))
do i=1,n_obs
read(1,*)obs_pt(1:2,i)
enddo
close(1)
endif

end subroutine openfiles

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output_loop(jt)
use param,only:output,dt,c_count,S_FLAG,SCAL_init
use sim_param,only:path,u,v,w,dudz,dudx,p,&
     RHSx,RHSy,RHSz,theta, txx, txy, txz, tyy, tyz, tzz
use sgsmodule,only:Cs_opt2!,Cs_opt2_avg  !--not essential, save mem.
use scalars_module2,only:scalar_slice
use immersedbc, only : fx, fy, fz
implicit none
integer,intent(in)::jt
!real(kind=rprec),dimension(ld,ny,nz)::usp,vsp  !--not used
real(kind=rprec),dimension(nz)::u_ndim
character(len=20)::req

character (64) :: fname, temp

integer::jx,jy,jz

!jt_total=jt_total+1  !--moved into main program

if (io_mean) call calculate_mean ()

if (use_avgslice) then
  if (mod(jt,c_count)==0) then
    call avgslice(jt)
    if ((S_FLAG) .and. (jt_total.GE.SCAL_init)) then ! Output scalar variables
      call scalar_slice(jt) ! Uses file unit numbers (36-47)
    end if
  end if
end if

if (output) then

  if (modulo (jt_total, n_vel_slice_write) == 0) call vel_slice ()

  if (mod(jt_total,base)==0) then

    write (fname, '(a,i6.6,a)') path // 'output/vel', jt_total, '.out'
    $if ($MPI)
      write (temp, '(".c",i0)') coord
      fname = trim (fname) // temp
    $endif

    open(1,file=fname,form='unformatted')

    call checkpoint (1)

    close(1)

    if (write_txz) then
      !--needed for validation (in overall stress balance)
      write (fname, '(a,i6.6,a)') path // 'output/txz', jt_total, '.out'
      $if ($MPI)
        write (temp, '(".c",i0)') coord
        fname = trim (fname) // temp
      $endif
      open (1, file=fname, form='unformatted')
      write (1) txz(:, :, 1:nz)
      close (1)
    end if

    if (write_f) then
      !--mainly for debugging
      write (fname, '(a,i6.6,a)') path // 'output/f', jt_total, '.out'
      $if ($MPI)
        write (temp, '(".c",i0)') coord
        fname = trim (fname) // temp
      $endif
      open (1, file=fname, form='unformatted')
      write (1) fx(:, :, 1:nz), fy(:, :, 1:nz), fz(:, :, 1:nz)
      close (1)
    end if

  end if

  if (mod(jt_total,nwrite)==0) then
    !--this does not make sense to me, instead added p, theta to above
    !--also Cs_opt2_avg commented out of sgsmodule
    !write(fname,'(A13,i6.6,A4)')path//'output/vel-',jt_total,'.out'
    !open(1,file=fname,form='unformatted')
    !if(S_FLAG)THEN
    !write(1) real(u),real(v),real(w),real(Cs_opt2_avg),real(p),real(theta)
    !else
    !write(1) real(u),real(v),real(w),real(Cs_opt2_avg),real(p)
    !end if
    !Cs_opt2_avg(:,:,:)=0._rprec
    !close(1)
    if(io_mean)call io_mean_out
    !if(io_lambda2)call io_lambda2_out
    if(io_spec)call post_spec(jt)
  end if
  if(time_spec.gt.0)call timeseries_spec
end if
end subroutine output_loop

$if (0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_1d ()
implicit none
use sim_param, only : u, v, w

!---------------------------------------------------------------------

if (use_write_1dx) then

end if

if (use_write_1dy) then

end if

end subroutine write_1d
$endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine avgslice (jt)
use sim_param,only:path,u,v,w,dudz,dvdz, txx, txz, tyy, tyz, tzz, p
use param,only:dz,p_count,c_count
use sgsmodule,only:Cs_opt2
implicit none
integer::i,j,k, jt
real(kind=rprec),dimension(nx,nz)::ap,au,av,aw,p2,u2,v2,w2,auw,avw,acs
real(kind=rprec),dimension(nx,nz)::adudz,advdz
real(kind=rprec),dimension(nx,nz)::atxx,atxz,atyy,atyz,atzz
real(kind=rprec)::tu1,tv1,tw1,ttxx,ttxz,ttyy,ttyz,ttzz,tdudz,tdvdz,&
     tu2,tv2,tw2,tp1,tp2,tuw,tvw,tCs,fr, arg1, arg2     


fr=(1._rprec/real(p_count,kind=rprec))*real(c_count,kind=rprec)
do k=1,Nz
do i=1,Nx
   tu1=0._rprec;tv1=0._rprec;tw1=0._rprec;tp1=0._rprec
   ttxx=0._rprec;ttxz=0._rprec;ttyy=0._rprec;ttyz=0._rprec
   ttzz=0._rprec;tdudz=0._rprec;tdvdz=0._rprec;tu2=0._rprec
   tv2=0._rprec;tw2=0._rprec;tp2=0._rprec;tuw=0._rprec;tvw=0._rprec
   tCs=0._rprec

   do j=1,Ny
      tu1=tu1+u(i,j,k)
      tv1=tv1+v(i,j,k)
      tw1=tw1+w(i,j,k)
      tp1=tp1+p(i,j,k)
      ttxx=ttxx+txx(i,j,k)
      ttxz=ttxz+txz(i,j,k)
      ttyy=ttyy+tyy(i,j,k)
      ttyz=ttyz+tyz(i,j,k)
      ttzz=ttzz+tzz(i,j,k)
      tdudz=tdudz+dudz(i,j,k)
      tdvdz=tdvdz+dvdz(i,j,k)
      tu2=tu2+u(i,j,k)*u(i,j,k)
      tv2=tv2+v(i,j,k)*v(i,j,k)
      tw2=tw2+w(i,j,k)*w(i,j,k)
      tp2=tp2+p(i,j,k)*p(i,j,k)
      tCs=tCs+Cs_opt2(i,j,k)
      if (k.gt.1) then
         arg1=(u(i,j,k)+u(i,j,k-1))/2.
         arg2=(v(i,j,k)+v(i,j,k-1))/2.
      else
         arg1=0._rprec
         arg2=0._rprec
      end if
      tuw=tuw+w(i,j,k)*arg1
      tvw=tvw+w(i,j,k)*arg2
   end do
   au(i,k)=au(i,k)+(fr)*tu1/Ny
   av(i,k)=av(i,k)+(fr)*tv1/Ny
   aw(i,k)=aw(i,k)+(fr)*tw1/Ny
   ap(i,k)=ap(i,k)+(fr)*tp1/Ny
   adudz(i,k)=adudz(i,k)+(fr)*tdudz/Ny
   advdz(i,k)=advdz(i,k)+(fr)*tdvdz/Ny
   u2(i,k)=u2(i,k)+(fr)*tu2/Ny
   v2(i,k)=v2(i,k)+(fr)*tv2/Ny
   w2(i,k)=w2(i,k)+(fr)*tw2/Ny
   atxx(i,k)=atxx(i,k)+(fr)*ttxx/Ny
   atxz(i,k)=atxz(i,k)+(fr)*ttxz/Ny
   atyy(i,k)=atyy(i,k)+(fr)*ttyy/Ny
   atyz(i,k)=atyz(i,k)+(fr)*ttyz/Ny
   atzz(i,k)=atzz(i,k)+(fr)*ttzz/Ny
   p2(i,k)=p2(i,k)+fr*tp2/Ny
   aCs(i,k)=aCs(i,k)+(fr)*tCs/Ny
   auw(i,k)=auw(i,k)+(fr)*tuw/Ny
   avw(i,k)=avw(i,k)+(fr)*tvw/Ny
end do
end do

if(mod(jt,p_count)==0.and.jt_total.gt.2000) then
  open(20,file=path//'output/aver_u.out',status="unknown",position="append")
  open(21,file=path//'output/aver_v.out',status="unknown",position="append")
  open(22,file=path//'output/aver_w.out',status="unknown",position="append")
  open(23,file=path//'output/aver_p.out',status="unknown",position="append")
  open(24,file=path//'output/aver_u2.out',status="unknown",position="append")
  open(25,file=path//'output/aver_v2.out',status="unknown",position="append")
  open(26,file=path//'output/aver_w2.out',status="unknown",position="append")
  open(27,file=path//'output/aver_txx.out',status="unknown",position="append")
  open(28,file=path//'output/aver_txz.out',status="unknown",position="append")
  open(29,file=path//'output/aver_tyy.out',status="unknown",position="append")
  open(30,file=path//'output/aver_tyz.out',status="unknown",position="append")
  open(31,file=path//'output/aver_tzz.out',status="unknown",position="append")
  open(32,file=path//'output/aver_p2.out',status="unknown",position="append")
  open(33,file=path//'output/aver_uw.out',status="unknown",position="append")
  open(34,file=path//'output/aver_vw.out',status="unknown",position="append")
  open(35,file=path//'output/aver_Cs2.out',status="unknown",position="append")
  open(36,file=path//'output/aver_dudz.out',status="unknown",position="append") 
  open(37,file=path//'output/aver_dvdz.out',status="unknown",position="append")
  do k=1,nz
      write(20,5168) real(jt_total*dz),(real(au(i,k)),i=1,nx)
      write(21,5168) real(jt_total*dz),(real(av(i,k)),i=1,nx)
      write(22,5168) real(jt_total*dz),(real(aw(i,k)),i=1,nx)
      write(23,5168) real(jt_total*dz),(real(ap(i,k)),i=1,nx)
      write(24,5168) real(jt_total*dz),(real(u2(i,k)),i=1,nx)
      write(25,5168) real(jt_total*dz),(real(v2(i,k)),i=1,nx)
      write(26,5168) real(jt_total*dz),(real(w2(i,k)),i=1,nx)
      write(27,5168) real(jt_total*dz),(real(atxx(i,k)),i=1,nx)
      write(28,5168) real(jt_total*dz),(real(atxz(i,k)),i=1,nx)
      write(29,5168) real(jt_total*dz),(real(atyy(i,k)),i=1,nx)
      write(30,5168) real(jt_total*dz),(real(atyz(i,k)),i=1,nx)
      write(31,5168) real(jt_total*dz),(real(atzz(i,k)),i=1,nx)
      write(32,5168) real(jt_total*dz),(real(p2(i,k)),i=1,nx)
      write(33,5168) real(jt_total*dz),(real(auw(i,k)),i=1,nx)
      write(34,5168) real(jt_total*dz),(real(avw(i,k)),i=1,nx)
      write(35,5168) real(jt_total*dz),(real(aCs(i,k)),i=1,nx)
      write(36,5168) real(jt_total*dz),(real(adudz(i,k)),i=1,nx)
      write(37,5168) real(jt_total*dz),(real(advdz(i,k)),i=1,nx)
  end do

!VK Zero out the outputted averages !!
   au=0._rprec;av=0._rprec;aw=0._rprec;ap=0._rprec;u2=0._rprec;v2=0._rprec
   w2=0._rprec;atxx=0._rprec;atxz=0._rprec;atyy=0._rprec;atyz=0._rprec
   atzz=0._rprec;p2=0._rprec;auw=0._rprec;avw=0._rprec;aCs=0._rprec
   adudz=0._rprec;advdz=0._rprec

! Close the files as soon as the writing is done as one needs to open them
! only when mod(jt,p_count)==0 !!
   close(19);close(20);close(22);close(23);close(24);close(25)
   close(26);close(27);close(28);close(29);close(30)
   close(31);close(33);close(34);close(35);close(36);close(37)

end if
 5168     format(1400(E14.5))
end subroutine avgslice

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
implicit none

integer,intent(in)::jt
integer, intent (in), optional :: lun_opt  !--if present, write to unit lun

integer, parameter :: lun_default = 11

integer::jx,jy,jz
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

if ((cumulative_time) .and. (lun == lun_default)) then

  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
    !--only do this for true final output, not intermediate recording
    open (1, file=fcumulative_time)
    write (1, *) jt_total
    close (1)
  end if

end if

end subroutine output_final

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine io_lambda2_out
!use sim_param,only:path
!implicit none
!character(len=24)::fname
!call lambda2()
!write(fname,'(A13,i6.6,A4)')path//'output/lam-',jt_total,'.out'
!open(1,file=fname,form='unformatted')
!write(1)nx,ny,nz
!write(1)real(lam2)
!close(1)
!end subroutine io_lambda2_out

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine lambda2()
!use types,only:rprec
!use sim_param,only:u,v,w,&
!     dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
!use param,only:dx,dy,dz
!implicit none
!!TSreal(kind=rprec),dimension(nx)::dpdx,dpdy,dpdz,aeu,awu,anu,asu,afu,abu
!!TSreal(kind=rprec)::fne,fnn,fnf,S11,S22,S33,S12,O12,S13,O13,S23,O23
!real(kind=rprec)::S11,S22,S33,S12,O12,S13,O13,S23,O23,&
!     ux,uy,uz,vx,vy,vz,wx,wy,wz
!integer::jx,jy,jz
!! following used for eispack call...
!integer::neis,nmeis,matzeis,ierreis,iv1eis(3)
!double precision,dimension(3,3)::aeis,zeis
!double precision,dimension(3)::wreis,wieis,fv1eis
!double precision::ave
!
!! assignments for eispack call
!neis=3
!nmeis=3
!matzeis=0
!ierreis=0
!!TSfne=1._rprec/dx
!!TSfnn=1._rprec/dy
!!TSfnf=1._rprec/dz
!!TSprint*,fne,fnn,fnf
!lam2=0._rprec
!
!! at level z=dz/2.  (Where dudz and dvdz are on UVP nodes)
!jz=1
!do jy=1,ny
!do jx=1,nx              
!   ux=dudx(jx,jy,1)  ! uvp-node
!   uy=dudy(jx,jy,1)  ! uvp-node
!   uz=dudz(jx,jy,1)  ! uvp-node
!   vx=dvdx(jx,jy,1)  ! uvp-node
!   vy=dvdy(jx,jy,1)  ! uvp-node
!   vz=dvdz(jx,jy,1)  ! uvp-node 
!! special case
!   wx=0.5_rprec*(dwdx(jx,jy,1)+dwdx(jx,jy,2))  ! uvp-node
!   wy=0.5_rprec*(dwdy(jx,jy,1)+dwdy(jx,jy,2))  ! uvp-node
!   wz=dwdz(jx,jy,1)  ! uvp-node
!   S11=ux          ! uvp-node
!   S12=0.5_rprec*(uy+vx) ! uvp-node
!! taken care of with wall stress routine
!   S13=0.5_rprec*(uz+wx) ! uvp
!   O12=0.5_rprec*(uy-vx) ! w-node
!   O13=0.5_rprec*(uz-wx) ! w-node
!   S22=vy          ! uvp-node
!! taken care of with wall stress routine 
!   S23=0.5_rprec*(vz+wy) ! uvp
!   O23=0.5_rprec*(vz-wy) ! w-node
!   S33=wz          ! uvp-node
!   aeis(1,1)=s11*s11+s12*s12+s13*s13-O12*O12-O13*O13
!   aeis(1,2)=s11*s12+s12*s22+s13*s23-O13*O23
!   aeis(1,3)=s11*s13+s12*s23+s13*s33+O12*O23
!   aeis(2,2)=s12*s12+s22*s22+s23*s23-O12*O12-O23*O23
!   aeis(2,3)=s12*s13+s22*s23+s23*s33-O12*O13
!   aeis(3,3)=s13*s13+s23*s23+s33*s33-O13*O13-O23*O23
!   aeis(2,1)=aeis(1,2)
!   aeis(3,1)=aeis(1,3)
!   aeis(3,2)=aeis(2,3)
!  write (*, *) 'rg temporarily removed, sorry'; stop
!  !call rg(nmeis,neis,aeis,wreis,wieis,matzeis,zeis,iv1eis,fv1eis,ierreis)
!   if(wreis(1).ge.wreis(2).and.wreis(1).le.wreis(3)) then
!      lam2(jx,jy,jz)=real(wreis(1),kind=rprec)
!   elseif(wreis(2).ge.wreis(1).and.wreis(2).le.wreis(3)) then
!      lam2(jx,jy,jz)=real(wreis(2),kind=rprec)
!   else
!      lam2(jx,jy,jz)=real(wreis(3),kind=rprec)
!   endif
!end do
!end do
!! calculate derivatives/strain on w-nodes
!do jz=2,nz-1  
!do jy=1,ny
!do jx=1,nx              
!   ux=0.5_rprec*(dudx(jx,jy,jz) + dudx(jx,jy,jz-1))  ! w-node
!   uy=0.5_rprec*(dudy(jx,jy,jz) + dudy(jx,jy,jz-1))  ! w-node
!   uz=dudz(jx,jy,jz)  ! w-node
!   vx=0.5_rprec*(dvdx(jx,jy,jz) + dvdx(jx,jy,jz-1))  ! w-node
!   vy=0.5_rprec*(dvdy(jx,jy,jz) + dvdy(jx,jy,jz-1))  ! w-node
!   vz=dvdz(jx,jy,jz)  ! w-node
!   wx=dwdx(jx,jy,jz)  ! w-node
!   wy=dwdy(jx,jy,jz)  ! w-node
!   wz=0.5_rprec*(dwdz(jx,jy,jz) + dwdz(jx,jy,jz-1))  ! w-node
!   S11=ux          ! w-node
!   S12=0.5_rprec*(uy+vx) ! w-node
!   S13=0.5_rprec*(uz+wx) ! w-node
!   O12=0.5_rprec*(uy-vx) ! w-node
!   O13=0.5_rprec*(uz-wx) ! w-node
!   S22=vy          ! w-node
!   S23=0.5_rprec*(vz+wy) ! w-node
!   O23=0.5_rprec*(vz-wy) ! w-node
!   S33=wz          ! w-node
!   aeis(1,1)=s11*s11+s12*s12+s13*s13-O12*O12-O13*O13
!   aeis(1,2)=s11*s12+s12*s22+s13*s23-O13*O23
!   aeis(1,3)=s11*s13+s12*s23+s13*s33+O12*O23
!   aeis(2,2)=s12*s12+s22*s22+s23*s23-O12*O12-O23*O23
!   aeis(2,3)=s12*s13+s22*s23+s23*s33-O12*O13
!   aeis(3,3)=s13*s13+s23*s23+s33*s33-O13*O13-O23*O23
!   aeis(2,1)=aeis(1,2)
!   aeis(3,1)=aeis(1,3)
!   aeis(3,2)=aeis(2,3)
!   !call rg(nmeis,neis,aeis,wreis,wieis,matzeis,zeis,iv1eis,fv1eis,ierreis)
!   if(wreis(1).ge.wreis(2).and.wreis(1).le.wreis(3)) then
!      lam2(jx,jy,jz)=real(wreis(1),kind=rprec)
!   elseif(wreis(2).ge.wreis(1).and.wreis(2).le.wreis(3)) then
!      lam2(jx,jy,jz)=real(wreis(2),kind=rprec)
!   else
!      lam2(jx,jy,jz)=real(wreis(3),kind=rprec)
!   endif
!end do
!end do
!end do
!
!!do k=2,nz-1
!!do j=2,ny-1
!!do i=2,nx-1
!!   dpdx(i) = u(i+1,j,k) - u(i-1,j,k)
!!   dpdy(i) = u(i,j+1,k) - u(i,j-1,k)
!!   dpdz(i) = u(i,j,k+1) - u(i,j,k-1)
!!   aeu(i)  = v(i+1,j,k) - v(i-1,j,k)
!!   awu(i)  = v(i,j+1,k) - v(i,j-1,k)
!!   anu(i)  = v(i,j,k+1) - v(i,j,k-1)
!!   asu(i)  = w(i+1,j,k) - w(i-1,j,k)
!!   afu(i)  = w(i,j+1,k) - w(i,j-1,k)
!!   abu(i)  = w(i,j,k+1) - w(i,j,k-1)
!
!!   S11=fne*dpdx(i)
!!   S22=fnn*awu(i)
!!   S33=fnf*abu(i)
!!   S12=0.5_rprec*(fnn*dpdy(i)+fne*aeu(i))
!!   O12=0.5_rprec*(fnn*dpdy(i)-fne*aeu(i))
!!   S13=0.5_rprec*(fnf*dpdz(i)+fne*asu(i))
!!   O13=0.5_rprec*(fnf*dpdz(i)-fne*asu(i))
!!   S23=0.5_rprec*(fnf*anu(i)+fnn*afu(i))
!!   O23=0.5_rprec*(fnf*anu(i)-fnn*afu(i))
!!   aeis(1,1)=s11*s11+s12*s12+s13*s13-O12*O12-O13*O13
!!   aeis(1,2)=s11*s12+s12*s22+s13*s23-O13*O23
!!   aeis(1,3)=s11*s13+s12*s23+s13*s33+O12*O23
!!   aeis(2,2)=s12*s12+s22*s22+s23*s23-O12*O12-O23*O23
!!   aeis(2,3)=s12*s13+s22*s23+s23*s33-O12*O13
!!   aeis(3,3)=s13*s13+s23*s23+s33*s33-O13*O13-O23*O23
!!   aeis(2,1)=aeis(1,2)
!!   aeis(3,1)=aeis(1,3)
!!   aeis(3,2)=aeis(2,3)
!!   call rg(nmeis,neis,aeis,wreis,wieis,matzeis,zeis,iv1eis,fv1eis,ierreis)
!!   if(wreis(1).ge.wreis(2).and.wreis(1).le.wreis(3)) then
!!      lam2(i,j,k)=real(wreis(1),kind=rprec)
!!   elseif(wreis(2).ge.wreis(1).and.wreis(2).le.wreis(3)) then
!!      lam2(i,j,k)=real(wreis(2),kind=rprec)
!!   else
!!      lam2(i,j,k)=real(wreis(3),kind=rprec)
!!   endif
!!enddo
!!enddo
!!enddo
!print*,'minmax',minval(lam2),maxval(lam2)
!end subroutine lambda2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine io_mean_out
implicit none
write(51)real(mean_u),real(mean_u2),real(mean_v),real(mean_v2),&
     real(mean_w),real(mean_w2)
!--the nz/4*3 stuff has to go
mean_u(jx_pls:jx_ple,jy_pls:jy_ple,1:nz/4*3)=0._rprec
mean_u2(jx_pls:jx_ple,jy_pls:jy_ple,1:nz/4*3)=0._rprec
mean_v(jx_pls:jx_ple,jy_pls:jy_ple,1:nz/4*3)=0._rprec
mean_v2(jx_pls:jx_ple,jy_pls:jy_ple,1:nz/4*3)=0._rprec
mean_w(jx_pls:jx_ple,jy_pls:jy_ple,1:nz/4*3)=0._rprec
mean_w2(jx_pls:jx_ple,jy_pls:jy_ple,1:nz/4*3)=0._rprec
end subroutine io_mean_out

subroutine calculate_mean
use sim_param,only:u,v,w
use sgsmodule,only:Cs_opt2!,Cs_opt2_avg
                          !--Cs_opt2_avg commented to save mem.
implicit none

!--commented to save mem.
!Cs_opt2_avg(:,:,:)=Cs_opt2_avg(:,:,:)+Cs_opt2(:,:,:)/nwrite

!TS
mean_u(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)=&
     mean_u(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)+&
     u(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)/nwrite
mean_u2(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)=&
     mean_u2(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)+&
     u(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)**2/nwrite
mean_v(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)=&
     mean_v(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)+&
     v(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)/nwrite
mean_v2(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)=&
     mean_v2(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)+&
     v(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)**2/nwrite
mean_w(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)=&
     mean_w(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)+&
     w(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)/nwrite
mean_w2(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)=&
     mean_w2(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)+&
     w(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)**2/nwrite
end subroutine calculate_mean

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine timeseries_spec
use sim_param,only:u,v,w,theta
implicit none
integer::jx,jy,jz,i
if(mod(jt_total,time_spec)==0.and.jt_total.gt.2000)then
jx=NX/8
jy=NY/2+1
jz=NZ/2
!TSwrite(15)real(u(jx+NX/24*2,jy,jz:jz+3)),real(u(jx+NX/24*4,jy,jz:jz+3)),&
!TS     real(u(jx+NX/24*6,jy,jz:jz+3)),&
!TS     real(v(jx+NX/24*2,jy,jz:jz+3)),real(v(jx+NX/24*4,jy,jz:jz+3)),&
!TS     real(v(jx+NX/24*6,jy,jz:jz+3)),&
!TS     real(w(jx+NX/24*2,jy,jz:jz+3)),real(w(jx+NX/24*4,jy,jz:jz+3)),&
!TS     real(w(jx+NX/24*6,jy,jz:jz+3))
endif
end subroutine timeseries_spec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine post_spec(jt)
use sim_param,only:path,u,v,w
use param,only:dz,z_i,pi,L_x
use fft,only:FFTW_REAL_TO_COMPLEX,FFTW_ESTIMATE
implicit none
real(kind=rprec),dimension(nx/2,nz)::spectra_u,spectra_v,spectra_w,&
     spectra_theta
integer,intent(in)::jt
integer::k,jz,z
character(len=25)::fname1,fname2,fname3,fname4
if(jt.eq.nwrite)then
write(fname1,'(A15,A4)')path//'output/spec_x','.dat'
open(1,file=fname1,form='formatted')
do jz=1,nz
   z=(jz-0.5)*dz*z_i
   write(1,*) (real(2*pi/L_x*k/z_i*z),k=1,nx/2-1)
enddo
close(1)
endif
!TSaspec_theta=0._rprec
do k=1,nz
   call spectrum(u(:, :, k), spectra_u(:,k))
!          print *,'U spectrum calculated'
   call spectrum(v(:, :, k), spectra_v(:,k))
!          print *,'V spectrum calculated'
   call spectrum(w(:, :, k), spectra_w(:,k))
!          print *,'W spectrum calculated'
!TS   call spectrum(theta,k,spectra_theta(:,k),plan)
!          print *,'THETA spectrum calculated'
end do
print *,'spectrum calculated'
write(fname1,'(A15,i6.6,A4)')path//'output/spec_u',jt_total,'.dat'
open(1,file=fname1,form='formatted')
write(fname2,'(A15,i6.6,A4)')path//'output/spec_v',jt_total,'.dat'
open(2,file=fname2,form='formatted')
write(fname3,'(A15,i6.6,A4)')path//'output/spec_w',jt_total,'.dat'
open(3,file=fname3,form='formatted')
!TSwrite(fname4,'(A14,i6.6,A4)')path//'output/spec_t',jt_total,'.dat'
!TSopen(4,file=fname4,form='formatted')
do jz=1,nz
   z=(jz-0.5)*dz*z_i
   write(1,*)real(spectra_u(2:nx/2,jz))
   write(2,*)real(spectra_v(2:nx/2,jz))
   write(3,*)real(spectra_w(2:nx/2,jz))
!TS   write(4,*)real(spectra_theta(2:nx/2,jz))
enddo
close(1);close(2);close(3);close(4)
end subroutine post_spec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine spectrum(u, spec)
use fft
implicit none      
real(kind=rprec),dimension(ld,ny),intent(in)::u
real(kind=rprec),dimension(nx/2),intent(out)::spec  !--assumes Nyquist is 0

integer::jy,jz,k
real(kind=rprec),dimension(nx)::vel_r,vel_c

integer*8, save :: plan
logical, save :: init = .false.

if (.not. init) then
  call rfftw_f77_create_plan(plan,nx,FFTW_REAL_TO_COMPLEX,FFTW_MEASURE)
  init = .true.
end if

! initialize
spec(:)=0._rprec
do jy=1,ny
   vel_r(:)= u(1:nx,jy)/real(nx,kind=rprec)
! check this normaliztion-part of forward
! call the fft
   call rfftw_f77_one(plan,vel_r,vel_c)
! compute magnitudes
! the 0.5 is the 1/2, all others are taken care of! (except maybe Nyquist)
   spec(1)=spec(1)+0.5*vel_c(1)*vel_c(1)
   do k=2,nx/2
      spec(k)=spec(k)+vel_c(k)*vel_c(k)+vel_c(nx+2-k)*vel_c(nx+2-k)
!        print *,'k,vel,spec',k,vel_c(k),spec(k)
   end do

   !--assume Nyquist is 0
   !spec(nx/2+1)=spec(nx/2+1)+vel_c(nx/2+1)*vel_c(nx/2+1)
end do
spec(:)=spec(:)/real(Ny,kind=rprec) ! for average over Ny
end subroutine spectrum

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine avg_stats ()
use param
use sim_param, only : u, v, w, txz
use fft, only : kx
implicit none

!--choose naming convention that does not conflict with qpost
character (*), parameter :: fubar_avg = 'output/ubar-avg_stats.dat'
character (*), parameter :: fupr2bar_avg = 'output/upr2bar-avg_stats.dat'
character (*), parameter :: fstressbar_avg = 'output/stressbar-avg_stats.dat'
character (*), parameter :: fEozbar_avg = 'output/Eozbar-avg_stats.dat'

integer, parameter :: hdr_len = 256

logical, parameter :: DEBUG = .false.

character (hdr_len) :: Eozbar_hdr

$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif
$if ($MPI)
  integer :: recvcounts(nproc)
  integer :: displs(nproc)
$endif
integer, save :: n_ubar_avg
integer, save :: n_upr2bar_avg
integer, save :: n_stressbar_avg
integer, save :: n_Eozbar_avg
integer :: jz

logical, save :: init = .false.

real (rprec) :: z
real (rprec) :: zu(1, nz_tot-1)
real (rprec) :: kz_z(2, nx/2)
real (rprec), save :: ubar_avg(1, nz_tot-1)  !--1 is <u>
real (rprec), save :: upr2bar_avg(3, nz_tot-1)  !--<u'^2>, <v'^2>, <w'^2>
real (rprec), save :: stressbar_avg(3, nz_tot-1)
                      !--1 is <u'w'>, 2 is <txz>, 3 is <u'w'> + <txz>
real (rprec), save :: Eozbar_avg(1, nx/2, nz_tot-1)  !--E11(k1,z)/z
!--tot is a temp for current stats at nz_tot size
real (rprec), save :: ubar_tot(1, nz_tot-1)  !--1 is <u>
real (rprec), save :: upr2bar_tot(3, nz_tot-1)  !--<u'^2>, <v'^2>, <w'^2>
real (rprec), save :: stressbar_tot(3, nz_tot-1)
                      !--1 is <u'w'>, 2 is <txz>, 3 is <u'w'> + <txz>
real (rprec), save :: Eozbar_tot(1, nx/2, nz_tot-1)  !--E11(k1,z)/z
real (rprec) :: upr(nx, ny), vpr(nx, ny), wpr(nx, ny)
real (rprec) :: ubar(nz-1), vbar(nz-1), wbar(nz-1)
real (rprec) :: upr2bar(3, nz-1)
real (rprec) :: stressbar(3, nz-1)
real (rprec) :: Eozbar(nx/2, nz-1)

!---------------------------------------------------------------------

if (.not. use_avg_stats) goto 001  !--exit cleanly

!--check whether or not to actually do anything
!--motivation for doing this way is that it cleans up interface in main
if (modulo (jt, n_avg_stats) /= 0) goto 001  !--do nothing, exit cleanly

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then  
  if (.not. init) then  !--initialization

    call init_avg (fubar_avg, 1, ubar_avg, n_ubar_avg)
    call init_avg (fupr2bar_avg, 1, upr2bar_avg, n_upr2bar_avg)
    call init_avg (fstressbar_avg, 1, stressbar_avg, n_stressbar_avg) 
    do jz = 1, nz-2
      call init_avg (fEozbar_avg, 2, Eozbar_avg(:, :, jz), n_Eozbar_avg,  &
                     leaveopn='yes')
    end do
    call init_avg (fEozbar_avg, 2, Eozbar_avg(:, :, nz-1), n_Eozbar_avg)

    init = .true.

  end if
end if

!--calculation of current stats
do jz = $lbz, nz-1

  ubar(jz) = sum (u(1:nx, 1:ny, jz)) / (nx * ny)
  vbar(jz) = sum (v(1:nx, 1:ny, jz)) / (nx * ny)
  wbar(jz) = sum (w(1:nx, 1:ny, jz)) / (nx * ny)

  if ( ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) .and.  &
       (jz == 1) ) then
    upr = 0._rprec
    vpr = 0._rprec
    wpr = 0._rprec
  else
    !--see qpost for u/w-node interpolation
    !--convention will be to put up, vp, wp on w-nodes
    upr = 0.5_rprec * (u(1:nx, 1:ny, jz) - ubar(jz) +  &
                       u(1:nx, 1:ny, jz-1) - ubar(jz-1))
    vpr = 0.5_rprec * (v(1:nx, 1:ny, jz) - vbar(jz) +  &
                       v(1:nx, 1:ny, jz-1) - vbar(jz-1))
    wpr = w(1:nx, 1:ny, jz) - wbar(jz)
  end if
 
  upr2bar(1, jz) = sum (upr**2) / (nx * ny)
  upr2bar(2, jz) = sum (vpr**2) / (nx * ny)
  upr2bar(3, jz) = sum (wpr**2) / (nx * ny)

  stressbar(1, jz) = sum (upr * wpr) / (nx * ny) 
  stressbar(2, jz) = sum (txz(1:nx, 1:ny, jz)) / (nx * ny)
  stressbar(3, jz) = sum (stressbar(1:2, jz))

  !--energy spectra
  call spectrum (u(:, :, jz), Eozbar(:, jz))  !--not /z yet
  z = (jz - 0.5_rprec) * dz
  Eozbar(:, jz) = Eozbar(:, jz) / z

end do

!--collect current stats into nz_tot sized arrays
$if ($MPI)

  if (DEBUG) then
    write (*, *) coord, ': ubar(1) = ', ubar(1)
  end if

  recvcounts = size (ubar)
  displs = coord_of_rank * recvcounts 
  call mpi_gatherv (ubar(1), size (ubar), 
  ,                &
                    ubar_tot(1, 1), recvcounts, displs, MPI_RPREC,  &
                    rank_of_coord(0), comm, ierr)

  recvcounts = size (upr2bar)
  displs = coord_of_rank * recvcounts
  call mpi_gatherv (upr2bar(1, 1), size (upr2bar), MPI_RPREC,          &
                    upr2bar_tot(1, 1), recvcounts, displs, MPI_RPREC,  &
                    rank_of_coord(0), comm, ierr)  

  recvcounts = size (stressbar)
  displs = coord_of_rank * recvcounts
  call mpi_gatherv (stressbar(1, 1), size (stressbar), MPI_RPREC,        &
                    stressbar_tot(1, 1), recvcounts, displs, MPI_RPREC,  &
                    rank_of_coord(0), comm, ierr)

  recvcounts = size (Eozbar)
  displs = coord_of_rank * recvcounts
  call mpi_gatherv (Eozbar(1, 1), size (Eozbar), MPI_RPREC,              &
                    Eozbar_tot(1, 1, 1), recvcounts, displs, MPI_RPREC,  &
                    rank_of_coord(0), comm, ierr)

$else

  ubar_tot(1, :) = ubar

  upr2bar_tot = upr2bar

  stressbar_tot = stressbar

  Eozbar_tot(1, :, :) = Eozbar

$endif

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
  !--calculation of cumulative average stats
  ubar_avg = (n_ubar_avg * ubar_avg + ubar_tot) / (n_ubar_avg + 1)
  n_ubar_avg = n_ubar_avg + 1

  upr2bar_avg = (n_upr2bar_avg * upr2bar_avg + upr2bar_tot) /  &
                (n_upr2bar_avg + 1)
  n_upr2bar_avg = n_upr2bar_avg + 1

  stressbar_avg = (n_stressbar_avg * stressbar_avg + stressbar_tot) /  &
                  (n_stressbar_avg + 1)
  n_stressbar_avg = n_stressbar_avg + 1

  Eozbar_avg = (n_Eozbar_avg * Eozbar_avg + Eozbar_tot) / (n_Eozbar_avg + 1)
  n_Eozbar_avg = n_Eozbar_avg + 1

  !--prepare list of z-coordinates
  forall (jz=1:nz_tot-1) zu(1, jz) = (jz - 0.5_rprec) * dz

  !--prepare  header, optional

  !--write out to file
  call write_avg (fubar_avg, n_ubar_avg, zu, ubar_avg)
  call write_avg (fupr2bar_avg, n_upr2bar_avg, zu, upr2bar_avg)
  call write_avg (fstressbar_avg, n_stressbar_avg, zu, stressbar_avg)

  !--this is a bit awkward: maybe just have another routine to do it right
  Eozbar_hdr = 'zone' !--this is for tecplot... 
  kz_z(1, :) = kx(1:nx/2, 1) * zu(1, 1)
  kz_z(2, :) = zu(1, 1)
  call write_avg (fEozbar_avg, n_Eozbar_avg, kz_z, Eozbar_avg(:, :, 1),  &
                  hdr=Eozbar_hdr) 

  do jz = 2, nz_tot - 1

    kz_z(1, :) = kx(1:nx/2, 1) * zu(1, jz)
    kz_z(2, :) = zu(1, jz)
  
    call write_avg (fEozbar_avg, n_Eozbar_avg, kz_z, Eozbar_avg(:, :, jz),  &
                    hdr=Eozbar_hdr, position='append') 

  end do

end if

001 continue  !--exit cleanly

end subroutine avg_stats

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine init_avg (file_avg, n_ccol, a_avg, n_avg, leaveopn)
implicit none

character (*), intent (in) :: file_avg
integer, intent (in) :: n_ccol  !--num. coord columns: x, y, etc.

real (rprec), intent (out) :: a_avg(:, :)
integer, intent (out) :: n_avg

character (*), optional, intent (in) :: leaveopn

character (128) :: buff

logical :: exst, opn

integer :: j

real (rprec) :: z(n_ccol)

!---------------------------------------------------------------------

inquire (file=file_avg, exist=exst, opened=opn)

if (exst) then

  if (.not. opn) then
    open (1, file=file_avg)

    read (1, '(a)') buff

    if (buff(1:1) == '#') then
      read (buff(2:), *) n_avg
    else
      write (*, *) 'avg_stats: error'
      write (*, *) trim (file_avg), ' does not have expected format on line 1'
      stop  !--need to replace all stops with nice mpi exits
    end if
  end if

  !--skip data header lines here
  do
    read (1, '(a)') buff
    if (trim (buff) == trim (end_hdr_avg)) exit
  end do

  do j = 1, size (a_avg, 2)
    read (1, *) z, a_avg(:, j)  !--z is just placeholder here
  end do

  if (present (leaveopn)) then
    if (leaveopn /= 'yes') close (1)  !--case sensitive here
  else
    close (1)
  end if

else

  n_avg = 0
  a_avg = 0._rprec

end if

end subroutine init_avg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_avg (file_avg, n_avg, x, a_avg, hdr, position)
implicit none

character (*), intent (in) :: file_avg
integer, intent (in) :: n_avg
real (rprec), intent (in) :: x(:, :)  !--coord columns for x, y, etc
real (rprec), intent (in) :: a_avg(:, :)

character (*), optional, intent (in) :: hdr
character (*), optional, intent (in) :: position

character (64) :: r_fmt, fmt
character (32) :: posn

integer :: j

!---------------------------------------------------------------------

!--check sizes compatible
if (size (x, 2) /= size (a_avg, 2)) then
  write (*, *) 'write_avg: error with sizes of x, a_avg'
  stop
end if

if (present (position)) then
  posn = position
else
  posn = 'rewind'
end if

open (1, file=file_avg, position=posn)

if (trim (posn) /= 'append') then  !--case sensitive
  write (1, '(a,i0)') '# ', n_avg  !--not needed when appending
end if

if (present (hdr)) then
  !--write data header, if present
  write (1, '(a)') trim (hdr)
end if

!--write something to indicate end of header, always do this
write (1, '(a)') end_hdr_avg

!--define output format
write (r_fmt, '(2(a,i0))') 'es', precision (1._rprec) + 7,  &
                           '.', precision (1._rprec)
write (fmt, '(a,i0,3a)') '(', size (x, 1) + size (a_avg, 1),  &
                         '(1x,', trim (r_fmt), '))'

!--write to file
do j = 1, size (a_avg, 2)
  write (1, fmt) x(:, j), a_avg(:, j)
end do

close (1)

end subroutine write_avg

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

logical, parameter :: DEBUG = .false.

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
  if ( DEBUG ) write (*, *) sub // ': wrote record ', rec
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

logical, parameter :: DEBUG = .false.

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

if ( DEBUG ) write (*, *) sub // ' : read record ', rec
    
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
