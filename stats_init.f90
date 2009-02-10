!***************************************************************
subroutine stats_init ()
!***************************************************************
!  This subroutine allocates the memory for arrays
!  used for statistical calculations 

use stat_defs
use param, only : nx,ny,nz,nsteps
implicit none

character(15) :: ci,cj,ck
character(120) :: fname
integer :: fid, i

! ------ These values are not to be changed ------
!  Don't change from false
taver_t%calc = .false.
taver_t%started = .false.
ui_t%started = .false.
!  Initialize with non used integer 
ui_t%ijk=-1
! ------ These values are not to be changed ------

!  Master switch for turning on or off all statistics
!  including instantaneous recordings and final timestep
!  output
stats_t%calc = .true.

!  Sub switches for statistics and output
 !  Turns Reynolds stresses calculations on or off 
rs_t%calc = .false.
!  Turns instantaneous velocity recording on or off
ui_t%calc = .true.
ui_t%global = .true.
!  Turns temporal averaged quantities on or off
aver_calc = .false.

!  All nstart and nend values are based
!  on jt and not jt_total
taver_t%nstart = 1
taver_t%nend = nsteps

ui_t%nstart = 1
ui_t%nend = nsteps
ui_t%nloc = 1
ui_t%ijk(:,1) = (/ nx/2+1, ny/2+1, nz/2+1 /)
ui_t%ijk(:,2) = (/ nx/2+1, ny/2+1, 1 /)
ui_t%ijk(:,3) = (/ nx/2+1, ny/2+1, nz /)

ui_t%nskip = 100

!  Set time summation calculations based on
!  dependants. Don't touch, depends on above
!  information
if(rs_t%calc) taver_t%calc = .true.

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
if(taver_t%calc) then 
  allocate(taver_t%u(nx, ny, $lbz:nz))
  allocate(taver_t%v(nx, ny, $lbz:nz))
  allocate(taver_t%w(nx, ny, $lbz:nz))
  allocate(taver_t%u2(nx, ny, $lbz:nz))
  allocate(taver_t%v2(nx, ny, $lbz:nz))
  allocate(taver_t%w2(nx, ny, $lbz:nz))
  allocate(taver_t%uw(nx, ny, $lbz:nz))
  allocate(taver_t%vw(nx, ny, $lbz:nz))
  allocate(taver_t%uv(nx, ny, $lbz:nz))
  allocate(taver_t%dudz(nz, ny, $lbz:nz))
  !  Initialize arrays
  taver_t%u=0.
  taver_t%v=0.
  taver_t%w=0.
  taver_t%u2=0.
  taver_t%v2=0.
  taver_t%w2=0.
  taver_t%uw=0.
  taver_t%vw=0.
  taver_t%uv=0.
  taver_t%dudz=0.  
endif

if(rs_t%calc) then
  allocate(rs_t%up2(nx, ny, $lbz:nz))
  allocate(rs_t%vp2(nx, ny, $lbz:nz))
  allocate(rs_t%wp2(nx, ny, $lbz:nz))
  allocate(rs_t%upwp(nx, ny, $lbz:nz))
  allocate(rs_t%vpwp(nx, ny, $lbz:nz))
  allocate(rs_t%upvp(nx, ny, $lbz:nz))
  rs_t%up2=0.
  rs_t%vp2=0.
  rs_t%wp2=0.
  rs_t%upwp=0.
  rs_t%vpwp=0.
  rs_t%upvp=0.
endif

!  Open files for instantaneous writing
if(ui_t%calc) then
  do i=1,ui_t%nloc
    write(ci,*) ui_t%ijk(1,i)
	write(cj,*) ui_t%ijk(2,i)
	write(ck,*) ui_t%ijk(3,i)
    write (fname,*) 'output/inst_uvw-',trim(adjustl(ci)),'-',  &
      trim(adjustl(cj)),'-', trim(adjustl(ck)),'.out'
	fname=trim(adjustl(fname))
	fid=3000*i
	open(unit = fid,file = fname,status="unknown",position="rewind")
	write(fid,*) 'variables= "t (s)", "u", "v", "w"'
  enddo
endif

return
end subroutine stats_init
