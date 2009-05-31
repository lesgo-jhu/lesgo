!***************************************************************
subroutine stats_init ()
!***************************************************************
!  This subroutine allocates the memory for arrays
!  used for statistical calculations 

use stat_defs
use param, only : dy,dz,nx,ny,nz,nsteps
implicit none

character(15) :: ci,cj,ck
character(120) :: fname
integer :: fid, i, j,k, jy, kz

! ------ These values are not to be changed ------
!  Don't change from false
tavg_t%calc     = .false.
tavg_t%started  = .false.
ui_pnt_t%started = .false.
ui_gbl_t%started = .false.
!  Initialize with non used integer 
ui_pnt_t%ijk=-1
! ------ These values are not to be changed ------

!  Master switch for turning on or off all statistics
!  including instantaneous recordings 
stats_t%calc = .false.

!  Sub switches for statistics and output
!  Turns Reynolds stresses calculations on or off 
rs_t%calc = .false.
!  Turns instantaneous velocity recording on or off
ui_pnt_t%calc = .false.
ui_gbl_t%calc = .false.
!  Turns temporal averaged quantities on or off
avg_calc = .false.

!  All nstart and nend values are based
!  on jt and not jt_total
tavg_t%nstart = 1
tavg_t%nend = nsteps

ui_pnt_t%nstart = 1
ui_pnt_t%nend   = nsteps
ui_pnt_t%nskip = 1

ui_pnt_t%nloc = 3
ui_pnt_t%ijk(:,1) = (/ nx/2+1, ny/2+1, nz/2+1 /)
ui_pnt_t%ijk(:,2) = (/ nx/2+1, ny/2+1, 1 /)
ui_pnt_t%ijk(:,3) = (/ nx/2+1, ny/2+1, nz /)

!ui_gbl_t%nstart = ui_pnt_t%nstart
ui_gbl_t%nstart = 80000
ui_gbl_t%nend   = nsteps
ui_gbl_t%nskip = 100

!  y-plane stats/data
yplane_t%avg=.false.
yplane_t%nstart = 1
yplane_t%nend   = nsteps
yplane_t%na     = 1
yplane_t%la(1)  = 2.0

!  z-plane stats/data
zplane_t%avg=.false.
zplane_t%nstart = 1
zplane_t%nend   = nsteps
zplane_t%na     = 1
zplane_t%la(1)  = 0.5
zplane_t%la(2)	= 2.25



!  Set time summation calculations based on
!  dependants. Don't touch, depends on above
!  information
if(rs_t%calc) tavg_t%calc = .true.

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
if(tavg_t%calc) then 
  allocate(tavg_t%u(nx, ny, $lbz:nz))
  allocate(tavg_t%v(nx, ny, $lbz:nz))
  allocate(tavg_t%w(nx, ny, $lbz:nz))
  allocate(tavg_t%u2(nx, ny, $lbz:nz))
  allocate(tavg_t%v2(nx, ny, $lbz:nz))
  allocate(tavg_t%w2(nx, ny, $lbz:nz))
  allocate(tavg_t%uw(nx, ny, $lbz:nz))
  allocate(tavg_t%vw(nx, ny, $lbz:nz))
  allocate(tavg_t%uv(nx, ny, $lbz:nz))
  allocate(tavg_t%dudz(nz, ny, $lbz:nz))
  !  Initialize arrays
  tavg_t%u=0.
  tavg_t%v=0.
  tavg_t%w=0.
  tavg_t%u2=0.
  tavg_t%v2=0.
  tavg_t%w2=0.
  tavg_t%uw=0.
  tavg_t%vw=0.
  tavg_t%uv=0.
  tavg_t%dudz=0.  
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

! Initialize information for y-planar stats/data
if(yplane_t%avg) then
  allocate(yplane_t%ua(Nx,yplane_t%na,Nz))
  allocate(yplane_t%va(Nx,yplane_t%na,Nz))
  allocate(yplane_t%wa(Nx,yplane_t%na,Nz))
  
  yplane_t%fa = 1./(dble(yplane_t%nend - yplane_t%nstart + 1))
  yplane_t%ua = 0.
  yplane_t%va = 0.
  yplane_t%wa = 0.
  yplane_t%istart = -1
  yplane_t%ldiff = 0.
  
!  Compute istart and ldiff
  do jy=1,yplane_t%na
    isearch_j: do j=1,ny
	  if((j-1)*dy >= yplane_t%la(jy)) then
	    yplane_t%istart(jy) = j-1
		yplane_t%ldiff(jy) = yplane_t%la(jy) - (j-1-1)*dy
		exit isearch_j
	  endif
	enddo isearch_j
  enddo
    
endif

! Initialize information for y-planar stats/data
if(zplane_t%avg) then
  allocate(zplane_t%ua(Nx,Ny,yplane_t%na))
  allocate(zplane_t%va(Nx,Ny,yplane_t%na))
  allocate(zplane_t%wa(Nx,Ny,yplane_t%na))
  zplane_t%fa = 1./(dble(zplane_t%nend - zplane_t%nstart + 1))
  zplane_t%ua = 0.
  zplane_t%va = 0.
  zplane_t%wa = 0.
  zplane_t%istart = -1
  zplane_t%ldiff = 0.  
  
!  Compute istart and ldiff
  do kz=1,zplane_t%na
    isearch_k: do k=1,nz
	  if((k-1)*dz + dz/2. >= zplane_t%la(kz)) then
	    zplane_t%istart(kz) = k-1
		zplane_t%ldiff(kz) = zplane_t%la(kz) - ((k-1-1)*dz + dz/2.)
		exit isearch_k
	  endif
	enddo isearch_k
	write(*,*) 'zplane_t%istart(kz) = ', zplane_t%istart(kz)
  enddo  
  
endif

!  Open files for instantaneous writing
if(ui_pnt_t%calc) then
  do i=1,ui_pnt_t%nloc
    write(ci,*) ui_pnt_t%ijk(1,i)
	write(cj,*) ui_pnt_t%ijk(2,i)
	write(ck,*) ui_pnt_t%ijk(3,i)
    write (fname,*) 'output/uvw-inst-',trim(adjustl(ci)),'-',  &
      trim(adjustl(cj)),'-', trim(adjustl(ck)),'.out'
	fname=trim(adjustl(fname))
	fid=3000*i
	open(unit = fid,file = fname,status="unknown",position="rewind")
	write(fid,*) 'variables= "t (s)", "u", "v", "w"'
  enddo
endif

return
end subroutine stats_init
