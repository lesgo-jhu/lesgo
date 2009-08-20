!***************************************************************
subroutine stats_init ()
!***************************************************************
!  This subroutine allocates the memory for arrays
!  used for statistical calculations 

use param, only : dy,dz,nx,ny,nz,nsteps,coord
use stat_defs
use grid_defs
implicit none

character(120) :: cx,cy,cz
character(120) :: fname
integer :: fid, i,j,k

! ------ These values are not to be changed ------
!  Don't change from false
tavg_t%calc     = .false.
tavg_t%started  = .false.
point_t%started = .false.
domain_t%started = .false.
!  Initialize with non used integer 
point_t%xyz=-1.
! ------ These values are not to be changed ------

!  Master switch for turning on or off all statistics
!  including instantaneous recordings 
!stats_t%calc = .true.

!  Sub switches for statistics and output
!  Turns temporal averaged quantities on or off
!aver_calc = .false.

!  All nstart and nend values are based
!  on jt and not jt_total
tavg_t%calc = .true.
tavg_t%nstart = 1
tavg_t%nend = nsteps

!  Turns Reynolds stresses calculations on or off 
rs_t%calc = .false.

!  Turns instantaneous velocity recording on or off
point_t%calc = .false.
point_t%nstart = 1
point_t%nend   = nsteps
point_t%nskip = 1
point_t%nloc = 2
point_t%xyz(:,1) = (/3., 2., 0.5/)
point_t%xyz(:,2) = (/1., 2., 0.5/)

domain_t%calc = .false.
domain_t%nstart = 100
domain_t%nend   = nsteps
domain_t%nskip = 100

!  y-plane stats/data
yplane_t%calc=.false.
yplane_t%nstart = 1
yplane_t%nend   = nsteps
yplane_t%nloc     = 1
yplane_t%loc(1)  = 2.0

!  z-plane stats/data
zplane_t%calc=.false.
zplane_t%nstart = 1
zplane_t%nend   = nsteps
zplane_t%nloc   = 1 
zplane_t%loc(1)  = 0.5

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
  allocate(tavg_t%dudz(nx, ny, $lbz:nz))
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
if(yplane_t%calc) then
  allocate(yplane_t%ua(Nx,yplane_t%nloc,Nz))
  allocate(yplane_t%va(Nx,yplane_t%nloc,Nz))
  allocate(yplane_t%wa(Nx,yplane_t%nloc,Nz))
  
  yplane_t%fa = 1./(dble(yplane_t%nend - yplane_t%nstart + 1))
  yplane_t%ua = 0.
  yplane_t%va = 0.
  yplane_t%wa = 0.
  yplane_t%istart = -1
  yplane_t%ldiff = 0.
!  Not really needed
  yplane_t%coord = -1
  
!  Compute istart and ldiff
  do j=1,yplane_t%nloc
    call find_istart(y,ny,yplane_t%loc(j),yplane_t%istart(j), yplane_t%ldiff(j))
  enddo
    
endif

! Initialize information for y-planar stats/data
if(zplane_t%calc) then
  allocate(zplane_t%ua(Nx,Ny,zplane_t%nloc))
  allocate(zplane_t%va(Nx,Ny,zplane_t%nloc))
  allocate(zplane_t%wa(Nx,Ny,zplane_t%nloc))
  zplane_t%fa = 1./(dble(zplane_t%nend - zplane_t%nstart + 1))
  zplane_t%ua = 0.
  zplane_t%va = 0.
  zplane_t%wa = 0.

!  Initialize 
  zplane_t%istart = -1
  zplane_t%ldiff = 0. 
  zplane_t%coord=-1 
  
!  Compute istart and ldiff
  do k=1,zplane_t%nloc

    $if ($MPI)
    if(zplane_t%loc(k) >= z(1) .and. zplane_t%loc(k) < z(nz)) then
      zplane_t%coord(k) = coord
      call find_istart(z,nz,zplane_t%loc(k),zplane_t%istart(k), zplane_t%ldiff(k))
    endif
    $else
    zplane_t%coord(k) = 0
    call find_istart(z,nz,zplane_t%loc(k),zplane_t%istart(k), zplane_t%ldiff(k))
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
      call find_istart(x,nx,point_t%xyz(1,i),point_t%istart(i), point_t%xdiff(i))
      call find_istart(y,ny,point_t%xyz(2,i),point_t%jstart(i), point_t%ydiff(i))
      call find_istart(z,nz,point_t%xyz(3,i),point_t%kstart(i), point_t%zdiff(i))

      write(cx,'(F9.4)') point_t%xyz(1,i)
      write(cy,'(F9.4)') point_t%xyz(2,i)
      write(cz,'(F9.4)') point_t%xyz(3,i)
      write (fname,*) 'output/uvw_inst-',trim(adjustl(cx)),'-',  &
        trim(adjustl(cy)),'-', trim(adjustl(cz)),'.dat'
      fname=trim(adjustl(fname))
      fid=3000*i
      open(unit = fid,file = fname,status="unknown",position="rewind")
      write(fid,*) 'variables= "t (s)", "u", "v", "w"'
    endif
  $else
    point_t%coord(i) = 0
    call find_istart(x,nx,point_t%xyz(1,i),point_t%istart(i), point_t%xdiff(i))
    call find_istart(y,ny,point_t%xyz(2,i),point_t%jstart(i), point_t%ydiff(i))
    call find_istart(z,nz,point_t%xyz(3,i),point_t%kstart(i), point_t%zdiff(i))

    write(cx,'(F9.4)') point_t%xyz(1,i)
    write(cy,'(F9.4)') point_t%xyz(2,i)
    write(cz,'(F9.4)') point_t%xyz(3,i)
    write (fname,*) 'output/uvw_inst-',trim(adjustl(cx)),'-',  &
      trim(adjustl(cy)),'-', trim(adjustl(cz)),'.dat'
    fname=trim(adjustl(fname))
    fid=3000*i
    open(unit = fid,file = fname,status="unknown",position="rewind")
    write(fid,*) 'variables= "t (s)", "u", "v", "w"'
  $endif
  
  enddo
endif

return
end subroutine stats_init

!**********************************************************************
subroutine find_istart(x,nx,px,istart,xdiff)
!**********************************************************************
implicit none

integer, intent(IN) :: nx
double precision, dimension(nx), intent(IN) :: x
double precision, intent(IN) :: px
integer, intent(OUT) :: istart
double precision, intent(OUT) :: xdiff

integer :: i

isearch: do i=1,nx
  if(x(i) >= px) then
    istart = i-1
    xdiff = px - x(istart)
    exit isearch
  endif
enddo isearch

return
end subroutine find_istart
