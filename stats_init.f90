!***************************************************************
subroutine stats_init ()
!***************************************************************
!  This subroutine allocates the memory for arrays
!  used for statistical calculations 

use param, only : L_x,L_y,L_z,dx,dy,dz,nx,ny,nz,nsteps,coord
use stat_defs
use grid_defs
use functions, only : index_start
implicit none

character(120) :: cx,cy,cz
character(120) :: fname
integer :: fid, i,j,k

! ------ These values are not to be changed ------
!  Don't change from false
tsum_t%calc     = .false.
tsum_t%started  = .false.
point_t%started = .false.
domain_t%started = .false.
!  Initialize with non used integer 
point_t%xyz=-1.
! ------ These values are not to be changed ------

!  All nstart and nend values are based
!  on jt and not jt_total
tsum_t%calc = .true.
tsum_t%nstart = 100000
tsum_t%nend = nsteps

!  Turns instantaneous velocity recording on or off
point_t%calc = .false.
point_t%nstart = 50000
point_t%nend   = nsteps
point_t%nskip = 10
point_t%nloc = 2
point_t%xyz(:,1) = (/L_x/2., L_y/2., 1.5_rprec/)
point_t%xyz(:,2) = (/L_x/2., L_y/2., 2.5_rprec/)

domain_t%calc = .true.
domain_t%nstart = 1
domain_t%nend   = nsteps
domain_t%nskip = 50000

!  y-plane stats/data
yplane_t%calc   = .false.
yplane_t%nstart = 990000
yplane_t%nend   = nsteps
yplane_t%nskip  = 100
yplane_t%nloc   = 2
yplane_t%loc(1) = 1.0
yplane_t%loc(2) = 3.0

!  z-plane stats/data
zplane_t%calc   = .false.
zplane_t%nstart = 8000
zplane_t%nend   = nsteps
zplane_t%nskip  = 100
zplane_t%nloc   = 2
zplane_t%loc(1) = 0.5
zplane_t%loc(2) = 1.5

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
if(tsum_t%calc) then 
  allocate(tsum_t%u(nx, ny, $lbz:nz))
  allocate(tsum_t%v(nx, ny, $lbz:nz))
  allocate(tsum_t%w(nx, ny, $lbz:nz))
  allocate(tsum_t%u2(nx, ny, $lbz:nz))
  allocate(tsum_t%v2(nx, ny, $lbz:nz))
  allocate(tsum_t%w2(nx, ny, $lbz:nz))
  allocate(tsum_t%uw(nx, ny, $lbz:nz))
  allocate(tsum_t%vw(nx, ny, $lbz:nz))
  allocate(tsum_t%uv(nx, ny, $lbz:nz))
  allocate(tsum_t%dudz(nx, ny, $lbz:nz))
  !  Initialize arrays
  tsum_t%u=0.
  tsum_t%v=0.
  tsum_t%w=0.
  tsum_t%u2=0.
  tsum_t%v2=0.
  tsum_t%w2=0.
  tsum_t%uw=0.
  tsum_t%vw=0.
  tsum_t%uv=0.
  tsum_t%dudz=0.
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
	yplane_t%istart(j) = index_start('j',dy,yplane_t%loc(j))
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
	  
	  point_t%istart(i) = index_start('i',dx,point_t%xyz(1,i))
	  point_t%jstart(i) = index_start('j',dy,point_t%xyz(2,i))
	  point_t%kstart(i) = index_start('k',dz,point_t%xyz(3,i))
	  
	  point_t%xdiff(i) = x(point_t%istart(i)) - point_t%xyz(1,i)
	  point_t%ydiff(i) = y(point_t%jstart(i)) - point_t%xyz(2,i)
	  point_t%zdiff(i) = z(point_t%kstart(i)) - point_t%xyz(3,i)

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
	point_t%istart(i) = index_start('i',dx,point_t%xyz(1,i))
	point_t%jstart(i) = index_start('j',dy,point_t%xyz(2,i))
	point_t%kstart(i) = index_start('k',dz,point_t%xyz(3,i))
	
	point_t%xdiff(i) = x(point_t%istart(i)) - point_t%xyz(1,i)
	point_t%ydiff(i) = y(point_t%jstart(i)) - point_t%xyz(2,i)
	point_t%zdiff(i) = z(point_t%kstart(i)) - point_t%xyz(3,i)

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

