program qpost_f
implicit none

integer, parameter :: rp = kind (1.d0)

integer,parameter:: nx=64,ny=64,nz=65
integer,parameter:: nx2=3*nx/2,ny2=3*ny/2
integer,parameter:: lh=nx/2+1,ld=2*lh

character (*), parameter :: fmt = '(7(1x,es17.10))'

!--merge output files by separate MPI processes
logical, parameter :: merge_MPI = .false.
integer, parameter :: np = 1
character (*), parameter :: MPI_suffix = '.c'  !--c for coord, then a number
character (64) :: MPI_fname(0:np-1)
integer :: lbz, ubz

real (rp),parameter::pi=3.1415926535897932384626433_rp
real (rp), parameter :: BOGUS = -1234567890._rp

real (rp),parameter::L_x=1._rp, L_y= 1._rp
real (rp),parameter::z_i=1._rp,L_z=1._rp * z_i !1._rp - 1._rp / nz

logical, parameter :: write_xslices = .true.

real (rp),parameter::dz=L_z/z_i/(nz-1)
real (rp),parameter::dx=L_x/nx,dy=L_y/ny

real (rp), save, dimension(ld,ny,nz)::fx, fy, fz
real (rp), save, dimension(ld,ny,nz)::fx_avg, fy_avg, fz_avg
!real (rp), save :: phi(ld, ny, 0:nz)  !--level set function
                               !--0-level not used for non-MPI
real (rp) :: x,y,z

!character (32) :: fmt

integer :: ip
integer :: jx,jy,jz,jt
integer :: i, j, k
integer :: ios

character (*), parameter :: digfmt = 'i6.6'  !--formatting digits in filenames
character (*), parameter :: fnamefmt = '(a,' // digfmt // ',a)'
character (*), parameter :: fbasefmt = '(a,' // digfmt // ')'

character(len=128)::fname, fbase
logical :: exst

integer :: nstart, nstop, nstep
integer :: jx_min, jx_max, jy_min, jy_max, jz_min, jz_max
integer :: nfiles
integer :: ist, jst

!---------------------------------------------------------------------

fx_avg = 0._rp
fy_avg = 0._rp
fz_avg = 0._rp

!--read in level-set function, if it is there, else put 0
!write (fname, '(a)') '../phi.out'
!
!if (merge_MPI) then
!
!  do ip = 0, np-1
!  
!    write (MPI_fname(ip), '(a,a,i0)') trim (fname), MPI_suffix, ip
!    
!    inquire (file=MPI_fname(ip), exist=exst)
!    
!    if (exst) then
!      write (*, *) 'attempting to read phi for ip = ', ip
!      open (1, file=MPI_fname(ip), form='unformatted')
!      lbz = ip * (nz-1) / np  !--these will overlap
!      ubz = lbz + (nz-1) / np + 1
!      read (1) phi(:, :, lbz:ubz)
!      close (1)
!    else
!      write (*, *) 'all or part of phi missing, setting phi to zero'
!      phi = 0._rp
!      exit
!    end if
!
!  end do
!
!else
!
!  inquire (file=fname, exist=exst)
!  if (exst) then
!    write (*, *) 'attempting to read phi'
!    write (*, *) 'expected size=', size (phi, 1), size (phi, 2),  &
!                 size (phi, 3) - 1
!
!    open (1, file=fname, form='unformatted')
!    read (1) phi(:, :, 1:nz)
!    close (1)
!
!    phi(:, :, 0) = BOGUS
!    
!  else
!    write (*, *) 'phi missing, setting phi to zero'
!    phi = 0._rp
!  end if
!
!end if

!--read nstart, nstop, nstep from command line
!--best usage would probably be e.g. $> echo "1000 2000 20" | ./qpost  
!--or $> echo "1000 2000 20" > file; ./qpost < file
write (*, *) 'Enter nstart, nstop, nstep:'

read (*, *) nstart, nstop, nstep

write (*, *) 'nstart = ', nstart
write (*, *) 'nstop = ', nstop
write (*, *) 'nstep = ', nstep

!--do some basic checks on nstart, nstop, nskip
if (nstart > nstop) then
  write (*, *) 'nstart > nstop, nothing to do'
  stop
end if

do jt=nstart,nstop,nstep

  write (fname, fnamefmt) 'fxyz.', jt, '.dat'

  if (merge_MPI) then

    do ip = 0, np-1
      write (MPI_fname(ip), '(a,a,i0)') trim (fname), MPI_suffix, ip
    end do

    do ip = 0, np-1
      inquire (file=MPI_fname(ip), exist=exst)
      if (.not. exst) then
        write (*, '(1x,a)') 'file ' // trim (MPI_fname(ip)) //  &
                            ' does not exist'
        write (*, '(1x,a)') 'stopping here (no averaging performed)'
        stop
      end if

      open (1, file=MPI_fname(ip), form='unformatted')
      !--note the small overlap: gives us an opportunity to check the
      !  overlapping is OK
      lbz = ip * (nz-1) / np + 1
      ubz = lbz + (nz-1) / np
      read (1) fx(:, :, lbz:ubz), fy(:, :, lbz:ubz), fz(:, :, lbz:ubz)
      close (1)

    end do

  else

    inquire (file=fname, exist=exst)
    if (.not. exst) then
      write (*, '(1x,a)') 'file ' // trim (fname) // ' does not exist'
      write (*, '(1x,a)') 'stopping here (no averaging performed)'
      stop
    end if

    open(1,file=fname,form='unformatted')
    read (1) fx, fy, fz
    close(1)

  end if

  fx_avg = fx_avg + fx
  fy_avg = fy_avg + fy
  fz_avg = fz_avg + fz

  write (fbase, fbasefmt) 'fxyz.', jt
  write (6, '(a)') fbase

  if (write_xslices) then
    call write_vel_xslices (fbase, fx, fy, fz)
  else
    call write_vel (fbase, fx, fy, fz)
  end if

end do

close (21)

nfiles = (nstop - nstart + nstep) / nstep  !--watch truncation here
write (*, *) 'nfiles = ', nfiles

!--normalize all averages
fx_avg = fx_avg / nfiles
fy_avg = fy_avg / nfiles
fz_avg = fz_avg / nfiles

if (write_xslices) then
  call write_vel_xslices ('fxyz.avg', fx_avg, fy_avg, fz_avg)
else
  call write_vel ('fxyz.avg', fx_avg, fy_avg, fz_avg)
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_vel (fbase, u, v, w)
implicit none

character (*), intent (in) :: fbase
real (rp), intent (in) :: u(ld, ny, nz), v(ld, ny, nz), w(ld, ny, nz)
!real (rp), intent (in) :: phi(ld, ny, 0:nz)

character (128) :: fname

integer :: jx, jy, jz

real (rp) :: x, y, z
real (rp) :: ww

!---------------------------------------------------------------------

!--formatted avg output
write (fname, '(a)') trim (fbase) // '.dat'
open (1, file=fname)
write(1, '(a)') 'variables="x" "y" "z" "fx" "fy" "fz"'
write(1, *) 'zone f=point, i=', nx, ', j=', ny, ', k=', nz-1

do jz=1,nz-1
  do jy=1,ny
    do jx=1,nx

      x=(jx-1)*dx
      y=(jy-1)*dy
      z=(jz-0.5_rp)*dz

      ww = 0.5_rp * (w(jx, jy, jz) + w(jx, jy, jz+1))

!      write (1, fmt) x, y, z, u(jx, jy, jz), v(jx, jy, jz), ww,  &
!                     phi(jx, jy, jz)

      write (1, fmt) x, y, z, u(jx, jy, jz), v(jx, jy, jz), ww

    end do
  end do
end do

close(1)

end subroutine write_vel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--write each x-plane to a separate file (convenient for big MPI files)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_vel_xslices (fbase, u, v, w)
implicit none

character (*), intent (in) :: fbase
real (rp), intent (in) :: u(ld, ny, nz), v(ld, ny, nz), w(ld, ny, nz)
!real (rp), intent (in) :: phi(ld, ny, 0:nz)

integer, parameter :: jxskip = nx / 4
integer, parameter :: jxstart = 1

character (128) :: fname

integer :: jx, jy, jz

real (rp) :: x, y, z
real (rp) :: ww

!---------------------------------------------------------------------

do jx = jxstart, nx, jxskip

  x = (jx - 1) * dx

  !--formatted avg output
  write (fname, '(a,i0,a)') trim (fbase) // '.i', jx,  '.dat'
  open (1, file=fname)
  write (1, '(a,es13.6)') '# x = ', x
  write(1, '(a)') 'variables="y" "z" "fx" "fy" "fz"'
  write(1, *) 'zone f=point, i=', ny, ', j=', nz-1

  ! write velocity vectors interpolated to uvp-nodes
  do jz = 1, nz - 1
    do jy = 1, ny

      y=(jy-1)*dy
      z=(jz-0.5_rp)*dz

      ww = 0.5_rp * (w(jx, jy, jz) + w(jx, jy, jz+1))

      !--watch number of fields in format
      !write (1, fmt) y, z, u(jx, jy, jz), v(jx, jy, jz), ww,  &
      !               phi(jx, jy, jz)
      write (1, fmt) y, z, u(jx, jy, jz), v(jx, jy, jz), ww

    end do
  end do

  close(1)

end do

end subroutine write_vel_xslices

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end program qpost_f
