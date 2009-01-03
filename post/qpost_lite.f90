program qpost_lite
implicit none

integer, parameter :: rp = kind (1.d0)

integer,parameter:: nx=128,ny=64,nz=65
integer,parameter:: nx2=3*nx/2,ny2=3*ny/2
integer,parameter:: lh=nx/2+1,ld=2*lh

character (*), parameter :: fmt = '(9(1x,es17.10))'

!--merge output files by separate MPI processes
logical, parameter :: merge_MPI = .false.
integer, parameter :: np = 1
character (*), parameter :: MPI_suffix = '.c'  !--c for coord, then a number
character (64) :: MPI_fname(0:np-1)
integer :: lbz, ubz

real (rp),parameter::pi=3.1415926535897932384626433_rp
real (rp), parameter :: BOGUS = -1234567890._rp

real (rp),parameter::L_x=2._rp, L_y= 1._rp
real (rp),parameter::z_i=1._rp,L_z=1._rp * z_i !1._rp - 1._rp / nz

logical, parameter :: use_mean_p_force = .false.
real (rp), parameter :: mean_p_force = 1._rp * z_i / L_z
                                       !--usually just z_i / L_z

logical, parameter :: write_xslices = .false.

real (rp),parameter::dz=L_z/z_i/(nz-1)
real (rp),parameter::dx=L_x/nx,dy=L_y/ny

real (rp), parameter :: z0 = 0.0001_rp  !--nondimensional already

real (rp), save, dimension(ld,ny,nz)::u,v,w,Cs_opt2_avg,Cs_opt2 !,F_LM,F_MM
real (rp), save, dimension(ld,ny,nz)::u_avg,v_avg,w_avg
real (rp), save :: phi(ld, ny, 0:nz)  !--level set function
                               !--0-level not used for non-MPI
real (rp) :: x,y,z
real (rp) :: fp

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
logical :: have_txz

integer :: nstart, nstop, nstep
integer :: jx_min, jx_max, jy_min, jy_max, jz_min, jz_max
integer :: nfiles
integer :: ist, jst

real (rp) :: U_infty_avg, U_infty
real (rp) :: u2, ww

!---------------------------------------------------------------------

u_avg=0._rp; v_avg=0._rp; w_avg=0._rp; Cs_opt2_avg=0._rp

U_infty_avg = 0._rp

if (use_mean_p_force) then
  fp = z_i / L_z
else
  fp = 0._rp
end if

!--read in level-set function, if it is there, else put 0
write (fname, '(a)') '../phi.out'

if (merge_MPI) then

  do ip = 0, np-1
  
    write (MPI_fname(ip), '(a,a,i0)') trim (fname), MPI_suffix, ip
    
    inquire (file=MPI_fname(ip), exist=exst)
    
    if (exst) then
      write (*, *) 'attempting to read phi for ip = ', ip
      open (1, file=MPI_fname(ip), form='unformatted')
      lbz = ip * (nz-1) / np  !--these will overlap
      ubz = lbz + (nz-1) / np + 1
      read (1) phi(:, :, lbz:ubz)
      close (1)
    else
      write (*, *) 'all or part of phi missing, setting phi to zero'
      phi = 0._rp
      exit
    end if

  end do

else

  inquire (file=fname, exist=exst)
  if (exst) then
    write (*, *) 'attempting to read phi'
    write (*, *) 'expected size=', size (phi, 1), size (phi, 2),  &
                 size (phi, 3) - 1

    open (1, file=fname, form='unformatted')
    read (1) phi(:, :, 1:nz)
    close (1)

    phi(:, :, 0) = BOGUS
    
  else
    write (*, *) 'phi missing, setting phi to zero'
    phi = 0._rp
  end if

end if

open (21, file='U_infty.dat')

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

!--initially assume that txz files are present until we find otherwise
have_txz = .true.

do jt=nstart,nstop,nstep

  write (fname, fnamefmt) 'vel', jt, '.out'
  !write (fname, fnamefmt) 'vel-', jt, '.out'

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
      !--problem here  with pressure: 0 layer???
      lbz = ip * (nz-1) / np + 1
      ubz = lbz + (nz-1) / np
      read (1) u(:, :, lbz:ubz),  &
               v(:, :, lbz:ubz),  &
               w(:, :, lbz:ubz),  &
               Cs_opt2(:, :, lbz:ubz),  &
               Cs_opt2(:, :, lbz:ubz),  &
               Cs_opt2(:, :, lbz:ubz),  &
               Cs_opt2(:, :, lbz:ubz)
               !p(:, :, lbz-1:ubz)
               !--first 3 Cs_opt2s eat up RHSx, RHSy, RHSz
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
    read (1) u, v, w, Cs_opt2, Cs_opt2, Cs_opt2, Cs_opt2 !, p
                      !--first 3 Cs_opt2 just to eat up RHSx, RHSy, RHSz
                      !--next 3
    !read(1) u,v,w,Cs_opt2_avg,p
    close(1)

  end if

  !--reference velocity: avg u over inflow plane
  U_infty = 0._rp
  U_infty = sum (u(1, 1:ny, 1:nz-1)) / size (u(1, 1:ny, 1:nz-1))

  write (*, *) 'U_infty = ', U_infty
  write (21, *) jt, U_infty

  U_infty_avg = U_infty_avg + U_infty

  u_avg = u_avg + u
  v_avg = v_avg + v
  w_avg = w_avg + w
  Cs_opt2_avg = Cs_opt2_avg + Cs_opt2

  write (fbase, fbasefmt) 'vel', jt
  write (6, '(a)') fbase

  if (write_xslices) then
    call write_vel_xslices (fbase, u, v, w, Cs_opt2, phi)
  else
    call write_vel (fbase, u, v, w, Cs_opt2, phi)
  end if

end do

close (21)

nfiles = (nstop - nstart + nstep) / nstep  !--watch truncation here
write (*, *) 'nfiles = ', nfiles

U_infty_avg = U_infty_avg / nfiles
write (*, *) 'U_infty_avg = ', U_infty_avg

!--normalize all averages
u_avg = u_avg / nfiles
v_avg = v_avg / nfiles
w_avg = w_avg / nfiles
Cs_opt2_avg = Cs_opt2_avg / nfiles

!--write unformatted avg file (different vars than below)
open (1, file='vel-avg.out', form='unformatted')
write (1) u_avg, v_avg, w_avg
close (1)

if (write_xslices) then
  call write_vel_xslices ('vel-avg', u_avg, v_avg, w_avg, Cs_opt2_avg, phi)
else
  call write_vel ('vel-avg', u_avg, v_avg, w_avg, Cs_opt2_avg, phi)
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_vel (fbase, u, v, w, Cs2, phi)
implicit none

character (*), intent (in) :: fbase
real (rp), intent (in) :: u(ld, ny, nz), v(ld, ny, nz), w(ld, ny, nz)
real (rp), intent (in) :: Cs2(ld, ny, nz), phi(ld, ny, 0:nz)

character (128) :: fname

integer :: jx, jy, jz

real (rp) :: x, y, z
real (rp) :: ww
real (rp) :: u2
real (rp) :: Cs2_tmp

!---------------------------------------------------------------------

!--formatted avg output
write (fname, '(a)') trim (fbase) // '.dat'
open (1, file=fname)
write(1, '(a)') 'variables="x" "y" "z" "u" "v" "w" "|u|" "Cs^2" "phi"'
write(1, *) 'zone f=point, i=', nx, ', j=', ny, ', k=', nz-1

! write velocity vectors interpolated to uvp-nodes
do jz=1,nz-1
  do jy=1,ny
    do jx=1,nx

      x=(jx-1)*dx
      y=(jy-1)*dy
      z=(jz-0.5_rp)*dz

      ww = 0.5_rp * (w(jx, jy, jz) + w(jx, jy, jz+1))
      u2 = u(jx, jy, jz)**2 + v(jx, jy, jz)**2 + ww**2

      if (jz == 1) then
        Cs2_tmp = Cs2(jx, jy, jz)
      else
        Cs2_tmp = 0.5_rp * (Cs2(jx,jy,jz)+Cs2(jx,jy,jz+1))
      end if
      
      write (1, fmt) x, y, z, u(jx, jy, jz), v(jx, jy, jz), ww,  &
                     sqrt (u2), Cs2_tmp, phi(jx, jy, jz)

    end do
  end do
end do

close(1)

end subroutine write_vel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--write each x-plane to a separate file (convenient for big MPI files)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_vel_xslices (fbase, u, v, w, Cs2, phi)
implicit none

character (*), intent (in) :: fbase
real (rp), intent (in) :: u(ld, ny, nz), v(ld, ny, nz), w(ld, ny, nz)
real (rp), intent (in) :: Cs2(ld, ny, nz), phi(ld, ny, 0:nz)

integer, parameter :: jxskip = 4
integer, parameter :: jxstart = 1

character (128) :: fname

integer :: jx, jy, jz

real (rp) :: x, y, z
real (rp) :: ww
real (rp) :: u2
real (rp) :: Cs2_tmp

!---------------------------------------------------------------------

do jx = jxstart, nx, jxskip

  x = (jx - 1) * dx

  !--formatted avg output
  write (fname, '(a,i0,a)') trim (fbase) // '.i', jx,  '.dat'
  open (1, file=fname)
  write (1, '(a,es13.6)') '# x = ', x
  write(1, '(a)') 'variables="y" "z" "u" "v" "w" "|u|" "Cs^2" "phi"'
  write(1, *) 'zone f=point, i=', ny, ', j=', nz-1

  ! write velocity vectors interpolated to uvp-nodes
  do jz = 1, nz - 1
    do jy = 1, ny

      y=(jy-1)*dy
      z=(jz-0.5_rp)*dz

      ww = 0.5_rp * (w(jx, jy, jz) + w(jx, jy, jz+1))
      u2 = u(jx, jy, jz)**2 + v(jx, jy, jz)**2 + ww**2

      if (jz == 1) then
        Cs2_tmp = Cs2(jx, jy, jz)
      else
        Cs2_tmp = 0.5_rp * (Cs2(jx,jy,jz)+Cs2(jx,jy,jz+1))
      end if

      !--fmt has 9 fields here
      write (1, fmt) y, z, u(jx, jy, jz), v(jx, jy, jz), ww,  &
                     sqrt (u2), Cs2_tmp, phi(jx, jy, jz)

    end do
  end do

  close(1)

end do

end subroutine write_vel_xslices

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end program qpost_lite
