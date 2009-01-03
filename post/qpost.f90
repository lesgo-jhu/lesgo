program qpost
use types
implicit none

integer,parameter:: nx=64,ny=64,nz=65
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

real (rp),parameter::L_x=1._rp, L_y= 1._rp
real (rp),parameter::z_i=1._rp,L_z=1._rp * z_i !1._rp - 1._rp / nz

logical, parameter :: use_mean_p_force = .false.
real (rp), parameter :: mean_p_force = 1._rp * z_i / L_z
                                       !--usually just z_i / L_z

logical, parameter :: write_xslices = .false.

logical, parameter :: plot_cylinder_stuff = .true.
real (rp), parameter :: diam = 0.0625_rp
real (rp), parameter :: xoffset = 0.3125_rp
real (rp), parameter :: yoffset = 0.5_rp
real (rp), parameter :: xstation = 1._rp  !--in cylinder-centered coords!
real (rp), parameter :: ystation = 0._rp  !--e.g. xstation = 0 == center
                                          !       xstation = 0.5 == surface
                                          !--normalized by diam

real (rp),parameter::dz=L_z/z_i/(nz-1)
real (rp),parameter::dx=L_x/nx,dy=L_y/ny

real (rp), parameter :: z0 = 0.0001_rp  !--nondimensional already

real (rp), dimension(ld,ny,nz)::u,v,w,Cs_opt2_avg,Cs_opt2 !,F_LM,F_MM
real (rp), dimension(ld,ny,nz)::u_avg,v_avg,w_avg,Cs_opt2_avg_avg
real (rp), dimension(ld, ny, nz) :: u2_avg, v2_avg, w2_avg  !--for rms calcs
real (rp), dimension(ld, ny, nz) :: urms, vrms, wrms
real (rp) :: phi(ld, ny, 0:nz)  !--level set function
real (rp) :: x,y,z
real (rp) :: fp

!character (32) :: fmt

integer :: ip
integer :: jx,jy,jz,jt
integer :: i, j, k
integer :: n_bldg
integer, allocatable :: bldg_pts(:,:)
integer :: px, py, lx, ly, lz
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

real (rp) :: p_infty, p_tmp
real (rp) :: U_infty_avg, U_infty
real (rp) :: u2, ww
real (rp) :: ubar(nz), vbar(nz), wbar(nz)
real (rp) :: ubar_avg(nz)!, vbar_avg(nz), wbar_avg(nz)
real (rp) :: up2(nz), vp2(nz), wp2(nz)
real (rp) :: up2_avg(nz), vp2_avg(nz), wp2_avg(nz)
real (rp) :: upwp(nz)
real (rp) :: upwp_avg(nz)
real (rp) :: txz_3d(ld, ny, nz), txz(nz)
real (rp) :: txz_avg(nz)
real (rp) :: Eu(nx/2+1, nz), Eu_avg(nx/2+1, nz)
real (rp) :: Csbar_avg(nz)
real (rp) :: ubarown, vbarown, wbarown
real (rp) :: up2barown_avg(nz), vp2barown_avg(nz), wp2barown_avg(nz)

!---------------------------------------------------------------------

u_avg=0._rp; v_avg=0._rp; w_avg=0._rp; Cs_opt2_avg_avg=0._rp

u2_avg = 0._rp
v2_avg = 0._rp
w2_avg = 0._rp

up2_avg = 0._rp; vp2_avg = 0._rp; wp2_avg = 0._rp
upwp_avg = 0._rp
txz_avg = 0._rp
ubar_avg = 0._rp
Eu_avg = 0._rp

Csbar_avg = 0._rp

up2barown_avg = 0._rp
vp2barown_avg = 0._rp
wp2barown_avg = 0._rp

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
    write (*, *) 'expected size=', size (phi, 1), size (phi, 2), size (phi, 3)

    open (1, file=fname, form='unformatted')
    read (1) phi(:, :, 1:nz)
    close (1)
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
  !Cs_opt2_avg_avg=Cs_opt2_avg_avg+Cs_opt2_avg
  !p_avg = p_avg + p

  u2_avg = u2_avg + u**2
  v2_avg = v2_avg + v**2
  w2_avg = w2_avg + w**2

  write (fbase, fbasefmt) 'vel', jt
  write (6, '(a)') fbase

  call write_vel (fbase, u, v, w, Cs_opt2, phi)
  if (write_xslices) call write_vel_xslices (fbase, u, v, w, Cs_opt2, phi)

  !--separate calculations of u'^2, v'^2, w'^2 on their own nodes
  !--first calc avg
  do jz = 1, nz
    ubarown = sum (u(1:nx, 1:ny, jz)) / (nx * ny)
    vbarown = sum (v(1:nx, 1:ny, jz)) / (nx * ny)
    wbarown = sum (w(1:nx, 1:ny, jz)) / (nx * ny)

    up2barown_avg(jz) = up2barown_avg(jz) +                                  &
                        sum ((u(1:nx, 1:ny, jz) - ubarown)**2) / (nx * ny)
   
    vp2barown_avg(jz) = vp2barown_avg(jz) +                                  &
                        sum ((v(1:nx, 1:ny, jz) - vbarown)**2) / (nx * ny)

    wp2barown_avg(jz) = wp2barown_avg(jz) +                                  &
                    sum ((w(1:nx, 1:ny, jz) - wbarown)**2) / (nx * ny)

  end do

  !--put everything on w-nodes here for stress balance
  do jz = 2, nz
    
    ubar(jz) = sum ( 0.5_rp * (u(1:nx, 1:ny, jz) + u(1:nx, 1:ny, jz-1)) )  &
               / (nx * ny)
    vbar(jz) = sum ( 0.5_rp * (v(1:nx, 1:ny, jz) + v(1:nx, 1:ny, jz-1)) )  &
               / (nx * ny)
    wbar(jz) = sum (w(1:nx, 1:ny, jz)) / (nx * ny)

    up2(jz) = sum ( (0.5_rp * (u(1:nx, 1:ny, jz) + u(1:nx, 1:ny, jz-1))  &
                     - ubar(jz))**2 ) / (nx * ny)
    vp2(jz) = sum ( (0.5_rp * (v(1:nx, 1:ny, jz) + v(1:nx, 1:ny, jz-1))  &
                     - vbar(jz))**2 ) / (nx * ny)
    wp2(jz) = sum ( (w(1:nx, 1:ny, jz) - wbar(jz))**2 ) / (nx * ny)

    upwp(jz) = sum ( (0.5_rp * (u(1:nx, 1:ny, jz) + u(1:nx, 1:ny, jz-1))  &
                      - ubar(jz)) *                                       &
                     (w(1:nx, 1:ny, jz) - wbar(jz)) ) / (nx * ny)

  end do

  !--at wall, these are zero
  wp2(1) = 0._rp
  upwp(1) = 0._rp

  !--at wall, these are unknown (so put zero)
  up2(1) = 0._rp
  vp2(1) = 0._rp

  up2_avg = up2_avg + up2
  vp2_avg = vp2_avg + vp2
  wp2_avg = wp2_avg + wp2

  upwp_avg = upwp_avg + upwp

  if (have_txz) then
    !--txz for stress balance
    !--if one file is missing then all the txz stuff is skipped
    write (fname, fnamefmt) 'txz', jt, '.out'

    if (merge_MPI) then

      do ip = 0, np-1
    
        write (MPI_fname(ip), '(a,a,i0)') trim (fname), MPI_suffix, ip

        inquire (file=MPI_fname(ip), exist=exst)

        if (exst) then
          open(1, file=MPI_fname(ip), form='unformatted')
          read (1) txz_3d(:, :, ip*nz/np + 1 : (ip+1)*nz/np)
          close(1)
        else
          write (*, *) 'warning: ' // trim (MPI_fname(ip)) //  &
                       ' not found, skipping txz'
          have_txz = .false.
        end if

      end do

    else

      inquire (file=fname, exist=exst)

      if (exst) then
        open (1, file=fname, form='unformatted')
        read (1) txz_3d
        close(1)
      else
        write (*, *) 'warning: ' // trim (fname) // ' not found, skipping txz'
        have_txz = .false.
      end if

    end if

    do jz = 1, nz
      txz(jz) = sum (txz_3d(1:nx, 1:ny, jz)) / (nx * ny)
    end do

    txz_avg = txz_avg + txz

  end if

  !--recalculate ubar for u-profile (on u nodes this time)
  ubar_avg = ubar_avg + sum (sum (u(1:nx, 1:ny, :), dim=1), dim=1) / (nx*ny)

  !--spectra
  do jz = 1, nz
    call spectrum (nx, ny, u(1:nx, 1:ny, jz), Eu(:, jz))
    Eu_avg(:, jz) = Eu_avg(:, jz) + Eu(:, jz)
  end do

  !--Cs
  do jz = 1, nz
    Csbar_avg(jz) = Csbar_avg(jz) +                                  &
                    sum ( sqrt (Cs_opt2(1:nx, 1:ny, jz)) ) / (nx * ny)
  end do
  
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

u2_avg = u2_avg / nfiles
v2_avg = v2_avg / nfiles
w2_avg = w2_avg / nfiles

up2_avg = up2_avg / nfiles
vp2_avg = vp2_avg / nfiles
wp2_avg = wp2_avg / nfiles

upwp_avg = upwp_avg / nfiles

txz_avg = txz_avg / nfiles  !--not needed if .not. have_txz

ubar_avg = ubar_avg / nfiles

Eu_avg = Eu_avg / nfiles

Csbar_avg = Csbar_avg / nfiles

up2barown_avg = up2barown_avg / nfiles
vp2barown_avg = vp2barown_avg / nfiles
wp2barown_avg = wp2barown_avg / nfiles

!--write unformatted avg file (different vars than below)
open (1, file='vel-avg.out', form='unformatted')
write (1) u_avg, v_avg, w_avg
close (1)

call write_vel ('vel-avg', u_avg, v_avg, w_avg, Cs_opt2_avg, phi)

if (write_xslices) then
  call write_vel_xslices ('vel-avg', u_avg, v_avg, w_avg, Cs_opt2_avg, phi)
end if

!--spatial rms distribution
urms = sqrt (u2_avg - u_avg**2)
vrms = sqrt (v2_avg - v_avg**2)
wrms = sqrt (w2_avg - w_avg**2)

open (1, file='urms.dat')
write (1, *) 'variables = "x" "y" "z", "urms"'
write (1, *) 'zone, f=point, i =', nx, ',j=', ny, ',k=', nz-1

do k = 1, nz-1

  z = (k - 0.5_rp) * dz
  
  do j = 1, ny
  
    y = (j-1) * dy
    
    do i = 1, nx

      x = (i-1) * dx

      write (1, '(4(es12.5,1x))') x, y, z, urms(i, j, k)
      
    end do
  end do
end do

close (1)

open (1, file='vrms.dat')
write (1, *) 'variables = "x" "y" "z", "vrms"'
write (1, *) 'zone, f=point, i =', nx, ',j=', ny, ',k=', nz-1

do k = 1, nz-1

  z = (k - 0.5_rp) * dz
  
  do j = 1, ny
  
    y = (j-1) * dy
    
    do i = 1, nx

      x = (i-1) * dx

      write (1, '(4(es12.5,1x))') x, y, z, vrms(i, j, k)
      
    end do
  end do
end do

close (1)

open (1, file='wrms.dat')
write (1, *) 'variables = "x" "y" "z", "wrms"'
write (1, *) 'zone, f=point, i =', nx, ',j=', ny, ',k=', nz-1

do k = 1, nz-1

  z = (k-1) * dz
  
  do j = 1, ny
  
    y = (j-1) * dy
    
    do i = 1, nx

      x = (i-1) * dx

      write (1, '(4(es12.5,1x))') x, y, z, wrms(i, j, k)
      
    end do
  end do
end do

close (1)

!--up2 profiles
open (1, file='up2-profile.dat')
write (1, *) 'variables = "z" "up2" "vp2" "wp2"'
write (1, *) 'zone, f=point, i=', nz

do jz = 1, nz
  write (1, *) (jz-0.5_rp)*dz, up2_avg(jz), vp2_avg(jz), wp2_avg(jz)
end do

close(1)

!--stress profile
if (have_txz) then

  open (1, file='stress-profile.dat')
  write (1, *) 'variables = "z" "total" "upwp" "txz"'
  write (1, *) 'zone, f=point, i=', nz

  do jz = 1, nz
    write (1, *) (jz-1)*dz, txz_avg(jz) + upwp_avg(jz), upwp_avg(jz),  &
                 txz_avg(jz)
  end do

  close(1)

end if

!--mean velocity profile
open (1, file='u-profile.dat')
write (1, *) 'variables = "z/H" "u/u<sub>*</sub>" "log-law"'
write (1, *) 'zone, f=point, i=', nz

do jz = 1, nz

  z = (jz-0.5_rp) * dz  !--actually this is (z/z_i)...

  write (1, *) z, ubar_avg(jz), 2.5_rp * log (z/z0)
  
end do

close (1)

!--energy spectra
open (1, file='E-spectra.dat')
write (1, '(a)') 'variables = "k<sub>1</sub>z" "z/H" ' //                  &
             '"E<sub>11</sub>(k<sub>1</sub>,z)/u<sub>*</sub><sup>2</sup>/z"'
write (1, *) 'zone, f=point, i=', nx/2, ', j=', nz

do jz = 1, nz
  z = (jz-0.5_rp)*dz
  do jx = 1, nx/2  !--assume Nyquist nx/2+1 is zero
    write (1, *) 2.*pi/L_x * (jx-1) * z, z, Eu_avg(jx, jz)/z
  end do
end do

close (1)

!--Cs
open (1, file='Cs.dat')
write (1, '(a)') 'variables = "Cs" "z/H"'
write (1, *) 'zone, f=point, i=', nz - 1

do jz = 1, nz - 1

  if (jz == 1) then
    z = (jz - 0.5_rp) * dz
  else
    z = (jz - 1) * dz
  end if
  
  write (1, *) Csbar_avg(jz), z / L_z 

end do

close (1)

!--up2 on own nodes
open (1, file='up2own.dat')
write (1, '(a)') 'variables = "z" "up2"'
write (1, *) 'zone, f=point, i=', nz - 1

do jz = 1, nz-1

  z = (jz - 0.5_rp) * dz
  write (1, *) z / L_z, up2barown_avg(jz)

end do

close (1)

!--vp2 on own nodes
open (1, file='vp2own.dat')
write (1, '(a)') 'variables = "z" "vp2"'
write (1, *) 'zone, f=point, i=', nz - 1

do jz = 1, nz-1

  z = (jz - 0.5_rp) * dz
  write (1, *) z / L_z, vp2barown_avg(jz)

end do

close (1)

!--wp2 on own nodes
open (1, file='wp2own.dat')
write (1, '(a)') 'variables = "z" "wp2"'
write (1, *) 'zone, f=point, i=', nz - 1

do jz = 1, nz-1

  z = (jz - 1) * dz
  write (1, *) z / L_z, wp2barown_avg(jz)

end do

close (1)

!--calculate (kz/u_*) (dUdz)
!--MAY NEED to use ustar = sqrt (L_z * mean_p_force)
open (1, file='dUdz-norm.dat')
write (1, *) 'variables = "(kz/u*)(dU/dz)" "z/H"'
write (1, *) 'zone, f=point, i=', nz - 1

do jz = 2, nz

  z = (jz - 1) *dz
  write (1, *) (0.4 * z / 1.) * (ubar_avg(jz) - ubar_avg(jz - 1)) / dz, z/L_z

end do

close (1)

!--writes x-profile at ystation and y-profile at xstation
!--performs interpolation
if (plot_cylinder_stuff) then

  call write_xyslices (u_avg, 'u_avg')
  call write_xyslices (v_avg, 'v_avg')
  call write_xyslices (urms, 'urms')
  call write_xyslices (vrms, 'vrms')

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
subroutine write_xyslices (var, tag)
  implicit none

  real (rp), intent (in) :: var(ld, ny, nz)
  character (*), intent (in) :: tag

  character (128) :: xfile, yfile

  integer :: jxst, jyst

  real (rp) :: xst, yst
  real (rp) :: ust, ust1, wgt, uinterp

  !-------------------------------------------------------------------

  xfile = trim (tag) // '--x-profile.dat'
  yfile = trim (tag) // '--y-profile.dat'

  !--x-profile

  yst = diam * ystation + yoffset  !--"dimensional" coordinate
  jyst = floor (yst / dy + 1._rp)
  wgt = ((jyst * dy) - yst) / dy  !--interpolation weight

  if (jyst == ny) then
    write (*, *) 'error: ystation is outside of domain'
    stop
  end if

  open (1, file=xfile)

  do jx = 1, nx
  
    x = (jx - 1) * dx
    x = (x - xoffset) / diam
    
    ust = sum (var(jx, jyst, :)) / nz
    ust1 = sum (var(jx, jyst+1, :)) / nz

    uinterp = wgt * ust + (1._rp - wgt) * ust1
    write (1, '(2(1x,es12.5))') x, uinterp

  end do
  
  close (1)

  !--y-profile

  xst = diam * xstation + xoffset
  jxst = floor (xst / dx + 1._rp)
  wgt = ((jxst * dx) - xst) / dx  !--interpolation weight

  if (jxst == nx) then
    write (*, *) 'error: xstation is outside of domain'
    stop
  end if

  open (1, file=yfile)

  do jy = 1, ny

    y = (jy - 1) * dy
    y = (y - yoffset) / diam

    ust = sum (var(jxst, jy, :)) / nz
    ust1 = sum (var(jxst+1, jy, :)) / nz

    uinterp = wgt * ust + (1._rp - wgt) * ust1
    write (1, '(2(1x,es12.5))') y, ust

  end do
  
  close (1)

  end subroutine write_xyslices

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end program qpost

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine spectrum (nx, ny, u, spec)
use types
implicit none

include 'fftw_f77.i'

integer, intent (in) :: nx, ny
real (rp),dimension(nx,ny),intent(in)::u
real (rp),dimension(nx/2+1),intent(out)::spec

integer :: jy, k
integer*8, save :: plan

logical, save :: init = .false.

real (rp),dimension(nx)::vel_r,vel_c

!---------------------------------------------------------------------

write (*, *) 'spectrum disabled for 64-bit mode experiment'

!if (.not. init) then
!  call rfftw_f77_create_plan(plan,nx,FFTW_REAL_TO_COMPLEX,FFTW_ESTIMATE)
!  init = .true.
!end if
!
!! initialize
!spec(:)=0._rp
!do jy=1,ny
!   vel_r(:)= u(:,jy) / nx
!! check this normaliztion-part of forward
!! call the fft
!   call rfftw_f77_one(plan,vel_r,vel_c)
!! compute magnitudes
!! the 0.5 is the 1/2, all others are taken care of! (except maybe Nyquist)
!   spec(1)=spec(1)+0.5_rp*vel_c(1)*vel_c(1)
!   do k=2,nx/2
!      spec(k)=spec(k)+vel_c(k)*vel_c(k)+vel_c(nx+2-k)*vel_c(nx+2-k)
!!        print *,'k,vel,spec',k,vel_c(k),spec(k)
!   end do
!   spec(nx/2+1)=spec(nx/2+1)+vel_c(nx/2+1)*vel_c(nx/2+1)
!end do
!spec(:)=spec(:) / ny ! for average over Ny

end subroutine spectrum
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
