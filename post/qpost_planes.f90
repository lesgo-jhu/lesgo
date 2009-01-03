!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--see helper routines at bottom of file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module qpost_planes_param
implicit none

save
public

integer, parameter :: rp = kind (1.d0)

integer,parameter:: nx=64,ny=64,nz=65
integer,parameter:: nx2=3*nx/2,ny2=3*ny/2
integer,parameter:: lh=nx/2+1,ld=2*lh

real (rp),parameter::L_x=1._rp, L_y= 1._rp
real (rp),parameter::z_i=1._rp,L_z=1._rp * z_i !1._rp - 1._rp / nz
real (rp),parameter::dz=L_z/z_i/(nz-1)
real (rp),parameter::dx=L_x/nx,dy=L_y/ny

end module qpost_planes_param

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program qpost_planes
use qpost_planes_param
implicit none

character (*), parameter :: digfmt = 'i6.6'  !--formatting digits in filenames
character (*), parameter :: fnamefmt = '(a,' // digfmt // ',a)'
character (*), parameter :: fbasefmt = '(a,' // digfmt // ')'

!--could have user enter plane info through stdin, if desired
integer, parameter :: nxpl = 1
integer, parameter :: nypl = 1
integer, parameter :: nzpl = 1
integer, parameter :: xpl(nxpl) = (/ nx/2 /)
integer, parameter :: ypl(nypl) = (/ ny/2 /)
integer, parameter :: zpl(nzpl) = (/ 1 /)

logical, parameter :: calc_avg = .false.

!--merge output files by separate MPI processes
logical, parameter :: merge_MPI = .false.
integer, parameter :: np = 1
character (*), parameter :: MPI_suffix = '.c'  !--c for coord, then a number
character (64) :: MPI_fname(0:np-1)
integer :: lbz, ubz

real (rp),parameter::pi=3.1415926535897932384626433_rp
real (rp), parameter :: BOGUS = -1234567890._rp

logical, parameter :: use_mean_p_force = .false.
real (rp), parameter :: mean_p_force = 1._rp * z_i / L_z
                                       !--usually just z_i / L_z

real (rp), parameter :: z0 = 0.0001_rp  !--nondimensional already

character(len=128)::fname, fbase

integer :: ip
integer :: jt
!integer :: jx,jy,jz
!integer :: i, j, k
integer :: ios
integer :: nstart, nstop, nstep
integer :: jx_min, jx_max, jy_min, jy_max, jz_min, jz_max
integer :: nfiles
integer :: ist, jst

logical :: exst

real (rp), save, dimension(ld,ny,nz)::u,v,w,Cs_opt2 !,F_LM,F_MM
real (rp), save, dimension(ld,ny,nz)::u_avg,v_avg,w_avg,Cs_opt2_avg
real (rp), save :: phi(ld, ny, 0:nz)  !--level set function
                               !--0-level not used for non-MPI
real (rp) :: U_infty_avg, U_infty
real (rp) :: u2, ww
!real (rp) :: x,y,z
!real (rp) :: fp

!---------------------------------------------------------------------

u_avg=0._rp; v_avg=0._rp; w_avg=0._rp; Cs_opt2_avg=0._rp

U_infty_avg = 0._rp

!if (use_mean_p_force) then
!  fp = z_i / L_z
!else
!  fp = 0._rp
!end if

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
  write (*, '(a)') fbase

  call write_xpl ( fbase, nxpl, xpl, u, v, w, Cs_opt2, phi )
  call write_ypl ( fbase, nypl, ypl, u, v, w, Cs_opt2, phi )
  call write_zpl ( fbase, nzpl, zpl, u, v, w, Cs_opt2, phi )

end do

close (21)

if ( calc_avg ) then
  nfiles = (nstop - nstart + nstep) / nstep  !--watch truncation here
  write (*, *) 'nfiles = ', nfiles
  
  U_infty_avg = U_infty_avg / nfiles
  write (*, *) 'U_infty_avg = ', U_infty_avg
  
  !--normalize all averages
  u_avg = u_avg / nfiles
  v_avg = v_avg / nfiles
  w_avg = w_avg / nfiles
  Cs_opt2_avg = Cs_opt2_avg / nfiles

  call write_xpl ( 'vel-avg', nxpl, xpl, u_avg, v_avg, w_avg, Cs_opt2_avg, phi )
  call write_ypl ( 'vel-avg', nypl, ypl, u_avg, v_avg, w_avg, Cs_opt2_avg, phi )
  call write_zpl ( 'vel-avg', nzpl, zpl, u_avg, v_avg, w_avg, Cs_opt2_avg, phi )

end if

end program qpost_planes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--helper routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_xpl ( pre, nxpl, xpl, u, v, w, Cs2, phi )
use qpost_planes_param
implicit none

character (*), intent (in) :: pre

integer, intent (in) :: nxpl
integer, intent (in) :: xpl(nxpl)

real (rp), intent (in) :: u(ld, ny, nz), v(ld, ny, nz), w(ld, ny, nz)
real (rp), intent (in) :: Cs2(ld, ny, nz), phi(ld, ny, 0:nz)

character (128) :: fname

integer :: i
integer :: d

!---------------------------------------------------------------------

!write (*, *) 'in write_xpl'

d = 1  !--x-planes
do i = 1, nxpl
    write ( fname, '(a,i0,a)' )  trim (pre) // '.xpl', xpl(i), '.dat'
    !write (*, *) fname
    call write_pl ( fname, d, xpl(i), u, v, w, Cs2, phi ) 
end do

end subroutine write_xpl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_ypl ( pre, nypl, ypl, u, v, w, Cs2, phi )
use qpost_planes_param
implicit none

character (*), intent (in) :: pre

integer, intent (in) :: nypl
integer, intent (in) :: ypl(nypl)

real (rp), intent (in) :: u(ld, ny, nz), v(ld, ny, nz), w(ld, ny, nz)
real (rp), intent (in) :: Cs2(ld, ny, nz), phi(ld, ny, 0:nz)

character (128) :: fname

integer :: i
integer :: d

!---------------------------------------------------------------------

d = 2  !--y-planes
do i = 1, nypl
    write ( fname, '(a,i0,a)' )  trim (pre) // '.ypl', ypl(i), '.dat'
    call write_pl ( fname, d, ypl(i), u, v, w, Cs2, phi ) 
end do

end subroutine write_ypl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_zpl ( pre, nzpl, zpl, u, v, w, Cs2, phi )
use qpost_planes_param
implicit none

character (*), intent (in) :: pre

integer, intent (in) :: nzpl
integer, intent (in) :: zpl(nzpl)

real (rp), intent (in) :: u(ld, ny, nz), v(ld, ny, nz), w(ld, ny, nz)
real (rp), intent (in) :: Cs2(ld, ny, nz), phi(ld, ny, 0:nz)

character (128) :: fname

integer :: i
integer :: d

!---------------------------------------------------------------------

d = 3  !--y-planes
do i = 1, nzpl
    write ( fname, '(a,i0,a)' )  trim (pre) // '.zpl', zpl(i), '.dat'
    call write_pl ( fname, d, zpl(i), u, v, w, Cs2, phi ) 
end do

end subroutine write_zpl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_pl ( fname, d, pl, u, v, w, Cs2, phi )
use qpost_planes_param
implicit none

character (*), intent (in) :: fname

integer, intent (in) :: d, pl

real (rp), intent (in) :: u(ld, ny, nz), v(ld, ny, nz), w(ld, ny, nz)
real (rp), intent (in) :: Cs2(ld, ny, nz), phi(ld, ny, 0:nz)

character (*), parameter :: fmt = '(9(es13.6,1x))'

integer, parameter :: lun = 1

integer :: i, j, k
integer :: imin, imax
integer :: jmin, jmax
integer :: kmin, kmax

logical :: exst, opn

real (rp) :: x, y, z
real (rp) :: ww
real (rp) :: u2
real (rp) :: Cs2_tmp

!---------------------------------------------------------------------

if ( d < 1 .or. d > 3 ) stop "write_pl: d out of range"

if ( d /= 1 ) then
    imin = 1
    imax = nx
else
    if ( pl < 1 .or. pl > nx ) stop "write_pl: pl out of range"
    imin = pl
    imax = pl
end if

if ( d /= 2 ) then
    jmin = 1
    jmax = ny
else
    if ( pl < 1 .or. pl > ny ) stop "write_pl: pl out of range"
    jmin = pl
    jmax = pl
end if

if ( d /= 3 ) then
    kmin = 1
    kmax = nz - 1
else
    if ( pl < 1 .or. pl > nz - 1 ) stop "write_pl: pl out of range"
    kmin = pl
    kmax = pl
end if

inquire ( lun, exist=exst, opened=opn )
if ( opn .or. .not. exst ) stop "write_pl: problem with lun"

open ( lun, file=fname, action='write', form='formatted', position='rewind' )
write ( lun, '(a)' ) 'variables="x" "y" "z" "u" "v" "w" "|u|" "Cs^2" "phi"'
write ( lun, '(3(a,i0))' ) 'zone f=point, i=', imax - imin + 1,  &
                           ', j=',             jmax - jmin + 1,  &
                           ', k=',             kmax - kmin + 1

do k = kmin, kmax
    do j = jmin, jmax
        do i = imin, imax

            x = (i - 1) * dx
            y = (j - 1) * dy
            z = (k - 0.5_rp) * dz

            ww = 0.5_rp * ( w(i, j, k) + w(i, j, k + 1) )
            u2 = u(i, j, k)**2 + v(i, j, k)**2 + ww**2

            if (k == 1) then
                Cs2_tmp = Cs2(i, j, k)
            else
                Cs2_tmp = 0.5_rp * ( Cs2(i, j, k) + Cs2(i, j, k + 1) )
            end if

            write (lun, fmt) x, y, z, u(i, j, k), v(i, j, k), ww,  &
                             sqrt (u2), Cs2_tmp, phi(i, j, k)

        end do
    end do
end do

close ( lun )

end subroutine write_pl
