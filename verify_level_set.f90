!!
!!  Copyright (C) 2009-2013  Johns Hopkins University
!!
!!  This file is part of lesgo.
!!
!!  lesgo is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License as published by
!!  the Free Software Foundation, either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  lesgo is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with lesgo.  If not, see <http://www.gnu.org/licenses/>.
!!

!--driver program in this file, after the module
module verify_level_set
use types, only : rp => rprec
implicit none

save
public

character (*), parameter :: mod_name = 'verify_level_set'

integer, parameter :: nd = 3

integer, parameter :: nxc(nd) = (/ 64, 32, 17 /)
integer, parameter :: nxf(nd) = (/ 96, 48, 25 /)

real (rp), parameter :: Lx(nd) = (/ 2._rp, 1._rp, 1._rp /)
                        !--Lx should be same for coarse, fine otherwise
                        !  why are we trying to compare them?

real (rp), parameter :: dxc(nd) = (/ Lx(1) / nxc(1),       &
                                     Lx(2) / nxc(2),       &
                                     Lx(3) / (nxc(3) - 1) /)
real (rp), parameter :: dxf(nd) = (/ Lx(1) / nxf(1),       &
                                     Lx(2) / nxf(2),       &
                                     Lx(3) / (nxf(3) - 1) /)

real (rp), parameter :: cutoff = dxc(1)  !--this is variable

real (rp) :: phic(nxc(1)+2, nxc(2), nxc(3))
real (rp) :: phif(nxf(1)+2, nxf(2), nxf(3))  !--not needed right now

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine read_phi (file, nx, phi)
implicit none

character (*), intent (in) :: file
integer, intent (in) :: nx(nd)

real (rp), intent (out) :: phi(nx(1)+2, nx(2), nx(3))

integer, parameter :: lun = 1

logical :: exst, opn

!---------------------------------------------------------------------

inquire (lun, exist=exst, opened=opn)

if ((.not. exst) .or. opn) then
  write (*, *) 'read_phi: problem with lun'
  stop
end if

inquire (file=file, exist=exst, opened=opn)

if ((.not. exst) .or. opn) then
  write (*, *) 'read_phi: problem with file ' // file
end if

open (lun, file=file, form='unformatted')

read (lun) phi

close (lun)

end subroutine read_phi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine read_uvw (file, nx, u, v, w)
implicit none

character (*), intent (in) :: file
integer, intent (in) :: nx(nd)

real (rp), intent (out) :: u(nx(1)+2, nx(2), nx(3)),  &
                           v(nx(1)+2, nx(2), nx(3)),  &
                           w(nx(1)+2, nx(2), nx(3))

integer, parameter :: lun = 1

logical :: exst, opn

!---------------------------------------------------------------------

inquire (lun, exist=exst, opened=opn)

if ((.not. exst) .or. opn) then
  write (*, *) 'read_uvw: problem with lun'
  stop
end if

inquire (file=file, exist=exst, opened=opn)

if ((.not. exst) .or. opn) then
  write (*, *) 'read_uvw: problem with file ' // file
end if

open (lun, file=file, form='unformatted')

read (lun) u, v, w  !--ignore rest of stuff in this file

close (lun)

end subroutine read_uvw

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine linf_error (uc, ufi, linf)
implicit none

real (rp), intent (in) :: uc(nxc(1)+2, nxc(2), nxc(3))
real (rp), intent (in) :: ufi(nxc(1)+2, nxc(2), nxc(3))
                          !--has beeninterpolated to coarse grid

real (rp), intent (out) :: linf

integer :: ic, jc, kc

!---------------------------------------------------------------------

linf = 0._rp

do kc = 1, nxc(3) - 1
  do jc = 1, nxc(2)
    do ic = 1, nxc(1)

      if (phic(ic, jc, kc) > cutoff) then
        linf = max ( linf, abs (uc(ic, jc, kc) - ufi(ic, jc, kc)) )
      end if

    end do
  end do
end do

end subroutine linf_error

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--a more sophisticated summing strategy may be required here 
!--may want to restrict this to regions with phi > 0
!--by default, the lp norm is unnormalized
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine lp_error (p, uc, ufi, lp, normalized)
implicit none

integer, intent (in) :: p
real (rp), intent (in) :: uc(nxc(1)+2, nxc(2), nxc(3))
real (rp), intent (in) :: ufi(nxc(1)+2, nxc(2), nxc(3))

real (rp), intent (out) :: lp

character (*), intent (in), optional :: normalized

integer :: counter
integer :: ic, jc, kc

!---------------------------------------------------------------------

if (p <= 0) then
  write (*, *) 'lp_error: p should be > 0'
  stop
end if

lp = 0._rp
counter = 0

do kc = 1, nxc(3) - 1
  do jc = 1, nxc(2)
    do ic = 1, nxc(1)

      if (phic(ic, jc, kc) > cutoff) then
        lp = lp + abs (uc(ic, jc, kc) - ufi(ic, jc, kc))**p
        counter = counter + 1
      end if

    end do
  end do
end do

if (present (normalized)) then

  select case (normalized)
    case ('yes')
      lp = lp / counter
    case ('no')
      !--default is unnormalized: do nothing
    case default
      write (*, *) 'lp_error: invalid normalized = ' // normalized
  end select

end if

lp = lp**(1._rp / p)

end subroutine lp_error

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--measures error in uc relative to uf (interpd to coarse grid)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine err_ptwise  (ufi, uc, ferr, node)
implicit none

real (rp), intent (in) :: ufi(nxc(1)+2, nxc(2), nxc(3))
real (rp), intent (in) :: uc(nxc(1)+2, nxc(2), nxc(3))

character (*), intent (in) :: ferr
character (*), intent (in), optional :: node

integer :: i, j, k

real (rp) :: err
real (rp) :: s
real (rp) :: x, y, z

!---------------------------------------------------------------------

if (present (node)) then

  select case (node)
    case ('u'); s = 0.5_rp
    case ('w'); s = 1._rp
    case default
      write (*, *) 'err_ptwise: invalid node=', node
      stop
  end select

else

  s = 0.5_rp  !--default is u-nodes    

end if

open (1, file=ferr)
write (1, '(a)') 'variables = "x" "y" "z" "err"'
write (1, '(3(a,i0))') 'zone, f=point, i=', nxc(1), ', j=', nxc(2),  &
                       ', k=', nxc(3) - 1

do k = 1, nxc(3) - 1

  z = (k - s) * dxc(3)

  do j = 1, nxc(2)

    y = (j - 1) * dxc(2)
    
    do i = 1, nxc(1)

      x = (i - 1) * dxc(1)
      
      err = uc(i, j, k) - ufi(i, j, k)

      write (1, '(4(es12.5,1x))') x, y, z, err
      
    end do
    
  end do

end do

close (1)

end subroutine err_ptwise

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine box_filter2 (u, u_filt, node)
use integrate, only : integrate_1d
implicit none

real (rp), intent (in) :: u(nxf(1)+2, nxf(2), nxf(3))
real (rp), intent (out) :: u_filt(nxc(1)+2, nxc(2), nxc(3))
character (*), intent (in), optional :: node  !--'u' or 'w'

integer :: ic, jc, kc
integer :: i, j, k
integer :: ii, jj
integer :: imin, imax, jmin, jmax, kmin, kmax
integer :: npts

real (rp) :: delta(nd)  !--filter size
real (rp) :: filtval
real (rp) :: s
real (rp) :: x, y, z
real (rp) :: xarr(nxf(1)), yarr(nxf(2)), zarr(nxf(3))
real (rp) :: xmin, xmax, ymin, ymax, zmin, zmax
real (rp) :: u_xfilt (nxc(1)+2, nxf(2), nxf(3))  !--yes, its mixed size
real (rp) :: u_yfilt (nxc(1)+2, nxc(2), nxf(3))  !--yes, its mixed size

!---------------------------------------------------------------------

if (present (node)) then

  select case (node)
    case ('u'); s = 0.5_rp
    case ('w'); s = 1._rp
    case default; write (*, *) 'interp: invalid node'; stop
  end select
  
else

  s = 0.5_rp  !--default is u-node

end if

delta = dxc

!--x-filtering: fine grid
write (*, *) 'box_filter2: starting x-filtering'

xarr = (/ ( (i - 1) * dxf(1), i=1, nxf(1) ) /)

do k = 1, nxf(3) - 1
  do j = 1, nxf(2)
    do i = 1, nxc(1)  !--yes, its coarse
    
      x = (i - 1) * dxc(1)

      xmin = x - 0.5_rp * delta(1)
      xmax = x + 0.5_rp * delta(1)

      u_xfilt(i, j, k) = ( integrate_1d (.true., xmin, xmax, dxf(1),  &
                                         xarr, u(1:nxf(1), j, k))   &
                         ) / delta(1)
      
    end do
  end do
end do

!--y-filtering: fine grid
write (*, *) 'box_filter2: starting y-filtering'

yarr = (/ ( (j - 1) * dxf(2), j=1, nxf(2) ) /)

do k = 1, nxf(3) - 1
  do j = 1, nxc(2)  !--yes, its coarse

    y = (j - 1) * dxc(2)
    ymin = y - 0.5_rp * delta(2)
    ymax = y + 0.5_rp * delta(2)

    do i = 1, nxc(1)  !--yes, its coarse

      u_yfilt(i, j, k) = ( integrate_1d (.true., ymin, ymax, dxf(2),  &
                                         yarr, u_xfilt(i, :, k))      &
                         ) / delta(2)

    end do

  end do
end do

!--z-filtering: also now convert to coarse grid
write (*, *) 'box_filter2: starting z-filtering'

!--may need to be more careful near boundary
zarr = (/ ( (k - s) * dxf(3), k=1, nxf(3) ) /)

do k = 2, nxc(3) - 1

  z = (k - s) * dxc(3)
  zmin = z - 0.5_rp * delta(3)
  zmax = z + 0.5_rp * delta(3)

  do j = 1, nxc(2)
    do i = 1, nxc(1)

      u_filt(i, j, k) = ( integrate_1d (.false., zmin, zmax, dxf(3),  &
                                        zarr, u_yfilt(i, j, :))       &
                        ) / delta(3)
                                       
    end do
  end do
  
end do

end subroutine box_filter2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine box_filter (u, u_filt, node)
implicit none

real (rp), intent (in) :: u(nxf(1)+2, nxf(2), nxf(3))
real (rp), intent (out) :: u_filt(nxc(1)+2, nxc(2), nxc(3))
character (*), intent (in), optional :: node  !--'u' or 'w'

integer :: ic, jc, kc
integer :: i, j, k
integer :: ii, jj
integer :: imin, imax, jmin, jmax, kmin, kmax
integer :: npts

real (rp) :: delta(nd)  !--filter size
real (rp) :: filtval
real (rp) :: s
real (rp) :: x, y, z
real (rp) :: xmin, xmax, ymin, ymax, zmin, zmax

!---------------------------------------------------------------------

if (present (node)) then

  select case (node)
    case ('u'); s = 0.5_rp
    case ('w'); s = 1._rp
    case default; write (*, *) 'interp: invalid node'; stop
  end select
  
else

  s = 0.5_rp  !--default is u-node

end if

delta = dxc

do kc = 1, nxc(3) - 1

  z = (kc - s) * dxc(3)

  zmin = z - 0.5_rp * delta(3)
  zmax = z + 0.5_rp * delta(3)

  kmin = ceiling (zmin / dxf(3) + s)
  kmax = floor (zmax / dxf(3) + s)

  do jc = 1, nxc(2)

    y = (jc - 1) * dxc(2)

    ymin = y - 0.5_rp * delta(2)
    ymax = y + 0.5_rp * delta(2)

    jmin = ceiling (ymin / dxf(2) + 1._rp)
    jmax = floor (ymax / dxf(2) + 1._rp)
    
    do ic = 1, nxc(1)

      x = (ic - 1) * dxc(1)

      xmin = x - 0.5_rp * delta(1)
      xmax = x + 0.5_rp * delta(1)

      !--this ceiling floor convention slightly decreases filter size
      !--this needs to be made more precise (2nd order)
      imin = ceiling (xmin / dxf(1) + 1._rp)
      imax = floor (xmax / dxf(1) + 1._rp)
      
      npts = 0
      filtval = 0._rp

      do k = max (1, kmin), min (nxf(3) - 1, kmax)

        do j = jmin, jmax

          jj = modulo (j - 1, nxf(2)) + 1
          
          do i = imin, imax
          
            ii = modulo (i - 1, nxf(1)) + 1 
            
            filtval = filtval + u(ii, jj, k)
            npts = npts + 1

          end do
        end do
      end do

      if (npts > 0) filtval = filtval / npts

      u_filt(ic, jc, kc) = filtval
      
    end do
  end do
end do

end subroutine box_filter

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine interp (u, u_interp, node)
implicit none

real (rp), intent (in) :: u(nxf(1)+2, nxf(2), nxf(3))
real (rp), intent (out) :: u_interp(nxc(1)+2, nxc(2), nxc(3))
character (*), intent (in), optional :: node  !--'u' or 'w'

integer :: i_c, j_c, k_c
integer :: i_f, j_f, k_f
integer :: i_f1, j_f1, k_f1

real (rp) :: x, y, z
real (rp) :: wx, wy, wz
real (rp) :: u_x00, u_x10, u_xy0, u_x01, u_x11, u_xy1
real (rp) :: s

!---------------------------------------------------------------------

if (present (node)) then

  select case (node)
    case ('u'); s = 0.5_rp
    case ('w'); s = 1._rp
    case default; write (*, *) 'interp: invalid node'; stop
  end select
  
else

  s = 0.5_rp  !--default is u-node

end if

do k_c = 1, nxc(3) - 1  !--skip last point

  z = (k_c - s) * dxc(3)

  k_f = floor (z / dxf(3) + s)
  if ( (k_f < 1) .or. (k_f > nxf(3)) ) then
    write (*, *) 'k_f out of bounds'
    stop
  end if
  
  wz = ((k_f - s) * dxf(3) - z) / dxf(3)

  do j_c = 1, nxc(2)

    y = (j_c - 1) * dxc(2)

    j_f = floor (y / dxf(2) + 1._rp)
    if ( (j_f < 1) .or. (j_f > nxf(2)) ) then
      write (*, *) 'j_f out of bounds'
      stop
    end if

    j_f1 = modulo (j_f, nxf(2)) + 1  !--wraps boundaries

    wy = ((j_f - 1) * dxf(2) - y) / dxf(2)

    do i_c = 1, nxc(1)

      x = (i_c - 1) * dxc(1)
      
      i_f = floor (x / dxf(1) + 1._rp)
      if ( (i_f < 0) .or. (i_f > nxf(1)) ) then
        write (*, *) 'i_f out of bounds'
        stop
      end if

      i_f1 = modulo (i_f, nxf(1)) + 1  !--wraps boundaries

      wx = ((i_f - 1) * dxf(1) - x) / dxf(1)
      
      u_x00 = u(i_f, j_f , k_f ) * (1._rp - wx) + u(i_f1, j_f , k_f ) * wx
      u_x10 = u(i_f, j_f1, k_f ) * (1._rp - wx) + u(i_f1, j_f1, k_f ) * wx
      u_x01 = u(i_f, j_f , k_f1) * (1._rp - wx) + u(i_f1, j_f , k_f1) * wx
      u_x11 = u(i_f, j_f1, k_f1) * (1._rp - wx) + u(i_f1, j_f1, k_f1) * wx

      u_xy0 = u_x00 * (1._rp - wy) + u_x10 * wy
      u_xy1 = u_x01 * (1._rp - wy) + u_xy1 * wy
     
      u_interp(i_c, j_c, k_c) = u_xy0 * (1._rp - wz) + u_xy1 * wz
      
    end do
  end do
end do

end subroutine interp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module verify_level_set

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program main
use verify_level_set
implicit none

character (*), parameter :: prec = '64x32x17'
character (*), parameter :: pref = '96x48x25'

character (*), parameter :: pathc = prec // '/'
character (*), parameter :: pathf = pref // '/'

character (*), parameter :: fphic = pathc // 'phi.out'
character (*), parameter :: fphif = pathf // 'phi.out'

character (*), parameter :: fuvwc = pathc // 'output/vel-avg.out'
character (*), parameter :: fuvwf = pathf // 'output/vel-avg.out'

character (*), parameter :: fptwise_int = 'err-ptwise-int--' //       &
                                          prec // '-' // pref // '.dat'
character (*), parameter :: fptwise_filt = 'err-ptwise-filt--' //      &
                                           prec // '-' // pref // '.dat'
character (*), parameter :: fptwise_filt2 = 'err-ptwise-filt2--' //      &
                                            prec // '-' // pref // '.dat'

real (rp) :: l1_err, l2_err, linf_err
real (rp) :: uc(nxc(1)+2, nxc(2), nxc(3)),  &
             vc(nxc(1)+2, nxc(2), nxc(3)),  &
             wc(nxc(1)+2, nxc(2), nxc(3))
real (rp) :: uf(nxf(1)+2, nxf(2), nxf(3)),  &
             vf(nxf(1)+2, nxf(2), nxf(3)),  &
             wf(nxf(1)+2, nxf(2), nxf(3))

real (rp) :: uf_int(nxc(1)+2, nxc(2), nxc(3))
real (rp) :: uf_filt(nxc(1)+2, nxc(2), nxc(3))

!---------------------------------------------------------------------

!--read coarse stuff
call read_phi (fphic, nxc, phic)
call read_uvw (fuvwc, nxc, uc, vc, wc)

!--read fine stuff
call read_phi (fphif, nxf, phif)
call read_uvw (fuvwf, nxf, uf, vf, wf)

!--interpolate fine stuff onto coarse grid
call interp (uf, uf_int)

call lp_error (1, uc, uf_int, l1_err, normalized='yes')
call lp_error (2, uc, uf_int, l2_err, normalized='yes')
call linf_error (uc, uf_int, linf_err)

write (*, *) 'interpolation:'
write (*, *) 'l_1 error = ', l1_err
write (*, *) 'l_2 error = ', l2_err
write (*, *) 'l_inf  = ', linf_err

!--point-wise error
call err_ptwise (uf_int, uc, fptwise_int)

!--filter fine stuff and put on coarse grid
call box_filter (uf, uf_filt)

call lp_error (1, uc, uf_filt, l1_err, normalized='yes')
call lp_error (2, uc, uf_filt, l2_err, normalized='yes')
call linf_error (uc, uf_filt, linf_err)

write (*, *) 'filtering:'
write (*, *) 'l_1 error = ', l1_err
write (*, *) 'l_2 error = ', l2_err
write (*, *) 'l_inf  = ', linf_err

!--point-wise error
call err_ptwise (uf_filt, uc, fptwise_filt)

!--experimental box filter
call box_filter2 (uf, uf_filt)

call lp_error (1, uc, uf_filt, l1_err, normalized='yes')
call lp_error (2, uc, uf_filt, l2_err, normalized='yes')
call linf_error (uc, uf_filt, linf_err)

write (*, *) 'filtering:'
write (*, *) 'l_1 error = ', l1_err
write (*, *) 'l_2 error = ', l2_err
write (*, *) 'l_inf  = ', linf_err

!--point-wise error
call err_ptwise (uf_filt, uc, fptwise_filt2)

end program main
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
