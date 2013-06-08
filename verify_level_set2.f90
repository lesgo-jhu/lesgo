!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--driver program in this file, after the module
!--the "2" version only work for r = 2 (i.e. grid doubling)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module verify_level_set2
implicit none

save
public

character (*), parameter :: mod_name = 'verify_level_set2'

integer, parameter :: rp = kind (1.d0)
integer, parameter :: nd = 3

integer, parameter :: nxc(nd) = (/ 32, 32, 9 /)
integer, parameter :: nxf(nd) = (/ 64, 64, 17/)

real (rp), parameter :: Lx(nd) = (/ 1._rp, 1._rp, 0.25_rp /)
                        !--Lx should be same for coarse, fine otherwise
                        !  why are we trying to compare them?

real (rp), parameter :: dxc(nd) = (/ Lx(1) / nxc(1),       &
                                     Lx(2) / nxc(2),       &
                                     Lx(3) / (nxc(3) - 1) /)
real (rp), parameter :: dxf(nd) = (/ Lx(1) / nxf(1),       &
                                     Lx(2) / nxf(2),       &
                                     Lx(3) / (nxf(3) - 1) /)

real (rp), parameter :: cutoff = 2 * dxc(1)  !--this is variable

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
!--this assumes r = 2
!--direct sampling, no filtering (except necessary interpolation)
!--interpolation for u-node quantities
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine sample (u, u_samp, node)
implicit none

real (rp), intent (in) :: u(nxf(1)+2, nxf(2), nxf(3))
real (rp), intent (out) :: u_samp(nxc(1)+2, nxc(2), nxc(3))
character (*), intent (in), optional :: node  !--'u' or 'w'

integer, parameter  :: r = 2

integer :: i, j, k
integer :: ic, jc, kc
integer :: s

!---------------------------------------------------------------------

!--check node is valid, is it present
if (present (node)) then

  select case (node)
    case ('u'); s =  1
    case ('w'); s = 0
    case default; write (*, *) 'sample: invalid node'; stop
  end select
  
else  !--assume u-node

  s = 1

end if

!--check r = 2
if ((nxf(1) / r) /= nxc(1)) then
  write (*, *) 'sample: problem with nxf(1), nxc(1)'
  stop
else if ((nxf(2) / r) /= nxc(2)) then
  write (*, *) 'sample: problem with nxf(2), nxc(2)'
  stop
else if ( ((nxf(3) - 1) / r) /= (nxc(3) - 1) ) then
  write (*, *) 'sample: problem with nxf(3), nxc(3)'
  stop
end if

!--interpolation required in z-direction for u-nodes
!--for w-nodes, no interp required, so k + s = k
do kc = 1, nxc(3) - 1
  do jc = 1, nxc(2)
    do ic = 1, nxc(1)

      i = r * ic - 1
      j = r * jc - 1
      k = r * kc - 1  !--careful

      u_samp(ic, jc, kc) = 0.5_rp * (u(i, j, k) + u(i, j, k + s))

    end do
  end do
end do

end subroutine sample
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module verify_level_set2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program main
use verify_level_set2
implicit none

character (*), parameter :: prec = '32x32x9-del2'
character (*), parameter :: pref = '64x64x17-del4'

character (*), parameter :: pathc = prec // '/run/'
character (*), parameter :: pathf = pref // '/run/'

character (*), parameter :: fphic = pathc // 'phi.out'
character (*), parameter :: fphif = pathf // 'phi.out'

character (*), parameter :: fuvwc = pathc // 'output/vel-avg.out'
character (*), parameter :: fuvwf = pathf // 'output/vel-avg.out'

character (*), parameter :: fptwise_samp = 'err-ptwise-samp--' //      &
                                           prec // '-' // pref // '.dat'
character (*), parameter :: fptwise_nondim_samp =                         &
                'err-ptwise-nondim-samp--' // prec // '-' // pref // '.dat'

real (rp) :: l1_err, l2_err, linf_err
real (rp) :: Uinf_c, Uinf_f
real (rp) :: uc(nxc(1)+2, nxc(2), nxc(3)),  &
             vc(nxc(1)+2, nxc(2), nxc(3)),  &
             wc(nxc(1)+2, nxc(2), nxc(3))
real (rp) :: uf(nxf(1)+2, nxf(2), nxf(3)),  &
             vf(nxf(1)+2, nxf(2), nxf(3)),  &
             wf(nxf(1)+2, nxf(2), nxf(3))
real (rp) :: uf_samp(nxc(1)+2, nxc(2), nxc(3))

!---------------------------------------------------------------------

!--read coarse stuff
call read_phi (fphic, nxc, phic)
call read_uvw (fuvwc, nxc, uc, vc, wc)

!--read fine stuff
call read_phi (fphif, nxf, phif)
call read_uvw (fuvwf, nxf, uf, vf, wf)

!--calculate Uinf, for normalization(s)
!Uinf_c = sum (uc(1, :, :)) / (nxc(2) * nxc(3))
!Uinf_f = sum (uf(1, :, :)) / (nxf(2) * nxf(3))
Uinf_c = sum (uc(1, 1, :)) / nxc(3)
Uinf_f = sum (uf(1, 1, :)) / nxf(3)

write (*, *) 'Uinf_c = ', Uinf_c
write (*, *) 'Uinf_f = ', Uinf_f

!--sample fine stuff directly: r = 2 here
call sample (uf, uf_samp)

call lp_error (1, uc, uf_samp, l1_err, normalized='yes')
call lp_error (2, uc, uf_samp, l2_err, normalized='yes')
call linf_error (uc, uf_samp, linf_err)

!--normalize errors by Uinf_f
write (*, *) 'sampling:'
write (*, *) 'e_1/Uinf_f = ', l1_err / Uinf_f
write (*, *) 'e_2/Uinf_f = ', l2_err / Uinf_f
write (*, *) 'e_inf/Uinf_f  = ', linf_err / Uinf_f

!--this has not been normalized by Uinf_f
call err_ptwise (uf_samp, uc, fptwise_samp)

!--repeat the analysis, but with u/Uinf_c and u/Uinf_f
!--we are using the exact Uinf_f, not the sampled one...does it matter?
uc = uc / Uinf_c
uf_samp = uf_samp / Uinf_f

call lp_error (1, uc, uf_samp, l1_err, normalized='yes')
call lp_error (2, uc, uf_samp, l2_err, normalized='yes')
call linf_error (uc, uf_samp, linf_err)

!--normalize errors by Uinf_f
write (*, *) 'sampling:'
write (*, *) 'e*_1 = ', l1_err
write (*, *) 'e*_2 = ', l2_err
write (*, *) 'e*_inf = ', linf_err

!--this has not been normalized by Uinf_f
call err_ptwise (uf_samp, uc, fptwise_nondim_samp)

end program main
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
