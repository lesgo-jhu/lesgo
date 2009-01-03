program cylinder_cd
implicit none

integer, parameter :: rp = kind (1.d0)

character (len=*), parameter :: prog_name = 'cylinder_cd'

character (len=*), parameter :: in_avg = 'Cd_cyl_avg.dat'
character (len=*), parameter :: in_local = 'Cd_cyl_local.dat'

character (len=*), parameter :: out_avg = 'post_Cd_cyl_avg.dat'
character (len=*), parameter :: out_local = 'post_Cd_cyl_local.dat'

integer, parameter :: nd = 3
integer, parameter :: n_method = 3

integer :: iost
integer :: ncalls
integer :: counter

real (rp) :: f(nd)

!---------------------------------------------------------------------

call process_avg ()

call process_local ()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine process_local ()
implicit none

character (len=*), parameter :: sub_name = prog_name // '.process_local'

character (len=20) :: fmt_str

integer :: i_pt, n_pt

real (rp), allocatable :: Cd(:, :)
real (rp), allocatable :: Cd_avg(:, :)

!---------------------------------------------------------------------

open (1, file = in_local)

! figure out number of pts in the cylinder by counting number of lines
! until first blank line or number of lines until ncalls changes
! or until first index (i_pt) goes back to 1
! third approach

n_pt = 0

count_npt_loop: do

  read (1, *, iostat=iost) i_pt ! forget the other stuff

  if (iost /= 0) then
    write (*, *) sub_name // ': file ended before n_pt determined'
    stop
  end if

  if (i_pt < n_pt) exit count_npt_loop

  n_pt = i_pt
  
end do count_npt_loop

if (n_pt < 1) then
  write (*, *) sub_name // ': n_pt < 1'
  stop
end if

write (*, *) sub_name // ': determined n_pt = ', n_pt

allocate (Cd(n_method, n_pt))
allocate (Cd_avg(n_method, n_pt))

Cd_avg = 0._rp
counter = 0

rewind (1)

read_local_loop: do

  do i_pt = 1, n_pt

    read (1, *, iostat=iost) i_pt, ncalls, f, Cd(:, i_pt)

    if (iost /= 0) then

      if (i_pt == 1) then
        exit read_local_loop
      else
        write (*, *) sub_name // ': file ended in mid vertical sweep'
        write (*, *) sub_name // ': i_pt = ', i_pt
        stop
      end if

    end if

  end do

  Cd_avg = Cd_avg + Cd
  counter = counter + 1
  
end do read_local_loop

close (1)

Cd_avg = Cd_avg / counter

! correction for data error: A_p was n_pt times too big
! really A_p <- A_p / n_pt, so Cd <- Cd * n_pt
Cd_avg = Cd_avg * n_pt

open (1, file = out_local)

write (1, '(a)') '# Time avg Cd_cyl_local:'
write (1, '(a)') '# k, Cd(1), Cd(2), ...'

write (fmt_str, '(a,i1,a)') '(i2,', n_method, '(x, f12.5))'

do i_pt = 1, n_pt

  write (1, fmt_str) i_pt, Cd_avg(:, i_pt)

end do

close (1)

end subroutine process_local

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine process_avg ()
implicit none

integer :: i

real (rp) :: Cd(n_method)
real (rp) :: Cd_avg(n_method)

!---------------------------------------------------------------------

open (1, file = in_avg)

counter = 0

Cd_avg = 0._rp

read_avg_loop: do

  read (1, *, iostat=iost) ncalls, f, Cd
  if (iost /= 0) exit read_avg_loop

  counter = counter + 1
  Cd_avg = Cd_avg + Cd

end do read_avg_loop

close (1)

if (counter > 0) then
  Cd_avg = Cd_avg / counter
else
  write (*, *) prog_name // ': error, counter = ', counter
  stop
end if

open (1, file = out_avg)
write (1, '(a)') 'Time avgs of Cd_cyl_avg:'

do i = 1, n_method
  write (1, '(a,i1,a,f12.5)') 'Method ', i, ': ', Cd_avg(i)
end do

close (1)
end subroutine process_avg
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end program cylinder_cd
