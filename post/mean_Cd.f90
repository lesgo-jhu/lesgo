program mean_Cd
implicit none

integer, parameter :: rp = kind (1.d0)

character (*), parameter :: file_in = 'Cd_dynamic.dat'

integer :: counter
integer :: junk
integer :: iost

real (rp) :: Cd, Cd_avg

!---------------------------------------------------------------------

open (1, file=file_in)

Cd_avg = 0._rp
counter = 0

read_file_loop: do

  read (1, *, iostat=iost) junk, Cd
  if (iost /= 0) exit read_file_loop

  Cd_avg = Cd_avg + Cd
  counter = counter + 1

end do read_file_loop

close (1)

if (counter > 0) then
  Cd_avg = Cd_avg / counter
  write (*, *) 'Cd_avg = ', Cd_avg
  write (*, *) 'counter = ', counter
else
  write (*, *) 'error, counter = 0'
  stop
end if

end program mean_Cd
