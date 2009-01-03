program compare

integer, parameter :: rp = kind (1.d0)

character (*), parameter :: file1 = 'vel000100.dat'
character (*), parameter :: file2 = '../../nompi-run/output/vel000100.dat'
character (*), parameter :: outfile = 'compare.fortran.out'

character (128) :: hdr2

real (rp) :: x1, y1, z1, u1, v1, w1
real (rp) :: x2, y2, z2, u2, v2, w2

open (21, file=file1)
!--advance past 2 line header
read (21, *)
read (21, '(a)') hdr2

open (22, file=file2)
read (22, *)
read (22, *)

open (1, file=outfile)
write (1, '(a)') 'variables = "x" "y" "z" "du" "dv" "dw"'
write (1, '(a)') hdr2

do

  read (21, *, iostat=ios1) x1, y1, z1, u1, v1, w1
  read (22, *, iostat=ios2) x2, y2, z2, u2, v2, w2
  if ((ios1 /= 0) .or. (ios2 /= 0)) then
    exit
  end if
  write (1, '(6(es12.5,1x))') x1, y1, z1, u1-u2, v1-v2, w1-w2
  
end do

end program compare
