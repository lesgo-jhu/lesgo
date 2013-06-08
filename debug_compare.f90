program debug_compare
implicit none

integer, parameter :: rp = kind (1.d0)

character (*), parameter :: outfile = 'debug_compare.dat'

character (128) :: file1, file2
character (128) :: line1, line2
character (128) :: path
character (128) :: hdr

integer :: ios1, ios2
integer :: nlines
integer :: ni, nj, nk
integer :: i1, j1, k1
integer :: i2, j2, k2

real (rp) :: max_abs_du
real (rp) :: eps
real (rp) :: du
real (rp) :: x1, y1, z1, u1
real (rp) :: x2, y2, z2, u2

write (*, *) 'Enter non-mpi file name:'
read (*, '(a)') file1
write (*, *) 'Enter relative path to MPI file:'
read (*, '(a)') path
write (*, *) 'Enter absolute tolerance:'
read (*, *) eps

file2 = trim (path) // '/' // trim (file1) // '.MPI'

write (*, *) 'using MPI file ' // file2
write (*, *) 'using tolerance eps =', eps

open (21, file=file1)
read (21, '(a)') hdr  !--of the form #debug_mod nx ny nz

!--extract ni, nj, nk
!--11 since 10 chars in '#debug_mod'
read (hdr(11:), *) ni, nj, nk
write (*, '(a,3(i0,1x))') 'ni, nj, nk = ', ni, nj, nk

open (22, file=file2)
read (22, *)  !--skip header

open (1, file=outfile)
write (1, '(a)') '#produced by debug_compare'
write (1, '(a)') '#non-mpi file:' // trim (file1)
write (1, '(a)') '#mpi-file:' // trim (file2)
write (1, '(a)') 'variables = "i" "j" "k" "diff"'
write (1, '(3(a,i0))') 'zone, f=point, i=', ni, ',j=', nj, ',k=', nk 

nlines = 0
max_abs_du = 0._rp

do

  read (21, '(a)', iostat=ios1) line1
  read (22, '(a)', iostat=ios2) line2
  if ((ios1 /= 0) .or. (ios2 /= 0)) then
    exit
  end if

  !--remove colon from line1, line2
  line1 = line1(1:index(line1, ':')-1) // line1(index(line1, ':')+1:)
  line2 = line2(1:index(line2, ':')-1) // line2(index(line2, ':')+1:)

  read (line1, *) i1, j1, k1, u1
  read (line2, *) i2, j2, k2, u2

  du = u1 - u2
  max_abs_du = max (du, max_abs_du)

  if (abs (du) > eps) then
    write (*, '(a,3(i0,1x))') 'tolerance exceeded at i, j, k =', i1, j1, k1
  end if
  write (1, '(3(i0,1x),es17.10)') i1, j1, k1, du
  nlines = nlines + 1
  
end do

write (*, '(a,i0,a)') 'wrote ', nlines, ' lines'
write (*, *) 'largest difference = ', max_abs_du

end program debug_compare
