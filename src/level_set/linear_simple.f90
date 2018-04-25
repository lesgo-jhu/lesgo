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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! to solve small, nonsingular  linear systems
! reference numerical recipes (fortran 90/95)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module linear_simple
!use precision
use param, only : rp => rprec
use messages
implicit none

private

public :: solve_linear

character (len=*), parameter :: mod_name = 'linear_simple'

interface assert_eq
  module procedure assert_eq2, assert_eq3, assert_eq4
end interface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! solves Ax = b for x
! does NOT destroy A or b
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine solve_linear (A, b, x)
implicit none

real (rp), intent (in) :: A(:, :), b(:)
real (rp), intent (out) :: x(size (b))

character (len=*), parameter :: sub_name = mod_name // '.solve_linear'

integer :: indx(size (A, 1))

real (rp) :: A_tmp(size (A, 1), size (A, 2))
real (rp) :: d

!----------------------------------------------------------------------

A_tmp = A
call ludcmp (A_tmp, indx, d)

x = b
call lubksb (A_tmp, indx, x)

end subroutine solve_linear

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! index of maxloc on an array
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!FUNCTION imaxloc (arr)
!INTEGER :: imaxloc
!
!REAL (RP), DIMENSION(:), INTENT(IN) :: arr
!
!!----------------------------------------------------------------------
!
!imaxloc = maxloc (arr(:), dim = 1)
!
!END FUNCTION imaxloc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Report and die if integers not all equal (used for size checking)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION assert_eq2(n1,n2,string)
CHARACTER(LEN=*), INTENT(IN) :: string
INTEGER, INTENT(IN) :: n1,n2
INTEGER :: assert_eq2
CHARACTER(LEN=1024) :: msg

character (len=*), parameter :: sub_name = mod_name // '.assert_eq2'

!----------------------------------------------------------------------

if (n1 == n2) then
  assert_eq2=n1
else
  write (msg, *) 'failed assert with tag:', n_l, string
  call error (sub_name, msg)
end if

END FUNCTION assert_eq2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION assert_eq3(n1,n2,n3,string)
CHARACTER(LEN=*), INTENT(IN) :: string
INTEGER, INTENT(IN) :: n1,n2,n3
INTEGER :: assert_eq3
CHARACTER(LEN=1024) :: msg

character (len=*), parameter :: sub_name = mod_name // '.assert_eq3'

!----------------------------------------------------------------------

if (n1 == n2 .and. n2 == n3) then
  assert_eq3=n1
else
  write (msg, *) 'failed assert with tag:', n_l, string
  call error (sub_name, msg)
end if

END FUNCTION assert_eq3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION assert_eq4(n1,n2,n3,n4,string)
CHARACTER(LEN=*), INTENT(IN) :: string
INTEGER, INTENT(IN) :: n1,n2,n3,n4
INTEGER :: assert_eq4
CHARACTER(LEN=1024) :: msg

character (len=*), parameter :: sub_name = mod_name // '.assert_eq4'

!----------------------------------------------------------------------

if (n1 == n2 .and. n2 == n3 .and. n3 == n4) then
  assert_eq4=n1
else
  write (msg, *) 'failed assert with tag:', n_l, string
  call error (sub_name, msg)
end if

END FUNCTION assert_eq4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Given an N × N input matrix a, this routine replaces it by the LU
! decomposition of a rowwise permutation of itself. On output, a is
! arranged as in equation (2.3.14); indx is an output vector of
! length N that records the row permutation effected by the partial
! pivoting: d is output as ±1 depending on whether the number of row
! interchanges was even or odd, respectively. This routine is used in
! combination with lubksb to solve linear equations or invert a matrix.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ludcmp (a, indx, d)
implicit none

real (rp), intent (inout) :: a(:, :)
integer, intent (out) :: indx(:)
real (rp), intent (out) :: d

character (len=*), parameter :: sub_name = mod_name // '.ludcmp'

real (rp), parameter :: TINY = 1.0e-20_rp  ! a small number

integer :: j,n,imax

real (rp) :: vv(size (a, 1))
             !--stores the implicit scaling of each row

!----------------------------------------------------------------------

n = assert_eq (size (a, 1), size(a, 2), size(indx), 'ludcmp')

d = 1.0_rp  ! no row interchanges yet

vv = maxval (abs (a), dim=2)  ! loop over rows to get implicit scaling info

if (any (vv == 0.0_rp)) call error (sub_name, 'singular matrix')
    !--There is a row of zeros

vv = 1.0_rp / vv  ! save the scaling

do j = 1, n

  ! find pivot row
  imax = (j - 1) + maxloc ( vv(j:n) * abs (a(j:n, j)), dim = 1 )

  ! do we need to interchange rows?
  if (j /= imax) then  ! yes

    ! interchange row, change d, and interchange scale factor
    call swap (a(imax, :), a(j, :))
    d = -d
    vv(imax) = vv(j)

  end if

  indx(j)=imax

  ! If the pivot element is zero the matrix is singular (at least
  ! to the precision of the algorithm). For some applications on
  ! singular matrices, it is desirable to substitute TINY for zero.
  if (a(j, j) == 0.0_rp) a(j, j) = TINY

  ! divide by the pivot element
  a(j+1:n, j) = a(j+1:n, j) / a(j, j)

  ! reduce remaining submatrix
  a(j+1:n, j+1:n) = a(j+1:n, j+1:n) -                  &
                    outerprod (a(j+1:n, j), a(j, j+1:n))

end do

end subroutine ludcmp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Solves the set of N linear equations AX = B. Here the N × N matrix a
! is input, not as the original matrix A, but rather as its LU
! decomposition, determined by the routine ludcmp.
! indx is input as the permutation vector of length N returned by ludcmp
! b is input as the right-hand-side vector B, also of length N, and
! returns with the solution vector X.
! a and indx are not modified by this routine and can be left in place
! for successive calls with di
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE lubksb(a,indx,b)
IMPLICIT NONE

REAL(RP), DIMENSION(:,:), INTENT(IN) :: a
INTEGER, DIMENSION(:), INTENT(IN) :: indx
REAL(RP), DIMENSION(:), INTENT(INOUT) :: b

INTEGER :: i,n,ii,ll
REAL(RP) :: summ

!----------------------------------------------------------------------

n=assert_eq(size(a,1),size(a,2),size(indx),'lubksb')

ii=0
! When ii is set to a positive value, it will become the index of the
! first nonvanishing element of b. We now do the forward substitution,
! equation (2.3.6). The only new wrinkle is to unscramble the
! permutation as we go.

do i=1,n

  ll=indx(i)
  summ=b(ll)
  b(ll)=b(i)

  if (ii /= 0) then

    summ = summ - dot_product(a(i,ii:i-1),b(ii:i-1))

  else if (summ /= 0.0) then

    ii=i
    ! a nonzero element was encountered, so from now on we will have
    ! to do the dot product above

  end if

  b(i)=summ

end do

do i=n,1,-1

  ! now we do the backsubstitution, equation (2.3.7)
  b(i) = (b(i) - dot_product(a(i,i+1:n),b(i+1:n))) / a(i,i)

end do

END SUBROUTINE lubksb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION outerprod(a,b)
REAL(RP), DIMENSION(:), INTENT(IN) :: a,b
REAL(RP), DIMENSION(size(a),size(b)) :: outerprod

outerprod = spread(a,dim=2,ncopies=size(b)) *  &
            spread(b,dim=1,ncopies=size(a))

END FUNCTION outerprod

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE swap(a,b)
REAL(RP), DIMENSION(:), INTENT(INOUT) :: a,b
REAL(RP), DIMENSION(SIZE(a)) :: dum

dum=a
a=b
b=dum

END SUBROUTINE swap
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module linear_simple
