!!
!!  Copyright (C) 2016  Johns Hopkins University
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

module util
    
use types, only : rprec
use param, only : pi, CHAR_BUFF_LENGTH
use messages
use string_util
implicit none

public

interface softplus
    module procedure :: softplus_scalar
    module procedure :: softplus_array
end interface softplus

interface logistic
    module procedure :: logistic_scalar
    module procedure :: logistic_array
end interface logistic

interface gaussian
    module procedure :: gaussian_scalar
    module procedure :: gaussian_array
end interface gaussian

interface interpolate
    module procedure :: interpolate_scalar
    module procedure :: interpolate_array
end interface interpolate

contains

! softplus function of the form sp(x) = ln(1 + exp(x-s)) with scalar input
function softplus_scalar(s, x) result(sp)
    implicit none
    real(rprec), intent(in)     :: x, s
    real(rprec)                 :: sp
    real(rprec), parameter      :: threshold = 100

    if (x - s > threshold) then
        sp = x - s
    else
        sp = log(1 + exp(x - s))
    endif
end function softplus_scalar

! softplus function of the form ln(1 + exp(x-s)) with array input
function softplus_array(s, x) result(sp)
    implicit none
    real(rprec), dimension(:), intent(in)   :: x
    real(rprec), intent(in)                 :: s
    real(rprec), dimension(:), allocatable  :: sp
    integer                                 :: i
    
    allocate( sp(size(x)) )
    do i = 1, size(x)
        sp(i) = softplus(s, x(i))
    end do
end function softplus_array

! logistic function of the form 1/(1 + exp(-(x-s)) with scalar input
function logistic_scalar (s, x) result(l)
    implicit none
    real(rprec), intent(in)     :: x, s
    real(rprec)                 :: l
    real(rprec), parameter      :: threshold = 100

    if (x - s > threshold) then
        l = 1.0
    else
        l = 1.0 / ( 1.0 + exp( -(x - s) ) )
    endif
end function logistic_scalar

! logistic function of the form 1/(1 + exp(-(x-s)) with array input
function logistic_array(s, x) result(l)
    implicit none
    real(rprec), dimension(:), intent(in)   :: x
    real(rprec), intent(in)                 :: s
    real(rprec), dimension(:), allocatable  :: l
    integer                                 :: i
    
    allocate( l(size(x)) )
    do i = 1, size(x)
        l(i) = logistic(s, x(i))
    end do
end function logistic_array

! normalized Gaussian with scalar input
function gaussian_scalar(x, x0, Delta) result(g)
    implicit none
    real(rprec), intent(in) :: x, x0, Delta
    real(rprec)             :: g
    
    g = 1.0 / (Delta * sqrt(2.0 *pi)) * exp(-0.5 * (x - x0) * (x - x0) / Delta / Delta);
    ! If near precision limit, set to zero 
    ! NOTE: remove this after check
    if (abs(g) < sqrt(epsilon(g))) g = 0.0
    
end function gaussian_scalar

! normalized Gaussian with array input
function gaussian_array(x, x0, Delta) result(g)
    implicit none
    real(rprec), dimension(:), intent(in)   :: x
    real(rprec), intent(in)                 :: x0, Delta
    real(rprec), dimension(:), allocatable  :: g
    integer                                 :: i
    
    allocate( g(size(x)) )
    do i = 1, size(x)
        g(i) = gaussian(x(i), x0, Delta)
    end do
end function gaussian_array

! Third-order upwind biased differencing for d/dx operator.
! Uses second order upwind at third point and last point and first order upwind at second point
! Assumes boundary condition is already applied at first point of u
function ddx_upwind3biased(u, dx) result(dudx)
    implicit none
    real(rprec), dimension(:), intent(in)   :: u
    real(rprec), intent(in)                 :: dx
    real(rprec), dimension(:), allocatable  :: dudx
    integer                                 :: i, N
    
    N  = size(u)
    allocate(dudx(N))
    dudx(1) = 0.d0
    dudx(2) = (u(2) - u(1))/dx
    ! NOTE: Change remove this index and make do go from 3 to N-1
    dudx(3) = (1.5 * u(3) - 2.0 * u(2) - 0.5 * u(1))/dx
    do i = 4, N-1
        dudx(i) = ((1.0 / 3.0) * u(i+1) + 0.5 * u(i) - u(i-1) + (1.0 / 6.0) * u(i-2))/dx
    end do
    dudx(N) = (1.5 * u(N) - 2 * u(N-1) + 0.5 * u(N-2))/dx
end function ddx_upwind3biased

! Third-order downwind biased differencing for d/dx operator.
! Uses second order upwind at third point and last point and first order upwind at second point
! Assumes boundary condition is already applied at last point of u
function ddx_downwind3biased(u, dx) result(dudx)
    implicit none
    real(rprec), dimension(:), intent(in)   :: u
    real(rprec), intent(in)                 :: dx
    real(rprec), dimension(:), allocatable  :: dudx
    integer                                 :: i, N
    
    N  = size(u)
    allocate(dudx(N))
    dudx(N)   = 0.d0
    dudx(N-1) = (-u(N-1) + u(N))/dx;
    dudx(N-2) = (-1.5*u(N-2) + 2.0*u(N-1) - 0.5*u(N))/dx;
    ! NOTE: Change remove this index and make do go from N-2 to 2
    do i = N-3, 2, -1
        dudx(i) = ((-1.0 / 3.0) * u(i-1) - 0.5 * u(i) + u(i+1) - (1.0 / 6.0) * u(i+2))/dx;
    end do
    dudx(1) = (-1.5*u(1) + 2.0*u(2) - 0.5*u(3))/dx;
end function ddx_downwind3biased

! First-order upwind differencing for d/dx operator.
! Assumes boundary condition is already applied at first point of v
function ddx_upwind1(u, dx) result(dudx)
    implicit none
    real(rprec), dimension(:), intent(in)   :: u
    real(rprec), intent(in)                 :: dx
    real(rprec), dimension(:), allocatable  :: dudx
    integer                                 :: i, N
    
    N  = size(u)
    allocate(dudx(N))
    dudx(1) = 0.d0
    do i = 2, N
        dudx(i) =  ( u(i) - u(i-1) ) / dx
    end do
end function ddx_upwind1

! First-order downwind differencing for d/dx operator.
! Assumes boundary condition is already applied at last point of v
function ddx_downwind1(u, dx) result(dudx)
    implicit none
    real(rprec), dimension(:), intent(in)   :: u
    real(rprec), intent(in)                 :: dx
    real(rprec), dimension(:), allocatable  :: dudx
    integer                                 :: i, N
    
    N = size(u)
    allocate(dudx(N))
    dudx(N) = 0.d0
    do i = N-1, 1, -1
        dudx(i) =  ( -u(i) + u(i+1) ) / dx
    end do
end function ddx_downwind1

! Test function for minimization routines
subroutine rosenbrock(x, f, g)
    implicit none
    real(rprec), dimension(:), intent(in)                 :: x
    real(rprec), intent(inout)                            :: f
    real(rprec), dimension(:), allocatable, intent(inout) :: g
    
    if (size(x) /= 2) then
        call error('util.rosenbrock','Rosenbrock requires an array of size 2.')
    end if
    if (.not. allocated(g)) then
        allocate(g(2))
    end if

    ! Function value
    f = (1.0 - x(1))**2 + 100.0 * (x(2) - x(1)**2)**2;
    
    ! Gradient
    g(1) = -2 * (1.0 - x(1)) - 4 * x(1) * 100.0 * (x(2) - x(1)**2);
    g(2) = 2 * 100.0 * (x(2) - x(1)**2);
    
end subroutine rosenbrock

! linear interpolation between 1D-arrays
subroutine interpolate_scalar(x, y, xi, yi)
    implicit none
    real(rprec), dimension(:), intent(in)  :: x, y
    real(rprec), intent(in)                :: xi
    real(rprec), intent(out)               :: yi
    real(rprec), dimension(:), allocatable :: xi_array, yi_array
    
    allocate( xi_array(1) )
    allocate( yi_array(1) )
    xi_array(1) = xi
    call interpolate(x, y, xi_array, yi_array)
    yi = yi_array(1)
    
end subroutine interpolate_scalar

! linear interpolation between 1D-arrays
subroutine interpolate_array(x, y, xi, yi)
    implicit none
    real(rprec), dimension(:), intent(in)  :: x, y, xi
    real(rprec), dimension(:), intent(out) :: yi
    integer     :: i, j   
    real(rprec) :: dx, t 
    
    if ( size(x) /= size(y) .or. size(xi) /= size(yi)) then
        call error('util.inteprolate','Interpolation pairs must be of equal size.')
    end if
    if ( size(x) /= size(y) ) then
        call error('util.inteprolate','Interpolation pairs must be of equal size.')
    end if
    do i = 2, size(x)
        if ( x(i-1) >= x(i) ) then
            call error('util.interpolate','array x must be monotonically increasing')
        end if
    end do
    do i = 2, size(xi)
        if ( xi(i-1) >= xi(i) ) then
            call error('util.interpolate','array xi must be monotonically increasing')
        end if
    end do

    j = 1
    do i = 1, size(xi)
        if ( xi(i) <= x(1) ) then
            yi(i) = y(1)
        else if ( xi(i) >= x( size(x) ) ) then
            yi(i) = y( size(x) )
        else
            do while ( xi(i) > x(j+1) )
                j = j + 1
            end do
            dx = x(j+1) - x(j)
            t = ( xi(i) - x(j) ) / dx
            yi(i) = (1-t) * y(j) + t * y(j+1)
        end if
    end do
    
end subroutine interpolate_array

! Matrix inverse. This is a wrapper for Lapack. A must be a square matrix.
function inverse(A) result(Ainv)
    implicit none
    real(rprec), dimension(:,:), intent(in)    :: A
    real(rprec), dimension(:,:), allocatable   :: Ainv
    integer                                    :: N, INFO,  LWORK
    integer, dimension(:), allocatable         :: IPIV
    real(rprec), dimension(:), allocatable     :: WORK
    character(CHAR_BUFF_LENGTH)                :: err_str
    
    ! Require a square matrix
    if ( size(A, 1) /= size(A, 2) ) then
        call error('util.inv','Matrix must be square')
    end if
    
    ! Allocate arrays
    N = size(A, 1)
    allocate( Ainv(N, N) )
    Ainv = A
    allocate( WORK(N) )
    allocate( IPIV(N) )
    
    ! Do LU decomposition
    call DGETRF( N, N, Ainv, N, IPIV, INFO )
    if ( INFO < 0 ) then
        call string_splice(err_str, 'Lapak error: The ', INFO,                           &
        '-th argument had an illegal value')
        call error('util.inv', err_str)
    else if ( INFO > 0 ) then
        call string_splice(err_str, 'Lapack error: U(', INFO, ',', INFO, ') is '         &
        // 'exactly zero. The factorization has been completed, but the factor U is '    &
        // 'exactly singular, and division by zero will occur if it is used to solve a ' &
        // 'system  of equations.')
        call error('util.inv', err_str)
    end if
    
    ! Query optimal workspace size
    call DGETRI( N, Ainv, N, IPIV, WORK, -1, INFO )
    if ( INFO < 0 ) then
        call string_splice(err_str, 'Lapak error: The ', INFO,                           &
        '-th argument had an illegal value')
        call error('util.inv', err_str)
    else if ( INFO > 0 ) then
        call string_splice(err_str, 'Lapack error: U(', INFO, ',', INFO, ') is '         &
        // 'exactly zero; the matrix is singular and its inverse could not be computed.')
        call error('util.inv', err_str)
    end if
    LWORK = INT(WORK(1))
    deallocate(WORK)
    allocate( WORK(LWORK) )
    
    ! Invert
    call DGETRI( N, Ainv, N, IPIV, WORK, LWORK, INFO )
    if ( INFO < 0 ) then
        call string_splice(err_str, 'Lapak error: The ', INFO,                           &
        '-th argument had an illegal value')
        call error('util.inv', err_str)
    else if ( INFO > 0 ) then
        call string_splice(err_str, 'Lapack error: U(', INFO, ',', INFO, ') is '         &
        // 'exactly zero; the matrix is singular and its inverse could not be computed.')
        call error('util.inv', err_str)
    end if

end function inverse


! Cholesky factorization. This is a wrapper for Lapack. A must be a symmetric and positive 
! definite
function chol(A) result(Ac)
    implicit none
    real(rprec), dimension(:,:), intent(in)  :: A
    real(rprec), dimension(:,:), allocatable :: Ac
    integer                                  :: N, INFO
    character(CHAR_BUFF_LENGTH)              :: err_str

    N = size(A,1)
    allocate ( Ac(N, N) )
    Ac = A
    
    call DPOTRF( 'U', N, Ac, N, INFO )
    
    if ( INFO < 0 ) then
        call string_splice(err_str, 'Lapak error: The ', INFO,                           &
        '-th argument had an illegal value')
        call error('util.chol', err_str)
    else if ( INFO > 0 ) then
        call string_splice(err_str, 'Lapack error: the leading minor of order ', INFO,   &
        'is not positive definite, and the factorization could not be completed.')
        call error('util.inv', err_str)
    end if
    
end function chol

! Compute outer product between two vectors
function outer_product(v1, v2) result(M)
    implicit none
    real(rprec), intent(in), dimension(:)    :: v1, v2
    real(rprec), dimension(:,:), allocatable :: M
    
    allocate( M(size(v1), size(v2)) )
    
    M = spread(v1, DIM=2, NCOPIES=size(v2) ) * &
        spread(v2, DIM=1, NCOPIES=size(v1) )
end function outer_product

! Random number generator seeding subroutine from gcc
! https://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html#RANDOM_005fSEED
subroutine init_random_seed()
    use iso_fortran_env, only: int64
    implicit none
    integer, allocatable :: seed(:)
    integer :: i, n, un, istat, dt(8), pid
    integer(int64) :: t

    call random_seed(size = n)
    allocate(seed(n))
    ! First try if the OS provides a random number generator
    open(newunit=un, file="/dev/urandom", access="stream", &
         form="unformatted", action="read", status="old", iostat=istat)
    if (istat == 0) then
       read(un) seed
       close(un)
    else
       ! Fallback to XOR:ing the current time and pid. The PID is
       ! useful in case one launches multiple instances of the same
       ! program in parallel.
       call system_clock(t)
       if (t == 0) then
          call date_and_time(values=dt)
          t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
               + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
               + dt(3) * 24_int64 * 60 * 60 * 1000 &
               + dt(5) * 60 * 60 * 1000 &
               + dt(6) * 60 * 1000 + dt(7) * 1000 &
               + dt(8)
       end if
       pid = getpid()
       t = ieor(t, int(pid, kind(t)))
       do i = 1, n
          seed(i) = lcg(t)
       end do
    end if
    call random_seed(put=seed)
contains
    ! This simple PRNG might not be good enough for real work, but is
    ! sufficient for seeding a better PRNG.
    function lcg(s)
      integer :: lcg
      integer(int64) :: s
      if (s == 0) then
         s = 104729
      else
         s = mod(s, 4294967296_int64)
      end if
      s = mod(s * 279470273_int64, 4294967291_int64)
      lcg = int(mod(s, int(huge(0), int64)), kind(0))
    end function lcg
end subroutine init_random_seed

! This is taken from http://www.netlib.org/random/random.f90
! Adapted from the following Fortran 77 code
!      ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.
!      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
!      VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.

!  The function random_normal() returns a normally distributed pseudo-random
!  number with zero mean and unit variance.

!  The algorithm uses the ratio of uniforms method of A.J. Kinderman
!  and J.F. Monahan augmented with quadratic bounding curves.
FUNCTION random_normal() result(nr)
    implicit none
    real(rprec) :: nr
    real(rprec), parameter  :: s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472
    real(rprec), parameter  :: r1 = 0.27597, r2 = 0.27846
    real(rprec)             :: u, v, x, y, q

    ! Generate P = (u,v) uniform in rectangle enclosing acceptance region
    do 
        call random_number(u)
        call random_number(v)
        v = 1.7156 * (v - 0.5)

        ! Evaluate the quadratic form
        x = u - s
        y = abs(v) - t
        q = x**2 + y*(a*y - b*x)

        ! Accept P if inside inner ellipse
        if (q < r1) exit
        ! Reject P if outside outer ellipse
        if (q > r2) cycle
        ! Reject P if outside acceptance region
        if (v**2 < -4.0 * log(u) * u**2) exit
    end do

    ! Return ratio of P's coordinates as the normal deviate
    nr = v/u
end function random_normal

end module util