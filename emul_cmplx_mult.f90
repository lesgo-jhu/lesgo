!**********************************************************************
module emul_cmplx_mult
!**********************************************************************
implicit none

save

public

contains

!**********************************************************************
subroutine emul_cmplx_mult_rci( a, a_c, b )
!**********************************************************************
!  This subroutine emulates the multiplication of two complex scalars
!  by emulating the input real vector (a) as a complex type. This 
!  subroutine ignores the real part of a_c (e.g. would use this when
!  real(a_c) = 0). The results from the multplication are output as a.
!
!  Input:
!  
!    a (real,size(2,1))  - input real vector
!    a_c (real)          - input imaginary part of complex scalar
!
!  Output:
!
!    b (real, size(2,1)) - output real vector

use types, only : rprec
implicit none

real(rprec), dimension(2), intent(in) :: a
real(rprec), intent(in) :: a_c
real(rprec), dimension(2), intent(out) :: b

!  Cached variables
real(rprec) :: a_c_i
  
!  Cache multi-usage variables
a_c_i = a_c

!  Perform multiplication
b(1) = - a(2) * a_c_i
b(2) =  a(1) * a_c_i

return

end subroutine emul_cmplx_mult_rci

!**********************************************************************
subroutine emul_cmplx_mult_inpl_rci( a, a_c )
!**********************************************************************
!  This subroutine emulates the multiplication of two complex scalars
!  by emulating the input real vector (a) as a complex type. This 
!  subroutine ignores the real part of a_c (e.g. would use this when
!  real(a_c) = 0). The results from the multplication are output as a.
!
!  Input:
!  
!    a (real,size(2,1))           - input/outpu real array
!    a_c (real)         - input imaginary part of complex array 
!

use types, only : rprec
implicit none

real(rprec), dimension(2), intent(inout) :: a
real(rprec), intent(in) :: a_c

!  Cached variables
real(rprec) :: a_r, a_i, a_c_i
  
!  Cache multi-usage variables
a_r = a(1)
a_i = a(2)
a_c_i = a_c

!  Perform multiplication
a(1) = - a_i * a_c_i
a(2) =  a_r * a_c_i

return

end subroutine emul_cmplx_mult_inpl_rci

!**********************************************************************
subroutine emul_cmplx_mult_rc_2D( a, a_c, nx_r, nx_c, ny, b )
!**********************************************************************
!  This subroutine emulates the multiplication of two complex 2D array
!  by emulating the input real array (a) as a complex type.
!
!  Input:
!  
!    a (real,size(nx_r,ny))     - input real array
!    a_c (complex,size(nx_c,ny))- input complex array 
!    nx_r (integer)             - size of real x dimension
!    nx_c (integer)             - size of complex x dimension (nx_c must be nx_r/2)
!    ny (integer)               - size of y dimension
!
!  Output:
! 
!    b (real,size(nx_r,ny))     - output real array
!

use types, only : rprec
implicit none

integer, intent(in) :: nx_r, nx_c, ny
real(rprec), dimension( nx_r, ny ), intent(in) :: a
complex(rprec), dimension( nx_c, ny), intent(in) :: a_c

real(rprec), dimension(nx_r, ny), intent(out) :: b

!  Cached variables
real(rprec) ::  a_r, a_i, a_c_r, a_c_i

integer :: i,j,ir,ii

!  Emulate complex multiplication
!  Using outer loop to get contingious memory access
do j=1, ny
  do i=1,nx_c

    !  Real and imaginary indicies of a
    ii = 2*i
    ir = ii-1
  
    !  Cache multi-usage variables
    a_r = a(ir,j)
    a_i = a(ii,j)
    a_c_r = real(a_c(i,j),kind=rprec)
    $if($DBLPREC)
    a_c_i = dimag(a_c(i,j))
    $else
    a_c_i = aimag(a_c(i,j))
    $endif
    
    !  Perform multiplication
    b(ir,j) = a_r * a_c_r - a_i * a_c_i
    b(ii,j) = a_r * a_c_i + a_i * a_c_r

  enddo
enddo

return

end subroutine emul_cmplx_mult_rc_2D

!**********************************************************************
subroutine emul_cmplx_mult_inpl_rc_2D( a, a_c, nx_r, nx_c, ny )
!**********************************************************************
!
!  This subroutine emulates the multiplication of two complex 2D array
!  by emulating the input real array (a) as a complex type. The resulting
!  multiplication is returned as array a
!
!  Input:
!  
!    a (real,size(nx_r,ny))     - input/output real array
!    a_c (complex,size(nx_c,ny))- input complex array 
!    nx_r (integer)             - size of real x dimension
!    nx_c (integer)             - size of complex x dimension (nx_c must be nx_r/2)
!    ny (integer)               - size of y dimension
!

use types, only : rprec
implicit none

integer, intent(in) :: nx_r, nx_c, ny
real(rprec), dimension( nx_r, ny ), intent(inout) :: a
complex(rprec), dimension( nx_c, ny), intent(in) :: a_c

!  Cached variables
real(rprec) ::  a_r, a_i, a_c_r, a_c_i

integer :: i,j,ir,ii

!  Emulate complex multiplication
!  Using outer loop to get contingious memory access
do j=1, ny
  do i=1,nx_c

    !  Real and imaginary indicies of a
    ii = 2*i
    ir = ii-1
  
    !  Cache multi-usage variables
    a_r = a(ir,j)
    a_i = a(ii,j)

    a_c_r = real(a_c(i,j),kind=rprec)
    $if($DBLPREC)
    a_c_i = dimag(a_c(i,j))
    $else
    a_c_i = aimag(a_c(i,j))
    $endif
    
    !  Perform multiplication
    a(ir,j) = a_r * a_c_r - a_i * a_c_i
    a(ii,j) = a_r * a_c_i + a_i * a_c_r

  enddo
enddo

return

end subroutine emul_cmplx_mult_inpl_rc_2D

!**********************************************************************
subroutine emul_cmplx_mult_rci_2D( a, a_c, nx_r, nx_c, ny, b )
!**********************************************************************
!  This subroutine emulates the multiplication of two complex 2D array
!  by emulating the input real array (a) as a complex type. This 
!  subroutine ignores the real part of a_c (e.g. would use this when
!  real(a_c) = 0)
!
!  Input:
!  
!    a (real,size(nx_r,ny))     - input real array
!    a_c (real,size(nx_c,ny))   - input imaginary part of complex array 
!    nx_r (integer)             - size of real x dimension
!    nx_c (integer)             - size of complex x dimension (nx_c must be nx_r/2)
!    ny (integer)               - size of y dimension
!
!  Output:
! 
!    b (real,size(nx_r,ny))     - output real array
!

use types, only : rprec
implicit none

integer, intent(in) :: nx_r, nx_c, ny
real(rprec), dimension(nx_r, ny), intent(in) :: a
real(rprec), dimension(nx_c, ny), intent(in) :: a_c

real(rprec), dimension(nx_r, ny), intent(out) :: b

!  Cached variables
real(rprec) ::  a_c_i

integer :: i,j,ii,ir

!  Emulate complex multiplication
do j=1, ny !  Using outer loop to get contingious memory access
  do i=1,nx_c

    !  Real and imaginary indicies of a
    ii = 2*i
    ir = ii-1
  
    !  Cache multi-usage variables
    a_c_i = a_c(i,j)

    !  Perform multiplication
    b(ir,j) = - a(ii,j) * a_c_i
    b(ii,j) =  a(ir,j) * a_c_i

  enddo
enddo

return

end subroutine emul_cmplx_mult_rci_2D

!**********************************************************************
subroutine emul_cmplx_mult_inpl_rci_2D( a, a_c, nx_r, nx_c, ny )
!**********************************************************************
!  This subroutine emulates the multiplication of two complex 2D array
!  by emulating the input real array (a) as a complex type. This 
!  subroutine ignores the real part of a_c (e.g. would use this when
!  real(a_c) = 0). The results from the multplication are output as a.
!
!  Input:
!  
!    a (real,size(nx_r,ny))     - input/outpu real array
!    a_c (real,size(nx_c,ny))   - input imaginary part of complex array 
!    nx_r (integer)             - size of real x dimension
!    nx_c (integer)             - size of complex x dimension (nx_c must be nx_r/2)
!    ny (integer)               - size of y dimension
!


use types, only : rprec
implicit none

integer, intent(in) :: nx_r, nx_c, ny
real(rprec), dimension(nx_r, ny), intent(inout) :: a
real(rprec), dimension(nx_c, ny), intent(in) :: a_c

!  Cached variables
real(rprec) :: a_r, a_i, a_c_i

integer :: i,j,ii,ir

!  Emulate complex multiplication
do j=1, ny !  Using outer loop to get contingious memory access
  do i=1,nx_c

    !  Real and imaginary indicies of a
    ii = 2*i
    ir = ii-1
  
    !  Cache multi-usage variables
    a_r = a(ir,j)
    a_i = a(ii,j)
    a_c_i = a_c(i,j)

    !  Perform multiplication
    a(ir,j) = - a_i * a_c_i
    a(ii,j) =  a_r * a_c_i

  enddo
enddo

return

end subroutine emul_cmplx_mult_inpl_rci_2D

!**********************************************************************
subroutine emul_cmplx_mult_rcr_2D( a, a_c, nx_r, nx_c, ny, b )
!**********************************************************************
!  This subroutine emulates the multiplication of two complex 2D array
!  by emulating the input real array (a) as a complex type. This 
!  subroutine ignores the imaginary part of a_c (e.g. would use this when
!  imag(a_c) = 0). 
!
!  Input:
!  
!    a (real,size(nx_r,ny))     - input/output real array
!    a_c (real,size(nx_c,ny))   - input real part of complex array 
!    nx_r (integer)             - size of real x dimension
!    nx_c (integer)             - size of complex x dimension (nx_c must be nx_r/2)
!    ny (integer)               - size of y dimension
!
!  Output:
!
!    b (real, size(nx_r,ny))    - output real array
!

use types, only : rprec
implicit none

integer, intent(in) :: nx_r, nx_c, ny
real(rprec), dimension(nx_r, ny), intent(inout) :: a
real(rprec), dimension(nx_c, ny), intent(in) :: a_c

real(rprec), dimension(nx_r, ny), intent(out) :: b

integer :: i,j,ii,ir

!  Emulate complex multiplication
do j=1, ny !  Using outer loop to get contingious memory access
  do i=1,nx_c

    !  Real and imaginary indicies of a
    ii = 2*i
    ir = ii-1

    !  Perform multiplication
    b(ir:ii,j) = a_c(i,j)*a(ir:ii,j)

  enddo
enddo

return

end subroutine emul_cmplx_mult_rcr_2D

!**********************************************************************
subroutine emul_cmplx_mult_inpl_rcr_2D( a, a_c, nx_r, nx_c, ny )
!**********************************************************************
!  This subroutine emulates the multiplication of two complex 2D array
!  by emulating the input real array (a) as a complex type. This 
!  subroutine ignores the imaginary part of a_c (e.g. would use this when
!  imag(a_c) = 0). The results from the multplication are output as a.
!
!  Input:
!  
!    a (real,size(nx_r,ny))     - input/output real array
!    a_c (real,size(nx_c,ny))   - input real part of complex array 
!    nx_r (integer)             - size of real x dimension
!    nx_c (integer)             - size of complex x dimension (nx_c must be nx_r/2)
!    ny (integer)               - size of y dimension
!

use types, only : rprec
implicit none

integer, intent(in) :: nx_r, nx_c, ny
real(rprec), dimension(nx_r, ny), intent(inout) :: a
real(rprec), dimension(nx_c, ny), intent(in) :: a_c

integer :: i,j,ii,ir

!  Emulate complex multiplication
do j=1, ny !  Using outer loop to get contingious memory access
  do i=1,nx_c

    !  Real and imaginary indicies of a
    ii = 2*i
    ir = ii-1

    !  Perform multiplication
    a(ir:ii,j) = a_c(i,j)*a(ir:ii,j)

  enddo
enddo

return

end subroutine emul_cmplx_mult_inpl_rcr_2D
end module emul_cmplx_mult
