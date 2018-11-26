module sim_var_2d_m

use iso_c_binding
! use defs
use types, only : rprec, cprec, MPI_RPREC, MPI_CPREC
use grid_m
implicit none
include 'fftw3.f03'

private
public :: sim_var_2d_t

logical, parameter :: UV_GRID = .true.
logical, parameter :: W_GRID = .false.

type :: sim_var_2d_t
    integer(c_intptr_t), private :: alloc_size
    type(c_ptr), private :: fftw_data, forward_plan, backward_plan
    real(rprec), pointer, dimension(:,:), private :: rdata
    real(rprec), pointer, dimension(:), private :: rdata_1d
    complex(cprec), pointer, dimension(:,:), private :: cdata
    complex(cprec), pointer, dimension(:), private :: cdata_1d
    real(rprec), pointer, dimension(:,:), public :: real
    complex(cprec), pointer, dimension(:,:), public :: cmplx
    type(grid_t), pointer, private :: grid
contains
    procedure, public :: forward => forward_inplace, forward_outofplace
    procedure, public :: backward
    procedure, public :: zero_nyquist
    procedure, public :: ddx, ddy, ddxy, filt_ddxy
    procedure, public :: test_filter
    procedure, public :: padd, unpadd
    procedure, public :: free
end type sim_var_2d_t

interface sim_var_2d_t
    module procedure constructor
end interface sim_var_2d_t

contains

!*******************************************************************************
function constructor(grid) result(this)
!*******************************************************************************
! Constructor for sim_var_2d_t
!
type(sim_var_2d_t) :: this
! integer(c_intptr_t) :: Nx_cint, Nkx_cint, Ny_cint
type(grid_t), target :: grid

! Set input parameters
this%grid => grid

! Number of complex values in x
! this%grid%Nkx = this%grid%Nkx

! Convert to c_intptr_t for FFTW calls
! Nx_cint = int(this%grid%Nx)
! Nkx_cint = int(this%grid%Nkx)
! Ny_tot_cint = int(this%grid%Ny_tot)

! Find allocation size
this%alloc_size = this%grid%Nkx*this%grid%Ny
!
! ! write(*,*) this%alloc_size, Ny_tot_cint, Ny_cint
!
! C pointer for FFTW (where the data is actually held)
this%fftw_data = fftw_alloc_complex(this%alloc_size)

! Fortran points to data that is properly typed
call c_f_pointer(this%fftw_data, this%rdata_1d, [2*this%alloc_size])
call c_f_pointer(this%fftw_data, this%cdata_1d, [this%alloc_size])

! Initialize to 0
this%rdata_1d(:) = 0.0_rprec

! Pointers in 2D
this%rdata(1:2*this%grid%Nkx,1:this%alloc_size/this%grid%Nkx) => this%rdata_1d
this%cdata(1:this%grid%Nkx,1:this%alloc_size/this%grid%Nkx) => this%cdata_1d

! Public access pointers that hide padding
this%real(1:,1:) => this%rdata(1:this%grid%Nx, 1:this%grid%Ny)
this%cmplx(1:,1:) => this%cdata(1:this%grid%Nkx, 1:this%grid%Ny)

! Create plans
this%forward_plan = fftw_plan_dft_r2c_2d(this%grid%Ny, this%grid%Nx,             &
    this%rdata(:,:), this%cdata(:,:), FFTW_MEASURE)
this%backward_plan = fftw_plan_dft_c2r_2d(this%grid%Ny, this%grid%Nx,            &
    this%cdata(:,:), this%rdata(:,:), FFTW_MEASURE)

end function constructor

!*******************************************************************************
subroutine forward_inplace(this)
!*******************************************************************************
class(sim_var_2d_t) :: this

call fftw_execute_dft_r2c(this%forward_plan, this%rdata(:,:), this%cdata(:,:))

end subroutine forward_inplace

!*******************************************************************************
subroutine forward_outofplace(this, that)
!*******************************************************************************
class(sim_var_2d_t), intent(inout) :: this, that

call fftw_execute_dft_r2c(this%forward_plan, this%rdata(:,:), that%cdata(:,:))

end subroutine forward_outofplace

!*******************************************************************************
subroutine backward(this)
!*******************************************************************************
class(sim_var_2d_t) :: this

call fftw_execute_dft_c2r(this%backward_plan, this%cdata(:,:), this%rdata(:,:))

this%rdata = this%rdata/this%grid%Ny/this%grid%Nx

end subroutine backward

!*******************************************************************************
subroutine zero_nyquist(this)
!*******************************************************************************
class(sim_var_2d_t) :: this
integer :: i

do i = 1, size(this%grid%nyquist(:,1))
    this%cmplx(this%grid%nyquist(i,1), this%grid%nyquist(i,2)) = 0._cprec
end do

end subroutine zero_nyquist

!*******************************************************************************
subroutine ddx(this, dudx)
!*******************************************************************************
class(sim_var_2d_t), intent(inout) :: this
class(sim_var_2d_t), intent(inout) :: dudx
integer :: i

! Check that derivative array is associated with the same grid
if ( .not. associated(this%grid, target=dudx%grid) ) then
    write(*,*) "derivative is not associated with the same grid"
    stop
end if

! Copy complex values to dudx
dudx%real = this%real

! Perform derivative calulation in spectral space
! \hat{dv/dx} = i*kx*\hat{v}
call dudx%forward()
do i = 1, this%grid%Nkx
    dudx%cmplx(i,:) = dudx%cmplx(i,:) * cmplx(0, this%grid%kx(i), kind=cprec)
end do
call dudx%zero_nyquist()
call dudx%backward()

end subroutine ddx

!*******************************************************************************
subroutine ddy(this, dudy)
!*******************************************************************************
class(sim_var_2d_t), intent(inout) :: this
class(sim_var_2d_t), intent(inout) :: dudy

integer :: i

! Check that derivative array is associated with the same grid
if ( .not. associated(this%grid, target=dudy%grid) ) then
    write(*,*) "derivative is not associated with the same grid"
    stop
end if

! Copy complex values to dudy
dudy%real = this%real

! Perform derivative calulation in spectral space
! \hat{dv/dy} = i*ky*\hat{v}
call dudy%forward()
do i = 1, this%grid%Ny
    dudy%cmplx(:,i) = dudy%cmplx(:,i) * cmplx(0, this%grid%ky(i), kind=cprec)
end do
call dudy%zero_nyquist()
call dudy%backward()

end subroutine ddy

!*******************************************************************************
subroutine ddxy(this, dudx, dudy)
!*******************************************************************************
class(sim_var_2d_t), intent(inout) :: this
class(sim_var_2d_t), intent(inout) :: dudx, dudy

integer :: i

! Check that derivative array is associated with the same grid
if ( .not. associated(this%grid, target=dudy%grid) ) then
    write(*,*) "derivative is not associated with the same grid"
    stop
end if

! Copy complex values to dudx
dudx%real = this%real

! Perform derivative calulations in spectral space
call dudx%forward()
do i = 1, this%grid%Ny
    dudy%cmplx(:,i) = dudx%cmplx(:,i) * cmplx(0, this%grid%ky(i), kind=cprec)
end do
do i = 1, this%grid%Nkx
    dudx%cmplx(i,:) = dudx%cmplx(i,:) * cmplx(0, this%grid%kx(i), kind=cprec)
end do
call dudx%zero_nyquist()
call dudy%zero_nyquist()
call dudx%backward()
call dudy%backward()

end subroutine ddxy

!*******************************************************************************
subroutine filt_ddxy(this, dudx, dudy)
!*******************************************************************************
class(sim_var_2d_t), intent(inout) :: this
class(sim_var_2d_t), intent(inout) :: dudx, dudy

integer :: i

! Check that derivative array is associated with the same grid
if ( .not. associated(this%grid, target=dudy%grid) ) then
    write(*,*) "derivative is not associated with the same grid"
    stop
end if

! Filter variable in spectral space
call this%forward()
call this%zero_nyquist()

! Perform derivative calulations in spectral space
do i = 1, this%grid%Ny
    dudy%cmplx(:,i) = this%cmplx(:,i) * cmplx(0, this%grid%ky(i), kind=cprec)
end do
do i = 1, this%grid%Nkx
    dudx%cmplx(i,:) = this%cmplx(i,:) * cmplx(0, this%grid%kx(i), kind=cprec)
end do

! Return variables to physical space
call this%backward()
call dudx%backward()
call dudy%backward()

end subroutine filt_ddxy

!*******************************************************************************
subroutine test_filter(this, f)
!*******************************************************************************
class(sim_var_2d_t), intent(inout) :: this
class(sim_var_2d_t), intent(inout) :: f

! Check that derivative array is associated with the same grid
if ( .not. associated(this%grid, target=f%grid) ) then
    write(*,*) "test filter variable is not associated with the same grid"
    stop
end if

! Copy complex values to dudx
f%real = this%real

! Perform test filter calulations in spectral space
call f%forward()
f%cmplx(:,:) = f%cmplx(:,:)*this%grid%G_test
call f%backward()

end subroutine test_filter

!*******************************************************************************
subroutine padd(this, this_big)
!*******************************************************************************
class(sim_var_2d_t), intent(inout) :: this, this_big
integer :: nx_h, ny_h, j_s, j_big_s
real(rprec) :: const

! Go to spectral space
call this%forward()

! Make sure big variable is zeroed
this_big%cmplx(:,:) = cmplx(0._rprec, 0._rprec, cprec)

! First set of wavenumbers
nx_h = this%grid%Nkx
ny_h = this%grid%Ny/2
const = (this_big%grid%Ny*this_big%grid%Nx)
const = const/(this%grid%Ny*this%grid%Nx)
this_big%cmplx(1:nx_h,1:ny_h) = const*this%cmplx(1:nx_h,1:ny_h)

! Second set of wavenumbers
j_s = ny_h + 2
j_big_s = this_big%grid%Ny - ny_h + 2
this_big%cmplx(1:nx_h,j_big_s:this_big%grid%Ny) =                            &
    const*this%cmplx(1:nx_h,j_s:this%grid%Ny)

! Make sure nyquist is not included
call this_big%zero_nyquist()

! Go back to physical space
call this%backward
call this_big%backward

end subroutine padd

!*******************************************************************************
subroutine unpadd(this_big, this)
!*******************************************************************************
class(sim_var_2d_t), intent(inout) :: this_big, this
integer :: nx_h, ny_h, j_s, j_big_s
real(rprec) :: const

! Go to spectral space
call this_big%forward()

! Make sure  variable is zeroed
this%cmplx(:,:) = cmplx(0._rprec, 0._rprec, cprec)

! First set of wavenumbers
nx_h = this%grid%Nkx
ny_h = this%grid%Ny/2
const = (this%grid%Ny*this%grid%Nx)
const = const/(this_big%grid%Ny*this_big%grid%Nx)
this%cmplx(1:nx_h,1:ny_h) = const*this_big%cmplx(1:nx_h,1:ny_h)

! Second set of wavenumbers
j_s = ny_h + 2
j_big_s = this_big%grid%Ny - ny_h + 2
this%cmplx(1:nx_h,j_s:this%grid%Ny) =                                          &
    const*this_big%cmplx(1:nx_h,j_big_s:this_big%grid%Ny)

! Make sure nyquist is not included
call this%zero_nyquist()

! Go back to physical space
call this%backward
call this_big%backward

end subroutine unpadd

!*******************************************************************************
subroutine free(this)
!*******************************************************************************
class(sim_var_2d_t) :: this

! Clean up plans, pointers, and data
call fftw_destroy_plan(this%forward_plan)
call fftw_destroy_plan(this%backward_plan)
nullify(this%real)
nullify(this%cmplx)
nullify(this%rdata)
nullify(this%cdata)
nullify(this%rdata_1d)
nullify(this%cdata_1d)
nullify(this%grid)
call fftw_free(this%fftw_data)

end subroutine free

end module sim_var_2d_m
