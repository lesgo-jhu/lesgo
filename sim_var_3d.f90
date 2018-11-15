module sim_var_3d_m

use iso_c_binding
! use defs
use types, only : rprec, cprec, MPI_RPREC, MPI_CPREC
use grid_m
implicit none
include 'fftw3.f03'

private
public :: sim_var_3d_t, UV_GRID, W_GRID

logical, parameter :: UV_GRID = .true.
logical, parameter :: W_GRID = .false.

type :: sim_var_3d_t
    integer(c_intptr_t), private :: alloc_size
    type(c_ptr), private :: fftw_data, forward_plan, backward_plan
    real(rprec), pointer, dimension(:,:,:), private :: rdata
    real(rprec), pointer, dimension(:), private :: rdata_1d
    complex(cprec), pointer, dimension(:,:,:), private :: cdata
    complex(cprec), pointer, dimension(:), private :: cdata_1d
    real(rprec), pointer, dimension(:,:,:), public :: real
    complex(cprec), pointer, dimension(:,:,:), public :: cmplx
    type(grid_t), pointer, private :: grid
    logical, public :: uvw_grid
contains
    procedure, public :: forward => forward_inplace, forward_outofplace
    ! procedure, public :: forward_inplace, forward_outofplace
    procedure, public :: backward
    procedure, public :: zero_nyquist
    procedure, public :: ddx, ddy, ddxy, filt_ddxy, ddz
    procedure, public :: test_filter
    procedure, public :: padd, unpadd
    procedure, public :: sync_down, sync_up, sync_downup
    final :: destructor
end type sim_var_3d_t

interface sim_var_3d_t
    module procedure constructor
end interface sim_var_3d_t

contains

!*******************************************************************************
function constructor(grid, uvw_grid) result(this)
!*******************************************************************************
! Constructor for sim_var_3d_t
!
type(sim_var_3d_t) :: this
logical, intent(in) :: uvw_grid
! integer(c_intptr_t) :: Nx_cint, Nkx_cint, Ny_cint
type(grid_t), target :: grid

! Set input parameters
this%grid => grid
this%uvw_grid = uvw_grid

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
this%fftw_data = fftw_alloc_complex(this%alloc_size*(this%grid%Nz+1))

! Fortran points to data that is properly typed
call c_f_pointer(this%fftw_data, this%rdata_1d,                                &
    [2*this%alloc_size*(this%grid%Nz+1)])
call c_f_pointer(this%fftw_data, this%cdata_1d,                                &
    [this%alloc_size*(this%grid%Nz+1)])

! Initialize to 0
this%rdata_1d(:) = 0.0_rprec

! Pointers in 2D
this%rdata(1:2*this%grid%Nkx,1:this%alloc_size/this%grid%Nkx, 0:this%grid%Nz)     &
    => this%rdata_1d
this%cdata(1:this%grid%Nkx,1:this%alloc_size/this%grid%Nkx, 0:this%grid%Nz)       &
    => this%cdata_1d

! Public access pointers that hide padding
this%real(1:,1:,0:) => this%rdata(1:this%grid%Nx, 1:this%grid%Ny, 0:this%grid%Nz)
this%cmplx(1:,1:,0:) => this%cdata(1:this%grid%Nkx, 1:this%grid%Ny, 0:this%grid%Nz)

! Create plans
this%forward_plan = fftw_plan_dft_r2c_2d(this%grid%Ny, this%grid%Nx,             &
    this%rdata(:,:,0), this%cdata(:,:,0), FFTW_MEASURE)
this%backward_plan = fftw_plan_dft_c2r_2d(this%grid%Ny, this%grid%Nx,            &
    this%cdata(:,:,0), this%rdata(:,:,0), FFTW_MEASURE)

end function constructor

!*******************************************************************************
subroutine forward_inplace(this)
!*******************************************************************************
class(sim_var_3d_t) :: this
integer :: i

do i = 0, this%grid%Nz
    call fftw_execute_dft_r2c(this%forward_plan, this%rdata(:,:,i),        &
        this%cdata(:,:,i))
end do

end subroutine forward_inplace

!*******************************************************************************
subroutine forward_outofplace(this, that)
!*******************************************************************************
class(sim_var_3d_t), intent(inout) :: this, that
integer :: i

do i = 0, this%grid%Nz
    call fftw_execute_dft_r2c(this%forward_plan, this%rdata(:,:,i),        &
        that%cdata(:,:,i))
end do

end subroutine forward_outofplace

!*******************************************************************************
subroutine backward(this)
!*******************************************************************************
class(sim_var_3d_t) :: this
integer :: i

do i = 0, this%grid%Nz
    call fftw_execute_dft_c2r(this%backward_plan, this%cdata(:,:,i),       &
        this%rdata(:,:,i))
end do
this%rdata = this%rdata/this%grid%Ny/this%grid%Nx

end subroutine backward

!*******************************************************************************
subroutine zero_nyquist(this)
!*******************************************************************************
class(sim_var_3d_t) :: this
integer :: i

do i = 1, size(this%grid%nyquist(:,1))
    this%cmplx(this%grid%nyquist(i,1), this%grid%nyquist(i,2),:) = 0._cprec
end do

end subroutine zero_nyquist

!*******************************************************************************
subroutine ddx(this, dudx)
!*******************************************************************************
class(sim_var_3d_t), intent(inout) :: this
class(sim_var_3d_t), intent(inout) :: dudx
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
    dudx%cmplx(i,:,:) = dudx%cmplx(i,:,:) * cmplx(0, this%grid%kx(i), kind=cprec)
end do
call dudx%zero_nyquist()
call dudx%backward()

end subroutine ddx

!*******************************************************************************
subroutine ddy(this, dudy)
!*******************************************************************************
class(sim_var_3d_t), intent(inout) :: this
class(sim_var_3d_t), intent(inout) :: dudy

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
    dudy%cmplx(:,i,:) = dudy%cmplx(:,i,:) * cmplx(0, this%grid%ky(i), kind=cprec)
end do
call dudy%zero_nyquist()
call dudy%backward()

end subroutine ddy

!*******************************************************************************
subroutine ddxy(this, dudx, dudy)
!*******************************************************************************
class(sim_var_3d_t), intent(inout) :: this
class(sim_var_3d_t), intent(inout) :: dudx, dudy

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
    dudy%cmplx(:,i,:) = dudx%cmplx(:,i,:) * cmplx(0, this%grid%ky(i), kind=cprec)
end do
do i = 1, this%grid%Nkx
    dudx%cmplx(i,:,:) = dudx%cmplx(i,:,:) * cmplx(0, this%grid%kx(i), kind=cprec)
end do
call dudx%zero_nyquist()
call dudy%zero_nyquist()
call dudx%backward()
call dudy%backward()

end subroutine ddxy

!*******************************************************************************
subroutine filt_ddxy(this, dudx, dudy)
!*******************************************************************************
class(sim_var_3d_t), intent(inout) :: this
class(sim_var_3d_t), intent(inout) :: dudx, dudy

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
    dudy%cmplx(:,i,:) = this%cmplx(:,i,:) * cmplx(0, this%grid%ky(i), kind=cprec)
end do
do i = 1, this%grid%Nkx
    dudx%cmplx(i,:,:) = this%cmplx(i,:,:) * cmplx(0, this%grid%kx(i), kind=cprec)
end do

! Return variables to physical space
call this%backward()
call dudx%backward()
call dudy%backward()

end subroutine filt_ddxy

!*******************************************************************************
subroutine ddz(this, dudz)
!*******************************************************************************
!
! This subroutine computes the partial derivative of f with respect to z using
! 2nd order finite differencing.
!
! If variable is on the uv grid, then derivative is on w grid. On bottom
! processor, provides 2:nz. On top processor provides 1:nz-1. All others provide
! 1:nz.

! If variable is on the w grid, then derivative is on the uv grid. On the top
! processor, provides dudz(1:nz-1). On all other processors, provides
! dudz(0:nz-1)
!
use param, only : BOGUS
class(sim_var_3d_t), intent(inout) :: this
class(sim_var_3d_t), intent(inout) :: dudz
integer :: i

! Check that input arguments are on opposing grids
if (this%uvw_grid .eqv. dudz%uvw_grid) then
    write(*,*) "cannot calculate ddz because variables are on the same grid"
    stop
end if

if (this%uvw_grid .eqv. UV_GRID) then

#if defined(PPMPI) && defined(PPSAFETYMODE)
    dudz%real(:,:,0) = BOGUS
#endif

    ! Calculate derivative.
    ! The ghost node information is available here
    ! if coord == 0, dudz(1) will be set in wallstress
    do i = 1, this%grid%nz
        dudz%real(:,:,i) = (this%real(:,:,i) - this%real(:,:,i-1))/this%grid%dz
    end do

#ifdef PPSAFETYMODE
    ! Not necessarily accurate at top and bottom boundary
    ! Set to BOGUS just to be safe
    if (this%grid%coord == 0) then
        dudz%real(:,:,1) = BOGUS
    end if
    if (this%grid%coord == this%grid%nproc-1) then
        dudz%real(:,:,this%grid%nz) = BOGUS
    end if
#endif

else

    do i = 0, this%grid%nz-1
        dudz%real(:,:,i) = (this%real(:,:,i+1) - this%real(:,:,i))/this%grid%dz
    end do

#ifdef PPSAFETYMODE
    ! bottom process cannot calculate dudz(0)
    if (this%grid%coord == 0) then
        dudz%real(:,:,0) = BOGUS
    endif

    ! All processes cannot calculate dudz(nz)
    dudz%real(:,:,this%grid%nz) = BOGUS
#endif

end if

end subroutine ddz

!*******************************************************************************
subroutine test_filter(this, f)
!*******************************************************************************
class(sim_var_3d_t), intent(inout) :: this
class(sim_var_3d_t), intent(inout) :: f

integer :: i

! Check that derivative array is associated with the same grid
if ( .not. associated(this%grid, target=f%grid) ) then
    write(*,*) "test filter variable is not associated with the same grid"
    stop
end if

! Copy complex values to dudx
f%real = this%real

! Perform test filter calulations in spectral space
call f%forward()
do i = 0, this%grid%Nz
    f%cmplx(:,:,i) = f%cmplx(:,:,i)*this%grid%G_test
end do
call f%backward()

end subroutine test_filter

! !*******************************************************************************
! subroutine padd (u_big,u)
! !*******************************************************************************
! ! puts arrays into larger, zero-padded arrays
! ! automatically zeroes the oddballs
! use types, only : rprec
! use param, only : ld,ld_big,nx,ny,ny2
! implicit none
!
! !  u and u_big are interleaved as complex arrays
! real(rprec), dimension(ld,ny), intent(in) :: u
! real(rprec), dimension(ld_big,ny2), intent(out) :: u_big
!
! integer :: ny_h, j_s, j_big_s
!
! ny_h = ny/2
!
! ! make sure the big array is zeroed!
! u_big(:,:) = 0._rprec
!
! ! note: split access in an attempt to maintain locality
! u_big(:nx,:ny_h) = u(:nx,:ny_h)
!
! ! Compute starting j locations for second transfer
! j_s = ny_h + 2
! j_big_s = ny2 - ny_h + 2
!
! u_big(:nx,j_big_s:ny2) = u(:nx,j_s:ny)
!
! end subroutine padd

!*******************************************************************************
subroutine padd(this, this_big)
!*******************************************************************************
class(sim_var_3d_t), intent(inout) :: this, this_big
integer :: nx_h, ny_h, j_s, j_big_s
real(rprec) :: const

! Go to spectral space
call this%forward()

! Make sure big variable is zeroed
this_big%cmplx(:,:,:) = cmplx(0._rprec, 0._rprec, cprec)

! First set of wavenumbers
nx_h = this%grid%Nkx
ny_h = this%grid%Ny/2
const = (this_big%grid%Ny*this_big%grid%Nx)
const = const/(this%grid%Ny*this%grid%Nx)
this_big%cmplx(1:nx_h,1:ny_h,:) = const*this%cmplx(1:nx_h,1:ny_h,:)

! Second set of wavenumbers
j_s = ny_h + 2
j_big_s = this_big%grid%Ny - ny_h + 2
this_big%cmplx(1:nx_h,j_big_s:this_big%grid%Ny,:) =                            &
    const*this%cmplx(1:nx_h,j_s:this%grid%Ny,:)

! Make sure nyquist is not included
call this_big%zero_nyquist()

! Go back to physical space
call this%backward
call this_big%backward

end subroutine padd

!*******************************************************************************
subroutine unpadd(this_big, this)
!*******************************************************************************
class(sim_var_3d_t), intent(inout) :: this_big, this
integer :: nx_h, ny_h, j_s, j_big_s
real(rprec) :: const

! Go to spectral space
call this_big%forward()

! Make sure  variable is zeroed
this%cmplx(:,:,:) = cmplx(0._rprec, 0._rprec, cprec)

! First set of wavenumbers
nx_h = this%grid%Nkx
ny_h = this%grid%Ny/2
const = (this%grid%Ny*this%grid%Nx)
const = const/(this_big%grid%Ny*this_big%grid%Nx)
this%cmplx(1:nx_h,1:ny_h,:) = const*this_big%cmplx(1:nx_h,1:ny_h,:)

! Second set of wavenumbers
j_s = ny_h + 2
j_big_s = this_big%grid%Ny - ny_h + 2
this%cmplx(1:nx_h,j_s:this%grid%Ny,:) =                                       &
    const*this_big%cmplx(1:nx_h,j_big_s:this_big%grid%Ny,:)

! Make sure nyquist is not included
call this%zero_nyquist()

! Go back to physical space
call this%backward
call this_big%backward

end subroutine unpadd


! !*******************************************************************************
! subroutine unpadd(cc,cc_big)
! !*******************************************************************************
! use types, only : rprec
! use param, only : ld,nx,ny,ny2,ld_big
! implicit none
!
! !  cc and cc_big are interleaved as complex arrays
! real(rprec), dimension( ld, ny ) :: cc
! real(rprec), dimension( ld_big, ny2 ) :: cc_big
!
! integer :: ny_h, j_s, j_big_s
!
! ny_h = ny/2
!
! cc(:nx,:ny_h) = cc_big(:nx,:ny_h)
!
! ! oddballs
! cc(ld-1:ld,:) = 0._rprec
! cc(:,ny_h+1) = 0._rprec
!
! ! Compute starting j locations for second transfer
! j_s = ny_h + 2
! j_big_s = ny2 - ny_h + 2
! cc(:nx,j_s:ny) = cc_big(:nx,j_big_s:ny2)
!
! end subroutine unpadd

! !*******************************************************************************
! subroutine spectral_project(this, this_big)
! !*******************************************************************************
! class(sim_var_3d_t), intent(inout) :: this, this_big
!
! ! Go to spectral space
! call this%forward()
!
! ! Transfer onto larger grid
! this_big%cmplx(:,:,:) = cmplx(0._rprec, cprec)
! this_big%cmplx(1:this%grid%Nkx, 1:this%grid%Ny, :) =                           &
!     (this_big%grid%Ny*this_big%grid%Nx)*this%cmplx(1:this%grid%Nkx, 1:this%grid%Ny, :)/(this%grid%Ny*this%grid%Nx)
!
! ! Go back to physical space
! call this%backward()
! call this_big%backward()
!
! end subroutine spectral_project

!*******************************************************************************
subroutine sync_down(this)
!*******************************************************************************
class(sim_var_3d_t), intent(inout) :: this
integer :: ierr, status
! Very Kludgy---Forced gfortran not to do an optimization
if (this%grid%coord == -1) write(*,*) shape(this%cmplx(:,:,1))

call mpi_sendrecv(this%cmplx(:,:,1), this%alloc_size, MPI_CPREC,      &
    this%grid%down, 1, this%cmplx(:,:,this%grid%Nz), this%alloc_size, &
    MPI_CPREC, this%grid%up, 1, this%grid%comm, status, ierr)


end subroutine sync_down

!*******************************************************************************
subroutine sync_up(this)
!*******************************************************************************
class(sim_var_3d_t), intent(inout) :: this
integer :: ierr, status
! Very Kludgy---Forced gfortran not to do an optimization
if (this%grid%coord == -1) write(*,*) shape(this%cmplx(:,:,1))

call mpi_sendrecv(this%cmplx(:,:,this%grid%Nz-1), this%alloc_size,    &
    MPI_CPREC, this%grid%up, 2, this%cmplx(:,:,0), this%alloc_size, &
    MPI_CPREC, this%grid%down, 2, this%grid%comm, status, ierr)

end subroutine sync_up

!*******************************************************************************
subroutine sync_downup(this)
!*******************************************************************************
class(sim_var_3d_t), intent(inout) :: this

call this%sync_down()
call this%sync_up()

end subroutine sync_downup

!*******************************************************************************
subroutine destructor(this)
!*******************************************************************************
type(sim_var_3d_t) :: this

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

end subroutine destructor

end module sim_var_3d_m
