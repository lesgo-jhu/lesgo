!!
!!  Copyright (C) 2009-2016  Johns Hopkins University
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

module rh_control

use types, only : rprec
use param, only : pi
use wake_model_class
use wake_model_adjoint_class
use minimized_class
use messages
implicit none

private
public MinimizedFarm

type, extends(Minimized) :: MinimizedFarm
    type(WakeModel)         :: w       ! wake (turbine)
    type(WakeModel)         :: iw      ! wake initial condition (turbine)
    type(WakeModelAdjoint)  :: wstar   ! adjoint wake (turbine)
    type(WakeModelAdjoint)  :: iwstar  ! adjoint wake initial condition  (turbine)
    real(rprec), dimension(:,:), allocatable     :: Ctp     ! local thrust coefficient (turbine, time)
    real(rprec), dimension(:), allocatable       :: t       ! times
    real(rprec), dimension(:), allocatable       :: Pref    ! reference signal
    real(rprec), dimension(:), allocatable       :: Pfarm   ! farm power
    real(rprec)                                  :: cost = 0._rprec    ! cost
    ! Gradient is indexed with time, but only 2:end counts since Ctp is fixed at 1
    real(rprec), dimension(:,:), allocatable     :: grad    ! gradient (turbine, time)
    real(rprec), dimension(:,:), allocatable     :: fdgrad  ! finite-difference gradient (turbine, time)
    real(rprec) :: Ctp0 = 1.33_rprec
    real(rprec) :: gamma = 0.0_rprec! 0.05_rprec
    real(rprec) :: eta = 0.0_rprec! 0.005_rprec
    real(rprec) :: POWER = 0._rprec
    real(rprec) :: cfl, dt              ! Constant time step
    integer     :: N, Nt
    logical     :: isDimensionless = .false.
contains
    procedure, public :: eval
    procedure, public :: get_Ctp_vector
    procedure, public :: initialize
    procedure, public :: makeDimensionless
    procedure, public :: makeDimensional
    procedure, public :: finiteDifferenceGradient
    procedure, public :: run => run_input
    procedure, private :: run_noinput
end type MinimizedFarm

interface MinimizedFarm
    procedure :: constructor
end interface MinimizedFarm

contains

function constructor(i_wm, i_t0, i_T, i_cfl, i_time, i_Pref) result(this)
    implicit none
    type(MinimizedFarm)                         :: this
    class(WakeModel), intent(in)                :: i_wm
    real(rprec), dimension(:), intent(in)       :: i_time, i_Pref
    real(rprec), intent(in)                     :: i_t0, i_T, i_cfl
    
    call this%initialize(i_wm, i_t0, i_T, i_cfl, i_time, i_Pref) 
end function constructor

subroutine initialize(this, i_wm, i_t0, i_T, i_cfl, i_time, i_Pref)
    use util, only : interpolate
    implicit none
    class(MinimizedFarm)                        :: this
    class(WakeModel), intent(in)                :: i_wm
    real(rprec), dimension(:), intent(in)       :: i_time, i_Pref
    real(rprec), intent(in)                     :: i_t0, i_T, i_cfl
    integer                                     :: i
        
    ! Set initial condition for wake model
    this%iw = i_wm
    call this%iw%makeDimensional
    this%w = this%iw
    this%N = this%w%N

    ! Create adjoint wake models
    this%iwstar = WakeModelAdjoint(this%w%s, this%w%U_infty, this%w%Delta, this%w%k, this%w%Dia, this%w%Nx)
    this%wstar = this%iwstar
    
    ! Dimension for power
    this%POWER = this%w%VELOCITY**3 
    
    ! Create time for the time horizon
    this%cfl      = i_cfl
    this%dt       = this%cfl * this%w%dx / this%w%U_infty
    this%Nt       = ceiling(i_T / this%dt)
    allocate( this%t(this%Nt) )
    do i = 1, this%Nt
        this%t(i) = i_t0 + this%dt * (i - 1)
    end do
    
    ! Interpolate the reference signal
    allocate( this%Pref(this%Nt) )
    allocate( this%Pfarm(this%Nt) )
    call interpolate(i_time, i_Pref, this%t, this%Pref)
    write(*,*) i_time
    write(*,*) i_Pref
    write(*,*) this%t
    write(*,*) this%Pref
    write(*,*) this%POWER
    
    ! Allocate other variables
    allocate( this%Ctp(this%N, this%Nt)    )
    allocate( this%grad(this%N, this%Nt)   )
    allocate( this%fdgrad(this%N, this%Nt) )
    
    
end subroutine initialize

subroutine run_input(this, i_t, i_Ctp)
    use util, only : interpolate
    implicit none
    class(MinimizedFarm), intent(inout)        :: this
    real(rprec), dimension(:), intent(in)      :: i_t
    real(rprec), dimension(:,:), intent(in)    :: i_Ctp
    integer                                    :: n
    
    ! Make dimensionless
    call this%makeDimensionless
    
    ! Interpolate input onto object
    do n = 1, this%N
        call interpolate(i_t, i_Ctp(n,:), this%t(2:) *  this%w%TIME, this%Ctp(n,2:))
    end do
    this%Ctp(:,1) = this%Ctp0
    
    call this%run_noinput
end subroutine run_input

subroutine run_noinput(this)
    use util, only : interpolate
    implicit none
    class(MinimizedFarm), intent(inout)        :: this
    real(rprec), dimension(:), allocatable     :: u
    integer                                    :: i, n, k
    real(rprec), dimension(:), allocatable     :: du_super
    real(rprec), dimension(:), allocatable     :: uhatstar, ustar ! adjoint forcing terms
    real(rprec), dimension(:,:,:), allocatable :: fstar           ! adjoint forcing (turbine, time, space)
    real(rprec), dimension(:), allocatable     :: g

    ! Allocate some variables
    allocate( du_super(this%w%Nx) )
    allocate( uhatstar(this%N) )
    allocate( ustar(this%w%Nx) )
    allocate( fstar(this%N, this%Nt, this%w%Nx) )
    
    ! Make dimensionless
    call this%makeDimensionless
        
    ! Reset cost and gradient
    this%cost = 0
    this%grad = 0
    this%Pfarm = 0
    fstar = 0.d0
    
    ! Assign initial conditions
    this%w = this%iw
    this%wstar = this%iwstar

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Forward Simulation
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! Iterate in time
    do k = 2, this%Nt
        ! Calculate next step
        call this%w%advance(this%Ctp(:, k), this%dt)
        
        ! Calculate farm power
        this%Pfarm(k) = sum(this%w%Phat)
        
        ! Calculate contribution to cost function
        this%cost = this%cost + (this%Pfarm(k) - this%Pref(k))**2 * this%dt
        
        ! Calculate contribution to the gradient (dJ_1/dCt'_n) at time k
        do n = 1, this%N
          this%grad(n, k) = 2.d0 * (this%Pfarm(k) - this%Pref(k)) * this%w%uhat(n)**3 * this%dt
        end do
        
        ! Calculate adjoint forcing
        uhatstar = -6.d0 * (this%Pfarm(k) - this%Pref(k)) * this%Ctp(:,k) * this%w%uhat**2
        ustar = 0
        do n = 1, this%N
            ustar = ustar + this%w%G(n,:) * uhatstar(n)
        end do
        du_super = this%w%U_infty - this%w%u
        do n = 1, this%N
            fstar(n, k, :) = - ustar * this%w%du(n,:) / du_super
            do i = 1, this%w%Nx
                if ( du_super(i) <= 1E-10 )   fstar(n, k, i) = 0._rprec
            end do
        end do
        
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Backward (Adjoint) Simulation
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate( g(this%N) )
    do k = this%Nt-1, 1, -1
        g = 0._rprec
        call this%wstar%retract(fstar(:,k+1,:), this%dt, g)
        do n = 1, this%N
            this%grad(n, k) = this%grad(n, k) - g(n) * 8.d0 * this%w%U_infty**2 &
            / ( ( 4.d0 + this%Ctp(n,k))**2 ) * this%dt
        end do
    end do
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Regularizations
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Regularization of derivative
    ! NOTE: This gamma is not the same as the C++ version anymore....find out what the right value should be!
    do n = 1, this%N
        do k = 2, this%Nt
            this%cost = this%cost + this%gamma * (this%Ctp(n,k) - this%Ctp(n,k-1))**2 / this%dt
            this%grad(n,k) = this%grad(n,k) + 2.d0 * this%gamma * (this%Ctp(n,k) - this%Ctp(n,k-1)) / this%dt
        end do
        do k = 2, this%Nt - 1
            this%grad(n,k) = this%grad(n,k) - 2.d0 * this%gamma * (this%Ctp(n,k+1) - this%Ctp(n,k)) / this%dt
        end do
    end do
    
    ! Regularization of difference to reference Ctp
    do n = 1, this%N
        do k = 2, this%Nt
            ! NOTE: This was in here, and I didn't notice!
            ! NOTE: probably should be removed pow((double)n / (p.Nt() - 1), 2.0) * pow(Ctp[i][n] - p.Ctp_ref(), 2.0);
            this%cost = this%cost + &! ((k * 1.d0) / (this%Nt - 1.d0))**2 * &
             this%eta * (this%Ctp(n,k) - this%Ctp0)**2 * this%dt
            this%grad(n,k) = this%grad(n,k) + &! ((k * 1.d0) / (this%Nt - 1.d0))**2 * &
             2.d0 * this%eta * (this%Ctp(n,k) - this%Ctp0) * this%dt
        end do
    end do
    
    ! Set first step of gradient to zeros since it is fixed
    this%grad(:,1) = 0.d0
    
end subroutine run_noinput

subroutine makeDimensionless(this)
    implicit none
    class(MinimizedFarm), intent(inout)  :: this
    
    if (.not.this%isDimensionless) then
        this%isDimensionless = .true.
        this%Pref = this%Pref / this%POWER
        this%t    = this%t    / this%w%TIME
        this%dt   = this%dt   / this%w%TIME
        call this%w%makeDimensionless
        call this%wstar%makeDimensionless
        call this%iw%makeDimensionless
        call this%iwstar%makeDimensionless
    end if
end subroutine makeDimensionless

subroutine makeDimensional(this)
    implicit none
    class(MinimizedFarm), intent(inout)  :: this
    
    if (this%isDimensionless) then
        this%isDimensionless = .false.
        this%Pref = this%Pref * this%POWER
        this%t    = this%t    * this%w%TIME
        this%dt   = this%dt   * this%w%TIME
        call this%w%makeDimensional
        call this%wstar%makeDimensional
        call this%iw%makeDimensional
        call this%iwstar%makeDimensional
    end if
end subroutine makeDimensional

subroutine eval(this, x, f, g)
    implicit none
    class(MinimizedFarm), intent(inout)                   :: this
    real(rprec), dimension(:), intent(in)                 :: x
    real(rprec), intent(inout)                            :: f
    real(rprec), dimension(:),  intent(inout)             :: g
    integer                         :: k

    ! Place x in this%Ctp 
    do k = 1, this%Nt-1
        this%Ctp(:,k + 1) = x((k - 1) * this%N + 1 : this%N * k)
    end do
    this%Ctp(:,1) = this%Ctp0

    ! Run model
    call this%run_noinput
    
    ! Return cost function
    f = this%cost
    
    ! Return gradient as vector
    g = 0.d0
    do k = 1, this%Nt-1
        g((k - 1) * this%N + 1 : this%N * k) = this%grad(:,k + 1)
    end do
end subroutine eval

subroutine finiteDifferenceGradient(this)
    implicit none
    class(MinimizedFarm), intent(inout)     :: this
    type(MinimizedFarm)                     :: mf
    integer :: n, k
    real(rprec) :: dCtp
    
    this%fdgrad = 0.d0
    dCtp = sqrt( epsilon( this%Ctp(n,k) ) )
    do n = 1, this%N
        do k = 2, this%Nt
            mf = this
            mf%Ctp(n,k) = mf%Ctp(n,k) + dCtp
            call mf%run_noinput
            this%fdgrad(n,k) = (mf%cost - this%cost) / dCtp
        end do
    end do
end subroutine finiteDifferenceGradient

function get_Ctp_vector(this) result(Ctp_vec)
    implicit none
    class(MinimizedFarm) :: this
    real(rprec), dimension(:),allocatable   :: Ctp_vec
    integer :: k
    
    allocate(Ctp_vec(size(this%Ctp) - this%N))
    
    do k = 1, this%Nt-1
        Ctp_vec((k - 1) * this%N + 1 : this%N * k) = this%Ctp(:,k + 1)
    end do
    
end function get_Ctp_vector
    

end module rh_control