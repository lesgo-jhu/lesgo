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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                   Wake Model Class
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module wake_model_class
use types, only : rprec
use util,  only : logistic, softplus, gaussian
use messages
implicit none

private
public WakeModel

type :: WakeModel
    real(rprec), dimension(:),   allocatable :: s    ! turbine location
    real(rprec), dimension(:),   allocatable :: k    ! wake expansion coefficient
    real(rprec), dimension(:),   allocatable :: x    ! streamwise coordinate
    real(rprec), dimension(:,:), allocatable :: G    ! Gaussian forcing function (turbine, space)
    real(rprec), dimension(:,:), allocatable :: d    ! dimensionless wake diameter (turbine, space)
    real(rprec), dimension(:,:), allocatable :: dp   ! d/dx of d (turbine, space)
    real(rprec), dimension(:,:), allocatable :: w    ! wake expansion function (turbine, space)
    real(rprec), dimension(:,:), allocatable :: fp   ! forcing prefactor f = fp*Ctp/(4+Ctp)  (turbine, space)
    real(rprec), dimension(:,:), allocatable :: du   ! velocity deficit (turbine, space)
    real(rprec), dimension(:), allocatable   :: u    ! superimposed velocity (space)
    real(rprec), dimension(:), allocatable   :: uhat ! estimated local turbine velocity (turbine)
    real(rprec), dimension(:), allocatable   :: Phat ! estimated turbine power (turbine)
    real(rprec), dimension(:), allocatable   :: Ctp  ! local thrust coefficient (turbine)
    real(rprec) :: U_infty = 0                       ! inlet velocity
    real(rprec) :: Delta   = 0                       ! Gaussian forcing width
    real(rprec) :: Dia     = 0                       ! rotor diameter
    real(rprec) :: dx      = 0                       ! rotor diameter
    integer     :: Nx      = 0                       ! Number of streamwise points
    integer     :: N       = 0                       ! Number of turbines
    logical     :: isDimensionless = .false.
!     logical     :: useRungeKutta = .false.
    real(rprec) :: LENGTH=0, VELOCITY=0, TIME=0, FORCE=0
contains
    procedure, private :: initialize_val
    procedure, private :: initialize_file
    procedure, public  :: print
    procedure, public  :: write_to_file
    procedure, public  :: makeDimensionless
    procedure, public  :: makeDimensional
    procedure, public  :: advance
    procedure, public  :: computeWakeExpansionFunctions
    procedure, private :: rhs
    
end type WakeModel

interface WakeModel
    module procedure :: constructor_val
    module procedure :: constructor_file
end interface WakeModel

contains

! Constructor for wake model with values given
function constructor_val(i_s, i_U_infty, i_Delta, i_k, i_Dia, i_Nx) result(this)
    implicit none
    type(WakeModel)                       :: this
    real(rprec), intent(in)               :: i_U_infty, i_Delta, i_Dia
    real(rprec), dimension(:), intent(in) :: i_s, i_k
    integer, intent(in)                   :: i_Nx
    
    call this%initialize_val(i_s, i_U_infty, i_Delta, i_k, i_Dia, i_Nx)
end function constructor_val

! Constructor for wake model that reads from file
function constructor_file(fstring) result(this)
    use open_file_fid_mod
    use param, only : CHAR_BUFF_LENGTH
    implicit none
    
    type(WakeModel)          :: this
    character(*), intent(in) :: fstring
    
    call this%initialize_file(fstring)

end function constructor_file

subroutine initialize_val(this, i_s, i_U_infty, i_Delta, i_k, i_Dia, i_Nx)
    implicit none
    class(WakeModel), intent(inout)       :: this
    real(rprec), intent(in)               :: i_U_infty, i_Delta, i_Dia
    real(rprec), dimension(:), intent(in) :: i_s, i_k
    integer, intent(in)                   :: i_Nx
    integer                               :: i
    
    ! Allocate based on number of turbines
    this%N = size(i_s)
    if ( size(i_k) /= this%N ) then
        call error('WakeModel.initialize','s and k must be the same size')
    end if
    allocate( this%s(this%N)    )
    allocate( this%k(this%N)    )
    allocate( this%uhat(this%N) )
    allocate( this%Phat(this%N) )
    allocate( this%Ctp(this%N)  )

    ! Allocate based on number of gridpoints
    this%Nx = i_Nx
    allocate( this%x(this%Nx) )
    allocate( this%u(this%Nx) )
    allocate( this%G(this%N,  this%Nx) )
    allocate( this%d(this%N,  this%Nx) )
    allocate( this%dp(this%N, this%Nx) )
    allocate( this%w(this%N,  this%Nx) )
    allocate( this%fp(this%N, this%Nx) )
    allocate( this%du(this%N, this%Nx) )

    ! Assign input arguments
    this%s          = i_s
    this%U_infty    = i_U_infty
    this%Delta      = i_Delta
    this%k          = i_k
    this%Dia        = i_Dia
    
    ! Normalization constants
    this%VELOCITY = i_U_infty
    this%LENGTH   = i_Dia
    this%TIME     = this%LENGTH / this%VELOCITY
    this%FORCE    = this%VELOCITY / this%TIME
    
    ! Calculate other variables
    this%dx = ( this%s(1) + this%s(this%N) ) / this%Nx
    this%x(1) = this%dx/2
    do i = 2, this%Nx
        this%x(i) = this%x(i-1) + this%dx
    end do
    do i = 1, this%N
        this%G(i,:) = gaussian(this%x, this%s(i), this%Delta)
    end do

    call this%computeWakeExpansionFunctions

    this%du(:,:) = 0.d0
    this%u(:)    = this%U_infty
    this%uhat(:) = this%U_infty
    this%Phat(:) = 0.d0
    this%Ctp(:)  = 0.d0
    
end subroutine initialize_val

subroutine initialize_file(this, fstring)
    use open_file_fid_mod
    use param, only : CHAR_BUFF_LENGTH
    implicit none

    class(WakeModel), intent(inout) :: this
    character(*), intent(in)        :: fstring
    integer :: i, fid
    
    !  Open vel.out (lun_default in io) for final output
    fid = open_file_fid(fstring, 'rewind', 'unformatted')
    read(fid) this%N, this%Nx, this%dx, this%Dia, this%Delta,                           &
              this%U_infty, this%isDimensionless
    read(fid) this%LENGTH, this%VELOCITY, this%TIME, this%FORCE
    
    allocate( this%s(this%N)    )
    allocate( this%k(this%N)    )
    allocate( this%uhat(this%N) )
    allocate( this%Phat(this%N) )
    allocate( this%Ctp(this%N)  )
    allocate( this%x(this%Nx) )
    allocate( this%u(this%Nx) )
    allocate( this%G(this%N,  this%Nx) )
    allocate( this%d(this%N,  this%Nx) )
    allocate( this%dp(this%N, this%Nx) )
    allocate( this%w(this%N,  this%Nx) )
    allocate( this%fp(this%N, this%Nx) )
    allocate( this%du(this%N, this%Nx) )    
    
    read(fid) this%s
    read(fid) this%k
    read(fid) this%x
    read(fid) this%du
    read(fid) this%u
    read(fid) this%uhat
    read(fid) this%Phat
    read(fid) this%Ctp
    close(fid)
    
    do i = 1, this%N
        this%G(i,:) = gaussian(this%x, this%s(i), this%Delta)
    end do
    call this%computeWakeExpansionFunctions
    
end subroutine initialize_file

subroutine computeWakeExpansionFunctions(this)
    implicit none
    class(WakeModel), intent(inout) :: this
    integer                         :: i
    
    do i = 1, this%N
        this%d(i,:)  = 1.0 + this%k(i) * softplus(2.0 *                                  &
            (this%s(i) + 2.0 * this%Delta)/this%Dia, 2.0*this%x/this%Dia)
        this%dp(i,:) = 2.0 * this%k(i) * logistic(2.0 *                                  &
            (this%s(i) + 2.0 * this%Delta)/this%Dia, 2.0*this%x/this%Dia)/this%Dia
        this%w(i,:)  = 2.0 * this%U_infty * this%dp(i,:) / this%d(i,:)
        this%fp(i,:) = 2.0 * this%U_infty**2 * this%G(i,:) / ( this%d(i,:)**2 )   
    end do

end subroutine 

subroutine makeDimensionless(this)
    implicit none
    class(WakeModel), intent(inout)  :: this
    
    if (.not.this%isDimensionless) then
        this%isDimensionless = .true.
        this%s       = this%s / this%LENGTH
        this%U_infty = this%U_infty / this%VELOCITY
        this%Delta   = this%Delta / this%LENGTH
        this%Dia     = this%Dia / this%LENGTH
        this%dx      = this%dx / this%LENGTH
        this%x       = this%x / this%LENGTH
        this%G       = this%G * this%LENGTH         ! G has units 1/length
        this%dp      = this%dp * this%LENGTH        ! dp has units 1/length
        this%w       = this%w * this%TIME           ! w has units 1/time     
        this%du      = this%du / this%VELOCITY
        this%u       = this%u / this%VELOCITY
        this%uhat    = this%uhat / this%VELOCITY
        this%Phat    = this%uhat / this%VELOCITY**3
        this%fp      = this%fp / this%FORCE
    end if
end subroutine makeDimensionless

subroutine makeDimensional(this)
    implicit none
    class(WakeModel), intent(inout)  :: this
    
    if (this%isDimensionless) then
        this%isDimensionless = .false.
        this%s       = this%s * this%LENGTH
        this%U_infty = this%U_infty * this%VELOCITY
        this%Delta   = this%Delta * this%LENGTH
        this%Dia     = this%Dia * this%LENGTH
        this%dx      = this%dx * this%LENGTH
        this%x       = this%x * this%LENGTH
        this%G       = this%G / this%LENGTH         ! G has units 1/length
        this%dp      = this%dp / this%LENGTH        ! dp has units 1/length
        this%w       = this%w / this%TIME           ! w has units 1/time     
        this%du      = this%du * this%VELOCITY
        this%u       = this%u * this%VELOCITY
        this%uhat    = this%uhat * this%VELOCITY
        this%Phat    = this%uhat * this%VELOCITY**3
        this%fp      = this%fp * this%FORCE
    end if
end subroutine makeDimensional

! Writes object to file
subroutine write_to_file(this, fstring)
    use open_file_fid_mod
    use param, only : CHAR_BUFF_LENGTH
    implicit none
    class(WakeModel), intent(in) :: this
    character(CHAR_BUFF_LENGTH), intent(in)  :: fstring
    integer :: i, fid
    
    !  Open vel.out (lun_default in io) for final output
    fid = open_file_fid(fstring, 'rewind', 'unformatted')
    write(fid) this%N, this%Nx, this%dx, this%Dia, this%Delta,                           &
               this%U_infty, this%isDimensionless
    write(fid) this%LENGTH, this%VELOCITY, this%TIME, this%FORCE
    write(fid) this%s
    write(fid) this%k
    write(fid) this%x
    write(fid) this%du
    write(fid) this%u
    write(fid) this%uhat
    write(fid) this%Phat
    write(fid) this%Ctp
    close(fid)

end subroutine write_to_file

! Prints all variables of the class to standard output
subroutine print(this)
    implicit none
    class(WakeModel), intent(in) :: this
    integer :: i

    write(*,*) ' U_infty         = ', this%U_infty
    write(*,*) ' Delta           = ', this%Delta
    write(*,*) ' Dia             = ', this%Dia
    write(*,*) ' Nx              = ', this%Nx
    write(*,*) ' x               = ', this%x
    write(*,*) ' dx              = ', this%dx
    write(*,*) ' isDimensionless = ', this%isDimensionless
    write(*,*) ' LENGTH          = ', this%LENGTH
    write(*,*) ' VELOCITY        = ', this%VELOCITY
    write(*,*) ' TIME            = ', this%TIME
    write(*,*) ' FORCE           = ', this%FORCE
    write(*,*) ' u               = ', this%u
    do i = 1, this%N
        write(*,*) ' Wake', i,':'
        write(*,*) '  uhat    = ', this%uhat(i)
        write(*,*) '  s       = ', this%s(i)
        write(*,*) '  k       = ', this%k(i)
        write(*,*) '  G       = ', this%G(i,:)
        write(*,*) '  d       = ', this%d(i,:)
        write(*,*) '  dp      = ', this%dp(i,:)
        write(*,*) '  w       = ', this%w(i,:)
    end do
end subroutine print

subroutine advance(this, Ctp, dt)
    implicit none
    class(WakeModel), intent(inout)        :: this
    real(rprec), intent(in)                :: dt
    real(rprec), dimension(:), intent(in)  :: Ctp
    integer                                :: i
    real(rprec), dimension(:), allocatable :: du_superimposed
    
    if (size(Ctp) /= this%N) then
        call error('WakeModel.advance','Ctp must be size N')
    end if

    ! Compute new wake deficit and superimpose wakes
    allocate(du_superimposed(this%Nx))
    du_superimposed = 0.0
    do i = 1, this%N
        this%du(i,:) = this%du(i,:) +                                                    &
          dt * this%rhs(this%du(i,:), this%fp(i,:) * this%Ctp(i) / (4.0 + this%Ctp(i)), i)
        du_superimposed = du_superimposed + this%du(i,:)**2
    end do
    du_superimposed = sqrt(du_superimposed)
    
    ! Find the velocity field
    this%u = this%U_infty - du_superimposed
    
    ! Find estimated velocities
    do i = 1, this%N
        this%Ctp(i) = Ctp(i)
        this%uhat(i) = sum(this%G(i,:) * this%u * this%dx)
        this%Phat(i) = this%Ctp(i) * this%uhat(i)**3
    end do
    
end subroutine advance

! Evaluates RHS of wake equation
function rhs(this, du, f, i) result(ddudt)
    use util, only : ddx_upwind1
    implicit none
    class(WakeModel), intent(in)           :: this
    real(rprec), dimension(:), intent(in)  :: f, du
    integer, intent(in)                    :: i
    real(rprec), dimension(:), allocatable :: ddudt, ddudx

    allocate(ddudt(this%Nx))
    allocate(ddudx(this%Nx))
    ddudx = ddx_upwind1(du, this%dx)

    ddudt = -this%U_infty * ddudx - this%w(i,:) * du + f
end function rhs

end module wake_model_class

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                   Wake Model Estimator Class
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module wake_model_estimator_class
use types, only : rprec
use messages
use wake_model_class
implicit none

private
public WakeModelEstimator

type :: WakeModelEstimator
    type(WakeModel), dimension(:), allocatable :: ensemble
    type(WakeModel)                            :: wm
    real(rprec)                                :: sigma_du, sigma_k, sigma_Phat
    integer                                    :: Ne, Nm, Ns
    real(rprec), dimension(:,:), allocatable   :: A, Aprime, Ahat, Ahatprime, E, D, Dprime
    real(rprec), dimension(:), allocatable     :: Abar, Ahatbar
    real(rprec) :: tau_U_infty = 300
contains
    procedure, private :: initialize_val
    procedure, private :: initialize_file
    procedure, public  :: write_to_file
    procedure, public  :: generateInitialEnsemble
    procedure, public  :: advance
    procedure, private :: advanceEnsemble
end type WakeModelEstimator

interface WakeModelEstimator
    module procedure :: constructor_val
    module procedure :: constructor_file
end interface WakeModelEstimator

contains 

! Constructor for wake model
function constructor_val(i_s, i_U_infty, i_Delta, i_k, i_Dia, i_Nx, i_Ne,                &
                         i_sigma_du, i_sigma_k, i_sigma_Phat) result(this)
    implicit none
    type(WakeModelEstimator)              :: this
    real(rprec), intent(in)               :: i_U_infty, i_Delta, i_Dia
    real(rprec), dimension(:), intent(in) :: i_s, i_k
    integer, intent(in)                   :: i_Nx
    integer, intent(in)                   :: i_Ne
    real(rprec), intent(in)               :: i_sigma_du, i_sigma_k, i_sigma_Phat
        
    call this%initialize_val(i_s, i_U_infty, i_Delta, i_k, i_Dia, i_Nx, i_Ne,            &
                             i_sigma_du, i_sigma_k, i_sigma_Phat)
end function constructor_val

! Constructor for wake model that reads from file
function constructor_file(fpath) result(this)
    use open_file_fid_mod
    use param, only : CHAR_BUFF_LENGTH
    implicit none
    
    type(WakeModelEstimator) :: this
    character(*), intent(in) :: fpath
    
    call this%initialize_file(fpath)

end function constructor_file

subroutine initialize_val(this, i_s, i_U_infty, i_Delta, i_k, i_Dia, i_Nx, i_Ne,         &
                          i_sigma_du, i_sigma_k, i_sigma_Phat)
    implicit none
    class(WakeModelEstimator), intent(inout)    :: this
    real(rprec), intent(in)                     :: i_sigma_du, i_U_infty, i_Delta, i_Dia
    real(rprec), dimension(:), intent(in)       :: i_s, i_k
    integer, intent(in)                         :: i_Nx
    integer, intent(in)                         :: i_Ne
    real(rprec), intent(in)                     :: i_sigma_k, i_sigma_Phat
    integer                                     :: i
    
    ! Set std deviations for noise
    this%sigma_du   = i_sigma_du
    this%sigma_k    = i_sigma_k
    this%sigma_Phat = i_sigma_Phat
    
    ! Create ensemble members
    this%Ne = i_Ne
    allocate( this%ensemble(this%Ne) )
    do i = 1, this%Ne
        this%ensemble(i) = WakeModel(i_s, i_U_infty, i_Delta, i_k, i_Dia, i_Nx)
    end do
    
    ! Create wake model estimate
    this%wm = WakeModel(i_s, i_U_infty, i_Delta, i_k, i_Dia, i_Nx)
    
    ! Allocate filter matrices
    this%Nm = this%wm%N                              ! Number of measurements
    this%Ns = this%wm%N * this%wm%Nx + this%wm%N - 1 ! Number of states
    allocate( this%Abar(this%Ns) )
    allocate( this%Ahatbar(this%Nm) )
    allocate( this%A(this%Ns, this%Ne) )
    allocate( this%Aprime(this%Ns, this%Ne) )
    allocate( this%Ahat(this%Nm, this%Ne) )
    allocate( this%Ahatprime(this%Nm, this%Ne) )
    allocate( this%E(this%Nm, this%Ne) )
    allocate( this%D(this%Nm, this%Ne) )
    allocate( this%Dprime(this%Nm, this%Ne) )
    
end subroutine initialize_val

subroutine initialize_file(this, fpath)
    use open_file_fid_mod
    use param, only : CHAR_BUFF_LENGTH
    use string_util, only : string_splice
    implicit none

    class(WakeModelEstimator)   :: this
    character(*), intent(in)    :: fpath
    character(CHAR_BUFF_LENGTH) :: fstring, str
    integer                     :: i, fid    

    fstring = fpath // '/wm_est.dat'
    fid = open_file_fid(fstring, 'rewind', 'unformatted')
    read(fid) this%Ne, this%Nm, this%Ns
    read(fid) this%sigma_du, this%sigma_k, this%sigma_Phat
    
    allocate( this%Abar(this%Ns) )
    allocate( this%Ahatbar(this%Nm) )
    allocate( this%A(this%Ns, this%Ne) )
    allocate( this%Aprime(this%Ns, this%Ne) )
    allocate( this%Ahat(this%Nm, this%Ne) )
    allocate( this%Ahatprime(this%Nm, this%Ne) )
    allocate( this%E(this%Nm, this%Ne) )
    allocate( this%D(this%Nm, this%Ne) )
    allocate( this%Dprime(this%Nm, this%Ne) )
    
    read(fid) this%A
    read(fid) this%Aprime
    read(fid) this%Ahat
    read(fid) this%Ahatprime
    read(fid) this%E
    read(fid) this%D
    read(fid) this%Dprime
    read(fid) this%Abar
    read(fid) this%Ahatbar
    close(fid)
    
    fstring = fpath // '/wm.dat'
    this%wm = WakeModel(fstring)
    
    allocate( this%ensemble(this%Ne) )
    do i = 1, this%Ne
        call string_splice( fstring, fpath // '/ensemble_', i, '.dat' )
        this%ensemble(i) = WakeModel(fstring)
    end do
    
    
end subroutine initialize_file

subroutine write_to_file(this, fpath)
    use open_file_fid_mod
    use param, only : CHAR_BUFF_LENGTH
    use string_util, only : string_splice
    implicit none
    
    class(WakeModelEstimator)   :: this
    character(*), intent(in)    :: fpath
    character(CHAR_BUFF_LENGTH) :: fstring, str
    integer                     :: i, fid    

    call system('mkdir -vp ' // fpath)
    fstring = fpath // '/wm_est.dat'
    fid = open_file_fid(fstring, 'rewind', 'unformatted')
    write(fid) this%Ne, this%Nm, this%Ns
    write(fid) this%sigma_du, this%sigma_k, this%sigma_Phat
    write(fid) this%A
    write(fid) this%Aprime
    write(fid) this%Ahat
    write(fid) this%Ahatprime
    write(fid) this%E
    write(fid) this%D
    write(fid) this%Dprime
    write(fid) this%Abar
    write(fid) this%Ahatbar
    close(fid)
    
    fstring = fpath // '/wm.dat'
    call this%wm%write_to_file(fstring)
    
    do i = 1, this%Ne
        call string_splice( fstring, fpath // '/ensemble_', i, '.dat' )
        call this%ensemble(i)%write_to_file(fstring)
    end do

end subroutine write_to_file

subroutine generateInitialEnsemble(this, Ctp)
    use util, only : random_normal, init_random_seed
    implicit none
    class(WakeModelEstimator), intent(inout)    :: this
    real(rprec), dimension(:)                   :: Ctp
    real(rprec)                                 :: dt
    real(rprec), parameter                      :: cfl = 0.99
    real(rprec)                                 :: FFT
    integer                                     :: i, j, N, Nx, jstart, jend
    
    if (size(Ctp) /= this%ensemble(1)%N) then
        call error('WakeModelEstimator.generateInitialEnsemble','Ctp must be of size N')
    end if
    
    write(*,*) 'Generating initial wake model filter ensemble...'

    ! Initialize random number generator
    call init_random_seed
    
    ! Calculate safe dt.
    dt          = cfl * this%wm%dx / this%wm%U_infty
    FFT         = this%wm%x(this%wm%Nx) / this%wm%U_Infty

    ! Do at least 1 FFT of simulations to get good ensemble statistics  
    N = this%wm%N  
    do i = 1, floor(FFT / dt)
        call this%advanceEnsemble(Ctp, dt)
    end do
    
    ! Place ensemble into a matrix with each member in a column
    this%Abar    = 0
    this%Ahatbar = 0
    N  = this%wm%N
    Nx = this%wm%Nx
    do i = 1, this%Ne
        do j = 1, N
            jstart = (j-1)*Nx+1
            jend   = j*Nx
            this%A(jstart:jend,i) = this%ensemble(i)%du(j,:)
        end do
        this%A(Nx*N+1:,i) = this%ensemble(i)%k(1:N-1)
        this%Ahat(:,i)    = this%ensemble(i)%Phat
        this%Abar         = this%Abar + this%A(:,i) / this%Ne
        this%Ahatbar      = this%Ahatbar + this%Ahat(:,i) / this%Ne
    end do
    do j = 1, this%Ne
        this%Aprime(:,j) = this%A(:,j) - this%Abar
        this%Ahatprime(:,j) = this%Ahat(:,j) - this%Ahatbar
    end do
    
end subroutine generateInitialEnsemble

subroutine advance(this, dt, Pm, Ctp)
    use util, only : random_normal, inverse
    implicit none
    class(WakeModelEstimator), intent(inout)    :: this
    real(rprec), intent(in)                     :: dt
    real(rprec), dimension(:), intent(in)       :: Pm, Ctp
    integer                                     :: i, j, N, Nx, jstart, jend
    real(rprec)                                 :: alpha, Uinftyi, r1_if

    N = this%wm%N
    Nx = this%wm%Nx

    if (size(Pm) /= N ) then
        call error('WakeModelEstimator.advance','Pm must be of size N')
    end if
    if (size(Ctp) /= N ) then
        call error('WakeModelEstimator.advance','Ctp must be of size N')
    end if

    ! Calculate noisy measurements
    this%E = 0
    do i = 1, N
        do j = 1, this%Ne
            this%E(i, j) = this%sigma_Phat * random_normal()
            this%D(i, j) = Pm(i) + this%E(i,j)
        end do
    end do
    this%Dprime = this%D - this%Ahat
    
    ! Update Anew = A + A'*Ahat'^T * (Ahat'*Ahat'^T + E*E^T)^-1 * D'
    ! Since the dimension is small, we don't bother doing the SVD. If the matrix becomes
    ! singular, then this should be considered as in section 4.3.2 of Everson(2003)
    this%A = this%A + matmul( matmul(this%Aprime, transpose(this%Ahatprime)),            &
    matmul(inverse(matmul(this%Ahatprime, transpose(this%Ahatprime)) +                   &
    matmul(this%E, transpose(this%E))), this%Dprime))
        
    ! Compute mean
    this%Abar = 0
    do i = 1, this%Ne
        this%Abar = this%Abar + this%A(:,i) / this%Ne;
    end do

    ! Filter U_infty
    alpha = dt / (this%tau_U_infty + dt)
    r1_if = this%wm%Ctp(1) / ( 4.d0 + this%wm%Ctp(1) )
    Uinftyi = ( Pm(1) / this%wm%Ctp(1) )**(1.d0/3.d0) / (1.d0 - r1_if)
    this%wm%U_infty = alpha * Uinftyi + (1 - alpha) * this%wm%U_infty
        
    ! Fill into objects
    do i = 1, this%Ne
        do j = 1, N
            jstart = (j-1)*Nx+1
            jend   = j*Nx
            this%ensemble(i)%du(j,:)  = this%A(jstart:jend,i)
        end do
        this%ensemble(i)%U_infty  = this%wm%U_infty
        this%ensemble(i)%k(1:N-1) = this%A(Nx*N+1:,i)
        this%ensemble(i)%k(N)     = this%ensemble(i)%k(N-1)
    end do
    do j = 1, N
        jstart = (j-1)*Nx+1
        jend   = j*Nx
        this%wm%du(j,:)  = this%Abar(jstart:jend)
    end do
    this%wm%k(1:N-1) = this%Abar(Nx*N+1:)
    this%wm%k(N)     = this%wm%k(N-1)


    ! Advance ensemble and mean estimate
    call this%advanceEnsemble(Ctp, dt)
    
    ! Place ensemble into a matrix with each member in a column
    this%Abar    = 0
    this%Ahatbar = 0
    N  = this%wm%N
    Nx = this%wm%Nx
    do i = 1, this%Ne
        do j = 1, N
            jstart = (j-1)*Nx+1
            jend   = j*Nx
            this%A(jstart:jend,i) = this%ensemble(i)%du(j,:)
        end do
        this%A(Nx*N+1:,i) = this%ensemble(i)%k(1:N-1)
        this%Ahat(:,i)    = this%ensemble(i)%Phat
        this%Abar         = this%Abar + this%A(:,i) / this%Ne
        this%Ahatbar      = this%Ahatbar + this%Ahat(:,i) / this%Ne
    end do
    do j = 1, this%Ne
        this%Aprime(:,j) = this%A(:,j) - this%Abar
        this%Ahatprime(:,j) = this%Ahat(:,j) - this%Ahatbar
    end do
    
end subroutine 

subroutine advanceEnsemble(this, Ctp, dt)
    use param, only : pi
    use util, only  : random_normal
    implicit none
    class(WakeModelEstimator), intent(inout) :: this
    real(rprec), dimension(:), intent(in)    :: Ctp
    real(rprec), intent(in)                  :: dt
    integer                                  :: i, j, N
    
    ! Advance step for objects
    N = this%wm%N
    do i = 1, this%Ne
        do j = 1,N
            this%ensemble(i)%k(j) = this%ensemble(i)%k(j)                                &
                                  + sqrt(dt) * this%sigma_k * random_normal()
            this%ensemble(i)%du(j,:) = this%ensemble(i)%du(j,:)                          &
                + sqrt(dt) * this%sigma_du * random_normal()                             &
                * this%ensemble(i)%G(j,:) * sqrt(2*pi) * this%ensemble(i)%Delta
        end do
        this%ensemble(i)%k(N) = this%ensemble(i)%k(N-1)
        call this%ensemble(i)%computeWakeExpansionFunctions
        call this%ensemble(i)%advance(Ctp, dt)
    end do
    call this%wm%computeWakeExpansionFunctions
    call this%wm%advance(Ctp, dt)
    
end subroutine advanceEnsemble

end module wake_model_estimator_class

