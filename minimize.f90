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

module minimized_class
use types, only : rprec
use messages
implicit none

private
public Minimized

type :: Minimized
    ! pointer to function to minimize. Not used if extending this base class
    procedure(MinimizedFunction), pointer, nopass :: fun => NULL()  
contains
    ! This procedure is used by the minimizer to find the minimum
    procedure, public :: eval
end type Minimized

interface Minimized
    module procedure :: constructor
end interface Minimized

!! Interface for the scalar-valued function to be minimized. This interface must evaluate 
!! the function at x and return the function value f and gradient g at this location.
abstract interface
    subroutine MinimizedFunction(x, f, g)
        use types, only : rprec
        implicit none
        real(rprec), dimension(:), intent(in)                 :: x  ! Point to evaluate
        real(rprec), intent(inout)                            :: f  ! Function value (scalar)
        real(rprec), dimension(:), intent(inout)              :: g  ! Function gradient
    end subroutine MinimizedFunction
end interface

contains

!! Constructor for Minimized class. This constructor takes a procedure pointer matching
!! the MinimizedFunction interface. This constructor can be overloaded for a class that
!! extends this base class.
function constructor(i_fun) result(this)
    implicit none
    procedure(MinimizedFunction)  :: i_fun
    type(Minimized)               :: this
    
    this%fun => i_fun
end function  constructor

!! Evaluates the function pointer. Overload this procedure if extending this base class
subroutine eval(this, x, f, g)
    implicit none
    class(Minimized), intent(inout)             :: this
    real(rprec), dimension(:), intent(in)       :: x
    real(rprec), intent(inout)                  :: f
    real(rprec), dimension(:), intent(inout)    :: g
    
    call this%fun(x, f, g)
end subroutine eval

end module minimized_class

! http://src.gnu-darwin.org/ports/science/mpb/work/mpb-1.4.2/src/matrices/minpack2-linmin.c
! published under that link using GPL 2
module line_search_class
use types, only : rprec
use messages
use minimized_class
implicit none

private
public LineSearch

type :: LineSearch
    class(Minimized), pointer :: mini => NULL()
    ! NOTE: All should be const! --- I don't know how to do this in fortran.
    real(rprec)     :: ftol = 0.1               ! Sufficient decrease condition parameter
    real(rprec)     :: gtol = 0.4               ! Curvature condition parameter
    real(rprec)     :: xtol = 1E-6              ! Maximum relative width of uncertainty
    real(rprec)     :: minStep = 0.0            ! Minimum step size
    real(rprec)     :: maxStep = 10000          ! Maximum step size
    integer         :: maxFunEval = 100         ! Maximum number of function evaluations
    ! NOTE: Should not be const 
    logical         :: bracketed = .false.      ! The minimum is current bracketed
    real(rprec)     :: curMinStep  = 0.0        ! Minimum step for current trial within interval of uncertainty
    real(rprec)     :: curMaxStep   = 0.0       ! Maximum step for current trial within interval of uncertainty
    character (11)  :: class_name = 'LineSearch'
contains
    procedure, public  :: search
    procedure, private :: take_step
end type LineSearch

interface LineSearch
    module procedure :: constructor
end interface LineSearch

contains

function constructor(i_mini, i_ftol, i_gtol, i_xtol, i_minStep, i_maxStep, i_maxFunEval) &
 result(this)
    implicit none
    type(LineSearch) :: this
    class(Minimized), target :: i_mini
    real(rprec), intent(in), optional   :: i_ftol, i_gtol, i_xtol, i_minStep, i_maxStep
    integer, intent(in), optional       :: i_maxFunEval

    ! Assign input arguments
    ! NOTE: It would be nice to check arguments here, but I can't seem to call errors 
    !       from this function
    this%mini => i_mini
    if ( present(i_ftol) ) then
        this%ftol = i_ftol
    end if
    if ( present(i_gtol) ) then
        this%gtol = i_gtol
    end if
    if ( present(i_xtol) ) then
        this%xtol = i_xtol
    end if
    if ( present(i_minStep) ) then
        this%minStep = i_minStep
    end if
    if ( present(i_maxStep) ) then
        this%maxStep = i_maxStep 
    end if
    if ( present(i_maxFunEval) ) then
        this%maxFunEval = i_maxFunEval 
    end if
end function constructor

subroutine search(this, x, f, g, s, stp, o_numFunEval)
    implicit none
    class(LineSearch), intent(inout)                        :: this
    real(rprec), dimension(:), intent(inout)                :: x
    real(rprec), dimension(:), intent(inout)                :: g
    real(rprec), intent(inout)                              :: f, stp
    real(rprec), dimension(:), intent(in)                   :: s
    integer, intent(out), optional                          :: o_numFunEval
    integer :: numFunEval
    character (100) :: sub_name
    real(rprec), dimension(:), allocatable  :: x0
    logical     :: take_step_success = .true.
    logical     :: stage1 = .true.
    real(rprec) :: f0, dd0, ddtest, width, width1 
    real(rprec) :: stx, fx, ddx   
    real(rprec) :: sty, fy, ddy
    real(rprec) :: dg, ftest1
    real(rprec) :: fm, fxm, fym, dgm, ddxm, ddym
    
    ! Initialization
    numFunEval = 0
    this%bracketed = .false.
    sub_name = this%class_name // '.search'
    
    ! Check input arguments
    if (size(x) == 0) then
        call error(sub_name, 'Initial guess must be non-empty.')
    end if
    if (size(x) /= size(g)) then
        call error(sub_name, 'gradient and input must be the same size')
    end if
    if (size(x) /= size(s)) then
        call error(sub_name, 'x and s must be the same size.')
    end if
    if (stp <= 0.0) then
        call error(sub_name, 'Initial step size guess must be > 0.')
    end if
    if (this%minStep < 0.0) then
        call error(sub_name, 'Maximum step size must be greater than minimum step size.')
    end if
    if (this%maxStep < this%minStep) then
        call error(sub_name, 'Minimum step size must be non-negative.')
    end if
    if (this%ftol < 0.0 .or. this%gtol < 0.0 .or. this%xtol < 0.0) then
        call error(sub_name, 'Search tolerances must be >= 0.')
    end if
    
    ! Initialization point x0 and associated function evaluation f0 and 
    ! directional derivative dd0
    allocate(x0(size(x)))
    x0 = x
    f0 = f
    dd0 = sum(g*s)
    if (dd0 >= 0.0) call error(sub_name, 'Gradient is non-decreasing')
    
    ! Initialize local variables.
    ddtest = this%ftol * dd0
    width  = this%maxStep -  this%minStep
    width1 = 2*width

    ! step, function, and directional derivative at the best step.
    stx = 0.0
    fx = f0
    ddx = dd0
    ! step, function, and derivative at the other endpoint of the interval of uncertainty.
    sty = 0.0
    fy = f0
    ddy = dd0
    ! Start of iteration
    do
        ! Interval of uncertainty
        if (this%bracketed) then
            this%curMinStep = min(stx, sty)
            this%curMaxStep = max(stx, sty)
        else
            this%curMinStep = stx
            this%curMaxStep = stp + 4.0 * (stp-stx)       
        end if
    
        ! Force the step to be within the bounds maxStep and minStep.
        stp = max(stp, this%minStep);
        stp = min(stp, this%maxStep);
        ! If an unusual termination is to occur then let stp be the lowest point 
        ! obtained so far.
        if ( (this%bracketed .and. (stp <= this%curMinStep .or. stp >= this%curMaxStep)) &
         .or. numFunEval >= this%maxFunEval - 1 .or. .not.take_step_success .or.         &
         (this%bracketed .and. this%curMaxStep - this%curMinStep <=                      &
         this%xtol * this%curMaxStep) )  then
            stp = stx
        end if

        ! Evaluate the function and gradient at stp and compute the directional derivative.
        x = x0 + stp * s
        call this%mini%eval(x, f, g)
        numFunEval = numFunEval + 1
        dg = sum(g*s)
        ftest1 = f0 + stp * ddtest

        ! Tests for convergence issues. Strange behaviors will return normally and 
        ! print a warning to std output
        if ((this%bracketed .and. (stp <= this%curMinStep .or. stp >= this%curMaxStep))  &
        .or. .not.take_step_success) then
            call mesg(sub_name, 'Rounding errors prevent further progress. There may'   &
             // ' not be a step which satisfies the sufficient decrease and curvature'   &
             // ' conditions. Tolerances may be too small')
            o_numFunEval = numFunEval
            return
        end if
        if (stp == this%maxStep .and. f <= ftest1 .and. dg < ddtest) then
            call mesg(sub_name,'The step is at the upper bound.')
            o_numFunEval = numFunEval
            return
        end if
        if (stp == this%minStep .and. (f > ftest1 .or. dg >= ddtest)) then
            call mesg(sub_name,'The step is at the lower bound.')
            o_numFunEval = numFunEval
            return
        end if
        if (numFunEval >= this%maxFunEval) then
            call mesg(sub_name,'Maximum number of function evaluations has been reached.')
            o_numFunEval = numFunEval
            return
        end if
        if (this%bracketed .and.                                                         &
         this%curMaxStep - this%curMinStep <= this%xtol * this%curMaxStep) then
            call mesg(sub_name,'Relative width of the interval of uncertainty has '     &
             // 'reached tolerance.')
            o_numFunEval = numFunEval
            return
        end if
        
        ! Test for convergence. i.e. the sufficient decrease and curvature conditions hold.
        if ( f <= ftest1 .and. abs(dg) <= this%gtol * (-dd0) ) then
           o_numFunEval = numFunEval
           return
        end if

        ! In the first stage we seek a step for which the modified function has a 
        ! nonpositive value and nonnegative
        if ( stage1 .and. f <= ftest1 .and. dg >= min(this%ftol, this%gtol) * dd0 ) then
            stage1 = .false.
        end if
        
        ! A modified function is used to predict the step only if we have not obtained a 
        ! step for which the modified function has a nonpositive function value and 
        ! nonnegative derivative, and if a lower function value has been obtained but the 
        ! decrease is not sufficient.
        if (stage1 .and. f <= fx .and. f > ftest1) then
            ! Define the modified function and derivative values.
            fm = f - stp * ddtest
            fxm = fx - stx * ddtest
            fym = fy - sty * ddtest
            dgm = dg - ddtest
            ddxm = ddx - ddtest
            ddym = ddy - ddtest

            ! Call cstep to update the interval of uncertainty and to compute the new step.
            call this%take_step(stx, fxm, ddxm, sty, fym, ddym, stp, fm, dgm, take_step_success)

            ! Reset the function and gradient values for f.
            fx = fxm + stx * ddtest
            fy = fym + sty * ddtest
            ddx = ddxm + ddtest
            ddy = ddym + ddtest
        else
            ! Call cstep to update the interval of uncertainty and to compute the new step.
            call this%take_step(stx, fx, ddx, sty, fy, ddy, stp, f, dg, take_step_success)
        end if
        
        ! Force a sufficient decrease in the size of the interval of uncertainty.
        if (this%bracketed) then
            if (abs(sty-stx) >= 2.0/3.0 * width1) then
                stp = stx + 0.5 * (sty - stx)
            end if
            width1 = width
            width = abs(sty - stx)
        end if
    
        continue
    end do
end subroutine search

subroutine take_step(this, stx, fx, dx, sty, fy, dy, stp, fp, dp, take_step_success)
    implicit none
    class(LineSearch), intent(inout)    :: this
    real(rprec), intent(inout)          :: stx, fx, dx, sty, fy, dy, stp, fp, dp
    logical, intent(inout)              :: take_step_success
    integer         :: info = 0
    real(rprec)     :: sgnd
    character (100) :: sub_name
    logical         :: bound
    real(rprec)     :: theta, s, gamma
    real(rprec)     :: p, q, r, stpc, stpq, stpf
    
    sub_name = this%class_name // '.search'
    
    ! Check the input parameters for errors.
    if ( ( this%bracketed .and. ( stp <= min(stx, sty) .or. stp >= max(stx, sty) ) )     &
     .or. dx * (stp - stx) >= 0.0 .or. this%curMaxStep < this%curMinStep ) then
        call error(sub_name, 'Invalid input argument')
    end if
 
    ! Determine if the derivatives have opposite sign.
    sgnd = dp * ( dx / abs(dx) );
    
    ! First case. A higher function value. The minimum is bracketed. If the cubic step is 
    ! closer to stx than the quadratic step, the cubic step is taken, else the average of 
    ! the cubic and quadratic steps is taken.
    if (fp > fx) then
        info = 1
        bound = .true.
        theta = 3.0 * (fx - fp) / (stp - stx) + dx + dp
        s = max(abs(theta), max(abs(dx),abs(dp)))
        gamma = s * sqrt((theta / s)**2 - (dx / s) * (dp / s))
        if (stp < stx)      gamma = -gamma
        p = (gamma - dx) + theta
        q = ( (gamma - dx) + gamma ) + dp
        r = p / q
        stpc = stx + r * (stp - stx);
        stpq = stx + ((dx / ((fx - fp) / (stp - stx)+dx)) / 2.0) * (stp - stx)
        if (abs(stpc - stx) < abs(stpq - stx)) then
            stpf = stpc
        else    
            stpf = stpc + (stpq - stpc) / 2.0;
        end if
        this%bracketed = .true.
    ! Second case. A lower function value and derivatives of opposite sign. The minimum is
    ! bracketed. If the cubic step is closer to stx than the quadratic (secant) step, the 
    ! cubic step is taken, else the quadratic step is taken.
    else if (sgnd < 0.0) then
        info = 2
        bound = .false.
        theta = 3.0 * (fx - fp) / (stp - stx) + dx + dp
        s = max(abs(theta), max(abs(dx), abs(dp)))
        gamma = s * sqrt((theta / s)**2 - (dx / s) * (dp / s));
        if (stp > stx) gamma = -gamma
        p = (gamma - dp) + theta
        q = ((gamma - dp) + gamma) + dx
        r = p/q
        stpc = stp + r * (stx - stp)
        stpq = stp + (dp / (dp - dx)) * (stx - stp)
        if (abs(stpc - stp) > abs(stpq - stp)) then
            stpf = stpc
        else
            stpf = stpq
        end if
        this%bracketed = .true.
    ! Third case. A lower function value, derivatives of the same sign, and the magnitude 
    ! of the derivative decreases. The cubic step is only used if the cubic tends to 
    ! inf0y in the direction of the step or if the minimum of the cubic is beyond stp. 
    ! Otherwise the cubic step is defined to be either curMinStep or curMaxStep. The 
    ! quadratic (secant) step is also computed and if the minimum is bracketed then the 
    ! the step closest to stx is taken, else the step farthest away is taken.
    else if(abs(dp) <abs(dx)) then
        info = 3
        bound = .true.
        theta = 3.0 * (fx - fp) / (stp - stx) + dx + dp
        s = max(abs(theta), max(abs(dx), abs(dp)))
        ! The case gamma = 0 only arises if the cubic does not tend to inf0y in the 
        ! direction of the step.
        gamma = s*sqrt(max(0.0,(theta / s)**2 - (dx / s) * (dp / s)))
        if (stp > stx)  gamma = -gamma
        p = (gamma - dp) + theta
        q = (gamma + (dx - dp)) + gamma
        r = p/q
        if (r < 0.0 .and. gamma /= 0.0) then
            stpc = stp + r*(stx - stp)
        else if (stp > stx) then
            stpc = this%curMaxStep
        else
            stpc = this%curMinStep
        end if

        stpq = stp + (dp/(dp-dx))*(stx - stp)
        if (this%bracketed) then
            if (abs(stp-stpc) < abs(stp-stpq)) then
                stpf = stpc
            else
                stpf = stpq
            end if
        else
            if (abs(stp-stpc) > abs(stp-stpq)) then
                stpf = stpc
            else
                stpf = stpq
            end if
        end if
    ! Fourth case. A lower function value, derivatives of the same sign, and the magnitude
    ! of the derivative does not decrease. If the minimum is not bracketed, the step is 
    ! either curMinStep or curMaxStep, else the cubic step is taken.
    else
        info = 4
        bound = .false.
        if (this%bracketed) then
            theta = 3*(fp - fy)/(sty - stp) + dy + dp
            s = max(abs(theta), max(abs(dy),abs(dp)))
            gamma = s*sqrt((theta/s)**2 - (dy/s)*(dp/s))
            if (stp > sty) gamma = -gamma;
            p = (gamma - dp) + theta
            q = ((gamma - dp) + gamma) + dy
            r = p/q
            stpc = stp + r*(sty - stp)
            stpf = stpc
        else if (stp > stx) then
            stpf = this%curMaxStep
        else
            stpf = this%curMinStep
        end if
    end if

    ! Update the interval of uncertainty. This update does not depend on the new step or the case analysis above.
    if (fp > fx) then
        sty = stp
        fy = fp
        dy = dp
    else
        if (sgnd < 0.0) then
            sty = stx
            fy = fx
            dy = dx
        end if
        stx = stp
        fx = fp
        dx = dp
    end if
  
    ! Compute the new step and safeguard it.
    stpf = min(this%curMaxStep, stpf)
    stpf = max(this%curMinStep, stpf)
    stp  = stpf
    if (this%bracketed .and. bound) then
        if (sty > stx) then
            stp = min(stx + (2.0/3.0) * (sty-stx), stp)
        else
            stp = max(stx + (2.0/3.0) * (sty-stx), stp)
        end if
    end if
    
    take_step_success = info > 0
    
end subroutine take_step

end module line_search_class

module lbfgsb_class
use types, only : rprec
use minimized_class
use messages

private
public lbfgsb

type :: lbfgsb
    class(Minimized), pointer :: mini => NULL()
    integer :: maxIter = 10000      ! maximum number of iterations
    real(rprec) :: lb = -10000000   ! -Infinity
    real(rprec) :: ub = 10000000    ! Infinity
    real(rprec) :: tol = 1E-6       ! convergence level
contains
   procedure, public :: minimize
   procedure, private :: minimize_priv
end type lbfgsb

interface lbfgsb
    module procedure :: constructor
end interface lbfgsb

contains

function constructor(i_mini, i_maxIter, i_lb, i_ub, i_tol) result(this)
    implicit none
    type(lbfgsb) :: this
    class(Minimized), target :: i_mini
    integer, intent(in), optional :: i_maxIter 
    real(rprec), intent(in), optional :: i_lb, i_ub, i_tol

    ! Assign input arguments
    if ( present(i_maxIter) )   this%maxIter = i_maxIter
    if ( present(i_tol) )       this%tol     = i_tol
    if ( present(i_lb) )        this%lb      = i_lb
    if ( present(i_ub) )        this%ub      = i_ub
    this%mini => i_mini
end function constructor

subroutine minimize(this, i_x, o_x)
    implicit none
    class(lbfgsb), intent(inout) :: this
    real(rprec), dimension(:), intent(in) :: i_x
    real(rprec), dimension(:), intent(out), optional :: o_x
    real(rprec), dimension(:), allocatable :: x_work
    integer :: n, m
    
    ! Set sizes of arrays
    n = size(i_x)
    m = 20
    
    ! Allocate work array
    allocate(x_work(size(i_x)))
    x_work = i_x
 
    ! Call private method
    call this%minimize_priv(x_work, n, m)
    
    ! Set output if present
    if ( present(o_x) ) o_x = x_work
    
end subroutine minimize

subroutine minimize_priv(this, x, n, m)
    implicit none
    class(lbfgsb), intent(inout) :: this
    integer, intent(in) :: n, m
    real(rprec), dimension(n), intent(inout) :: x
    real(rprec), dimension(n) :: l, u, g
    real(rprec), dimension((2*m + 5)*n + 11*m*m + 8*m) :: wa
    integer, dimension(3*n) :: iwa
    integer, dimension(n) :: nbd
    real(rprec) :: f, factr, pgtol
    character*60 :: task, csave
    logical, dimension(4) :: lsave
    integer, dimension(44) :: isave
    real(rprec), dimension(29) :: dsave
    integer :: iprint
    integer :: iter
    
    l = this%lb
    u = this%ub
    nbd = 2
    factr = this%tol/epsilon(1._rprec)
    pgtol = 0._rprec
    task = 'START'
    iprint = -1
    
    iter = 0
    do while (iter < this%maxIter)
        call setulb(n, m, x, l, u, nbd, f, g, factr, pgtol, wa, iwa, task,     &
            iprint, csave, lsave, isave, dsave)
        
        ! On a return with task(1:2)='FG', the user must evaluate the
        !   function f and gradient g at the returned value of x.
        if (task(1:2) == 'FG') then
            call this%mini%eval(x, f, g)
        end if
        
        ! On a return with task(1:5)='NEW_X', an iteration of the
        !   algorithm has concluded, and f and g contain f(x) and g(x)
        !   respectively.  The user can decide whether to continue or stop
        !   the iteration. 
        if (task(1:5) == 'NEW_X') then
            iter = iter + 1
            cycle
        end if
        
        ! When (1:4)='CONV', the termination test in L-BFGS-B has been 
        !   satisfied;
        if (task(1:4) == 'CONV') then
            exit
        end if
        
        ! When task(1:4)='ABNO', the routine has terminated abnormally
        !   without being able to satisfy the termination conditions,
        !   x contains the best approximation found,
        !   f and g contain f(x) and g(x) respectively;
        if (task(1:4) == 'ABNO') then
            call mesg('lbfgs:minimize','routine terminated abnormally')
            exit
        end if
        
        ! When task(1:5)='ERROR', the routine has detected an error in the
        !   input parameters;
        if (task(1:4) == 'ABNO') then
            call mesg('lbfgs:minimize','routine detected an error')
            exit
        end if
        
    end do

    ! Print result
    write(*,*) 'L-BFGS-B terminated after ', iter, 'iterations. Minimum f = ',f
    
end subroutine minimize_priv

end module lbfgsb_class

module conjugate_gradient_class
use types, only : rprec
use line_search_class
use minimized_class
use messages
implicit none

private
public ConjugateGradient

type :: ConjugateGradient
    class(Minimized), pointer :: mini => NULL()
    real(rprec)         :: eps = 1E-10      ! A small number
    integer             :: maxIter = 10000  ! maximum number of CG iterations
    real(rprec)         :: tol = 1E-6       ! convergence level
    real(rprec)         :: f = 0.0          ! current function evaluation
    real(rprec)         :: fp = 0.0         ! previous function evaluation
    integer             :: fnev = 0         ! number of function evaluations
    real(rprec)         :: gamma = 0.0      ! conjugate direction step
    type(LineSearch)    :: ls               ! line search class
    real(rprec), dimension(:), allocatable :: gd    ! conjugate direction
    real(rprec), dimension(:), allocatable :: x, xp ! current and previous location
    real(rprec), dimension(:), allocatable :: g, gp ! current and previous gradients
    real(rprec), dimension(:), allocatable :: h, hp ! current and previous search direction
    real(rprec) :: lb = -10000000!-Infinity
    real(rprec) :: ub = 10000000!Infinity
contains
   procedure, public    :: minimize
   procedure, private   :: evaluate_gamma
end type ConjugateGradient

interface ConjugateGradient
    module procedure :: constructor
end interface ConjugateGradient

contains

function constructor(i_mini, i_maxIter, i_lb, i_ub, i_tol) result(this)
    implicit none
    type(ConjugateGradient) :: this
    class(Minimized), target :: i_mini
    integer, intent(in), optional       :: i_maxIter 
    real(rprec), intent(in), optional   :: i_lb, i_ub, i_tol

    ! Assign input arguments
    if ( present(i_maxIter) )   this%maxIter = i_maxIter
    if ( present(i_tol) )       this%tol     = i_tol
    if ( present(i_lb) )        this%lb      = i_lb
    if ( present(i_ub) )        this%ub      = i_ub
    this%mini => i_mini
    this%ls = LineSearch(i_mini)
end function constructor

subroutine evaluate_gamma(this)
    implicit none
    class(ConjugateGradient), intent(inout) :: this
    
    this%gd = this%g - this%gp
    this%gamma = sum(this%g * this%gd) / sum(this%gp * this%gp)
end subroutine evaluate_gamma

subroutine minimize(this, i_x, o_x)
    implicit none
    class(ConjugateGradient), intent(inout) :: this
    real(rprec), dimension(:), intent(in)   :: i_x
    real(rprec), dimension(:), intent(out), optional :: o_x
    real(rprec) :: d, delta_f, stp
    integer     :: i, j
    integer     :: dummy = 0
    
    ! Allocate arrays
    allocate(this%xp(size(i_x)))
    allocate(this%gp(size(i_x)))
    allocate(this%hp(size(i_x)))
    allocate(this%x(size(i_x)))
    allocate(this%g(size(i_x)))
    allocate(this%h(size(i_x)))
    allocate(this%gd(size(i_x)))
    
    ! Evaluate gradient and function at starting position
    this%xp = i_x
    call this%mini%eval(this%xp, this%fp, this%gp)
    this%fnev = this%fnev + 1
    this%hp = -this%gp
    
    ! Assign other arrays
    this%f  = this%fp
    this%g  = this%gp
    this%h  = this%hp
    this%x  = this%xp
    this%gd = this%g
    
    d = sum(this%g * this%h)
    delta_f = -0.5 * d
    
    do i = 1, this%maxIter
        ! Compute derivative with respect to stp at origin
        d = sum(this%g * this%h)

        ! If derivative is positive, reset at the gradient
        if (d > 0) then
            this%h = -this%g
            d = sum(this%g * this%h)
            delta_f = -0.5 * d
        end if

        ! Search in the direction hp
        stp = min(1.0, -2.0 * delta_f / d)
        call this%ls%search(this%x, this%f, this%g, this%h, stp, dummy)
        this%fnev = this%fnev + dummy
        
        ! Safeguard step against bounds
        do j = 1, size(this%x)
            this%x(j) = min( max(this%x(j), this%lb), this%ub)
        end do
        call this%mini%eval(this%x, this%f, this%g)

        ! Check for convergence
        if ( 2.0 * abs(this%fp - this%f) <= this%tol * ( abs(this%fp) + abs(this%f) +    &
         this%eps ) .or. sum(this%g*this%g) == 0 ) then
            ! Set output if present
            if ( present(o_x) ) o_x = this%x
            
            ! Evaluate minimization at current point
            call this%mini%eval(this%x, this%f, this%g)
            
            ! Print result
            write(*,*) 'Conjugate gradient terminated after ', i,              &
                'iterations. Minimum f = ', this%f
            
            return
        else if (this%f > this%fp) then
            ! Reset to previous point
            this%x = this%xp;
            this%f = this%fp;
            this%g = this%gp;
            
            ! Set output if present
            if ( present(o_x) ) o_x = this%x
            
            ! Evaluate minimization at current point
            call this%mini%eval(this%x, this%f, this%g)
            
            ! Print result
            write(*,*) 'Conjugate gradient terminated after ', i,              &
                'iterations. Minimum f = ', this%f
            return
        end if

        ! Compute next search direction
        call this%evaluate_gamma
        this%h = -this%g + this%gamma * this%hp

        ! Swap arrays
        delta_f = max(this%fp - this%f, 10 * epsilon(this%f))
        
        this%gp = this%g
        this%hp = this%h
        this%xp = this%x
        this%fp = this%f
    end do
    
    ! Evaluate minimization at current point
    call this%mini%eval(this%x, this%f, this%g)
    
    ! Set output if present
    if ( present(o_x) ) o_x = this%x
    
    ! Print result
            write(*,*) 'Conjugate gradient terminated after ', this%maxIter,   &
                'iterations. Minimum f = ', this%f
end subroutine minimize

end module conjugate_gradient_class
