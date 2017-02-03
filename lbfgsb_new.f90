module lbfgsb_class
use types, only : rprec
use minimize
use messages

private
public lbfgsb

type :: lbfgsb
    class(minimize_t), pointer :: mini => NULL()
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
    class(minimize_t), target :: i_mini
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