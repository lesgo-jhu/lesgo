!*******************************************************************************
subroutine stability(L, zo_s, phi_m, phi_h, psi_m, psi_h)
!*******************************************************************************
use types, only : rprec
use param, only : pi, dz, zo
! Stability functions from Hydrology (2005) Brutsaert
real(rprec), intent(in) :: L, zo_s
real(rprec), intent(out) :: phi_m, phi_h, psi_m, psi_h
! Constants for stable region
real(rprec), parameter :: am=6.1_rprec, bm=2.5_rprec, ah=5.3_rprec, bh=1.1_rprec
! Constants for unstable region
real(rprec), parameter :: a=0.33_rprec, b=0.41_rprec, c=0.33_rprec,            &
    d=0.057_rprec, n=0.78_rprec

! Neutral
if (abs(L) > 1/epsilon(0._rprec)) then
    phi_m = 1._rprec
    phi_h = 1._rprec
    psi_m = 0._rprec
    psi_h = 0._rprec
! Stable or unstable
else
    psi_m = -calc_psi_m(dz*0.5_rprec/L) + calc_psi_m(zo/L)
    ! write(*,*) dz*0.5_rprec/L, zo/L, psi_m
    ! stop
    psi_h = -calc_psi_h(dz*0.5_rprec/L) + calc_psi_h(zo_s/L)
    call calc_phi(dz*0.5_rprec/L, phi_m, phi_h)
end if

contains

!*******************************************************************************
subroutine calc_phi(zeta, o_phi_m, o_phi_h)
!*******************************************************************************
real(rprec), intent(in) :: zeta
real(rprec), intent(out) :: o_phi_m, o_phi_h
real(rprec) :: y

if ( zeta < 0._rprec) then
    y = -zeta
    if (y > b**(-3) ) then
        o_phi_m = 1._rprec
    else
        o_phi_m = ( a + (b*(y**(4._rprec/3._rprec))) )/(a+y)
    end if

    o_phi_h = (c + (d*(y**n)))/(c+(y**n))
else
    o_phi_m = 1._rprec + am*(zeta + zeta**bm                                     &
        * ((1._rprec + zeta**bm)**(-1._rprec + 1._rprec/bm)))                  &
        / (zeta +((1._rprec + zeta**bm)**(1._rprec/bm)))
    o_phi_h = 1._rprec + ah*(zeta + zeta**bh                                     &
        * ((1._rprec + zeta**bh)**(-1._rprec + 1._rprec/bh)))                  &
        / (zeta +((1._rprec + zeta**bh)**(1._rprec/bh)))
end if

end subroutine calc_phi

!*******************************************************************************
function calc_psi_m(zeta)
!*******************************************************************************
real(rprec), intent(in) :: zeta
real(rprec) :: calc_psi_m
real(rprec) :: y, yy, xx, psi_zero

if (zeta < 0._rprec) then
    y = -zeta

    if (y > b**(-3) ) then
        yy =  (b**(-3._rprec))
        xx = (yy/a)**(1._rprec/3._rprec)
        psi_zero = -log(a) + ((3._rprec**(1._rprec/2._rprec))                  &
            * b*(a**(1._rprec/3._rprec))*(pi/6._rprec))
        calc_psi_m = log(a+yy) - (3._rprec*b*(yy**(1._rprec/3._rprec)))             &
            + (((b*(a**(1._rprec/3._rprec)))/2._rprec)*log(((1+xx)**2)         &
            / (1 - xx + (xx**2)))) + ((3._rprec**(1._rprec/2._rprec))          &
            * b*(a**(1._rprec/3._rprec)))*(atan(((2._rprec*xx)-1)              &
            / (3._rprec**(1._rprec/2._rprec)))) + psi_zero
    else
        xx = (y/a)**(1._rprec/3._rprec)
        psi_zero = -log(a) + ((3._rprec**(1._rprec/2._rprec))                  &
            * b*(a**(1._rprec/3._rprec))*(pi/6._rprec))
        calc_psi_m = log(a+y) - (3._rprec*b*(y**(1._rprec/3._rprec)))               &
            + (((b*(a**(1._rprec/3._rprec)))/2._rprec)                         &
            * log(((1+xx)**2)/(1 - xx + (xx**2)))) + sqrt(3._rprec)            &
            * b*(a**(1._rprec/3._rprec))*(atan(((2._rprec*xx)-1)               &
            / (3._rprec**(1._rprec/2._rprec)))) + psi_zero
    end if
else
    calc_psi_m = -am*log(zeta + ((1._rprec + zeta**bm)**(1._rprec/bm)))
end if

end function  calc_psi_m

!*******************************************************************************
function calc_psi_h(zeta)
!*******************************************************************************
real(rprec), intent(in) :: zeta
real(rprec) :: calc_psi_h
real(rprec) :: y

if (zeta < 0._rprec) then
    y = -zeta
    calc_psi_h = ((1-d)/n)*log((c+(y**n))/c)
else
    calc_psi_h = -ah*log(zeta + ((1._rprec + zeta**bh)**(1._rprec/bh)))
end if

end function calc_psi_h

end subroutine stability
