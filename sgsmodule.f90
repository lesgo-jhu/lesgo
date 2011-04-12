module sgsmodule
use types,only:rprec
use param,only:ld,nx,ny,nz
implicit none
private ld,ny,nz

$if ($MPI)
    $define $lbz 0
$else
    $define $lbz 1
$endif

! The following are for dynamic Lagranrian SGS models (model=4,5) 
real(kind=rprec),parameter::opftime=1.5_rprec   ! (Meneveau, Lund, Cabot; JFM 1996)
real(kind=rprec),dimension(ld,ny,nz)::F_LM,F_MM,F_QN,F_NN,Beta

! The following are for dynamically updating T, the timescale for Lagrangian averaging
!   F_ee2 is the running average of (eij*eij)^2
!   F_deedt2 is the running average of [d(eij*eij)/dt]^2
!   ee_past is the array (eij*eij) for the past timestep
$if ($DYN_TN)
real(kind=rprec),dimension(ld,ny,$lbz:nz) :: F_ee2, F_deedt2
real(kind=rprec),dimension(ld,ny,$lbz:nz) :: ee_past
$endif

!real(kind=rprec),dimension(ld,ny,nz)::Betaclip  !--not used

real(rprec), dimension(ld, ny, nz) :: Nu_t      ! eddy viscosity
real(kind=rprec),dimension(ld,ny,$lbz:nz)::u_lag,v_lag,w_lag
integer ::jt_count
real(kind=rprec),dimension(ld,ny,nz)::Cs_opt2   ! (C_s)^2, Dynamic Smag coeff
integer :: count_clip, count_all

!**********************************************************************                 
contains
!**********************************************************************                 

real(kind=rprec) function rtnewt(A, jz)
use types,only:rprec
integer,parameter :: jmax=100
real(kind=rprec) :: x1,x2,xacc
integer :: j, jz
real(kind=rprec) :: df,dx,f
real(kind=rprec), dimension(0:5) :: A
x1 = 0._rprec
x2 = 15._rprec  ! try to find the largest root first....hmm
xacc = 0.001_rprec ! doesn't need to be that accurate...
rtnewt = 0.5_rprec*(x1+x2)
do j=1,jmax
   f = A(0)+rtnewt*(A(1)+rtnewt*(A(2)+rtnewt*(A(3)+rtnewt*(A(4)+rtnewt*A(5)))))
   df = A(1) + rtnewt*(2._rprec*A(2) + rtnewt*(3._rprec*A(3) +&
        rtnewt*(4._rprec*A(4) + rtnewt*(5._rprec*A(5)))))
   dx=f/df
   rtnewt = rtnewt - dx
!        if ((x1-rtnewt)*(rtnewt-x2) < 0.) STOP 'rtnewt out of bounds'
   if (abs(dx) < xacc) return
end do
rtnewt = 1._rprec  ! if dont converge fast enough
write(6,*) 'using beta=1 at jz= ', jz
end function rtnewt

end module sgsmodule
