!**********************************************************************
module sgsmodule
!**********************************************************************
use types,only:rprec
use param2,only:ld,ny,nz
implicit none
private ld,ny,nz
!TS In genreal, ofttime is 1.5 (Meneveau et al., 1996)
real(kind=rprec),parameter::opftime=1.5_rprec
real(kind=rprec),allocatable, dimension(:,:,:) :: F_LM,F_MM,F_QN,F_NN,Beta
!real(kind=rprec),dimension(ld,ny,nz)::Betaclip  !--not used
!xxxx----- Added by Vij - 04/14/04--xxxx----------------
! Nu_t is needed for scalar sgs
! For more details look into scalars_module.f90
$if ($MPI)
  real (rprec), allocatable, dimension(:,:,:) :: Nu_t
$else
  real(kind=rprec), allocatable, dimension(:,:,:)::Nu_t
$endif
!xxxx---- Vij change ends here --xxxx-------------
real(kind=rprec),allocatable, dimension(:,:,:)::u_lag,v_lag,w_lag
integer ::jt_count
real(kind=rprec),allocatable, dimension(:,:,:)::Cs_opt2!,Cs_opt2_avg
                 !--Cs_opt2_avg commented to save mem.
contains

!**********************************************************************
subroutine alloc_sgsmodule()
!**********************************************************************
implicit none

allocate(F_LM(ld,ny,nz))
allocate(F_MM(ld,ny,nz))
allocate(F_QN(ld,ny,nz))
allocate(F_NN(ld,ny,nz))
allocate(Beta(ld,ny,nz))
allocate(Nu_t(ld,ny,nz))

$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif
allocate(u_lag(ld,ny,$lbz:nz))
allocate(v_lag(ld,ny,$lbz:nz))
allocate(w_lag(ld,ny,$lbz:nz))

allocate(Cs_opt2(ld,ny,nz))

return
end subroutine alloc_sgsmodule

!**********************************************************************
real(kind=rprec) function rtnewt(A, jz)
!**********************************************************************
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
