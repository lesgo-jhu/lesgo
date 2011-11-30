!**********************************************************************
module sgs_param
!**********************************************************************
use types, only : rprec

save
private rprec
public

! For all sgs models
    real(rprec) :: delta, nu
    real(rprec), dimension (:,:,:), allocatable :: S11, S12, S22, S33, S13, S23
    real(rprec), dimension(:,:,:), allocatable :: Nu_t      ! eddy viscosity
    integer ::jt_count
    real(rprec), dimension(:,:,:), allocatable ::Cs_opt2   ! (C_s)^2, Dynamic Smag coeff
    integer :: count_clip, count_all

! For Lagrangian models (4,5)
    real(rprec), parameter :: opftime = 1.5_rprec   ! (Meneveau, Lund, Cabot; JFM 1996)
    real(rprec), dimension(:,:,:), allocatable :: F_LM, F_MM, F_QN, F_NN, Beta
    real(rprec) :: lagran_dt = 0._rprec

! The following are for dynamically updating T, the timescale for Lagrangian averaging
!   F_ee2 is the running average of (eij*eij)^2
!   F_deedt2 is the running average of [d(eij*eij)/dt]^2
!   ee_past is the array (eij*eij) for the past timestep
    $if ($DYN_TN)
    real(rprec), dimension(:,:,:), allocatable :: F_ee2, F_deedt2
    real(rprec), dimension(:,:,:), allocatable :: ee_past
    $endif

contains

!**********************************************************************
subroutine sgs_param_init ()
!**********************************************************************
use param, only : ld,ny,nz,lbz,molec,nu_molec,u_star,z_i,dx,dy,dz
use test_filtermodule,only:filter_size

implicit none

! Allocate arrays

    ! For all sgs models:
    allocate ( S11(ld,ny,nz), &
         S12(ld,ny,nz), &
         S13(ld,ny,nz), &
         S22(ld,ny,nz), &
         S23(ld,ny,nz), &
         S33(ld,ny,nz) )

    allocate ( Nu_t(ld,ny,nz), Cs_opt2(ld,ny,nz) )

        S11 = 0.0_rprec
        S12 = 0.0_rprec
        S13 = 0.0_rprec
        S22 = 0.0_rprec
        S23 = 0.0_rprec
        S33 = 0.0_rprec

        Nu_t = 0.0_rprec
        Cs_opt2 = 0.0_rprec

    ! For Lagrangian models:
    allocate ( F_LM(ld,ny,lbz:nz), F_MM(ld,ny,lbz:nz), &
               F_QN(ld,ny,lbz:nz), F_NN(ld,ny,lbz:nz), &
               Beta(ld,ny,lbz:nz) )

        F_LM = 0.0_rprec
        F_MM = 0.0_rprec
        F_QN = 0.0_rprec
        F_NN = 0.0_rprec
        Beta = 0.0_rprec

    ! Lagrangian zero-crossing time scale variables
    $if ($DYN_TN)
    allocate ( F_ee2(ld,ny,lbz:nz), F_deedt2(ld,ny,lbz:nz), &
               ee_past(ld,ny,lbz:nz) }

        F_ee2 = 0.0_rprec
        F_deedt2 = 0.0_rprec
        ee_past = 0.0_rprec
    $endif

! Set constants
    delta=filter_size*(dx*dy*dz)**(1._rprec/3._rprec) ! nondimensional

    if (molec) then
        nu = (nu_molec/(u_star*z_i))    ! dimensionless
    else
        nu = 0._rprec
    end if   

return
end subroutine sgs_param_init

end module sgs_param
