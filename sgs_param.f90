!**********************************************************************
module sgs_param
!**********************************************************************
use types, only : rprec

save
private rprec
public

! For all sgs models
    real(rprec) :: delta, nu
    real(rprec), dimension(:,:,:), allocatable :: S11, S12, S22, S33, S13, S23
    real(rprec), dimension(:,:,:), allocatable :: Nu_t      ! eddy viscosity
    integer ::jt_count
    real(rprec), dimension(:,:,:), allocatable ::Cs_opt2   ! (C_s)^2, Dynamic Smag coeff
    integer :: count_clip, count_all

    real(rprec), dimension(:,:), allocatable :: L11,L12,L13,L22,L23,L33
    real(rprec), dimension(:,:), allocatable :: M11,M12,M13,M22,M23,M33

    real(rprec), dimension(:,:), allocatable :: S_bar,S11_bar,S12_bar,S13_bar, &
                                     S22_bar,S23_bar,S33_bar,S_S11_bar,S_S12_bar,S_S13_bar, &
                                     S_S22_bar, S_S23_bar, S_S33_bar
    real(rprec), dimension(:,:), allocatable :: S, u_bar,v_bar,w_bar

! For Lagrangian models (4,5)
    real(rprec), parameter :: opftime = 1.5_rprec   ! (Meneveau, Lund, Cabot; JFM 1996)
    real(rprec), dimension(:,:,:), allocatable :: F_LM, F_MM, F_QN, F_NN, Beta
    real(rprec) :: lagran_dt = 0._rprec

! For scale dependent models (3,5)
    real(rprec), dimension(:,:), allocatable :: Q11,Q12,Q13,Q22,Q23,Q33     ! these are only for model 5
    real(rprec), dimension(:,:), allocatable :: N11,N12,N13,N22,N23,N33     ! these are only for model 5

    real(rprec), dimension(:,:), allocatable :: S_hat,S11_hat,S12_hat,S13_hat, &
                                     S22_hat,S23_hat,S33_hat,S_S11_hat,S_S12_hat,S_S13_hat, &
                                     S_S22_hat, S_S23_hat, S_S33_hat
    real(rprec), dimension(:,:), allocatable :: u_hat,v_hat,w_hat

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
use param, only : ld,ny,nz,lbz,molec,nu_molec,u_star,z_i,dx,dy,dz,sgs_model
use test_filtermodule,only:filter_size

implicit none

! Allocate arrays

    ! For all sgs models:
    allocate ( S11(ld,ny,nz), S12(ld,ny,nz), S13(ld,ny,nz), &
        S22(ld,ny,nz), S23(ld,ny,nz), S33(ld,ny,nz) )

    allocate ( Nu_t(ld,ny,nz), Cs_opt2(ld,ny,nz) )

    allocate ( L11(ld,ny), L12(ld,ny), L13(ld,ny), &
        L22(ld,ny), L23(ld,ny), L33(ld,ny) )

    allocate ( M11(ld,ny), M12(ld,ny), M13(ld,ny), &
        M22(ld,ny), M23(ld,ny), M33(ld,ny) )

    allocate ( S_bar(ld,ny), S11_bar(ld,ny), S12_bar(ld,ny), &
        S13_bar(ld,ny), S22_bar(ld,ny), S23_bar(ld,ny), &
        S33_bar(ld,ny), S_S11_bar(ld,ny), S_S12_bar(ld,ny), &
        S_S13_bar(ld,ny), S_S22_bar(ld,ny), S_S23_bar(ld,ny), &
        S_S33_bar(ld,ny) )

    allocate ( S(ld,ny), u_bar(ld,ny), v_bar(ld,ny), w_bar(ld,ny)  )

        S11 = 0.0_rprec; S12 = 0.0_rprec; S13 = 0.0_rprec
        S22 = 0.0_rprec; S23 = 0.0_rprec; S33 = 0.0_rprec

        Nu_t = 0.0_rprec; Cs_opt2 = 0.0_rprec

        L11 = 0.0_rprec; L12 = 0.0_rprec; L13 = 0.0_rprec
        L22 = 0.0_rprec; L23 = 0.0_rprec; L33 = 0.0_rprec

        M11 = 0.0_rprec; M12 = 0.0_rprec; M13 = 0.0_rprec
        M22 = 0.0_rprec; M23 = 0.0_rprec; M33 = 0.0_rprec

        S_bar(ld,ny) = 0.0_rprec; S11_bar(ld,ny) = 0.0_rprec
        S12_bar(ld,ny) = 0.0_rprec; S13_bar(ld,ny) = 0.0_rprec
        S22_bar(ld,ny) = 0.0_rprec; S23_bar(ld,ny) = 0.0_rprec
        S33_bar(ld,ny) = 0.0_rprec; S_S11_bar(ld,ny) = 0.0_rprec
        S_S12_bar(ld,ny) = 0.0_rprec; S_S13_bar(ld,ny) = 0.0_rprec
        S_S22_bar(ld,ny) = 0.0_rprec; S_S23_bar(ld,ny) = 0.0_rprec
        S_S33_bar(ld,ny) = 0.0_rprec;

        S = 0.0_rprec
        u_bar = 0.0_rprec; v_bar = 0.0_rprec; w_bar = 0.0_rprec

    ! For Lagrangian models:
    if ((sgs_model .eq. 4).or.(sgs_model .eq. 5)) then
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

    endif

    ! For scale dependent models:
    if ((sgs_model .eq. 3).or.(sgs_model .eq. 5)) then
        allocate ( S_hat(ld,ny), S11_hat(ld,ny), S12_hat(ld,ny), &
            S13_hat(ld,ny), S22_hat(ld,ny), S23_hat(ld,ny), &
            S33_hat(ld,ny), S_S11_hat(ld,ny), S_S12_hat(ld,ny), &
            S_S13_hat(ld,ny), S_S22_hat(ld,ny), S_S23_hat(ld,ny), &
            S_S33_hat(ld,ny) )

        allocate ( u_hat(ld,ny), v_hat(ld,ny), w_hat(ld,ny)  )

            S_hat(ld,ny) = 0.0_rprec; S11_hat(ld,ny) = 0.0_rprec
            S12_hat(ld,ny) = 0.0_rprec; S13_hat(ld,ny) = 0.0_rprec
            S22_hat(ld,ny) = 0.0_rprec; S23_hat(ld,ny) = 0.0_rprec
            S33_hat(ld,ny) = 0.0_rprec; S_S11_hat(ld,ny) = 0.0_rprec
            S_S12_hat(ld,ny) = 0.0_rprec; S_S13_hat(ld,ny) = 0.0_rprec
            S_S22_hat(ld,ny) = 0.0_rprec; S_S23_hat(ld,ny) = 0.0_rprec
            S_S33_hat(ld,ny) = 0.0_rprec;

            u_hat = 0.0_rprec; v_hat = 0.0_rprec; w_hat = 0.0_rprec

        if (sgs_model .eq. 5) then
            allocate ( Q11(ld,ny), Q12(ld,ny), Q13(ld,ny), &
                Q22(ld,ny), Q23(ld,ny), Q33(ld,ny) )

            allocate ( N11(ld,ny), N12(ld,ny), N13(ld,ny), &
                N22(ld,ny), N23(ld,ny), N33(ld,ny) )

                Q11 = 0.0_rprec; Q12 = 0.0_rprec; Q13 = 0.0_rprec
                Q22 = 0.0_rprec; Q23 = 0.0_rprec; Q33 = 0.0_rprec

                N11 = 0.0_rprec; N12 = 0.0_rprec; N13 = 0.0_rprec
                N22 = 0.0_rprec; N23 = 0.0_rprec; N33 = 0.0_rprec
        endif

    endif

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
