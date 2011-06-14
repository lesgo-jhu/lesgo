$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif

! this is the w-node version
!--provides Cs_opt2 1:nz
!--MPI: required u,v on 0:nz, except bottom node 1:nz

subroutine lagrange_Sdep(S11,S12,S13,S22,S23,S33)
! standard dynamic model to calculate the Smagorinsky coefficient
! this is done layer-by-layer to save memory
! everything is done to be on uv-nodes
! -note: we need to calculate |S| here, too.
! stuff is done on uv-nodes
! can save more mem if necessary. mem requirement ~ n^2, not n^3
use types,only:rprec
use param
use sim_param,only:u,v,w
use sgsmodule,only:F_LM,F_MM,F_QN,F_NN,beta,Cs_opt2,opftime,lagran_dt
use test_filtermodule
use immersedbc,only:n_bldg,bldg_pts,building_interp
$if ($DYN_TN)
use sgsmodule, only:F_ee2,F_deedt2,ee_past
$endif
$if($LVLSET)
use level_set, only : level_set_Cs_lag_dyn
$endif
implicit none

$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif

integer :: jx,jy,jz
integer :: i, px, py, lx, ly, lz
integer :: counter1,counter2,counter3,counter4,counter5
integer :: istart, iend

real(kind=rprec):: tf1,tf2,tf1_2,tf2_2 ! Size if the second test filter
real(kind=rprec) :: fractus
real(kind=rprec), dimension(ld,ny,nz) :: S11,S12,S13,S22,S23,S33
real(kind=rprec) :: Betaclip  !--scalar to save mem., otherwise (ld,ny,nz)
real(kind=rprec), dimension(ld,ny) :: Cs_opt2_2d,Cs_opt2_4d

real(kind=rprec), dimension(ld,ny) :: S,tempos
real(kind=rprec), dimension(ld,ny) :: L11,L12,L13,L22,L23,L33
real(kind=rprec), dimension(ld,ny) :: Q11,Q12,Q13,Q22,Q23,Q33
real(kind=rprec), dimension(ld,ny) :: M11,M12,M13,M22,M23,M33
real(kind=rprec), dimension(ld,ny) :: N11,N12,N13,N22,N23,N33

real(kind=rprec), dimension(nz) :: LMvert,MMvert,QNvert,NNvert
real(kind=rprec), dimension(ld,ny) :: LM,MM,QN,NN,Tn,epsi,dumfac
real(kind=rprec), dimension(ld,ny) :: ee_now

real(kind=rprec), dimension(ld,ny) :: S_bar,S11_bar,S12_bar,&
     S13_bar,S22_bar,S23_bar,S33_bar,S_S11_bar, S_S12_bar,&
     S_S13_bar, S_S22_bar, S_S23_bar, S_S33_bar
real(kind=rprec), dimension(ld,ny) :: S_hat,S11_hat,S12_hat,&
     S13_hat,S22_hat,S23_hat,S33_hat,S_S11_hat, S_S12_hat,&
     S_S13_hat, S_S22_hat, S_S23_hat, S_S33_hat

real(kind=rprec), dimension(ld,ny) :: u_bar,v_bar,w_bar
real(kind=rprec), dimension(ld,ny) :: u_hat,v_hat,w_hat

real(kind=rprec) :: delta,const
real(kind=rprec) :: opftdelta,powcoeff

real(kind=rprec), parameter :: zero=1.e-24_rprec ! zero = infimum(0)

logical, save :: F_LM_MM_init = .false.
logical, save :: F_QN_NN_init = .false.

!---------------------------------------------------------------------
$if ($VERBOSE)
write (*, *) 'started lagrange_Sdep'
$endif

! Set coefficients
    delta = filter_size*(dx*dy*dz)**(1._rprec/3._rprec)
    opftdelta = opftime*delta
    powcoeff = -1._rprec/8._rprec
    fractus= 1._rprec/real(ny*nx,kind=rprec)
    const = 2._rprec*(delta**2)
    tf1=2._rprec
    tf2=4._rprec
    tf1_2=tf1**2
    tf2_2=tf2**2

! "Rearrange" F_* (running averages) so that their new positions (i,j,k) 
!   correspond to the current (i,j,k) particle
call interpolag_Sdep()

! For each horizontal level, calculate Lij(:,:), Qij(:,:), Mij(:,:), and Nij(:,:).  
!   Then update the running averages, F_*(:,:,jz), which are used to 
!   calculate Cs_opt2(:,:,jz).
do jz = 1,nz
    ! Calculate Lij
        ! Interp u,v,w onto w-nodes and store result as u_bar,v_bar,w_bar
        ! (except for very first level which should be on uvp-nodes)
        if ( ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) .and.  &
            (jz == 1) ) then  ! uvp-nodes
            u_bar(:,:) = u(:,:,1)
            v_bar(:,:) = v(:,:,1)
            w_bar(:,:) = .25_rprec*w(:,:,2)
        else  ! w-nodes
            u_bar(:,:) = .5_rprec*(u(:,:,jz) + u(:,:,jz-1)) 
            v_bar(:,:) = .5_rprec*(v(:,:,jz) + v(:,:,jz-1))  
            w_bar(:,:) = w(:,:,jz)
        end if
        u_hat = u_bar
        v_hat = v_bar
        w_hat = w_bar  

        ! First term before filtering (not the final value)
        L11=u_bar*u_bar
        L12=u_bar*v_bar
        L13=u_bar*w_bar
        L23=v_bar*w_bar
        L22=v_bar*v_bar
        L33=w_bar*w_bar              

        Q11 = u_bar*u_bar
        Q12 = u_bar*v_bar
        Q13 = u_bar*w_bar
        Q22 = v_bar*v_bar
        Q23 = v_bar*w_bar
        Q33 = w_bar*w_bar

        ! Filter first term and add the second term to get the final value
        call test_filter(u_bar,G_test)   ! in-place filtering
        call test_filter(v_bar,G_test)
        call test_filter(w_bar,G_test)
        call test_filter(L11,G_test)  
        L11 = L11 - u_bar*u_bar  
        call test_filter(L12,G_test)
        L12 = L12 - u_bar*v_bar
        call test_filter(L13,G_test)
        L13 = L13 - u_bar*w_bar
        call test_filter(L22,G_test)
        L22 = L22 - v_bar*v_bar
        call test_filter(L23,G_test)
        L23 = L23 - v_bar*w_bar
        call test_filter(L33,G_test)
        L33 = L33 - w_bar*w_bar       
        
        call test_filter(u_hat,G_test_test)
        call test_filter(v_hat,G_test_test)
        call test_filter(w_hat,G_test_test)
        call test_filter(Q11,G_test_test)
        Q11 = Q11 - u_hat*u_hat
        call test_filter(Q12,G_test_test)
        Q12 = Q12 - u_hat*v_hat
        call test_filter(Q13,G_test_test)
        Q13 = Q13 - u_hat*w_hat
        call test_filter(Q22,G_test_test)
        Q22 = Q22 - v_hat*v_hat
        call test_filter(Q23,G_test_test)
        Q23 = Q23 - v_hat*w_hat
        call test_filter(Q33,G_test_test)
        Q33 = Q33 - w_hat*w_hat

    ! Calculate |S|
        S(:,:) = sqrt(2._rprec*(S11(:,:,jz)**2+S22(:,:,jz)**2+S33(:,:,jz)**2+&
            2._rprec*(S12(:,:,jz)**2+S13(:,:,jz)**2+S23(:,:,jz)**2)))        
            
    ! Select Sij for this level for test-filtering, saving results as Sij_bar
    !   note: Sij is already on w-nodes
        S11_bar(:,:) = S11(:,:,jz)  
        S12_bar(:,:) = S12(:,:,jz)  
        S13_bar(:,:) = S13(:,:,jz)  
        S22_bar(:,:) = S22(:,:,jz)  
        S23_bar(:,:) = S23(:,:,jz)  
        S33_bar(:,:) = S33(:,:,jz)
       
        S11_hat = S11_bar
        S12_hat = S12_bar
        S13_hat = S13_bar
        S22_hat = S22_bar
        S23_hat = S23_bar
        S33_hat = S33_bar       

        call test_filter(S11_bar,G_test)
        call test_filter(S12_bar,G_test)
        call test_filter(S13_bar,G_test)
        call test_filter(S22_bar,G_test)
        call test_filter(S23_bar,G_test)
        call test_filter(S33_bar,G_test)

        call test_filter(S11_hat,G_test_test)
        call test_filter(S12_hat,G_test_test)
        call test_filter(S13_hat,G_test_test)
        call test_filter(S22_hat,G_test_test)
        call test_filter(S23_hat,G_test_test)
        call test_filter(S33_hat,G_test_test)
        
    ! Calculate |S_bar| (the test-filtered Sij)      
        S_bar = sqrt(2._rprec*(S11_bar**2 + S22_bar**2 + S33_bar**2 +&
            2._rprec*(S12_bar**2 + S13_bar**2 + S23_bar**2)))

    ! Calculate |S_hat| (the test-test-filtered Sij)     
        S_hat = sqrt(2._rprec*(S11_hat**2 + S22_hat**2 + S33_hat**2 +&
            2._rprec*(S12_hat**2 + S13_hat**2 + S23_hat**2)))

    ! Calculate |S|Sij then test-filter this quantity 
        S_S11_bar(:,:) = S(:,:)*S11(:,:,jz)
        S_S12_bar(:,:) = S(:,:)*S12(:,:,jz)
        S_S13_bar(:,:) = S(:,:)*S13(:,:,jz)
        S_S22_bar(:,:) = S(:,:)*S22(:,:,jz)
        S_S23_bar(:,:) = S(:,:)*S23(:,:,jz)
        S_S33_bar(:,:) = S(:,:)*S33(:,:,jz)

        S_S11_hat(:,:) = S_S11_bar(:,:)
        S_S12_hat(:,:) = S_S12_bar(:,:)
        S_S13_hat(:,:) = S_S13_bar(:,:)
        S_S22_hat(:,:) = S_S22_bar(:,:)
        S_S23_hat(:,:) = S_S23_bar(:,:)
        S_S33_hat(:,:) = S_S33_bar(:,:)

        call test_filter(S_S11_bar,G_test)
        call test_filter(S_S12_bar,G_test)
        call test_filter(S_S13_bar,G_test)
        call test_filter(S_S22_bar,G_test)
        call test_filter(S_S23_bar,G_test)
        call test_filter(S_S33_bar,G_test)     

        call test_filter(S_S11_hat,G_test_test)
        call test_filter(S_S12_hat,G_test_test)
        call test_filter(S_S13_hat,G_test_test)
        call test_filter(S_S22_hat,G_test_test)
        call test_filter(S_S23_hat,G_test_test)
        call test_filter(S_S33_hat,G_test_test)  

    ! Calculate Mij and Nij          
        M11 = const*(S_S11_bar - tf1_2*S_bar*S11_bar)
        M12 = const*(S_S12_bar - tf1_2*S_bar*S12_bar)
        M13 = const*(S_S13_bar - tf1_2*S_bar*S13_bar)
        M22 = const*(S_S22_bar - tf1_2*S_bar*S22_bar)
        M23 = const*(S_S23_bar - tf1_2*S_bar*S23_bar)
        M33 = const*(S_S33_bar - tf1_2*S_bar*S33_bar)

        N11 = const*(S_S11_hat - tf2_2*S_hat*S11_hat)
        N12 = const*(S_S12_hat - tf2_2*S_hat*S12_hat)
        N13 = const*(S_S13_hat - tf2_2*S_hat*S13_hat)
        N22 = const*(S_S22_hat - tf2_2*S_hat*S22_hat)
        N23 = const*(S_S23_hat - tf2_2*S_hat*S23_hat)
        N33 = const*(S_S33_hat - tf2_2*S_hat*S33_hat)

    ! Calculate LijMij, MijMij, etc for each point in the plane        
        LM = L11*M11+L22*M22+L33*M33+2._rprec*(L12*M12+L13*M13+L23*M23)
        MM = M11**2+M22**2+M33**2+2._rprec*(M12**2+M13**2+M23**2)
        QN = Q11*N11+Q22*N22+Q33*N33+2._rprec*(Q12*N12+Q13*N13+Q23*N23)
        NN = N11**2+N22**2+N33**2+2._rprec*(N12**2+N13**2+N23**2)

    ! Calculate ee_now (the current value of eij*eij)
            $if ($DYN_TN)       
            ee_now = L11**2+L22**2+L33**2+2._rprec*(L12**2+L13**2+L23**2) &
                    -2._rprec*LM*Cs_opt2(:,:,jz) + MM*Cs_opt2(:,:,jz)**2
            $endif           
        
    ! Initialize (???)        
        if (inilag) then
          if ((.not. F_LM_MM_init) .and. (jt == cs_count .or. jt == DYN_init)) then
             print *,'F_MM and F_LM initialized' 
             F_MM (:,:,jz) = MM
             F_LM (:,:,jz) = 0.03_rprec*MM
             F_MM(ld-1:ld,:,jz)=1._rprec
             F_LM(ld-1:ld,:,jz)=1._rprec

             if (jz == nz) F_LM_MM_init = .true.
          end if
        end if

    ! Inflow (???)        
        if(inflow)then
           iend = floor (buff_end * nx + 1._rprec)
           iend = modulo (iend - 1, nx) + 1
           istart = floor ((buff_end - buff_len) * nx + 1._rprec)
           istart = modulo (istart - 1, nx) + 1
           
           Tn=merge(.1_rprec*const*S**2,MM,MM.le..1_rprec*const*S**2)
           MM=Tn
           LM(istart + 1:iend, 1:ny) = 0._rprec
           F_LM(istart + 1:iend, 1:ny, jz) = 0._rprec
           Tn=merge(.1_rprec*const*S**2,NN,NN.le..1_rprec*const*S**2)
           NN=Tn
           QN(istart + 1:iend, 1:ny) = 0._rprec
           F_QN(istart + 1:iend, 1:ny, jz) = 0._rprec
        endif

    ! Update running averages (F_LM, F_MM)
        ! Determine averaging timescale (for 2-delta filter)
            $if ($DYN_TN)   ! based on Taylor timescale
                Tn = 4._rprec*pi*sqrt(F_ee2(:,:,jz)/F_deedt2(:,:,jz))   
            $else           ! based on Meneveau, Lund, and Cabot paper (JFM 1996)
                Tn = max (F_LM(:,:,jz) * F_MM(:,:,jz), real(zero))
                Tn = opftdelta*(Tn**powcoeff)
                ! Clip, if necessary
                Tn(:,:) = max(real(zero),real(Tn(:,:)))                
            $endif           
            
        ! Calculate new running average = old*(1-epsi) + instantaneous*epsi                 
            dumfac = lagran_dt/Tn 
            epsi = dumfac / (1._rprec+dumfac)

            F_LM(:,:,jz)=(epsi*LM + (1._rprec-epsi)*F_LM(:,:,jz))
            F_MM(:,:,jz)=(epsi*MM + (1._rprec-epsi)*F_MM(:,:,jz))
            ! Clip, if necessary
            F_LM(:,:,jz)= max(real(zero),real(F_LM(:,:,jz)))

    ! Calculate Cs_opt2 (for 2-delta filter)
        ! Add +zero in demomenator to avoid division by identically zero
        Cs_opt2_2d(:,:) = F_LM(:,:,jz)/(F_MM(:,:,jz) + zero)
        Cs_opt2_2d(ld,:) = zero
        Cs_opt2_2d(ld-1,:) = zero
        ! Clip, if necessary
        Cs_opt2_2d(:,:)=max(real(zero),real(Cs_opt2_2d(:,:)))

    ! Initialize (???)           
        if (inilag) then
          if ((.not. F_QN_NN_init) .and. (jt == cs_count .or. jt == DYN_init)) then
            print *,'F_NN and F_QN initialized'
            F_NN (:,:,jz) = NN
            F_QN (:,:,jz) = 0.03_rprec*NN
            F_NN(ld-1:ld,:,jz)=1._rprec
            F_QN(ld-1:ld,:,jz)=1._rprec

            if (jz == nz) F_QN_NN_init = .true.
          end if
        end if

    ! Update running averages (F_QN, F_NN)
        ! Determine averaging timescale (for 4-delta filter)
            $if ($DYN_TN)   ! based on Taylor timescale
                ! Keep the same as 2-delta filter 
            $else           ! based on Meneveau, Cabot, Lund paper (JFM 1996)
                Tn =max(real(F_QN(:,:,jz)*F_NN(:,:,jz)),real(zero))
                Tn=opftdelta*(Tn**powcoeff)
                ! Clip, if necessary
                Tn(:,:) = max(real(zero),real(Tn(:,:)))                     
            $endif   

        ! Calculate new running average = old*(1-epsi) + instantaneous*epsi                
            dumfac = lagran_dt/Tn
            epsi = dumfac / (1._rprec+dumfac)

            F_QN(:,:,jz)=(epsi*QN + (1._rprec-epsi)*F_QN(:,:,jz))
            F_NN(:,:,jz)=(epsi*NN + (1._rprec-epsi)*F_NN(:,:,jz))
            ! Clip, if necessary
            F_QN(:,:,jz)= max(real(zero),real(F_QN(:,:,jz)))

    ! Calculate Cs_opt2 (for 4-delta filter)
       ! Add +zero in demomenator to avoid division by identically zero
        Cs_opt2_4d(:,:) = F_QN(:,:,jz)/(F_NN(:,:,jz) + zero)
        Cs_opt2_4d(ld,:) = zero
        Cs_opt2_4d(ld-1,:) = zero
        ! Clip, if necessary
        Cs_opt2_4d(:,:)=max(real(zero),real(Cs_opt2_4d(:,:)))

    ! Calculate Beta and count how many are clipped
        ! Add +zero in demomenator to avoid division by identically zero
        Beta(:,:,jz)=&
             (Cs_opt2_4d(:,:)/(Cs_opt2_2d(:,:) + zero))**(log(tf1)/(log(tf2)-log(tf1)))
        counter1=0      
        counter2=0
        
        do jx=1,Nx
        do jy=1,Ny
           if (Beta(jx,jy,jz).LE.1/(tf1*tf2)) then
              counter1=counter1+1
           end if
        end do
        end do

        !--MPI: this is not valid
        if ( ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) .and.  &
             (jz == nz) ) then
          Beta(:,:,jz)=1._rprec
        end if
        
    ! Clip Beta and set Cs_opt2 for each point in the plane
        do jy = 1, ny
          do jx = 1, ld  !--perhaps only nx is needed
            Betaclip=max(real(Beta(jx,jy,jz)),real(1._rprec/(tf1*tf2)))
            Cs_opt2(jx,jy,jz)=Cs_opt2_2d(jx,jy)/Betaclip
          end do
        end do
        Cs_opt2(ld,:,jz) = zero
        Cs_opt2(ld-1,:,jz) = zero
        ! Clip, if necessary
        Cs_opt2(:,:,jz)=max(real(zero),real(Cs_opt2(:,:,jz)))

    ! Save planar averages every 200 timesteps (not currently written anywhere?)
    !    if (mod(jt,200) == 0) then
    !       LMvert(jz) = sum(sum(LM,DIM=1),DIM=1)/ny/nx
    !       MMvert(jz) = sum(sum(MM,DIM=1),DIM=1)/ny/nx
    !       QNvert(jz) = sum(sum(QN,DIM=1),DIM=1)/ny/nx
    !       NNvert(jz) = sum(sum(NN,DIM=1),DIM=1)/ny/nx
    !    end if
         
end do
! this ends the main jz=1,nz loop     -----------------------now repeat for other horiz slices

$if ($LVLSET)
  ! Zero Cs_opt2 inside objects
  call level_set_Cs_lag_dyn ()
$endif

$if ($VERBOSE)
write (*, *) 'finished lagrange_Sdep'
$endif

end subroutine lagrange_Sdep
