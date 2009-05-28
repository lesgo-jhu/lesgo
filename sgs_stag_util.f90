! put everything onto w-nodes, follow original version
!--provides txx, txy, tyy, tzz for jz=1:nz-1; txz, tyz for 1:nz
subroutine sgs_stag ()
use types,only:rprec
use param
use param2
use sim_param,only: u,v,w,dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,  &
                    txx, txy, txz, tyy, tyz, tzz
use sgsmodule,only:u_lag,v_lag,w_lag,Cs_opt2,Nu_t
!------------------------- Vij Comment begins---------------
! 04/14/2004 - Added Nu_t to the list of variables from sgsmodule
! 05/19/2004 - Replaced all references to visc by Nu_t; deleted local 
!              declaration of visc
! 05/19/2004 - Replace all .5 by 0.5
!-------------------------Vij Comment ends ------------------
use bottombc,only:zo
use immersedbc,only:building_mask,building_interp
use test_filtermodule,only:filter_size
use messages
use debug_mod
$if ($LVLSET)
  use level_set, only : level_set_BC, level_set_Cs, level_set_smooth_vel
$endif
implicit none

character (*), parameter :: sub_name = 'sgs_stag'

logical, parameter :: DEBUG = .false.

real(kind=rprec),dimension(nz)::l,ziko,zz
!real (rprec), dimension (ld, ny, nz) :: S11, S12, S22, S33, S13, S23
real(kind=rprec),dimension(ld,ny,nz):: dissip
real(kind=rprec),dimension(ld,ny) :: txzp, tyzp,S
real(kind=rprec) :: delta, nu, const
$if ($LVLSET)
  !--experimental: holds derivs of smoothed velocity field
  !real (rprec), dimension (ld, ny, nz) :: u_s, v_s, w_s
  !real (rprec), dimension (ld, ny, nz) :: dudx_s, dudy_s, dudz_s,  &
  !                                        dvdx_s, dvdy_s, dvdz_s,  &
  !                                        dwdx_s, dwdy_s, dwdz_s
$endif 

integer::jx,jy,jz
integer :: jz_min

if (VERBOSE) call enter_sub (sub_name)

delta=filter_size*(dx*dy*dz)**(1._rprec/3._rprec) ! nondimensional
! Cs is Smagorinsky's constant. l is a filter size (non-dim.)  

if (DEBUG) then
  call DEBUG_write (dudx(:, :, 1:nz), 'sgs_stag.dudx.a')
  call DEBUG_write (dudy(:, :, 1:nz), 'sgs_stag.dudy.a')
  call DEBUG_write (dudz(:, :, 1:nz), 'sgs_stag.dudz.a')
  call DEBUG_write (dvdx(:, :, 1:nz), 'sgs_stag.dvdx.a')
  call DEBUG_write (dvdy(:, :, 1:nz), 'sgs_stag.dvdy.a')
  call DEBUG_write (dvdz(:, :, 1:nz), 'sgs_stag.dvdz.a')
  call DEBUG_write (dwdx(:, :, 1:nz), 'sgs_stag.dwdx.a')
  call DEBUG_write (dwdy(:, :, 1:nz), 'sgs_stag.dwdy.a')
  call DEBUG_write (dwdz(:, :, 1:nz), 'sgs_stag.dwdz.a')
end if

!$if (LVLSET)
  !--experimental: smooth velocity, recalculate derivatives
  !--this is Expensive: ram and cpu, and IMHO not a good idea
  !u_s = u
  !v_s = v
  !w_s = w
  !call level_set_smooth_vel (u_s, v_s, w_s)
  
  !call filt_da (u_s, dudx_s, dudy_s)
  !call filt_da (v_s, dvdx_s, dvdy_s)
  !call filt_da (w_s, dwdx_s, dwdy_s)

  !call ddz_uv (dudz_s, u_s)
  !call ddz_uv (dvdz_s, v_s)

  !call ddz_w (dwdz_s, w_s)
  
  !call calc_Sij (dudx_s, dudy_s, dudz_s,  &
  !               dvdx_s, dvdy_s, dvdz_s,  &
  !               dwdx_s, dwdy_s, dwdz_s)

!$else
  call calc_Sij ()
 
				 
!$endif

if (DEBUG) then
  call DEBUG_write (S11(:, :, 1:nz), 'sgs_stag.S11.b')
  call DEBUG_write (S12(:, :, 1:nz), 'sgs_stag.S12.b')
  call DEBUG_write (S13(:, :, 1:nz), 'sgs_stag.S13.b')
  call DEBUG_write (S22(:, :, 1:nz), 'sgs_stag.S22.b')
  call DEBUG_write (S23(:, :, 1:nz), 'sgs_stag.S23.b')
  call DEBUG_write (S33(:, :, 1:nz), 'sgs_stag.S33.b')
end if

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
  ! save the wall level
  txzp(:,:)=txz(:,:,1)
  tyzp(:,:)=tyz(:,:,1)
end if

! This part computes the average velocity during cs_count times steps
! This is used with the lagrangian model only
if (model == 4 .OR. model==5) then
   u_lag = u_lag+u
   v_lag = v_lag+v
   w_lag = w_lag+w
end if

if (sgs) then!ref01
   if((model == 1))then  !For traditional Smagorinsky   ref02


     $if ($LVLSET)
 
       l = delta
       call level_set_Cs (delta)
 
     $else

       ! Define parameters (Co and nn) for wallfunction
       Cs_opt2 = Co**2  ! constant coefficient

       if (lbc_mom == 'stress free') then

         l = delta
         
       else
       
         if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
           do jz = 1, nz

             if (jz == 1) then
                zz(jz) = (jz - 0.5_rprec) * dz
             else  ! w-nodes
                zz(jz) = (jz - 1) * dz
             end if

             ! z's nondimensional, l here is on w-nodes, except at bottom
             l(jz) = ( Co**(nnn)*(vonk*zz(jz))**(-nnn) +  &
                       (delta)**(-nnn) )**(-1._rprec/nnn)
           end do
         else
           do jz = 1, nz
             zz(jz) = ((jz - 1) + coord * (nz - 1)) * dz  !--w-nodes

             ! z's nondimensional, l here is on w-nodes, except at bottom
             l(jz) = ( Co**(nnn)*(vonk*zz(jz))**(-nnn) +  &
                       (delta)**(-nnn) )**(-1._rprec/nnn)

           end do 
         end if

       end if

     $endif

   else ! for dynamic procedures: initialization  !cont ref02
   
     l = delta  ! constant equal to delta
     if ((jt == 1) .and. (inilag)) then !ref05
       print *,'CS_opt2 initialiazed'
       Cs_opt2 = 0.03_rprec

    ! make sure that the "node conventions" in these dynamic models
    ! match with those done here
    elseif ( ((jt.GE.DYN_init).OR.(initu)) .AND.  &
             (mod(jt,cs_count)==0) ) then!cont ref05
             
      if (jt ==DYN_init) print *,'running dynamic model = ',model
      if (model == 2) then  ! Standard dynamic model !ref06
        !ziko = Cs_opt2(5,5,:)  !--what does this do?
        !--provides ziko 1:nz
        call std_dynamic(ziko,S11,S12,S13,S22,S23,S33)
        forall (jz = 1:nz) Cs_opt2(:, :, jz) = ziko(jz)
      else if (model==3) then !cont ref06
        ! Plane average dynamic  continue ref 05
        !ziko = Cs_opt2(5,5,:)  !--what does this do?
        call scaledep_dynamic(ziko,S11,S12,S13,S22,S23,S33)
        do jz = 1, nz
          Cs_opt2(:, :, jz) = ziko(jz)
        end do
      else if (model==4.) then ! Lagrangian SS !cont ref06
        call lagrange_Ssim(S11,S12,S13,S22,S23,S33)
      elseif (model==5) then ! LagSD !cont ref06
        call lagrange_Sdep(S11,S12,S13,S22,S23,S33)
      end if !end ref06
      
    end if !end ref05
    ! for if jt == 1 and 2 else ifs end ref 5
  end if !(for model =1 or else) end ref02   
end if ! for sgs-if end ref01

if (DEBUG) then
  call DEBUG_write (Cs_opt2, 'sgs_stag.Cs_opt2')
end if

! define |S| and viscosity on w-nodes (on uvp node for jz=1)
!$comp parallel do default(shared) private(jx,jy,jz)
!--MPI: going up to nz saves communication
do jz = 1, nz
do jy=1,ny
do jx=1,nx
   S(jx,jy) = sqrt(2._rprec*(S11(jx,jy,jz)**2 + S22(jx,jy,jz)**2 +&
        S33(jx,jy,jz)**2 + 2._rprec*(S12(jx,jy,jz)**2 +&
        S13(jx,jy,jz)**2 + S23(jx,jy,jz)**2)))
   Nu_t(jx,jy,jz)=S(jx,jy)*Cs_opt2(jx,jy,jz)*(l(jz)**2)
   dissip(jx,jy,jz)=-(S(jx,jy)**3)*Cs_opt2(jx,jy,jz)*(l(jz)**2)
end do
end do

!--should have option to turn off
!if (mod(jt,100) == 0) then
!      if(jz==1.or.jz==nz/4.or.jz==nz/2) then
!         write(90)real(jt*dt),real(jz),real(Cs_opt2(1:NX,1:NY,jz))
!         write(93)real(jt*dt),real(jz),real(Nu_t(1:NX,1:NY,jz))
!         write(94)real(jt*dt),real(jz),real(dissip(1:NX,1:NY,jz))
!      end if
!end if
end do

!	if (mod(jt,50)==0) pause
!	pause
!$comp end parallel do

! evaluate nu_t= c_s^2 l^2 |S|
if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
  ! bottom

  select case (lbc_mom)

    case ('wall')
      !$comp parallel do default(shared) private(jx,jy,const)
      do jy=1,ny
      do jx=1,nx
         nu = 0._rprec
         const = 0._rprec
         if (sgs) then
            const = Nu_t(jx,jy,1)
         end if
         if (molec) then
            nu = (nu_molec/(u_star*z_i)) ! dim/less 
         end if
      ! everything on right node here
         txx(jx,jy,1) = -2._rprec*(const+nu)*S11(jx,jy,1)
         txy(jx,jy,1) = -2._rprec*(const+nu)*S12(jx,jy,1)
         tyy(jx,jy,1) = -2._rprec*(const+nu)*S22(jx,jy,1)
         tzz(jx,jy,1) = -2._rprec*(const+nu)*S33(jx,jy,1)
      end do
      end do
      !$comp end parallel do

    case ('stress free')

      do jy=1,ny
      do jx=1,nx
         nu=0._rprec
         const=0._rprec
         if (sgs) then
            const=0.5_rprec*(Nu_t(jx,jy,1) + Nu_t(jx,jy,2))
         end if
         if (molec) then
            nu=(nu_molec/(u_star*z_i)) ! dim/less
         end if
         txx(jx,jy,1)=-2._rprec*(const+nu)*&
               0.5_rprec*(S11(jx,jy,1) + S11(jx,jy,2))
         txy(jx,jy,1)=-2._rprec*(const+nu)*&
              0.5_rprec*(S12(jx,jy,1) + S12(jx,jy,2))
         tyy(jx,jy,1)=-2._rprec*(const+nu)*&
              0.5_rprec*(S22(jx,jy,1) + S22(jx,jy,2))
         tzz(jx,jy,1)=-2._rprec*(const+nu)*&
              0.5_rprec*(S33(jx,jy,1) + S33(jx,jy,2))
      end do
      end do

  end select
  
  jz_min = 2  !--wall level set by wallstress routine

else

  jz_min = 1

end if

! middle
!$comp parallel do default(shared) private(jx,jy,jz,const)	
!--note that nz level is available for the interpolations at (jz+1) here
!--txx, tyy, tzz, txy not needed at nz
do jz=jz_min, nz - 1
do jy=1,ny
do jx=1,nx
   nu=0._rprec
   const=0._rprec
   if (sgs) then
      const=0.5_rprec*(Nu_t(jx,jy,jz) + Nu_t(jx,jy,jz+1))
   end if
   if (molec) then
      nu=(nu_molec/(u_star*z_i)) ! dim/less
   end if
   txx(jx,jy,jz)=-2._rprec*(const+nu)*&
        0.5_rprec*(S11(jx,jy,jz) + S11(jx,jy,jz+1))
   txy(jx,jy,jz)=-2._rprec*(const+nu)*&
        0.5_rprec*(S12(jx,jy,jz) + S12(jx,jy,jz+1))
   tyy(jx,jy,jz)=-2._rprec*(const+nu)*&
        0.5_rprec*(S22(jx,jy,jz) + S22(jx,jy,jz+1))
   tzz(jx,jy,jz)=-2._rprec*(const+nu)*&
        0.5_rprec*(S33(jx,jy,jz) + S33(jx,jy,jz+1))
end do
end do
end do
!$comp end parallel do

!$comp parallel do default(shared) private(jx,jy,jz,const)	
do jz = jz_min, nz - 1
do jy=1,ny
do jx=1,nx
   nu = 0._rprec
   const = 0._rprec
   if (sgs) const=Nu_t(jx,jy,jz)
   if (molec) nu=(nu_molec/(u_star*z_i)) ! dim/less
   txz(jx,jy,jz)=-2._rprec*(const + nu) * S13(jx,jy,jz)
   tyz(jx,jy,jz)=-2._rprec*(const + nu) * S23(jx,jy,jz)
end do
end do
end do
!$comp end parallel do

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
  ! don't need to do this, unless we accidentally mess with the wall level
  txz(:,:,1)=txzp(:,:)
  tyz(:,:,1)=tyzp(:,:)
end if

if (DEBUG) then
  call DEBUG_write (txx(:, :, 1:nz), 'sgs_stag.txx.a')
  call DEBUG_write (txy(:, :, 1:nz), 'sgs_stag.txy.a')
  call DEBUG_write (txz(:, :, 1:nz), 'sgs_stag.txz.a')
end if

$if ($LVLSET)
  !--at this point tij are only set for 1:nz-1
  !--at this point u, v, w are set for 0:nz, except bottom process
  !  is 1:nz
  !--some MPI synchronizing may be done in here, but this will be kept
  !  separate from the rest of the code (at the risk of some redundancy)
  call level_set_BC ()
$endif

if (DEBUG) then
  call DEBUG_write (txx(:, :, 1:nz), 'sgs_stag.txx.b')
  call DEBUG_write (txy(:, :, 1:nz), 'sgs_stag.txy.b')
  call DEBUG_write (txz(:, :, 1:nz), 'sgs_stag.txz.b')
end if

$if ($MPI)
  !--recv information for top nodes: txy, txz only
  !--other components not needed, since no equation is solved there
  !  (values are just copied)
  call mpi_sendrecv (txz(1, 1, 1), ld*ny, MPI_RPREC, down, 3,  &
                     txz(1, 1, nz), ld*ny, MPI_RPREC, up, 3,   &
                     comm, status, ierr)
  call mpi_sendrecv (tyz(1, 1, 1), ld*ny, MPI_RPREC, down, 4,  &
                     tyz(1, 1, nz), ld*ny, MPI_RPREC, up, 4,   &
                     comm, status, ierr)

  !--also set 0-layer to bogus values
  txx(:, :, 0) = BOGUS
  txy(:, :, 0) = BOGUS
  txz(:, :, 0) = BOGUS
  tyy(:, :, 0) = BOGUS
  tyz(:, :, 0) = BOGUS
  tzz(:, :, 0) = BOGUS  !--tzz is updated in main (move to here)

$endif

txx(:, :, nz) = BOGUS
txy(:, :, nz) = BOGUS
tyy(:, :, nz) = BOGUS
tzz(:, :, nz) = BOGUS
  
if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
  txz(:,:,nz)=0._rprec
  tyz(:,:,nz)=0._rprec
end if

if (VERBOSE) call exit_sub (sub_name)

return
end subroutine sgs_stag
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_Sij ()
use Sij_defs 
use types, only : rprec
use sim_param,only:dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
implicit none

integer::jx,jy,jz
integer :: jz_min

$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif

!real (rprec), dimension (ld, ny, $lbz:nz) :: dudx, dudy, dudz,  &
!                                             dvdx, dvdy, dvdz,  &
!                                             dwdx, dwdy, dwdz

real (rprec) :: ux, uy, uz, vx, vy, vz, wx, wy, wz
!---------------------------------------------------------------------                                     
if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then

  ! calculate Sij on w-nodes
  ! calculate |S| on w-nodes

  select case (lbc_mom)

    case ('wall')  !--u-node for first point
    
      ! at level z=dz/2.  (Where dudz and dvdz are on UVP nodes)
      do jy=1,ny
      do jx=1,nx              
         ux=dudx(jx,jy,1)  ! uvp-node
         uy=dudy(jx,jy,1)  ! uvp-node
         uz=dudz(jx,jy,1)  ! uvp-node
         vx=dvdx(jx,jy,1)  ! uvp-node
         vy=dvdy(jx,jy,1)  ! uvp-node
         vz=dvdz(jx,jy,1)  ! uvp-node 
         ! special case
         wx=0.5_rprec*(dwdx(jx,jy,1)+dwdx(jx,jy,2))  ! uvp-node
         wy=0.5_rprec*(dwdy(jx,jy,1)+dwdy(jx,jy,2))  ! uvp-node
         wz=dwdz(jx,jy,1)  ! uvp-node
         S11(jx,jy,1)=ux          ! uvp-node
         S12(jx,jy,1)=0.5_rprec*(uy+vx) ! uvp-node
         ! taken care of with wall stress routine
         S13(jx,jy,1)=0.5_rprec*(uz+wx) ! uvp
         S22(jx,jy,1)=vy          ! uvp-node
         ! taken care of with wall stress routine 
         S23(jx,jy,1)=0.5_rprec*(vz+wy) ! uvp
         S33(jx,jy,1)=wz          ! uvp-node
      end do
      end do
  
    case ('stress free')  !--w-nodes like usual, using BC assumptions

      do jy=1,ny
      do jx=1,nx              
         ux=0.5_rprec*(dudx(jx,jy,1) + dudx(jx,jy,1))  ! w-node
         uy=0.5_rprec*(dudy(jx,jy,1) + dudy(jx,jy,1))  ! w-node
         uz=dudz(jx,jy,1)  ! w-node, 0
         vx=0.5_rprec*(dvdx(jx,jy,1) + dvdx(jx,jy,1))  ! w-node
         vy=0.5_rprec*(dvdy(jx,jy,1) + dvdy(jx,jy,1))  ! w-node
         vz=dvdz(jx,jy,1)  ! w-node, 0
         wx=dwdx(jx,jy,1)  ! w-node, 0
         wy=dwdy(jx,jy,1)  ! w-node, 0
         !--not sure best way to estimate dwdz (w-node)
         !--perhaps use odd BC on w?
         !--same problem occurs at top bdry
         !--for now, this is consistent with way its done at top
         wz=0.5_rprec*(dwdz(jx,jy,1) + 0._rprec) ! w-node
                                      !--this is what is done at top,
                                      !  since ddz_w () implicity does this
         S11(jx,jy,1)=ux          ! w-node
         S12(jx,jy,1)=0.5_rprec*(uy+vx) ! w-node
         S13(jx,jy,1)=0.5_rprec*(uz+wx) ! w-node
         S22(jx,jy,1)=vy          ! w-node
         S23(jx,jy,1)=0.5_rprec*(vz+wy) ! w-node
         S33(jx,jy,1)=wz          ! w-node
      end do
      end do

  end select
  
  jz_min = 2

else

  jz_min = 1

end if

$if ($MPI)
  !--this is only required b/c of the unnatural placement of all strains
  !  onto w-nodes, be careful not to overwrite nz on top process with garbage
  !--dwdz(jz=0) is already known, except at bottom process (OK)
  !call mpi_sendrecv (dwdz(1, 1, nz-1), ld*ny, MPI_RPREC, up, 1,  &
  !                   dwdz(1, 1, 0), ld*ny, MPI_RPREC, down, 1,   &
  !                   comm, status, ierr)
  call mpi_sendrecv (dwdz(1, 1, 1), ld*ny, MPI_RPREC, down, 2,  &
                     dwdz(1, 1, nz), ld*ny, MPI_RPREC, up, 2,   &
                     comm, status, ierr)
$endif

! calculate derivatives/strain on w-nodes
!$vvohmygod parallel do default(shared) private(jx,jy,jz)	
!--in MPI version, calculating up to nz saves some interprocess exchange
!  later but note that dwdz is not provided w/o some communication
!  (unless its the top process) 
do jz=jz_min, nz
do jy=1,ny
do jx=1,nx              
   ux=0.5_rprec*(dudx(jx,jy,jz) + dudx(jx,jy,jz-1))  ! w-node
   uy=0.5_rprec*(dudy(jx,jy,jz) + dudy(jx,jy,jz-1))  ! w-node
   uz=dudz(jx,jy,jz)  ! w-node
   vx=0.5_rprec*(dvdx(jx,jy,jz) + dvdx(jx,jy,jz-1))  ! w-node
   vy=0.5_rprec*(dvdy(jx,jy,jz) + dvdy(jx,jy,jz-1))  ! w-node
   vz=dvdz(jx,jy,jz)  ! w-node
   wx=dwdx(jx,jy,jz)  ! w-node
   wy=dwdy(jx,jy,jz)  ! w-node
   wz=0.5_rprec*(dwdz(jx,jy,jz) + dwdz(jx,jy,jz-1))  ! w-node
   S11(jx,jy,jz)=ux          ! w-node
   S12(jx,jy,jz)=0.5_rprec*(uy+vx) ! w-node
   S13(jx,jy,jz)=0.5_rprec*(uz+wx) ! w-node
   S22(jx,jy,jz)=vy          ! w-node
   S23(jx,jy,jz)=0.5_rprec*(vz+wy) ! w-node
   S33(jx,jy,jz)=wz          ! w-node
end do
end do
end do
!$ffohmygod end parallel do

return
end subroutine calc_Sij
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

