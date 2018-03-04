!! Calculate the sub-grid scale diffusivity of scalar quantities
subroutine scalars_sgs(scalar,dsdx,dsdy,dsdz,Ds_opt2,kappa_t,s_Tn_all,s_Beta,I_LM,I_MM,I_QN,I_NN,ds2_clips,I_LM_MM_init,I_QN_NN_init,sigma_s)

use types, only: rprec
use param
use sgs_stag_util, only: calc_Sij
use sgs_param, only: S11,S12,S13,S22,S23,S33,S
use sgs_param, only: delta,s_lagran_dt
use s_lagrange_Sdep, only: scalars_lagrange_Sdep

implicit none

integer :: jx,jy,jz
real(rprec), dimension(ld,ny,lbz:nz), intent(in) :: scalar 
real(rprec), dimension(ld,ny,lbz:nz), intent(in) :: dsdx 
real(rprec), dimension(ld,ny,lbz:nz), intent(in) :: dsdy 
real(rprec), dimension(ld,ny,lbz:nz), intent(in) :: dsdz 
real(rprec), dimension(nz) :: Ds_1D
real(rprec), dimension(ld,ny,nz), intent(inout) :: Ds_opt2
real(rprec), dimension(ld,ny,nz), intent(out) :: kappa_t
real(rprec), dimension(ld,ny,lbz:nz), intent(inout) :: I_LM, I_MM, I_QN, I_NN
real(rprec), dimension(ld,ny,lbz:nz), intent(inout) :: s_Tn_all, s_Beta, ds2_clips
logical, intent(inout) :: I_LM_MM_init
logical, intent(inout) :: I_QN_NN_init
real(rprec), dimension(lbz:nz), intent(inout) :: sigma_s

$if ($CFL_DT)
    if (scalars_sgs_model == 5) then
       if ( (jt.ge.ssgs_init-cs_count+1) ) then
             s_lagran_dt = lagrand_dt_s + dt
       end if
    end if
$else
    s_lagran_dt = cs_count*dt
$endif

!Re-calculate Sij components
call calc_Sij

if ( (jt.ge.ssgs_init).and.(mod(jt_total,cs_count)==0) ) then
!   if (jt == ssgs_init ) print*, 'running dynamics sgs model for scalars = ', scalars_sgs_model            
   if (scalars_sgs_model == 3) then
      call scalars_scaledep_dynamic(jt,scalar,dsdx,dsdy,dsdz,Ds_1D)
      do jz=1,nz
         Ds_opt2(:,:,jz) = Ds_1D(jz)
      end do
   else if (scalars_sgs_model == 5) then 
           call scalars_lagrange_Sdep(scalar,dsdx,dsdy,dsdz,Ds_opt2,s_Tn_all,s_Beta,I_LM,I_MM,I_QN,I_NN,ds2_clips,I_LM_MM_init,I_QN_NN_init,sigma_s)
   end if 
end if

!Calculate eddy diffusivity kappa
!Stored on w-nodes for entire domain except
!on uvp node for jz=1 and 'wall' BC
do jz=1,nz
   do jy=1,ny
      do jx=1,nx
         S(jx,jy) = sqrt( 2._rprec*(S11(jx,jy,jz)**2 +           S22(jx,jy,jz)**2 +&
                                    S33(jx,jy,jz)**2 + 2._rprec*(S12(jx,jy,jz)**2 +&
                                    S13(jx,jy,jz)**2 +           S23(jx,jy,jz)**2 )))
         kappa_t(jx,jy,jz)=S(jx,jy)*Ds_opt2(jx,jy,jz)*delta**2
      end do
   end do
end do          

end subroutine scalars_sgs
