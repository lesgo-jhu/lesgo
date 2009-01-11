module scalars_module2
use types,only:rprec
use param2, jt_global => jt  !--rename jt to avoid naming conflicts here
                            !--only needed since jt added to params
use bottombc,only:T_s,q_s,q_mix,zo_avg,phi_m,psi_m,phi_h,psi_h
!makes obukhov functions available
use sim_param,only:u,v,w,theta,q
implicit none

! Preserve averaged variable outputs in memory for scalar_slice subroutine
!real(kind=rprec),dimension(nx,nz):: atheta,aq,t2,q2,asgs_t3,asgs_q3,awt,awq,adTdz,adqdz,abeta,anu_t

!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
! Part II of the scalar files .. also look at scalars_module.f90
!contains subroutines:
! ic_scal --- initialize velocity fields and scalars
! patch_or_remote - initialize surface boundary conditions
! scalar_in - read surface data from an external data file (remote-sensed)
! append_zeros_string - append zeros in front of a string - called by scalar_in (NOT USED !!)
! scalar_slice - outputs time-averaged x-z slices of scalar variables
! controlled by c_count & p_count from param; Uses file unit numbers from (36-47)
! obukhov_slice - outputs the obukhov variables (phi,psi, L,wstar);
! called fromn scalar_slice and toggled by parameter obukhov_output (specified in scalars_module.f)
! 
! Authored by Vijayant Kumar
! Last modified - April 25, 2004
!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX

contains

!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
subroutine ic_scal()
!subroutine ic_scal(u,v,w,theta,q,sign_wt_s,theta_mean)
!c...Log profile that is modified to flatten at z=z_i
!c .. Modified by Vijayant to put scalars back
! Last modified April 14, 2004
use types,only:rprec
use sim_param,only:u,v,w,theta,q
use bottombc
implicit none
real(kind=rprec),dimension(nz)::ubar
real(kind=rprec)::rms, noise, arg, arg2,theta_mean
real(kind=rprec)::z,w_star,T_star,q_star,ran3
integer::jx,jy,jz,seed

!      if (patch_flag .eq. 1) theta_mean=theta_s1
!      if (patch_flag .eq. 1) theta_mean=T_init
       theta_mean=T_init
      
      if (wt_s .eq. 0.0) then
! Compute the values of w_star etc. using the default value of
! wt_s = 0.06
      w_star=(9.81/theta_mean*0.06*z_i)**(1./3.)
! w_star is of O(1) with z_i=500 and wt_s=0.06
      T_star=0.06/w_star
      q_star=T_star
      else
      w_star=sign((9.81/theta_mean*abs(wt_s)*z_i)**(1./3.),wt_s)
      T_star=wt_s/w_star
      q_star=T_star
      end if
      seed=-112 ! Should be part of param actually.. change later
      print *,'w_star,T_star,q_star,seed',w_star,T_star,q_star,seed
        
!TS      zo_avg=sum(zo)/(nx*ny)

      print *,'Modified Log Profile for IC'
       do jz=1,nz
        z=(real(jz)-0.5)*dz*z_i
!c IC in equilibrium with rough surface (rough dominates in effective zo)
        arg2=z/zo_avg
        arg=(1./vonk)*log(arg2)!-1./(2.*vonk*z_i*z_i)*z*z
        if (coriolis_forcing) then
!        ubar(jz)=ug/u_star
        ubar(jz)=arg/30
        else
        ubar(jz)=arg
        end if
!C sc: I changed around some parenthesis here
        if (z.gt.(0.6*z_i)) then
!        print *, 'Statement executed for the scalars'
        ubar(jz)=ubar(jz-1)
        end if
       end do

!      do jz=1,nz
!       print *,'k, ubar:',jz,ubar(jz)
!      end do

      rms = 3.
      do jz=1,nz
        do jy=1,ny
          do jx=1,nx
!c...Ran3 returns uniform RV between 0 and 1. (sigma_rv=0.289)
!c...Taking std dev of vel as 1 at all heights

!cVK Note that if you put wt_s = 0 symbolizing neutal conditions
!c u should also put L_z=z_i i.e. the inversion layer height should
!c be equal to the height of the domain and in that case the second
!c part of the subsequent if loop will never execute. This is
!c ensured by putting an OR statement in the if clause, which makes 
!c sure that only the first part of if block is executed and not the
!c block after else

            z=(real(jz)-0.5)*dz*z_i
!          if ((z.le.z_i) .OR. (wt_s .eq. 0)) then
          if (z.le.z_i) then
           noise=rms/0.289*(ran3(seed)-0.5)
           u(jx,jy,jz)=noise*(1.-z/z_i)*w_star/u_star+ubar(jz)
           noise=rms/0.289*(ran3(seed)-0.5)
           v(jx,jy,jz)=noise*(1.-z/z_i)*w_star/u_star !noise
           noise=rms/0.289*(ran3(seed)-0.5)
           w(jx,jy,jz)=noise*(1.-z/z_i)*w_star/u_star
           noise=rms/0.289*(ran3(seed)-0.5)
           theta(jx,jy,jz)=(theta_mean+10.*noise*(1-z/z_i)*T_star)/T_scale
           noise=rms/0.289*(ran3(seed)-0.5)
           q(jx,jy,jz)=q_mix+50.*noise*(1-z/z_i)*q_star
          else
           noise=rms/0.289*(ran3(seed)-0.5)
           u(jx,jy,jz)=noise*w_star/u_star*0.01+ubar(jz)
           noise=rms/0.289*(ran3(seed)-0.5)
           v(jx,jy,jz)=noise*w_star/u_star*0.01
           noise=rms/0.289*(ran3(seed)-0.5)
           w(jx,jy,jz)=noise*w_star/u_star*0.01
! is this right for the neutral case? can't we put just theta=0., q=0. ?
           noise=rms/0.289*(ran3(seed)-0.5)
           theta(jx,jy,jz)=(theta_mean+(z-z_i)*inv_strength)/T_scale

           noise=rms/0.289*(ran3(seed)-0.5)
           q(jx,jy,jz)=q_mix+10.*noise*(1-z/z_i)*q_star
          end if
          end do
        end do
!       print *,'k, T_bar:',jz,sum(theta(:,:,jz))/(nx*ny)
      end do
      do jz=1,nz
!       print *,'k, ubar:',jz,sum(u(:,:,jz))/(nx*ny)
      end do

!c...BC for W
      do jy=1,ny
        do jx=1,nx
          w(jx,jy,1)=0._rprec
          w(jx,jy,nz)=0._rprec
        end do
      end do

!c...BC for U, V
      do jy=1,ny
        do jx=1,nx
          u(jx,jy,nz)=u(jx,jy,nz-1)
          v(jx,jy,nz)=v(jx,jy,nz-1)
        end do
      end do

do jz=1,nz
     write(6,7781) jz,sum(u(:,:,jz))/(nx*ny),sum(v(:,:,jz))/&
     (nx*ny),sum(w(:,:,jz))/(nx*ny),sum(theta(:,:,jz))/(nx*ny)
end do
7781 format('k, ubar, vbar, wbar,T_bar:',(1x,I3,1x,F9.4,1x,F9.4,1x,F9.4,1x,F9.4))

end subroutine ic_scal

!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
subroutine patch_or_remote()
use param
use param2
use bottombc
use scalars_module
implicit none
!z_os already defined in scalars_module
!zo,T_s and q_s defined in bottombc
! April 20, 2004 - so far contains both patch_or_remote
! and scalar_in subroutine

if (patch_flag .eq. 1) then
print *, 'Assigning values to patches'
!call patches(zo,T_s,q_s2d,patch,patchnum)
call patches()

!z_os already defined in scalars_module
!zo,T_s and q_s defined in bottombc
z_os(:,:)=(1./10.)*zo(:,:)

!c sets temperature field and roughness field at the bottom , x-y plane
!c Added by Vijayant
else if (remote_flag .eq. 1) then
print *, 'Assigning remote-sensed values to the surface'
   call scalar_in() ! Updated T_s and zo loaded from bottombc
!   call scalar_in(T_s,zo,crap2,crap3)
!   T_s=T_s+273 ! Convert to Kelvin
!   T_s_init_mean=sum(T_s)/(nx*ny)

!   if (remote_homog_flag .eq. 1) then
!print *,'Homogeinizing the remote b.c.s'
!T_s_init_mean=sum(T_s)/(nx*ny)
!zo=sum(zo)/(nx*ny)
!   end if
! Non-dimensionalize T
   T_s=(T_s)/T_scale
   zo=zo/z_i  
! z_os is the scalar roughness length. Divide momentum
! roughness length by 10 to get scalar roughness length
! data (look in the scalar_in routine above)
z_os(:,:)=(1./10.)*zo(:,:)
end if

end subroutine patch_or_remote


!!!xxxxxxxx--------VIJ---------XXXXXXXXXXXX-------------
!!!xxxxxxxx--------VIJ---------XXXXXXXXXXXX-------------
subroutine scalar_in()
!subroutine scalar_in(T_s,zo)
!subroutine scalar_in(T_s,zo,suffix2,suffix3)
!c This reads in the scalar input from a interpolated scalar
!c file and assigns it to the x-y plane at z=0 i.e. at the ground
!c This for the moment is used to read the momentum roughness, z0
!c and temperature from the USDA remote-sensed data set. The data 
!c can be interpolated using bilinear/cubic/spline interpolation (MATLAB)
!c for the grid size(nx,ny)
!c Authored by Vijayant Kumar
!c Last modified on April 11th, 2004
use param
use param2
use bottombc,only:T_s,zo !Load the variables from bottombc and update in here
       implicit none
       integer:: ii,jj
       character(len=6):: suffix2,suffix3
      
       write(suffix2,'(i6.6)') nx ! create a string from nx   
!       call append_zeros_string(suffix2,nx) ! Append leading zeros to string
       
      if (coarse_grain_flag) then
      write(suffix3,'(i6.6)') stencil_pts    
!      call append_zeros_string(suffix3,stencil_pts) 

 
open(unit=77,file='../interp_data/coarse_grained/interp_temp_cg_'&
//suffix2(4:6)//'pts_'//suffix3(4:6)//'.out',status='unknown')
open(unit=78,file='../interp_data/coarse_grained/interp_z_m_cg_'&
//suffix2(4:6)//'pts_'//suffix3(4:6)//'.out',status='unknown')
       
     
print *,'interp_temp_cg_'//suffix2(4:6)//'pts_'&
//suffix3(4:6)//'.out loaded from scalar_in.f'
       
      do jj=1,ny
          read(77,5169) (T_s(ii,jj),ii=1,nx)
          read(78,5169) (zo(ii,jj),ii=1,nx)
      enddo
      T_s=T_s+273.15 ! Convert T to Kelvin
      close(77)
      close(78)

       else

!open(unit=77,file='../interp_data/interp_temp_'//suffix2(4:6)&
!//'X'//suffix2(4:6)//'_cubic.out',status='unknown')
!open(unit=78,file='../interp_data/interp_z_m_'//suffix2(4:6)&
!//'X'//suffix2(4:6)//'_cubic.out',status='unknown')
!print *,'interp_temp_'//suffix2(4:6)//'X'//suffix2(4:6)//'_cubic.out'

open(unit=77,file='./interp_data/interp_temp_'//suffix2(4:6)&
//'X'//suffix2(4:6)//'_cubic.out',status='unknown')
open(unit=78,file='./interp_data/interp_z_m_'//suffix2(4:6)&
//'X'//suffix2(4:6)//'_cubic.out',status='unknown')
print *,'interp_temp_'//suffix2(4:6)//'X'//suffix2(4:6)//'_cubic.out'
       

      do jj=1,ny
!          do ii=1,nx
          read(77,5169) (T_s(ii,jj),ii=1,nx)
          read(78,5169) (zo(ii,jj),ii=1,nx)
      enddo
      print *,'Mean T_s = ',sum(T_s)/(nx*ny)
      close(77)
      close(78)
          T_s=T_s+273.15 ! Convert T to Kelvin
      end if
5169     format(1400(E16.10))
end subroutine scalar_in

!!!xxxxxxxxxxxx--------VIJ-------------XXXXXXXXXXXXXXXXXXX
!!!xxxxxxxxxxxx--------VIJ-------------XXXXXXXXXXXXXXXXXXX
!!!xxx-appends character zeros to a character string--XXXX
!!! NOT NEEDED ANYMORE - May 5, 2004

!!subroutine append_zeros_string(string_in,number_in)
!!integer :: number_in
!!character(len=6) string_in

!!if (number_in<100) then
!!    if (number_in<10) then
!!        string_in(4:6)='00'//string_in(6:6)
!!    else
!!        string_in(4:6)='0'//string_in(5:6)
!!    end if
!!end if

!!end subroutine append_zeros_string

!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
!!!xxxxxxxxxx-------scalar output subroutine-----XXXXXXXXXXXXXXXXXXXXX

subroutine scalar_slice(jt)
!subroutine scalar_slice(w,t,q,sgs_t3,sgs_q3,dTdz,beta_scal,Nu_t,jt)
!c This is exactly the same like the subroutine avgslice with the
!c only difference being that it averages the scalar variables
!c to find the y-averaged instantaneous x-z slices of variables
!c t,q,sgs_t3,sgs_q3 and their variances such as t2, q2.
!c It also outputs the average covariance between wt and wq
use sgsmodule,only: Nu_t 
use scalars_module 
implicit none
integer :: jt
integer:: i, j, k
real(kind=rprec),dimension(nx,nz),save:: atheta,aq,t2,q2,asgs_t3,asgs_q3,awt,awq
real(kind=rprec),dimension(nx,nz),save:: adTdz,adqdz,abeta,anu_t
real(kind=rprec):: ttheta1,tq1,tt2,tq2,tsgst,tsgsq,twt,twq,tdTdz,tdqdz,arg1,arg2,fr
real(kind=rprec):: tbeta,tnu_t

fr=(1./p_count)*(c_count)

if (jt .EQ. SCAL_init) then
atheta=0._rprec;aq=0._rprec;t2=0._rprec;q2=0._rprec;asgs_t3=0._rprec;
asgs_q3=0._rprec;awt=0._rprec;awq=0._rprec;adTdz=0._rprec;adqdz=0._rprec;
abeta=0._rprec;anu_t=0._rprec
end if

do k=1,nz
do i=1,nx
ttheta1=0._rprec;tq1=0._rprec;tt2=0._rprec;tq2=0._rprec;tsgst=0._rprec;
tsgsq=0._rprec;twt=0._rprec;twq=0._rprec;tdTdz=0._rprec;tdqdz=0._rprec;
tbeta=0._rprec;tnu_t=0._rprec;

do j=1,ny  
    ttheta1=ttheta1+theta(i,j,k)
    tq1=tq1+q(i,j,k)
    tt2=tt2+theta(i,j,k)*theta(i,j,k)
    tq2=tq2+q(i,j,k)*q(i,j,k)
    tsgst=tsgst+sgs_t3(i,j,k)
    tsgsq=tsgsq+sgs_q3(i,j,k)
    tdTdz=tdTdz+dTdz(i,j,k)
    tdqdz=tdqdz+dqdz(i,j,k)
    tbeta=tbeta+beta_scal(i,j,k)
    tnu_t=tnu_t+Nu_t(i,j,k)
if (k.gt.1) then
arg1=(theta(i,j,k)+theta(i,j,k-1))/2.
arg2=(q(i,j,k)+q(i,j,k-1))/2.
else
arg1=0._rprec
arg2=0._rprec
end if
twt=twt+w(i,j,k)*arg1
twq=twq+w(i,j,k)*arg2
end do

atheta(i,k)=atheta(i,k)+(fr)*ttheta1/ny
aq(i,k)=aq(i,k)+(fr)*tq1/ny
t2(i,k)=t2(i,k)+(fr)*tt2/ny
q2(i,k)=q2(i,k)+(fr)*tq2/ny
asgs_t3(i,k)=asgs_t3(i,k)+(fr)*tsgst/ny
asgs_q3(i,k)=asgs_q3(i,k)+(fr)*tsgsq/ny
awt(i,k)=awt(i,k)+(fr)*twt/ny
awq(i,k)=awq(i,k)+(fr)*twq/ny
adTdz(i,k)=adTdz(i,k)+(fr)*tdTdz/ny
adqdz(i,k)=adqdz(i,k)+(fr)*tdqdz/ny
abeta(i,k)=abeta(i,k)+(fr)*tbeta/ny
anu_t(i,k)=anu_t(i,k)+(fr)*tnu_t/ny
end do
end do
      
if (mod(jt,p_count)==0) then
open(36,file=path//'output/aver_theta.out',status="unknown",position="append")
open(37,file=path//'output/aver_q.out',status="unknown",position="append")
open(38,file=path//'output/aver_t2.out',status="unknown",position="append")
open(39,file=path//'output/aver_q2.out',status="unknown",position="append")
open(40,file=path//'output/aver_sgs_t3.out',status="unknown",position="append")
open(41,file=path//'output/aver_sgs_q3.out',status="unknown",position="append")
open(42,file=path//'output/aver_wt.out', status="unknown",position="append")
open(43,file=path//'output/aver_wq.out',status="unknown",position="append")
open(44,file=path//'output/aver_dtdz.out',status="unknown",position="append")
open(45,file=path//'output/aver_dqdz.out',status="unknown",position="append")
open(46,file=path//'output/aver_beta_scalar.out',status="unknown",position="append")
open(47,file=path//'output/aver_Nu_t.out',status="unknown",position="append")

do k=1,nz
write(36,5168) k*dz,(atheta(i,k),i=1,nx)
write(37,5168) k*dz,(aq(i,k),i=1,nx)
write(38,5168) k*dz,(t2(i,k),i=1,nx)
write(39,5168) k*dz,(q2(i,k),i=1,nx)
write(40,5168) k*dz,(asgs_t3(i,k),i=1,nx)
write(41,5168) k*dz,(asgs_q3(i,k),i=1,nx)
write(42,5168) k*dz,(awt(i,k),i=1,nx)
write(43,5168) k*dz,(awq(i,k),i=1,nx)
write(44,5168) k*dz,(adTdz(i,k),i=1,nx)
write(45,5168) k*dz,(adqdz(i,k),i=1,nx)
write(46,5168) k*dz,(abeta(i,k),i=1,nx)
write(47,5168) k*dz,(anu_t(i,k),i=1,nx)
end do

atheta=0._rprec;aq=0._rprec;t2=0._rprec;q2=0._rprec;asgs_t3=0._rprec;
asgs_q3=0._rprec;awt=0._rprec;awq=0._rprec;adTdz=0._rprec;adqdz=0._rprec;
abeta=0._rprec;anu_t=0._rprec

close(36);close(37);close(38);close(39);close(40);close(41);close(42);close(43)
close(44);close(45);close(46);close(47)
end if

if (obukhov_output==1) then
    call obukhov_slice(jt) !output obukhov variables ...
end if


5168     format(1400(E14.5))

end subroutine scalar_slice


!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX

subroutine obukhov_slice(jt)
!subroutine obukhov_slice(phi_m,psi_m,phi_h,psi_h,L,wstar,jt)
use scalars_module,only: L,wstar
use sim_param,only: path
implicit none
integer:: i, k, jt
!real(kind=rprec),dimension(nx,nz):: phi_m,psi_m,phi_h,psi_h
!real(kind=rprec),dimension(nx,nz):: L_fin,wstar_fin
real(kind=rprec),dimension(nx,nz),save:: aphi_m,apsi_m,aphi_h,apsi_h,aL,awstar
real(kind=rprec):: fr


fr=(1./(p_count))*(c_count)
!c ------------------------VK----------------------------------
if (jt .EQ. SCAL_init) then
aphi_m=0._rprec;apsi_m=0._rprec;aphi_h=0._rprec;apsi_h=0._rprec
aL=0._rprec;awstar=0._rprec
end if
!c ------------------------VK----------------------------------

do k=1,ny
do i=1,nx
aphi_m(i,k)=aphi_m(i,k)+(fr)*phi_m(i,k)
apsi_m(i,k)=apsi_m(i,k)+(fr)*psi_m(i,k)
aphi_h(i,k)=aphi_h(i,k)+(fr)*phi_h(i,k)
apsi_h(i,k)=apsi_h(i,k)+(fr)*psi_h(i,k)
aL(i,k)=aL(i,k)+(fr)*L(i,k)
awstar(i,k)=awstar(i,k)+(fr)*wstar(i,k)
end do
end do

if (mod(jt,p_count)==0) then
!cprint *, 'outputting averaged data to files'
open(48,file=path//'output/aver_phi_m.out',&
status="unknown",position="append")
open(49,file=path//'output/aver_psi_m.out',&
status="unknown",position="append")
open(50,file=path//'output/aver_phi_h.out',&
status="unknown",position="append")
open(52,file=path//'output/aver_psi_h.out',&
status="unknown",position="append")
open(53,file=path//'output/aver_L.out',&
status="unknown",position="append")
open(54,file=path//'output/aver_wstar.out',&
status="unknown",position="append")

do k=1,ny
write(48,5168) jt*dz,(aphi_m(i,k),i=1,nx)
write(49,5168) jt*dz,(apsi_m(i,k),i=1,nx)
write(50,5168) jt*dz,(aphi_h(i,k),i=1,nx)
write(52,5168) jt*dz,(apsi_h(i,k),i=1,nx)
write(53,5168) jt*dz,(aL(i,k),i=1,nx)
write(54,5168) jt*dz,(awstar(i,k),i=1,nx)
end do

aphi_m=0._rprec;apsi_m=0._rprec;aphi_h=0._rprec;apsi_h=0._rprec
aL=0._rprec;awstar=0._rprec
end if
!TS
close(48);close(49);close(50);close(52);close(53);close(54)

 5168     format(1400(E14.5))

end subroutine obukhov_slice

!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX

end module scalars_module2
