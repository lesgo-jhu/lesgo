module scalars_module
! HUMIDITY subroutines in place but not yet turned on !!
use types,only:rprec
use param, jt_global => jt  !--rename to avoid name clashes
                            !--could also modify all routines to access jt
                            !  from param module, not argument list
use sim_param,only:u,v,w,theta,q,path
use bottombc,only:zo,num_patch,zot,T_s,q_s,phi_m,psi_m,phi_h,psi_h,&
                  avgpatch
!Includes patches subroutine
use sgsmodule,only:Nu_t
implicit none
!!!!!!--------------------------------------------
! Part I of the scalar files - contains the basic subroutines
! Also look at scalars_module2.f90 for other subroutines !! 
! CONTAINS subroutines :
! theta_all_in_one - Performs the various steps for theta calculation
! humidity_all_in_one - Performs the various steps for humidity
! scalar_RHS_calc - computes the RHS side of the scalar evolution equation
! calcbeta - computes the buoyancy term for temperature
! step_scalar - time evolves the scalar
! obukhov - computes the obukhov similarity terms for use in scalars,wallstress and derivwall
! Authored by Vijayant Kumar
! Last modified - April 24, 2004
!!!!!!--------------------------------------------
!real(kind=rprec),dimension(ld,ny,nz):: scalar
!real(kind=rprec),dimension(ld,ny,nz):: theta,q ! theta and q specified in sim_param
real(kind=rprec),dimension(ld,ny,nz):: beta_scal,Pr_
!real(kind=rprec),dimension(ld,ny,nz):: dsdx,dsdy,dsdz
real(kind=rprec),dimension(ld,ny,nz):: dTdz,dqdz ! Only ones needed for output
! Might need to add x and y derivatives here in case they need to be outputted
! Right now they are in the "scalar"_all_in_one routines below !!
real(kind=rprec),dimension(ld,ny,nz):: RHS_Tf,RHS_T,RHS_qf,RHS_q
real(kind=rprec), dimension(ld,ny,nz):: sgs_t3,sgs_q3 !defines the surface sgs flux

!real(kind=rprec),dimension(:,:,:),allocatable:: scalar,q
!real(kind=rprec),parameter::g=9.81,inv_strength=0.008 !inversion strength (K/m)
real(kind=rprec),dimension(nx,ny)::L,wstar !defines obukhov length and convective vel scale, w_star
real(kind=rprec),dimension(nx,ny)::z_os !T_s and q_s are defined in bottombc.f90
real(kind=rprec),dimension(nx,ny)::ustar_avg ! Defines the local u_star as calculated in obukhov
integer, parameter:: obukhov_output=1 !Controls whether obukhov variables are outputted by scalar_slice

!integer,parameter:: patch_flag=1,remote_flag=0
!integer,dimension(nx,ny)::patch
!integer,dimension(num_patch)::patchnum

!if ((S_FLAG) .and. (jt .eq. SCAL_init)) then
!        if (allocated(scalar)) deallocate(scalar)
!        allocate(scalar(ld,ny,nz))
!        if (allocated(q)) deallocate(q)
!        allocate(q(ld,ny,nz))
!else if((jt .lt. SCAL_init))
!        break
!end if
contains
!----------------------XXXXXXXXXXXVIJXXXXXXXXXXX-----------------------
!----------------------XXXXXXXXXXXVIJXXXXXXXXXXX-----------------------
subroutine theta_all_in_one(jt)
!subroutine theta_all_in_one(theta)
!subroutine scalar_all_in_one(theta)
!subroutine scalar_all_in_one
!use scalars_module
!use scalars_module
use immersedbc,only:building_interp_one
implicit none
real:: wt_s_current
integer :: jt
real(kind=rprec),dimension(ld,ny,nz)::dTdx,dTdy
!real(kind=rprec),dimension(ld,ny,nz)::RHS_T, RHS_Tf
!use scalars_module,theta => scalar

!z_os=zo/10._rprec !Assume - Scalar roughness length = momentum roughness length / 10
wt_s_current=wt_s

Pr_=Pr !Right now set the Prandtl num matrix equal to a constant Prandtl
! number as specified in param. could use the dynamic model ideal to compute Pr as well !!
! The plan-averaged dynamic Prandtl number model is already coded. just need to put it in !!
if(use_bldg)call building_interp_one(theta,.04_rprec,3)
call filt_da(theta,dTdx,dTdy)
call ddz_uv (dTdz,theta)
dTdz(:,:,Nz)=inv_strength/T_scale*z_i ! Valid for temperature
!print *,'DTDZ remote in theta_all_in_one',dTdz(5,5,:)
if (S_FLAG) then
 RHS_Tf=RHS_T
!call scalar_RHS_calc(T_s,z_os,RHS_T,sgs_t3,jt,psi_h,phi_h,Pr_,wt_s_current)
call scalar_RHS_calc(theta,dTdx,dTdy,dTdz,T_s,z_os,RHS_T,sgs_t3,Pr_,wt_s_current)
!TS
if(sflux_flag)then
beta_scal=0._rprec
else
call calcbeta(theta)
endif
! Calculates the buoyancy term which gets added to the vertical momentum equation
!call calcbeta(theta,beta_scal)
if (jt==1. .and. (.not. initu)) then
RHS_Tf=RHS_T
end if

call step_scalar(theta,RHS_T,RHS_Tf)
end if
end subroutine theta_all_in_one


!----------------------XXXXXXXXXXXVIJXXXXXXXXXXX-----------------------
!----------------------XXXXXXXXXXXVIJXXXXXXXXXXX-----------------------
subroutine humidity_all_in_one(jt)
! Not correct now.. DONT USE ....
use sim_param
!use scalars_module
implicit none
real(kind=rprec),dimension(ld,ny,nz)::dqdx,dqdy
real:: wt_s_current
integer:: jt

wt_s_current=wt_s

Pr_=Pr !Right now set the Prandtl no matrix equal to a constant Prandtl
! number as specified in param. could use the dynamic model ideal to compute Pr as well !!
! The plan-averaged dynamic Prandtl number model is already coded. just need to put it in !!

call filt_da(q,dqdx,dqdy)
call ddz_uv (dqdz,q) ! Calculate vertical derivatives !! 
dqdz(:,:,Nz)=0 ! Valid for humidity and other passive scalars

if (S_FLAG) then
 RHS_q=RHS_q
! q_s is the surface humidity - specified using patch_or_remote()
call scalar_RHS_calc(q,dqdx,dqdy,dqdz,q_s,z_os,RHS_q,sgs_q3,Pr_,wt_s_current)

if (jt==1. .and. (.not. initu)) then
RHS_qf=RHS_q
end if

call step_scalar(q,RHS_q,RHS_qf)
end if

end subroutine humidity_all_in_one


!----------------------XXXXXXXXXXXVIJXXXXXXXXXXX-----------------------
!----------------------XXXXXXXXXXXVIJXXXXXXXXXXX-----------------------
subroutine calcbeta (scalar)
!subroutine calcbeta (scalar, beta_scal)
! This calculates the buoyancy term (beta_scal) to be added to the vertical
! momentum equation for temperature
! Authored by Vijayant Kumar
! Last updated April 14, 2004
implicit none
integer::i, j, k
!real(kind=rprec),dimension(ld,ny,nz),intent(out)::beta_scal
real(kind=rprec),dimension(ld,ny,nz),intent(in)::scalar
real(kind=rprec),dimension(nz)::scalar_bar
real(kind=rprec)::g_hat,above, below
!..Non-dimensionalize gravity
g_hat=g*(z_i/(u_star**2))

!..Note Beta is stored on W nodes, but Theta is on UVP nodes
!....We do not time-advance the ground nodes, so start at k=2
! VK: Inserted the averaging code inside this file itself
! rather than doing it in prof
do k=1,nz
scalar_bar(k)=0.0    
   do j=1,ny
      do i=1,nx
        scalar_bar(k)=scalar_bar(k)+scalar(i,j,k)
      end do
   end do
scalar_bar(k)=scalar_bar(k)/(nx*ny)
end do
!....We do not time-advance the ground nodes, so start at k=2
do k=2,Nz
      do j=1,Ny
             do i=1,nx
                above=(scalar(i,j,k)-scalar_bar(k))/scalar_bar(k)
                below=(scalar(i,j,k-1)-scalar_bar(k-1))/scalar_bar(k-1)
                beta_scal(i,j,k)=g_hat*(above + below)/2.
             end do
      end do
end do
return
end subroutine calcbeta
!----------------------XXXXXXXXXXXVIJXXXXXXXXXXX-----------------------
!----------------------XXXXXXXXXXXVIJXXXXXXXXXXX-----------------------


subroutine scalar_RHS_calc(scalar,dsdx,dsdy,dsdz,S_Surf,z_os,RHS,sgs_vert,Pr_,surf_flux_current)
!subroutine scalar_RHS_calc(S_Surf,z_os,RHS,sgs_vert,jt,psi_h,phi_h,Pr_,wt_s_current)
!subroutine scalar_RHS_calc(s,dsdx,dsdy,dsdz,u,v,w,Nu_t,txz,tyz,&
!S_Surf,z_os,patch,patchnum,RHS,sgs_vert,jt,psi_h,phi_h,Pr_,&
!wt_s_current)
!cVK - s, here is the scalar inputted into the subroutine
use immersedbc,only:n_bldg,bldg_pts
implicit none
integer:: i, j, k,ni
!integer:: patch(nx,ny),patchnum(types)
real:: surf_flux_current
!real(kind=rprec),dimension(ld,ny,nz):: u,v,w - !No need as already invoked using sim_param
real(kind=rprec),dimension(ld,ny,nz):: dsdx,dsdy,dsdz
real(kind=rprec),dimension(ld,ny,nz):: RHS,temp
real(kind=rprec),dimension(ld,ny,nz):: scalar
real(kind=rprec),dimension(ld,ny,nz):: dtemp,txz,tyz,sgs_vert,Pr_
!real,dimension(ld,ny,nz):: dtemp,s,txz,tyz,sgs_vert,Pr_,Nu_t
real(kind=rprec),dimension(ld_big,ny2,nz):: u_m,v_m,w_m,dsdx_m,dsdy_m,dsdz_m
 real(kind=rprec),dimension(ld_big,ny2,nz):: RHS_m
!cVK - changed the dimensions for RHS_m,u_m etc. to ld_big
!cVK as otherwise it causes segmentation errors
real(kind=rprec),dimension(nx,ny):: ustar_local,S_Surf,surf_flux,z_os
!real(kind=rprec),dimension(nx,ny):: psi_h,phi_h
!real(kind=rprec),dimension(nx,ny):: psi_h,phi_h,ustar2
!TS CHECK WHAT I CHANGE
real(kind=rprec)::wt_s2
!TS
integer::px,py,lx,ly,lz
!  call ddz_uv (dsdz,scalar) 
!  if ((scalar(2,2,2)<2.and. nproc*L_z>z_i).AND.(wt_s.ne.0.0)) then
!  dsdz(:,:,Nz)=inv_strength/T_scale*z_i ! Valid for temperature
!  else if (wt_s .eq. 0.0) then
!  dsdz(:,:,Nz)=0 ! Valid if heat flux is set as zero - for neutral test with scalars
!  else
!  dsdz(:,:,Nz)=0 ! Valid for humidity and other passive scalars
!  end if

if (patch_flag==1) then
 if (sflux_flag) then !PASSIVE SCALAR IF BLOCK
    surf_flux=0._rprec
     k=num_patch
!Gaussian distribution sigma=1 P=1/(sigma*sqrt(2pi))*exp(-(x-nu)^2/(2sigma^2))
     do j=int(zot(k,4)),int(zot(k,5))
     do i=int(zot(k,2)),int(zot(k,3))
     wt_s2=-(&
      (real(i,kind=rprec)-.5_rprec*(zot(k,2)+zot(k,3)))**2+&
      (real(j,kind=rprec)-.5_rprec*(zot(k,4)+zot(k,5)))**2)/sqrt(2._rprec)
     surf_flux(i,j)&
             =surf_flux_current*exp(wt_s2)/T_scale/u_star
     enddo
     enddo
     ustar_local=ustar_avg !set ustar as value computed in obukhov
  else    ! PASSIVE SCALAR FLUX IF BLOCK CONTINUES 
    call avgpatch(ustar_avg,ustar_local)
  if (lbc==1.and.scalar(1,1,1)<2) then
     do k=1,num_patch
        wt_s2=(-1.)**(k+1)*surf_flux_current 
        surf_flux(int(zot(k,2)):int(zot(k,3)),int(zot(k,4)):int(zot(k,5)))&
             =wt_s2/T_scale/u_star
     enddo
  else if (lbc==0.and.scalar(1,1,1)<2) then
    do j=1,ny
    do i=1,nx 
        surf_flux(i,j)=(S_Surf(i,j)-scalar(i,j,1))*vonk*ustar_local(i,j)&
                /(phi_h(i,j)*log(dz/(2._rprec*z_os(i,j))))
    enddo
    enddo
  endif
 endif
  dsdz(1:nx,1:ny,1) =-phi_h(1:nx,1:ny)*surf_flux(1:nx,1:ny)&
       /(ustar_local(1:nx,1:ny)*vonk*DZ*0.5_rprec)
!TSend if !YUHENG FLUX IF BLOCK ENDS HERE
elseif (remote_flag==1) then
!     ustar_local=(txz(1:nx,1:ny,1)**2+tyz(1:nx,1:ny,1)**2)**.25
     ustar_local=ustar_avg
!ustar_local_mean=SUM(ustar_local(1:16,1:16))/float(nx*ny)
!print *,'ustar_avg',ustar_local_mean
    do j=1,ny
    do i=1,nx
!     if (ustar_local(i,j)==0.) then
!     ustar_local(i,j)=ustar_local_mean
!     print *,'replaced ustar by avg value',ustar_local(i,j)
!     end if
surf_flux(i,j)=(S_Surf(i,j)-scalar(i,j,1))*vonk*ustar_local(i,j)&
/(phi_h(i,j)*log(dz/(2.*z_os(i,j))))
!surf_flux(i,j)=(S_Surf(i,j)-s(i,j,1))*vonk*ustar_local(i,j)
!     &/(log(dz/(2.*z_os(i,j)))-psi_h(i,j))
dsdz(i,j,1) =-phi_h(i,j)*surf_flux(i,j)/(ustar_local(i,j)*vonk*DZ/2.)
    end do
    end do
 end if

!c...Setting the Convective term in the RHS matrix.
!c..........dsdz and w stored on half-nodes...so we average them
!c..........at DZ/2.  dsdz on UVP node, w=(w(DZ)-0)/2.
!c..........at lid w=0 and dsdz=0.
!c Note that dsdz(k=1) is on u,v,p node(dz/2)
!c while dsdz(k>1) is on w nodes.. This is reflected in the computations
!c of RHS for k=1.

!CVK - How is dsdz = 0 at the lid (e.g. temp . we impose
!CVK a gradient of inv_strength from dimen.h
!TTS
if(use_bldg)then
do ni=1,n_bldg
   px=bldg_pts(1,ni)
   py=bldg_pts(2,ni)
   lx=bldg_pts(3,ni)
   ly=bldg_pts(4,ni)
   lz=bldg_pts(5,ni)
   do k=1,lz
   do j=py,py+ly
   do i=px,px+lx
   dsdx(i,j,k)=0._rprec;dsdy(i,j,k)=0._rprec;dsdz(i,j,k)=0._rprec
   enddo
   enddo
   enddo
enddo
endif

 call dealias1(u,u_m)
 call dealias1(v,v_m)
 call dealias1(w,w_m)
 call dealias1(dsdx,dsdx_m)
 call dealias1(dsdy,dsdy_m)
 call dealias1(dsdz,dsdz_m)

! Now compute the RHS term of the filtered scalar equation. 
! Note that this is the advection term with the scalar as 
! the diffusion term has been thrown away. This is done step 
! by step for each of the expanded arrays from dealias1 separately
! for the first node & last node AND the rest of the nodes. 
do k=2,Nz-1
  do j=1,Ny2
    do i=1,Nx2
    RHS_m(i,j,k)=u_m(i,j,k)*dsdx_m(i,j,k)+v_m(i,j,k)*dsdy_m(i,j,k)&
    +(w_m(i,j,k)*dsdz_m(i,j,k)+w_m(i,j,k+1)*dsdz_m(i,j,k+1))*0.5_rprec
     end do
   end do
 end do
 do j=1,Ny2
   do i=1,Nx2
    RHS_m(i,j,1)=u_m(i,j,1)*dsdx_m(i,j,1)+v_m(i,j,1)*dsdy_m(i,j,1)&
    +(0.5_rprec*w_m(i,j,2))*dsdz_m(i,j,2)
    RHS_m(i,j,Nz)=u_m(i,j,Nz)*dsdx_m(i,j,Nz)+v_m(i,j,Nz)*dsdy_m(i,j,Nz)
   end do
 end do
 call dealias2(RHS,RHS_m)
!c...Now building the SGS part of the RHS.
! Here the sgs_term for scalars is built up using Nu_t from sgs_stag_W.f
! and dividing it by the turbulent Prandtl # specified in dimen.h
!c....Note: Since we bring the Convective term to RHS its sign changes.
!c....Below "Temp" is used for SGS flux; its divergence is added to RHS

!VK.. Should we apply dealiasing also to these multiplication of 
!..   arrays as was done for RHS earlier
!TS
do k=1,Nz
do j=1,Ny
do i=1,Nx
   temp(i,j,k)=(1._rprec/Pr_(i,j,k))*Nu_t(i,j,k)*dsdx(i,j,k)
end do
end do
end do
call DDX (dtemp, temp)

do k=1,Nz
do j=1,Ny
do i=1,Nx
   RHS(i,j,k) = (-1._rprec*RHS(i,j,k) + dtemp(i,j,k))
   temp(i,j,k)=(1._rprec/Pr_(i,j,k))*Nu_t(i,j,k)*dsdy(i,j,k)
end do
end do
end do
call DDY (dtemp, temp)   

!c...Use MO flux at wall for the scalar sgs term !
!c Note that the total contribution to the scalar sgs term at
!c the first node comes from the surface flux computed above from
!c the specified heat flux, wt_s

do j=1,Ny
do i=1,Nx
   RHS(i,j,1) = RHS(i,j,1) + dtemp(i,j,1)
   temp(i,j,1) = -1._rprec*surf_flux(i,j)
   sgs_vert(i,j,1) =surf_flux(i,j)
end do
end do
!c...And the values computed from Nu_t (using the sgs_model for momentum)
!c...in interior
!c.....Note Nu_t on UVP, dsdz on W nodes....interpolate Nu_t onto W nodes     
! Note sgs_vert is -1*temp because sgs_vert is modeled as -Nu_t*dsdz/Pr
! while temp is the same term but w/o the minus sign due to the additional
! minus outside the scalar fluctuation flux term in RHS

do k=2,Nz
do j=1,Ny
do i=1,Nx
    RHS(i,j,k) = RHS(i,j,k) + dtemp(i,j,k)
    temp(i,j,k)=(1._rprec/Pr_(i,j,k))*0.5_rprec*&
         (Nu_t(i,j,k)+Nu_t(i,j,k-1))*dsdz(i,j,k)
    sgs_vert(i,j,k)=-1._rprec*temp(i,j,k)
end do
   end do
 end do
!c...The SGS_z flux is on the W nodes, but DDZ_W will put it back on UVP nodes! 
!c Also note that sgs_vert(i,j,k) influences the computations in 
! OBUKHOV.f and is not involved in any computations in this routine.
! sgs_t3(i,j,1) (<w'theta'> is used for computing wt at the surface in OBUKHOV)

  call DDZ_w (dtemp, temp)
  do k=1,Nz
    Do j=1,Ny
    do i=1,Nx
    RHS(i,j,k) = RHS(i,j,k) + dtemp(i,j,k)
    end do
    end do
  end do

!5167  format (110(1x,e14.5)) 

!return 

end subroutine scalar_RHS_calc
!----------------------XXXXXXXXXXXVIJXXXXXXXXXXX-----------------------
!----------------------XXXXXXXXXXXVIJXXXXXXXXXXX-----------------------



subroutine step_scalar(scalar,RHS_pre,RHS_post)
!subroutine step_scalar(scalar,RHS_T,RHS_Tf,wt_s_current)
use immersedbc,only:n_bldg,bldg_pts
implicit none
integer:: i,j,k
real(kind=rprec),dimension(ld,ny,nz)::scalar, RHS_pre, RHS_post
!TS
integer::px,py,lx,ly,lz,ni
!real(kind=rprec)::wt_s_current
!cVK - This routine moves the scalar field (scalar in this case)
!cVK - forward in time using the scalar from previous time step
!cVK - and the RHS terms from the previous two time steps 
!cVK - using second order Adams-Bashforth scheme

do k=1,nz
do j=1,ny
do i=1,nx
   scalar(i,j,k)= scalar(i,j,k)+&
        dt*(1.5_rprec*RHS_pre(i,j,k)-0.5_rprec*RHS_post(i,j,k))
end do
end do
end do
!TS BUILDING CASES
if(use_bldg)then
do k=1,n_bldg
   px=bldg_pts(1,k)
   py=bldg_pts(2,k)
   lx=bldg_pts(3,k)
   ly=bldg_pts(4,k)
   lz=bldg_pts(5,k)
   scalar(px:px+lx,py:py+ly,1:lz)=0._rprec
enddo
endif
!VK The following OR clause to the if statement was added on February 4th, 2003
!VK to account for the case when we put wt_s=0 as this served as a test for the
!VK Coriolis with both momentum and scalars running but decoupled.. as is the case 
!VK for neutral conditions (wt_s=0)

!if ((nproc*L_z .eq. z_i) .OR. (wt_s .eq. 0)) then
!scalar(:,:,nz)=scalar(:,:,nz-1)
!else if ((scalar(2,2,2)<2.and. nproc*L_z>z_i) .AND. (wt_s .ne. 0)) then
!scalar(:,:,nz)=scalar(:,:,nz-1)+0.003/T_scale*z_i*dz !ubc
!end if

!VK Note that this subroutine was designed to be working with a set of scalars (incl.
!VK temperature and humidity and so, these boundary conditions as given below should
!VK be interpreted in the right context and not just for temperature
!VK For example, the first if block refers to a condition with humidity while the
!VK second and third statements are focussed to temperature

if (scalar(2,2,2)<2.and. nproc*L_z>1.) then ! for temperature and non-neutral case
     scalar(:,:,Nz)=scalar(:,:,Nz-1)+inv_strength/T_scale*z_i*dz !ubc 
! inv_strength - refers to the slope of the inversion ~ 0.003K/Km for temperature
else ! for everything else - neutral and passive scalars (may be modified depending on need)
     scalar(:,:,Nz)=scalar(:,:,Nz-1)
end if
!if (scalar(2,2,2)>2.and. nproc*L_z>z_i) then
!     scalar(:,:,Nz)=scalar(:,:,Nz-1)
!else if ((wt_s_current .eq. 0.0).OR.(nproc*L_z.eq.z_i)) then
!     scalar(:,:,Nz)=scalar(:,:,Nz-1)
!else
!     scalar(:,:,Nz)=scalar(:,:,Nz-1)+inv_strength/T_scale*z_i*dz !ubc
!end if

!return
end subroutine step_scalar


!----------------------XXXXXXXXXXXVIJXXXXXXXXXXX-----------------------
!----------------------XXXXXXXXXXXVIJXXXXXXXXXXX-----------------------
subroutine obukhov (jt)  
use types,only:rprec
use test_filtermodule
!use scalars_module
implicit none
integer,intent(in)::jt
integer:: jx,jy
!real(kind=rprec), dimension(ld,ny):: sgs_t3,sgs_q3,u,v,theta
!real(kind=rprec), dimension(nx,ny):: ustar_avg
real(kind=rprec), dimension(ld,ny):: wt_avg,wq_avg,theta_avg,u1,v1
!real(kind=rprec), dimension(nx,ny):: psi_m,psi_h,phi_m,phi_h,zeta,zo
real(kind=rprec), dimension(nx,ny):: x,zeta ! wstar, L already defined above
!real(kind=rprec), dimension(nx,ny):: wstar,L,x,zeta
real(kind=rprec) g_,wt_,wq_,ustar_,theta_,L_,zo_,u_avg(nx,ny),fr
! parameter(g=9.81) ! defined in scalars_module.f90
!double precision, dimension(lh,ny) :: G_test,G_test_test
real(kind=rprec) obuk_L,obuk_ustar,obuk_phi_m,obuk_phi_h,obuk_psi_m,obuk_psi_h,obuk_zo   
!TS ADD sflux_flag
if (jt.LT.SCAL_init.OR.(.NOT.S_FLAG)) then  
    !print *,'obukhov variables = 1 & return before S_FLAG'  
    psi_m=0.  
    phi_m=1.  
    psi_h=0.  
    phi_h=1. 
    L=0
    wstar=0 
!    print *,'return before S_FLAG'  
    return
    print *,'return before S_FLAG'  
end if  

!  nondimensionalize g
!  if (coriolis_forcing) then g_=g/((ug/30.)**2/z_i)!   else
   g_=g/(u_star**2/z_i)  ! end if

theta_avg=theta(:,:,1) 
wt_avg=sgs_t3(:,:,1) ! We need only the surface flux - defined by sgs
wq_avg=sgs_q3(:,:,1) ! We need only the surface flux - defined by sgs
zo_=exp(sum(log(zo(1:nx,1:ny)))/float(nx*ny))
! averages over x-y plane @ z = 1
wt_=sum(sgs_t3(1:nx,1:ny,1))/float(nx*ny)
wq_=sum(sgs_q3(1:nx,1:ny,1))/float(nx*ny)
ustar_=sum(sqrt(u(1:nx,1:ny,1)**2+v(1:nx,1:ny,1)**2))/float(nx*ny)*vonK/&
(log(0.5*dz/zo_)-sum(psi_m(1:nx,1:ny))/float(nx*ny))
theta_=sum(theta_avg(1:nx,1:ny))/float(nx*ny)

!if (remote_flag==0 .and. num_patch==1) then
if (patch_flag==1 .and. num_patch==1) then  
  do jx=1,nx
    do jy=1,ny
wt_avg(jx,jy)=wt_
wq_avg(jx,jy)=wq_
ustar_avg(jx,jy)=ustar_
theta_avg(jx,jy)=theta_
    end do
  end do
else
  u1=u(:,:,1)
  v1=v(:,:,1)
  call test_filter(u1,G_test)
  call test_filter(v1,G_test)
  call test_filter(theta_avg,G_test)
  call test_filter(wt_avg,G_test)
  call test_filter(wq_avg,G_test)
  do jx=1,nx
    do jy=1,ny
     u_avg(jx,jy)=sqrt(u1(jx,jy)**2+v1(jx,jy)**2)
    end do
  end do
  ustar_avg(1:nx,:)=u_avg(:,:)*vonK/(log(0.5_rprec*dz/zo(:,:))-psi_m(:,:))
end if
!TS
if(sflux_flag)return
! Compute Obukhov Length
    
   do jx=1,ny
   do jy=1,nx
    L(jx,jy)=-ustar_avg(jx,jy)**3/(vonk*g_/theta_avg(jx,jy)*wt_avg(jx,jy))
    wstar(jx,jy)=sign((g_/theta_avg(jx,jy)*abs(wt_avg(jx,jy))*z_i)**(1./3.),wt_avg(jx,jy))
!c  for unstable conditions
      if ((L(jx,jy)<0.) .and. (wt_avg(jx,jy) .ne. 0.)) then
             x(jx,jy)=(1.-16.*dz/2./L(jx,jy))**.25
             psi_m(jx,jy)=2.*log((1.+x(jx,jy))/2.)+&
             log((1.+x(jx,jy)**2)/2.)-2.*atan(x(jx,jy))+pi/2.
             psi_h(jx,jy)=2.*log((1.+x(jx,jy)**2)/2.)
             phi_m(jx,jy)=x(jx,jy)**(-1)
             phi_h(jx,jy)=x(jx,jy)**(-2)
      else if ((L(jx,jy)>0.).and.(wt_avg(jx,jy).ne. 0.)) then
! Implementing new formulations for phi and psi for stable case
! using Cheng & Brutsaert (2004): source - Brutsaert's book from
! Marc's Hydrology course
             zeta(jx,jy)=0.5*dz/L(jx,jy)
             phi_m(jx,jy)=1.+6.1*(zeta(jx,jy)+zeta(jx,jy)**2.5*&
             ((1+zeta(jx,jy)**2.5)**(-1.0+1/2.5)))/(zeta(jx,jy)+&
             ((1+zeta(jx,jy)**2.5)**(1/2.5)))
             phi_h(jx,jy)=phi_m(jx,jy) ! Reynolds analogy
             psi_m(jx,jy)=-1.*6.1*log(zeta(jx,jy)+((1+zeta(jx,jy)**2.5)**(1/2.5)))
             psi_h(jx,jy)=psi_m(jx,jy) ! Reynolds analogy
             
!             phi_m(jx,jy)=1.+5.*(zeta(jx,jy))
!             phi_h(jx,jy)=phi_m(jx,jy) ! Reynolds analogy
!             psi_m(jx,jy)=-1.*6.1*log(zeta(jx,jy)+((1+zeta(jx,jy)**2.5)**(1/2.5)))
!             psi_h(jx,jy)=psi_m(jx,jy) ! Reynolds analogy
      else
             psi_m(jx,jy)=0.
             psi_h(jx,jy)=0.
             phi_m(jx,jy)=1.
             phi_h(jx,jy)=1.
      end if ! (Loop5 ends)

  end do
  end do

  L_=-ustar_**3/(vonk*g_/theta_*wt_)
  

!-------------------- OUTPUT ------------------------------
if (mod(jt,c_count)==0) then
  fr=float(c_count)/float(p_count)
  obuk_L=obuk_L+fr*L_
  obuk_ustar=obuk_ustar+fr*ustar_
  obuk_phi_m=obuk_phi_m+fr*sum(phi_m(:,:))/float(nx*ny)
  obuk_psi_m=obuk_psi_m+fr*sum(psi_m(:,:))/float(nx*ny)
  obuk_phi_h=obuk_phi_h+fr*sum(phi_h(:,:))/float(nx*ny)
  obuk_psi_h=obuk_psi_h+fr*sum(psi_h(:,:))/float(nx*ny)
  obuk_zo=obuk_zo+fr*zo_
end if

if (mod(jt,p_count)==0) then
    open (unit=47,file=path//'output/mo.out',status="unknown",position="append")
  write(47,5168) jt*dt,obuk_L,obuk_ustar,obuk_phi_m,obuk_psi_m,obuk_phi_h,obuk_psi_h,obuk_zo
  close(47)

    obuk_L=0.
    obuk_ustar=0.
    obuk_phi_m=0.
    obuk_phi_h=0.
    obuk_psi_m=0.
    obuk_psi_h=0.
    obuk_zo=0.
end if

 5168     format(1400(E14.5))
  return
end subroutine obukhov 

!----------------------XXXXXXXXXXXVIJXXXXXXXXXXX-----------------------
!----------------------XXXXXXXXXXXVIJXXXXXXXXXXX-----------------------

end module scalars_module


