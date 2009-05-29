module immersedbc
use types,only:rprec
use param2,only:ld,ny,nz
implicit none
private
public n_bldg,bldg_pts,fx,fy,fz,u_des,v_des,w_des&
     ,building_mask,building_interp,building_mask_one,building_interp_one&
     ,wallstress_building,walldudx_building
integer::n_bldg
integer,allocatable::bldg_pts(:,:)
!real(kind=rprec), dimension(ld,ny,nz)::fx,fy,fz,u_des,v_des,w_des
real(kind=rprec), allocatable, dimension(:,:,:)::fx,fy,fz,u_des,v_des,w_des
logical,parameter::callag=.false.
!--is this for the building or what?
real(kind=rprec),parameter::zo_avg=2._rprec

contains

subroutine building_mask(u,v,w)
!TSuse sim_param,only:u,v,w
implicit none
integer::px,py,lx,ly,lz
integer::i
!real(kind=rprec),dimension(ld,ny,nz),intent(inout)::u,v,w

! this sets pressure grad's inside bldg's
do i=1,n_bldg
   px=bldg_pts(1,i)
   py=bldg_pts(2,i)
   lx=bldg_pts(3,i)
   ly=bldg_pts(4,i)
   lz=bldg_pts(5,i)
   u(px:px+lx,py:py+ly,1:lz)=0._rprec
   v(px:px+lx,py:py+ly,1:lz)=0._rprec
   w(px:px+lx,py:py+ly,1:lz)=0._rprec
enddo
end subroutine building_mask

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--this appears to apply SOR (or underrelaxation) to 2d laplace 
!  equation (horz. planes), presumably to smooth the solution within
!  the buildings
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine building_interp(u,v,w,weight,iternum)
!TSuse sim_param,only:u,v,w
implicit none
integer::px,py,lx,ly,lz
integer::jx,jy,jz,i,j
integer,intent(in)::iternum
real(kind=rprec),dimension(ld,ny,nz),intent(inout)::u,v,w
real(kind=rprec),intent(in)::weight
real(kind=rprec)::weightp
weightp=(1-weight)*.25_rprec
do i=1,n_bldg
   px=bldg_pts(1,i)
   py=bldg_pts(2,i)
   lx=bldg_pts(3,i)
   ly=bldg_pts(4,i)
   lz=bldg_pts(5,i)
   do j=1,lx*ly*iternum
      do jz=1,lz
      do jy=py,py+ly
      do jx=px,px+lx
         u(jx,jy,jz)=u(jx,jy,jz)*weight+&
              u(jx-1,jy,jz)*weightp+u(jx+1,jy,jz)*weightp+&
              u(jx,jy-1,jz)*weightp+u(jx,jy+1,jz)*weightp
         v(jx,jy,jz)=v(jx,jy,jz)*weight+&
              v(jx-1,jy,jz)*weightp+v(jx+1,jy,jz)*weightp+&
              v(jx,jy-1,jz)*weightp+v(jx,jy+1,jz)*weightp
         w(jx,jy,jz)=w(jx,jy,jz)*weight+&
              w(jx-1,jy,jz)*weightp+w(jx+1,jy,jz)*weightp+&
              w(jx,jy-1,jz)*weightp+w(jx,jy+1,jz)*weightp
      enddo
      enddo
      enddo
   enddo
end do
end subroutine building_interp

subroutine building_interp_one(u,weight,iternum)
!TSuse sim_param,only:u,v,w
implicit none
integer::px,py,lx,ly,lz
integer::jx,jy,jz,i,j
integer,intent(in)::iternum
real(kind=rprec),dimension(ld,ny,nz),intent(inout)::u
real(kind=rprec),intent(in)::weight
real(kind=rprec)::weightp

weightp=(1-weight)*.25_rprec

do i=1,n_bldg
   px=bldg_pts(1,i)
   py=bldg_pts(2,i)
   lx=bldg_pts(3,i)
   ly=bldg_pts(4,i)
   lz=bldg_pts(5,i)
   do j=1,lx*ly*iternum
      do jz=1,lz
      do jy=py,py+ly
      do jx=px,px+lx
         u(jx,jy,jz)=u(jx,jy,jz)*weight+&
              u(jx-1,jy,jz)*weightp+u(jx+1,jy,jz)*weightp+&
              u(jx,jy-1,jz)*weightp+u(jx,jy+1,jz)*weightp
      enddo
      enddo
      enddo
   enddo
end do
end subroutine building_interp_one

subroutine building_mask_one(stress,value)
implicit none
integer::px,py,lx,ly,lz
integer::i
real(kind=rprec)::value
real(kind=rprec),dimension(ld,ny,nz),intent(inout)::stress
! this sets pressure grad's inside bldg's
do i=1,n_bldg
   px=bldg_pts(1,i)
   py=bldg_pts(2,i)
   lx=bldg_pts(3,i)
   ly=bldg_pts(4,i)
   lz=bldg_pts(5,i)
   stress(px:px+lx,py:py+ly,1:lz)=value
end do
end subroutine building_mask_one

subroutine walldudx_building
use types,only:rprec
use param,only:vonk
use param2,only:dy,dx,dz,ld,lh,nx,ny,nz,z_i
use sim_param,only:u,v,w,dudx,dvdx,dwdx,dudy,dvdy,dwdy,dudz,dvdz,dwdz
implicit none
integer::jx,jy,jz,i,j
integer::px,py,lx,ly,lz
real(kind=rprec)::u1,v1,ustar,u_avg,u_aver
real(kind=rprec)::const,dz1,zo,zo_s
zo=zo_avg/z_i
zo_s=zo_avg/z_i
do i=1,n_bldg
   px=bldg_pts(1,i)
   py=bldg_pts(2,i)
   lx=bldg_pts(3,i)
   ly=bldg_pts(4,i)
   lz=bldg_pts(5,i)
   dz1=.5_rprec*dz
   jz=lz+1
   do jy=py,py+ly
   do jx=px,px+lx
    u1=u(jx,jy,jz)*.5_rprec
    v1=v(jx,jy,jz)*.5_rprec
    u_avg=sqrt(u1**2+v1**2)
    ustar=u_avg*vonk/log(dz1/zo)
    const=-(ustar**2)/u_avg
    !this is as in Moeng 84
    dudz(jx,jy,jz)=ustar/(dz1*vonK)*u1/u_avg
    dvdz(jx,jy,jz)=ustar/(dz1*vonK)*v1/u_avg
   end do
   end do
   jy=py+ly+1
   dz1=dy
   do jz=2,lz
   do jx=px,px+lx
    u1=u(jx,jy,jz)
    v1=(w(jx,jy,jz)+w(jx,jy,jz+1))*.5_rprec
    u_avg=sqrt(u1**2+v1**2)
    ustar=u_avg*vonk/log(dz1/zo_s)
    const=-(ustar**2)/u_avg
    !this is as in Moeng 84
    dudy(jx,jy,jz)=ustar/(dz1*vonK)*u1/u_avg
    dwdy(jx,jy,jz)=ustar/(dz1*vonK)*v1/u_avg
   end do
   end do
   dwdy(px:px+lx,jy,1)=0._rprec
   do jz=lz,2,-1
    dwdy(px:px+lx,jy,jz)=(dwdy(px:px+lx,jy,jz)+dwdy(px:px+lx,jy,jz-1))*.5_rprec
   enddo
   jy=py-1
   do jz=2,lz
   do jx=px,px+lx
    u1=u(jx,jy,jz)
    v1=(w(jx,jy,jz)+w(jx,jy,jz+1))*.5_rprec
    u_avg=sqrt(u1**2+v1**2)
    ustar=u_avg*vonk/log(dz1/zo_s)
    const=-(ustar**2)/u_avg
    !this is as in Moeng 84
    dudy(jx,jy,jz)=-ustar/(dz1*vonK)*u1/u_avg
    dwdy(jx,jy,jz)=-ustar/(dz1*vonK)*v1/u_avg
   end do
   end do
   dwdy(px:px+lx,jy,1)=0._rprec
   do jz=lz,2,-1
    dwdy(px:px+lx,jy,jz)=(dwdy(px:px+lx,jy,jz)+dwdy(px:px+lx,jy,jz-1))*.5_rprec
   enddo
   dz1=dx
   jx=px+lx+1
   do jz=2,lz
   do jy=py,py+ly
    u1=v(jx,jy,jz)
    v1=(w(jx,jy,jz)+w(jx,jy,jz+1))*.5_rprec
    u_avg=sqrt(u1**2+v1**2)
    ustar=u_avg*vonk/log(dz1/zo_s)
    const=-(ustar**2)/u_avg
!    !this is as in Moeng 84
    dvdx(jx,jy,jz)=-ustar/(dz1*vonK)*u1/u_avg
    dwdx(jx,jy,jz)=-ustar/(dz1*vonK)*v1/u_avg
   end do
   end do
   dwdx(jx,py:py+ly,1)=0._rprec
   do jz=lz,2,-1
    dwdx(jx,py:py+ly,jz)=(dwdx(jx,py:py+ly,jz)+dwdx(jx,py:py+ly,jz-1))*.5_rprec
   enddo
   jx=px-1
   do jz=2,lz
   do jy=py,py+ly
    u1=v(jx,jy,jz)
    v1=(w(jx,jy,jz)+w(jx,jy,jz+1))*.5_rprec
    u_avg=sqrt(u1**2+v1**2)
    ustar=u_avg*vonk/log(dz1/zo_s)
    const=-(ustar**2)/u_avg
    !this is as in Moeng 84
    dvdx(jx,jy,jz)=ustar/(dz1*vonK)*u1/u_avg
    dwdx(jx,jy,jz)=ustar/(dz1*vonK)*v1/u_avg
   end do
   end do
   dwdx(jx,py:py+ly,1)=0._rprec
   do jz=lz,2,-1
    dwdx(jx,py:py+ly,jz)=(dwdx(jx,py:py+ly,jz)+dwdx(jx,py:py+ly,jz-1))*.5_rprec
   enddo
   do j=1,lx*ly*3
      do jz=1,lz
      do jy=py,py+ly
      do jx=px,px+lx
   dvdx(jx,jy,jz)=dvdx(jx,jy,jz)*.2_rprec+&
        dvdx(jx-1,jy,jz)*.2_rprec+dvdx(jx+1,jy,jz)*.2_rprec+&
        dvdx(jx,jy-1,jz)*.2_rprec+dvdx(jx,jy+1,jz)*.2_rprec
   dwdx(jx,jy,jz)=dwdx(jx,jy,jz)*.2_rprec+&
        dwdx(jx-1,jy,jz)*.2_rprec+dwdx(jx+1,jy,jz)*.2_rprec+&
        dwdx(jx,jy-1,jz)*.2_rprec+dwdx(jx,jy+1,jz)*.2_rprec
   dudy(jx,jy,jz)=dudy(jx,jy,jz)*.2_rprec+&
        dudy(jx-1,jy,jz)*.2_rprec+dudy(jx+1,jy,jz)*.2_rprec+&
        dudy(jx,jy-1,jz)*.2_rprec+dudy(jx,jy+1,jz)*.2_rprec
   dwdy(jx,jy,jz)=dwdy(jx,jy,jz)*.2_rprec+&
        dwdy(jx-1,jy,jz)*.2_rprec+dwdy(jx+1,jy,jz)*.2_rprec+&
        dwdy(jx,jy-1,jz)*.2_rprec+dwdy(jx,jy+1,jz)*.2_rprec
   dudz(jx,jy,jz)=dudz(jx,jy,jz)*.2_rprec+&
        dudz(jx-1,jy,jz)*.2_rprec+dudz(jx+1,jy,jz)*.2_rprec+&
        dudz(jx,jy-1,jz)*.2_rprec+dudz(jx,jy+1,jz)*.2_rprec
   dvdz(jx,jy,jz)=dvdz(jx,jy,jz)*.2_rprec+&
        dvdz(jx-1,jy,jz)*.2_rprec+dvdz(jx+1,jy,jz)*.2_rprec+&
        dvdz(jx,jy-1,jz)*.2_rprec+dvdz(jx,jy+1,jz)*.2_rprec
   enddo
   enddo
   enddo
   enddo
enddo
end subroutine walldudx_building

subroutine wallstress_building(txy,txz,tyz)
use types,only:rprec
use param,only:vonk
use param2,only:dy,dx,dz,ld,lh,nx,ny,nz,z_i
use sim_param,only:u,v,w
implicit none
integer::jx,jy,jz,i,j
integer::px,py,lx,ly,lz
real(kind=rprec),dimension(ld,ny,nz),intent(out)::txz,tyz,txy
real(kind=rprec)::u1,v1,ustar,u_avg,u_aver
real(kind=rprec)::const,dz1,zo,zo_s
zo=zo_avg/z_i
zo_s=zo_avg/z_i
do i=1,n_bldg
   px=bldg_pts(1,i)
   py=bldg_pts(2,i)
   lx=bldg_pts(3,i)
   ly=bldg_pts(4,i)
   lz=bldg_pts(5,i)
   dz1=.5_rprec*dz
   jz=lz+1
   do jy=py,py+ly
   do jx=px,px+lx
    u1=u(jx,jy,jz)*.5_rprec
    v1=v(jx,jy,jz)*.5_rprec
    u_avg=sqrt(u1**2+v1**2)
    ustar=u_avg*vonk/log(dz1/zo)
    const=-(ustar**2)/u_avg
    txz(jx,jy,jz)=const*u1
    tyz(jx,jy,jz)=const*v1
   end do
   end do
   jy=py+ly+1
   dz1=dy
   do jz=2,lz
   do jx=px,px+lx
    u1=u(jx,jy,jz)
    v1=(w(jx,jy,jz)+w(jx,jy,jz+1))*.5_rprec
    u_avg=sqrt(u1**2+v1**2)
    ustar=u_avg*vonk/log(dz1/zo_s)
    const=-(ustar**2)/u_avg
    txy(jx,jy,jz)=const*u1
    tyz(jx,jy,jz)=const*v1
   end do
   end do
   do jz=lz,2,-1
    tyz(px:px+lx,jy,jz)=(tyz(px:px+lx,jy,jz)+tyz(px:px+lx,jy,jz-1))*.5_rprec
   enddo
   jy=py-1
   do jz=2,lz
   do jx=px,px+lx
    u1=u(jx,jy,jz)
    v1=(w(jx,jy,jz)+w(jx,jy,jz+1))*.5_rprec
    u_avg=sqrt(u1**2+v1**2)
    ustar=u_avg*vonk/log(dz1/zo_s)
    const=-(ustar**2)/u_avg
    txy(jx,jy,jz)=-const*u1
    tyz(jx,jy,jz)=-const*v1
   end do
   end do
   do jz=lz,2,-1
    tyz(px:px+lx,jy,jz)=(tyz(px:px+lx,jy,jz)+tyz(px:px+lx,jy,jz-1))*.5_rprec
   enddo
   dz1=dx
   jx=px+lx+1
   do jz=2,lz
   do jy=py,py+ly
    u1=v(jx,jy,jz)
    v1=(w(jx,jy,jz)+w(jx,jy,jz+1))*.5_rprec
    u_avg=sqrt(u1**2+v1**2)
    ustar=u_avg*vonk/log(dz1/zo_s)
    const=-(ustar**2)/u_avg
    txy(jx,jy,jz)=-const*u1
    txz(jx,jy,jz)=-const*v1
   end do
   end do
   do jz=lz,2,-1
    txz(jx,py:py+ly,jz)=(txz(jx,py:py+ly,jz)+txz(jx,py:py+ly,jz-1))*.5_rprec
   enddo
   jx=px-1
   do jz=2,lz
   do jy=py,py+ly
    u1=v(jx,jy,jz)
    v1=(w(jx,jy,jz)+w(jx,jy,jz+1))*.5_rprec
    u_avg=sqrt(u1**2+v1**2)
    ustar=u_avg*vonk/log(dz1/zo_s)
    const=-(ustar**2)/u_avg
    txy(jx,jy,jz)=const*u1
    txz(jx,jy,jz)=const*v1
   end do
   end do
   do jz=lz,2,-1
    txz(jx,py:py+ly,jz)=(txz(jx,py:py+ly,jz)+txz(jx,py:py+ly,jz-1))*.5_rprec
   enddo
   do j=1,lx*ly*2
      do jz=1,lz
      do jy=py,py+ly
      do jx=px,px+lx
   txy(jx,jy,jz)=txy(jx,jy,jz)*.04_rprec+&
        txy(jx-1,jy,jz)*.24_rprec+txy(jx+1,jy,jz)*.24_rprec+&
        txy(jx,jy-1,jz)*.24_rprec+txy(jx,jy+1,jz)*.24_rprec
   txz(jx,jy,jz)=txz(jx,jy,jz)*.04_rprec+&
        txz(jx-1,jy,jz)*.24_rprec+txz(jx+1,jy,jz)*.24_rprec+&
        txz(jx,jy-1,jz)*.24_rprec+txz(jx,jy+1,jz)*.24_rprec
   tyz(jx,jy,jz)=tyz(jx,jy,jz)*.04_rprec+&
        tyz(jx-1,jy,jz)*.24_rprec+tyz(jx+1,jy,jz)*.24_rprec+&
        tyz(jx,jy-1,jz)*.24_rprec+tyz(jx,jy+1,jz)*.24_rprec
   enddo
   enddo
   enddo
   enddo
enddo
end subroutine wallstress_building

!TS SUBROUTINES FOR CUBIC SPLINE INTERPOLATION
SUBROUTINE bcuint(y,y1,y2,y12,x1dif,x2dif,M,N,x1,x2,ansy,jz)
implicit none
INTEGER::M,N,jz
REAL(KIND=rprec),INTENT(IN)::x1dif,x2dif
REAL(KIND=rprec),DIMENSION(M),INTENT(IN)::x1
REAL(KIND=rprec),DIMENSION(N),INTENT(IN)::x2
REAL(KIND=rprec),DIMENSION(M,N),INTENT(INOUT)::ansy
REAL(KIND=rprec),INTENT(IN)::y(4),y1(4),y12(4),y2(4)
INTEGER::i,j,k
REAL(KIND=rprec)::t,u
REAL(KIND=rprec),DIMENSION(4,4)::c
call bcucof(y,y1,y2,y12,x1dif,x2dif,c)
!if(jz.eq.16)then
!   print*,'y1',real(y),x1dif,x2dif
!   print*,'x12',x1,x2
!endif
if(x1dif.eq.0._rprec.or.x2dif.eq.0._rprec)stop 'bad input in bcuint'
do j=2,N-1
do i=2,M-1
t=x1(i)/x1dif
u=x2(j)/x2dif
ansy(i,j)=0._rprec
do k=4,1,-1
   ansy(i,j)=t*ansy(i,j)+((c(k,4)*u+c(k,3))*u+c(k,2))*u+c(k,1)
enddo
enddo
enddo
END SUBROUTINE bcuint

SUBROUTINE bcucof(y,y1,y2,y12,d1,d2,c)
IMPLICIT NONE
REAL(KIND=rprec), INTENT(IN) :: d1,d2
REAL(KIND=rprec), DIMENSION(4), INTENT(IN) :: y,y1,y2,y12
REAL(KIND=rprec), DIMENSION(4,4), INTENT(OUT) :: c
REAL(KIND=rprec), DIMENSION(16) :: x
REAL(KIND=rprec), DIMENSION(16,16) :: wt
DATA wt /1,0,-3,2,4*0,-3,0,9,-6,2,0,-6,4,&
     8*0,3,0,-9,6,-2,0,6,-4,10*0,9,-6,2*0,-6,4,2*0,3,-2,6*0,-9,6,&
     2*0,6,-4,4*0,1,0,-3,2,-2,0,6,-4,1,0,-3,2,8*0,-1,0,3,-2,1,0,-3,&
     2,10*0,-3,2,2*0,3,-2,6*0,3,-2,2*0,-6,4,2*0,3,-2,0,1,-2,1,5*0,&
     -3,6,-3,0,2,-4,2,9*0,3,-6,3,0,-2,4,-2,10*0,-3,3,2*0,2,-2,2*0,&
     -1,1,6*0,3,-3,2*0,-2,2,5*0,1,-2,1,0,-2,4,-2,0,1,-2,1,9*0,-1,2,&
     -1,0,1,-2,1,10*0,1,-1,2*0,-1,1,6*0,-1,1,2*0,2,-2,2*0,-1,1/
x(1:4)=y
x(5:8)=y1*d1
x(9:12)=y2*d2
x(13:16)=y12*d1*d2
x=matmul(wt,x)
c=reshape(x,(/4,4/),order=(/2,1/))
END SUBROUTINE bcucof
end module immersedbc
