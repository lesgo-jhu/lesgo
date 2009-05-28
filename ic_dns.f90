subroutine ic_dns()
use types,only:rprec
use param
use param2
use sim_param,only:u,v,w
implicit none
real(kind=rprec),dimension(nz)::ubar
real(kind=rprec)::rms,temp,ran3
integer::jx,jy,jz,seed,z

if (inflow) then

  ! uniform flow case:
  ubar = face_avg

else

  seed=-112

  ! calculate height of first uvp point in wall units
  ! lets do a laminar case (?)
  do jz=1,nz
  
     z=(real(jz)-.5_rprec)*dz ! non-dimensional
     ubar(jz)=(u_star*z_i/nu_molec)*z*(1._rprec-.5_rprec*z) ! non-dimensional
  !         ubar(jz)=0.
  end do  
end if

do jz=1,nz
  print *,'jz, ubar:',jz,ubar(jz)
end do
! rms=0.0001 seems to work in some cases
! the "default" rms of a unif variable is 0.289
rms=0.2_rprec
do jz=1,nz
  do jy=1,ny
     do jx=1,nx
       u(jx,jy,jz)=ubar(jz)+(rms/.289_rprec)*(ran3(seed)-.5_rprec)/u_star
       v(jx,jy,jz)=0._rprec+(rms/.289_rprec)*(ran3(seed)-.5_rprec)/u_star
       w(jx,jy,jz)=0._rprec+(rms/.289_rprec)*(ran3(seed)-.5_rprec)/u_star
    end do
  end do
end do

! make sure w-mean is 0
temp=0._rprec
do jz=1,nz
   do jy=1,ny
      do jx=1,nx
         temp=temp+w(jx,jy,jz)
      end do
   end do
end do
temp=temp/(nx*ny*nz)

do jz=1,nz
   do jy=1,ny
      do jx=1,nx
         w(jx,jy,jz)=w(jx,jy,jz)-temp
      end do
   end do
end do
      
w(:,:,1)=0._rprec
w(:,:,nz)=0._rprec
u(:,:,nz)=u(:,:,nz-1)
v(:,:,nz)=v(:,:,nz-1)
end subroutine ic_dns
