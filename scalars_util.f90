module scalars_util

implicit none

contains

!************************************************************************
subroutine RHS_scalar_solver (s,dsdx,dsdy,dsdz,RHS,surface_flux,kappa_t)
!************************************************************************
use types, only: rprec
use param, only: nx, ny, nz, ld, ld_big, z_i, dz, vonk, &
                 initu, nx2, ny2, lbz, zo, jt,          &
                 coord, nproc, dz, ml_depth
use sim_param, only: u, v, w
use derivatives, only: ddx, ddy, ddz_w, ddz_uv
use test_filtermodule

implicit none

integer :: i, j, k, n, jz, jy, jx, jz_min

real(rprec), dimension(ld,ny,lbz:nz),intent(in) :: s
real(rprec), dimension(ld,ny,lbz:nz), intent(in) :: dsdx, dsdy, dsdz
real(rprec), dimension(ld,ny), intent(in) :: surface_flux
real(rprec), dimension(ld,ny,lbz:nz) :: RHS
real(rprec), dimension(nz) :: z
integer :: k_ml_global, k_abs
real(rprec), allocatable, dimension(:,:,:) :: u_m,v_m,w_m,dsdx_m,dsdy_m,dsdz_m,RHS_m
real(rprec), allocatable, dimension(:,:,:) :: temp,dtemp
real(rprec), dimension(ld,ny)::u1,v1
real(rprec), dimension(ld,ny,nz), intent(in) :: kappa_t

allocate(u_m(ld_big,ny2,lbz:nz))
allocate(v_m(ld_big,ny2,lbz:nz))
allocate(w_m(ld_big,ny2,lbz:nz))
allocate(dsdx_m(ld_big,ny2,lbz:nz))
allocate(dsdy_m(ld_big,ny2,lbz:nz))
allocate(dsdz_m(ld_big,ny2,lbz:nz))
allocate(RHS_m(ld_big,ny2,lbz:nz))
allocate(temp(ld,ny,lbz:nz))
allocate(dtemp(ld,ny,lbz:nz))

call dealias1(u,u_m,lbz)
call dealias1(v,v_m,lbz)
call dealias1(w,w_m,lbz)
call dealias1(dsdx,dsdx_m,lbz)
call dealias1(dsdy,dsdy_m,lbz)
call dealias1(dsdz,dsdz_m,lbz)

$if($MPI)
   if (coord ==0) then
      jz_min = 2
   else 
      jz_min = 1
   end if
$else
   jz_min = 2
$endif

!Compute RHS of scalar transport equation
!Advection term

do k=jz_min, Nz-1
   do j=1,Ny2
      do i=1,Nx2
         RHS_m(i,j,k) = u_m(i,j,k)*dsdx_m(i,j,k) + v_m(i,j,k)*dsdy_m(i,j,k) &
                      + (w_m(i,j,k)*dsdz_m(i,j,k) + w_m(i,j,k+1)*dsdz_m(i,j,k+1))*0.5_rprec
      end do
   end do
end do

$if($MPI)
    if (coord ==0) then
       do j=1,Ny2
          do i=1,Nx2
             RHS_m(i,j,1)  = u_m(i,j,1)*dsdx_m(i,j,1) + v_m(i,j,1)*dsdy_m(i,j,1) &
                             + (0.5_rprec*w_m(i,j,2))*dsdz_m(i,j,2)
          end do
       end do
    end if
$else
    do j=1,Ny2
       do i=1,Nx2
          RHS_m(i,j,1)  = u_m(i,j,1)*dsdx_m(i,j,1) + v_m(i,j,1)*dsdy_m(i,j,1) &
                          + (0.5_rprec*w_m(i,j,2))*dsdz_m(i,j,2)
       end do
    end do
$endif

$if($MPI)
    if (coord == nproc-1) then
       do j=1,Ny2
          do i=1,Nx2
             RHS_m(i,j,Nz) = u_m(i,j,Nz)*dsdx_m(i,j,Nz) + v_m(i,j,Nz)*dsdy_m(i,j,Nz)
          end do
       end do
    end if
$else
    do j=1,Ny2
       do i=1,Nx2
          RHS_m(i,j,Nz) = u_m(i,j,Nz)*dsdx_m(i,j,Nz) + v_m(i,j,Nz)*dsdy_m(i,j,Nz)
       end do
    end do
$endif

call dealias2(RHS,RHS_m,lbz)

$if($MPI)
   if (coord == 0) then
       do j=1,Ny
          do i=1,Nx
             temp(i,j,1)=kappa_t(i,j,1)*dsdx(i,j,1)
          end do
       end do
   end if
$else
    do j=1,Ny
       do i=1,Nx
          temp(i,j,1)=kappa_t(i,j,1)*dsdx(i,j,1)
       end do
    end do
$endif

do k=jz_min,Nz-1
      do j=1,Ny
         do i=1,Nx
            temp(i,j,k) = 0.5_rprec*(kappa_t(i,j,k)+kappa_t(i,j,k+1))*dsdx(i,j,k)
         end do
      end do
end do

call ddx(temp,dtemp,lbz)

$if($MPI)
   if (coord == 0) then
      do j=1,Ny
         do i=1,Nx
            RHS(i,j,1) = (-1._rprec*RHS(i,j,1)) + dtemp(i,j,1)
            temp(i,j,1) = kappa_t(i,j,1)*dsdy(i,j,1)
         end do
      end do
   end if
$else
    do j=1,Ny
       do i=1,Nx
          RHS(i,j,1) = (-1._rprec*RHS(i,j,1)) + dtemp(i,j,1)
          temp(i,j,1) = kappa_t(i,j,1)*dsdy(i,j,1)
       end do
    end do
$endif


do k=jz_min,Nz-1
   do j=1,Ny
      do i=1,Nx
         RHS(i,j,k) = (-1._rprec*RHS(i,j,k) + dtemp(i,j,k))
      end do 
   end do
end do

do k=jz_min,Nz-1
      do j=1,Ny
         do i=1,Nx
            temp(i,j,k) = 0.5_rprec*(kappa_t(i,j,k)+kappa_t(i,j,k+1))*dsdy(i,j,k)
         end do 
      end do
end do

call ddy(temp,dtemp,lbz)

$if($MPI)
   if (coord == 0) then 
      do j=1,Ny
         do i=1,Nx
            RHS(i,j,1) = RHS(i,j,1) + dtemp(i,j,1)
            ! Get surface flux from similarity relationship
            temp(i,j,1) = -1.0_rprec*surface_flux(i,j)
         end do
      end do
   end if
$else
    do j=1,Ny
       do i=1,Nx
          RHS(i,j,1) = RHS(i,j,1) + dtemp(i,j,1)
          ! Get surface flux from similarity relationship
          temp(i,j,1) = -1.0_rprec*surface_flux(i,j)  
       end do
    end do
$endif

do k=jz_min,Nz-1
   do j=1,Ny
      do i=1,Nx
         RHS(i,j,k) = RHS(i,j,k) + dtemp(i,j,k)
      end do
   end do
end do


do k=jz_min,Nz
      do j=1,Ny
         do i=1,Nx
            temp(i,j,k) = kappa_t(i,j,k)*dsdz(i,j,k)
         end do
      end do
end do

call ddz_w(temp,dtemp,lbz)

do k=1,Nz-1
   do j=1,Ny
      do i=1,Nx
         RHS(i,j,k) = RHS(i,j,k) + dtemp(i,j,k)
      end do
   end do
end do

!if (coord==0) then
!   do jz=1,nz
!      write(*,*) 'jt, jz, kappa_t(25,25,jz)', jt, jz, kappa_t(25,25,jz)
!   end do
!end if

end subroutine RHS_scalar_solver

!************************************************************************
subroutine step_scalar (s, RHS_pre, RHS_post, top_bc)
!************************************************************************
use types, only: rprec
use param, only: nx, ny, nz, ld, l_z, z_i, dz, lbz, dt, &
                 coord, nproc, jt
 
implicit none

integer :: i, j, k, jz

real(rprec), dimension(ld,ny,lbz:nz), intent(inout) :: s 
real(rprec), dimension(ld,ny,lbz:nz), intent(in) :: RHS_pre, RHS_post
real(rprec), intent(in) :: top_bc

do k=1,nz
   do j=1,ny
      do i=1,nx
         s(i,j,k) = s(i,j,k) + &
                         dt*(1.5_rprec*RHS_pre(i,j,k) - 0.5_rprec*RHS_post(i,j,k))
      end do
   end do
end do

$if($MPI)
    if (coord == nproc-1) then
       s(:,:,Nz) = s(:,:,Nz-1) + top_bc*z_i*dz
    end if
$else
      s(:,:,Nz) = s(:,:,Nz-1) + top_bc*z_i*dz
$endif

!if (coord==nproc-1) write(*,*) 'jt, s(25,25,nz-1), s(25,25,nz)', jt, s(25,25,nz-1), s(25,25,nz)

!if (coord==nproc-1) write(*,*) 'top_bc=', top_bc


end subroutine step_scalar


end module
