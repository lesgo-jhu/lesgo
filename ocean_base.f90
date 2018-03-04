module ocean_base

use types, only: rprec
$if ($MPI)
use mpi_defs
$endif

implicit none

!Ocean parameters
logical, parameter :: OCEAN_FLAG=.true.
logical, parameter :: sponge_damping =.true.

contains


!**********************************************************************
subroutine inst_mixed_layer_depth(jt, theta)              
!**********************************************************************

use types, only: rprec
use param, only: ld, ny, lbz, nz, jzmax, nx, ny, coord, dz, z_i, nproc, &
                 ierr, comm, MPI_RPREC, MPI_INTEGER, down, up, status
use sim_param, only: w, wpthetap_zplane, thetap2_zplane
use functions, only: interp_to_uv_grid
use mpi


implicit none 

integer :: i,j,k, jz_min, k_ml
integer :: k_ml_global 
integer, intent(in) :: jt
real(rprec) :: favg, z, wpthetap_min 
real(rprec) ::  wpthetap_min_global
real(rprec), dimension(ld,ny,lbz:nz), intent(in) :: theta
real(rprec), dimension(ld,ny,lbz:nz) :: w_uv
real(rprec), dimension(ld,ny,lbz:nz) :: wtheta, theta2
real(rprec), dimension(lbz:nz) ::  w_zplane, theta_zplane, &
                                 wtheta_zplane, theta2_zplane

!Total number of grid points
favg = real(nx*ny,kind=rprec)

wpthetap_min = 0._rprec
wpthetap_min_global = 0._rprec
k_ml =123456789._rprec

wtheta(1:nx,1:ny,lbz:nz) = 0._rprec
theta2(1:nx,1:ny,lbz:nz) = 0._rprec

!Interpolate theta to w-grid to calculate heat flux
w_uv(1:nx,1:ny,lbz:nz) = interp_to_uv_grid( w(1:nx,1:ny,lbz:nz), lbz ) 

!Calculate the horizontal average of the heat flux wT
do k=lbz,jzmax
   w_zplane(k) = 0._rprec
   theta_zplane(k) = 0._rprec
   wtheta_zplane(k) = 0._rprec
   theta2_zplane(k) = 0._rprec
   do j=1,ny
      do i=1,nx
         w_zplane(k) = w_zplane(k) + w_uv(i,j,k)
         theta_zplane(k) = theta_zplane(k) + theta(i,j,k)
         wtheta(i,j,k) = w_uv(i,j,k) * theta(i,j,k)
         wtheta_zplane(k) = wtheta_zplane(k) + wtheta(i,j,k)
         theta2(i,j,k) = theta(i,j,k)*theta(i,j,k)
         theta2_zplane(k) = theta2_zplane(k) + theta2(i,j,k)
      enddo
   enddo
   w_zplane(k) = w_zplane(k)/favg
   theta_zplane(k) = theta_zplane(k)/favg
   wtheta_zplane(k) = wtheta_zplane(k)/favg
   theta2_zplane(k) = theta2_zplane(k)/favg
enddo

!Calculate the horizontal average of the fluctuating part of the heat flux w'T; 
do k=lbz,jzmax
   wpthetap_zplane(k) = wtheta_zplane(k) - w_zplane(k)*theta_zplane(k)
   thetap2_zplane(k) = theta2_zplane(k) - theta_zplane(k)*theta_zplane(k)
enddo


$if($MPI) 
   call mpi_sendrecv (wpthetap_zplane(1), 1, MPI_RPREC, down, 1, &
                      wpthetap_zplane(nz), 1, MPI_RPREC, up, 1, &
                      comm, status, ierr)
   call mpi_sendrecv (wpthetap_zplane(nz-1), 1, MPI_RPREC, up, 2, &
                      wpthetap_zplane(0), 1, MPI_RPREC, down, 2, &
                     comm, status, ierr)
$endif


$if($MPI) 
   call mpi_sendrecv (thetap2_zplane(1), 1, MPI_RPREC, down, 3, &
                      thetap2_zplane(nz), 1, MPI_RPREC, up, 3, &
                      comm, status, ierr)
   call mpi_sendrecv (thetap2_zplane(nz-1), 1, MPI_RPREC, up, 4, &
                      thetap2_zplane(0), 1, MPI_RPREC, down, 4, &
                     comm, status, ierr)
$endif

$if($MPI)
   if (coord ==0) then
      jz_min = 2
   else
      jz_min = 1
   end if
$else
   jz_min = 2
$endif

!Find the minimum of horizontally-averaged w'T' on each processor
do k=jz_min, jzmax
   if (wpthetap_zplane(k).lt.wpthetap_min) then
      wpthetap_min = wpthetap_zplane(k)
   endif
enddo

!Find the global minimum of horizontally-averaged w'T'
call mpi_allreduce (wpthetap_min, wpthetap_min_global, 1, MPI_RPREC, MPI_MIN, comm, ierr)

!Find the z-level of the  global minimum of the horizontally-averaged w'T'
do k=jz_min, jzmax
   if (wpthetap_zplane(k) == wpthetap_min_global) then
      k_ml = (coord) *(nz-1) + k 
      !print*, 'A. jt, coord, k, k_ml, wpthetap_min_global', jt, coord, k, k_ml, wpthetap_min_global 
   endif
enddo


!If the global minimum z-level is found on more than one processor choose the processor
!with the lowest minimum z-level
call mpi_allreduce (k_ml, k_ml_global, 1, MPI_INTEGER, MPI_MIN, comm, ierr)

!if ( (coord ==0 ).and. (mod(jt,100).eq.0) ) then 
!   print*, 'jt, k_ml_global, wpthetap_min_global', jt, k_ml_global, wpthetap_min_global 
!end if 

end subroutine inst_mixed_layer_depth


!**********************************************************************
subroutine thermocline_depth(k_ml_global)              
!**********************************************************************

use param, only : coord, nproc, nz, z_i, dz, ml_depth, lbz
use types, only : rprec

implicit none

integer, intent(out) :: k_ml_global
real(rprec) ::  z
integer :: k, k_abs, jz_min

$if($MPI)
   if (coord ==0) then
      jz_min = 2
   else 
      jz_min = 1
   end if
$else
   jz_min = 2
$endif


!Find the vertical level of the end of the mixed layer
!$if($MPI)
do k=jz_min, (nproc*(nz-1))+1
   k_abs = coord*(nz-1) + k 
   z = (k - 0.5)*dz*z_i
   if ((z.gt.ml_depth).and.((z-(dz*z_i)).le.ml_depth)) then
      k_ml_global = k-1
      exit
   end if
end do
!$endif

end subroutine thermocline_depth


!**********************************************************************
subroutine prandtl_profile(kappa_tt)              
!**********************************************************************

use param, only : coord, nproc, ld, nz, ny, Pr_lam, Pr_turb, z_i, &
                  dz, jt
use types, only : rprec
use sgs_param, only: Nu_t

implicit none

integer :: k_ml_global
real(rprec), dimension(nz) :: Pr
real(rprec), dimension(ld,ny,nz), intent(out)::kappa_tt
integer :: k, k_abs

call thermocline_depth(k_ml_global) 

do k=1, nz
   k_abs = coord*(nz-1) + k
   if (k_ml_global.ge.k_abs) then
      Pr(k) = Pr_turb
   else
      Pr(k) = Pr_lam
   end if
      kappa_tt(:,:,k) = Nu_t(:,:,k)/Pr(k)
end do 

end subroutine prandtl_profile


!**********************************************************************
subroutine schmidt_profile(kappa_ts)              
!**********************************************************************

use param, only : coord, nproc, ld, nz, ny, Sc_lam, Sc_turb, z_i, &
                  dz, jt
use types, only : rprec
use sgs_param, only: Nu_t

implicit none

integer :: k_ml_global
real(rprec), dimension(nz) :: Sc
real(rprec), dimension(ld,ny,nz), intent(out)::kappa_ts
integer :: k, k_abs

call thermocline_depth(k_ml_global) 

do k=1, nz
   k_abs = coord*(nz-1) + k
   if (k_ml_global.ge.k_abs) then
      Sc(k) = Sc_turb
   else
      Sc(k) = Sc_lam
   end if
      kappa_ts(:,:,k) = Nu_t(:,:,k)/Sc(k)
end do 

end subroutine schmidt_profile


!**********************************************************************
subroutine setsponge()              
!**********************************************************************
!Sets relaxation term to vertical momentum equation in top quarter
!of domain. Relaxation time scale is 50s with a factor of 5. 
! Reference: Nieuwstadt et al. (1991)

use param, only : coord, nproc, nz, lbz, z_i, nz_tot, u_star
use types, only : rprec
use sim_param, only : sponge

real(rprec) :: factor, sponge_top
integer :: k_abs, k

sponge_top = z_i / (50._rprec*u_star)

$if($MPI)
   factor = 9._rprec / (nz_tot - 3*nz_tot/4 +1)
   do k=1,nz
        k_abs = coord*(nz-1) + k
        if (k_abs.gt. 3*nz_tot/4 +1) then 
           sponge(k) = sponge_top*5._rprec**((k_abs-nz_tot)*factor)
        end if
   end do
$else
   factor = 9._rprec / (nz - 3*nz/4 +1)
   sponge(nz) = z_i / (50._rprec*u_star)
   do k=nz-1,3*nz/4+1,-1 
      sponge(k) = sponge(k+1)/5._rprec**factor
   end do
$endif

end subroutine setsponge


!**********************************************************************
subroutine ml_base_on_grid(z_ml)              
!**********************************************************************

use types, only: rprec
use param, only: ml_depth, coord, nz, dz, z_i, comm, ierr, MPI_RPREC
use mpi

integer :: k
real(rprec) :: z, z_ml_temp
real(rprec), intent(out) :: z_ml

do k=1,nz
   $if ($MPI)
   z = (coord * (nz-1) + k - 0.5_rprec) * dz * z_i
   $else
   z = (k - 0.5_rprec) * dz * z_i
   $endif
   if (z.le.ml_depth) then
      z_ml_temp = z
   end if
end do 


call mpi_allreduce(z_ml_temp, z_ml, 1, MPI_RPREC, MPI_MAX, comm, ierr)



end subroutine ml_base_on_grid

end module ocean_base
