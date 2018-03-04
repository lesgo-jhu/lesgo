!!
!!  Copyright (C) 2009-2013  Johns Hopkins University
!!
!!  This file is part of lesgo.
!!
!!  lesgo is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License as published by
!!  the Free Software Foundation, either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  lesgo is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with lesgo.  If not, see <http://www.gnu.org/licenses/>.
!!

subroutine ic()
  use types,only:rprec
  use param
  use sim_param,only:u,v,w
  use messages, only : error
  $if ($TURBINES)
  use turbines, only: turbine_vel_init
  $endif 
  
  implicit none

  character(*), parameter :: sub_name = 'ic'

  $if ($DEBUG)
  logical, parameter :: DEBUG = .false.
  $endif 

  $if ($TURBINES)
    real(rprec) :: zo_turbines
  $endif    
  
  integer::jx,jy,jz,seed
  integer :: jz_abs

  real(kind=rprec),dimension(nz)::ubar
  real(kind=rprec)::rms, noise, sigma_rv, arg, arg2
  real(kind=rprec)::z,w_star

  $if ($TURBINES)
  zo_turbines = 0._rprec
  $endif

  $if( $CPS ) 

  call boundary_layer_ic()

  $else
  
  if ( inflow ) then  !--no turbulence
     call uniform_ic()
  else
     call boundary_layer_ic()
  end if

  $endif

  !VK Display the mean vertical profiles of the initialized variables on the
  !screen
  !do jz=1,nz
  !   $if ($MPI)
  !   z = (coord*(nz-1) + jz - 0.5_rprec) * dz
  !   $else
  !   z = (jz - 0.5_rprec) * dz
  !   $endif
  !   write(6,7780) jz,z,sum(u(1:nx,:,jz))/float(nx*ny),sum(v(1:nx,:,jz))/&
  !        float(nx*ny),sum(w(1:nx,:,jz))/float(nx*ny)
  !end do
!7780 format('jz,z,ubar,vbar,wbar:',(1x,I3,4(1x,F9.4)))

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine uniform_ic()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

u = inflow_velocity 
v = 0.05_rprec * inflow_velocity
w = 0._rprec

return
end subroutine uniform_ic

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine boundary_layer_ic()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Log profile that is modified to flatten at z=z_i
! This is a hot mess (JSG 20111221)

use param, only: coriol, u_star, z_i
use mpi_defs
$if($CYL_SKEW_LS)
use cyl_skew_base_ls
$endif
implicit none

real(rprec), dimension(lbz:nz) :: ui_avg, vi_avg, wi_avg
real(rprec) :: z_hat
real(rprec), parameter :: nu = 4e-3_rprec   !Eshwan - eddy viscosity for McWilliams simulation
real(rprec) :: nu_bar                       !Eshwan - non-dimensionalize eddy viscosity
real(rprec) :: C, D                         !Constants for initial velocity profile from McWilliams et al et. 3.7
real(rprec) :: delta

interface
   function ran3(idum)
     integer :: idum
     real(8) :: ran3
   end function ran3
end interface

!w_star=(9.81_rprec/T_init*wt_s*z_i)**(1._rprec/3._rprec)
w_star = u_star

!if( coord == 0 ) write(*,*) '------> Creating modified log profile for IC'
do jz=1,nz

   $if ($MPI)
   z = (coord*(nz-1) + jz - 0.5_rprec) * dz
   $else
   z=(jz-.5_rprec)*dz
   $endif

   ! Another kludge for creating a channel profile. For now ignoring location of actual surface.
   $if($CYL_SKEW_LS)
   if( use_top_surf ) then
      if( z > L_z / 2 ) z = L_z - z
   endif
   $endif

   ! IC in equilibrium with rough surface (rough dominates in effective zo)
   arg2=z/zo
   arg=(1._rprec/vonk)*log(arg2)!-1./(2.*vonk*z_i*z_i)*z*z

   $if($LVLSET)
   ! Kludge to adjust magnitude of velocity profile
   ! Not critical - may delete
   arg = 0.357*arg
   $endif

   $if ($TURBINES)
   call turbine_vel_init (zo_turbines)
   arg2=z/zo_turbines
   arg=(1._rprec/vonk)*log(arg2)!-1./(2.*vonk*z_i*z_i)*z*z          
   $endif        

   !ubar(jz)=arg

   ! Added by VK for making the u less than 1...need to change this
   ! initialization routine

   if (coriolis_forcing) then
      ubar(jz)=arg/30._rprec
   else
      ubar(jz)=arg

   end if

   if ((coriolis_forcing).and.(z.gt.(.5_rprec))) ubar(jz)=ug

end do

!To generate an initial random velocity profile for each simulation
call init_random_seed()
call random_number(u)
call random_number(v)
call random_number(w)

!Center random number about 0 and rescale
rms = 1._rprec
sigma_rv = 0.289_rprec

u = (rms/sigma_rv)*(u-0.5_rprec)
v = (rms/sigma_rv)*(v-0.5_rprec)
w = (rms/sigma_rv)*(w-0.5_rprec)


!Constants for initial velocity profile - McWilliams et al. eq. 3.7 
nu_bar = nu/(u_star*z_i)
C = 1._rprec/(sqrt(2._rprec*coriol*nu_bar))
D = (1._rprec/sqrt(2._rprec)) * sqrt(coriol/nu_bar)

delta = sqrt(2._rprec*nu_bar/coriol)

do jz=1,nz
   $if ($MPI)
       jz_abs = coord * (nz-1) + jz
       z = (coord * (nz-1) + jz - 0.5_rprec) * dz * z_i          !dimensions in meters
       z_hat = -1._rprec*(z/z_i)                                 !non-dimensional - inverted z for ocean simulation 
   $else
       jz_abs = jz
       z = (jz_abs - .5_rprec) * dz * z_i                        !dimensions in meters
       z_hat = -1._rprec * (z/z_i)                               ! non-dimensional - inverted z for ocean simulation 
   $endif
   seed = -80 - jz_abs
   call random_seed(seed)
   do jy=1,ny
      do jx=1,nx
         !McWilliams simulation - Eshwan
         !Generate the same random initial profile for each simulation
         !if (z.le.ml_depth) then
         !   call random_number(noise)
         !   u(jx,jy,jz)=(rms/sigma_rv)*(noise-0.5_rprec)*(1._rprec-z/z_i) + C*exp(D*z_hat)*(cos(D*z_hat) + sin(D*z_hat)) 
         !   call random_number(noise)
         !   v(jx,jy,jz)=(rms/sigma_rv)*(noise-0.5_rprec)*(1._rprec-z/z_i) + C*exp(D*z_hat)*(-1._rprec*cos(D*z_hat) + sin(D*z_hat)) 
         !   call random_number(noise)
         !   w(jx,jy,jz)=(rms/sigma_rv)*(noise-0.5_rprec)*(1._rprec-z/z_i)
         !else
         !   call random_number(noise)
         !   u(jx,jy,jz)=(rms/sigma_rv)*(noise-0.5_rprec)*0.01_rprec + C*exp(D*z_hat)*(cos(D*z_hat) + sin(D*z_hat)) 
         !   call random_number(noise)
         !   v(jx,jy,jz)=(rms/sigma_rv)*(noise-0.5_rprec)*0.01_rprec + C*exp(D*z_hat)*(-1._rprec*cos(D*z_hat) + sin(D*z_hat)) 
         !   call random_number(noise)
         !   w(jx,jy,jz)=(rms/sigma_rv)*(noise-0.5_rprec)*0.01_rprec
         !end if

         !My simulations
         !if (z.le.ml_depth) then
         !   call random_number(noise)
         !   u(jx,jy,jz)=(rms/sigma_rv)*(noise-0.5_rprec)*(1._rprec-z/z_i) + exp(z_hat/delta)*cos(z_hat/delta) 
         !   !u(jx,jy,jz)=exp(z_hat/delta)*cos(z_hat/delta) 
         !   !u(jx,jy,jz)=(rms/sigma_rv)*(noise-0.5_rprec)*(1._rprec-z/z_i)  
         !   call random_number(noise)
         !   v(jx,jy,jz)=(rms/sigma_rv)*(noise-0.5_rprec)*(1._rprec-z/z_i) + exp(z_hat/delta)*sin(z_hat/delta) 
         !   !v(jx,jy,jz)=exp(z_hat/delta)*sin(z_hat/delta) 
         !   !v(jx,jy,jz)=(rms/sigma_rv)*(noise-0.5_rprec)*(1._rprec-z/z_i)  
         !   call random_number(noise)
         !   w(jx,jy,jz)=(rms/sigma_rv)*(noise-0.5_rprec)*(1._rprec-z/z_i)
         !else
         !   call random_number(noise)
         !   u(jx,jy,jz)=(rms/sigma_rv)*(noise-0.5_rprec)*0.01_rprec + exp(z_hat/delta)*cos(z_hat/delta) 
         !   !u(jx,jy,jz)=exp(z_hat/delta)*cos(z_hat/delta) 
         !   !u(jx,jy,jz)=(rms/sigma_rv)*(noise-0.5_rprec)*0.01_rprec  
         !   call random_number(noise)
         !   v(jx,jy,jz)=(rms/sigma_rv)*(noise-0.5_rprec)*0.01_rprec + exp(z_hat/delta)*sin(z_hat/delta) 
         !   !v(jx,jy,jz)=exp(z_hat/delta)*sin(z_hat/delta) 
         !   !v(jx,jy,jz)=(rms/sigma_rv)*(noise-0.5_rprec)*0.01_rprec  
         !   call random_number(noise)
         !   w(jx,jy,jz)=(rms/sigma_rv)*(noise-0.5_rprec)*0.01_rprec
         !end if
        
         !Generate a random initial profile for each simulation
         if (z.le.ml_depth) then
            u(jx,jy,jz)=u(jx,jy,jz)*(1._rprec-z/z_i) + exp(z_hat/delta)*cos(z_hat/delta) 
            v(jx,jy,jz)=v(jx,jy,jz)*(1._rprec-z/z_i) + exp(z_hat/delta)*sin(z_hat/delta) 
            w(jx,jy,jz)=w(jx,jy,jz)*(1._rprec-z/z_i) 
         else
            u(jx,jy,jz)=u(jx,jy,jz)*0.01_rprec + exp(z_hat/delta)*cos(z_hat/delta) 
            v(jx,jy,jz)=v(jx,jy,jz)*0.01_rprec + exp(z_hat/delta)*sin(z_hat/delta) 
            w(jx,jy,jz)=w(jx,jy,jz)*0.01_rprec
         end if

          !Logarithmic profile
!         if (z.le.z_i) then
!            noise=rms/.289_rprec*(ran3(seed)-.5_rprec)
!            u(jx,jy,jz)=noise*(1._rprec-z/z_i)*w_star/u_star+ubar(jz)
!            !u(jx,jy,jz)=(1._rprec-z/z_i)*w_star/u_star+ubar(jz)  !Eshwan
!            noise=rms/.289_rprec*(ran3(seed)-0.5_rprec)
!            v(jx,jy,jz)=noise*(1._rprec-z/z_i)*w_star/u_star 
!            !v(jx,jy,jz)=(1._rprec-z/z_i)*w_star/u_star           !Eshwan
!            noise=rms/.289_rprec*(ran3(seed)-.5_rprec)
!            w(jx,jy,jz)=noise*(1._rprec-z/z_i)*w_star/u_star
!            !w(jx,jy,jz)=(1._rprec-z/z_i)*w_star/u_star           !Eshwan
!         else
!            noise=rms/.289_rprec*(ran3(seed)-.5_rprec)
!            u(jx,jy,jz)=noise*w_star/u_star*.01_rprec+ubar(jz)
!            !u(jx,jy,jz)=w_star/u_star*.01_rprec+ubar(jz)         !Eshwan
!            noise=rms/.289_rprec*(ran3(seed)-0.5_rprec)
!            v(jx,jy,jz)=noise*w_star/u_star*.01_rprec                 
!            !v(jx,jy,jz)=w_star/u_star*.01_rprec                  !Eshwan
!            noise=rms/.289_rprec*(ran3(seed)-0.5_rprec)
!            w(jx,jy,jz)=noise*w_star/u_star*.01_rprec
!            !w(jx,jy,jz)=w_star/u_star*.01_rprec                  !Eshwan
!         end if
      end do
   end do
  !write(*,*) z, noise, noise
end do

!BC for W
if (coord == 0) then
   w(1:nx, 1:ny, 1) = 0._rprec
end if

!Set upper boundary condition as zero gradient in u and v and no peneration in w
$if ($MPI)
if (coord == nproc-1) then
   $endif    
   w(1:nx, 1:ny, nz) = 0._rprec
   u(1:nx, 1:ny, nz) = u(1:nx, 1:ny, nz-1)
   v(1:nx, 1:ny, nz) = v(1:nx, 1:ny, nz-1)
   $if ($MPI)
end if
$endif

!Exchange ghost node information for u, v, and w
$if($MPI)
    call mpi_sync_real_array(u, lbz, MPI_SYNC_DOWNUP)
    call mpi_sync_real_array(v, lbz, MPI_SYNC_DOWNUP)
    call mpi_sync_real_array(w, lbz, MPI_SYNC_DOWNUP)
$endif

!Average initial velocity horizontally
do jz=1, Nz
   
   $if ($MPI)
   jz_abs = coord * (nz-1) + jz
   z = (coord * (nz-1) + jz - 0.5_rprec) * dz * z_i    !dimensions in meters
   $else
   jz_abs = jz
   z = (jz-.5_rprec) * dz * z_i                        !dimensions in meters
   $endif
   
   ui_avg(jz) = 0.0_rprec
   vi_avg(jz) = 0.0_rprec
   wi_avg(jz) = 0.0_rprec
   do jy=1,Ny
      do jx=1,Nx
         ui_avg(jz) = ui_avg(jz) + u(jx,jy,jz)
         vi_avg(jz) = vi_avg(jz) + v(jx,jy,jz)
         wi_avg(jz) = wi_avg(jz) + w(jx,jy,jz)
      end do
   end do
   ui_avg(jz) = ui_avg(jz)/real(Nx*Ny) 
   vi_avg(jz) = vi_avg(jz)/real(Nx*Ny) 
   wi_avg(jz) = wi_avg(jz)/real(Nx*Ny) 

!   print*, 'Horizontally averaged initial velocity'
    !write(*,*), z, ui_avg(jz), vi_avg(jz)
end do

!if (coord==0) write(*,*) 'dz, nz_tot, z_i, nu_bar', dz, nz_tot, z_i, nu_bar

return
end subroutine boundary_layer_ic

end subroutine ic
