!!
!!  Copyright (C) 2010-2013  Johns Hopkins University
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

!**********************************************************************
module derivatives
!**********************************************************************
! 
! This module contains all of the major subroutines used for computing
! derivatives.
!
implicit none

save
private

public ddx, &
     ddy, &
     ddxy, &
     ddz_uv, &
     ddz_w, &
     filt_da

contains

!**********************************************************************
subroutine ddx(f,dfdx,lbz)              
!**********************************************************************
!
!  This subroutine computes the partial derivative of f with respect to
!  x using spectral decomposition.
!  

use types,only:rprec
use param,only:ld,nx,ny,nz
use fft
use emul_complex, only : OPERATOR(.MULI.)
implicit none

integer::jz

integer, intent(in) :: lbz
real(rprec), dimension(:,:,lbz:), intent(in) :: f
real(rprec), dimension(:,:,lbz:), intent(inout) :: dfdx

real(rprec) :: const

const = 1._rprec / ( nx * ny )

! Loop through horizontal slices
do jz=lbz,nz

  !  Use dfdx to hold f; since we are doing IN_PLACE FFT's this is required
  dfdx(:,:,jz)=const*f(:,:,jz)
  call dfftw_execute_dft_r2c(forw, dfdx(:,:,jz),dfdx(:,:,jz))

  ! Zero padded region and Nyquist frequency
  !  dfdx_c(lh,:,jz)=0._rprec ! Complex version
  !  dfdx_c(:,ny/2+1,jz)=0._rprec !  Complex version
  dfdx(ld-1:ld,:,jz) = 0._rprec
  dfdx(:,ny/2+1,jz) = 0._rprec

  ! Compute coefficients for pseudospectral derivative calculation
  !  dfdx_c(:,:,jz)=eye*kx(:,:)*dfdx_c(:,:,jz) ! Complex version

  !  Use complex emulation of dfdx to perform complex multiplication
  !  Optimized version for real(eye*kx)=0; only passing imaginary part of eye*kx
  dfdx(:,:,jz) = dfdx(:,:,jz) .MULI. kx

  ! Perform inverse transform to get pseudospectral derivative
  call dfftw_execute_dft_c2r(back, dfdx(:,:,jz), dfdx(:,:,jz))

enddo

return

end subroutine ddx

!**********************************************************************
subroutine ddy(f,dfdy, lbz)              
!**********************************************************************
!
!  This subroutine computes the partial derivative of f with respect to
!  y using spectral decomposition.
!  
use types,only:rprec
use param,only:ld,nx,ny,nz
use fft
use emul_complex, only : OPERATOR(.MULI.)
implicit none      

integer::jz
  
integer, intent(in) :: lbz
real(rprec), dimension(:,:,lbz:), intent(in) :: f
real(rprec), dimension(:,:,lbz:), intent(inout) :: dfdy

real(rprec) :: const

const = 1._rprec / ( nx * ny )

! Loop through horizontal slices
do jz=lbz,nz    

  !  Use dfdy to hold f; since we are doing IN_PLACE FFT's this is required
  dfdy(:,:,jz)=const*f(:,:,jz)  
  call dfftw_execute_dft_r2c(forw, dfdy(:,:,jz), dfdy(:,:,jz))

  ! Zero padded region and Nyquist frequency
  !  dfdy_c(lh,:,jz)=0._rprec ! Complex version
  !  dfdy_c(:,ny/2+1,jz)=0._rprec !  Complex version
  dfdy(ld-1:ld,:,jz) = 0._rprec
  dfdy(:,ny/2+1,jz) = 0._rprec   

  ! Compute coefficients for pseudospectral derivative calculation
  !  dfdy_c(:,:,jz)=eye*ky(:,:)*dfdy_c(:,:,jz) ! Complex version

  !  Use complex emulation of dfdy to perform complex multiplication
  !  Optimized version for real(eye*ky)=0; only passing imaginary part of eye*ky
  dfdy(:,:,jz) = dfdy(:,:,jz) .MULI. ky

  ! Perform inverse transform to get pseudospectral derivative
  call dfftw_execute_dft_c2r(back, dfdy(:,:,jz), dfdy(:,:,jz))
end do

return

end subroutine ddy

!**********************************************************************
subroutine ddxy (f, dfdx, dfdy, lbz)              
!**********************************************************************
use types,only:rprec
use param,only:ld,nx,ny,nz
use fft
use emul_complex, only : OPERATOR(.MULI.)
implicit none
integer::jz

integer, intent(in) :: lbz
! only need complex treatment
real(rprec), dimension(:,:,lbz:), intent(in) :: f
real(rprec), dimension(:,:,lbz:), intent(inout) :: dfdx,dfdy

real(rprec) :: const

const = 1._rprec / ( nx * ny )

!...Loop through horizontal slices
do jz=lbz,nz
   ! temporay storage in dfdx_c, this was don't mess up f_c
   dfdx(:,:,jz) = const*f(:,:,jz)   
   call dfftw_execute_dft_r2c(forw, dfdx(:,:,jz), dfdx(:,:,jz))

   !dfdx_c(lh,:,jz)=0._rprec
   !dfdx_c(:,ny/2+1,jz)=0._rprec
   dfdx(ld-1:ld,:,jz)=0._rprec
   dfdx(:,ny/2+1,jz)=0._rprec

! derivatives: must to y's first here, because we're using dfdx as storage
   !dfdy_c(:,:,jz)=eye*ky(:,:)*dfdx_c(:,:,jz)
   !dfdx_c(:,:,jz)=eye*kx(:,:)*dfdx_c(:,:,jz)
   dfdy(:,:,jz) = dfdx(:,:,jz) .MULI. ky
   dfdx(:,:,jz) = dfdx(:,:,jz) .MULI. kx

! the oddballs for derivatives should already be dead, since they are for f

! inverse transform 
   call dfftw_execute_dft_c2r(back, dfdx(:,:,jz), dfdx(:,:,jz))
   call dfftw_execute_dft_c2r(back, dfdy(:,:,jz), dfdy(:,:,jz))
 
end do
end subroutine ddxy

!**********************************************************************
subroutine ddz_uv(f, dfdz, lbz)
!**********************************************************************
!
! first derivative in z direction for boundary layer (2nd order numerics)
!--provides dfdz(:, :, 2:nz), value at jz=1 is not touched
!--MPI: provides dfdz(:, :, 1:nz) where f is on uvp-node, except at
!  bottom process it only supplies 2:nz
!
use types,only:rprec
use param,only:nx,ny,nz,dz,BOGUS
#ifdef PPMPI
use param,only:nproc,coord
#endif
implicit none

integer, intent(in) :: lbz
real (rprec), dimension (:, :, lbz:), intent (in) :: f
real (rprec), dimension(:, :, lbz:), intent (inout) :: dfdz

integer::jx,jy,jz
real (rprec) :: const

const=1._rprec/dz

#ifdef PPMPI

#ifdef PPSAFETYMODE
  dfdz(:, :, 0) = BOGUS
#endif

  if (coord > 0) then
    !--ghost node information is available here
    !--watch the z-dimensioning!
    !--if coord == 0, dudz(1) will be set in wallstress
    do jy=1,ny
    do jx=1,nx    
       dfdz(jx,jy,1)=const*(f(jx,jy,1)-f(jx,jy,0))
    end do
    end do
  end if

#endif

do jz=2,nz-1
do jy=1,ny
do jx=1,nx    
   dfdz(jx,jy,jz)=const*(f(jx,jy,jz)-f(jx,jy,jz-1))
end do
end do
end do

!--should integrate this into the above loop, explicit zeroing is not
!  needed, since dudz, dvdz are forced to zero by copying the u, v fields
!--also, what happens when called with tzz? 
#ifdef PPMPI 

  if (coord == nproc-1) then
    dfdz(:,:,nz)=0._rprec  !--do not need to do this...
  else
    do jy=1,ny
    do jx=1,nx
       dfdz(jx,jy,nz)=const*(f(jx,jy,nz)-f(jx,jy,nz-1))
    end do
    end do
  endif

#else

  dfdz(:,:,nz)=0._rprec  !--do not need to do this...

#endif

return
end subroutine ddz_uv

!**********************************************************************
subroutine ddz_w(f, dfdz, lbz)
!**********************************************************************
!
!  First deriv in z direction for boundary layer (2nd order numerics)
!  F is on w nodes and dFdz is on uvp nodes
!    -23 January 1996
!--MPI: provides 0:nz-1, except top process has 0:nz, and bottom process
!  has 1:nz-1
!
use types,only:rprec
use param,only:nx,ny,nz,dz,BOGUS
#ifdef PPMPI
use param,only:nproc,coord
#endif
implicit none

integer, intent(in) :: lbz
real (rprec), dimension (:, :, lbz:), intent (in) :: f
real(kind=rprec),dimension(:, :, lbz:), intent (inout) :: dfdz

real(kind=rprec)::const
integer::jx,jy,jz

const=1._rprec/dz
do jz=lbz,nz-1
do jy=1,ny
do jx=1,nx
   dfdz(jx,jy,jz)=const*(f(jx,jy,jz+1)-f(jx,jy,jz))
end do
end do
end do

#ifdef PPMPI
    if (coord == 0) then
      !--bottom process cannot calculate dfdz(jz=0)
#ifdef PPSAFETYMODE
      dfdz(:, :, lbz) = BOGUS
#endif      
    endif
    if (coord == nproc-1) then
      dfdz(:,:,nz)=0._rprec !dfdz(:,:,Nz-1) ! c? any better ideas for sponge?
    else
#ifdef PPSAFETYMODE
      dfdz(:, :, nz) = BOGUS
#endif      
    end if
#else
  dfdz(:,:,nz)=0._rprec !dfdz(:,:,Nz-1) ! c? any better ideas for sponge?
#endif


end subroutine ddz_w

!**********************************************************************
subroutine filt_da(f,dfdx,dfdy, lbz)
!**********************************************************************
!
! kills the oddball components in f, and calculates in horizontal derivatives
!--supplies results for jz=$lbz:nz
!--MPI: all these arrays are 0:nz
!--MPI: on exit, u,dudx,dudy,v,dvdx,dvdy,w,dwdx,dwdy valid on jz=0:nz,
!  except on bottom process (0 level set to BOGUS, starts at 1)
!
use types,only:rprec
use param,only:ld,nx,ny,nz
use fft
use emul_complex, only : OPERATOR(.MULI.)
implicit none
integer::jz

integer, intent(in) :: lbz
real(rprec), dimension(:, :, lbz:), intent(inout) :: f
real(rprec), dimension(:, :, lbz:), intent(inout) :: dfdx, dfdy

real(rprec) :: const

const = 1._rprec/(nx*ny)

! loop through horizontal slices
do jz=lbz,nz

  f(:,:,jz)=const*f(:,:,jz)  
  call dfftw_execute_dft_r2c(forw, f(:,:,jz), f(:,:,jz))

  ! what exactly is the filter doing here? in other words, why is routine
  ! called filt_da? not filtering anything

  ! Zero padded region and Nyquist frequency
  !  f_c(lh,:,jz)=0._rprec ! Complex version
  !  f_c(:,ny/2+1,jz)=0._rprec !  Complex version
  f(ld-1:ld,:,jz) = 0._rprec
  f(:,ny/2+1,jz) = 0._rprec   

  !  Compute in-plane derivatives
  !  dfdy_c(:,:,jz)=eye*ky(:,:)*f_c(:,:,jz) !  complex version
  !  dfdx_c(:,:,jz)=eye*kx(:,:)*f_c(:,:,jz) !  complex version
  dfdx(:,:,jz) = f(:,:,jz) .MULI. kx
  dfdy(:,:,jz) = f(:,:,jz) .MULI. ky

  ! the oddballs for derivatives should already be dead, since they are for f
  ! inverse transform 
  call dfftw_execute_dft_c2r(back, f(:,:,jz), f(:,:,jz))
  call dfftw_execute_dft_c2r(back, dfdx(:,:,jz), dfdx(:,:,jz))
  call dfftw_execute_dft_c2r(back, dfdy(:,:,jz), dfdy(:,:,jz))

end do

return
end subroutine filt_da

end module derivatives


