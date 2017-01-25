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
     filt_da, &
     !ddx_direct, &
     !ddy_only, &
     !dft_direct_forw_2d, &
     !dft_direct_back_2d, &
     !filt_da_direct, &
     !dft_direct_forw_2d_new, &
     !dft_direct_back_2d_new, &
     dft_direct_forw_2d_n, &
     dft_direct_back_2d_n, &
     dft_direct_forw_2d_n_yonly, &
     dft_direct_back_2d_n_yonly, &
     dft_direct_forw_2d_n_yonlyC, &
     dft_direct_back_2d_n_yonlyC, &
     dft_direct_forw_2d_n_yonlyC_big, &
     dft_direct_back_2d_n_yonlyC_big, &
     dft_direct_forw_2d_n_big, &
     dft_direct_back_2d_n_big, &
     !filt_da_direct_n, &
     filt_da_kxspace, &
     !ddx_n, &
     !ddy_n, &
     phys2wave, &
     phys2wave_pr, &
     wave2phys, &
     wave2phys_pr, &
     phys2waveF, &
     wave2physF, &
     phys2waveF_pr, &
     wave2physF_pr, &
     ddx_kxspace, &
     ddy_kxspace, &
     convolve, &
     convolve2, &
     convolve_rnl, &
     ky_spectra_calc_jb

contains

!**********************************************************************
subroutine ddx(f,dfdx,lbz)              
!**********************************************************************
!  This subroutine computes the partial derivative of f with respect to
!  x using spectral decomposition.

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
  $if ($FFTW3)
  call dfftw_execute_dft_r2c(forw, dfdx(:,:,jz),dfdx(:,:,jz))
  $else
  call rfftwnd_f77_one_real_to_complex(forw,dfdx(:,:,jz),fftwNull_p)       
  $endif

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
  $if ($FFTW3)
  call dfftw_execute_dft_c2r(back, dfdx(:,:,jz), dfdx(:,:,jz))
  $else
  call rfftwnd_f77_one_complex_to_real(back,dfdx(:,:,jz),fftwNull_p)
  $endif

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
  $if ($FFTW3)
  call dfftw_execute_dft_r2c(forw, dfdy(:,:,jz), dfdy(:,:,jz))
  $else
  call rfftwnd_f77_one_real_to_complex(forw,dfdy(:,:,jz),fftwNull_p)     
  $endif

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
  $if ($FFTW3)
  call dfftw_execute_dft_c2r(back, dfdy(:,:,jz), dfdy(:,:,jz))
  $else
  call rfftwnd_f77_one_complex_to_real(back,dfdy(:,:,jz),fftwNull_p)     
  $endif
end do

return

end subroutine ddy

!!$!**********************************************************************
!!$subroutine ddx_n(f,dfdx,lbz)              
!!$!**********************************************************************
!!$!  This subroutine computes the partial derivative of f with respect to
!!$!  x using spectral decomposition.
!!$
!!$use types,only:rprec
!!$use param,only:ld,nx,ny,nz
!!$use fft
!!$use emul_complex, only : OPERATOR(.MULI.)
!!$implicit none
!!$
!!$integer::jz
!!$
!!$integer, intent(in) :: lbz
!!$real(rprec), dimension(:,:,lbz:), intent(in) :: f
!!$real(rprec), dimension(:,:,lbz:), intent(inout) :: dfdx
!!$
!!$real(rprec) :: const
!!$
!!$const = 1._rprec / ( nx * ny )
!!$
!!$! Loop through horizontal slices
!!$do jz=lbz,nz
!!$
!!$  !  Use dfdx to hold f; since we are doing IN_PLACE FFT's this is required
!!$  dfdx(:,:,jz)=const*f(:,:,jz)
!!$  $if ($FFTW3)
!!$  !!call dfftw_execute_dft_r2c(forw, dfdx(:,:,jz),dfdx(:,:,jz))
!!$  call dft_direct_forw_2d_n( dfdx(:,:,jz)  )      !!jb
!!$  $else
!!$  call rfftwnd_f77_one_real_to_complex(forw,dfdx(:,:,jz),fftwNull_p)       
!!$  $endif
!!$
!!$  ! Zero padded region and Nyquist frequency
!!$  !  dfdx_c(lh,:,jz)=0._rprec ! Complex version
!!$  !  dfdx_c(:,ny/2+1,jz)=0._rprec !  Complex version
!!$  dfdx(ld-1:ld,:,jz) = 0._rprec
!!$  dfdx(:,ny/2+1,jz) = 0._rprec
!!$
!!$  ! Compute coefficients for pseudospectral derivative calculation
!!$  !  dfdx_c(:,:,jz)=eye*kx(:,:)*dfdx_c(:,:,jz) ! Complex version
!!$
!!$  !  Use complex emulation of dfdx to perform complex multiplication
!!$  !  Optimized version for real(eye*kx)=0; only passing imaginary part of eye*kx
!!$  dfdx(:,:,jz) = dfdx(:,:,jz) .MULI. kx
!!$
!!$  ! Perform inverse transform to get pseudospectral derivative
!!$  $if ($FFTW3)
!!$  !!call dfftw_execute_dft_c2r(back, dfdx(:,:,jz), dfdx(:,:,jz))
!!$  call dft_direct_back_2d_n( dfdx(:,:,jz)  )      !!jb
!!$  $else
!!$  call rfftwnd_f77_one_complex_to_real(back,dfdx(:,:,jz),fftwNull_p)
!!$  $endif
!!$
!!$enddo
!!$
!!$return
!!$
!!$end subroutine ddx_n

!!$!**********************************************************************
!!$subroutine ddy_n(f,dfdy, lbz)              
!!$!**********************************************************************
!!$!
!!$!  This subroutine computes the partial derivative of f with respect to
!!$!  y using spectral decomposition.
!!$!  
!!$use types,only:rprec
!!$use param,only:ld,nx,ny,nz
!!$use fft
!!$use emul_complex, only : OPERATOR(.MULI.)
!!$implicit none      
!!$
!!$integer::jz
!!$  
!!$integer, intent(in) :: lbz
!!$real(rprec), dimension(:,:,lbz:), intent(in) :: f
!!$real(rprec), dimension(:,:,lbz:), intent(inout) :: dfdy
!!$
!!$real(rprec) :: const
!!$
!!$const = 1._rprec / ( nx * ny )
!!$
!!$! Loop through horizontal slices
!!$do jz=lbz,nz    
!!$
!!$  !  Use dfdy to hold f; since we are doing IN_PLACE FFT's this is required
!!$  dfdy(:,:,jz)=const*f(:,:,jz)  
!!$  $if ($FFTW3)
!!$  !!call dfftw_execute_dft_r2c(forw, dfdy(:,:,jz), dfdy(:,:,jz))
!!$  call dft_direct_forw_2d_n( dfdy(:,:,jz)  )      !!jb
!!$  $else
!!$  call rfftwnd_f77_one_real_to_complex(forw,dfdy(:,:,jz),fftwNull_p)     
!!$  $endif
!!$
!!$  ! Zero padded region and Nyquist frequency
!!$  !  dfdy_c(lh,:,jz)=0._rprec ! Complex version
!!$  !  dfdy_c(:,ny/2+1,jz)=0._rprec !  Complex version
!!$  dfdy(ld-1:ld,:,jz) = 0._rprec
!!$  dfdy(:,ny/2+1,jz) = 0._rprec   
!!$
!!$  ! Compute coefficients for pseudospectral derivative calculation
!!$  !  dfdy_c(:,:,jz)=eye*ky(:,:)*dfdy_c(:,:,jz) ! Complex version
!!$
!!$  !  Use complex emulation of dfdy to perform complex multiplication
!!$  !  Optimized version for real(eye*ky)=0; only passing imaginary part of eye*ky
!!$  dfdy(:,:,jz) = dfdy(:,:,jz) .MULI. ky
!!$
!!$  ! Perform inverse transform to get pseudospectral derivative
!!$  $if ($FFTW3)
!!$  !!call dfftw_execute_dft_c2r(back, dfdy(:,:,jz), dfdy(:,:,jz))
!!$  call dft_direct_back_2d_n( dfdy(:,:,jz)  )      !!jb
!!$  $else
!!$  call rfftwnd_f77_one_complex_to_real(back,dfdy(:,:,jz),fftwNull_p)     
!!$  $endif
!!$end do
!!$
!!$return
!!$
!!$end subroutine ddy_n

!**********************************************************************
subroutine ddx_kxspace(f,dfdx,lbz)              
!**********************************************************************
!  This subroutine computes the partial derivative of f with respect to
!  x using spectral decomposition.

use types,only:rprec
use param,only:ld,nx,ny,nz
use fft
use emul_complex, only : OPERATOR(.MULI.)
implicit none

integer::jz

integer, intent(in) :: lbz
real(rprec), dimension(:,:,lbz:), intent(in) :: f
real(rprec), dimension(:,:,lbz:), intent(inout) :: dfdx

! Loop through horizontal slices
do jz=lbz,nz

  dfdx(:,:,jz) = f(:,:,jz)

  ! Zero padded region and Nyquist frequency
  !  dfdx_c(lh,:,jz)=0._rprec ! Complex version
  !  dfdx_c(:,ny/2+1,jz)=0._rprec !  Complex version
  dfdx(ld-1:ld,:,jz) = 0._rprec
  dfdx(:,ny/2+1,jz) = 0._rprec

  !  Use complex emulation of dfdx to perform complex multiplication
  !  Optimized version for real(eye*kx)=0; only passing imaginary part of eye*kx
  dfdx(:,:,jz) = dfdx(:,:,jz) .MULI. kx

enddo

return
end subroutine ddx_kxspace

!**********************************************************************
subroutine ddy_kxspace(f,dfdy, lbz)              
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

! Loop through horizontal slices
do jz=lbz,nz    

  dfdy(:,:,jz) = f(:,:,jz)

  ! Zero padded region and Nyquist frequency
  !  dfdy_c(lh,:,jz)=0._rprec ! Complex version
  !  dfdy_c(:,ny/2+1,jz)=0._rprec !  Complex version
  dfdy(ld-1:ld,:,jz) = 0._rprec
  dfdy(:,ny/2+1,jz) = 0._rprec   

  !  Use complex emulation of dfdy to perform complex multiplication
  !  Optimized version for real(eye*ky)=0; only passing imaginary part of eye*ky
  dfdy(:,:,jz) = dfdy(:,:,jz) .MULI. ky

end do

return
end subroutine ddy_kxspace

!**********************************************************************
subroutine ddxy (f, dfdx, dfdy, lbz)              
!**********************************************************************
use types,only:rprec
use param,only:nx,ny,nz
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
   $if ($FFTW3)
   call dfftw_execute_dft_r2c(forw, dfdx(:,:,jz), dfdx(:,:,jz))
   $else
   call rfftwnd_f77_one_real_to_complex(forw,dfdx(:,:,jz),fftwNull_p)
   $endif

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
   $if ($FFTW3)
   call dfftw_execute_dft_c2r(back, dfdx(:,:,jz), dfdx(:,:,jz))
   call dfftw_execute_dft_c2r(back, dfdy(:,:,jz), dfdy(:,:,jz))
   $else
   call rfftwnd_f77_one_complex_to_real(back,dfdx(:,:,jz),fftwNull_p)
   call rfftwnd_f77_one_complex_to_real(back,dfdy(:,:,jz),fftwNull_p)
   $endif
 
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
use param,only:nx,ny,nz,dz,coord,nproc,BOGUS,ld
implicit none

integer, intent(in) :: lbz
real (rprec), dimension (:, :, lbz:), intent (in) :: f
real (rprec), dimension(:, :, lbz:), intent (inout) :: dfdz

integer::jx,jy,jz
real (rprec) :: const

const=1._rprec/dz

$if ($MPI)

$if ($SAFETYMODE)
  dfdz(:, :, 0) = BOGUS
$endif

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

$endif

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
$if ($MPI) 

  if (coord == nproc-1) then
    dfdz(:,:,nz)=0._rprec  !--do not need to do this...
  else
    do jy=1,ny
    do jx=1,ld !nx    !!jb   !! fourier   !! fourier
       dfdz(jx,jy,nz)=const*(f(jx,jy,nz)-f(jx,jy,nz-1))
    end do
    end do
  endif

$else

  dfdz(:,:,nz)=0._rprec  !--do not need to do this...

$endif

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
use param,only:nx,ny,nz,dz,coord,nproc,BOGUS, ld
implicit none

integer, intent(in) :: lbz
real (rprec), dimension (:, :, lbz:), intent (in) :: f
real(kind=rprec),dimension(:, :, lbz:), intent (inout) :: dfdz

real(kind=rprec)::const
integer::jx,jy,jz

const=1._rprec/dz
do jz=lbz,nz-1
do jy=1,ny
do jx=1,ld !nx    !!jb   !! fourier   !! fourier
   dfdz(jx,jy,jz)=const*(f(jx,jy,jz+1)-f(jx,jy,jz))
end do
end do
end do

$if ($MPI)
    if (coord == 0) then
      !--bottom process cannot calculate dfdz(jz=0)
$if ($SAFETYMODE)
      dfdz(:, :, lbz) = BOGUS
$endif      
    endif
    if (coord == nproc-1) then
      dfdz(:,:,nz)=0._rprec !dfdz(:,:,Nz-1) ! c? any better ideas for sponge?
    else
$if ($SAFETYMODE)
      dfdz(:, :, nz) = BOGUS
$endif      
    end if
$else
  dfdz(:,:,nz)=0._rprec !dfdz(:,:,Nz-1) ! c? any better ideas for sponge?
$endif


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
use param,only:ld,nx,ny,nz,coord
use fft
use emul_complex, only : OPERATOR(.MULI.)
implicit none
integer::jz,jy,jx

integer, intent(in) :: lbz
real(rprec), dimension(:, :, lbz:), intent(inout) :: f
real(rprec), dimension(:, :, lbz:), intent(inout) :: dfdx, dfdy

real(rprec) :: const

const = 1._rprec/(nx*ny)

! loop through horizontal slices
do jz=lbz,nz

!!$if(coord==0 .and. jz==1) then
!!$     write(*,*) 'filt_da 1'
!!$     do jx=1,nx+2
!!$        do jy=1,ny
!!$           write(*,*) jx,jy,f(jx,jy,jz)
!!$        enddo
!!$     enddo
!!$  endif



  f(:,:,jz)=const*f(:,:,jz)  

  $if ($FFTW3)
  call dfftw_execute_dft_r2c(forw, f(:,:,jz), f(:,:,jz))
  $else
  call rfftwnd_f77_one_real_to_complex(forw,f(:,:,jz),fftwNull_p)
  $endif

!!$  if(coord==0 .and. jz==1) then
!!$     write(*,*) 'filt_da 2'
!!$     do jx=1,nx+2
!!$        do jy=1,ny
!!$           write(*,*) jx,jy,f(jx,jy,jz)
!!$        enddo
!!$     enddo
!!$  endif


  ! what exactly is the filter doing here? in other words, why is routine
  ! called filt_da? not filtering anything

  ! Zero padded region and Nyquist frequency
  !  f_c(lh,:,jz)=0._rprec ! Complex version
  !  f_c(:,ny/2+1,jz)=0._rprec !  Complex version
  f(ld-1:ld,:,jz) = 0._rprec
  f(:,ny/2+1,jz) = 0._rprec   


!if(coord==0 .and. jz==1) then      !!jb
!   write(*,*) 'filt_da 3'
!   write(*,*) f(:,:,jz)
!endif


  !  Compute in-plane derivatives
  !  dfdy_c(:,:,jz)=eye*ky(:,:)*f_c(:,:,jz) !  complex version
  !  dfdx_c(:,:,jz)=eye*kx(:,:)*f_c(:,:,jz) !  complex version
  dfdx(:,:,jz) = f(:,:,jz) .MULI. kx
  dfdy(:,:,jz) = f(:,:,jz) .MULI. ky

  ! the oddballs for derivatives should already be dead, since they are for f
  ! inverse transform 
  $if ($FFTW3)
  call dfftw_execute_dft_c2r(back, f(:,:,jz), f(:,:,jz))
  call dfftw_execute_dft_c2r(back, dfdx(:,:,jz), dfdx(:,:,jz))
  call dfftw_execute_dft_c2r(back, dfdy(:,:,jz), dfdy(:,:,jz))
  $else    
  call rfftwnd_f77_one_complex_to_real(back,   f(:,:,jz),fftwNull_p)
  call rfftwnd_f77_one_complex_to_real(back,dfdx(:,:,jz),fftwNull_p)
  call rfftwnd_f77_one_complex_to_real(back,dfdy(:,:,jz),fftwNull_p)
  $endif  

!!$if(coord==0 .and. jz==1) then
!!$     write(*,*) 'filt_da 4'
!!$     do jx=1,nx+2
!!$        do jy=1,ny
!!$           write(*,*) jx,jy,f(jx,jy,jz)
!!$        enddo
!!$     enddo
!!$  endif

end do

return
end subroutine filt_da

!!**********************************************************************
!subroutine ddx_direct(f, dfdx, lbz)
!!**********************************************************************
!!
!! Computes x derivative spectrally using direct summation instead of
!! FFTs. For use with RNL.
!!
!use types,only:rprec
!use param,only:ld,nx,ny,nz,kx_num,kx_vec,pi,coord,L_x,nproc
!implicit none

!integer :: jx, jy, jz, k
!integer, intent(in) :: lbz
!real(rprec), dimension(:, :, lbz:), intent(in) :: f
!real(rprec), dimension(:, :, lbz:), intent(out) :: dfdx
!real(rprec) :: const

!real(rprec), dimension(kx_num) :: pre
!complex(rprec), dimension(kx_num) :: fhat
!complex(rprec) :: imag = (0.0, 1.0)
!complex(rprec) :: const1, const2

!pre(:) = 2._rprec
!pre(1) = 1._rprec
!const1 = -L_x/nx * imag
!const2 = -1 * const1

!!!$if (coord == 0) then
!!!$   print*, 'pre: ', pre
!!!$   print*, 'const1: ', const1
!!!$   print*, 'const2: ', const2
!!!$endif

!dfdx(:,:,:) = 0._rprec

!do jy=1,ny
!do jz=lbz,nz

   !fhat(:)  = ( 0._rprec, 0.0_rprec  )

   !do  k=1,kx_num
      !do jx=1,nx
         !fhat(k) = fhat(k) + f(jx,jy,jz) * exp(const1 * kx_vec(k) * real(jx-1,rprec) )
      !enddo
   !enddo

!!!$   if (coord == nproc-1 .and. jz==1 .and. jy == 1) then
!!!$   print*, 'Before deriv: '
!!!$   do k=1,kx_num
!!!$      write(*,*) 'k, fhat(k): ', k, fhat(k)
!!!$   enddo
!!!$   endif

   !do k=1,kx_num
      !fhat(k) = imag * kx_vec(k) * fhat(k)
   !enddo

!!!$   if (coord == nproc-1 .and. jz==1 .and. jy == 1) then
!!!$   print*, 'After deriv: '
!!!$   do k=1,kx_num
!!!$      write(*,*) 'k, fhat(k): ', k, fhat(k)
!!!$   enddo
!!!$   endif
   !!!fhat(5) = (0._rprec, 0._rprec)
   !do jx=1,nx
      !do   k=1,kx_num
         !dfdx(jx,jy,jz) = dfdx(jx,jy,jz) + pre(k)*fhat(k)*exp(const2 * kx_vec(k) * real(jx-1,rprec))
      !enddo
   !enddo
   
!enddo
!enddo

!dfdx = dfdx / real(nx,rprec)

!return
!end subroutine ddx_direct

!!**********************************************************************
!subroutine ddy_only(f,dfdy, lbz)              
!!**********************************************************************
!!
!!  This subroutine computes the partial derivative of f with respect to
!!  y using spectral decomposition (using 1d FFT instead of 2d)
!!  
!use types,only:rprec
!use param,only:ld,nx,ny,nz,coord
!use fft
!implicit none      

!integer :: jx, jy, jz
!integer, intent(in) :: lbz
!real(rprec), dimension(:,:,lbz:), intent(in) :: f
!real(rprec), dimension(:,:,lbz:), intent(inout) :: dfdy
!real(rprec) :: const
!complex(rprec) :: imag = (0.0, 1.0)
!complex(rprec), dimension(ny/2+1) :: fhat

!const = 1._rprec / ny

!! Loop through horizontal slices
!do jz=lbz,nz    
!do jx=1,nx

  !!  Use dfdy to hold f; since we are doing IN_PLACE FFT's this is required
  !dfdy(jx,:,jz) = const * f(jx,:,jz)  
  !$if ($FFTW3)

  !call dfftw_execute_dft_r2c(forw_y, dfdy(jx,:,jz), fhat(:) )
  !$else
  !write(*,*) 'ERROR - no DFT available!'
  !$endif

  !fhat(ny/2+1) = (0._rprec, 0._rprec)

  !if (coord==0) then
     !if (jx==1 .and. jz==1) then
        !write(*,*) 'y fft'
        !write(*,*) fhat(:)
     !endif
  !endif

  !do jy=1,ny/2 !+1
     !fhat(jy) = fhat(jy) * ky_vec(jy) * imag
  !enddo

  !! Perform inverse transform to get pseudospectral derivative
  !$if ($FFTW3)
  !call dfftw_execute_dft_c2r(back_y, fhat(:), dfdy(jx,:,jz))
  !$else
  !write(*,*) 'ERROR - no DFT available!'
  !$endif

!enddo
!enddo

!return

!end subroutine ddy_only

!!**********************************************************************
!subroutine dft_direct_forw_2d(f)
!!**********************************************************************
!!
!! Computes 2d DFT directly, no FFT.
!!
!use types,only:rprec
!use param,only:ld,nx,ny,nz,kx_num,kx_vec,pi,coord,L_x,nproc,L_y
!use fft,only:ky_vec
!implicit none

!integer :: jx, jy, i, j, ii, ir, jy_r
!real(rprec), dimension(:, :), intent(inout) :: f
!integer :: ky_num

!!complex(rprec), dimension(kx_num, ny/2+1) :: fhat
!complex(rprec), dimension(nx, ny) :: fhat
!complex(rprec) :: imag = (0.0_rprec, 1.0_rprec)
!complex(rprec) :: c1, c2

!ky_num = ny/2+1

!c1 = -L_x/real(nx,rprec) * imag
!c2 = -L_y/real(ny,rprec) * imag

!fhat(:,:)  = ( 0._rprec, 0._rprec  )

!do jx=1,kx_num
!do jy=1,ny
!do i=1,nx
!do j=1,ny
!fhat(jx,jy)=fhat(jx,jy)+f(i,j)*exp(c1*kx_vec(jx)*real(i-1,rprec) + c2*ky_vec(jy)*real(j-1,rprec));
!end do
!end do
!end do
!end do

!!fhat = fhat / real(nx*ny,rprec)

!!!$if (coord == 0) then
!!!$   write(*,*) 'dft forw'
!!!$do jx=1,kx_num
!!!$do jy=1,ny
!!!$   write(*,*) jx,jy,fhat(jx,jy)
!!!$enddo
!!!$enddo
!!!$endif

!!!$if (coord==0) then
!!!$do jx=1,kx_num
!!!$do jy=1,ky_num
!!!$write(*,*) 'forw', jx, jy, fhat(jx,jy)
!!!$enddo
!!!$enddo
!!!$endif

!f(:,:) = 0._rprec

!do jx=1, kx_num
   !ii = 2*jx     ! imag index
   !ir = ii-1     ! real index
!do jy=1, ny
   !f(ir,jy) = real ( fhat(jx,jy) ,rprec )
   !f(ii,jy) = aimag( fhat(jx,jy) )
!enddo
!enddo


!!!$do jx=1, kx_num
!!!$   ii = 2*jx     ! imag index
!!!$   ir = ii-1     ! real index
!!!$do jy=2, 4
!!!$   jy_r = ny - jy + 2;   !! reflected y index
!!!$   f(ir,jy) = real ( fhat(jx,jy) ,rprec )
!!!$   f(ii,jy) = aimag( fhat(jx,jy) )
!!!$   f(ir,jy_r) = f(ir, jy)
!!!$   f(ii,jy_r) = -f(ii, jy)
!!!$enddo
!!!$enddo


!return
!end subroutine dft_direct_forw_2d

!!**********************************************************************
!subroutine dft_direct_back_2d(f)
!!**********************************************************************
!!
!! Computes 2d inverse DFT directly, no FFT.
!!
!use types,only:rprec
!use param,only:ld,nx,ny,nz,kx_num,pi,coord,L_x,nproc,L_y,kx_vec
!use fft,only:ky_vec
!implicit none

!integer :: jx, jy, i, j, ii, ir
!real(rprec), dimension(:, :), intent(inout) :: f
!integer :: ky_num

!complex(rprec), dimension(nx,ny) :: fhat
!complex(rprec) :: imag = (0.0_rprec, 1.0_rprec)
!complex(rprec) :: c1, c2

!real(rprec), dimension(kx_num,ny) :: pre

!pre(:,:) = 1._rprec
!pre(2:kx_num-1,:) = 2._rprec

!ky_num = ny/2+1
!fhat(:,:)  = ( 0._rprec, 0._rprec  )

!do jx=1,kx_num
   !ii = 2*jx     ! imag index
   !ir = ii-1     ! real index
!do jy=1,ny
   !fhat(jx,jy) = cmplx( f(ir,jy), f(ii,jy), rprec )
!enddo
!enddo

!!!$if (coord == 0) then
!!!$   write(*,*) 'dft back'
!!!$do jx=1,kx_num
!!!$do jy=1,ny
!!!$   write(*,*) jx,jy,fhat(jx,jy)
!!!$enddo
!!!$enddo
!!!$endif

!c1 = L_x/real(nx,rprec) * imag
!c2 = L_y/real(ny,rprec) * imag

!f(:,:) = 0._rprec

!do i=1,Nx
!do j=1,Ny
!do jx=1,nx/2+1
!do jy=1,ny
!f(i,j)=f(i,j)+pre(jx,jy)*fhat(jx,jy)*exp(c1*kx_vec(jx)*real(i-1,rprec) + &
                                         !c2*ky_vec(jy)*real(j-1,rprec));
!end do
!end do
!end do
!end do

!!f = f / real(nx*ny,rprec)

!return
!end subroutine dft_direct_back_2d

!!**********************************************************************
!subroutine dft_direct_forw_2d_new(f)
!!**********************************************************************
!!
!! Computes 2d DFT directly, no FFT.
!!
!use types,only:rprec
!use param,only:ld,nx,ny,nz,kx_num,kx_vec,pi,coord,L_x,nproc,L_y
!use fft !,only:ky_vec
!implicit none

!integer :: jx, jy, i, j, ii, ir, jy_r
!real(rprec), dimension(:, :), intent(inout) :: f
!integer :: ky_num
!real(rprec) :: const

!!complex(rprec), dimension(kx_num, ny/2+1) :: fhat
!complex(rprec), dimension(nx, ny) :: fhat
!complex(rprec), dimension(nx, ny) :: fhat2
!complex(rprec), dimension(nx, ny/2+1) :: fhaty
!complex(rprec) :: imag = (0.0_rprec, 1.0_rprec)
!complex(rprec) :: c1, c2

!ky_num = ny/2+1

!c1 = -L_x/real(nx,rprec) * imag
!c2 = -L_y/real(ny,rprec) * imag

!fhat(:,:)  = ( 0._rprec, 0._rprec  )
!fhat2(:,:)  = ( 0._rprec, 0._rprec  )
!fhaty(:,:)  = ( 0._rprec, 0._rprec  )

!!const = 1._rprec / ny
!!f = const * f  

!do jx=1,nx
  !call dfftw_execute_dft_r2c(forw_y, f(jx,:), fhaty(jx,:) )
!enddo

!if (coord == 0) then
!write(*,*) 'fhaty FORW ++++++++++++++++++++++++++++++++++++'
!do jx=1,kx_num
!do jy=1,ny
   !write(*,*) 'jx, jy, fhaty: ', jx,jy,fhaty(jx,jy)
!enddo
!enddo
!endif

!fhat2( :, 1:ny/2+1 ) = fhaty( :, 1:ny/2+1 )

!do jy=2,ny/2
   !jy_r = ny - jy + 2;   !! reflected y index
   !fhat2(:,jy_r) = conjg(fhaty(:,jy))
!enddo

!!if (coord == 0) then
!!write(*,*) 'fhat2 FORW ++++++++++++++++++++++++++++++++++++'
!!do jx=1,nx
!!do jy=1,ny
!!   write(*,*) 'jx, jy, fhat2: ', jx,jy,fhat2(jx,jy)
!!enddo
!!enddo
!!endif

!do jx=1,kx_num
!do i=1,nx
!fhat(jx,:)=fhat(jx,:)+fhat2(i,:)*exp(c1*kx_vec(jx)*real(i-1,rprec) );
!end do
!end do

!if (coord == 0) then
!write(*,*) 'fhat FORW ++++++++++++++++++++++++++++++++++++'
!do jx=1,kx_num
!do jy=1,ny
   !write(*,*) 'jx, jy, fhat: ', jx,jy,fhat(jx,jy)
!enddo
!enddo
!endif

!!fhat = fhat / real(nx*ny,rprec)

!f(:,:) = 0._rprec

!do jx=1, kx_num
   !ii = 2*jx     ! imag index
   !ir = ii-1     ! real index
!do jy=1, ny
   !f(ir,jy) = real ( fhat(jx,jy) ,rprec )
   !f(ii,jy) = aimag( fhat(jx,jy) )
!enddo
!enddo

!return
!end subroutine dft_direct_forw_2d_new

!!**********************************************************************
!subroutine dft_direct_back_2d_new(f)
!!**********************************************************************
!!
!! Computes 2d inverse DFT directly, no FFT.
!!
!use types,only:rprec
!use param,only:ld,nx,ny,nz,kx_num,pi,coord,L_x,nproc,L_y,kx_vec
!use fft !,only:ky_vec
!implicit none

!integer :: jx, jy, i, j, ii, ir
!real(rprec), dimension(:, :), intent(inout) :: f

!complex(rprec), dimension(nx,ny) :: fhat, fhat2
!complex(rprec), dimension(nx, ny/2+1) :: fhaty
!complex(rprec) :: imag = (0.0_rprec, 1.0_rprec)
!complex(rprec) :: c1, c2

!!!!real(rprec), dimension(kx_num) :: pre
!real(rprec), dimension(kx_num*2 - 2) :: kx_vec_long

!kx_vec_long(1:kx_num) = kx_vec(1:kx_num)

!do jx=2,nx/2
   !ii = nx - jx + 2;   !! reflected x index
   !kx_vec_long(ii) = kx_vec(jx) * (-1._rprec)
!enddo

!if (coord == 0) then
   !write(*,*) 'kx_vec_long: ', kx_vec_long
!endif

!!!pre(:) = 1._rprec
!!!pre(2:kx_num-1) = 2._rprec

!fhat(:,:)  = ( 0._rprec, 0._rprec  )
!fhat2(:,:)  = ( 0._rprec, 0._rprec  )
!fhaty(:,:)  = ( 0._rprec, 0._rprec  )

!do jx=1,kx_num
   !ii = 2*jx     ! imag index
   !ir = ii-1     ! real index
!do jy=1,ny
   !fhat(jx,jy) = cmplx( f(ir,jy), f(ii,jy), rprec )
!enddo
!enddo

!if (coord == 0) then
!write(*,*) 'fhat BACK ++++++++++++++++++++++++++++++++++++'
!do jx=1,kx_num
!do jy=1,ny
   !write(*,*) 'jx, jy, fhat: ', jx,jy,fhat(jx,jy)
!enddo
!enddo
!endif



!c1 = L_x/real(nx,rprec) * imag

!f(:,:) = 0._rprec

!fhat(6,1) = conjg(fhat(4,1))
!fhat(7,1) = conjg(fhat(3,1))
!fhat(8,1) = conjg(fhat(2,1))

!fhat(6,2) = conjg(fhat(4,8))
!fhat(7,2) = conjg(fhat(3,8))
!fhat(8,2) = conjg(fhat(2,8))
!fhat(6,3) = conjg(fhat(4,7))
!fhat(7,3) = conjg(fhat(3,7))
!fhat(8,3) = conjg(fhat(2,7))
!fhat(6,4) = conjg(fhat(4,6))
!fhat(7,4) = conjg(fhat(3,6))
!fhat(8,4) = conjg(fhat(2,6))
!fhat(6,5) = conjg(fhat(4,5))
!fhat(7,5) = conjg(fhat(3,5))
!fhat(8,5) = conjg(fhat(2,5))

!do i=1,nx
!do jx=1,8 !kx_num
!fhat2(i,:)=fhat2(i,:)+fhat(jx,:)*exp(c1*kx_vec_long(jx)*real(i-1,rprec) );
!end do
!end do

!if (coord == 0) then
!write(*,*) 'fhat2 BACK ++++++++++++++++++++++++++++++++++++'
!do jx=1,nx
!do jy=1,ny
   !write(*,*) 'jx, jy, fhat2: ', jx,jy,fhat2(jx,jy)
!enddo
!enddo
!endif

!do jx=1,nx
  !call dfftw_execute_dft_c2r(back_y, fhat2(jx,1:ny/2+1), f(jx,:) )
!enddo

!!!$do i=1,Nx
!!!$do j=1,Ny
!!!$do jx=1,nx/2+1
!!!$do jy=1,ny
!!!$f(i,j)=f(i,j)+pre(jx,jy)*fhat(jx,jy)*exp(c1*kx_vec(jx)*real(i-1,rprec) + &
!!!$                                         c2*ky_vec(jy)*real(j-1,rprec));
!!!$end do
!!!$end do
!!!$end do
!!!$end do

!!f = f / real(nx*ny,rprec)

!return
!end subroutine dft_direct_back_2d_new

!**********************************************************************
subroutine dft_direct_forw_2d_n(f)
!**********************************************************************
!
! Computes 2d DFT directly, no FFT in x-direction.
!
use types,only:rprec
use param,only:ld,nx,ny,nz,kx_num,kx_vec,pi,coord,L_x,nproc,L_y
use fft !,only:ky_vec
implicit none

integer :: jx, jy, i, j, ii, ir, jy_r, jx_s
real(rprec), dimension(:, :), intent(inout) :: f
real(rprec) :: const

complex(rprec), dimension(nx, ny) :: fhat
complex(rprec), dimension(nx, ny) :: fhat2
complex(rprec), dimension(nx, ny/2+1) :: fhaty

fhat(:,:)  = ( 0._rprec, 0._rprec  )
fhat2(:,:)  = ( 0._rprec, 0._rprec  )
fhaty(:,:)  = ( 0._rprec, 0._rprec  )

!!$if (coord == 0) then
!!$   print*, 'grg0 >>>>>>>>>>>>>>>'
!!$   do jy=1,ny
!!$      write(*,*) jy, f(1,jy)
!!$   enddo
!!$endif

!! 1d FFTs in the spanwise direction
do jx=1,nx
  call dfftw_execute_dft_r2c(forw_y, f(jx,:), fhaty(jx,:) )
enddo

!!$if (coord == 0) then
!!$   print*, 'grg1 >>>>>>>>>>>>>>>'
!!$   do jy=1,ny/2+1
!!$      write(*,*) jy, fhaty(1,jy)
!!$   enddo
!!$endif

!! copy to larger array
fhat2( :, 1:ny/2+1 ) = fhaty( :, 1:ny/2+1 )

!! due to symmetry, populate negative ky with conj of positive ky
do jy=2,ny/2
   jy_r = ny - jy + 2;   !! reflected y index
   fhat2(:,jy_r) = conjg(fhaty(:,jy))
enddo

!!$if (coord == 0) then
!!$   print*, 'grg2 >>>>>>>>>>>>>>>'
!!$   do jy=1,ny
!!$      write(*,*) jy, fhat2(1,jy)
!!$   enddo
!!$endif

!! now DFT in streamwise direction
do jx=1,kx_num
jx_s = kx_veci ( jx )  !! index/storage location
do i=1,nx
  fhat(jx_s,:)=fhat(jx_s,:)+fhat2(i,:)*expn(jx,i)   !exp(c1*kx_vec(jx)*real(i-1,rprec) );
end do
end do

!!$if (coord == 0) then
!!$   print*, 'grg3 >>>>>>>>>>>>>>>'
!!$   do jy=1,ny
!!$      write(*,*) jy, fhat(1,jy)
!!$   enddo
!!$endif

f(:,:) = 0._rprec

!! interleave the complex array into a real array
do jx=1, kx_num
   jx_s = kx_veci( jx )   !!---jb
   ii = 2*jx_s     ! imag index
   ir = ii-1       ! real index
do jy=1, ny
   f(ir,jy) = real ( fhat(jx_s,jy) ,rprec )
   f(ii,jy) = real(aimag( fhat(jx_s,jy) ),rprec)   !! or dimag?
enddo
enddo

!!$if (coord == 0) then
!!$   print*, 'grg4 >>>>>>>>>>>>>>>'
!!$   do jy=1,ny
!!$      write(*,*) jy, f(1,jy)
!!$   enddo
!!$endif

return
end subroutine dft_direct_forw_2d_n

!**********************************************************************
subroutine dft_direct_back_2d_n(f)
!**********************************************************************
!
! Computes 2d inverse DFT directly, no FFT in x-direction.
!
use types,only:rprec
use param,only:ld,nx,ny,nz,kx_num,pi,coord,L_x,nproc,L_y,kx_vec
use fft !,only:ky_vec
implicit none

integer :: jx, jy, i, j, ii, ir, jx_r, jy_r, jx_s, jx2
real(rprec), dimension(:, :), intent(inout) :: f

complex(rprec), dimension(nx,ny) :: fhat, fhat2
complex(rprec), dimension(nx, ny/2+1) :: fhaty

fhat(:,:)  = ( 0._rprec, 0._rprec  )
fhat2(:,:)  = ( 0._rprec, 0._rprec  )
fhaty(:,:)  = ( 0._rprec, 0._rprec  )

!!$if (coord == 0) then
!!$   print*, 'aaa4 >>>>>>>>>>>>>>>'
!!$   do jy=1,ny
!!$      write(*,*) jy, f(1,jy)
!!$   enddo
!!$endif

!! un-interleave the real array into a complex array
do jx=1,kx_num
   jx_s = kx_veci( jx )
   ii = 2*jx_s     ! imag index
   ir = ii-1     ! real index
do jy=1,ny
   fhat(jx_s,jy) = cmplx( f(ir,jy), f(ii,jy), rprec )
enddo
enddo

!!$if (coord == 0) then
!!$   print*, 'aaa3 >>>>>>>>>>>>>>>'
!!$   do jy=1,ny
!!$      write(*,*) jy, fhat(1,jy)
!!$   enddo
!!$endif

f(:,:) = 0._rprec

!! populate the negative kx with the conj of positive kx
do i= 2, nx/2     !!jb - I believe the +1 not needed here... test
  jx_r = nx - i + 2
  fhat(jx_r,1) = conjg( fhat(i,1) )
enddo

!!$if (coord == 0) then
!!$   print*, 'aaa2 >>>>>>>>>>>>>>>'
!!$   do jy=1,ny
!!$      write(*,*) jy, fhat(1,jy)
!!$   enddo
!!$endif

!! populate (-kx,ky) with conj of (kx,-ky)
do j = 2     , ny/2+1
  jy_r = ny - j + 2
do i = nx/2+2, nx
  jx_r = nx - i + 2
  fhat(i,j) = conjg(fhat(jx_r ,jy_r))
enddo
enddo

!!$if (coord == 0) then
!!$   print*, 'aaa1 >>>>>>>>>>>>>>>'
!!$   do jy=1,ny
!!$      write(*,*) jy, fhat(1,jy)
!!$   enddo
!!$endif

!! Nyquist mode ( kx=nx/2 ) is not taken into account
!!$do i=1,nx
!!$do jx=2,kx_num
!!$jx_s = kx_veci( jx )
!!$jx2 = nx - jx_s + 2
!!$fhat2(i,:)=fhat2(i,:)+fhat(jx_s,:) * expp(jx,i)  !exp( c1*kx_vec(jx)*real(i-1,rprec) )
!!$fhat2(i,:)=fhat2(i,:)+fhat(jx2,:) * expn(jx,i)    !exp(-c1*kx_vec(jx)*real(i-1,rprec) )
!!$end do
!!$end do

!!$do i=1,nx
!!$do jx=2,kx_num-1
!!$jx_s = kx_veci( jx )
!!$jx2 = nx - jx_s + 2
!!$fhat2(i,:)=fhat2(i,:)+fhat(jx_s,:) * expp(jx,i)  !exp( c1*kx_vec(jx)*real(i-1,rprec) )
!!$fhat2(i,:)=fhat2(i,:)+conjg(fhat(jx_s,:)) * expn(jx,i)    !exp(-c1*kx_vec(jx)*real(i-1,rprec) )
!!$end do
!!$end do

do i =1,nx
do jx=1,kx_num+kx_num-2
   jx_s = kx_veci_n( jx )
   fhat2(i,:) = fhat2(i,:) + fhat(jx_s,:) * expp_n(jx_s,i)
end do
end do


!!$if (coord == 0) then
!!$   print*, 'aaa0 >>>>>>>>>>>>>>>'
!!$   do jy=1,ny
!!$      write(*,*) jy, fhat2(1,jy)
!!$   enddo
!!$endif

!!$!! add in kx=0
!!$do i=1,nx
!!$  fhat2(i,:) = fhat2(i,:) + fhat(1,:)
!!$enddo

!!$if (coord == 0) then
!!$   print*, 'aaa-1 >>>>>>>>>>>>>>>'
!!$   do jy=1,ny
!!$      write(*,*) jy, fhat2(1,jy)
!!$   enddo
!!$endif

!do i=1,nx
!do jx=1,nx
!jx_s = kx_vec( jx )
!jx2 = nx - jx_s + 2
!fhat2(i,:)=fhat2(i,:)+fhat(jx,:)*exp(c1*kx_vec(jx)*real(i-1,rprec) );
!fhat2(i,:)=fhat2(i,:)+fhat(jx2,:) *exp(-c1*kx_vec(jx)*real(i-1,rprec) );
!end do
!end do

do jx=1,nx
  call dfftw_execute_dft_c2r(back_y, fhat2(jx,1:ny/2+1), f(jx,:) )
enddo

!!$if (coord == 0) then
!!$   print*, 'aaa-2 >>>>>>>>>>>>>>>'
!!$   do jy=1,ny
!!$      write(*,*) jy, f(1,jy)
!!$   enddo
!!$endif

return
end subroutine dft_direct_back_2d_n


!!!!!!!!!!!!!!**********!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!**********************************************************************
subroutine dft_direct_forw_2d_n_yonly(f)
!**********************************************************************
!
! Computes 2d DFT directly, no FFT in x-direction.
!
use types,only:rprec
use param,only:ld,nx,ny,nz,kx_num,kx_vec,pi,coord,L_x,nproc,L_y
use fft !,only:ky_vec
implicit none

integer :: jx, jy, i, j, ii, ir, jy_r, jx_s
real(rprec), dimension(:, :), intent(inout) :: f
real(rprec) :: const

complex(rprec), dimension(nx, ny) :: fhat
complex(rprec), dimension(nx, ny) :: fhat2
complex(rprec), dimension(nx, ny/2+1) :: fhaty

fhat(:,:)  = ( 0._rprec, 0._rprec  )
fhat2(:,:)  = ( 0._rprec, 0._rprec  )
fhaty(:,:)  = ( 0._rprec, 0._rprec  )

!! 1d FFTs in the spanwise direction
do jx=1,nx
  call dfftw_execute_dft_r2c(forw_y, f(jx,:), fhaty(jx,:) )
enddo

!! copy to larger array
fhat2( :, 1:ny/2+1 ) = fhaty( :, 1:ny/2+1 )

!! due to symmetry, populate negative ky with conj of positive ky
do jy=2,ny/2
   jy_r = ny - jy + 2;   !! reflected y index
   fhat2(:,jy_r) = conjg(fhaty(:,jy))
enddo

!!$!! now DFT in streamwise direction
!!$do jx=1,kx_num
!!$jx_s = kx_veci ( jx )  !! index/storage location
!!$do i=1,nx
!!$  fhat(jx_s,:)=fhat(jx_s,:)+fhat2(i,:)*expn(jx,i)   !exp(c1*kx_vec(jx)*real(i-1,rprec) );
!!$end do
!!$end do

!if (coord == 0) then
! write(*,*) 'FHAT >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
! do jx=1,nx
!   write(*,*) 'jx,fhat(jx,1): ', jx, fhat(jx,1)
! enddo
!endif

f(:,:) = 0._rprec

!! interleave the complex array into a real array
do jx=1, kx_num
   jx_s = kx_veci( jx )   !!---jb
   ii = 2*jx_s     ! imag index
   ir = ii-1       ! real index
do jy=1, ny
   f(ir,jy) = real ( fhat2(jx_s,jy) ,rprec )
   f(ii,jy) = aimag( fhat2(jx_s,jy) )   !! or dimag?
enddo
enddo

return
end subroutine dft_direct_forw_2d_n_yonly

!**********************************************************************
subroutine dft_direct_back_2d_n_yonly(f)
!**********************************************************************
!
! Computes 2d inverse DFT directly, no FFT in x-direction.
!
use types,only:rprec
use param,only:ld,nx,ny,nz,kx_num,pi,coord,L_x,nproc,L_y,kx_vec
use fft !,only:ky_vec
implicit none

integer :: jx, jy, i, j, ii, ir, jx_r, jy_r, jx_s, jx2
real(rprec), dimension(:, :), intent(inout) :: f

complex(rprec), dimension(nx,ny) :: fhat, fhat2
complex(rprec), dimension(nx, ny/2+1) :: fhaty

fhat(:,:)  = ( 0._rprec, 0._rprec  )
fhat2(:,:)  = ( 0._rprec, 0._rprec  )
fhaty(:,:)  = ( 0._rprec, 0._rprec  )

!! un-interleave the real array into a complex array
do jx=1,kx_num
   jx_s = kx_veci( jx )
   ii = 2*jx_s     ! imag index
   ir = ii-1     ! real index
do jy=1,ny
   fhat(jx_s,jy) = cmplx( f(ir,jy), f(ii,jy), rprec )
enddo
enddo

f(:,:) = 0._rprec

!!$!! populate the negative kx with the conj of positive kx
!!$do i= 2, nx/2     !!jb - I believe the +1 not needed here... test
!!$  jx_r = nx - i + 2
!!$  fhat(jx_r,1) = conjg( fhat(i,1) )
!!$enddo

!!$!! populate (-kx,ky) with conj of (kx,-ky)
!!$do j = 2     , ny/2+1
!!$  jy_r = ny - j + 2
!!$do i = nx/2+2, nx
!!$  jx_r = nx - i + 2
!!$  fhat(i,j) = conjg(fhat(jx_r ,jy_r))
!!$enddo
!!$enddo

!!$!! Nyquist mode ( kx=nx/2 ) is not taken into account
!!$do i=1,nx
!!$do jx=2,kx_num
!!$jx_s = kx_veci( jx )
!!$jx2 = nx - jx_s + 2
!!$fhat2(i,:)=fhat2(i,:)+fhat(jx_s,:) * expp(jx,i)  !exp( c1*kx_vec(jx)*real(i-1,rprec) )
!!$fhat2(i,:)=fhat2(i,:)+fhat(jx2,:) * expn(jx,i)    !exp(-c1*kx_vec(jx)*real(i-1,rprec) )
!!$end do
!!$end do

!!$!! add in kx=0
!!$do i=1,nx
!!$  fhat2(i,:) = fhat2(i,:) + fhat(1,:)
!!$enddo

!do i=1,nx
!do jx=1,nx
!jx_s = kx_vec( jx )
!jx2 = nx - jx_s + 2
!fhat2(i,:)=fhat2(i,:)+fhat(jx,:)*exp(c1*kx_vec(jx)*real(i-1,rprec) );
!fhat2(i,:)=fhat2(i,:)+fhat(jx2,:) *exp(-c1*kx_vec(jx)*real(i-1,rprec) );
!end do
!end do

do jx=1,nx
  call dfftw_execute_dft_c2r(back_y, fhat(jx,1:ny/2+1), f(jx,:) )
enddo

return
end subroutine dft_direct_back_2d_n_yonly

!**********************************************************************
subroutine dft_direct_forw_2d_n_yonlyC(f)
!**********************************************************************
!
! Computes 2d DFT directly, no FFT in x-direction.
!
use types,only:rprec
use param,only:ld,nx,ny,nz,kx_num,kx_vec,pi,coord,L_x,nproc,L_y
use fft !,only:ky_vec
implicit none

integer :: jx, jy, i, j, ii, ir, jy_r, jx_s
real(rprec), dimension(:, :), intent(inout) :: f
real(rprec) :: const

complex(rprec), dimension(nx, ny) :: fhat
!complex(rprec), dimension(nx, ny) :: fhat2
!complex(rprec), dimension(nx, ny/2+1) :: fhaty

const = 1._rprec / real(ny,rprec)

fhat(:,:)  = ( 0._rprec, 0._rprec  )
!fhat2(:,:)  = ( 0._rprec, 0._rprec  )
!fhaty(:,:)  = ( 0._rprec, 0._rprec  )

f(:,:) = const * f(:,:)

!! un-interleave the real array into a complex array
do jx=1,kx_num
   jx_s = kx_veci( jx )
   ii = 2*jx_s     ! imag index
   ir = ii-1     ! real index
do jy=1,ny
   fhat(jx_s,jy) = cmplx( f(ir,jy), f(ii,jy), rprec )
enddo
enddo

f(:,:) = 0._rprec

! remains in complex space
do jx=1,nx
  call dfftw_execute_dft(ycomp_forw, fhat(jx,1:ny), fhat(jx,1:ny) )
enddo

!! interleave the complex array into a real array
do jx=1, kx_num
   jx_s = kx_veci( jx )   !!---jb
   ii = 2*jx_s     ! imag index
   ir = ii-1       ! real index
do jy=1, ny
   f(ir,jy) = real ( fhat(jx_s,jy) ,rprec )
   f(ii,jy) = aimag( fhat(jx_s,jy) )   !! or dimag?
enddo
enddo

return
end subroutine dft_direct_forw_2d_n_yonlyC

!**********************************************************************
subroutine dft_direct_back_2d_n_yonlyC(f)
!**********************************************************************
!
! Computes 2d inverse DFT directly, no FFT in x-direction.
!
use types,only:rprec
use param,only:ld,nx,ny,nz,kx_num,pi,coord,L_x,nproc,L_y,kx_vec
use fft !,only:ky_vec
implicit none

integer :: jx, jy, i, j, ii, ir, jx_r, jy_r, jx_s, jx2
real(rprec), dimension(:, :), intent(inout) :: f

complex(rprec), dimension(nx,ny) :: fhat !, fhat2
!complex(rprec), dimension(nx, ny/2+1) :: fhaty

fhat(:,:)  = ( 0._rprec, 0._rprec  )
!fhat2(:,:)  = ( 0._rprec, 0._rprec  )
!fhaty(:,:)  = ( 0._rprec, 0._rprec  )

!! un-interleave the real array into a complex array
do jx=1,kx_num
   jx_s = kx_veci( jx )
   ii = 2*jx_s     ! imag index
   ir = ii-1     ! real index
do jy=1,ny
   fhat(jx_s,jy) = cmplx( f(ir,jy), f(ii,jy), rprec )
enddo
enddo

f(:,:) = 0._rprec

! remains in complex space
do jx=1,nx
  call dfftw_execute_dft(ycomp_back, fhat(jx,1:ny), fhat(jx,1:ny) )
enddo

!! interleave the complex array into a real array
do jx=1, kx_num
   jx_s = kx_veci( jx )   !!---jb
   ii = 2*jx_s     ! imag index
   ir = ii-1       ! real index
do jy=1, ny
   f(ir,jy) = real ( fhat(jx_s,jy) ,rprec )
   f(ii,jy) = aimag( fhat(jx_s,jy) )   !! or dimag?
enddo
enddo

return
end subroutine dft_direct_back_2d_n_yonlyC

! >>>>>>>>>>>>>>>> spectra_jb >>>>>>>>>>>>>>>>>>>

!**********************************************************************
function ky_spectra_calc_jb( f1, g1, rc_case ) result(fg_out)
!**********************************************************************
!
! Computes ky spectra
!
use types,only:rprec
use param,only:ld,nx,ny,nz,coord,nproc,fourier
use fft !,only:ky_vec
implicit none

real(rprec), dimension(:, :), intent(in) :: f1, g1
integer, intent(in) :: rc_case
integer :: jx, jy, i, j, ii, ir, jy_r, jx_s
real(rprec) :: const

real(rprec), dimension(nx, ny) :: f, g
complex(rprec), dimension(nx, ny/2+1) :: fhat, ghat, fghat
real(rprec), dimension(ld, ny) :: fg_out

f(:,:) = 0._rprec
g(:,:) = 0._rprec
fg_out(:,:) = 0._rprec

const = 1._rprec / real(ny,rprec)
f(:,:) = const * f1(1:nx,1:ny)
g(:,:) = const * g1(1:nx,1:ny)
fhat(:,:)  = ( 0._rprec, 0._rprec  )
ghat(:,:)  = ( 0._rprec, 0._rprec  )
fghat(:,:)  = ( 0._rprec, 0._rprec  )

! remains in complex space
do jx=1,nx
  call dfftw_execute_dft_r2c(forw_span_spectra, f(jx,1:ny), fhat(jx,1:ny/2+1) )
  call dfftw_execute_dft_r2c(forw_span_spectra, g(jx,1:ny), ghat(jx,1:ny/2+1) )
enddo

do jx=1,nx
do jy=1,ny/2+1
   fghat(jx,jy) = fhat(jx,jy) * conjg( ghat(jx,jy) )
enddo
enddo

! note: Fortran intrinsic function 'real' works fine here for GNU compiler, but this
! operation might require the instrinsic 'realpart' function with other compilers
select case (rc_case)
case(0)  !! f1 = g1, thus fg_out is purely real
   do jy=1,ny/2+1
      fg_out(1:nx,jy) = real( fghat(1:nx,jy) )
   enddo
case(1)  !! f1 /= g1, thus fg_out is NOT purely real so we take its magnitude
   do jy=1,ny/2+1
      fg_out(1:nx,jy) = real( abs(fghat(1:nx,jy)) )
   enddo
!!$   if (coord == 0) then
!!$      do jy=1,ny/2+1
!!$         write(*,*) 'fghat(1,jy): ', jy, fghat(1,jy),fg_out(1,jy)
!!$      enddo
!!$   endif
end select

return
end function ky_spectra_calc_jb

! <<<<<<<<<<<< spectra_jb <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


!**********************************************************************
subroutine dft_direct_forw_2d_n_yonlyC_big(f)
!**********************************************************************
!
! Computes 2d DFT directly, no FFT in x-direction.
!
use types,only:rprec
use param,only:ld,nx,ny2,nz,kx_num,kx_vec,pi,coord,L_x,nproc,L_y
use fft !,only:ky_vec
implicit none

integer :: jx, jy, i, j, ii, ir, jy_r, jx_s
real(rprec), dimension(:, :), intent(inout) :: f
real(rprec) :: const

complex(rprec), dimension(nx, ny2) :: fhat

fhat(:,:)  = ( 0._rprec, 0._rprec  )

const = 1._rprec / real(ny2,rprec)

f(:,:) = const * f(:,:)

!! un-interleave the real array into a complex array
do jx=1,kx_num
   jx_s = kx_veci( jx )
   ii = 2*jx_s     ! imag index
   ir = ii-1     ! real index
do jy=1,ny2
   fhat(jx_s,jy) = cmplx( f(ir,jy), f(ii,jy), rprec )
enddo
enddo

f(:,:) = 0._rprec

! remains in complex space
do jx=1,nx
  call dfftw_execute_dft(ycomp_forw_big, fhat(jx,1:ny2), fhat(jx,1:ny2) )
enddo

!! interleave the complex array into a real array
do jx=1, kx_num
   jx_s = kx_veci( jx )   !!---jb
   ii = 2*jx_s     ! imag index
   ir = ii-1       ! real index
do jy=1, ny2
   f(ir,jy) = real ( fhat(jx_s,jy) ,rprec )
   f(ii,jy) = aimag( fhat(jx_s,jy) )   !! or dimag?
enddo
enddo

return
end subroutine dft_direct_forw_2d_n_yonlyC_big

!**********************************************************************
subroutine dft_direct_back_2d_n_yonlyC_big(f)
!**********************************************************************
!
! Computes 2d inverse DFT directly, no FFT in x-direction.
!
use types,only:rprec
use param,only:ld,nx,ny2,nz,kx_num,pi,coord,L_x,nproc,L_y,kx_vec,nx2
use fft !,only:ky_vec
implicit none

integer :: jx, jy, i, j, ii, ir, jx_r, jy_r, jx_s, jx2
real(rprec), dimension(:, :), intent(inout) :: f

complex(rprec), dimension(nx2,ny2) :: fhat

fhat(:,:)  = ( 0._rprec, 0._rprec  )

!! un-interleave the real array into a complex array
do jx=1,kx_num
   jx_s = kx_veci( jx )
   ii = 2*jx_s     ! imag index
   ir = ii-1     ! real index
do jy=1,ny2
   fhat(jx_s,jy) = cmplx( f(ir,jy), f(ii,jy), rprec )
enddo
enddo

f(:,:) = 0._rprec

! remains in complex space
do jx=1,nx2
  call dfftw_execute_dft(ycomp_back_big, fhat(jx,1:ny2), fhat(jx,1:ny2) )
enddo

!! interleave the complex array into a real array
do jx=1, kx_num
   jx_s = kx_veci( jx )   !!---jb
   ii = 2*jx_s     ! imag index
   ir = ii-1       ! real index
do jy=1, ny2
   f(ir,jy) = real ( fhat(jx_s,jy) ,rprec )
   f(ii,jy) = aimag( fhat(jx_s,jy) )   !! or dimag?
enddo
enddo

return
end subroutine dft_direct_back_2d_n_yonlyC_big

!!**********************************************************************
!subroutine filt_da_direct(f,dfdx,dfdy, lbz)
!!**********************************************************************
!!
!! kills the oddball components in f, and calculates in horizontal derivatives
!!--supplies results for jz=$lbz:nz
!!--MPI: all these arrays are 0:nz
!!--MPI: on exit, u,dudx,dudy,v,dvdx,dvdy,w,dwdx,dwdy valid on jz=0:nz,
!!  except on bottom process (0 level set to BOGUS, starts at 1)
!!
!use types,only:rprec
!use param,only:ld,nx,ny,nz,coord
!use fft
!use emul_complex, only : OPERATOR(.MULI.)
!implicit none
!integer::jz,jy,jx

!integer, intent(in) :: lbz
!real(rprec), dimension(:, :, lbz:), intent(inout) :: f
!real(rprec), dimension(:, :, lbz:), intent(inout) :: dfdx, dfdy

!real(rprec) :: const

!const = 1._rprec/(nx*ny)

!! loop through horizontal slices
!do jz=lbz,nz

  !f(:,:,jz)=const*f(:,:,jz)  

  !!$if ($FFTW3)
  !call dft_direct_forw_2d( f(:,:,jz) )
  !!$else
  !!write(*,*) 'WARNING - no subroutine available'
  !!$endif

  !! Zero padded region and Nyquist frequency
  !f(ld-1:ld,:,jz) = 0._rprec
  !f(:,ny/2+1,jz) = 0._rprec   

  !!  Compute in-plane derivatives
  !dfdx(:,:,jz) = f(:,:,jz) .MULI. kx
  !dfdy(:,:,jz) = f(:,:,jz) .MULI. ky

  !! the oddballs for derivatives should already be dead, since they are for f
  !! inverse transform 
  !!$if ($FFTW3)
  !call dft_direct_back_2d( f(:,:,jz) )
  !call dft_direct_back_2d( dfdx(:,:,jz) )
  !call dft_direct_back_2d( dfdy(:,:,jz) )
  !!$else
  !!write(*,*) 'WARNING - no subroutine available'
  !!$endif  

!end do

!return
!end subroutine filt_da_direct

!!$!**********************************************************************
!!$subroutine filt_da_direct_n(f,dfdx,dfdy, lbz)
!!$!**********************************************************************
!!$!
!!$! kills the oddball components in f, and calculates in horizontal derivatives
!!$!--supplies results for jz=$lbz:nz
!!$!--MPI: all these arrays are 0:nz
!!$!--MPI: on exit, u,dudx,dudy,v,dvdx,dvdy,w,dwdx,dwdy valid on jz=0:nz,
!!$!  except on bottom process (0 level set to BOGUS, starts at 1)
!!$!
!!$use types,only:rprec
!!$use param,only:ld,nx,ny,nz,coord
!!$use fft
!!$use emul_complex, only : OPERATOR(.MULI.)
!!$implicit none
!!$integer::jz,jy,jx
!!$
!!$integer, intent(in) :: lbz
!!$real(rprec), dimension(:, :, lbz:), intent(inout) :: f
!!$real(rprec), dimension(:, :, lbz:), intent(inout) :: dfdx, dfdy
!!$
!!$real(rprec) :: const
!!$
!!$const = 1._rprec/(nx*ny)
!!$
!!$! loop through horizontal slices
!!$do jz=lbz,nz
!!$
!!$  f(:,:,jz)=const*f(:,:,jz)  
!!$
!!$  !$if ($FFTW3)
!!$  call dft_direct_forw_2d_n( f(:,:,jz) )
!!$  !call dft_direct_forw_2d_n_yonlyC( f(:,:,jz) )
!!$  !$else
!!$  !write(*,*) 'WARNING - no subroutine available'
!!$  !$endif
!!$
!!$  ! Zero padded region and Nyquist frequency
!!$  f(ld-1:ld,:,jz) = 0._rprec
!!$  f(:,ny/2+1,jz) = 0._rprec   
!!$
!!$  !  Compute in-plane derivatives
!!$  dfdx(:,:,jz) = f(:,:,jz) .MULI. kx
!!$  dfdy(:,:,jz) = f(:,:,jz) .MULI. ky
!!$
!!$  ! the oddballs for derivatives should already be dead, since they are for f
!!$  ! inverse transform 
!!$  !$if ($FFTW3)
!!$  call dft_direct_back_2d_n( f(:,:,jz) )
!!$  call dft_direct_back_2d_n( dfdx(:,:,jz) )
!!$  call dft_direct_back_2d_n( dfdy(:,:,jz) )
!!$  !call dft_direct_back_2d_n_yonlyC( f(:,:,jz) )
!!$  !call dft_direct_back_2d_n_yonlyC( dfdx(:,:,jz) )
!!$  !call dft_direct_back_2d_n_yonlyC( dfdy(:,:,jz) )
!!$  !$else
!!$  !write(*,*) 'WARNING - no subroutine available'
!!$  !$endif  
!!$
!!$end do
!!$
!!$return
!!$end subroutine filt_da_direct_n

!**********************************************************************
subroutine filt_da_kxspace(f, dfdx, dfdy, lbz)
!**********************************************************************
!
! kills the oddball components in f, and calculates in horizontal derivatives
!--supplies results for jz=$lbz:nz
!--MPI: all these arrays are 0:nz
!--MPI: on exit, u,dudx,dudy,v,dvdx,dvdy,w,dwdx,dwdy valid on jz=0:nz,
!  except on bottom process (0 level set to BOGUS, starts at 1)
!
use types, only: rprec
use param, only: ld, nz
use fft
use emul_complex, only : OPERATOR(.MULI.)
implicit none

integer :: jz
integer, intent(in) :: lbz
real(rprec), dimension(:, :, lbz:), intent(inout) :: f
real(rprec), dimension(:, :, lbz:), intent(inout) :: dfdx, dfdy

! loop through horizontal slices
do jz=lbz,nz

  ! Zero padded region and Nyquist frequency
  f(ld-1:ld,:,jz) = 0._rprec
  f(:,ny/2+1,jz) = 0._rprec   

  !  Compute in-plane derivatives
  dfdx(:,:,jz) = f(:,:,jz) .MULI. kx
  dfdy(:,:,jz) = f(:,:,jz) .MULI. ky

end do

return
end subroutine filt_da_kxspace

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!**********************************************************************
subroutine dft_direct_forw_2d_n_big(f)
!**********************************************************************
!
! Computes 2d DFT directly, no FFT.
!
use types,only:rprec
use param,only:ld,nx,ny,nz,kx_num,kx_vec,pi,coord,L_x,nproc,L_y,nx2,ny2
use fft !,only:ky_vec
implicit none

integer :: jx, jy, i, j, ii, ir, jy_r, jx_s
real(rprec), dimension(:, :), intent(inout) :: f
real(rprec) :: const

complex(rprec), dimension(nx2, ny2) :: fhat
complex(rprec), dimension(nx2, ny2) :: fhat2
complex(rprec), dimension(nx2, ny2/2+1) :: fhaty

fhat(:,:)  = ( 0._rprec, 0._rprec  )
fhat2(:,:)  = ( 0._rprec, 0._rprec  )
fhaty(:,:)  = ( 0._rprec, 0._rprec  )

!! 1d FFTs in the spanwise direction
do jx=1,nx2
  call dfftw_execute_dft_r2c(forw_y_big, f(jx,:), fhaty(jx,:) )
enddo

!! copy to larger array
fhat2( :, 1:ny2/2+1 ) = fhaty( :, 1:ny2/2+1 )

!! due to symmetry
do jy=2,ny2/2
   jy_r = ny2 - jy + 2;   !! reflected y index
   fhat2(:,jy_r) = conjg(fhaty(:,jy))
enddo

!! now DFT in streamwise direction
do jx=1,kx_num
jx_s = kx_veci ( jx )  !! index/storage location
do i=1,nx2
  fhat(jx_s,:)=fhat(jx_s,:)+fhat2(i,:)*expn_big(jx,i)!exp(c1*kx_vec(jx)*real(i-1,rprec));
end do
end do

!if (coord == 0) then
! write(*,*) 'FHAT >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
! do jx=1,nx
!   write(*,*) 'jx,fhat(jx,1): ', jx, fhat(jx,1)
! enddo
!endif

f(:,:) = 0._rprec

do jx=1, kx_num
   jx_s = kx_veci( jx )   !!---jb
   ii = 2*jx_s     ! imag index
   ir = ii-1     ! real index
do jy=1, ny2
   f(ir,jy) = real ( fhat(jx_s,jy) ,rprec )
   f(ii,jy) = aimag( fhat(jx_s,jy) )
enddo
enddo

return
end subroutine dft_direct_forw_2d_n_big

!**********************************************************************
subroutine dft_direct_back_2d_n_big(f)
!**********************************************************************
!
! Computes 2d inverse DFT directly, no FFT.
!
use types,only:rprec
use param,only:ld,nx,ny,nz,kx_num,pi,coord,L_x,nproc,L_y,kx_vec,nx2,ny2
use fft !,only:ky_vec
implicit none

integer :: jx, jy, i, j, ii, ir, jx_r, jy_r, jx_s, jx2
real(rprec), dimension(:, :), intent(inout) :: f

complex(rprec), dimension(nx2,ny2) :: fhat, fhat2
complex(rprec), dimension(nx2, ny2/2+1) :: fhaty

fhat(:,:)  = ( 0._rprec, 0._rprec  )
fhat2(:,:)  = ( 0._rprec, 0._rprec  )
fhaty(:,:)  = ( 0._rprec, 0._rprec  )

do jx=1,kx_num
   jx_s = kx_veci( jx )
   ii = 2*jx_s     ! imag index
   ir = ii-1     ! real index
do jy=1,ny2
   fhat(jx_s,jy) = cmplx( f(ir,jy), f(ii,jy), rprec )
enddo
enddo

f(:,:) = 0._rprec

do i= 2, nx2/2     !!jb - I believe the +1 not needed here... test
  jx_r = nx2 - i + 2
  fhat(jx_r,1) = conjg( fhat(i,1) )
enddo

do j = 2     , ny2/2+1
  jy_r = ny2 - j + 2    !!!!
do i = nx2/2+2, nx2
  jx_r = nx2 - i + 2
  fhat(i,j) = conjg(fhat(jx_r ,jy_r))
enddo
enddo

!!$!! Nyquist mode ( kx=nx/2 ) is not taken into account
!!$do i=1,nx2
!!$do jx=2,kx_num
!!$jx_s = kx_veci( jx )
!!$jx2 = nx2 - jx_s + 2
!!$fhat2(i,:)=fhat2(i,:)+fhat(jx_s,:) * expp_big(jx,i) !exp( c1*kx_vec(jx)*real(i-1,rprec) )
!!$fhat2(i,:)=fhat2(i,:)+fhat(jx2,:) * expn_big(jx,i)  !exp(-c1*kx_vec(jx)*real(i-1,rprec) )
!!$end do
!!$end do
!!$
!!$!! add in kx=0
!!$do i=1,nx2
!!$  fhat2(i,:) = fhat2(i,:) + fhat(1,:)
!!$enddo

!! Nyquist mode ( kx=nx/2 ) is not taken into account
do i =1,nx2
do jx=1,kx_num+kx_num-2
   jx_s = kx_veci_n_big( jx )
   fhat2(i,:) = fhat2(i,:) + fhat(jx_s,:) * expp_big_n(jx_s,i)
end do
end do

!do i=1,nx
!do jx=1,nx
!jx_s = kx_vec( jx )
!jx2 = nx - jx_s + 2
!fhat2(i,:)=fhat2(i,:)+fhat(jx,:)*exp(c1*kx_vec(jx)*real(i-1,rprec) );
!fhat2(i,:)=fhat2(i,:)+fhat(jx2,:) *exp(-c1*kx_vec(jx)*real(i-1,rprec) );
!end do
!end do

do jx=1,nx2    !!!!
  call dfftw_execute_dft_c2r(back_y_big, fhat2(jx,1:ny2/2+1), f(jx,:) )  !!!!
enddo

return
end subroutine dft_direct_back_2d_n_big

!**********************************************************************
subroutine phys2wave(f)   !!jb
!**********************************************************************
!
! Transform field from physical space to wavenumber space.
! (x,y) --> (kx,ky)
!
use types,only:rprec
use param,only:nx,ny,nz,lbz
use fft
implicit none

integer :: jz
real(rprec), dimension(:, :, lbz:), intent(inout) :: f
real(rprec) :: const

const = 1._rprec/(nx*ny)

do jz=lbz, nz
   f(:,:,jz) = const * f(:,:,jz)  !! normalization
   call dfftw_execute_dft_r2c(forw, f(:,:,jz), f(:,:,jz) )
   !call dft_direct_forw_2d_n( f(:,:,jz) )
enddo

return
end subroutine phys2wave

!**********************************************************************
subroutine wave2phys(f)   !!jb
!**********************************************************************
!
! Transform field from wavenumber space to physical space.
! (kx, ky) --> (x,y)
!
use types,only:rprec
use param,only:nz,lbz
use fft
implicit none

integer :: jz
real(rprec), dimension(:, :, lbz:), intent(inout) :: f

do jz=lbz, nz
  call dfftw_execute_dft_c2r(back, f(:,:,jz), f(:,:,jz) )
  !call dft_direct_back_2d_n( f(:,:,jz) )
enddo

return
end subroutine wave2phys

!**********************************************************************
subroutine phys2wave_pr(f)   !!jb  !! for pressure gradients only
!**********************************************************************
!
! Transform field from physical space to wavenumber space.
! (x,y) --> (kx,ky)
!
! for pressure gradients only due to different array size in z-direction
use types,only:rprec
use param,only:nx,ny,nz
use fft
implicit none

integer :: jz
real(rprec), dimension(:, :, :), intent(inout) :: f
real(rprec) :: const

const = 1._rprec/(nx*ny)

do jz = 1, nz
   f(:,:,jz) = const * f(:,:,jz)  !! normalization
   call dfftw_execute_dft_r2c(forw, f(:,:,jz), f(:,:,jz) )
   !call dft_direct_forw_2d_n( f(:,:,jz) )
enddo

return
end subroutine phys2wave_pr

!**********************************************************************
subroutine wave2phys_pr(f)   !!jb
!**********************************************************************
!
! Transform field from wavenumber space to physical space.
! (kx, ky) --> (x,y)
!
! for pressure gradients only due to different array size in z-dir
use types,only:rprec
use param,only:nz
use fft
implicit none

integer :: jz
real(rprec), dimension(:, :, :), intent(inout) :: f

do jz=1, nz
  call dfftw_execute_dft_c2r(back, f(:,:,jz), f(:,:,jz) )
  !call dft_direct_back_2d_n( f(:,:,jz) )
enddo

return
end subroutine wave2phys_pr

!**********************************************************************
subroutine wave2physF(uhat, u)   !!jb
!**********************************************************************
!
! Transform field from wavenumber space to physical space.
! (kx, ky) --> (x,y)
!
use types,only:rprec
use param,only:nx,ny,nz,lbz,nxp
use fft
implicit none

integer :: jx, jy, jz, ii, ir, ii_, ir_
real(rprec), dimension(ld, ny, lbz:nz), intent(in) :: uhat
real(rprec), dimension(nxp+2, ny, lbz:nz), intent(out) :: u

u(:,:,:) = 0._rprec

!write(*,*) 'yt1: ', u(1:2,1,1)

do jx=1, kx_num
   ii_ = 2*jx
   ir_ = ii_ - 1

   ii = 2 * kxs_in(jx) + 2
   ir = ii-1

   u(ir:ii, :, :) = uhat(ir_:ii_, :, :)
enddo
!write(*,*) 'yt2: ', u(1:2,1,1)

!!$print*, 'mid >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
!!$do jx=1,nxp+2
!!$   do jy=1,ny
!!$      write(*,*) jx, jy, u(jx,jy,3)
!!$   enddo
!!$enddo


do jz=lbz, nz
  call dfftw_execute_dft_c2r(back_fourier, u(1:nxp+2,1:ny,jz), u(1:nxp+2,1:ny,jz) )
enddo

!write(*,*) 'yt3: ', u(1:2,1,1)

return
end subroutine wave2physF

!**********************************************************************
subroutine phys2waveF(u, uhat)   !!jb
!**********************************************************************
!
! Transform field from physical space to wavenumber space.
! (x,y) --> (kx, ky)
!
use types,only:rprec
use param,only:nz,lbz,nxp
use fft
implicit none

integer :: jx, jz, ii, ir, ii_, ir_
real(rprec), dimension(ld, ny, lbz:nz), intent(out) :: uhat
real(rprec), dimension(nxp+2, ny, lbz:nz), intent(inout) :: u  !!should be 'in' only, but using 'inout' for convenvience since we modify it

real(rprec) :: const
const = 1._rprec/(nxp*ny)

uhat = 0._rprec

do jz=lbz, nz
   u(:,:,jz) = const * u(:,:,jz)  !! normalization
   call dfftw_execute_dft_r2c(forw_fourier, u(:,:,jz), u(:,:,jz) )
enddo

do jx=1, kx_num
   ii_ = 2*jx
   ir_ = ii_ - 1

   ii = 2 * kxs_in(jx) + 2
   ir = ii-1

   uhat(ir_:ii_, :, :) = u(ir:ii, :, :)
enddo

return
end subroutine phys2waveF

!**********************************************************************
subroutine wave2physF_pr(uhat, u)   !!jb
!**********************************************************************
!
! Transform field from wavenumber space to physical space.
! (kx, ky) --> (x,y)
!
use types,only:rprec
use param,only:nx,ny,nz,nxp
use fft
implicit none

integer :: jx, jy, jz, ii, ir, ii_, ir_
real(rprec), dimension(ld, ny, nz), intent(in) :: uhat
real(rprec), dimension(nxp+2, ny, nz), intent(out) :: u

u(:,:,:) = 0._rprec

do jx=1, kx_num
   ii_ = 2*jx
   ir_ = ii_ - 1

   ii = 2 * kxs_in(jx) + 2
   ir = ii-1

   u(ir:ii, :, :) = uhat(ir_:ii_, :, :)
enddo

!!$print*, 'mid >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
!!$do jx=1,nxp+2
!!$   do jy=1,ny
!!$      write(*,*) jx, jy, u(jx,jy,3)
!!$   enddo
!!$enddo


do jz=1, nz
  call dfftw_execute_dft_c2r(back_fourier, u(:,:,jz), u(:,:,jz) )
enddo


return
end subroutine wave2physF_pr

!**********************************************************************
subroutine phys2waveF_pr(u, uhat)   !!jb
!**********************************************************************
!
! Transform field from physical space to wavenumber space.
! (x,y) --> (kx, ky)
!
use types,only:rprec
use param,only:nz,nxp
use fft
implicit none

integer :: jx, jz, ii, ir, ii_, ir_
real(rprec), dimension(ld, ny, nz), intent(out) :: uhat
real(rprec), dimension(nxp+2, ny, nz), intent(inout) :: u  !!should be 'in' only, but using 'inout' for convenvience since we modify it

real(rprec) :: const
const = 1._rprec/(nxp*ny)

uhat = 0._rprec

do jz=1, nz
   u(:,:,jz) = const * u(:,:,jz)  !! normalization
   call dfftw_execute_dft_r2c(forw_fourier, u(:,:,jz), u(:,:,jz) )
enddo

do jx=1, kx_num
   ii_ = 2*jx
   ir_ = ii_ - 1

   ii = 2 * kxs_in(jx) + 2
   ir = ii-1

   uhat(ir_:ii_, :, :) = u(ir:ii, :, :)
enddo

return
end subroutine phys2waveF_pr

!**********************************************************************
function convolve(f,g) result(out)              
!**********************************************************************
!  This function computes the convolution of f and g, which are both
!  in Fourier space. The f and g arrays contain complex numbers stored as 
!  interleaved real arrays.

use types, only: rprec
use param, only: ld,nx,ny,nz,kx_num,coord,nx2
use fft
use functions, only: interleave_r2c, interleave_c2r
implicit none

real(rprec), dimension(:,:), intent(in) :: f, g
!real(rprec), dimension(ld_big,ny2) :: out
real(rprec), dimension(ld,ny2) :: out

!!complex(rprec), dimension(nx+nx-1,ny2) :: fc, gc, outc
complex(rprec), dimension(nx2+nx2-1,ny2) :: fc, gc, outc
!complex(rprec), dimension(nx+nx-1,ny2) :: fc, gc, outc
integer :: jx,jy,jx_s,ii,ir,k,i,j,st,jx1,jx2,s

fc(:,:)  = ( 0._rprec, 0._rprec  )
gc(:,:)  = ( 0._rprec, 0._rprec  )
outc(:,:)  = ( 0._rprec, 0._rprec  )
out(:,:)  = 0._rprec

!! f, g are in kx space, structured as real arrays
!! now re-structure as complex arrays
fc(1:nx2,:) = interleave_r2c( f(:,:) )
gc(1:nx2,:) = interleave_r2c( g(:,:) )
!fc(1:nx,:) = interleave_r2c( f(:,:) )
!gc(1:nx,:) = interleave_r2c( g(:,:) )

!!$if (coord==0) then
!!$print*, 'fc, gc: '
!!$   do jx=1,nx2+nx2-1
!!$      write(*,*) jx, fc(jx,1), gc(jx,1)
!!$   enddo
!!$endif

!! get the -kx modes from the +kx modes
do jx=2,nx2/2
   fc(nx2-jx+2,:) = conjg( fc(jx,:) )
   gc(nx2-jx+2,:) = conjg( gc(jx,:) )
enddo
!!$do jx=2,nx/2
!!$   fc(nx-jx+2,:) = conjg( fc(jx,:) )
!!$   gc(nx-jx+2,:) = conjg( gc(jx,:) )
!!$enddo

!!$if (coord==0) then
!!$print*, 'fc, gc: '
!!$   do jx=1,nx2+nx2-1
!!$      write(*,*) jx, fc(jx,1), gc(jx,1)
!!$   enddo
!!$endif

do jx1 = 1, (nx2 + nx2 - 1)
do jx2 = 1, nx2
   s = jx1 - jx2 + 1
   if(s .gt. 0 .and. s .le. nx2) then
      outc(jx1,:) = outc(jx1,:) + fc(jx2,:) * gc(s,:)
   endif
enddo
enddo
!!$do jx1 = 1, (nx + nx - 1)
!!$do jx2 = 1, nx
!!$   s = jx1 - jx2 + 1
!!$   if(s .gt. 0 .and. s .le. nx2) then
!!$      outc(jx1,:) = outc(jx1,:) + fc(jx2,:) * gc(s,:)
!!$   endif
!!$enddo
!!$enddo

!!$if (coord==0) then
!!$print*, 'outc: '
!!$   do jx=1,nx2+nx2-1
!!$      write(*,*) jx, outc(jx,1)
!!$   enddo
!!$endif

outc(1:nx2-1,:) = outc(1:nx2-1,:) + outc(nx2+1:nx2+nx2-1,:)
!outc(1:nx-1,:) = outc(1:nx-1,:) + outc(nx+1:nx+nx-1,:)

!! re-structure the complex array into a real array
out(:,:) = interleave_c2r( outc(1:nx,:) )
!out(9:10,:) = 0._rprec

return

end function convolve

!**********************************************************************
function convolve2(f,g) result(out)              
!**********************************************************************
!  This function computes the convolution of f and g, which are both
!  in Fourier space. The f and g arrays contain complex numbers stored as 
!  interleaved real arrays.

use types, only: rprec
use param, only: ld,nx,ny,nz,kx_num,coord,nx2
use fft
use functions, only: interleave_r2c, interleave_c2r
implicit none

real(rprec), dimension(:,:), intent(in) :: f, g
real(rprec), dimension(ld,ny) :: out

complex(rprec), dimension(nx+nx-1,ny) :: fc, gc, outc
integer :: jx,jy,jx_s,ii,ir,k,i,j,st,jx1,jx2,s

fc(:,:)  = ( 0._rprec, 0._rprec  )
gc(:,:)  = ( 0._rprec, 0._rprec  )
outc(:,:)  = ( 0._rprec, 0._rprec  )
out(:,:)  = 0._rprec

!! f, g are in kx space, structured as real arrays
!! now re-structure as complex arrays
fc(1:nx,:) = interleave_r2c( f(:,:) )
gc(1:nx,:) = interleave_r2c( g(:,:) )

do jx=2,nx/2
   fc(nx-jx+2,:) = conjg( fc(jx,:) )
   gc(nx-jx+2,:) = conjg( gc(jx,:) )
enddo

do jx1 = 1, (nx + nx - 1)
do jx2 = 1, nx
   s = jx1 - jx2 + 1
   if(s .gt. 0 .and. s .le. nx2) then
      outc(jx1,:) = outc(jx1,:) + fc(jx2,:) * gc(s,:)
   endif
enddo
enddo

outc(1:nx-1,:) = outc(1:nx-1,:) + outc(nx+1:nx+nx-1,:)

!! re-structure the complex array into a real array
out(:,:) = interleave_c2r( outc(1:nx,:) )

return

end function convolve2



!! this version (below) works
!**********************************************************************
function convolve_rnl(f,g) result(out)              
!**********************************************************************
!  This function computes the convolution of f and g, which are both
!  in Fourier space. The f and g arrays contain complex numbers stored as 
!  interleaved real arrays.

use types, only: rprec
use param, only: ld,nx,ny,nz,kx_num,coord,nx2
use fft
use functions, only: interleave_r2c, interleave_c2r
implicit none

real(rprec), dimension(:,:), intent(in) :: f, g
!real(rprec), dimension(ld_big,ny2) :: out
real(rprec), allocatable, dimension(:,:) :: out

!complex(rprec), dimension(nx2+nx2-1,ny2) :: fc, gc, outc
complex(rprec), allocatable, dimension(:,:) :: fc, gc, outc
integer :: jx,jy,jx_s,ii,ir,k,i,j,st,jx1,jx2,s

integer :: ld_, nx_, ny_

ld_ = size(f,1)    !! either ld or ld_big
nx_ = ld_ - 2      !! either nx or nx2
ny_ = size(f,2)    !! either ny or ny2

!!write(*,*) 'convolve: ', ld_, nx_, ny_

allocate(  out(ld_, ny_) )
allocate(   fc(nx_ + nx_ - 1, ny_) )
allocate(   gc(nx_ + nx_ - 1, ny_) )
allocate( outc(nx_ + nx_ - 1, ny_) )

fc(:,:)  = ( 0._rprec, 0._rprec  )
gc(:,:)  = ( 0._rprec, 0._rprec  )
outc(:,:)  = ( 0._rprec, 0._rprec  )
out(:,:)  = 0._rprec

!! f, g are in kx space, structured as real arrays
!! now re-structure as complex arrays
fc(1:nx,:) = interleave_r2c( f(:,:) )
gc(1:nx,:) = interleave_r2c( g(:,:) )

!! get the -kx modes from the +kx modes
do jx=2,nx_/2
   fc(nx_-jx+2,:) = conjg( fc(jx,:) )
   gc(nx_-jx+2,:) = conjg( gc(jx,:) )
enddo

!! UV part
outc(1,:) = outc(1,:) + fc(1,:) * gc(1,:)

!!Uv, vU parts
do jx1 = 2, (nx_+nx_-1)
!do jx2 = 1, nx2
   !s = jx1 - jx2 + 1
   !if(s .gt. 0) then
      outc(jx1,:) = outc(jx1,:) + fc(1,:) * gc(jx1,:)
      outc(jx1,:) = outc(jx1,:) + fc(jx1,:) * gc(1,:)
   !endif
!enddo
enddo

!! < uv > part
do jx1 = 2, nx_/2 !(nx2+nx2-1)
   !if( nx2 + 1 - jx1 .gt. 0 ) then
      outc(1,:) = outc(1,:) + fc(jx1,:) * conjg( gc(jx1,:) ) 
      outc(1,:) = outc(1,:) + conjg(fc(jx1,:))*gc(jx1,:)
   !endif
enddo

outc(1:nx_-1,:) = outc(1:nx_-1,:) + outc(nx_+1:nx_+nx_-1,:)

!! re-structure the complex array into a real array
out(:,:) = interleave_c2r( outc(1:nx,:) )

return

end function convolve_rnl


!!$!**********************************************************************
!!$function convolve_rnl(f,g) result(out)              
!!$!**********************************************************************
!!$!  This function computes the convolution of f and g, which are both
!!$!  in Fourier space. The f and g arrays contain complex numbers stored as 
!!$!  interleaved real arrays.
!!$
!!$use types, only: rprec
!!$use param, only: ld,nx,ny,nz,kx_num,coord,nx2,kx_veci_n_big
!!$use fft
!!$use functions, only: interleave_r2c, interleave_c2r
!!$implicit none
!!$
!!$real(rprec), dimension(:,:), intent(in) :: f, g
!!$real(rprec), dimension(ld_big,ny2) :: out
!!$
!!$complex(rprec), dimension(nx2+nx2-1,ny2) :: fc, gc, outc
!!$integer :: jx,jy,jx_s,ii,ir,k,i,j,st,jx1,jx2,s
!!$
!!$fc(:,:)  = ( 0._rprec, 0._rprec  )
!!$gc(:,:)  = ( 0._rprec, 0._rprec  )
!!$outc(:,:)  = ( 0._rprec, 0._rprec  )
!!$out(:,:)  = 0._rprec
!!$
!!$!! f, g are in kx space, structured as real arrays
!!$!! now re-structure as complex arrays
!!$fc(1:nx,:) = interleave_r2c( f(:,:) )
!!$gc(1:nx,:) = interleave_r2c( g(:,:) )
!!$
!!$!! get the -kx modes from the +kx modes
!!$do jx=2,nx2/2
!!$   fc(nx2-jx+2,:) = conjg( fc(jx,:) )
!!$   gc(nx2-jx+2,:) = conjg( gc(jx,:) )
!!$enddo
!!$
!!$!! UV part
!!$outc(1,:) = outc(1,:) + fc(1,:) * gc(1,:)
!!$
!!$!!Uv, vU parts
!!$do jx1 = 2, (nx2+nx2-1)
!!$!do jx2 = 1, nx2
!!$   !s = jx1 - jx2 + 1
!!$   !if(s .gt. 0) then
!!$      outc(jx1,:) = outc(jx1,:) + fc(1,:) * gc(jx1,:)
!!$      outc(jx1,:) = outc(jx1,:) + fc(jx1,:) * gc(1,:)
!!$   !endif
!!$!enddo
!!$enddo
!!$
!!$!! < uv > part
!!$do jx1 = 2, nx/2 !(nx2+nx2-1)
!!$   !if( nx2 + 1 - jx1 .gt. 0 ) then
!!$      outc(1,:) = outc(1,:) + fc(jx1,:) * conjg( gc(jx1,:) ) 
!!$      outc(1,:) = outc(1,:) + conjg(fc(jx1,:))*gc(jx1,:)
!!$   !endif
!!$enddo
!!$
!!$outc(1:nx2-1,:) = outc(1:nx2-1,:) + outc(nx2+1:nx2+nx2-1,:)
!!$
!!$!! re-structure the complex array into a real array
!!$out(:,:) = interleave_c2r( outc(1:nx,:) )
!!$
!!$return
!!$
!!$end function convolve_rnl



!!$subroutine padd_kxspace (u_big,u)
!!$use types,only:rprec
!!$use param,only:ld,ld_big,nx,ny,ny2
!!$implicit none
!!$
!!$!  u and u_big are interleaved as complex arrays
!!$real(rprec), dimension(ld,ny), intent(in) :: u
!!$real(rprec), dimension(ld,ny2), intent(out) :: u_big
!!$
!!$integer :: ny_h, j_s, j_big_s
!!$
!!$ny_h = ny/2
!!$
!!$! make sure the big array is zeroed!
!!$u_big(:,:) = 0._rprec
!!$
!!$! note: split access in an attempt to maintain locality
!!$u_big(:nx,:ny_h) = u(:nx,:ny_h)
!!$
!!$! Compute starting j locations for second transfer
!!$j_s = ny_h + 2
!!$j_big_s = ny2 - ny_h + 2
!!$
!!$u_big(:nx,j_big_s:ny2) = u(:nx,j_s:ny)
!!$end subroutine padd_kxspace

end module derivatives
