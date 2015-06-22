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
     ddx_direct, &
     ddy_only, &
     dft_direct_forw_2d, &
     dft_direct_back_2d, &
     filt_da_direct, &
     dft_direct_forw_2d_new, &
     dft_direct_back_2d_new, &
     dft_direct_forw_2d_n, &
     dft_direct_back_2d_n, &
     filt_da_direct_n, &
     ddx_n, &
     ddy_n

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

!**********************************************************************
subroutine ddx_n(f,dfdx,lbz)              
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
  !!call dfftw_execute_dft_r2c(forw, dfdx(:,:,jz),dfdx(:,:,jz))
  call dft_direct_forw_2d_n( dfdx(:,:,jz)  )      !!jb
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
  !!call dfftw_execute_dft_c2r(back, dfdx(:,:,jz), dfdx(:,:,jz))
  call dft_direct_back_2d_n( dfdx(:,:,jz)  )      !!jb
  $else
  call rfftwnd_f77_one_complex_to_real(back,dfdx(:,:,jz),fftwNull_p)
  $endif

enddo

return

end subroutine ddx_n

!**********************************************************************
subroutine ddy_n(f,dfdy, lbz)              
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
  !!call dfftw_execute_dft_r2c(forw, dfdy(:,:,jz), dfdy(:,:,jz))
  call dft_direct_forw_2d_n( dfdy(:,:,jz)  )      !!jb
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
  !!call dfftw_execute_dft_c2r(back, dfdy(:,:,jz), dfdy(:,:,jz))
  call dft_direct_back_2d_n( dfdy(:,:,jz)  )      !!jb
  $else
  call rfftwnd_f77_one_complex_to_real(back,dfdy(:,:,jz),fftwNull_p)     
  $endif
end do

return

end subroutine ddy_n

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
use param,only:nx,ny,nz,dz,coord,nproc,BOGUS
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
    do jx=1,nx
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
use param,only:nx,ny,nz,dz,coord,nproc,BOGUS
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

$if ($USE_RNL)
!**********************************************************************
subroutine ddx_direct(f, dfdx, lbz)
!**********************************************************************
!
! Computes x derivative spectrally using direct summation instead of
! FFTs. For use with RNL.
!
use types,only:rprec
use param,only:ld,nx,ny,nz,kx_num,kx_vec,pi,coord,L_x,nproc
implicit none

integer :: jx, jy, jz, k
integer, intent(in) :: lbz
real(rprec), dimension(:, :, lbz:), intent(in) :: f
real(rprec), dimension(:, :, lbz:), intent(out) :: dfdx
real(rprec) :: const

real(rprec), dimension(kx_num) :: pre
complex(rprec), dimension(kx_num) :: fhat
complex(rprec) :: imag = (0.0, 1.0)
complex(rprec) :: const1, const2

pre(:) = 2._rprec
pre(1) = 1._rprec
const1 = -L_x/nx * imag
const2 = -1 * const1

!!$if (coord == 0) then
!!$   print*, 'pre: ', pre
!!$   print*, 'const1: ', const1
!!$   print*, 'const2: ', const2
!!$endif

dfdx(:,:,:) = 0._rprec

do jy=1,ny
do jz=lbz,nz

   fhat(:)  = ( 0._rprec, 0.0_rprec  )

   do  k=1,kx_num
      do jx=1,nx
         fhat(k) = fhat(k) + f(jx,jy,jz) * exp(const1 * kx_vec(k) * real(jx-1,rprec) )
      enddo
   enddo

!!$   if (coord == nproc-1 .and. jz==1 .and. jy == 1) then
!!$   print*, 'Before deriv: '
!!$   do k=1,kx_num
!!$      write(*,*) 'k, fhat(k): ', k, fhat(k)
!!$   enddo
!!$   endif

   do k=1,kx_num
      fhat(k) = imag * kx_vec(k) * fhat(k)
   enddo

!!$   if (coord == nproc-1 .and. jz==1 .and. jy == 1) then
!!$   print*, 'After deriv: '
!!$   do k=1,kx_num
!!$      write(*,*) 'k, fhat(k): ', k, fhat(k)
!!$   enddo
!!$   endif
   !!fhat(5) = (0._rprec, 0._rprec)
   do jx=1,nx
      do   k=1,kx_num
         dfdx(jx,jy,jz) = dfdx(jx,jy,jz) + pre(k)*fhat(k)*exp(const2 * kx_vec(k) * real(jx-1,rprec))
      enddo
   enddo
   
enddo
enddo

dfdx = dfdx / real(nx,rprec)

return
end subroutine ddx_direct

!**********************************************************************
subroutine ddy_only(f,dfdy, lbz)              
!**********************************************************************
!
!  This subroutine computes the partial derivative of f with respect to
!  y using spectral decomposition (using 1d FFT instead of 2d)
!  
use types,only:rprec
use param,only:ld,nx,ny,nz,coord
use fft
implicit none      

integer :: jx, jy, jz
integer, intent(in) :: lbz
real(rprec), dimension(:,:,lbz:), intent(in) :: f
real(rprec), dimension(:,:,lbz:), intent(inout) :: dfdy
real(rprec) :: const
complex(rprec) :: imag = (0.0, 1.0)
complex(rprec), dimension(ny/2+1) :: fhat

const = 1._rprec / ny

! Loop through horizontal slices
do jz=lbz,nz    
do jx=1,nx

  !  Use dfdy to hold f; since we are doing IN_PLACE FFT's this is required
  dfdy(jx,:,jz) = const * f(jx,:,jz)  
  $if ($FFTW3)

  call dfftw_execute_dft_r2c(forw_y, dfdy(jx,:,jz), fhat(:) )
  $else
  write(*,*) 'ERROR - no DFT available!'
  $endif

  fhat(ny/2+1) = (0._rprec, 0._rprec)

  if (coord==0) then
     if (jx==1 .and. jz==1) then
        write(*,*) 'y fft'
        write(*,*) fhat(:)
     endif
  endif

  do jy=1,ny/2 !+1
     fhat(jy) = fhat(jy) * ky_vec(jy) * imag
  enddo

  ! Perform inverse transform to get pseudospectral derivative
  $if ($FFTW3)
  call dfftw_execute_dft_c2r(back_y, fhat(:), dfdy(jx,:,jz))
  $else
  write(*,*) 'ERROR - no DFT available!'
  $endif

enddo
enddo

return

end subroutine ddy_only

!**********************************************************************
subroutine dft_direct_forw_2d(f)
!**********************************************************************
!
! Computes 2d DFT directly, no FFT.
!
use types,only:rprec
use param,only:ld,nx,ny,nz,kx_num,kx_vec,pi,coord,L_x,nproc,L_y
use fft,only:ky_vec
implicit none

integer :: jx, jy, i, j, ii, ir, jy_r
real(rprec), dimension(:, :), intent(inout) :: f
integer :: ky_num

!complex(rprec), dimension(kx_num, ny/2+1) :: fhat
complex(rprec), dimension(nx, ny) :: fhat
complex(rprec) :: imag = (0.0_rprec, 1.0_rprec)
complex(rprec) :: c1, c2

ky_num = ny/2+1

c1 = -L_x/real(nx,rprec) * imag
c2 = -L_y/real(ny,rprec) * imag

fhat(:,:)  = ( 0._rprec, 0._rprec  )

do jx=1,kx_num
do jy=1,ny
do i=1,nx
do j=1,ny
fhat(jx,jy)=fhat(jx,jy)+f(i,j)*exp(c1*kx_vec(jx)*real(i-1,rprec) + c2*ky_vec(jy)*real(j-1,rprec));
end do
end do
end do
end do

!fhat = fhat / real(nx*ny,rprec)

!!$if (coord == 0) then
!!$   write(*,*) 'dft forw'
!!$do jx=1,kx_num
!!$do jy=1,ny
!!$   write(*,*) jx,jy,fhat(jx,jy)
!!$enddo
!!$enddo
!!$endif

!!$if (coord==0) then
!!$do jx=1,kx_num
!!$do jy=1,ky_num
!!$write(*,*) 'forw', jx, jy, fhat(jx,jy)
!!$enddo
!!$enddo
!!$endif

f(:,:) = 0._rprec

do jx=1, kx_num
   ii = 2*jx     ! imag index
   ir = ii-1     ! real index
do jy=1, ny
   f(ir,jy) = real ( fhat(jx,jy) ,rprec )
   f(ii,jy) = aimag( fhat(jx,jy) )
enddo
enddo


!!$do jx=1, kx_num
!!$   ii = 2*jx     ! imag index
!!$   ir = ii-1     ! real index
!!$do jy=2, 4
!!$   jy_r = ny - jy + 2;   !! reflected y index
!!$   f(ir,jy) = real ( fhat(jx,jy) ,rprec )
!!$   f(ii,jy) = aimag( fhat(jx,jy) )
!!$   f(ir,jy_r) = f(ir, jy)
!!$   f(ii,jy_r) = -f(ii, jy)
!!$enddo
!!$enddo


return
end subroutine dft_direct_forw_2d

!**********************************************************************
subroutine dft_direct_back_2d(f)
!**********************************************************************
!
! Computes 2d inverse DFT directly, no FFT.
!
use types,only:rprec
use param,only:ld,nx,ny,nz,kx_num,pi,coord,L_x,nproc,L_y,kx_vec
use fft,only:ky_vec
implicit none

integer :: jx, jy, i, j, ii, ir
real(rprec), dimension(:, :), intent(inout) :: f
integer :: ky_num

complex(rprec), dimension(nx,ny) :: fhat
complex(rprec) :: imag = (0.0_rprec, 1.0_rprec)
complex(rprec) :: c1, c2

real(rprec), dimension(kx_num,ny) :: pre

pre(:,:) = 1._rprec
pre(2:kx_num-1,:) = 2._rprec

ky_num = ny/2+1
fhat(:,:)  = ( 0._rprec, 0._rprec  )

do jx=1,kx_num
   ii = 2*jx     ! imag index
   ir = ii-1     ! real index
do jy=1,ny
   fhat(jx,jy) = cmplx( f(ir,jy), f(ii,jy), rprec )
enddo
enddo

!!$if (coord == 0) then
!!$   write(*,*) 'dft back'
!!$do jx=1,kx_num
!!$do jy=1,ny
!!$   write(*,*) jx,jy,fhat(jx,jy)
!!$enddo
!!$enddo
!!$endif

c1 = L_x/real(nx,rprec) * imag
c2 = L_y/real(ny,rprec) * imag

f(:,:) = 0._rprec

do i=1,Nx
do j=1,Ny
do jx=1,nx/2+1
do jy=1,ny
f(i,j)=f(i,j)+pre(jx,jy)*fhat(jx,jy)*exp(c1*kx_vec(jx)*real(i-1,rprec) + &
                                         c2*ky_vec(jy)*real(j-1,rprec));
end do
end do
end do
end do

!f = f / real(nx*ny,rprec)

return
end subroutine dft_direct_back_2d

!**********************************************************************
subroutine dft_direct_forw_2d_new(f)
!**********************************************************************
!
! Computes 2d DFT directly, no FFT.
!
use types,only:rprec
use param,only:ld,nx,ny,nz,kx_num,kx_vec,pi,coord,L_x,nproc,L_y
use fft !,only:ky_vec
implicit none

integer :: jx, jy, i, j, ii, ir, jy_r
real(rprec), dimension(:, :), intent(inout) :: f
integer :: ky_num
real(rprec) :: const

!complex(rprec), dimension(kx_num, ny/2+1) :: fhat
complex(rprec), dimension(nx, ny) :: fhat
complex(rprec), dimension(nx, ny) :: fhat2
complex(rprec), dimension(nx, ny/2+1) :: fhaty
complex(rprec) :: imag = (0.0_rprec, 1.0_rprec)
complex(rprec) :: c1, c2

ky_num = ny/2+1

c1 = -L_x/real(nx,rprec) * imag
c2 = -L_y/real(ny,rprec) * imag

fhat(:,:)  = ( 0._rprec, 0._rprec  )
fhat2(:,:)  = ( 0._rprec, 0._rprec  )
fhaty(:,:)  = ( 0._rprec, 0._rprec  )

!const = 1._rprec / ny
!f = const * f  

do jx=1,nx
  call dfftw_execute_dft_r2c(forw_y, f(jx,:), fhaty(jx,:) )
enddo

if (coord == 0) then
write(*,*) 'fhaty FORW ++++++++++++++++++++++++++++++++++++'
do jx=1,kx_num
do jy=1,ny
   write(*,*) 'jx, jy, fhaty: ', jx,jy,fhaty(jx,jy)
enddo
enddo
endif

fhat2( :, 1:ny/2+1 ) = fhaty( :, 1:ny/2+1 )

do jy=2,ny/2
   jy_r = ny - jy + 2;   !! reflected y index
   fhat2(:,jy_r) = conjg(fhaty(:,jy))
enddo

!if (coord == 0) then
!write(*,*) 'fhat2 FORW ++++++++++++++++++++++++++++++++++++'
!do jx=1,nx
!do jy=1,ny
!   write(*,*) 'jx, jy, fhat2: ', jx,jy,fhat2(jx,jy)
!enddo
!enddo
!endif

do jx=1,kx_num
do i=1,nx
fhat(jx,:)=fhat(jx,:)+fhat2(i,:)*exp(c1*kx_vec(jx)*real(i-1,rprec) );
end do
end do

if (coord == 0) then
write(*,*) 'fhat FORW ++++++++++++++++++++++++++++++++++++'
do jx=1,kx_num
do jy=1,ny
   write(*,*) 'jx, jy, fhat: ', jx,jy,fhat(jx,jy)
enddo
enddo
endif

!fhat = fhat / real(nx*ny,rprec)

f(:,:) = 0._rprec

do jx=1, kx_num
   ii = 2*jx     ! imag index
   ir = ii-1     ! real index
do jy=1, ny
   f(ir,jy) = real ( fhat(jx,jy) ,rprec )
   f(ii,jy) = aimag( fhat(jx,jy) )
enddo
enddo

return
end subroutine dft_direct_forw_2d_new

!**********************************************************************
subroutine dft_direct_back_2d_new(f)
!**********************************************************************
!
! Computes 2d inverse DFT directly, no FFT.
!
use types,only:rprec
use param,only:ld,nx,ny,nz,kx_num,pi,coord,L_x,nproc,L_y,kx_vec
use fft !,only:ky_vec
implicit none

integer :: jx, jy, i, j, ii, ir
real(rprec), dimension(:, :), intent(inout) :: f

complex(rprec), dimension(nx,ny) :: fhat, fhat2
complex(rprec), dimension(nx, ny/2+1) :: fhaty
complex(rprec) :: imag = (0.0_rprec, 1.0_rprec)
complex(rprec) :: c1, c2

!!!real(rprec), dimension(kx_num) :: pre
real(rprec), dimension(kx_num*2 - 2) :: kx_vec_long

kx_vec_long(1:kx_num) = kx_vec(1:kx_num)

do jx=2,nx/2
   ii = nx - jx + 2;   !! reflected x index
   kx_vec_long(ii) = kx_vec(jx) * (-1._rprec)
enddo

if (coord == 0) then
   write(*,*) 'kx_vec_long: ', kx_vec_long
endif

!!pre(:) = 1._rprec
!!pre(2:kx_num-1) = 2._rprec

fhat(:,:)  = ( 0._rprec, 0._rprec  )
fhat2(:,:)  = ( 0._rprec, 0._rprec  )
fhaty(:,:)  = ( 0._rprec, 0._rprec  )

do jx=1,kx_num
   ii = 2*jx     ! imag index
   ir = ii-1     ! real index
do jy=1,ny
   fhat(jx,jy) = cmplx( f(ir,jy), f(ii,jy), rprec )
enddo
enddo

if (coord == 0) then
write(*,*) 'fhat BACK ++++++++++++++++++++++++++++++++++++'
do jx=1,kx_num
do jy=1,ny
   write(*,*) 'jx, jy, fhat: ', jx,jy,fhat(jx,jy)
enddo
enddo
endif



c1 = L_x/real(nx,rprec) * imag

f(:,:) = 0._rprec

fhat(6,1) = conjg(fhat(4,1))
fhat(7,1) = conjg(fhat(3,1))
fhat(8,1) = conjg(fhat(2,1))

fhat(6,2) = conjg(fhat(4,8))
fhat(7,2) = conjg(fhat(3,8))
fhat(8,2) = conjg(fhat(2,8))
fhat(6,3) = conjg(fhat(4,7))
fhat(7,3) = conjg(fhat(3,7))
fhat(8,3) = conjg(fhat(2,7))
fhat(6,4) = conjg(fhat(4,6))
fhat(7,4) = conjg(fhat(3,6))
fhat(8,4) = conjg(fhat(2,6))
fhat(6,5) = conjg(fhat(4,5))
fhat(7,5) = conjg(fhat(3,5))
fhat(8,5) = conjg(fhat(2,5))

do i=1,nx
do jx=1,8 !kx_num
fhat2(i,:)=fhat2(i,:)+fhat(jx,:)*exp(c1*kx_vec_long(jx)*real(i-1,rprec) );
end do
end do

if (coord == 0) then
write(*,*) 'fhat2 BACK ++++++++++++++++++++++++++++++++++++'
do jx=1,nx
do jy=1,ny
   write(*,*) 'jx, jy, fhat2: ', jx,jy,fhat2(jx,jy)
enddo
enddo
endif

do jx=1,nx
  call dfftw_execute_dft_c2r(back_y, fhat2(jx,1:ny/2+1), f(jx,:) )
enddo

!!$do i=1,Nx
!!$do j=1,Ny
!!$do jx=1,nx/2+1
!!$do jy=1,ny
!!$f(i,j)=f(i,j)+pre(jx,jy)*fhat(jx,jy)*exp(c1*kx_vec(jx)*real(i-1,rprec) + &
!!$                                         c2*ky_vec(jy)*real(j-1,rprec));
!!$end do
!!$end do
!!$end do
!!$end do

!f = f / real(nx*ny,rprec)

return
end subroutine dft_direct_back_2d_new

!**********************************************************************
subroutine dft_direct_forw_2d_n(f)
!**********************************************************************
!
! Computes 2d DFT directly, no FFT.
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

do jx=1,nx
  call dfftw_execute_dft_r2c(forw_y, f(jx,:), fhaty(jx,:) )
enddo

fhat2( :, 1:ny/2+1 ) = fhaty( :, 1:ny/2+1 )

do jy=2,ny/2
   jy_r = ny - jy + 2;   !! reflected y index
   fhat2(:,jy_r) = conjg(fhaty(:,jy))
enddo

do jx=1,kx_num
jx_s = kx_veci ( jx )  !! index/storage location
do i=1,nx
  fhat(jx_s,:)=fhat(jx_s,:)+fhat2(i,:)*expn(jx,i)   !exp(c1*kx_vec(jx)*real(i-1,rprec) );
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
do jy=1, ny
   f(ir,jy) = real ( fhat(jx_s,jy) ,rprec )
   f(ii,jy) = aimag( fhat(jx_s,jy) )
enddo
enddo

return
end subroutine dft_direct_forw_2d_n

!**********************************************************************
subroutine dft_direct_back_2d_n(f)
!**********************************************************************
!
! Computes 2d inverse DFT directly, no FFT.
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

do jx=1,kx_num
   jx_s = kx_veci( jx )
   ii = 2*jx_s     ! imag index
   ir = ii-1     ! real index
do jy=1,ny
   fhat(jx_s,jy) = cmplx( f(ir,jy), f(ii,jy), rprec )
enddo
enddo

f(:,:) = 0._rprec

do i= 2, nx/2     !!jb - I believe the +1 not needed here... test
  jx_r = nx - i + 2
  fhat(jx_r,1) = conjg( fhat(i,1) )
enddo

do j = 2     , ny/2+1
  jy_r = ny - j + 2
do i = nx/2+2, nx
  jx_r = nx - i + 2
  fhat(i,j) = conjg(fhat(jx_r ,jy_r))
enddo
enddo

!! Nyquist mode ( kx=nx/2 ) is not taken into account
do i=1,nx
do jx=2,kx_num
jx_s = kx_veci( jx )
jx2 = nx - jx_s + 2
fhat2(i,:)=fhat2(i,:)+fhat(jx_s,:) * expp(jx,i)  !exp( c1*kx_vec(jx)*real(i-1,rprec) )
fhat2(i,:)=fhat2(i,:)+fhat(jx2,:) * expn(jx,i)    !exp(-c1*kx_vec(jx)*real(i-1,rprec) )
end do
end do

!! add in kx=0
do i=1,nx
  fhat2(i,:) = fhat2(i,:) + fhat(1,:)
enddo

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

return
end subroutine dft_direct_back_2d_n


!**********************************************************************
subroutine filt_da_direct(f,dfdx,dfdy, lbz)
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

  f(:,:,jz)=const*f(:,:,jz)  

  !$if ($FFTW3)
  call dft_direct_forw_2d( f(:,:,jz) )
  !$else
  !write(*,*) 'WARNING - no subroutine available'
  !$endif

  ! Zero padded region and Nyquist frequency
  f(ld-1:ld,:,jz) = 0._rprec
  f(:,ny/2+1,jz) = 0._rprec   

  !  Compute in-plane derivatives
  dfdx(:,:,jz) = f(:,:,jz) .MULI. kx
  dfdy(:,:,jz) = f(:,:,jz) .MULI. ky

  ! the oddballs for derivatives should already be dead, since they are for f
  ! inverse transform 
  !$if ($FFTW3)
  call dft_direct_back_2d( f(:,:,jz) )
  call dft_direct_back_2d( dfdx(:,:,jz) )
  call dft_direct_back_2d( dfdy(:,:,jz) )
  !$else
  !write(*,*) 'WARNING - no subroutine available'
  !$endif  

end do

return
end subroutine filt_da_direct

!**********************************************************************
subroutine filt_da_direct_n(f,dfdx,dfdy, lbz)
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

  f(:,:,jz)=const*f(:,:,jz)  

  !$if ($FFTW3)
  call dft_direct_forw_2d_n( f(:,:,jz) )
  !$else
  !write(*,*) 'WARNING - no subroutine available'
  !$endif

  ! Zero padded region and Nyquist frequency
  f(ld-1:ld,:,jz) = 0._rprec
  f(:,ny/2+1,jz) = 0._rprec   

  !  Compute in-plane derivatives
  dfdx(:,:,jz) = f(:,:,jz) .MULI. kx
  dfdy(:,:,jz) = f(:,:,jz) .MULI. ky

  ! the oddballs for derivatives should already be dead, since they are for f
  ! inverse transform 
  !$if ($FFTW3)
  call dft_direct_back_2d_n( f(:,:,jz) )
  call dft_direct_back_2d_n( dfdx(:,:,jz) )
  call dft_direct_back_2d_n( dfdy(:,:,jz) )
  !$else
  !write(*,*) 'WARNING - no subroutine available'
  !$endif  

end do

return
end subroutine filt_da_direct_n
$endif

end module derivatives
