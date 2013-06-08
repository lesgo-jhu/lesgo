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

module test_filtermodule
use types,only:rprec
use param,only:lh,ny

private lh,ny

    integer, parameter::filter_size=1   ! the implicit filter (1=grid size)
    real(rprec) :: alpha_test = 2.0_rprec * filter_size ! alpha is ratio of test filter to grid filter widths
    real(rprec) :: alpha_test_test = 4.0_rprec * filter_size

    real(rprec), dimension(:,:), allocatable :: G_test, G_test_test

contains

!**********************************************************************
subroutine test_filter_init()
!**********************************************************************
! Creates the kernels which will be used for filtering the field

use types,only:rprec
use param,only:lh,nx,ny,dx,dy,pi,ifilter,sgs_model
use fft

implicit none

real(rprec):: delta_test, kc2_test, delta_test_test, kc2_test_test

! Allocate the arrays
    allocate ( G_test(lh,ny) )

! Include the normalization
    G_test = 1._rprec/(nx*ny)  ! normalization for the forward FFT

! Filter characteristic width
    delta_test = alpha_test * sqrt(dx*dy)  ! "2d-delta", not full 3d one

! Calculate the kernel
    if(ifilter==1) then ! spectral cutoff filter
       if (sgs_model==6.OR.sgs_model==7) then
          print *, 'Use Gaussian or Top-hat filter for mixed models'
          stop
       endif
       kc2_test = (pi/(delta_test))**2
       where (real(k2, rprec) >= kc2_test) G_test = 0._rprec
    else if(ifilter==2) then ! Gaussian filter
       G_test=exp(-(delta_test)**2*k2/(4._rprec*6._rprec))*G_test   
    else if(ifilter==3) then ! Top-hat (Box) filter
       G_test= (sin(kx*delta_test/2._rprec)*sin(ky*delta_test/2._rprec)+1E-8)/&
            (kx*delta_test/2._rprec*ky*delta_test/2._rprec+1E-8)*G_test
    endif

! since our k2 has zero at Nyquist, we have to do this by hand
   G_test(lh,:) = 0._rprec
   G_test(:,ny/2+1) = 0._rprec


! Second test filter, if necessary
if ((sgs_model == 3) .or. (sgs_model == 5)) then  !--scale dependent dynamic

    ! Allocate the arrays
        allocate ( G_test_test(lh,ny) )

    ! Include the normalization
        G_test_test = 1._rprec/(nx*ny)    

    ! Filter characteristic width
        delta_test_test = alpha_test_test * sqrt(dx*dy)

    ! Calculate the kernel
        if(ifilter==1) then ! spectral cutoff filter
           if (sgs_model==6.OR.sgs_model==7) then
              print *, 'Use Gaussian or Top-hat filter for mixed models'
              stop
           endif
           kc2_test_test = (pi/(delta_test_test))**2
           where (real(k2, rprec) >= kc2_test_test) G_test_test = 0._rprec
        else if(ifilter==2) then ! Gaussian filter
           G_test_test=exp(-(delta_test_test)**2*k2/(4._rprec*6._rprec))*G_test_test  
        else if(ifilter==3) then ! Top-hat (Box) filter
           G_test_test= (sin(kx*delta_test_test/2._rprec)*sin(ky*delta_test_test/2._rprec)+1E-8)/&
                (kx*delta_test_test/2._rprec*ky*delta_test_test/2._rprec+1E-8)*G_test_test
        endif

    ! since our k2 has zero at Nyquist, we have to do this by hand
       G_test_test(lh,:) = 0._rprec
       G_test_test(:,ny/2+1) = 0._rprec

endif

return   
end subroutine test_filter_init
!**********************************************************************

!**********************************************************************
subroutine test_filter ( f )
!**********************************************************************
! note: this filters in-place, so input is ruined
use types,only:rprec
use fft
use emul_complex, only : OPERATOR(.MULR.)
implicit none

  real(rprec), dimension(:,:), intent(inout) :: f

!  Perform in-place FFT
  call rfftwnd_f77_one_real_to_complex(forw,f,fftwNull_p)

! Perform f = G_test*f, emulating f as complex
! Nyquist frequency and normalization is taken care of with G_test
  f = f .MULR. G_test

  call rfftwnd_f77_one_complex_to_real(back,f,fftwNull_p)

return
end subroutine test_filter

!**********************************************************************
subroutine test_test_filter ( f )
!**********************************************************************
! note: this filters in-place, so input is ruined
use types,only:rprec
use fft
use emul_complex, only : OPERATOR(.MULR.)
implicit none

  real(rprec), dimension(:,:), intent(inout) :: f

!  Perform in-place FFT
  call rfftwnd_f77_one_real_to_complex(forw,f,fftwNull_p)

! Perform f = G_test*f, emulating f as complex
! Nyquist frequency and normalization is taken care of with G_test_test
  f = f .MULR. G_test_test

  call rfftwnd_f77_one_complex_to_real(back,f,fftwNull_p)

return
end subroutine test_test_filter

end module test_filtermodule
