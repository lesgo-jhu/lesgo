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

!**********************************************************************
program main
!**********************************************************************
!
! Main file for lesgo solver 
! Contains main time-loop
! 

use types, only : rprec
use clocks 
use param
use sim_param
use grid_defs, only : grid_build
use io, only : energy, output_loop, output_final, jt_total
use io, only : write_tau_wall, energy_kx_spectral_complex, energy_kx_spectral_complex_fourier
use fft
use derivatives, only : filt_da, ddz_uv, ddz_w
use derivatives, only : ddx, ddy
use derivatives, only : dft_direct_forw_2d_n, dft_direct_back_2d_n  !!jb
use derivatives, only : dft_direct_forw_2d_n_yonly, dft_direct_back_2d_n_yonly  !!jb
use derivatives, only : dft_direct_forw_2d_n_yonlyC, dft_direct_back_2d_n_yonlyC  !!jb
use derivatives, only : dft_direct_forw_2d_n_yonlyC_big, dft_direct_back_2d_n_yonlyC_big  !!jb
use derivatives, only : wave2phys, phys2wave, filt_da_kxspace, convolve, convolve2  !!jb
use derivatives, only : wave2phys_pr, phys2wave_pr
use derivatives, only : wave2physF, phys2waveF
use derivatives, only : dft_direct_back_2d_n_big
use test_filtermodule
use cfl_util
use sgs_hist
use sgs_stag_util, only : sgs_stag, sgs_stag_fourier
use forcing
use functions, only: x_avg, get_tau_wall, interleave_r2c, interleave_c2r
!!use emul_complex, only : OPERATOR(.MULC.)    !! jb

$if ($MPI)
use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWN
$endif

$if ($LVLSET)
use level_set, only : level_set_global_CA, level_set_vel_err
use level_set_base, only : global_CA_calc
  
$if ($RNS_LS)
use rns_ls, only : rns_elem_force_ls
$endif
$endif

$if ($TURBINES)
use turbines, only : turbines_forcing, turbine_vel_init
$endif

$if ($DEBUG)
use debug_mod
$endif

use messages

implicit none

character (*), parameter :: prog_name = 'main'

$if ($DEBUG)
logical, parameter :: DEBUG = .false.
$endif

integer :: jt_step, nstart,jx,jy,jz  !!jb jx,jy
integer :: qq  !!jb
real(kind=rprec) rmsdivvel,ke, maxcfl
real (rprec):: tt
real (rprec) :: triggerFactor    !!jb
real (rprec) :: jtime1, jtime2   !!jb

!!complex(rprec), dimension(8) :: uc, vc, wc   !!jb
!!complex(rprec), dimension(8,8) :: uc
!!real(rprec), dimension(16) :: ur, vr, wr   !!jb
!!real(rprec), dimension(16,32,32) :: u2, v2
!!real(rprec), dimension(20,9,9) :: u2_big, v2_big
!!real(rprec), dimension(10,12,9) :: q
!!real(rprec), dimension(10,8,9) :: ut

type(clock_t) :: clock, clock_total

$if($MPI)
! Buffers used for MPI communication
real(rprec) :: rbuffer
$endif

! Start the clocks, both local and total
call clock_start( clock )

! Initialize time variable
tt = 0
jt = 0
jt_total = 0

! Initialize all data
call initialize()

if(coord == 0) then
   call clock_stop( clock )
   $if($MPI)
   write(*,'(1a,E15.7)') 'Initialization wall time: ', clock % time
   $else
   write(*,'(1a,E15.7)') 'Initialization cpu time: ', clock % time
   $endif
   
   !$if ($USE_RNL)
   !!kx_vec = kx_vec * 2._rprec * pi / L_x   !! aspect ratio change
   
   if ( coord == 0 .and. fourier) then
      write(*,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write(*,*) 'SIMULATING IN FOURIER SPACE !!!'
      write(*,*) '<<<<<<<<<<<<<<<<<<<<<<<<<<'
   endif

   if ( coord == 0 .and. fourier) then   
      write(*,*) '=================================='
      write(*,*) 'RNL modes >>>> '
      write(*,*) 'kx_num: ', kx_num
      write(*,*) 'kx_vec: ', kx_vec
      write(*,*) 'kx_veci: ', kx_veci
      write(*,*) 'kx_veci_n: ', kx_veci_n
      write(*,*) 'kx_veci_n_big: ', kx_veci_n_big
      write(*,*) 'L_x, L_y: ', L_x, L_y
      write(*,*) '=================================='
   endif
   !$endif
   
endif

call clock_start( clock_total )

! Initialize starting loop index 
! If new simulation jt_total=0 by definition, if restarting jt_total
! provided by total_time.dat
nstart = jt_total+1

! BEGIN TIME LOOP
time_loop: do jt_step = nstart, nsteps   
   ! Get the starting time for the iteration
   call clock_start( clock )
   call cpu_time(jtime1)   !!jb

!!$   if (coord == 0) then
!!$      !do jx=1,nx
!!$      !do jy=1,ny
!!$      !u(jx,jy,:) = 2.2 + 2.1*sin(L_x/(nx)*(jx-1)*1.0) + 2.3*sin(L_y/(ny)*(jy-1)*1.0)
!!$      !u(jx,jy,:) = u(jx,jy,:) + 2.7*sin(L_x/(nx)*(jx-1)*2.0) + 2.4*sin(L_y/(ny)*(jy-1)*2.0)
!!$      !u(jx,jy,:) = u(jx,jy,:) + 3.7*sin(L_x/(nx)*(jx-1)*5.0) + 3.4*sin(L_y/(ny)*(jy-1)*3.0)
!!$      !u(jx,jy,:) = u(jx,jy,:) + 4.7*sin(L_x/(nx)*(jx-1)*6.0) + 1.4*sin(L_y/(ny)*(jy-1)*2.0)
!!$      !enddo
!!$      !enddo
!!$      v = u
!!$      w = u
!!$      call wave2phys( v )
!!$      print*, '>>      u:', u(:,1,1)
!!$      print*, '>>      v:',  sum( v(1:nx,1,1) ) / nx
!!$      print*, '>>    v^2:', (sum( v(1:nx,1,1) ) / nx)**2
!!$      print*, '>>    v^3:', (sum( v(1:nx,1,1) ) / nx)**3
!!$      call dft_direct_back_2d_n_yonlyC( u(:,:,1) )
!!$      print*, '>>   u(1):', u(1,1,1)
!!$      print*, '>> u(1)^2:', u(1,1,1)**2
!!$      print*, '>> u(1)^3:', u(1,1,1)**3
!!$   endif

!!$   if ( fourier ) then    !!jb
!!$      if ( coord == 0 ) then
!!$         v = u
!!$         call wave2phys( u )
!!$         call wave2physF( v , vF )
!!$         print*, 'u: ',  u(1:nx,1,1)
!!$         print*, 'v: ', vF(1:nx,1,1)
!!$      endif
!!$   endif


!!$  if (coord == 0) then
!!$     if ( fourier ) then    !!jb
!!$        w(:,:,:) = u(:,:,:)
!!$        call wave2phys( u )
!!$        v = u**2
!!$
!!$        call dft_direct_back_2d_n_yonlyC( w(:,:,3) )
!!$        u(:,:,3) = convolve2( w(:,:,3), w(:,:,3) )
!!$        call dft_direct_forw_2d_n_yonlyC( u(:,:,3) )
!!$        call wave2phys( u )
!!$     endif
!!$     print*, 'test >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
!!$     do jx=1,nx
!!$        do jy=1,ny
!!$           write(*,*) v(jx,jy,3), u(jx,jy,3)
!!$        enddo
!!$     enddo
!!$  endif


!!$  if (coord == 0) then
!!$     if ( fourier ) then    !!jb
!!$        call wave2phys( u )
!!$        call wave2phys( v )
!!$        call wave2phys( w )
!!$     endif
!!$     print*, 'test >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
!!$     do jx=1,nx
!!$        do jy=1,ny
!!$           write(*,*) u(jx,jy,1), v(jx,jy,3), w(jx,jy,nz-1)
!!$        enddo
!!$     enddo
!!$  endif


!!$  if ( fourier ) then    !!jb
!!$     call phys2wave( u )
!!$     call phys2wave( v )
!!$     call phys2wave( w )
!!$  endif
!!$
!!$  if (coord == 0) then
!!$      print*, '4 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
!!$      do jx=1,3
!!$         do jy=1,3
!!$            write(*,*) v(jx,jy,0),v(jx,jy,1),v(jx,jy,nz-1),v(jx,jy,3)
!!$         enddo
!!$      enddo
!!$   endif
!!$
!!$  if ( fourier ) then    !!jb
!!$     call wave2phys( u )
!!$     call wave2phys( v )
!!$     call wave2phys( w )
!!$  endif
!!$
!!$  if (coord == 0) then
!!$      print*, '5 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
!!$      do jx=1,3
!!$         do jy=1,3
!!$            write(*,*) v(jx,jy,0),v(jx,jy,1),v(jx,jy,nz-1),v(jx,jy,3)
!!$         enddo
!!$      enddo
!!$   endif

   if( use_cfl_dt ) then

      if ( fourier ) then   !!jb
         call wave2physF( u, uF )
         call wave2physF( v, vF )
         call wave2physF( w, wF )
      endif

      dt_f = dt
      dt = get_cfl_dt()
      dt_dim = dt * z_i / u_star
    
      tadv1 = 1._rprec + 0.5_rprec * dt / dt_f
      tadv2 = 1._rprec - tadv1

   endif

!!$   if (coord == 0) then
!!$     if ( fourier ) then    !!jb
!!$        call wave2phys( u )
!!$        call wave2phys( v )
!!$        call wave2phys( w )
!!$     endif
!!$     print*, 'test >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
!!$     do jx=1,3
!!$        do jy=1,3
!!$           write(*,*) u(jx,jy,nz), v(jx,jy,nz-1), w(jx,jy,1)
!!$        enddo
!!$     enddo
!!$  endif


   ! Advance time
   jt_total = jt_step
   jt = jt + 1
   total_time = total_time + dt
   total_time_dim = total_time_dim + dt_dim
   tt=tt+dt
  
    $if ($DEBUG)
        if (DEBUG) write (*, *) $str($context_doc), ' reached line ', $line_num
    $endif

    ! Save previous time's right-hand-sides for Adams-Bashforth Integration
    ! NOTE: RHS does not contain the pressure gradient
    RHSx_f = RHSx
    RHSy_f = RHSy
    RHSz_f = RHSz

    $if ($DEBUG)
    if (DEBUG) then
        call DEBUG_write (u(:, :, 1:nz), 'main.p.u')
        call DEBUG_write (v(:, :, 1:nz), 'main.p.v')
        call DEBUG_write (w(:, :, 1:nz), 'main.p.w')
    end if
    $endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!$if (coord == 0) then
!!$print*, 'TEST DFTs >>>>>>>>>>>>>>>>>>>>>>>>'
!!$
!!$if (fourier) then
!!$   call wave2phys(u)
!!$endif
!!$
!!$do jx=1,nx
!!$do jy=1,ny
!!$do jz=1,nz
!!$!u(jx,jy,jz) = 2.2 + 2.1*sin(L_x/(nx)*(jx-1)*1.0)
!!$!u(jx,jy,jz) = 2.7 + 2.1*sin(L_x/(nx)*(jx-1)*1.0) + 1.7*sin(L_x/(nx)*(jx-1)*2.0) + 1.3*sin(L_x/(nx)*(jx-1)*3.0) + 1.6*sin(L_x/(nx)*(jx-1)*4.0)
!!$v(jx,jy,jz) = u(jx,jy,jz)
!!$w(jx,jy,jz) = u(jx,jy,jz)
!!$enddo
!!$enddo
!!$enddo
!!$
!!$call dfftw_execute_dft_r2c(forw, u(:,:,4), u(:,:,4))
!!$!call dfftw_execute_dft_r2c(forw, v(:,:,4), v(:,:,4))
!!$call dft_direct_forw_2d_n( v(:,:,4) )
!!$call dfftw_execute_dft_r2c(forw, w(:,:,4), w(:,:,4))
!!$call padd( u_big(:,:,4), u(:,:,4) )
!!$call padd( v_big(:,:,4), v(:,:,4) )
!!$call padd( w_big(:,:,4), w(:,:,4) )
!!$
!!$!call dfftw_execute_dft_c2r(back_big, u_big(:,:,4), u_big(:,:,4))
!!$!call dfftw_execute_dft_c2r(back_big, v_big(:,:,4), v_big(:,:,4))
!!$!call dfftw_execute_dft_c2r(back_big, w_big(:,:,4), w_big(:,:,4))
!!$
!!$print*, '----------------------------------------'
!!$
!!$!call dft_direct_back_2d_n_yonlyC_big( v_big(:,:,4) )
!!$!call dft_direct_forw_2d_n_yonlyC_big( v_big(:,:,4) )
!!$
!!$call dfftw_execute_dft_c2r( back_big, u_big(:,:,4), u_big(:,:,4) )
!!$call dft_direct_back_2d_n_big( v_big(:,:,4) )
!!$call dfftw_execute_dft_c2r( back_big, w_big(:,:,4), w_big(:,:,4) )
!!$
!!$!call dfftw_execute_dft_r2c( forw_big, u_big(:,:,4), u_big(:,:,4) )
!!$!call dfftw_execute_dft_r2c( forw_big, v_big(:,:,4), v_big(:,:,4) )
!!$!call dfftw_execute_dft_r2c( forw_big, w_big(:,:,4), w_big(:,:,4) )
!!$
!!$!call dfftw_execute_dft_c2r(back, u(:,:,4), u(:,:,4))
!!$!call dft_direct_back_2d_n( v(:,:,4) )
!!$
!!$print*, 'after >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
!!$do jx=1,ld_big
!!$do jy=1,ny2
!!$   write(*,*) jx, jy, u_big(jx,jy,4), v_big(jx,jy,4), w_big(jx,jy,4)
!!$   !!write(*,*) jx, jy, u(jx,jy,4), v(jx,jy,4), w(jx,jy,4)*real(nx*ny,rprec)
!!$enddo
!!$enddo
!!$print*, 'end test >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
!!$endif



!!$!! testing convolve >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!!$if (coord == 0) then 
!!$print*, 'TEST CONVOLVE >>>>>>>>>>>>>>'
!!$do jx=1,nx
!!$do jy=1,ny
!!$do jz=1,nz
!!$u(jx,jy,jz) = 2.7 + 2.1*sin(L_x/(nx)*(jx-1)*1.0) + 1.7*sin(L_x/(nx)*(jx-1)*2.0) + 1.3*sin(L_x/(nx)*(jx-1)*3.0) + 1.6*sin(L_x/(nx)*(jx-1)*4.0) !+ 1.6*sin(L_y/(ny)*(jy-1)*1.0)
!!$v(jx,jy,jz) = 1.2 + 1.2*sin(L_x/(nx)*(jx-1)*1.0) + 6.3*sin(L_x/(nx)*(jx-1)*2.0) + 3.1*sin(L_x/(nx)*(jx-1)*3.0)+ 1.2*sin(L_x/(nx)*(jx-1)*4.0) !+ 3.7*sin(L_y/(ny)*(jy-1)*1.0)
!!$!u(jx,jy,jz) = 2.2 + 2.1*sin(L_x/(nx)*(jx-1)*1.0)
!!$!v(jx,jy,jz) = 1.4 + 1.3*sin(L_x/(nx)*(jx-1)*1.0)
!!$!v(jx,jy,jz) = u(jx,jy,jz)
!!$w(jx,jy,jz) = u(jx,jy,jz)
!!$enddo
!!$enddo
!!$enddo
!!$
!!$u_big = 0._rprec
!!$v_big = 0._rprec
!!$w_big = 0._rprec
!!$u(:,:,4) = u(:,:,4) / (nx*ny)
!!$v(:,:,4) = v(:,:,4) / (nx*ny)
!!$call dfftw_execute_dft_r2c(forw, u(:,:,4), u(:,:,4))
!!$call dfftw_execute_dft_r2c(forw, v(:,:,4), v(:,:,4))
!!$call padd( u_big(:,:,4), u(:,:,4) )
!!$call padd( v_big(:,:,4), v(:,:,4) )
!!$call dfftw_execute_dft_r2c(back_big, u_big(:,:,4), u_big(:,:,4))
!!$call dfftw_execute_dft_r2c(back_big, v_big(:,:,4), v_big(:,:,4))
!!$
!!$w_big(:,:,3) = u_big(:,:,4) * v_big(:,:,4)
!!$call dfftw_execute_dft_r2c(forw_big, w_big(:,:,3), w_big(:,:,3))
!!$call unpadd( w(:,:,3), w_big(:,:,3) )
!!$
!!$
!!$call dfftw_execute_dft_r2c(forw_big, u_big(:,:,4), u_big(:,:,4))
!!$call dfftw_execute_dft_r2c(forw_big, v_big(:,:,4), v_big(:,:,4))
!!$!call dft_direct_forw_2d_n( v(:,:,jz) )
!!$!call dft_direct_forw_2d_n_yonlyC( w(:,:,jz) )
!!$
!!$print*, 'u, v in kx, ky space: >>>>>>>>>>>>>>>>>>'
!!$do jx=1,ld_big
!!$   write(*,*) jx, u_big(jx,1,4), v_big(jx,1,4)
!!$enddo
!!$
!!$print*, 'CONVOLVED >>>>>>>>>>>>>>'
!!$w_big(:,:,4) = convolve( u_big(:,:,4), v_big(:,:,4) )
!!$print*, 'AFTER CONVOLVED >>>>>>>>>>>>>>'
!!$print*, '         convolved,         physical: '
!!$do jx=1,ld_big
!!$   write(*,*) jx, w_big(jx,1,4)/(nx*ny*3/2*3/2), w(jx,1,3)
!!$enddo
!!$endif
!!$!! testing convolve <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!!$!! testing yonly >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!!$if (coord == 0) then 
!!$print*, 'BEFORE >>>>>>>>>>>>>>'
!!$
!!$do jx=1,nx
!!$do jy=1,ny
!!$do jz=1,nz
!!$u(jx,jy,jz) = 2.7 + 4.5*(sin(L_x/(nx)*(jx-1)*1.0) + sin(L_y/(ny)*(jy-1)*2.0))
!!$v(jx,jy,jz) = u(jx,jy,jz)
!!$w(jx,jy,jz) = u(jx,jy,jz)
!!$enddo
!!$enddo
!!$enddo
!!$
!!$print*, 'u in (x, y) space: >>>>>>>>>>>>>>>>>>'
!!$do jx=1,ld
!!$   write(*,*) jx, u(jx,:,4)
!!$enddo
!!$
!!$call dfftw_execute_dft_r2c(forw, u(:,:,4), u(:,:,4))
!!$
!!$print*, 'u in (kx, ky) space: >>>>>>>>>>>>>>>>>>'
!!$do jx=1,ld
!!$   write(*,*) jx, u(jx,:,4)/(nx*ny)
!!$enddo
!!$
!!$call dft_direct_back_2d_n_yonlyC( u(:,:,4) )
!!$
!!$print*, 'u in (kx, y) space: >>>>>>>>>>>>>>>>>>'
!!$do jx=1,ld
!!$   write(*,*) jx, u(jx,:,4)/(nx*ny)
!!$enddo
!!$
!!$
!!$endif
!!$!! testing yonly <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


!!$!! testing yonly >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!!$if (coord == 0) then 
!!$print*, 'BEFORE >>>>>>>>>>>>>>'
!!$
!!$do jx=1,nx
!!$do jy=1,ny
!!$do jz=1,nz
!!$u(jx,jy,jz) = 2.7 + 4.5*(sin(L_x/(nx)*(jx-1)*1.0) + sin(L_y/(ny)*(jy-1)*1.0))
!!$!v(jx,jy,jz) = u(jx,jy,jz)
!!$!w(jx,jy,jz) = u(jx,jy,jz)
!!$enddo
!!$enddo
!!$enddo
!!$
!!$print*, 'u in (x, y) space: >>>>>>>>>>>>>>>>>>'
!!$do jx=1,ld
!!$   write(*,*) jx, u(jx,:,4)
!!$enddo
!!$call dfftw_execute_dft_r2c(forw, u(:,:,4), u(:,:,4))
!!$u = u / (nx*ny)
!!$call padd(u_big(:,:,4), u(:,:,4) )
!!$print*, 'u_big in (kx, ky) space, padded: >>>>>>>>>>>>>>>>>>'
!!$do jx=1,ld_big
!!$   write(*,*) jx, u_big(jx,:,4)
!!$enddo
!!$call dfftw_execute_dft_r2c(back_big, u_big(:,:,4), u_big(:,:,4))
!!$print*, 'u_big in (x, y) space: >>>>>>>>>>>>>>>>>>'
!!$do jx=1,ld_big
!!$   write(*,*) jx, u_big(jx,:,4)
!!$enddo
!!$call dfftw_execute_dft_r2c(forw_big, u_big(:,:,4), u_big(:,:,4))
!!$u_big = u_big / (nx2*ny2)
!!$print*, 'u_big in (kx, ky) space: >>>>>>>>>>>>>>>>>>'
!!$do jx=1,ld_big
!!$   write(*,*) jx, u_big(jx,:,4)
!!$enddo
!!$call dft_direct_back_2d_n_yonlyC_big( u_big(:,:,4) )
!!$print*, 'u_big in (kx, y) space: >>>>>>>>>>>>>>>>>>'
!!$do jx=1,ld_big
!!$   write(*,*) jx, u_big(jx,:,4)
!!$enddo
!!$call dft_direct_forw_2d_n_yonlyC_big( u_big(:,:,4) )
!!$print*, 'u_big in (kx, ky) space: >>>>>>>>>>>>>>>>>>'
!!$do jx=1,ld_big
!!$   write(*,*) jx, u_big(jx,:,4)
!!$enddo
!!$call unpadd( u(:,:,4), u_big(:,:,4) )
!!$call dfftw_execute_dft_r2c(back, u(:,:,4), u(:,:,4))
!!$print*, 'u in (x, y) space: >>>>>>>>>>>>>>>>>>'
!!$do jx=1,ld
!!$   write(*,*) jx, u(jx,:,4)
!!$enddo
!!$
!!$endif
!!$!! testing yonly <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


!!$do jz=1,nz
!!$   call dfftw_execute_dft_r2c(back, u(:,:,jz), u(:,:,jz))
!!$   call dft_direct_back_2d_n( v(:,:,jz) )
!!$   call dft_direct_back_2d_n_yonlyC( w(:,:,jz) )
!!$enddo

!!$
!!$if (coord == 0) then
!!$ do jx=1,ld
!!$ do jy=1,ny
!!$    write(*,*) jx, jy, u(jx,jy,4)/(nx*ny), v(jx,jy,4)/(nx*ny), w(jx,jy,4)/ny
!!$ enddo
!!$ enddo
!!$endif
!!$
!!$if (coord == 0) then 
!!$  print*, 'DERIVS >>>>>>>>>>>>>>'
!!$endif
!!$
!!$do jx=1,nx
!!$do jy=1,ny
!!$do jz=1,nz
!!$u(jx,jy,jz) = 1.0 + sin(L_y/(ny)*(jy-1)*1.0) + sin(L_x/(nx)*(jx-1)*2.0) 
!!$v(jx,jy,jz) = u(jx,jy,jz)
!!$w(jx,jy,jz) = u(jx,jy,jz)
!!$enddo
!!$enddo
!!$enddo
!!$
!!$if (coord == 0) then
!!$ do jx=1,ld
!!$ do jy=1,ny
!!$    write(*,*) jx, jy, u(jx,jy,1), v(jx,jy,1)
!!$ enddo
!!$ enddo
!!$endif
!!$
!!$
!!$call filt_da (u, dudx, dudy, lbz)
!!$call filt_da_direct_n (v, dvdx, dvdy, lbz)    
!!$
!!$if (coord == 0) then
!!$ do jx=1,ld
!!$ do jy=1,ny
!!$    write(*,*) jx, jy, dudx(jx,jy,1), dvdx(jx,jy,1)
!!$ enddo
!!$ enddo
!!$endif
!!$
!!$if (coord == 0) then
!!$ do jx=1,ld
!!$ do jy=1,ny
!!$    write(*,*) jx, jy, dudy(jx,jy,1), dvdy(jx,jy,1) !*ny
!!$ enddo
!!$ enddo
!!$endif
!!$
!!$if (coord == 0) then 
!!$  print*, 'AFTER >>>>>>>>>>>>>>'
!!$endif

!!$if (coord == 0) then
!!$ do jx=1,ld
!!$ do jy=1,ny
!!$    write(*,*) jx, jy, u(jx,jy,4)/(nx*ny), v(jx,jy,4)/(nx*ny), w(jx,jy,4)/(ny)
!!$ enddo
!!$ enddo
!!$endif



!!$uc(1) = cmplx(1, 2, rprec)
!!$uc(2) = cmplx(3, 4, rprec)
!!$uc(3) = cmplx(5, 6, rprec)
!!$uc(4) = cmplx(7, 8, rprec)
!!$uc(5) = cmplx(9, 10, rprec)
!!$uc(6) = cmplx(11, 12, rprec)
!!$uc(7) = cmplx(13, 14, rprec)
!!$uc(8) = cmplx(15, 16, rprec)
!!$
!!$
!!$do jy = 1, 16
!!$ ur(jy) = jy
!!$enddo
!!$
!!$if (coord == 0) then
!!$ do jy=1,8
!!$  write(*,*) jy, uc(jy)
!!$ enddo
!!$endif
!!$
!!$if (coord == 0) then
!!$ do jy=1,16
!!$  write(*,*) jy, ur(jy)
!!$ enddo
!!$endif
!!$
!!$call dfftw_execute_dft(ccomp_forw, uc(:), uc(:))
!!$call dfftw_execute_dft(rcomp_forw, ur(:), ur(:))

!!$if (coord == 0) then 
!!$  print*, 'AFTER >>>>>>>>>>>>>>'
!!$endif

!!$if (coord == 0) then
!!$ do jy=1,8
!!$  write(*,*) jy, uc(jy)
!!$ enddo
!!$endif
!!$
!!$if (coord == 0) then
!!$ do jy=1,16
!!$  write(*,*) jy, ur(jy)
!!$ enddo
!!$endif
!!$
!!$call dfftw_execute_dft(ccomp_back, uc(:), uc(:))
!!$call dfftw_execute_dft(rcomp_back, ur(:), ur(:))

!!$if (coord == 0) then 
!!$  print*, 'AFTER BACK >>>>>>>>>>>>>>'
!!$endif

!!$if (coord == 0) then
!!$ do jy=1,8
!!$  write(*,*) jy, uc(jy) / 8
!!$ enddo
!!$endif
!!$
!!$if (coord == 0) then
!!$ do jy=1,16
!!$  write(*,*) jy, ur(jy) / 8
!!$ enddo
!!$endif


!!$if (coord == 0) then 
!!$  print*, 'BEFORE >>>>>>>>>>>>>>'
!!$endif
!!$
!!$do jx=1,nx
!!$do jy=1,ny
!!$do jz=1,nz
!!$u(jx,jy,jz) = 2.0+sin(L_x/(nx)*(jx-1)*1.0)+sin(L_y/(ny)*(jy-1)*2.0)
!!$v(jx,jy,jz) = u(jx,jy,jz)
!!$w(jx,jy,jz) = u(jx,jy,jz)
!!$enddo
!!$enddo
!!$enddo
!!$!!sin( L_y/(ny-2)*(jy-1) * 3.0_rprec )
!!$
!!$if (coord == 0) then
!!$ do jx=1,ld
!!$ do jy=1,ny
!!$    write(*,*) jx, jy, u(jx,jy,1), v(jx,jy,1), w(jx,jy,1)
!!$ enddo
!!$ enddo
!!$endif
!!$
!!$do jz=1,nz
!!$  call dfftw_execute_dft_r2c(forw, u(:,:,jz), u(:,:,jz))
!!$  call dft_direct_forw_2d_n( v(:,:,jz) )  
!!$  call dft_direct_forw_2d_n_yonlyC( w(:,:,jz) )
!!$enddo
!!$
!!$if (coord == 0) then 
!!$  print*, 'MIDDLE >>>>>>>>>>>>>>>>>>'
!!$endif
!!$
!!$
!!$if (coord == 0) then
!!$ do jx=1,ld
!!$ do jy=1,ny
!!$    write(*,*) jx, jy, u(jx,jy,1), v(jx,jy,1), w(jx,jy,1)
!!$ enddo
!!$ enddo
!!$endif
!!$
!!$
!do jy=1,ny
!do jz=1,nz
!   call dfftw_execute_dft_r2c(forw_span_spectra2, v(:,jy,jz), v(:,jy,jz))
!enddo
!enddo
!!$
!!$do jz=1,nz
!!$   call dfftw_execute_dft_r2c(back, u(:,:,jz), u(:,:,jz))
!!$   call dft_direct_back_2d_n( v(:,:,jz) )
!!$   call dft_direct_back_2d_n_yonlyC( w(:,:,jz) )
!!$enddo
!!$
!!$if (coord == 0) then 
!!$  print*, 'AFTER >>>>>>>>>>>>>>'
!!$endif
!!$
!!$if (coord == 0) then
!!$ do jx=1,ld
!!$ do jy=1,ny
!!$    write(*,*) jx, jy, u(jx,jy,1)/(nx*ny), v(jx,jy,1)/(nx*ny), w(jx,jy,1)
!!$ enddo
!!$ enddo
!!$endif

!!$do jy=1,ny  ! ld = nx+2
!!$  if (coord .eq. 0) then
!!$    write(*,*) jy, u(1,jy,1)
!!$  endif
!!$enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    triggerFactor = 5.0_rprec    !!jb
    if (trigger) then
       if (jt_total == trig_on ) nu_molec = nu_molec / triggerFactor
       if (jt_total == trig_off) nu_molec = nu_molec * triggerFactor
    endif

    ! Calculate velocity derivatives
    ! Calculate dudx, dudy, dvdx, dvdy, dwdx, dwdy (in Fourier space)
!!$  do jx=1,nx
!!$  do jy=1,ny
!!$    u(jx,jy,:) = sin( L_x/nx*(jx-1) * 3.0_rprec ) + sin( L_x/nx*(jx-1) * 1.0_rprec ) + .56
!!$  enddo
!!$  enddo

    if (fourier) then
       call filt_da_kxspace (u, dudx, dudy, lbz)    
       call filt_da_kxspace (v, dvdx, dvdy, lbz)    
       call filt_da_kxspace (w, dwdx, dwdy, lbz)    
    else
       call filt_da (u, dudx, dudy, lbz)
       call filt_da (v, dvdx, dvdy, lbz)
       call filt_da (w, dwdx, dwdy, lbz)
    endif

    !call mpi_barrier(comm,ierr)
!!$    if (coord == 0) then
!!$       if (fourier) then
!!$          call wave2phys(dudx)
!!$          call wave2phys(dudy)
!!$          call wave2phys(dvdx)
!!$          call wave2phys(dvdy)
!!$          call wave2phys(dwdx)
!!$          call wave2phys(dwdy)
!!$          call wave2phys(u)
!!$          call wave2phys(v)
!!$          call wave2phys(w)
!!$       endif
!!$       print*, 'test 1 >>>>>>>>>>>>>>>>>>>>>>>'
!!$       do jx=1,3
!!$       do jy=1,3
!!$          write(*,*) dudx(jx,jy,3),dudy(jx,jy,3),dvdx(jx,jy,3)
!!$       enddo
!!$       enddo
!!$       print*, 'test 2 >>>>>>>>>>>>>>>>>>>>>>>'
!!$       do jx=1,3
!!$       do jy=1,3
!!$          write(*,*) dvdy(jx,jy,nz-1),dwdx(jx,jy,nz),dwdy(jx,jy,3)
!!$       enddo
!!$       enddo
!!$       print*, 'test 3 >>>>>>>>>>>>>>>>>>>>>>>'
!!$       do jx=1,3
!!$       do jy=1,3
!!$          write(*,*) u(jx,jy,3),v(jx,jy,1),w(jx,jy,2)
!!$       enddo
!!$       enddo
!!$    endif

    ! Calculate dudz, dvdz using finite differences (for 1:nz on uv-nodes)
    !  except bottom coord, only 2:nz
    call ddz_uv(u, dudz, lbz)
    call ddz_uv(v, dvdz, lbz)
       
    ! Calculate dwdz using finite differences (for 0:nz-1 on w-nodes)
    !  except bottom coord, only 1:nz-1
    call ddz_w(w, dwdz, lbz)

!!$    if (coord == 0) then
!!$       if (fourier) then
!!$          call wave2phys(dudz)
!!$          call wave2phys(dvdz)
!!$          call wave2phys(dwdz)
!!$       endif
!!$       print*, 'here: >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
!!$       do jx=1,3
!!$       do jy=1,3
!!$          write(*,*) dudz(jx,jy,2), dvdz(jx,jy,2), dwdz(jx,jy,2)
!!$       enddo
!!$       enddo
!!$    endif

    $if ($DEBUG)
    if (DEBUG) then
        call DEBUG_write (u(:, :, 1:nz), 'main.q.u')
        call DEBUG_write (v(:, :, 1:nz), 'main.q.v')
        call DEBUG_write (w(:, :, 1:nz), 'main.q.w')
        call DEBUG_write (dudx(:, :, 1:nz), 'main.q.dudx')
        call DEBUG_write (dudy(:, :, 1:nz), 'main.q.dudy')
        call DEBUG_write (dudz(:, :, 1:nz), 'main.q.dudz')
        call DEBUG_write (dvdx(:, :, 1:nz), 'main.q.dvdx')
        call DEBUG_write (dvdy(:, :, 1:nz), 'main.q.dvdy')
        call DEBUG_write (dvdz(:, :, 1:nz), 'main.q.dvdz')
        call DEBUG_write (dwdx(:, :, 1:nz), 'main.q.dwdx')
        call DEBUG_write (dwdy(:, :, 1:nz), 'main.q.dwdy')
        call DEBUG_write (dwdz(:, :, 1:nz), 'main.q.dwdz')
    end if
    $endif

!!$    if (coord == 0) then
!!$       if (fourier) then
!!$          call wave2phys( tyz )
!!$          call wave2phys( dudz )
!!$          call wave2phys( u )
!!$       endif
!!$       print*, 'here wall 1: >>>>>>>>>>>>>>>>>>'
!!$       do jx=1,3
!!$       do jy=1,3
!!$          write(*,*) dudz(jx,jy,1:3)
!!$       enddo
!!$       enddo
!!$       if (fourier) then
!!$          call phys2wave( tyz )
!!$          call phys2wave( dudz )
!!$          call phys2wave( u )
!!$       endif
!!$    endif

    
!!$    if (coord == 0) then
!!$       do jx=1,nx
!!$       do jy=1,ny
!!$          write(*,*) u(jx,jy,1), v(jx,jy,1), dudz(jx,jy,1)
!!$       enddo
!!$       enddo
!!$       print*, 'yoy >>>>>>>>>>>>>>>>>>>>>>>>'
!!$       do jx=1,nx
!!$       do jy=1,ny
!!$          write(*,*) dvdz(jx,jy,1), txz(jx,jy,1), tyz(jx,jy,1)
!!$       enddo
!!$       enddo
!!$    endif

    ! Calculate wall stress and derivatives at the wall (txz, tyz, dudz, dvdz at jz=1)
    !   using the velocity log-law
    !   MPI: bottom process only
    if (dns_bc) then
        if (coord == 0) then
            call wallstress_dns ()
        end if
    else    ! "impose" wall stress 
        if (coord == 0) then
           !if (fourier) then
           !   call wallstress_fourier ()
           !else
              call wallstress ()
           !endif
        end if
    end if

!!$    if (coord == 0) then
!!$       if (.not. fourier ) then
!!$          call phys2wave( dudz )
!!$          call phys2wave( dvdz )
!!$          call phys2wave( txz )
!!$          call phys2wave( tyz )
!!$          call phys2wave( u )
!!$          call phys2wave( v )
!!$       endif
!!$       print*, 'here wall 1: >>>>>>>>>>>>>>>>>>'
!!$       do jx=1,ld
!!$       do jy=1,ny
!!$          write(*,*) txz(jx,jy,1), tyz(jx,jy,1) 
!!$       enddo
!!$       enddo
!!$       print*, 'here wall 2: >>>>>>>>>>>>>>>>>>'
!!$       do jx=1,ld
!!$       do jy=1,ny
!!$          write(*,*) dudz(jx,jy,1), dvdz(jx,jy,1)
!!$       enddo
!!$       enddo
!!$    endif

    ! Calculate turbulent (subgrid) stress for entire domain
    !   using the model specified in param.f90 (Smag, LASD, etc)
    !   MPI: txx, txy, tyy, tzz at 1:nz-1; txz, tyz at 1:nz (stress-free lid)
    if (dns_bc .and. molec) then
!!$       $if ($MPI)    !--jb, copied from sgs_stag_util.f90
!!$       ! dudz calculated for 0:nz-1 (on w-nodes) except bottom process
!!$       ! (only 1:nz-1) exchange information between processors to set
!!$       ! values at nz from jz=1 above to jz=nz below
!!$       call mpi_sync_real_array( dwdz(:,:,1:), 1, MPI_SYNC_DOWN )
!!$       $endif
!!$
!!$       call dns_stress(txx,txy,txz,tyy,tyz,tzz)
!!$
!!$       $if ($MPI)    !--jb
!!$       call mpi_sync_real_array( txz, 0, MPI_SYNC_DOWN )
!!$       call mpi_sync_real_array( tyz, 0, MPI_SYNC_DOWN )
!!$       $endif
       call sgs_stag() 
    else     
!!$       if (fourier) then
!!$         call wave2phys( dudx )
!!$         call wave2phys( dudy )
!!$         call wave2phys( dudz )
!!$         call wave2phys( dvdx )
!!$         call wave2phys( dvdy )
!!$         call wave2phys( dvdz )
!!$         call wave2phys( dwdx )
!!$         call wave2phys( dwdy )
!!$         call wave2phys( dwdz )
!!$         
!!$         call wave2phys( txx )
!!$         call wave2phys( txy )
!!$         call wave2phys( txz )
!!$         call wave2phys( tyy )
!!$         call wave2phys( tyz )
!!$         call wave2phys( tzz )
!!$       endif
       
       !if (fourier) then
       !   call sgs_stag_fourier()
       !else
          call sgs_stag()    !! not updated yet for fourier for sgs=true
       !endif

!!$       if (fourier) then
!!$         call phys2wave( dudx )
!!$         call phys2wave( dudy )
!!$         call phys2wave( dudz )
!!$         call phys2wave( dvdx )
!!$         call phys2wave( dvdy )
!!$         call phys2wave( dvdz )
!!$         call phys2wave( dwdx )
!!$         call phys2wave( dwdy )
!!$         call phys2wave( dwdz )
!!$
!!$         call phys2wave( txx )
!!$         call phys2wave( txy )
!!$         call phys2wave( txz )
!!$         call phys2wave( tyy )
!!$         call phys2wave( tyz )
!!$         call phys2wave( tzz )
!!$       endif

     end if                 !! also should be updated to not calc Nu_t if sgs=false

!!$ if (coord == 0) then
!!$    if ( .not. fourier ) then
!!$       call phys2wave( txx )
!!$       call phys2wave( txy )
!!$       call phys2wave( txz )
!!$       call phys2wave( tyy )
!!$       call phys2wave( tyz )
!!$       call phys2wave( tzz )
!!$       !call phys2wave( u )
!!$       !call phys2wave( v )
!!$       !call phys2wave( w )
!!$       !call phys2wave( txz )
!!$    endif
!!$    print*, 'here wall 1: !>>>>>>>>>>>>>>>>>>'
!!$    do jx=1,nx
!!$    do jy=1,ny
!!$       write(*,*) txx(jx,jy,1), txy(jx,jy,1), txz(jx,jy,1)
!!$       !write(*,*) u(jx,jy,1), v(jx,jy,1), w(jx,jy,1)
!!$    enddo
!!$    enddo
!!$    print*, 'here wall 2: !>>>>>>>>>>>>>>>>>>'
!!$    do jx=1,nx
!!$    do jy=1,ny
!!$       write(*,*) tyy(jx,jy,1), tyz(jx,jy,1), tzz(jx,jy,1)
!!$       !write(*,*) u(jx,jy,3), v(jx,jy,3), w(jx,jy,3)
!!$    enddo
!!$    enddo
!!$endif

    ! Exchange ghost node information (since coords overlap) for tau_zz
    !   send info up (from nz-1 below to 0 above)
    $if ($MPI)
        call mpi_sendrecv (tzz(:, :, nz-1), ld*ny, MPI_RPREC, up, 6,   &
                           tzz(:, :, 0), ld*ny, MPI_RPREC, down, 6,  &
                           comm, status, ierr)
    $endif

    $if ($DEBUG)
    if (DEBUG) then
        call DEBUG_write (divtx(:, :, 1:nz), 'main.r.divtx')
        call DEBUG_write (divty(:, :, 1:nz), 'main.r.divty')
        call DEBUG_write (divtz(:, :, 1:nz), 'main.r.divtz')
        call DEBUG_write (txx(:, :, 1:nz), 'main.r.txx')
        call DEBUG_write (txy(:, :, 1:nz), 'main.r.txy')
        call DEBUG_write (txz(:, :, 1:nz), 'main.r.txz')
        call DEBUG_write (tyy(:, :, 1:nz), 'main.r.tyy')
        call DEBUG_write (tyz(:, :, 1:nz), 'main.r.tyz')
        call DEBUG_write (tzz(:, :, 1:nz), 'main.r.tzz')
    end if
    $endif    

    !!if (fourier) call wave2phys(txz)
    !!if (fourier) call wave2phys(tyz)

    ! Compute divergence of SGS shear stresses     
    !   the divt's and the diagonal elements of t are not equivalenced in this version
    !   provides divtz 1:nz-1, except 1:nz at top process
    call divstress_uv (divtx, divty, txx, txy, txz, tyy, tyz) ! saves one FFT with previous version

!!$    if (coord == 0) then
!!$       if (fourier) then
!!$          call wave2phys(divtx)
!!$          call wave2phys(divty)
!!$       endif
!!$       do jx=1,nx
!!$       do jy=1,ny
!!$          write(*,*) jx,jy,divtx(jx,jy,3),divty(jx,jy,3)
!!$       enddo
!!$       enddo
!!$    endif

    call divstress_w(divtz, txz, tyz, tzz)

!!$    if (coord == nproc-1) then
!!$       if (.not. fourier) then
!!$          !call wave2phys(divtx)
!!$          !call wave2phys(divty)
!!$          !call wave2phys(divtz)
!!$          call phys2wave(divtx)
!!$          call phys2wave(divty)
!!$          call phys2wave(divtz)
!!$          call phys2wave(txz)
!!$       endif
!!$       print*, 'compare divs: >>>>>>>>>>>>>>>>'
!!$       do jx=1,ld
!!$       do jy=1,ny
!!$          write(*,*) divtx(jx,jy,nz-1),divty(jx,jy,nz-1),divtz(jx,jy,nz-1)
!!$          !write(*,*) jx, jy, txz(jx,jy,1:2)
!!$       enddo
!!$       enddo
!!$    endif


    $if ($DEBUG)
    if (DEBUG) then
        call DEBUG_write (divtx(:, :, 1:nz), 'main.s.divtx')
        call DEBUG_write (divty(:, :, 1:nz), 'main.s.divty')
        call DEBUG_write (divtz(:, :, 1:nz), 'main.s.divtz')
        call DEBUG_write (RHSx(:, :, 1:nz), 'main.preconvec.RHSx')
    end if
    $endif

!!$    if (coord == nproc-1) then
!!$       if (fourier) then
!!$          call wave2phys(u)
!!$          call wave2phys(v)
!!$          call wave2phys(w)
!!$          call wave2phys(dudy)
!!$          call wave2phys(dudz)
!!$          call wave2phys(dvdx)
!!$          call wave2phys(dvdz)
!!$          call wave2phys(dwdx)
!!$          call wave2phys(dwdy)
!!$          call wave2phys(RHSx)
!!$          call wave2phys(RHSy)
!!$          call wave2phys(RHSz)
!!$       endif
!!$       print*, 'yot: >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
!!$       qq = 1
!!$       do jx=1,nx
!!$       do jy=1,ny
!!$          write(*,*) u(jx,jy,qq), v(jx,jy,qq), w(jx,jy,qq)
!!$          write(*,*) dudy(jx,jy,qq), dudz(jx,jy,qq), dvdx(jx,jy,qq)
!!$          write(*,*) dvdz(jx,jy,qq), dwdx(jx,jy,qq), dwdy(jx,jy,qq)
!!$          write(*,*) RHSx(jx,jy,qq), RHSy(jx,jy,qq), RHSz(jx,jy,qq)
!!$          print*, '---------------------------------------------'
!!$       enddo
!!$       enddo
!!$    endif


    ! Calculates u x (omega) term in physical space. Uses 3/2 rule for
    ! dealiasing. Stores this term in RHS (right hand side) variable
    $if ($USE_RNL)  
    call convec(u,v,w,dudy,dudz,dvdx,dvdz,dwdx,dwdy,RHSx,RHSy,RHSz)

    if (.not. fourier) then
       u_rnl = u - x_avg(u)
       v_rnl = v - x_avg(v)
       w_rnl = w - x_avg(w)
       dudy_rnl = dudy - x_avg(dudy)
       dudz_rnl = dudz - x_avg(dudz)
       dvdx_rnl = dvdx - x_avg(dvdx)
       dvdz_rnl = dvdz - x_avg(dvdz)
       dwdx_rnl = dwdx - x_avg(dwdx)
       dwdy_rnl = dwdy - x_avg(dwdy)

       call convec(u_rnl,v_rnl,w_rnl,dudy_rnl,dudz_rnl,dvdx_rnl,dvdz_rnl,dwdx_rnl,dwdy_rnl,RHSx_rnl,RHSy_rnl,RHSz_rnl)

       RHSx_rnl = RHSx_rnl - x_avg(RHSx_rnl)
       RHSy_rnl = RHSy_rnl - x_avg(RHSy_rnl)
       RHSz_rnl = RHSz_rnl - x_avg(RHSz_rnl)

       RHSx = RHSx - RHSx_rnl
       RHSy = RHSy - RHSy_rnl
       RHSz = RHSz - RHSz_rnl
    endif

!!$    if (coord == 0) then
!!$       if (.not. fourier) then
!!$          call phys2wave( RHSx )
!!$          call phys2wave( RHSy )
!!$          call phys2wave( RHSz )
!!$          !call phys2wave( u )
!!$          !call phys2wave( v )
!!$          !call phys2wave( w )
!!$       endif
!!$       print*, 'here wall: >>>>>>>>>>>>>>>>>>'
!!$       do jx=1,ld
!!$          do jy=1,ny
!!$             write(*,*) RHSx(jx,jy,3), RHSy(jx,jy,3), RHSz(jx,jy,3)
!!$             !write(*,*) u(jx,jy,3), v(jx,jy,3), w(jx,jy,3)
!!$          enddo
!!$       enddo
!!$    endif


!!$    if (coord == 0) then
!!$       if (fourier) then
!!$          print*, 'fourier rhs: >>>>>>>>>>>'
!!$          call wave2phys(RHSx)
!!$          call wave2phys(RHSy)
!!$          call wave2phys(RHSz)
!!$
!!$          do jx=1,nx
!!$          do jy=1,ny
!!$             write(*,*) jx, jy, RHSx(jx,jy,nz-1), RHSy(jx,jy,nz-1), RHSz(jx,jy,1)
!!$          enddo
!!$          enddo
!!$       else
!!$        !call dfftw_execute_dft_r2c(forw, RHSx(:,:,5), RHSx(:,:,5))
!!$        !call dfftw_execute_dft_r2c(forw, RHSy(:,:,5), RHSy(:,:,5))
!!$        !call dfftw_execute_dft_r2c(forw, RHSz(:,:,5), RHSz(:,:,5))
!!$        print*, 'phys rhs: >>>>>>>>>>>>>>>>'
!!$        do jx=1,nx
!!$        do jy=1,ny
!!$           write(*,*) jx, jy, RHSx(jx,jy,nz-1), RHSy(jx,jy,nz-1), RHSz(jx,jy,1)
!!$        enddo
!!$        enddo
!!$
!!$     endif
!!$
!!$    endif

    $else
    call convec(u,v,w,dudy,dudz,dvdx,dvdz,dwdx,dwdy,RHSx,RHSy,RHSz)
    $endif

!!$    if (coord == 0) then
!!$       if (fourier) then
!!$          call wave2phys(u)
!!$          call wave2phys(v)
!!$          call wave2phys(w)
!!$          call wave2phys(dudy)
!!$          call wave2phys(dudz)
!!$          call wave2phys(dvdx)
!!$          call wave2phys(dvdz)
!!$          call wave2phys(dwdx)
!!$          call wave2phys(dwdy)
!!$          call wave2phys(RHSx)
!!$          call wave2phys(RHSy)
!!$          call wave2phys(RHSz)
!!$       endif
!!$       print*, 'yot: >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
!!$       qq = 2
!!$       do jx=1,nx
!!$       do jy=1,ny
!!$          write(*,*) u(jx,jy,qq), v(jx,jy,qq), w(jx,jy,qq)
!!$          write(*,*) dudy(jx,jy,qq), dudz(jx,jy,qq), dvdx(jx,jy,qq)
!!$          write(*,*) dvdz(jx,jy,qq), dwdx(jx,jy,qq), dwdy(jx,jy,qq)
!!$          write(*,*) RHSx(jx,jy,qq), RHSy(jx,jy,qq), RHSz(jx,jy,qq)
!!$          print*, '---------------------------------------------'
!!$       enddo
!!$       enddo
!!$    endif

!!$    if (coord == 0) then
!!$       print*, 'TEST INTERLEAVE >>>>>>>>>>>>>>>>>>>>>>>'
!!$       call wave2phys(u)
!!$       call dfftw_execute_dft_r2c(forw, u(:,:,3), u(:,:,3) )
!!$       v = u
!!$       uc(:,:) = interleave_r2c( u(:,:,3) )
!!$       do jy=1,ny
!!$       do jx=1,ld
!!$          write(*,*) jx, jy, uc(jx,jy)
!!$       enddo
!!$       enddo
!!$       print*, '---------------------------------------'
!!$       u(:,:,3) = interleave_c2r( uc(:,:) )
!!$       do jy=1,ny
!!$       do jx=1,ld
!!$          write(*,*) jx, jy, v(jx,jy,3), u(jx,jy,3)
!!$       enddo
!!$       enddo
!!$       print*, 'END TEST INTERLEAVE <<<<<<<<<<<<<<<<<<<<'
!!$    endif

    $if ($DEBUG)
    if (DEBUG) then
        call DEBUG_write (RHSx(:, :, 1:nz), 'main.postconvec.RHSx')
    end if
    $endif

!!$    if (coord == nproc-1) then
!!$       print*, 'yot >>>>>>>>>>>'
!!$       if (fourier) then
!!$          call wave2phys(RHSx)
!!$          call wave2phys(RHSy)
!!$          call wave2phys(RHSz)
!!$          call wave2phys( w )
!!$       endif
!!$       do jx=1,ld
!!$       do jy=1,1
!!$          write(*,*) RHSy(jx,jy,1:2), RHSz(jx,jy,1:2)
!!$       enddo
!!$       enddo
!!$    endif


    ! Add div-tau term to RHS variable 
    !   this will be used for pressure calculation
    RHSx(:, :, 1:nz-1) = -RHSx(:, :, 1:nz-1) - divtx(:, :, 1:nz-1)
    RHSy(:, :, 1:nz-1) = -RHSy(:, :, 1:nz-1) - divty(:, :, 1:nz-1)
    RHSz(:, :, 1:nz-1) = -RHSz(:, :, 1:nz-1) - divtz(:, :, 1:nz-1)

!!$     if (coord == 0) then
!!$       print*, 'here >> !'
!!$       if (fourier) then
!!$          call wave2phys(RHSx)
!!$          call wave2phys(RHSy)
!!$          call wave2phys(RHSz)
!!$       endif
!!$       do jx=1,nx
!!$       do jy=1,ny
!!$          write(*,*) RHSx(jx,jy,2), RHSy(jx,jy,2), RHSz(jx,jy,2)
!!$       enddo
!!$       enddo
!!$    endif

    ! Coriolis: add forcing to RHS
    if (coriolis_forcing) then
        ! This is to put in the coriolis forcing using coriol,ug and vg as
        ! precribed in param.f90. (ug,vg) specfies the geostrophic wind vector
        ! Note that ug and vg are non-dimensional (using u_star in param.f90)
        RHSx(:, :, 1:nz-1) = RHSx(:, :, 1:nz-1) +                 &
                             coriol * v(:, :, 1:nz-1) - coriol * vg
        RHSy(:, :, 1:nz-1) = RHSy(:, :, 1:nz-1) -                 &
                             coriol * u(:, :, 1:nz-1) + coriol * ug
    end if

    !--calculate u^(*) (intermediate vel field)
    !  at this stage, p, dpdx_i are from previous time step
    !  (assumes old dpdx has NOT been added to RHSx_f, etc)
    !  we add force (mean press forcing) here so that u^(*) is as close
    !  to the final velocity as possible
    if (use_mean_p_force) then
       if (fourier) then
          !! only add to the mean (kx=0) mode
          !! no need to transform mean_p_force to kx space
          RHSx(1, 1, 1:nz-1) = RHSx(1, 1, 1:nz-1) + mean_p_force !/real(nx*ny,rprec)
       else
          RHSx(:, :, 1:nz-1) = RHSx(:, :, 1:nz-1) + mean_p_force
       endif
    end if

!!$     if (coord == 0) then
!!$       if (.not. fourier) then
!!$          call phys2wave(RHSx)
!!$          call phys2wave(RHSy)
!!$          call phys2wave(RHSz)
!!$       endif
!!$       print*, 'here >><<>><<>><<>>'
!!$       do jx=1,nx
!!$       do jy=1,1
!!$          write(*,*) jx,RHSx(jx,jy,3), RHSy(jx,jy,3), RHSz(jx,jy,3)
!!$       enddo
!!$       enddo
!!$    endif


!!$    if (coord == 0) then
!!$    ut(:,:,:) = 0._rprec
!!$    ut(1:8,1:8,:) = 1._rprec
!!$    !call dfftw_execute_dft_r2c(forw, ut, ut)
!!$    call phys2wave(ut)
!!$    print*, 'ut 1 >>>>>>>>>>>>>>>>>>>>>>'
!!$    do jx=1,8
!!$       do jy=1,8
!!$          write(*,*) jx,jy,ut(jx,jy,3)
!!$       enddo
!!$    enddo
!!$    call wave2phys(ut)
!!$    print*, 'ut 2 >>>>>>>>>>>>>>>>>>>>>>'
!!$    do jx=1,8
!!$       do jy=1,8
!!$          write(*,*) jx,jy,ut(jx,jy,3)
!!$       enddo
!!$    enddo
!!$ endif


    $if ($DEBUG)
    if (DEBUG) then
        call DEBUG_write (u(:, :, 1:nz), 'main.a.u')
        call DEBUG_write (v(:, :, 1:nz), 'main.a.v')
        call DEBUG_write (w(:, :, 1:nz), 'main.a.w')
        call DEBUG_write (RHSx(:, :, 1:nz), 'main.a.RHSx')
        call DEBUG_write (RHSy(:, :, 1:nz), 'main.a.RHSy')
        call DEBUG_write (RHSz(:, :, 1:nz), 'main.a.RHSz')
        call DEBUG_write (RHSz_f(:, :, 1:nz), 'main.a.RHSx_f')
        call DEBUG_write (RHSz_f(:, :, 1:nz), 'main.a.RHSy_f')
        call DEBUG_write (RHSz_f(:, :, 1:nz), 'main.a.RHSz_f')
        call DEBUG_write (dpdx(:, :, 1:nz), 'main.a.dpdx')
        call DEBUG_write (dpdy(:, :, 1:nz), 'main.a.dpdy')
        call DEBUG_write (dpdz(:, :, 1:nz), 'main.a.dpdz')
        call DEBUG_write (fxa(:, :, 1:nz), 'main.a.fxa')
        call DEBUG_write (force, 'main.a.force')
    end if
    $endif

    !//////////////////////////////////////////////////////
    !/// APPLIED FORCES                                 ///
    !//////////////////////////////////////////////////////
    !  In order to save memory the arrays fxa, fya, and fza are now only defined when needed.
    !  For Levelset RNS all three arrays are assigned. 
    !  For turbines at the moment only fxa is assigned.
    !  Look in forcing_applied for calculation of forces.
    !  Look in sim_param.f90 for the assignment of the arrays.
        
    !  Applied forcing (forces are added to RHS{x,y,z})
    
    call forcing_applied()

    !  Update RHS with applied forcing
    $if ($TURBINES and not ($LVLSET and $RNS_LS))
    RHSx(:,:,1:nz-1) = RHSx(:,:,1:nz-1) + fxa(:,:,1:nz-1)
    $elseif ($LVLSET and $RNS_LS)
    RHSx(:,:,1:nz-1) = RHSx(:,:,1:nz-1) + fxa(:,:,1:nz-1)
    RHSy(:,:,1:nz-1) = RHSy(:,:,1:nz-1) + fya(:,:,1:nz-1)
    RHSz(:,:,1:nz-1) = RHSz(:,:,1:nz-1) + fza(:,:,1:nz-1)    
    $endif

    !//////////////////////////////////////////////////////
    !/// EULER INTEGRATION CHECK                        ///
    !////////////////////////////////////////////////////// 
    ! Set RHS*_f if necessary (first timestep) 
    if ((jt_total == 1) .and. (.not. initu)) then
      ! if initu, then this is read from the initialization file
      ! else for the first step put RHS_f=RHS
      !--i.e. at first step, take an Euler step
      RHSx_f=RHSx
      RHSy_f=RHSy
      RHSz_f=RHSz
    end if    

    !//////////////////////////////////////////////////////
    !/// INTERMEDIATE VELOCITY                          ///
    !//////////////////////////////////////////////////////     
    ! Calculate intermediate velocity field
    !   only 1:nz-1 are valid
    u(:, :, 1:nz-1) = u(:, :, 1:nz-1) +                   &
                     dt * ( tadv1 * RHSx(:, :, 1:nz-1) +  &
                            tadv2 * RHSx_f(:, :, 1:nz-1) )
    v(:, :, 1:nz-1) = v(:, :, 1:nz-1) +                   &
                     dt * ( tadv1 * RHSy(:, :, 1:nz-1) +  &
                            tadv2 * RHSy_f(:, :, 1:nz-1) )
    w(:, :, 1:nz-1) = w(:, :, 1:nz-1) +                   &
                     dt * ( tadv1 * RHSz(:, :, 1:nz-1) +  &
                            tadv2 * RHSz_f(:, :, 1:nz-1) )

!!$    if (coord == 0) then
!!$       if (.not. fourier) then
!!$          call phys2wave(RHSx)
!!$          call phys2wave(RHSy)
!!$          call phys2wave(RHSz)
!!$       endif
!!$       print*, 'here a >>>>>>>>>>>'
!!$       do jx=1,ld
!!$       do jy=1,ny
!!$          write(*,*) jx,jy, RHSx(jx,jy,1), RHSy(jx,jy,3), RHSz(jx,jy,nz-1)
!!$       enddo
!!$       enddo
!!$    endif

    ! Set unused values to BOGUS so unintended uses will be noticable
    $if ($SAFETYMODE)
    $if ($MPI)
        u(:, :, 0) = BOGUS
        v(:, :, 0) = BOGUS
        w(:, :, 0) = BOGUS
    $endif
    !--this is an experiment
    !--u, v, w at jz = nz are not useful either, except possibly w(nz), but that
    !  is supposed to zero anyway?
    !--this has to do with what bc are imposed on intermediate velocity    
    u(:, :, nz) = BOGUS
    v(:, :, nz) = BOGUS
    w(:, :, nz) = BOGUS
    $endif
    
    $if ($DEBUG)
    if (DEBUG) then
        call DEBUG_write (u(:, :, 1:nz), 'main.b.u')
        call DEBUG_write (v(:, :, 1:nz), 'main.b.v')
        call DEBUG_write (w(:, :, 1:nz), 'main.b.w')
    end if
    $endif

    !//////////////////////////////////////////////////////
    !/// PRESSURE SOLUTION                              ///
    !//////////////////////////////////////////////////////
    ! Solve Poisson equation for pressure
    !   div of momentum eqn + continuity (div-vel=0) yields Poisson eqn
    !   do not need to store p --> only need gradient
    !   provides p, dpdx, dpdy at 0:nz-1 and dpdz at 1:nz-1
    
!!$    if (fourier) then
!!$      call wave2phys(u)
!!$      call wave2phys(v)
!!$      call wave2phys(w)
!!$      call wave2phys(divtz)
!!$    endif

!!$    if (coord == 0) then
!!$       print*, 'u >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
!!$       do jx=1,ld
!!$       do jy=1,ny
!!$          write(*,*) jx,jy, u(jx,jy,nz-1:nz)
!!$       enddo
!!$       enddo
!!$       print*, 'v >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
!!$       do jx=1,ld
!!$       do jy=1,ny
!!$          write(*,*) jx,jy, v(jx,jy,nz-1:nz)
!!$       enddo
!!$       enddo
!!$       print*, 'w >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
!!$       do jx=1,ld
!!$       do jy=1,ny
!!$          write(*,*) jx,jy, w(jx,jy,nz-1:nz)
!!$       enddo
!!$       enddo
!!$       print*, 'divtz >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
!!$       do jx=1,ld
!!$       do jy=1,ny
!!$          write(*,*) jx,jy, divtz(jx,jy,nz-1:nz)
!!$       enddo
!!$       enddo
!!$
!!$    endif

    if (fourier) then
       if (coord == 0) then
          u(3:4,:,:) = 0._rprec
          v(3:4,:,:) = 0._rprec
          w(3:4,:,:) = 0._rprec
       endif
    endif

    call press_stag_array()

!!$    if (fourier) then
!!$      call phys2wave(u)
!!$      call phys2wave(v)
!!$      call phys2wave(w)
!!$      call phys2wave(divtz)
!!$      call phys2wave_pr(dpdx)
!!$      call phys2wave_pr(dpdy)
!!$      call phys2wave_pr(dpdz)
!!$    endif


!!$    if (coord == nproc-1) then
!!$       print*, 'yot >>>>>>>>>>>'
!!$       if (fourier) then
!!$          call wave2phys(RHSy)
!!$          call wave2phys(RHSz)
!!$       endif
!!$       do jx=1,ld
!!$       do jy=1,ny
!!$          write(*,*) RHSy(jx,jy,1:2), RHSz(jx,jy,1:2)
!!$       enddo
!!$       enddo
!!$    endif


    ! Add pressure gradients to RHS variables (for next time step)
    !   could avoid storing pressure gradients - add directly to RHS
    RHSx(:, :, 1:nz-1) = RHSx(:, :, 1:nz-1) - dpdx(:, :, 1:nz-1)
    RHSy(:, :, 1:nz-1) = RHSy(:, :, 1:nz-1) - dpdy(:, :, 1:nz-1)
    RHSz(:, :, 1:nz-1) = RHSz(:, :, 1:nz-1) - dpdz(:, :, 1:nz-1)

!!$    if (coord == 0) then
!!$       if (.not. fourier) then
!!$          call phys2wave_pr(dpdx)
!!$          call phys2wave_pr(dpdy)
!!$          call phys2wave_pr(dpdz)
!!$          !call phys2wave(RHSx)
!!$          !call phys2wave(RHSy)
!!$          !call phys2wave(RHSz)
!!$          call phys2wave(u)
!!$          call phys2wave(v)
!!$          call phys2wave(w)
!!$       endif
!!$       do jx=1,ld
!!$       do jy=1,ny
!!$          write(*,*) jx,jy, dpdx(jx,jy,1), dpdy(jx,jy,1), dpdz(jx,jy,3)
!!$       enddo
!!$       enddo
!!$       print*, '^^^ pressure ^^^ >>>> RHS >>>>'
!!$       do jx=1,ld
!!$       do jy=1,ny
!!$          !write(*,*) jx,jy, RHSx(jx,jy,3), RHSy(jx,jy,1), RHSz(jx,jy,nz-1)
!!$          write(*,*) jx,jy, u(jx,jy,3), v(jx,jy,1), w(jx,jy,nz-1)
!!$       enddo
!!$       enddo
!!$    endif

    $if ($DEBUG)
    if (DEBUG) then
        call DEBUG_write (dpdx(:, :, 1:nz), 'main.dpdx')
        call DEBUG_write (dpdy(:, :, 1:nz), 'main.dpdy')
        call DEBUG_write (dpdz(:, :, 1:nz), 'main.dpdz')
    end if
    $endif

    $if ($DEBUG)
    if (DEBUG) then
        call DEBUG_write (u(:, :, 1:nz), 'main.d.u')
        call DEBUG_write (v(:, :, 1:nz), 'main.d.v')
        call DEBUG_write (w(:, :, 1:nz), 'main.d.w')
    end if
    $endif

    !//////////////////////////////////////////////////////
    !/// INDUCED FORCES                                 ///
    !//////////////////////////////////////////////////////    
    ! Calculate external forces induced forces. These are
    ! stored in fx,fy,fz arrays. We are calling induced 
    ! forces before applied forces as some of the applied
    ! forces (RNS) depend on the induced forces and the 
    ! two are assumed independent
    call forcing_induced()

    !//////////////////////////////////////////////////////
    !/// PROJECTION STEP                                ///
    !//////////////////////////////////////////////////////   
    ! Projection method provides u,v,w for jz=1:nz
    !   uses fx,fy,fz calculated above
    !   for MPI: syncs 1 -> Nz and Nz-1 -> 0 nodes info for u,v,w    
    call project ()

!!$    if (coord == 0) then
!!$       if (.not. fourier) then
!!$          call phys2wave(u)
!!$          call phys2wave(v)
!!$          call phys2wave(w)
!!$       endif
!!$       print*, 'here a >>>>>>>>>>>'
!!$       do jx=1,ld
!!$       do jy=1,ny
!!$          write(*,*) jx,jy, u(jx,jy,1), v(jx,jy,3), w(jx,jy,nz-1)
!!$       enddo
!!$       enddo
!!$    endif


    $if($LVLSET and $RNS_LS)
    !  Compute the relavent force information ( include reference quantities, CD, etc.)
    !  of the RNS elements using the IBM force; No modification to f{x,y,z} is
    !  made here.
    call rns_elem_force_ls()
    $endif
   
    ! Write ke to file
    if (modulo (jt_total, nenergy) == 0) call energy (ke)

    $if ($LVLSET)
      if( global_CA_calc ) call level_set_global_CA ()
    $endif


!!$u = 0._rprec
!!$do jx=1,nx
!!$do jy=1,ny
!!$do jz=1,nz
!!$u(jx,jy,jz) = 2.7         + 2.1*sin(L_x/nx*(jx-1)*1.0)
!!$u(jx,jy,jz) = u(jx,jy,jz) + 1.1*sin(L_y/ny*(jy-1)*2.0)
!!$u(jx,jy,jz) = u(jx,jy,jz) + 1.7*sin(2*3.1415926535897*L_z/(nz)*(jz-1)*3.0)
!!$enddo
!!$enddo
!!$enddo
!!$
!!$call phys2wave( u )
    
    ! Write output files
    !! Currently transforms u,v,w,txz to physical space for recording tavg
    !! Stats of other quantities not currently available
    !! Separately transforms u,v,w for instantaneous recording in inst_write subroutine
    !! Also handles checkpointing

    call output_loop()  

    ! Check the total time of the simulation up to this point on the master node and send this to all
   
    if (modulo (jt_total, wbase) == 0) then

       if ( fourier ) then   !!jb
         call wave2physF( u, uF )
         call wave2physF( v, vF )
         call wave2physF( w, wF )
         call wave2physF( txz, txzF )
         call wave2physF( tyz, tyzF )
!!$         call wave2phys( u )
!!$         call wave2phys( v )
!!$         call wave2phys( w )
!!$         call wave2phys( txz )
!!$         call wave2phys( tyz )
!!$         uF = u
!!$         vF = v
!!$         wF = w
!!$         txzF = txz  !! for tau_wall recording
!!$         tyzF = tyz
!!$         call phys2wave( u )
!!$         call phys2wave( v )
!!$         call phys2wave( w )
!!$         call phys2wave( txz )
!!$         call phys2wave( tyz )
       endif
       
       ! Get the ending time for the iteration
       call clock_stop( clock )
       call clock_stop( clock_total )

       ! Calculate rms divergence of velocity
       ! only written to screen, not used otherwise
       call rmsdiv (rmsdivvel)
       maxcfl = get_max_cfl()
       if (coord == 0) then
          write(*,*)
          write(*,'(a)') '========================================'
          write(*,'(a)') 'Time step information:'
          write(*,'(a,i9)') '  Iteration: ', jt_total
          write(*,'(a,E15.7)') '  Time step: ', dt
          write(*,'(a,E15.7)') '  CFL: ', maxcfl
          write(*,'(a,2E15.7)') '  AB2 TADV1, TADV2: ', tadv1, tadv2
          write(*,*) 
          write(*,'(a)') 'Flow field information:'          
          write(*,'(a,E15.7)') '  Velocity divergence metric: ', rmsdivvel
          write(*,'(a,E15.7)') '  Kinetic energy: ', ke
          write(*,'(a,E15.7)') '  Wall stress: ', get_tau_wall()
          write(*,*)
          $if($MPI)
          write(*,'(1a)') 'Simulation wall times (s): '
          $else
          write(*,'(1a)') 'Simulation cpu times (s): '
          $endif
          write(*,'(1a,E15.7)') '  Iteration: ', clock % time
          write(*,'(1a,E15.7)') '  Cumulative: ', clock_total % time
          write(*,'(a)') '========================================'
       end if
       if (fourier) then
          if (coord == 0) then            !!jb
             write(*,'(a)') '===================== BOTTOM WALL =================='
             write(*,*) 'u: ', uF(1,1,1:2)
             write(*,*) 'v: ', vF(1,1,1:2)
             write(*,*) 'w: ', wF(1,1,1:2)
             write(*,'(a)') '===================================================='
          endif
          call mpi_barrier(comm, ierr)
          if (coord == nproc-1) then            !!jb
             write(*,'(a)') '=============== TOP ===================='
             write(*,*) 'u: ', uF(1,1,nz-1)
             write(*,*) 'v: ', vF(1,1,nz-1)
             write(*,*) 'w: ', wF(1,1,nz-1), wF(1,1,nz)
             write(*,'(a)') '======================================='
          endif
          call mpi_barrier(comm, ierr)
       else
          if (coord == 0) then            !!jb
             write(*,'(a)') '===================== BOTTOM WALL =================='
             write(*,*) 'u: ', u(1,1,1:2)
             write(*,*) 'v: ', v(1,1,1:2)
             write(*,*) 'w: ', w(1,1,1:2)
             write(*,'(a)') '===================================================='
          endif
          call mpi_barrier(comm, ierr)
          if (coord == nproc-1) then            !!jb
             write(*,'(a)') '=============== TOP ===================='
             write(*,*) 'u: ', u(1,1,nz-1)
             write(*,*) 'v: ', v(1,1,nz-1)
             write(*,*) 'w: ', w(1,1,nz-1), w(1,1,nz)
             write(*,'(a)') '======================================='
          endif
          call mpi_barrier(comm, ierr)
       endif
       if (coord == 0) call write_tau_wall()     !!jb
       if (fourier) then
          call energy_kx_spectral_complex_fourier(uF,vF,wF)    !!jb
       else
          call energy_kx_spectral_complex(u,v,w)    !!jb
       endif
       ! Check if we are to check the allowable runtime
       if( runtime > 0 ) then

          ! Determine the processor that has used most time and communicate this.
          ! Needed to make sure that all processors stop at the same time and not just some of them
          $if($MPI)
          call mpi_allreduce(clock_total % time, rbuffer, 1, MPI_RPREC, MPI_MAX, MPI_COMM_WORLD, ierr)
          clock_total % time = rbuffer
          $endif
       
          ! If maximum time is surpassed go to the end of the program
          if ( clock_total % time >= real(runtime,rprec) ) then
             call mesg( prog_name, 'Specified runtime exceeded. Exiting simulation.')
             exit time_loop
          endif

       endif

!!$       if ( fourier ) then   !!jb
!!$         call phys2waveF(uF)
!!$         call phys2waveF(vF)
!!$         call phys2waveF(wF)
!!$         call phys2waveF(txzF)
!!$         call phys2waveF(tyzF, tyz)
!!$         !call phys2wave(dudx)
!!$         !call phys2wave(dvdy)
!!$         !call phys2wave(dwdz)
!!$      endif

   end if

end do time_loop
! END TIME LOOP
call cpu_time(jtime2)    !!jb

! Finalize
close(2)
    
! Write total_time.dat and tavg files
! calls checkpoint, which already handle phys2wave transforms
call output_final()

! Stop wall clock
call clock_stop( clock_total )
$if($MPI)
if( coord == 0 )  write(*,"(a,e15.7)") 'Simulation wall time (s) : ', clock_total % time
$else
if( coord == 0 )  write(*,"(a,e15.7)") 'Simulation cpu time (s) : ', clock_total % time
$endif

call finalize()

if(coord == 0) write(*,*) 'jT: ', jtime2-jtime1
if(coord == 0) write(*,'(a)') 'Simulation complete'

end program main
