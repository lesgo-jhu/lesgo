subroutine press_stag(p_hat,dfdx,dfdy)   
! p_hat contains the physical space pressure on exit
!-------------------    
! Boundary Layer version with 4th order derivs in vertical.
!  04 December 1995
!	Mods.
!	12/6: added upper and lower boundary conditions.	
!    12/8: Corrected sign on imag. x and y deriv of RHS.
!	12/17 Forcing pressure equal zero at wall 
!			Removed forcing of <P>=0 for all z.
!		    Will need to change BC when hetero. surface stress.
!	12/17 Added ficticious node below wall (Note solver on Nz+1 x Nz+1
!          prev. version of this subroutine saved as p_bl4.old.12.17
!    12/18 Revised 2st deriv stencil at wall (avg of deriv at -1 and at 1)
!    12/21 Redid FDD to 2nd order accurate.
!    12/22 Now value of P(wall) diagnosed from prev P, to match gradient BC
!....1/13: major changes.
!		Broke out mean pressure for separate solution
!    1/20 back to tridag for natrix solution (same sol'n as LUDCMP...A Keeper!)
!....1/23 Staggered solution
!.........Using Nz+1 levels for computing P, but tossing out level below ground
!....4/1 Changed sign on Div T_iz at wall and lid (five places)
!-------------------          
use types,only:rprec
use param
use sim_param,only:u,v,w,RHSx,RHSy,RHSz,RHSx_f,RHSy_f,RHSz_f, divtz
use fft
use immersedbc,only:fx,fy,fz  ! only for forcing experiment
use debug_mod
!$undefine $MPI
$if ($MPI)
  use mpi_transpose_mod
$endif
implicit none      
complex(kind=rprec),dimension(lh,ny,0:nz)::p_hat
real(kind=rprec),dimension(ld,ny,nz)::rH_x,rH_y,rH_z
complex(kind=rprec),dimension(lh,ny,nz)::H_x,H_y,H_z
equivalence (rH_x,H_x),(rH_y,H_y),(rH_z,H_z)
real(kind=rprec),dimension(ld,ny)::rtopw, rbottomw
complex(kind=rprec),dimension(lh,ny)::topw,bottomw
equivalence (rtopw,topw),(rbottomw,bottomw)
complex(kind=rprec),dimension(lh,ny,nz),intent(out)::dfdx,dfdy
real(kind=rprec)::const,ignore_me
! remove this stuff!
integer::jx,jy,jz,k

character (64) :: fname
logical, parameter :: DEBUG = .true.

$if ($MPI)
  !--for now, we will do the transpose stuff out-of-place,
  !  but it should be possible to do in place (perhaps use equivalences)
  !--transpose of complex NOT just transpose of its equivalenced
  !  real array
  integer, parameter :: nz_t = (nz-1) * nproc + 1
  integer, parameter :: nx_t = nx / nproc
  !integer, save :: rank_of_coord(0:nproc-1) = -1
  integer :: ip
  integer :: tag

  logical, save :: init = .false.

  complex (rprec) :: temp_t(ny, nx_t/2), temp(nx_t/2, ny)
  complex (rprec), dimension (nz_t, ny, nx_t/2) :: H_x_t, H_y_t, H_z_t
                              !--local transposes                              
                              !--problem: sizes are not known till run time
                              !--only know total size of array
                              !--no nyquist in 3rd dim (x)
  complex (rprec), dimension (0:nz_t, ny, nx_t/2) :: p_hat_t
  complex (rprec), dimension (ny, nx_t/2) :: bottomw_t, topw_t 

  complex (rprec), dimension (nz_t+1) :: RHS_col, p_col
  real (rprec), dimension (nz_t+1) :: a1, b1, c1

  real (rprec) :: kx_t, ky_t
$else
  complex(kind=rprec),dimension(nz+1)::RHS_col,p_col
  real(kind=rprec),dimension(nz+1)::a1,b1,c1
$endif

if (DEBUG) write (*, *) $str($context_doc), coord, ' started press_stag'

!$if ($MPI)
!  if (.not. init) then
!    !--initialization of rank_of_coord
!    do ip = 0, nproc-1
!      call mpi_cart_rank (comm, (/ ip /), rank_of_coord(ip), ierr)
!    end do
!
!    init = .true.
!  end if
!$endif

p_hat(:,:,0)=0._rprec ! dont know for sure if this is right

$if ($MPI)
  p_hat_t(0, :, :) = 0._rprec
$endif

!==========================================================================
! Get the right hand side ready 
! Loop over levels    
const=1._rprec/(nx*ny)
do jz=1,nz
! temp storage for sum of RHS terms.  normalized for fft
! sc: recall that the old timestep guys already contain the pressure
!   term

   ! no forces
   rH_x(:, :, jz) = const / tadv1 * (u(:, :, jz) / dt)
   rH_y(:, :, jz) = const / tadv1 * (v(:, :, jz) / dt)
   rH_z(:, :, jz) = const / tadv1 * (w(:, :, jz) / dt)

   call rfftwnd_f77_one_real_to_complex(forw,rH_x(:,:,jz),ignore_me)
   call rfftwnd_f77_one_real_to_complex(forw,rH_y(:,:,jz),ignore_me)
   call rfftwnd_f77_one_real_to_complex(forw,rH_z(:,:,jz),ignore_me)     
end do

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
  ! sc: could do out of place transform if wanted to...
  rbottomw(:,:)=const*divtz(:,:,1)
  call rfftwnd_f77_one_real_to_complex(forw,rbottomw(:,:),ignore_me)
end if
if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
  rtopw(:,:)=const*divtz(:,:,nz)
  call rfftwnd_f77_one_real_to_complex(forw,rtopw(:,:),ignore_me)
end if

if (DEBUG) then
  write (*, *) 'divtz:'
  write (*, *) '(1, 1, 1): ', divtz(1, 1, 1)
  write (*, *) '(1, 1, 2): ', divtz(1, 1, 2)
  write (*, *) '(2, 1, 1): ', divtz(2, 1, 1)
  write (fname, '(a,i0)') 'debug.divtz.', jt
  $if ($MPI)
    write (fname, '(a)') trim (fname) // $str(.MPI)
  $endif
  open (1, file=fname)
  do jz = 1, nz
    do jy = 1, ny
      do jx = 1, lh
        write (1, *) jx, jy, jz, divtz(jx, jy, jz)
      end do
    end do
  end do
  close (1)
end if

! set oddballs to 0
! probably can get rid of this if we're more careful below
H_x(lh,:,:)=0._rprec
H_y(lh,:,:)=0._rprec
H_z(lh,:,:)=0._rprec
H_x(:,ny/2+1,:)=0._rprec
H_y(:,ny/2+1,:)=0._rprec
H_z(:,ny/2+1,:)=0._rprec
!--with MPI; topw and bottomw are only on top & bottom processes
topw(lh,:)=0._rprec
topw(:,ny/2+1)=0._rprec
bottomw(lh,:)=0._rprec
bottomw(:,ny/2+1)=0._rprec

if (DEBUG) then
  write (*, *) 'bottomw:'
  write (*, *) '(1, 1): ', bottomw(1, 1)
  write (*, *) '(2, 2): ', bottomw(2, 2)
  write (*, *) '(3, 2): ', bottomw(3, 2)
  write (fname, '(a,i0)') 'debug.bottomw.', jt
  $if ($MPI)
    write (fname, '(a)') trim (fname) // $str(.MPI)
  $endif
  open (1, file=fname)
  do jy = 1, ny
    do jx = 1, lh
      write (1, *) jx, jy, bottomw(jx, jy)
    end do
  end do
  close (1)

  write (*, *) 'topw:'
  write (*, *) '(1, 1): ', topw(1, 1)
  write (*, *) '(2, 2): ', topw(2, 2)
  write (*, *) '(3, 2): ', topw(3, 2)
  write (fname, '(a,i0)') 'debug.topw.', jt
  $if ($MPI)
    write (fname, '(a)') trim (fname) // $str(.MPI)
  $endif
  open (1, file=fname)
  do jy = 1, ny
    do jx = 1, lh
      write (1, *) jx, jy, topw(jx, jy)
    end do
  end do
  close (1)
end if

if (DEBUG) then
  write (*, *) 'H_x:'
  write (*, *) '(1, 1, 1): ', H_x(1, 1, 1)
  write (*, *) '(2, 1, 1): ', H_x(2, 1, 1)
  write (*, *) '(1, 1, 2): ', H_x(1, 1, 2)
  write (fname, '(a,i0)') 'debug.H_x_hat.', jt
  $if ($MPI)
    write (fname, '(a)') trim (fname) // $str(.MPI)
  $endif
  open (1, file=fname)
  do jz = 1, nz
    do jy = 1, ny
      do jx = 1, lh
        write (1, *) jx, jy, jz, H_x(jx, jy, jz)
      end do
    end do
  end do
  close (1)

  write (*, *) 'H_y:'
  write (*, *) '(1, 1, 1): ', H_y(1, 1, 1)
  write (*, *) '(2, 1, 1): ', H_y(2, 1, 1)
  write (*, *) '(1, 1, 2): ', H_y(1, 1, 2)
  write (fname, '(a,i0)') 'debug.H_y_hat.', jt
  $if ($MPI)
    write (fname, '(a)') trim (fname) // $str(.MPI)
  $endif
  open (1, file=fname)
  do jz = 1, nz
    do jy = 1, ny
      do jx = 1, lh
        write (1, *) jx, jy, jz, H_y(jx, jy, jz)
      end do
    end do
  end do
  close (1)

  write (*, *) 'H_z:'
  write (*, *) '(1, 1, 1): ', H_z(1, 1, 1)
  write (*, *) '(2, 1, 1): ', H_z(2, 1, 1)
  write (*, *) '(1, 1, 2): ', H_z(1, 1, 2)
  write (fname, '(a,i0)') 'debug.H_z_hat.', jt
  $if ($MPI)
    write (fname, '(a)') trim (fname) // $str(.MPI)
  $endif
  open (1, file=fname)
  do jz = 1, nz
    do jy = 1, ny
      do jx = 1, lh
        write (1, *) jx, jy, jz, H_z(jx, jy, jz)
      end do
    end do
  end do
  close (1)
end if


if (DEBUG) write (*, *) $str($context_doc), coord, ' before mpi-if'

$if ($MPI)
  !--I think this is where we do the transpose
  !--transpose will change alot (e.g. which indices are looped over),
  !  so its probably better to keep MPI & non-MPI versions separate
  !  basically from here down

 
  !--transpose
  !--could perhaps do in place, but that would require that we reshape
  !  the array being transposed, and need to be careful of aliasing
  !--the transposed arrays do not contain the x Nyquist freqs
  !--nz-1 here:  the only process that needs to transfer nz information
  !  is the top process, and this is done separately
  !--H_{x,y,z}_t is dimensioned (nproc*(nz-1) + 1) X ny X (nx/2/nproc), the
  !  transpose will only fill nproc*(nz-1) X ny X (nx/2/nproc)  
  call mpi_transpose (nx/2, ny, nz-1, H_x, H_x_t)
  call mpi_transpose (nx/2, ny, nz-1, H_y, H_y_t)
  call mpi_transpose (nx/2, ny, nz-1, H_z, H_z_t)
  !--transfer nz line from top process
  call mpi_transpose_top_line (H_x, H_x_t)
  call mpi_transpose_top_line (H_y, H_y_t)
  call mpi_transpose_top_line (H_z, H_z_t)

  if (DEBUG) write (*, *) $str($context_doc), coord, ' after mpi_transpose H'

  !--transpose topw, bottomw (different since they are not "full" arrays)
  !--only coordinates 0, nproc-1 take part in sending, all procs recv
  !--really should use scatterv, gatherv here
  if (coord == 0) then  !--send bottomw

    do ip = 1, nproc-1
      !--use  bottomw_t as a buffer for chunk ip
      bottomw_t = transpose (bottomw(ip*nx_t/2 + 1 : (ip+1)*nx_t/2, :))
      !--send chunk ip to proc with coord ip
      call mpi_send (bottomw_t(1, 1), ny*nx_t/2, MPI_CPREC,  &
                     rank_of_coord(ip), ip, comm, ierr)
    end do

    !--now transpose chunk 0 (local)
    bottomw_t = transpose (bottomw(1:nx_t/2, :))

  else  !--recv chunk of bottomw_t from proc with coord 0   

    call mpi_recv (bottomw_t(1, 1), ny*nx_t/2, MPI_CPREC,     &
                   rank_of_coord(0), coord, comm, status, ierr)

  end if

  if (DEBUG) write (*, *) $str($context_doc), coord, ' after bottomw'

  if (DEBUG) then
    write (*, *) $str($context_doc), coord, 'bottomw_t(1,1) = ', bottomw_t(1, 1)
  end if

  if (coord == nproc-1) then  !--send topw

    do ip = 0, nproc-2
      !--use topw_t as a buffer for chunk ip
      topw_t = transpose (topw(ip*nx_t/2 + 1 : (ip+1)*nx_t/2, :))
      !--send chunk to proc with coord ip
      call mpi_send (topw_t(1, 1), ny*nx_t/2, MPI_CPREC,  &
                     rank_of_coord(ip), ip, comm, ierr)
    end do

    !--now transpose chunk nproc-1
    topw_t = transpose (topw((nproc-1)*nx_t/2 + 1 : nproc*nx_t/2, :))

  else  !--recv chunk of topw from proc with coord nproc-1

    call mpi_recv (topw_t(1, 1), ny*nx_t/2, MPI_CPREC,              &
                   rank_of_coord(nproc-1), coord, comm, status, ierr)  

  end if

  if (DEBUG) write (*, *) $str($context_doc), coord, ' after topw'

  do jx = 1, nx_t/2  !--jx is now third dimension

    kx_t = (2._rprec * pi / L_x) * ((jx-1) + coord * nx_t/2)

    !--the jx=lh part is skipped (all amplitudes zero)
 
    do jy = 1, ny

      if (jy == ny/2 + 1) then
        p_hat_t(:, jy, jx) = 0._rprec
        cycle  !--jy-loop
      end if

      ky_t = (2._rprec * pi / L_y) * (modulo (jy - 1 + ny/2, ny) - ny/2)

      if (jx * jy == 1) then

        !--zero wavenumber solution 
        p_hat_t(1, jy, jx) = p_hat_t(0, jy, jx) - dz * bottomw_t(jy, jx)
        do jz = 2, nz_t
          p_hat_t(jz, jy, jx) = p_hat_t(jz-1, jy, jx) + H_z_t(jz, jy, jx) * dz
        end do

        if (DEBUG) then
          write (*, *) '0-wave ', coord, ' p_hat_t(1) =', p_hat_t(1, 1, 1)
          write (*, *) '0-wave ', coord, ' p_hat_t(0) =', p_hat_t(0, 1, 1)
          write (*, *) '0-wave ', coord, ' jx, jy =', jx, jy
          write (*, *) '0-wave ', coord, ' dz =', dz
        end if

      else

        !--set up tridiagonal matrix solution for pressure Fourier coeffs
        !--near wall nodes
        a1(1) = 0._rprec
        b1(1) = -1._rprec
        c1(1) = 1._rprec
        RHS_col(1) = -dz * bottomw_t(jy, jx) ! complex 

        !--interior nodes
        !--plan to store this rhs, then try to do only one "forward" transpose
        do jz = 2, nz_t
          RHS_col(jz) = eye * (kx_t * H_x_t(jz-1, jy, jx) +          &
                               ky_t * H_y_t(jz-1, jy, jx)) +         &
                        (H_z_t(jz, jy, jx) - H_z_t(jz-1, jy, jx)) / dz
          a1(jz) = 1._rprec / (dz**2)
          b1(jz) = -(kx_t**2 + ky_t**2 + 2._rprec / (dz**2))
          c1(jz) = 1._rprec / (dz**2)
        end do

        !--top nodes 
        RHS_col(nz_t+1) = -dz * topw_t(jy, jx)
        a1(nz_t+1)= -1._rprec
        b1(nz_t+1) = 1._rprec
        c1(nz_t+1) = 0._rprec

        !if (DEBUG) then
        !  write (*, *) 'jx, jy = ', jx, jy
        !  write (*, *) 'a1 = ', a1
        !  write (*, *) 'b1 = ', b1
        !  write (*, *) 'c1 = ', c1
        !  write (*, *) 'RHS_col = ', RHS_col
        !end if

        !--tridiagonal solver
        !--modified to handle complex
        call tridag (a1, b1, c1, RHS_col, p_col, nz_t+1)

        !--copy p_col back to p_hat array
        !--since jz indices occupy consective storage in p_hat_t
        !  we could just call tridag with p_hat_t(:, jy, jx) and avoid this
        !  (except for that shift of indices...)
        do jz = 0, nz_t
          p_hat_t(jz,jy,jx) = p_col(jz+1)
        end do

      end if  !--end if from zero wavenumber solution

    end do  !--end jy loop
  end do  !--end jx loop

  if (DEBUG) write (*, *) $str($context_doc), coord, ' after jy,jx loops'

  !--transpose p_hat back
  !--do not need to transpose back H or bottomw or topw
  !--this will NOT do p_hat_t(0,1:ny,1:nx_t/2) -> p_hat(1:nx/2, 1:ny, 0),
  !  so if this is not zero, then it must also be transmitted 
  !--this will NOT do p_hat_t(nz_t, :, :) -> p_hat(:, :, nz) (top proc only)
  call mpi_transpose (nz_t-1, ny, nx_t/2, p_hat_t(1:nz_t, 1:ny, 1:nx_t/2),  &
                      p_hat(1:nx/2, 1:ny, 1:nz))

  !--now transfer p_hat_t(0,1:ny,1:nx_t/2) -> p_hat(1:nx/2, 1:ny, 0)
  !  @ bottom process only
  if (coord == 0) then

    !--probably much better ways to do this
    !  nonblocking?, also move local transpose to sender?
    do ip = 1, nproc-1
      !--recv chunk of size ny*nx_t/2 from coord ip
      call mpi_recv (temp_t(1, 1), ny*nx_t/2, MPI_CPREC, rank_of_coord(ip),  &
                     ip, comm, status, ierr)
      !--copy chunk into position
      p_hat(ip*nx_t/2 + 1 : (ip+1)*nx_t/2, 1:ny, 0) = transpose (temp_t)
    end do

    !--now transpose chunk 0
    p_hat(1:nx_t/2, 1:ny, 0) = transpose (p_hat_t(0, 1:ny, 1:nx_t/2))

  else

    !--send chunks of size ny*nx_t/2 to coord 0
    temp_t = p_hat_t(0, :, :)
    call mpi_send (temp_t(1, 1), ny*nx_t/2, MPI_CPREC, rank_of_coord(0),  &
                   coord, comm, ierr)

  end if

  !--transfer p_hat_t(nz_t, :, :) -> p_hat(1:, :, nz) @ top process
  !--not sure if its better to have sender to local transpose
  if (coord == nproc-1) then

    do ip = 0, nproc-2
      call mpi_recv (temp(1, 1), ny*nx_t/2, MPI_CPREC,        &
                     rank_of_coord(ip), ip, comm, status, ierr)
      p_hat(ip*nx_t/2 + 1:(ip+1)*nx_t/2, 1:ny, nz) = temp
    end do

    !--now transpose chunk nproc-1
    p_hat((nproc-1)*nx_t/2 + 1:nproc*nx_t/2, 1:ny, nz) =              &
                              transpose (p_hat_t(nz_t, 1:ny, 1:nx_t/2))

  else  !--send chunk to coord nproc-1

    temp = transpose (p_hat_t(nz_t, 1:ny, 1:nx_t/2))
    call mpi_send (temp(1, 1), ny*nx_t/2, MPI_CPREC,        &
                   rank_of_coord(nproc-1), coord, comm, ierr)
                   
  end if

  if (DEBUG) then
    write (*, *) coord, 'pre: p_hat(3,2,0)=', p_hat(3,2,0)
    write (*, *) coord, 'pre: p_hat(3,2,1)=', p_hat(3,2,1)
    write (*, *) coord, 'pre: p_hat(3,2,nz-1)=', p_hat(3,2,nz-1)
    write (*, *) coord, 'pre: p_hat(3,2,nz)=', p_hat(3,2,nz)
  end if

  !--reintroduce overlap at nz-1 <-> 0, nz <-> 1 here
  call mpi_sendrecv (p_hat(1, 1, nz-1), lh*ny, MPI_CPREC, up, 1,  &
                     p_hat(1, 1, 0), lh*ny, MPI_CPREC, down, 1,   &
                     comm, status, ierr)
  call mpi_sendrecv (p_hat(1, 1, 1), lh*ny, MPI_CPREC, down, 2, &
                     p_hat(1, 1, nz), lh*ny, MPI_CPREC, up, 2,  &
                     comm, status, ierr)

  if (DEBUG) then
    write (*, *) coord, 'post: p_hat(3,2,0)=', p_hat(3,2,0)
    write (*, *) coord, 'post: p_hat(3,2,1)=', p_hat(3,2,1)
    write (*, *) coord, 'post: p_hat(3,2,nz-1)=', p_hat(3,2,nz-1)
    write (*, *) coord, 'post: p_hat(3,2,nz)=', p_hat(3,2,nz)
    if (coord == 0) then
      write (*, *) coord, 'compare: p_hat_t(:, 2, 3) =', p_hat_t(:, 2, 3)
    end if
  end if

  if (DEBUG) then
    write (*, *) 'p_hat_t:'
    write (*, *) '(0, 1, 1): ', p_hat_t(0, 1, 1)
    write (*, *) '(0, 2, 2): ', p_hat_t(0, 2, 2)
    write (*, *) '(1, 1, 1): ', p_hat_t(1, 1, 1)
    write (*, *) '(2, 1, 1): ', p_hat_t(2, 1, 1)
    write (*, *) '(1, 1, 2): ', p_hat_t(1, 1, 2)
  end if

  !--need to make sure x-Nyquist freqs are zero here
  p_hat(lh, :, :) = 0._rprec

  if (DEBUG) write (*, *) $str($context_doc), coord, ' after transpose p_hat'

$else
  !==========================================================================
  ! Loop over (Kx,Ky) to solve for Pressure amplitudes
  do jy=1,ny
     if (jy.eq.ny/2+1) then
        p_hat(:,jy,:)=0._rprec
        goto 5111
     end if
  ! since this is all in Fourier space, don't we need to skip ny/2+1?
     do jx=1,lh-1
  ! the lh part should be zero, right?
  ! this is the k=0 part
  ! Zero wavenumber problem is non-unique for matrix solution
  ! So just integrate up from wall with arbitrary P=const starting point.
  ! why is the jz index the innermost loop?
        if (jx*jy.eq.1) then          
  ! what is in p_hat(:,:,0) here?
  ! does the (:,:,1) entry need to be real as in original?
  ! p_hat(jx,jy,0) should be 0. here
          p_hat(jx,jy,1)=p_hat(jx,jy,0)-dz*bottomw(jx,jy)
          if (DEBUG) then
            write (*, *) '0-wave: p_hat(1, 1, 1) = ', p_hat(1, 1, 1)
            write (*, *) '0-wave: p_hat(1, 1, 0) = ', p_hat(1, 1, 0)
            write (*, *) '0-wave: jx, jy =', jx, jy
            write (*, *) '0-wave: dz =', dz
          end if
          do jz=2,nz
            p_hat(jx,jy,jz)=p_hat(jx,jy,jz-1)+H_z(jx,jy,jz)*dz
          end do   
       else
  ! See end of wavenumber loops for endif
  ! Set up matrix solution for pressure Fourier coeff's
  ! Near Wall Nodes
           a1(1)=0._rprec
           b1(1)=-1._rprec
           c1(1)=1._rprec
           RHS_col(1)=-dz*bottomw(jx,jy) ! complex 
  ! Interior Nodes
           do jz=2,nz
              RHS_col(jz)=eye*(kx(jx,jy)*H_x(jx,jy,jz-1)+&
                   ky(jx,jy)*H_y(jx,jy,jz-1))&
                   +(H_z(jx,jy,jz)-H_z(jx,jy,jz-1))/dz
  ! this part OK
              a1(jz)=1._rprec/(dz**2)
              b1(jz)=-(kx(jx,jy)**2+ky(jx,jy)**2+2._rprec/(dz**2))
              c1(jz)=1._rprec/(dz**2)
           end do

  ! top nodes 
           RHS_col(nz+1)=-topw(jx,jy)*dz
  ! note: check to make sure the topw's aren't automatically 0
  ! this would save some time
           a1(nz+1)= -1._rprec
           b1(nz+1)=1._rprec
           c1(nz+1)=0._rprec

           !if (DEBUG) then
           !  write (*, *) 'jx, jy = ', jx, jy
           !  write (*, *) 'a1 = ', a1
           !  write (*, *) 'b1 = ', b1
           !  write (*, *) 'c1 = ', c1
           !  write (*, *) 'RHS_col = ', RHS_col
           !end if

  ! Call Solver
  ! modified to handle complex
           call tridag (a1,b1,c1,RHS_col,p_col,nz+1)

  ! Put Pressure amplitudes in Matrix
           do jz=0,nz
  ! if called above directly with a row of p_hat, we'd probably encounter
  ! some serious cache thrashing??? try it
              p_hat(jx,jy,jz)=p_col(jz+1)
           end do   
  ! Ending IF from zero wavenumber solution
        end if
  ! End loop over kx,ky
     end do
   p_hat(lh,jy,:)=0._rprec
  5111 continue
  end do
$endif

if (DEBUG) call DEBUG_write (p_hat, 'press_stag.d.p_hat')

!=========================================================================== 
!...Now need to get p_hat(wave,level) to physical p(jx,jy,jz)   
!.....Loop over height levels     

if (DEBUG) write (*, *) $str($context_doc), ' reached line ', $line_num

call rfftwnd_f77_one_complex_to_real(back,p_hat(:,:,0),ignore_me)
do jz=1,nz
do jy=1,ny
do jx=1,lh
! complex
   dfdx(jx,jy,jz)=eye*kx(jx,jy)*p_hat(jx,jy,jz)
   dfdy(jx,jy,jz)=eye*ky(jx,jy)*p_hat(jx,jy,jz)
! note the oddballs of p_hat are already 0, so we should be OK here
end do
end do
call rfftwnd_f77_one_complex_to_real(back,dfdx(:,:,jz),ignore_me)
call rfftwnd_f77_one_complex_to_real(back,dfdy(:,:,jz),ignore_me)
call rfftwnd_f77_one_complex_to_real(back,p_hat(:,:,jz),ignore_me)
end do

$if ($MPI)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  contains
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_transpose_top_line (a, a_t)
  implicit none

  complex (rprec), intent (in) :: a(lh, ny, nz)
  complex (rprec), intent (inout) :: a_t(nz_t, ny, nx_t/2)

  complex (rprec) :: temp1(ny, nx_t/2)

  !-------------------------------------------------------------------

  if (coord == nproc-1) then

    do ip = 0, nproc-2 
      !--use temp as a buffer for chunk ip
      temp1 = transpose (a(ip*nx_t/2 + 1 : (ip+1)*nx_t/2, :, nz))
      !--send chunk to proc with coord ip
      call mpi_send (temp1(1, 1), ny*nx_t/2, MPI_CPREC,  &
                     rank_of_coord(ip), ip, comm, ierr)
    end do

    !--now transpose chunk nproc-1
    a_t(nz_t, :, 1:nx_t/2) = transpose ( a((nproc-1)*nx_t/2 + 1 :  &
                                            nproc*nx_t/2, :, nz) )

  else  !--recv chunk of a from proc with coord nproc-1

    call mpi_recv (temp1(1, 1), ny*nx_t/2, MPI_CPREC,                &
                   rank_of_coord(nproc-1), coord, comm, status, ierr)  
    a_t(nz_t, :, 1:nx_t/2) = temp1  !note not contiguous

  end if

  end subroutine mpi_transpose_top_line
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
$endif
end subroutine press_stag
