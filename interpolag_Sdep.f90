! this is the w-node version
subroutine interpolag_Sdep()
! This subroutine computes the values of F_LM and F_MM 
! at positions (x-u*dt) for use in the subroutine lagrangian
use types,only:rprec
!use param,only:ld,nx,ny,nz,dx,dy,dz,dt,cs_count,c_count
use param
use sgsmodule
implicit none
integer :: jx, jy, jz 
integer :: jjx, jjy, jjz
integer :: jz_min
$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif

!--saves here: force heap storage
!real(kind=rprec), save, dimension(ld,ny,$lbz:nz) :: xp, yp, zp
                   !--use u_lag, v_lag, w_lag instead of xp, yp ,zp
!real(kind=rprec), save, dimension(ld,ny,$lbz:nz) :: u_temp, v_temp
                   !--removed u_temp, v_temp to save mem.

real(kind=rprec) :: frac_x,frac_y,frac_z
real(kind=rprec) :: comp_x,comp_y,comp_z
real(kind=rprec), save, dimension(nx+2,ny+2,nz+2) :: FF_LM,FF_MM
real(kind=rprec), save, dimension(nx+2,ny+2,nz+2) :: FF_QN,FF_NN

integer :: addx, addy, addz
integer :: jxaddx, jyaddy, jzaddz
real(kind=rprec), save, dimension(nx+2,ny+2,nz+2) :: Beta_t

!---------------------------------------------------------------------
$if ($VERBOSE)
write (*, *) 'started interpolag_Sdep'
$endif

!     creates dummy arrays FF_LM and FF_MM to use in the subroutine
 
do jz=1,nz
do jy=1,ny
do jx=1,nx
   FF_LM(jx+1,jy+1,jz+1) = F_LM(jx,jy,jz)
   FF_MM(jx+1,jy+1,jz+1) = F_MM(jx,jy,jz)
   FF_QN(jx+1,jy+1,jz+1) = F_QN(jx,jy,jz)
   FF_NN(jx+1,jy+1,jz+1) = F_NN(jx,jy,jz)
end do
end do
end do

!This is a bit like witch craft but the following lines do take care
!of all the edges including the corners
FF_LM(1,:,:) = FF_LM(nx+1,:,:)
FF_LM(nx+2,:,:) = FF_LM(2,:,:)

FF_LM(:,1,:) = FF_LM(:,ny+1,:) 
FF_LM(:,ny+2,:) = FF_LM(:,2,:) 

$if ($MPI)
  !--send F_LM @ jz=nz-1 to F_LM @ jz=0'
  !  i.e. FF_LM @ jz=nz to FF_LM @ jz=1'
  call mpi_sendrecv (FF_LM(1, 1, nz), (nx+2)*(ny+2), MPI_RPREC, up, 1,   &
                     FF_LM(1, 1, 1), (nx+2)*(ny+2), MPI_RPREC, down, 1,  &
                     comm, status, ierr)
  !--F_LM @ jz=nz and F_LM @ jz=1' should already be in sync (test?)
  !  i.e. FF_LM @ jz=nz+1 and FF_LM @ jz=2'
$endif

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
  FF_LM(:,:,1) = FF_LM(:,:,2)
end if

$if ($MPI)
  !--send F_LM @ jz=2 to F_LM @ jz=nz+1'
  !  i.e. FF_LM @ jz=3 to FF_LM @ jz=nz+2'
  call mpi_sendrecv (FF_LM(1, 1, 3), (nx+2)*(ny+2), MPI_RPREC, down, 2,   &
                     FF_LM(1, 1, nz+2), (nx+2)*(ny+2), MPI_RPREC, up, 2,  &
                     comm, status, ierr)
$endif

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
  FF_LM(:,:,nz+2) = FF_LM(:,:,nz+1)
end if

FF_MM(1,:,:) = FF_MM(nx+1,:,:)
FF_MM(nx+2,:,:) = FF_MM(2,:,:)

FF_MM(:,1,:) = FF_MM(:,ny+1,:) 
FF_MM(:,ny+2,:) = FF_MM(:,2,:) 

$if ($MPI)
  !--send F_MM @ jz=nz-1 to F_MM @ jz=0'
  !  i.e. FF_MM @ jz=nz to FF_MM @ jz=1'
  call mpi_sendrecv (FF_MM(1, 1, nz), (nx+2)*(ny+2), MPI_RPREC, up, 3,   &
                     FF_MM(1, 1, 1), (nx+2)*(ny+2), MPI_RPREC, down, 3,  &
                     comm, status, ierr)
  !--F_MM @ jz=nz and F_MM @ jz=1' should already be in sync (test?)
  !  i.e. FF_MM @ jz=nz+1 and FF_MM @ jz=2'
$endif

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
  FF_MM(:,:,1) = FF_MM(:,:,2)
end if

$if ($MPI)
  !--send F_MM @ jz=2 to F_MM @ jz=nz+1'
  !  i.e. FF_MM @ jz=3 to FF_MM @ jz=nz+2'
  call mpi_sendrecv (FF_MM(1, 1, 3), (nx+2)*(ny+2), MPI_RPREC, down, 4,   &
                     FF_MM(1, 1, nz+2), (nx+2)*(ny+2), MPI_RPREC, up, 4,  &
                     comm, status, ierr)
$endif

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
  FF_MM(:,:,nz+2) = FF_MM(:,:,nz+1) 
end if

FF_QN(1,:,:) = FF_QN(nx+1,:,:)
FF_QN(nx+2,:,:) = FF_QN(2,:,:)

FF_QN(:,1,:) = FF_QN(:,ny+1,:) 
FF_QN(:,ny+2,:) = FF_QN(:,2,:) 

$if ($MPI)
  !--send F_QN @ jz=nz-1 to F_QN @ jz=0'
  !  i.e. FF_QN @ jz=nz to FF_QN @ jz=1'
  call mpi_sendrecv (FF_QN(1, 1, nz), (nx+2)*(ny+2), MPI_RPREC, up, 5,   &
                     FF_QN(1, 1, 1), (nx+2)*(ny+2), MPI_RPREC, down, 5,  &
                     comm, status, ierr)
  !--F_QN @ jz=nz and F_QN @ jz=1' should already be in sync (test?)
  !  i.e. FF_QN @ jz=nz+1 and FF_QN @ jz=2'
$endif

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
  FF_QN(:,:,1) = FF_QN(:,:,2)
end if

$if ($MPI)
  !--send F_QN @ jz=2 to F_QN @ jz=nz+1'
  !  i.e. FF_QN @ jz=3 to FF_QN @ jz=nz+2'
  call mpi_sendrecv (FF_QN(1, 1, 3), (nx+2)*(ny+2), MPI_RPREC, down, 6,   &
                     FF_QN(1, 1, nz+2), (nx+2)*(ny+2), MPI_RPREC, up, 6,  &
                     comm, status, ierr)
$endif

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
  FF_QN(:,:,nz+2) = FF_QN(:,:,nz+1) 
end if

FF_NN(1,:,:) = FF_NN(nx+1,:,:)
FF_NN(nx+2,:,:) = FF_NN(2,:,:)

FF_NN(:,1,:) = FF_NN(:,ny+1,:) 
FF_NN(:,ny+2,:) = FF_NN(:,2,:) 

$if ($MPI)
  !--send F_NN @ jz=nz-1 to F_NN @ jz=0'
  !  i.e. FF_NN @ jz=nz to FF_NN @ jz=1'
  call mpi_sendrecv (FF_NN(1, 1, nz), (nx+2)*(ny+2), MPI_RPREC, up, 7,   &
                     FF_NN(1, 1, 1), (nx+2)*(ny+2), MPI_RPREC, down, 7,  &
                     comm, status, ierr)
  !--F_NN @ jz=nz and F_NN @ jz=1' should already be in sync (test?)
  !  i.e. FF_NN @ jz=nz+1 and FF_NN @ jz=2'
$endif

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
  FF_NN(:,:,1) = FF_NN(:,:,2)
end if

$if ($MPI)
  !--send F_NN @ jz=2 to F_NN @ jz=nz+1'
  !  i.e. FF_NN @ jz=3 to FF_NN @ jz=nz+2'
  call mpi_sendrecv (FF_NN(1, 1, 3), (nx+2)*(ny+2), MPI_RPREC, down, 8,   &
                     FF_NN(1, 1, nz+2), (nx+2)*(ny+2), MPI_RPREC, up, 8,  &
                     comm, status, ierr)
$endif

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
  FF_NN(:,:,nz+2) = FF_NN(:,:,nz+1) 
end if
! end of witch craft

!--this is re-written below: save mem by removing u_temp, v_temp
!!	puts u_lag and and v_lag on cs nodes
!u_temp = u_lag/real(cs_count,kind=rprec)
!v_temp = v_lag/real(cs_count,kind=rprec)
!
!$if ($MPI)
!  !--not sure if u_lag, u_temp at 0 are needed yet
!  u_lag(:, :, $lbz) = BOGUS
!  v_lag(:, :, $lbz) = BOGUS
!$endif
!
!if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
!
!  u_lag (:,:,1) = u_temp(:,:,1)
!  v_lag (:,:,1) = v_temp(:,:,1)
!
!  jz_min = 2
!
!else
!
!  jz_min = 1
!end if
!
!do jz = jz_min, nz
!   u_lag (:,:,jz) = 0.5_rprec*(u_temp(:,:,jz) + u_temp(:,:,jz-1))
!   v_lag (:,:,jz) = 0.5_rprec*(v_temp(:,:,jz) + v_temp(:,:,jz-1))
!end do

!	puts u_lag and and v_lag on cs nodes
u_lag = u_lag/real(cs_count,kind=rprec)
v_lag = v_lag/real(cs_count,kind=rprec)

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
  jz_min = 2
else
  jz_min = 1
end if

!--must do backwards due to data depenency on plane below current one
do jz = nz, jz_min, -1
  u_lag(:, :, jz) = 0.5_rprec * (u_lag(:, :, jz) + u_lag(:, :, jz-1))
  v_lag(:, :, jz) = 0.5_rprec * (v_lag(:, :, jz) + v_lag(:, :, jz-1))
end do

w_lag = w_lag/real(cs_count,kind=rprec)

$if ($MPI)
  !--not sure if 0-level used
  w_lag(:, :, $lbz) = BOGUS
$endif

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
  w_lag (:,:,1) = 0.25_rprec*w_lag (:,:,2)
end if

!     computes the 3-D inverse displacement arrays that describe
!     the location where the point was at the previous step
!xp = -u_lag*dt*real(cs_count,kind=rprec)/dx    !! is u already normalized
!yp = -v_lag*dt*real(cs_count,kind=rprec)/dy   
!zp = -w_lag*dt*real(cs_count,kind=rprec)/dz
!--use u_lag to store xp, etc.
u_lag = -u_lag*dt*real(cs_count,kind=rprec)/dx    !! is u already normalized
v_lag = -v_lag*dt*real(cs_count,kind=rprec)/dy   
w_lag = -w_lag*dt*real(cs_count,kind=rprec)/dz

$if ($MPI)
  !--perhaps can remove 0-level altogether
  !xp(:, :, $lbz) = BOGUS
  !yp(:, :, $lbz) = BOGUS
  !zp(:, :, $lbz) = BOGUS
  u_lag(:, :, $lbz) = BOGUS
  v_lag(:, :, $lbz) = BOGUS
  w_lag(:, :, $lbz) = BOGUS
$endif

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then

  !	because first plane is on u,v,p nodes
  !	this corrects for the fact that the first cell
  !	in the z direction has height dz/2
  !	it doubles the zp fraction if this fraction relates to the cell
  !     beneath it
  do jy=1,ny
  do jx=1,nx
     !zp(jx,jy,2) = zp(jx,jy,2) + min (zp(jx,jy,2), 0._rprec)
     !zp(jx,jy,1) = zp(jx,jy,1) + max (zp(jx,jy,1), 0._rprec)
     w_lag(jx,jy,2) = w_lag(jx,jy,2) + min (w_lag(jx,jy,2), 0._rprec)
     w_lag(jx,jy,1) = w_lag(jx,jy,1) + max (w_lag(jx,jy,1), 0._rprec)
  end do
  end do

  w_lag(:,:,2) = min (1._rprec, w_lag(:,:,2))
  w_lag(:,:,1) = min (1._rprec, w_lag(:,:,2))
  w_lag(:,:,2) = max (-1._rprec, w_lag(:,:,2))
  w_lag(:,:,1) = max (-1._rprec, w_lag(:,:,2))

end if

!if(mod(jt,c_count).eq.0) print*,'Lagrangian CFL condition= ',  &
!                         maxval ( abs (xp(1:nx, :, 1:nz)) )
if(mod(jt,c_count).eq.0) print*,'Lagrangian CFL condition= ',  &
                         maxval ( abs (u_lag(1:nx, :, 1:nz)) )


do jz=1,nz
 jjz = jz+1
 do jy=1,ny
  jjy = jy+1
   do jx=1,nx
      jjx = jx+1
!     the are the values to add to the indices jx, jy, and jz
!     the are +1 or -1 depending on what cube should be used for interpolation
   !addx = int(sign(1._rprec,xp(jx,jy,jz)))
   !addy = int(sign(1._rprec,yp(jx,jy,jz)))
   !addz = int(sign(1._rprec,zp(jx,jy,jz)))
   addx = int(sign(1._rprec,u_lag(jx,jy,jz)))
   addy = int(sign(1._rprec,v_lag(jx,jy,jz)))
   addz = int(sign(1._rprec,w_lag(jx,jy,jz)))
   jxaddx = jjx + addx 
   jyaddy = jjy + addy
   jzaddz = jjz + addz
!     computes the relative weights given to F_** in the cube depending on point location
   !comp_x = abs(xp(jx,jy,jz))
   !comp_y = abs(yp(jx,jy,jz))
   !comp_z = abs(zp(jx,jy,jz)) 
   comp_x = abs(u_lag(jx,jy,jz))
   comp_y = abs(v_lag(jx,jy,jz))
   comp_z = abs(w_lag(jx,jy,jz)) 
   frac_x = 1._rprec - comp_x
   frac_y = 1._rprec - comp_y
   frac_z = 1._rprec - comp_z

!     computes interpolated F_LM

   F_LM(jx,jy,jz)=frac_x*frac_y*&
        (FF_LM(jjx,jjy,jjz)*frac_z+FF_LM(jjx,jjy,jzaddz)*comp_z)&
       + frac_x*comp_y*&
      (FF_LM(jjx,jyaddy,jjz)*frac_z+FF_LM(jjx,jyaddy,jzaddz)*comp_z)&
       + comp_x*frac_y*&
      (FF_LM(jxaddx,jjy,jjz)*frac_z+FF_LM(jxaddx,jjy,jzaddz)*comp_z)&
       + comp_x*comp_y*&
      (FF_LM(jxaddx,jyaddy,jjz)*frac_z+FF_LM(jxaddx,jyaddy,jzaddz)*comp_z)
 
!     computes interpolated F_MM
   F_MM(jx,jy,jz)=frac_x*frac_y*&
        (FF_MM(jjx,jjy,jjz)*frac_z+FF_MM(jjx,jjy,jzaddz)*comp_z)&
        + frac_x*comp_y*&
        (FF_MM(jjx,jyaddy,jjz)*frac_z+FF_MM(jjx,jyaddy,jzaddz)*comp_z)&
        + comp_x*frac_y*&
        (FF_MM(jxaddx,jjy,jjz)*frac_z+FF_MM(jxaddx,jjy,jzaddz)*comp_z)&
        + comp_x*comp_y*&
        (FF_MM(jxaddx,jyaddy,jjz)*frac_z+FF_MM(jxaddx,jyaddy,jzaddz)*comp_z)

!     computes interpolated F_QN
   F_QN(jx,jy,jz)=frac_x*frac_y*&
        (FF_QN(jjx,jjy,jjz)*frac_z+FF_QN(jjx,jjy,jzaddz)*comp_z)&
        + frac_x*comp_y*&
        (FF_QN(jjx,jyaddy,jjz)*frac_z+FF_QN(jjx,jyaddy,jzaddz)*comp_z)&
        + comp_x*frac_y*&
        (FF_QN(jxaddx,jjy,jjz)*frac_z+FF_QN(jxaddx,jjy,jzaddz)*comp_z)&
        + comp_x*comp_y*&
        (FF_QN(jxaddx,jyaddy,jjz)*frac_z+FF_QN(jxaddx,jyaddy,jzaddz)*comp_z)
 
!     computes interpolated F_NN
   F_NN(jx,jy,jz)=frac_x*frac_y*&
       (FF_NN(jjx,jjy,jjz)*frac_z+FF_NN(jjx,jjy,jzaddz)*comp_z)&
      + frac_x*comp_y*&
       (FF_NN(jjx,jyaddy,jjz)*frac_z+FF_NN(jjx,jyaddy,jzaddz)*comp_z)&
      + comp_x*frac_y*&
       (FF_NN(jxaddx,jjy,jjz)*frac_z+FF_NN(jxaddx,jjy,jzaddz)*comp_z)&
      + comp_x*comp_y*&
       (FF_NN(jxaddx,jyaddy,jjz)*frac_z+FF_NN(jxaddx,jyaddy,jzaddz)*comp_z)
end do
end do
end do

u_lag = 0._rprec
v_lag = 0._rprec
w_lag = 0._rprec

$if ($VERBOSE)
write (*, *) 'finished interpolag_Sdep'
$endif

end subroutine interpolag_Sdep
