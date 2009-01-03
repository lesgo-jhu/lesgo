      program interp
! interpolate initial conditions--in physical space
! this is not intended to be a subroutine...
! only for 2d-stuff!!!!!
! could automatically determine the size of the smaller ones from size of C  the file.
      implicit none
      integer, parameter :: nx_s=64, ny_s=64, nz_s=64
      integer, parameter :: nx_b=120, ny_b=120 , nz_b=120
      character(len=9) :: path='./'
      double precision, dimension(nx_s+2,ny_s,nz_s) :: u_s, v_s, w_s,
     &  RHSx_s, RHSy_s , RHSz_s, cs_s, FLM_s, FMM_s
      double precision, dimension(nx_b+2,ny_b,nz_b) :: u_b, v_b, w_b,
     &  RHSx_b, RHSy_b ,RHSz_b, cs_b, FLM_b, FMM_b
            
      write(6,*) 'Going to interpolate "'//path//'vel.out".'
      write(6,*) 'start size is ', nx_s,' X ',nx_s,' X ', nz_s
      write(6,*) 'new size is ', nx_s,' X ',nx_s,' X ', nz_s
  
      open(1,file='vel.out',form='unformatted')
      read(1) u_s,v_s,w_s,RHSx_s,RHSy_s,RHSz_s,cs_s,FLM_s,FMM_s
      close(1)

	call interpolate(u_s,nx_s,ny_s,nz_s,u_b,nx_b,ny_b,nz_b)
	call interpolate(v_s,nx_s,ny_s,nz_s,v_b,nx_b,ny_b,nz_b)
	call interpolate(w_s,nx_s,ny_s,nz_s,w_b,nx_b,ny_b,nz_b)
	call interpolate(RHSx_s,nx_s,ny_s,nz_s,RHSx_b,nx_b,ny_b,nz_b)
	call interpolate(RHSy_s,nx_s,ny_s,nz_s,RHSy_b,nx_b,ny_b,nz_b)
	call interpolate(RHSz_s,nx_s,ny_s,nz_s,RHSz_b,nx_b,ny_b,nz_b)
	call interpolate(cs_s,nx_s,ny_s,nz_s,cs_b,nx_b,ny_b,nz_b)
	call interpolate(FLM_s,nx_s,ny_s,nz_s,FLM_b,nx_b,ny_b,nz_b)
	call interpolate(FMM_s,nx_s,ny_s,nz_s,FMM_b,nx_b,ny_b,nz_b)


      open(1,file=path//'vel-interp.out',form='unformatted')
      write(1) u_b,v_b,w_b,RHSx_b,RHSy_b,RHSz_b,cs_b,FLM_b,FMM_b
      close(1)

      end program interp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
	subroutine interpolate(u_s,nx_s,ny_s,nz_s,u_b,nx_b,ny_b,nz_b)

      implicit none
      integer :: nx_s, ny_s, nz_s, nx_b, ny_b, nz_b
      double precision, dimension(nx_b+2,ny_b,nz_b) :: u_b
      double precision, dimension(nx_s+2,ny_s,nz_s) :: u_s
      double precision :: rx, ry, rz, P1, P2, P3, P4, P5, P6
	double precision :: dx, dy, dz
      integer :: jx, jy, jz, i, j, k, i_wrap, j_wrap				!!????????

      rx = real(nx_s)/real(nx_b)
	ry = real(ny_s)/real(ny_b)
      rz = real(nz_s-1)/real(nz_b-1)

      do jz=1,nz_b-1  ! we already know the last point!
        k = floor((jz-1)*rz) + 1
        dz = (jz-1)*rz-(k-1)
		do jy=1,ny_b
          j = floor((jy-1)*ry) + 1
          j_wrap = modulo(j+1-1,ny_s) + 1
		dy = (jy-1)*ry-(j-1)
			do jx=1,nx_b
			i = floor((jx-1)*rx) + 1
			i_wrap = modulo(i+1-1,nx_s) + 1
			dx = (jx-1)*rx-(i-1)

		P1 = (1.-dx)*u_s(i,j,k) + dx*u_s(i_wrap,j,k)
		P2 = (1.-dx)*u_s(i,j,k+1) + dx*u_s(i_wrap,j,k+1)
		P3 = (1.-dx)*u_s(i,j_wrap,k) + dx*u_s(i_wrap,j_wrap,k)
		P4 = (1.-dx)*u_s(i,j_wrap,k+1) + dx*u_s(i_wrap,j_wrap,k+1)
		P5 = (1.-dy)*P1 + dy*P3 
		P6 = (1.-dy)*P2 + dy*P4

			u_b(jx,jy,jz) = (1.-dz)*P5 + dz*P6

			end do
		end do
      end do

! still have to do nz_b level
	do jy=1,ny_b
        j = floor((jy-1)*ry) + 1
        j_wrap = modulo(j+1-1,ny_s) + 1
	dy = (jy-1)*ry-(j-1)
		do jx=1,nx_b
		i = 1 + floor((jx-1)*rx)
		i_wrap = modulo(i+1-1,nx_s) + 1
		dx = (jx-1)*rx-(i-1)

		P1 = (1.-dx)*u_s(i,j,nz_s) + dx*u_s(i_wrap,j,nz_s)
		P2 = (1.-dx)*u_s(i,j_wrap,nz_s) + dx*u_s(i_wrap,j_wrap,nz_s)

		u_b(jx,jy,nz_b) = (1.-dy)*P1 + dy*P2

		end do

      end do
      u_b(nx_b+1:nx_b+2,:,:) = 0.

      end subroutine interpolate



 
