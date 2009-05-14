program cylinder_skew
use error
use root
implicit none

type cs
     double precision :: x,y,z, x0,y0,z0, phi
end type cs

type(cs), allocatable, dimension(:,:,:) :: gcs_t
type(cs) :: lcs_t, scs_t  ! global, local, skewed coordinate system respectivly
type(cs) :: tcs_t

integer, parameter :: Nx=128, Ny=128, Nz=128

double precision, parameter :: skew_angle=1.*3.14/180. !  In radians
double precision, parameter :: crad = 0.1 !  Cylinder radius
double precision, parameter :: clength = 0.5 !  Cylinder length
double precision, parameter :: thresh = 1.e-9
double precision :: echeck, aa,bb,cc,dd,ee, dist
double complex, dimension(4) :: rx,ry
integer :: lcase,i,j,k,icount,ir,nf

double precision, parameter :: xmin=0., xmax=1., dx=(xmax-xmin)/(Nx-1)
double precision, parameter :: ymin=0., ymax=1., dy=(ymax-ymin)/(Ny-1)
double precision, parameter :: zmin=0., zmax=1., dz=(zmax-zmin)/(Nz-1)
double precision, parameter :: a=crad/cos(skew_angle), b=crad

!  Allocate x,y,z for all coordinate systems
allocate(gcs_t(nx,ny,nz))

!  Create grid in the global coordinate system
do k=1,Nz
  do j=1,ny
    do i=1,nx
      gcs_t(i,j,k)%x=(i-1)*dx
      gcs_t(i,j,k)%y=(j-1)*dy
      gcs_t(i,j,k)%z=(k-1)*dz
    enddo
  enddo
enddo

!  Specify lcs origin in terms of gcs
lcs_t%x0=0.333
lcs_t%y0=0.5
lcs_t%z0=0.25

!  Loop over all global coordinates
do k=1,Nz
  do j=1,ny
    do i=1,nx
!  Convert global coordinate to local coordinate system (no rotation)
      lcs_t%x=gcs_t(i,j,k)%x - lcs_t%x0
      lcs_t%y=gcs_t(i,j,k)%y - lcs_t%y0
      lcs_t%z=gcs_t(i,j,k)%z - lcs_t%z0

!  Transform local coordinates to skewed (local) coordinates
      scs_t%x = lcs_t%x - lcs_t%z*dtan(skew_angle)
      scs_t%y = lcs_t%y
      scs_t%z = lcs_t%z/dcos(skew_angle)

!  Set up coefficients for minimum dist calculations
      aa=b*b*(b*b-a*a)
      bb=2.*a**2*b**2*(b*b-a*a)*scs_t%x
      cc=a**4*b**2*scs_t%x**2+a**2*b**4*scs_t%y**2 - a**2*b**2*(b**2-a**2)
      dd=-2.*a**4*b**2*(b**2-a**2)*scs_t%x
      ee=-a**6*b**2*scs_t%x**2

!  Check where the point lies wrt to the scs.
      echeck = (scs_t%x/a)**2 + (scs_t%y/b)**2
      if(scs_t%z > 0 .and. scs_t%z < clength) then
        if(echeck > 1) then
          lcase=1
        else
          lcase=4
        endif
      elseif(scs_t%z < 0) then
        if(echeck > 1) then
          lcase=3
        else
          lcase=6
        endif
      else
         if(echeck > 1) then
          lcase=2
        else
          lcase=5
        endif       
      endif
      if(skew_angle <= thresh) then !  Non skewed cylinder
!  Check which case is applied
        if(lcase == 1) then
          gcs_t(i,j,k)%phi = dsqrt(lcs_t%x**2 + lcs_t%y**2) - crad
        elseif(lcase == 2) then
          gcs_t(i,j,k)%phi = dsqrt((dsqrt(lcs_t%x**2 + lcs_t%y**2) - crad)**2 + (lcs_t%z - clength)**2)
        elseif(lcase == 3) then    
          gcs_t(i,j,k)%phi = dsqrt((dsqrt(lcs_t%x**2 + lcs_t%y**2) - crad)**2 + lcs_t%z**2)
        elseif(lcase == 4) then
          gcs_t(i,j,k)%phi = -dabs((dsqrt(lcs_t%x**2 + lcs_t%y**2) - crad))
        elseif(lcase == 5) then
          gcs_t(i,j,k)%phi = dabs(lcs_t%z - clength)
        else
          gcs_t(i,j,k)%phi = dabs(lcs_t%z)
        endif
      else !  Skewed cylinder
        gcs_t(i,j,k)%phi=huge(1.)
        if(lcase == 1 .or. lcase == 4) then
          call RootPol(bb/aa,cc/aa,dd/aa,ee/aa,rx(1),rx(2),rx(3),rx(4))
          icount=0
          do ir=1,4
            ry(ir) = scs_t%y/(1 - (rx(ir) - scs_t%x)/rx(ir)*a**2/b**2)
            if(dimag(rx(ir)) <= thresh) then
              icount=icount+1
              tcs_t%x = real(rx(ir)) + scs_t%z*sin(skew_angle)
              tcs_t%y = real(ry(ir))
              tcs_t%z = scs_t%z*dsin(skew_angle)
!  Compute distance
              dist = sqrt((lcs_t%x - tcs_t%x)**2 + (lcs_t%y - tcs_t%y)**2 + (lcs_t%z - tcs_t%z)**2)
              if(dist < gcs_t(i,j,k)%phi) then
                gcs_t(i,j,k)%phi = dist !  Keep it      
              endif
            endif
          enddo
          if(lcase == 4) gcs_t(i,j,k)%phi = -gcs_t(i,j,k)%phi
          if(icount == 0) then
            write(*,*) 'No exclusively real roots found!'
            stop
          endif
        elseif(lcase == 2 .or. lcase == 3) then
          call RootPol(bb/aa,cc/aa,dd/aa,ee/aa,rx(1),rx(2),rx(3),rx(4))
          icount=0
          do ir=1,4
            ry(ir) = scs_t%y/(1 - (rx(ir) - scs_t%x)/rx(ir)*a**2/b**2)
            if(dimag(rx(ir)) <= thresh) then
              icount=icount+1
              if(lcase == 2) then
                tcs_t%x = real(rx(ir)) + clength*dsin(skew_angle)
                tcs_t%y = real(ry(ir))
                tcs_t%z = clength*dsin(skew_angle)
              elseif(lcase == 3) then
                tcs_t%x = real(rx(ir))
                tcs_t%y = real(ry(ir))
                tcs_t%z = 0.
              endif
!  Compute distance
              dist = dsqrt((lcs_t%x - tcs_t%x)**2 + (lcs_t%y - tcs_t%y)**2 + (lcs_t%z - tcs_t%z)**2)
              if(dist < gcs_t(i,j,k)%phi) then
                gcs_t(i,j,k)%phi = dist !  Keep it      
              endif
            endif
          enddo   
          if(icount == 0) then
            write(*,*) 'No exclusively real roots found!'
            stop
          endif
        elseif(lcase == 5) then
          gcs_t(i,j,k)%phi = dabs(scs_t%z - clength)
          write(*,*) 'scs_t%z - clength = ', scs_t%z - clength
        elseif(lcase == 6) then
          gcs_t(i,j,k)%phi = dabs(scs_t%z)
        endif
      endif  

      if(gcs_t(i,j,k)%phi > 0. .and. gcs_t(i,j,k)%phi <= thresh) gcs_t(i,j,k)%phi = 0.  

    enddo
  enddo
enddo
nf=1
!  Create tecplot formatted velocity field file  
    open (unit = 2,file = 'cylinder_skew.dat', status='unknown',form='formatted', &
      action='write',position='rewind')
    write(2,*) 'variables = "x", "y", "z", "phi"'; 
    write(2,"(1a,i9,1a,i3,1a,i3,1a,i3,1a,i3)") 'ZONE T="', &
      nf,'", DATAPACKING=POINT, i=', Nx,', j=',Ny, ', k=', Nz
    write(2,"(1a)") ''//adjustl('DT=(DOUBLE DOUBLE DOUBLE DOUBLE)')//''

    do k=1,nz
      do j=1,ny
        do i=1,nx
          write(2,*) gcs_t(i,j,k)%x, gcs_t(i,j,k)%y, gcs_t(i,j,k)%z, gcs_t(i,j,k)%phi
        enddo
      enddo
    enddo
    close(2)
!  Create plt file the hard way
!    write(fpreplt,*) 'preplot ',ftec,' && rm -v ', ftec
!    call system(fpreplt);


stop
end program cylinder_skew
