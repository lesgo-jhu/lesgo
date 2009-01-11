!************************************************************
subroutine read_inputs ()
!************************************************************
!
!  This subroutine reads the default input file, lesgo.inp. 
!  Several parameters are also set which use the information
!  from the input file.
!

use param2
implicit none

open(88,file='lesgo.inp',status='old')
!  Read namelists from input file
read(88,NML=geometry)
read(88,NML=tparam)
read(88,NML=coriolis)
read(88,NML=ic)
read(88,NML=ic)
read(88,NML=pparam)
read(88,NML=les)
read(88,NML=bc)
close(88)

!  Compute parameters dependant on input file information
nz = (nz -1)/nproc + 1
nz_tot = (nz - 1) * nproc + 1
nx2 = 3*nx/2
ny2 = 3*ny/2
L_z = 2._rprec*z_1/nproc
lh=nx/2+1
ld=2*lh
lh_big=nx2/2+1
ld_big=2*lh_big
! set the aspect ratio of the box, already nondimensional
dz=L_z/z_i/(nz-1)
dx=L_x/nx
dy=L_y/n

!--Coriolis stuff
! coriol=non-dim coriolis parameter,
! ug=horiz geostrophic vel, vg=transverse geostrophic vel
! u_star=0.45 if coriolis_forcing=.FALSE. and =ug if coriolis_forcing=.TRUE
coriol=9.125E-05*z_i/u_star
ug=u_star/u_star
vg=0._rprec/u_star

mean_p_force = 1._rprec * z_i/L_z/nproc !--usually just z_i/L_

return
end subroutine read_inputs
