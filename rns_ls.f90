!**********************************************************************
module rns_ls
!**********************************************************************
use types, only : rprec
use param
use rns_base_ls
implicit none

save
private

public :: rns_init_ls

contains

!**********************************************************************
subroutine rns_init_ls()
!**********************************************************************

implicit none

character(64) :: fname, temp

!!!  Open file which to write global data
!!fname = path // 'cylinder_skew_gen_ls.out'
!!$if ($MPI)
!!  write (temp, '(".c",i0)') coord
!!  fname = trim (fname) // temp
!!$endif

!!!  Read in cylinder_skew_gen.dat file
!!open (unit = 2,file = fname, status='old',form='formatted', &
!!  action='read',position='rewind')

!!allocate(igen(ngen))
!!allocate(kbottom_inside(ngen))
!!allocate(kbottom(ngen))
!!allocate(dz_bottom(ngen))
!!allocate(ktop_inside(ngen))
!!allocate(ktop(ngen))
!!allocate(dz_top(ngen))
!!allocate(itype(nx+2,ny,nz))

!!do ng=1,ngen
!!  read(2,*) igen(ng), kbottom_inside(ng), kbottom(ng), &
!!    dz_bottom(ng), ktop_inside(ng), ktop(ng), &
!!    dz_top(ng)
!!enddo
!!close(2)

!!!  Check 1st generation only for need ground association
!!if(igen(1) /= -1) then
!!  !  Open file which to write global data
!!  fname = path // 'cylinder_skew_point_ls.out'
!!  $if ($MPI)
!!    write (temp, '(".c",i0)') coord
!!    fname = trim (fname) // temp
!!  $endif

!!  !  Read in cylinder_skew_gen.dat file
!!  open (unit = 2,file = fname, status='old',form='formatted', &
!!    action='read',position='rewind')
!!  do k=1,nz
!!    do j = 1,ny
!!      do i = 1,nx+2
!!        read(2,*) itype(i,j,k)
!!      enddo
!!    enddo
!!  enddo
!!   
!!  close(2)
!!endif

return
end subroutine rns_init_ls

end module rns_ls
