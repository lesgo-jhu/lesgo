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

module trees_io_ls
use trees_base_ls
implicit none

save
private

public :: draw_tree_array, write_ta_data, read_ta_data

character (*), parameter :: mod_name = 'trees_io_ls'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--assumes file with unit no. lun is opened and positioned correctly
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
recursive subroutine draw_branches (branch, lun)
implicit none

type (branch_type), intent (in) :: branch

integer, intent (in) :: lun

character (*), parameter :: sub_name = mod_name // '.draw_branches'

integer, parameter :: n_face = 4
integer, parameter :: n_height = 2  !--must be > 1

integer :: i, j
integer :: n_tot, n_cap, n_base

logical :: add_base_this  !--local, changeable version of add_base

real (rp) :: pi
real (rp) :: h, dh, theta, dtheta, theta0
real (rp) :: r, r0, r_top, l
real (rp) :: x_local, y_local, z_local
real (rp) :: x_abs, y_abs, z_abs

!----------------------------------------------------------------------
pi = acos (-1._rp)

if (add_cap) then

  select case (cap_shape)
    case ('hemispherical')
      n_cap = 2
    case ('rectangular')
      n_cap = 2
    case default
      call error (sub_name, 'invalid cap_shape')
  end select

else

  n_cap = 1  !--draw lid

end if

if (add_base .and. (branch % gen > 0)) then
  add_base_this = add_base
else  !--do not draw base for trunks
  add_base_this = .false.
end if

if (add_base_this) then

  select case (base_shape)
    case ('hemispherical')
      n_base = 2
    case ('rectangular')
      n_base = 2
    case default
      call error (sub_name, 'invalid cap_shape')
  end select

else

  n_base = 1  !--draw bottom

end if

n_tot = n_height + n_cap + n_base

if (use_tecplot) call write_tecplot_zone_hdr ()

!--define n_face X n_height mesh to represent the branch
dtheta = 2._rp * pi / n_face

l = branch % l

dh = l / (n_height - 1)

select case (branch_cross_section)
  case ('circular')

    r0 = (branch % d) / 2._rp
    theta0 = branch % twist

  case ('square')

    if (n_face /= 4) then
      call error (sub_name, 'use n_face = 4 with square X-section')
    end if

    r0 = (branch % d) / sqrt (2._rp)
    theta0 = pi / 4._rp + (branch % twist)

  case ('square+plate')

    call mesg (sub_name, 'cannot draw square+plate cross-section')
    goto 001  !--do nothing, exit this nicely
    
  case default

    call error (sub_name, 'invalid branch_cross_section')

end select

r_top = r0 * (1._rp - (branch % taper))

do j = 1, n_tot

  theta = theta0 

  call set_h ()

  call set_radius ()  !--set the radius used for drawing

  call draw_cross_section ()

end do

if (associated (branch % sub_branch)) then  !--this is safer, perhaps

  do i = 1, branch % n_sub_branch

    call draw_branches (branch % sub_branch(i), lun)  !--recursion

  end do

end if

001 continue

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine draw_cross_section ()
  implicit none

  do i = 1, n_face + 1  !--(+1) needed to complete circle

    x_local = r * cos (theta)
    y_local = r * sin (theta) 
    z_local = h

    x_abs = (branch % x0(1)) +               &
            x_local * (branch % x_hat(1)) +  &
            y_local * (branch % y_hat(1)) +  &
            z_local * (branch % z_hat(1))

    y_abs = (branch % x0(2)) +               &
            x_local * (branch % x_hat(2)) +  &
            y_local * (branch % y_hat(2)) +  &
            z_local * (branch % z_hat(2))

    z_abs = (branch % x0(3)) +               &
            x_local * (branch % x_hat(3)) +  &
            y_local * (branch % y_hat(3)) +  &
            z_local * (branch % z_hat(3))

    write (lun, '(3(f12.5,1x))') x_abs, y_abs, z_abs

    theta = theta + dtheta

  end do

  end subroutine draw_cross_section

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine set_h ()
  implicit none

  if (j <= n_base) then

    if (add_base_this) then

      select case (base_shape)
        case ('hemispherical')  !--n_base = 2

          if (j == 1) then
            h = -(branch % d) / 2._rp
          else  !--j = 2
            h = -((branch % d) / 2._rp) / sqrt (2._rp)
          end if

        case ('rectangular')  !--n_base = 2

          h = -(branch % d) / 2._rp

      end select

    else  !--n_base = 1, so j = 1

      h = 0._rp

    end if

  else if ( j > n_height + n_base) then

    if (add_cap) then

      select case (cap_shape)
        case ('hemispherical')  !--n_cap = 2
          
          if (j == n_tot - 1) then
            h = l + (1._rp - branch % taper) * ((branch % d) / 2._rp)  &
                    / sqrt (2._rp)
          else  !--j = n_tot
            h = l + (1._rp - branch % taper) * (branch % d) / 2._rp
          end if
    
        case ('rectangular')  !--c_cap = 2

          h = l + (1._rp - branch % taper) * (branch % d) / 2._rp

      end select

    else  !--n_cap = 1, so j = n_tot

      h = l
      
    end if

  else

    h = (j - n_base - 1) * dh

  end if

  end subroutine set_h
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !--not that taper is not included with caps
  !--depends on h: call AFTER set_h
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine set_radius ()
  implicit none

  if (j <= n_base) then

    if (add_base_this) then

      select case (base_shape)
        case ('hemispherical')  !--n_base is 2

          if (j == 1) then
            r = 0._rp
          else
            r = r0 / sqrt (2._rp)
          end if

        case ('rectangular')  !--n_base is 2

          if (j == 1) then
            r = 0._rp
          else
            r = r0
          end if

      end select

    else  !--n_base is 1

      r = 0._rp

    end if
    
  else if (j > n_height + n_base) then

    if (add_cap) then
    
      select case (cap_shape)
        case ('hemispherical')  !--n_cap is 2

          if (j == n_tot - 1) then
            r = r_top / sqrt (2._rp)
          else
            r = 0._rp
          end if

        case ('rectangular')  !--n_cap is 2

          if (j == n_tot - 1) then
            r = r_top
          else
            r = 0._rp
          end if

      end select

    else  !--n_cap is 1

      r = 0._rp

    end if

  else

    r = r0 * (1._rp - (h / l) * (branch % taper))  !--add the taper

  end if

  end subroutine set_radius

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine write_tecplot_zone_hdr ()
  implicit none

  write (lun, '(a,i0,a,i0,a)') 'zone, f=point, i = ', n_face + 1,  &
                               ', j = ', n_tot, ', k = 1'
                                                           
  end subroutine write_tecplot_zone_hdr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine draw_branches

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine draw_tree_array (name)
implicit none

character (*), intent (in) :: name
character(1024) :: msg
character (*), parameter :: sub_name = mod_name //      &
                                           '.draw_tree_array'

integer, parameter :: lun = 1

integer :: i

logical :: opn

!----------------------------------------------------------------------
inquire (unit = lun, opened = opn)
if (opn) then
  write (msg, '(a,i0,a)') 'unit ', lun, ' is already open'
  call error (sub_name, msg)
else
  open (lun, file = name)
end if

if (use_tecplot) then  ! write header
  write (lun, *) 'variables = "x" "y" "z"'
end if

do i = 1, n_tree
  call draw_branches (tree_array(i) % trunk, lun)    
end do

close (lun)

end subroutine draw_tree_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--read in branch data from file
!--assumes the tree data structure already exists, only some things
!  need to be read from file, e.g. velscale
!--this must be kept in sync with write_br_data
!--does not care about order
!--we make allowance for some spaces here (not really needed)
!--some error checking in here to let us know in case of error
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
recursive subroutine read_br_data (br, lun)
implicit none

type (branch_type), intent (inout) :: br

integer, intent (in) :: lun

character (*), parameter :: sub = mod_name // '.read_br_data'
                            !--shortened name since used alot here

character (256) :: line

integer :: i
integer :: ieq
integer :: itmp

logical :: ltmp

!----------------------------------------------------------------------

read (lun, '(a)') line
if (index (line, 'begin branch') == 0) then
  call error (sub, "no 'begin branch' statement")
end if

do

  read (lun, '(a)') line

  !--is this the end of this branch?
  if (index (line, 'end branch') > 0) exit

  ieq = index (line, '=') 
  !--is there an equals sign?
  !  if not, then this line does not contain useful data (except
  !  maybe for error checking), so skip to next line
  if (ieq == 0) cycle  !--read next line
 
  !--the commented reads here indicate that the value should be known
  !  from trees_setup anyway, so that case is a do-nothing
  select case ( trim (line(1:ieq-1)) )
  case ('ident')
    read (line(ieq+1:), *) itmp
    if (itmp /= br % ident) call error (sub, 'ident mismatch')
  case ('itree')
    read (line(ieq+1:), *) itmp
    if (itmp /= br % itree) call error (sub, 'itree mismatch')
  case ('gen')
    read (line(ieq+1:), *) itmp
    if (itmp /= br % gen) call error (sub, 'gen mismatch')
  case ('n_sub_branch')
    read (line(ieq+1:), *) itmp
    if (itmp /= br % n_sub_branch) call error (sub, 'n_sub_branch mismatch')
  case ('zone')
    read (line(ieq+1:), *) itmp
    if (itmp /= br % zone) call warn (sub, 'zone mismatch')
                                !--difference zones allowed
  case ('nbboxpt')
    !read (line(ieq+1:), *) br % nbboxpt
  case ('nrespt')
    !read (line(ieq+1:), *) br % nrespt
  case ('resolved')
    read (line(ieq+1:), *) ltmp
    if (ltmp .neqv. br % resolved) call error (sub, 'resolved mismatch')
  case ('l')
    !read (line(ieq+1:), *) br % l
  case ('d')
    !read (line(ieq+1:), *) br % d
  case ('taper')
    !read (line(ieq+1:), *) br % taper
  case ('twist')
    !read (line(ieq+1:), *) br % twist
  case ('width_bbox')
    !read (line(ieq+1:), *) br % width_bbox
  case ('height_bbox')
    !read (line(ieq+1:), *) br % height_bbox
  case ('root_height')
    !read (line(ieq+1:), *) br % root_height
  case ('fcoeff')
    !read (line(ieq+1:), *) br % fcoeff
  case ('fdist_coeff')
    !read (line(ieq+1:), *) br % fdist_coeff
  case ('A')
    read (line(ieq+1:), *) br % A
  case ('fcoeff_dir')
    !read (line(ieq+1:), *) br % fcoeff_dir
  case ('fnorm')
    !read (line(ieq+1:), *) br % fnorm
  case ('Mdyn')
    !read (line(ieq+1:), *) br % Mdyn
  case ('x0')
    !read (line(ieq+1:), *) br % x0
  case ('abs_dir')
    !read (line(ieq+1:), *) br % abs_dir
  case ('rel_dir')
    !read (line(ieq+1:), *) br % rel_dir
  case ('resf')
    read (line(ieq+1:), *) br % resf
  case ('unresf')
    read (line(ieq+1:), *) br % unresf
  case ('resftot')
    read (line(ieq+1:), *) br % resftot
  case ('unresftot')
    read (line(ieq+1:), *) br % unresftot
  case ('ftot')
    read (line(ieq+1:), *) br % ftot
  case ('velscale')
    read (line(ieq+1:), *) br % velscale
  case ('x_hat')
    !read (line(ieq+1:), *) br % x_hat
  case ('y_hat')
    !read (line(ieq+1:), *) br % y_hat
  case ('z_hat')
    !read (line(ieq+1:), *) br % z_hat
  case default
    call error ( sub, 'invalid data line' // n_l // trim (line) )
  end select

end do

if (associated (br % sub_branch)) then

  do i = 1, br % n_sub_branch

    call read_br_data (br % sub_branch(i), lun)

  end do

end if

end subroutine read_br_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine read_tree_data (tree, lun)
implicit none

type (tree_type), intent (inout) :: tree

character (*), parameter :: sub_name = mod_name // '.read_tree_data'

integer, intent (in) :: lun

character (256) :: line

!----------------------------------------------------------------------

read (lun, '(a)') line

if (index (line, 'begin tree') == 0) call error (sub_name,             &
                                                 "missing 'begin tree'")

call read_br_data (tree % trunk, lun)

read (lun, '(a)') line

if (index (line, 'end tree') == 0) call error (sub_name,           &
                                               "missing 'end tree'")

end subroutine read_tree_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine read_ta_data (name)
implicit none

character (*), intent (in) :: name
character(1024) :: msg
character (*), parameter :: sub_name = mod_name // '.read_ta_data'

integer, parameter :: lun = 1

character (256) :: line

integer :: i

logical :: opn, exst

!----------------------------------------------------------------------

if (.not. tree_array_initialized) then
  call error (sub_name, 'tree_array not initialized')
end if

inquire (unit=lun, opened=opn)
if (opn) then
  write (msg, '(a,i0,a)') 'unit ', lun, ' is already open'
  call error (sub_name, msg)
end if

inquire (file=name, exist=exst, opened=opn)
if (.not. exst) then
  call error (sub_name, 'file ' // trim (name) // ' does not exist')
else if (opn) then
  call error (sub_name, 'file ' // trim (name) // 'is already open')
end if

open (lun, file=name, action='read')

do i = 1, n_tree

  read (lun, '(a)') line
  
  !--this check is not precise
  if ((index (line, 'begin tree') == 0) .or.                &
      (index (line, 'of') < index (line, 'begin tree'))) then
    call error (sub_name, 'invalid begin tree i of n_tree line' // n_l //  &
                          'line =' // line)
  end if  

  call read_tree_data (tree_array(i), lun)

end do

close (lun)

end subroutine read_ta_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--dumps all branch data to file
!--may want to make format nicer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
recursive subroutine write_br_data (br, lun)
implicit none

type (branch_type), intent (in) :: br

integer, intent (in) :: lun

character (*), parameter:: ifmt = 'i0'
character (*), parameter :: rfmt = 'es13.6'

character (64) :: fmt

integer :: i

!----------------------------------------------------------------------

write (lun, '(a)') 'begin branch'

fmt = '(a,' // ifmt // ')'

write (lun, fmt) 'ident = ', br % ident
write (lun, fmt) 'itree = ', br % itree
write (lun, fmt) 'gen = ', br % gen
write (lun, fmt) 'n_sub_branch = ', br % n_sub_branch
write (lun, fmt) 'zone = ', br % zone
write (lun, fmt) 'nbboxpt = ', br % nbboxpt
write (lun, fmt) 'nrespt = ', br % nrespt

if (associated (br % bboxpt)) then
  write (lun, '(a)') 'branch % bboxpt is associated'
  
  !write (lun, '(a)') 'br % bboxpt:'
  !fmt = '(4(' // ifmt // ',1x))'
  !do i = 1, br % nbboxpt
  !  write (lun, fmt) i, br % bboxpt(:, i)
  !end do
else
  write (lun, '(a)') 'br % bboxpt is NOT associated'
end if

write (lun, '(a,l1)') 'resolved = ', br % resolved

fmt = '(a,' // rfmt // ')'  !--for reals

write (lun, fmt) 'l = ', br % l
write (lun, fmt) 'd = ', br % d
write (lun, fmt) 'taper = ', br % taper
write (lun, fmt) 'twist = ', br % twist
write (lun, fmt) 'width_bbox = ', br % width_bbox
write (lun, fmt) 'height_bbox = ', br % height_bbox
write (lun, fmt) 'root_height = ', br % root_height

write (fmt, '(a,i0,a)') '(a,', nfcoeff, '(' // rfmt // ',1x))'

write (lun, fmt) 'fcoeff = ', br % fcoeff
write (lun, fmt) 'fdist_coeff = ', br % fdist_coeff
write (lun, fmt) 'A = ', br % A

write (fmt, '(a,i0,a)') '(a,', nd * nfcoeff, '(' // rfmt // ',1x))'

write (lun, fmt) 'fcoeff_dir = ', br % fcoeff_dir
write (lun, fmt) 'fnorm = ', br % fnorm

fmt = '(a,3(' // rfmt // ',1x))'  !--could write nd instead of 3

write (lun, fmt) 'Mdyn = ', br % Mdyn
write (lun, fmt) 'x0 = ', br % x0
write (lun, fmt) 'abs_dir = ', br % abs_dir
write (lun, fmt) 'rel_dir = ', br % rel_dir
write (lun, fmt) 'resf = ', br % resf
write (lun, fmt) 'unresf = ', br % unresf
write (lun, fmt) 'resftot = ', br % resftot
write (lun, fmt) 'unresftot = ', br % unresftot
write (lun, fmt) 'ftot = ', br % ftot
write (lun, fmt) 'velscale = ', br % velscale
write (lun, fmt) 'x_hat = ', br % x_hat
write (lun, fmt) 'y_hat = ', br % y_hat
write (lun, fmt) 'z_hat = ', br % z_hat

if (associated (br % sub_branch)) then
  write (lun, '(a)') 'sub_branch is associated'
else
  write (lun, '(a)') 'sub_branch is NOT associated'
end if

if (associated (br % parent_branch)) then
  write (lun, '(a)') 'parent_branch is associated'
else
  write (lun, '(a)') 'parent_branch is NOT associated'
end if

write (lun, '(a)') 'end branch'

if (associated (br % sub_branch)) then

  do i = 1, br % n_sub_branch

    call write_br_data (br % sub_branch(i), lun)

  end do

end if

end subroutine write_br_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_tree_data (tree, lun)
implicit none

type (tree_type), intent (in) :: tree

character (*), parameter :: sub_name = mod_name // '.write_tree_data'

integer, intent (in) :: lun

!----------------------------------------------------------------------

write (lun, '(a)') 'begin tree'

call write_br_data (tree % trunk, lun)

write (lun, '(a)') 'end tree'

end subroutine write_tree_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_ta_data (name)
implicit none

character (*), intent (in) :: name
character(1024) :: msg
character (*), parameter :: sub_name = mod_name // '.write_ta_data'

integer, parameter :: lun = 1

integer :: i

logical :: opn

!----------------------------------------------------------------------

if (.not. tree_array_initialized) then
  call error (sub_name, 'tree_array not initialized')
end if

inquire (unit = lun, opened = opn)
if (opn) then
  write (msg, '(a,i0,a)') 'unit ', lun, ' is already open'
  call error (sub_name, msg)
else
  open (lun, file = name)
end if

do i = 1, n_tree

  write (lun, '(a,i0,a,i0)') 'begin tree ', i, ' of ', n_tree 

  call write_tree_data (tree_array(i), lun)

end do

close (lun)

end subroutine write_ta_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--really only for debugging: lets us know when forcing is not working
!--commented out to remove sim_param dependence from this module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine output_tree_array_grid_vel (name)
!implicit none
!
!character (*), intent (in) :: name
!
!character (*), parameter :: sub_name = mod_name //                 &
!                                       '.output_tree_array_grid_vel'
!
!integer, parameter :: lun = 1
!
!integer :: i
!
!logical :: opn
!
!!----------------------------------------------------------------------
!
!inquire (unit = lun, opened = opn)
!if (opn) then
!  write (msg, '(a,i0,a)') 'unit ', lun, ' is already open'
!  call error (sub_name, msg)
!else
!  open (lun, file = name)
!end if
!
!do i = 1, n_tree
!
!  write (lun, *) 'begin tree ', i, ' of ', n_tree
!
!  call output_tree_grid_vel (tree_array(i), lun)
!
!end do
!
!close (lun)
!
!end subroutine output_tree_array_grid_vel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine output_tree_grid_vel (tree, lun)
!implicit none
!
!type (tree_type), intent (in) :: tree
!
!integer, intent (in) :: lun
!
!character (*), parameter :: sub_name = mod_name //           &
!                                       '.output_tree_grid_vel'
!
!!----------------------------------------------------------------------
!
!write (lun, *) 'begin tree'
!
!call warn (sub_name, 'output_branch_grid_vel has been disabled')
!!call output_branch_grid_vel (tree % trunk, lun)
!
!end subroutine output_tree_grid_vel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! outputs grid velocity at the nodes of the branches
! really only for debugging, lets us know if our forcing is messed up
!--careful: adding this will add sim_param dependence to this file,
!  which is not good for trees_pre
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!recursive subroutine output_branch_grid_vel (branch, lun)
!use sim_param, only : u, v, w
!
!type (branch_type), intent (in) :: branch
!
!integer, intent (in) :: lun
!
!character (*), parameter :: sub_name = mod_name //             &
!                                       '.output_branch_grid_vel'
!
!integer :: i, j, k
!integer :: i_pt
!
!!----------------------------------------------------------------------
!
!if (branch % n_pt < 1) call error (sub_name, 'n_pt < 1')
!
!write (lun, *) 'begin branch'
!
!do i_pt = 1, branch % n_pt
!
!  i = branch % jx(1, i_pt)
!  j = branch % jx(2, i_pt)
!  k = branch % jx(3, i_pt)
!
!  write (lun, *) i, j, k, u(i, j, k), v(i, j, k),   &
!                 0.5_rp * (w(i, j, k) + w(i, j, k+1))
!
!end do
!
!if (associated (branch % sub_branch)) then
!
!  do i = 1, branch % n_sub_branch
!    call output_branch_grid_vel (branch % sub_branch(i), lun)
!  end do
!
!end if
!
!end subroutine output_branch_grid_vel
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module trees_io_ls
