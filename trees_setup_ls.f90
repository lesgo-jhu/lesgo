module trees_setup_ls
use trees_base_ls
implicit none

save
private

public ::  sdistfcn_tree_array, fill_tree_array

character (*), parameter :: mod_name = 'trees_setup'

integer :: ident = 0 !--global index used to set branch % ident 

real (rp), parameter :: zero_vector(nd) = 0._rp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--calculates signed distance function for whole tree array
!--could use fmm to speed things up
!--when brident is present, it stores the branch % ident for pts
!  inside the branches and is -1 otherwise
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine sdistfcn_tree_array (ip, mx, phi, brident)
implicit none

integer, intent (in) :: ip  !--specifies chunk (or process)
                            !  local chunk k=0 <-> ip*(mx(3)-1)
                            !  local chunk k=nz <-> ip*(mx(3)-1)+mx(3)
integer, intent (in) :: mx(nd)

real (rp), intent (out) :: phi(mx(1), mx(2), 0:mx(3))  !--note 0 here
integer, intent (out), optional :: brident(mx(1), mx(2), mx(3))

character (*), parameter :: sub = mod_name // '.sdistfcn_tree_array'

real (rp), parameter :: phi_init = huge (1._rp)
real (rp), parameter :: phi_thresh = 10._rp * epsilon (0._rp)

integer :: i

!---------------------------------------------------------------------
write(*,*) 'epsilon(0) = ', epsilon(0.);
if (.not. grid % initialized) then
  call error (sub, 'grid must be initialized')
end if

!--check size of phi array vs. size of grid
do i = 1, nd-1
  if (mx(i) < grid % nx(i)) then
    call error (sub, 'phi is too small in dimension', i)
  end if
  !--not sure how to check third dimension now with the chunks
end do

!--init phi to large values
phi = phi_init

if (present (brident)) brident = -1

do i = 1, n_tree

  call mesg (sub, 'trunk % abs_dir =', tree_array(i) % trunk % abs_dir)
  call mesg (sub, 'trunk % x_hat =', tree_array(i) % trunk % x_hat)
  call mesg (sub, 'trunk % y_hat =', tree_array(i) % trunk % y_hat)
  call mesg (sub, 'trunk % z_hat =', tree_array(i) % trunk % z_hat)

  if (present (brident)) then
    call sdistfcn_branch (tree_array(i) % trunk, ip, mx, phi, brident)
  else
    call sdistfcn_branch (tree_array(i) % trunk, ip, mx, phi)
  end if
  
end do

!--where phi is 0 < phi < thresh, set phi to 0
!--this is just experimental--a better way must exist
!where ((0._rp < phi) .and. (phi < phi_thresh)) phi = 0._rp
where ((abs (phi) < phi_thresh)) phi = 0._rp

end subroutine sdistfcn_tree_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--assumes grid is initialized, and is not larger than phi
!--puts phi on level_set_node
!--this is slow
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
recursive subroutine sdistfcn_branch (br, ip, mx, phi, brident)
implicit none

type (branch_type), intent (in) :: br
integer, intent (in) :: ip
integer, intent (in) :: mx(3)
real (rp), intent (inout) :: phi(mx(1), mx(2), 0:mx(3))
                             !--only 1:nx, 1:ny used
                             !--note 0
integer, intent (inout), optional :: brident(mx(1), mx(2), mx(3))  

character (*), parameter :: sub_name = mod_name // '.sdistfcn_branch'

integer, parameter :: phi_node = 1  !--u-nodes

real (rp), parameter :: thresh = 10._rp * epsilon (0._rp)

integer :: i, j, k, ktot

real (rp) :: cosine, sine
real (rp) :: d_para, d_perp
real (rp) :: dist
real (rp) :: dist_sq
real (rp) :: l
real (rp) :: rb0, rb, rbl
real (rp) :: t, tgt
real (rp) :: x_para(nd), x_perp(nd)
real (rp) :: x(nd)

!---------------------------------------------------------------------

if (.not. br % resolved) return

l = br % l
rb0 = 0.5_rp * (br % d)  !--base radius
t = br % taper
rbl = rb0 * (1._rp - t)  !--top radius
tgt = rb0 * t / l
sine = rb0 * t / sqrt ((rb0 * t)**2 + l**2)
cosine  = l / sqrt ((rb0 * t)**2 + l**2)

!--not sure how openmp will handle the internal sub variables and
!  present (brident)
!$omp parallel do                                             &
!$omp private(x,d_para,x_para,d_perp,x_perp,rb,dist_sq,dist)  &
!$omp shared(phi, brident)
do k = 0, mx(3)  !--mx(3) should be (grid % nx(3) - 1) / (np) + 1

  ktot = ip * (mx(3) - 1) + k

  call mesg (sub_name, 'processing k =', ktot)

  x(3) = pt_of_grid (ktot, 3, phi_node) - br % x0(3)

  do j = 1, mx(2)  !--mx(2) should be grid % nx(2)

    x(2) = pt_of_grid (j, 2, phi_node) - br % x0(2)

    do i = 1, mx(1)  !--mx(3) should be grid % nx(1)

      x(1) = pt_of_grid (i, 1, phi_node) - br % x0(1)

      !--this part depends on cross section
      d_para = dot_product (x, br % abs_dir)
      x_para = d_para * (br % abs_dir)

      x_perp = x - x_para
      d_perp = mag (x_perp)

      !--set rb with taper? or explicitly code in?
      rb = rb0 * (1._rp - (d_para / l) * t)  !--careful when d_para < 0 or
                                             !  when d_para > l

      select case (branch_cross_section)
        case ('circular')
          call dist_circle ()
        case ('square')
          call dist_square ()
        case ('square+plate')
          call dist_square ()
          dist_sq = dist
          call dist_plate ()
          dist = min (dist_sq, dist)
        case default
          call error (sub_name, 'invalid branch_cross_section')
      end select

      !--need to take sign into account here
      if (phi(i, j, k) >= 0._rp) then  !--cannot modify interior values...

        if (dist < phi(i, j, k)) then
        
          phi(i, j, k) = dist

          if (present (brident) .and. (dist <= 0._rp)) then
            if (k > 0) brident(i, j, k) = br % ident
                !--current convention is weird, phi from 0, brident from 1
                !  however, this is assumed in the main code so must change
                !  all or none
          end if

        end if

      end if

    end do

  end do

end do
!$omp end parallel do

if (associated (br % sub_branch)) then

  do i = 1, br % n_sub_branch
    if (present (brident)) then
      call sdistfcn_branch (br % sub_branch(i), ip, mx, phi, brident)
    else
      call sdistfcn_branch (br % sub_branch(i), ip, mx, phi)
    end if
  end do

end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !--puts in a flat plate shaped sort of like a stealth bomber:
  !          *
  !        *****
  !      *********
  !    *************
  !  *****************
  !    *****   *****
  !      *       *
  !--this is supposed to simulate the case where the SGS branches
  !  basically fill in the space they occupy
  !--only does something when br % gen = n_gen
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine dist_plate ()
  implicit none

  real (rp) :: d1, d2, d3
  real (rp) :: xi_part(nd)
  real (rp) :: u(nd), u1(nd), u2(nd), u3(nd)  !--temp variables
  real (rp) :: xp(nd), xi(nd)
  real (rp) :: xi_hat(nd), eta_hat(nd), zeta_hat(nd)
  real (rp) :: hdepth, s  !--dimensions of plate shape
  real (rp) :: ratio
  
  !-------------------------------------------------------------------

  if (br % gen /= tree_array(br % itree) % n_gen) return

  ratio = tree_array(br % itree) % ratio
  
  hdepth = ratio * rb0  !--half-depth
  s = (rb0 + (l + rb0) * ratio / (1._rp - ratio)) / sqrt (2._rp)
  !--s is the half length of one of the plate sides

  !--assumes x is rel. to x0 already
  !--make xp the branch-local coords
  xp(1) = dot_product (x, br % x_hat)
  xp(2) = dot_product (x, br % y_hat)
  xp(3) = dot_product (x, br % z_hat)

  !--now make rel. to center of plate shape
  xp = xp - l * (/ 0._rp, 0._rp, 1._rp /)  !--along the branch-local z-axis

  !--these are relative to branch directions (x_hat, y_hat, z_hat)
  xi_hat = (/ 1._rp, 0._rp, 0._rp /)  !--already normalized

  eta_hat = (/ 0._rp, 1._rp, -1._rp /)
  eta_hat = eta_hat / mag (eta_hat)

  zeta_hat = (/ 0._rp, 1._rp, 1._rp /)
  zeta_hat = zeta_hat / mag (zeta_hat)
  
  !--rotate again to get 45 degree clockwise rotation in x_hat-z_hat plane
  xi(1) = dot_product (xp, xi_hat)
  xi(2) = dot_product (xp, eta_hat)
  xi(3) = dot_product (xp, zeta_hat)

  !--redefine xi_hat, etc. to be in their own coordinates now
  xi_hat = (/ 1._rp, 0._rp, 0._rp /)
  eta_hat = (/ 0._rp, 1._rp, 0._rp /)
  zeta_hat = (/ 0._rp, 0._rp, 1._rp /)

  !--to be subtracted from u vectors
  xi_part = sign (1._rp, xi(1)) * min (abs (xi(1)), hdepth) *  &
             xi_hat

  !--most of the inequalities here are "greedy": they include boundaries
  !  this ensures that points do not get missed, and is not a problem since
  !  the distance function is continuous across the boundaries anyway
  if (xi(2) >= s .and. (0._rp <= xi(3) .and. xi(3) <= s)) then
    !--case 2
    u = xi - s * eta_hat - xi_part
    dist = mag (u - dot_product (u, zeta_hat) * zeta_hat)

  else if ((-s <= xi(2) .and. xi(2) <= 0._rp) .and. xi(3) <= -s) then
    !--case 3
    u = xi + s * zeta_hat - xi_part
    dist = mag (u - dot_product (u, eta_hat) * eta_hat)

  else if ((0._rp <= xi(2) .and. xi(2) <= s) .and.  &
           (-s <= xi(3) .and. xi(3) <= 0._rp)) then
    !--case 4, 5
    if ( abs (xi(3)) >= abs (xi(2)) ) then
      !--case 4
      u = xi - xi_part
      dist = mag (u - dot_product (u, zeta_hat) * zeta_hat)
    else
      !--case 5
      u = xi - xi_part
      dist = mag (u - dot_product (u, eta_hat) * eta_hat)
    end if
    
  else if ((xi(2) >= 0._rp .and. xi(3) <= -s) .and.  &
           ( abs (xi(3)) >= abs (xi(2)) )) then
    !--case 6
    u = xi + s * zeta_hat - xi_part
    dist = mag (u)  !--no proj

  else if ((xi(2) >= s .and. xi(3) <= 0._rp) .and.  &
           ( abs (xi(3)) <= abs (xi(2)) )) then
    !--case 7
    u = xi - s * eta_hat - xi_part
    dist = mag (u)  !--no proj

  else if (xi(2) <= -s .and. xi(3) <= -s) then
    !--case 8
    u = xi + s * (eta_hat + zeta_hat) - xi_part
    dist = mag (u)  !--no proj

  else if (xi(2) >= s .and. xi(3) >= s) then
    !--case 9
    u = xi - s * (eta_hat + zeta_hat) - xi_part
    dist = mag (u)  !--no proj

  else if (xi(2) <= -s .and. abs (xi(3)) <= s) then
    !--case 10
    u = xi + s * (eta_hat + zeta_hat) - xi_part
    dist = mag (u - dot_product (u, zeta_hat) * zeta_hat)

  else if (abs (xi(2)) <= s .and. xi(3) >= s) then
    !--case 11
    u = xi - s * (eta_hat + zeta_hat) - xi_part
    dist = mag (u - dot_product (u, eta_hat) * eta_hat)

  else if (xi(2) <= -s .and. xi(3) >= s) then
    !-case 12
    u = xi - s * (-eta_hat + zeta_hat) - xi_part
    dist = mag (u)  !--no proj

  else  !--case 1

    if (abs (xi(1)) >= hdepth) then  !--outside

      dist = abs (xi(1)) - hdepth

    else  !--inside

      !--set dist after this if-block
      if ((0._rp <= xi(2) .and. xi(2) <= s) .and.  &
          (0._rp <= xi(3) .and. xi(3) <= s)) then
        !--quadrant 1
        u1 = xi - xi_part
        d1 = mag (u1 - dot_product (u1, eta_hat) * eta_hat)

        u2 = xi - s * zeta_hat - xi_part
        d2 = mag (u2 - dot_product (u2, eta_hat) * eta_hat)

        u3 = xi - s * eta_hat - xi_part
        d3 = mag (u3 - dot_product (u3, zeta_hat) * zeta_hat)

       else if ((-s <= xi(2) .and. xi(2) <= 0._rp) .and.  &
                (0._rp <= xi(3) .and. xi(3) <= s)) then
         !--quadrant 2
         u1 = xi - xi_part
         d1 = mag (u1)  !--no proj

         u2 = xi + s * eta_hat - xi_part
         d2 = mag (u2 - dot_product (u2, zeta_hat) * zeta_hat)

         u3 = xi - s * zeta_hat - xi_part
         d3 = mag (u3 - dot_product (u3, eta_hat) * eta_hat)

       else
          !--quadrant 3
          u1 = xi - xi_part
          d1 = mag (u1 - dot_product (u1, zeta_hat) * zeta_hat)

          u2 = xi + s * eta_hat - xi_part
          d2 = mag (u2 - dot_product (u2, zeta_hat) * zeta_hat)

          u3 = xi + s * zeta_hat - xi_part
          d3 = mag (u3 - dot_product (u3, eta_hat) * eta_hat)

        end if

        dist = -min ( d1, d2, d3, hdepth - abs (xi(1)) )

    end if
    
  end if

  end subroutine dist_plate

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine dist_circle ()
  implicit none

  if (add_base .and. add_cap) then
    call dist_circle_bc ()
  else if (add_cap) then  !--add_base is false
    call dist_circle_c ()
  else if (add_base) then  !--add_cap is false
    call dist_circle_b ()
  else  !--no bases/caps
    call dist_circle_nobc ()
  end if
  
  end subroutine dist_circle

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine dist_circle_b ()
  implicit none

  call error (sub_name, 'dist_circle_b not implemented')

  end subroutine dist_circle_b

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine dist_circle_c ()
  implicit none

  call error (sub_name, 'dist_circle_c not implemented')

  end subroutine dist_circle_c
 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine dist_circle_bc ()
  implicit none

  !integer :: code

  !-------------------------------------------------------------------

  if (d_para <= 0._rp) then

    if (base_shape == 'hemispherical') then

      dist = mag (x) - rb0
      
    else if (base_shape == 'rectangular') then

      !--special case: allowing so-called rectangular base with
      !  circular cross-section
      if (d_para > -rb0) then

        if (-d_para > d_perp) then
          dist = -(rb0 + d_para)
        else
          dist = d_perp - rb0
        end if

      else

        if (d_perp < rb0) then
          dist = -(rb0 + d_para)
        else
          dist = sqrt ((d_perp - rb0)**2 + (d_para + rb0)**2)
        end if
        
      end if

    else
      call error (sub_name, 'invalid base_shape')
    end if
    !code = 1
    
  else if ((d_para < (d_perp - rb0) * tgt) .and. (d_perp > rb0)) then
    !--know d_para > 0 here

    dist = sqrt ((d_perp - rb0)**2 + d_para**2)
    !code = 2
    
  else if ((d_para < l) .and. (d_perp < rb)) then
    !--know d_para > 0 here

    if ((d_perp < rbl) .and. (l - d_para < (rbl - d_perp) * tgt)) then

      dist = -sqrt ((l - d_para)**2 + (rbl - d_perp)**2)
      !code = 3
      
    else

      dist = - (rb - d_perp) * cosine
      !code = 4
      
    end if
    
  else if ((d_para >= l) .and. (d_perp**2 + (d_para - l)**2 <= rbl**2)) then

    dist = sqrt (d_perp**2 + (d_para - l)**2) - rbl
    !code = 5
    
  else if ((d_para >= l) .and. (d_para - l >= (d_perp - rbl) * tgt)) then

    dist = sqrt (d_perp**2 + (d_para - l)**2) - rbl
    !code = 6
    
  else

    dist = (d_perp - rb) * cosine
    !code = 7
    
  end if

 !write (*, '(1x,a,4(1x,i0),es12.5)') 'i, j, k, code, dist =',  &
 !                                    i, j, k, code, dist

  end subroutine dist_circle_bc

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine dist_circle_nobc ()
  implicit none

  call error (sub_name, 'dist_circle_nobc not implemented yet')

  end subroutine dist_circle_nobc
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine dist_square ()
  implicit none

  if (add_base .and. add_cap) then
    call dist_square_bc ()
  else if (add_cap) then
    call dist_square_c ()
  else if (add_base) then
    call dist_square_b ()
  else
    call dist_square_nobc ()
  end if

  end subroutine dist_square

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine dist_square_b ()
  implicit none

  call error (sub_name, 'dist_square_b not implemented yet')
  
  end subroutine dist_square_b

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine dist_square_c ()
  implicit none

  call error (sub_name, 'dist_square_c not implemented yet')
  
  end subroutine dist_square_c

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine dist_square_bc ()
  implicit none

  real (rp) :: c(nd), u(nd)
  real (rp) :: s1, s2
  real (rp) :: x1, x2
  
  !-------------------------------------------------------------------

  x1 = dot_product (x, br % x_hat)
  x2 = dot_product (x, br % y_hat)

  s1 = sign (1._rp, x1)
  s2 = sign (1._rp, x2)

  if (d_para < -rb0) then
    !--below the base
    if ((abs (x1) > rb0) .and. (abs (x2) > rb0)) then
      !--closest point is lower corners of base

      c = rb0 * (-(br % z_hat) + s1 * (br % x_hat) + s2 * (br % y_hat))

      dist = mag (x - c)

    else if ((abs (x1) <= rb0) .and. (abs (x2) > rb0)) then

       c = rb0 * (-(br % z_hat) + s2 * (br % y_hat))
       c = x - c
       c = c - (dot_product (br % x_hat, c)) * (br % x_hat)
       dist = mag (c)
       
    else if ((abs (x1) > rb0) .and. (abs (x2) <= rb0)) then

      c = rb0 * (-(br % z_hat) + s1 * (br % x_hat))
      c = x - c
      c = c - (dot_product (br % y_hat, c)) * (br % y_hat)
      dist = mag (c)

    else

      dist = -(d_para + rb0)

    end if

  else if (d_para < 0._rp) then  !--know -rb0 <= d
    !--level with the base
    if ((abs (x1) > rb0) .and. (abs (x2) > rb0)) then

      c = rb0 * (s1 * (br % x_hat) + s2 * (br % y_hat))
      dist = mag (x_perp - c)  !--x_perp, not x

    else if ((abs (x1) <= rb0) .and. (abs (x2) > rb0)) then

      dist = abs (x2) - rb0

    else if ((abs (x1) > rb0) .and. (abs (x2) <= rb0)) then

      dist = abs (x1) - rb0

    else  !--inside base

      if ( abs (d_para) > max (abs (x1), abs (x2)) ) then

        dist = -(rb0 - abs (d_para))

      else
      
        dist = max (abs (x1), abs (x2)) - rb0

      end if
      
    end if

  else if (d_para > l + rbl) then
    !--above the cap
    if ((abs (x1) > rbl) .and. (abs (x2) > rbl)) then

      c = (l + rbl) * (br % z_hat) +                  &
          rbl * (s1 * (br % x_hat) + s2 * (br % y_hat))
      dist = mag (x - c)

    else if ((abs (x1) <= rbl) .and. (abs (x2) > rbl)) then

      c = (l + rbl) * (br % z_hat) + rbl * s2 * (br % y_hat)
      c = x - c
      c = c - dot_product (c, br % x_hat) * (br % x_hat)
      dist = mag (c)

    else if ((abs (x1) > rbl) .and. (abs (x2) <= rbl)) then

      c = (l + rbl) * (br % z_hat) + rbl * s1 * (br % x_hat)
      c = x - c
      c = c - dot_product (c, br % y_hat) * (br % y_hat)
      dist = mag (c)

    else

      dist = d_para - (l + rbl)

    end if

   else if (d_para >= l) then  !--know d_para <= l + rbl

     if ((abs (x1) > rbl) .and. (abs (x2) > rbl)) then

       c = rbl * (s1 * (br % x_hat) + s2 * (br % y_hat))
       c = x_perp - c
       
       if (d_para - l > mag (c) * tgt) then  !--closer to vertical

         dist = mag (c)

       else  !--closer to slanted line

         !u = (/ -s1 * sine, -s2 * sine, cosine /)
         u = (-s1 * sine) * (br % x_hat) + (-s2 * sine) * (br % y_hat) +  &
             cosine * (br % z_hat)

         u = u / mag (u)  !--unit vector along slanted corner line
         c = rb0 * (s1 * (br % x_hat) + s2 * (br % y_hat))
             !--reference corner slanted line passes through
         c = x - c
         c = c - dot_product (c, u) * u
         dist = mag (c)

       end if

     else if ((abs (x1) <= rbl) .and. (abs (x2) > rbl)) then

       if ( d_para - l > (abs (x2) - rbl) * tgt ) then  !--closer to vertical
         
         dist = abs (x2) - rbl

       else  !--closer to slanted plane

         dist = (abs (x2) - rb) * cosine  !--it is rb, not rbl!

       end if

     else if ((abs (x1) > rbl) .and. (abs (x2) <= rbl)) then

       if ( d_para - l > (abs (x1) - rbl) * tgt ) then  !--closer to vertical

         dist = abs (x1) - rbl

       else  !--closer to slanted plane

         dist = (abs (x1) - rb) * cosine  !--it is rb, not rbl!

       end if

     else  !--inside the cap

       if ( d_para - l > max (abs (x1), abs (x2)) ) then

         dist = -((l + rbl) - d_para)

       else

         dist = max (abs (x1), abs (x2)) - rbl

       end if
  
    end if

  else  !--d_para >= 0 & d_para <= l

    if ((abs (x1) > rb) .and. (abs (x2) > rb)) then

      !--check if closer to lower corner
      c = rb0 * (s1 * (br % x_hat) + s2 * (br % y_hat))

      if (d_para < mag (x_perp - c) * tgt) then  !--corner is closer

        dist = mag (x - c)  !--x here, not x_perp

      else  !--slanted line is closer

        !--this is in branch-local coordinates!
        !u = (/ -s1 * sine, - s2 * sine, cosine /)

        u = (-s1 * sine) * (br % x_hat) + (-s2 * sine) * (br % y_hat) +  &
            cosine * (br % z_hat)
        u = u / mag (u)

        c = rb0 * (s1 * (br % x_hat) + s2 * (br % y_hat))  !--ref. corner
        c = x - c
        c = c - dot_product (c, u) * u
        dist = mag (c)

      end if

    else if ((abs (x1) <= rb) .and. (abs (x2) > rb)) then

      if (d_para < (abs (x2) - rb0) * tgt) then
        !--line is closer
        dist = sqrt ((abs (x2) - rb0)**2 + d_para**2)

      else  !--slanted plane is closer

        dist = (abs (x2) - rb) * cosine  !--rb, not rb0

      end if

    else if ((abs (x1) > rb) .and. (abs (x2) <= rb)) then

      if (d_para < (abs (x1) - rb0) * tgt) then
        !--line is closer
        dist = sqrt ((abs (x1) - rb0)**2 + d_para**2)

      else  !--slanted plane is closer

        dist = (abs (x1) - rb) * cosine  !--rb, not rb0

      end if

    else  !--inside stem

      dist = (max (abs (x1), abs (x2)) - rb) * cosine

    end if

  end if
 
  end subroutine dist_square_bc

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine dist_square_nobc ()
  implicit none

  call error (sub_name, 'dist_square_nobc not implemented yet')

  end subroutine dist_square_nobc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine sdistfcn_branch

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! fills the module tree array
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fill_tree_array ()
implicit none

character (*), parameter :: sub_name = mod_name // '.fill_tree_array'

integer :: i

! for now, these parameters go here
!real (rp) :: l(n_tree)  !--experimental
!real (rp) :: d(n_tree)  !--experimental
!real (rp) :: x0(nd, n_tree)  !--experimental
!real (rp) :: y0(n_tree)

!----------------------------------------------------------------------

!--read tree configuration file, allocates tree_array
call read_trees_conf ()  !--always returns n_tree > 0

!if (grid_bound) then
!
!  if (grid % initialized) then  ! specify in terms of grid pts
!
!    !--(40, 16) for r=1/2 case (128^3)
!    !--(40, 20) for r=1/2.5 case (128^3)
!    !--(45, 18) for r=1/3 case (128^3)
!    l = 15 * (grid % dx(3))  ! length along trunk_dir
!    d = 6 * (grid % dx(1))  ! assume dx(1) = dx(2), length perp trunk_dir
!
!    ! note: all three directions may be set here, but only two are
!    !       actually used
!    !       parallel to trunk_dir is not used here, that is set below
!
!    ! just for extended case (Lx=2)
!    !x0(1, 1) = ((grid % nx(1)) / 4) * (grid % dx(1))
!    x0(1, 1) = ((grid % nx(1)) / 2) * (grid % dx(1))
!    !x0(1, 1:2) = ((grid % nx(1)) / 4) * (grid % dx(1))
!    !x0(1, 3:4) = (3 * (grid % nx(1)) / 4) * (grid % dx(1))
!
!    x0(2, 1) = ((grid % nx(2)) / 2) * (grid % dx(2))
!    !x0(2, (/ 1, 3 /)) = ((grid % nx(2)) / 4) * (grid % dx(2))
!    !x0(2, (/ 2, 4 /)) = (3 * (grid % nx(2)) / 4) * (grid % dx(2))
!
!    !x0(3, 1) = ((grid % nx(3)) / 2) * (grid % dx(3))
!    x0(3, 1) = ((grid % nx(3)) / 2 - 0.5_rp) * (grid % dx(3))
!               !--experiment: want to have x0 on u-node
!
!  else
!
!    call error (sub_name, 'for grid_bound, initialize grid first')
!
!  end if
!
!else  ! specify in terms of real values
!
!  l = 0.2_rp
!  d = 0.1_rp * l
!
!  !x0(1:2) = 0.25_rp
!  !x0(3:4) = 0.75_rp
!
!  !y0((/ 1, 3 /)) = 0.25_rp
!  !y0((/ 2, 4 /)) = 0.75_rp
!
!end if

! now make sure we start at edge of grid
! this overrides what is specified above
! old comment (still applies, but generalized to other directions now):
!! assume trees has its root at z = 0
!!tree % trunk % x0 = (/ x0, y0, 0._rp /)
!! setting x0(3) = 0. is problematic since the trunks horizontal children
!! will not lie on grid nodes.
!! experiment: shift tree up by dx(3)/2 (starts at u-node)
!! this means the trunks children should also lie on u-nodes

!call warn (sub_name, 'base of tree height may not be set correctly')

!i = maxloc (abs (trunk_dir), dim=1)
!
!if (grid_bound .and. grid % initialized) then
!
!  x0(i, :) = grid % x_min(i, tree_node)
!
!else
!
!  call warn (sub_name, 'grid not initialized, ' //       &
!                       'setting base of tree at height 0')
!  x0(i, :) = 0._rp
!
!end if

do i = 1, n_tree

  call init_tree (i, tree_array(i))

end do

end subroutine fill_tree_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--provides a valid n_tree > 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine read_trees_conf ()
use messages
use string_util, only : eat_whtspc
implicit none

character (*), parameter :: sub = 'read_trees_conf'

character (*), parameter :: ftrees_conf = 'trees.conf'
character (*), parameter :: comment = '#'
character (*), parameter :: ldelim = '('  !--no whitespace allowed
character (*), parameter :: rdelim = ')'  !--no whitespace allowed
character (*), parameter :: begin_tree_block = '{'
character (*), parameter :: end_tree_block = '}'
character (*), parameter :: esyntax = 'syntax error at line'

integer, parameter :: lun = 1
integer, parameter :: BUFF_LEN = 256

logical, parameter :: DEBUG = .true.

character (BUFF_LEN) :: buff

integer :: eqpos
integer :: i
integer :: i_tree
integer :: ios
integer :: line
integer :: n_sub_branch
integer :: state

logical :: exst
logical :: have_n_tree, have_n_sub_branch, have_rel_dir, have_root_height

!---------------------------------------------------------------------

if (DEBUG) call enter_sub (sub)

inquire (file=ftrees_conf, exist=exst)

if (exst) then
  open (lun, file=ftrees_conf, action='read')
else
  call error (sub, 'file ' // ftrees_conf // ' does not exist')
end if

line = 0
have_n_tree = .false.
have_n_sub_branch = .false.
state = 0
i_tree = 0

!--format description
!  * 'comment' is comment token, must be first nonblank character
!  * blank lines are allowed
!  * <variable name> = <variable value> set the variables: e.g.
!    n_tree = 4
!  * trees are defined using tree blocks of the form
!    tree = {
!      l = ...
!      d = ...
!      x0 = ...
!      ...
!    }
!    the {,} delimiters are reqd to be in the lines as shown

do

  read (lun, '(a)', iostat=ios) buff
  if (ios /= 0) exit

  line = line + 1

  !--remove leading/intermediate whitespace
  call eat_whtspc (buff)

  if (verify (buff, ' ') == 0) cycle  !--drop blank lines
  
  if (buff (1:len (comment)) == comment) cycle  !--drop comment lines

  if (state /= 0) then
    !--only legal way is to be at an end of tree block
    !--check for end_tree_block at start of buff: have removed whitespace
    if (index (buff, end_tree_block) == 1) then
      call mesg (sub, 'exiting tree block')
      state = 0  !--reset into normal state
      cycle
    end if
  end if
  
  !--isolate up to the first equals sign
  eqpos = index (buff, '=')

  if (eqpos == 0) then  !--for now, invalid format if no equals sign
    call error (sub, 'no equals sign in line', line) 
  end if
  
  if (len_trim (buff) == eqpos) then  !--invalid if nothing after equals
    call error (sub, 'nothing after equals sign in line', line) 
  end if

  select case (buff(1:eqpos-1))

    case ('d')

      if (DEBUG) call mesg (sub, 'd case selected')
      read (buff(eqpos+1:), *) tree_array(i_tree) % d
      call mesg (sub, 'read d =', tree_array(i_tree) % d)

    case ('l')

      if (DEBUG) call mesg (sub, 'l case selected')
      read (buff(eqpos+1:), *) tree_array(i_tree) % l
      call mesg (sub, 'read l =', tree_array(i_tree) % l)
      
    case ('max_res_gen')

      if (DEBUG) call mesg (sub, 'max_res_gen case selected')
      read (buff(eqpos+1:), *) tree_array(i_tree) % max_res_gen
      call mesg (sub, 'read max_res_gen =', tree_array(i_tree) % max_res_gen)
      
    case ('n_gen')

      if (DEBUG) call mesg (sub, 'n_gen case selected')
      read (buff(eqpos+1:), *) tree_array(i_tree) % n_gen
      call mesg (sub, 'read n_gen =', tree_array(i_tree) % n_gen)
      
    case ('n_sub_branch')

      if (DEBUG) call mesg (sub, 'n_sub_branch case selected')
      read (buff(eqpos+1:), *) n_sub_branch
      call mesg (sub, 'read n_sub_branch =', n_sub_branch)

      tree_array(i_tree) % n_sub_branch = n_sub_branch
      
      !--now we must allocate certain arrays
      allocate ( tree_array(i_tree) % root_height(n_sub_branch) )
      allocate ( tree_array(i_tree) % rel_dir(nd, n_sub_branch) )
      allocate ( tree_array(i_tree) % twist(n_sub_branch) )
      !--set default twist here for now: need to move somewhere
      tree_array(i_tree) % twist = 0._rp

      have_n_sub_branch = .true.

    case ('n_tree')

      if (DEBUG) call mesg (sub, 'n_tree case selected')
      call case_n_tree ()

    case ('ratio')

      if (DEBUG) call mesg (sub, 'ratio case selected')
      read (buff(eqpos+1:), *) tree_array(i_tree) % ratio
      call mesg (sub, 'read ratio =', tree_array(i_tree) % ratio)
      
    case ('rel_dir')

      if (DEBUG) call mesg (sub, 'rel_dir case selected')
      if (.not. have_n_sub_branch) then
        call error (sub, 'n_sub_branch must be specified before rel_dir')
      end if
      read (buff(eqpos+1:), *) tree_array(i_tree) % rel_dir
	  write(*,*) 'tree_array(i_tree) % rel_dir = ', tree_array(i_tree) % rel_dir
      do i = 1, n_sub_branch
        call mesg (sub, 'rel_dir =', tree_array(i_tree) % rel_dir(:, i))
      end do

      !--normalize
      do i = 1, n_sub_branch
        tree_array(i_tree) % rel_dir(:, i) =                          &
                               tree_array(i_tree) % rel_dir(:, i) /   &
                               mag (tree_array(i_tree) % rel_dir(:, i))
      end do

      have_rel_dir = .true.
   
    case ('root_height')

      if (DEBUG) call mesg (sub, 'root_height case selected')
      !--make n_sub_branch is known
      if (.not. have_n_sub_branch) then
        call error (sub, 'n_sub_branch must specified before root_height')
      end if
      read (buff(eqpos+1:), *) tree_array(i_tree) % root_height
      call mesg (sub, 'read root_height =', tree_array(i_tree) % root_height)

      have_root_height = .true.

    case ('taper')

      if (DEBUG) call mesg (sub, 'taper case selected')
      read (buff(eqpos+1:), *) tree_array(i_tree) % taper
      call mesg (sub, 'read taper =', tree_array(i_tree) % taper)

    case ('tree')

      if (DEBUG) call mesg (sub, 'tree case selected')
      call case_tree ()

    case ('trunk_dir')

      if (DEBUG) call mesg (sub, 'trunk_dir case selected')
      read (buff(eqpos+1:), *) tree_array(i_tree) % trunk_dir
      call mesg (sub, 'read trunk_dir =', tree_array(i_tree) % trunk_dir)

    case ('trunk_twist')

      if (DEBUG) call mesg (sub, 'trunk_twist case selected')
      read (buff(eqpos+1:), *) tree_array(i_tree) % trunk_twist
      call mesg (sub, 'read trunk_twist =', tree_array(i_tree) % trunk_twist)
      call mesg (sub, 'converting trunk_twist to radians')
      !--assume input in degrees, convert to radians
      tree_array(i_tree) % trunk_twist = (acos (-1._rp) / 180) *        &
                                         tree_array(i_tree) % trunk_twist

    case ('twist')

      if (DEBUG) call mesg (sub, 'twist case selected')
      read (buff(eqpos+1:), *) tree_array(i_tree) % twist
      call mesg (sub, 'read twist =', tree_array(i_tree) % twist)
      call mesg (sub, 'converting twist to radians')
      !--assume input in degrees, convert to radians
      tree_array(i_tree) % twist = (acos (-1._rp) / 180) *  &
                                   tree_array(i_tree) % twist
 
    case ('x0')

      if (DEBUG) call mesg (sub, 'x0 case selected')
      read (buff(eqpos+1:), *) tree_array(i_tree) % x0
      call mesg (sub, 'read x0 =', tree_array(i_tree) % x0)
    case default
      call error (sub, 'invalid variable tag ' // buff(1:eqpos-1) //  &
                  ' at line', line)
  end select
  
end do

close (lun)

!--issue error if missing some required data
if (.not. (have_n_tree .and. have_n_sub_branch .and.  &
           have_rel_dir .and. have_root_height) ) then
  call error (sub, 'missing required trees data')
end if

if (DEBUG) call exit_sub (sub)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine case_n_tree ()
  implicit none

  read (buff(eqpos + 1:), *) n_tree
  if (n_tree <= 0) then
    call error (sub, 'n_tree must be > 0 at line', line)
  end if
  call mesg (sub, 'read n_tree =', n_tree)

  have_n_tree = .true.

  !--now allocate the tree_array we need to read in
  !--do not allocate the trunks yet!
  allocate (tree_array(n_tree))

  end subroutine case_n_tree

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine case_tree ()
  implicit none

  call mesg (sub, 'entering tree block')
  
  !--check for beginning of tree block right after equals (since we
  !  removed whitespace), and rest of line is blank
  if ((buff(eqpos + 1 : eqpos + len (begin_tree_block)) /=           &
       begin_tree_block) .or.                                        &
      (verify (buff(eqpos + len (begin_tree_block) + 1:), ' ') /= 0)) then
    call error (sub, 'invalid tree block syntax at line', line)
  end if

  state = 1  !--we alter the state when in a tree block
  i_tree = i_tree + 1
  
  !--now read in tree data until we read a line with a end_tree_block

  end subroutine case_tree
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine read_trees_conf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! careful with the pts here: careful memory leaks
! may need to do subroutines instead
! we are assuming nd = 3 here
!--fills some things that are not read in config file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine init_tree (itree, tree)
implicit none

integer, intent (in) :: itree
type (tree_type), intent (inout) :: tree

character (*), parameter :: sub_name = mod_name // '.init_tree'

!logical, parameter :: DEBUG = .true. 

real (rp), parameter :: pi = 3.14159265359_rp  ! precision not needed
real (rp), parameter :: thresh = 100._rp * epsilon (1._rp)

real (rp) :: x_hat(nd), y_hat(nd)

!---------------------------------------------------------------------

if (DEBUG) call enter_sub (sub_name)

!--allocate the trunk
allocate (tree % trunk)

ident = ident + 1  
tree % trunk % ident = ident

tree % trunk % itree = itree

tree % trunk % gen = 0

tree % trunk % n_sub_branch = tree % n_sub_branch

if (tree % trunk % gen <= tree % max_res_gen) then
  tree % trunk % resolved = .true.
else
  tree % trunk % resolved = .false.
end if

!--somehow determine if defaults are set or not?
tree % trunk % l = tree % l
tree % trunk % d = tree % d

tree % trunk % x0 = tree % x0

tree % trunk % root_height = -1._rp  ! bogus, since parent not defined

!--set taper
tree % trunk % taper = tree % taper

tree % trunk % twist = tree_array(itree) % trunk_twist

! direction of trunk: absolute coordinates
tree % trunk % abs_dir = real (tree % trunk_dir, rp)
!--normalize
tree % trunk % abs_dir = (tree % trunk % abs_dir) /  &
                         mag (tree % trunk % abs_dir)

if (DEBUG) then
  call mesg (sub_name, 'tree % trunk_dir =', tree % trunk_dir)
  call mesg (sub_name, 'tree % trunk % abs_dir =', tree % trunk % abs_dir)
end if

! relative to absolute coords (not really needed for trunk)
! 8/11/04 note: not sure if this is valid 
tree % trunk % rel_dir = real (tree % trunk_dir, rp)
tree % trunk % rel_dir = (tree % trunk % rel_dir) /  &
                         mag (tree % trunk % rel_dir)

!--local trunk coordinate directions in abs coordinates
!--this may be sort of weird for trees not aligned in z-dir
!--these are NOT normalized yet: wait to check if degenerate
tree % trunk % z_hat = tree % trunk % abs_dir
tree % trunk % x_hat = cross_product (tree % trunk % z_hat,    &
                                      (/ 0._rp, 0._rp, 1._rp /))
tree % trunk % y_hat = cross_product (tree % trunk % z_hat,  &
                                      tree % trunk % x_hat)

!--coords degenerate: just use absolute instead
if (maxval ( abs (tree % trunk % y_hat - zero_vector) ) < thresh) then
  tree % trunk % x_hat = (/ 1._rp, 0._rp, 0._rp /)
  tree % trunk % y_hat = (/ 0._rp, 1._rp, 0._rp /)
  tree % trunk % z_hat = (/ 0._rp, 0._rp, 1._rp /)
end if

!--normalize
tree % trunk % x_hat = (tree % trunk % x_hat) / mag (tree % trunk % x_hat)
tree % trunk % y_hat = (tree % trunk % y_hat) / mag (tree % trunk % y_hat)
tree % trunk % z_hat = (tree % trunk % z_hat) / mag (tree % trunk % z_hat)

!--apply trunk twist
x_hat = cos (tree % trunk % twist) * (tree % trunk % x_hat) +  &
        sin (tree % trunk % twist) * (tree % trunk % y_hat)
y_hat = cos (tree % trunk % twist) * (tree % trunk % y_hat) -  &
        sin (tree % trunk % twist) * (tree % trunk % x_hat)
tree % trunk % x_hat = x_hat
tree % trunk % y_hat = y_hat
!--z_hat is untouched

!--in the cases where we do want the x_hat, y_hat, z_hat along cartesian
!  directions, but we have applied a rotation, we correct be making
!  epsilon changes
call cartesian_correction (tree % trunk % x_hat, tree % trunk % y_hat,  &
                           tree % trunk % z_hat)

!--normalize again due to errors in sin/cos? or cartesian correction
tree % trunk % x_hat = (tree % trunk % x_hat) / mag (tree % trunk % x_hat)
tree % trunk % y_hat = (tree % trunk % y_hat) / mag (tree % trunk % y_hat)
tree % trunk % z_hat = (tree % trunk % z_hat) / mag (tree % trunk % z_hat)

! note: we leave (tree % trunk % parent_branch) nullified

call heightwidth_bbox_br (tree % trunk)

!call Ap_bbox_br (tree % trunk) 

! use abs dir to determine zone
! actually, only smallest resolved and unresolved need zones
call set_zone_branch (tree % trunk)

! relative coordinates are not defined for trunk: it has no parent

! forces not needed yet

! construct sub branches: need rules for constructing branches
! for now, we are interested in self-similar trees, so we have
! to come up with some recursive way to do this.
if (tree % n_gen > 0) then
  call add_sub_branches (tree % trunk)
end if

!--option sanity check
if (tree % n_gen == tree % max_res_gen) then

  !if (use_unres_f_test) then
  !  call error (sub_name, 'cannot have use_unres_f_test = .true.' // n_l //  &
  !                        'there are no unresolved branches')
  !end if
  
else if (tree % n_gen < tree % max_res_gen) then

  call error (sub_name, 'n_gen < max_res_gen does not make sense')

end if

if (DEBUG) call exit_sub (sub_name)

end subroutine init_tree

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--for components of direction vectors close to 0, set exactly to 0
!--this routine DOES NOT normalize the corrected vectors
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cartesian_correction (x_hat, y_hat, z_hat)
implicit none

real (rp), intent (inout) :: x_hat(nd), y_hat(nd), z_hat(nd)

character (*), parameter :: sub_name = mod_name // '.cartesian_correction'

integer :: i

!---------------------------------------------------------------------

if (VERBOSE) call enter_sub (sub_name)

do i = 1, nd
  if ( abs (x_hat(i)) < epsilon (x_hat(i)) ) x_hat(i) = 0.0_rp
  if ( abs (y_hat(i)) < epsilon (y_hat(i)) ) y_hat(i) = 0.0_rp
  if ( abs (z_hat(i)) < epsilon (z_hat(i)) ) z_hat(i) = 0.0_rp
end do

if (VERBOSE) call exit_sub (sub_name)

end subroutine cartesian_correction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! adds a sub-branch to branch
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
recursive subroutine add_sub_branches (p)
implicit none

!type (branch_type), intent (inout) :: p  ! parent
type (branch_type), pointer :: p  ! parent

character (*), parameter :: sub_name = mod_name // '.add_sub_branches'

real (rp), parameter :: thresh = 100._rp * epsilon (1._rp)
real (rp), parameter :: pi = 3.14159265359_rp

integer :: i

real (rp) :: d
real (rp) :: dir_disp(nd)
real (rp) :: disp(nd)
real (rp) :: h_p
real (rp) :: x_hat(nd), y_hat(nd)

type (branch_type), pointer :: t  ! this

!---------------------------------------------------------------------

if (DEBUG) call enter_sub (sub_name)

if ((p % n_sub_branch) < 1) then
  if (DEBUG) call exit_sub (sub_name)
  return  ! do nothing
end if

allocate (p % sub_branch(p % n_sub_branch))

! careful: n_sub_branch is variable, cannot change while in this loop
do i = 1, p % n_sub_branch

  t => p % sub_branch(i)

  t % parent_branch => p  ! set pointer back to parent
  
  ident = ident + 1
  t % ident = ident

  t % itree = p % itree

  t % gen = (p % gen) + 1

  if (DEBUG) then
    if (t % gen > tree_array (t % itree) % n_gen) then
      call error (sub_name, 'invalid gen')
    end if
  end if

  if (DEBUG) then
    call mesg (sub_name, 'adding branch ', i ,         &
                         ' on level ', t % gen)
  end if

  t % n_sub_branch = tree_array (t % itree) % n_sub_branch

  if (t % gen <= tree_array(t % itree) % max_res_gen) then
    t % resolved = .true.
  else
    t % resolved = .false.
  end if
  
  t % l = tree_array (t % itree) % ratio * (p % l)
  t % d = tree_array (t % itree) % ratio * (p % d)

  ! areas set after abs_dir

  t % root_height = tree_array(t % itree) % root_height(i)

  t % taper = tree_array(t % itree) % taper

  t % twist = tree_array(t % itree) % twist(i)

  !--calculation of x0 moved after that of abs_dir, to facilitate
  !  putting sub-branches on the outside of the parent
  !t % x0 = (p % x0) + (t % root_height) * (p % l) * (p % abs_dir)

  !if (DEBUG) call mesg (sub_name, 't % x0(3) = ', t % x0(3))

  ! direction relative to parents coords
  t % rel_dir(:) = tree_array(t % itree) % rel_dir(:, i)

  t % abs_dir(1) = (t % rel_dir(1)) * (p % x_hat(1)) +  &
                   (t % rel_dir(2)) * (p % y_hat(1)) +  &
                   (t % rel_dir(3)) * (p % z_hat(1))

  t % abs_dir(2) = (t % rel_dir(1)) * (p % x_hat(2)) +  &
                   (t % rel_dir(2)) * (p % y_hat(2)) +  &
                   (t % rel_dir(3)) * (p % z_hat(2))

  t % abs_dir(3) = (t % rel_dir(1)) * (p % x_hat(3)) +  &
                   (t % rel_dir(2)) * (p % y_hat(3)) +  &
                   (t % rel_dir(3)) * (p % z_hat(3))

  !--calculation of x0
  h_p = (t % root_height) * (p % l)
  t % x0 = (p % x0) + h_p * (p % abs_dir)

  if (sub_branches_outside) then
    !--assumes we are adding cap of side d/2 to top, otherwise
    !  the d/2 part is not needed at the top
    !--preserves root height (when its < 1) when moving to outside since
    !  |cross_product| = sin(ang) for unit vectors
    !--when root height = 1, root height is NOT preserved, and branch is
    !  based d/2 in direction of branch
    d = (p % d) * (1._rp - (p % taper) * h_p / (p % l))

    !--here, we basically ignore squared edges, and always displace by d/2
    if (abs (t % root_height - 1._rp) < epsilon (0._rp)) then

      if (add_cap) then
        disp = 0.5_rp * d * (t % abs_dir)
      else
        disp = 0.0_rp
      end if

    else

      dir_disp = (t % abs_dir) -                                      &
                 (p % abs_dir) * dot_product (t % abs_dir, p % abs_dir)
      if ( mag ( dir_disp ) > epsilon (0.0_rp) ) then
        dir_disp = dir_disp / mag ( dir_disp )
          !--comment out for compatibility with older versions of code
      else
        dir_disp = 0.0_rp
      end if
      disp = 0.5_rp * d * dir_disp

    end if
    
    t % x0 = (t % x0) + disp

  end if

  if (DEBUG) call mesg (sub_name, 't % x0(3) = ', t % x0(3))

  ! now we have the absolute direction of this sub-branch,
  ! we can calculate the x_hat, y_hat, z_hat of this sub-branch
  ! i.e. local sub-branch coordinate unit vectors, rel. to abs frame
  ! this is a particular convention to define the local coord. system
  ! these need to be unit vectors

  ! original convention
  !t % z_hat = t % abs_dir
  !t % y_hat = cross_product (t % z_hat, p % abs_dir)
  !t % x_hat = cross_product (t % y_hat, t % z_hat)

  ! convention that keeps y unit vector in absolute y-z plane
  ! for grid-bound trees (I think)
  ! mainly useful for 2-d cruciform trees
  t % z_hat = t % abs_dir
  t % x_hat = cross_product (t % z_hat, p % abs_dir)
  t % y_hat = cross_product (t % z_hat, t % x_hat)

  ! if this coordinate system is degenerate (e.g. when p, t parallel)
  ! then t inherits p local coordinate system
  if (maxval (abs (t % y_hat - zero_vector)) < thresh) then

    if (DEBUG) then
      call mesg (sub_name, "t's coords degenerate, using p's")
    end if

    t % x_hat = p % x_hat
    t % y_hat = p % y_hat
    t % z_hat = p % z_hat

  else  !--normalize

    t % x_hat = t % x_hat / mag (t % x_hat)
    t % y_hat = t % y_hat / mag (t % y_hat)
    t % z_hat = t % z_hat / mag (t % z_hat)

  end if

  !--apply twist
  x_hat = cos (t % twist) * (t % x_hat) + sin (t % twist) * (t % y_hat)
  y_hat = cos (t % twist) * (t % y_hat) - sin (t % twist) * (t % x_hat)
  t % x_hat = x_hat
  t % y_hat = y_hat
  
  call cartesian_correction (t % x_hat, t % y_hat, t % z_hat)

  !--normalize again (errors in twist/correction)
  t % x_hat = t % x_hat / mag (t % x_hat)
  t % y_hat = t % y_hat / mag (t % y_hat)
  t % z_hat = t % z_hat / mag (t % z_hat)

  if (DEBUG) then
    call mesg (sub_name, 't % x_hat = ', t % x_hat)
    call mesg (sub_name, 't % y_hat = ', t % y_hat)
    call mesg (sub_name, 't % z_hat = ', t % z_hat)
  end if

  !--evaluate bounding box height, width
  call heightwidth_bbox_br (t)

  if (DEBUG) then
    call mesg (sub_name, 't % ident = ', t % ident)
    call mesg (sub_name, 't % height_bbox = ', t % height_bbox)
    call mesg (sub_name, 't % width_bbox = ', t % width_bbox)
  end if

  !--evaluate bounding box projected area
  !call Ap_bbox_br (t)

  !--now we have the abs_dir, can set the zone
  call set_zone_branch (t)

  ! check if any more sub-branches
  if (t % gen < tree_array (t % itree) % n_gen) then
    ! recursion
    call add_sub_branches (t)
  end if

end do

if (DEBUG) call exit_sub (sub_name)

end subroutine add_sub_branches

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine heightwidth_bbox_br (br)
use param, only : L_x, L_y, L_z, nproc  !--added to detect bbox out of domain
implicit none

type (branch_type), intent (inout) :: br

character (*), parameter :: sub_name = mod_name // '.heightwidth_bbox_br'

logical, parameter :: Y_experiment = .false.
logical, parameter :: warn_only = .true.

real (rp) :: r, l, d
real (rp) :: height, width

!---------------------------------------------------------------------

r = tree_array(br % itree) % ratio
l = br % l
d = br % d
!l = tree_array(br % itree) % l
!d = tree_array(br % itree) % d

if (Y_experiment) then  !--for right-angled Y-shapes

  call warn (sub_name, 'Y_experiment is active')
  
  height = l + (d / 2._rp + r * l) * (1._rp + sqrt (2._rp) * r) /  &
               (1._rp - r**2) / sqrt (2._rp)
  width = sqrt (2._rp) * (d / 2._rp + r * l) *      &
          (1._rp + sqrt (2._rp) * r) / (1._rp - r**2)


else  !--default conservative bounding box (branches perpendicular)

  height = l
  if (add_cap) height = height + 0.5_rp * d
  height = height / (1._rp - r)

  width = d + 2._rp * r * height  !--the 1/1-r is already in height
                                  !--add_cap case: already in height
                                  !--add option for sub_branches_outside?
end if

br % height_bbox = height
br % width_bbox = width

!--assumes trunk is upright (aligned with z)
if (height > nproc * L_z) then
    if ( warn_only ) then
        call warn (sub_name, 'bbox height > nproc * L_z')
    else
        call error (sub_name, 'bbox height > nproc * L_z')
    end if
else if (width > min (L_x, L_y)) then
    if ( warn_only ) then
        call warn (sub_name, 'bbox width > L_x or L_y')
    else
        call error (sub_name, 'bbox width > L_x or L_y')
    end if
end if

end subroutine heightwidth_bbox_br

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!--calculates the projected area of bbox as seen from flow_dir
!!--br % height_bbox, br % width_bbox must already be set
!!--br % x_hat, br % y_hat, br % z_hat must already be set
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine Ap_bbox_br (br)
!implicit none
!
!type (branch_type), intent (inout) :: br
!
!character (*), parameter :: sub_name = mod_name // '.Ap_bbox_br'
!
!!---------------------------------------------------------------------
!
!br % Ap_bbox = (br % width_bbox) * (br % height_bbox) *      &
!               (abs (dot_product (flow_dir, br % x_hat)) +   &
!                abs (dot_product (flow_dir, br % y_hat))) +  &
!               (br % width_bbox) * (br % width_bbox) *       &
!               abs (dot_product (flow_dir, br % z_hat))
!
!end subroutine Ap_bbox_br

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_zone_branch (branch)
implicit none

type (branch_type), intent (inout) :: branch

character (*), parameter :: sub_name = mod_name // '.set_zone_branch'

character (*), parameter :: zone_selection = 'single'
                            !--selects how zones defined
                            !  'b.f', 'b.f|b.p', 'absolute', 'single'
                            !  'hv', 'b.f+b.p'

integer, parameter :: zone_bogus = -1  ! bogus init value
integer, parameter :: e1(nd) = (/ 1, 0, 0 /)
integer, parameter :: e2(nd) = (/ 0, 1, 0 /)
integer, parameter :: e3(nd) = (/ 0, 0, 1 /)

integer :: int_dir(nd)
integer :: int_b_dot_f
integer :: int_b_dot_p

!----------------------------------------------------------------------

if (DEBUG) call enter_sub (sub_name)

if (nzone <= 0) then
  write (msg, *) 'invalid nzone = ', nzone
  call error (sub_name, msg)
end if

!--may want to change zone selection methods

! for now, set zone depending on abs_dir
int_dir = nint (branch % abs_dir)

! this can be -1, 0, +1
int_b_dot_f = nint (dot_product (branch % abs_dir, flow_dir))

select case (zone_selection)  ! plan on making internal subroutines

  case ('absolute')

    if (nzone /= 4) then
      write (msg, *) 'nzone = ', nzone, n_l,           &
                     'expecting 4 with ', zone_selection
      call error (sub_name, msg)
    end if

    ! can alias (sort of equivalence) zones here
    if (all (int_dir == e1)) then
      branch % zone = 1
    else if (all (int_dir == -e1)) then
      branch % zone = 2
    else if (all (int_dir == e2)) then
      branch % zone = 3
    else if (all (int_dir == -e2)) then
      branch % zone = 3
    else if (all (int_dir == e3)) then
      branch % zone = 4
    else if (all (int_dir == -e3)) then
      branch % zone = 4
    else
      write (msg, *) 'not sure which zone for ind_dir = ', int_dir
      call error (sub_name, msg)
    end if

  case ('b.f')

    if (nzone /= 3) then
      write (msg, *) 'nzone = ', nzone, n_l,           &
                     'expecting 3 with ', zone_selection
      call error (sub_name, msg)
    end if

    select case (int_b_dot_f)
      case (-1)
        branch % zone = 1
      case (0)
        branch % zone = 2
      case (1)
        branch % zone = 3
      case default
        write (msg, *) 'not which zone for int_b_dot_f = ', int_b_dot_f
        call error (sub_name, msg)
    end select

  case ('b.f+b.p')

    if (nzone /= 4) then
      write (msg, *) 'nzone = ', nzone, n_l,           &
                     'expecting 4 with ', zone_selection
      call error (sub_name, msg)
    end if

    if (associated (branch % parent_branch)) then
      !--for 3d-cross geometry, this can only be 0, 1
      int_b_dot_p = nint (dot_product (branch % abs_dir,                &
                                       branch % parent_branch % abs_dir))
    else

      if (branch % gen == 0) then  !--trunk
        int_b_dot_p = 1  !--treat as if straight out of parent
                         !--this is mainly to avoid errors caused by
                         !  assigning a bogus zone here
      else
        call error (sub_name, 'expecting associated parent_branch')
      end if

    end if
    
    select case (int_b_dot_f)
      case (0)
      
        if (int_b_dot_p == 0) then
          branch % zone = 1  !--cross-stream/bent
        else
          branch % zone = 2  !--cross-stream/straight
        end if

      case (-1)
        branch % zone = 3  !--upstream
      case (1)
        branch % zone = 4  !--downstream
      case default
        write (msg, *) 'not which zone for int_b_dot_f = ', int_b_dot_f
        call error (sub_name, msg)
    end select

  case ('b.f|b.p')

    if (nzone /= 6) then
      write (msg, *) 'nzone = ', nzone, n_l,           &
                     'expecting 6 with ', zone_selection
      call error (sub_name, msg)
    end if

    if (associated (branch % parent_branch)) then

      ! for current geometries, this can only be 0, +1
      int_b_dot_p = nint (dot_product (branch % abs_dir,                &
                                       branch % parent_branch % abs_dir))

      if (int_b_dot_p == 0) then

        if (int_b_dot_f == -1) then
          branch % zone = 1
        else if (int_b_dot_f == 0) then
          branch % zone = 2
        else if (int_b_dot_f == 1) then
          branch % zone = 3
        else
          write (msg, *) 'unexpected value of int_b_dot_f = ', int_b_dot_f
          call error (sub_name, msg)
        end if

      else if (int_b_dot_p == 1) then

        if (int_b_dot_f == -1) then
          branch % zone = 4
        else if (int_b_dot_f == 0) then
          branch % zone = 5
        else if (int_b_dot_f == 1) then
          branch % zone = 6
        else
          write (msg, *) 'unexpected value of int_b_dot_f = ', int_b_dot_f
          call error (sub_name, msg)
        end if

      else

        write (msg, *) 'unexpected value of int_b_dot_p = ', int_b_dot_p
        call error (sub_name, msg)

      end if

    else  ! this should be the trunk, which has no parent

      if (branch % gen == 0) then
        branch % zone = zone_bogus  ! for now
        ! NOTE: in future, make it so that only branches that need zones
        !       actually have valid ones (easier to catch mistakes?)
      else
        call error (sub_name, 'expecting gen 0')
      end if
    
    end if

  case ('hv')

    if (nzone /= 2) then
      write (msg, *) 'nzone = ', nzone, n_l,           &
                     'expecting 2 with ', zone_selection
      call error (sub_name, msg)
    end if

    if (int_dir(3) == 0) then
      branch % zone = 1
    else
      branch % zone = 2
    end if

  case ('single')

    if (nzone /= 1) then
      write (msg, *) 'nzone = ', nzone, n_l,           &
                     'expecting 1 with ', zone_selection
      call error (sub_name, msg)
    end if

    branch % zone = 1

  case default

    write (msg, *) 'invalid zone_selection = ', zone_selection
    call error (sub_name, msg)

end select
  
if (DEBUG) call exit_sub (sub_name)

end subroutine set_zone_branch
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module trees_setup_ls
