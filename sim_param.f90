module sim_param
use types, only : rprec
use param, only : ld, ny, nz, lbz
implicit none

save
public

logical :: sim_param_initialized = .false.

!--still testing allocatable array implementation:
!  ifc 7.1 segfaults
!  ifort 8.1 ok
!  xlf segfaults when in MPI mode 256^3/32 cpu (need to test other combos)
    
real (rprec), dimension (:, :, :), allocatable :: u, v, w
real (rprec), dimension (:, :, :), allocatable :: dudx, dudy, dudz,  &
                                                  dvdx, dvdy, dvdz,  &
                                                  dwdx, dwdy, dwdz,  &
                                                  RHSx, RHSy, RHSz,  &
                                                  RHSx_f, RHSy_f, RHSz_f

real (rprec), dimension (:, :, :), allocatable :: dpdx, dpdy, dpdz

real (rprec), dimension (:, :, :), allocatable :: txx, txy, tyy
real (rprec), dimension (:, :, :), allocatable :: txz, tyz, tzz

real (rprec), target, dimension (:, :, :), allocatable :: p

real (rprec), dimension (:, :, :), allocatable :: divtx, divty, divtz

real (rprec), dimension (:, :, :), allocatable :: fx, fy, fz, &
                                                  fxa, fya, fza

real (rprec), dimension (:, :, :), allocatable :: theta, q
    !--Added for scalars


contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine sim_param_init ( array_list_opt )
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
!--this is needed to make some post-processing code compilable, since
!  is only allocates the space to be used
!--initialized all data to zero
!
implicit none

character (*), intent (in), optional :: array_list_opt
    !--comma separated list of arrays to init

integer, parameter :: narray_max = 128
integer, parameter :: narray_name_len = 8  !--can increase if needed

character (*), parameter :: array_list_def = 'u, v, w,' //                 &
                                             'dudx, dudy, dudz,' //        &
                                             'dvdx, dvdy, dvdz,' //        &
                                             'dwdx, dwdy, dwdz,' //        &
                                             'RHSx, RHSy, RHSz,' //        &
                                             'RHSx_f, RHSy_f, RHSz_f,' //  &
                                             'dpdx, dpdy, dpdz,' //        &
                                             'txx, txy, tyy,' //           &
                                             'txz, tyz, tzz,' //           &
                                             'p,' //                       &
                                             'divtx, divty, divtz,' //     &
                                             'fx, fy, fz,' //              &
                                             'fxa, fya, fza,' //           &
                                             'theta, q'

$if ($DEBUG)
logical, parameter :: DEBUG = .true.
$endif

character (narray_name_len * narray_max) :: array_list
character (narray_name_len) :: array(narray_max)
character (narray_name_len) :: alloced_array(narray_max)


integer :: i
integer :: ios

!---------------------------------------------------------------------


!if ( .not. present ( arrays_to_init ) ) then
!    !--initialize all arrays
!    allocate ( u(ld, ny, lbz:nz),  &
!               v(ld, ny, lbz:nz),  &
!               w(ld, ny, lbz:nz) )
!    
!    u = 0.0_rprec
!    v = 0.0_rprec
!    w = 0.0_rprec
!    
!    allocate( dudx(ld, ny, lbz:nz),  &
!              dudy(ld, ny, lbz:nz),  &
!              dudz(ld, ny, lbz:nz),  &
!              dvdx(ld, ny, lbz:nz),  &
!              dvdy(ld, ny, lbz:nz),  &
!              dvdz(ld, ny, lbz:nz),  &
!              dwdx(ld, ny, lbz:nz),  &
!              dwdy(ld, ny, lbz:nz),  &
!              dwdz(ld, ny, lbz:nz),  &
!              RHSx(ld, ny, lbz:nz),  &
!              RHSy(ld, ny, lbz:nz),  &
!              RHSz(ld, ny, lbz:nz),  &
!              RHSx_f(ld, ny, lbz:nz),  &
!              RHSy_f(ld, ny, lbz:nz),  &
!              RHSz_f(ld, ny, lbz:nz) )
!    
!    dudx = 0.0_rprec
!    dudy = 0.0_rprec
!    dudz = 0.0_rprec
!    dvdx = 0.0_rprec
!    dvdy = 0.0_rprec
!    dvdz = 0.0_rprec
!    dwdx = 0.0_rprec
!    dwdy = 0.0_rprec
!    dwdz = 0.0_rprec
!    RHSx = 0.0_rprec
!    RHSy = 0.0_rprec
!    RHSz = 0.0_rprec
!    RHSx_f = 0.0_rprec
!    RHSy_f = 0.0_rprec
!    RHSz_f = 0.0_rprec
!    
!    allocate ( dpdx(ld, ny, nz),  &
!               dpdy(ld, ny, nz),  &
!               dpdz(ld, ny, nz) )
!    
!    dpdx = 0.0_rprec
!    dpdy = 0.0_rprec
!    dpdz = 0.0_rprec
!
!    allocate ( txx(ld, ny, lbz:nz),  &
!               txy(ld, ny, lbz:nz),  &
!               tyy(ld, ny, lbz:nz),  &
!               txz(ld, ny, lbz:nz),  &
!               tyz(ld, ny, lbz:nz),  &
!               tzz(ld, ny, lbz:nz) )
!
!    txx = 0.0_rprec
!    txy = 0.0_rprec
!    tyy = 0.0_rprec
!    txz = 0.0_rprec
!    tyz = 0.0_rprec
!    tzz = 0.0_rprec
!    
!    allocate ( p(ld, ny, 0:nz) )
!    
!    p = 0.0_rprec
!    
!    allocate ( divtx(ld, ny, lbz:nz),  &
!               divty(ld, ny, lbz:nz),  &
!               divtz(ld, ny, lbz:nz) )
!    
!    divtx = 0.0_rprec
!    divty = 0.0_rprec
!    divtz = 0.0_rprec
!    
!    allocate ( theta(ld, ny, nz),  &
!               q(ld, ny, nz) )
!    
!    theta = 0.0_rprec
!    q = 0.0_rprec
!
!else

array = ''
alloced_array = ''

if ( present ( array_list_opt ) ) then
    array_list = array_list_opt
else  !--use default list: all arrays allocated
    array_list = array_list_def
end if

read ( array_list, *, iostat=ios ) array
!--comma or space separated

!--ios < 0 indicates array_list is too small  (or too large) to fit
!  into array, which for now is OK (at least when its too small)
if ( ios > 0 ) then
    write (*, *) 'sim_param_init: error reading array list'
    stop
end if

do i = 1, size ( array )
    
    if ( trim (array(i)) == '' ) exit
        !--assume rest of array(i:) is blank as well

    select case ( trim ( array(i) ) )
    case ( 'u' )
        allocate ( u(ld, ny, lbz:nz) )
        u = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'v' )
        allocate ( v(ld, ny, lbz:nz) )
        v = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'w' )
        allocate ( w(ld, ny, lbz:nz) )
        w = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'dudx' )
        allocate( dudx(ld, ny, lbz:nz) )
        dudx = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'dudy' )
        allocate( dudy(ld, ny, lbz:nz) )
        dudy = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'dudz' )
        allocate( dudz(ld, ny, lbz:nz) )
        dudz = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'dvdx' )
        allocate( dvdx(ld, ny, lbz:nz) )
        dvdx = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'dvdy' )
        allocate( dvdy(ld, ny, lbz:nz) )
        dvdy = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'dvdz' )
        allocate( dvdz(ld, ny, lbz:nz) )
        dvdz = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'dwdx' )
        allocate( dwdx(ld, ny, lbz:nz) )
        dwdx = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'dwdy' )
        allocate( dwdy(ld, ny, lbz:nz) )
        dwdy = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'dwdz' )
        allocate( dwdz(ld, ny, lbz:nz) )
        dwdz = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'RHSx' )
        allocate( RHSx(ld, ny, lbz:nz) )
        RHSx = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'RHSy' )
        allocate( RHSy(ld, ny, lbz:nz) )
        RHSy = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'RHSz' )
        allocate( RHSz(ld, ny, lbz:nz) )
        RHSz = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'RHSx_f' )
        allocate( RHSx_f(ld, ny, lbz:nz) )
        RHSx_f = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'RHSy_f' )
        allocate( RHSy_f(ld, ny, lbz:nz) )
        RHSy_f = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'RHSz_f' )
        allocate( RHSz_f(ld, ny, lbz:nz) )
        RHSz_f = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'dpdx' )
        allocate ( dpdx(ld, ny, nz) )
        dpdx = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'dpdy' )
        allocate ( dpdy(ld, ny, nz) )
        dpdy = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'dpdz' )
        allocate ( dpdz(ld, ny, nz) )
        dpdz = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'txx' )
        allocate ( txx(ld, ny, lbz:nz) )
        txx = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'txy' )
        allocate ( txy(ld, ny, lbz:nz) )
        txy = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'tyy' )
        allocate ( tyy(ld, ny, lbz:nz) )
        tyy = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'txz' )
        allocate ( txz(ld, ny, lbz:nz) )
        txz = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'tyz' )
        allocate ( tyz(ld, ny, lbz:nz) )
        tyz = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'tzz' )
        allocate ( tzz(ld, ny, lbz:nz) )
        tzz = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'p' )
        allocate ( p(ld, ny, 0:nz) )
        p = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'divtx' )
        allocate ( divtx(ld, ny, lbz:nz) )
        divtx = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'divty' )
        allocate ( divty(ld, ny, lbz:nz) )
        divty = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'divtz' )
        allocate ( divtz(ld, ny, lbz:nz) )
        divtz = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'fx' )
        allocate ( fx(ld, ny, nz) )
        fx = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'fy' )
        allocate ( fy(ld, ny, nz) )
        fy = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )  
    case ( 'fz' )
        allocate ( fz(ld, ny, nz) )
        fz = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'fxa' )
        allocate ( fxa(ld, ny, nz) )
        fxa = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'fya' )
        allocate ( fya(ld, ny, nz) )
        fya = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )  
    case ( 'fza' )
        allocate ( fza(ld, ny, nz) )
        fza = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'theta' )
        allocate ( theta(ld, ny, nz) )
        theta = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case ( 'q' )
        allocate ( q(ld, ny, nz) )
        q = 0.0_rprec
        write ( alloced_array(i), '(a)' ) trim ( array(i) )
    case default
        write (*, *) 'sim_param_init: invalid array ' // trim ( array(i) )
        stop
    end select

end do

write(*,*) 'Setting sim_param_initialized true'
sim_param_initialized = .true.
    !--beware this does NOT guarentee that a particular set of arrays
    !  was initialized!

$if ($DEBUG)
if ( DEBUG ) then

    write (*, *) 'sim_param_init: the following arrays were initialized'

    do i = 1, size (alloced_array)
        if ( trim ( alloced_array(i) ) == '' ) cycle
        write ( *, '(3a)' ) 'array=<', trim (alloced_array(i)), '>'
    end do
    
end if
$endif

end subroutine sim_param_init
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module sim_param
