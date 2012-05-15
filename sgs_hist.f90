module sgs_hist

use types, only: rprec
use param
use stat_defs, only: HISTcs2_t, HISTtn_t, HISTnu_t, HISTee_t

implicit none

save
private
public sgs_hist_init, sgs_hist_update_vals, sgs_hist_finalize

! Used for reading and writing data to file (for cumulative runs)
    character(*), parameter :: fname_cs2_base = path // 'output/hist_cs2.z-'
    character(*), parameter :: fname_tn_base = path // 'output/hist_tn.z-'
    character(*), parameter :: fname_nu_base = path // 'output/hist_nu.z-'
    character(*), parameter :: fname_ee_base = path // 'output/hist_ee.z-'
    character(*), parameter :: fname_histp_base = path // 'hist_param.out'
    character(*), parameter :: fname_histv_base = path // 'hist_vals.out'
    character(64) :: cl, fname_cs2, fname_tn, fname_nu, fname_ee
    character(64) :: fname_histp, fname_histv

!*********************************************************************
contains
!*********************************************************************
!---------------------------------------------------------------------
subroutine sgs_hist_init()
!  Allocate, determine bin values, and initialize the counts to zero
!  If sgs_hist_cumulative is true, then these values will be read 
!    from hist_{param,vals}.out rather than lesgo.conf.

use grid_defs, only : grid_t
use functions, only : cell_indx
use string_util, only : string_concat
implicit none

integer :: k, n
logical :: exstp = .false.
logical :: exstv = .false.
integer :: nproc_test, nz_tot_test
real(rprec) :: L_z_test

! If sgs_hist_cumulative is true, try to read data from files
  if (sgs_hist_cumulative) then
    ! Check if files exist
      ! Build file name      
        fname_histp = ''; call string_concat( fname_histp, fname_histp_base // '.c', coord )
        fname_histv = ''; call string_concat( fname_histv, fname_histv_base // '.c', coord )

      ! Check for existence
        inquire (file=fname_histp, exist=exstp)
        inquire (file=fname_histv, exist=exstv)

    if ( exstp .and. exstv) then
      ! Read data from file (overwrites values from lesgo.conf or param.f90)
        write(*,*) 'Reading from files:', fname_histp, fname_histv

        open(unit=1,file=fname_histp,position='rewind',action='read',form='formatted')
          read(1,*) nproc_test, nz_tot_test, L_z_test
          if ( (nproc_test.ne.nproc) .or. (nz_tot_test.ne.nz_tot) .or. (L_z_test.ne.L_z) ) then
            write(*,*) 'Sgs-hist cumulative time error: nproc, nz_tot, and L_z must match'
            write(*,*) 'Starting sgs-histograms from scratch'
            exstp = .false.
          else
            read(1,*) sgs_hist_nloc
            deallocate ( sgs_hist_loc )
            allocate ( sgs_hist_loc(sgs_hist_nloc) )
            read(1,*) sgs_hist_loc
            read(1,*) cs2_bmin, cs2_bmax, cs2_nbins
            read(1,*) tn_bmin, tn_bmax, tn_nbins
            read(1,*) nu_bmin, nu_bmax, nu_nbins
            read(1,*) ee_bmin, ee_bmax, ee_nbins
          endif
        close(1)

    else
      ! Start from scratch
        if (.not.exstp) write(*,*) 'File not found:', fname_histp
        if (.not.exstv) write(*,*) 'File not found:', fname_histv
        write(*,*) 'Starting sgs-histograms from scratch'

    endif

  else
    write(*,*) 'Starting sgs-histograms from scratch'

  endif


! Allocate space to define plane locations - done for all procs
    allocate ( HISTcs2_t % istart (sgs_hist_nloc) )
    allocate ( HISTtn_t % istart (sgs_hist_nloc) )
    allocate ( HISTnu_t % istart (sgs_hist_nloc) )
    allocate ( HISTee_t % istart (sgs_hist_nloc) )

    allocate ( HISTcs2_t % coord (sgs_hist_nloc) )
    allocate ( HISTtn_t % coord (sgs_hist_nloc) )
    allocate ( HISTnu_t % coord (sgs_hist_nloc) )
    allocate ( HISTee_t % coord (sgs_hist_nloc) )

    allocate ( HISTcs2_t % ldiff (sgs_hist_nloc) )
    allocate ( HISTtn_t % ldiff (sgs_hist_nloc) )
    allocate ( HISTnu_t % ldiff (sgs_hist_nloc) )
    allocate ( HISTee_t % ldiff (sgs_hist_nloc) )

    allocate ( HISTcs2_t % hist_t (sgs_hist_nloc) )
    allocate ( HISTtn_t % hist_t (sgs_hist_nloc) )
    allocate ( HISTnu_t % hist_t (sgs_hist_nloc) )
    allocate ( HISTee_t % hist_t (sgs_hist_nloc) )


! Determine which coords have each location
    ! Initialize 
    HISTcs2_t % istart(:) = -1
    HISTcs2_t % coord(:) = -1
    HISTcs2_t % ldiff(:) = 0._rprec

    HISTtn_t % istart(:) = -1
    HISTtn_t % coord(:) = -1
    HISTtn_t % ldiff(:) = 0._rprec

    HISTnu_t % istart(:) = -1
    HISTnu_t % coord(:) = -1
    HISTnu_t % ldiff(:) = 0._rprec

    HISTee_t % istart(:) = -1
    HISTee_t % coord(:) = -1
    HISTee_t % ldiff(:) = 0._rprec

    ! For each zplane, find the coord where it belongs and set
    !   variables in that coord only   
    do k=1,sgs_hist_nloc

        ! If in this processor
        if(sgs_hist_loc(k) >= grid_t%z(1) .and. sgs_hist_loc(k) < grid_t%z(nz)) then

            ! Determine istart and ldiff (and set coord for future reference)
            HISTcs2_t % coord(k) = coord    ! should be -1 (default) if not using MPI
            HISTcs2_t % istart(k) = cell_indx('k',dz,sgs_hist_loc(k))
            HISTcs2_t % ldiff(k) = sgs_hist_loc(k) - grid_t%z(HISTcs2_t % istart(k))

            HISTtn_t % coord(k) = coord
            HISTtn_t % istart(k) = cell_indx('k',dz,sgs_hist_loc(k))
            HISTtn_t % ldiff(k) = sgs_hist_loc(k) - grid_t%z(HISTtn_t % istart(k))

            HISTnu_t % coord(k) = coord
            HISTnu_t % istart(k) = cell_indx('k',dz,sgs_hist_loc(k))
            HISTnu_t % ldiff(k) = sgs_hist_loc(k) - grid_t%z(HISTnu_t % istart(k))

            HISTee_t % coord(k) = coord
            HISTee_t % istart(k) = cell_indx('k',dz,sgs_hist_loc(k))
            HISTee_t % ldiff(k) = sgs_hist_loc(k) - grid_t%z(HISTee_t % istart(k))


            ! ALLOCATE
            ! For reference: bin_index given by:
            !   bin_index = min( ceil( max(value-bmin,-0.5_rprec) /db ), nbins+1 )
            ! Also:
            !   0 catches values below bmin
            !   nbins+1 catches values above bmax

            allocate ( HISTcs2_t % hist_t(k) % bins (0:cs2_nbins+1) )
            allocate ( HISTcs2_t % hist_t(k) % vals (0:cs2_nbins+1) )

            allocate ( HISTtn_t % hist_t(k) % bins (0:tn_nbins+1) )
            allocate ( HISTtn_t % hist_t(k) % vals (0:tn_nbins+1) )

            allocate ( HISTnu_t % hist_t(k) % bins (0:nu_nbins+1) )
            allocate ( HISTnu_t % hist_t(k) % vals (0:nu_nbins+1) )

            allocate ( HISTee_t % hist_t(k) % bins (0:ee_nbins+1) )
            allocate ( HISTee_t % hist_t(k) % vals (0:ee_nbins+1) )


            ! Define bins range for each variable
            !   For now this assumes all levels have the same range

                ! Number of bins
                HISTcs2_t % hist_t(k) % nbins = cs2_nbins
                HISTtn_t  % hist_t(k) % nbins = tn_nbins
                HISTnu_t  % hist_t(k) % nbins = nu_nbins
                HISTee_t  % hist_t(k) % nbins = ee_nbins

                ! Bin minimum
                HISTcs2_t % hist_t(k) % bmin = cs2_bmin
                HISTtn_t  % hist_t(k) % bmin = tn_bmin
                HISTnu_t  % hist_t(k) % bmin = nu_bmin
                HISTee_t  % hist_t(k) % bmin = ee_bmin

                ! Bin maximum
                HISTcs2_t % hist_t(k) % bmax = cs2_bmax
                HISTtn_t  % hist_t(k) % bmax = tn_bmax
                HISTnu_t  % hist_t(k) % bmax = nu_bmax
                HISTee_t  % hist_t(k) % bmax = ee_bmax

                ! Set spacing between bins
                HISTcs2_t % hist_t(k) % db = (cs2_bmax - cs2_bmin) / cs2_nbins
                HISTtn_t  % hist_t(k) % db = (tn_bmax - tn_bmin) / tn_nbins
                HISTnu_t  % hist_t(k) % db = (nu_bmax - nu_bmin) / nu_nbins
                HISTee_t  % hist_t(k) % db = (ee_bmax - ee_bmin) / ee_nbins

            ! Set bins variable to center of each bin range
                forall(n=0:cs2_nbins+1) HISTcs2_t % hist_t(k) % bins(n) = & 
                                        cs2_bmin + (n-0.5_rprec) * HISTcs2_t % hist_t(k) % db

                forall(n=0:tn_nbins+1) HISTtn_t % hist_t(k) % bins(n) = & 
                                       tn_bmin + (n-0.5_rprec) * HISTtn_t % hist_t(k) % db

                forall(n=0:nu_nbins+1) HISTnu_t % hist_t(k) % bins(n) = & 
                                       nu_bmin + (n-0.5_rprec) * HISTnu_t % hist_t(k) % db

                forall(n=0:ee_nbins+1) HISTee_t % hist_t(k) % bins(n) = & 
                                       ee_bmin + (n-0.5_rprec) * HISTee_t % hist_t(k) % db

            ! Initialize vals to zero
                HISTcs2_t % hist_t(k) % vals = 0.0_rprec
                HISTtn_t  % hist_t(k) % vals = 0.0_rprec
                HISTnu_t  % hist_t(k) % vals = 0.0_rprec
                HISTee_t  % hist_t(k) % vals = 0.0_rprec

        endif

    enddo
 
! If reading from file, do so now for bins and vals
  if ( exstp .and. exstv ) then

    open(unit=1,file=fname_histv,position='rewind',action='read',form='unformatted')
      do k=1,sgs_hist_nloc
        read(1) HISTcs2_t%hist_t(k)%bins, HISTcs2_t%hist_t(k)%vals
        read(1) HISTtn_t%hist_t(k)%bins, HISTtn_t%hist_t(k)%vals
        read(1) HISTnu_t%hist_t(k)%bins, HISTnu_t%hist_t(k)%vals
        read(1) HISTee_t%hist_t(k)%bins, HISTee_t%hist_t(k)%vals
      enddo
   close(1)

  endif

end subroutine sgs_hist_init
!---------------------------------------------------------------------

!---------------------------------------------------------------------
subroutine sgs_hist_update_vals()

use stat_defs, only: hist_binit
use functions, only: linear_interp
use sgs_param, only: Cs_opt2, Tn_all, Nu_t, ee_now
$if ($LVLSET)
use level_set_base, only: phi
$endif

implicit none

integer :: i, j, k
$if ($LVLSET)
real(rprec), dimension(nx,ny) :: phi_plane
$endif
real(rprec), dimension(nx,ny) :: cs2_plane, tn_plane, nu_plane, ee_plane

! For each z-plane location
    do k=1,sgs_hist_nloc

        if ( HISTcs2_t % coord(k) == coord) then

            ! Interpolate variables to z-plane location
            do j=1,ny
            do i=1,nx
                $if ($LVLSET)
                phi_plane(i,j) = linear_interp( phi(i,j,HISTcs2_t % istart(k)), &
                                 phi(i,j,HISTcs2_t % istart(k)+1), dz, HISTcs2_t % ldiff(k) )
                $endif
                cs2_plane(i,j) = linear_interp( Cs_opt2(i,j,HISTcs2_t % istart(k)), &
                                 Cs_opt2(i,j,HISTcs2_t % istart(k)+1), dz, HISTcs2_t % ldiff(k) )
            enddo  
            enddo

            ! Bin values
            $if ($LVLSET)
                call hist_binit( HISTcs2_t % hist_t(k), cs2_plane, phi_plane )
            $else
                call hist_binit( HISTcs2_t % hist_t(k), cs2_plane )
            $endif

        !endif

        !if ( HISTtn_t % coord(k) == coord) then
         if (sgs_model.gt.3) then
            ! Interpolate variables to z-plane location
            do j=1,ny
            do i=1,nx
                tn_plane(i,j) = linear_interp( Tn_all(i,j,HISTtn_t % istart(k)), &
                                 Tn_all(i,j,HISTtn_t % istart(k)+1), dz, HISTtn_t % ldiff(k) )
            enddo  
            enddo            

            ! Bin values
            $if ($LVLSET)
                ! For now, use same z-planes for all variables so phi's are the same
                call hist_binit( HISTtn_t % hist_t(k), tn_plane, phi_plane )
            $else
                call hist_binit( HISTtn_t % hist_t(k), tn_plane )
            $endif

         endif
        !endif 

        !if ( HISTnu_t % coord(k) == coord) then

            ! Interpolate variables to z-plane location
            do j=1,ny
            do i=1,nx
                nu_plane(i,j) = linear_interp( Nu_t(i,j,HISTnu_t % istart(k)), &
                                 Nu_t(i,j,HISTnu_t % istart(k)+1), dz, HISTnu_t % ldiff(k) )
            enddo  
            enddo

            ! Bin values
            $if ($LVLSET)
                ! For now, use same z-planes for all variables so phi's are the same
                call hist_binit( HISTnu_t % hist_t(k), nu_plane, phi_plane )
            $else
                call hist_binit( HISTnu_t % hist_t(k), nu_plane )
            $endif

        !endif 

        !if ( HISTee_t % coord(k) == coord) then
         if (sgs_model.gt.1) then
            ! Interpolate variables to z-plane location
            do j=1,ny
            do i=1,nx
                ee_plane(i,j) = linear_interp( ee_now(i,j,HISTee_t % istart(k)), &
                                 ee_now(i,j,HISTee_t % istart(k)+1), dz, HISTee_t % ldiff(k) )
            enddo  
            enddo

            ! Bin values
            $if ($LVLSET)
                ! For now, use same z-planes for all variables so phi's are the same
                call hist_binit( HISTee_t % hist_t(k), ee_plane, phi_plane )
            $else
                call hist_binit( HISTee_t % hist_t(k), ee_plane )
            $endif
          endif

        endif

    enddo

end subroutine sgs_hist_update_vals
!---------------------------------------------------------------------

!---------------------------------------------------------------------
subroutine sgs_hist_finalize()
! Output data and deallocate arrays

use string_util
implicit none
include 'tecryte.h'

    integer :: k

! Write structures to file (to continue analysis, if necessary)
    fname_histp=''; call string_concat( fname_histp, fname_histp_base // '.c', coord )
    fname_histv=''; call string_concat( fname_histv, fname_histv_base // '.c', coord )

    open(unit=1,file=fname_histp,position='rewind',action='write',form='formatted')
        write(1,*) nproc, nz_tot, L_z
        write(1,*) sgs_hist_nloc
        write(1,*) sgs_hist_loc
        write(1,*) cs2_bmin, cs2_bmax, cs2_nbins
        write(1,*) tn_bmin, tn_bmax, tn_nbins
        write(1,*) nu_bmin, nu_bmax, nu_nbins
        write(1,*) ee_bmin, ee_bmax, ee_nbins
    close(1)

    open(unit=1,file=fname_histv,position='rewind',action='write',form='unformatted')
        do k=1,sgs_hist_nloc
            write(1) HISTcs2_t%hist_t(k)%bins, HISTcs2_t%hist_t(k)%vals
            write(1) HISTtn_t%hist_t(k)%bins, HISTtn_t%hist_t(k)%vals
            write(1) HISTnu_t%hist_t(k)%bins, HISTnu_t%hist_t(k)%vals
            write(1) HISTee_t%hist_t(k)%bins, HISTee_t%hist_t(k)%vals
        enddo
    close(1)

! Normalize and write final curves to file
do k=1,sgs_hist_nloc

    write(cl,'(F9.4)') sgs_hist_loc(k)

    if ( HISTcs2_t % coord(k) == coord) then

        ! Normalize (integrates to one) - do not include junk bins at each end
            HISTcs2_t % hist_t(k) % vals = HISTcs2_t % hist_t(k) % vals / &
                ( sum( HISTcs2_t % hist_t(k) % vals(1:HISTcs2_t%hist_t(k)%nbins) ) * HISTcs2_t % hist_t(k) % db )

        ! Write to file
            ! Update file names to include z-plane location
            write(fname_cs2,*) fname_cs2_base,trim(adjustl(cl)),'.dat'
            fname_cs2=trim(adjustl(fname_cs2))

            ! Write bins and vals to file
            call write_tecplot_header_ND(fname_cs2, 'rewind', 2, (/ HISTcs2_t % hist_t(k) % nbins+2 /), &
                '"Cs2", "pdf"', numtostr(k,6), 2)
            call write_real_data_1D(fname_cs2, 'append', 'formatted', 1, HISTcs2_t % hist_t(k) % nbins+2, &
                (/ HISTcs2_t % hist_t(k) % vals /), 0, HISTcs2_t % hist_t(k) % bins)

    !endif

    !if ( HISTtn_t % coord(k) == coord) then
      if (sgs_model.gt.3) then
        ! Normalize (integrates to one)
            HISTtn_t % hist_t(k) % vals = HISTtn_t % hist_t(k) % vals / &
                ( sum( HISTtn_t % hist_t(k) % vals(1:HISTtn_t%hist_t(k)%nbins) ) * HISTtn_t % hist_t(k) % db )

        ! Write to file
            ! Update file names to include z-plane location
            write(fname_tn,*) fname_tn_base,trim(adjustl(cl)),'.dat'
            fname_tn=trim(adjustl(fname_tn))

            ! Write bins and vals to file
            call write_tecplot_header_ND(fname_tn, 'rewind', 2, (/ HISTtn_t % hist_t(k) % nbins+2 /), &
                '"T", "pdf"', numtostr(k,6), 2)
            call write_real_data_1D(fname_tn, 'append', 'formatted', 1, HISTtn_t % hist_t(k) % nbins+2, &
                (/ HISTtn_t % hist_t(k) % vals /), 0, HISTtn_t % hist_t(k) % bins)
      endif
    !endif

    !if ( HISTnu_t % coord(k) == coord) then

        ! Normalize (integrates to one)
            HISTnu_t % hist_t(k) % vals = HISTnu_t % hist_t(k) % vals / &
                ( sum( HISTnu_t % hist_t(k) % vals(1:HISTnu_t%hist_t(k)%nbins) ) * HISTnu_t % hist_t(k) % db )

        ! Write to file
            ! Update file names to include z-plane location
            write(fname_nu,*) fname_nu_base,trim(adjustl(cl)),'.dat'
            fname_nu=trim(adjustl(fname_nu))

            ! Write bins and vals to file
            call write_tecplot_header_ND(fname_nu, 'rewind', 2, (/ HISTnu_t % hist_t(k) % nbins+2 /), &
                '"Nu_t", "pdf"', numtostr(k,6), 2)
            call write_real_data_1D(fname_nu, 'append', 'formatted', 1, HISTnu_t % hist_t(k) % nbins+2, &
                (/ HISTnu_t % hist_t(k) % vals /), 0, HISTnu_t % hist_t(k) % bins)

    !endif

    !if ( HISTee_t % coord(k) == coord) then
      if (sgs_model.gt.1) then
        ! Normalize (integrates to one)
            HISTee_t % hist_t(k) % vals = HISTee_t % hist_t(k) % vals / &
                ( sum( HISTee_t % hist_t(k) % vals(1:HISTee_t%hist_t(k)%nbins) ) * HISTee_t % hist_t(k) % db )

        ! Write to file
            ! Update file names to include z-plane location
            write(fname_ee,*) fname_ee_base,trim(adjustl(cl)),'.dat'
            fname_ee=trim(adjustl(fname_ee))

            ! Write bins and vals to file
            call write_tecplot_header_ND(fname_ee, 'rewind', 2, (/ HISTee_t % hist_t(k) % nbins+2 /), &
                '"ee", "pdf"', numtostr(k,6), 2)
            call write_real_data_1D(fname_ee, 'append', 'formatted', 1, HISTee_t % hist_t(k) % nbins+2, &
                (/ HISTee_t % hist_t(k) % vals /), 0, HISTee_t % hist_t(k) % bins)
      endif
    endif

enddo 

! Deallocate
    deallocate ( HISTcs2_t % istart )
    deallocate ( HISTtn_t % istart )
    deallocate ( HISTnu_t % istart )
    deallocate ( HISTee_t % istart )

    deallocate ( HISTcs2_t % coord )
    deallocate ( HISTtn_t % coord )
    deallocate ( HISTnu_t % coord )
    deallocate ( HISTee_t % coord )

    deallocate ( HISTcs2_t % ldiff )
    deallocate ( HISTtn_t % ldiff )
    deallocate ( HISTnu_t % ldiff )
    deallocate ( HISTee_t % ldiff )

    deallocate ( HISTcs2_t % hist_t )
    deallocate ( HISTtn_t % hist_t )
    deallocate ( HISTnu_t % hist_t )
    deallocate ( HISTee_t % hist_t )  

end subroutine sgs_hist_finalize
!---------------------------------------------------------------------

end module sgs_hist
