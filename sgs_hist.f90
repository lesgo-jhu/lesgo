module sgs_hist

use types, only: rprec
use param
use stat_defs, only: HISTcs2, HISTtn, HISTnu, HISTee

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

use grid_defs, only : grid
use functions, only : cell_indx
use string_util, only : string_concat
implicit none
include 'tecryte.h'

integer :: k, n
logical :: exstp = .false.
logical :: exstv = .false.
integer :: nproc_test, nz_tot_test
real(rprec) :: L_z_test
integer :: fid1, fid2

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

        fid1 = open_file ( fname_histp, 'rewind', 'formatted' )
          read(fid1,*) nproc_test, nz_tot_test, L_z_test
          if ( (nproc_test.ne.nproc) .or. (nz_tot_test.ne.nz_tot) .or. (L_z_test.ne.L_z) ) then
            write(*,*) 'Sgs-hist cumulative time error: nproc, nz_tot, and L_z must match'
            write(*,*) 'Starting sgs-histograms from scratch'
            exstp = .false.
          else
            read(fid1,*) sgs_hist_nloc
            deallocate ( sgs_hist_loc )
            allocate ( sgs_hist_loc(sgs_hist_nloc) )
            read(fid1,*) sgs_hist_loc
            read(fid1,*) cs2_bmin, cs2_bmax, cs2_nbins
            read(fid1,*) tn_bmin, tn_bmax, tn_nbins
            read(fid1,*) nu_bmin, nu_bmax, nu_nbins
            read(fid1,*) ee_bmin, ee_bmax, ee_nbins
          endif
        close(fid1)

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
    allocate ( HISTcs2 % istart (sgs_hist_nloc) )
    allocate ( HISTtn % istart (sgs_hist_nloc) )
    allocate ( HISTnu % istart (sgs_hist_nloc) )
    allocate ( HISTee % istart (sgs_hist_nloc) )

    allocate ( HISTcs2 % coord (sgs_hist_nloc) )
    allocate ( HISTtn % coord (sgs_hist_nloc) )
    allocate ( HISTnu % coord (sgs_hist_nloc) )
    allocate ( HISTee % coord (sgs_hist_nloc) )

    allocate ( HISTcs2 % ldiff (sgs_hist_nloc) )
    allocate ( HISTtn % ldiff (sgs_hist_nloc) )
    allocate ( HISTnu % ldiff (sgs_hist_nloc) )
    allocate ( HISTee % ldiff (sgs_hist_nloc) )

    allocate ( HISTcs2 % hist (sgs_hist_nloc) )
    allocate ( HISTtn % hist (sgs_hist_nloc) )
    allocate ( HISTnu % hist (sgs_hist_nloc) )
    allocate ( HISTee % hist (sgs_hist_nloc) )


! Determine which coords have each location
    ! Initialize 
    HISTcs2 % istart(:) = -1
    HISTcs2 % coord(:) = -1
    HISTcs2 % ldiff(:) = 0._rprec

    HISTtn % istart(:) = -1
    HISTtn % coord(:) = -1
    HISTtn % ldiff(:) = 0._rprec

    HISTnu % istart(:) = -1
    HISTnu % coord(:) = -1
    HISTnu % ldiff(:) = 0._rprec

    HISTee % istart(:) = -1
    HISTee % coord(:) = -1
    HISTee % ldiff(:) = 0._rprec

    ! For each zplane, find the coord where it belongs and set
    !   variables in that coord only   
    do k=1,sgs_hist_nloc

        ! If in this processor
        if(sgs_hist_loc(k) >= grid%z(1) .and. sgs_hist_loc(k) < grid%z(nz)) then

            ! Determine istart and ldiff (and set coord for future reference)
            HISTcs2 % coord(k) = coord    ! should be -1 (default) if not using MPI
            HISTcs2 % istart(k) = cell_indx('k',dz,sgs_hist_loc(k))
            HISTcs2 % ldiff(k) = sgs_hist_loc(k) - grid%z(HISTcs2 % istart(k))

            HISTtn % coord(k) = coord
            HISTtn % istart(k) = cell_indx('k',dz,sgs_hist_loc(k))
            HISTtn % ldiff(k) = sgs_hist_loc(k) - grid%z(HISTtn % istart(k))

            HISTnu % coord(k) = coord
            HISTnu % istart(k) = cell_indx('k',dz,sgs_hist_loc(k))
            HISTnu % ldiff(k) = sgs_hist_loc(k) - grid%z(HISTnu % istart(k))

            HISTee % coord(k) = coord
            HISTee % istart(k) = cell_indx('k',dz,sgs_hist_loc(k))
            HISTee % ldiff(k) = sgs_hist_loc(k) - grid%z(HISTee % istart(k))


            ! ALLOCATE
            ! For reference: bin_index given by:
            !   bin_index = min( ceil( max(value-bmin,-0.5_rprec) /db ), nbins+1 )
            ! Also:
            !   0 catches values below bmin
            !   nbins+1 catches values above bmax

            allocate ( HISTcs2 % hist(k) % bins (0:cs2_nbins+1) )
            allocate ( HISTcs2 % hist(k) % vals (0:cs2_nbins+1) )

            allocate ( HISTtn % hist(k) % bins (0:tn_nbins+1) )
            allocate ( HISTtn % hist(k) % vals (0:tn_nbins+1) )

            allocate ( HISTnu % hist(k) % bins (0:nu_nbins+1) )
            allocate ( HISTnu % hist(k) % vals (0:nu_nbins+1) )

            allocate ( HISTee % hist(k) % bins (0:ee_nbins+1) )
            allocate ( HISTee % hist(k) % vals (0:ee_nbins+1) )


            ! Define bins range for each variable
            !   For now this assumes all levels have the same range

                ! Number of bins
                HISTcs2 % hist(k) % nbins = cs2_nbins
                HISTtn  % hist(k) % nbins = tn_nbins
                HISTnu  % hist(k) % nbins = nu_nbins
                HISTee  % hist(k) % nbins = ee_nbins

                ! Bin minimum
                HISTcs2 % hist(k) % bmin = cs2_bmin
                HISTtn  % hist(k) % bmin = tn_bmin
                HISTnu  % hist(k) % bmin = nu_bmin
                HISTee  % hist(k) % bmin = ee_bmin

                ! Bin maximum
                HISTcs2 % hist(k) % bmax = cs2_bmax
                HISTtn  % hist(k) % bmax = tn_bmax
                HISTnu  % hist(k) % bmax = nu_bmax
                HISTee  % hist(k) % bmax = ee_bmax

                ! Set spacing between bins
                HISTcs2 % hist(k) % db = (cs2_bmax - cs2_bmin) / cs2_nbins
                HISTtn  % hist(k) % db = (tn_bmax - tn_bmin) / tn_nbins
                HISTnu  % hist(k) % db = (nu_bmax - nu_bmin) / nu_nbins
                HISTee  % hist(k) % db = (ee_bmax - ee_bmin) / ee_nbins

            ! Set bins variable to center of each bin range
                forall(n=0:cs2_nbins+1) HISTcs2 % hist(k) % bins(n) = & 
                                        cs2_bmin + (n-0.5_rprec) * HISTcs2 % hist(k) % db

                forall(n=0:tn_nbins+1) HISTtn % hist(k) % bins(n) = & 
                                       tn_bmin + (n-0.5_rprec) * HISTtn % hist(k) % db

                forall(n=0:nu_nbins+1) HISTnu % hist(k) % bins(n) = & 
                                       nu_bmin + (n-0.5_rprec) * HISTnu % hist(k) % db

                forall(n=0:ee_nbins+1) HISTee % hist(k) % bins(n) = & 
                                       ee_bmin + (n-0.5_rprec) * HISTee % hist(k) % db

            ! Initialize vals to zero
                HISTcs2 % hist(k) % vals = 0.0_rprec
                HISTtn  % hist(k) % vals = 0.0_rprec
                HISTnu  % hist(k) % vals = 0.0_rprec
                HISTee  % hist(k) % vals = 0.0_rprec

        endif

    enddo
 
! If reading from file, do so now for bins and vals
  if ( exstp .and. exstv ) then

    fid2 = open_file ( fname_histv, 'rewind', 'unformatted' )
      do k=1,sgs_hist_nloc
        read(fid2) HISTcs2%hist(k)%bins, HISTcs2%hist(k)%vals
        read(fid2) HISTtn%hist(k)%bins, HISTtn%hist(k)%vals
        read(fid2) HISTnu%hist(k)%bins, HISTnu%hist(k)%vals
        read(fid2) HISTee%hist(k)%bins, HISTee%hist(k)%vals
      enddo
   close(fid2)

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

        if ( HISTcs2 % coord(k) == coord) then

            ! Interpolate variables to z-plane location
            do j=1,ny
            do i=1,nx
                $if ($LVLSET)
                phi_plane(i,j) = linear_interp( phi(i,j,HISTcs2 % istart(k)), &
                                 phi(i,j,HISTcs2 % istart(k)+1), dz, HISTcs2 % ldiff(k) )
                $endif
                cs2_plane(i,j) = linear_interp( Cs_opt2(i,j,HISTcs2 % istart(k)), &
                                 Cs_opt2(i,j,HISTcs2 % istart(k)+1), dz, HISTcs2 % ldiff(k) )
            enddo  
            enddo

            ! Bin values
            $if ($LVLSET)
                call hist_binit( HISTcs2 % hist(k), cs2_plane, phi_plane )
            $else
                call hist_binit( HISTcs2 % hist(k), cs2_plane )
            $endif

        !endif

        !if ( HISTtn % coord(k) == coord) then
         if (sgs_model.gt.3) then
            ! Interpolate variables to z-plane location
            do j=1,ny
            do i=1,nx
                tn_plane(i,j) = linear_interp( Tn_all(i,j,HISTtn % istart(k)), &
                                 Tn_all(i,j,HISTtn % istart(k)+1), dz, HISTtn % ldiff(k) )
            enddo  
            enddo            

            ! Bin values
            $if ($LVLSET)
                ! For now, use same z-planes for all variables so phi's are the same
                call hist_binit( HISTtn % hist(k), tn_plane, phi_plane )
            $else
                call hist_binit( HISTtn % hist(k), tn_plane )
            $endif

         endif
        !endif 

        !if ( HISTnu % coord(k) == coord) then

            ! Interpolate variables to z-plane location
            do j=1,ny
            do i=1,nx
                nu_plane(i,j) = linear_interp( Nu_t(i,j,HISTnu % istart(k)), &
                                 Nu_t(i,j,HISTnu % istart(k)+1), dz, HISTnu % ldiff(k) )
            enddo  
            enddo

            ! Bin values
            $if ($LVLSET)
                ! For now, use same z-planes for all variables so phi's are the same
                call hist_binit( HISTnu % hist(k), nu_plane, phi_plane )
            $else
                call hist_binit( HISTnu % hist(k), nu_plane )
            $endif

        !endif 

        !if ( HISTee % coord(k) == coord) then
         if (sgs_model.gt.1) then
            ! Interpolate variables to z-plane location
            do j=1,ny
            do i=1,nx
                ee_plane(i,j) = linear_interp( ee_now(i,j,HISTee % istart(k)), &
                                 ee_now(i,j,HISTee % istart(k)+1), dz, HISTee % ldiff(k) )
            enddo  
            enddo

            ! Bin values
            $if ($LVLSET)
                ! For now, use same z-planes for all variables so phi's are the same
                call hist_binit( HISTee % hist(k), ee_plane, phi_plane )
            $else
                call hist_binit( HISTee % hist(k), ee_plane )
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
    integer :: fid1, fid2

! Write structures to file (to continue analysis, if necessary)
    fname_histp=''; call string_concat( fname_histp, fname_histp_base // '.c', coord )
    fname_histv=''; call string_concat( fname_histv, fname_histv_base // '.c', coord )

    fid1 = open_file ( fname_histp, 'rewind', 'formatted' )
        write(fid1,*) nproc, nz_tot, L_z
        write(fid1,*) sgs_hist_nloc
        write(fid1,*) sgs_hist_loc
        write(fid1,*) cs2_bmin, cs2_bmax, cs2_nbins
        write(fid1,*) tn_bmin, tn_bmax, tn_nbins
        write(fid1,*) nu_bmin, nu_bmax, nu_nbins
        write(fid1,*) ee_bmin, ee_bmax, ee_nbins
    close(fid1)

    fid2 = open_file ( fname_histv, 'rewind', 'unformatted' )
        do k=1,sgs_hist_nloc
            if (allocated(HISTcs2%hist(k)%bins)) write(fid2) HISTcs2%hist(k)%bins, HISTcs2%hist(k)%vals
            if (allocated(HISTtn%hist(k)%bins)) write(fid2) HISTtn%hist(k)%bins, HISTtn%hist(k)%vals
            if (allocated(HISTnu%hist(k)%bins)) write(fid2) HISTnu%hist(k)%bins, HISTnu%hist(k)%vals
            if (allocated(HISTee%hist(k)%bins)) write(fid2) HISTee%hist(k)%bins, HISTee%hist(k)%vals
        enddo
    close(fid2)


! Normalize and write final curves to file
do k=1,sgs_hist_nloc

    write(cl,'(F9.4)') sgs_hist_loc(k)

    if ( HISTcs2 % coord(k) == coord) then

        ! Normalize (integrates to one) - do not include junk bins at each end
            HISTcs2 % hist(k) % vals = HISTcs2 % hist(k) % vals / &
                ( sum( HISTcs2 % hist(k) % vals(1:HISTcs2%hist(k)%nbins) ) * HISTcs2 % hist(k) % db )

        ! Write to file
            ! Update file names to include z-plane location
            write(fname_cs2,*) fname_cs2_base,trim(adjustl(cl)),'.dat'
            fname_cs2=trim(adjustl(fname_cs2))

            ! Write bins and vals to file
            call write_tecplot_header_ND(fname_cs2, 'rewind', 2, (/ HISTcs2 % hist(k) % nbins+2 /), &
                '"Cs2", "pdf"', numtostr(k,6), 2)
            call write_real_data_1D(fname_cs2, 'append', 'formatted', 1, HISTcs2 % hist(k) % nbins+2, &
                (/ HISTcs2 % hist(k) % vals /), 0, HISTcs2 % hist(k) % bins)

    !endif

    !if ( HISTtn % coord(k) == coord) then
      if (sgs_model.gt.3) then
        ! Normalize (integrates to one)
            HISTtn % hist(k) % vals = HISTtn % hist(k) % vals / &
                ( sum( HISTtn % hist(k) % vals(1:HISTtn%hist(k)%nbins) ) * HISTtn % hist(k) % db )

        ! Write to file
            ! Update file names to include z-plane location
            write(fname_tn,*) fname_tn_base,trim(adjustl(cl)),'.dat'
            fname_tn=trim(adjustl(fname_tn))

            ! Write bins and vals to file
            call write_tecplot_header_ND(fname_tn, 'rewind', 2, (/ HISTtn % hist(k) % nbins+2 /), &
                '"T", "pdf"', numtostr(k,6), 2)
            call write_real_data_1D(fname_tn, 'append', 'formatted', 1, HISTtn % hist(k) % nbins+2, &
                (/ HISTtn % hist(k) % vals /), 0, HISTtn % hist(k) % bins)
      endif
    !endif

    !if ( HISTnu % coord(k) == coord) then

        ! Normalize (integrates to one)
            HISTnu % hist(k) % vals = HISTnu % hist(k) % vals / &
                ( sum( HISTnu % hist(k) % vals(1:HISTnu%hist(k)%nbins) ) * HISTnu % hist(k) % db )

        ! Write to file
            ! Update file names to include z-plane location
            write(fname_nu,*) fname_nu_base,trim(adjustl(cl)),'.dat'
            fname_nu=trim(adjustl(fname_nu))

            ! Write bins and vals to file
            call write_tecplot_header_ND(fname_nu, 'rewind', 2, (/ HISTnu % hist(k) % nbins+2 /), &
                '"Nu_t", "pdf"', numtostr(k,6), 2)
            call write_real_data_1D(fname_nu, 'append', 'formatted', 1, HISTnu % hist(k) % nbins+2, &
                (/ HISTnu % hist(k) % vals /), 0, HISTnu % hist(k) % bins)

    !endif

    !if ( HISTee % coord(k) == coord) then
      if (sgs_model.gt.1) then
        ! Normalize (integrates to one)
            HISTee % hist(k) % vals = HISTee % hist(k) % vals / &
                ( sum( HISTee % hist(k) % vals(1:HISTee%hist(k)%nbins) ) * HISTee % hist(k) % db )

        ! Write to file
            ! Update file names to include z-plane location
            write(fname_ee,*) fname_ee_base,trim(adjustl(cl)),'.dat'
            fname_ee=trim(adjustl(fname_ee))

            ! Write bins and vals to file
            call write_tecplot_header_ND(fname_ee, 'rewind', 2, (/ HISTee % hist(k) % nbins+2 /), &
                '"ee", "pdf"', numtostr(k,6), 2)
            call write_real_data_1D(fname_ee, 'append', 'formatted', 1, HISTee % hist(k) % nbins+2, &
                (/ HISTee % hist(k) % vals /), 0, HISTee % hist(k) % bins)
      endif
    endif

enddo 

! Deallocate
    deallocate ( HISTcs2 % istart )
    deallocate ( HISTtn % istart )
    deallocate ( HISTnu % istart )
    deallocate ( HISTee % istart )

    deallocate ( HISTcs2 % coord )
    deallocate ( HISTtn % coord )
    deallocate ( HISTnu % coord )
    deallocate ( HISTee % coord )

    deallocate ( HISTcs2 % ldiff )
    deallocate ( HISTtn % ldiff )
    deallocate ( HISTnu % ldiff )
    deallocate ( HISTee % ldiff )

    deallocate ( HISTcs2 % hist )
    deallocate ( HISTtn % hist )
    deallocate ( HISTnu % hist )
    deallocate ( HISTee % hist )  

end subroutine sgs_hist_finalize
!---------------------------------------------------------------------

end module sgs_hist
