!///////////////////////////////////////////////////////////////////////
module hdf5_io
!///////////////////////////////////////////////////////////////////////
! This module contains all the subroutines for writing HDF5 
! 

  $if($MPI)
    use mpi
  $endif   

  use types,only:rprec
  use HDF5 ! HDF5 data library  
  use h5lt ! Easier to use hdf5 implementation
  use grid_defs, only : grid   ! Domain coordinates
  use param, only : nx, ny, nz, ierr, comm, coord, jt_total, & 
  total_time, path
  use sim_param, only : u,v,w ! velocity field
  use string_util
  use messages

  implicit none

  ! String to save the name of the HDF5 output file
  character (64) :: fname_hdf5_domain
  character(LEN=40) :: grid_file  ! Name of file to write out

  ! Subroutines to write HDF5 data and XDMF description 
  contains


subroutine get_file_name()
!***********************************************************************
!
! This subroutine stores the name of the hdf5 file to be used
! It is intended to be called everytime fo that the name of the file
! gets updated every time-step
  character(len=5) :: coord_str ! The string character of comm variable

  ! Erase the value stored in the variable
  fname_hdf5_domain=''
  
  ! Grid name
  write(coord_str, '(i4.4)') coord
  grid_file='./output/grid.c'//trim(coord_str)//'.h5'
  
  ! Create the file name used for the output
  ! .hdf5 is for hdf5 format
  $if( $BINARY )
  call string_splice( fname_hdf5_domain, path // 'output/binary.', jt_total,'.h5')
  $else
  call string_splice( fname_hdf5_domain, path // 'output/data.', jt_total,'.h5')
  $endif

  $if ($MPI)
  call string_concat( fname_hdf5_domain, '.c', coord)
  $endif
  
  write(*,*) fname_hdf5_domain

end subroutine get_file_name  

  
subroutine write_grid_hdf5()
!***********************************************************************
!
! This subroutine writes the the velocity as an HDF5 file
! Writes velocity vector in domain section provided by the processor

  integer(HSIZE_T), dimension(3) :: data_dims ! Data size
  integer(HID_T) :: file_id   ! Error and file identifier
  integer(HID_T) :: group_id      ! Group identifier
  integer i,j,k
  real(rprec), dimension(4) :: data_coord(nx,ny,nz,3)
  real(rprec), pointer, dimension(:) :: x,y,z
  logical :: file_exists

  ! This will create the names for the the grid and data files
  call get_file_name()
    
  ! Find out if the grid file exists, if not create it
  inquire(file=trim(grid_file), exist=file_exists)

  ! This will create the first instance of the file
  ! This subroutine only runs if the grid file does not exist
  if(.not. file_exists) then
    ! Name of the grid file
  
    ! Set pointers
    x => grid % x
    y => grid % y
    z => grid % z

    ! Initialize values
    data_dims = (/nx,ny,nz/)

    do k=1, nz
      do j=1, ny
        do i=1, nx
          data_coord(i,j,k,1) = x(i)
          data_coord(i,j,k,2) = y(j)
          data_coord(i,j,k,3) = z(k)
        enddo
      enddo
    enddo
  
    !  write(*,*) 'Not working',x,y,z

    ! Initialize FORTRAN interface.
    call h5open_f(ierr) 

    ! Create a new file using default properties.
    call h5fcreate_f(trim(grid_file), H5F_ACC_TRUNC_F, file_id, ierr) 

    ! Create the group
    call h5gcreate_f(file_id, 'Coordinates', group_id, ierr)

    ! Write coordinates 
    call h5LTmake_dataset_double_f(group_id, 'X', 3, data_dims, & 
    data_coord(1:nx,1:ny,1:nz,1) ,ierr)
    call h5LTmake_dataset_double_f(group_id, 'Y', 3, data_dims, & 
    data_coord(1:nx,1:ny,1:nz,2) ,ierr)
    call h5LTmake_dataset_double_f(group_id, 'Z', 3, data_dims, & 
    data_coord(1:nx,1:ny,1:nz,3) ,ierr)

    ! Close the group.
    call h5gclose_f(group_id, ierr)

    ! Close the file.
    call h5fclose_f(file_id, ierr)

    ! Close FORTRAN interface.
    call h5close_f(ierr)

  endif
    
    $if($MPI)
      ! Ensure that all processes finish before attempting to write 
      ! additional files. Otherwise it may flood the system with 
      ! too many I/O requests and crash the process 
      call mpi_barrier( comm, ierr )
    $endif
end subroutine write_grid_hdf5

subroutine write_velocity_domain_hdf5()
!***********************************************************************
!
! This subroutine writes the the velocity as an HDF5 file
! Writes velocity vector in domain section provided by the processor

   integer(HSIZE_T), dimension(3) :: data_dims !(/nx,ny,nz/)
   integer(HID_T) :: file_id   ! Error and file identifier
   integer(HID_T) :: group_id      ! Group identifier
   
  ! Initialize values
  data_dims = (/nx,ny,nz/)

    
  ! Initialize FORTRAN interface.
  call h5open_f(ierr) 

  ! Create a new file using default properties.
  call h5fcreate_f(trim(fname_hdf5_domain), H5F_ACC_TRUNC_F, file_id, ierr) 

  ! Write the dataset.
  call h5gcreate_f(file_id, 'velocity', group_id, ierr)
  
  call h5LTmake_dataset_double_f(group_id, 'u', 3, data_dims, &
  u(1:nx,1:ny,1:nz) ,ierr)
  call h5LTmake_dataset_double_f(group_id, 'v', 3, data_dims, &
  v(1:nx,1:ny,1:nz) ,ierr)
  call h5LTmake_dataset_double_f(group_id, 'w', 3, data_dims, &
  w(1:nx,1:ny,1:nz) ,ierr)

  ! Close the group.
  call h5gclose_f(group_id, ierr)

  ! Close the file.
  call h5fclose_f(file_id, ierr)

  ! Close FORTRAN interface.
  call h5close_f(ierr)
    
  $if($MPI)
    ! Ensure that all processes finish before attempting to write 
    ! additional files. Otherwise it may flood the system with 
    ! too many I/O requests and crash the process 
    call mpi_barrier( comm, ierr )
  $endif
 call writexdmf()

end subroutine write_velocity_domain_hdf5

subroutine write_domain_hdf5() ! 
!***********************************************************************
!
! This subroutine calls the subroutine for writing all the fields in
! the domain (velocity, pressure, etc)

  call write_velocity_domain_hdf5()
  
end subroutine write_domain_hdf5 

subroutine writexdmf()
!***********************************************************************
!
! This subroutine writes the xml file describing the HDF5 data
! It is used by Paraview to read the data

  character(40) :: file_xmf ! Name of xmf file
  integer :: file_id
  logical :: file_exists
   
!  file_xmf='./output/data_output.c'//coord_str//'.',jt_total,'.xmf'
  call string_splice(file_xmf,'./output/data_output.c',coord,'.',jt_total,'.xmf')
!  inquire(file=trim(file_xmf), exist=file_exists)

  ! This will create the first instance of the file
!  if(.not. file_exists) then
  open(file_id,file=trim(file_xmf))
  write(file_id,10) total_time, nz,ny,nx
!  close(file_id)
  
!  else

  ! Open the existing xmf file 
  open(file_id,file=trim(file_xmf), ACCESS = 'APPEND',STATUS = 'old')
  ! The variable name grid_file(10:) is used to eliminate 
  ! the path "./output/" from the beginning
  write(file_id,20) nz,ny,nx, trim(grid_file(10:))//':/Coordinates/X'
  write(file_id,20) nz,ny,nx, trim(grid_file(10:))//':/Coordinates/Y'
  write(file_id,20) nz,ny,nx, trim(grid_file(10:))//':/Coordinates/Z'
  write(file_id,30)

  write(file_id,40) 'X', nz,ny,nx, trim(fname_hdf5_domain(10:))//':/velocity/u'
  write(file_id,40) 'Y', nz,ny,nx, trim(fname_hdf5_domain(10:))//':/velocity/v'
  write(file_id,40) 'Z', nz,ny,nx, trim(fname_hdf5_domain(10:))//':/velocity/w'

  write(file_id,50)
!  endif 

  close(file_id)

! Formating to write xdmf
10  format( '<?xml version="1.0" ?>',/,&
            '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" [] >',/,&
            '<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="2">',/,&
            '  <Domain Name="les-go">',/,&
            '  <Grid GridType="Collection" CollectionType="Temporal">',/,&
            '    <Grid Name="grid" GridType="Uniform">',/,&
            '      <Time Type="Single" Value=" ',  f15.12 , ' "/>',/,&
            '      <Topology Type="3DSMesh" NumberOfElements=" ', i9, i9, i9, ' "/> ',/,&
            '         <Geometry GeometryType="X_Y_Z"> ')
            
20  format( '        <DataItem Dimensions="',i9,i9,i9, ' " NumberType="Double" Precision="8" Format="HDF">',/,&
            '          'A,/,& 
            '          </DataItem>')

30  format( '         </Geometry>')


40  format( '         <Attribute Name="Velocity',A,'" AttributeType="Scalar" Center="Node">',/,&
            '           <DataItem Dimensions="',i9,i9, i9, ' " NumberType="Double" Precision="8" Format="HDF"> ',/,&
            '          'A,/,&
            '           </DataItem>',/,&
            '           </Attribute>')

50  format( '    </Grid>',/,&
            '    </Grid>',/,&
            '  </Domain>',/,&
            '</Xdmf>')
end subroutine writexdmf




end module hdf5_io
