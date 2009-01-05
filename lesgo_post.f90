!
!    ################################################################################
!    ################################################################################
!    ##                                                                            ## 
!    ##                            postprocessing.f90                              ##
!    ##                                                                            ##
!    ##                               Developed by                                 ##
!    ##                             Marcelo Chamecki                               ##
!    ##                                                                            ##
!    ##                        The Johns Hopkins University                        ##
!    ##                                                                            ##
!    ##                                03/30/2006                                  ##
!    ##                                                                            ##
!    ################################################################################
!    ################################################################################
!
!    ################################################################################
!
!    PURPOSE:  This program calculates average profiles from LES outputs.
!
!    ################################################################################
!
MODULE input_defs

! Parameters
  INTEGER,PARAMETER          :: nx=40           ! Number of x-nodes
  INTEGER,PARAMETER          :: ny=40           ! Number of y-nodes
  INTEGER,PARAMETER          :: nz=40           ! Number of z-nodes
  REAL(KIND=8),PARAMETER     :: kappa=0.4D0     ! Von Karman constant

  INTEGER,PARAMETER          :: Lz=1         ! Height of domain
  REAL(KIND=8),PARAMETER     :: dt=2.e-4    ! Timestep

  CHARACTER(len=120)         :: run='output'    ! Run to be processed
  INTEGER,PARAMETER          :: nt=1            ! Number of time outputs
  INTEGER,PARAMETER          :: tini=1     ! Initial timestep for averaging
  INTEGER,PARAMETER          :: tend=1000     ! Final timestep for averaging

namelist/parameters/&
nx,ny,nz,kappa,Lz,dt,run,nt,tini,tend

END MODULE input_defs




PROGRAM POSTPROCESSING

! No implicit variables
  IMPLICIT NONE
  use input_defs
!! Parameters
!  INTEGER,PARAMETER          :: nx=40           ! Number of x-nodes
!  INTEGER,PARAMETER          :: ny=40           ! Number of y-nodes
!  INTEGER,PARAMETER          :: nz=40           ! Number of z-nodes
 ! REAL(KIND=8),PARAMETER     :: kappa=0.4D0     ! Von Karman constant

!  INTEGER,PARAMETER          :: Lz=1         ! Height of domain
!  REAL(KIND=8),PARAMETER     :: dt=2.e-4    ! Timestep

!  CHARACTER(len=120)         :: run='output'    ! Run to be processed
!  INTEGER,PARAMETER          :: nt=1            ! Number of time outputs
!  INTEGER,PARAMETER          :: tini=1     ! Initial timestep for averaging
!  INTEGER,PARAMETER          :: tend=1000     ! Final timestep for averaging

! Main variables
  CHARACTER(len=120)         :: file            ! File name
  CHARACTER(len=120)         :: dir             ! Folder
  REAL(KIND=8),ALLOCATABLE   :: z(:,:)          ! Vertical coordinates
  REAL(KIND=8),ALLOCATABLE   :: avgtx(:,:)      ! Averages

! Variables to determine ustar
  INTEGER                    :: num             ! Number of points for time averaging
  REAL(KIND=8),ALLOCATABLE   :: t(:)            ! Time
  REAL(KIND=8),ALLOCATABLE   :: dat(:,:,:)      ! Data
  REAL(KIND=8),ALLOCATABLE   :: ustar(:)        ! Data averaged in x

! Auxiliar variables
  INTEGER                    :: i, j, k         ! Counters

!
! BEGINNING CODE   
!
!  Read input information
  open(88,file='postprocessing.inp',status='old')
  read(88,NML=parameters)
  close(88)

  ! Memory allocation
  ALLOCATE(z(3,nz),avgtx(20,nz))
  ALLOCATE(t(nt),dat(nt,3,nx),ustar(nt))
  
  ! Folder with output of run to be processed
!  dir='../'//TRIM(run)//'/'
  dir='/home/anderson/LES/LES/ANDERSON_2/output/'

  ! Read data from file
  file=TRIM(dir)//'aver_txz.out'
  OPEN(10,FILE='aver_txz.out') !,STATUS='OLD',ACTION='READ') 
  file=TRIM(dir)//'/aver_tyz.out'
  OPEN(11,FILE='aver_tyz.out') !,STATUS='OLD',ACTION='READ')    
  ! Read only first vertical level
  DO i=1,nt
    READ(10,*)t(i),(dat(i,1,k),k=1,nx)
    READ(11,*)t(i),(dat(i,2,k),k=1,nx)
    DO j=2,nz
      READ(10,*)t(i),(dat(i,3,k),k=1,nx)
      READ(11,*)t(i),(dat(i,3,k),k=1,nx)
    END DO
  END DO
  CLOSE(10)
  CLOSE(11)
  
  ! Average in x
  ustar(:)=SQRT((SUM(dat(:,1,:),2)/nx)**2+(SUM(dat(:,2,:),2)/nx)**2)
  ustar=SQRT(ustar)
  
  PRINT*,SUM(ustar)/nt

  ! Average data from all files
  file=TRIM(dir)//'aver_u.out'
  CALL AVG_X_T(nx,nz,nt,dt,tini,tend,file,ustar,avgtx(1,:))

  file=TRIM(dir)//'aver_v.out'
  CALL AVG_X_T(nx,nz,nt,dt,tini,tend,file,ustar,avgtx(2,:))

  file=TRIM(dir)//'aver_w.out'
  CALL AVG_X_T(nx,nz,nt,dt,tini,tend,file,ustar,avgtx(3,:))

  file=TRIM(dir)//'aver_dudz.out'
  CALL AVG_X_T(nx,nz,nt,dt,tini,tend,file,ustar,avgtx(4,:))

  file=TRIM(dir)//'aver_dvdz.out'
  CALL AVG_X_T(nx,nz,nt,dt,tini,tend,file,ustar,avgtx(5,:))


  ustar=ustar*ustar
  

  file=TRIM(dir)//'aver_u2.out'
  CALL AVG_X_T(nx,nz,nt,dt,tini,tend,file,ustar,avgtx(6,:))

  file=TRIM(dir)//'aver_v2.out'
  CALL AVG_X_T(nx,nz,nt,dt,tini,tend,file,ustar,avgtx(7,:))

  file=TRIM(dir)//'aver_w2.out'
  CALL AVG_X_T(nx,nz,nt,dt,tini,tend,file,ustar,avgtx(8,:))

  file=TRIM(dir)//'aver_uw.out'
  CALL AVG_X_T(nx,nz,nt,dt,tini,tend,file,ustar,avgtx(9,:))

  file=TRIM(dir)//'aver_vw.out'
  CALL AVG_X_T(nx,nz,nt,dt,tini,tend,file,ustar,avgtx(10,:))

  file=TRIM(dir)//'aver_txx.out'
  CALL AVG_X_T(nx,nz,nt,dt,tini,tend,file,ustar,avgtx(11,:))

  file=TRIM(dir)//'aver_tyy.out'
  CALL AVG_X_T(nx,nz,nt,dt,tini,tend,file,ustar,avgtx(12,:))

  file=TRIM(dir)//'aver_tzz.out'
  CALL AVG_X_T(nx,nz,nt,dt,tini,tend,file,ustar,avgtx(13,:))

  file=TRIM(dir)//'aver_txz.out'
  CALL AVG_X_T(nx,nz,nt,dt,tini,tend,file,ustar,avgtx(14,:))

  file=TRIM(dir)//'aver_tyz.out'
  CALL AVG_X_T(nx,nz,nt,dt,tini,tend,file,ustar,avgtx(15,:))
  
  ustar=1.D0

  file=TRIM(dir)//'aver_Cs.out'
  CALL AVG_X_T(nx,nz,nt,dt,tini,tend,file,ustar,avgtx(16,:))

  file=TRIM(dir)//'aver_beta_sgs.out'
  CALL AVG_X_T(nx,nz,nt,dt,tini,tend,file,ustar,avgtx(17,:))

  file=TRIM(dir)//'aver_betaclip_sgs.out'
  CALL AVG_X_T(nx,nz,nt,dt,tini,tend,file,ustar,avgtx(18,:))

  file=TRIM(dir)//'aver_Cs_Ssim.out'
  CALL AVG_X_T(nx,nz,nt,dt,tini,tend,file,ustar,avgtx(19,:))



  ! Average <u> - dimensional
  file=TRIM(dir)//'aver_u.out'
  ustar=1.D0
  CALL AVG_X_T(nx,nz,nt,dt,tini,tend,file,ustar,avgtx(20,:))

  
  
  ! Some calculations for output
  
  ! Determine z-coordinates
  DO i=1,nz
    
    ! uv-nodes
    z(1,i)=(i-0.5D0)/nz
    
    ! w-nodes
    z(2,i)=(1.D0*i)/nz
    
  END DO

  ! cs-nodes
  z(3,1)=z(1,1)
  z(3,2:nz)=z(2,1:nz-1)  
  
  ! Variances
  avgtx(6,:)=avgtx(6,:)-avgtx(1,:)*avgtx(1,:)
  avgtx(7,:)=avgtx(7,:)-avgtx(2,:)*avgtx(2,:)
  avgtx(8,:)=avgtx(8,:)-avgtx(3,:)*avgtx(3,:)
  
  ! Covariances
  avgtx(9,:)=avgtx(9,:)-avgtx(1,:)*avgtx(3,:)
  avgtx(10,:)=avgtx(10,:)-avgtx(2,:)*avgtx(3,:)

  ! Vertical gradient
  avgtx(4,:)=(kappa*z(3,:))*avgtx(4,:)
  ! Correction for the first grid node (see Porte-Agel (2000) appendix)
!  avgtx(4,1)=avgtx(4,1)/LOG(3.D0)
  
  ! Output results
  file='vertical_profiles_'//TRIM(run)//'.txt'
  OPEN(UNIT=10,FILE=file,STATUS='UNKNOWN',ACTION='WRITE')    
  DO i=1,nz
    WRITE(10,'(23F14.8)')z(1,i),       &   ! z
                         z(2,i),       &   ! z
                         z(3,i),       &   ! z
                         avgtx(1,i),   &   ! <u>
                         avgtx(2,i),   &   ! <v>
                         avgtx(3,i),   &   ! <w>
                         avgtx(4,i),   &   ! <du/dz>
                         avgtx(5,i),   &   ! <dv/dz>
                         avgtx(6,i),   &   ! var{u}
                         avgtx(7,i),   &   ! var{v}
                         avgtx(8,i),   &   ! var{w}
                         avgtx(9,i),   &   ! covar{uw}
                         avgtx(10,i),  &   ! covar{vw}
                         avgtx(11,i),  &   ! txx
                         avgtx(12,i),  &   ! tyy
                         avgtx(13,i),  &   ! tzz
                         avgtx(14,i),  &   ! txz
                         avgtx(15,i),  &   ! tyz
                         avgtx(16,i),  &   ! cs
                         avgtx(17,i),  &   ! beta
                         avgtx(18,i),  &   ! beta_clip
                         avgtx(19,i),  &   ! cs_rns
                         avgtx(20,i)       ! dimensional <u>
  END DO
  CLOSE(10) 
  
  
STOP 
END PROGRAM POSTPROCESSING



!
!    ################################################################################
!    ################################################################################
!    ##                                                                            ## 
!    ##                                  AVG_X_T                                   ##
!    ##                                                                            ##
!    ##                               Developed by                                 ##
!    ##                             Marcelo Chamecki                               ##
!    ##                                                                            ##
!    ##                        The Johns Hopkins University                        ##
!    ##                                                                            ##
!    ##                                03/30/2006                                  ##
!    ##                                                                            ##
!    ################################################################################
!    ################################################################################
!
!    ################################################################################
!
!    PURPOSE:  This routine average LES output in x and time.
!
!    ################################################################################
!
!    INPUTS:
!      n         -> number of points in the serie
!      dat       -> serie
!
!    OUTPUTS:
!      avg       -> average value
!      rms       -> root-mean-square
!
!    ################################################################################
!

SUBROUTINE AVG_X_T(nx,nz,nt,dt,tini,tend,file,ustar,avg)

! No implicit variables
  IMPLICIT none

! Main variables
  INTEGER                    :: nx              ! Number of x-nodes
  INTEGER                    :: nz              ! Number of z-nodes
  INTEGER                    :: nt              ! Number of time outputs
  INTEGER                    :: num             ! Number of points for time averaging
  INTEGER       	     :: tini		! Initial time for averaging
  INTEGER       	     :: tend		! Final time for averaging
  REAL(KIND=8)               :: dt              ! Time step
  REAL(KIND=8)               :: time            ! Time
  CHARACTER(len=120)         :: file            ! File name
  INTEGER,DIMENSION(nt)            :: t         ! Time
  REAL(KIND=8),DIMENSION(nt,nz,nx) :: dat       ! Data
  REAL(KIND=8),DIMENSION(nt,nz)    :: avgx      ! Data averaged in x
  REAL(KIND=8),DIMENSION(nz)       :: avg       ! Data averaged in x and t
  REAL(KIND=8),DIMENSION(nt)       :: ustar     ! u*(t)

! Auxiliar variables
  INTEGER                    :: i, j, k         ! Counters

!
! BEGINNING CODE   
!

  ! Read data from file
  OPEN(UNIT=10,FILE=file,STATUS='OLD',ACTION='READ')    
  DO i=1,nt
    DO j=1,nz
      READ(10,*)time,(dat(i,j,k),k=1,nx)
    END DO
    
    ! Time to iteration number
    t(i)=INT(time/dt)
    
    ! Proper normalization
    dat(i,:,:)=dat(i,:,:)/ustar(i)
    
  END DO
  CLOSE(10)
  
  ! Average in x
  avgx(:,:)=SUM(dat(:,:,:),3)/nx
  
  ! Average in time
  DO i=1,nz
    num=COUNT(t(:)>=tini .AND. t(:)<=tend)
    avg(i)=SUM(avgx(:,i),MASK=(t(:)>=tini .AND. t(:)<=tend))/num
  END DO
  

RETURN
END SUBROUTINE AVG_X_T
