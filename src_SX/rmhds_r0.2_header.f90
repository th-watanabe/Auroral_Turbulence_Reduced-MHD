MODULE RMHDS_header
!-------------------------------------------------------------------------------
!
!    Nonlinear fluxtube reduced MHD spectral code RMHDS
!
!      Header for general use in the ITG_fluxtube code
!
!      RMHDS r0.1 ( T.-H.Watanabe, July 2016)
!
!      Coding style is based on GKV-plus r2.1
!
!-------------------------------------------------------------------------------

  implicit none

  public

  integer, parameter :: DP = selected_real_kind(14)

!--------------------------------------
!  Dimension size (grid numbers)
!--------------------------------------

  integer, parameter :: nxw = 16, nyw = 16
  integer, parameter :: global_nz = 32
  integer, parameter :: nb = 2

!--------------------------------------
!  Data distribution for MPI
!--------------------------------------

  integer, parameter :: nprocz = 4

!--------------------------------------
!  Number of OMP threads
!--------------------------------------

  integer, parameter :: nomp = 2

!--------------------------------------
!  Parameters for variable sizes
!--------------------------------------

  integer, parameter :: nx = (nxw/3)*2, ny = (nyw/3)*2  

  integer, parameter :: nz = global_nz / nprocz

  integer, parameter :: nxyz  = (2*nx+1)*(ny+1)*(2*nz),    &
                        nxy   = (2*nx+1)*(ny+1)     

  integer, parameter :: nnx = nxw*2, nny = nyw*2

!--------------------------------------
!  Constants
!--------------------------------------

  real(kind=DP),    parameter :: pi  = 3.141592653589793_DP, &
                                 twopi = pi * 2._DP,         &
                                 eps = 0.0000000001_DP
  complex(kind=DP), parameter :: ui  = ( 0._DP, 1._DP )


!--------------------------------------
!  Parameters for time
!--------------------------------------

  real(kind=DP) :: e_limit                           ! elapsed time limit of a job
  real(kind=DP) :: tend                              ! end time
  real(kind=DP) :: dtout_mag, dtout_ion, dtout_eng   ! time-spacing for output


!--------------------------------------
!  Configuration parameters to be 
!    initialized in init subroutine
!--------------------------------------

  real(kind=DP), dimension(-nx:nx,0:ny,-nz:nz-1) :: ksq, ksqi

  real(kind=DP), dimension(-nx:nx)          :: kx
  real(kind=DP), dimension(0:ny)            :: ky
  real(kind=DP), dimension(-nz:nz-1)        :: zz
  real(kind=DP), dimension(-nz:nz-1)        :: dpara, valf, b0

  complex(kind=DP), dimension(0:ny)       :: ck
  integer, dimension(0:ny)                :: dj

! --- box size
  real(kind=DP) :: lx, ly, lz

! --- basic parameters for reduced MHD
  real(kind=DP) :: nu,   &    ! viscosity
                   eta,  &    ! resitivity
                   s_hat      ! shear parameter

! --- parameters for M-I coupling
  real(kind=DP) :: e0,   &    ! convection electric field
                   idns0,&    ! ionospheric number density
                   mp,   &    ! Pedersen mobility normalized by the ExB mobility
                   mh,   &    ! Hall mobility normalized by the ExB mobility
                   dperp,&    ! diffusion coefficient
                   alpha      ! normalized recombination rate : alpha * n0

  real(kind=DP) :: dt


!--------------------------------------
!  Type of calculation
!--------------------------------------

  character(9)  :: calc_type, fd_type, fd_filt
  character(5)  :: bndry_type

!--------------------------------------
!  Parameters for numerical settings
!--------------------------------------

  integer :: inum
  integer :: loop

! --- unit numbers for I/O
  integer, parameter :: inml =  9, &
                        olog = 10, &
                        icnt = 20, &
                        omag = 30, &
                        oion = 31, &
                        ocnt = 50, &
                        omph = 60, &
                        omps = 61, &
                        oidn = 62, &
                        oicr = 63, &
                        oiog = 64


END MODULE RMHDS_header
