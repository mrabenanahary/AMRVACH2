!> This module contains definitions of global parameters and variables and some
!> generic functions/subroutines used in AMRVAC.
!>
!> \todo Move the parameters to the relevant (physics) modules
module mod_global_parameters
  use mod_physicaldata
  use mod_connectivity
  use mpi
  use mod_constants
  use mod_variables
  use mod_basic_types

  implicit none
  public

  ! Parameters

  character(len=*), parameter :: undefined = 'undefined'

  !> @todo Move mpi related variables to e.g. mod_comm

  !> The number of MPI tasks
  integer :: npe

  !> The rank of the current MPI task
  integer :: mype

  !> The MPI communicator
  integer :: icomm

  !> A global MPI error return code
  !> @todo Make local
  integer :: ierrmpi

  !> MPI file handle for logfile
  integer :: log_fh
  !> MPI type for block including ghost cells and its size
  integer :: type_block, size_block
  !> MPI type for block coarsened by 2, and for its children blocks
  integer :: type_coarse_block, type_sub_block(2,2)
  !> MPI type for IO: block excluding ghost cells
  integer :: type_block_io, size_block_io
  !> MPI type for IO: transformed data excluding ghost cells
  integer :: type_block_io_tf, size_block_io_tf
  !> MPI type for IO: cell corner (xc) or cell center (xcc) coordinates
  integer :: type_block_xc_io,type_block_xcc_io
  !> MPI type for IO: cell corner (wc) or cell center (wcc) variables
  integer :: type_block_wc_io,type_block_wcc_io

  integer :: itag, irecv, isend
  integer, dimension(:), allocatable :: recvrequest, sendrequest
  integer, dimension(:,:), allocatable :: recvstatus, sendstatus

  ! geometry and domain setups

  !> the mesh range (within a block with ghost cells)
  integer :: ixMlo1,ixMlo2,ixMhi1,ixMhi2


  !> Indices for  cartesian coordinates
  integer :: x_ = -1
  integer :: y_ = -1
  !> Indices for cylindrical coordinates FOR TESTS, negative value when not used:
  integer :: r_ = -1
  integer :: phi_ = -1
  integer :: z_ = -1,theta_=-1

  logical :: r_in       = .false.
  logical :: phi_in     = .false.
  logical :: theta_in   = .false.
  logical :: x_in       = .false.
  logical :: y_in       = .false.
  logical :: z_in       = .false.

  !> Number of spatial dimensions for grid variables
  integer, parameter :: ndim=2

  !> Number of spatial dimensions (components) for vector variables
  integer  :: ndir=ndim

  !> minimum and maximum domain boundaries for each dimension
  real(kind=dp)     :: xprobmin1,xprobmin2,xprobmax1,xprobmax2
  real(kind=dp)     :: box_limit(2,ndim)

  !> Cartesian geometry or not
  logical :: slab

  !> number of grid blocks in domain per dimension, in array over levels
  integer, dimension(:), allocatable :: ng1,ng2
  !> extent of grid blocks in domain per dimension, in array over levels
  real(kind=dp)   , dimension(:), allocatable :: dg1,dg2

  !> number of cells for each dimension in level-one mesh
  integer :: domain_nx1,domain_nx2

  !> number of cells for each dimension in grid block excluding ghostcells
  integer :: block_nx1,block_nx2

  !> Lower index of grid block arrays (always 1)
  integer, parameter :: ixGlo1 = 1, ixGlo2 = 1

  !> Upper index of grid block arrays
  integer :: ixGhi1,ixGhi2

  !> Number of ghost cells surrounding a grid
  integer :: nghostcells

  !> Switch to use stretched grid
  logical :: stretched_grid
  !> Stretched Cartesian geometry or not
  logical :: slab_stretched
  !> Switch to set stretched dimension
  logical :: stretched_dim(ndim)
  !> Switch to set symmetrically stretched dimension
  logical :: stretched_symm_dim(ndim)
  !> stretch factor between cells at AMR level 1, per dimension
  real(kind=dp)    ::  qstretch_baselevel(ndim)
  !> (even) number of (symmetrically) stretched
  !> blocks at AMR level 1, per dimension
  integer ::  nstretchedblocks_baselevel(ndim)
  !> (even) number of (symmetrically) stretched blocks per level and dimension
  integer, allocatable ::  nstretchedblocks(:,:)
  !> physical extent of stretched border in symmetric stretching
  real(kind=dp)    :: xstretch1,xstretch2
  !> Stretching factors and first cell size for each AMR level and dimension
  real(kind=dp)   , allocatable :: qstretch(:,:), dxfirst(:,:),  dxfirst_1mq(:,&
     :), dxmid(:,:)

  !> grid hierarchy info (level and grid indices)
  integer, parameter :: nodehi=2+1
  integer, parameter :: plevel_=1
  integer, parameter :: pig1_=plevel_+1,pig2_=plevel_+2

  integer, allocatable :: node(:,:)
  integer, allocatable :: node_sub(:,:)

  !> grid location info (corner coordinates and grid spacing)
  integer, parameter :: rnodehi=3*2
  integer, parameter :: rpxmin0_=0
  integer, parameter :: rpxmin1_=rpxmin0_+1,rpxmin2_=rpxmin0_+2
  integer, parameter :: rpxmax0_=2
  integer, parameter :: rpxmax1_=rpxmax0_+1,rpxmax2_=rpxmax0_+2
  integer, parameter :: rpdx1_=2*2+1,rpdx2_=2*2+2

  !> Corner coordinates
  real(kind=dp)   , allocatable :: rnode(:,:)
  real(kind=dp)   , allocatable :: rnode_sub(:,:)

  real(kind=dp)   , allocatable :: dx(:,:)
  real(kind=dp)    :: dxlevel(ndim)

  ! IO related quantities

  !> Maximum number of saves that can be defined by tsave or itsave
  integer, parameter :: nsavehi=100

  !> Number of output methods
  integer, parameter :: nfile         = 5

  !> Names of the output methods
  character(len=40), parameter  :: output_names(nfile) = ['log      ',&
      'normal   ', 'slice    ', 'collapsed', 'analysis ']

  !> If collapse(DIM) is true, generate output integrated over DIM
  logical :: collapse(ndim)

  !> Save output of type N on times tsave(:, N)
  real(kind=dp)    :: tsave(nsavehi,nfile)

  !> \todo Move tsavelast to amrvac.t
  real(kind=dp)    :: tsavelast(nfile)

  !> Repeatedly save output of type N when dtsave(N) simulation time has passed
  real(kind=dp)    :: dtsave(nfile)

  !> Save output of type N on iterations itsave(:, N)
  integer :: itsave(nsavehi,nfile)

  !> \todo remove itsavelast?
  integer :: itsavelast(nfile)

  !> Repeatedly save output of type N when ditsave(N) time steps have passed
  integer :: ditsave(nfile)

  !> \todo Move to amrvac.t
  integer :: isavet(nfile)

  !> \todo Move to amrvac.t
  integer :: isaveit(nfile)

  !> The level at which to produce line-integrated / collapsed output
  integer :: collapseLevel

  !> Number of saved files of each type
  !> \todo Move to mod_input_output
  integer :: n_saves(1:nfile)

  !> to monitor timeintegration loop at given wall-clock time intervals
  real(kind=dp)    :: time_between_print

  !> accumulated wall-clock time spent on boundary conditions
  real(kind=dp)    :: time_bc

  !> IO: snapshot and collapsed views output numbers/labels
  integer :: snapshotnext, collapseNext, icollapse

  !> Constant indicating log output
  integer, parameter :: filelog_      = 1

  !> Constant indicating regular output
  integer, parameter :: fileout_      = 2

  !> Constant indicating slice output
  integer, parameter :: fileslice_    = 3

  !> Constant indicating collapsed output
  integer, parameter :: filecollapse_ = 4

  !> Constant indicating analysis output (see @ref analysis.md)
  integer, parameter :: fileanalysis_ = 5

  !> Unit for standard input
  integer, parameter :: unitstdin=5

  !> Unit for standard output
  integer, parameter :: unitterm=6

  !> Unit for error messages
  integer, parameter :: uniterr=6

  !> \todo Move to mod_input_output
  integer, parameter :: unitpar=9
  integer, parameter :: unitconvert=10
  integer, parameter :: unitslice=11
  integer, parameter :: unitsnapshot=12
  integer, parameter :: unitcollapse=13
  integer, parameter :: unitanalysis=14

  !> Number of auxiliary variables that are only included in output
  integer :: nwauxio

  !> IO switches for conversion
  logical              :: nocartesian
  logical, allocatable :: w_write(:)
  logical, allocatable :: writelevel(:)
  real(kind=dp)        :: writespshift(ndim,2)
  integer              :: level_io, level_io_min, level_io_max

  !> Which par files are used as input
  character(len=std_len), allocatable :: par_files(:)

  !> Base file name for simulation output, which will be followed by a number
  character(len=std_len) :: base_filename

  !> If not 'unavailable', resume from snapshot with this base file name
  character(len=std_len) :: restart_from_file

  !> Which type of log to write: 'normal', 'special', 'regression_test'
  character(len=std_len) :: typefilelog

  !> Resume from the snapshot with this index
  integer :: snapshotini

  !> If true and restart_from_file is given, convert snapshots to
  !> other file formats
  logical                :: convert

  !> If true, already convert to output format during the run
  logical                :: autoconvert

  !> If true, convert from conservative to primitive variables in output
  logical                :: saveprim

  !> Which format to use when converting
  !>
  !> Options are: tecplot, tecplotCC, vtu, vtuCC, vtuB, vtuBCC,
  !> tecplotmpi, tecplotCCmpi, vtumpi, vtuCCmpi, vtuBmpi, vtuBCCmpi, pvtumpi, pvtuCCmpi,
  !> pvtuBmpi, pvtuBCCmpi, tecline, teclinempi, onegrid
  character(len=std_len) :: convert_type

  character(len=std_len) :: collapse_type

  !> Conversion factors the primitive variables
  real(kind=dp)   , allocatable :: w_convert_factor(:)

  real(kind=dp)    :: length_convert_factor

  !> Conversion factor for time unit
  real(kind=dp)          :: time_convert_factor

  integer                :: saveigrid

  !> Stores the memory and load imbalance, used in printlog
  real(kind=dp)          :: Xload, Xmemory

  !> Save a snapshot before crash a run met unphysical values
  logical :: crash=.false.

  ! Physics factors

  !> Physical scaling factor for length
  real(kind=dp)    :: unit_length=1.d0

  !> Physical scaling factor for time
  real(kind=dp)    :: unit_time=1.d0

  !> Physical scaling factor for density
  real(kind=dp)    :: unit_density=1.d0

  !> Physical scaling factor for velocity
  real(kind=dp)    :: unit_velocity=0.d0

  !> Physical scaling factor for temperature
  real(kind=dp)    :: unit_temperature=1.d0

  !> Physical scaling factor for pressure
  real(kind=dp)    :: unit_pressure=1.d0

  !> Physical scaling factor for magnetic field
  real(kind=dp)    :: unit_magneticfield=1.d0

  !> Physical scaling factor for number density
  real(kind=dp)    :: unit_numberdensity=1.d0

  !> error handling
  real(kind=dp)    :: small_temperature,small_pressure,small_density

  !> amplitude of background dipolar, quadrupolar, octupolar, user's field
  real(kind=dp)    :: Bdip=0.d0
  real(kind=dp)    :: Bquad=0.d0
  real(kind=dp)    :: Boct=0.d0
  real(kind=dp)    :: Busr=0.d0

  !> amplitude of background time depending magnetic moment
  real(kind=dp)    :: mu0dip         =0.0_dp
  real(kind=dp)    :: mu0theta       =0.0_dp
  real(kind=dp)    :: mu0phi         =0.0_dp
  real(kind=dp)    :: mutheta_period =0.0_dp
  real(kind=dp)    :: muphi_period   =0.0_dp


  !> check and optionally fix unphysical small values (density, gas pressure)
  logical :: check_small_values=.false.

  !> split potential or linear force-free magnetic field as background B0 field
  logical :: B0field=.false.

  !> split potential or linear force-free time depending magnetic field as background B0 field
  logical :: B0field_reset =.false.

  !> Use SI units (.true.) or use cgs units (.false.)
  logical :: SI_unit=.false.

  !> Solve energy equation or not
  logical :: phys_energy=.true.

  !> Solve polytropic process instead of solving total energy
  logical :: solve_internal_e=.false.

  !> Enable to strictly conserve the angular momentum
  !> (works both in cylindrical and spherical coordinates)
  logical :: angmomfix=.false.

  !> Use particles module or not
  logical :: use_particles=.false.

  ! AMR switches

  !> The maximum number of grid blocks in a processor
  integer :: max_blocks

  !> The maximum number of levels in the grid refinement
  integer, parameter :: nlevelshi = 20

  !> Maximal number of AMR levels
  integer :: refine_max_level

  !> Weights of variables used to calculate error for mesh refinement
  real(kind=dp)   , allocatable :: w_refine_weight(:)

  !> Fix the AMR grid after this time
  real(kind=dp)    :: tfixgrid

  !> Fix the AMR grid after this many time steps
  integer :: itfixgrid

  !> Reconstruct the AMR grid once every ditregrid iteration(s)
  integer :: ditregrid

  !> refinement: lohner estimate wavefilter setting
  real(kind=dp)   , allocatable :: amr_wavefilter(:)

  integer                       :: refine_criterion
  logical                       :: prolongprimitive
  logical                       :: coarsenprimitive

  !> Error tolerance for refinement decision
  real(kind=dp)   , allocatable :: refine_threshold(:)
  real(kind=dp)   , allocatable :: derefine_ratio(:)

  !> If true, rebuild the AMR grid upon restarting
  logical :: reset_grid

  !> Number of cells as buffer zone
  !> \todo is it necessary?
  integer :: nbufferx1,nbufferx2

  integer :: levmin
  integer :: levmax
  integer :: levmax_sub

  ! Miscellaneous

  !> problem switch allowing different setups in same usr_mod.t
  integer           :: iprob

  !> Kronecker delta tensor
  integer :: kr(3,3)

  !> Levi-Civita tensor
  integer :: lvc(3,3,3)

  ! Time integration aspects

  real(kind=dp)    :: dt
  real(kind=dp)   , allocatable :: dt_grid(:)

  logical :: time_advance

  !> The Courant (CFL) number used for the simulation
  real(kind=dp)    :: courantpar

  logical          :: small_getdt_avergare =.false.
  real(kind=dp)    :: small_dt_coef        =1.0_dp
  !> How to compute the CFL-limited time step.
  !>
  !> Options are 'maxsum': max(sum(c/dx)); 'summax': sum(max(c/dx)) and
  !> 'minimum: max(c/dx), where the summations loop over the grid dimensions and
  !> c is the velocity. The default 'maxsum' is the conventiontal way of
  !> computing CFL-limited time steps.
  character(len=std_len) :: typecourant

  !> If dtpar is positive, it sets the timestep dt, otherwise courantpar is used
  !> to limit the time step based on the Courant condition.
  real(kind=dp)    :: dtpar

  !> For resistive MHD, the time step is also limited by the diffusion time:
  !> \f$ dt < dtdiffpar \times dx^2/eta \f$
  real(kind=dp)    :: dtdiffpar

  !> The global simulation time
  real(kind=dp)    :: global_time

  !> Start time for the simulation
  real(kind=dp)    :: time_init

  !> End time for the simulation
  real(kind=dp)    :: time_max

  !> Stop the simulation when the time step becomes smaller than this value
  real(kind=dp)    :: dtmin

  !> If true, reset iteration count and global_time to original values, and
  !> start writing snapshots at index 0
  logical :: reset_time

  !> If true, reset iteration count to 0
  logical :: reset_it

  !> If true, call initonegrid_usr upon restarting
  logical :: firstprocess

  !> Number of time steps taken
  integer :: it

  !> Stop the simulation after this many time steps have been taken
  integer :: it_max

  !> initial iteration count
  integer :: it_init

  !> If > 1, then in the first slowsteps-1 time steps dt is reduced
  !> by a factor \f$ 1 - (1- step/slowsteps)^2 \f$
  integer :: slowsteps

  ! Method switches

  !> Index of the sub-step in a multi-step time integrator
  integer :: istep

  !> How many sub-steps the time integrator takes
  integer :: nstep

  !> Which time integrator to use
  character(len=std_len) :: time_integrator

  !> What should be used as a basis for the limiting in TVD methods. Options are
  !> 'original', 'previous' and 'predictor'.
  !>
  !> By default, the original value is used in 1D and for dimensional splitting,
  !> while for dimensionally unsplit multidimensional case (dimsplit=F), TVDLF
  !> and TVD-MUSCL uses the previous value from wold for limiting.
  character(len=std_len) :: typelimited

  !> How to apply dimensional splitting to the source terms, see
  !> @ref disretization.md
  character(len=std_len) :: typesourcesplit

  !> Which spatial discretization to use (per grid level)
  character(len=std_len), allocatable :: flux_scheme(:)

  !> The spatial dicretization for the predictor step when using a two
  !> step method
  character(len=std_len), allocatable :: typepred1(:)

  !> Type of slope limiter used for reconstructing variables on cell edges
  integer, allocatable :: type_limiter(:)

  !> Type of slope limiter used for computing gradients or divergences, when
  !> typegrad or typediv are set to 'limited'
  integer, allocatable :: type_gradient_limiter(:)

  !> \todo Remove / replace with limiter
  integer :: typelimiter

  !> \todo Remove / replace with gradient_limiter
  integer :: typegradlimiter

  !> Limiter used for prolongation to refined grids and ghost cells
  character(len=std_len) :: typeprolonglimit

  !> Which type of entropy fix to use with Riemann-type solvers
  character(len=std_len), allocatable :: typeentropy(:)

  !> Which type of TVD method to use
  character(len=std_len) :: typetvd

  !> Which type of TVDLF method to use
  character(len=std_len) :: typeboundspeed

  character(len=std_len) :: typeaverage
  character(len=std_len) :: typedimsplit
  character(len=std_len) :: typeaxial='default'
  character(len=std_len) :: typepoly

  integer                       :: nxdiffusehllc
  real(kind=dp)   , allocatable :: entropycoef(:)
  real(kind=dp)                 :: tvdlfeps
  logical, allocatable          :: loglimit(:), logflag(:)
  logical                       :: flathllc,flatcd,flatsh
  !> Use split or unsplit way to add user's source terms, default: unsplit
  logical                       :: source_split_usr
  logical                       :: dimsplit

  character(len=std_len) :: typediv,typegrad,typecurl

  !> global fastest wave speed needed in fd scheme and glm method
  real(kind=dp)    :: cmax_global

  !> global fastest flow speed needed in glm method
  real(kind=dp)    :: vmax_global

  !> need global maximal wave speed
  logical :: need_global_cmax=.false.

  !> need global maximal flow speed
  logical :: need_global_vmax=.false.

  ! Boundary region parameters

  !> True for dimensions with periodic boundaries
  logical :: periodB(ndim)

  !> Indicates whether there is a pole at a boundary
  logical :: poleB(2,ndim)

  !> True for dimensions with aperiodic boundaries
  logical :: aperiodB(ndim)

  !> True for save physical boundary cells in dat files
  logical :: save_physical_boundary

  !> Array indicating the type of boundary condition per variable and per
  !> physical boundary
  character(len=std_len), allocatable :: typeboundary(:, :)

  character(len=std_len) :: typeghostfill='linear',prolongation_method
  logical :: internalboundary

  ! parameters for bc_phys
  integer, parameter :: ismin1=-1+2*1,ismin2=-1+2*2
  integer, parameter :: ismax1=2*1,ismax2=2*2

  type fluxalloc
     real(kind=dp)   , dimension(:,:,:), pointer:: flux => null()
  end type fluxalloc
  !> store flux to fix conservation
  type(fluxalloc), dimension(:,:,:), allocatable :: pflux

  logical, allocatable :: phyboundblock(:)

  !> set logical false and true array fro the grid cell
  logical, allocatable, target, save    :: grid_logical_false(:,:)
  logical, allocatable, target, save    :: grid_logical_true(:,:)
  !$OMP THREADPRIVATE(dxlevel)
  !$OMP THREADPRIVATE(saveigrid)
  !$OMP THREADPRIVATE(typelimiter,typegradlimiter)

contains

  !> Cross product of two vectors
  pure subroutine cross_product(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,a,b,axb)
    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2
    real(kind=dp)   , intent(in) :: a(ixImin1:ixImax1,ixImin2:ixImax2,3),&
        b(ixImin1:ixImax1,ixImin2:ixImax2,3)
    real(kind=dp)   , intent(out) :: axb(ixImin1:ixImax1,ixImin2:ixImax2,3)
    !-------------------------------------------------------------------------

    axb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)=a(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       2)*b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,3)-a(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,3)*b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)
    axb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)=a(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       3)*b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)-a(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,1)*b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,3)
    axb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,3)=a(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       1)*b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)-a(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,2)*b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)
  end subroutine cross_product

end module mod_global_parameters
