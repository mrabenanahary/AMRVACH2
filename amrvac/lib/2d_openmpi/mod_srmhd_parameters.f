!> Paramter module for Special Relativistic Magneto-hydrodynamics
module mod_srmhd_parameters
  ! made by Z. MELIANI 14/02/2018
  use mod_global_parameters, only: std_len
  use mod_constants
  use mod_physics
  implicit none

    !> Whether thermal conduction is used
  logical, public                  :: srmhd_thermal_conduction = .false.

  !> Whether radiative cooling is added
  logical, public                  :: srmhd_radiative_cooling = .false.

   !> Whether sole energy equation
  logical, public                  :: srmhd_energy = .true.
  !> Whether synge eos  is added
  logical, public                  :: srmhd_eos = .false.

  !> Whether viscosity is added
  logical, public                  :: srmhd_viscosity = .false.

  !> Whether gravity is added
  logical, public                  :: srmhd_gravity = .false.

  !> Whether Hall-MHD is used
  logical, public                  :: srmhd_Hall = .false.

  !> Whether particles module is added
  logical, public                  :: srmhd_particles = .false.

  !> Whether magnetofriction is added
  logical, public                  :: srmhd_magnetofriction = .false.

  !> Whether GLM-MHD is used
  logical, public                  :: srmhd_glm = .false.

  !> dust not implimented in srmhd
  logical, public                  :: srmhd_dust = .false.

  !> Whether divB cleaning sources are added splitting from fluid solver
  logical, public                  :: source_split_divb = .false.

  !> GLM-MHD parameter: ratio of the diffusive and advective time scales for div b
  !> taking values within [0, 1]
  real(kind=dp)   , public                :: srmhd_glm_alpha = 0.5d0

  !> MHD fourth order
  logical, public                  :: srmhd_4th_order = .false.

  !> Number of tracer species
  integer, public                  :: srmhd_n_tracer = 0

  !> Index of the moving frame density (in the w array)
  integer, public                  :: rho_

   !> Index of the lab frame density (in the w array)
  integer, public                  :: d_
  !> Indices of the momentum density
  integer, allocatable, public     :: mom(:)

  !> Index of the energy density (-1 if not present)
  integer, public                  :: e_

  !> Index of the gas pressure (-1 if not present) should equal e_
  integer, public                  :: p_


  !> Indices of the magnetic field
  integer, allocatable, public     :: mag(:)

  !> Indices of the GLM psi
  integer, public     :: psi_
  !> Indices of the Lorentz factor
  integer, public     :: lfac_

  !> Indices of the inertia
  integer, public     :: xi_



  !> Indices of the tracers
  integer, allocatable, public     :: tracer(:)

  !> The adiabatic index
  real(kind=dp)   , public                :: srmhd_gamma = 5.d0/3.0d0

  !> The adiabatic constant
  real(kind=dp)   , public                :: srmhd_adiab = 1.0d0

  !> The MHD resistivity
  real(kind=dp)   , public                :: srmhd_eta = 0.0d0

  !> The MHD hyper-resistivity
  real(kind=dp)   , public                :: srmhd_eta_hyper = 0.0d0

  !> TODO: what is this?
  real(kind=dp)   , public                :: srmhd_etah = 0.0d0


  logical,          public         :: srmhd_checkNR   = .true.
  real(kind=dp)   , public         :: srmhd_absaccnr  = 1.0d-8
  real(kind=dp)   , public         :: srmhd_tolerNR   = 1.0d-9
  real(kind=dp)   , public         :: srmhd_maxdspeed = 1.0d-7
  real(kind=dp)   , public         :: srmhd_maxspeed  = 0.9999
  real(kind=dp)   , public         :: srmhd_maxspeed2 = 0.99998
  integer, public                  :: srmhd_maxiterationNR=100
  integer, public                  :: srmhd_maxiter_nr    =10000
  real(kind=dp)   , public         :: srmhd_limitvalue    =2.0d-8

  !> The small_est allowed energy
  real(kind=dp),public             :: small_e

  !> The small_est allowed pressure
  !real(kind=dp),public             :: small_pressure = 0.0_dp

  !> The small_est allowed inertia
  real(kind=dp),public             :: small_xi

  ! smaller values for speed
  real(kind=dp)   , public         :: small_vec2  = 0.0
  !> The number of waves
  integer                          :: nwwave=8

  !> Method type to clean divergence of B
  character(len=std_len), public     :: typedivbfix  = 'linde'

  !> Method type in a integer for good performance
  integer :: type_divb

  !> Coefficient of diffusive divB cleaning
  real(kind=dp)    :: divbdiff     = 0.8d0

  !> Update all equations due to divB cleaning
  character(len=std_len) ::    typedivbdiff = 'all'

  !> Use a compact way to add resistivity
  logical :: compactres   = .false.

  !> Add divB wave in Roe solver
  logical, public :: divbwave     = .true.

  !> Helium abundance over Hydrogen
  real(kind=dp)   , public      :: He_abundance=0.0d0!0.1d0

  !> To control divB=0 fix for boundary
  logical, public     :: boundary_divbfix(2*2)=.true.

  !> To skip * layer of ghost cells during divB=0 fix for boundary
  integer, public     :: boundary_divbfix_skip(2*2)=0

  !> B0 field is force-free
  logical, public     :: B0field_forcefree=.true.
  ! DivB cleaning methods
  integer, parameter :: divb_none          = 0
  integer, parameter :: divb_glm1          = 1
  integer, parameter :: divb_glm2          = 2
  integer, parameter :: divb_powel         = 3
  integer, parameter :: divb_janhunen      = 4
  integer, parameter :: divb_linde         = 5
  integer, parameter :: divb_lindejanhunen = 6
  integer, parameter :: divb_lindepowel    = 7
  integer, parameter :: divb_lindeglm      = 8

  logical,allocatable,target                     :: srmhd_iw_average(:)
  type (phys_variables_indices), target, public  :: srmhd_ind
  type(physconfig),target,public                 :: srmhd_config

  type funcd_srmhd
   procedure(gene_funcd),pointer,nopass :: procfuncd
  end type funcd_srmhd


abstract interface
  subroutine gene_funcd(xi,F,dF,lfac,d,ssqr,tau,taud,bsqr,sdotb,sb2,E2,ierror)
   use mod_global_parameters

   double precision, intent(in)  :: xi,d,ssqr,tau,taud,bsqr,sdotb,sb2,E2
   double precision, intent(out) :: F,dF,lfac
   integer, intent(inout)        :: ierror
   ! .. local ..
   double precision  :: dlfacdxi,rhoh,dhdxi,rho,drhodxi,xibsqr
   double precision  :: vsqr,p,dpdxi,invlfac

  end subroutine gene_funcd
end  interface
!===================================================================================
type phys_srmhd_parameters

    !> Whether thermal conduction is used
  logical, public                  :: thermal_conduction

  !> Whether radiative cooling is added
  logical, public                  :: radiative_cooling

   !> Whether sole energy equation
  logical, public                  :: energy
  !> Whether synge eos  is added
  logical, public                  :: eos

  !> Whether viscosity is added
  logical, public                  :: viscosity

  !> Whether gravity is added
  logical, public                  :: gravity

  !> Whether Hall-MHD is used
  logical, public                  :: Hall

  !> Whether particles module is added
  logical, public                  :: particles

  !> Whether magnetofriction is added
  logical, public                  :: magnetofriction

  !> Whether GLM-MHD is used
  logical, public                  :: glm

  !> dust not implimented in srmhd
  logical, public                  :: dust

  !> Whether divB cleaning sources are added splitting from fluid solver
  logical, public                  :: source_split_divb

  !> GLM-MHD parameter: ratio of the diffusive and advective time scales for div b
  !> taking values within [0, 1]
  real(kind=dp)   , public                :: glm_alpha

  !> MHD fourth order
  logical, public                  :: fourth_order

  !> Number of tracer species
  integer, public                  :: n_tracer

  !> Index of the moving frame density (in the w array)
  integer, public                  :: rho_

   !> Index of the lab frame density (in the w array)
  integer, public                  :: d_
  !> Indices of the momentum density
  integer, allocatable, public     :: mom(:)

  !> Index of the energy density (-1 if not present)
  integer, public                  :: e_

  !> Index of the gas pressure (-1 if not present) should equal e_
  integer, public                  :: p_


  !> Indices of the magnetic field
  integer, allocatable, public     :: mag(:)

  !> Indices of the GLM psi
  integer, public     :: psi_
  !> Indices of the Lorentz factor
  integer, public     :: lfac_

  !> Indices of the inertia
  integer, public     :: xi_



  !> Indices of the tracers
  integer, allocatable, public     :: tracer(:)

  !> The adiabatic index
  real(kind=dp)   , public                :: gamma

  !> The adiabatic constant
  real(kind=dp)   , public                :: adiab

  !> The MHD resistivity
  real(kind=dp)   , public                :: eta

  !> The MHD hyper-resistivity
  real(kind=dp)   , public                :: eta_hyper

  !> TODO: what is this?
  real(kind=dp)   , public                :: etah


  logical,          public         :: checkNR
  real(kind=dp)   , public         :: absaccnr
  real(kind=dp)   , public         :: tolerNR
  real(kind=dp)   , public         :: maxdspeed
  real(kind=dp)   , public         :: maxspeed
  real(kind=dp)   , public         :: maxspeed2
  integer, public                  :: maxiterationNR
  integer, public                  :: maxiter_nr
  real(kind=dp)   , public         :: limitvalue

  !> The small_est allowed energy
  real(kind=dp),public             :: small_e

  !> The small_est allowed pressure
  !real(kind=dp),public             :: small_pressure = 0.0_dp

  !> The small_est allowed inertia
  real(kind=dp),public             :: small_xi

  ! smaller values for speed
  real(kind=dp)   , public         :: small_vec2
  !> The number of waves
  integer                          :: nwwave=8

  !> Method type to clean divergence of B
  character(len=std_len), public     :: typedivbfix

  !> Method type in a integer for good performance
  integer :: type_divb

  !> Coefficient of diffusive divB cleaning
  real(kind=dp)    :: divbdiff

  !> Update all equations due to divB cleaning
  character(len=std_len) ::    typedivbdiff

  !> Use a compact way to add resistivity
  logical :: compactres

  !> Add divB wave in Roe solver
  logical, public :: divbwave

  !> Helium abundance over Hydrogen
  real(kind=dp)   , public      :: He_abundance!0.1d0

  !> To control divB=0 fix for boundary
  logical, public     :: boundary_divbfix(2*2)

  !> To skip * layer of ghost cells during divB=0 fix for boundary
  integer, public     :: boundary_divbfix_skip(2*2)

  !> B0 field is force-free
  logical, public     :: B0field_forcefree
  ! DivB cleaning methods
  integer :: divb_none
  integer :: divb_glm1
  integer :: divb_glm2
  integer :: divb_powel
  integer :: divb_janhunen
  integer :: divb_linde
  integer :: divb_lindejanhunen
  integer :: divb_lindepowel
  integer :: divb_lindeglm
end type phys_srmhd_parameters
end module mod_srmhd_parameters
