!> This module defines the procedures of a physics module. It contains function
!> pointers for the various supported routines. An actual physics module has to
!> set these pointers to its implementation of these routines.
module mod_physics
  use mod_global_parameters, only: name_len,ndir,ndim,std_len
  use mod_physics_hllc
  use mod_physics_roe
  use mod_physics_ppm
  use mod_constants

  implicit none
  public

  !> String describing the physics type of the simulation
  character(len=name_len) :: physics_type = ""

  !> To use wider stencils in flux calculations. A value of 1 will extend it by
  !> one cell in both directions, in any dimension
  integer :: phys_wider_stencil = 0

  !> Whether the physics routines require diagonal ghost cells, for example for
  !> computing a curl.
  logical :: phys_req_diagonal = .true.

  !> Array per direction per variable, which can be used to specify that certain
  !> fluxes have to be treated differently
  integer, allocatable :: flux_type(:, :)

  !> Indicates a normal flux
  integer, parameter   :: flux_default        = 0
  !> Indicates the flux should be treated with tvdlf
  integer, parameter   :: flux_tvdlf          = 1
  !> Indicates dissipation should be omitted
  integer, parameter   :: flux_no_dissipation = 2

  logical,pointer :: phys_iw_average(:)
  type phys_variables_indices
    integer              :: rho_
    integer,allocatable  :: mom(:)
    integer,allocatable  :: mag(:)
    integer              :: e_
    integer              :: pressure_

    integer              :: lfac_
    integer              :: xi_
    integer              :: psi_
    integer, allocatable :: tracer(:)

    integer, allocatable :: dust_rho(:)
    integer, allocatable :: dust_mom(:,:)
    integer, allocatable :: chemical_element(:)

    integer, allocatable :: epsinf_rho(:)
    integer, allocatable :: epsinf_rho0(:)
    integer, allocatable :: epsinf_Gmax(:)
    integer, allocatable :: epsinf_Gmaxadiab(:)
    integer              :: epsinf_B2turb_

    integer              :: Lcool1_
    integer              :: dtcool1_
    integer              :: mup_
  end type phys_variables_indices

  type(phys_variables_indices),pointer  :: phys_ind


  type physconfig
    logical                   :: ismhd
    logical                   :: isrel
    logical                   :: SI_unit

    integer                   :: nw
    integer                   :: nwfluxbc
    integer                   :: nwflux
    integer                   :: nwhllc

    logical                   :: dust_on
    integer                   :: dust_n_species
    real(dp)                  :: dust_small_density

    logical                   :: energy
    real(dp)                  :: gamma
    real(dp)                  :: adiab
    real(dp)                  :: temperature_isotherm
    logical                   :: eos
    logical                   :: thermal_conduction
    logical                   :: radiative_cooling
    logical                   :: viscosity
    logical                   :: particles
    character(len=30)         :: chemical_gas_type
    real(dp)                  :: He_abundance
    real(dp)                  :: mean_mass
    real(dp)                  :: mean_mup
    real(dp)                  :: mean_ne_to_nH
    real(dp)                  :: mean_nall_to_nH
    logical                   :: mean_mup_on

    logical                   :: cool_saveL
    logical                   :: cool_savedT
    logical                   :: cool_savemu

    real(dp)                  :: cool_tlow
    real(dp)                  :: cool_cfrac
    real(dp)                  :: cool_He_abundance
    integer                   :: cool_npoint
    !> Name of cooling curve
    character(len=std_len)    :: cool_curve

    !> Name of cooling method
    character(len=std_len)    :: cool_method

    !> Fixed temperature not lower than tlow
    logical                   :: cool_Tfix

    !> Add cooling source in a split way (.true.) or un-split way (.false.)
    logical                   :: cool_rc_split

    !> inonisation rate
    real(kind=dp)             :: cool_fn


    logical                   :: tracer_on
    integer                   :: n_tracer
    real(dp)                  :: unit_velocity
    real(dp)                  :: small_density
    real(dp)                  :: small_pressure
    real(dp)                  :: small_energy
    real(dp)                  :: small_xi
    real(dp)                  :: small_v2

    integer                   :: nwwave

    logical                   :: is4th_order
    logical                   :: compactres
    logical                   :: divbwave
    !> B0 field is force-free
    logical                   :: B0field_forcefree
    character(len=std_len)    :: typedivbfix
    integer                   :: type_divb
    character(len=std_len)    :: typedivbdiff
    logical                   :: source_split_divb
    !> To control divB=0 fix for boundary
    logical                   :: boundary_divbfix(2*ndim)
    !> To skip * layer of ghost cells during divB=0 fix for boundary
    integer                   :: boundary_divbfix_skip(2*ndim)

    logical                   :: magnetofriction
    logical                   :: hall
    real(kind=dp)             :: eta
    real(kind=dp)             :: eta_hyper
    real(kind=dp)             :: etah
    logical                   :: glm
    real(kind=dp)             :: glm_alpha


    integer                   :: maxiterationNR
    real(dp)                  :: absaccNR
    real(dp)                  :: tolerNr
    logical                   :: checkNR
    real(dp)                  :: maxdspeed
    real(dp)                  :: maxspeed
    real(dp)                  :: maxspeed2

    logical                   :: gravity
    logical                   :: chemical_on
    integer                   :: chemical_n_species
    real(dp)                  :: chemical_small_density

    logical                   :: epsinf_on
    integer                   :: epsinf_n_species
    real(dp)                  :: epsinf_small_Gammamax
    real(dp)                  :: epsinf_small_density
    real(dp)                  :: epsinf_small_Bfield_turb
  end type physconfig

  type(physconfig), pointer   :: phys_config

  procedure(sub_activate), pointer        :: phys_activate             => null()
  procedure(sub_check_params), pointer    :: phys_check_params         => null()
  procedure(sub_convert), pointer         :: phys_to_conserved         => null()
  procedure(sub_convert), pointer         :: phys_to_primitive         => null()
  procedure(sub_modify_wLR), pointer      :: phys_modify_wLR           => null()
  procedure(sub_get_cmax), pointer        :: phys_get_cmax             => null()
  procedure(sub_get_cbounds), pointer     :: phys_get_cbounds          => null()
  procedure(sub_get_flux), pointer        :: phys_get_flux             => null()
  procedure(sub_get_v_idim), pointer      :: phys_get_v_idim           => null()
  procedure(sub_get_dt), pointer          :: phys_get_dt               => null()
  procedure(sub_add_source_geom), pointer :: phys_add_source_geom      => null()
  procedure(sub_add_source), pointer      :: phys_add_source           => null()
  procedure(sub_get_aux), pointer         :: phys_get_aux              => null()
  procedure(sub_check_w), pointer         :: phys_check_w              => null()
  procedure(sub_get_pthermal), pointer    :: phys_get_pthermal         => null()
  procedure(sub_get_csound2), pointer     :: phys_get_csound2          => null()
  procedure(sub_get_comove_B2), pointer   :: phys_get_comove_B2        => null()
  procedure(sub_boundary_adjust), pointer :: phys_boundary_adjust      => null()
  procedure(sub_write_info), pointer      :: phys_write_info           => null()
  procedure(sub_angmomfix), pointer       :: phys_angmomfix            => null()
  procedure(sub_small_values), pointer    :: phys_handle_small_values  => null()
  procedure(sub_get_temperature), pointer :: phys_get_temperature      => null()
  procedure(sub_fill_chemical_ionisation), pointer :: phys_fill_chemical_ionisation => null()
  procedure(sub_get_4u_from_3v), pointer  :: phys_get_4u_from_3v       => null()
  procedure(sub_get_3v_from_4u), pointer  :: phys_get_3v_from_4u       => null()
  abstract interface


     subroutine sub_activate
     end subroutine sub_activate

     subroutine sub_check_params
     end subroutine sub_check_params

     subroutine sub_boundary_adjust
       use mod_global_parameters
     end subroutine sub_boundary_adjust

     subroutine sub_convert(ixI^L, ixO^L, w, x)
       use mod_global_parameters
       integer, intent(in)                :: ixI^L, ixO^L
       real(kind=dp)      , intent(inout) :: w(ixI^S, nw)
       real(kind=dp)      , intent(in)    :: x(ixI^S, 1:^ND)
     end subroutine sub_convert

     subroutine sub_modify_wLR(wLC, wRC, ixI^L, ixO^L, idir)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L, idir
       real(kind=dp)      , intent(inout) :: wLC(ixI^S,1:nw), wRC(ixI^S,1:nw)
     end subroutine sub_modify_wLR

     subroutine sub_get_cmax(w, x, ixI^L, ixO^L, idim, cmax)
       use mod_global_parameters
       integer, intent(in)                :: ixI^L, ixO^L, idim
       real(kind=dp)      , intent(in)    :: w(ixI^S, nw), x(ixI^S, 1:^ND)
       real(kind=dp)      , intent(inout) :: cmax(ixI^S)
     end subroutine sub_get_cmax

     subroutine sub_get_v_idim(w,x,ixI^L,ixO^L,idim,v)
       use mod_global_parameters

       integer, intent(in)              :: ixI^L, ixO^L, idim
       real(kind=dp)      , intent(in)  :: w(ixI^S,nw), x(ixI^S,1:^ND)
       real(kind=dp)      , intent(out) :: v(ixI^S)

     end subroutine sub_get_v_idim

     subroutine sub_get_cbounds(wLC, wRC, wLp, wRp, x, ixI^L, ixO^L, idim, cmax, cmin)
       use mod_global_parameters
       integer, intent(in)                :: ixI^L, ixO^L, idim
       real(kind=dp)      , intent(in)    :: wLC(ixI^S, nw), wRC(ixI^S, nw)
       real(kind=dp)      , intent(in)    :: wLp(ixI^S, nw), wRp(ixI^S, nw)
       real(kind=dp)      , intent(in)    :: x(ixI^S, 1:^ND)
       real(kind=dp)      , intent(inout) :: cmax(ixI^S)
       real(kind=dp)      , intent(inout), optional :: cmin(ixI^S)
     end subroutine sub_get_cbounds

     subroutine sub_get_flux(wC, w, x, ixI^L, ixO^L, idim, f)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L, idim
       real(kind=dp)      , intent(in)    :: wC(ixI^S, 1:nw)
       real(kind=dp)      , intent(in)    :: w(ixI^S, 1:nw)
       real(kind=dp)      , intent(in)    :: x(ixI^S, 1:^ND)
       real(kind=dp)      , intent(out)   :: f(ixI^S, nwflux)
     end subroutine sub_get_flux

     subroutine sub_add_source_geom(qdt, ixI^L, ixO^L, wCT, w, x)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L
       real(kind=dp)      , intent(in)    :: qdt, x(ixI^S, 1:^ND)
       real(kind=dp)      , intent(inout) :: wCT(ixI^S, 1:nw), w(ixI^S, 1:nw)
     end subroutine sub_add_source_geom

     subroutine sub_add_source(qdt, ixI^L, ixO^L, wCT, w, x, &
          qsourcesplit, active)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L
       real(kind=dp)      , intent(in)    :: qdt
       real(kind=dp)      , intent(in)    :: wCT(ixI^S, 1:nw), x(ixI^S, 1:ndim)
       real(kind=dp)      , intent(inout) :: w(ixI^S, 1:nw)
       logical, intent(in)             :: qsourcesplit
       logical, intent(inout)          :: active !< Needs to be set to true when active
     end subroutine sub_add_source

     subroutine sub_get_dt(w, ixI^L, ixO^L, dtnew, dx^D, x)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L
       real(kind=dp)      , intent(in)    :: dx^D, x(ixI^S, 1:^ND)
       real(kind=dp)      , intent(in)    :: w(ixI^S, 1:nw)
       real(kind=dp)      , intent(inout) :: dtnew
     end subroutine sub_get_dt

     subroutine sub_get_aux(clipping,w,x,ixI^L,ixO^L,subname,use_oldaux,sqrB,SdotB)
       use mod_global_parameters
       integer, intent(in)                :: ixI^L, ixO^L
       real(kind=dp)      , intent(in)    :: x(ixI^S,1:ndim)
       real(kind=dp)      , intent(inout) :: w(ixI^S,nw)
       logical, intent(in)                :: clipping
       character(len=*), intent(in)       :: subname
       logical         , intent(in), target, optional :: use_oldaux(ixI^S)
       real(kind=dp)   , intent(in), target, optional :: sqrB(ixI^S),SdotB(ixI^S)
     end subroutine sub_get_aux

     subroutine sub_check_w(primitive, ixI^L, ixO^L, w, w_flag)
       use mod_global_parameters
       logical, intent(in)          :: primitive
       integer, intent(in)          :: ixI^L, ixO^L
       real(kind=dp)      , intent(in) :: w(ixI^S,nw)
       integer, intent(inout)       :: w_flag(ixI^S)
     end subroutine sub_check_w

     subroutine sub_get_pthermal(w,x,ixI^L,ixO^L,pth)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L
       real(kind=dp)      , intent(in)    :: w(ixI^S,nw)
       real(kind=dp)      , intent(in)    :: x(ixI^S,1:ndim)
       real(kind=dp)      , intent(out)   :: pth(ixI^S)
     end subroutine sub_get_pthermal

     subroutine sub_get_csound2(w,x,ixI^L,ixO^L,csound2)
       use mod_global_parameters
       integer, intent(in)                :: ixI^L, ixO^L
       real(kind=dp)      , intent(in)    :: w(ixI^S,nw)
       real(kind=dp)      , intent(in)    :: x(ixI^S,1:ndim)
       real(kind=dp)      , intent(out)   :: csound2(ixI^S)
     end subroutine sub_get_csound2

     subroutine sub_get_comove_B2(ixI^L,ixO^L,x,w,sqrB)
       use mod_global_parameters
       integer, intent(in)                :: ixI^L, ixO^L
       real(kind=dp)      , intent(in)    :: w(ixI^S,nw)
       real(kind=dp)      , intent(in)    :: x(ixI^S,1:ndim)
       real(kind=dp)      , intent(out)   :: sqrB(ixO^S)
     end subroutine sub_get_comove_B2

     !> Calculate temperature within ixO^L
     subroutine sub_get_temperature( ixI^L, ixO^L,w, x, temperature)
       use mod_global_parameters

       integer, intent(in)          :: ixI^L, ixO^L
       real(dp), intent(in)         :: w(ixI^S, nw)
       real(dp), intent(in)         :: x(ixI^S, 1:ndim)
       real(dp), intent(out)        :: temperature(ixI^S)
       real(dp)                     :: pth(ixI^S)
       !----------------------------------------------------

     end subroutine sub_get_temperature

     !> set the chemical parameters
     subroutine sub_fill_chemical_ionisation(He_abundance_sub,gas_type, &
        mean_nall_to_nH,mean_mass,mean_mup,mean_ne_to_nH)
        use mod_global_parameters
        implicit none
        real(kind=dp), intent(in)    :: He_abundance_sub
        character(len=*), intent(in) :: gas_type
        real(kind=dp), intent(out)   :: mean_nall_to_nH,mean_mass,mean_mup,mean_ne_to_nH
        !----------------------------------------------------
     end   subroutine sub_fill_chemical_ionisation


     subroutine sub_write_info(file_handle)
       integer, intent(in) :: file_handle
     end subroutine sub_write_info

     subroutine sub_angmomfix(fC,x,wnew,ixI^L,ixO^L,idim)
       use mod_global_parameters
       integer, intent(in)                :: ixI^L, ixO^L
       real(kind=dp)      , intent(in)       :: x(ixI^S,1:ndim)
       real(kind=dp)      , intent(inout)    :: fC(ixI^S,1:nwflux,1:ndim), wnew(ixI^S,1:nw)
       integer, intent(in)                :: idim
     end subroutine sub_angmomfix

     subroutine sub_small_values(primitive, w, x, ixI^L, ixO^L, subname,flag_error)
       use mod_global_parameters
       use mod_small_values
       logical, intent(in)             :: primitive
       integer, intent(in)             :: ixI^L,ixO^L
       real(kind=dp)      , intent(inout) :: w(ixI^S,1:nw)
       real(kind=dp)      , intent(in)    :: x(ixI^S,1:ndim)
       character(len=*), intent(in)    :: subname
       integer, optional, intent(in)   :: flag_error(ixI^S)
     end subroutine sub_small_values

     subroutine sub_get_4u_from_3v(ixI^L,ixO^L,vtou,lfac)
       ! made by Z. Meliani 13/02/2018
       use mod_global_parameters
       implicit none
       integer, intent(in)             :: ixI^L,ixO^L
       real(kind=dp)   , intent(inout) :: vtou(ixI^S,1:ndir)
       real(kind=dp)   , intent(inout) :: lfac(ixI^S)
     end    subroutine sub_get_4u_from_3v

     subroutine sub_get_3v_from_4u(ixI^L,ixO^L,lfac,vtou)
  !   made by Z. Meliani 13/02/2018
      use mod_global_parameters
      implicit none
      integer, intent(in)             :: ixI^L,ixO^L
      real(kind=dp)   , intent(in)    :: lfac(ixI^S)
      real(kind=dp)   , intent(inout) :: vtou(ixI^S,1:ndir)
     end subroutine sub_get_3v_from_4u
  end interface

contains

  subroutine phys_default_config

        phys_config%ismhd                 = .false.
        phys_config%isrel                 = .false.
        phys_config%SI_unit               = .false.

        phys_config%dust_on               = .false.
        phys_config%dust_n_species        = 0
        phys_config%dust_small_density    = -1.0_dp

        phys_config%energy               = .true.
        phys_config%gamma                = 5.d0/3.0d0
        phys_config%adiab                = 1.0_dp
        phys_config%temperature_isotherm = -1.0_dp
        phys_config%eos                  = .false.
        phys_config%thermal_conduction   = .false.
        phys_config%radiative_cooling    = .false.
        phys_config%viscosity            = .false.
        phys_config%particles            = .false.
        phys_config%He_abundance         = 0.1_dp


        phys_config%cool_saveL           = .false.
        phys_config%cool_savedT          = .false.
        phys_config%cool_npoint          = 4000
        phys_config%cool_curve           = 'JCcorona'
        phys_config%cool_method          = 'exact'
        phys_config%cool_cfrac           = 0.1_dp
        phys_config%cool_tlow            = bigdouble
        phys_config%cool_Tfix            = .false.
        phys_config%cool_saveL           = .false.
        phys_config%cool_savedT          = .false.
        phys_config%cool_fn              = 1.0d-3
        phys_config%tracer_on            = .false.
        phys_config%n_tracer             = 0
        phys_config%unit_velocity        = 1.0_dp

        phys_config%small_density        = 0.0_dp
        phys_config%small_pressure       = 0.0_dp
        phys_config%small_energy         = 0.0_dp
        phys_config%small_v2             = 0.0_dp
        phys_config%small_xi             = 0.0_dp

        phys_config%nwwave               = 5
        phys_config%is4th_order          = .false.
        phys_config%compactres           = .false.
        phys_config%B0field_forcefree    = .false.
        phys_config%typedivbfix          = 'linde'
        phys_config%type_divb            = 1
        phys_config%typedivbdiff         = 'all'
        phys_config%source_split_divb    = .false.
        phys_config%divbwave             = .true.
        phys_config%glm                  = .false.

        phys_config%Hall                 = .false.
        phys_config%eta                  = 0.0_dp
        phys_config%eta_hyper            = 0.0_dp
        phys_config%etah                 = 0.0_dp
        phys_config%glm_alpha            = 0.5_dp

        phys_config%magnetofriction      = .false.

        phys_config%maxiterationNR       = 1000
        phys_config%absaccNR             = 1.0d-12
        phys_config%tolerNr              = 1.0d-9
        phys_config%checkNR              = .true.
        phys_config%maxdspeed            = 1.0d-9
        phys_config%maxspeed             = 1.0_dp-phys_config%maxdspeed
        phys_config%boundary_divbfix     = .true.
        phys_config%boundary_divbfix_skip= 0

        phys_config%chemical_on          = .false.
        phys_config%chemical_n_species   = 0
        phys_config%chemical_small_density = 0.0_dp

        phys_config%gravity              = .false.
  end subroutine phys_default_config



  subroutine phys_check()
    use mod_global_parameters, only: nw, ndir

    use mod_physics_hllc, only: phys_hllc_check
    use mod_physics_roe, only: phys_roe_check
    use mod_physics_ppm, only: phys_ppm_check

    if (physics_type == "") call mpistop("Error: no physics module loaded")

    call phys_hllc_check()
    call phys_roe_check()
    call phys_ppm_check()

    ! Checks whether the required physics methods have been defined
    if (.not. associated(phys_check_params)) &
         phys_check_params => dummy_check_params

    if (.not. associated(phys_to_conserved)) &
         call mpistop("Error: phys_to_conserved not defined")

    if (.not. associated(phys_to_primitive)) &
         call mpistop("Error: phys_to_primitive not defined")

    if (.not. associated(phys_modify_wLR)) &
         phys_modify_wLR => dummy_modify_wLR

    if (.not. associated(phys_get_cmax)) &
         call mpistop("Error: no phys_get_cmax not defined")

    if (.not. associated(phys_get_cbounds)) &
         call mpistop("Error: no phys_get_cbounds not defined")

    if (.not. associated(phys_get_flux)) &
         call mpistop("Error: no phys_get_flux not defined")

    if (.not. associated(phys_get_dt)) &
         call mpistop("Error: no phys_get_dt not defined")

    if (.not. associated(phys_add_source_geom)) &
         phys_add_source_geom => dummy_add_source_geom

    if (.not. associated(phys_add_source)) &
         phys_add_source => dummy_add_source

    if (.not. associated(phys_get_aux)) &
         phys_get_aux => dummy_get_aux

    if (.not. associated(phys_check_w)) &
         phys_check_w => dummy_check_w

    if (.not. associated(phys_get_pthermal)) &
         phys_get_pthermal => dummy_get_pthermal

    if (.not. associated(phys_get_csound2)) &
         phys_get_csound2=> dummy_get_csound2

    if (.not. associated(phys_get_comove_B2)) &
         phys_get_comove_B2 => dummy_get_comove_B2

    if (.not. associated(phys_boundary_adjust)) &
         phys_boundary_adjust => dummy_boundary_adjust

    if (.not. associated(phys_write_info)) &
         phys_write_info => dummy_write_info

    if (.not. associated(phys_angmomfix)) &
         phys_angmomfix => dummy_angmomfix

    if (.not. associated(phys_handle_small_values)) &
         phys_handle_small_values => dummy_small_values
    if (.not. associated(phys_get_temperature)) &
         phys_get_pthermal => dummy_get_temperature
    if (.not. associated(phys_get_4u_from_3v))   &
          phys_get_4u_from_3v=>dummy_get_4u_from_3v
    if (.not. associated(phys_get_3v_from_4u))   &
          phys_get_3v_from_4u=>dummy_get_3v_from_4u
  end subroutine phys_check

  subroutine dummy_init_params
  end subroutine dummy_init_params

  subroutine dummy_check_params
  end subroutine dummy_check_params

  subroutine dummy_modify_wLR(wLC, wRC, ixI^L, ixO^L, idir)
    use mod_global_parameters
    integer, intent(in)                :: ixI^L, ixO^L, idir
    real(kind=dp)      , intent(inout)    :: wLC(ixI^S,1:nw), wRC(ixI^S,1:nw)
  end subroutine dummy_modify_wLR

  subroutine dummy_add_source_geom(qdt, ixI^L, ixO^L, wCT, w, x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    real(kind=dp)      , intent(in)    :: qdt, x(ixI^S, 1:^ND)
    real(kind=dp)      , intent(inout) :: wCT(ixI^S, 1:nw), w(ixI^S, 1:nw)
  end subroutine dummy_add_source_geom

  subroutine dummy_add_source(qdt, ixI^L, ixO^L, wCT, w, x, &
       qsourcesplit, active)
    use mod_global_parameters
    integer, intent(in)                :: ixI^L, ixO^L
    real(kind=dp)      , intent(in)    :: qdt
    real(kind=dp)      , intent(in)    :: wCT(ixI^S, 1:nw), x(ixI^S, 1:ndim)
    real(kind=dp)      , intent(inout) :: w(ixI^S, 1:nw)
    logical, intent(in)                :: qsourcesplit
    logical, intent(inout)             :: active
    ! Don't have to set active, since it starts as .false.
  end subroutine dummy_add_source

  subroutine dummy_get_aux(clipping,w,x,ixI^L,ixO^L,subname,use_oldaux,sqrB,SdotB)
    use mod_global_parameters
    integer, intent(in)                :: ixI^L, ixO^L
    real(kind=dp)      , intent(in)    :: x(ixI^S,1:ndim)
    real(kind=dp)      , intent(inout) :: w(ixI^S,nw)
    logical, intent(in)                :: clipping
    character(len=*),intent(in)        :: subname
    logical         , intent(in), target, optional :: use_oldaux(ixI^S)
    real(kind=dp)   , intent(in), target, optional :: sqrB(ixI^S),SdotB(ixI^S)
  end subroutine dummy_get_aux

  subroutine dummy_check_w(primitive, ixI^L, ixO^L, w, w_flag)
    use mod_global_parameters
    logical, intent(in)          :: primitive
    integer, intent(in)          :: ixI^L, ixO^L
    real(kind=dp)   , intent(in) :: w(ixI^S,nw)
    integer, intent(inout)       :: w_flag(ixG^T)

    w_flag(ixO^S) = 0             ! All okay
  end subroutine dummy_check_w

  subroutine dummy_get_pthermal(w, x, ixI^L, ixO^L, pth)
    use mod_global_parameters

    integer, intent(in)                :: ixI^L, ixO^L
    real(kind=dp)      , intent(in)    :: w(ixI^S, nw)
    real(kind=dp)      , intent(in)    :: x(ixI^S, 1:ndim)
    real(kind=dp)      , intent(out)   :: pth(ixI^S)

    pth(ixO^S) = 0
    call mpistop("No get_pthermal method specified")
  end subroutine dummy_get_pthermal


  subroutine dummy_get_temperature(w, x, ixI^L, ixO^L, temperature)
    use mod_global_parameters

    integer, intent(in)                :: ixI^L, ixO^L
    real(kind=dp)      , intent(in)    :: w(ixI^S, nw)
    real(kind=dp)      , intent(in)    :: x(ixI^S, 1:ndim)
    real(kind=dp)      , intent(out)   :: temperature(ixI^S)

    temperature(ixO^S) = 0
    call mpistop("No get_ temperature method specified")
  end subroutine dummy_get_temperature
  subroutine dummy_get_csound2(w, x, ixI^L, ixO^L, csound2)
    use mod_global_parameters

    integer, intent(in)     :: ixI^L, ixO^L
    real(dp), intent(in)    :: w(ixI^S,nw)
    real(dp), intent(in)    :: x(ixI^S,1:ndim)
    real(dp), intent(out)   :: csound2(ixI^S)

    csound2(ixO^S) = 0
    call mpistop("No get_csound2 method specified")
  end subroutine dummy_get_csound2

  subroutine dummy_get_comove_B2( ixI^L, ixO^L, x,w,sqrB)
    use mod_global_parameters

    integer, intent(in)                :: ixI^L, ixO^L
    real(kind=dp)      , intent(in)    :: w(ixI^S, nw)
    real(kind=dp)      , intent(in)    :: x(ixI^S, 1:ndim)
    real(kind=dp)      , intent(out)   :: sqrB(ixO^S)

    sqrB(ixO^S) = 0
    call mpistop("No get_comove_B2 method specified")
  end subroutine dummy_get_comove_B2

  subroutine dummy_boundary_adjust
  end subroutine dummy_boundary_adjust

  subroutine dummy_write_info(fh)
    use mod_global_parameters
    integer, intent(in)                 :: fh !< File handle
    integer, dimension(MPI_STATUS_SIZE) :: st
    integer                             :: er

    ! Number of physics parameters
    integer, parameter                  :: n_par = 0

    call MPI_FILE_WRITE(fh, n_par, 1, MPI_INTEGER, st, er)
  end subroutine dummy_write_info

  subroutine dummy_angmomfix(fC,x,wnew,ixI^L,ixO^L,idim)
    use mod_global_parameters
    real(kind=dp)      , intent(in)       :: x(ixI^S,1:ndim)
    real(kind=dp)      , intent(inout)    :: fC(ixI^S,1:nwflux,1:ndim), wnew(ixI^S,1:nw)
    integer, intent(in)                :: ixI^L, ixO^L
    integer, intent(in)                :: idim
  end subroutine dummy_angmomfix

  subroutine dummy_small_values(primitive, w, x, ixI^L, ixO^L, subname,flag_error)
    use mod_global_parameters
    use mod_small_values
    logical, intent(in)                :: primitive
    integer, intent(in)                :: ixI^L,ixO^L
    real(kind=dp)      , intent(inout) :: w(ixI^S,1:nw)
    real(kind=dp)      , intent(in)    :: x(ixI^S,1:ndim)
    character(len=*), intent(in)       :: subname
    integer, optional, intent(in)      :: flag_error(ixI^S)
  end subroutine dummy_small_values

  subroutine dummy_get_4u_from_3v(ixI^L,ixO^L,vtou,lfac)
    ! made by Z. Meliani 13/02/2018
    use mod_global_parameters
    implicit none
    integer, intent(in)             :: ixI^L,ixO^L
    real(kind=dp)   , intent(inout) :: vtou(ixI^S,1:ndir)
    real(kind=dp)   , intent(inout) :: lfac(ixI^S)

    lfac(ixO^S)        = 1.0_dp
    vtou(ixO^S,1:ndir) = vtou(ixO^S,1:ndir)
  end    subroutine dummy_get_4u_from_3v

  subroutine dummy_get_3v_from_4u(ixI^L,ixO^L,lfac,vtou)
    ! made by Z. Meliani 13/02/2018
    use mod_global_parameters
    implicit none
    integer, intent(in)             :: ixI^L,ixO^L
    real(kind=dp)   , intent(in)    :: lfac(ixI^S)
    real(kind=dp)   , intent(inout) :: vtou(ixI^S,1:ndir)
     vtou(ixO^S,1:ndir) = vtou(ixO^S,1:ndir)
  end subroutine dummy_get_3v_from_4u
end module mod_physics
