!> Hydrodynamics physics module
module mod_hd_phys
  use mod_global_parameters, only: std_len
  use mod_constants
  use mod_physics
  implicit none
  private

  !> Whether an energy equation is used
  logical, public, protected              :: hd_energy = .true.

  !> Whether thermal conduction is added
  logical, public, protected              :: hd_thermal_conduction = .false.

  !> Whether radiative cooling is added
  logical, public, protected              :: hd_radiative_cooling = .false.

  !> Whether dust is added
  logical, public, protected              :: hd_dust = .false.

  !> Whether dust is added
  logical, public, protected              :: hd_chemical = .false.

  !> Whether viscosity is added
  logical, public, protected              :: hd_viscosity = .false.

  !> Whether grackle chemistry is added
  logical                   :: hd_use_grackle
  integer                   :: hd_gr_primordial_chemistry
  logical                   :: hd_gr_metal_cooling
  logical                   :: hd_gr_dust_chemistry
  logical                   :: hd_gr_with_radiative_cooling

  !> Whether gravity is added
  logical, public, protected              :: hd_gravity = .false.

  !> Whether gravity hydrostatic equilibrum improvement is enabled whereas gravity is enabled
  logical, public, protected              :: hd_gravity_hse = .false.
  character(len=30), public, protected    :: hd_gravity_hse_scheme = 'twoordernohse'
  logical, public, protected              :: hd_use_gravity_g = .false.
  integer, public, protected              :: hd_use_gravity_g_split = 0

  !> Whether particles module is added
  logical, public, protected              :: hd_particles = .false.

  !> Number of tracer species
  integer, public, protected              :: hd_n_tracer = 0

  real(dp), public                        :: hd_unit_velocity = 1.0_dp
  real(dp), public                        :: hd_unit_temperature = 1.0_dp

  !> Index of the density (in the w array)
  integer, public, protected              :: rho_

  !> Indices of the momentum density
  integer, allocatable, public, protected :: mom(:)

  !> Indices of grackle chemistry species
  integer, allocatable, public, protected :: gr_species(:)

  !> Indices of the tracers
  integer, allocatable, public, protected :: tracer(:)

  !> Index of the energy density (-1 if not present)
  integer, public, protected              :: e_

  !> Index of the gas pressure (-1 if not present) should equal e_
  integer, public, protected              :: p_

  integer, public, protected              :: gamma_
  integer, public, protected              :: temperature_

  !> Index of the gravity potential (-1 if not present)
  integer, public, protected              :: grav_phi_

 integer, public, protected              :: rhoHI_
  integer, public, protected              :: rhoHII_
  integer, public, protected              :: rhoHeI_
  integer, public, protected              :: rhoHeII_
  integer, public, protected              :: rhoHeIII_
  integer, public, protected              :: rhoe_
  integer, public, protected              :: rhoH2I_
  integer, public, protected              :: rhoHM_
  integer, public, protected              :: rhoH2II_
  integer, public, protected              :: rhoDI_
  integer, public, protected              :: rhoDII_
  integer, public, protected              :: rhoHDI_
  integer, public, protected              :: rhometal_
  integer, public, protected              :: rhodust_


  !> The adiabatic index
  real(dp), public                :: hd_gamma = 5.d0/3.0d0

  !> The adiabatic constant
  real(dp), public                :: hd_adiab = 1.0_dp

  logical, public                :: hd_isotherm_on = .false.
  !> The isotherm temperature constant
  real(dp), public                :: hd_temperature_isotherm = 1.0_dp

  !> The small_est allowed energy
  real(dp), protected             :: small_e
  !> The small_est allowed density
  real(dp), public, protected     :: hd_small_density = 0.0_dp
  !> The small_est allowed pressure
  real(dp), public, protected     :: hd_small_pressure= 0.0_dp

  !> Helium abundance over Hydrogen
  real(dp), public, protected  :: He_abundance=0.1d0

  !> Type of gas
  character(len=30), public, protected :: hd_chemical_gas_type='fullyionised'
  !> The number of waves
  integer :: nwwave = 4
  !> use mean mup to calculate the temperature
  logical           :: hd_mean_mup_on
  logical, allocatable,target    :: hd_iw_average(:)

  type(physconfig),target,public                 :: hd_config
  type(physconfig),target,public                 :: pre_hd_config
  type(phys_variables_indices),target,public     :: hd_ind
  ! for local use
  real(kind=dp), private         :: gamma_1
  ! Public methods
  public :: hd_phys_init
  public :: hd_kin_en
  public :: hd_get_pthermal
  public :: hd_get_temperature
  public :: hd_get_gamma
  public :: hd_get_mmw
  public :: hd_to_conserved
  public :: hd_to_primitive
  public :: hd_small_values_floor
contains



  !> Read this module s parameters from a file
  subroutine hd_read_params(files)
    use mod_global_parameters
    character(len=*), intent(in) :: files(:)
    integer                      :: i_file,i_reason
    character(len=70)            :: error_message

    namelist /hd_list/ hd_energy, hd_n_tracer, hd_unit_velocity, hd_unit_temperature, &
    hd_gamma, hd_adiab, &
    hd_dust, hd_thermal_conduction, hd_radiative_cooling, hd_viscosity, &
    hd_gravity, hd_gravity_hse, hd_gravity_hse_scheme,hd_use_gravity_g, &
    He_abundance, SI_unit, hd_particles,hd_small_density,hd_small_pressure,&
    hd_chemical,hd_chemical_gas_type,hd_mean_mup_on,hd_temperature_isotherm,&
    hd_use_grackle, hd_gr_primordial_chemistry, hd_gr_metal_cooling,&
    hd_gr_with_radiative_cooling,&
    hd_gr_dust_chemistry,hd_isotherm_on

    !---------------------------------------------------------
    error_message = 'At '//' mod_hd_phys.t'//'  in the procedure : hd_read_params'
    Loop_iparfile : do i_file = 1, size(files)
      open(unitpar, file=trim(files(i_file)), status="old")
      read(unitpar, hd_list, iostat=i_reason)
      cond_ierror : if(i_reason>0)then
       write(*,*)' Error in reading the parameters file : ',trim(files(i_file))
       write(*,*)' Error at namelist: ', 'hd_list'
       write(*,*)' Error number = ',i_reason
       write(*,*)' The code stops now '
       call mpistop(trim(error_message))
      elseif(i_reason<0)then cond_ierror
       write(*,*)' Reache the end of the file  : ',trim(files(i_file))
       write(*,*)' Error at namelist: hd_list'
       write(*,*)' Error number = ',i_reason
       write(*,*)' The code stops now '
       call mpistop(trim(error_message))
      else cond_ierror
       write(*,*)' End of reading of the hd_list'
      end if cond_ierror
      close(unitpar)

    end do Loop_iparfile


  end subroutine hd_read_params





  !> Write this module s parameters to a snapshot
  subroutine hd_write_info(fh)
    use mod_global_parameters
    integer, intent(in)                 :: fh
    integer, parameter                  :: n_par = 1
    real(dp)                    :: values(n_par)
    character(len=name_len)             :: names(n_par)
    integer, dimension(MPI_STATUS_SIZE) :: st
    integer                             :: er

    call MPI_FILE_WRITE(fh, n_par, 1, MPI_INTEGER, st, er)

    names(1) = "gamma"
    values(1) = hd_gamma
    call MPI_FILE_WRITE(fh, values, n_par, MPI_DOUBLE_PRECISION, st, er)
    call MPI_FILE_WRITE(fh, names, n_par * name_len, MPI_CHARACTER, st, er)
  end subroutine hd_write_info

  !> Add fluxes in an angular momentum conserving way
  subroutine hd_angmomfix(fC,x,wnew,ixI^L,ixO^L,idim)
    use mod_global_parameters
    real(dp), intent(in)       :: x(ixI^S,1:ndim)
    real(dp), intent(inout)    :: fC(ixI^S,1:nwflux,1:ndim),  wnew(ixI^S,1:nw)
    integer, intent(in)                :: ixI^L, ixO^L
    integer, intent(in)                :: idim
    integer                            :: hxO^L, kxC^L, iw
    real(dp)                           :: inv_volume(ixI^S)

    ! shifted indexes
    hxO^L=ixO^L-kr(idim,^D);
    ! all the indexes
    kxCmin^D=hxOmin^D;
    kxCmax^D=ixOmax^D;

    inv_volume = 1.0_dp/block%dvolume(ixO^S)

    select case(typeaxial)
    case ("cylindrical")
      do iw=1,nwflux
        if (idim==r_ .and. iw==iw_mom(phi_)) then
          fC(kxC^S,iw,idim)= fC(kxC^S,iw,idim)*(x(kxC^S,r_)+half*block%dx(kxC^S,idim))
          wnew(ixO^S,iw)=wnew(ixO^S,iw) + (fC(ixO^S,iw,idim)-fC(hxO^S,iw,idim)) * &
               (inv_volume/x(ixO^S,idim))
        else
          wnew(ixO^S,iw)=wnew(ixO^S,iw) + (fC(ixO^S,iw,idim)-fC(hxO^S,iw,idim)) * &
                inv_volume
        endif
      enddo
    case ("spherical")
      do iw=1,nwflux
        if     (idim==r_ .and. (iw==iw_mom(2) .or. iw==iw_mom(phi_))) then
          fC(kxC^S,iw,idim)= fC(kxC^S,iw,idim)*(x(kxC^S,idim)+half*block%dx(kxC^S,idim))
          wnew(ixO^S,iw)=wnew(ixO^S,iw) + (fC(ixO^S,iw,idim)-fC(hxO^S,iw,idim)) * &
               (inv_volume/x(ixO^S,idim))
        elseif (idim==2  .and. iw==iw_mom(phi_)) then
          fC(kxC^S,iw,idim)=fC(kxC^S,iw,idim)*sin(x(kxC^S,idim)+half*block%dx(kxC^S,idim)) ! (x(4,3,1)-x(3,3,1)))
          wnew(ixO^S,iw)=wnew(ixO^S,iw) + (fC(ixO^S,iw,idim)-fC(hxO^S,iw,idim)) * &
               (inv_volume/sin(x(ixO^S,idim)))
        else
          wnew(ixO^S,iw)=wnew(ixO^S,iw) + (fC(ixO^S,iw,idim)-fC(hxO^S,iw,idim)) * &
                inv_volume
        endif
      enddo

    end select

  end subroutine hd_angmomfix


  subroutine hd_phys_set_config
    use mod_global_parameters
    use mod_physics
    implicit none


    use_particles                  = hd_particles

    hd_config%dust_on              = hd_dust
    hd_config%energy               = hd_energy
    hd_config%gravity              = hd_gravity
    !> Whether gravity hydrostatic equilibrum improvement is enabled whereas gravity is enabled
    hd_config%gravity_hse          = hd_gravity_hse
    hd_config%gravity_hse_scheme   = hd_gravity_hse_scheme
    hd_config%use_gravity_g        = hd_use_gravity_g
    hd_config%chemical_on          = hd_chemical
    hd_config%chemical_gas_type    = hd_chemical_gas_type
    hd_config%He_abundance         = He_abundance
    hd_config%mean_mup_on          = hd_mean_mup_on
    hd_config%dust_n_species       = 0
    hd_config%ismhd                = .false.
    hd_config%isrel                = .false.
    hd_config%n_tracer             = hd_n_tracer
    hd_config%unit_velocity        = hd_unit_velocity
    hd_config%unit_temperature     = hd_unit_temperature
    hd_config%gamma                = hd_gamma
    hd_config%adiab                = hd_adiab
    hd_config%isotherm_on          = hd_isotherm_on
    hd_config%temperature_isotherm = hd_temperature_isotherm

    gamma_1 = hd_config%gamma   - 1.0_dp
    hd_config%small_density        = small_density
    hd_config%small_pressure       = small_pressure
    if(hd_config%energy) then
      hd_config%small_energy       = hd_config%small_pressure/gamma_1
    end if
    hd_config%radiative_cooling    = hd_radiative_cooling
    hd_config%use_grackle          = hd_use_grackle
    hd_config%gr_primordial_chemistry = hd_gr_primordial_chemistry
    hd_config%gr_metal_cooling     = hd_gr_metal_cooling
    hd_config%gr_dust_chemistry    = hd_gr_dust_chemistry
    hd_config%gr_with_radiative_cooling = hd_gr_with_radiative_cooling
    hd_config%thermal_conduction   = hd_thermal_conduction
    hd_config%viscosity            = hd_viscosity
    hd_config%particles            = hd_particles

    phys_energy                    = hd_config%energy

    if(hd_config%isotherm_on)then
      if(hd_config%energy)then
        write(*,*) 'hd_energy=T & hd_isotherm_on=T,...'
        call mpistop('.. but isothermal case incompatible with energy equation')
      end if
    end if
    if(hd_config%isotherm_on)then
      if(dabs(hd_config%gamma-1.0_dp)>smalldouble)then
        call mpistop('Isothermal gas but hd_gamma/=1')
      end if
    end if

  end subroutine hd_phys_set_config



  subroutine hd_phys_to_phys

    use mod_physics

    phys_config     => hd_config
    phys_ind        => hd_ind
    !allocate(hd_iw_average(1:nw))
    !hd_iw_average = .true.

    phys_get_dt                   => hd_get_dt
    phys_get_cmax                 => hd_get_cmax
    phys_get_aux                 => hd_get_aux
    phys_get_cbounds              => hd_get_cbounds
    phys_get_flux                 => hd_get_flux
    phys_get_v_idim               => hd_get_v
    phys_add_source_geom          => hd_add_source_geom
    phys_add_source               => hd_add_source
    phys_to_conserved             => hd_to_conserved
    phys_to_primitive             => hd_to_primitive
    phys_check_params             => hd_check_params
    phys_check_w                  => hd_check_w
    phys_get_pthermal             => hd_get_pthermal
    phys_get_gamma                => hd_get_gamma
    phys_get_mup                  => hd_get_mmw
    phys_get_kin_en               => hd_kin_en
    phys_get_temperature          => hd_get_temperature
    phys_get_csound2              => hd_get_csound2
    phys_write_info               => hd_write_info
    phys_handle_small_values      => hd_handle_small_values
    phys_angmomfix                => hd_angmomfix
    phys_fill_chemical_ionisation =>hd_fill_chemical_ionisation
    ! Whether diagonal ghost cells are required for the physics
    phys_req_diagonal = .true. !.false.
  end subroutine hd_phys_to_phys

  !> Initialize the module
  subroutine hd_phys_init()
    use mod_global_parameters
    use mod_thermal_conduction
    use mod_radiative_cooling
    use mod_chemical, only: chemical_init
    use mod_dust, only: dust_init,dust_n_species
    use mod_viscosity, only: viscosity_init
    use mod_gravity, only: gravity_init
    use mod_particles, only: particles_init
    use mod_physics

    integer :: itr, idir
    real(kind=dp)  :: mp,kB
    !------------------------------------------------------
    physics_type = "hd"
    call phys_default_config(hd_config)
    call hd_read_params(par_files)
    call hd_phys_set_config()
    call hd_fill_phys_indices()

    call hd_phys_to_phys()




    ! set the mean molecular mass and electron density
    call hd_fill_chemical_ionisation(hd_config%He_abundance,hd_config%chemical_gas_type, &
       hd_config%mean_nall_to_nH,hd_config%mean_mass,&
      hd_config%mean_mup,hd_config%mean_ne_to_nH)

    ! derive units from basic units
    call hd_physical_units
    !isotherm case
    if(hd_config%isotherm_on)then
      if(SI_unit) then
          mp=mp_SI
          kB=kB_SI
      else
          mp=mp_cgs
          kB=kB_cgs
      end if

      if(hd_config%temperature_isotherm>0)then
        ! isotherm case : need to compute adiab = csound**2
        hd_config%adiab = kB/(hd_config%mean_mup*mp)*hd_config%temperature_isotherm
      elseif(hd_config%adiab>=0)then
        ! isotherm case : need to compute isotherm temperature T
        hd_config%temperature_isotherm = (hd_config%adiab/kB)*hd_config%mean_mup*mp
      end if

      ! normalise value for isotherm temperature
      hd_config%temperature_isotherm=hd_config%temperature_isotherm/unit_temperature
      ! normalise value for adiab
      hd_config%adiab  = hd_config%adiab/(unit_velocity*unit_velocity)

      ! save in local hd paramteres
      hd_adiab= hd_config%adiab
      hd_temperature_isotherm = hd_config%temperature_isotherm

      if(hd_config%energy)then
        write(*,*) 'hd_config%isotherm_on = T but hd_config%energy = T'
        call mpistop('Isothermal gas : must have hd_config%isotherm_on=T,hd_config%energy=F')
      end if

      if(dabs(pre_hd_config%gamma-1.0_dp)>smalldouble)then
        call mpistop('Isothermal gas but hd_gamma/=1')
      else
        write(*,*) '*********************************************************'
        write(*,*) 'MPI-AMRVAC runs HD adiabatic isothermal run (no energy solved)'
        write(*,*) 'with temperature = ',hd_config%temperature_isotherm*unit_temperature
        write(*,*) 'and c_sound**2 = ',hd_config%adiab*(unit_velocity*unit_velocity)
        write(*,*) '*********************************************************'
      end if

    else
      if(.not.hd_config%energy)then
        write(*,*) 'hd_config%isotherm_on = F but hd_config%energy = F'
        call mpistop('Isothermal gas : must have hd_config%isotherm_on=T,hd_config%energy=F')
      else
        write(*,*) '*********************************************************'
        write(*,*) 'MPI-AMRVAC runs HD run (gamma/=1) with energy solved'
        write(*,*) 'temperature can be variable, and we have'
        write(*,*) 'gamma =',hd_config%gamma
        write(*,*) '*********************************************************'
      end if
    end if

    if(hd_config%adiab<0)then
      write(*,*) " mod_usr.t : non-physical hd_adiab < 0 set "
      call mpistop(" .par parameter file : hd_list parameters uncorrect")
    end if

    if (hd_config%dust_on) call dust_init(hd_ind,hd_config,rho_, mom(:), e_)

    if(hd_config%chemical_on)then
      call chemical_init(hd_ind,hd_config,rho_, mom(:), e_)
    end if
    call hd_fill_convert_factor


    ! initialize thermal conduction module
    if (hd_config%thermal_conduction) then
      if (.not. hd_config%energy) &
           call mpistop("thermal conduction needs hd_config%energy=T")
      call thermal_conduction_init(hd_config%gamma)
    end if

    ! Initialize radiative cooling module
    if (hd_config%radiative_cooling) then
      if (.not. hd_config%energy) call mpistop("radiative cooling needs hd_config%energy=T")
      call radiative_cooling_init(hd_ind,hd_config,hd_config%gamma,He_abundance)
    end if

    ! Initialize viscosity module
    if (hd_config%viscosity) call viscosity_init(phys_wider_stencil,phys_req_diagonal)

    ! Initialize gravity module
    if (hd_config%gravity) call gravity_init()

    ! Initialize particles module
    if (hd_config%particles) then
       call particles_init()
       phys_req_diagonal = .true.
    end if

    ! Check whether custom flux types have been defined
    if (.not. allocated(flux_type)) then
       allocate(flux_type(ndir, nw))
       flux_type = flux_default
    else if (any(shape(flux_type) /= [ndir, nw])) then
       call mpistop("phys_check error: flux_type has wrong shape")
    end if



    nvector      = 1 ! No. vector vars
    allocate(iw_vector(nvector))
    iw_vector(1) = mom(1) - 1

    allocate(hd_iw_average(1:nw))
    hd_iw_average = .true.

    !automatically_
    if(hd_config%gravity)then
      if(hd_config%gravity_hse)then
            phys_modify_wLR => improve_gravity_hse
      end if
    end if


    phys_iw_average => hd_iw_average

    hd_config%energy = hd_config%energy

    phys_iw_average => hd_iw_average
    phys_config     => hd_config
    !phys_ind        => hd_ind
    ! refill again to add
    call hd_fill_convert_factor
  end subroutine hd_phys_init

  subroutine improve_gravity_hse(wLp, wRp, wCT, x, ixI^L, ixO^L, idir,&
  qdt, qtC, qt, dx^D, method, wnew, wold, local_block, wLC, wRC)
    use mod_global_parameters
    use mod_usr_methods, only: usr_reset_solver,usr_gravity_potential
    use mod_finite_volume, only: reconstruct_LR, reconstruct_LR3
    character(len=*), intent(in)                         :: method
    integer, intent(in)             :: ixI^L, ixO^L, idir
    real(kind=dp)   , dimension(ixI^S,1:ndim), intent(in) ::  x
    real(kind=dp)   , intent(in)                         :: qdt, qtC, qt, dx^D
    real(kind=dp)      , intent(inout) :: wLp(ixI^S,1:nw), wRp(ixI^S,1:nw), wCT(ixI^S,1:nw)
    type(walloc)       , intent(inout)    :: local_block
    real(kind=dp)   , intent(inout) :: wLC(ixI^S,1:nw), wRC(ixI^S,1:nw)
    real(kind=dp)   , dimension(ixI^S,1:nw)               :: wnew, wold
    ! ..........local...........
    ! primitive w at cell center
    real(kind=dp)   , dimension(ixI^S,1:nw) :: wprim,wpert,wpertLp,wpertRp
    real(kind=dp)   , dimension(ixI^S,1:nw) :: wpertLC, wpertRC
    integer                                 :: ixOfinite^L
    real(kind=dp)   , dimension(1:ndim)     :: dxinv
    character(len=len(method))              :: method_loc
    integer :: iw, ix^L, hxO^L, ixC^L, ixCR^L, kxL^L, kxC^L, kxR^L
    real(kind=dp)   , dimension(ixI^S,1:ndim) ::  xhalfR,xhalfL
    real(kind=dp)   , dimension(ixI^S) ::  gravity_potential
    real(kind=dp)   , dimension(ixI^S) ::  gravity_potentialR, gravity_potentialL

    !Recompute output part of the indexes from mod_finite_volume.t:finite_volume
    ixOfinitemax^D=ixOmax^D;
    ixOfinitemin^D=ixOmin^D+kr(idir,^D);



    wprim=wCT
    call phys_to_primitive(ixI^L,ixI^L,wprim,x)


    ^D&dxinv(^D)=-qdt/dx^D;

    ! initialise the flux scheme
    method_loc=trim(method)

    !allow the user to change localy the numerical method
    !Mialy : this should be at the end of the finite volume treatment
    !after cooling source term and then L1 is being properly computed
    !so that solver is changed for the next time step
    call usr_reset_solver(ixI^L,ixI^L,idir,qt,wprim,x,method,&
      type_limiter(node(plevel_,saveigrid)),method_loc,typelimiter)


    ! use interface value of w0 at idir=finite_volume:idim
    local_block%iw0=idir

    hxO^L=ixOfinite^L-kr(idir,^D);
    ! ixC is centered index in the idir=finite_volume:idim
    ! direction from ixOmin-1/2 to ixOmax+1/2
    ixCmax^D=ixOfinitemax^D; ixCmin^D=hxOmin^D;

    kxCmin^D=ixImin^D; kxCmax^D=ixImax^D-kr(idir,^D);!i or j
    kxR^L=kxC^L+kr(idir,^D);! i+1 or j+1
    kxL^L=kxC^L-kr(idir,^D);! i-1 or j-1



    ! wRp and wLp are defined at the same locations, and will correspond to
    ! the left and right reconstructed values at a cell face. Their indexing
    ! is similar to cell-centered values, but in direction idir=finite_volume:idim they are
    ! shifted half a cell towards the 'lower' direction.
    !reconstruction or not of the left and right state variables
    !at the cell interfaces and their conservative variables
     wRp(kxC^S,1:nw)=wprim(kxR^S,1:nw)
     wLp(kxC^S,1:nw)=wprim(kxC^S,1:nw)
     !Recompute pressure for hydrostatic equilibrum:
     xhalfR(ixI^S,1:ndim)=x(ixI^S,1:ndim)
     xhalfL(ixI^S,1:ndim)=x(ixI^S,1:ndim)
     xhalfR(ixI^S,idir)=x(ixI^S,idir)-half*local_block%dx(ixI^S,idir)
     xhalfL(ixI^S,idir)=x(ixI^S,idir)+half*local_block%dx(ixI^S,idir)

     !where(DABS(xhalfR(ixI^S,r_))<=smalldouble)
       !xhalfR(ixI^S,r_)=x(ixI^S,idir)-half*local_block%dx(ixI^S,idir)+smalldouble
     !end where
     !where(DABS(xhalfL(ixI^S,r_))<=smalldouble)
       !xhalfL(ixI^S,r_)=x(ixI^S,idir)+half*local_block%dx(ixI^S,idir)-smalldouble
     !end where

     call usr_gravity_potential(ixI^L,ixO^L,wCT,xhalfR,gravity_potentialR)
     call usr_gravity_potential(ixI^L,ixO^L,wCT,x,gravity_potential)
     call usr_gravity_potential(ixI^L,ixO^L,wCT,xhalfL,gravity_potentialL)


     wRp(kxC^S,iw_e)=wprim(kxR^S,iw_e)-wprim(kxR^S,iw_rho)*&
                     (gravity_potentialR(kxR^S)-gravity_potential(kxR^S))
     wLp(kxC^S,iw_e)=wprim(kxC^S,iw_e)-wprim(kxC^S,iw_rho)*&
                     (gravity_potentialL(kxC^S)-gravity_potential(kxC^S))


    ! Determine stencil size
    ! phys_wider_stencil from mod_physics and set in each physics mod_..._phys.t module
    {ixCRmin^D = ixCmin^D - phys_wider_stencil\}
    {ixCRmax^D = ixCmax^D + phys_wider_stencil\}



    ! apply limited reconstruction for left and right status at cell interfaces
    select case (typelimited)
    case ('previous')
    call reconstruct_LR(ixI^L,ixCR^L,ixCR^L,idir,wold,wprim,wLC,wRC,wLp,wRp,x,.true.)
    case ('predictor')
    call reconstruct_LR(ixI^L,ixCR^L,ixCR^L,idir,wprim,wprim,wLC,wRC,wLp,wRp,x,.false.)
    case default
    call mpistop("Error in reconstruction: no such base for limiter")
    end select

    !write(*,*) 'phys_modify_wLR called until end !!'

  end subroutine improve_gravity_hse

  subroutine improve_gravity_hse_exact(wLp, wRp, wCT, x, ixI^L, ixO^L, idir,&
  qdt, qtC, qt, dx^D, method, wnew, wold, local_block, wLC, wRC)
    use mod_global_parameters
    use mod_usr_methods, only: usr_reset_solver, usr_gravity_fpotential
    use mod_finite_volume, only: reconstruct_LR, reconstruct_LR3
    character(len=*), intent(in)                         :: method
    integer, intent(in)             :: ixI^L, ixO^L, idir
    real(kind=dp)   , dimension(ixI^S,1:ndim), intent(in) ::  x
    real(kind=dp)   , intent(in)                         :: qdt, qtC, qt, dx^D
    real(kind=dp)      , intent(inout) :: wLp(ixI^S,1:nw), wRp(ixI^S,1:nw), wCT(ixI^S,1:nw)
    type(walloc)       , intent(inout)    :: local_block
    real(kind=dp)   , intent(inout) :: wLC(ixI^S,1:nw), wRC(ixI^S,1:nw)
    real(kind=dp)   , dimension(ixI^S,1:nw)               :: wnew, wold
    ! ..........local...........
    ! primitive w at cell center
    real(kind=dp)   , dimension(ixI^S,1:nw) :: wprim,wpert,wpertLp,wpertRp
    real(kind=dp)   , dimension(ixI^S,1:nw) :: wpertLC, wpertRC
    real(kind=dp)   , dimension(ixI^S)      :: wphikR,wphikC,wphikL
    real(kind=dp)   , dimension(ixI^S,1:nw) :: xhalfR,xhalfL
    integer                                 :: ixOfinite^L
    real(kind=dp)   , dimension(1:ndim)     :: dxinv
    character(len=len(method))              :: method_loc
    integer :: iw, ix^L, hxO^L, ixC^L, ixCR^L, kxL^L, kxC^L, kxR^L, myidim



    write(*,*) 'improve_gravity_hse_exact in mod_hd_phys.t:'
    write(*,*) 'use_gravity_g is true but'
    write(*,*) 'but this part is not yet implemented'
    call mpistop('Code part not yet implemented')


    !write(*,*) 'phys_modify_wLR called until end !!'




  end subroutine improve_gravity_hse_exact


  subroutine hd_fill_chemical_ionisation(He_abundance_sub,hd_chemical_gas_type,mean_nall_to_nH &
                                        ,mean_mass,mean_mup,mean_ne_to_nH)
    use mod_global_parameters
    implicit none
    real(kind=dp), intent(in)    :: He_abundance_sub
    character(len=*), intent(in) :: hd_chemical_gas_type
    real(kind=dp), intent(out)   :: mean_nall_to_nH,mean_mass,mean_mup,mean_ne_to_nH
    !----------------------------------------------------------
    if(He_abundance_sub>0.0_dp)then
      mean_nall_to_nH = (1.0_dp+2.0_dp*He_abundance_sub)
      mean_mass       = (1.0_dp+4.0_dp*He_abundance_sub)
      select case(trim(hd_chemical_gas_type))
        case('fullymolecular')
          mean_mup = mean_mass/(0.5_dp+He_abundance_sub)
          mean_ne_to_nH =0.0_dp  ! is the value used in implimentation of low temperarure cooling table DM2
        case('fullyatomic')
          mean_mup = mean_mass/(1.0_dp+He_abundance_sub)
          mean_ne_to_nH =0.0_dp ! is the value used in implimentation of low temperarure cooling table DM2
        case('ionised')
          mean_mup = mean_mass/(2.0_dp+2.0_dp*He_abundance_sub)
          mean_ne_to_nH = (1.0_dp+He_abundance_sub)
        case('fullyionised')
          mean_mup = mean_mass/(2.0_dp+3.0_dp*He_abundance_sub)
          mean_ne_to_nH = (1.0_dp+2.0_dp*He_abundance_sub)
       case default
         write(*,*) 'The chemical gas type : ', trim(hd_chemical_gas_type)
          write (*,*) "Undefined gas chemical type entered in mod_hd_phys.t "
          call mpistop('The stops at hd_fill_chemical_ionisation in src/hd/mod_hd_phys.t')
      end select
      !hd_config%mean_mup          = (2.0_dp+3.0_dp*He_abundance_sub)

    else if(dabs(He_abundance_sub)<=smalldouble)then
      mean_mass         = 1.0_dp
      mean_mup          = 1.0_dp
      mean_ne_to_nH     = 1.0_dp
    end if
  end   subroutine hd_fill_chemical_ionisation

  subroutine hd_fill_phys_indices
    use mod_global_parameters
    implicit none
    ! .. local ..
    integer :: itr, idir
    !-------------------------------------

    ! Determine flux variables
    rho_ = var_set_rho()
    allocate(mom(ndir))
    mom(:) = var_set_momentum(ndir)

    ! Set index of energy variable
    if (hd_config%energy) then
      nwwave = 4
      e_     = var_set_energy() ! energy density
      p_     = e_               ! gas pressure
    else
      nwwave = 3
      e_     = -1
      p_     = -1
    end if


    !>Grackle:
    iw_HI_density = -1
    iw_HII_density = -1
    iw_HeI_density = -1
    iw_HeII_density = -1
    iw_HeIII_density = -1
    iw_e_density = -1
    iw_H2I_density = -1
    iw_HM_density = -1
    iw_H2II_density = -1
    iw_DI_density = -1
    iw_DII_density = -1
    iw_HDI_density = -1
    iw_metal_density = -1
    iw_dust_density = -1

    using_grackle : if(hd_config%use_grackle)then
      prim_chem : if(hd_config%gr_primordial_chemistry>=0.and.hd_config%gr_primordial_chemistry<=3)then
        !six species: H, H+, He, He+, He++, e- (in this order)
        iw_HI_density = var_set_fluxvar('HI', 'HI', need_bc=.true.)
        iw_HII_density = var_set_fluxvar('HII', 'HII', need_bc=.true.)
        iw_HeI_density = var_set_fluxvar('HeI', 'HeI', need_bc=.true.)
        iw_HeII_density = var_set_fluxvar('HeII', 'HeII', need_bc=.true.)
        iw_HeIII_density = var_set_fluxvar('HeIII', 'HeIII', need_bc=.true.)
        iw_e_density = var_set_fluxvar('electron_density', 'electron_density', need_bc=.true.)
        ! + 3 species : H2, H-, H2+ (in this order)
        iw_H2I_density = var_set_fluxvar('H2I', 'H2I', need_bc=.true.)
        iw_HM_density = var_set_fluxvar('HM', 'HM', need_bc=.true.)
        iw_H2II_density = var_set_fluxvar('H2II', 'H2II', need_bc=.true.)
        ! + 3 species : D, D+, HD.
        iw_DI_density = var_set_fluxvar('DI', 'DI', need_bc=.true.)
        iw_DII_density = var_set_fluxvar('DII', 'DII', need_bc=.true.)
        iw_HDI_density = var_set_fluxvar('HDI', 'HDI', need_bc=.true.)
      else
        write(*,*) 'Unknown primordial chemistry network'
        call mpistop('You must set hd_gr_primordial_chemistry to  0,1,2 or 3!')
      end if prim_chem
      !metal cooling
      iw_metal_density = var_set_fluxvar('metal_density', 'metal_density', need_bc=.true.)

      !dust chemistry
      iw_dust_density = var_set_fluxvar('dust_density', 'dust_density', need_bc=.true.)

      !six species: H, H+, He, He+, He++, e- (in this order)
      hd_ind%HI_density_ = iw_HI_density
      hd_ind%HII_density_ = iw_HII_density
      hd_ind%HeI_density_ = iw_HeI_density
      hd_ind%HeII_density_ = iw_HeII_density
      hd_ind%HeIII_density_ = iw_HeIII_density
      hd_ind%e_density_ = iw_e_density
      ! + 3 species : H2, H-, H2+ (in this order)
      hd_ind%H2I_density_ = iw_H2I_density
      hd_ind%HM_density_ = iw_HM_density
      hd_ind%H2II_density_ = iw_H2II_density
      ! + 3 species : D, D+, HD.
      hd_ind%DI_density_ = iw_DI_density
      hd_ind%DII_density_ = iw_DII_density
      hd_ind%HDI_density_ = iw_HDI_density
      !metal cooling
      hd_ind%metal_density_ = iw_metal_density
      !dust chemistry
      hd_ind%dust_density_ = iw_dust_density
      rhoHI_ = iw_HI_density
      rhoHII_ = iw_HII_density
      rhoHeI_ = iw_HeI_density
      rhoHeII_ = iw_HeII_density
      rhoHeIII_ = iw_HeIII_density
      rhoe_ = iw_e_density
      rhoH2I_ = iw_H2I_density
      rhoH2II_ = iw_H2II_density
      rhoHM_ = iw_HM_density
      rhoDI_ = iw_DI_density
      rhoDII_ = iw_DII_density
      rhoHDI_ = iw_HDI_density
      rhometal_ = iw_metal_density
      rhodust_ = iw_dust_density

    end if using_grackle


    if(hd_config%n_tracer>0) then
     allocate(tracer(hd_config%n_tracer))

     ! Set starting index of tracers
     do itr = 1, hd_config%n_tracer
      tracer(itr) = var_set_fluxvar("trc", "trp", itr, need_bc=.true.)
     end do
    end if

    ! Set index(es) of gravity variable(s)
    if (hd_config%gravity) then
      grav_phi_=var_set_grav_phi('Gphi')
    else
      grav_phi_=-1
    end if



    allocate(hd_ind%mom(ndir))
    hd_ind%rho_        =rho_
    hd_ind%mom(:)      =mom(:)
    !hd_ind%mag(:)=-1
    hd_ind%e_          =e_
    hd_ind%pressure_   =p_
    hd_ind%grav_phi_   =grav_phi_
    hd_ind%lfac_       =-1
    hd_ind%xi_         =-1
    hd_ind%psi_        =-1



    if(hd_config%mean_mup_on)then
      hd_ind%mup_ = var_set_extravar('mup', 'mup')
    end if

    hd_ind%mythetafield_ = var_set_extravar('theta', 'theta')
    hd_ind%mytheta_zero_ = var_set_extravar('thetazero', 'thetazero')


    !>Grackle:



    iw_rhoX = -1
    iw_rhoY = -1

    iw_gamma = -1
    iw_temperature = -1

    if(hd_config%use_grackle)then

      iw_rhoX = var_set_extravar('rhoX', 'rhoX')
      iw_rhoY = var_set_extravar('rhoY', 'rhoY')

      iw_gamma = var_set_extravar('gamma', 'gamma')
      iw_temperature = var_set_extravar('temperature', 'temperature')

      hd_ind%rhoX_=iw_rhoX
      hd_ind%rhoY_=iw_rhoY

      hd_ind%gamma_ =iw_gamma
      gamma_ = iw_gamma
      hd_ind%temperature_ = iw_temperature
      temperature_ = iw_temperature

    end if




    !End of Grackle part


    if(hd_config%n_tracer>0)then
      allocate(hd_ind%tracer(hd_config%n_tracer))
      hd_ind%tracer(:) = tracer(:)
    end if
    hd_config%nw       = nw
    hd_config%nwflux   = nwflux
    hd_config%nwhllc   = nwflux
    hd_config%nwfluxbc = nwfluxbc

  end subroutine hd_fill_phys_indices

  subroutine hd_get_aux(clipping,w,x,ixI^L,ixO^L,subname,use_oldaux,sqrB,SdotB)
    use mod_global_parameters
    use mod_usr_methods
    integer, intent(in)                :: ixI^L, ixO^L
    real(kind=dp)      , intent(in)    :: x(ixI^S,1:ndim)
    real(kind=dp)      , intent(inout) :: w(ixI^S,nw)
    logical, intent(in)                :: clipping
    character(len=*), intent(in)       :: subname
    logical         , intent(in), target, optional :: use_oldaux(ixI^S)
    real(kind=dp)   , intent(in), target, optional :: sqrB(ixI^S),SdotB(ixI^S)

    !---local------
    double precision, allocatable :: gravity_potential(:^D&)
    integer :: iw_g_dir,idim
    !---------------------------------------------------------------------------
    if (.not. associated(usr_gravity_potential)) then
      write(*,*) "mod_usr.t: please point usr_gravity_potential to a subroutine"
      write(*,*) "like the phys_gravity_potential in mod_usr_methods.t"
      call mpistop("hd_get_aux: usr_gravity_potential not defined")
    else
    if(.not.allocated(gravity_potential))allocate(gravity_potential(ixI^S))
      gravity_potential(ixO^S)=1.0_dp
      call usr_gravity_potential(ixI^L,ixO^L,w,x,gravity_potential)
      !==========================================
      w(ixO^S,grav_phi_)=gravity_potential(ixO^S)
      !==========================================
    if(allocated(gravity_potential))deallocate(gravity_potential)
  end if



end subroutine hd_get_aux

  subroutine hd_check_params
    use mod_global_parameters
    use mod_dust, only: dust_check_params

    if (.not. hd_config%energy) then
       if (hd_config%gamma <= 0.0_dp) call mpistop ("Error: hd_gamma <= 0")
       if (hd_config%adiab < 0.0_dp) call mpistop ("Error: hd_adiab < 0")
       hd_config%small_pressure= hd_config%adiab*hd_config%small_density**hd_config%gamma
    else
       if (hd_config%gamma <= 0.0_dp .or. hd_config%gamma == 1.0_dp) &
            call mpistop ("Error: hd_gamma <= 0 or hd_gamma == 1.0")
       small_e = hd_config%small_pressure/gamma_1
    end if

    if (hd_config%dust_on) call dust_check_params()

  end subroutine hd_check_params


  subroutine hd_physical_units
    use mod_global_parameters
    use mod_dust, only : dust_physical_units
    real(dp) :: mp,kB
    ! Derive scaling units


    if(SI_unit) then
      mp=mp_SI
      kB=kB_SI
    else
      mp=mp_cgs
      kB=kB_cgs
    end if

    if(unit_velocity<smalldouble) then
      unit_density=phys_config%mean_mass*mp*unit_numberdensity
      unit_pressure= unit_density/(phys_config%mean_mup*mp)*kB*unit_temperature
      unit_velocity=dsqrt(unit_pressure/unit_density)
      unit_time=unit_length/unit_velocity
    else
      unit_density=phys_config%mean_mass*mp*unit_numberdensity
      unit_pressure=unit_density*unit_velocity**2.0_dp
      unit_temperature=phys_config%mean_mup*mp*unit_pressure/(unit_density*kB)
      unit_time=unit_length/unit_velocity




    end if

    hd_config%unit_velocity = unit_velocity
    hd_config%unit_temperature = unit_temperature




  end subroutine hd_physical_units

  subroutine hd_fill_convert_factor
    use mod_global_parameters
    use mod_dust, only : dust_physical_units
    use mod_radiative_cooling, only: rad_cooling_physical_units
    !---------------------------------------------------
    if(allocated(w_convert_factor))deallocate(w_convert_factor)

     allocate(w_convert_factor(nw))
     ! initialisation of the convert factor for the array
     w_convert_factor           = 1.0_dp
     ! set the convertion factor for the density
     w_convert_factor(rho_)     = unit_density
     ! set the convertion factor for the energy and pressure
     if(hd_config%energy)w_convert_factor(e_)       = unit_density*unit_velocity**2.0_dp

     ! set the convertion factor for the speed and the moment
     if(saveprim)then
      w_convert_factor(mom(:))  = unit_velocity
     else
      w_convert_factor(mom(:))  = unit_density*unit_velocity
     end if


     if(hd_config%use_grackle)then
      !HI
      w_convert_factor(rhoHI_)=unit_density
      !HII
      w_convert_factor(rhoHII_)=unit_density
      !HeI
      w_convert_factor(rhoHeI_)=unit_density
      !HeII
      w_convert_factor(rhoHeII_)=unit_density
      !HeIII
      w_convert_factor(rhoHeIII_)=unit_density
      !electrons
      w_convert_factor(rhoe_)=unit_density
      !HM
      w_convert_factor(rhoe_)=unit_density
      !H2I
      w_convert_factor(rhoH2I_)=unit_density
      !H2II
      w_convert_factor(rhoH2II_)=unit_density
      !DI
      w_convert_factor(rhoDI_)=unit_density
      !DII
      w_convert_factor(rhoDII_)=unit_density
      !HDI
      w_convert_factor(rhoHDI_)=unit_density
      !metal
      w_convert_factor(rhometal_)=unit_density
      !dust
      w_convert_factor(rhodust_)=unit_density

      w_convert_factor(iw_rhoX)=unit_density
      w_convert_factor(iw_rhoY)=unit_density
     end if


     if (hd_config%dust_on) call dust_physical_units
     if (hd_config%radiative_cooling) call rad_cooling_physical_units
     if(hd_config%mean_mup_on)w_convert_factor(hd_ind%mup_) = hd_config%mean_mup
     w_convert_factor(hd_ind%gamma_) = 1.0_dp
     w_convert_factor(hd_ind%temperature_) = unit_temperature
     time_convert_factor   = unit_time
     length_convert_factor = unit_length


  end   subroutine hd_fill_convert_factor
  !> Returns 0 in argument flag where values are ok
  subroutine hd_check_w(primitive, ixI^L, ixO^L, w, flag)
    use mod_global_parameters
    use mod_dust, only     : dust_check_w
    use mod_chemical, only : chemical_check_w
    implicit none


    logical, intent(in)          :: primitive
    integer, intent(in)          :: ixI^L, ixO^L
    real(dp), intent(in)         :: w(ixI^S, nw)
    integer, intent(inout)       :: flag(ixI^S)
    ! .. local ..
    real(dp)                     :: tmp(ixI^S),gamm(ixI^S)
    !-----------------------
    flag(ixO^S) = 0

    if(hd_config%dust_on)call dust_check_w(primitive, ixI^L, ixO^L, flag, w)
    if(hd_config%chemical_on)call chemical_check_w(primitive, ixI^L, ixO^L, flag, w)

    if(hd_config%use_grackle)then
      where(w(ixO^S, rho_) < hd_config%small_density) flag(ixO^S) = rho_
      where(w(ixO^S, rhoHI_) < hd_config%small_density) flag(ixO^S) = rhoHI_
      where(w(ixO^S, rhoHII_) < hd_config%small_density) flag(ixO^S) = rhoHII_
      where(w(ixO^S, rhoHM_) < hd_config%small_density) flag(ixO^S) = rhoHM_
      where(w(ixO^S, rhoHeI_) < hd_config%small_density) flag(ixO^S) = rhoHeI_
      where(w(ixO^S, rhoHeII_) < hd_config%small_density) flag(ixO^S) = rhoHeII_
      where(w(ixO^S, rhoHeIII_) < hd_config%small_density) flag(ixO^S) = rhoHeIII_
      where(w(ixO^S, rhoH2I_) < hd_config%small_density) flag(ixO^S) = rhoH2I_
      where(w(ixO^S, rhoH2II_) < hd_config%small_density) flag(ixO^S) = rhoH2II_
      where(w(ixO^S, rhoDI_) < hd_config%small_density) flag(ixO^S) = rhoDI_
      where(w(ixO^S, rhoDII_) < hd_config%small_density) flag(ixO^S) = rhoDII_
      where(w(ixO^S, rhoHDI_) < hd_config%small_density) flag(ixO^S) = rhoHDI_
      where(w(ixO^S, rhoe_) < hd_config%small_density) flag(ixO^S) = rhoe_
      where(w(ixO^S, rhometal_) < hd_config%small_density) flag(ixO^S) = rhometal_
      where(w(ixO^S, rhodust_) < hd_config%small_density) flag(ixO^S) = rhodust_
    end if

    if (hd_config%energy) then
       if (primitive) then
          where(w(ixO^S, e_) < hd_config%small_pressure) flag(ixO^S) = e_
       else
         call hd_get_gamma(w, ixI^L, ixO^L, gamm)
         where(w(ixO^S, rho_) > hd_config%small_density)
          !old : tmp(ixO^S) = gamma_1*(w(ixO^S, e_) - &
          tmp(ixO^S) = (gamm(ixO^S)-1.0_dp)*(w(ixO^S, e_) - &
               0.5_dp * sum(w(ixO^S, mom(:))**2.0_dp, dim=ndim+1) / w(ixO^S, rho_))
          elsewhere
            tmp(ixO^S) =   w(ixO^S, e_)
         end where
          where(tmp(ixO^S) < hd_config%small_pressure) flag(ixO^S) = e_
       endif
    end if




  end subroutine hd_check_w


  !> Returns 0 in argument flag where values are ok
  subroutine hd_check_pth(ixI^L, ixO^L, pth, flag)
    use mod_global_parameters
    use mod_dust, only : dust_check_w

    integer, intent(in)          :: ixI^L, ixO^L
    real(dp), intent(in)         :: pth(ixI^S)
    integer, intent(inout)       :: flag(ixI^S)
    !-----------------------
    flag(ixO^S) = 0
    if (hd_config%energy) then
      where(pth(ixO^S) < hd_config%small_pressure) flag(ixO^S) = e_
    end if
  end subroutine hd_check_pth

  !> Transform primitive variables into conservative ones
  subroutine hd_to_conserved(ixI^L, ixO^L, w, x)
    use mod_global_parameters
    use mod_dust, only: dust_to_conserved
    integer, intent(in)             :: ixI^L, ixO^L
    real(dp), intent(inout) :: w(ixI^S, nw)
    real(dp), intent(in)    :: x(ixI^S, 1:ndim)
  !  real(dp)                :: invgam
    real(dp)                :: gamm(ixI^S)
    integer                         :: idir, itr

    if (hd_config%energy) then
    !   invgam = 1.0_dp/(hd_config%gamma - 1.0_dp)
       ! Calculate total energy from pressure and kinetic energy
       ! old :
       !w(ixO^S, e_) = w(ixO^S, e_)/gamma_1  + &
       !      0.5_dp * sum(w(ixO^S, mom(:))**2.0_dp, dim=ndim+1) * w(ixO^S, rho_)
       ! new :
       call hd_get_gamma(w, ixI^L, ixO^L, gamm)
       w(ixO^S, e_) = w(ixO^S, e_)/(gamm(ixO^S)-1.0_dp)  + &
             0.5_dp * sum(w(ixO^S, mom(:))**2.0_dp, dim=ndim+1) * w(ixO^S, rho_)
    end if

    ! Convert velocity to momentum
    do idir = 1, ndir
       w(ixO^S, mom(idir)) = w(ixO^S, rho_) * w(ixO^S, mom(idir))
    end do

    if (hd_config%dust_on) then
      call dust_to_conserved(ixI^L, ixO^L, w, x)
    end if

    if (check_small_values) &
     call hd_handle_small_values(.false., w, x, ixI^L, ixO^L, 'hd_to_conserved')

  end subroutine hd_to_conserved

  !> Transform conservative variables into primitive ones
  subroutine hd_to_primitive(ixI^L, ixO^L, w, x)
    use mod_global_parameters
    use mod_dust, only: dust_to_primitive
    implicit none

    integer, intent(in)             :: ixI^L, ixO^L
    real(dp), intent(inout)         :: w(ixI^S, nw)
    real(dp), intent(in)            :: x(ixI^S, 1:ndim)
    ! .. local ..
    integer                         :: itr, idir
    real(dp)                        :: inv_rho(ixO^S)
    real(dp)                        :: gamm(ixI^S)
    !--------------------------------
    if (check_small_values) call hd_handle_small_values(.true., w, x, &
                                              ixI^L, ixO^L, 'hd_to_primitive')

    inv_rho = 1.0_dp / w(ixO^S, rho_)

    if (hd_config%energy) then
       ! Compute pressure
       ! old :
       !w(ixO^S, e_) = gamma_1 * (w(ixO^S, e_) - &
       !     hd_kin_en(w, ixI^L, ixO^L, inv_rho))
       call hd_get_gamma(w,ixI^L, ixO^L, gamm)
       w(ixO^S, e_) = (gamm(ixO^S)-1.0_dp)* (w(ixO^S, e_) - &
            hd_kin_en(w, ixI^L, ixO^L, inv_rho))
    end if

    ! Convert momentum to velocity
    do idir = 1, ndir
       w(ixO^S, mom(idir)) = w(ixO^S, mom(idir)) * inv_rho
    end do

    ! Convert dust momentum to dust velocity
    if (hd_config%dust_on) then
      call dust_to_primitive(ixI^L, ixO^L, w, x)
    end if

    if (check_small_values) call hd_handle_small_values(.true., w, x, ixI^L, ixO^L, 'hd_to_primitive')

  end subroutine hd_to_primitive

  subroutine e_to_rhos(ixI^L, ixO^L, w, x)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    real(dp)             :: w(ixI^S, nw)
    real(dp), intent(in) :: x(ixI^S, 1:ndim)
    real(dp), dimension(ixI^S)   :: gamm

    if (hd_config%energy) then
       call hd_get_gamma(w, ixI^L, ixO^L, gamm)
       !w(ixO^S, e_) = gamma_1 * &
       w(ixO^S, e_) = (gamm(ixO^S)-1.0_dp) * &
       w(ixO^S, rho_)**(-(gamm(ixO^S)-1.0_dp)) * &
            (w(ixO^S, e_) - hd_kin_en(w, ixI^L, ixO^L))
    else
       call mpistop("energy from entropy can not be used with -eos = iso !")
    end if
  end subroutine e_to_rhos

  subroutine rhos_to_e(ixI^L, ixO^L, w, x)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    real(dp)             :: w(ixI^S, nw)
    real(dp), intent(in) :: x(ixI^S, 1:ndim)
    real(dp), dimension(ixI^S)   :: gamm

    if (hd_config%energy) then
       call hd_get_gamma(w, ixI^L, ixO^L, gamm)
       !w(ixO^S, e_) = w(ixO^S, rho_)**gamma_1 * w(ixO^S, e_) / &
       !    gamma_1 + hd_kin_en(w, ixI^L, ixO^L)
       w(ixO^S, e_) = w(ixO^S, rho_)**(gamm(ixO^S)-1.0_dp) * w(ixO^S, e_) / &
            (gamm(ixO^S)-1.0_dp) + hd_kin_en(w, ixI^L, ixO^L)
    else
       call mpistop("entropy from energy can not be used with -eos = iso !")
    end if
  end subroutine rhos_to_e

  !> Calculate v_i = m_i / rho within ixO^L
  subroutine hd_get_v(w, x, ixI^L, ixO^L, idim, v)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L, idim
    real(dp), intent(in)  :: w(ixI^S, nw), x(ixI^S, 1:ndim)
    real(dp), intent(out) :: v(ixI^S)

    v(ixO^S) = w(ixO^S, mom(idim)) / w(ixO^S, rho_)
  end subroutine hd_get_v

  !> Calculate cmax_idim = csound + abs(v_idim) within ixO^L
  subroutine hd_get_cmax(w, x, ixI^L, ixO^L, idim, cmax)
    use mod_global_parameters
    use mod_dust, only: dust_get_cmax

    integer, intent(in)                       :: ixI^L, ixO^L, idim
    real(dp), intent(in)              :: w(ixI^S, nw), x(ixI^S, 1:ndim)
    real(dp), intent(inout)           :: cmax(ixI^S)
    real(dp)                          :: csound(ixI^S)
    real(dp)                          :: v(ixI^S)

    call hd_get_v(w, x, ixI^L, ixO^L, idim, v)
    call hd_get_csound2(w,x,ixI^L,ixO^L,csound)
    csound(ixO^S) = sqrt(csound(ixO^S))

    !write(*,*) idim,') ','hd_get_cmax csound of iteration it = ', it,&
    !' = ', csound(ixO^S)


    !write(*,*) idim,') ','hd_get_cmax v(ixO^S) of iteration it = ', it,&
    !' = ', v(ixO^S)
    cmax(ixO^S) = abs(v(ixO^S))+csound(ixO^S)


    !write(*,*) idim,') ','hd_get_cmax cmax(ixO^S) of iteration it = ', it,&
    !' = ', cmax(ixO^S)

    if (hd_config%dust_on) then
      call dust_get_cmax(w, x, ixI^L, ixO^L, idim, cmax)
    end if
  end subroutine hd_get_cmax

  !> Calculate cmax_idim = csound + abs(v_idim) within ixO^L
  subroutine hd_get_cbounds(wLC, wRC, wLp, wRp, x, ixI^L, ixO^L, idim, cmax, cmin)
    use mod_global_parameters
    use mod_dust, only: dust_get_cmax

    integer, intent(in)               :: ixI^L, ixO^L, idim
    ! conservative left and right status
    real(dp), intent(in)              :: wLC(ixI^S, nw), wRC(ixI^S, nw)
    ! primitive left and right status
    real(dp), intent(in)              :: wLp(ixI^S, nw), wRp(ixI^S, nw)
    real(dp), intent(in)              :: x(ixI^S, 1:ndim)
    real(dp), intent(inout)           :: cmax(ixI^S)
    real(dp), intent(inout), optional :: cmin(ixI^S)
    ! .. local ..
    real(dp)                          :: wmean(ixI^S,nw)
    real(dp), dimension(ixI^S)        :: umean, dmean, csoundL, csoundR,&
                                         tmp1,tmp2,tmp3
    real(dp), dimension(ixI^S)   :: gammLp,gammRp
    !---------------------------------------------------
    if (typeboundspeed/='cmaxmean') then
      ! This implements formula (10.52) from "Riemann Solvers and Numerical
      ! Methods for Fluid Dynamics" by Toro.

      tmp1(ixO^S)=sqrt(wLp(ixO^S,rho_))
      tmp2(ixO^S)=sqrt(wRp(ixO^S,rho_))
      tmp3(ixO^S)=1.0_dp/(sqrt(wLp(ixO^S,rho_))+sqrt(wRp(ixO^S,rho_)))
      umean(ixO^S)=(wLp(ixO^S,mom(idim))*tmp1(ixO^S)+wRp(ixO^S,mom(idim))*tmp2(ixO^S))*tmp3(ixO^S)

      if(hd_config%energy) then
        !csoundL(ixO^S)=hd_config%gamma*wLp(ixO^S,p_)/wLp(ixO^S,rho_)
        !csoundR(ixO^S)=hd_config%gamma*wRp(ixO^S,p_)/wRp(ixO^S,rho_)
        call hd_get_gamma(wLp, ixI^L, ixO^L, gammLp)
        csoundL(ixO^S)=gammLp(ixO^S)*wLp(ixO^S,p_)/wLp(ixO^S,rho_)
        call hd_get_gamma(wRp, ixI^L, ixO^L, gammRp)
        csoundR(ixO^S)=gammRp(ixO^S)*wRp(ixO^S,p_)/wRp(ixO^S,rho_)
      else
        csoundL(ixO^S)=hd_config%gamma*hd_config%adiab*wLp(ixO^S,rho_)**gamma_1
        csoundR(ixO^S)=hd_config%gamma*hd_config%adiab*wRp(ixO^S,rho_)**gamma_1
      end if

      dmean(ixO^S) = (tmp1(ixO^S)*csoundL(ixO^S)+tmp2(ixO^S)*csoundR(ixO^S)) * &
           tmp3(ixO^S) + 0.5_dp*tmp1(ixO^S)*tmp2(ixO^S)*tmp3(ixO^S)**2 * &
           (wRp(ixO^S,mom(idim))-wLp(ixO^S,mom(idim)))**2

      dmean(ixO^S)=sqrt(dmean(ixO^S))
      if(present(cmin)) then
        cmin(ixO^S)=umean(ixO^S)-dmean(ixO^S)
        cmax(ixO^S)=umean(ixO^S)+dmean(ixO^S)
      else
        cmax(ixO^S)=abs(umean(ixO^S))+dmean(ixO^S)
      end if

      if (hd_config%dust_on) then
        wmean(ixO^S,1:nwflux)=0.5d0*(wLC(ixO^S,1:nwflux)+wRC(ixO^S,1:nwflux))
        call dust_get_cmax(wmean, x, ixI^L, ixO^L, idim, cmax, cmin)
      end if

    else

      wmean(ixO^S,1:nwflux)=0.5d0*(wLC(ixO^S,1:nwflux)+wRC(ixO^S,1:nwflux))
      tmp1(ixO^S)=wmean(ixO^S,mom(idim))/wmean(ixO^S,rho_)
      call hd_get_csound2(wmean,x,ixI^L,ixO^L,csoundR)
      csoundR(ixO^S) = sqrt(csoundR(ixO^S))

      if(present(cmin)) then
        cmax(ixO^S)=max(tmp1(ixO^S)+csoundR(ixO^S),zero)
        cmin(ixO^S)=min(tmp1(ixO^S)-csoundR(ixO^S),zero)
      else
        cmax(ixO^S)=abs(tmp1(ixO^S))+csoundR(ixO^S)
      end if

      if (hd_config%dust_on) then
        call dust_get_cmax(wmean, x, ixI^L, ixO^L, idim, cmax, cmin)
      end if
    end if

  end subroutine hd_get_cbounds

  !> Calculate the square of the thermal sound speed csound2 within ixO^L.
  !> csound2=gamma*p/rho
  subroutine hd_get_csound2(w,x,ixI^L,ixO^L,csound2)
    use mod_global_parameters
    integer, intent(in)     :: ixI^L, ixO^L
    real(dp), intent(in)    :: w(ixI^S,nw)
    real(dp), intent(in)    :: x(ixI^S,1:ndim)
    real(dp), intent(out)   :: csound2(ixI^S)
    real(dp), dimension(ixI^S)   :: gamm
    !------------------------------------------
    if(hd_config%energy) then

    !write(*,*) 'hd_get_csound2 w(ixO^S,rho_) of iteration it = ', it,&
    !' = ', w(ixO^S,rho_)
      call hd_get_pthermal(w,x,ixI^L,ixO^L,csound2, 'hd_get_csound2')
    !write(*,*) 'hd_get_csound2 p(csound2) of iteration it = ', it,&
    !' = ', csound2(ixO^S)
      !csound2(ixO^S)=hd_config%gamma*csound2(ixO^S)/w(ixO^S,rho_)
      call hd_get_gamma(w, ixI^L, ixO^L, gamm)
      !write(*,*) 'hd_get_csound2 gamm(csound2) of iteration it = ', it,&
      !' = ', gamm(ixO^S)
      csound2(ixO^S)=gamm(ixO^S)*csound2(ixO^S)/w(ixO^S,rho_)
    else
      if(dabs(gamma_1)<smalldouble) then
        csound2(ixO^S)=hd_config%adiab
      else
        csound2(ixO^S)=hd_config%gamma*hd_config%adiab*w(ixO^S,rho_)**gamma_1
      end if
    end if
  end subroutine hd_get_csound2

  !> Calculate thermal pressure=(gamma-1)*(e-0.5*m**2/rho) within ixO^L
  subroutine hd_get_pthermal(w, x, ixI^L, ixO^L, pth,subname)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    real(dp), intent(in)         :: w(ixI^S, nw)
    real(dp), intent(in)         :: x(ixI^S, 1:ndim)
    real(dp), intent(out)        :: pth(ixI^S)
    character(len=*), optional     :: subname
    real(dp)                     :: gamm(ixI^S),sgmerr(1:1)
    !----------------------------------------------------
    if(present(subname)) then
      !write(*,*) 'hd_get_pthermal called from ', trim(subname)
    else
     !write(*,*) 'hd_get_pthermal called from outside'
    end if

    if (hd_config%energy) then
      call hd_get_gamma(w,ixI^L, ixO^L, gamm)
      !write(*,*) 'hd_get_pthermal gamma(csound2) of iteration it = ', it,&
      !' = ', gamm(ixO^S)
      !write(*,*) 'hd_get_pthermal w(ixO^S, e_) of iteration it = ', it,&
      !' = ', w(ixO^S, e_)*w_convert_factor(e_)
      pth(ixO^S) = (gamm(ixO^S)-1.0_dp) * (w(ixO^S, e_) - &
           hd_kin_en(w, ixI^L, ixO^L))
      !pth(ixO^S) = gamma_1 * (w(ixO^S, e_) - &
      !    hd_kin_en(w, ixI^L, ixO^L))
    else
       pth(ixO^S) = hd_config%adiab * w(ixO^S, rho_)**hd_config%gamma
    end if

    !if (check_small_values)call hd_handle_small_values_pressure(ixI^L,ixO^L,x,&
    !                                                            'ptherm',pth,w)
  end subroutine hd_get_pthermal

  !> Calculate adiabatic index
  subroutine hd_get_gamma(w,ixI^L, ixO^L, gammaeff)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    real(dp), intent(in)         :: w(ixI^S, nw)
    real(dp), intent(out)        :: gammaeff(ixI^S)
    real(kind=dp), dimension(ixI^S ):: nH2, nother
    real(kind=dp), dimension(ixI^S ):: gammaoth
    real(kind=dp), dimension(ixI^S ):: inv_gammaH2_1, inv_gammaoth_1
    real(kind=dp)                :: tvar
    !real(kind=dp), dimension(ixI^S ):: tvar
    integer :: imesh^D
    !----------------------------------------------------
    !write(*,*) 'temperature field = ', w(ixO^S,temperature_)*w_convert_factor(temperature_)
    if (hd_config%use_grackle) then
      !write(*,*) 'hd_get_gamma H2I',w(ixO^S,rhoH2I_)
      !write(*,*) 'hd_get_gamma H2II',w(ixO^S,rhoH2II_)
      !write(*,*) 'hd_get_gamma HeI',w(ixO^S,rhoHeI_)
      !write(*,*) 'hd_get_gamma HeII',w(ixO^S,rhoHeII_)
      !write(*,*) 'hd_get_gamma HeIII',w(ixO^S,rhoHeIII_)
      !write(*,*) 'hd_get_gamma HI',w(ixO^S,rhoHI_)
      !write(*,*) 'hd_get_gamma HII',w(ixO^S,rhoHII_)
      !write(*,*) 'hd_get_gamma e-',w(ixO^S,rhoe_)
      nH2(ixO^S) = 0.5d0*(w(ixO^S,rhoH2I_) + &
      w(ixO^S,rhoH2II_))
      nother(ixO^S) = (w(ixO^S,rhoHeI_) + w(ixO^S,rhoHeII_) +&
                   w(ixO^S,rhoHeIII_))/4.0d0 +&
                   w(ixO^S,rhoHI_) + w(ixO^S,rhoHII_) +&
                   w(ixO^S,rhoe_)
      {do imesh^D=ixOmin^D,ixOmax^D|\}
      if(nH2(imesh^D)/nother(imesh^D) > 1.0d-3)then
        tvar = 6100.0d0/(w(imesh^D, temperature_)*&
        w_convert_factor(temperature_))
        if (tvar > 10.0)then
          inv_gammaH2_1(imesh^D) = 0.5d0*5.0d0
        else
          inv_gammaH2_1(imesh^D) = 0.5d0*(5.0d0 +&
                             2.0d0 * tvar**2.0d0 *&
                          DEXP(tvar)/(DEXP(tvar)-1.0d0)**2.0d0)
        end if
      else
        inv_gammaH2_1(imesh^D) = 2.5d0
      end if
      {end do^D&|\}
      !tvar(ixO^S)=6100.0d0/(w(ixO^S, temperature_)*&
      !w_convert_factor(temperature_))
      !inv_gammaH2_1(ixO^S)=merge(merge(0.5d0*5.0d0,&
      !                0.5d0*(5.0d0 +&
      !                2.0d0*tvar(ixO^S)**2.0d0*&
      !                DEXP(tvar(ixO^S))/&
      !                (DEXP(tvar(ixO^S))-1.0d0)**2.0d0),&
      !tvar(ixO^S)>10.0d0),&
      !2.5d0,nH2(ixO^S)/nother(ixO^S)>1.0d-3)

      gammaoth(ixO^S) = hd_config%gamma
      inv_gammaoth_1(ixO^S) = 1.0d0/(gammaoth(ixO^S)-1.0d0)

      gammaeff(ixO^S) = 1.0d0 + (nH2(ixO^S) + nother(ixO^S))/&
                   (nH2(ixO^S)*inv_gammaH2_1(ixO^S) +&
                nother(ixO^S)*inv_gammaoth_1(ixO^S))
    else
       gammaeff(ixO^S) = hd_config%gamma
    end if

    !if (check_small_values)call hd_handle_small_values_pressure(ixI^L,ixO^L,x,&
    !                                                            'gammaeff',gammaeff,w)
  end subroutine hd_get_gamma

  !> Calculate mean molecular weight
  subroutine hd_get_mmw(w, x, ixI^L, ixO^L, mmw)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    real(dp), intent(in)         :: w(ixI^S, nw)
    real(dp), intent(in)         :: x(ixI^S, 1:ndim)
    real(dp), intent(out)        :: mmw(ixI^S)
    !----------------------------------------------------
    if (hd_config%use_grackle) then
      mmw(ixO^S) = w(ixO^S,rho_)/&
        (w(ixO^S,rhoHI_)+w(ixO^S,rhoHII_)+&
        w(ixO^S,rhoHM_)+&
        0.25d0*(w(ixO^S,rhoHeI_)+w(ixO^S,rhoHeII_)+&
        w(ixO^S,rhoHeIII_))+&
        0.5d0*(w(ixO^S,rhoH2I_)+w(ixO^S,rhoH2II_))+&
        w(ixO^S,rhometal_)/16.0)

       !mmw(ixO^S) = w(ixO^S,hd_ind%H2I_density_)
    else
       mmw(ixO^S) = hd_config%mean_mup
    end if

    !if (check_small_values)call hd_handle_small_values_pressure(ixI^L,ixO^L,x,&
    !                                                            'mmw',mmw,w)
  end subroutine hd_get_mmw

  !> Calculate temperature within ixO^L
  subroutine hd_get_temperature( ixI^L, ixO^L,w, x, temperature)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    real(dp), intent(in)         :: w(ixI^S, nw)
    real(dp), intent(in)         :: x(ixI^S, 1:ndim)
    real(dp), intent(out)        :: temperature(ixI^S)
    real(dp), dimension(ixI^S)   :: pth
    real(dp), dimension(ixI^S)   :: meanmw
    logical , dimension(ixI^S)   :: patch_mult_mup
    real(kind=dp)  :: mp,kB
    !----------------------------------------------------

    if(SI_unit) then
        mp=mp_SI
        kB=kB_SI
    else
        mp=mp_cgs
        kB=kB_cgs
    end if


    if (hd_config%energy) then
      call hd_get_pthermal(w, x, ixI^L, ixO^L, pth,'hd_get_temperature')
      temperature(ixO^S) =&
      mp*w_convert_factor(p_)*pth(ixO^S)/&
      (kB*w_convert_factor(rho_)*w(ixO^S,rho_))
      if(hd_config%mean_mup_on)then
        call hd_get_mmw(w, x, ixI^L, ixO^L, meanmw)
        temperature(ixO^S) =temperature(ixO^S) *&
        meanmw(ixO^S)
      end if
      temperature(ixO^S) =temperature(ixO^S)/unit_temperature
    else
      Temperature(ixO^S) = hd_config%temperature_isotherm !necessarilly isotherm
      if(hd_config%mean_mup_on)then
        patch_mult_mup(ixO^S) = dabs(w(ixO^S, phys_ind%mup_)-1.0_dp)>smalldouble &
             .or. w(ixO^S, phys_ind%mup_)>smalldouble
        if(any(patch_mult_mup(ixO^S)))then
          where(patch_mult_mup(ixO^S))
            temperature(ixO^S) =temperature(ixO^S)*w(ixO^S,phys_ind%mup_)
          end where
        end if
      end if
    end if

    return

    if (hd_config%energy) then
      call hd_get_pthermal(w, x, ixI^L, ixO^L, pth)
      temperature(ixO^S) = pth(ixO^S)/w(ixO^S, rho_)
      if(hd_config%mean_mup_on)then
        patch_mult_mup(ixO^S) = dabs(w(ixO^S, phys_ind%mup_)-1.0_dp)>smalldouble &
             .or. w(ixO^S, phys_ind%mup_)>smalldouble
        if(any(patch_mult_mup(ixO^S)))then
          where(patch_mult_mup(ixO^S))
            temperature(ixO^S) =temperature(ixO^S)*w(ixO^S,phys_ind%mup_)
          end where
        end if
      end if
    else
      Temperature(ixO^S) = hd_config%temperature_isotherm !necessarilly isotherm
      if(hd_config%mean_mup_on)then
        patch_mult_mup(ixO^S) = dabs(w(ixO^S, phys_ind%mup_)-1.0_dp)>smalldouble &
             .or. w(ixO^S, phys_ind%mup_)>smalldouble
        if(any(patch_mult_mup(ixO^S)))then
          where(patch_mult_mup(ixO^S))
            temperature(ixO^S) =temperature(ixO^S)*w(ixO^S,phys_ind%mup_)
          end where
        end if
      end if
    end if
  end subroutine hd_get_temperature


  ! Calculate flux f_idim[iw]
  subroutine hd_get_flux_cons(w, x, ixI^L, ixO^L, idim, f)
    use mod_global_parameters
    use mod_dust, only: dust_get_flux
    use mod_chemical, only: chemical_get_flux
    integer, intent(in)             :: ixI^L, ixO^L, idim
    real(dp), intent(in)            :: w(ixI^S, 1:nw), x(ixI^S, 1:ndim)
    real(dp), intent(out)           :: f(ixI^S, nwflux)
    real(dp)                        :: pth(ixI^S), v(ixI^S)
    integer                         :: idir, itr

    call hd_get_pthermal(w, x, ixI^L, ixO^L, pth,'hd_get_flux_cons')
    call hd_get_v(w, x, ixI^L, ixO^L, idim, v)

    f(ixO^S, rho_) = v(ixO^S) * w(ixO^S, rho_)

    !Mialy s addings:

    if(hd_config%use_grackle)then
      f(ixO^S, rhoHI_) = v(ixO^S) * w(ixO^S, rhoHI_)
      f(ixO^S, rhoHII_) = v(ixO^S) * w(ixO^S, rhoHII_)
      f(ixO^S, rhoHM_) = v(ixO^S) * w(ixO^S, rhoHM_)
      f(ixO^S, rhoHeI_) = v(ixO^S) * w(ixO^S, rhoHeI_)
      f(ixO^S, rhoHeII_) = v(ixO^S) * w(ixO^S, rhoHeII_)
      f(ixO^S, rhoHeIII_) = v(ixO^S) * w(ixO^S, rhoHeIII_)
      f(ixO^S, rhoH2I_) = v(ixO^S) * w(ixO^S, rhoH2I_)
      f(ixO^S, rhoH2II_) = v(ixO^S) * w(ixO^S, rhoH2II_)
      f(ixO^S, rhoe_) = v(ixO^S) * w(ixO^S, rhoe_)
      f(ixO^S, rhoDI_) = v(ixO^S) * w(ixO^S, rhoDI_)
      f(ixO^S, rhoDII_) = v(ixO^S) * w(ixO^S, rhoDII_)
      f(ixO^S, rhoHDI_) = v(ixO^S) * w(ixO^S, rhoHDI_)
      f(ixO^S, rhometal_) = v(ixO^S) * w(ixO^S, rhometal_)
      f(ixO^S, rhodust_) = v(ixO^S) * w(ixO^S, rhodust_)
    end if

    ! Momentum flux is v_i*m_i, +p in direction idim
    do idir = 1, ndir
      f(ixO^S, mom(idir)) = v(ixO^S) * w(ixO^S, mom(idir))
    end do

    f(ixO^S, mom(idim)) = f(ixO^S, mom(idim)) + pth(ixO^S)

    if(hd_config%energy) then
      ! Energy flux is v_i*e + v*p ! Check? m_i/rho*p
      f(ixO^S, e_) = v(ixO^S) * (w(ixO^S, e_) + pth(ixO^S))
    end if

    do itr = 1, hd_config%n_tracer
       f(ixO^S, tracer(itr)) = v(ixO^S) * w(ixO^S, tracer(itr))
    end do

    ! Dust fluxes
    if (hd_config%dust_on) then
      call dust_get_flux(w, x, ixI^L, ixO^L, idim, f)
    end if

    ! chemical fluxes

    if (hd_config%chemical_on) then
      call chemical_get_flux(w, x, ixI^L, ixO^L, idim, f)
    end if
  end subroutine hd_get_flux_cons

  ! Calculate flux f_idim[iw]
  subroutine hd_get_flux(wC, w, x, ixI^L, ixO^L, idim, f)
    use mod_global_parameters
    use mod_dust, only: dust_get_flux_prim
    use mod_chemical, only: chemical_get_flux_prim
    use mod_viscosity, only: visc_get_flux_prim ! viscInDiv

    integer, intent(in)             :: ixI^L, ixO^L, idim
    ! conservative w
    real(dp), intent(in)            :: wC(ixI^S, 1:nw)
    ! primitive w
    real(dp), intent(in)            :: w(ixI^S, 1:nw)
    real(dp), intent(in)            :: x(ixI^S, 1:ndim)
    real(dp), intent(out)           :: f(ixI^S, nwflux)
    real(dp)                        :: pth(ixO^S)
    integer                         :: idir, itr

    if (hd_config%energy) then
       pth(ixO^S) = w(ixO^S,p_)
    else
       pth(ixO^S) = hd_config%adiab * w(ixO^S, rho_)**hd_config%gamma
    end if

    f(ixO^S, rho_) = w(ixO^S,mom(idim)) * w(ixO^S, rho_)

    if(hd_config%use_grackle)then
      f(ixO^S, rhoHI_) = w(ixO^S,mom(idim)) * w(ixO^S, rhoHI_)
      f(ixO^S, rhoHII_) = w(ixO^S,mom(idim)) * w(ixO^S, rhoHII_)
      f(ixO^S, rhoHM_) = w(ixO^S,mom(idim)) * w(ixO^S, rhoHM_)
      f(ixO^S, rhoHeI_) = w(ixO^S,mom(idim)) * w(ixO^S, rhoHeI_)
      f(ixO^S, rhoHeII_) = w(ixO^S,mom(idim)) * w(ixO^S, rhoHeII_)
      f(ixO^S, rhoHeIII_) = w(ixO^S,mom(idim)) * w(ixO^S, rhoHeIII_)
      f(ixO^S, rhoH2I_) = w(ixO^S,mom(idim)) * w(ixO^S, rhoH2I_)
      f(ixO^S, rhoH2II_) = w(ixO^S,mom(idim)) * w(ixO^S, rhoH2II_)
      f(ixO^S, rhoe_) = w(ixO^S,mom(idim)) * w(ixO^S, rhoe_)
      f(ixO^S, rhoDI_) = w(ixO^S,mom(idim)) * w(ixO^S, rhoDI_)
      f(ixO^S, rhoDII_) = w(ixO^S,mom(idim)) * w(ixO^S, rhoDII_)
      f(ixO^S, rhoHDI_) = w(ixO^S,mom(idim)) * w(ixO^S, rhoHDI_)
      f(ixO^S, rhometal_) = w(ixO^S,mom(idim)) * w(ixO^S, rhometal_)
      f(ixO^S, rhodust_) = w(ixO^S,mom(idim)) * w(ixO^S, rhodust_)
    end if

    ! Momentum flux is v_i*m_i, +p in direction idim
    do idir = 1, ndir
      f(ixO^S, mom(idir)) = w(ixO^S,mom(idim)) * wC(ixO^S, mom(idir))
    end do

    f(ixO^S, mom(idim)) = f(ixO^S, mom(idim)) + pth(ixO^S)

    if(hd_config%energy) then
      ! Energy flux is v_i*e + v*p ! Check? m_i/rho*p
      f(ixO^S, e_) = w(ixO^S,mom(idim)) * (wC(ixO^S, e_) + w(ixO^S,p_))
    end if

    do itr = 1, hd_config%n_tracer
       f(ixO^S, tracer(itr)) = w(ixO^S,mom(idim)) * w(ixO^S, tracer(itr))
    end do

    ! Dust fluxes
    if (hd_config%dust_on) then
      call dust_get_flux_prim(w, x, ixI^L, ixO^L, idim, f)
    end if

    ! Viscosity fluxes - viscInDiv
    if (hd_config%viscosity) then
      call visc_get_flux_prim(w, x, ixI^L, ixO^L, idim, f, hd_config%energy)
    endif

    ! chemical fluxes
    if (hd_config%chemical_on) then
      call chemical_get_flux_prim(w, x, ixI^L, ixO^L, idim, f)
    end if
  end subroutine hd_get_flux

  !> Add geometrical source terms to w
  !>
  !> Notice that the expressions of the geometrical terms depend only on ndir,
  !> not ndim. Eg, they are the same in 2.5D and in 3D, for any geometry.
  !>
  !> Ileyk : to do :
  !>     - address the source term for the dust
  subroutine hd_add_source_geom(qdt, ixI^L, ixO^L, wCT, w, x)
    use mod_global_parameters
    use mod_viscosity, only: visc_add_source_geom ! viscInDiv

    integer, intent(in)             :: ixI^L, ixO^L
    real(dp), intent(in)    :: qdt, x(ixI^S, 1:ndim)
    real(dp), intent(inout) :: wCT(ixI^S, 1:nw), w(ixI^S, 1:nw)
    ! to change and to set as a parameter in the parfile once the possibility to
    ! solve the equations in an angular momentum conserving form has been
    ! implemented (change tvdlf.t eg)
    real(dp) :: tmp(ixI^S),tmp1(ixI^S)
    integer                         :: iw,idir, h1x^L{^NOONED, h2x^L}
    integer :: mr_,mphi_ ! Polar var. names

    mr_=mom(1); mphi_=mom(1)-1+phi_ ! Polar var. names

    select case (typeaxial)
    case ("cylindrical")
       ! s[mr]=(pthermal+mphi**2/rho)/radius
       call hd_get_pthermal(wCT,x,ixI^L,ixO^L,tmp,'hd_add_source_geom 1')
       if(phi_>0) then
         tmp(ixO^S)=tmp(ixO^S)+wCT(ixO^S,mphi_)**2/wCT(ixO^S,rho_)
         w(ixO^S,mr_)=w(ixO^S,mr_)+qdt*tmp(ixO^S)/x(ixO^S,1)
         ! s[mphi]=(-mphi*mr/rho)/radius
         ! Ileyk : beware the index permutation : mphi=2 if -phi=2 (2.5D
         ! (r,theta) grids) BUT mphi=3 if -phi=3 (for 2.5D (r,z) grids)
         if(.not. angmomfix) then
           tmp(ixO^S)=-wCT(ixO^S,mphi_)*wCT(ixO^S,mr_)/wCT(ixO^S,rho_)
           w(ixO^S,mphi_)=w(ixO^S,mphi_)+qdt*tmp(ixO^S)/x(ixO^S,1)
         end if
       else
         ! s[mr]=2pthermal/radius
         w(ixO^S,mr_)=w(ixO^S,mr_)+qdt*tmp(ixO^S)/x(ixO^S,1)
       end if
    case ("spherical")
       h1x^L=ixO^L-kr(1,^D); {^NOONED h2x^L=ixO^L-kr(2,^D);}
       ! s[mr]=((mtheta**2+mphi**2)/rho+2*p)/r
       call hd_get_pthermal(wCT,x,ixI^L,ixO^L,tmp1,'hd_add_source_geom 2')
       tmp(ixO^S)=tmp1(ixO^S)*x(ixO^S,1) &
            *(block%surfaceC(ixO^S,1)-block%surfaceC(h1x^S,1)) &
            /block%dvolume(ixO^S)
       if(ndir>1) then
         do idir=2,ndir
           tmp(ixO^S)=tmp(ixO^S)+wCT(ixO^S,mom(idir))**2/wCT(ixO^S,rho_)
         end do
       end if
       w(ixO^S,mr_)=w(ixO^S,mr_)+qdt*tmp(ixO^S)/x(ixO^S,1)

       {^NOONED
       ! s[mtheta]=-(mr*mtheta/rho)/r+cot(theta)*(mphi**2/rho+p)/r
       tmp(ixO^S)=tmp1(ixO^S)*x(ixO^S,1) &
            *(block%surfaceC(ixO^S,2)-block%surfaceC(h2x^S,2)) &
            /block%dvolume(ixO^S)
       if(ndir==3) tmp(ixO^S)=tmp(ixO^S)+(wCT(ixO^S,mom(3))**2/wCT(ixO^S,rho_))/tan(x(ixO^S,2))
       if (.not. angmomfix) tmp(ixO^S)=tmp(ixO^S)-(wCT(ixO^S,mom(2))*wCT(ixO^S,mr_))/wCT(ixO^S,rho_)
       w(ixO^S,mom(2))=w(ixO^S,mom(2))+qdt*tmp(ixO^S)/x(ixO^S,1)

       if(ndir==3) then
         ! s[mphi]=-(mphi*mr/rho)/r-cot(theta)*(mtheta*mphi/rho)/r
         if(.not. angmomfix) then
           tmp(ixO^S)=-(wCT(ixO^S,mom(3))*wCT(ixO^S,mr_))/wCT(ixO^S,rho_)&
                      -(wCT(ixO^S,mom(2))*wCT(ixO^S,mom(3)))/wCT(ixO^S,rho_)/tan(x(ixO^S,2))
           w(ixO^S,mom(3))=w(ixO^S,mom(3))+qdt*tmp(ixO^S)/x(ixO^S,1)
         end if
       end if
       }
    end select

    if (hd_viscosity) call visc_add_source_geom(qdt,ixI^L,ixO^L,wCT,w,x)

  end subroutine hd_add_source_geom

  ! w[iw]= w[iw]+qdt*S[wCT, qtC, x] where S is the source based on wCT within ixO
  subroutine hd_add_source(qdt,ixI^L,ixO^L,wCT,w,x,qsourcesplit,active)
    use mod_global_parameters
    use mod_radiative_cooling, only: radiative_cooling_add_source
    use mod_chemical, only: chemical_add_source
    use mod_dust, only: dust_add_source
    use mod_viscosity, only: viscosity_add_source
    use mod_gravity, only: gravity_add_source

    integer, intent(in)     :: ixI^L, ixO^L
    real(dp), intent(in)    :: qdt
    real(dp), intent(in)    :: wCT(ixI^S, 1:nw), x(ixI^S, 1:ndim)
    real(dp), intent(inout) :: w(ixI^S, 1:nw)
    logical, intent(in)     :: qsourcesplit
    logical, intent(inout)  :: active

    if(hd_config%dust_on) then
      call dust_add_source(qdt,ixI^L,ixO^L,wCT,w,x,qsourcesplit,active)
    end if

    if(hd_config%radiative_cooling) then
      call radiative_cooling_add_source(qdt,ixI^L,ixO^L,wCT,w,x,&
           qsourcesplit,active)
    end if

    if(hd_config%viscosity) then
      call viscosity_add_source(qdt,ixI^L,ixO^L,wCT,w,x,&
           hd_config%energy,qsourcesplit,active)
    end if

    if(hd_config%gravity) then
      call gravity_add_source(qdt,ixI^L,ixO^L,wCT,w,x,&
           hd_config%energy,qsourcesplit,active)
    end if

    if(hd_config%chemical_on) then
      call chemical_add_source(qdt,ixI^L,ixO^L,wCT,w,x,&
               qsourcesplit,active)
    end if
    if (check_small_values) call hd_handle_small_values(.false., w, x, &
                                 ixI^L, ixO^L, 'hd_add_source')

  end subroutine hd_add_source

  subroutine hd_get_dt(w, ixI^L, ixO^L, dtnew, dx^D, x)
    use mod_global_parameters
    use mod_dust, only: dust_get_dt
    use mod_radiative_cooling, only: cooling_get_dt
    use mod_chemical, only: chemical_get_dt
    use mod_viscosity, only: viscosity_get_dt
    use mod_gravity, only: gravity_get_dt

    integer, intent(in)             :: ixI^L, ixO^L
    real(dp), intent(in)    :: dx^D, x(ixI^S, 1:^ND)
    real(dp), intent(in)    :: w(ixI^S, 1:nw)
    real(dp), intent(inout) :: dtnew

    dtnew = bigdouble

    if(hd_config%dust_on) then
      call dust_get_dt(w, ixI^L, ixO^L, dtnew, dx^D, x)
    end if

    if(hd_radiative_cooling) then
      call cooling_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
    end if

    if(hd_viscosity) then
      call viscosity_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
    end if

    if(hd_config%gravity) then
      call gravity_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
    end if
    if(hd_config%chemical_on) then
      call chemical_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
    end if
  end subroutine hd_get_dt


  !beware : here, input wCT called only from mod_finite_volume.t contains
  ! conservative variables. It must be first set to primitives at the beginning
  ! of this subroutine, and wCT must be conservative again at the end of this
  !subroutine. Or create a local variable to set to primitives, so that wCT
  !is not affected.
  ! from mod_finite_volume.t:
  ! * here idir = idim for idim=idim^LIM=idimmin^D:idimmax^D
  ! which is input in mod_finite_volume.t:finite_volume(...) itself
  ! * ixI^L= finite_volume:ixI^L,
  ! and ixO^L = finite_volume:ixC^L from mod_finite_volume.t:finite_volume(...) itself,
  ! where
  ! finite_volume:ixC is centered index in the finite_volume:idim direction from ixOmin-1/2 to ixOmax+1/2
  ! ==> finite_volume:ixCmax^D=finite_volume:ixOmax^D; finite_volume:ixCmin^D=finite_volume:hxOmin^D;
  ! where finite_volume:hxO^L=finite_volume:ixO^L-kr(idim,^D)
  ! so that if you want to retrieve back finite_volume:ixO^L, you must do
  ! ** ixOfinitemax^D=ixOmax^D=finite_volume:ixCmax^D=finite_volume:ixOmax^D ;
  ! ** ixOfinitemin^D=ixOmin^D+kr(idim,^D)=finite_volume:ixCmin^D+kr(idim,^D)
  !                  =finite_volume:hxOmin^D+kr(idim,^D)
  !                  =finite_volume:ixOmin^D-kr(idim,^D)+kr(idim,^D)
  !                  =finite_volume:ixOmin^D
  !/



  function hd_kin_en(w, ixI^L, ixO^L, inv_rho) result(ke)
    use mod_global_parameters, only: nw, ndim
    integer, intent(in)            :: ixI^L, ixO^L
    real(dp), intent(in)           :: w(ixI^S, nw)
    real(dp)                       :: ke(ixO^S)
    real(dp), intent(in), optional :: inv_rho(ixO^S)

    if (present(inv_rho)) then
       ke = 0.5_dp * sum(w(ixO^S, mom(:))**2.0_dp, dim=ndim+1) * inv_rho
    else
       ke = 0.5_dp * sum(w(ixO^S, mom(:))**2.0_dp, dim=ndim+1) / w(ixO^S, rho_)
    end if
  end function hd_kin_en



  function hd_inv_rho(w, ixI^L, ixO^L) result(inv_rho)
    use mod_global_parameters, only: nw, ndim
    integer, intent(in)           :: ixI^L, ixO^L
    real(dp), intent(in)  :: w(ixI^S, nw)
    real(dp)              :: inv_rho(ixO^S)

    ! Can make this more robust
    inv_rho = 1.0_dp / w(ixO^S, rho_)
  end function hd_inv_rho

  subroutine hd_handle_small_values(primitive, w, x, ixI^L, ixO^L, subname,&
                                   flag_error)
    use mod_global_parameters
    use mod_small_values
    use mod_dust, only : dust_set_floor
    implicit none

    logical, intent(in)             :: primitive
    integer, intent(in)             :: ixI^L,ixO^L
    real(dp), intent(inout)         :: w(ixI^S,1:nw)
    real(dp), intent(in)            :: x(ixI^S,1:ndim)
    character(len=*), intent(in)    :: subname
    integer, optional, intent(in)   :: flag_error(ixI^S)

    real(dp) :: smallone
    integer :: idir, ierror,flag(ixI^S)


    if (small_values_method == "ignore") return

    call hd_check_w(primitive, ixI^L, ixO^L, w, flag)

    if (any(flag(ixO^S) /= 0)) then
      select case (small_values_method)
      case ("replace")
        call hd_small_values_floor(primitive, w, x, ixI^L, ixO^L, subname,&
                                     flag)
        if(hd_config%dust_on)call dust_set_floor(primitive,ixI^L,ixO^L,flag,x,w)
      case ("average")
        call small_values_average(ixI^L, ixO^L,subname, w, x, flag,ierror)

        if(small_values_force_floor.and.ierror/=0)then
          call hd_small_values_floor(primitive, w, x, ixI^L, ixO^L, subname,&
                                     flag)
         end if
         cond_dust : if(hd_config%dust_on)then
          cond_dust_error : if(any(flag(ixO^S)<0)) then
           call dust_set_floor(primitive,ixI^L,ixO^L,flag,x,w)
          end if cond_dust_error
         end if cond_dust
      case default
        call small_values_error(w, x, ixI^L, ixO^L, flag, subname)
      end select
    end if

  end subroutine hd_handle_small_values
  !-----------------------------------------------------------------
  !> reset floor variables
  subroutine hd_small_values_floor(primitive, w, x, ixI^L, ixO^L, subname,&
                                   flag)
    use mod_global_parameters
    use mod_small_values
    use mod_dust, only : dust_set_floor
    implicit none

    logical, intent(in)             :: primitive
    integer, intent(in)             :: ixI^L,ixO^L
    real(dp), intent(inout)         :: w(ixI^S,1:nw)
    real(dp), intent(in)            :: x(ixI^S,1:ndim)
    character(len=*), intent(in)    :: subname
    integer, intent(in)             :: flag(ixI^S)
    ! .. local ..
    real(dp)                        :: smallone
    integer                         :: idir
     if(all(flag(ixO^S) <= 0)) return
        where(flag(ixO^S) > 0) w(ixO^S,rho_) = hd_config%small_density

        if(hd_config%use_grackle)then
          where(flag(ixO^S) > 0) w(ixO^S, rhoHI_) = hd_config%small_density/14.0_dp
          where(flag(ixO^S) > 0) w(ixO^S, rhoHII_) = hd_config%small_density/14.0_dp
          where(flag(ixO^S) > 0) w(ixO^S, rhoHM_) = hd_config%small_density/14.0_dp
          where(flag(ixO^S) > 0) w(ixO^S, rhoHeI_) = hd_config%small_density/14.0_dp
          where(flag(ixO^S) > 0) w(ixO^S, rhoHeII_) = hd_config%small_density/14.0_dp
          where(flag(ixO^S) > 0) w(ixO^S, rhoHeIII_) = hd_config%small_density/14.0_dp
          where(flag(ixO^S) > 0) w(ixO^S, rhoH2I_) = hd_config%small_density/14.0_dp
          where(flag(ixO^S) > 0) w(ixO^S, rhoH2II_) = hd_config%small_density/14.0_dp
          where(flag(ixO^S) > 0) w(ixO^S, rhoDI_) = hd_config%small_density/14.0_dp
          where(flag(ixO^S) > 0) w(ixO^S, rhoDII_) = hd_config%small_density/14.0_dp
          where(flag(ixO^S) > 0) w(ixO^S, rhoHDI_) = hd_config%small_density/14.0_dp
          where(flag(ixO^S) > 0) w(ixO^S, rhoe_) = hd_config%small_density/14.0_dp
          where(flag(ixO^S) > 0) w(ixO^S, rhometal_) = hd_config%small_density/14.0_dp
          where(flag(ixO^S) > 0) w(ixO^S, rhodust_) = hd_config%small_density/14.0_dp
        end if

        do idir = 1, ndir
          where(flag(ixO^S) > 0) w(ixO^S, mom(idir)) = 0.0_dp
        end do

        if (hd_config%energy) then
          if(primitive) then
            smallone = hd_config%small_pressure
          else
            smallone = small_e
          end if
          where(flag(ixO^S) > 0) w(ixO^S,e_) = smallone
        end if
        if(hd_config%dust_on)call dust_set_floor(primitive,ixI^L,ixO^L,flag,x,w)
  end subroutine hd_small_values_floor


  !-----------------------------------------------------------------
  !> handle small pressure
  subroutine hd_handle_small_values_pressure(ixI^L,ixO^L,x,subname,pth,w)
    use mod_global_parameters
    use mod_small_values
    implicit none
    integer, intent(in)             :: ixI^L,ixO^L
    real(dp), intent(in)            :: x(ixI^S, 1:ndim)
    real(dp), intent(inout)         :: pth(ixI^S)
    real(dp), intent(in)            :: w(ixI^S,1:nw)
    character(len=*), intent(in)    :: subname
    !.. local ..
    real(dp)                        :: smallone
    integer                         :: idir, flag(ixI^S)
    !-------------------------------------------------------
    if (small_values_method == "ignore") return

    call hd_check_pth( ixI^L, ixO^L, pth, flag)

    if (any(flag(ixO^S) /= 0)) then
      select case (small_values_method)
      case ("replace")
        where(flag(ixO^S) > 0) pth(ixO^S) = hd_config%small_pressure
      case ("average")
        call small_var_average(ixI^L, ixO^L, pth, x, flag)
      case default
        call small_values_error(w, x, ixI^L, ixO^L, flag, subname)
      end select
    end if
  end subroutine hd_handle_small_values_pressure
end module mod_hd_phys
