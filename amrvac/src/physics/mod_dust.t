!> Module for including dust species, which interact with the gas through a drag
!> force
module mod_dust
  use mod_global_parameters, only: std_len,dp
  use mod_physics

  implicit none
  private

  !> The number of dust species
  integer, public, protected      :: dust_n_species = 0

  integer, protected              :: gas_rho_ = -1
  integer, allocatable, protected :: gas_mom(:)
  integer, protected              :: gas_e_   = -1

  !> Mean molecular weight of gas molecules
  real(kind=dp), protected, public :: gas_mu = -huge(1.0d0)

  !> Indices of the dust densities
  integer, allocatable, public, protected :: dust_rho(:)

  !> Indices of the dust momentum densities
  integer, allocatable, public, protected :: dust_mom(:, :)

  !> Size of each dust species
  real(kind=dp), allocatable, public :: dust_size(:)

  !> Internal density of each dust species
  real(kind=dp), allocatable, public :: dust_density(:)

  !> Dust temperature (if dust_temperature_type is constant)
  real(kind=dp) :: dust_temperature = -1.0d0

  !> If dust_temperature_type is stellar, it will be calculated according to Tielens (2005),
  !> eqn. 5.44 using a stellar luminosity in solar luminosities
  real(kind=dp) :: dust_stellar_luminosity = -1.0d0

  !> Set small dust densities to zero to avoid numerical problems
  logical :: dust_small_to_zero = .false.

  !> Minimum dust density
  real(kind=dp) :: dust_min_rho = -1.0d0

  !> TODO: 1. Introduce this generically in advance, 2: document
  logical :: dust_source_split = .false.

  !> What type of dust drag force to use. Can be 'Kwok', 'sticking', 'linear',or 'none'.
  character(len=std_len) :: dust_method = 'Kwok'

  !> Can be 'graphite' or 'silicate', affects the dust temperature
  character(len=std_len) :: dust_species = 'graphite'

  !> Determines the dust temperature, can be 'constant', 'ism', or 'stellar'
  character(len=std_len) :: dust_temperature_type = 'constant'
  !>
  integer                :: dust_it_diff=-1
  real(kind=dp)          :: dust_dt_small=1.0d-18


  type dust_parameters
    !> The number of dust species
    integer      :: n_species

    !> Mean molecular weight of gas molecules
    real(kind=dp) :: gas_mu


      !> Dust temperature (if dust_temperature_type is constant)
      real(kind=dp) :: temperature

      !> If dust_temperature_type is stellar, it will be calculated according to Tielens (2005),
      !> eqn. 5.44 using a stellar luminosity in solar luminosities
      real(kind=dp) :: stellar_luminosity

      !> Set small dust densities to zero to avoid numerical problems
      logical :: small_to_zero

      !> Minimum dust density
      real(kind=dp) :: min_rho

      !> TODO: 1. Introduce this generically in advance, 2: document
      logical :: source_split

      !> What type of dust drag force to use. Can be 'Kwok', 'sticking', 'linear',or 'none'.
      character(len=std_len) :: method

      !> Can be 'graphite' or 'silicate', affects the dust temperature
      character(len=std_len) :: species

      !> Determines the dust temperature, can be 'constant', 'ism', or 'stellar'
      character(len=std_len) :: temperature_type
      !>
      integer                :: it_diff
      real(kind=dp)          :: dt_small
  end type dust_parameters

  type dust_phys
     !> Size of each dust species
     real(kind=dp), allocatable :: dsize(:)

     !> Internal density of each dust species
     real(kind=dp), allocatable:: density(:)
    type(dust_parameters)  :: myconfig
  end type dust_phys
  type(dust_phys)          :: dust_inuse
  type(dust_parameters)    :: dust_config


  real(kind=dp), private   :: dust_temperature_norma_coef
  ! Public methods
  public :: dust_init
  public :: dust_physical_units
  public :: dust_get_dt
  public :: dust_get_flux
  public :: dust_get_cmax
  public :: dust_get_flux_prim
  public :: dust_get_cmax_prim
  public :: dust_add_source
  public :: dust_to_conserved
  public :: dust_to_primitive
  public :: dust_check_params
  public :: dust_check_w
  public :: dust_set_floor
  public :: dust_average_drag_force
  public :: dust_average_dspeed
  public :: dust_average_dustdensity
  public :: get_3d_dragforce
contains

  subroutine dust_init(phys_indices_inuse,phys_config_inuse,g_rho, g_mom, g_energy)
    use mod_global_parameters
    use mod_physics
    type(phys_variables_indices)    :: phys_indices_inuse
    type(physconfig)                :: phys_config_inuse
    integer, intent(in)             :: g_rho
    integer, intent(in)             :: g_mom(ndir)
    integer, intent(in)             :: g_energy ! Negative value if not present
    integer                         :: idust, idir
    character(len=2)                :: dim
    !---------------------------------------------------------

    call dust_set_default(dust_inuse)
    call dust_read_params(par_files)



    allocate(gas_mom(ndir))
    gas_rho_ = g_rho
    gas_mom  = g_mom
    gas_e_   = g_energy

    dust_temperature_norma_coef = (unit_velocity**2.0_dp*mH_cgs) / &
         (kB_cgs)

    allocate(dust_size(dust_inuse%myconfig%n_species))
    allocate(dust_density(dust_inuse%myconfig%n_species))
    dust_size(:)    = -1.0d0
    dust_density(:) = -1.0d0

    allocate(dust_rho(dust_inuse%myconfig%n_species))
    allocate(dust_mom(ndir, dust_inuse%myconfig%n_species))

    ! Set index of dust densities
    Loop_idust1 : do idust = 1, dust_inuse%myconfig%n_species
      dust_rho(idust) = var_set_fluxvar("rhod", "rhod", idust)
    end do Loop_idust1

    ! Dust momentum
    Loop_idir1 : do idir = 1, ndir
      write(dim, "(I0,A)") idir, "d"
      Loop_idust2 : do idust = 1, dust_n_species
        dust_mom(idir, idust) = var_set_fluxvar("m"//dim, "v"//dim, idust)
      end do Loop_idust2
    end do Loop_idir1



    phys_config_inuse%dust_n_species     = dust_inuse%myconfig%n_species
    phys_config_inuse%dust_small_density = dust_inuse%myconfig%min_rho

    allocate(phys_indices_inuse%dust_rho(dust_n_species),&
             phys_indices_inuse%dust_mom(ndir, dust_n_species))
    phys_indices_inuse%dust_rho(:)   =  dust_rho(:)
    phys_indices_inuse%dust_mom(:,:) =  dust_mom(:,:)
  end subroutine dust_init

  !> set default values for dust configuration parameters
  subroutine dust_set_default(self)
    type(dust_phys)  :: self

        self%myconfig%n_species          = 0
        self%myconfig%min_rho            = -1.0d0
        self%myconfig%small_to_zero      = .false.
        self%myconfig%temperature        = -1.0d0
        self%myconfig%temperature_type   = 'constant'
        self%myconfig%species            = 'graphite'
        self%myconfig%source_split       = .false.
        self%myconfig%method             = 'Kwok'
        self%myconfig%stellar_luminosity = -1.0d0
        self%myconfig%it_diff            = -1
        self%myconfig%dt_small           = 0.0_dp
        self%myconfig%gas_mu             = huge(1.0_dp)
  end subroutine dust_set_default

  !> Read this module s parameters from a file
  subroutine dust_read_params(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: i_file

    namelist /dust_list/ dust_n_species, dust_min_rho, gas_mu, dust_method, &
         dust_small_to_zero, dust_source_split, dust_temperature, &
         dust_temperature_type,dust_it_diff,dust_dt_small

    do i_file = 1, size(files)
      open(unitpar, file=trim(files(i_file)), status="old")
      read(unitpar, dust_list, end=111)
111   close(unitpar)
    end do
        dust_inuse%myconfig%n_species          = dust_n_species
        dust_inuse%myconfig%min_rho            = dust_min_rho
        dust_inuse%myconfig%small_to_zero      = dust_small_to_zero
        dust_inuse%myconfig%temperature        = dust_temperature
        dust_inuse%myconfig%temperature_type   = dust_temperature_type
        dust_inuse%myconfig%species            = dust_species
        dust_inuse%myconfig%source_split       = dust_source_split
        dust_inuse%myconfig%method             = dust_method
        dust_inuse%myconfig%stellar_luminosity = dust_stellar_luminosity
        dust_inuse%myconfig%it_diff            = dust_it_diff
        dust_inuse%myconfig%dt_small           = dust_dt_small
        dust_inuse%myconfig%gas_mu             = gas_mu
  end subroutine dust_read_params


  subroutine dust_check_params()
    if (dust_inuse%myconfig%gas_mu <= 0.0d0) call mpistop ("Dust error: gas_mu (molecular weight)"//&
         "negative or not set")

    if (dust_inuse%myconfig%temperature_type == "constant") then
       if (dust_inuse%myconfig%temperature < 0.0d0) then
          call mpistop("Dust error: dust_temperature < 0 is not set")
       end if
    else if (dust_inuse%myconfig%temperature_type == "stellar") then
       if (dust_inuse%myconfig%stellar_luminosity < 0.0d0) then
          call mpistop("Dust error: dust_stellar_luminosity < 0 is not set")
       end if
    end if

    if (any(dust_size < 0.0d0))then
         write(*,*)' in dust_check_params At the mod_dust, the dust_size<0 '
         write(*,*)' dust_size = ',dust_size
         write(*,*)' The code will stop'
         call mpistop("Dust error: any(dust_size < 0) is not set")
    end if
    if (any(dust_density < 0.0d0)) &
         call mpistop("Dust error: any(dust_density < 0) is not set")
  end subroutine dust_check_params


  subroutine dust_physical_units
    use mod_global_parameters
    use mod_physics
    integer   :: idust
    Loop_idust0 : do idust = 1, dust_inuse%myconfig%n_species
     w_convert_factor(dust_rho(idust))     = unit_density
     if(saveprim)then
      w_convert_factor(dust_mom(:,idust))  = unit_velocity
     else
      w_convert_factor(dust_mom(:,idust))  = unit_density*unit_velocity
     end if
    end do  Loop_idust0
  end subroutine dust_physical_units

  subroutine dust_to_conserved(ixI^L, ixO^L, w, x)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    real(kind=dp), intent(inout)    :: w(ixI^S, nw)
    real(kind=dp), intent(in)       :: x(ixI^S, 1:ndim)
    integer                         :: idust, idir

    do idust = 1, dust_inuse%myconfig%n_species
      ! Convert velocity to momentum
      do idir = 1, ndir
        w(ixO^S, dust_mom(idir, idust)) = w(ixO^S, dust_rho(idust)) * &
             w(ixO^S, dust_mom(idir, idust))
      end do
    end do
  end subroutine dust_to_conserved

  subroutine dust_to_primitive(ixI^L, ixO^L, w, x)
    use mod_global_parameters

    integer, intent(in)              :: ixI^L, ixO^L
    real(kind=dp), intent(inout)     :: w(ixI^S, nw)
    real(kind=dp), intent(in)        :: x(ixI^S, 1:ndim)
    integer                          :: idust, idir

    Loop_idust : do idust = 1, dust_inuse%myconfig%n_species
      ! Convert momentum to velocity
      Loop_idir: do idir = 1, ndir
        where (w(ixO^S, dust_rho(idust)) > dust_inuse%myconfig%min_rho)
          w(ixO^S, dust_mom(idir, idust)) = w(ixO^S, dust_mom(idir, idust)) / &
               w(ixO^S, dust_rho(idust))
        elsewhere
          w(ixO^S, dust_mom(idir, idust)) = 0
        end where
      end do Loop_idir
    end do Loop_idust
  end subroutine dust_to_primitive

  subroutine dust_get_flux(w, x, ixI^L, ixO^L, idim, f)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L, idim
    real(kind=dp), intent(in)       :: w(ixI^S, 1:nw), x(ixI^S, 1:^ND)
    real(kind=dp), intent(inout)    :: f(ixI^S, nwflux)
    integer                         :: idust, idir

    Loop_idust : do idust = 1, dust_inuse%myconfig%n_species
      where (w(ixO^S, dust_rho(idust)) > dust_inuse%myconfig%min_rho)
        f(ixO^S, dust_rho(idust)) = w(ixO^S, dust_mom(idim, idust))
      elsewhere             ! TODO: remove?
        f(ixO^S, dust_rho(idust)) = 0.0d0
      end where

      Loop_idir: do idir = 1, ndir
        f(ixO^S, dust_mom(idir, idust)) = w(ixO^S, dust_mom(idir, idust)) * &
             get_vdust(w, ixI^L, ixO^L, idim, idust)
      end do Loop_idir
    end do Loop_idust
  end subroutine dust_get_flux

  subroutine dust_get_flux_prim(w, x, ixI^L, ixO^L, idim, f)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L, idim
    real(kind=dp), intent(in)       :: w(ixI^S, 1:nw), x(ixI^S, 1:^ND)
    real(kind=dp), intent(inout)    :: f(ixI^S, nwflux)
    integer                         :: idust, idir

    do idust = 1, dust_inuse%myconfig%n_species
      where (w(ixO^S, dust_rho(idust)) > dust_inuse%myconfig%min_rho)
        f(ixO^S, dust_rho(idust)) = w(ixO^S, dust_mom(idim, idust))*w(ixO^S, dust_rho(idust))
      elsewhere             ! TODO: remove?
        f(ixO^S, dust_rho(idust)) = 0.0d0
      end where

      do idir = 1, ndir
        f(ixO^S, dust_mom(idir, idust)) = w(ixO^S, dust_mom(idir, idust)) * &
        w(ixO^S, dust_rho(idust)) * get_vdust_prim(w, ixI^L, ixO^L, idim, idust)
      end do
    end do
  end subroutine dust_get_flux_prim

  function get_vdust(w, ixI^L, ixO^L, idim, idust) result(vdust)
    use mod_global_parameters, only: nw
    integer, intent(in)          :: ixI^L, ixO^L, idim, idust
    real(kind=dp), intent(in)    :: w(ixI^S, nw)
    real(kind=dp)                :: vdust(ixO^S)

    where (w(ixO^S, dust_rho(idust)) > dust_inuse%myconfig%min_rho)
      vdust = w(ixO^S, dust_mom(idim, idust)) / w(ixO^S, dust_rho(idust))
    elsewhere
      vdust                       = 0.0d0
    end where
  end function get_vdust

  function get_vdust_prim(w, ixI^L, ixO^L, idim, idust) result(vdust)
    use mod_global_parameters, only: nw
    integer, intent(in)           :: ixI^L, ixO^L, idim, idust
    real(kind=dp), intent(in)     :: w(ixI^S, nw)
    real(kind=dp)                 :: vdust(ixO^S)

    where (w(ixO^S, dust_rho(idust)) > dust_inuse%myconfig%min_rho)
      vdust = w(ixO^S, dust_mom(idim, idust))
    elsewhere
      vdust = 0.0d0
    end where
  end function get_vdust_prim

  ! Force dust density to zero if dust_rho <= dust_min_rho
  subroutine set_dusttozero(qdt, ixI^L, ixO^L,  wCT,  w, x)
    use mod_global_parameters

    integer, intent(in)           :: ixI^L, ixO^L
    real(kind=dp), intent(in)     :: qdt
    real(kind=dp), intent(in)     :: wCT(ixI^S, 1:nw), x(ixI^S, 1:ndim)
    real(kind=dp), intent(inout)  :: w(ixI^S, 1:nw)

    integer                       :: n, idir
    integer                       :: flag(ixI^S)
    flag(ixO^S) = 0
    call dust_check_w(.false., ixI^L, ixO^L, flag, w)
    call dust_set_floor(.false.,ixI^L,ixO^L,flag,x,w)
  end subroutine set_dusttozero

  ! Calculate drag force based on Epstein's or Stokes' law
  ! From Kwok 1975, page 584 (between eqn 8 and 9)
  subroutine get_3d_dragforce(ixI^L, ixO^L, w, x, fdrag, ptherm, vgas,fd_flag)
    use mod_global_parameters

    integer, intent(in)              :: ixI^L, ixO^L
    real(kind=dp), intent(in)        :: x(ixI^S, 1:ndim)
    real(kind=dp), intent(in)        :: w(ixI^S, 1:nw)
    real(kind=dp), intent(out)      :: &
         fdrag(ixI^S, 1:ndir, 1:dust_inuse%myconfig%n_species)
    real(kind=dp), intent(in)        :: ptherm(ixI^S), vgas(ixI^S, ndir)
    logical, intent(out)            :: fd_flag(ixI^S)
    real(kind=dp), dimension(ixI^S) :: vt2, deltav, fd, vdust
    real(kind=dp)                   :: alpha_T(ixI^S, 1:dust_inuse%myconfig%n_species)
    integer                         :: idust, idir, it_diff
    real(kind=dp)                   :: K,coef
    logical, save                   :: it_dustreset
    integer, save                   :: it_dustsave =-100000
    !---------------------------------
    vt2(ixO^S) = 3.0d0*ptherm(ixO^S)/w(ixO^S, gas_rho_)

    select case( TRIM(dust_inuse%myconfig%method) )
    case ('Kwok') ! assume sticking coefficient equals 0.25
     if(dt>dust_inuse%myconfig%dt_small.and.it>it_dustsave+10) then
       it_diff       = dust_inuse%myconfig%it_diff
       it_dustreset  = .true.
       it_dustsave   = -100
      else
        if(it_dustreset) then
         it_dustsave  = it
         it_dustreset = .false.
       end if
       it_diff=100000
      end if
      if(dust_inuse%myconfig%it_diff>0)coef=(real(it,kind=dp)/real(it_diff,kind=dp))**4.0_dp


      Loop_idir1 : do idir = 1, ndir

        Loop_idust1: do idust = 1, dust_inuse%myconfig%n_species
          where(w(ixO^S, dust_rho(idust)) > dust_inuse%myconfig%min_rho)
            vdust(ixO^S)  = w(ixO^S, dust_mom(idir, idust)) / w(ixO^S, dust_rho(idust))
            deltav(ixO^S) = (vgas(ixO^S, idir)-vdust(ixO^S))

             ! 0.75 from sticking coefficient
             fd(ixO^S)     = 0.75d0*w(ixO^S, dust_rho(idust))*w(ixO^S, gas_rho_)*deltav(ixO^S) &
                  / (dust_density(idust) * dust_size(idust))

             ! 0.75 from spherical grainvolume
             fd(ixO^S)     = -fd(ixO^S)*0.75d0*dsqrt(vt2(ixO^S) + deltav(ixO^S)**2.0_dp)

          elsewhere
            fd(ixO^S) = 0.0d0
          end where
          call  dust_handle_largevariation_dragforce(.false.,ixI^L,ixO^L,&
                                                         fd,fd_flag,w)

         !  test_todelete : if(it<=it_diff)then
        !    fd(ixO^S) = coef*fd(ixO^S)
         ! ! elseif(it<=it_diff.and.it>it_diff/10.0)then
         ! !   where(.not.fd_flag(ixO^S))
         ! !     fd(ixO^S) = coef *fd(ixO^S)
         ! !   elsewhere
         ! !      fd(ixO^S) = coef**0.99*fd(ixO^S)
         ! !   end where
         !  end if   test_todelete
          !call dust_handle_largevariation_dragforce(.false.,ixI^L,ixO^L,fd,fd_flag,w)
          ! testttt

          if(it<dust_inuse%myconfig%it_diff)fd(ixO^S) = coef*fd(ixO^S)
          fdrag(ixO^S, idir, idust) = fd(ixO^S)


        end do Loop_idust1

      end do Loop_idir1

    case ('sticking') ! Calculate sticking coefficient based on the gas and dust temperatures
      !  Equation from Decin et al. 2006
      if (gas_e_ < 0) call mpistop("dust sticking requires gas energy")

      call dust_get_sticking(w, x, ixI^L, ixO^L, alpha_T, ptherm)

      Loop_idir2 : do idir = 1, ndir
        do idust = 1, dust_inuse%myconfig%n_species
          where(w(ixO^S, dust_rho(idust))>dust_inuse%myconfig%min_rho)
            vdust(ixO^S)  = w(ixO^S,dust_mom(idir, idust)) / w(ixO^S, dust_rho(idust))
            deltav(ixO^S) = (vgas(ixO^S, idir)-vdust(ixO^S))
            fd(ixO^S)     = (one-alpha_T(ixO^S,idust)) * w(ixO^S, dust_rho(idust))*w(ixO^S, gas_rho_) * &
                 deltav(ixO^S) / (dust_density(idust)*dust_size(idust))
            fd(ixO^S)     = -fd(ixO^S) * 0.75d0 * dsqrt(vt2(ixO^S) + deltav(ixO^S)**2)
          else where
            fd(ixO^S) = 0.0d0
          end where
          fdrag(ixO^S, idir,idust) = fd(ixO^S)
        end do
      end do Loop_idir2
    case('linear') !linear with Deltav, for testing (see Laibe & Price 2011)
      K = 3.4d5 / dust_inuse%myconfig%n_species
      do idir = 1, ndir
        do idust = 1, dust_inuse%myconfig%n_species
          where(w(ixO^S, dust_rho(idust))>dust_inuse%myconfig%min_rho)
            vdust(ixO^S)  = w(ixO^S,dust_mom(idir, idust))/w(ixO^S, dust_rho(idust))
            deltav(ixO^S) = (vgas(ixO^S, idir)-vdust(ixO^S))

            fd(ixO^S)     = -K*deltav(ixO^S)
          else where
            fd(ixO^S) = 0.0d0
          end where

          fdrag(ixO^S, idir,idust) = fd(ixO^S)
        end do
      end do
    case('none')
      fdrag(ixO^S, :, :) = 0.0d0
    case default
      call mpistop( "=== This dust method has not been implemented===" )
    end select



  end subroutine get_3d_dragforce

  !> Get sticking coefficient
  !>
  !> Assume cgs units, and use of convert factors for conversion
  !> Equation from Decin et al. 2006
  subroutine dust_get_sticking(w, x, ixI^L, ixO^L, alpha_T, ptherm)
    use mod_global_parameters
    integer, intent(in)            :: ixI^L, ixO^L
    real(kind=dp), intent(in)      :: x(ixI^S, 1:ndim)
    real(kind=dp), intent(in)      :: w(ixI^S, 1:nw)
    real(kind=dp), intent(out)     :: alpha_T(ixI^S, 1:dust_inuse%myconfig%n_species)
    real(kind=dp), intent(in)      :: ptherm(ixI^S)
    real(kind=dp)                  :: Tgas(ixI^S)
    integer                        :: idust

    call get_tdust(w, x, ixI^L, ixO^L, alpha_T)

    Tgas(ixO^S) = (ptherm(ixO^S) / w(ixO^S, gas_rho_))*&
        dust_temperature_norma_coef

    do idust = 1, dust_inuse%myconfig%n_species
      alpha_T(ixO^S,idust) =  max(0.35d0 * exp(-sqrt((Tgas(ixO^S) + &
           alpha_T(ixO^S,idust))/5.0d2))+0.1d0, smalldouble)
    end do
  end subroutine dust_get_sticking

  !> Returns dust temperature (in K), either as constant or based on equ. 5.41,
  !> 5.42 and 5.44 from Tielens (2005)
  !>
  !> Note that this calculation assumes cgs!!!! with conversion between physical
  !> and scaled quantities done through the convert factors!!!!
  !>
  !> It takes as input the stellar luminosoity in solar units and/or a fixed
  !> dust temperature in Kelvin
  subroutine get_tdust(w, x, ixI^L, ixO^L, Td)
    use mod_global_parameters

    integer, intent(in)         :: ixI^L, ixO^L
    real(kind=dp), intent(in)   :: x(ixI^S, 1:ndim)
    real(kind=dp), intent(in)   :: w(ixI^S, 1:nw)
    real(kind=dp), intent(out)  :: Td(ixI^S, 1:dust_inuse%myconfig%n_species)
    real(kind=dp)               :: G0(ixO^S)
    integer                     :: idust
    !-------------------------------------------------------------------

    select case( trim(dust_inuse%myconfig%temperature_type) )
    case( 'constant' )
      Td(ixO^S, :) = dust_inuse%myconfig%temperature
    case( 'ism' )

      select case( trim(dust_species) )
      case( 'graphite' )
        do idust = 1, dust_inuse%myconfig%n_species
          Td(ixO^S, idust) = 15.8d0*((0.0001d0/(dust_size(idust)*length_convert_factor))**0.06d0)
        end do
      case( 'silicate' )
        do idust = 1, dust_inuse%myconfig%n_species
          Td(ixO^S, idust) = 13.6d0*((0.0001d0/(dust_size(idust)*length_convert_factor))**0.06d0)
        end do
      case default
        call mpistop( "=== Dust species undetermined===" )
      end select

    case( 'stellar' )

      select case( trim(typeaxial) )
      case( 'spherical' )
        G0(ixO^S) = max(x(ixO^S, 1)*length_convert_factor, smalldouble)
      case( 'cylindrical' )
        G0(ixO^S) = max(dsqrt(sum(x(ixO^S,:)**2,dim=ndim+1))*length_convert_factor, smalldouble)
      end select

      G0(ixO^S) = 2.1d4*(dust_inuse%myconfig%stellar_luminosity/1.0d8)*((3.0857d17/G0(ixO^S))**2)

      select case( trim(dust_species) )
      case( 'graphite' )
        do idust = 1, dust_inuse%myconfig%n_species
          Td(ixO^S, idust) = 61.0d0*((0.0001d0/(dust_size(idust)*length_convert_factor))**0.06d0) &
               *(G0(ixO^S)**(one/5.8d0))
        end do
      case( 'silicate' )
        do idust = 1, dust_inuse%myconfig%n_species
          Td(ixO^S, idust) = 50.0d0*((0.0001d0/(dust_size(idust)*length_convert_factor))**0.06d0) &
               *(G0(ixO^S)**(one/6.0d0))
        end do
      case default
        call mpistop( "=== Dust species undetermined===" )
      end select

    case default
      write(*,*)' The dust temperature type is : ' ,trim(dust_inuse%myconfig%temperature_type)
      call mpistop( "=== Dust temperature undetermined===" )
    end select

  end subroutine get_tdust

  !> w[iw]= w[iw]+qdt*S[wCT,  x] where S is the source based on wCT within ixO
  subroutine dust_add_source(qdt, ixI^L, ixO^L, wCT,w, x, qsourcesplit, active)
    use mod_global_parameters

    integer, intent(in)              :: ixI^L, ixO^L
    real(kind=dp), intent(in)        :: qdt
    real(kind=dp), intent(in)        :: wCT(ixI^S, 1:nw), x(ixI^S, 1:ndim)
    real(kind=dp), intent(inout)     :: w(ixI^S, 1:nw)
    logical, intent(in)              :: qsourcesplit
    logical, intent(inout)           :: active

    ! .. local ..
    logical          :: fd_flag(ixI^S)
    real(kind=dp)    :: ptherm(ixI^S), vgas(ixI^S, ndir),vdust(ixI^S),dv(ixI^S)
    real(kind=dp)    :: fdrag(ixI^S, ndir, dust_inuse%myconfig%n_species)
    integer          :: idust,idir, it_diff
    !------------------------------------------------------------
    select case( TRIM(dust_inuse%myconfig%method) )
    case( 'none' )
      !do nothing here
    case default !all regular dust methods here
      if (qsourcesplit .eqv. dust_inuse%myconfig%source_split) then
        active = .true.

        call phys_get_pthermal(wCT, x, ixI^L, ixO^L, ptherm)
        do idir=1,ndir
          vgas(ixO^S,idir)=wCT(ixO^S,gas_mom(idir))/wCT(ixO^S,gas_rho_)
        end do

        call get_3d_dragforce(ixI^L, ixO^L, wCT, x, fdrag, ptherm, vgas,fd_flag)

        fdrag(ixO^S,:, :) = fdrag(ixO^S,:, :) * qdt


        Loop_idir : do idir = 1, ndir

          Loop_idust : do idust = 1, dust_inuse%myconfig%n_species
            w(ixO^S, gas_mom(idir))  = w(ixO^S, gas_mom(idir)) + &
                 fdrag(ixO^S, idir, idust)

            if (gas_e_ > 0) then
              w(ixO^S, gas_e_) = w(ixO^S, gas_e_) + (wCT(ixO^S, gas_mom(idir)) / &
                   wCT(ixO^S, gas_rho_)) * fdrag(ixO^S, idir, idust)
            end if

            w(ixO^S, dust_mom(idir, idust)) = w(ixO^S, dust_mom(idir, idust)) - &
                 fdrag(ixO^S, idir, idust)
          end do Loop_idust
        end do Loop_idir
        ! it_diff=2000000
        ! test_todelete : if(it<=it_diff)then
        !
        !   where(w(ixO^S,dust_rho(idust))>dust_min_rho)
        !     vdust(ixO^S) = w(ixO^S, dust_mom(idir,idust))/w(ixO^S,dust_rho(idust))
        !     dv(ixO^S) = vgas(ixO^S,idir)-vdust(ixO^S)
        !     w(ixO^S,dust_mom(idir,idust))  = (w(ixO^S, dust_mom(idir,idust))  +&
        !            (dv(ixO^S))*(it-it_diff)/it_diff)*w(ixO^S,dust_rho(idust))
        !   end where
        ! end if test_todelete

        if (dust_small_to_zero) then
          call set_dusttozero(qdt, ixI^L, ixO^L,  wCT,  w, x)
        end if
      endif
    end select

  end subroutine dust_add_source

  !> Get dt related to dust and gas stopping time (Laibe 2011)
  subroutine dust_get_dt(w, ixI^L, ixO^L, dtnew, dx^D, x)
    use mod_global_parameters

    integer, intent(in)              :: ixI^L, ixO^L
    real(kind=dp), intent(in)        :: dx^D, x(ixI^S, 1:^ND)
    real(kind=dp), intent(in)        :: w(ixI^S, 1:nw)
    real(kind=dp), intent(inout)     :: dtnew
    ! .. local ..
    real(kind=dp)                              :: vgas(ixI^S, ndir)
    real(kind=dp), dimension(1:dust_inuse%myconfig%n_species):: dtdust
    real(kind=dp), dimension(ixI^S)            :: vt2, deltav, tstop, vdust,ptherm
    real(kind=dp), dimension(ixI^S,1:dust_inuse%myconfig%n_species) :: alpha_T
    real(kind=dp)                              :: K
    integer                                    :: idust, idir
    !----------------------------------------------
    call phys_get_pthermal(w, x, ixI^L, ixO^L, ptherm)
    do idir = 1, ndir
      vgas(ixO^S,idir)=w(ixO^S,gas_mom(idir))/w(ixO^S,gas_rho_)
    end do
    select case( TRIM(dust_inuse%myconfig%method) )

    case( 'Kwok' ) ! assume sticking coefficient equals 0.25
      dtdust(:) = bigdouble

      vt2(ixO^S) = 3.0d0*ptherm(ixO^S)/w(ixO^S, gas_rho_)

      ! Tgas, mu = mean molecular weight
      ptherm(ixO^S) = ( ptherm(ixO^S) * w_convert_factor(gas_e_) * &
           mH_cgs*dust_inuse%myconfig%gas_mu) / (w(ixO^S, gas_rho_) * &
           w_convert_factor(gas_rho_)*kB_cgs)

      do idir = 1, ndir
        do idust = 1, dust_inuse%myconfig%n_species
          where(w(ixO^S, dust_rho(idust))>dust_inuse%myconfig%min_rho)
            vdust(ixO^S)  = w(ixO^S,dust_mom(idir, idust))/w(ixO^S, dust_rho(idust))
            deltav(ixO^S) = (vgas(ixO^S, idir)-vdust(ixO^S))
            tstop(ixO^S)  = 4.0d0*(dust_density(idust)*dust_size(idust))/ &
                 (3.0d0*(0.75d0)*dsqrt(vt2(ixO^S) + &
                 deltav(ixO^S)**2)*(w(ixO^S, dust_rho(idust)) + &
                 w(ixO^S, gas_rho_)))
          else where
            tstop(ixO^S) = bigdouble
          end where
        end do
      end do

      dtnew = min(minval(dtdiffpar*dtdust(:)), dtnew)

    case( 'sticking' ) ! Calculate sticking coefficient based on the gas temperature
      dtdust(:) = bigdouble

      vt2(ixO^S) = 3.0d0*ptherm(ixO^S)/w(ixO^S, gas_rho_)

      ! Sticking coefficient
      call dust_get_sticking(w, x, ixI^L, ixO^L, alpha_T, ptherm)

      ! Tgas, mu = mean molecular weight
      ptherm(ixO^S) = ( ptherm(ixO^S)*w_convert_factor(gas_e_) * &
           mH_cgs*dust_inuse%myconfig%gas_mu) / (w(ixO^S, gas_rho_) * &
           w_convert_factor(gas_rho_)*kB_cgs)

      do idir = 1, ndir
        do idust = 1, dust_inuse%myconfig%n_species
          where(w(ixO^S, dust_rho(idust))>dust_inuse%myconfig%min_rho)
            vdust(ixO^S)  = w(ixO^S,dust_mom(idir, idust))/w(ixO^S, dust_rho(idust))
            deltav(ixO^S) = (vgas(ixO^S, idir)-vdust(ixO^S))
            tstop(ixO^S)  = 4.0d0*(dust_density(idust)*dust_size(idust))/ &
                 (3.0d0*(one-alpha_T(ixO^S,idust))*dsqrt(vt2(ixO^S) + &
                 deltav(ixO^S)**2)*(w(ixO^S, dust_rho(idust)) + &
                 w(ixO^S, gas_rho_)))
          else where
            tstop(ixO^S) = bigdouble
          end where

          dtdust(idust) = min(minval(tstop(ixO^S)), dtdust(idust))
        end do
      end do

      dtnew = min(minval(dtdiffpar*dtdust(:)), dtnew)

    case('linear') !linear with Deltav, for testing (see Laibe & Price 2011)
      K = 3.4d5/dust_inuse%myconfig%n_species
      dtdust(:) = bigdouble

      do idust = 1, dust_inuse%myconfig%n_species
        where(w(ixO^S, dust_rho(idust))>dust_inuse%myconfig%min_rho)
          tstop(ixO^S)  = (w(ixO^S, dust_rho(idust))*w(ixO^S, gas_rho_))/ &
               (K*(w(ixO^S, dust_rho(idust)) + w(ixO^S, gas_rho_)))
        else where
          tstop(ixO^S) = bigdouble
        end where

        dtdust(idust) = min(minval(tstop(ixO^S)), dtdust(idust))
      end do

      dtnew = min(minval(dtdiffpar*dtdust(:)), dtnew)
    case('none')
      ! no dust timestep
    case default
      call mpistop( "=== This dust method has not been implemented===" )
    end select

    if (dtnew < dtmin) then

      write(unitterm,*)"-------------------------------------"
      write(unitterm,*)"Warning: found DUST related time step too small! dtnew=", dtnew
      write(unitterm,*)"on grid with index:", saveigrid," grid level=", node(plevel_, saveigrid)
      write(unitterm,*)"grid corners are=",{^D&rnode(rpxmin^D_, saveigrid), rnode(rpxmax^D_, saveigrid)}
      write(unitterm,*)" dtdust =", dtdust(:)
      write(unitterm,*)"on processor:", mype
      write(unitterm,*)"-------------------------------------"
    endif

  end subroutine dust_get_dt

  ! Note that cmax and cmin are assumed to be initialized
  subroutine dust_get_cmax(w, x, ixI^L, ixO^L, idim, cmax, cmin)
    use mod_global_parameters

    integer, intent(in)                     :: ixI^L, ixO^L, idim
    real(kind=dp), intent(in)               :: w(ixI^S, nw), x(ixI^S, 1:^ND)
    real(kind=dp), intent(inout)            :: cmax(ixI^S)
    real(kind=dp), intent(inout), optional  :: cmin(ixI^S)
    real(kind=dp)                           :: vdust(ixI^S)
    integer                                 :: idust

    do idust = 1, dust_inuse%myconfig%n_species
      vdust(ixO^S) = get_vdust(w, ixI^L, ixO^L, idim, idust)

      if (present(cmin)) then
        cmin(ixO^S) = min(cmin(ixO^S), vdust(ixO^S))
        cmax(ixO^S) = max(cmax(ixO^S), vdust(ixO^S))
      else
        cmax(ixO^S) = max(cmax(ixO^S), abs(vdust(ixO^S)))
      end if
    end do
  end subroutine dust_get_cmax

  ! Note that cmax and cmin are assumed to be initialized
  subroutine dust_get_cmax_prim(w, x, ixI^L, ixO^L, idim, cmax, cmin)
    use mod_global_parameters

    integer, intent(in)                     :: ixI^L, ixO^L, idim
    real(kind=dp), intent(in)               :: w(ixI^S, nw), x(ixI^S, 1:^ND)
    real(kind=dp), intent(inout)            :: cmax(ixI^S)
    real(kind=dp), intent(inout), optional  :: cmin(ixI^S)
    real(kind=dp)                           :: vdust(ixI^S)
    integer                                 :: idust

    do idust = 1, dust_inuse%myconfig%n_species
      vdust(ixO^S) = get_vdust_prim(w, ixI^L, ixO^L, idim, idust)
      if (present(cmin)) then
        cmin(ixO^S) = min(cmin(ixO^S), vdust(ixO^S))
        cmax(ixO^S) = max(cmax(ixO^S), vdust(ixO^S))
      else
        cmax(ixO^S) = max(cmax(ixO^S), abs(vdust(ixO^S)))
      end if
    end do
  end subroutine dust_get_cmax_prim



  !--------------------------------------------------------------------
  !> subroutine check error in dust
  subroutine dust_check_w(primitive, ixI^L, ixO^L, flag, w)
    use mod_global_parameters
    implicit none
    logical, intent(in)        :: primitive
    integer, intent(in)        :: ixI^L, ixO^L
    integer, intent(inout)     :: flag(ixI^S)
    real(dp), intent(in)       :: w(ixI^S,1:nw)
    ! .. local ..
    integer                    :: idust
    !------------------------------------
    Loop_dust :  do idust=1,dust_inuse%myconfig%n_species
      where(w(ixO^S,dust_rho(idust))<=dust_inuse%myconfig%min_rho)
        flag(ixO^S) = - dust_rho(idust)
      end where
    end do Loop_dust
  end subroutine dust_check_w
  !--------------------------------------------------------------------
  !> subroutine to handel small values
  subroutine dust_set_floor(primitive,ixI^L,ixO^L,flag,x,w)
    use mod_global_parameters
    implicit none
    logical, intent(in)        :: primitive
    integer, intent(in)        :: ixI^L, ixO^L
    integer, intent(in)        :: flag(ixI^S)
    real(dp), intent(in)       :: x(ixI^S,1:ndim)
    real(dp), intent(inout)    :: w(ixI^S,1:nw)
    ! .. local ..
    integer                    :: idust,idir
    logical, dimension(ixI^S)  :: patch_correct
    !------------------------------------
    Loop_dust :  do idust=1,dust_inuse%myconfig%n_species
      patch_correct(ixO^S) = (flag(ixO^S)==-dust_rho(idust)).or.(flag(ixO^S)>0)
      where(patch_correct(ixO^S))
        w(ixO^S,dust_rho(idust)) = 0.0_dp
      end where
      Loop_idir :  do idir=1,ndir
        where(patch_correct(ixO^S))
          w(ixO^S,dust_mom(idir,idust)) = 0.0_dp
        end where
      end do Loop_idir
    end do Loop_dust
  end subroutine dust_set_floor
  !------------------------------------------------------------------------------------------
  !> subroutine to compute large variation in drag force
  subroutine dust_mean_dragforce(ixI^L,ixO^L,forcedrag,mean_force,mean_flag)
    use mod_global_parameters
    use mod_small_values
    implicit none
    integer, intent(in)              :: ixI^L, ixO^L
    real(kind=dp), intent(in)        :: forcedrag(ixI^S)
    real(kind=dp), intent(out)       :: mean_force(ixI^S)
    logical, intent(out)             :: mean_flag(ixI^S)
    ! .. local ..
    integer                          :: ix^D,kxO^L
    real(kind=dp)                    :: dust_max_s_dragforce,sum_force
    !----------------------------------------------------------------
    dust_max_s_dragforce = 1.2d0
    {do ix^DB= ixO^LIM^DB\}
      {kxOmin^D= max(ix^D-small_values_daverage, ixOmin^D);
      kxOmax^D= min(ix^D+small_values_daverage, ixOmax^D);\}
      SUM_force=sum(forcedrag(kxO^S))
      mean_force(ix^D)=(sum_force-forcedrag(ix^D))&
        /({^D&(kxOmax^D-kxOmin^D+1)*}-1)
      if(dabs(forcedrag(ix^D))>smalldouble&
         .and.dabs(mean_force(ix^D))>smalldouble)then
          mean_flag(ix^D)=dabs((forcedrag(ix^D)-mean_force(ix^D))/mean_force(ix^D))<dust_max_s_dragforce
      else
          mean_flag(ix^D)=.true.
      end if

    {end do^D&\}

  end subroutine dust_mean_dragforce
  !------------------------------------------------------------------------------------------
  !> subroutine to handel large variation in drag force
  subroutine dust_handle_largevariation_dragforce(primitive,ixI^L,ixO^L,&
                                                 forcedrag,mean_flag,w)
   use mod_global_parameters
   implicit none
   integer, intent(in)              :: ixI^L, ixO^L
   logical, intent(in)              :: primitive
   real(kind=dp), intent(inout)     :: forcedrag(ixI^S)
   logical, intent(inout)           :: mean_flag(ixI^S)
   real(kind=dp), intent(in)        :: w(ixI^S,1:nw)

   !.. local ..
    real(kind=dp)                   :: mean_force(ixI^S)

   !----------------------
   call dust_mean_dragforce(ixI^L,ixO^L,forcedrag,mean_force,mean_flag)
   !call dust_average_drag_force(ixI^L,ixO^L,mean_flag,w)

   where(.not.mean_flag(ixO^S))
    forcedrag(ixO^S) = mean_force(ixO^S)
   end where
  end subroutine dust_handle_largevariation_dragforce

  !------------------------------------------------------------------------------------------
  !> subroutine to handel large variation in drag force
  subroutine dust_average_drag_force(ixI^L,ixO^L,idust,w_flag,w)

      use mod_global_parameters
      use mod_small_values
    integer, intent(in)              :: ixI^L, ixO^L,idust
    logical, intent(in)              :: w_flag(ixI^S)
    real(kind=dp), intent(inout)     :: w(ixI^S, 1:nw)
    integer                          :: iw, kxO^L, ix^D, i,ivar



    {do ix^DB= ixO^LIM^DB\}

    ! point with local failure identified by w_flag
    if (.not.w_flag(ix^D)) then
      ! verify in cube with border width small_values_daverage the presence of
      ! cells where all went ok
      do i = 1, max(small_values_daverage, 1)
        {kxOmin^D= max(ix^D-i, ixOmin^D);
        kxOmax^D= min(ix^D+i, ixOmax^D);\}

        ! in case cells are fine within smaller cube than
        ! the userset small_values_daverage: use that smaller cube
        if (any(w_flag(kxO^S))) exit
      end do

      if (any(w_flag(kxO^S))) then
        ! within surrounding cube, cells without problem were found

        ! faulty cells are corrected by averaging here
        ! only average those which were ok and replace faulty cells
        do ivar = 1,ndir+1
           if(ivar==1) then
             iw = dust_rho(idust)
           else
             iw =  dust_mom(ivar-1,idust)
           end if
          w(ix^D, iw) = sum(w(kxO^S, iw), w_flag(kxO^S))&
               / count(w_flag(kxO^S))
        end do



      else
        write(*,*) "no cells without error were found in cube of size", &
             small_values_daverage
        write(*,*) "at index:", ix^D
        write(*,*) "w_flag(ix^D):", w_flag(ix^D)
        write(*,*)  "w = ",w(ix^D, 1:nw)
        write(*,*) 'iteration',it
        write(*,*) "Saving status at the previous time step"
        crash=.true.
        call mpistop('is wrong')
      end if
    end if
    {enddo^D&\}
  end  subroutine dust_average_drag_force

  !------------------------------------------------------------------------------------------
  !> subroutine to handel large variation in drag force
  subroutine dust_average_dspeed(ixI^L,ixO^L,idust,w_flag,w)

      use mod_global_parameters
      use mod_small_values
    integer, intent(in)             :: ixI^L, ixO^L,idust
    logical, intent(in)             :: w_flag(ixI^S)
    real(kind=dp), intent(inout)    :: w(ixI^S, 1:nw)
    integer                         :: iw, kxO^L, ix^D, i,ivar



    {do ix^DB= ixO^LIM^DB\}

    ! point with local failure identified by w_flag
    if (.not.w_flag(ix^D)) then
      ! verify in cube with border width small_values_daverage the presence of
      ! cells where all went ok
      do i = 1, max(small_values_daverage, 1)
        {kxOmin^D= max(ix^D-i, ixOmin^D);
        kxOmax^D= min(ix^D+i, ixOmax^D);\}

        ! in case cells are fine within smaller cube than
        ! the userset small_values_daverage: use that smaller cube
        if (any(w_flag(kxO^S))) exit
      end do

      if (any(w_flag(kxO^S))) then
        ! within surrounding cube, cells without problem were found

        ! faulty cells are corrected by averaging here
        ! only average those which were ok and replace faulty cells
        do ivar = 1,ndim+1
           if(ivar==1) then
             cycle
             iw = dust_rho(idust)

             w(ix^D, iw) = sum(w(kxO^S, iw), mask =w_flag(kxO^S))&
                             / count(w_flag(kxO^S))
           else

             iw =  dust_mom(ivar-1,idust)
             if(dabs(w(ix^D, iw))>smalldouble)then
              w(ix^D, iw) = sum(w(kxO^S, iw)-w(kxO^S,gas_mom(ivar-1)),&
                               mask=w_flag(kxO^S))&
                             / count(w_flag(kxO^S))+w(ix^D, gas_mom(ivar-1))
             end if
           end if



        end do



      else
        write(*,*) "no cells without error were found in cube of size", &
             small_values_daverage
        write(*,*) "at index:", ix^D
        write(*,*) "w_flag(ix^D):", w_flag(ix^D)
        write(*,*)  "w = ",w(ix^D, 1:nw)
        write(*,*) 'iteration',it
        write(*,*) "Saving status at the previous time step"
        crash=.true.
        call mpistop('is wrong')
      end if
    end if
    {enddo^D&\}
  end  subroutine dust_average_dspeed

    !------------------------------------------------------------------------------------------
    !> subroutine to handel large variation in drag force
    subroutine dust_average_dustdensity(ixI^L,ixO^L,idust,w_flag,w)

        use mod_global_parameters
        use mod_small_values
      integer, intent(in)             :: ixI^L, ixO^L,idust
      logical, intent(in)             :: w_flag(ixI^S)
      real(kind=dp), intent(inout)    :: w(ixI^S, 1:nw)
      integer                         :: kxO^L, ix^D



      {do ix^DB= ixO^LIM^DB\}

      ! point with local failure identified by w_flag
      if (.not.w_flag(ix^D)) then
        ! verify in cube with border width small_values_daverage the presence of
        ! cells where all went ok

          {kxOmin^D= max(ix^D-small_values_daverage, ixOmin^D);
          kxOmax^D= min(ix^D+small_values_daverage, ixOmax^D);\}
        if (any(w_flag(kxO^S))) then
          ! within surrounding cube, cells without problem were found

          ! faulty cells are corrected by averaging here
          ! only average those which were ok and replace faulty cells

               if(dabs(w(ix^D, dust_rho(idust)))>smalldouble)then
                w(ix^D, dust_rho(idust)) = sum(w(kxO^S, dust_rho(idust))*w(kxO^S,gas_rho_),&
                                 mask=w_flag(kxO^S))&
                               / count(w_flag(kxO^S))/w(ix^D, gas_rho_)
               end if

        else
          write(*,*) "no cells without error were found in cube of size", &
               small_values_daverage
          write(*,*) "at index:", ix^D
          write(*,*) "w_flag(ix^D):", w_flag(ix^D)
          write(*,*)  "w = ",w(ix^D, 1:nw)
          write(*,*) 'iteration',it
          write(*,*) "Saving status at the previous time step"
          crash=.true.
          call mpistop('is wrong')
        end if
      end if
      {enddo^D&\}
    end  subroutine dust_average_dustdensity
end module mod_dust
