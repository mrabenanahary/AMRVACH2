!> Special Relativistic Magneto-hydrodynamics module
module mod_srhd_phys
  use mod_global_parameters, only: std_len
  use mod_srhd_parameters
  use mod_srhd_eos

  implicit none
  private

!  !> Whether an energy equation is used
!  logical, public, protected              :: srhd_energy = .true.

!  !> Whether thermal conduction is used
!  logical, public, protected              :: srhd_thermal_conduction = .false.

!  !> Whether radiative cooling is added
!  logical, public, protected              :: srhd_radiative_cooling = .false.


!  !> Whether synge eos  is added
!  logical, public, protected              :: srhd_eos = .false.

!  !> Whether viscosity is added
!  logical, public, protected              :: srhd_viscosity = .false.

!  !> Whether gravity is added
!  logical, public, protected              :: srhd_gravity = .false.

!  !> Whether Hall-MHD is used
!  logical, public, protected              :: srhd_Hall = .false.

!  !> Whether particles module is added
!  logical, public, protected              :: srhd_particles = .false.

!  !> Whether magnetofriction is added
!  logical, public, protected              :: srhd_magnetofriction = .false.

!  !> Whether GLM-MHD is used
!  logical, public, protected              :: srhd_glm = .false.

!  !> Whether divB cleaning sources are added splitting from fluid solver
!  logical, public, protected              :: source_split_divb = .false.

!  !> GLM-MHD parameter: ratio of the diffusive and advective time scales for div b
!  !> taking values within [0, 1]
!  double precision, public                :: srhd_glm_alpha = 0.5d0

!  !> MHD fourth order
!  logical, public, protected              :: srhd_4th_order = .false.

!  !> Number of tracer species
!  integer, public, protected              :: srhd_n_tracer = 0

!  !> Index of the density (in the w array)
!  integer, public, protected              :: rho_
!
!  !> Indices of the momentum density
!  integer, allocatable, public, protected :: mom(:)
!
!  !> Index of the energy density (-1 if not present)
!  integer, public, protected              :: e_

!  !> Index of the gas pressure (-1 if not present) should equal e_
!  integer, public, protected              :: p_

!  !> Indices of the magnetic field
!  integer, allocatable, public, protected :: mag(:)

!  !> Indices of the GLM psi
!  integer, public, protected :: psi_

!  !> Indices of the tracers
!  integer, allocatable, public, protected :: tracer(:)
!
!  !> The adiabatic index
!  double precision, public                :: srhd_gamma = 5.d0/3.0d0
!
!  !> The adiabatic constant
!  double precision, public                :: srhd_adiab = 1.0d0
!
!  !> The MHD resistivity
!  double precision, public                :: srhd_eta = 0.0d0

!  !> The MHD hyper-resistivity
!  double precision, public                :: srhd_eta_hyper = 0.0d0

!  !> TODO: what is this?
!  double precision, public                :: srhd_etah = 0.0d0


!  logical,          public, protected     :: srhd_checkNR = .true.
!  double precision, public, protected     :: srhd_absaccnr=1.0d-8
!  double precision, public, protected     :: srhd_tolerNR =1.0d-9

!  !> The small_est allowed energy
!  double precision, protected             :: small_e


!  !> The small_est allowed inertia
!  double precision, protected             :: small_xi

!  ! smaller values for speed
!  double precision, public, protected             :: small_vec2  = 0.0
!  !> The number of waves
!  integer :: nwwave=8

!  !> Method type to clean divergence of B
!  character(len=std_len), public, protected :: typedivbfix  = 'linde'
!
!  !> Method type in a integer for good performance
!  integer :: type_divb
!
!  !> Coefficient of diffusive divB cleaning
!  double precision :: divbdiff     = 0.8d0
!
!  !> Update all equations due to divB cleaning
!  character(len=std_len) ::    typedivbdiff = 'all'

!  !> Use a compact way to add resistivity
!  logical :: compactres   = .false.

!  !> Add divB wave in Roe solver
!  logical, public :: divbwave     = .true.

!  !> Helium abundance over Hydrogen
!  double precision, public, protected  :: He_abundance=0.1d0

!  !> To control divB=0 fix for boundary
!  logical, public, protected :: boundary_divbfix(2*^ND)=.true.

!  !> To skip * layer of ghost cells during divB=0 fix for boundary
!  integer, public, protected :: boundary_divbfix_skip(2*^ND)=0

!  !> B0 field is force-free
!  logical, public, protected :: B0field_forcefree=.true.


!  ! DivB cleaning methods
!  integer, parameter :: divb_none          = 0
!  integer, parameter :: divb_glm1          = 1
!  integer, parameter :: divb_glm2          = 2
!  integer, parameter :: divb_powel         = 3
!  integer, parameter :: divb_janhunen      = 4
!  integer, parameter :: divb_linde         = 5
!  integer, parameter :: divb_lindejanhunen = 6
!  integer, parameter :: divb_lindepowel    = 7
!  integer, parameter :: divb_lindeglm      = 8

  ! Public methods
  public :: srhd_phys_init
  public :: srhd_kin_en_primitive
  public :: srhd_get_pthermal
  public :: srhd_get_p_mag
  public :: srhd_get_p_total
  public :: srhd_get_v
  public :: srhd_get_v_idim
  public :: srhd_to_conserved
  public :: srhd_to_primitive
  public :: srhd_get_csound2
  public :: srhd_get_csound_prim
  public :: get_divb
  public :: get_current
  public :: get_normalized_divb
  public :: srhd_get_4u_from_3v
contains

  !> Read this module"s parameters from a file
  subroutine srhd_read_params(files)
    ! made by Z. Meliani 20/02/2018
    use mod_global_parameters
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /srhd_list/ srhd_energy, srhd_eos,srhd_n_tracer, srhd_gamma, &
                          srhd_adiab, srhd_eta, srhd_eta_hyper, &
                          srhd_etah, srhd_glm_alpha, srhd_magnetofriction,&
                          srhd_thermal_conduction, srhd_radiative_cooling, &
                          srhd_Hall, srhd_gravity, srhd_viscosity, &
                          srhd_4th_order, typedivbfix, source_split_divb, &
                          divbdiff, typedivbdiff, compactres, divbwave, srhd_glm,&
                          He_abundance, SI_unit, B0field,&
                          B0field_forcefree, Bdip, Bquad, Boct, Busr, &
                          srhd_particles, boundary_divbfix, &
                          boundary_divbfix_skip, &
                          srhd_maxiterationNR, srhd_absaccNR,srhd_tolerNr,&
                          srhd_checkNR,srhd_maxdspeed,small_vec2
    write(*,*)'Reading srhd_list'
    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, srhd_list, end=111)
111    close(unitpar)
    end do

  end subroutine srhd_read_params

  !> Write this module's parameters to a snapsoht
  subroutine srhd_write_info(fh)
    ! made by Z. Meliani 10/02/2018
    use mod_global_parameters
    integer, intent(in)                 :: fh
    integer, parameter                  :: n_par = 1
    double precision                    :: values(n_par)
    character(len=name_len)             :: names(n_par)
    integer, dimension(MPI_STATUS_SIZE) :: st
    integer                             :: er

    call MPI_FILE_WRITE(fh, n_par, 1, MPI_INTEGER, st, er)

    names(1) = "gamma"
    values(1) = srhd_gamma
    call MPI_FILE_WRITE(fh, values, n_par, MPI_DOUBLE_PRECISION, st, er)
    call MPI_FILE_WRITE(fh, names, n_par * name_len, MPI_CHARACTER, st, er)
  end subroutine srhd_write_info

  subroutine srhd_angmomfix(fC,x,wnew,ixI^L,ixO^L,idim)
    use mod_global_parameters
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision, intent(inout)    :: fC(ixI^S,1:nwflux,1:ndim),  wnew(ixI^S,1:nw)
    integer, intent(in)                :: ixI^L, ixO^L
    integer, intent(in)                :: idim
    integer                            :: hxO^L, kxC^L, iw
    double precision                   :: inv_volume(ixI^S)

    call mpistop("to do")
    ! ! shifted indexes
    ! hxO^L=ixO^L-kr(idim,^D);
    ! ! all the indexes
    ! kxCmin^D=hxOmin^D;
    ! kxCmax^D=ixOmax^D;
    !
    ! inv_volume = 1.0d0/block%dvolume(ixO^S)
    !
    ! select case(typeaxial)
    ! case ("cylindrical")
    !   do iw=1,nwflux
    !     if (idim==r_ .and. iw==iw_mom(phi_)) then
    !       fC(kxC^S,iw,idim)= fC(kxC^S,iw,idim)*(x(kxC^S,r_)+half*block%dx(kxC^S,r_))
    !       wnew(ixO^S,iw)=wnew(ixO^S,iw) + (fC(ixO^S,iw,idim)-fC(hxO^S,iw,idim)) * &
    !            (inv_volume/x(ixO^S,r_))
    !     else
    !       wnew(ixO^S,iw)=wnew(ixO^S,iw) + (fC(ixO^S,iw,idim)-fC(hxO^S,iw,idim)) * &
    !             inv_volume
    !     endif
    !   enddo
    ! case ("spherical")
    !   do iw=1,nwflux
    !     if     (idim==r_ .and. (iw==iw_mom(2) .or. iw==iw_mom(phi_))) then
    !       fC(kxC^S,iw,idim)= fC(kxC^S,iw,idim)*(x(kxC^S,r_)+half*block%dx(kxC^S,r_))
    !       wnew(ixO^S,iw)=wnew(ixO^S,iw) + (fC(ixO^S,iw,idim)-fC(hxO^S,iw,idim)) * &
    !            (inv_volume/x(ixO^S,r_))
    !     elseif (idim==2  .and. iw==iw_mom(phi_)) then
    !       fC(kxC^S,iw,idim)=fC(kxC^S,iw,idim)*dsin(x(kxC^S,2)+half*block%dx(kxC^S,2)) ! (x(4,3,1)-x(3,3,1)))
    !       wnew(ixO^S,iw)=wnew(ixO^S,iw) + (fC(ixO^S,iw,idim)-fC(hxO^S,iw,idim)) * &
    !            (inv_volume/dsin(x(ixO^S,2)))
    !     else
    !       wnew(ixO^S,iw)=wnew(ixO^S,iw) + (fC(ixO^S,iw,idim)-fC(hxO^S,iw,idim)) * &
    !             inv_volume
    !     endif
    !   enddo
    !
    !   ! if (idim==r_) then
    !   !   fC(kxC^S,iw_mom(phi_),idim)= fC(kxC^S,iw_mom(phi_),idim)*(x(kxC^S,r_)+half*block%dx(kxC^S,r_))
    !   !   fC(kxC^S,iw_mom(phi_),idim)= fC(kxC^S,iw_mom(phi_),idim)*(x(kxC^S,r_)+half*block%dx(kxC^S,r_))
    !   !   wnew(ixO^S,iw_mom(phi_))=wnew(ixO^S,iw_mom(phi_)) + (fC(ixO^S,iw_mom(phi_),idim)-fC(hxO^S,iw_mom(phi_),idim)) * &
    !   !        (inv_volume/x(ixO^S,r_))
    !   !
    !   ! elseif (idim==2) then
    !   !   fC(hxOmin1:hxOmax1,hxOmin2,hxOmin3:hxOmax3,iw,idim)=fC(hxOmin1:hxOmax1,hxOmin2,hxOmin3:hxOmax3,iw,idim)*dsin(x(hxOmin1:hxOmax1,hxOmin2,hxOmin3:hxOmax3,2)+half*block%dx(hxOmin1:hxOmax1,hxOmin2,hxOmin3:hxOmax3,2)) ! (x(4,3,1)-x(3,3,1)))
    !   !   fC(ixO^S,iw,idim)=fC(ixO^S,iw,idim)*dsin(x(ixO^S,2)+half*block%dx(ixO^S,2)) ! (x(4,3,1)-x(3,3,1)))
    !   !   wnew(ixO^S,iw)=wnew(ixO^S,iw) + (fC(ixO^S,iw,idim)-fC(hxO^S,iw,idim)) * &
    !   !        (inv_volume/dsin(x(ixO^S,2)))
    !   ! endif
    !
    ! end select

  end subroutine srhd_angmomfix

  subroutine srhd_phys_init()
    ! made by Z. Meliani 10/02/2018
    use mod_global_parameters
    use mod_thermal_conduction
    use mod_radiative_cooling
    use mod_viscosity, only: viscosity_init
    use mod_gravity, only: gravity_init
    use mod_particles, only: particles_init
    use mod_magnetofriction, only: magnetofriction_init
    use mod_physics

    integer :: itr, idir
    unit_velocity=const_c
    call srhd_read_params(par_files)

    physics_type = "srmhd"
    phys_energy=srhd_energy
    ! set default gamma for polytropic/isothermal process
    if(.not.srhd_energy) srhd_gamma=1.d0
    use_particles=srhd_particles


    ! Determine flux variables
    rho_ = var_set_rho()
    d_=rho_
    allocate(mom(ndir))
    mom(:) = var_set_momentum(ndir)

    ! Set index of energy variable
    if (srhd_energy) then
      nwwave = 8
      e_     = var_set_energy() ! energy density
      p_     = e_               ! gas pressure
    else
      nwwave = 7
      e_     = -1
      p_     = -1
    end if



    allocate(tracer(srhd_n_tracer))

    ! Set starting index of tracers
    do itr = 1, srhd_n_tracer
      tracer(itr) = var_set_fluxvar("trc", "trp", itr, need_bc=.false.)
    end do

    ! Set index for auxiliary variables
    pp_  = var_set_auxvar('pp')
    lfac_=var_set_auxvar('lfac_')
    !nwfluxbc=nwfluxbc+2

    nvector      = 1 ! No. vector vars
    allocate(iw_vector(nvector))
    iw_vector(1) = mom(1) - 1   ! TODO: why like this?


    ! Check whether custom flux types have been defined
    if (.not. allocated(flux_type)) then
       allocate(flux_type(ndir, nw))
       flux_type = flux_default
    else if (any(shape(flux_type) /= [ndir, nw])) then
       call mpistop("phys_check error: flux_type has wrong shape")
    end if
    do idir=1,ndir
      if(ndim>1) flux_type(idir,mag(idir))=flux_tvdlf
    end do


    srhd_maxspeed=1.0-srhd_maxdspeed

    allocate(srhd_iw_average(1:nw))
    srhd_iw_average = .true.

    phys_get_dt              => srhd_get_dt
    phys_get_cmax            => srhd_get_cmax
    phys_get_cbounds         => srhd_get_cbounds
    phys_get_flux            => srhd_get_flux
    phys_get_v_idim          => srhd_get_v_idim
    phys_add_source_geom     => srhd_add_source_geom
    phys_add_source          => srhd_add_source
    phys_to_conserved        => srhd_to_conserved
    phys_to_primitive        => srhd_to_primitive
    phys_get_aux             => srhd_get_auxiliary
    phys_check_params        => srhd_check_params
    phys_check_w             => srhd_check_w
    phys_get_pthermal        => srhd_get_pthermal
    phys_boundary_adjust     => srhd_boundary_adjust
    phys_write_info          => srhd_write_info
    phys_angmomfix           => srhd_angmomfix
    phys_handle_small_values => srhd_handle_small_values


    phys_iw_average =>srhd_iw_average
    ! Whether diagonal ghost cells are required for the physics
    if(type_divb < divb_linde) phys_req_diagonal = .false.

    ! derive units from basic units
    call srhd_physical_units()

    if(.not. srhd_energy .and. srhd_thermal_conduction) then
      call mpistop("thermal conduction needs srhd_energy=T")
    end if
    if(.not. srhd_energy .and. srhd_radiative_cooling) then
      call mpistop("radiative cooling needs srhd_energy=T")
    end if

    ! initialize thermal conduction module
    !if (srhd_thermal_conduction) then
    !  phys_req_diagonal = .true.
    !  call thermal_conduction_init(srhd_gamma)
    !end if

    ! Initialize radiative cooling module
    !if (srhd_radiative_cooling) then
    !  call radiative_cooling_init(srhd_gamma,He_abundance)
    !end if

    ! Initialize viscosity module
    !if (srhd_viscosity) call viscosity_init(phys_wider_stencil,phys_req_diagonal)

    ! Initialize gravity module
    if(srhd_gravity) then
      call gravity_init()
    end if

    ! Initialize particles module
    if(srhd_particles) then
      call particles_init()
      phys_req_diagonal = .true.
    end if

    srhd_config%dust_on        = .false.
    srhd_config%dust_n_species = 0
    srhd_config%energy         = srhd_energy
    srhd_config%gamma          = srhd_gamma
    phys_iw_average =>srhd_iw_average
    phys_config     => srhd_config
  end subroutine srhd_phys_init

  subroutine srhd_check_params
    use mod_global_parameters

    ! after user parameter setting
    gamma_1=srhd_gamma-1.0_dp
    inv_gamma_1=1.d0/gamma_1
    gamma_to_gamma_1=srhd_gamma/gamma_1
    if (.not. srhd_energy) then
       if (srhd_gamma <= 0.0d0) call mpistop ("Error: srhd_gamma <= 0")
       if (srhd_adiab < 0.0d0) call mpistop ("Error: srhd_adiab < 0")
       small_pressure = srhd_adiab*small_density**srhd_gamma
    else
       if (srhd_gamma <= 0.0d0 .or. srhd_gamma == 1.0d0) &
            call mpistop ("Error: srhd_gamma <= 0 or srhd_gamma == 1")
       small_e = small_pressure * inv_gamma_1
    end if
    small_xi=small_density+gamma_to_gamma_1*small_pressure

  end subroutine srhd_check_params

  subroutine srhd_physical_units()
    use mod_global_parameters
    double precision :: mp,kB,miu0
    ! Derive scaling units
    if(SI_unit) then
      mp=mp_SI
      kB=kB_SI
      miu0=miu0_SI
    else
      mp=mp_cgs
      kB=kB_cgs
      miu0=4.d0*dpi
    end if
    if(unit_velocity==0) then
      call mpistop ("Error: in srmhd the unit_velocity=c")
    else
      unit_density=(1.d0+4.d0*He_abundance)*mp*unit_numberdensity
      unit_pressure=unit_density*unit_velocity**2
      unit_temperature=unit_pressure/((2.d0+3.d0*He_abundance)*unit_numberdensity*kB)
      unit_time=unit_length/unit_velocity
    end if

  end subroutine srhd_physical_units

  subroutine srhd_check_w(primitive,ixI^L,ixO^L,w,flag)
    ! made by Z. Meliani 10/02/2018
    use mod_global_parameters

    logical, intent(in)            :: primitive
    integer, intent(in)            :: ixI^L, ixO^L
    double precision, intent(in)   :: w(ixI^S,nw)
    integer, intent(inout)         :: flag(ixI^S)
    double precision :: tmp(ixI^S)

    flag(ixO^S)=0
    where(w(ixO^S, rho_)/w(ixO^S,lfac_) < small_density) flag(ixO^S) = rho_

    cond_energy : if (srhd_energy) then
       cond_einternal : if (block%e_is_internal) then
          where(w(ixO^S, e_) < small_pressure*inv_gamma_1) flag(ixO^S) = e_
       else cond_einternal
         is_prim : if (primitive)then
           where(w(ixO^S, p_) < small_pressure) flag(ixO^S) = p_ ! p_=e_
         else is_prim
           ! Calculate pressure=(gamma-1)*(e-0.5*(2ek+2eb))
           ! in srmhd we check if xi-p-d is positive
             ! do be done
             where(w(ixO^S, e_) < small_e) flag(ixO^S) = e_
         end if is_prim
       end if cond_einternal
    end if cond_energy
  end subroutine srhd_check_w

  !> Transform primitive variables into conservative ones
  subroutine srhd_to_conserved(ixI^L,ixO^L,w,x)
    ! made by Z. Meliani 10/02/2018
    use mod_global_parameters
    integer, intent(in)                :: ixI^L, ixO^L
    double precision, intent(inout)    :: w(ixI^S, nw)
    double precision, intent(in)       :: x(ixI^S, 1:ndim)
    integer                            :: idir, itr
    double precision, dimension(ixO^S) :: sqrU,rhoh,xi
    integer, dimension(ixO^S)          :: flag_error
    character(len=30)                  :: subname_loc


    subname_loc= 'srhd_to_conserved'
    flag_error(ixO^S)=0
    where(w(ixO^S,rho_)<small_density)flag_error(ixO^S)=rho_


    sqrU(ixO^S)    = sum(w(ixO^S, mom(:))**2, dim=ndim+1)
    w(ixO^S,lfac_) = dsqrt(1.0_dp+sqrU(ixO^S))
    sqrV=sqrU/w(ixO^S,lfac_)
    ! fill the auxiliary variable xi and density D
    call srhd_get_enthalpy(ixO^L,w(ixO^S,rho_),w(ixO^S,p_),rhoh)
    ! with enthalpy w: xi= lfac^2 rhoh
    xi(ixO^S) = w(ixO^S,lfac_)**2.0D0*rhoh(ixO^S)
    ! density: d = lfac * rho
    w(ixO^S,rho_)=w(ixO^S,lfac_)*w(ixO^S,rho_)


    ! Convert velocity to momentum
    ! s= (xi + B^2) * v - (v.B) * B
    ! re-use rhoh as rhoh=lfac*rhoh
    rhoh(ixO^S)=  w(ixO^S,lfac_)* rhoh(ixO^S)
    Loop_idirmom : do idir = 1, ndir
       w(ixO^S, mom(idir)) = rhoh(ixO^S)* w(ixO^S, mom(idir))
    end do Loop_idirmom

    cond_energy : if (srhd_energy) then
       where(sqrV(ixO^S)<0.0_dp) flag_error(ixO^S)=e_
       ! E = xi - p - D
       w(ixO^S,e_)=xi(ixO^S) - w(ixO^S,p_) - w(ixO^S,d_)
    end if cond_energy
    if (check_small_values) call srhd_handle_small_values(.false., &
                                  w, x, ixI^L, ixO^L,trim(subname_loc),&
                                  flag_error=flag_error)
  end subroutine srhd_to_conserved

  !> Transform conservative variables into primitive ones
  subroutine srhd_to_primitive(ixI^L,ixO^L,w,x)
    ! made by Z. Meliani 10/02/2018
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    double precision                :: inv_rho(ixI^S)
    integer                         :: itr, idir
    character(len=30)               :: subname_loc

    subname_loc='srhd_to_primitive'

    ! get auxiliary variables
    call srhd_get_auxiliary(.true.,w,x,ixI^L,ixO^L,subname_loc)
    w(ixO^S,rho_) = w(ixO^S,d_)/w(ixO^S,lfac_)
    if (srhd_energy) then

    ! re-use inv_rho to store inverse of the density
    inv_rho = 1.0_dp / w(ixO^S, rho_)


    Loop_idir : do idir=1,ndir
     w(ixO^S,mom(idir)) = w(ixO^S,lfac_)*w(ixO^S,mom(idir))
    end do Loop_idir

    if (check_small_values) call srhd_handle_small_values(.true., w, x&
                                 , ixI^L, ixO^L,trim(subname_loc))
  end subroutine srhd_to_primitive

  subroutine srhd_handle_small_values(primitive, w, x, ixI^L, ixO^L, &
                                       subname,flag_error)
    ! made by Z. Meliani 10/02/2018
    use mod_global_parameters
    use mod_small_values
    logical, intent(in)             :: primitive
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    character(len=*), intent(in)    :: subname
    integer, optional, intent(in)   :: flag_error(ixO^S)

    ! .. local ..
    double precision :: smallone
    integer :: idir, ierror,flag(ixI^S)
    !---------------------------------------------------------
    if (small_values_method == "ignore") return
    if(present(flag_error)) then
     flag(ixO^S) = flag_error(ixO^S)
    else
     call srhd_check_w(primitive, ixI^L, ixO^L, w, flag)
    end if
    if (any(flag(ixO^S) /= 0)) then
       select case (small_values_method)
       case ("replace")
          where(flag(ixO^S) /= 0) w(ixO^S,rho_) = small_density

          do idir = 1, ndir
             where(flag(ixO^S) /= 0) w(ixO^S, mom(idir)) = 0.0d0
          end do

          if (srhd_energy) then
             if(primitive) then
               smallone = small_pressure
             else
               smallone = small_e
             end if
             where(flag(ixO^S) /= 0) w(ixO^S,e_) = smallone
          end if
       case ("average")
          call small_values_average(ixI^L, ixO^L, w, x, flag,ierror)
       case default
          call small_values_error(w, x, ixI^L, ixO^L, flag, subname)
       end select
    end if
  end subroutine srhd_handle_small_values




 !> Calculate thermal pressure for enthalpy and density
  subroutine srhd_get_pressure_fromprimitive(ixI^L,ixO^L,w,pth)
    ! made by Z. Meliani 10/02/2018
    use mod_global_parameters
    implicit none
    integer, intent(in)            :: ixI^L, ixO^L
    double precision, intent(in)   :: w(ixI^S,nw)
    double precision, intent(out)  :: pth(ixI^S)

    double precision               :: rhoh(ixO^S)

    ptherm(ixO^S)=w(ixO^S,pp_)

  end subroutine srhd_get_pressure_fromprimitive


  !> Calculate thermal pressure=(gamma-1)*(e-0.5*m**2/rho-b**2/2) within ixO^L
  subroutine srhd_get_pthermal(w,x,ixI^L,ixO^L,pth)
    ! made by Z. Meliani 10/02/2018
    use mod_global_parameters
    implicit none

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: w(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: pth(ixI^S)

    double precision                :: rho(ixO^S),rhoh(ixO^S)


    call srhd_get_auxiliary(.true.,w,x,ixI^L,ixO^L,'srhd_get_pthermal')
    ptherm(ixO^S)=w(ixO^S,p_)


  end subroutine srhd_get_pthermal


 !> Calculate the square of the thermal sound speed csound2 within ixO^L.
  !> csound2=gamm*p/rho
  subroutine srhd_get_csound2_prim(w,x,ixI^L,ixO^L,csound2)
    ! made by Z. Meliani 13/02/2018
    use mod_global_parameters
    implicit none
    integer, intent(in)               :: ixI^L, ixO^L
    double precision, intent(in)      :: w(ixI^S,nw)
    double precision, intent(in)      :: x(ixI^S,1:ndim)
    double precision, intent(out)     :: csound2(ixO^S)

    double precision,dimension(ixO^S) :: rhoh


    rhoh=w(ixO^S,xi_)/w(ixO^S,lfac_)**2.0d0

    call srhd_get_csound2_prim_eos(ixI^L,ixO^L,x,w(ixO^S,rho_),&
                                    rhoh,w(ixO^S,p_),csound2)
  end subroutine srhd_get_csound2_prim

  !> Convert energy to entropy
  subroutine e_to_rhos(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision,intent(inout)  :: w(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)

    if (srhd_energy) then
      ! Calculate pressure and Lorentz factor from conservative variables only
      ! these are put in lfac_ and p_ auxiliaries
      call getaux(.true.,w,ixI^L,ixO^L,patchtrue,'e_to_rhos')
      w(ixO^S,rhos_)=w(ixO^S,d_)*w(ixO^S,p_)*(w(ixO^S,lfac_)/w(ixO^S,d_))**eqpar(gamma_)
    else
      call mpistop("e_to_rhos can not be used without energy equation!")
    end if
  end subroutine e_to_rhos

  !> Convert entropy to energy
  subroutine rhos_to_e(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in) :: ixI^L, ixO^L
    double precision :: w(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)

    if (srhd_energy) then
      ! Calculate pressure and Lorentz factor from conservative variables BUT
      ! with tau replaced by DS, these are put in lfac_ and p_ auxiliaries
      call getaux2(.true.,w,ixI^L,ixO^L,'rhos_to_e')

      xi(ixO^S)=w(ixO^S,lfac_)*w(ixO^S,d_)+w(ixO^S,lfac_)**2*w(ixO^S,p_)*eqpar(gamma_)/(eqpar(gamma_)-one)
      w(ixO^S,e_)=xi(ixO^S)-w(ixO^S,p_)-w(ixO^S,d_)

       if(.not.block%e_is_internal) &
         w(ixO^S, e_) =w(ixO^S, e_) + srhd_kin_en_primitive(w, ixI^L, ixO^L) + &
            srhd_mag_en_primitive(w, ixI^L, ixO^L)
    else
       call mpistop("rhos_to_e can not be used without energy equation!")
    end if
  end subroutine rhos_to_e

  !> Calculate v vector
  subroutine srhd_get_v(w,x,ixI^L,ixO^L,v,xi)
    ! made by Z. Meliani 13/02/2018
    use mod_global_parameters

    integer, intent(in)                  :: ixI^L, ixO^L
    double precision, intent(in)         :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(out)        :: v(ixI^S,1:ndir)
    double precision, intent(in), option :: xi(ixI^S)
    ! .. local ..
    double precision                     :: xi_loc(ixI^S)
    integer :: idir
    !----------------------------------------------------------
     if(present(xi)) then
      Loop_idir_v2_a: do idir=1,ndir
       v(ixO^S,idir) = w(ixO^S, mom(idir))/xi(ixO^S)
     end do Loop_idir_v2_a
     else
      xi_loc(ixO^S) = w(ixO^S,e_)+w(ixO^S,d_)+w(ixO^S,p_)
      Loop_idir_v2: do idir=1,ndir
       v(ixO^S,idir) = w(ixO^S, mom(idir))/xi_loc(ixO^S)
      end do Loop_idir_v2
    end if
  end subroutine srhd_get_v

  !> Calculate v vector
  subroutine srhd_get_v_prim(w,x,ixI^L,ixO^L,v)
    ! made by Z. Meliani 13/02/2018
    use mod_global_parameters

    integer, intent(in)                  :: ixI^L, ixO^L
    double precision, intent(in)         :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(out)        :: v(ixI^S,1:ndir)
    ! .. local ..
    integer :: idir
    !----------------------------------------------------------
      Loop_idir_v2_a: do idir=1,ndir
       v(ixO^S,idir) = w(ixO^S, mom(idir))/w(ixO^S,lfac_)
     end do Loop_idir_v2_a
  end subroutine srhd_get_v_prim

  !> Calculate v component
  subroutine srhd_get_v_idim(w,x,ixI^L,ixO^L,idim,v_idim)
    ! made by Z. Meliani 10/02/2018
    use mod_global_parameters

    integer, intent(in)                    :: ixI^L, ixO^L, idim
    double precision, intent(in)           :: w(ixI^S,nw), x(ixI^S,1:ndim)
    double precision, intent(out)          :: v_idim(ixI^S)
    ! .. local ..
    double precision                       :: xi(ixO^S)

     xi(ixO^S)     = w(ixO^S,e_)+w(ixO^S,d_)+w(ixO^S,p_)
     v_idim(ixO^S) = w(ixO^S, mom(idim)) /xi(ixO^S,xi_)
  end subroutine srhd_get_v_idim

  !> Calculate v component
  subroutine srhd_get_v_idim_loc(w,x,ixI^L,ixO^L,idim,v_idim)
    ! made by Z. Meliani 13/02/2018
    use mod_global_parameters

    integer, intent(in)                    :: ixI^L, ixO^L, idim
    double precision, intent(in)           :: w(ixI^S,nw), x(ixI^S,1:ndim)
    double precision, intent(out)          :: v_idim(ixI^S)
    ! .. local ..
    double precision                       :: xi(ixO^S)


    xi(ixO^S)     = w(ixO^S,e_)+w(ixO^S,d_)+w(ixO^S,p_)
    v_idim(ixO^S) = w(ixO^S, mom(idim)) /xi(ixO^S,xi_)

  end subroutine srhd_get_v_idim_loc

 !> Calculate v^2
  subroutine srhd_get_v2(w,x,ixI^L,ixO^L,v2)
    ! made by Z. Meliani 13/02/2018
    use mod_global_parameters

    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S,nw), x(ixI^S,1:ndim)
    double precision, intent(out) :: v2(ixO^S)
    ! .. local ..
    double precision                        :: xi(ixO^S)


     xi(ixO^S)     = w(ixO^S,e_)+w(ixO^S,d_)+w(ixO^S,p_)
     v2(ixO^S) = sum(w(ixO^S, mom(1):mom(ndir))**2.0 ,dim=ndim+1)&
                      /xi(ixO^S,xi_)**2.0

  end subroutine srhd_get_v2

  subroutine srhd_get_cspeed(ixI^L,ixO^L,idim,x,w,rhoh,vidim,csound2,cmax&
                                   ,from_cbound,cmin)
    ! made by Z. Meliani 10/02/2018
    use mod_global_parameters
    implicit none
    integer, intent(in)          :: ixI^L, ixO^L, idim
    double precision, intent(in) :: w(ixI^S, nw), x(ixI^S,1:ndim)
    double precision, intent(in), dimension(ixO^S) :: rhoh, csound2
    double precision, intent(inout)          :: cmax(ixI^S)
    double precision, optional, intent(inout):: cmin(ixI^S)
    logical, optional, intent(in)            :: from_cbound



    double precision, dimension(ixO^S):: vidim,v2


    if(.not.present(from_cbound).or..not.present(cmin))&
                                        vidim(ixO^S)=dabs(vidim(ixO^S))

    cond_onedir: if(ndir==1)then
     cmax(ixO^S)=(vidim(ixO^S)+dsqrt(csound2(ixO^S)))&
                 /(1.0d0+dsqrt(csound2(ixO^S)*vidim(ixO^S)**2.0d0))
     if(present(cmin))cmin(ixO^S)=(vidim(ixO^S)-dsqrt(csound2(ixO^S)))&
               /(1.0d0+dsqrt(csound2(ixO^S)*vidim(ixO^S)**2.0d0))
    else cond_onedir

     call srhd_get_v2(w,x,ixI^L,ixO^L,v2)
     cmax(ixO^S)=(vidim(ixO^S)*(1.0d0-csound2(ixO^S))&
                  +dsqrt(csound2(ixO^S)*(1.0d0-v2(ixO^S))*((1.0d0-v2(ixO^S)*csound2(ixO^S))&
                         -vidim(ixO^S)**2.0d0*(1.0d0-csound2(ixO^S)))))&
                 /( 1.0d0-v2(ixO^S)*csound2(ixO^S))
     if(present(cmin))cmin(ixO^S)=(vidim(ixO^S)*(1.0d0-csound2(ixO^S))&
                  -dsqrt(csound2(ixO^S)*(1.0d0-v2(ixO^S))*(1.0d0-v2(ixO^S)*csound2(ixO^S)&
                         -vidim(ixO^S)**2.0d0*(1.0d0-csound2(ixO^S)))))&
                 /( 1.0d0-v2(ixO^S)*csound2(ixO^S))
    end if cond_onedir



 end subroutine srhd_get_cspeed


  !> Calculate cmax_idim using gammie method within ixO^L
  subroutine srhd_get_cmax(w,x,ixI^L,ixO^L,idim,cmax)
    ! made by Z. Meliani 13/02/2018
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L, idim
    double precision, intent(in) :: w(ixI^S, nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: cmax(ixI^S)

    double precision, dimension(ixO^S):: xi,vidim
    double precision, dimension(ixO^S):: v2,csound2

    xi(ixO^S)=w(ixO^S,e_)+w(ixO^S,d_)+w(ixO^S,p_)
    call srhd_get_csound2(w,x,ixI^L,ixO^L,rhoh,csound2)
    call srhd_get_v_idim_loc(w,x,ixI^L,ixO^L,idim,vidim)
    call srhd_get_cspeed(ixI^L,ixO^L,idim,x,w,rhoh=rhoh,vidim=vidim,cmax=cmax)

  end subroutine srhd_get_cmax

  !> Estimating bounds for the minimum and maximum signal velocities
  subroutine srhd_get_cbounds(wLC,wRC,wLp,wRp,x,ixI^L,ixO^L,idim,cmax,cmin)
    ! made by Z. Meliani 13/02/2018
    use mod_global_parameters

    integer, intent(in)                       :: ixI^L, ixO^L, idim
    double precision, intent(in)              :: wLC(ixI^S, nw), wRC(ixI^S, nw)
    double precision, intent(in)              :: wLp(ixI^S, nw), wRp(ixI^S, nw)
    double precision, intent(in)              :: x(ixI^S,1:ndim)
    double precision, intent(inout)           :: cmax(ixI^S)
    double precision, intent(inout), optional :: cmin(ixI^S)

    double precision                   :: wmean(ixI^S,nw)
    double precision, dimension(ixO^S) :: umean, dmean, csound2Lp, &
                                          csound2Rp, tmp1,tmp2,tmp3, &
                                          rhohLp,rhohRp, &
                                          cmaxL,cmaxR,csound2,rhoh,&
                                          vidim,vidimLp,vidimRp
    character(len=30)                  :: subname_loc

    subname_loc='srhd_get_cbounds'
    if (typeboundspeed/='cmaxmean') then
      ! This implements formula (10.52) from "Riemann Solvers and Numerical
      ! Methods for Fluid Dynamics" by Toro.
      rhohLp=wLp(ixO^S,xi_)/wLp(ixO^S,lfac_)**2.0d0
      rhohRp=wRp(ixO^S,xi_)/wRp(ixO^S,lfac_)**2.0d0

      call srhd_get_csound2(wLp,x,ixI^L,ixO^L,rhohLp,csound2Lp)
      call srhd_get_csound2(wRp,x,ixI^L,ixO^L,rhohRp,csound2Rp)




      call srhd_get_v_idim_loc(wLp,x,ixI^L,ixO^L,idim,vidimLp)
      call srhd_get_v_idim_loc(wRp,x,ixI^L,ixO^L,idim,vidimRp)

      tmp1(ixO^S)=sqrt(xiLp(ixO^S))
      tmp2(ixO^S)=sqrt(xiRp(ixO^S))
      tmp3(ixO^S)=1.d0/(sqrt(xiLp(ixO^S))&
                  +sqrt(wRp(ixO^S,rho_)))

      umean(ixO^S)=(wLp(ixO^S,mom(idim))*tmp1(ixO^S)&
                    +wRp(ixO^S,mom(idim))*tmp2(ixO^S))*tmp3(ixO^S)

     call srhd_get_cspeed(ixI^L,ixO^L,idim,x,wLp,rhoh=rhohLp,&
                               vidim=vidimLp,cmax=cmaxL)

     call srhd_get_cspeed(ixI^L,ixO^L,idim,x,wRp,rhoh=rhohRp,&
                               vidim=vidimRp,cmax=cmaxR)


      dmean(ixO^S)=(tmp1(ixO^S)*cmaxL(ixO^S)**2&
                    +tmp2(ixO^S)*cmaxR(ixO^S)**2)*tmp3(ixO^S)+&
                    0.5d0*tmp1(ixO^S)*tmp2(ixO^S)*tmp3(ixO^S)**2*&
                    (wRp(ixO^S,mom(idim))-wLp(ixO^S,mom(idim)))**2
      dmean(ixO^S)=sqrt(dmean(ixO^S))
      if(present(cmin)) then
        cmin(ixO^S)=umean(ixO^S)-dmean(ixO^S)
        cmax(ixO^S)=umean(ixO^S)+dmean(ixO^S)
      else
        cmax(ixO^S)=abs(umean(ixO^S))+dmean(ixO^S)
      end if
    else
      wmean(ixO^S,1:nwflux)=0.5d0*(wLC(ixO^S,1:nwflux)+wRC(ixO^S,1:nwflux))
      ! get auxiliary variables
      call srhd_get_auxiliary(.true.,wmean,x,ixI^L,ixO^L,subname_loc)
      rhoh=wmean(ixO^S,xi_)/wmean(ixO^S,lfac_)**2.0d0

      call srhd_get_csound2(wmean,x,ixI^L,ixO^L,rhoh,csound2)

      call srhd_get_v_idim_loc(wmean,x,ixI^L,ixO^L,idim,vidim)
      call srhd_get_cspeed(ixI^L,ixO^L,idim,x,wmean,rhoh=rhoh,vidim=vidim,cmax=cmax&
                                 ,from_cbound=.true.,cmin=cmin)
      if(present(cmin)) then
        cmax(ixO^S)=min(max(cmax(ixO^S),0.0D0),1.0d0)
        cmin(ixO^S)=max(min(cmin(ixO^S),0.0D0),-1.0d0)
      else
        cmax(ixO^S)=min(cmax(ixO^S),1.0d0)
      end if
    end if
  end subroutine srhd_get_cbounds
    ! made by Z. Meliani 10/02/2018


  !> Calculate fast magnetosonic wave speed
  subroutine srhd_get_csound_prim(w,x,ixI^L,ixO^L,idim,csound)
    ! made by Z. Meliani 10/02/2018
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L, idim
    double precision, intent(in) :: w(ixI^S, nw), x(ixI^S,1:ndim)
    double precision, intent(out):: csound(ixI^S)
    double precision             :: inv_rho(ixO^S)

    inv_rho=1.d0/w(ixO^S,rho_)

    if(srhd_energy) then
      csound(ixO^S)=dsqrt(srhd_gamma*w(ixO^S,p_)*inv_rho)
    else
      csound(ixO^S)=dsqrt(srhd_gamma*srhd_adiab*w(ixO^S,rho_)**gamma_1)
    end if

  end subroutine srhd_get_csound_prim




  !> Calculate the square of the thermal sound speed csound2 within ixO^L.
  !> csound2=gamma*p/rho
  subroutine srhd_get_csound2(w,x,ixI^L,ixO^L,csound2)
    ! made by Z. Meliani 13/02/2018
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: w(ixI^S,nw),rhoh(ixO^S)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: csound2(ixO^S)

    double precision                :: rho(ixO^S),rhoh(ixO^S)

    if(srhd_energy) then
      rho  = w(ixO^S,d_)/w(ixO^S,lfac_)
      rhoh = (w(ixO^S,e_)+w(ixO^S,d_)+w(ixO^S,p_))/w(ixO^S,lfac_)
      call srhd_get_csound2_eos(ixI^L,ixO^L,x,rho,w(ixI^S,p_),rhoh,csound2)
    else
      csound2(ixO^S)=srhd_gamma*srhd_adiab*rho(ixO^S)**gamma_1
    end if
  end subroutine srhd_get_csound2



  !> Calculate total pressure within ixO^L
  subroutine srhd_get_p_total(w,x,ixI^L,ixO^L,p)
    ! made by Z. Meliani 10/02/2018
    use mod_global_parameters
    implicit none

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: w(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: p(ixI^S)

    call srhd_get_pthermal(w,x,ixI^L,ixO^L,p)


  end subroutine srhd_get_p_total






  !> Calculate fluxes within ixO^L.
  subroutine srhd_get_flux(wC,wP,x,ixI^L,ixO^L,idim,f)
    ! made by Z. Meliani 10/02/2018
    use mod_global_parameters
    implicit none
    integer, intent(in)          :: ixI^L, ixO^L, idim
    ! conservative w
    double precision, intent(in) :: wC(ixI^S,nw)
    ! primitive w
    double precision, intent(in) :: wP(ixI^S,nw)
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision,intent(out) :: f(ixI^S,nwflux)

    double precision             :: v(ixI^S,1:ndir)
    integer                      :: idirmin, iw, idir



    call srhd_get_v_prim(wP,x,ixI^L,ixO^L,v)
    ! Get flux of density
    f(ixO^S,rho_)=wP(ixO^S,mom(idim))*wP(ixO^S,rho_)

    ! Get flux of tracer
    do iw=1,srhd_n_tracer
      f(ixO^S,tracer(iw))=wP(ixO^S,mom(idim))*wP(ixO^S,tracer(iw))
    end do

    ! Get flux of momentum
    ! f_i[m_k]=v_i*m_k [+ptotal if i==k]
    Loop_idir_mom : do idir=1,ndir
      f(ixO^S,mom(idir))= v(ixO^S,idir)*wC(ixO^S,mom(idim))
    end do Loop_idir_mom
    f(ixO^S,mom(idim))=wP(ixO^S,pp_)+f(ixO^S,mom(idim))

    ! Get flux of energy
    ! f_i[e]=v_i*e+v_i*ptotal
    is_notiso : if (srhd_energy) then
       is_internal : if (block%e_is_internal) then
          f(ixO^S,e_)=v(ixO^S,idim)*wP(ixO^S,pp_)
       else is_internal
          f(ixO^S,e_)=v(ixO^S,idim)*(wC(ixO^S,e_) + wP(ixO^S,pp_))
       end if is_internal

    end if is_notiso



  end subroutine srhd_get_flux

  !> w[iws]=w[iws]+qdt*S[iws,wCT] where S is the source based on wCT within ixO
  subroutine srhd_add_source(qdt,ixI^L,ixO^L,wCT,w,x,qsourcesplit,active)
    ! made by Z. Meliani 10/02/2018
    use mod_global_parameters
    use mod_radiative_cooling, only: radiative_cooling_add_source
    use mod_viscosity, only: viscosity_add_source
    use mod_gravity, only: gravity_add_source

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    logical, intent(in)             :: qsourcesplit
    logical, intent(inout)            :: active

    if (.not. qsourcesplit) then
      ! Source for solving internal energy
      if (srhd_energy .and. block%e_is_internal) then
        active = .true.
        call internal_energy_add_source(qdt,ixI^L,ixO^L,wCT,w,x)
      endif


    end if


!    if(srhd_radiative_cooling) then
!      call radiative_cooling_add_source(qdt,ixI^L,ixO^L,wCT,&
!           w,x,qsourcesplit,active)
!    end if

!    if(srhd_viscosity) then
!      call viscosity_add_source(qdt,ixI^L,ixO^L,wCT,&
!           w,x,srhd_energy,qsourcesplit,active)
!    end if

    if(srhd_gravity) then
      call gravity_add_source(qdt,ixI^L,ixO^L,wCT,&
           w,x,srhd_energy,qsourcesplit,active)
    end if
    call srhd_get_auxiliary(.true.,w,x,ixI^L,ixO^L,'srhd_add_source')

  end subroutine srhd_add_source

  subroutine internal_energy_add_source(qdt,ixI^L,ixO^L,wCT,w,x)
    ! made by Z. Meliani 10/02/2018
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision                :: pth(ixI^S),v(ixI^S,1:ndir),divv(ixI^S)

    call srhd_get_v(wCT,x,ixI^L,ixI^L,v)
    call divvector(v,ixI^L,ixO^L,divv)
    call srhd_get_pthermal(wCT,x,ixI^L,ixO^L,pth)
    w(ixO^S,e_)=w(ixO^S,e_)-qdt*pth(ixO^S)*divv(ixO^S)

  end subroutine internal_energy_add_source







  !> If resistivity is not zero, check diffusion time limit for dt
  subroutine srhd_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
    use mod_global_parameters
    use mod_usr_methods
    use mod_radiative_cooling, only: cooling_get_dt
    use mod_viscosity, only: viscosity_get_dt
    use mod_gravity, only: gravity_get_dt

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: dtnew
    double precision, intent(in)    :: dx^D
    double precision, intent(in)    :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)

    integer                       :: idirmin,idim
    double precision              :: dxarr(ndim)
    double precision              :: current(ixI^S,7-2*ndir:3),eta(ixI^S)

    dtnew = bigdouble

    ^D&dxarr(^D)=dx^D;


    if(srhd_gravity) then
      call gravity_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
    end if

  end subroutine srhd_get_dt

  ! Add geometrical source terms to w
  subroutine srhd_add_source_geom(qdt,ixI^L,ixO^L,wCT,w,x)
    ! made by Z. Meliani 10/02/2018
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: wCT(ixI^S,1:nw), w(ixI^S,1:nw)

    integer          :: iw,idir, h1x^L{^NOONED, h2x^L}
    double precision,dimension(ixI^S) :: tmp,ptot,tmp2,xi
    integer          :: mr_,mphi_ ! Polar var. names

    character(len=30)  :: subname_loc

    cond_noslab: if(typeaxial /= 'slab')then
     subname_loc='srhd_add_geom'
     ! get auxiliary variables
     call srhd_get_auxiliary(.true.,wCT,x,ixI^L,ixO^L,subname_loc)

     call srhd_get_pthermal(wCT,x,ixI^L,ixO^L,ptot)


     mr_=mom(1); mphi_=mom(1)-1+phi_  ! Polar var. names
     xi(ixO^S) = wCT(ixO^S,e_)+w(ixO^S,p_)+w(ixO^S,d_)


    select case (typeaxial)
    case ('cylindrical')
      if (angmomfix) then
        call mpistop("angmomfix not implemented yet in MHD")
      endif
       if(phi_>0) then
         w(ixO^S,mr_)=w(ixO^S,mr_)+qdt/x(ixO^S,1)*(ptot(ixO^S)+&
                   (wCT(ixO^S,mphi_)**2)&
                   /xi(ixO^S))

         w(ixO^S,mphi_)=w(ixO^S,mphi_)+qdt/x(ixO^S,1)*(&
                  -(wCT(ixO^S,mphi_)*wCT(ixO^S,mr_))&
                  /xi(ixO^S))


       else
         w(ixO^S,mr_)=w(ixO^S,mr_)+qdt/x(ixO^S,1)*ptot(ixO^S)
       end if
    case ('spherical')
       h1x^L=ixO^L-kr(1,^D); {^NOONED h2x^L=ixO^L-kr(2,^D);}

       ! m1
       tmp(ixO^S)=ptot(ixO^S)*x(ixO^S,1) &
         *(block%surfaceC(ixO^S,1)-block%surfaceC(h1x^S,1))/block%dvolume(ixO^S)
       if(ndir>1) then
         do idir=2,ndir
           tmp(ixO^S)=tmp(ixO^S)+wCT(ixO^S,mom(idir))**2.0&
                                /xi(ixO^S)
         end do
       end if
       w(ixO^S,mom(1))=w(ixO^S,mom(1))+qdt*tmp(ixO^S)/x(ixO^S,1)


       {^NOONED
       ! m2
       !tmp(ixO^S)=tmp1(ixO^S)
       ! This will make hydrostatic p=const an exact solution
       w(ixO^S,mom(2))=w(ixO^S,mom(2))+qdt*ptot(ixO^S) &
            *(block%surfaceC(ixO^S,2)-block%surfaceC(h2x^S,2)) &
            /block%dvolume(ixO^S)

       tmp(ixO^S)=-wCT(ixO^S,mom(1))*wCT(ixO^S,mom(2))&
                    /xi(ixO^S)

       if(ndir==3) then
         tmp(ixO^S)=tmp(ixO^S)+wCT(ixO^S,mom(3))**2/xi(ixO^S) &
              *dcos(x(ixO^S,2))/dsin(x(ixO^S,2))
       end if
       w(ixO^S,mom(2))=w(ixO^S,mom(2))+qdt*tmp(ixO^S)/x(ixO^S,1)

       }

       if(ndir==3) then
         ! m3
         if(.not.angmomfix) then
           tmp(ixO^S)=-(wCT(ixO^S,mom(3))*wCT(ixO^S,mom(1))&
                       /xi(ixO^S)) {^NOONED &
                -(wCT(ixO^S,mom(2))*wCT(ixO^S,mom(3))/xi(ixO^S)) &
                *dcos(x(ixO^S,2))/dsin(x(ixO^S,2)) }
           w(ixO^S,mom(3))=w(ixO^S,mom(3))+qdt*tmp(ixO^S)/x(ixO^S,1)
         else
           call mpistop("angmomfix not implemented yet in SRHD")
         end if

       end if
    end select
  end if  cond_noslab
  end subroutine srhd_add_source_geom



  !> compute kinetic energy from primitive
  function srhd_kin_en_primitive(w, ixI^L, ixO^L, inv_rho) result(ke)
    ! made by Z. Meliani 10/02/2018
    use mod_global_parameters, only: nw, ndim
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S, nw)
    double precision              :: ke(ixO^S)
    double precision, intent(in), optional :: inv_rho(ixO^S)

    if (present(inv_rho)) then
       ke = (w(ixO^S,xi_)-1.0d0) * inv_rho
    else
       ke = (w(ixO^S,xi_)-1.0d0) / w(ixO^S, rho_)
    end if
  end function srhd_kin_en_primitive





   !> subtroutine srhd_get_auxiliary calcule using srhd_con2prim to calculate the enthalpy and the lorentz factor
   subroutine srhd_get_auxiliary(clipping,w,x,ixI^L,ixO^L,subname)
    ! made by Z. Meliani 10/02/2018
    use mod_global_parameters
    use mod_srhd_con2prim
    implicit none

    logical, intent(in)             :: clipping
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    character(len=*), intent(in)    :: subname

    integer                             :: ix^D,ierror
    integer                             :: flag_error(ixO^S)
    character(len=len(subname)+30):: subname_loc
  ! print*,', is yooou ',subname,saveigrid
   {do ix^DB=ixOmin^DB,ixOmax^DB\}
    call srhd_con2prim(w(ix^D,d_),w(ix^D,mom(1):mom(ndir)),w(ix^D,e_),&
             w(ix^D,mag(1):mag(ndir)),w(ix^D,lfac_),w(ix^D,xi_),ierror)
    if(check_small_values)then
     if(ierror/=0) then
       flag_error(ix^D) = ierror
     else
       flag_error(ix^D) = 0
     end if

    end if
   {enddo^D&\}
   is_check_small : if(check_small_values)then
    is_flag_on : if(any(flag_error(ixO^S)/=0))then
     subname_loc='srhd_get_auxiliary from -> '//trim(subname)
     call srhd_handle_small_values(.false., &
                                    w, x, &
                                    ixI^L, ixO^L,subname_loc,&
                                    flag_error=flag_error)
    end if is_flag_on
   end if is_check_small
   end subroutine srhd_get_auxiliary

   subroutine srhd_get_4u_from_3v(ixI^L,ixO^L,vtou)
    ! made by Z. Meliani 13/02/2018
    use mod_global_parameters
    implicit none
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(inout) :: vtou(ixI^S,1:ndir)

    double precision                :: lfac(ixO^S)
    integer                         :: idir


    lfac= 1.0d0/dsqrt(1.0d0-sum(vtou(ixO^S,1:ndir)**2.0,dim=ndim+1))
    Loop_idir: do idir=1,ndir
     vtou(ixO^S,idir)=lfac(ixO^S)*vtou(ixO^S,idir)
    end do Loop_idir
   end subroutine srhd_get_4u_from_3v
end module mod_srhd_phys
