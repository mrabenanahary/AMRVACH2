!> Special Relativistic Magneto-hydrodynamics module
module mod_srmhd_phys
  use mod_global_parameters, only: std_len
  use mod_srmhd_parameters
  use mod_srmhd_eos

  implicit none
  private



  ! Public methods
  public :: srmhd_phys_init
  public :: srmhd_phys_clean
  public :: srmhd_kin_en_primitive
  public :: srmhd_get_pthermal
  public :: srmhd_get_p_mag
  public :: srmhd_get_p_total
  public :: srmhd_get_v
  public :: srmhd_get_v_idim
  public :: srmhd_to_conserved
  public :: srmhd_to_primitive
  public :: srmhd_get_csound2
  public :: srmhd_get_csound_prim
  public :: get_divb
  public :: get_current
  public :: get_normalized_divb
  public :: srmhd_get_4u_from_3v
  public :: srmhd_get_3v_from_4u
contains

  !> Read this module"s parameters from a file
  subroutine srmhd_read_params(files)
    ! made by Z. Meliani 20/02/2018
    use mod_global_parameters
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /srmhd_list/ srmhd_energy, srmhd_eos,srmhd_n_tracer, srmhd_gamma,&
        srmhd_adiab, srmhd_eta, srmhd_eta_hyper, srmhd_etah, srmhd_glm_alpha,&
        srmhd_magnetofriction,srmhd_thermal_conduction,&
        srmhd_radiative_cooling, srmhd_Hall, srmhd_gravity, srmhd_viscosity,&
        srmhd_4th_order, typedivbfix, source_split_divb, divbdiff,&
        typedivbdiff, compactres, divbwave, srmhd_glm,He_abundance, SI_unit,&
        B0field,B0field_forcefree, Bdip, Bquad, Boct, Busr, srmhd_particles,&
        boundary_divbfix, boundary_divbfix_skip, srmhd_maxiterationNR,&
        srmhd_absaccNR,srmhd_tolerNr,srmhd_checkNR,srmhd_maxdspeed,small_vec2
    if(mype==0)write(*,*)'Reading srmhd_list'
    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, srmhd_list)
       close(unitpar)
    end do


  end subroutine srmhd_read_params

  !> Write this module's parameters to a snapsoht
  subroutine srmhd_write_info(fh)
    ! made by Z. Meliani 10/02/2018
    use mod_global_parameters
    integer, intent(in)                 :: fh
    integer, parameter                  :: n_par = 1
    real(kind=dp)                       :: values(n_par)
    character(len=name_len)             :: names(n_par)
    integer, dimension(MPI_STATUS_SIZE) :: st
    integer                             :: er

    call MPI_FILE_WRITE(fh, n_par, 1, MPI_INTEGER, st, er)

    names(1) = "gamma"
    values(1) = srmhd_gamma
    call MPI_FILE_WRITE(fh, values, n_par, MPI_DOUBLE_PRECISION, st, er)
    call MPI_FILE_WRITE(fh, names, n_par * name_len, MPI_CHARACTER, st, er)
  end subroutine srmhd_write_info

  subroutine srmhd_angmomfix(fC,x,wnew,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,idim)
    use mod_global_parameters
    real(kind=dp)   , intent(in)       :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    real(kind=dp)   , intent(inout)    :: fC(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nwflux,1:ndim),  wnew(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    integer, intent(in)                :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    integer, intent(in)                :: idim
    integer                            :: hxOmin1,hxOmin2,hxOmax1,hxOmax2,&
        kxCmin1,kxCmin2,kxCmax1,kxCmax2, iw
    real(kind=dp)                      :: inv_volume(ixImin1:ixImax1,&
       ixImin2:ixImax2)

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

  end subroutine srmhd_angmomfix

  subroutine srmhd_phys_init()
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
    call srmhd_read_params(par_files)

    physics_type                = "srmhd"
    phys_energy                 = srmhd_energy
    srmhd_config%dust_on        = .false.
    srmhd_config%dust_n_species = 0
    srmhd_config%ismhd          = .true.
    srmhd_config%isrel          = .true.
    srmhd_config%He_abundance   = He_abundance
    srmhd_config%n_tracer       = srmhd_n_tracer
    ! set default gamma for polytropic/isothermal process
    if(.not.srmhd_energy) srmhd_gamma=1.d0
    use_particles=srmhd_particles
    if(ndim==1) typedivbfix='none'
    select case (typedivbfix)
    case ('none')
      type_divb = divb_none
    case ('glm1')
      srmhd_glm          = .true.
      need_global_cmax = .true.
      type_divb        = divb_glm1
    case ('glm2')
      srmhd_glm          = .true.
      need_global_cmax = .true.
      need_global_vmax = .true.
      type_divb        = divb_glm2
    case ('powel', 'powell')
      type_divb = divb_powel
    case ('janhunen')
      type_divb = divb_janhunen
    case ('linde')
      type_divb = divb_linde
    case ('lindejanhunen')
      type_divb = divb_lindejanhunen
    case ('lindepowel')
      type_divb = divb_lindepowel
    case ('lindeglm')
      srmhd_glm          = .true.
      need_global_cmax = .true.
      need_global_vmax = .true.
      type_divb        = divb_lindeglm
    case default
      call mpistop('Unknown divB fix')
    end select


    !nwfluxbc=nwfluxbc+2
    call srmhd_fill_phys_indices



    nvector      = 2 ! No. vector vars
    allocate(iw_vector(nvector))
    iw_vector(1) = mom(1) - 1   ! TODO: why like this?
    iw_vector(2) = mag(1) - 1   ! TODO: why like this?

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
    if(srmhd_glm .and. ndim>1) flux_type(:,psi_)=flux_tvdlf

    srmhd_maxspeed=1.0-srmhd_maxdspeed
    srmhd_maxspeed2=srmhd_maxspeed**2.0_dp
    allocate(srmhd_iw_average(1:nw))
    srmhd_iw_average         = .true.
    srmhd_iw_average(mag(:)) = .false.

    phys_get_dt              => srmhd_get_dt
    phys_get_cmax            => srmhd_get_cmax
    phys_get_cbounds         => srmhd_get_cbounds
    phys_get_flux            => srmhd_get_flux
    phys_get_v_idim          => srmhd_get_v_idim
    phys_add_source_geom     => srmhd_add_source_geom
    phys_add_source          => srmhd_add_source
    phys_to_conserved        => srmhd_to_conserved
    phys_to_primitive        => srmhd_to_primitive
    phys_get_aux             => srmhd_get_auxiliary
    phys_check_params        => srmhd_check_params
    phys_check_w             => srmhd_check_w
    phys_get_pthermal        => srmhd_get_pthermal
    phys_boundary_adjust     => srmhd_boundary_adjust
    phys_write_info          => srmhd_write_info
    phys_angmomfix           => srmhd_angmomfix
    phys_handle_small_values => srmhd_handle_small_values
    phys_ind                  =>srmhd_ind



    ! Whether diagonal ghost cells are required for the physics
    !if(type_divb < divb_linde) phys_req_diagonal = .false.

    ! derive units from basic units
    call srmhd_physical_units()

    if(.not. srmhd_energy .and. srmhd_thermal_conduction) then
      call mpistop("thermal conduction needs srmhd_energy=T")
    end if
    if(.not. srmhd_energy .and. srmhd_radiative_cooling) then
      call mpistop("radiative cooling needs srmhd_energy=T")
    end if


    ! initialize thermal conduction module
    !if (srmhd_thermal_conduction) then
    !  phys_req_diagonal = .true.
    !  call thermal_conduction_init(srmhd_gamma)
    !end if

    ! Initialize radiative cooling module
    !if (srmhd_radiative_cooling) then
    !  call radiative_cooling_init(srmhd_gamma,He_abundance)
    !end if

    ! Initialize viscosity module
    !if (srmhd_viscosity) call viscosity_init(phys_wider_stencil,phys_req_diagonal)

    ! Initialize gravity module
    if(srmhd_gravity) then
      call gravity_init()
    end if

    ! Initialize particles module
    if(srmhd_particles) then
      call particles_init()
      phys_req_diagonal = .true.
    end if

    ! initialize magnetofriction module
    !if(srmhd_magnetofriction) then
    !  phys_req_diagonal = .true.
    !  call magnetofriction_init()
    !end if

    if(type_divb==divb_glm1) then
      ! Solve the Riemann problem for the linear 2x2 system for normal
      ! B-field and GLM_Psi according to Dedner 2002:
      phys_modify_wLR => glmSolve
    end if


    srmhd_config%dust_on        = .false.
    srmhd_config%dust_n_species = 0
    srmhd_config%energy         = srmhd_energy

    phys_iw_average => srmhd_iw_average
    phys_config     => srmhd_config
    phys_ind        => srmhd_ind
  end subroutine srmhd_phys_init
  !> deallocate array used in srmhd
  subroutine srmhd_phys_clean
  !  if(allocated(srmhd_logical_false))deallocate(srmhd_logical_false)
  !  if(allocated(srmhd_logical_true))deallocate(srmhd_logical_true)
  if(allocated(tracer)) deallocate(tracer)

  end subroutine srmhd_phys_clean

  subroutine srmhd_fill_phys_indices
    use mod_global_parameters
    implicit none
    ! .. local ..
    integer :: itr, idir
    !-------------------------------------

    ! Determine flux variables
    rho_ = var_set_rho()
    d_=rho_
    allocate(mom(ndir))
    mom(:) = var_set_momentum(ndir)

    ! Set index of energy variable
    if (srmhd_energy) then
      nwwave = 8
      e_     = var_set_energy() ! energy density
      p_     = e_               ! gas pressure
    else
      nwwave = 7
      e_     = -1
      p_     = -1
    end if

    allocate(mag(ndir))
    mag(:) = var_set_bfield(ndir)
    if (srmhd_glm) then
      psi_ = var_set_fluxvar('psi', 'psi', need_bc=.false.)
    else
      psi_ = -1
    end if

    allocate(tracer(srmhd_n_tracer))

    ! Set starting index of tracers
    do itr = 1, srmhd_n_tracer
      tracer(itr) = var_set_fluxvar("trc", "trp", itr, need_bc=.false.)
    end do

    ! Set index for auxiliary variables
    xi_  = var_set_auxvar('xi')
    lfac_=var_set_auxvar('lfac')


    allocate(srmhd_ind%mom(ndir),srmhd_ind%mag(ndir))
    srmhd_ind%rho_  =rho_
    srmhd_ind%mom(:)=mom(:)
    srmhd_ind%mag(:)=mag(:)
    srmhd_ind%e_    =e_
    srmhd_ind%pressure_   =p_
    srmhd_ind%lfac_ =lfac_
    srmhd_ind%xi_   =xi_
    srmhd_ind%psi_  =psi_
    if(srmhd_n_tracer>0)then
      allocate(srmhd_ind%tracer(srmhd_n_tracer))
      srmhd_ind%tracer(:) = tracer(:)
    end if
  end subroutine srmhd_fill_phys_indices



  !> check paramters
  subroutine srmhd_check_params
    use mod_global_parameters
    implicit none

    ! after user parameter setting
    srmhd_gamma_1          = srmhd_gamma-1.0_dp
    inv_srmhd_gamma_1      = 1.d0/srmhd_gamma_1
    gamma_to_srmhd_gamma_1 = srmhd_gamma/srmhd_gamma_1
    if (.not. srmhd_energy) then
       if (srmhd_gamma <= 0.0d0) call mpistop ("Error: srmhd_gamma <= 0")
       if (srmhd_adiab < 0.0d0) call mpistop ("Error: srmhd_adiab < 0")
       small_pressure = srmhd_adiab*small_density**srmhd_gamma
    else
       if (srmhd_gamma <= 0.0d0 .or. srmhd_gamma == 1.0d0) call mpistop &
          ("Error: srmhd_gamma <= 0 or srmhd_gamma == 1")
       small_e = small_pressure * inv_srmhd_gamma_1
    end if
    small_xi=small_density+gamma_to_srmhd_gamma_1*small_pressure

  end subroutine srmhd_check_params

  subroutine srmhd_physical_units()
    use mod_global_parameters
    real(kind=dp)    :: mp,kB,miu0
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
      if(unit_density/=1.0)then
        unit_numberdensity= unit_density/((1.d0+4.d0*He_abundance)*mp)
      else
        unit_density=(1.d0+4.d0*He_abundance)*mp*unit_numberdensity
      end if
      unit_pressure=unit_density*unit_velocity**2
      unit_temperature=unit_pressure/((2.d0+&
         3.d0*He_abundance)*unit_numberdensity*kB)
      unit_magneticfield=sqrt(miu0*unit_pressure)
      unit_time=unit_length/unit_velocity
    end if

  end subroutine srmhd_physical_units

  subroutine srmhd_check_w(primitive,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,w,flag)
    ! made by Z. Meliani 10/02/2018
    use mod_global_parameters

    logical, intent(in)            :: primitive
    integer, intent(in)            :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    real(kind=dp)   , intent(in)   :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw)
    integer, intent(inout)         :: flag(ixImin1:ixImax1,ixImin2:ixImax2)

    ! ---------------------------------------------

    flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=0
    if(primitive)then
     where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         rho_) < small_density) flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = rho_
    else
     where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, d_)/w(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2,lfac_) < small_density) flag(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2) = rho_
    end if
    cond_energy : if (srmhd_energy) then
       cond_einternal : if (block%e_is_internal) then
          where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              e_) < small_pressure*inv_srmhd_gamma_1) flag(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2) = e_
       else cond_einternal
         is_prim : if (primitive)then
           where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               p_) < small_pressure) flag(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2) = p_ !p_=e_
         else is_prim
           where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               e_) < small_e) flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = e_
           ! Calculate pressure=(gamma-1)*(e-0.5*(2ek+2eb))
           ! in srmhd we check if xi-p-d is positive
             ! do be done
             !where(tmp(ixO^S) < small_pressure) flag(ixO^S) = e_
         end if is_prim
       end if cond_einternal
    end if cond_energy

  end subroutine srmhd_check_w

  !> Transform primitive variables into conservative ones
  subroutine srmhd_to_conserved(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,w,x)
    ! made by Z. Meliani 10/02/2018
    use mod_global_parameters
    implicit none

    integer, intent(in)                :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    real(kind=dp)   , intent(inout)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
        nw)
    real(kind=dp)   , intent(in)       :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:ndim)
    ! .. local ..
    integer                            :: idir, itr
    real(kind=dp)   , dimension(ixOmin1:ixOmax1,ixOmin2:ixOmax2) :: sqrU,B2,&
       VdotB,rhoh,sqrV
    integer, dimension(ixOmin1:ixOmax1,ixOmin2:ixOmax2)          :: flag_error
    character(len=30)                  :: subname_loc
    !------------------------------------------

    subname_loc= 'srmhd_to_conserved'
    flag_error(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=0
    where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       rho_)<small_density)flag_error(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=rho_


    sqrU(ixOmin1:ixOmax1,ixOmin2:ixOmax2)    = sum(w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2, mom(:))**2, dim=ndim+1)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,lfac_) = &
       dsqrt(1.0_dp+sqrU(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
    sqrV=sqrU/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,lfac_)**2.0_dp
    ! fill the auxiliary variable xi and density D
    call srmhd_get_enthalpy(ixOmin1,ixOmin2,ixOmax1,ixOmax2,w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,rho_),w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,p_),rhoh)
    ! with enthalpy w: xi= lfac^2 rhoh

    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,xi_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       lfac_)**2.0D0*rhoh(ixOmin1:ixOmax1,ixOmin2:ixOmax2)


    ! density: d = lfac * rho
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       lfac_)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)

    ! B2 and VdotB
    call srmhd_get_B2andVdotB(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,w,.false.,B2=B2,VdotB=VdotB)
    ! Convert velocity to momentum
    ! s= (xi + B^2) * v - (v.B) * B
    ! re-use rhoh as rhoh=(xi+B^2)/lfac
    rhoh(ixOmin1:ixOmax1,ixOmin2:ixOmax2)= (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       xi_)+B2(ixOmin1:ixOmax1,ixOmin2:ixOmax2))/w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,lfac_)
    Loop_idirmom : do idir = 1, ndir
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mom(idir)) = rhoh(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2) * w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           mom(idir))-VdotB(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,mag(idir))
    end do Loop_idirmom

    cond_energy : if (srmhd_energy) then
       ! re-use sqrU=v^2 B^2 - (v.B)^2
       sqrU(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = B2(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)*sqrV(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)-VdotB(ixOmin1:ixOmax1,ixOmin2:ixOmax2)**2.0_dp
       ! sqrU should positive
       where(sqrU(ixOmin1:ixOmax1,ixOmin2:ixOmax2)<0.0_dp) &
          flag_error(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=e_
       ! re-use sqrV=xi - p -D
       sqrV(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,xi_) - w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          p_) - w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,d_)
       where(sqrV(ixOmin1:ixOmax1,ixOmin2:ixOmax2)<0.0_dp) &
          flag_error(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=e_
       ! E = xi - p +(B^2+v^2 B^2 - (v.B)^2)/2- D
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)=sqrV(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2) +0.5_dp*(B2(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2) + sqrU(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
      !if(.not.block%e_is_internal) w(ixO^S,e_)=w(ixO^S,e_) + x
      if(type_divb==divb_glm2) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_) + 0.5d0*w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,psi_)**2
    end if cond_energy


    if (check_small_values) call srmhd_handle_small_values(.false., w, x,&
        ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
       trim(subname_loc),flag_error=flag_error)


  end subroutine srmhd_to_conserved

  !> Transform conservative variables into primitive ones
  subroutine srmhd_to_primitive(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,w,x)
    ! made by Z. Meliani 10/02/2018
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    real(kind=dp)   , intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw)
    real(kind=dp)   , intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:ndim)
    ! .. local ..
    real(kind=dp)                   :: inv_rho(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2),B2(ixOmin1:ixOmax1,ixOmin2:ixOmax2),&
       VdotB(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    integer                         :: itr, idir
    character(len=30)               :: subname_loc
    !--------------------------------------------------------
    subname_loc='srmhd_to_primitive'

    ! get auxiliary variables

    call srmhd_get_auxiliary(.true.,w,x,ixImin1,ixImin2,ixImax1,ixImax2,&
       ixOmin1,ixOmin2,ixOmax1,ixOmax2,subname_loc)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_) = w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,d_)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,lfac_)

    if (srmhd_energy) then
      ! Calculate pressure = (gamma-1) * (e-ek-eb)
      if(.not.block%e_is_internal) then
        call srmhd_get_pressure_fromprimitive(ixImin1,ixImin2,ixImax1,ixImax2,&
           ixOmin1,ixOmin2,ixOmax1,ixOmax2,w,inv_rho)
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,p_)   = inv_rho(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2)
      else
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          p_)   = srmhd_gamma_1*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, e_)
      end if
    end if
    ! re-use inv_rho to store inverse of the density
    inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = 1.0_dp / w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2, rho_)

    call srmhd_get_B2andVdotB(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,w,.true.,B2=B2,VdotB=VdotB)
    Loop_idir : do idir=1,ndir
     w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(idir)) = w(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2,lfac_)*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        mom(idir))+VdotB*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        mag(idir)))/(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,xi_)+B2)
    end do Loop_idir

    if (check_small_values) call srmhd_handle_small_values(.true., w, x,&
        ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
       trim(subname_loc))

  end subroutine srmhd_to_primitive

  subroutine srmhd_handle_small_values(primitive, w, x, ixImin1,ixImin2,&
     ixImax1,ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, subname,flag_error)
    ! made by Z. Meliani 10/02/2018
    use mod_global_parameters
    use mod_small_values
    logical, intent(in)             :: primitive
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    real(kind=dp)   , intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    real(kind=dp)   , intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    character(len=*), intent(in)    :: subname
    integer, optional, intent(in)   :: flag_error(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)
    ! .. local ..
    real(kind=dp)    :: smallone
    integer          :: idir, ierror,flag(ixImin1:ixImax1,ixImin2:ixImax2)
    ! ---------------------------------------------------


    if (small_values_method == "ignore") return
    if(present(flag_error)) then
     flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = flag_error(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2)
    else
     call srmhd_check_w(primitive, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
        ixOmin2,ixOmax1,ixOmax2, w, flag)
    end if

    if (any(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2) /= 0)) then
       select case (small_values_method)
       case ("replace")
          where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2) /= 0) w(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,rho_) = small_density

          do idir = 1, ndir
             where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2) /= 0) &
                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mom(idir)) = 0.0d0
          end do

          if (srmhd_energy) then
             if(primitive) then
               smallone = small_pressure
             else
               smallone = small_e
             end if
             where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2) /= 0) &
                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_) = smallone
          end if
       case ("average")
          call small_values_average(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
             ixOmin2,ixOmax1,ixOmax2, subname, w, x, flag,ierror)
          if(small_values_force_floor.and.crash)then
            call srmhd_small_values_floor(primitive, w, x, ixImin1,ixImin2,&
               ixImax1,ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, flag)
          end if
       case default
          call small_values_error(w, x, ixImin1,ixImin2,ixImax1,ixImax2,&
              ixOmin1,ixOmin2,ixOmax1,ixOmax2, flag, subname)
       end select
    end if
  end subroutine srmhd_handle_small_values


  subroutine srmhd_small_values_floor(primitive, w, x, ixImin1,ixImin2,ixImax1,&
     ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, flag)
    ! made by Z. Meliani 10/02/2018
    use mod_global_parameters
    use mod_small_values
    logical, intent(in)             :: primitive
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    real(kind=dp)   , intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    real(kind=dp)   , intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    integer, optional, intent(in)   :: flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    ! .. local ..
    integer                         :: idir
    real(kind=dp)                   :: smallone
    !-------------------------------------------------------------
    where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2) /= 0) w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,rho_) = small_density

    do idir = 1, ndir
             where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2) /= 0) &
                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mom(idir)) = 0.0d0
    end do

    if (srmhd_energy) then
      if(primitive) then
         smallone = small_pressure
      else
        smallone = small_e
      end if
      where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2) /= 0) w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,e_) = smallone
    end if
  end subroutine srmhd_small_values_floor


 !> Calculate thermal pressure for enthalpy and density
  subroutine srmhd_get_pressure_fromprimitive(ixImin1,ixImin2,ixImax1,ixImax2,&
     ixOmin1,ixOmin2,ixOmax1,ixOmax2,w,pth)
    ! made by Z. Meliani 10/02/2018
    use mod_global_parameters
    implicit none
    integer, intent(in)            :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    real(kind=dp)   , intent(in)   :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw)
    real(kind=dp)   , intent(out)  :: pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    real(kind=dp)                  :: rhoh(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    !--------------------------

    rhoh(ixOmin1:ixOmax1,ixOmin2:ixOmax2)  = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       xi_)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,lfac_)**2.0_dp
    call srmhd_get_pressure_primitive_eos(ixImin1,ixImin2,ixImax1,ixImax2,&
       ixOmin1,ixOmin2,ixOmax1,ixOmax2,w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_),&
       rhoh,pth)

  end subroutine srmhd_get_pressure_fromprimitive


  !> Calculate thermal pressure=(gamma-1)*(e-0.5*m**2/rho-b**2/2) within ixO^L
  subroutine srmhd_get_pthermal(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,pth)
    ! made by Z. Meliani 10/02/2018
    use mod_global_parameters
    implicit none

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    real(kind=dp)   , intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw)
    real(kind=dp)   , intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    real(kind=dp)   , intent(out)   :: pth(ixImin1:ixImax1,ixImin2:ixImax2)

    real(kind=dp)                   :: rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2),&
       rhoh(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    rho        = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,d_)/w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,lfac_)
    rhoh       = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,xi_)/w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,lfac_)**2.0d0
    call srmhd_get_pthermal_eos(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2,x,rho,rhoh,w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       e_),pth)
  end subroutine srmhd_get_pthermal


 !> Calculate the square of the thermal sound speed csound2 within ixO^L.
  !> csound2=gamm*p/rho
  subroutine srmhd_get_csound2_prim(w,x,ixImin1,ixImin2,ixImax1,ixImax2,&
     ixOmin1,ixOmin2,ixOmax1,ixOmax2,csound2)
    ! made by Z. Meliani 13/02/2018
    use mod_global_parameters
    implicit none
    integer, intent(in)               :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    real(kind=dp)   , intent(in)      :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw)
    real(kind=dp)   , intent(in)      :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    real(kind=dp)   , intent(out)     :: csound2(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)

    real(kind=dp)   ,dimension(ixOmin1:ixOmax1,ixOmin2:ixOmax2) :: rhoh


    rhoh=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,xi_)/w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,lfac_)**2.0d0

    call srmhd_get_csound2_prim_eos(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2,x,w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_),rhoh,&
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,p_),csound2)
  end subroutine srmhd_get_csound2_prim

  !> Convert energy to entropy
  subroutine e_to_rhos(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
     ixOmax2,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    real(kind=dp)   ,intent(inout)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw)
    real(kind=dp)   , intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)

    if (srmhd_energy) then
      if(.not.block%e_is_internal) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          e_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          e_) - srmhd_kin_en_primitive(w, ixImin1,ixImin2,ixImax1,ixImax2,&
          ixOmin1,ixOmin2,ixOmax1,ixOmax2) - srmhd_mag_en_primitive(w, ixImin1,&
         ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2)
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, e_) = srmhd_gamma_1* &
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          rho_)**(1.0d0 - srmhd_gamma) * w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          e_)
    else
      call mpistop("e_to_rhos can not be used without energy equation!")
    end if
  end subroutine e_to_rhos

  !> Convert entropy to energy
  subroutine rhos_to_e(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
     ixOmax2,w,x)
    use mod_global_parameters
    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2
    real(kind=dp)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw)
    real(kind=dp)   , intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)

    if (srmhd_energy) then
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, e_) = w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2, rho_)**srmhd_gamma_1 * w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2, e_) * inv_srmhd_gamma_1
       if(.not.block%e_is_internal) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           e_) =w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           e_) + srmhd_kin_en_primitive(w, ixImin1,ixImin2,ixImax1,ixImax2,&
           ixOmin1,ixOmin2,ixOmax1,ixOmax2) + srmhd_mag_en_primitive(w,&
           ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2)
    else
       call mpistop("rhos_to_e can not be used without energy equation!")
    end if
  end subroutine rhos_to_e

  !> Calculate v vector
  subroutine srmhd_get_v(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,v,B2,VdotB)
    ! made by Z. Meliani 13/02/2018
    use mod_global_parameters

    integer, intent(in)           :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    real(kind=dp)   , intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
        x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    real(kind=dp)   , intent(out) :: v(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir)

    real(kind=dp)   , optional, intent(in)  :: VdotB(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2),B2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    real(kind=dp)                           :: sub_VdotB(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2),sub_B2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    integer :: idir



    is_B2_in : if(present(B2)) then
     Loop_idir_b_v1 : do idir=1,ndir
      v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir) = (w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2, mom(idir)) + VdotB*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         mag(idir)))/(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,xi_)+B2)
     end do Loop_idir_b_v1
    else  is_B2_in
     call srmhd_get_B2andVdotB(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
        ixOmax1,ixOmax2,w,.true.,B2=sub_B2,VdotB=sub_VdotB)
     Loop_idir_b_v2: do idir=1,ndir
      v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir) = (w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2, mom(idir))+sub_VdotB*w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,mag(idir)))/(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         xi_)+sub_B2)
     end do Loop_idir_b_v2
    end if is_B2_in
  end subroutine srmhd_get_v



  !> Calculate v component
  subroutine srmhd_get_v_idim(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,idim,v_idim)
    ! made by Z. Meliani 10/02/2018
    use mod_global_parameters

    integer, intent(in)                    :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2, idim
    real(kind=dp)   , intent(in)           :: w(ixImin1:ixImax1,&
       ixImin2:ixImax2,nw), x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    real(kind=dp)   , intent(out)          :: v_idim(ixImin1:ixImax1,&
       ixImin2:ixImax2)

    real(kind=dp)                          :: sub_VdotB(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2),sub_B2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

     call srmhd_get_B2andVdotB(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
        ixOmax1,ixOmax2,w,.true.,B2=sub_B2,VdotB=sub_VdotB)
     v_idim(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = (w(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2, mom(idim)) + sub_VdotB*w(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2,mag(idim)))/(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        xi_)+sub_B2)
  end subroutine srmhd_get_v_idim

  !> Calculate v component
  subroutine srmhd_get_v_idim_loc(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,idim,v_idim,B2,VdotB)
    ! made by Z. Meliani 13/02/2018
    use mod_global_parameters

    integer, intent(in)                    :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2, idim
    real(kind=dp)   , intent(in)           :: w(ixImin1:ixImax1,&
       ixImin2:ixImax2,nw), x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    real(kind=dp)   , intent(out)          :: v_idim(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)
    real(kind=dp)   , optional, intent(in) :: VdotB(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2),B2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    real(kind=dp)                          :: sub_VdotB(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2),sub_B2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    if(present(B2))then
     v_idim(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = (w(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2, mom(idim)) + VdotB*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        mag(idim)))/(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,xi_)+B2)
    else
     call srmhd_get_B2andVdotB(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
        ixOmax1,ixOmax2,w,.true.,B2=sub_B2,VdotB=sub_VdotB)
     v_idim(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = (w(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2, mom(idim)) + sub_VdotB*w(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2,mag(idim)))/(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        xi_)+sub_B2)
    end if
  end subroutine srmhd_get_v_idim_loc

 !> Calculate v^2
  subroutine srmhd_get_v2(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,v2,B2,VdotB)
    ! made by Z. Meliani 13/02/2018
    use mod_global_parameters

    integer, intent(in)           :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    real(kind=dp)   , intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw),&
        x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    real(kind=dp)   , intent(out) :: v2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    real(kind=dp)   , optional, intent(in)  :: VdotB(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2),B2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    real(kind=dp)                           :: sub_VdotB(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2),sub_B2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    real(kind=dp)                           :: v_num(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,1:ndir)
    integer :: idir

    if(present(B2)) then
     Loop_idir_B2: do idir=1,ndir
       v_num(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir)=w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2, mom(idir)) + VdotB*w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,mag(idir))
     end do  Loop_idir_B2
     v2(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = sum(v_num(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2, 1:ndir)**2.0 ,dim=ndim+1)/(w(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2,xi_)+B2)**2.0
    else
     call srmhd_get_B2andVdotB(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
        ixOmax1,ixOmax2,w,.true.,B2=sub_B2,VdotB=sub_VdotB)
     Loop_idir_noB2: do idir=1,ndir
       v_num(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir)=w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2, mom(idir)) + sub_VdotB*w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,mag(idir))
     end do  Loop_idir_noB2

     v2(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = sum(v_num(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2, :)**2.0 ,dim=ndim+1)/(w(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2,xi_)+sub_B2)**2.0
    end if
  end subroutine srmhd_get_v2
  !========================================================================
  !> subroutine that calculate the characteristic speed
  subroutine srmhd_get_cmax_gammie(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,idim,x,w,rhoh,B2,VdotB,calfven,vidim,csound2,cmax,&
     patch_gammie,from_cbound,cmin)
    ! made by Z. Meliani 10/02/2018
    use mod_global_parameters
    implicit none
    integer, intent(in)          :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, idim
    real(kind=dp)   , intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw),&
        x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    real(kind=dp)   , intent(in), dimension(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2) :: rhoh,B2,VdotB&
                                                     ,vidim,calfven, csound2
    real(kind=dp)   , intent(inout)          :: cmax(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    real(kind=dp)   , optional, intent(inout):: cmin(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    logical, optional, intent(in)            :: from_cbound
    logical, optional, intent(in)            :: patch_gammie(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)
    ! .. local ..
    real(kind=dp)   , dimension(ixOmin1:ixOmax1,ixOmin2:ixOmax2):: A,B
    real(kind=dp)   , dimension(ixOmin1:ixOmax1,ixOmin2:ixOmax2):: v2
    !-----------------------------------------------------------


   is_patch : if(present(patch_gammie))then
    where(patch_gammie(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
     A(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = csound2(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2)+calfven(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2)-csound2(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2)*calfven(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
     B(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=vidim(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2)**2.0d0
    end where
    !if(.not.present(from_cbound))&
    !       where(patch_gammie(ixO^S))vidim(ixO^S)=dabs(vidim(ixO^S))

    cond_onedir_patch: if(ndir==1)then
     where(patch_gammie(ixOmin1:ixOmax1,ixOmin2:ixOmax2))cmax(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2)=(vidim(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2)+dsqrt(A))/(1.0d0+dsqrt(A*B))
     if(present(cmin))where(patch_gammie(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2))cmin(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2)=(vidim(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2)-dsqrt(A))/(1.0d0+dsqrt(A*B))
    else cond_onedir_patch

     call srmhd_get_v2(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
        ixOmax1,ixOmax2,v2,B2=B2,VdotB=VdotB)
     where(patch_gammie(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
      cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(vidim(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)*(1.0d0-A)+dsqrt(A*(1.0d0-v2(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2))*((1.0d0-v2(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)*A)-B(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)*(1.0d0-A))))/( 1.0d0-v2(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)*A)
     end where
     if(present(cmin))where(patch_gammie(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2))cmin(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2)=(vidim(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2)*(1.0d0-A)-dsqrt(A*(1.0d0-v2(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2))*(1.0d0-v2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*A&
                         -B(ixOmin1:ixOmax1,&
                            ixOmin2:ixOmax2)*(1.0d0-A))))/( &
                            1.0d0-v2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*A)
    end if cond_onedir_patch

   else is_patch


    A = csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)+calfven(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)-csound2(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)*calfven(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    ! use B to save vidim**2
  !  if(.not.present(from_cbound).or..not.present(cmin))&
  !                                      vidim(ixO^S)=dabs(vidim(ixO^S))

    B(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=vidim(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)**2.0d0

    cond_onedir: if(ndir==1)then
     cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(vidim(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2)+dsqrt(A))/(1.0d0+dsqrt(A*B))
     if(present(cmin))cmin(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2)=(vidim(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2)-dsqrt(A))/(1.0d0+dsqrt(A*B))
    else cond_onedir

     call srmhd_get_v2(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
        ixOmax1,ixOmax2,v2,B2=B2,VdotB=VdotB)

     cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(vidim(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2)*(1.0_dp-A)+dsqrt(A*(1.0_dp-v2)*((1.0_dp-v2*A)-&
        B*(1.0_dp-A))))/( 1.0_dp-v2*A)

     if(present(cmin))cmin(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2)=(vidim(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2)*(1.0d0-A)-dsqrt(A*(1.0d0-v2(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2))*(1.0d0-v2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*A&
                         -B(ixOmin1:ixOmax1,&
                            ixOmin2:ixOmax2)*(1.0d0-A))))/( &
                            1.0d0-v2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*A)
    end if cond_onedir

   end if is_patch

  end subroutine srmhd_get_cmax_gammie


  !> Calculate cmax_idim using gammie method within ixO^L
  subroutine srmhd_get_cmax(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,idim,cmax)
    ! made by Z. Meliani 13/02/2018
    use mod_global_parameters

    integer, intent(in)               :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2, idim
    real(kind=dp)   , intent(in)      :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
        nw), x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    real(kind=dp)   , intent(inout)   :: cmax(ixImin1:ixImax1,ixImin2:ixImax2)
    ! .. local ..
    real(kind=dp)   , dimension(ixOmin1:ixOmax1,ixOmin2:ixOmax2):: rhoh,B2,&
       VdotB,calfven2,vidim
    real(kind=dp)   , dimension(ixOmin1:ixOmax1,ixOmin2:ixOmax2):: v2,csound2
    !-------------------------------------------------------


    rhoh=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,xi_)/w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,lfac_)**2.0d0
    call srmhd_get_csound2(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,rhoh,csound2)
    call srmhd_get_B2andVdotB(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,w,.true.,B2=B2,VdotB=VdotB)
    call srmhd_get_calfven2(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,x,w,.true.,calfven2,rhoh=rhoh,B2=B2)
    call srmhd_get_v_idim_loc(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2,idim,vidim,B2=B2,VdotB=VdotB)

    call srmhd_get_cmax_gammie(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,idim,x,w,rhoh,B2,VdotB,calfven2,dabs(vidim),csound2,&
       cmax)


  end subroutine srmhd_get_cmax

  !> Estimating bounds for the minimum and maximum signal velocities
  subroutine srmhd_get_cbounds(wLC,wRC,wLp,wRp,x,ixImin1,ixImin2,ixImax1,&
     ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,idim,cmax,cmin)
    ! made by Z. Meliani 13/02/2018
    use mod_global_parameters
    implicit none
    integer, intent(in)                       :: ixImin1,ixImin2,ixImax1,&
       ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, idim
    real(kind=dp)   , intent(in)              :: wLC(ixImin1:ixImax1,&
       ixImin2:ixImax2, nw), wRC(ixImin1:ixImax1,ixImin2:ixImax2, nw)
    real(kind=dp)   , intent(in)              :: wLp(ixImin1:ixImax1,&
       ixImin2:ixImax2, nw), wRp(ixImin1:ixImax1,ixImin2:ixImax2, nw)
    real(kind=dp)   , intent(in)              :: x(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:ndim)
    real(kind=dp)   , intent(inout)           :: cmax(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    real(kind=dp)   , intent(inout), optional :: cmin(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    ! .. local ..
    real(kind=dp)                      :: wmean(ixImin1:ixImax1,&
       ixImin2:ixImax2,nw)
    real(kind=dp)   , dimension(ixOmin1:ixOmax1,ixOmin2:ixOmax2) :: umean,&
        dmean, csound2Lp, csound2Rp, tmp1,tmp2,tmp3, B2Lp,B2Rp,VdotBLp,VdotBRp,&
        rhohLp,rhohRp, cmaxL,cmaxR,csound2,rhoh,B2,VdotB, vidim,vidimLp,&
       vidimRp,calfvenLp,calfvenRp,calfven
    character(len=30)                  :: subname_loc

    subname_loc='srmhd_get_cbounds'
    if (typeboundspeed/='cmaxmean') then
      ! This implements formula (10.52) from "Riemann Solvers and Numerical
      ! Methods for Fluid Dynamics" by Toro.
      rhohLp=wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,xi_)/wLp(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,lfac_)**2.0d0
      rhohRp=wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,xi_)/wRp(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,lfac_)**2.0d0

      call srmhd_get_csound2(wLp,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2,rhohLp,csound2Lp)
      call srmhd_get_csound2(wRp,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2,rhohRp,csound2Rp)
      call srmhd_get_B2andVdotB(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2,wLp,.false.,B2=B2Lp,VdotB=VdotBLp)
      call srmhd_get_B2andVdotB(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2,wRp,.false.,B2=B2Rp,VdotB=VdotBRp)

      call srmhd_get_calfven2(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,x,wLp,.false.,calfvenLp,rhoh=rhohLp,B2=B2Lp)
      call srmhd_get_calfven2(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,x,wRp,.false.,calfvenRp,rhoh=rhohRp,B2=B2Rp)

      call srmhd_get_v_idim_loc(wLp,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2,idim,vidimLp,B2=B2Lp,VdotB=VdotBLp)
      call srmhd_get_v_idim_loc(wRp,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2,idim,vidimRp,B2=B2Rp,VdotB=VdotBRp)

      tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=sqrt(wLp(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,xi_)+B2Lp(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
      tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=sqrt(wRp(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,xi_)+B2Rp(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
      tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=1.d0/(sqrt(wLp(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,xi_)+B2Lp(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2))+sqrt(wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         rho_)+B2Rp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)))

      umean(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(wLp(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,mom(idim))*tmp1(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)+wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         mom(idim))*tmp2(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2))*tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

     call srmhd_get_cmax_gammie(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
        ixOmin2,ixOmax1,ixOmax2,idim,x,wLp,rhoh=rhohLp,B2=B2Lp,VdotB=VdotBLp,&
        calfven=calfvenLp,vidim=dabs(vidimLp),csound2=csound2Lp,cmax=cmaxL)

     call srmhd_get_cmax_gammie(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
        ixOmin2,ixOmax1,ixOmax2,idim,x,wRp,rhoh=rhohRp,B2=B2Rp,VdotB=VdotBRp,&
        calfven=calfvenRp,vidim=dabs(vidimRp),csound2=csound2Rp,cmax=cmaxR)


      dmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(tmp1(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)*cmaxL(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)**2+tmp2(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)*cmaxR(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)**2)*tmp3(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)+0.5d0*tmp1(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)*tmp2(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)*tmp3(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)**2*(wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         mom(idim))-wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(idim)))**2
      dmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=sqrt(dmean(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2))
      if(present(cmin)) then
        cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=umean(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2)-dmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=umean(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2)+dmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
      else
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=abs(umean(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2))+dmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
      end if
    else
      wmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nw)=0.5_dp*(wLC(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,1:nw)+wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nw))
      ! get auxiliary variables

      call srmhd_get_auxiliary(.true.,wmean,x,ixImin1,ixImin2,ixImax1,ixImax2,&
         ixOmin1,ixOmin2,ixOmax1,ixOmax2,subname_loc)
      rhoh=wmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,xi_)/wmean(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,lfac_)**2.0d0

      call srmhd_get_csound2(wmean,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2,rhoh,csound2)
      call srmhd_get_B2andVdotB(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2,wmean,.true.,B2=B2,VdotB=VdotB)
      call srmhd_get_calfven2(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,x,wmean,.true.,calfven,rhoh=rhoh,B2=B2)
      call srmhd_get_v_idim_loc(wmean,x,ixImin1,ixImin2,ixImax1,ixImax2,&
         ixOmin1,ixOmin2,ixOmax1,ixOmax2,idim,vidim,B2=B2,VdotB=VdotB)
      call srmhd_get_cmax_gammie(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2,idim,x,wmean,rhoh=rhoh,B2=B2,VdotB=VdotB,&
         calfven=calfven,vidim=vidim,csound2=csound2,cmax=cmax,&
         from_cbound=.true.,cmin=cmin)


      if(present(cmin)) then
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=min(max(cmax(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2),0.0D0),1.0d0)
        cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=max(min(cmin(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2),0.0D0),-1.0d0)
      else
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=min(cmax(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2),1.0d0)
      end if
    end if
  end subroutine srmhd_get_cbounds
    ! made by Z. Meliani 10/02/2018
  !> Calculate the Aflven speed
  subroutine srmhd_get_calfven2(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,x,w,conserve,calfven,rhoh,B2)
   use mod_global_parameters
   implicit none
   integer, intent(in)           :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
      ixOmin2,ixOmax1,ixOmax2
   real(kind=dp)   , intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2, 1/nw),&
       x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
   logical,          intent(in)  :: conserve
   real(kind=dp)   , intent(out) :: calfven(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
   real(kind=dp)   , optional    :: rhoh(ixOmin1:ixOmax1,ixOmin2:ixOmax2),&
      B2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

   real(kind=dp)                 :: sub_B2(ixOmin1:ixOmax1,ixOmin2:ixOmax2),&
      sub_rhoh(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

   if(present(B2))then
    calfven=B2/(B2+rhoh)
   else
    sub_rhoh=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,xi_)/w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,lfac_)**2.0d0
    call srmhd_get_B2andVdotB(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,w,conserve,B2=sub_B2)
    calfven=sub_B2/(sub_B2+sub_rhoh)
   end if
  end subroutine srmhd_get_calfven2


  !> Calculate fast magnetosonic wave speed
!  subroutine srmhd_get_cfast(w,x,ixI^L,ixO^L,idim,csound)
!    use mod_global_parameters

!    integer, intent(in)          :: ixI^L, ixO^L, idim
!    real(kind=dp)   , intent(in) :: w(ixI^S, nw), x(ixI^S,1:ndim)
!    real(kind=dp)   , intent(out):: csound(ixI^S)
!    real(kind=dp)    :: cfast2(ixI^S), AvMinCs2(ixI^S), b2(ixI^S), kmax
!    real(kind=dp)    :: inv_rho(ixO^S)

!    inv_rho=1.d0/w(ixO^S,rho_)

!    call srmhd_get_csound2(w,x,ixI^L,ixO^L,csound)
!    ! store |B|^2 in v
!    b2(ixO^S)        = srmhd_mag_en_all(w,ixI^L,ixO^L)
!    cfast2(ixO^S)   = b2(ixO^S) * inv_rho+csound(ixO^S)
!    AvMinCs2(ixO^S) = cfast2(ixO^S)**2-4.0d0*csound(ixO^S) &
!         * srmhd_mag_i_all(w,ixI^L,ixO^L,idim)**2 &
!         * inv_rho

!    where(AvMinCs2(ixO^S)<zero)
!       AvMinCs2(ixO^S)=zero
!    end where

!    AvMinCs2(ixO^S)=sqrt(AvMinCs2(ixO^S))
!
!    if (.not. MHD_Hall) then
!       csound(ixO^S) = sqrt(half*(cfast2(ixO^S)+AvMinCs2(ixO^S)))
!    else
!       ! take the Hall velocity into account:
!       ! most simple estimate, high k limit:
!       ! largest wavenumber supported by grid: Nyquist (in practise can reduce by some factor)
!       kmax = dpi/min({dxlevel(^D)},bigdouble)*half
!       csound(ixO^S) = max(sqrt(half*(cfast2(ixO^S)+AvMinCs2(ixO^S))), &
!            srmhd_etah * sqrt(b2(ixO^S))*inv_rho*kmax)
!    end if
!
!  end subroutine srmhd_get_csound

  !> Calculate fast magnetosonic wave speed
  subroutine srmhd_get_csound_prim(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,idim,csound)
    ! made by Z. Meliani 10/02/2018
    use mod_global_parameters

    integer, intent(in)          :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, idim
    real(kind=dp)   , intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw),&
        x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    real(kind=dp)   , intent(out):: csound(ixImin1:ixImax1,ixImin2:ixImax2)
    real(kind=dp)    :: cfast2(ixImin1:ixImax1,ixImin2:ixImax2),&
        AvMinCs2(ixImin1:ixImax1,ixImin2:ixImax2), b2(ixImin1:ixImax1,&
       ixImin2:ixImax2), kmax
    real(kind=dp)    :: inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    inv_rho=1.d0/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)

    if(srmhd_energy) then
      csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=dsqrt(srmhd_gamma*w(&
         ixOmin1:ixOmax1,ixOmin2:ixOmax2,p_)*inv_rho)
    else
      csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=dsqrt(srmhd_gamma*srmhd_adiab*w(&
         ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)**srmhd_gamma_1)
    end if

  end subroutine srmhd_get_csound_prim


 !> Calculate magnetic pressure within ixO^L including magnetic pressure
  subroutine srmhd_get_p_mag(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,pmag)
    ! made by Z. Meliani 10/02/2018
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    real(kind=dp)   , intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw)
    real(kind=dp)   , intent(out)   :: pmag(ixImin1:ixImax1,ixImin2:ixImax2)

    pmag(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = 0.5d0 * (sum(w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2, mag(:))**2, dim=ndim+1)+sum(w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2, mag(:))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mom(:)),&
        dim=ndim+1)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       lfac_))/ w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,lfac_)

  end subroutine srmhd_get_p_mag

  !> Calculate the square of the thermal sound speed csound2 within ixO^L.
  !> csound2=gamma*p/rho
  subroutine srmhd_get_csound2(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,rhoh,csound2)
    ! made by Z. Meliani 13/02/2018
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    real(kind=dp)   , intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw),&
       rhoh(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    real(kind=dp)   , intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    real(kind=dp)   , intent(out)   :: csound2(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)

    real(kind=dp)                   :: rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    if(srmhd_energy) then
      rho=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,d_)/w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,lfac_)
      call srmhd_get_csound2_eos(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2,x,rho,rhoh,csound2)
    else
      csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=srmhd_gamma*srmhd_adiab*rho(&
         ixOmin1:ixOmax1,ixOmin2:ixOmax2)**srmhd_gamma_1
    end if
  end subroutine srmhd_get_csound2



  !> Calculate total pressure within ixO^L including magnetic pressure
  subroutine srmhd_get_p_total(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,p,B2,VdotB)
    ! made by Z. Meliani 10/02/2018
    use mod_global_parameters
    implicit none

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    real(kind=dp)   , intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw)
    real(kind=dp)   , intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    real(kind=dp)   , intent(out)   :: p(ixImin1:ixImax1,ixImin2:ixImax2)
    real(kind=dp)   , optional, intent(in) :: B2(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2),VdotB(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    real(kind=dp)                   :: pmag(ixImin1:ixImax1,ixImin2:ixImax2)
    call srmhd_get_pthermal(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2,p)
    call srmhd_get_pmag(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,pmag,B2=B2,VdotB=VdotB)
    p(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = p(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2) + pmag(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
  end subroutine srmhd_get_p_total

  subroutine srmhd_get_pmag(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,pmag,B2,VdotB)
    ! made by Z. Meliani 10/02/2018
   use mod_global_parameters
   implicit none
    integer, intent(in)                    :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    real(kind=dp)   , intent(in)           :: w(ixImin1:ixImax1,&
       ixImin2:ixImax2,nw)
    real(kind=dp)   , intent(in)           :: x(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:ndim)
    real(kind=dp)   , intent(out)          :: pmag(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    real(kind=dp)   , optional, intent(in) :: B2(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2),VdotB(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    real(kind=dp)   , dimension(ixOmin1:ixOmax1,ixOmin2:ixOmax2)     :: sub_B2,&
       sub_VdotB

    is_presentB2 : if(present(B2)) then
     pmag(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = 0.5d0*(VdotB(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2)**2.0d0  + B2(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,lfac_)**2.0d0)
    else  is_presentB2
     call srmhd_get_B2andVdotB(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
        ixOmax1,ixOmax2,w,.true.,B2=sub_B2,VdotB=sub_VdotB)
     pmag(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = 0.5d0*(sub_VdotB(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2)**2.0d0  + sub_B2(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,lfac_)**2.0d0)
    end if is_presentB2
  end subroutine srmhd_get_pmag

  !> Calculate B2 and/ot VdotB in ixO^S from conserve or primitive variables
  subroutine srmhd_get_B2andVdotB(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,w,wconserve,B2,VdotB)
    ! made by Z. Meliani 13/02/2018
    use mod_global_parameters
    implicit none
    integer, intent(in)                      :: ixImin1,ixImin2,ixImax1,&
       ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2
    ! conservative or primitive w
    real(kind=dp)   , intent(in)             :: w(ixImin1:ixImax1,&
       ixImin2:ixImax2,nw)
    logical, intent(in)                      :: wconserve
    real(kind=dp)   ,optional, intent(inout) :: B2(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2),VdotB(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    is_VdotB_present : if(present(VdotB))then
     is_conserve : if(wconserve)then
       VdotB = sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(:))*w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,mom(:)),dim=ndim+1)/w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,xi_)
     else is_conserve
       VdotB = sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(:))*w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,mom(:)),dim=ndim+1)/w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,lfac_)
     end if is_conserve
    end if is_VdotB_present
    if(present(B2))B2    = sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(:))**2.0,&
       dim=ndim+1)

  end subroutine srmhd_get_B2andVdotB


  !> Calculate fluxes within ixO^L.
  subroutine srmhd_get_flux(wC,wP,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,idim,f)
    ! made by Z. Meliani 10/02/2018
    use mod_global_parameters
    implicit none
    integer, intent(in)          :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, idim
    ! conservative w
    real(kind=dp)   , intent(in) :: wC(ixImin1:ixImax1,ixImin2:ixImax2,nw)
    ! primitive w
    real(kind=dp)   , intent(in) :: wP(ixImin1:ixImax1,ixImin2:ixImax2,nw)
    real(kind=dp)   , intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    real(kind=dp)   ,intent(out) :: f(ixImin1:ixImax1,ixImin2:ixImax2,nwflux)

    real(kind=dp)                :: ptotal(ixImin1:ixImax1,ixImin2:ixImax2)
    real(kind=dp)                :: B2(ixOmin1:ixOmax1,ixOmin2:ixOmax2),&
       VdotB(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    real(kind=dp)                :: v(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir)
    integer                      :: idirmin, iw, idir


    call srmhd_get_B2andVdotB(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,wp,.false.,B2=B2,VdotB=VdotB)
    call srmhd_get_p_total(wC,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2,ptotal,B2=B2,VdotB=VdotB)
    call srmhd_get_v(wC,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,v,B2=B2,VdotB=VdotB)
    ! Get flux of density
    !>@f[F(\rho) = u^{idims}\cdot D]@f
    f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)=wp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       mom(idim))*wp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)

    ! Get flux of tracer
    do iw=1,srmhd_n_tracer
      f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,tracer(iw))=wp(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,mom(idim))*wp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         tracer(iw))
    end do

    ! Get flux of momentum
    ! f_i[m_k]=v_i*m_k-B_k*b_i [+ptotal if i==k]
    Loop_idir_mom : do idir=1,ndir
      f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(idir))= v(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,idir)*wC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         mom(idim))-wP(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         mag(idir))*(VdotB(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*v(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,idim)+wP(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         mag(idim))/wP(ixOmin1:ixOmax1,ixOmin2:ixOmax2,lfac_)**2.0d0)

    end do Loop_idir_mom
    f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(idim))=ptotal(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)+f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(idim))

    ! Get flux of energy
    ! f_i[e]=v_i*e+v_i*ptotal-b_i*(b_k*v_k)
    is_notiso : if (srmhd_energy) then
       is_internal : if (block%e_is_internal) then
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)=v(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,idim)*wP(ixOmin1:ixOmax1,ixOmin2:ixOmax2,p_)
       else is_internal
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)=v(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,idim)*(wC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             e_) + ptotal(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2))- wP(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             mag(idim))*VdotB(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
          if(type_divb==divb_glm2) then
            f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_) = f(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,e_) + vmax_global*wP(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,psi_)*wP(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               mag(idim))
          end if

       end if is_internal
    end if is_notiso

    ! compute flux of magnetic field
    ! f_i[b_k]=v_i*b_k-v_k*b_i
    Loop_idir_Bflux : do idir=1,ndir
      is_idim_Bflux : if (idim==idir) then
        ! f_i[b_i] should be exactly 0, so we do not use the transport flux
        is_glm_Bflux : if (srmhd_glm) then
           if(type_divb==divb_glm1) then
             f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(idir))=wP(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,psi_)
           else
             f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                mag(idir))=vmax_global*wP(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                psi_)
           end if
        else is_glm_Bflux
           f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(idir))=zero
        end if is_glm_Bflux
      else is_idim_Bflux
        f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(idir))=v(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,idim)*wP(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           mag(idir))-wP(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           mag(idim))*v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir)

      end if is_idim_Bflux
    end do Loop_idir_Bflux

    if (srmhd_glm) then
      if(type_divb==divb_glm1) then
        !f_i[psi]=Ch^2*b_{i} Eq. 24e and Eq. 38c Dedner et al 2002 JCP, 175, 645
        f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           psi_)  = cmax_global**2*wP(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           mag(idim))
      else
        !f_i[psi]=Ch*b_{i} Eq. 3.16e Derigs et al 2018 JCP, 364, 420
        f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           psi_)  = vmax_global*wP(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(idim))
      end if
    end if

  end subroutine srmhd_get_flux

  !> w[iws]=w[iws]+qdt*S[iws,wCT] where S is the source based on wCT within ixO
  subroutine srmhd_add_source(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,wCT,w,x,qsourcesplit,active)
    ! made by Z. Meliani 10/02/2018
    use mod_global_parameters
    use mod_radiative_cooling, only: radiative_cooling_add_source
    use mod_viscosity, only: viscosity_add_source
    use mod_gravity, only: gravity_add_source

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    real(kind=dp)   , intent(in)    :: qdt
    real(kind=dp)   , intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    real(kind=dp)   , intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    logical, intent(in)             :: qsourcesplit
    logical, intent(inout)          :: active
    !-----------------------------------------------------


    if (.not. qsourcesplit) then
      ! Source for solving internal energy
      if (srmhd_energy .and. block%e_is_internal) then
        active = .true.
        call internal_energy_add_source(qdt,ixImin1,ixImin2,ixImax1,ixImax2,&
           ixOmin1,ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
      endif

      ! Source for B0 splitting
!      if (B0field) then
!        active = .true.
!        call add_source_B0split(qdt,ixI^L,ixO^L,wCT,w,x)
!      end if

      ! Sources for resistivity in eqs. for e, B1, B2 and B3
!      if (abs(srmhd_eta)>smalldouble)then
!        active = .true.
!        call add_source_res2(qdt,ixI^L,ixO^L,wCT,w,x)
!      end if

!      if (srmhd_eta_hyper>0.d0)then
!        active = .true.
!        call add_source_hyperres(qdt,ixI^L,ixO^L,wCT,w,x)
!      end if
    end if

      
    if(.not.source_split_divb .and. .not.qsourcesplit .and. istep==nstep) then
      ! Sources related to div B
      select case (type_divb)
      case (divb_none)
        ! Do nothing
      case (divb_glm1)
        active = .true.
        call add_source_glm1(dt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2,pw(saveigrid)%wold,w,x)
      case (divb_glm2)
        active = .true.
        call add_source_glm2(dt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2,pw(saveigrid)%wold,w,x)
      case (divb_powel)
        active = .true.
        call add_source_powel(dt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2,pw(saveigrid)%wold,w,x)
      case (divb_janhunen)
        active = .true.
        call add_source_janhunen(dt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2,pw(saveigrid)%wold,w,x)
      case (divb_linde)
        active = .true.
        call add_source_linde(dt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2,pw(saveigrid)%wold,w,x)
      case (divb_lindejanhunen)
        active = .true.
        call add_source_linde(dt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2,pw(saveigrid)%wold,w,x)
        call add_source_janhunen(dt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2,pw(saveigrid)%wold,w,x)
      case (divb_lindepowel)
        active = .true.
        call add_source_linde(dt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2,pw(saveigrid)%wold,w,x)
        call add_source_powel(dt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2,pw(saveigrid)%wold,w,x)
      case (divb_lindeglm)
        active = .true.
        call add_source_linde(dt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2,pw(saveigrid)%wold,w,x)
        call add_source_glm2(dt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2,pw(saveigrid)%wold,w,x)
      case default
        call mpistop('Unknown divB fix')
      end select
    else if(source_split_divb .and. qsourcesplit) then
      ! Sources related to div B
      select case (type_divb)
      case (divb_none)
        ! Do nothing
      case (divb_glm1)
        active = .true.
        call add_source_glm1(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
      case (divb_glm2)
        active = .true.
        call add_source_glm2(dt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2,pw(saveigrid)%wold,w,x)
      case (divb_powel)
        active = .true.
        call add_source_powel(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
      case (divb_janhunen)
        active = .true.
        call add_source_janhunen(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
      case (divb_linde)
        active = .true.
        call add_source_linde(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
      case (divb_lindejanhunen)
        active = .true.
        call add_source_linde(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
        call add_source_janhunen(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
      case (divb_lindepowel)
        active = .true.
        call add_source_linde(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
        call add_source_powel(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
      case (divb_lindeglm)
        active = .true.
        call add_source_linde(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
        call add_source_glm2(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
      case default
        call mpistop('Unknown divB fix')
      end select
    end if
   

!    if(srmhd_radiative_cooling) then
!      call radiative_cooling_add_source(qdt,ixI^L,ixO^L,wCT,&
!           w,x,qsourcesplit,active)
!    end if

!    if(srmhd_viscosity) then
!      call viscosity_add_source(qdt,ixI^L,ixO^L,wCT,&
!           w,x,srmhd_energy,qsourcesplit,active)
!    end if

    if(srmhd_gravity) then
      call gravity_add_source(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2,wCT,w,x,srmhd_energy,qsourcesplit,active)
    end if

    call srmhd_get_auxiliary(.true.,w,x,ixImin1,ixImin2,ixImax1,ixImax2,&
       ixOmin1,ixOmin2,ixOmax1,ixOmax2,'srmhd_add_source')

  end subroutine srmhd_add_source

  subroutine internal_energy_add_source(qdt,ixImin1,ixImin2,ixImax1,ixImax2,&
     ixOmin1,ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
    ! made by Z. Meliani 10/02/2018
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    real(kind=dp)   , intent(in)    :: qdt
    real(kind=dp)   , intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    real(kind=dp)   , intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    real(kind=dp)                   :: pth(ixImin1:ixImax1,ixImin2:ixImax2),&
       v(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir),divv(ixImin1:ixImax1,&
       ixImin2:ixImax2)

    call srmhd_get_v(wCT,x,ixImin1,ixImin2,ixImax1,ixImax2,ixImin1,ixImin2,&
       ixImax1,ixImax2,v)
    call divvector(v,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
       ixOmax2,divv)
    call srmhd_get_pthermal(wCT,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2,pth)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       e_)-qdt*pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*divv(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)

  end subroutine internal_energy_add_source



  subroutine add_source_glm1(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
    ! Add divB related sources to w within ixO
    ! corresponding to Dedner JCP 2002, 175, 645 _equation 24_
    ! giving the EGLM-MHD scheme
    use mod_global_parameters
    use mod_geometry

    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2
    real(kind=dp)   , intent(in) :: qdt, wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    real(kind=dp)   , intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    real(kind=dp)   :: divb(ixImin1:ixImax1,ixImin2:ixImax2)
    integer          :: idim,idir
    real(kind=dp)    :: gradPsi(ixImin1:ixImax1,ixImin2:ixImax2)

    ! We calculate now div B
    call get_divb(wCT,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
       ixOmax2,divb)

    ! dPsi/dt =  - Ch^2/Cp^2 Psi
    if (srmhd_glm_alpha < zero) then
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,psi_) = &
         abs(srmhd_glm_alpha)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,psi_)
    else
      ! implicit update of Psi variable
      ! equation (27) in Mignone 2010 J. Com. Phys. 229, 2117
      if(slab) then
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           psi_) = dexp(-qdt*cmax_global*srmhd_glm_alpha/minval(dxlevel(:)))*w(&
           ixOmin1:ixOmax1,ixOmin2:ixOmax2,psi_)
      else
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           psi_) = dexp(-qdt*cmax_global*srmhd_glm_alpha/minval(block%ds(&
           ixOmin1:ixOmax1,ixOmin2:ixOmax2,:),dim=ndim+1))*w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,psi_)
      end if
    end if

    ! gradient of Psi
    do idim=1,ndim
       select case(typegrad)
       case("central")
          call gradient(wCT(ixImin1:ixImax1,ixImin2:ixImax2,psi_),ixImin1,&
             ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,idim,&
             gradPsi)
       case("limited")
          call gradientS(wCT(ixImin1:ixImax1,ixImin2:ixImax2,psi_),ixImin1,&
             ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,idim,&
             gradPsi)
       end select
       if (srmhd_energy .and. .not.block%e_is_internal) then
       ! e  = e  -qdt (b . grad(Psi))
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_) = w(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,e_)-qdt*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            mag(idim))*gradPsi(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
       end if
    end do

    ! m = m - qdt b div b
    do idir=1,ndir
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(idir))=w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,mom(idir))-qdt*srmhd_mag_i_all(w,ixImin1,ixImin2,&
         ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
         idir)*divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    end do

    if (check_small_values) call srmhd_handle_small_values(.false.,w,x,ixImin1,&
       ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
       'add_source_glm1')

  end subroutine add_source_glm1

  subroutine add_source_glm2(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
    ! Add divB related sources to w within ixO
    ! corresponding to Eq. 3.17 Derigs et al 2018 JCP, 364, 420
    use mod_global_parameters
    use mod_geometry

    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2
    real(kind=dp)   , intent(in) :: qdt, wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    real(kind=dp)   , intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    real(kind=dp)   :: divb(ixImin1:ixImax1,ixImin2:ixImax2),v(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:ndir)
    integer          :: idim,idir
    real(kind=dp)    :: gradPsi(ixImin1:ixImax1,ixImin2:ixImax2,ndim),&
        Bf(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir)

    ! calculate now div B
    call get_divb(wCT,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
       ixOmax2,divb)

    ! calculate velocity
    call srmhd_get_v(wCT,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,v)

    ! gradient of Psi
    do idim=1,ndim
       select case(typegrad)
       case("central")
         call gradient(wCT(ixImin1:ixImax1,ixImin2:ixImax2,psi_),ixImin1,&
            ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,idim,&
            gradPsi(ixImin1:ixImax1,ixImin2:ixImax2,idim))
       case("limited")
         call gradientS(wCT(ixImin1:ixImax1,ixImin2:ixImax2,psi_),ixImin1,&
            ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,idim,&
            gradPsi(ixImin1:ixImax1,ixImin2:ixImax2,idim))
       end select
    end do

    if(B0field) then
      Bf(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir)=wCT(ixImin1:ixImax1,&
         ixImin2:ixImax2,mag(1:ndir))+block%B0(ixImin1:ixImax1,ixImin2:ixImax2,&
         1:ndir,0)
    else
      Bf(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir)=wCT(ixImin1:ixImax1,&
         ixImin2:ixImax2,mag(1:ndir))
    end if

    if (srmhd_energy .and. .not.block%e_is_internal) then
       ! e = e - qdt ( (v . b) * div b + (grad psi . v) * psi)
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          e_) - qdt * (divb(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2) * sum(v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          :)*Bf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,:),&
          dim=ndim+1) + wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          psi_)*sum(v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          1:ndim)*gradPsi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:ndim),dim=ndim+1))
    end if

    ! b_i = b_i - qdt * v_i * div b
    do idir=1,ndir
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(idir))=w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,mag(idir))-qdt*v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         idir)*divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    end do

    ! m_i = m_i - qdt * b_i * div b
    do idir=1,ndir
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(idir))=w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,mom(idir))-qdt*Bf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         idir)*divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    end do

    ! psi = psi - qdt * (v . grad(psi) + alpha * psi)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,psi_) = w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,psi_)-qdt*(sum(v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       1:ndim)*gradPsi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:ndim),&
       dim=ndim+1)+vmax_global/0.18d0*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       psi_))

    if (check_small_values) call srmhd_handle_small_values(.false.,w,x,ixImin1,&
       ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
       'add_source_glm')

  end subroutine add_source_glm2

  !> Add divB related sources to w within ixO corresponding to Powel
  subroutine add_source_powel(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    real(kind=dp)   , intent(in)    :: qdt,   wCT(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    real(kind=dp)   , intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    real(kind=dp)                   :: divb(ixImin1:ixImax1,ixImin2:ixImax2),&
       v(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir)
    integer                         :: idir

    ! We calculate now div B
    call get_divb(wCT,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
       ixOmax2,divb)

    ! calculate velocity
    call srmhd_get_v(wCT,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,v)

    if (srmhd_energy .and. .not.block%e_is_internal) then
      ! e = e - qdt (v . b) * div b
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         e_)-qdt*sum(v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,:)*wCT(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,mag(:)),dim=ndim+1)*divb(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)
    end if

    ! b = b - qdt v * div b
    do idir=1,ndir
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(idir))=w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,mag(idir))-qdt*v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         idir)*divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    end do

    ! m = m - qdt b div b
    do idir=1,ndir
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(idir))=w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,mom(idir))-qdt*srmhd_mag_i_all(w,ixImin1,ixImin2,&
         ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
         idir)*divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    end do

    if (check_small_values) call srmhd_handle_small_values(.false.,w,x,ixImin1,&
       ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
       'add_source_powel')

  end subroutine add_source_powel

  subroutine add_source_janhunen(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
    ! Add divB related sources to w within ixO
    ! corresponding to Janhunen, just the term in the induction equation.
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    real(kind=dp)   , intent(in)    :: qdt,   wCT(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    real(kind=dp)   , intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    real(kind=dp)                   :: divb(ixImin1:ixImax1,ixImin2:ixImax2)
    integer                         :: idir

    ! We calculate now div B
    call get_divb(wCT,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
       ixOmax2,divb)

    ! b = b - qdt v * div b
    do idir=1,ndir
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(idir))=w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,mag(idir))-qdt*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         mom(idir))/wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         rho_)*divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    end do

    if (check_small_values) call srmhd_handle_small_values(.false.,w,x,ixImin1,&
       ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
       'add_source_janhunen')

  end subroutine add_source_janhunen

  subroutine add_source_linde(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
    ! Add Linde's divB related sources to wnew within ixO
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    real(kind=dp)   , intent(in)    :: qdt, wCT(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    real(kind=dp)   , intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    integer :: idim, idir, ixpmin1,ixpmin2,ixpmax1,ixpmax2, i1,i2, iside
    real(kind=dp)    :: divb(ixImin1:ixImax1,ixImin2:ixImax2),&
       graddivb(ixImin1:ixImax1,ixImin2:ixImax2)
    logical, dimension(-1:1,-1:1) :: leveljump

    ! Calculate div B
    ixpmin1=ixOmin1-1;ixpmin2=ixOmin2-1;ixpmax1=ixOmax1+1;ixpmax2=ixOmax2+1;
    call get_divb(wCT,ixImin1,ixImin2,ixImax1,ixImax2,ixpmin1,ixpmin2,ixpmax1,&
       ixpmax2,divb)

    ! for AMR stability, retreat one cell layer from the boarders of level jump
    do i2=-1,1
    do i1=-1,1
      if(i1==0.and.i2==0) cycle
      if(neighbor_type(i1,i2,saveigrid)==2 .or. neighbor_type(i1,i2,&
         saveigrid)==4) then
        leveljump(i1,i2)=.true.
      else
        leveljump(i1,i2)=.false.
      end if
    end do
    end do

    ixpmin1=ixOmin1;ixpmin2=ixOmin2;ixpmax1=ixOmax1;ixpmax2=ixOmax2;
    do idim=1,ndim
      select case(idim)
       case(1)
          do iside=1,2
            i1=kr(1,1)*(2*iside-3);i2=kr(2,1)*(2*iside-3);
            if (leveljump(i1,i2)) then
              if (iside==1) then
                ixpmin1=ixOmin1-i1
              else
                ixpmax1=ixOmax1-i1
              end if
            end if
          end do
       
       case(2)
          do iside=1,2
            i1=kr(1,2)*(2*iside-3);i2=kr(2,2)*(2*iside-3);
            if (leveljump(i1,i2)) then
              if (iside==1) then
                ixpmin2=ixOmin2-i2
              else
                ixpmax2=ixOmax2-i2
              end if
            end if
          end do
       
      end select
    end do

    ! Add Linde's diffusive terms
    do idim=1,ndim
       ! Calculate grad_idim(divb)
       select case(typegrad)
       case("central")
         call gradient(divb,ixImin1,ixImin2,ixImax1,ixImax2,ixpmin1,ixpmin2,&
            ixpmax1,ixpmax2,idim,graddivb)
       case("limited")
         call gradientS(divb,ixImin1,ixImin2,ixImax1,ixImax2,ixpmin1,ixpmin2,&
            ixpmax1,ixpmax2,idim,graddivb)
       end select

       ! Multiply by Linde's eta*dt = divbdiff*(c_max*dx)*dt = divbdiff*dx**2
       if (slab) then
          graddivb(ixpmin1:ixpmax1,ixpmin2:ixpmax2)=graddivb(ixpmin1:ixpmax1,&
             ixpmin2:ixpmax2)*divbdiff/(1.0d0/dxlevel(1)**2+&
             1.0d0/dxlevel(2)**2)
       else
          graddivb(ixpmin1:ixpmax1,ixpmin2:ixpmax2)=graddivb(ixpmin1:ixpmax1,&
             ixpmin2:ixpmax2)*divbdiff /(1.0d0/block%ds(ixpmin1:ixpmax1,&
             ixpmin2:ixpmax2,1)**2+1.0d0/block%ds(ixpmin1:ixpmax1,&
             ixpmin2:ixpmax2,2)**2)
       end if

       w(ixpmin1:ixpmax1,ixpmin2:ixpmax2,mag(idim))=w(ixpmin1:ixpmax1,&
          ixpmin2:ixpmax2,mag(idim))+graddivb(ixpmin1:ixpmax1,ixpmin2:ixpmax2)

       if (srmhd_energy .and. typedivbdiff=='all' .and. &
          .not.block%e_is_internal) then
         ! e += B_idim*eta*grad_idim(divb)
         w(ixpmin1:ixpmax1,ixpmin2:ixpmax2,e_)=w(ixpmin1:ixpmax1,&
            ixpmin2:ixpmax2,e_)+wCT(ixpmin1:ixpmax1,ixpmin2:ixpmax2,&
            mag(idim))*graddivb(ixpmin1:ixpmax1,ixpmin2:ixpmax2)
       end if
    end do

    if (check_small_values) call srmhd_handle_small_values(.false.,w,x,ixImin1,&
       ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
       'add_source_linde')

  end subroutine add_source_linde

  !> Calculate div B within ixO
  subroutine get_divb(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,divb)

    use mod_global_parameters
    use mod_geometry

    integer, intent(in)                :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    real(kind=dp)   , intent(in)       :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw)
    real(kind=dp)                      :: divb(ixImin1:ixImax1,&
       ixImin2:ixImax2)

    real(kind=dp)                      :: bvec(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndir)

    bvec(ixImin1:ixImax1,ixImin2:ixImax2,:)=w(ixImin1:ixImax1,ixImin2:ixImax2,&
       mag(:))

    select case(typediv)
    case("central")
      call divvector(bvec,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,divb)
    case("limited")
      call divvectorS(bvec,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,divb)
    end select
  end subroutine get_divb

  !> get dimensionless div B = |divB| * volume / area / |B|
  subroutine get_normalized_divb(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,divb)

    use mod_global_parameters

    integer, intent(in)                :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    real(kind=dp)   , intent(in)       :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw)
    real(kind=dp)                      :: divb(ixImin1:ixImax1,&
       ixImin2:ixImax2), dsurface(ixImin1:ixImax1,ixImin2:ixImax2)

    integer :: ixAmin1,ixAmin2,ixAmax1,ixAmax2,idims

    call get_divb(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
       ixOmax2,divb)
    if(slab) then
      divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=0.5d0*abs(divb(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2))/sqrt(srmhd_mag_en_all(w,ixImin1,ixImin2,ixImax1,&
         ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2))/sum(1.d0/dxlevel(:))
    else
      ixAmin1=ixOmin1-1;ixAmin2=ixOmin2-1;
      ixAmax1=ixOmax1-1;ixAmax2=ixOmax2-1;
      dsurface(ixOmin1:ixOmax1,ixOmin2:ixOmax2)= &
         sum(block%surfaceC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,:),dim=ndim+1)
      do idims=1,ndim
        ixAmin1=ixOmin1-kr(idims,1);ixAmin2=ixOmin2-kr(idims,2)
        ixAmax1=ixOmax1-kr(idims,1);ixAmax2=ixOmax2-kr(idims,2);
        dsurface(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=dsurface(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2)+block%surfaceC(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
           idims)
      end do
      divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=abs(divb(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2))/sqrt(srmhd_mag_en_all(w,ixImin1,ixImin2,ixImax1,&
         ixImax2,ixOmin1,ixOmin2,ixOmax1,&
         ixOmax2))*block%dvolume(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)/dsurface(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    end if

  end subroutine get_normalized_divb

  !> Calculate idirmin and the idirmin:3 components of the common current array
  !> make sure that dxlevel(^D) is set correctly.
  subroutine get_current(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,idirmin,current)
    use mod_global_parameters
    use mod_geometry

    integer :: idirmin0
    integer :: ixOmin1,ixOmin2,ixOmax1,ixOmax2, idirmin, ixImin1,ixImin2,&
       ixImax1,ixImax2
    real(kind=dp)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    integer :: idir

    ! For ndir=2 only 3rd component of J can exist, ndir=1 is impossible for MHD
    real(kind=dp)    :: current(ixImin1:ixImax1,ixImin2:ixImax2,7-2*ndir:3),&
       bvec(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir)

    idirmin0 = 7-2*ndir

    bvec(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir)=w(ixImin1:ixImax1,&
       ixImin2:ixImax2,mag(1:ndir))

    call curlvector(bvec,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,current,idirmin,idirmin0,ndir)

    if(B0field) current(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       idirmin0:3)=current(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       idirmin0:3)+block%J0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idirmin0:3)

  end subroutine get_current

  !> If resistivity is not zero, check diffusion time limit for dt
  subroutine srmhd_get_dt(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,dtnew,dx1,dx2,x)
    use mod_global_parameters
    use mod_usr_methods
    use mod_radiative_cooling, only: cooling_get_dt
    use mod_viscosity, only: viscosity_get_dt
    use mod_gravity, only: gravity_get_dt

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    real(kind=dp)   , intent(inout) :: dtnew
    real(kind=dp)   , intent(in)    :: dx1,dx2
    real(kind=dp)   , intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    real(kind=dp)   , intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)

    integer                       :: idirmin,idim
    real(kind=dp)                 :: dxarr(ndim)
    real(kind=dp)                 :: current(ixImin1:ixImax1,ixImin2:ixImax2,&
       7-2*ndir:3),eta(ixImin1:ixImax1,ixImin2:ixImax2)

    dtnew = bigdouble

    dxarr(1)=dx1;dxarr(2)=dx2;


    if(srmhd_gravity) then
      call gravity_get_dt(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,dtnew,dx1,dx2,x)
    end if

  end subroutine srmhd_get_dt

  ! Add geometrical source terms to w
  subroutine srmhd_add_source_geom(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
    ! made by Z. Meliani 10/02/2018
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    real(kind=dp)   , intent(in)    :: qdt, x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    real(kind=dp)   , intent(inout) :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw), w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    integer          :: iw,idir, h1xmin1,h1xmin2,h1xmax1,h1xmax2, h2xmin1,&
       h2xmin2,h2xmax1,h2xmax2
    real(kind=dp)    :: tmp(ixImin1:ixImax1,ixImin2:ixImax2),&
       ptot(ixImin1:ixImax1,ixImin2:ixImax2),tmp2(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    real(kind=dp)    :: B2(ixOmin1:ixOmax1,ixOmin2:ixOmax2),&
       VdotB(ixOmin1:ixOmax1,ixOmin2:ixOmax2),VdotB2(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2),xiplusB2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    integer :: mr_,mphi_ ! Polar var. names
    integer :: br_,bphi_
    character(len=30)  :: subname_loc

    if(typeaxial /= 'slab')then
     subname_loc='srmhd_add_geom'
     ! get auxiliary variables
     call srmhd_get_auxiliary(.true.,wCT,x,ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2,subname_loc)
     call srmhd_get_B2andVdotB(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
        ixOmax1,ixOmax2,wCT,.true.,B2=B2,VdotB=VdotB)
     call srmhd_get_p_total(wCT,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
        ixOmin2,ixOmax1,ixOmax2,ptot)
     VdotB2=VdotB**2.0
     xiplusB2(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = wCT(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2,xi_)**2.0+B2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

     mr_=mom(1); mphi_=mom(1)-1+phi_  ! Polar var. names
     br_=mag(1); bphi_=mag(1)-1+phi_
    end if

    select case (typeaxial)
    case ('cylindrical')
      if (angmomfix) then
        call mpistop("angmomfix not implemented yet in MHD")
      endif
       if(phi_>0) then
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mr_)=w(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,mr_)+qdt/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            1)*(ptot(ixOmin1:ixOmax1,ixOmin2:ixOmax2)+(wCT(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,mphi_)**2-VdotB2(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            bphi_)**2)/xiplusB2(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2)-(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            bphi_)/wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,lfac_))**2.0)

         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mphi_)=w(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,mphi_)+qdt/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            1)*(-(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            mphi_)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            mr_)-VdotB2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*wCT(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,br_)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            bphi_))/xiplusB2(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2) +wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            bphi_)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            br_)/wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,lfac_)**2)

         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,bphi_)=w(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,bphi_)+qdt/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            1)*(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,bphi_)*wCT(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,mr_) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            br_)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            mphi_)) /xiplusB2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
       else
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mr_)=w(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,mr_)+qdt/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            1)*ptot(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
       end if
       if(srmhd_glm) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,br_)=w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,br_)+qdt*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          psi_)/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)
    case ('spherical')
       h1xmin1=ixOmin1-kr(1,1);h1xmin2=ixOmin2-kr(1,2)
       h1xmax1=ixOmax1-kr(1,1);h1xmax2=ixOmax2-kr(1,2)
       h2xmin1=ixOmin1-kr(2,1);h2xmin2=ixOmin2-kr(2,2)
       h2xmax1=ixOmax1-kr(2,1);h2xmax2=ixOmax2-kr(2,2);

       ! m1
       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=ptot(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          1) *(block%surfaceC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          1)-block%surfaceC(h1xmin1:h1xmax1,h1xmin2:h1xmax2,&
          1))/block%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
       if(ndir>1) then
         do idir=2,ndir
           tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2)+(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              mom(idir))**2.0-(VdotB2(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2)*B2(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2)))/(xiplusB2(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2))-wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              mag(idir))**2/wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,lfac_)**2.0
         end do
       end if
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(1))=w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,mom(1))+qdt*tmp(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)
       ! b1
       if(srmhd_glm) then
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(1))=w(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,mag(1))+qdt/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            1)*2.0d0*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,psi_)
       end if

       
       ! m2
 !tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
       ! This will make hydrostatic p=const an exact solution
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(2))=w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,mom(2))+qdt*ptot(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2) *(block%surfaceC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          2)-block%surfaceC(h2xmin1:h2xmax1,h2xmin2:h2xmax2,&
          2)) /block%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=-((wCT(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,mom(1))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          mom(2))-VdotB2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*wCT(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,mag(1))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          mag(2)))/xiplusB2(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          mag(1))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          mag(2))/wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,lfac_)**2.0)

       if(ndir==3) then
         tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2)+((wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            mom(3))**2-VdotB2(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            mag(3))**2.0)/xiplusB2(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2) -(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            mag(3))/wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            lfac_))**2.0)*dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            2))/dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2))
       end if
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(2))=w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,mom(2))+qdt*tmp(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)
       ! b2
       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(wCT(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,mom(1))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          mag(2)) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          mom(2))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          mag(1)))/xiplusB2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

       if(srmhd_glm) then
         tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2) + dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            2))/dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2))*wCT(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,psi_)
       end if
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(2))=w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,mag(2))+qdt*tmp(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)
      

       if(ndir==3) then
         ! m3
         if(.not.angmomfix) then
           tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=-((wCT(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2,mom(3))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              mom(1))-VdotB2(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              mag(3))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              mag(1)))/xiplusB2(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              mag(3))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              mag(1))/wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              lfac_)**2.0d0)  -((wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              mom(2))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              mom(3))-VdotB2(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              mag(2))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              mag(3)))/xiplusB2(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              mag(2))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              mag(3))/wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              lfac_)**2.0d0) *dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              2))/dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2))
           w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(3))=w(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2,mom(3))+qdt*tmp(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2)/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)
         else
           call mpistop("angmomfix not implemented yet in MHD")
         end if
         ! b3
         tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(wCT(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,mom(1))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            mag(3)) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            mom(3))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            mag(1)))/xiplusB2(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2)  -(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            mom(3))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            mag(2)) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            mom(2))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            mag(3)))*dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            2)) /(xiplusB2(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2)*dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)))
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(3))=w(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,mag(3))+qdt*tmp(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2)/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)
       end if
    end select
  end subroutine srmhd_add_source_geom

  !> Compute 2 times total magnetic energy
  function srmhd_mag_en_all(w, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2) result(mge)
    use mod_global_parameters
    integer, intent(in)           :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    real(kind=dp)   , intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw)
    real(kind=dp)                 :: mge(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    if (B0field) then
      mge = sum((w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          mag(:))+block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,:,block%iw0))**2,&
          dim=ndim+1)
    else
      mge = sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mag(:))**2, dim=ndim+1)
    end if
  end function srmhd_mag_en_all

  !> Compute full magnetic field by direction
  function srmhd_mag_i_all(w, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,idir) result(mgf)
    use mod_global_parameters
    integer, intent(in)           :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, idir
    real(kind=dp)   , intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw)
    real(kind=dp)                 :: mgf(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    if (B0field) then
      mgf = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          mag(idir))+block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir,block%iw0)
    else
      mgf = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mag(idir))
    end if
  end function srmhd_mag_i_all

  !> Compute evolving magnetic energy from primitive
  function srmhd_mag_en_primitive(w, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2) result(mge)
    use mod_global_parameters, only: nw, ndim
    integer, intent(in)           :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    real(kind=dp)   , intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw)
    real(kind=dp)                 :: mge(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    mge = sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mag(:))**2_dp, dim=ndim+1)
    mge = 0.5_dp*(mge(ixOmin1:ixOmax1,ixOmin2:ixOmax2)+ mge(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)*sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mom(:))**2_dp,&
        dim=ndim+1)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       lfac_)+sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mag(:))*w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2, mom(:))))/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,lfac_)
  end function srmhd_mag_en_primitive

  !> compute kinetic energy from primitive
  function srmhd_kin_en_primitive(w, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, inv_rho) result(ke)
    ! made by Z. Meliani 10/02/2018
    use mod_global_parameters, only: nw, ndim
    integer, intent(in)           :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    real(kind=dp)   , intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw)
    real(kind=dp)                 :: ke(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    real(kind=dp)   , intent(in), optional :: inv_rho(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)

    if (present(inv_rho)) then
       ke = (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,xi_)-1.0d0) * inv_rho
    else
       ke = (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,xi_)-1.0d0) / w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2, rho_)
    end if
  end function srmhd_kin_en_primitive

!  !>get Hall speed
!  subroutine srmhd_getv_Hall(w,x,ixI^L,ixO^L,vHall)
!    use mod_global_parameters

!    integer, intent(in)             :: ixI^L, ixO^L
!    real(kind=dp)   , intent(in)    :: w(ixI^S,nw)
!    real(kind=dp)   , intent(in)    :: x(ixI^S,1:ndim)
!    real(kind=dp)   , intent(inout) :: vHall(ixI^S,1:3)

!    integer          :: idir, idirmin
!    real(kind=dp)    :: current(ixI^S,7-2*ndir:3)

!    ! Calculate current density and idirmin
!    call get_current(w,ixI^L,ixO^L,idirmin,current)
!    vHall(ixO^S,1:3) = zero
!    vHall(ixO^S,idirmin:3) = - srmhd_etah*current(ixO^S,idirmin:3)
!    do idir = idirmin, 3
!       vHall(ixO^S,idir) = vHall(ixO^S,idir)/w(ixO^S,rho_)
!    end do

!  end subroutine srmhd_getv_Hall

!  subroutine srmhd_getdt_Hall(w,x,ixI^L,ixO^L,dx^D,dthall)
!    use mod_global_parameters

!    integer, intent(in) :: ixI^L, ixO^L
!    real(kind=dp)   , intent(in)    :: dx^D
!    real(kind=dp)   , intent(in)    :: w(ixI^S,1:nw)
!    real(kind=dp)   , intent(in)    :: x(ixI^S,1:ndim)
!    real(kind=dp)   , intent(out)   :: dthall
!    !.. local ..
!    real(kind=dp)    :: dxarr(ndim)
!    real(kind=dp)    :: bmag(ixI^S)

!    dthall=bigdouble

!    ! because we have that in cmax now:
!    return

!    ^D&dxarr(^D)=dx^D;

!    if (.not. B0field) then
!       bmag(ixO^S)=sqrt(sum(w(ixO^S,mag(:))**2, dim=ndim+1))
!       bmag(ixO^S)=sqrt(sum((w(ixO^S,mag(:)) + block%B0(ixO^S,1:ndir,block%iw0))**2))
!    end if

!    if(slab) then
!      dthall=dtdiffpar*minval(dxarr(1:ndim))**2.0d0/(srmhd_etah*maxval(bmag(ixO^S)/w(ixO^S,rho_)))
!    else
!      dthall=dtdiffpar*minval(block%ds(ixO^S,1:ndim))**2.0d0/(srmhd_etah*maxval(bmag(ixO^S)/w(ixO^S,rho_)))
!    end if

!  end subroutine srmhd_getdt_Hall

  !> This implements eq. (42) in Dedner et al. 2002 JcP 175
  !> Gives the Riemann solution on the interface
  !> for the normal B component and Psi in the GLM-MHD system.
  !> 23/04/2013 Oliver Porth
  subroutine glmSolve(wLC,wRC,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,idir)
    use mod_global_parameters
    real(kind=dp)   , intent(inout) :: wLC(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw), wRC(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2, idir
    real(kind=dp)                   :: dB(ixImin1:ixImax1,ixImin2:ixImax2),&
        dPsi(ixImin1:ixImax1,ixImin2:ixImax2)

    dB(ixOmin1:ixOmax1,ixOmin2:ixOmax2)   = wRC(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,mag(idir)) - wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       mag(idir))
    dPsi(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = wRC(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,psi_) - wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,psi_)

    wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(idir))   = 0.5d0 * &
       (wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(idir)) + wLC(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,mag(idir))) - 0.5d0/cmax_global * dPsi(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)
    wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,psi_)       = 0.5d0 * &
       (wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,psi_) + wLC(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,psi_)) - 0.5d0*cmax_global * dB(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)

    wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(idir)) = wLC(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,mag(idir))
    wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,psi_) = wLC(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,psi_)

  end subroutine glmSolve

  subroutine srmhd_boundary_adjust
    use mod_global_parameters
    integer :: iB, idim, iside, iigrid, igrid
    integer :: ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixOmin1,ixOmin2,ixOmax1,&
       ixOmax2, i1,i2

    ixGmin1=ixGlo1;ixGmin2=ixGlo2;ixGmax1=ixGhi1;ixGmax2=ixGhi2;
     do iigrid=1,igridstail; igrid=igrids(iigrid);
        if(.not.phyboundblock(igrid)) cycle
        saveigrid=igrid
        block=>pw(igrid)
        dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid);
        do idim=1,ndim
           ! to avoid using as yet unknown corner info in more than 1D, we
           ! fill only interior mesh ranges of the ghost cell ranges at first,
           ! and progressively enlarge the ranges to include corners later
           do iside=1,2
              i1=kr(1,idim)*(2*iside-3);i2=kr(2,idim)*(2*iside-3);
              if (neighbor_type(i1,i2,igrid)/=1) cycle
              iB=(idim-1)*2+iside
              if(.not.boundary_divbfix(iB)) cycle
              if(any(typeboundary(:,iB)=="special")) then
                ! MF nonlinear force-free B field extrapolation and data driven
                ! require normal B of the first ghost cell layer to be untouched by
                ! fixdivB=0 process, set boundary_divbfix_skip(iB)=1 in par file
                select case (idim)
                case (1)
                   if (iside==2) then
                      ! maximal boundary
                      ixOmin1=ixGmax1+1-nghostcells+&
                         boundary_divbfix_skip(2*1);ixOmin2=ixGmin2;
                      ixOmax1=ixGmax1;ixOmax2=ixGmax2;
                   else
                      ! minimal boundary
                      ixOmin1=ixGmin1;ixOmin2=ixGmin2;
                      ixOmax1=ixGmin1-1+nghostcells-boundary_divbfix_skip(2*1-&
                         1);ixOmax2=ixGmax2;
                   end if 
                case (2)
                   if (iside==2) then
                      ! maximal boundary
                      ixOmin1=ixGmin1
                      ixOmin2=ixGmax2+1-nghostcells+&
                         boundary_divbfix_skip(2*2);
                      ixOmax1=ixGmax1;ixOmax2=ixGmax2;
                   else
                      ! minimal boundary
                      ixOmin1=ixGmin1;ixOmin2=ixGmin2;
                      ixOmax1=ixGmax1
                      ixOmax2=ixGmin2-1+nghostcells-boundary_divbfix_skip(2*2-&
                         1);
                   end if 
                end select
                call fixdivB_boundary(ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixOmin1,&
                   ixOmin2,ixOmax1,ixOmax2,pw(igrid)%wb,pw(igrid)%x,iB)
              end if
           end do
        end do
     end do

  end subroutine srmhd_boundary_adjust

  subroutine fixdivB_boundary(ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,w,x,iB)
    use mod_global_parameters

    integer, intent(in) :: ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,iB
    real(kind=dp)   , intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:nw)
    real(kind=dp)   , intent(in) :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:ndim)

    real(kind=dp)    :: dx1x2,dx1x3,dx2x1,dx2x3,dx3x1,dx3x2
    integer :: ix1,ix2,ixFmin1,ixFmin2,ixFmax1,ixFmax2

    select case(iB)
     case(1)
       ! 2nd order CD for divB=0 to set normal B component better
       if(srmhd_energy.and..not.block%e_is_internal) call &
          srmhd_to_primitive(ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixOmin1,ixOmin2,&
          ixOmax1,ixOmax2,w,x)
       
       ixFmin1=ixOmin1+1
       ixFmax1=ixOmax1+1
       ixFmin2=ixOmin2+1
       ixFmax2=ixOmax2-1
       if(slab) then
         dx1x2=dxlevel(1)/dxlevel(2)
         do ix1=ixFmax1,ixFmin1,-1
           w(ix1-1,ixFmin2:ixFmax2,mag(1))=w(ix1+1,ixFmin2:ixFmax2,&
              mag(1)) +dx1x2*(w(ix1,ixFmin2+1:ixFmax2+1,mag(2))-w(ix1,&
              ixFmin2-1:ixFmax2-1,mag(2)))
         enddo
       else
         do ix1=ixFmax1,ixFmin1,-1
           w(ix1-1,ixFmin2:ixFmax2,mag(1))=( (w(ix1+1,ixFmin2:ixFmax2,&
              mag(1))+w(ix1,ixFmin2:ixFmax2,mag(1)))*block%surfaceC(ix1,&
              ixFmin2:ixFmax2,1)+(w(ix1,ixFmin2+1:ixFmax2+1,mag(2))+w(ix1,&
              ixFmin2:ixFmax2,mag(2)))*block%surfaceC(ix1,ixFmin2:ixFmax2,&
              2)-(w(ix1,ixFmin2:ixFmax2,mag(2))+w(ix1,ixFmin2-1:ixFmax2-1,&
              mag(2)))*block%surfaceC(ix1,ixFmin2-1:ixFmax2-1,&
              2) )/block%surfaceC(ix1-1,ixFmin2:ixFmax2,1)-w(ix1,&
              ixFmin2:ixFmax2,mag(1))
         end do
       end if
      
       
       if(srmhd_energy.and..not.block%e_is_internal) call &
          srmhd_to_conserved(ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixOmin1,ixOmin2,&
          ixOmax1,ixOmax2,w,x)
     case(2)
       if(srmhd_energy.and..not.block%e_is_internal) call &
          srmhd_to_primitive(ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixOmin1,ixOmin2,&
          ixOmax1,ixOmax2,w,x)
       
       ixFmin1=ixOmin1-1
       ixFmax1=ixOmax1-1
       ixFmin2=ixOmin2+1
       ixFmax2=ixOmax2-1
       if(slab) then
         dx1x2=dxlevel(1)/dxlevel(2)
         do ix1=ixFmin1,ixFmax1
           w(ix1+1,ixFmin2:ixFmax2,mag(1))=w(ix1-1,ixFmin2:ixFmax2,&
              mag(1)) -dx1x2*(w(ix1,ixFmin2+1:ixFmax2+1,mag(2))-w(ix1,&
              ixFmin2-1:ixFmax2-1,mag(2)))
         enddo
       else
         do ix1=ixFmin1,ixFmax1
           w(ix1+1,ixFmin2:ixFmax2,mag(1))=( (w(ix1-1,ixFmin2:ixFmax2,&
              mag(1))+w(ix1,ixFmin2:ixFmax2,mag(1)))*block%surfaceC(ix1-1,&
              ixFmin2:ixFmax2,1)-(w(ix1,ixFmin2+1:ixFmax2+1,mag(2))+w(ix1,&
              ixFmin2:ixFmax2,mag(2)))*block%surfaceC(ix1,ixFmin2:ixFmax2,&
              2)+(w(ix1,ixFmin2:ixFmax2,mag(2))+w(ix1,ixFmin2-1:ixFmax2-1,&
              mag(2)))*block%surfaceC(ix1,ixFmin2-1:ixFmax2-1,&
              2) )/block%surfaceC(ix1,ixFmin2:ixFmax2,1)-w(ix1,ixFmin2:ixFmax2,&
              mag(1))
         end do
       end if
      
       
       if(srmhd_energy.and..not.block%e_is_internal) call &
          srmhd_to_conserved(ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixOmin1,ixOmin2,&
          ixOmax1,ixOmax2,w,x)
     case(3)
       if(srmhd_energy.and..not.block%e_is_internal) call &
          srmhd_to_primitive(ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixOmin1,ixOmin2,&
          ixOmax1,ixOmax2,w,x)
       
       ixFmin1=ixOmin1+1
       ixFmax1=ixOmax1-1
       ixFmin2=ixOmin2+1
       ixFmax2=ixOmax2+1
       if(slab) then
         dx2x1=dxlevel(2)/dxlevel(1)
         do ix2=ixFmax2,ixFmin2,-1
           w(ixFmin1:ixFmax1,ix2-1,mag(2))=w(ixFmin1:ixFmax1,ix2+1,&
              mag(2)) +dx2x1*(w(ixFmin1+1:ixFmax1+1,ix2,&
              mag(1))-w(ixFmin1-1:ixFmax1-1,ix2,mag(1)))
         enddo
       else
         do ix2=ixFmax2,ixFmin2,-1
           w(ixFmin1:ixFmax1,ix2-1,mag(2))=( (w(ixFmin1:ixFmax1,ix2+1,&
              mag(2))+w(ixFmin1:ixFmax1,ix2,&
              mag(2)))*block%surfaceC(ixFmin1:ixFmax1,ix2,&
              2)+(w(ixFmin1+1:ixFmax1+1,ix2,mag(1))+w(ixFmin1:ixFmax1,ix2,&
              mag(1)))*block%surfaceC(ixFmin1:ixFmax1,ix2,&
              1)-(w(ixFmin1:ixFmax1,ix2,mag(1))+w(ixFmin1-1:ixFmax1-1,ix2,&
              mag(1)))*block%surfaceC(ixFmin1-1:ixFmax1-1,ix2,&
              1) )/block%surfaceC(ixFmin1:ixFmax1,ix2-1,2)-w(ixFmin1:ixFmax1,&
              ix2,mag(2))
         end do
       end if
      
       
       if(srmhd_energy.and..not.block%e_is_internal) call &
          srmhd_to_conserved(ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixOmin1,ixOmin2,&
          ixOmax1,ixOmax2,w,x)
     case(4)
       if(srmhd_energy.and..not.block%e_is_internal) call &
          srmhd_to_primitive(ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixOmin1,ixOmin2,&
          ixOmax1,ixOmax2,w,x)
       
       ixFmin1=ixOmin1+1
       ixFmax1=ixOmax1-1
       ixFmin2=ixOmin2-1
       ixFmax2=ixOmax2-1
       if(slab) then
         dx2x1=dxlevel(2)/dxlevel(1)
         do ix2=ixFmin2,ixFmax2
           w(ixFmin1:ixFmax1,ix2+1,mag(2))=w(ixFmin1:ixFmax1,ix2-1,&
              mag(2)) -dx2x1*(w(ixFmin1+1:ixFmax1+1,ix2,&
              mag(1))-w(ixFmin1-1:ixFmax1-1,ix2,mag(1)))
         end do
       else
         do ix2=ixFmin2,ixFmax2
           w(ixFmin1:ixFmax1,ix2+1,mag(2))=( (w(ixFmin1:ixFmax1,ix2-1,&
              mag(2))+w(ixFmin1:ixFmax1,ix2,&
              mag(2)))*block%surfaceC(ixFmin1:ixFmax1,ix2-1,&
              2)-(w(ixFmin1+1:ixFmax1+1,ix2,mag(1))+w(ixFmin1:ixFmax1,ix2,&
              mag(1)))*block%surfaceC(ixFmin1:ixFmax1,ix2,&
              1)+(w(ixFmin1:ixFmax1,ix2,mag(1))+w(ixFmin1-1:ixFmax1-1,ix2,&
              mag(1)))*block%surfaceC(ixFmin1-1:ixFmax1-1,ix2,&
              1) )/block%surfaceC(ixFmin1:ixFmax1,ix2,2)-w(ixFmin1:ixFmax1,ix2,&
              mag(2))
         end do
       end if
      
       
       if(srmhd_energy.and..not.block%e_is_internal) call &
          srmhd_to_conserved(ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixOmin1,ixOmin2,&
          ixOmax1,ixOmax2,w,x)
     
     case default
       call mpistop("Special boundary is not defined for this region")
    end select

   end subroutine fixdivB_boundary

   !> subtroutine srmhd_get_auxiliary calcule using srmhd_con2prim to calculate the enthalpy and the lorentz factor
   subroutine srmhd_get_auxiliary(clipping,w,x,ixImin1,ixImin2,ixImax1,ixImax2,&
      ixOmin1,ixOmin2,ixOmax1,ixOmax2,subname,use_oldaux,sqrB,SdotB)
    ! made by Z. Meliani 10/02/2018
    use mod_global_parameters
    use mod_srmhd_con2prim
    implicit none

    logical, intent(in)                    :: clipping
    integer, intent(in)                    :: ixImin1,ixImin2,ixImax1,ixImax2,&
       ixOmin1,ixOmin2,ixOmax1,ixOmax2
    real(kind=dp)   , intent(in)           :: x(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:ndim)
    real(kind=dp)   , intent(inout)        :: w(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:nw)
    character(len=*), intent(in)           :: subname
    logical         , intent(in), target,&
        optional :: use_oldaux(ixImin1:ixImax1,ixImin2:ixImax2)
    real(kind=dp)   , intent(in), target, optional :: sqrB(ixImin1:ixImax1,&
       ixImin2:ixImax2),SdotB(ixImin1:ixImax1,ixImin2:ixImax2)
    ! .. local
    integer                               :: ix1,ix2,ierror,idir
    integer                               :: flag_error(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)
    character(len=len(subname)+30)        :: subname_loc
    logical, pointer                      :: use_oldaux_in(:,:)
    real(kind=dp), pointer                :: sqrB_in(:,:),SdotB_in(:,:)
    real(kind=dp),target                  :: SdotB_loc(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2),sqrB_loc(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    real(kind=dp), dimension(ixOmin1:ixOmax1,ixOmin2:ixOmax2)       :: sqrS,&
       SdotB2,E2
    real(kind=dp), dimension(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:ndir):: tmp,v
    !---------------------------------------------------------------------



where(dabs(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
   mag(3)))<smalldouble)w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(3))=0.0_dp
where(dabs(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
   mag(1)))<smalldouble)w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(1))=0.0_dp
where(dabs(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
   mag(2)))<smalldouble)w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(2))=0.0_dp
    if(present(use_oldaux))then
      use_oldaux_in=>use_oldaux
    else
      use_oldaux_in=>grid_logical_false
    end if

    if(present(sqrB))then
      sqrB_in=>sqrB
    else
      sqrB_loc=SUM(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(1):mag(ndir))**2.0_dp,&
         dim=ndim+1)
      sqrB_in=>sqrB_loc
    end if

    if(present(SdotB))then
      SdotB_in=>SdotB
    else
      SdotB_loc=SUM(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         mom(1):mom(ndir))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(1):mag(ndir)),&
         dim=ndim+1)
      SdotB_in=>SdotB_loc
    end if

    sqrS(ixOmin1:ixOmax1,ixOmin2:ixOmax2)  = sum(w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,mom(1):mom(ndir))**2.0_dp,dim=ndim+1)
    SdotB2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)= SdotB_in(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)**2.0_dp
    E2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)    = sqrB_in(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)*sqrS(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)-SdotB2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)




   do ix2=ixOmin2,ixOmax2
   do ix1=ixOmin1,ixOmax1
   !PRINT*,' is wrong ',ix^D,w(ix^D,:)
    call srmhd_con2prim(w(ix1,ix2,d_),w(ix1,ix2,e_),use_oldaux_in(ix1,ix2),&
       sqrB_in(ix1,ix2),sqrS(ix1,ix2),SdotB_in(ix1,ix2),SdotB2(ix1,ix2),E2(ix1,&
       ix2),w(ix1,ix2,lfac_),w(ix1,ix2,xi_),ierror)

    if(check_small_values)then
     if(ierror/=0) then
       flag_error(ix1,ix2) = ierror
     else
       flag_error(ix1,ix2) = 0
     end if

    end if

   enddo
   enddo

    ! compute v*(xi+B^2)
    cond_isflow : if(any(sqrS(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2) /= zero.and.flag_error(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)/=0))then
        where(sqrS(ixOmin1:ixOmax1,ixOmin2:ixOmax2)/= &
           zero.and.flag_error(ixOmin1:ixOmax1,ixOmin2:ixOmax2)/=0)
          E2(ixOmin1:ixOmax1,ixOmin2:ixOmax2) =sdotb_in(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,xi_)
        end where
        Loop_idir0 : do idir=1,ndir
         where(sqrS(ixOmin1:ixOmax1,ixOmin2:ixOmax2)/= &
            zero.and.flag_error(ixOmin1:ixOmax1,ixOmin2:ixOmax2)/=0)
           tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir)=w(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2,mom(idir)) + w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              mag(idir))*E2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
         end where
        end do Loop_idir0
        where(sqrS(ixOmin1:ixOmax1,ixOmin2:ixOmax2)/= &
           zero.and.flag_error(ixOmin1:ixOmax1,ixOmin2:ixOmax2)/=0)
        ! reuse  E2,SdotB2
         E2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)     = &
            dsqrt(sum(tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:ndir)**2.0_dp,&
            dim=ndim+1))
         SdotB2(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = &
            dsqrt(1.0_dp-1.0_dp/(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,lfac_)**2))
         E2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)     = SdotB2(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2)/E2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
        endwhere
        Loop_idir : do idir=1,ndir
         where(sqrS(ixOmin1:ixOmax1,ixOmin2:ixOmax2)/= &
            zero.and.flag_error(ixOmin1:ixOmax1,ixOmin2:ixOmax2)/=0)
           v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir)=E2(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2)*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir)
         end where
        end do  Loop_idir
        where(sqrS(ixOmin1:ixOmax1,ixOmin2:ixOmax2)/= &
           zero.and.flag_error(ixOmin1:ixOmax1,ixOmin2:ixOmax2)/=0)
         E2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=sum(v(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,1:ndir)**2.0_dp,dim=ndim+1)
        end where
        where(sqrS(ixOmin1:ixOmax1,ixOmin2:ixOmax2)/= &
           zero.and.flag_error(ixOmin1:ixOmax1,ixOmin2:ixOmax2)/=0)
         where(E2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)>one)flag_error(&
            ixOmin1:ixOmax1,ixOmin2:ixOmax2)=6
        end where
    end if cond_isflow

   is_check_small : if(check_small_values)then

    is_flag_on : if(any(flag_error(ixOmin1:ixOmax1,ixOmin2:ixOmax2)/=0))then
     subname_loc='srmhd_get_auxiliary from -> '//trim(subname)
     call srmhd_handle_small_values(.false., w, x, ixImin1,ixImin2,ixImax1,&
        ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2,subname_loc,&
        flag_error=flag_error)
    end if is_flag_on
   end if is_check_small
   nullify(SdotB_in,sqrB_in,use_oldaux_in)


 end subroutine srmhd_get_auxiliary


! !> subtroutine srmhd_get_auxiliary calcule using srmhd_con2prim to calculate the enthalpy and the lorentz factor
! subroutine srmhd_get_auxiliary_loc(clipping,w,x,ixI^L,ixO^L,subname,use_oldaux,&
!                                   sqrB,SdotB)
!  ! made by Z. Meliani 10/02/2018
!  use mod_global_parameters
!  use mod_srmhd_con2prim
!  implicit none
!
!  logical, intent(in)                    :: clipping
!  integer, intent(in)                    :: ixI^L,ixO^L
!  real(kind=dp)   , intent(in)           :: x(ixI^S,1:ndim)
!  real(kind=dp)   , intent(inout)        :: w(ixI^S,1:nw)
!  character(len=*), intent(in)           :: subname
!   logical         , intent(in), target, optional :: use_oldaux(ixI^S)
!   real(kind=dp)   , intent(in), target, optional :: sqrB(ixI^S),SdotB(ixI^S) ! .. local
!  integer                         :: ix^D,ierror
!  integer                         :: flag_error(ixO^S)
!  character(len=len(subname)+30)  :: subname_loc
!  logical, pointer                :: use_oldaux_in(:^D&)
!  real(kind=dp), pointer          :: sqrB_in(:^D&),SdotB_in(:^D&)
!  !---------------------------------------------------------------------
!     if(present(use_oldaux))then
!       use_oldaux_in=>use_oldaux
!     else
!       use_oldaux_in=>srmhd_logical_false
!     end if
!
!
! {do ix^DB=ixOmin^DB,ixOmax^DB\}
!  call srmhd_con2prim(w(ix^D,d_),w(ix^D,mom(1):mom(ndir)),w(ix^D,e_),&
!           w(ix^D,mag(1):mag(ndir)),w(ix^D,lfac_),w(ix^D,xi_),ierror)
!  if(check_small_values)then
!   if(ierror/=0) then
!     flag_error(ix^D) = ierror
!   else
!     flag_error(ix^D) = 0
!   end if
!
!  end if
! {enddo^D&\}
! is_check_small : if(check_small_values)then
!  is_flag_on : if(any(flag_error(ixO^S)/=0))then
!   subname_loc='srmhd_get_auxiliary from -> '//trim(subname)
!   call srmhd_handle_small_values(.false., &
!                                  w, x, &
!                                  ixI^L, ixO^L,subname_loc,&
!                                  flag_error=flag_error)
!  end if is_flag_on
! end if is_check_small
! end subroutine srmhd_get_auxiliary_loc

   subroutine srmhd_get_4u_from_3v(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
      ixOmin2,ixOmax1,ixOmax2,vtou,lfac)
    ! made by Z. Meliani 13/02/2018
    use mod_global_parameters
    implicit none
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    real(kind=dp)   , intent(inout) :: vtou(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndir)
    real(kind=dp)   , intent(inout) :: lfac(ixImin1:ixImax1,ixImin2:ixImax2)
    !.. local ..
    integer                         :: idir
    !------------------------------------------------------

    lfac(ixOmin1:ixOmax1,ixOmin2:ixOmax2)= &
       1.0d0/dsqrt(1.0d0-sum(vtou(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:ndir)**2.0,&
       dim=ndim+1))
    Loop_idir: do idir=1,ndir
     vtou(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir)=lfac(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2)*vtou(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir)
    end do Loop_idir
   end subroutine srmhd_get_4u_from_3v

    ! -------------------------------------------------------
    !> subroutine to compute 3 vector from 4 vector speed
   subroutine srmhd_get_3v_from_4u(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
      ixOmin2,ixOmax1,ixOmax2,lfac,vtou)
    ! made by Z. Meliani 13/02/2018
    use mod_global_parameters
    implicit none
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    real(kind=dp)   , intent(in)    :: lfac(ixImin1:ixImax1,ixImin2:ixImax2)
    real(kind=dp)   , intent(inout) :: vtou(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndir)

    !.. local ..
    integer                         :: idir
    !------------------------------------------------------

    Loop_idir: do idir=1,ndir
     vtou(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir)=vtou(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2,idir)/lfac(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    end do Loop_idir
  end subroutine srmhd_get_3v_from_4u
end module mod_srmhd_phys
