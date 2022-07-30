module mod_usr

  use mod_dust
  use mod_physics
  use mod_global_parameters
  use mod_obj_global_parameters
  use mod_obj_mat
  use mod_obj_cloud
  use mod_obj_ism
  use mod_obj_sn_remnant
  use mod_obj_usr_unit
  use mod_obj_cla_jet
  use mod_grackle_parameters
  use mod_grackle_chemistry
  implicit none
  save
  real(dp) :: theta, kx, ly, vc

  type usr_config
    logical           :: physunit_on
    logical           :: ism_on
    logical           :: profile_use_hse
    logical           :: cloud_on
    logical           :: jet_yso_on
    logical           :: grackle_chemistry_on
    logical           :: ism_list_diff
    logical           :: cloud_list_diff
    logical           :: reset_medium
    integer           :: cloud_number
    integer           :: ind_jet_associate_ism
    integer           :: ind_cloud_associate_ism
    integer           :: ism_number
    integer           :: jet_yso_number
    integer           :: cloud_structure
    character(len=30) :: coordinate_system
    character(len=30) :: filename
    logical           :: cloud_profile_on
    logical           :: cloud_profile_density_on
    logical           :: cloud_profile_pressure_on
    logical           :: cloud_profile_velocity_on

    logical           :: reset_flux_scheme_on
    logical           :: reset_limiter_on
    character(len=30) :: reset_flux_scheme_diffuse
    character(len=30) :: reset_flux_scheme_old_method
    real(kind=dp)     :: reset_flux_scheme_thresholdL1_max
    real(kind=dp)     :: reset_flux_scheme_thresholdL1_min

    real(kind=dp)     :: density_dusttogas_maxlimit
    real(kind=dp)     :: density_dusttogas_minlimit
    real(kind=dp)     :: temperature_max

    logical           :: phys_isotherm_on
    real(kind=dp)     :: phys_adiab
    real(kind=dp)     :: phys_temperature_isotherm
    character(len=20) :: phys_inuse


  end type usr_config
  type(usr_config)    :: usrconfig
  integer, parameter  :: n_dust_max = 20
  real(dp) :: SUM_MASS   = 0.0_dp
  real(dp) :: SUM_VOLUME = 0.0_dp

  ! add type for fluxes here
  ! Objects
  type (ISM),allocatable,target      :: ism_surround(:)
  type (cloud),allocatable,target    :: cloud_medium(:)
  type (cla_jet),allocatable,target  :: jet_yso(:)
  type(gr_objects)                   :: grackle_object
  type(grackle_type)                 :: grackle_structure
  !Default objects
  type (ISM),target                  :: ism_default
  type (cloud),target                :: cloud_default
  type (cla_jet),target              :: jet_yso_default
  type (physconfig),target           :: pre_phys_config_default
  type (dust),target                 :: dust_mialy
  type (dust),allocatable,target     :: the_dust_inuse(:)
  type(gr_objects)                   :: grackle_default
  !type(star) :: star_ms
  !type(star) :: sun

  type(usrphysical_unit) :: usr_physunit





contains
  subroutine usr_init
    ! order of call : 1,
    ! amrvac.t:read_arguments() ->  amrvac.t:usr_init() (THIS PROCEDURE)
    use mod_hd, only : hd_activate
    use mod_mhd, only : mhd_activate
    ! .. local ..
    integer :: i_cloud,i_ism
    !-------------------------------------------
    ! configuration of procedures to be used in this project
    usr_set_parameters  => initglobaldata_usr
    ! order of call : 2,
    ! amrvac.t:read_arguments() ->  amrvac.t:usr_init()
    ! -> amrvac.t:initialize_amrvac():
    ! phys_check() --> read_par_files() (in mod_initialize.t:initialize_amrvac)
    ! (read parameters from filelist,savelist,stoplist,methodlist,boundlist,meshlist,paramlist,emissionlist)
    ! --> initialize_vars() (Initialize (and allocate) simulation and grid variables)
    ! --> bc_data_init()+read_data_init() (Possibly load boundary condition data or initial data)
    ! --> usr_set_parameters (THIS PROCEDURE)
    usr_init_one_grid   => initonegrid_usr !> order of call : 3
    usr_special_bc      => specialbound_usr
    usr_aux_output      => specialvar_output
    usr_add_aux_names   => specialvarnames_output
    usr_source          => specialsource_usr
    usr_refine_grid     => specialrefine_usr
    usr_special_global  => usr_global_var
    usr_process_grid    => process_grid_usr
    usr_get_dt          => special_get_dt
    usr_internal_bc     => usr_special_internal_bc
    usr_reset_solver    => special_reset_solver

    usr_gravity_potential => special_pointmass_gravity_potential
    usr_gravity_fpotential => special_pointmass_gravity_fpotential
    call usr_set_default_parameters ! >mod_obj_usr_yso_jet.t



    call usr_physunit%set_default ! >mod_obj_usr_unit.t



    ! set default values for ISMs configuration
    ! this also set ism_default%myboundaries%myconfig with their default values
    call ism_default%set_default ! >mod_obj_ism.t


    ! set default values for clouds configuration
    ! this also set cloud_default%myboundaries%myconfig with their default values
    call cloud_default%set_default !>mod_obj_cloud.t

    ! set default values for jets configuration
    ! since jet_yso has no bc, no need here to
    ! also set jet_yso_default%myboundaries%myconfig which does not exist
    call jet_yso_default%set_default !>mod_obj_cla_jet.t

    ! set default values for grackle configuration
    call grackle_default%set_default_config !>mod_grackle_chemistry

    ! read .par parameters for ISMs, then, clouds, then jets
    call usr_params_read(par_files) ! >mod_obj_usr_yso_jet.t
    ! WHAT usr_params_read DOES:
    !0) Read user-defined user configuration usrconfig from &usr_list in .par
    !1) Read user-defined physical units in .par
    !2) If ism on, set all usr ISMs parameters to their default values
    ! and then read their user-defined parameters in .par
    ! and sets their boundary conditions types with the user-defined ism_config
    !3) If cloud on, set all usr clouds parameters to their default values
    ! and then read their user-defined parameters in .par
    !4) If jet on, set all usr jets parameters to their default values
    ! and then read their user-defined parameters in .par
    ! TO DO: 5) If computeflux is on, set all usr fluxes parameters to their default values
    ! and then read their user-defined parameters in .par

    if(usrconfig%profile_use_hse)then
      usr_gravity         => special_effective_gravity!special_pointmass_gravity
    else
      usr_gravity         => special_pointmass_gravity!special_effective_gravity
    end if

    ! complet all physical unit in use
    if(usrconfig%physunit_on) then
      call usr_physunit%set_complet(trim(usrconfig%phys_inuse))
    end if
    call usr_physical_unit

    ! set the code global coordinates system in use
    call set_coordinate_system(trim(usrconfig%coordinate_system))
    select case(trim(usrconfig%phys_inuse))
    case('hd')
     call hd_activate
    case('mhd')
     call mhd_activate
    end select

    call usr_check_conflict


  end subroutine usr_init
  !------------------------------------------------------------------

  subroutine special_effective_gravity(ixI^L,ixO^L,wCT,x,gravity_field)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    real(dp), intent(in)    :: x(ixI^S,1:ndim)
    real(dp), intent(in)    :: wCT(ixI^S,1:nw)
    real(dp), intent(out)   :: gravity_field(ixI^S,ndim)
    real(dp)                :: gravity_field_sphrc(ixI^S,ndim)
    real(dp)                        :: Ggrav, Mpoint, alpha,sqrtthree,threesqrtthree, halfdrr
    real(kind=dp), dimension(ixI^S) :: r_distance,radius
    real(dp)                        :: xd(ixI^S,1:ndim)
    real(dp)                        :: wprim(ixI^S,1:nw)
    real(kind=dp), dimension(ixI^S) :: rd_distance,rdradius
    real(kind=dp), dimension(ixI^S) :: theta_profile
    real(kind=dp), dimension(ixI^S) :: costta,sintta,rr
    real(kind=dp), dimension(ixI^S) :: fc,dfcdr,dfcdt
    logical                        :: xpatch(ixI^S)
    real(kind=dp), dimension(1:ndim) :: zero_dim,dx_local,dx_local_min
    integer                      :: i_ism,idim,level,ix^D
    !-----------------------------------

    wprim(ixO^S,1:nw)=wCT(ixO^S,1:nw)
    call phys_to_primitive(ixI^L,ixI^L,wprim,x)
    Ggrav = constusr%G

    gravity_field(ixO^S,1:ndim)=0.0_dp
    gravity_field_sphrc(ixO^S,1:ndim)=0.0_dp
    {zero_dim(^D)=0.0_dp\}!-dx(1,1)
    {^D&
    xd(ixO^S,^D)=x(ixO^S,^D)-kr(r_,^D)*ism_surround(i_ism)%myconfig%profile_rd
    }




    ! Here we set the graviationnal acceleration inside the whole domain
    ! from the state vector wCT
    ! here the called wCT contains conservative variables
      if(usrconfig%ism_on)then
        i_ism =0
        level = node(plevel_,saveigrid)
        ^D&dx_local(^D)=((xprobmax^D-xprobmin^D)/(domain_nx^D))/(2.0_dp**(level-1.0_dp));
        call usr_distance(ixI^L,ixO^L,typeaxial,&
                          zero_dim,x,r_distance)
        call usr_distance(ixI^L,ixO^L,typeaxial,&
                          zero_dim,xd,rd_distance)
        call usr_get_theta(ixI^L,ixO^L,x,theta_profile)
        Mpoint = ism_surround(i_ism)%myconfig%profile_Mstar
        radius(ixO^S)=r_distance(ixO^S)/ism_surround(i_ism)%myconfig%profile_rd
        rdradius(ixO^S)=rd_distance(ixO^S)/ism_surround(i_ism)%myconfig%profile_rd
        sqrtthree = DSQRT(3.0_dp)
        threesqrtthree = 3.0_dp * sqrtthree

        Loop_isms : do i_ism=0,usrconfig%ism_number-1
        select case(trim(ism_surround(i_ism)%myconfig%profile_density))
          case('Ulrich1976')

          ^D&dx_local_min(^D)=min(ism_surround(i_ism)%myconfig%escapencells(^D)*dx_local(^D),&
                          ism_surround(i_ism)%myconfig%escapencellsglobal(^D)*dx(^D,1));

          dx_local_min(:)=dx_local_min(:)/ism_surround(i_ism)%myconfig%profile_rd
          xpatch(ixI^S)=.true.


          halfdrr = 0.5_dp*DSQRT(SUM(dx_local_min**2.0_dp))

          rr(ixO^S)=radius(ixO^S)
          sintta(ixO^S)=x(ixO^S,r_)/DSQRT(x(ixO^S,r_)**2.0_dp+&
          x(ixO^S,z_)**2.0_dp)
          costta(ixO^S)=x(ixO^S,z_)/DSQRT(x(ixO^S,r_)**2.0_dp+&
                                 x(ixO^S,z_)**2.0_dp)



          !Case 2: r>1:
          where(rr(ixO^S)>1.0_dp)
            gravity_field_sphrc(ixO^S,1) = rr(ixO^S) ** ( 3.0_dp/2.0_dp ) * &
            ( ( 4.0_dp * ( rr(ixO^S) - 1.0_dp ) * DSINH( (1.0_dp/3.0_dp) * &
            DASINH( ( threesqrtthree * rr(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / &
            ( 2.0_dp * ( rr(ixO^S) - 1.0_dp ) ** (3.0_dp/2.0_dp) ) ) ) **2.0_dp - 1.0_dp ) / rr(ixO^S) &
            + 1.0_dp ) * DSQRT( ( sqrtthree * DCOS( theta_profile(ixO^S) ) * 1.0_dp / DSINH( (1.0_dp/3.0_dp) * &
            DASINH( ( threesqrtthree * rr(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / &
            ( 2.0_dp * ( rr(ixO^S) - 1.0_dp ) ** (3.0_dp/2.0_dp) ) ) ) ) / &
            ( 2.0_dp * DSQRT( rr(ixO^S) - 1.0_dp ) ) + 1.0_dp ) * &
            ( -3.0_dp / ( 2.0_dp * rr(ixO^S) ** (5.0_dp/2.0_dp) * &
            ( ( 4.0_dp * ( rr(ixO^S) - 1.0_dp ) * DSINH( (1.0_dp/3.0_dp) * &
            DASINH( ( threesqrtthree * rr(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / &
            ( 2.0_dp * ( rr(ixO^S) - 1.0_dp ) ** (3.0_dp/2.0_dp) ) ) ) ** 2.0_dp - &
            1.0_dp ) / rr(ixO^S) + 1.0_dp ) * DSQRT( ( sqrtthree * DCOS( theta_profile(ixO^S) ) * &
            ( 1.0_dp / DSINH( (1.0_dp/3.0_dp) * DASINH( ( threesqrtthree * rr(ixO^S) * &
            DCOS( theta_profile(ixO^S) ) ) / ( 2.0_dp * ( rr(ixO^S) - 1.0_dp ) ** (3.0_dp/2.0_dp) ) ) ) ) ) / &
            ( 2.0_dp * DSQRT( rr(ixO^S) - 1.0_dp ) ) + 1.0_dp ) ) - &
            ( ( ( 8.0_dp * ( rr(ixO^S) - 1.0_dp ) * ( ( threesqrtthree * DCOS( theta_profile(ixO^S) ) ) / &
            ( 2.0_dp * ( rr(ixO^S) - 1.0_dp ) ** (3.0_dp/2.0_dp) ) - ( 9.0_dp * sqrtthree * rr(ixO^S) * &
            DCOS( theta_profile(ixO^S) ) ) / ( 4.0_dp * ( rr(ixO^S) - 1.0_dp ) ** (5.0_dp/2.0_dp) ) ) * &
            DSINH( (1.0_dp/3.0_dp) * DASINH( ( threesqrtthree * rr(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / &
            ( 2 * ( rr(ixO^S) - 1.0_dp ) ** (3.0_dp/2.0_dp) ) ) ) * DCOSH( (1.0_dp/3.0_dp) * DASINH( ( threesqrtthree * &
            rr(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / ( 2.0_dp * ( rr(ixO^S) - 1.0_dp ) ** (3.0_dp/2.0_dp) ) ) ) ) / &
            ( 3.0_dp * DSQRT( ( 27.0_dp * rr(ixO^S) ** 2.0_dp * DCOS( theta_profile(ixO^S) ) ** 2.0_dp ) / &
            ( 4.0_dp * ( rr(ixO^S) - 1.0_dp ) ** 3.0_dp ) + 1.0_dp ) ) + &
            4.0_dp * DSINH( (1.0_dp/3.0_dp) * DASINH( ( threesqrtthree * rr(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / ( 2.0_dp * &
            ( rr(ixO^S) - 1.0_dp ) ** (3.0_dp/2.0_dp) ) ) ) ** 2.0_dp ) / rr(ixO^S) - &
            ( 4.0_dp * ( rr(ixO^S) - 1.0_dp ) * DSINH( (1.0_dp/3.0_dp) * &
            DASINH( ( threesqrtthree * rr(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / &
            ( 2.0_dp * ( rr(ixO^S) - 1.0_dp ) ** (3.0_dp/2.0_dp) ) ) ) **2.0_dp - 1.0_dp ) / &
            rr(ixO^S) ** 2.0_dp ) / ( rr(ixO^S) ** (3.0_dp/2.0_dp) * ( ( 4.0_dp * ( rr(ixO^S) - 1.0_dp ) * &
            DSINH( (1.0_dp/3.0_dp) * DASINH( ( threesqrtthree * rr(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / &
            ( 2.0_dp * ( rr(ixO^S) - 1.0_dp ) ** (3.0_dp/2.0_dp) ) ) ) ** 2.0_dp - 1.0_dp ) / rr(ixO^S) + 1.0_dp ) ** 2.0_dp * &
            DSQRT( ( sqrtthree * DCOS( theta_profile(ixO^S) ) * ( 1.0_dp / SINH( (1.0_dp/3.0_dp) * &
            DASINH( ( threesqrtthree * rr(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / &
            ( 2.0_dp * ( rr(ixO^S) - 1.0_dp ) ** (3.0_dp/2.0_dp) ) ) ) ) ) / &
            ( 2.0_dp * DSQRT( rr(ixO^S) - 1.0_dp ) ) + 1.0_dp ) ) - &
            (-( DCOS( theta_profile(ixO^S) ) * ( ( threesqrtthree * DCOS( theta_profile(ixO^S) ) ) / &
            ( 2.0_dp * ( rr(ixO^S) - 1.0_dp ) ** (3.0_dp/2.0_dp) ) - ( 9.0_dp * DSQRT(3.0_dp) * &
            rr(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / ( 4.0_dp * ( rr(ixO^S) - 1.0_dp ) ** (5.0_dp/2.0_dp) ) ) * &
            ( 1.0_dp / DTANH( (1.0_dp/3.0_dp) * DASINH( ( threesqrtthree * rr(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / &
            ( 2.0_dp * ( rr(ixO^S) - 1.0_dp ) ** (3.0_dp/2.0_dp) ) ) ) ) *&
            (1.0_dp / DSINH( (1.0_dp/3.0_dp) * DASINH( ( threesqrtthree * rr(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / &
            ( 2.0_dp * ( rr(ixO^S) - 1.0_dp ) ** (3.0_dp/2.0_dp) ) ) ) ) ) / &
            ( 2.0_dp * sqrtthree * DSQRT( rr(ixO^S) - 1.0_dp ) * DSQRT( ( 27.0_dp * rr(ixO^S) ** 2.0_dp * &
            DCOS( theta_profile(ixO^S) ) ** 2.0_dp ) / ( 4.0_dp * ( rr(ixO^S) - 1.0_dp ) ** 3.0_dp ) + 1.0_dp ) ) - &
            ( sqrtthree * DCOS ( theta_profile(ixO^S) ) * &
            ( 1.0_dp / DSINH( (1.0_dp/3.0_dp) * DASINH( ( threesqrtthree * rr(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / &
            ( 2.0_dp * ( rr(ixO^S) - 1.0_dp ) ** (3.0_dp/2.0_dp) ) ) ) ) ) / &
            ( 4.0_dp * ( rr(ixO^S) - 1.0_dp ) ** (3.0_dp/2.0_dp) ) ) / &
            ( 2.0_dp * rr(ixO^S) ** (3.0_dp/2.0_dp) * ( ( 4.0_dp * ( rr(ixO^S) - 1.0_dp ) * &
            DSINH( (1.0_dp/3.0_dp) * DASINH( ( threesqrtthree * rr(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / &
            ( 2.0_dp * ( rr(ixO^S) - 1.0_dp ) ** (3.0_dp/2.0_dp) ) ) ) ** 2.0_dp - 1.0_dp) / rr(ixO^S) + 1.0_dp ) * &
            ( ( sqrtthree * DCOS( theta_profile(ixO^S) ) * &
            ( 1.0_dp / DSINH( (1.0_dp/3.0_dp) * DASINH( ( threesqrtthree * rr(ixO^S) * &
            DCOS ( theta_profile(ixO^S) ) ) / ( 2.0_dp * ( rr(ixO^S) - 1.0_dp ) ** (3.0_dp/2.0_dp) ) ) ) ) ) / &
            ( 2.0_dp * DSQRT( rr(ixO^S) - 1.0_dp ) ) + 1.0_dp ) ** (3.0_dp/2.0_dp) ) )

            gravity_field_sphrc(ixO^S,2) = rr(ixO^S) ** (1.0_dp/2.0_dp) * &
            ( ( 4.0_dp * ( rr(ixO^S) - 1.0_dp ) * DSINH( (1.0_dp/3.0_dp) * &
            DASINH( ( threesqrtthree * rr(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / &
            ( 2.0_dp * ( rr(ixO^S) - 1.0_dp ) ** (3.0_dp/2.0_dp) ) ) ) ** 2.0_dp - 1.0_dp ) / &
            rr(ixO^S) + 1.0_dp ) * DSQRT( ( sqrtthree * DCOS( theta_profile(ixO^S) ) * &
            ( 1.0_dp / DSINH( (1.0_dp/3.0_dp) * DASINH( ( threesqrtthree * rr(ixO^S) * &
            DCOS( theta_profile(ixO^S) ) ) / ( 2.0_dp * ( rr(ixO^S) - 1.0_dp ) ** (3.0_dp/2.0_dp) ) ) ) ) ) / &
            ( 2.0_dp * DSQRT( rr(ixO^S) - 1.0_dp ) ) + 1.0_dp ) * ( ( 4.0_dp * sqrtthree * DSIN( theta_profile(ixO^S) ) * &
            DSINH( (1.0_dp/3.0_dp) * DASINH( ( threesqrtthree * rr(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / &
            ( 2.0_dp * ( rr(ixO^S) - 1.0_dp ) ** (3.0_dp/2.0_dp) ) ) ) * DCOSH( (1.0_dp/3.0_dp) * &
            DASINH( ( threesqrtthree * rr(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / &
            ( 2.0_dp * ( rr(ixO^S) - 1.0_dp ) ** (3.0_dp/2.0_dp) ) ) ) ) / &
            ( DSQRT( rr(ixO^S) - 1.0_dp ) * rr(ixO^S) ** (3.0_dp/2.0_dp) * &
            DSQRT( ( 27.0_dp * rr(ixO^S) ** 2.0_dp * DCOS( theta_profile(ixO^S) ) ** 2.0_dp ) / &
            ( 4.0_dp * ( rr(ixO^S) - 1.0_dp ) ** 3.0_dp ) + 1.0_dp ) * &
            ( ( 4.0_dp * ( rr(ixO^S) - 1.0_dp ) * DSINH( (1.0_dp/3.0_dp) * DASINH( ( threesqrtthree * rr(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / &
            ( 2.0_dp * ( rr(ixO^S) - 1.0_dp ) ** (3.0_dp/2.0_dp) ) ) ) ** 2.0_dp - 1.0_dp ) / rr(ixO^S)  + 1.0_dp ) ** 2.0_dp * &
            DSQRT( ( sqrtthree * DCOS( theta_profile(ixO^S) ) * &
            ( 1.0_dp / DSINH( (1.0_dp/3.0_dp) * DASINH( ( threesqrtthree * rr(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / &
            ( 2.0_dp * ( rr(ixO^S) - 1.0_dp ) ** (3.0_dp/2.0_dp) ) ) ) ) ) / &
            ( 2.0_dp * DSQRT( rr(ixO^S) - 1.0_dp ) ) + 1.0_dp ) ) - ( ( 3.0_dp * rr(ixO^S) * &
            DSIN( theta_profile(ixO^S) ) * DCOS( theta_profile(ixO^S) ) * &
            ( 1.0_dp / DTANH( (1.0_dp/3.0_dp) * DASINH( ( threesqrtthree * rr(ixO^S) * &
            DCOS( theta_profile(ixO^S) ) ) / ( 2.0_dp * ( rr(ixO^S) - 1.0_dp ) ** (3.0_dp/2.0_dp) ) ) ) ) * &
            ( 1.0_dp / DSINH( (1.0_dp/3.0_dp) * DASINH( ( threesqrtthree * rr(ixO^S) * &
            DCOS( theta_profile(ixO^S) ) ) / ( 2.0_dp * ( rr(ixO^S) - 1.0_dp ) ** (3.0_dp/2.0_dp) ) ) ) ) ) / &
            ( 4.0_dp * ( rr(ixO^S) - 1.0_dp ) ** 2.0_dp * DSQRT( ( 27.0_dp * rr(ixO^S) ** 2.0_dp * DCOS( theta_profile(ixO^S) ) ** 2.0_dp ) / &
            ( 4.0_dp * ( rr(ixO^S) - 1.0_dp ) ** 3.0_dp ) + 1.0_dp ) ) - &
            ( sqrtthree * DSIN( theta_profile(ixO^S) ) * &
            ( 1.0_dp / DSINH( (1.0_dp/3.0_dp) * DASINH( ( threesqrtthree * rr(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / &
            ( 2.0_dp * ( rr(ixO^S) - 1.0_dp ) ** (3.0_dp/2.0_dp) ) ) ) ) ) / &
            ( 2.0_dp * DSQRT( rr(ixO^S) - 1.0_dp ) ) ) / ( 2.0_dp * rr(ixO^S) ** (3.0_dp/2.0_dp) * &
            ( ( 4.0_dp * ( rr(ixO^S) - 1.0_dp ) * DSINH( (1.0_dp/3.0_dp) * DASINH( ( threesqrtthree * rr(ixO^S) * &
            DCOS( theta_profile(ixO^S) ) ) / ( 2.0_dp * ( rr(ixO^S) - 1.0_dp ) ** (3.0_dp/2.0_dp) ) ) ) ** 2.0_dp - 1.0_dp ) / &
            rr(ixO^S) + 1.0_dp ) * ( ( sqrtthree * DCOS( theta_profile(ixO^S) ) * &
            ( 1.0_dp / DSINH( (1.0_dp/3.0_dp) * DASINH( ( threesqrtthree * rr(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / &
            ( 2.0_dp * ( rr(ixO^S) - 1.0_dp ) ** (3.0_dp/2.0_dp) ) ) ) ) ) / &
            ( 2.0_dp * DSQRT( rr(ixO^S) - 1.0_dp) ) + 1.0_dp ) ** (3.0_dp/2.0_dp) ) )

          end where



          !Case 3: r<1 and xi>0:
          where((rr(ixO^S)<1.0_dp).and.((((rr(ixO^S)*costta(ixO^S))/2.0_dp)**2.0_dp)-&
          (((1.0_dp-rr(ixO^S))/3.0_dp)**3.0_dp)>=0.0_dp))

            gravity_field_sphrc(ixO^S,1) = rr(ixO^S) ** (3.0_dp/2.0_dp) * &
            ( ( 4.0_dp * (1.0_dp - rr(ixO^S) ) * DCOSH( (1.0_dp/3.0_dp) * &
            DACOSH( ( threesqrtthree * rr(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / &
            ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) ) ) ** 2.0_dp - 1.0_dp ) / &
            rr(ixO^S) + 1.0_dp ) * DSQRT( ( sqrtthree * DCOS( theta_profile(ixO^S) ) * &
            ( 1.0_dp / DCOSH( (1.0_dp/3.0_dp) * DACOSH( ( threesqrtthree * rr(ixO^S) * &
            DCOS( theta_profile(ixO^S) ) ) / ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) ) ) ) ) / &
            ( 2.0_dp * DSQRT( 1.0_dp - rr(ixO^S) ) ) + 1.0_dp ) * &
            (- ( ( ( 8.0_dp * ( 1.0_dp - rr(ixO^S) ) * ( ( 9.0_dp * sqrtthree * rr(ixO^S) * &
            DCOS( theta_profile(ixO^S) ) ) / ( 4.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (5.0_dp/2.0_dp) ) + &
            ( threesqrtthree * DCOS( theta_profile(ixO^S) ) ) / &
            ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) ) * &
            DCOSH( (1.0_dp/3.0_dp) * DACOSH( ( threesqrtthree * rr(ixO^S) * &
            DCOS( theta_profile(ixO^S) ) ) / ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) ) ) * &
            DSINH( (1.0_dp/3.0_dp) * DACOSH( ( threesqrtthree * rr(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / &
            ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) ) ) ) / &
            ( 3.0_dp * DSQRT( ( threesqrtthree * rr(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / &
            ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) - 1.0_dp ) * DSQRT( ( threesqrtthree * &
            rr(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) + &
            1.0_dp ) ) - 4.0_dp * DCOSH( (1.0_dp/3.0_dp) * DACOSH( ( threesqrtthree * rr(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / &
            ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) ) ) ** 2.0_dp ) / rr(ixO^S) - &
            ( 4.0_dp * ( 1.0_dp - rr(ixO^S) ) * DCOSH( (1.0_dp/3.0_dp) * &
            DACOSH( ( threesqrtthree * rr(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / &
            ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) ) ) ** 2.0_dp - 1.0_dp ) / rr(ixO^S) ** 2.0_dp ) / &
            ( rr(ixO^S) ** (3.0_dp/2.0_dp) * ( ( 4.0_dp * ( 1.0_dp - rr(ixO^S) ) * &
            DCOSH( (1.0_dp/3.0_dp) * DACOSH( ( threesqrtthree * rr(ixO^S) * &
            DCOS( theta_profile(ixO^S) ) ) / ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) ) ) ** 2.0_dp - 1.0_dp ) / &
            rr(ixO^S) + 1.0_dp ) ** 2.0_dp * DSQRT( ( sqrtthree * DCOS( theta_profile(ixO^S) ) * &
            ( 1.0_dp / DCOSH( (1.0_dp/3.0_dp) * DACOSH( ( threesqrtthree * &
            rr(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) ) ) ) ) / &
            ( 2.0_dp * DSQRT( 1.0_dp - rr(ixO^S) ) ) + 1.0_dp ) ) - &
            ( ( sqrtthree * DCOS( theta_profile(ixO^S) ) * &
            ( 1.0_dp / DCOSH( (1.0_dp/3.0_dp) * DACOSH( ( threesqrtthree * &
            rr(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) ) ) ) ) / &
            ( 4.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) - &
            ( DCOS( theta_profile(ixO^S) ) * ( ( 9.0_dp * sqrtthree * rr(ixO^S) * &
            DCOS( theta_profile(ixO^S) ) ) / ( 4.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (5.0_dp/2.0_dp) ) + &
            ( threesqrtthree * DCOS( theta_profile(ixO^S) ) ) / &
            ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) ) * &
            ( 1.0_dp / DCOSH( (1.0_dp/3.0_dp) * DACOSH( ( threesqrtthree * &
            rr(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) ) ) ) * &
            DTANH( (1.0_dp/3.0_dp) * DACOSH( ( threesqrtthree * rr(ixO^S) * &
            DCOS( theta_profile(ixO^S) ) ) / ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) ) ) ) / &
            ( 2.0_dp * sqrtthree * DSQRT( 1.0_dp - rr(ixO^S) ) * &
            DSQRT( ( threesqrtthree * rr(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / &
            ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) - 1.0_dp ) * &
            DSQRT( ( threesqrtthree * rr(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / &
            ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) + 1.0_dp ) ) ) / &
            ( 2.0_dp * rr(ixO^S) ** (3.0_dp/2.0_dp) * ( ( 4.0_dp * ( 1.0_dp - rr(ixO^S) ) * &
            DCOSH( (1.0_dp/3.0_dp) * DACOSH( ( threesqrtthree * rr(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / &
            ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) ) ) ** 2.0_dp - 1.0_dp ) / rr(ixO^S) + 1.0_dp ) * &
            ( ( sqrtthree * DCOS( theta_profile(ixO^S) ) * &
            ( 1.0_dp / COSH( (1.0/3.0_dp * DACOSH( ( threesqrtthree * rr(ixO^S) * &
            DCOS( theta_profile(ixO^S) ) ) / ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) ) ) ) ) / &
            ( 2.0_dp * DSQRT( 1.0_dp - rr(ixO^S) ) ) + 1.0_dp ) ** (3.0_dp/2.0_dp) ) ) - &
            3.0_dp / ( 2.0_dp * rr(ixO^S) ** (5.0_dp/2.0_dp) * ( ( 4.0_dp * ( 1.0_dp - rr(ixO^S) ) * &
            cosh( (1.0_dp/3.0_dp) * DACOSH( ( threesqrtthree * rr(ixO^S) * &
            DCOS( theta_profile(ixO^S) ) ) / ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) ) ) ** 2.0_dp  - &
            1.0_dp ) / rr(ixO^S) + 1.0_dp ) * DSQRT( ( sqrtthree * DCOS( theta_profile(ixO^S) ) * &
            ( 1.0_dp / DCOSH( (1.0_dp/3.0_dp) * DACOSH( ( threesqrtthree * rr(ixO^S) * &
            DCOS( theta_profile(ixO^S) ) ) / ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) ) ) ) ) / &
            ( 2.0_dp * DSQRT( 1.0_dp - rr(ixO^S) ) ) + 1.0_dp ) ) )

            gravity_field_sphrc(ixO^S,2) = DSQRT( rr(ixO^S) ) * &
            ( ( 4.0_dp * ( 1.0_dp - rr(ixO^S) ) * DCOSH( (1.0_dp/3.0_dp) * DACOSH( ( threesqrtthree * rr(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / &
            ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) ) ) ** 2.0_dp - 1.0_dp ) / rr(ixO^S) + 1.0_dp ) * &
            DSQRT( ( sqrtthree * DCOS( theta_profile(ixO^S) ) * &
            ( 1.0_dp / DCOSH( (1.0_dp/3.0_dp) * DACOSH( ( threesqrtthree * rr(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / &
            ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) ) ) ) ) / &
            ( 2.0_dp * DSQRT( 1.0_dp - rr(ixO^S) ) ) + 1.0_dp ) * &
            ( ( 4.0_dp * sqrtthree * DSIN( theta_profile(ixO^S) ) * &
            DCOSH( (1.0_dp/3.0_dp) * DACOSH( ( threesqrtthree * rr(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / &
            ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) ) ) * &
            DSINH( (1.0_dp/3.0_dp) * DACOSH( ( threesqrtthree * rr(ixO^S) * &
            DCOS( theta_profile(ixO^S) ) ) / ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) ) ) ) / &
            ( DSQRT( 1.0_dp - rr(ixO^S) ) * rr(ixO^S) ** (3.0_dp/2.0_dp) * &
            DSQRT( ( threesqrtthree * rr(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / &
            ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp ) ) - 1.0_dp ) * &
            DSQRT( ( threesqrtthree * rr(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / &
            ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp ) ) + 1.0_dp ) * &
            ( ( 4.0_dp * ( 1.0_dp - rr(ixO^S) ) * DCOSH( (1.0_dp/3.0_dp) * DACOSH( ( threesqrtthree * &
            rr(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / &
            ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) ) ) ** 2.0_dp - 1.0_dp ) / rr(ixO^S) &
            + 1.0_dp ) ** 2.0_dp * DSQRT( ( sqrtthree * DCOS( theta_profile(ixO^S) ) * &
            ( 1.0_dp / COSH( (1.0_dp/3.0_dp) * DACOSH( ( threesqrtthree * rr(ixO^S) * &
            DCOS( theta_profile(ixO^S) ) ) / ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) ) ) ) ) / &
            ( 2.0_dp * DSQRT( 1.0_dp - rr(ixO^S) ) ) + 1.0_dp ) ) &
            - &
            ( ( 3.0_dp * rr(ixO^S) * DSIN( theta_profile(ixO^S) ) * DCOS( theta_profile(ixO^S) ) * &
            DTANH( (1.0_dp/3.0_dp) * DACOSH( ( threesqrtthree * rr(ixO^S) * &
            DCOS( theta_profile(ixO^S) ) ) / ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) ) ) * &
            ( 1.0_dp / DCOSH( (1.0_dp/3.0_dp) * DACOSH( ( threesqrtthree * rr(ixO^S) * &
            DCOS( theta_profile(ixO^S) ) ) / ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) ) ) ) ) / &
            ( 4.0_dp * ( 1.0_dp - rr(ixO^S) ) ** 2.0_dp * DSQRT( ( threesqrtthree * rr(ixO^S) * &
            DCOS( theta_profile(ixO^S) ) ) / ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) - &
            1.0_dp ) * DSQRT( ( threesqrtthree * rr(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / &
            ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) + 1.0_dp ) ) - &
            ( sqrtthree * DSIN( theta_profile(ixO^S) ) * &
            ( 1.0_dp / DCOSH( (1.0_dp/3.0_dp) * DACOSH( ( threesqrtthree * rr(ixO^S) * &
            DCOS( theta_profile(ixO^S) ) ) / ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) ) ) ) ) / &
            ( 2.0_dp * DSQRT( 1.0_dp - rr(ixO^S) ) ) )  / &
            ( 2.0_dp * rr(ixO^S) ** (3.0_dp/2.0_dp) * ( ( 4.0_dp * ( 1.0_dp - rr(ixO^S) ) * &
            DCOSH( (1.0_dp/3.0_dp) * DACOSH( ( threesqrtthree * rr(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / &
            ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) ) ) ** 2.0_dp - 1.0_dp ) / rr(ixO^S) + 1.0_dp) * &
            ( ( sqrtthree * DCOS( theta_profile(ixO^S) ) * &
            ( 1.0_dp / DCOSH( (1.0_dp/3.0_dp) * DACOSH( ( threesqrtthree * rr(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / &
            ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) ) ) ) ) / &
            ( 2.0_dp * DSQRT( 1.0_dp - rr(ixO^S) ) ) + 1.0_dp ) ** (3.0_dp/2.0_dp) ) )

          end where



          !Case 4: r<1 and xi<0:
          where((rr(ixO^S)<1.0_dp).and.((((rr(ixO^S)*costta(ixO^S))/2.0_dp)**2.0_dp)-&
          (((1.0_dp-rr(ixO^S))/3.0_dp)**3.0_dp)<0.0_dp))

            gravity_field_sphrc(ixO^S,1) = rr(ixO^S) ** (3.0_dp/2.0_dp) * &
            ( ( 4.0_dp * ( 1.0_dp - rr(ixO^S) ) * DCOS( (1.0_dp/3.0) * DACOS( ( threesqrtthree * &
            rr(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) ) ) ** 2.0_dp - &
            1.0_dp ) / rr(ixO^S) + 1.0_dp ) * DSQRT( ( sqrtthree * DCOS( theta_profile(ixO^S) ) * &
            ( 1.0_dp / DCOS( (1.0_dp/3.0_dp) * DACOS( (threesqrtthree * rr(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / &
            ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) ) ) ) ) / &
            ( 2.0_dp * DSQRT( 1.0_dp - rr(ixO^S) ) ) + 1.0_dp ) * &
            ( - ( ( ( 8.0_dp * ( 1.0_dp - rr(ixO^S) ) * ( ( 9.0_dp * sqrtthree * rr(ixO^S) * &
            DCOS( theta_profile(ixO^S) ) ) / ( 4.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (5.0_dp/2.0_dp) ) + &
            ( threesqrtthree * DCOS( theta_profile(ixO^S) ) ) / ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) ) * &
            DCOS( (1.0_dp/3.0_dp) * DACOS( ( threesqrtthree * rr(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / &
            ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) ) ) * &
            DSIN( (1.0_dp/3.0_dp) * DACOS( ( threesqrtthree * rr(ixO^S) * &
            DCOS( theta_profile(ixO^S) ) ) / ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) ) ) ) / &
            ( 3.0_dp * DSQRT( 1.0_dp - ( 27.0_dp * rr(ixO^S) ** 2.0_dp * DCOS( theta_profile(ixO^S) ) ** 2.0_dp) / &
            ( 4.0_dp * ( 1.0_dp - rr(ixO^S) ) ** 3.0_dp ) ) ) - 4.0_dp * DCOS( (1.0_dp/3.0_dp) * DACOS( ( threesqrtthree * rr(ixO^S) * &
            DCOS( theta_profile(ixO^S) ) ) / ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) ) ) ** 2.0_dp ) / rr(ixO^S) - &
            ( 4.0_dp * ( 1.0_dp - rr(ixO^S) ) * DCOS( (1.0_dp/3.0_dp) * &
            DACOS( ( threesqrtthree * rr(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / &
            ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) ) ) ** 2.0_dp - 1.0_dp ) / &
            rr(ixO^S) ** 2.0_dp ) / ( rr(ixO^S) ** (3.0_dp/2.0_dp) * &
            ( ( 4.0_dp * ( 1.0_dp - rr(ixO^S) ) * DCOS( (1.0_dp/3.0_dp) * DACOS( ( threesqrtthree * &
            rr(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) ) ) ** 2.0_dp - &
            1.0_dp ) / rr(ixO^S) + 1.0_dp ) ** 2.0_dp * DSQRT( ( sqrtthree * DCOS( theta_profile(ixO^S) ) * &
            ( 1.0_dp / DCOS( (1.0_dp/3.0_dp) * DACOS( ( threesqrtthree * rr(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / &
            ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) ) ) ) ) / &
            ( 2.0_dp * DSQRT( 1.0_dp - rr(ixO^S) ) ) + 1.0_dp ) )  - &
            ( ( sqrtthree * DCOS( theta_profile(ixO^S) ) * &
            ( 1.0_dp / DCOS( (1.0_dp/3.0_dp) * DACOS( ( threesqrtthree * rr(ixO^S) * &
            DCOS( theta_profile(ixO^S) ) ) / ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) ) ) ) ) / &
            ( 4.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) - &
            ( DCOS( theta_profile(ixO^S) ) * ( ( 9.0_dp * sqrtthree * rr(ixO^S) * &
            DCOS( theta_profile(ixO^S) ) ) / ( 4.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (5.0_dp/2.0_dp) ) + &
            ( threesqrtthree * DCOS( theta_profile(ixO^S) ) ) / ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) ) * &
            ( 1.0_dp / DCOS( (1.0_dp/3.0_dp) * DACOS( ( threesqrtthree * rr(ixO^S) * &
            DCOS( theta_profile(ixO^S) ) ) / ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) ) ) ) * &
            DTAN( (1.0_dp/3.0_dp) * DACOS( ( threesqrtthree * rr(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / &
            ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) ) ) ) / ( 2.0_dp * sqrtthree * &
            DSQRT( 1.0_dp - rr(ixO^S) ) * DSQRT( 1.0_dp - ( 27.0_dp * rr(ixO^S) ** 2.0_dp * DCOS( theta_profile(ixO^S) ) ** 2.0_dp ) / &
            ( 4.0_dp * ( 1.0_dp - rr(ixO^S) ) ** 3.0_dp ) ) ) ) / ( 2.0_dp * rr(ixO^S) ** (3.0_dp/2.0_dp) * &
            ( ( 4.0_dp * ( 1.0_dp - rr(ixO^S) ) * DCOS( (1.0_dp/3.0_dp) * DACOS( ( threesqrtthree * rr(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / &
            ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) ) ) ** 2.0_dp - 1.0_dp ) / &
            rr(ixO^S) + 1.0_dp ) * ( ( sqrtthree * DCOS( theta_profile(ixO^S) ) * &
            ( 1.0_dp / DCOS( (1.0_dp/3.0_dp) * DACOS( ( threesqrtthree * rr(ixO^S) * &
            DCOS( theta_profile(ixO^S) ) ) / ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) ) ) ) ) / &
            ( 2.0_dp * DSQRT( 1.0_dp - rr(ixO^S) ) ) + 1.0_dp ) ** (3.0_dp/2.0_dp) ) - &
            3.0_dp/( 2.0_dp * rr(ixO^S) ** (5.0_dp/2.0_dp) * ( ( 4.0_dp * ( 1.0_dp - rr(ixO^S) ) * &
            DCOS( (1.0_dp/3.0_dp) * DACOS( ( threesqrtthree * rr(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / &
            ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) ) ) ** 2.0_dp - 1.0_dp ) / rr(ixO^S) + 1.0_dp ) * &
            DSQRT( ( sqrtthree * DCOS( theta_profile(ixO^S) ) * &
            1.0_dp / DCOS( (1.0_dp/3.0_dp) * DACOS( ( threesqrtthree * rr(ixO^S) * &
            DCOS( theta_profile(ixO^S) ) ) / ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) ) ) ) / &
            ( 2.0_dp * DSQRT( 1.0_dp - rr(ixO^S) ) ) + 1.0_dp ) ) )


            gravity_field_sphrc(ixO^S,2) = DSQRT( rr(ixO^S) ) * &
            ( ( 4.0_dp * ( 1.0_dp - rr(ixO^S) ) * DCOS( (1.0_dp/3.0_dp) * &
            DACOS( ( threesqrtthree * rr(ixO^S) * cos( theta_profile(ixO^S) ) ) / &
            ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) ) ) ** 2.0_dp - 1.0_dp) / &
            rr(ixO^S) + 1.0_dp ) * DSQRT( ( sqrtthree * DCOS( theta_profile(ixO^S) ) * &
            ( 1.0_dp / DCOS( (1.0_dp/3.0_dp) * DACOS( ( threesqrtthree * rr(ixO^S) * &
            DCOS( theta_profile(ixO^S) ) ) / &
            ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) ) ) ) ) / &
            ( 2.0_dp * DSQRT( 1.0_dp - rr(ixO^S) ) ) + 1.0_dp ) * &
            ( ( 4.0_dp * sqrtthree * DSIN( theta_profile(ixO^S) ) * &
            DCOS( (1.0_dp/3.0_dp) * DACOS( ( threesqrtthree * rr(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / &
            ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) ) ) * DSIN( (1.0_dp/3.0_dp) * &
            DACOS( ( threesqrtthree * rr(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / &
            ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) ) ) ) / &
            ( DSQRT( 1.0_dp - rr(ixO^S) ) * rr(ixO^S) ** (3.0_dp/2.0_dp) * &
            DSQRT( 1.0_dp - ( 27.0_dp * rr(ixO^S) ** 2.0_dp * DCOS( theta_profile(ixO^S) ) ** 2.0_dp ) / &
            ( 4.0_dp * ( 1.0_dp - rr(ixO^S) ) ** 3.0_dp ) ) * &
            ( ( 4.0_dp * ( 1.0_dp - rr(ixO^S) ) * DCOS( (1.0_dp/3.0_dp) * DACOS( ( threesqrtthree * &
            rr(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / &
            ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) ) ) ** 2.0_dp - 1.0_dp ) / &
            rr(ixO^S) + 1.0_dp ) ** 2.0_dp * DSQRT( ( sqrtthree * DCOS( theta_profile(ixO^S) ) * &
            ( 1.0_dp / DCOS( (1.0_dp/3.0_dp) * DACOS( ( threesqrtthree * rr(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / &
            ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) ) ) ) ) / &
            ( 2.0_dp * DSQRT( 1.0_dp - rr(ixO^S) ) ) + 1.0_dp ) ) - &
            ( ( 3.0_dp * rr(ixO^S) * DSIN( theta_profile(ixO^S) ) * DCOS( theta_profile(ixO^S) ) * &
            DTAN( (1.0_dp/3.0_dp) * DACOS( ( threesqrtthree * rr(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / &
            ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) ) ) * &
            ( 1.0_dp / DCOS( (1.0_dp/3.0_dp) * DACOS( ( threesqrtthree * rr(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / &
            ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) ) ) ) ) / &
            ( 4.0_dp * ( 1.0_dp - rr(ixO^S) ) ** 2.0_dp * DSQRT( 1.0_dp - ( 27.0_dp * &
            rr(ixO^S) ** 2.0_dp * DCOS( theta_profile(ixO^S) ) ** 2.0_dp ) / &
            ( 4.0_dp * ( 1.0_dp - rr(ixO^S) ) ** 3.0_dp ) ) ) - &
            ( sqrtthree * DSIN( theta_profile(ixO^S) ) * &
            ( 1.0_dp / DCOS( (1.0_dp/3.0_dp) * DACOS( ( threesqrtthree * rr(ixO^S) * &
            DCOS( theta_profile(ixO^S) ) ) / ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) ) ) ) ) / &
            ( 2.0_dp * DSQRT( 1.0_dp - rr(ixO^S) ) ) ) / ( 2.0_dp * rr(ixO^S) ** (3.0_dp/2.0_dp) * &
            ( ( 4.0_dp * ( 1.0_dp - rr(ixO^S) ) * DCOS( (1.0_dp/3.0_dp) * DACOS( ( threesqrtthree * &
            rr(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / &
            ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) ) ) ** 2.0_dp - 1.0_dp ) / rr(ixO^S) + 1.0_dp ) * &
            ( ( sqrtthree * DCOS( theta_profile(ixO^S) ) * &
            1.0_dp / DCOS( (1.0_dp/3.0_dp) * DACOS( ( threesqrtthree * rr(ixO^S) * &
            DCOS( theta_profile(ixO^S) ) ) / &
            ( 2.0_dp * ( 1.0_dp - rr(ixO^S) ) ** (3.0_dp/2.0_dp) ) ) ) ) / &
            ( 2.0_dp * DSQRT( 1.0_dp - rr(ixO^S) ) ) + 1.0_dp ) ** (3.0_dp/2.0_dp) ) )

          end where



          !Case 1: r==1:
          where(DABS(rr(ixO^S)-1.0_dp)<=halfdrr)

          gravity_field_sphrc(ixO^S,1) = (- 3.0_dp * DCOS( theta_profile(ixO^S) ) ** (2.0_dp/3.0_dp) - &
          3.0_dp * rr(ixO^S) + 1.0_dp ) / ( 2.0_dp * rr(ixO^S) * ( 3.0_dp * &
          DCOS( theta_profile(ixO^S) ) ** (2.0_dp/3.0_dp) + rr(ixO^S) - 1.0_dp ) )

          gravity_field_sphrc(ixO^S,2) = ( DSIN( theta_profile(ixO^S) ) * ( 9.0_dp * &
          DCOS( theta_profile(ixO^S) ) ** (2.0_dp/3.0_dp) + &
          rr(ixO^S) + 5.0_dp ) ) / ( 3.0_dp * ( DCOS( theta_profile(ixO^S) ) ** (2.0_dp/3.0_dp) + &
          1.0_dp ) * DCOS( theta_profile(ixO^S) ) ** (1.0_dp/3.0_dp) * &
          ( 3.0_dp * DCOS( theta_profile(ixO^S) ) ** (2.0_dp/3.0_dp) + rr(ixO^S) - 1.0_dp ) )

          end where

          where(x(ixO^S,z_)<0.0_dp)
            gravity_field_sphrc(ixO^S,2) = -gravity_field_sphrc(ixO^S,2)
          end where

         gravity_field(ixO^S,r_) = (gravity_field_sphrc(ixO^S,1)*&
         x(ixO^S,r_)/DSQRT(x(ixO^S,r_)**2.0_dp+x(ixO^S,z_)**2.0_dp))&
         +(gravity_field_sphrc(ixO^S,2)*&
         x(ixO^S,z_)/DSQRT(x(ixO^S,r_)**2.0_dp+x(ixO^S,z_)**2.0_dp))
         gravity_field(ixO^S,z_) = (gravity_field_sphrc(ixO^S,1)*&
         x(ixO^S,z_)/DSQRT(x(ixO^S,r_)**2.0_dp+x(ixO^S,z_)**2.0_dp))&
         -(gravity_field_sphrc(ixO^S,2)*&
         x(ixO^S,r_)/DSQRT(x(ixO^S,r_)**2.0_dp+x(ixO^S,z_)**2.0_dp))



         {^D&gravity_field(ixO^S,^D) = gravity_field(ixO^S,^D) *&
         ism_surround(i_ism)%myconfig%profile_Beta\}



         if(ism_surround(i_ism)%myconfig%profile_add_pmgrav_to_hse)then
           gravity_field(ixO^S,r_)=gravity_field(ixO^S,r_)+&
           (-Ggrav*Mpoint/&
           (r_distance(ixO^S)*r_distance(ixO^S)))*x(ixO^S,r_)/DSQRT(x(ixO^S,r_)**2.0_dp+x(ixO^S,z_)**2.0_dp)
           gravity_field(ixO^S,z_)=gravity_field(ixO^S,z_)+&
           (-Ggrav*Mpoint/&
           (r_distance(ixO^S)*r_distance(ixO^S)))*x(ixO^S,z_)/DSQRT(x(ixO^S,r_)**2.0_dp+x(ixO^S,z_)**2.0_dp)
          end if



          where(r_distance(ixO^S)<=0.5_dp*DSQRT(SUM(dx_local_min(1:ndim)**2.0_dp)))
             {^D&gravity_field(ixO^S,^D)=0.0_dp\}
          end where

         case('Lee2001')


          alpha = ism_surround(i_ism)%myconfig%profile_p
          gravity_field(ixO^S,r_) = - alpha * (x(ixO^S,r_)**2.0_dp-&
          x(ixO^S,z_)**2.0_dp)/(x(ixO^S,r_)*(x(ixO^S,r_)**2.0_dp+&
          x(ixO^S,z_)**2.0_dp))

          gravity_field(ixO^S,z_) = - 2.0_dp * alpha * x(ixO^S,z_) / &
          (x(ixO^S,r_)**2.0_dp+x(ixO^S,z_)**2.0_dp)
            !call mpistop('the code arrives here !!! ')

          level = node(plevel_,saveigrid)
          ^D&dx_local(^D)=((xprobmax^D-xprobmin^D)/(domain_nx^D))/(2.0_dp**(level-1));
          where(DABS(x(ixO^S,r_))<=min(ism_surround(i_ism)%myconfig%escapencells(r_)*&
                dx_local(r_),ism_surround(i_ism)%myconfig%escapencellsglobal(r_)*dx(r_,1)))
              ^D&gravity_field(ixO^S,^D)=0.0_dp\
          end where

          case default

            gravity_field(ixO^S,1:ndim)=0.0_dp

          end select
          end do Loop_isms
        end if

  end subroutine special_effective_gravity

  !-----------------------------------------------------------
  !> Calculate gravitational acceleration in each dimension
  subroutine special_pointmass_gravity(ixI^L,ixO^L,wCT,x,gravity_field)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    real(dp), intent(in)    :: x(ixI^S,1:ndim)
    real(dp), intent(in)    :: wCT(ixI^S,1:nw)
    real(dp), intent(out)   :: gravity_field(ixI^S,ndim)
    real(dp)                        :: Ggrav, Mpoint, alpha
    real(kind=dp), dimension(ixI^S) :: r_distance
    real(kind=dp), dimension(ixI^S) :: theta_profile
    real(kind=dp), dimension(1:ndim) :: zero_dim,dx_local
    integer                      :: i_ism,idim,level
    !-----------------------------------

    Ggrav = constusr%G

    gravity_field(ixO^S,1:ndim)=0.0_dp
    {zero_dim(^D)=0.0_dp\}!-dx(1,1)



    ! Here we set the graviationnal acceleration inside the whole domain
    ! from the state vector wCT
    ! here the called wCT contains conservative variables
      if(usrconfig%ism_on)then
        i_ism =0
        call usr_distance(ixI^L,ixO^L,typeaxial,&
                          zero_dim,x,r_distance)
        call usr_get_theta(ixI^L,ixO^L,x,theta_profile)
        Mpoint = ism_surround(i_ism)%myconfig%profile_Mstar

        Loop_isms : do i_ism=0,usrconfig%ism_number-1
        select case(trim(ism_surround(i_ism)%myconfig%profile_density))
          case('Ulrich1976')
          if(.not.phys_config%gravity_hse)then
            select case(typeaxial)
              case('spherical')

                gravity_field(ixO^S,r_)=-Ggrav*Mpoint/(r_distance(ixO^S)*r_distance(ixO^S))

              case('cylindrical')

                gravity_field(ixO^S,r_)=(-Ggrav*Mpoint/&
                (r_distance(ixO^S)*r_distance(ixO^S)))*x(ixO^S,r_)/DSQRT(x(ixO^S,r_)**2.0_dp+x(ixO^S,z_)**2.0_dp)
                gravity_field(ixO^S,z_)=(-Ggrav*Mpoint/&
                (r_distance(ixO^S)*r_distance(ixO^S)))*x(ixO^S,z_)/DSQRT(x(ixO^S,r_)**2.0_dp+x(ixO^S,z_)**2.0_dp)

              case('slab','slabstretch')


                if(ndim<3)then
                  gravity_field(ixO^S,x_)=(-Ggrav*Mpoint/&
                  (r_distance(ixO^S)*r_distance(ixO^S)))*x(ixO^S,r_)/DSQRT(x(ixO^S,r_)**2.0_dp+x(ixO^S,z_)**2.0_dp)
                  gravity_field(ixO^S,y_)=(-Ggrav*Mpoint/&
                  (r_distance(ixO^S)*r_distance(ixO^S)))*x(ixO^S,z_)/DSQRT(x(ixO^S,r_)**2.0_dp+x(ixO^S,z_)**2.0_dp)
                  if(ndir>2)then
                    gravity_field(ixO^S,z_) = 0.0_dp
                  end if
                else
                  gravity_field(ixO^S,x_)= 0.0_dp ! TO DO
                  gravity_field(ixO^S,y_)= 0.0_dp ! TO DO
                  if(ndir>2)then
                     gravity_field(ixO^S,z_) = 0.0_dp ! TO DO
                  end if
                end if

              case default
                 call mpistop('Unknown typeaxial')
            end select

          else


            if(phys_config%use_gravity_g)then
              write(*,*) 'special_pointmass_gravity in mod_usr.t:'
              write(*,*) 'use_gravity_g is true but'
              write(*,*) 'but this part is not yet implemented'
              call mpistop('Code part not yet implemented')
            else

              write(*,*) 'special_pointmass_gravity in mod_usr.t:'
              write(*,*) 'use_gravity_g is false but'
              write(*,*) 'but this part is not yet implemented'
              call mpistop('Code part not yet implemented')
            end if

          end if
          !gravity_field(ixO^S,1:ndim)=gravity_field(ixO^S,1:ndim)/&
          !((usr_physunit%myconfig%length**3.0_dp)/(usr_physunit%myconfig%mass*&
          !usr_physunit%myconfig%time**2.0_dp))

          case('Lee2001')


          alpha = ism_surround(i_ism)%myconfig%profile_p
          gravity_field(ixO^S,r_) = - alpha * (x(ixO^S,r_)**2.0_dp-&
          x(ixO^S,z_)**2.0_dp)/(x(ixO^S,r_)*(x(ixO^S,r_)**2.0_dp+&
          x(ixO^S,z_)**2.0_dp))

          gravity_field(ixO^S,z_) = - 2.0_dp * alpha * x(ixO^S,z_) / &
          (x(ixO^S,r_)**2.0_dp+x(ixO^S,z_)**2.0_dp)
            !call mpistop('the code arrives here !!! ')

          level = node(plevel_,saveigrid)
          ^D&dx_local(^D)=((xprobmax^D-xprobmin^D)/(domain_nx^D))/(2.0_dp**(level-1));
          where(DABS(x(ixO^S,r_))<=min(ism_surround(i_ism)%myconfig%escapencells(r_)*&
                dx_local(r_),ism_surround(i_ism)%myconfig%escapencellsglobal(r_)*dx(r_,1)))
              ^D&gravity_field(ixO^S,^D)=0.0_dp\
          end where

          case default

            gravity_field(ixO^S,1:ndim)=0.0_dp

          end select
          end do Loop_isms
        end if

  end subroutine special_pointmass_gravity


  subroutine special_pointmass_gravity_potential(ixI^L,ixO^L,wCT,x,gravity_field)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    real(dp), intent(in)    :: x(ixI^S,1:ndim)
    real(dp), intent(in)    :: wCT(ixI^S,1:nw)
    real(dp), intent(out)   :: gravity_field(ixI^S)
    real(dp)                        :: Ggrav, Mpoint, alpha
    real(kind=dp), dimension(ixI^S) :: r_distance
    real(kind=dp), dimension(ixI^S) :: theta_profile
    real(kind=dp), dimension(1:ndim) :: zero_dim
    integer                      :: i_ism,idim
    !-----------------------------------





    Ggrav = constusr%G

    gravity_field(ixO^S)=0.0_dp
    {zero_dim(^D)=0.0_dp\}!-dx(1,1)



    ! Here we set the graviationnal acceleration inside the whole domain
    ! from the state vector wCT
    ! here the called wCT contains conservative variables
      if(usrconfig%ism_on)then
        i_ism =0
        call usr_distance(ixI^L,ixO^L,typeaxial,&
                          zero_dim,x,r_distance)
        call usr_get_theta(ixI^L,ixO^L,x,theta_profile)
        Mpoint = ism_surround(i_ism)%myconfig%profile_Mstar

        Loop_isms : do i_ism=0,usrconfig%ism_number-1
        select case(trim(ism_surround(i_ism)%myconfig%profile_density))
          case('Ulrich1976')
            gravity_field(ixO^S)=-Ggrav*Mpoint/r_distance(ixO^S)

            !gravity_field(ixO^S,1:ndim)=gravity_field(ixO^S,1:ndim)/&
            !((usr_physunit%myconfig%length**3.0_dp)/(usr_physunit%myconfig%mass*usr_physunit%myconfig%time**2.0_dp))
          case('Lee2001')

            alpha = ism_surround(i_ism)%myconfig%profile_p
            gravity_field(ixO^S) = alpha * &
            DLOG((r_distance(ixO^S)/DSIN(theta_profile(ixO^S)))/&
            (ism_surround(i_ism)%myconfig%profile_rw/1.0_dp)) !theta_0=pi/2
            !gravity_field(ixO^S,phi_) = 0.0_dp
            !treat theta=0 where dp/dR=-rho*dPhi/dR=0

            !call mpistop('the code arrives here !!! ')

          case default

            gravity_field(ixO^S)=0.0_dp

          end select
        end do Loop_isms


      end if

      ! [p] = [p_i,j - rho_i,j*(phi_i+1-phi_i)/2] ==> renormalise to more usefulf pressure:
      ! [phi]= [p/rho]

      ! de-normalize
      gravity_field(ixO^S)=gravity_field(ixO^S)*&
      usr_physunit%myconfig%length*&
      usr_physunit%myconfig%mass/&
      (unit_density*(unit_length/unit_velocity)**(2.0_dp)*&
      (usr_physunit%myconfig%length*usr_physunit%myconfig%length))

      !re-normalize:

      gravity_field(ixO^S)=gravity_field(ixO^S)/&
      (usr_physunit%myconfig%pressure/usr_physunit%myconfig%density)

  end subroutine special_pointmass_gravity_potential


  function special_pointmass_gravity_fpotential(ixI^L,ixO^L,wCT,x) result(gravity_field)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    real(dp), intent(in)    :: x(ixI^S,1:ndim)
    real(dp), intent(in)    :: wCT(ixI^S,1:nw)
    real(dp)                :: gravity_field(ixO^S)
    real(dp)                        :: Ggrav, Mpoint
    real(kind=dp), dimension(ixI^S) :: r_distance
    real(kind=dp), dimension(ixI^S) :: theta_profile
    real(kind=dp), dimension(1:ndim) :: zero_dim
    integer                      :: i_ism,idim
    !-----------------------------------



    Ggrav = constusr%G

    gravity_field(ixO^S)=0.0_dp
    {zero_dim(^D)=0.0_dp\}!-dx(1,1)



    ! Here we set the graviationnal acceleration inside the whole domain
    ! from the state vector wCT
    ! here the called wCT contains conservative variables
      if(usrconfig%ism_on)then
        i_ism =0
        call usr_distance(ixI^L,ixO^L,typeaxial,&
                          zero_dim,x,r_distance)
        call usr_get_theta(ixI^L,ixO^L,x,theta_profile)
        Mpoint = ism_surround(i_ism)%myconfig%profile_Mstar

        gravity_field(ixO^S)=-Ggrav*Mpoint/r_distance(ixO^S)

        !gravity_field(ixO^S,1:ndim)=gravity_field(ixO^S,1:ndim)/&
        !((usr_physunit%myconfig%length**3.0_dp)/(usr_physunit%myconfig%mass*usr_physunit%myconfig%time**2.0_dp))
      end if

      ! [p] = [p_i,j - rho_i,j*(phi_i+1-phi_i)/2] ==> renormalise to more usefulf pressure:
      ! [phi]= [p/rho]

      ! de-normalize
      gravity_field(ixO^S)=gravity_field(ixO^S)*&
      usr_physunit%myconfig%length*&
      usr_physunit%myconfig%mass/&
      (unit_density*(unit_length/unit_velocity)**(2.0_dp)*&
      (usr_physunit%myconfig%length*usr_physunit%myconfig%length))

      !re-normalize:

      gravity_field(ixO^S)=gravity_field(ixO^S)/&
      (usr_physunit%myconfig%pressure/usr_physunit%myconfig%density)

  end function special_pointmass_gravity_fpotential

  !> default usr parameters from a file
  subroutine usr_set_default_parameters
    !-------------------------------------
    usrconfig%phys_inuse                    = 'hd'
    usrconfig%phys_isotherm_on             = .false.
    usrconfig%phys_adiab                    = -1.0_dp
    usrconfig%phys_temperature_isotherm    = -1.0_dp
    usrconfig%physunit_on                   = .false.
    usrconfig%ism_on                        = .false.
    usrconfig%profile_use_hse               = .false.
    usrconfig%cloud_on                      = .false.
    usrconfig%jet_yso_on                    = .false.
    usrconfig%grackle_chemistry_on          = .true.
    usrconfig%cloud_number                  = 1
    usrconfig%ism_number                    = 1
    usrconfig%jet_yso_number                = 1
    usrconfig%ind_jet_associate_ism         = 0
    usrconfig%ind_cloud_associate_ism       = -1
    usrconfig%ism_list_diff                 = .false.
    usrconfig%cloud_list_diff               = .false.
    usrconfig%reset_medium                  = .false.
    usrconfig%coordinate_system             = 'slab'
    usrconfig%cloud_structure               = 0
    usrconfig%filename                      = 'mod_usr_jet_yso_cla.t'


    usrconfig%reset_flux_scheme_on          = .false.
    usrconfig%reset_limiter_on              = .false.
    usrconfig%reset_flux_scheme_diffuse     = 'none'
    usrconfig%reset_flux_scheme_old_method  = 'none'
    usrconfig%reset_flux_scheme_thresholdL1_max = 0.0_dp
    usrconfig%reset_flux_scheme_thresholdL1_min = 0.0_dp


    usrconfig%density_dusttogas_maxlimit    = 1.0d6
    usrconfig%density_dusttogas_minlimit    = 1.0d-6


    usrconfig%temperature_max               = bigdouble

    usrconfig%cloud_profile_density_on      = .false.
    usrconfig%cloud_profile_pressure_on     = .false.
    usrconfig%cloud_profile_velocity_on     = .false.




  end subroutine usr_set_default_parameters
  !------------------------------------------------------------------
  !> Read this module s parameters from a file
  subroutine usr_params_read(files)
    use mod_hd, only : link_hd_pre_phys
    character(len=*), intent(in) :: files(:)
    ! .. local ..
    integer                      :: i_file,i_reason
    character(len=70)            :: error_message
    integer                      :: i_cloud,i_ism,i_jet_yso
    !-------------------------------------
    namelist /usr_list/ usrconfig
    !-------------------------------------

    error_message = 'At '//trim(usrconfig%filename)//'  in the procedure : usr_params_read'

    if(mype==0)write(*,*)'Reading usr_list'
    Loop_ifile : do i_file = 1, size(files)
       open(unitpar, file=trim(files(i_file)), status="old")
       read(unitpar, usr_list, iostat=i_reason)
       cond_ierror : if(i_reason>0)then
        write(*,*)' Error in reading the parameters file : ',trim(files(i_file))
        write(*,*)' Error at namelist: ', trim(usrconfig%filename)
        write(*,*)' The code stops now '
        call mpistop(trim(error_message))
       elseif(i_reason<0)then cond_ierror
        write(*,*)' Reache the end of the file  : ',trim(files(i_file))
        write(*,*)' Error at namelist: usr_list'
        write(*,*)' The code stops now '
        call mpistop(trim(error_message))
       else cond_ierror
        write(*,*)' End of reading of the usr_list'
       end if cond_ierror
       close(unitpar)
    end do Loop_ifile


    ! Must pre-default and rpre-ead mhd/hd_params
    ! first for Grackle configuration
    select case(trim(usrconfig%phys_inuse))
    case('hd')
     call link_hd_pre_phys()
    case('mhd')
     call link_hd_pre_phys()
    end select
    call phys_default_pre_config()
    call phys_pre_read_params(files)

    !link which treatment is made of energy and temperature
    usrconfig%phys_isotherm_on = pre_phys_config%isotherm_on
    usrconfig%phys_temperature_isotherm = pre_phys_config%temperature_isotherm
    usrconfig%phys_adiab = pre_phys_config%adiab
    !end : link which treatment is made of energy and temperature


    !>Link mod_hd use_grackle answer to usrconfig%grackle_chemistry_on:
    usrconfig%grackle_chemistry_on = pre_phys_config%use_grackle
    write(*,*) 'usrconfig%grackle_chemistry_on = ',usrconfig%grackle_chemistry_on


    !1) Read user-defined physical units
    if(usrconfig%physunit_on)then
      call usr_physunit%read_parameters(usr_physunit%myconfig,files)
    else
      call usr_unit_read(files)
      call usr_physunit%set_to_one
    end if


    !2) If ism on, set all usr ISMs parameters to their default values
    ! and then read their user-defined parameters
    if(usrconfig%ism_on)then
      allocate(ism_surround(0:usrconfig%ism_number-1))
      Loop_allism : do i_ism =0,usrconfig%ism_number-1
       ism_surround(i_ism)%myconfig        = ism_default%myconfig
       ism_surround(i_ism)%mydust%myconfig = ism_default%mydust%myconfig
       ism_surround(i_ism)%myboundaries%myconfig = ism_default%myboundaries%myconfig
       ism_surround(i_ism)%myconfig%myindice=i_ism
       call ism_surround(i_ism)%read_parameters(ism_surround(i_ism)%myconfig,files)
      end do Loop_allism
    else
      if(usrconfig%ism_number>0)then
        usrconfig%ism_number = 0
      end if
    end if

    !3) If cloud on, set all usr clouds parameters to their default values
    ! and then read their user-defined parameters
    if(usrconfig%cloud_on)then
      allocate(cloud_medium(0:usrconfig%cloud_number-1))
      Loop_allcloud : do i_cloud =0,usrconfig%cloud_number-1
       cloud_medium(i_cloud)%myconfig          = cloud_default%myconfig
       cloud_medium(i_cloud)%mydust%myconfig   = cloud_default%mydust%myconfig
       cloud_medium(i_cloud)%myconfig%myindice = i_cloud
       !cloud_medium(i_cloud)%myboundaries%myconfig = cloud_default%myboundary%myconfig
       call cloud_medium(i_cloud)%read_parameters(files,cloud_medium(i_cloud)%myconfig)
      end do Loop_allcloud
    else
      if(usrconfig%cloud_number>0)then
        usrconfig%cloud_number = 0
      end if
    end if

    !4) If jet on, set all usr jets parameters to their default values
    ! and then read their user-defined parameters
    if(usrconfig%jet_yso_on)then
      allocate(jet_yso(0:usrconfig%jet_yso_number-1))
      Loop_alljetagn : do i_jet_yso =0,usrconfig%jet_yso_number-1
       jet_yso(i_jet_yso)%myconfig          = jet_yso_default%myconfig
       jet_yso(i_jet_yso)%mydust%myconfig   = jet_yso_default%mydust%myconfig
       !jet_yso(i_jet_yso)%myboundaries%myconfig = jet_yso_default%myboundary%myconfig
       jet_yso(i_jet_yso)%myconfig%myindice = i_jet_yso

       call jet_yso(i_jet_yso)%read_parameters(files,jet_yso(i_jet_yso)%myconfig)
     end do Loop_alljetagn
   else
    if(usrconfig%jet_yso_number>0)then
     usrconfig%jet_yso_number = 0
    end if
  end if

  !5) If grackle_chemistry_on on, set all usr grackle config parameters to their default values
  ! and then read their user-defined parameters
  if(usrconfig%grackle_chemistry_on)then

    !i)
    grackle_object%myconfig%use_grackle = 0

    grackle_object%myconfig          = grackle_default%myconfig
    grackle_object%myparams          = grackle_default%myparams
    allocate(tstst(4))
    call grackle_object%read_parameters(grackle_object%myconfig,files)

    !>Link mod_hd answers (phys_config%gr_primordial_chemistry,phys_config%gr_metal_cooling
    ! phys_config%gr_dust_chemistry,phys_config%gr_with_radiative_cooling) to grackle_object%myconfig:

    !ii)
    write(*,*) 'pre_phys_config%gr_primordial_chemistry = ',pre_phys_config%gr_primordial_chemistry
    grackle_object%myconfig%gr_primordial_chemistry = pre_phys_config%gr_primordial_chemistry

    !iii)
    if(pre_phys_config%gr_with_radiative_cooling)then
      grackle_object%myconfig%gr_with_radiative_cooling = 1
    else
      grackle_object%myconfig%gr_with_radiative_cooling = 0
    end if

    !iv)
    if(pre_phys_config%gr_metal_cooling)then
      grackle_object%myconfig%gr_metal_cooling = 1
    else
      grackle_object%myconfig%gr_metal_cooling = 0
    end if

    !v)
    if(pre_phys_config%gr_dust_chemistry)then
      grackle_object%myconfig%gr_dust_chemistry = 1
    else
      grackle_object%myconfig%gr_dust_chemistry = 0
    end if


    ! Gamma of grackle defined from hd_list->hd_gamma:
    grackle_object%myconfig%gr_gamma = pre_phys_config%gamma
    write(*,*) 'grackle_object%myconfig%gr_gamma = ', grackle_object%myconfig%gr_gamma

    write(*,*) 'grackle_object%myconfig%use_grackle : ', grackle_object%myconfig%use_grackle
    write(*,*) 'grackle_object%myconfig%gr_primordial_chemistry : ', grackle_object%myconfig%gr_primordial_chemistry
    write(*,*) 'grackle_object%myconfig%gr_metal_cooling : ', grackle_object%myconfig%gr_metal_cooling
    write(*,*) 'grackle_object%myconfig%gr_dust_chemistry : ', grackle_object%myconfig%gr_dust_chemistry


    write(*,*) 'Test use_grackle == 0 : ',grackle_object%myconfig%use_grackle

    call grackle_object%allocate_fields_parameters(usrconfig%ism_number,usrconfig%jet_yso_number,usrconfig%cloud_number)
    write(*,*) 'Number of elements:'
    write(*,*) '> ISM:',usrconfig%ism_number
    write(*,*) '> Jet:',usrconfig%jet_yso_number
    write(*,*) '> Cloud:',usrconfig%cloud_number
    write(*,*) '> Total:',usrconfig%ism_number+usrconfig%jet_yso_number+usrconfig%cloud_number
    write(*,*) 'Size of grackle fields parameters arrays'
    write(*,*) '(grackle_object%myparams%number_of_objects):', grackle_object%myparams%number_of_objects

    call grackle_object%set_fields_default(usrconfig%ism_on,usrconfig%jet_yso_on,usrconfig%cloud_on)

    write(*,*) 'Test gr_H2I_density == 2.0 : ', gr_H2I_density
    write(*,*) 'Test gr_patches_name == ism,jet, or cloud : ', gr_patches_name

    call grackle_object%read_fields_parameters(files)

    write(*,*) 'Test gr_patches_name == ism,jet, or cloud : ', gr_patches_name
    write(*,*) 'Test gr_patches_name == ism_uniform,jet_uniform, or cloud_uniform : ', gr_density_method

    ! do the following in the source adding subroutine to save performance and memory:
    !iresult = set_default_chemistry_parameters(gr_struct%grackle_data)
    !write(*,*) 'Grackle structure configuration defaulting successfully done !'
   else
    grackle_object%myconfig%use_grackle = 0
    grackle_object%myconfig%gr_with_radiative_cooling = 0
    grackle_object%myconfig%gr_primordial_chemistry = 0
    grackle_object%myconfig%gr_metal_cooling = 0
    grackle_object%myconfig%gr_dust_chemistry = 0
   end if
  end subroutine usr_params_read

  !> subroutine to clean memory at the end
  subroutine usr_clean_memory_final
    if(usrconfig%ism_on)then
      if(allocated(ism_surround))deallocate(ism_surround)
    end if
    if(usrconfig%cloud_on)then
      if(allocated(cloud_medium))deallocate(cloud_medium)
    end if
    if(allocated(the_dust_inuse))deallocate(the_dust_inuse)
  end subroutine usr_clean_memory_final

!-------------------------------------------------------------------
!> subroutine read unit used in the code
  subroutine usr_unit_read(files)
   implicit none
   character(len=*), intent(in) :: files(:)
   integer                      :: i_file

   namelist /usr_unit_list/ unit_length , unit_time,unit_velocity,          &
                      unit_density, unit_numberdensity,                     &
                      unit_pressure,unit_temperature


  if(mype==0)write(*,*)'Reading usr_unit_list'
  Loop_read_usrfile : do i_file = 1, size(files)
         open(unitpar, file=trim(files(i_file)), status="old")
         read(unitpar, usr_unit_list, end=109)
  109    close(unitpar)
  end do Loop_read_usrfile
 end subroutine usr_unit_read
  !-----------------------------------------------------------
  !> subroutine to check configuration conflits
  subroutine usr_check_conflict
    implicit none
    ! .. local ..
    integer  :: i_ism,i_cloud,i_jet_yso
    !------------------------------

    cond_dust_on : if(.not.phys_config%dust_on)then
      if(usrconfig%ism_on)then
       Loop_isms : do i_ism=0,usrconfig%ism_number-1
        ism_surround(i_ism)%myconfig%dust_on =.false.
       end do Loop_isms
      end if
      if(usrconfig%cloud_on)then
       Loop_clouds : do i_cloud=0,usrconfig%cloud_number-1
        cloud_medium(i_cloud)%myconfig%dust_on =.false.
       end do Loop_clouds
      end if

      if(usrconfig%jet_yso_on)then
       Loop_jet_yso : do i_jet_yso=0,usrconfig%jet_yso_number-1
        jet_yso(i_jet_yso)%myconfig%dust_on =.false.
      end do Loop_jet_yso
      end if


    end if  cond_dust_on

    if(.not.usrconfig%reset_medium )then
      if(usrconfig%ism_on)then
        Loop_isms2 : do i_ism=0,usrconfig%ism_number-1
          ism_surround(i_ism)%myconfig%reset_coef=-1.0_dp
        end do  Loop_isms2
      end if
    end if

    usrconfig%cloud_profile_on=usrconfig%cloud_profile_density_on   .or.&
                               usrconfig%cloud_profile_pressure_on  .or.&
                               usrconfig%cloud_profile_velocity_on

  end   subroutine usr_check_conflict
  !-----------------------------------------------------------
  !> subroutine to normalize parameters in the code
  subroutine usr_normalise_parameters
   implicit none
   integer            :: idust



   constusr%G         = constusr%G*&
                      (unit_density*(unit_length/unit_velocity)**(2.0_dp))


   constusr%clight                      = constusr%clight/unit_velocity

   ! complet all physical unit in use
   if(usrconfig%physunit_on) then
      call usr_physunit%fillphysunit
    end if


    w_convert_factor(phys_ind%rho_)                       = unit_density
    if(phys_config%energy)w_convert_factor(phys_ind%e_)   = unit_density*unit_velocity**2.0
    if(saveprim)then
     w_convert_factor(phys_ind%mom(:))                    = unit_velocity
    else
     w_convert_factor(phys_ind%mom(:))                    = unit_density*unit_velocity
    end if
    time_convert_factor                                   = unit_time
    length_convert_factor                                 = unit_length


    if(trim(usrconfig%reset_flux_scheme_old_method)=='none')then
      usrconfig%reset_flux_scheme_old_method='hllc'
    end if
    if(usrconfig%reset_flux_scheme_on) then
      if(trim(usrconfig%reset_flux_scheme_diffuse)=='none')then
        usrconfig%reset_flux_scheme_diffuse='tvdlf'
      end if
    end if

    if(phys_config%dust_on)then
          w_convert_factor(phys_ind%dust_rho(:))              = unit_density
          if(saveprim)then
            Loop_dust1: do idust = 1, phys_config%dust_n_species
              w_convert_factor(phys_ind%dust_mom(:,idust))     = unit_velocity
            end do Loop_dust1
          else
            Loop_dust2: do idust = 1, phys_config%dust_n_species
              w_convert_factor(phys_ind%dust_mom(:,idust))     = unit_density*unit_velocity
            end do Loop_dust2
          end if
    end if
   usrconfig%reset_flux_scheme_thresholdL1_max= &
                                  usrconfig%reset_flux_scheme_thresholdL1_max &
                                 / usr_physunit%myconfig%luminosity
   usrconfig%reset_flux_scheme_thresholdL1_min= &
                                  usrconfig%reset_flux_scheme_thresholdL1_min &
                                 / usr_physunit%myconfig%luminosity
   usrconfig%temperature_max = usrconfig%temperature_max/usr_physunit%myconfig%temperature
   cond_iso : if(.not.phys_config%energy) then
      usrconfig%phys_adiab = usrconfig%phys_adiab/usr_physunit%myconfig%velocity**2.0_dp
      usrconfig%phys_temperature_isotherm = usrconfig%phys_temperature_isotherm / &
          usr_physunit%myconfig%temperature
   end if cond_iso
  end subroutine usr_normalise_parameters


!-------------------------------------------------------------------------
subroutine initglobaldata_usr
 use mod_variables
 implicit none
 ! .. local ..
 integer        :: i_cloud,i_ism,i_jet_yso,n_objects,n_object_w_dust
 real(kind=dp)  :: mp,kB
 type(gr_objects) :: gr_obj
 !------------------------------------
 !/



 call gr_obj%test_chemistry('mod_usr_yso_jet : initglobaldata_usr => usr_set_parameters')



  n_objects       = 0
  n_object_w_dust = 0
  itr             = 1


  if(SI_unit) then
      mp=mp_SI
      kB=kB_SI
  else
      mp=mp_cgs
      kB=kB_cgs
  end if

  ! complet usr configuration
  cond_isotherm_0 : if(usrconfig%phys_isotherm_on)then
    if(usrconfig%phys_temperature_isotherm>0)then
      ! isotherm case : need to compute adiab = csound**2
      usrconfig%phys_adiab = kB/(phys_config%mean_mup*mp)*&
      usrconfig%phys_temperature_isotherm
    elseif(usrconfig%phys_adiab>=0)then
      ! isotherm case : need to compute isotherm temperature T
      usrconfig%phys_temperature_isotherm = (usrconfig%phys_adiab/kB)*&
      phys_config%mean_mup*mp
    end if
    ! normalise value for isotherm temperature
    usrconfig%phys_temperature_isotherm = usrconfig%phys_temperature_isotherm/&
    unit_temperature
    ! normalise value for adiab
    usrconfig%phys_adiab  = usrconfig%phys_adiab/&
    (unit_velocity*unit_velocity)

    PRINT*,' is the same ',phys_config%adiab,usrconfig%phys_adiab
    if(dabs(phys_config%adiab*phys_config%unit_velocity**2.0_dp-&
    usrconfig%phys_adiab*unit_velocity**2.0_dp)>smalldouble)then
      write(*,*) 'the constant adiab is not well set in phys or at user site:'
      write(*,*) 'phys_config%adiab = ', phys_config%adiab*&
      usr_physunit%myconfig%velocity**2.0_dp
      write(*,*) 'usr_config%adiab = ', usrconfig%phys_adiab*&
      phys_config%unit_velocity**2.0_dp
      write(*,*) 'phys_config%temperature_isotherm = ', phys_config%temperature_isotherm*&
      phys_config%unit_temperature
      write(*,*) 'usrconfig%phys_temperature_isotherm = ', usrconfig%phys_temperature_isotherm*&
      usr_physunit%myconfig%temperature
      write(*,*) ' mod_usr.t : mean_mup = ', phys_config%mean_mup
    end if

    if(dabs(pre_phys_config%gamma-1.0_dp)>smalldouble)then
     write(*,*) 'pre_phys_config%gamma = ',pre_phys_config%gamma
     call mpistop('Isothermal perfect gas but hd_gamma/=1')
    else
     write(*,*) '*********************************************************'
     write(*,*) 'MPI-AMRVAC runs HD adiabatic isothermal run (no energy solved)'
     write(*,*) 'with temperature = ',usrconfig%phys_temperature_isotherm*unit_temperature
     write(*,*) 'and c_sound**2 = ',usrconfig%phys_adiab*(unit_velocity*unit_velocity)
     write(*,*) '*********************************************************'
    end if
  end if cond_isotherm_0



 ! complet ism parameters
 if(usrconfig%ism_on)then
  Loop_isms : do i_ism=0,usrconfig%ism_number-1
   ism_surround(i_ism)%myconfig%itr=itr
   if(ism_surround(i_ism)%myconfig%Mach_number_tomov_obj>0.0_dp) then
     if(usrconfig%jet_yso_on.and.usrconfig%ind_jet_associate_ism>=0)then
      ism_surround(i_ism)%myconfig%velocity_ofmov_obj=jet_yso(usrconfig%ind_jet_associate_ism)%myconfig%velocity
     else if(usrconfig%cloud_on.and.usrconfig%ind_cloud_associate_ism>=0)then
      ism_surround(i_ism)%myconfig%velocity_ofmov_obj=jet_yso(usrconfig%ind_cloud_associate_ism)%myconfig%velocity
     end if
   end if
   if(usrconfig%phys_isotherm_on) then
     if(ism_surround(i_ism)%myconfig%temperature>0.0_dp)then
       write(*,*) ' the ism temperature will reset to the imposed isotherm temperarture:', &
         usrconfig%phys_temperature_isotherm*unit_temperature, ' K'
     end if
     ism_surround(i_ism)%myconfig%temperature = usrconfig%phys_temperature_isotherm
   end if
   call ism_surround(i_ism)%set_complet
   call ism_surround(i_ism)%normalize(usr_physunit)
   if(ism_surround(i_ism)%myconfig%dust_on)n_object_w_dust=n_object_w_dust+1
  end do Loop_isms
  itr=ism_surround(usrconfig%ism_number-1)%myconfig%itr+1
  n_objects = n_objects + usrconfig%ism_number
 end if



 ! complet cloud parameters
 cond_cloud_on : if(usrconfig%cloud_on) then
   Loop_clouds : do i_cloud=0,usrconfig%cloud_number-1
    cloud_medium(i_cloud)%myconfig%itr=itr
    if(usrconfig%phys_isotherm_on) then
      if(cloud_medium(i_cloud)%myconfig%temperature>0.0_dp)then
       write(*,*) ' the cloud temperature will reset to the imposed isotherm temperarture:', &
         usrconfig%phys_temperature_isotherm, ' K'
      end if
      cloud_medium(i_cloud)%myconfig%temperature = usrconfig%phys_temperature_isotherm
    end if
    call cloud_medium(i_cloud)%set_complet
    call cloud_medium(i_cloud)%normalize(usr_physunit)
    if(cloud_medium(i_cloud)%myconfig%dust_on)n_object_w_dust=n_object_w_dust+1
   end do Loop_clouds
   itr=cloud_medium(usrconfig%cloud_number-1)%myconfig%itr+1
   n_objects = n_objects + usrconfig%cloud_number
 end if cond_cloud_on


 ! complet jet parameters
 if(usrconfig%jet_yso_on) then
   Loop_jet_yso : do i_jet_yso=0,usrconfig%jet_yso_number-1
    jet_yso(i_jet_yso)%myconfig%itr=itr
    jet_yso(i_jet_yso)%myconfig%pressure_associate_ism = ism_surround(0)%myconfig%pressure &
                                                         *usr_physunit%myconfig%pressure

    jet_yso(i_jet_yso)%myconfig%density_associate_ism  = ism_surround(0)%myconfig%density &
                                                        *usr_physunit%myconfig%density
    if(usrconfig%phys_isotherm_on) then
      if(jet_yso(i_jet_yso)%myconfig%temperature>0.0_dp)then
       write(*,*) ' the jet temperature will reset to the imposed isotherm temperarture:', &
         usrconfig%phys_temperature_isotherm, ' K'
      end if
      jet_yso(i_jet_yso)%myconfig%temperature = usrconfig%phys_temperature_isotherm
    end if
    call jet_yso(i_jet_yso)%set_complet
    call jet_yso(i_jet_yso)%normalize(usr_physunit)
    if(jet_yso(i_jet_yso)%myconfig%dust_on)n_object_w_dust=n_object_w_dust+1
   end do Loop_jet_yso
   itr=jet_yso(usrconfig%jet_yso_number-1)%myconfig%itr+1
   n_objects = n_objects + usrconfig%jet_yso_number
 else
   zjet_=min(z_,ndim)
   thetajet_ = z_
   rjet_=r_
 end if


 ! complet grackle units in gr_struct
 ! do the following in the source adding subroutine to save performance and memory:
 if(usrconfig%grackle_chemistry_on)then
   call grackle_object%set_complet
   write(*,*) 'Grackle usr configuration completing/consistency check successfully done!'

 end if

 ! normalize grackle units in gr_struct
 ! do the following in the source adding subroutine to save performance and memory:
 !if(usrconfig%grackle_chemistry_on)then
   !call ism_surround(i_ism)%normalize(usr_physunit)
   !call grackle_object%normalize(grackle_structure,usr_physunit)
   !write(*,*) 'Grackle structure normalization successfully done!'

 !end if
 write(*,*) 'End of non-user normalization/completing successfully done!'


 !Add here the set_complet of the fluxes !!!!

 if(phys_config%dust_on)allocate(the_dust_inuse(n_object_w_dust))

 call usr_normalise_parameters
 if(mype==0)call usr_write_setting



end subroutine initglobaldata_usr
  !> The initial conditions
  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    ! initialize one grid
    ! in the src code amrini.t, ixI^L = ixG^LL = ixGlo1,ixGlo2,ixGhi1,ixGhi2
    ! = range delimiting the whole domain (including boundary ghost cells), i.e. between integers
    ! 1 and the highest possible indices for the coordinates for the grid for each dimension
    ! and ixO^L = ixM^LL = ixMlo1,ixMlo2,ixMhi1,ixMhi2 =
    ! = range delimiting the mesh (i.e. the grid without boundary layers)
    implicit none

    integer, intent(in)     :: ixI^L,ixO^L
    real(dp), intent(in)    :: x(ixI^S,1:ndim)
    real(dp), intent(inout) :: w(ixI^S,1:nw)
    !.. local ..
    real(dp)      :: res
    integer       :: ix^D,na,flag(ixI^S)
    integer       :: i_cloud,i_ism,i_jet_yso,i_dust,i_start,i_end
    logical, save :: first=.true.
    logical       :: patch_all(ixI^S)
    logical       :: patch_inuse(ixI^S)
    type(dust)    :: dust_dummy
    integer       :: i_object_w_dust
    real(dp)      :: cloud_profile(ixI^S,1:nw)
    real(kind=dp), dimension(ixI^S) :: theta_profile,d_profile,r_normalized
    real(kind=dp), dimension(ixI^S) :: cos_theta_zero, theta_zero
    ! .. only test ..
    real(dp)      ::old_w(ixO^S,1:nw)
    !-----------------------------------------/
    patch_all(ixO^S) = .true.
    if(first)then
      if(mype==0) then
        write(*,*)'Jet start :-)'
      endif
      first=.false.
    endif

    !Let s initiate theta and cos_theta_zero
    call usr_get_theta(ixI^L,ixO^L,x,theta_profile)
    w(ixI^S,phys_ind%mythetafield_)=theta_profile(ixI^S)
    if(ism_surround(i_ism)%myconfig%profile_density_on)then
      if(trim(ism_surround(i_ism)%myconfig%profile_density)=='Ulrich1976')then

      call ism_surround(i_ism)%set_profile_distance(ixI^L,ixO^L,x,d_profile)
      r_normalized(ixI^S) = (d_profile(ixI^S)/ism_surround(i_ism)%myconfig%profile_rd)
      call usr_ulrich1976_costheta_zero(ixI^L, ixO^L,x,r_normalized,&
      theta_profile,cos_theta_zero)
      theta_zero(ixI^S) = DACOS(cos_theta_zero(ixI^S))
      w(ixI^S,phys_ind%mytheta_zero_)=DCOS(theta_zero(ixI^S))

      end if
    end if
    i_object_w_dust=0
    ! set the ism
    cond_ism_on : if(usrconfig%ism_on) then
      Loop_isms : do i_ism=0,usrconfig%ism_number-1
       ism_surround(i_ism)%subname='initonegrid_usr'
       call ism_surround(i_ism)%set_w(ixI^L,ixO^L,global_time,x,w)
       patch_all(ixO^S) =  patch_all(ixO^S) .and. .not.ism_surround(i_ism)%patch(ixO^S)
       if(ism_surround(i_ism)%myconfig%dust_on)then
         i_object_w_dust = i_object_w_dust +1
         if(.not.allocated(the_dust_inuse(i_object_w_dust)%patch))allocate(the_dust_inuse(i_object_w_dust)%patch(ixI^S))
         the_dust_inuse(i_object_w_dust)%myconfig    = ism_surround(i_ism)%mydust%myconfig
         if(.not.allocated(the_dust_inuse(i_object_w_dust)%the_ispecies))&
         allocate(the_dust_inuse(i_object_w_dust)%the_ispecies&
                   (the_dust_inuse(i_object_w_dust)%myconfig%idust_first:&
                    the_dust_inuse(i_object_w_dust)%myconfig%idust_last))
         the_dust_inuse(i_object_w_dust)%patch(ixO^S)=ism_surround(i_ism)%mydust%patch(ixO^S)
       end if
      end do Loop_isms
    end if cond_ism_on


    ! set jet agn
    cond_jet_on : if(usrconfig%jet_yso_on)then
      Loop_jet_yso : do i_jet_yso=0,usrconfig%jet_yso_number-1
       jet_yso(i_jet_yso)%subname='initonegrid_usr'
       call jet_yso(i_jet_yso)%set_w(ixI^L,ixO^L,global_time,x,w)

       patch_inuse(ixO^S) = jet_yso(i_jet_yso)%patch(ixO^S)
       cond_ism_onjet : if(usrconfig%ism_on)then
        Loop_isms2 : do i_ism=0,usrconfig%ism_number-1
          if(ism_surround(i_ism)%myconfig%tracer_on)then
           where(patch_inuse(ixO^S))
             w(ixO^S,phys_ind%tracer(ism_surround(i_ism)%myconfig%itr))=0.0_dp
           end where
          end if


          cond_ism_dust_onjet : if(ism_surround(i_ism)%myconfig%dust_on) then
            ism_surround(i_ism)%mydust%patch(ixO^S)=.not.patch_inuse(ixO^S)
            i_start= ism_surround(i_ism)%mydust%myconfig%idust_first
            i_end  = ism_surround(i_ism)%mydust%myconfig%idust_last
            Loop_ism_idustjet :  do i_dust=i_start,i_end
              ism_surround(i_ism)%mydust%the_ispecies(i_dust)%patch(ixO^S)=&
                                           .not.patch_inuse(ixO^S)
            end do Loop_ism_idustjet
          call ism_surround(i_ism)%mydust%set_w_zero(ixI^L,ixO^L,x,w)
          end if cond_ism_dust_onjet
        end do Loop_isms2
       end if cond_ism_onjet

       patch_all(ixO^S) =  patch_all(ixO^S) .and. .not.patch_inuse(ixO^S)
       if(jet_yso(i_jet_yso)%myconfig%dust_on)then
         i_object_w_dust = i_object_w_dust +1
         if(.not.allocated(the_dust_inuse(i_object_w_dust)%patch))&
                  allocate(the_dust_inuse(i_object_w_dust)%patch(ixI^S))
         the_dust_inuse(i_object_w_dust)%myconfig    = jet_yso(i_jet_yso)%mydust%myconfig
         if(.not.allocated(the_dust_inuse(i_object_w_dust)%the_ispecies))&
         allocate(the_dust_inuse(i_object_w_dust)%the_ispecies&
                   (the_dust_inuse(i_object_w_dust)%myconfig%idust_first:&
                    the_dust_inuse(i_object_w_dust)%myconfig%idust_last))
         the_dust_inuse(i_object_w_dust)%patch(ixO^S)= jet_yso(i_jet_yso)%patch(ixO^S)
       end if
     end do Loop_jet_yso
    end if cond_jet_on


    ! set one cloud
    cloud_on : if(usrconfig%cloud_on)then

      Loop_clouds : do i_cloud=0,usrconfig%cloud_number-1
       cloud_medium(i_cloud)%subname='initonegrid_usr'
       if(allocated(cloud_medium(i_cloud)%patch))deallocate(cloud_medium(i_cloud)%patch)
       allocate(cloud_medium(i_cloud)%patch(ixI^S))
       cloud_medium(i_cloud)%patch(ixO^S) = x(ixO^S,zjet_)>cloud_medium(i_cloud)%myconfig%extend(zjet_) +&
                                dtan(cloud_medium(i_cloud)%myconfig%extend(3)*&
                                     usr_physunit%myconfig%length*dpi/180.0_dp)*(x(ixO^S,1)-xprobmin1)


       cond_cloud_profile_on : if(usrconfig%cloud_profile_on)then
        call usr_set_profile_cloud(ixI^L,ixO^L,x,cloud_medium(i_cloud),cloud_profile)
       else cond_cloud_profile_on
        cloud_profile      = 1.0_dp
       end if cond_cloud_profile_on

       call cloud_medium(i_cloud)%set_w(ixI^L,ixO^L,global_time,x,w,&
                        usr_density_profile=cloud_profile(ixI^S,phys_ind%rho_),&
                        usr_pressure_profile=cloud_profile(ixI^S,phys_ind%pressure_),&
                        usr_velocity_profile=cloud_profile(ixI^S,phys_ind%mom(1):phys_ind%mom(ndir)))

       patch_inuse(ixO^S) = cloud_medium(i_cloud)%patch(ixO^S)
       ism_is_oncld : if(usrconfig%ism_on)then
        Loop_isms_cloud0 : do i_ism=0,usrconfig%ism_number-1
          cond_ism_tracer_oncld : if(ism_surround(i_ism)%myconfig%tracer_on)then
            where(patch_inuse(ixO^S))
             w(ixO^S,phys_ind%tracer(ism_surround(i_ism)%myconfig%itr))=0.0_dp
            end where
          end if cond_ism_tracer_oncld
          cond_ism_dust_oncld : if(ism_surround(i_ism)%myconfig%dust_on) then
            ism_surround(i_ism)%mydust%patch(ixO^S)=.not.patch_inuse(ixO^S)
            i_start= ism_surround(i_ism)%mydust%myconfig%idust_first
            i_end  = ism_surround(i_ism)%mydust%myconfig%idust_last
            Loop_ism_idustcld:  do i_dust=i_start,i_end
              ism_surround(i_ism)%mydust%the_ispecies(i_dust)%patch(ixO^S)=&
                                           .not.patch_inuse(ixO^S)
            end do Loop_ism_idustcld
            call ism_surround(i_ism)%mydust%set_w_zero(ixI^L,ixO^L,x,w)
          end if cond_ism_dust_oncld
        end do Loop_isms_cloud0
       end if ism_is_oncld

       jet_is_on : if(usrconfig%jet_yso_on)then
        Loop_jet_yso_clean : do i_jet_yso=0,usrconfig%jet_yso_number-1
         if(jet_yso(i_jet_yso)%myconfig%tracer_on)then

           where(patch_inuse(ixO^S))
             w(ixO^S,phys_ind%tracer(jet_yso(i_jet_yso)%myconfig%itr))=0.0_dp
           end where
         end if
        end do Loop_jet_yso_clean
       end if jet_is_on
       patch_all(ixO^S) =  patch_all(ixO^S) .and. .not.patch_inuse(ixO^S)
       if(cloud_medium(i_cloud)%myconfig%dust_on)then
         i_object_w_dust = i_object_w_dust +1
         if(.not.allocated(the_dust_inuse(i_object_w_dust)%patch))allocate(the_dust_inuse(i_object_w_dust)%patch(ixI^S))
         the_dust_inuse(i_object_w_dust)%myconfig     = cloud_medium(i_cloud)%mydust%myconfig
         if(.not.allocated(the_dust_inuse(i_object_w_dust)%the_ispecies))&
         allocate(the_dust_inuse(i_object_w_dust)%the_ispecies&
         (the_dust_inuse(i_object_w_dust)%myconfig%idust_first:&
          the_dust_inuse(i_object_w_dust)%myconfig%idust_last))
         call the_dust_inuse(i_object_w_dust)%set_patch(ixI^L,ixO^L, cloud_medium(i_cloud)%patch)
      !   the_dust_inuse(i_object_w_dust)%patch(ixO^S) = cloud_medium(i_cloud)%patch(ixO^S)
       end if
      end do Loop_clouds
    end if cloud_on





    if(any(patch_all(ixO^S)))then
      call usr_fill_empty_region(ixI^L,ixO^L,0.0_dp,patch_all,x,w)
    end if


  ! put dust to zero in all other zones
    cond_dust_on : if(phys_config%dust_on) then
      dust_dummy%myconfig%idust_first = 1
      dust_dummy%myconfig%idust_last  = phys_config%dust_n_species
      call dust_dummy%set_allpatch(ixI^L,ixO^L,the_dust_inuse)
      call dust_dummy%set_w_zero(ixI^L,ixO^L,x,w)
      call dust_dummy%clean_memory
      if(allocated(the_dust_inuse(i_object_w_dust)%patch))deallocate(the_dust_inuse(i_object_w_dust)%patch)
    end if cond_dust_on

    !Add grackle initialization here :

    ! check is if initial setting is correct
    call  phys_check_w(.true., ixI^L, ixO^L, w, flag)

    if(any(flag(ixO^S)>0)) PRINT*,' is error',maxval(flag(ixO^S)),minval(w(ixO^S,phys_ind%pressure_))






    ! get conserved variables to be used in the code
    patch_all=.true.
    call usr_check_w(ixI^L,ixO^L,patch_all,.true.,'initonegrid_usr',0.0_dp,x,w)
    call phys_to_conserved(ixI^L,ixO^L,w,x)

    call usr_clean_memory


  end subroutine initonegrid_usr
  !============================================================================
  subroutine specialsource_usr(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
    use mod_dust
    implicit none

    integer, intent(in)             :: ixI^L, ixO^L, iw^LIM
    real(dp), intent(in)            :: qdt, qtC, qt
    real(dp), intent(in)            :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    real(dp), intent(inout)         :: w(ixI^S,1:nw)
    ! .. local ..
    integer                         :: i_cloud,i_ism,i_jet_yso
    logical, dimension(ixI^S)       :: escape_patch
    real(kind=dp), dimension(ixI^S) :: source_filter,density_ratio,pressure_ism
    real(kind=dp), dimension(ixI^S) :: wjettracer,wismtracer
    real(kind=dp)                   :: usr_loc_tracer_small_density

    integer                         :: idir,idust,ix^D,ixL^D,ixR^L
    integer                         :: patch_back_cloud(ixI^S)
    !----------------------------------------------------------

    cond_reset : if(usrconfig%reset_medium) then
     escape_patch(ixI^S) =.false.
     cond_jet_on : if(usrconfig%jet_yso_on)then
       Loop_jet_yso : do i_jet_yso=0,usrconfig%jet_yso_number-1
         cond_jet_tracer_on : if(jet_yso(i_jet_yso)%myconfig%tracer_on)  then
          jet_yso(i_jet_yso)%subname='specialsource_usr'
          call jet_yso(i_jet_yso)%alloc_set_patch(ixI^L,ixI^L,qt,x,use_tracer=.false.,w=w)
          escape_patch(ixI^S) =  escape_patch(ixI^S).or.jet_yso(i_jet_yso)%patch(ixI^S)
          if(any(x(ixI^S,zjet_)>jet_yso(i_jet_yso)%myconfig%z_impos))then
            usr_loc_tracer_small_density = &
                jet_yso(i_jet_yso)%myconfig%tracer_small_density/10.0_dp
          else
            usr_loc_tracer_small_density = &
                jet_yso(i_jet_yso)%myconfig%tracer_small_density
          end if
          escape_patch(ixI^S) =  escape_patch(ixI^S).or.&
             w(ixI^S,phys_ind%tracer(jet_yso(i_jet_yso)%myconfig%itr))&
             >usr_loc_tracer_small_density
          escape_patch(ixI^S) =  escape_patch(ixI^S).or.&
           x(ixI^S,rjet_)<jet_yso(i_jet_yso)%myconfig%r_out_impos

          where(w(ixI^S,phys_ind%tracer(jet_yso(i_jet_yso)%myconfig%itr))<&
                    usr_loc_tracer_small_density)
            w(ixI^S,phys_ind%tracer(jet_yso(i_jet_yso)%myconfig%itr)) = 0.0_dp
          end where
          !escape_patch(ixO^S) = escape_patch(ixO^S).or. &
          !                      (x(ixO^S,x_)>1.1*jet_yso(i_jet_yso)%myconfig%r_in_init.and.&
          !                      x(ixO^S,x_)<1.1*jet_yso(i_jet_yso)%myconfig%r_out_impos)


        end if cond_jet_tracer_on

       end do Loop_jet_yso

     end if cond_jet_on


     ! set one cloud
     cond_usr_cloud : if(usrconfig%cloud_on)then
      Loop_clouds : do i_cloud=0,usrconfig%cloud_number-1
       cond_cloud_tracer : if(cloud_medium(i_cloud)%myconfig%tracer_on)  then
        cloud_medium(i_cloud)%subname='specialsource_usr'
        if(allocated(cloud_medium(i_cloud)%patch))deallocate(cloud_medium(i_cloud)%patch)
        call cloud_medium(i_cloud)%alloc_set_patch(ixI^L,ixI^L,qt,x,use_tracer=.true.,w=w)
        escape_patch(ixI^S) =  escape_patch(ixI^S).or.cloud_medium(i_cloud)%patch(ixI^S)
       end if cond_cloud_tracer
      end do Loop_clouds
     end if cond_usr_cloud


     ! add force in ISM
     cond_ism_on : if(usrconfig%ism_on) then
        if(usrconfig%jet_yso_on)then
!          Loop_isms0 : do i_ism=0,usrconfig%ism_number-1
!            where(w(ixO^S,phys_ind%rho_)<min(ism_surround(i_ism)%myconfig%density,&
!                                          minval(jet_yso(:)%myconfig%density))/10.0_dp)
!              escape_patch(ixO^S)=.false.
!            elsewhere(w(ixO^S,phys_ind%rho_)>ism_surround(i_ism)%myconfig%density.and.&
!              w(ixO^S,phys_ind%tracer(ism_surround(i_ism)%myconfig%itr))>5.0d3.and.&
!             x(ixO^S,zjet_)>2.0*maxval(jet_yso(:)%myconfig%z_impos))
!              escape_patch(ixO^S)=.true.
!            end where
!          end do Loop_isms0
        else
          escape_patch(ixI^S)=.false.
        end if

        cond_ism_noescape : if(.not.(all(escape_patch(ixI^S))))then
          Loop_isms : do i_ism=0,usrconfig%ism_number-1


            cond_reset_ism0 : if(ism_surround(i_ism)%myconfig%reset_on) then

              source_filter(ixI^S) = 0.5_dp*(1.0_dp-&
                tanh((x(ixI^S,zjet_)-ism_surround(i_ism)%myconfig%reset_distance(zjet_)/4.0)/&
                (2.0_dp*ism_surround(i_ism)%myconfig%reset_scale(zjet_)/4.0)))/2.0_dp


              cond_jeton_ism0: if(usrconfig%jet_yso_on)then
                Loop_jet_yso2 : do i_jet_yso=0,usrconfig%jet_yso_number-1
                  cond_jet_tracer_on2 : if(jet_yso(i_jet_yso)%myconfig%tracer_on)  then

                    where(escape_patch(ixI^S))
                      source_filter(ixI^S) = 0.0_dp
                    elsewhere
                     source_filter(ixI^S) = source_filter(ixI^S)*(1.0_dp-&
                            tanh((jet_yso(i_jet_yso)%myconfig%r_out_impos-x(ixI^S,rjet_))/&
                                       (0.1_dp*(jet_yso(i_jet_yso)%myconfig%r_out_impos))))

                     source_filter(ixI^S) = source_filter(ixI^S)&
                      *wCT(ixI^S,phys_ind%tracer(ism_surround(i_ism)%myconfig%itr))/&
                       merge(ism_surround(i_ism)%myconfig%tracer_init_density, &
                          wCT(ixI^S,phys_ind%rho_),                               &
                          ism_surround(i_ism)%myconfig%tracer_init_density>0.0_dp)

                    end where
                  end if cond_jet_tracer_on2
                end do Loop_jet_yso2
              end if cond_jeton_ism0
            end if cond_reset_ism0

            cond_tracer_ism_on : if((ism_surround(i_ism)%myconfig%reset_on &
                          .and.ism_surround(i_ism)%myconfig%tracer_on).or. &
                               ism_surround(i_ism)%myconfig%profile_force_on)then
              where(source_filter(ixI^S)<smalldouble)
                escape_patch(ixI^S)=.true.
              end where


              where(.not.escape_patch(ixI^S))&
              w(ixI^S,phys_ind%tracer(ism_surround(i_ism)%myconfig%itr))=&
               merge(ism_surround(i_ism)%myconfig%tracer_init_density, &
                 wCT(ixI^S,phys_ind%rho_),                               &
                 ism_surround(i_ism)%myconfig%tracer_init_density>0.0_dp)

              call ism_surround(i_ism)%add_source(ixI^L,ixO^L,iw^LIM,x,qdt,qtC,&
                                                  wCT,qt,w,use_tracer=.true.,&
                                                  escape_patch=escape_patch,&
                                                  source_filter=source_filter)
            end if cond_tracer_ism_on
          end do  Loop_isms
        end if cond_ism_noescape
      end if cond_ism_on

     call usr_clean_memory
    end if cond_reset


    cond_ismon_mup : if(usrconfig%ism_on) then
      cond_jeton_mup : if(usrconfig%jet_yso_on)then
        Loop_jet_yso_mup : do i_jet_yso=0,usrconfig%jet_yso_number-1
          Loop_isms_mup : do i_ism=0,usrconfig%ism_number-1

            cond_mup_on: if(phys_config%mean_mup_on.and.jet_yso(i_jet_yso)%myconfig%tracer_on&
                 .and. jet_yso(i_jet_yso)%myconfig%tracer_init_density>0.0_dp )then

              wjettracer(ixO^S)  = w(ixO^S,phys_ind%tracer(jet_yso(i_jet_yso)%myconfig%itr))/ &
                           jet_yso(i_jet_yso)%myconfig%tracer_init_density
              wismtracer(ixO^S)  = w(ixO^S,phys_ind%tracer(ism_surround(i_ism)%myconfig%itr))/ &
                           ism_surround(i_ism)%myconfig%tracer_init_density

              w(ixO^S,phys_ind%mup_)=(wjettracer(ixO^S)*jet_yso(i_jet_yso)%myconfig%mean_mup &
                + ism_surround(i_ism)%myconfig%mean_mup*wismtracer(ixO^S))/&
                (wjettracer(ixO^S)+wismtracer(ixO^S))

            end if cond_mup_on

          end do Loop_isms_mup
        end do Loop_jet_yso_mup
      end if cond_jeton_mup
    end if cond_ismon_mup

    cond_dust_on : if(phys_config%dust_on) then
    Loop_idust : do idust =1, phys_config%dust_n_species
      density_ratio(ixI^S)=w(ixI^S, phys_ind%dust_rho(idust))/w(ixI^S, phys_ind%rho_)

     Loop_idir2 : do idir = 1,ndir

     where(density_ratio(ixI^S)>smalldouble.and.&
          ((density_ratio(ixI^S)<usrconfig%density_dusttogas_minlimit&
           .or.density_ratio(ixI^S)>usrconfig%density_dusttogas_maxlimit).or.&
        (w(ixI^S, phys_ind%dust_mom(idir,idust))*w(ixI^S,phys_ind%mom(idir))<0.0_dp)))
       w(ixI^S, phys_ind%dust_mom(idir,idust)) = w(ixI^S,phys_ind%mom(idir))*&
                                                 density_ratio(ixI^S)



     end where
      where(density_ratio(ixI^S)>usrconfig%density_dusttogas_maxlimit)
        w(ixI^S, phys_ind%dust_rho(idust)) = &
          usrconfig%density_dusttogas_maxlimit*w(ixI^S, phys_ind%rho_)
        w(ixI^S, phys_ind%dust_mom(idir,idust)) = w(ixI^S,phys_ind%mom(idir))*&
                              usrconfig%density_dusttogas_maxlimit
      elsewhere(density_ratio(ixI^S)<usrconfig%density_dusttogas_minlimit**2.0_dp)
        w(ixI^S, phys_ind%dust_rho(idust)) = 0.0_dp
        w(ixI^S, phys_ind%dust_mom(idir,idust)) =  0.0_dp
      end where
     end do Loop_idir2
    end do Loop_idust
   end if cond_dust_on

!  call usr_check_w(ixI^L,ixO^L,.true.,'specialsource_usr',qt,x,w)
  end subroutine specialsource_usr
  !========================================================================
  subroutine process_grid_usr(igrid,level,ixI^L,ixO^L,qt,w,x)
    use mod_global_parameters
    implicit none
    integer, intent(in)             :: igrid,level,ixI^L,ixO^L
    real(kind=dp)   , intent(in)    :: qt
    real(kind=dp)   , intent(in)    :: x(ixI^S,1:ndim)
    real(kind=dp)   , intent(inout) :: w(ixI^S,1:nw)

    ! .. local ..
    integer                         :: idust,idir
    real(dp)                        :: small_dust_rho,coef
    real(dp), dimension(ixI^S)      :: temperature
    logical, dimension(ixI^S)       :: patch_correct,patch_slow
    integer       :: i_cloud,i_ism,i_jet_yso,i_dust,i_start,i_end
    logical, save :: first=.true.
    logical       :: patch_all(ixI^S)
    integer       :: i_object
    !---------------------------------------------------
    ! control the high temperature issues
    cond_no_isotherm : if(phys_config%energy)then
     cond_check_T_high : if(usrconfig%temperature_max<bigdouble)then
      call phys_get_temperature(ixI^L, ixI^L,w,x,Temperature)
      cond_Thigh : if(any(Temperature(ixI^S)>usrconfig%temperature_max))then
        call phys_to_primitive(ixI^L,ixI^L,w,x)
        where(Temperature(ixI^S)>usrconfig%temperature_max)
          w(ixI^S,phys_ind%pressure_) = w(ixI^S,phys_ind%pressure_)&
                                     *usrconfig%temperature_max/Temperature(ixI^S)
        end where
  !      call usr_check_w(ixI^L,ixO^L,.false.,'process_grid_usr',qt,x,w)
        call phys_to_conserved(ixI^L,ixI^L,w,x)
      end if cond_Thigh
     end if cond_check_T_high
    end if  cond_no_isotherm
return
     jet_is_on_var : if(usrconfig%jet_yso_on)then
       Loop_jet_yso_var : do i_jet_yso=0,usrconfig%jet_yso_number-1
         cond_var_on : if(jet_yso(i_jet_yso)%myconfig%variation_on)then
           call jet_yso(i_jet_yso)%process_grid(ixI^L,ixO^L,qt,x,w)
         end if cond_var_on
       end do Loop_jet_yso_var
     end if jet_is_on_var


    cond_usr_cloud : if(usrconfig%cloud_on)then
     Loop_clouds : do i_cloud=0,usrconfig%cloud_number-1
      cond_cloud_on : if(cloud_medium(i_cloud)%myconfig%time_cloud_on>0.0_dp.and.&
                            dabs(cloud_medium(i_cloud)%myconfig%time_cloud_on-qt)<smalldouble)  then

       cloud_medium(i_cloud)%subname='process_grid_usr'
       if(allocated(cloud_medium(i_cloud)%patch))deallocate(cloud_medium(i_cloud)%patch)
       allocate(cloud_medium(i_cloud)%patch(ixI^S))
       cloud_medium(i_cloud)%patch(ixO^S) = x(ixO^S,zjet_)>cloud_medium(i_cloud)%myconfig%extend(zjet_) +&
                                dtan(cloud_medium(i_cloud)%myconfig%extend(3)*&
                                     usr_physunit%myconfig%length*dpi/180.0_dp)*(x(ixO^S,1)-xprobmin1)
       call cloud_medium(i_cloud)%set_w(ixI^L,ixO^L,global_time,x,w)

       ism_is_on2 : if(usrconfig%ism_on)then
         if(ism_surround(0)%myconfig%tracer_on)then
           where(cloud_medium(i_cloud)%patch(ixO^S))
             w(ixO^S,phys_ind%tracer(ism_surround(0)%myconfig%itr))=0.0_dp
           end where
         end if
       end if ism_is_on2

       jet_is_on : if(usrconfig%jet_yso_on)then
        Loop_jet_yso_clean : do i_jet_yso=0,usrconfig%jet_yso_number-1
         if(jet_yso(i_jet_yso)%myconfig%tracer_on)then
           where(cloud_medium(i_cloud)%patch(ixO^S))
             w(ixO^S,phys_ind%tracer(jet_yso(0)%myconfig%itr))=0.0_dp
           end where
         end if
        end do Loop_jet_yso_clean
       end if jet_is_on

       patch_all(ixO^S) =  patch_all(ixO^S) .and. .not.cloud_medium(i_cloud)%patch(ixO^S)
       if(cloud_medium(i_cloud)%myconfig%dust_on)the_dust_inuse(i_object)=cloud_medium(i_cloud)%mydust
       i_object = i_object +1
     end if cond_cloud_on
    end do Loop_clouds
    end if cond_usr_cloud





    cond_dust_on : if(phys_config%dust_on) then
     small_dust_rho = 1.0d-4!ism_surround%mydust%myconfig%min_limit_rel

     call phys_to_primitive(ixI^L,ixI^L,w,x)
     ! handel small density dust
     Loop_idust : do idust =1, dust_n_species
      where(w(ixI^S, phys_ind%dust_rho(idust))<max(small_dust_rho*w(ixI^S,phys_ind%rho_),&
         ism_surround(0)%mydust%myconfig%min_limit_abs))
        w(ixI^S, phys_ind%dust_rho(idust))= 0.8* min(small_dust_rho*w(ixI^S,phys_ind%rho_),&
             ism_surround(0)%mydust%myconfig%min_limit_abs)
        patch_correct(ixI^S) = .true.
      elsewhere
        patch_correct(ixI^S) = .false.
      end where
      ! handel large density dust
      where(w(ixI^S,phys_ind%rho_)<0.9*ism_surround(0)%myconfig%density)
       where(w(ixI^S, phys_ind%dust_rho(idust))>ism_surround(0)%mydust%myconfig%max_limit_rel*w(ixI^S,phys_ind%rho_))
        w(ixI^S, phys_ind%dust_rho(idust))=0.8*ism_surround(0)%mydust%myconfig%max_limit_rel*w(ixI^S,phys_ind%rho_)
        patch_slow(ixI^S) = .true.
       elsewhere
        patch_slow(ixI^S) =.false.
       end where
      end where
      ! handel large cmax in rarefied region
      ! do idir = 1, ndim
      !   call phys_get_cmax(w,x,ixI^L,ixO^L,idir,cmax)
      !
      ! end do


    !  new_dvflag(ixI^S)=.true.
    !  new_dfflag(ixI^S)=.true.

    !  vt2(ixI^S) = 3.0d0*w(ixI^S,e_)/w(ixI^S, rho_)
      Loop_idir1 : do idir = 1,ndim
       where(patch_correct(ixI^S))
        w(ixI^S, phys_ind%dust_mom(idir,idust))=0.0_dp
       end where
       where(patch_slow(ixI^S))
               w(ixI^S, phys_ind%dust_mom(idir,idust))=w(ixI^S,phys_ind%mom(idir))
       end where
      end do   Loop_idir1
     end do Loop_idust
     call phys_to_conserved(ixI^L,ixI^L,w,x)
   end if cond_dust_on

  call usr_clean_memory


  end subroutine process_grid_usr
  !---------------------------------------------------------------------
  !-------------------------------------------------------------------------
  subroutine specialbound_usr(qt,ixI^L,ixO^L,iB,w,x)
    ! special boundary types, user defined
    integer, intent(in)     :: ixO^L, iB, ixI^L
    real(dp), intent(in)    :: qt, x(ixI^S,1:ndim)
    real(dp), intent(inout) :: w(ixI^S,1:nw)
    ! .. local ..
    integer                 :: flag(ixI^S)
    integer                 :: i_cloud,i_ism,i_jet_yso
    integer                 :: i_object_w_dust
    integer                 :: idims,iside,idims2,iside2,iw2
    integer                 :: i_start,i_end,i_dust
    logical                 :: patch_all(ixI^S),patch_inuse(ixI^S)
    real(dp)                :: cloud_profile(ixI^S,1:nw)
    real(dp),allocatable    :: wp(:^D&,:)
    logical                    :: isboundary,to_fix,bc_to_fix,some_unfixed
    character(len=30)          :: myboundary_cond
    logical,allocatable        :: mynotmixed_fixed_bound(:,:,:)
    !-------------------------------------
    ! * According to subroutine bc_phys in which it is called by mod_ghostcells_update.t twice:
    ! The input integers are ixI^L=ixG^L=ixG^LL=ixGlo1,ixGlo2,ixGhi1,ixGhi2
    ! and ixO^L=ixI^L defined from ixB^L in bc_phys as :
    ! > idims = 1, iside==2 (maximal boundary)
    ! ==> ixImin1=ixBmax1+1-nghostcells;ixImin2=ixBmin2;ixImax1=ixBmax1;ixImax2=ixBmax2;
    ! > idims = 1, iside==1 (minimal boundary)
    ! ==> ixImin1=ixBmin1;ixImin2=ixBmin2;ixImax1=ixBmin1-1+nghostcells;ixImax2=ixBmax2;
    ! > idims = 2, iside==2 (maximal boundary)
    ! ==> ixImin1=ixBmin1;ixImin2=ixBmax2+1-nghostcells;ixImax1=ixBmax1;ixImax2=ixBmax2;
    ! > idims = 2, iside==1 (minimal boundary)
    ! ==> ixImin1=ixBmin1;ixImin2=ixBmin2;ixImax1=ixBmax1;ixImax2=ixBmin2-1+nghostcells;
    ! to sum up, the range delimiting each boundary individually
    !TO MODIFY

    if(.not.allocated(mynotmixed_fixed_bound))&
     allocate(mynotmixed_fixed_bound(1:ndim,1:2,1:nwfluxbc))


    patch_all(ixO^S) = .true.
    i_object_w_dust = 1

    ! original code :
    idims = ceiling(real(iB,kind=dp)/2.0_dp)
    iside = iB-2*(idims-1)
    ! * iB=1=> idims = ceiling(real(1,dp)/2.0_dp) = ceiling(1.0/2.0) = ceiling(0.5)
    ! => idims = 1,iside = 1 - 2*(1-1) = 1
    ! * iB=2=> idims = ceiling(2.0/2.0) = ceiling(1.0)
    ! => idims = 1,iside = 2 - 2*(1-1) = 2
    ! * iB=3=> idims = ceiling(3.0/2.0) = ceiling(1.5)
    ! => idims = 2,iside = 3 - 2*(2-1) = 3 - 2 = 1
    ! * iB=4=> idims = ceiling(4.0/2.0) = ceiling(2.0)
    ! => idims = 2,iside = 4 - 2*(2-1) = 4 - 2 = 2
    ! * iB=5=> idims = ceiling(5.0/2.0) = ceiling(2.5)
    ! => idims = 3,iside = 5 - 2*(3-1) = 5 - 4 = 1
    ! * iB=6=> idims = ceiling(6.0/2.0) = ceiling(3.0)
    ! => idims = 3,iside = 6 - 2*(3-1) = 6 - 4 = 2



    !if(mype==0)then
    !  write(*,*) '==========User-defined special boundary conditions initialized=============='
    !end if
  ! set the ism
    cond_ism_on: if(usrconfig%ism_on)then
     Loop_isms : do i_ism=0,usrconfig%ism_number-1
      ism_surround(i_ism)%subname='specialbound_usr'
        call ism_surround(i_ism)%set_w(ixI^L,ixO^L,qt,x,w,isboundary_iB=(/idims,iside/))
        patch_all(ixO^S) =  patch_all(ixO^S) .and. .not.ism_surround(i_ism)%patch(ixO^S)

        !Here we set the user-defined boundary conditions already read from .par parameters files
        !call ism_surround(i_ism)%myboundaries%set_w(ixI^L,ixO^L,iB,idims,iside,&
        !                              ism_surround(i_ism)%patch,x,w)
        ! if does not work, maybe try with replacingism_surround(i_ism)%patch by .true. array everywhere
        !print*,'all(x(ixO^S,idims)<=xprobmin2)=',all(x(ixO^S,idims)<=xprobmin2)
        !print*,'boundary conditions: (idims,iside,iB)=',idims,iside,iB
        !print*, ism_surround(i_ism)%myboundaries%myconfig%boundary_type(idims,iside)

        !if(.not.allocated(ism_surround(i_ism)%patch)) then
         !allocate(ism_surround(i_ism)%patch(ixI^S))
         !ism_surround(i_ism)%patch              = .true.
        !end if

        myboundary_cond = ism_surround(i_ism)%myboundaries%myconfig%boundary_type(idims,iside)
        do idims2=1,ndim
          do iside2=1,2
            do iw2=1,nwfluxbc
              mynotmixed_fixed_bound(idims2,iside2,iw2)=.not.ism_surround(i_ism)%myboundaries%myconfig%mixed_fixed_bound(idims2,iside2,iw2)
            end do
          end do
        end do

        if(trim(myboundary_cond)/='fix'.or.any(mynotmixed_fixed_bound))then
          call ism_surround(i_ism)%myboundaries%set_w(ixI^L,ixO^L,iB,idims,iside,&
                                    ism_surround(i_ism)%patch,x,w)
        !else
        !  call ism_surround(i_ism)%set_w(ixI^L,ixO^L,qt,x,w,isboundary_iB=(/idims,iside/))
        end if
        !w(ixO^S,1:nw)=wp(ixO^S,1:nw)
        !deallocate(wp)

      !patch_all(ixO^S) =  patch_all(ixO^S) .and.(.not.ism_surround(i_ism)%patch(ixO^S))
      !if(ism_surround(i_ism)%myconfig%dust_on)the_dust_inuse(i_object_w_dust)=ism_surround(i_ism)%mydust
      i_object_w_dust = i_object_w_dust +1
     end do Loop_isms
   end if cond_ism_on
   !if(mype==0)then
   !   write(*,*) '=================================='
   !end if

  ! get conserved variables to be used in the code
  patch_all(ixO^S) = .true.
  call usr_check_w(ixI^L,ixO^L,patch_all,.false.,'specialbound_usr',qt,x,w)
  !call phys_to_conserved(ixI^L,ixO^L,w,x)

  call usr_clean_memory

  end subroutine specialbound_usr




  !-----------------------------------------------------------
  !> subroutine to check w
  subroutine usr_check_w(ixI^L,ixO^L,patch_check,is_conserved,subname,qt,x,w)
   !use mod_radiative_cooling, only : coolconfig
   implicit none
   integer, intent(in)             :: ixO^L,  ixI^L
   logical, intent(in)             :: patch_check(ixI^S)
   logical, intent(in)             :: is_conserved
   character(len=*), intent(in)    :: subname
   real(dp), intent(in)            :: qt
   real(dp), intent(in)            :: x(ixI^S,1:ndim)
   real(dp), intent(inout)         :: w(ixI^S,1:nw)
   ! .. local ..
   real(kind=dp), dimension(ixI^S) :: pressure, Temperature
   !------------------------------
   if(all(.not.patch_check(ixO^S)))return

   cond_cons : if(is_conserved) then

  !  pressure(ixO^S) = w(ixO^S,phys_ind%pressure_)/(phys_config%gamma-1.0_dp)
    call phys_get_pthermal(w, x, ixI^L, ixO^L, pressure)
    call phys_get_temperature(ixI^L, ixO^L,w,x,Temperature)
   if(any(pressure(ixO^S)<phys_config%small_pressure&
          .and. patch_check(ixO^S))) then
     write(*,*) 'The pressure set at usr procedure is smaller than small_pressure'
     write(*,*) 'is called from the subroutine : ', subname
     write(*,*) 'is conserved value'
     write(*,*) 'The minimum of pressure set is : ', minval(pressure(ixO^S))
     write(*,*) 'The small_pressure is : ', phys_config%small_pressure

     call mpistop ('The code stops at usr_check_w in mod_usr_yso_jet.t')
    end if
    if(any(w(ixO^S,phys_ind%rho_)<phys_config%small_density&
          .and. patch_check(ixO^S))) then
     write(*,*) 'The density set at usr procedure is smaller than small_densitiy'
     write(*,*) 'is called from the subroutine : ', subname
     write(*,*) 'is conserved value'
     write(*,*) 'The minimum of density set is : ', minval(w(ixO^S,phys_ind%rho_))
     write(*,*) 'The small_density is : ', phys_config%small_density
     call mpistop ('The code stops at usr_check_w in mod_usr_yso_jet.t')
    end if
    if(phys_config%radiative_cooling) then
      if(any(Temperature(ixO^S)>phys_config%cool_tlow&
            .and. patch_check(ixO^S))) then
       write(*,*) 'The temperature set at usr procedure is bigger than cooling temperature'
       write(*,*) 'The minimum of temperature set is : ', maxval(Temperature(ixO^S))
       write(*,*) 'The cool_tlow is : ', phys_config%cool_tlow
       write(*,*) 'is called from the subroutine : ', subname
       !write(*,*) 'is Tlow = ',phys_config%cool_tlow*unit_temperature, ' K ', &
         !' and the cells The temperature ', &
         !Temperature(ixO^S)*unit_temperature, ' K'
       write(*,*) 'is conserved value'
       !call mpistop ('The code stops at usr_check_w in mod_usr_yso_jet.t')
      end if
    end if
   else cond_cons
    if(phys_config%energy) then
      if(any(w(ixO^S,phys_ind%pressure_)<phys_config%small_pressure&
            .and. patch_check(ixO^S))) then
        write(*,*) 'The pressure set at usr procedure is smaller than small_pressure'
        write(*,*) 'is called from the subroutine : ', subname
        write(*,*) 'The minimum of pressure set is : ', minval(w(ixO^S,phys_ind%pressure_))
        write(*,*) 'The small_pressure is : ', phys_config%small_pressure
        call mpistop ('The code stops at usr_check_w in mod_usr_yso_jet.t')
      end if
    end if
    if(any(w(ixO^S,phys_ind%rho_)<phys_config%small_density&
          .and. patch_check(ixO^S))) then
     write(*,*) 'The density set at usr procedure is smaller than small_density'
     write(*,*) 'is called from the subroutine : ', subname
     write(*,*) 'The minimum of density set is : ', minval(w(ixO^S,phys_ind%rho_))
     write(*,*) 'The small_density is : ', phys_config%small_density
     call mpistop ('The code stops at usr_check_w in mod_usr_yso_jet.t')
    end if

   end if cond_cons
  end subroutine usr_check_w

!----------------------------------------------------------------
  subroutine usr_clean_memory
    implicit none
    ! .. local ..
    integer   :: i_cloud,i_ism,i_jet_yso
    !------------------------------
        if(usrconfig%ism_on)then
          Loop_isms : do i_ism=0,usrconfig%ism_number-1
           call ism_surround(i_ism)%clean_memory
          end do Loop_isms
        end if
        if(usrconfig%cloud_on)then
          Loop_clouds : do i_cloud=0,usrconfig%cloud_number-1
           call cloud_medium(i_cloud)%clean_memory
          end do Loop_clouds
        end if

        if(usrconfig%jet_yso_on)then
          Loop_jet_yso : do i_jet_yso=0,usrconfig%jet_yso_number-1
           jet_yso(i_jet_yso)%subname='usr_clean_memory'
           call jet_yso(i_jet_yso)%clean_memory
         end do Loop_jet_yso
        end if


  end subroutine usr_clean_memory


     !> Enforce additional refinement or coarsening
     !> One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.
     !> you must set consistent values for integers refine/coarsen:
     !> refine = -1 enforce to not refine
     !> refine =  0 doesn't enforce anything
     !> refine =  1 enforce refinement
     !> coarsen = -1 enforce to not coarsen
     !> coarsen =  0 doesn't enforce anything
     !> coarsen =  1 enforce coarsen
     !> e.g. refine for negative first coordinate x < 0 as
     !> if (any(x(ix^S,1) < zero)) refine=1
     subroutine specialrefine_usr(igrid,level,ixI^L,ixO^L,qt,w,x,refine,coarsen)
       use mod_global_parameters
       integer, intent(in)          :: igrid, level, ixI^L, ixO^L
       real(dp), intent(in)         :: qt, w(ixI^S,1:nw), x(ixI^S,1:ndim)
       integer, intent(inout)       :: refine, coarsen
       ! .. local ..
       integer                      :: level_min,level_max,level_need
       logical                      :: patch_cond
       real(dp)                     :: dx_loc(1:ndim)
       integer                      :: i_jet_yso
      !----------------------------------------
      cond_ysojet_on : if(usrconfig%jet_yso_on)then
        Loop_jet_yso : do i_jet_yso=0,usrconfig%jet_yso_number-1
          cond_init_jet: if(dabs(qt-jet_yso(i_jet_yso)%myconfig%time_cla_jet_on)&
                           <smalldouble) then
            jet_yso(i_jet_yso)%subname='specialrefine_usr'
            ^D&dx_loc(^D)=4.0_dp*rnode(rpdx^D_,igrid);
            call jet_yso(i_jet_yso)%set_patch(ixI^L,ixO^L,qt,x,&
                                              force_refine=1,dx_loc=dx_loc)

            if(any(jet_yso(i_jet_yso)%patch(ixO^S)))then
             level_need= nint(dlog((dabs(jet_yso(i_jet_yso)%myconfig%r_out_init&
                                      -jet_yso(i_jet_yso)%myconfig%r_in_init))/&
                               domain_nx1)/dlog(2.0_dp))
             level_min = max(jet_yso(i_jet_yso)%myconfig%refine_min_level&
                             ,level_need)
             level_max = jet_yso(i_jet_yso)%myconfig%refine_max_level

             else
             level_min = 1
             level_max = max(jet_yso(i_jet_yso)%myconfig%refine_max_level - 2,level_min)
            end if
            patch_cond=.true.
            call user_fixrefineregion(level,level_min,level_max,&
                                     patch_cond,refine,coarsen)
            call jet_yso(i_jet_yso)%clean_memory
          end if cond_init_jet
        end do Loop_jet_yso
      end if  cond_ysojet_on



     end subroutine specialrefine_usr
  !=====================================================================
     subroutine user_fixrefineregion(level,level_min,level_max,patch_cond,refine,coarsen)
     integer, intent(in)    :: level,level_min,level_max
     logical,intent(in)     :: patch_cond
     integer, intent(inout) :: refine, coarsen
     ! .. local ..
     !-------------------------------------------
     if(patch_cond)then
      if(level>level_max)then
        coarsen = 1
        refine  = -1
      else if(level==level_max)then
        coarsen = 0
        refine  = -1
      end if
      if(level<level_min)then
        coarsen = -1
        refine  =  1
      else if(level==level_min)then
        coarsen = -1
        refine  =  0
      end if
     end if
     end subroutine user_fixrefineregion

     subroutine special_get_dt(w,ixI^L,ixO^L,qt,dtnew,dx^D,x)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L
       double precision, intent(in)    :: dx^D,qt, x(ixI^S,1:ndim)
       double precision, intent(in)    :: w(ixI^S,1:nw)
       double precision, intent(inout) :: dtnew
       ! .. local ...
       integer                         :: i_cloud
       !--------------------------------------------------------------
        cond_usr_cloud : if(usrconfig%cloud_on)then
          Loop_clouds : do i_cloud=0,usrconfig%cloud_number-1
            cond_cloud_on : if(cloud_medium(i_cloud)%myconfig%time_cloud_on>0.0_dp.and.&
                              cloud_medium(i_cloud)%myconfig%time_cloud_on>qt)  then
              dtnew= cloud_medium(i_cloud)%myconfig%time_cloud_on-qt
            end if   cond_cloud_on
          end do Loop_clouds
        end if cond_usr_cloud

     end subroutine special_get_dt
   !> This subroutine can be used to artificially overwrite ALL conservative
   !> variables in a user-selected region of the mesh, and thereby act as
   !> an internal boundary region. It is called just before external (ghost cell)
   !> boundary regions will be set by the BC selection. Here, you could e.g.
   !> want to introduce an extra variable (nwextra, to be distinguished from nwaux)
   !> which can be used to identify the internal boundary region location.
   !> Its effect should always be local as it acts on the mesh.
   subroutine usr_special_internal_bc(level,qt,ixI^L,ixO^L,w,x)
     use mod_global_parameters
     integer, intent(in)             :: ixI^L,ixO^L,level
     real(kind=dp)   , intent(in)    :: qt
     real(kind=dp)   , intent(inout) :: w(ixI^S,1:nw)
     real(kind=dp)   , intent(in)    :: x(ixI^S,1:ndim)
     ! .. local ..
     integer                         :: i_jet_yso
     !--------------------------------------------------------------

    cond_jet_on : if(usrconfig%jet_yso_on)then

      Loop_jet_yso : do i_jet_yso=0,usrconfig%jet_yso_number-1
        cond_jet_start : if(qt>jet_yso(i_jet_yso)%myconfig%time_cla_jet_on)then
         cond_jet_impose : if(jet_yso(i_jet_yso)%myconfig%z_impos>box_limit(1,zjet_))then
          call jet_yso(i_jet_yso)%set_patch(ixI^L,ixO^L,qt,x)
          cond_insid_jet : if(any(jet_yso(i_jet_yso)%patch(ixI^S)))then
            call phys_to_primitive(ixI^L,ixO^L,w,x)
            call jet_yso(i_jet_yso)%set_w(ixI^L,ixO^L,qt,x,w)
            ! get conserved variables to be used in the code
            call phys_to_conserved(ixI^L,ixO^L,w,x)
          end if cond_insid_jet
          call usr_clean_memory
         end if cond_jet_impose
        end if cond_jet_start
      end do Loop_jet_yso

    end if cond_jet_on

   end subroutine usr_special_internal_bc
  !> special output
  subroutine specialvar_output(ixI^L,ixO^L,win,x,normconv)
  ! this subroutine can be used in convert, to add auxiliary variables to the
  ! converted output file, for further analysis using tecplot, paraview, ....
  ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
  ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
  ! corresponding normalization values (default value 1)
  !/
    use mod_physics
    use mod_dust
    implicit none
    integer, intent(in)        :: ixI^L,ixO^L
    real(dp), intent(in)       :: x(ixI^S,1:ndim)
    real(dp)                   :: win(ixI^S,nw+nwauxio)
    real(dp)                   :: normconv(0:nw+nwauxio)
    ! .. local ..
    real(dp)                   :: w(ixI^S,nw),w_init(ixI^S,nw)
    real(dp)                   :: error_var(ixM^T)
    integer                    :: iw, level,idir,idust
    integer, parameter         :: icsound_             = 1
    integer, parameter         :: itemperature_      = 2
    integer, parameter         :: ilevel_            = 3
    integer, parameter         :: indensity_         = 4
    integer, parameter         :: igravfield1         = 5
    integer, parameter         :: igravfield2         = 6
    integer, parameter         :: igravfield3         = 7
    integer, parameter         :: igravphi_         = 8
    integer, parameter         :: ierror_lohner_rho_ = 9
    integer, parameter         :: ierror_lohner_p_   = 10
    integer, parameter         :: ideltav_dust11_    = 11
    integer, parameter         :: ideltav_dust12_    = 12
    integer, parameter         :: ipos_rCC_           = 13
    integer, parameter         :: ipos_zCC_           = 14
    integer, parameter         :: irelerr_rho_           = 15
    integer, parameter         :: irelerr_p_           = 16
    integer, parameter         :: irelerr_v1_          = 17
    integer, parameter         :: irelerr_v2_          = 18
    integer, parameter         :: irelerr_v3_          = 19
    integer, parameter         :: irhod1torho_       = 20

    real(dp),dimension(ixI^S)                                  :: csound2,temperature
    real(dp),dimension(ixI^S,ndir)                             :: vgas
    real(dp),dimension(ixI^S,ndir,phys_config%dust_n_species)  :: vdust
    double precision :: gravity_field(ixI^S,ndim)
    double precision :: gravity_potential(ixI^S)


    !----------------------------------------------------


    w(ixI^S,1:nw) = win(ixI^S,1:nw)
    level = node(plevel_,saveigrid)

    call phys_handle_small_values(.false., w, x, ixI^L, ixO^L,'specialvar_output')
    Loop_iw :  do iw = 1,nwauxio
    select case(iw)
    case(icsound_)
      call phys_get_csound2(w,x,ixI^L,ixO^L,csound2)
      win(ixO^S,nw+icsound_) = csound2(ixO^S)*usr_physunit%myconfig%velocity**2.0_dp
          !dsqrt(SUM(w(ixO^S,phys_ind%mom(1):phys_ind%mom(ndir))**2.0_dp&
                                !,dim=ndim+1)/csound2(ixO^S))&
                                !/(w(ixO^S,phys_ind%rho_))
      !/
      case(itemperature_)
        call phys_get_temperature( ixI^L, ixO^L,w, x, temperature)
        win(ixO^S,nw+itemperature_) = temperature(ixO^S)*unit_temperature
        normconv(nw+itemperature_)  = 1.0_dp
      case(ilevel_)
        normconv(nw+ilevel_)     = 1.0_dp
        win(ixO^S,nw+ilevel_) = node(plevel_,saveigrid)

      case(indensity_)
        normconv(nw+indensity_)     = 1.0_dp
        win(ixO^S,nw+indensity_)    = w(ixO^S,phys_ind%rho_)*unit_density/(phys_config%mean_mass*mp_cgs)
      case(igravfield1,igravfield2,igravfield3)
        if (.not. associated(usr_gravity)) then
          write(*,*) "mod_usr.t: please point usr_gravity to a subroutine"
          write(*,*) "like the phys_gravity in mod_usr_methods.t"
          call mpistop("gravity_add_source: usr_gravity not defined")
        else
          call usr_gravity(ixI^L,ixO^L,w,x,gravity_field)
        end if

        normconv(nw+igravfield1)     = 1.0_dp
        win(ixO^S,nw+igravfield1)    = 0.0_dp
        normconv(nw+igravfield2)     = 1.0_dp
        win(ixO^S,nw+igravfield2)    = 0.0_dp
        normconv(nw+igravfield3)     = 1.0_dp
        win(ixO^S,nw+igravfield3)    = 0.0_dp
        {^D&win(ixO^S,nw+igravfield^D)    = gravity_field(ixO^S,^D)*&
        usr_physunit%myconfig%mass/&
        (unit_density*(unit_length/unit_velocity)**(2.0_dp)*&
        (usr_physunit%myconfig%length*usr_physunit%myconfig%length))
        }
      case(igravphi_)
        if (.not. associated(usr_gravity_potential)) then
          write(*,*) "mod_usr.t: please point usr_gravity to a subroutine"
          write(*,*) "like the phys_gravity in mod_usr_methods.t"
          call mpistop("gravity_add_source: usr_gravity not defined")
        else
          call usr_gravity_potential(ixI^L,ixO^L,w,x,gravity_potential)
        end if
        normconv(nw+igravphi_)     = 1.0_dp
        win(ixO^S,nw+igravphi_)    = 0.0_dp
        win(ixO^S,nw+igravphi_)    = gravity_potential(ixO^S) ! same thing than w(ixO^S,phys_ind%grav_phi_)
        !*&
        !usr_physunit%myconfig%length*&
        !usr_physunit%myconfig%mass/&
        !(unit_density*(unit_length/unit_velocity)**(2.0_dp)*&
        !(usr_physunit%myconfig%length*usr_physunit%myconfig%length))

      case(ipos_rCC_)
        normconv(nw+ipos_rCC_)     = 1.0_dp
        win(ixO^S,nw+ipos_rCC_)    = x(ixO^S,r_)*unit_length
      case(ipos_zCC_)
        normconv(nw+ipos_zCC_)     = 1.0_dp
        win(ixO^S,nw+ipos_zCC_)    = x(ixO^S,z_)*unit_length
      case(ierror_lohner_rho_)
        normconv(nw+iw)     = 1.0_dp
        win(ixG^T,nw+ierror_lohner_rho_) = 0.0
        call usr_mat_get_Lohner_error(ixI^L, ixM^LL,level,phys_ind%rho_,w,error_var)
        win(ixM^T,nw+ierror_lohner_rho_) = error_var(ixM^T)
      case(ierror_lohner_p_)
        normconv(nw+iw)     = 1.0_dp
        win(ixG^T,nw+ierror_lohner_p_) = 0.0_dp
        call usr_mat_get_Lohner_error(ixI^L, ixM^LL,level,phys_ind%pressure_,w,error_var)
        win(ixM^T,nw+ierror_lohner_p_) = error_var(ixM^T)

      case(ideltav_dust11_)
        normconv(nw+iw)     = 1.0_dp
        dust_on_deltav11 : if(phys_config%dust_on)then
            idir=1
            vgas(ixI^S,idir)=w(ixI^S,phys_ind%mom(idir))/w(ixI^S,phys_ind%rho_)
            idust = 1
            where(w(ixI^S,phys_ind%dust_rho(idust))>phys_config%dust_small_density)
              vdust(ixI^S,idir,idust)=w(ixI^S,phys_ind%dust_mom(idir, idust))/&
                                  w(ixI^S,phys_ind%dust_rho(idust))
              win(ixI^S,nw+ideltav_dust11_) = vgas(ixI^S,idir)  - vdust(ixI^S,idir,idust)
            elsewhere
              vdust(ixI^S,idir,idust)       = 0.0_dp
              win(ixI^S,nw+ideltav_dust11_) = 0.0_dp
            endwhere
        end if  dust_on_deltav11
        !/
      case(ideltav_dust12_)
        normconv(nw+iw)     = 1.0_dp
        dust_on_deltav12 : if(phys_config%dust_on)then
            idir=2
            vgas(ixI^S,idir)=w(ixI^S,phys_ind%mom(idir))/w(ixI^S,phys_ind%rho_)
            idust = 1
            where(w(ixI^S,phys_ind%dust_rho(idust))>phys_config%dust_small_density)
              vdust(ixI^S,idir,idust)=w(ixI^S,phys_ind%dust_mom(idir, idust))/&
                                  w(ixI^S,phys_ind%dust_rho(idust))
              win(ixI^S,nw+ideltav_dust12_) = vgas(ixI^S,idir)  - vdust(ixI^S,idir,idust)

            elsewhere
              vdust(ixI^S,idir,idust)       = 0.0_dp
              win(ixI^S,nw+ideltav_dust12_) = 0.0_dp
            endwhere
        end if  dust_on_deltav12
        !/
      case(irhod1torho_)
        normconv(nw+iw)     = 1.0_dp
        idust = 1
        where(w(ixI^S,phys_ind%dust_rho(idust))>phys_config%dust_small_density)
          win(ixI^S,nw+irhod1torho_) = w(ixI^S,phys_ind%dust_rho(idust)) /&
                                          w(ixI^S,phys_ind%rho_)
        elsewhere
          win(ixI^S,nw+irhod1torho_) = 0.0_dp
        end where

      case(irelerr_rho_,irelerr_p_,irelerr_v1_,irelerr_v2_,irelerr_v3_)
        call initonegrid_usr(ixI^L,ixO^L,w_init,x)
        select case(iw)
        case(irelerr_rho_)
          normconv(nw+irelerr_rho_)     = 1.0_dp
          win(ixO^S,nw+irelerr_rho_)    = 100.0_dp*((w(ixO^S,phys_ind%rho_)-&
                                          w_init(ixO^S,phys_ind%rho_))/&
                                          w_init(ixO^S,phys_ind%rho_))
        case(irelerr_p_)
          normconv(nw+irelerr_p_)     = 1.0_dp
          win(ixO^S,nw+irelerr_p_)    = 100.0_dp*((w(ixO^S,phys_ind%pressure_)-&
                                          w_init(ixO^S,phys_ind%pressure_))/&
                                          w_init(ixO^S,phys_ind%pressure_))
        case(irelerr_v1_)
          normconv(nw+irelerr_v1_)     = 1.0_dp
          win(ixO^S,nw+irelerr_v1_)  = 100.0_dp*((w(ixO^S,phys_ind%mom(1))-&
                                          w_init(ixO^S,phys_ind%mom(1)))/&
                                          w_init(ixO^S,phys_ind%mom(1)))
        case(irelerr_v2_)
          normconv(nw+irelerr_v2_)     = 1.0_dp
          if(ndir>=2)then
            win(ixO^S,nw+irelerr_v2_)  = 100.0_dp*((w(ixO^S,phys_ind%mom(2))-&
                                            w_init(ixO^S,phys_ind%mom(2)))/&
                                            w_init(ixO^S,phys_ind%mom(2)))
          else
            win(ixO^S,nw+irelerr_v2_)  = 0.0_dp
          end if
        case(irelerr_v3_)
          normconv(nw+irelerr_v3_)     = 1.0_dp
          if(ndir>=3)then
            win(ixO^S,nw+irelerr_v3_)  = 100.0_dp*((w(ixO^S,phys_ind%mom(3))-&
                                            w_init(ixO^S,phys_ind%mom(3)))/&
                                            w_init(ixO^S,phys_ind%mom(3)))
          else
            win(ixO^S,nw+irelerr_v3_)  = 0.0_dp
          end if
        end select
      case default
       write(*,*)'is not implimented at specialvar_output in mod_user'
       call mpistop('the code stops!')
    end select
  end do Loop_iw
    !----------------------------------------------------

   !w(ixI^S,1:nw)=win(ixI^S,1:nw)


  end subroutine specialvar_output


  !> this subroutine is ONLY to be used for computing auxiliary variables
  !> which happen to be non-local (like div v), and are in no way used for
  !> flux computations. As auxiliaries, they are also not advanced

  subroutine specialvarnames_output(varnames)
  ! newly added variables need to be concatenated with the w_names/primnames string
    character(len=*), intent(inout) :: varnames(:)
    ! .. local ..
    integer                    :: iw
    integer, parameter         :: icsound_             = 1
    integer, parameter         :: itemperature_      = 2
    integer, parameter         :: ilevel_            = 3
    integer, parameter         :: indensity_         = 4
    integer, parameter         :: igravfield1         = 5
    integer, parameter         :: igravfield2         = 6
    integer, parameter         :: igravfield3         = 7
    integer, parameter         :: igravphi_         = 8
    integer, parameter         :: ierror_lohner_rho_ = 9
    integer, parameter         :: ierror_lohner_p_   = 10
    integer, parameter         :: ideltav_dust11_    = 11
    integer, parameter         :: ideltav_dust12_    = 12
    integer, parameter         :: ipos_rCC_           = 13
    integer, parameter         :: ipos_zCC_           = 14
    integer, parameter         :: irelerr_rho_           = 15
    integer, parameter         :: irelerr_p_           = 16
    integer, parameter         :: irelerr_v1_          = 17
    integer, parameter         :: irelerr_v2_          = 18
    integer, parameter         :: irelerr_v3_          = 19
    integer, parameter         :: irhod1torho_       = 20
    character(len=30)          :: strint
    !----------------------------------------------------

    Loop_iw : do  iw = 1,nwauxio
      select case(iw)
      case(icsound_)
        varnames(icsound_)               = 'csound2'
      case(itemperature_)
        varnames(itemperature_)        ='temperature'
      case(ilevel_)
        varnames(ilevel_)              = 'level'
      case(indensity_)
        varnames(indensity_)           = 'number density'
      case(igravfield1,igravfield2,igravfield3)
        varnames(igravfield1)           = 'gravfield1'
        varnames(igravfield2)           = 'gravfield2'
        varnames(igravfield3)           = 'gravfield3'
      case(igravphi_)
        varnames(igravphi_)           = 'grav_phi'
      case(ierror_lohner_rho_)
        varnames(ierror_lohner_rho_)   = 'erroramrrho'
      case(ierror_lohner_p_)
        varnames(ierror_lohner_p_)     = 'erroamrp'
      case(ideltav_dust11_)
        varnames(ideltav_dust11_)      = 'deltav_dust11'
      case(ideltav_dust12_)
        varnames(ideltav_dust12_)      = 'deltav_dust12'
      case(irhod1torho_)
        varnames(irhod1torho_)         = 'rhod1torho'
      case(ipos_rCC_)
        varnames(ipos_rCC_)         = 'r_cell'
      case(ipos_zCC_)
        varnames(ipos_zCC_)         = 'z_cell'
      case(irelerr_rho_)
        varnames(irelerr_rho_)   = 'relerrorrho'
      case(irelerr_p_)
        varnames(irelerr_p_)   = 'relerrorp'
      case(irelerr_v1_)
        varnames(irelerr_v1_)   = 'relerrorv1'
      case(irelerr_v2_)
        varnames(irelerr_v2_)   = 'relerrorv2'
      case(irelerr_v3_)
        varnames(irelerr_v3_)   = 'relerrorv3'
      end select
    end do Loop_iw



  end subroutine specialvarnames_output
  !---------------------------------------------------------------------
  !> subroutine to fill the space regions that are not filled by the model
  subroutine usr_fill_empty_region(ixI^L,ixO^L,qt,patchw_empty,x,w)
    use mod_global_parameters
    implicit none
    integer, intent(in)         :: ixI^L,ixO^L
    real(kind=dp), intent(in)   :: qt
    logical, intent(in)         :: patchw_empty(ixI^S)
    real(kind=dp),intent(in)    :: x(ixI^S,1:ndir)
    real(kind=dp),intent(inout) :: w(ixI^S,1:nw)
    ! .. local ..
    integer                     :: idir
    !------------------------------------------------
    where(patchw_empty(ixO^S))
      w(ixO^S,phys_ind%rho_)      = 1.0_DP
    end where
    cond_no_isotherm : if(phys_config%energy)then
      where(patchw_empty(ixO^S))
        w(ixO^S,phys_ind%pressure_)= 1.0d-2
      end where
    end if cond_no_isotherm
    Loop_idir_v : do idir=1,ndir
     where(patchw_empty(ixO^S))
      w(ixO^S,phys_ind%mom(idir)) = 0.0_dp
     end where
     if(phys_config%ismhd) then
        where(patchw_empty(ixO^S))
          w(ixO^S,phys_ind%mag(idir)) = 0.0_dp
        end where
      end if
    end do Loop_idir_v
    if(phys_config%mean_mup_on) then
     where(patchw_empty(ixO^S))
        w(ixO^S,phys_ind%mup_) = phys_config%mean_mup
     end where
    end if
  end subroutine usr_fill_empty_region

!---------------------------------------------------------------------
  !> Initialize the method and limiter
  subroutine special_reset_solver(ixI^L,ixO^L,idims,qt,w,x,old_method,old_limiter,usr_method,usr_limiter)
    use mod_global_parameters
    use mod_limiter
    use mod_radiative_cooling
    !integer,parameter               :: limiter_config =12
    character(len=75)   :: filename_config
    !logical                         :: exist,bigtest,midtest
    integer, intent(in)             :: ixI^L, ixO^L,idims
    real(kind=dp), intent(in)       :: qt
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(in)    :: w(ixI^S,1:nw)
    character(len=*),intent(in)     :: old_method
    integer, intent(in)             :: old_limiter
    character(len=*),intent(inout)  :: usr_method
    integer, intent(inout)          :: usr_limiter
    ! .. local ..
    real(kind=dp), dimension(ixO^S) :: theta!,maxL1tab
    real(kind=dp)                   :: theta_min!,moddt,maxL1
    !------------------------------------------------------
    if(phys_config%radiative_cooling)then
      !filename_config=trim(base_filename)//'.limiter'
      !PRINT*,qt*usr_physunit%myconfig%time,dtsave(2)*usr_physunit%myconfig%time
      !moddt = MODULO(qt*usr_physunit%myconfig%time,dtsave(2)*usr_physunit%myconfig%time)
      !print*,moddt
      !if(moddt==0.0_dp)then
      !bigtest = any(w(ixI^S,phys_ind%Lcool1_)>usrconfig%reset_flux_scheme_thresholdL1_min)
      !midtest=any((w(ixI^S,phys_ind%Lcool1_)>0.0_dp).and.(w(ixI^S,phys_ind%Lcool1_)<usrconfig%reset_flux_scheme_thresholdL1_max))
      !bigtest=bigtest.and.midtest
      !if(bigtest)then
        !where(w(ixI^S,phys_ind%Lcool1_)<usrconfig%reset_flux_scheme_thresholdL1_max)
        !maxL1tab(ixO^S) = w(ixI^S,phys_ind%Lcool1_)
        !elsewhere
        !maxL1tab(ixO^S) = 0.0_dp
        !end where
        !maxL1 = MAXVAL(w(ixI^S,phys_ind%Lcool1_))!maxL1tab(ixO^S))
        !inquire(file=filename_config, exist=exist)
        !if (exist) then
        !  open(12, file=filename_config, status="old", position="append", action="write")
        !else
          !open(12, file=filename_config, status="new", action="write")
          !write(limiter_config,*)'reset_flux_scheme_thresholdL1_max unit = ',usr_physunit%myconfig%luminosity
        !end if
        !write(limiter_config,*)'============================================================'
        !write(limiter_config,*)'at time t(yrs) = ', qt*usr_physunit%myconfig%time/3.1557600e7
        !write(limiter_config,*)'reset_flux_scheme_thresholdL1_max (without unit) = ',usrconfig%reset_flux_scheme_thresholdL1_max
        !write(limiter_config,*)' reset_flux_scheme_thresholdL1_max (with unit) = ',usrconfig%reset_flux_scheme_thresholdL1_max*usr_physunit%myconfig%luminosity
        !write(limiter_config,*)'reset_flux_scheme_thresholdL1_min (without unit) = ',usrconfig%reset_flux_scheme_thresholdL1_min
        !write(limiter_config,*)' reset_flux_scheme_thresholdL1_min (with unit) = ',usrconfig%reset_flux_scheme_thresholdL1_min*usr_physunit%myconfig%luminosity
        !write(limiter_config,*)'L1<reset_flux_scheme_thresholL1_max, L1>reset_flux_scheme_thresholL1_min : ', &
        !any(w(ixI^S,phys_ind%Lcool1_)<usrconfig%reset_flux_scheme_thresholdL1_max), &
        !any(w(ixI^S,phys_ind%Lcool1_)>usrconfig%reset_flux_scheme_thresholdL1_min)
        !write(limiter_config,*)' min, max : ', &
        !MINVAL(w(ixI^S,phys_ind%Lcool1_))*usr_physunit%myconfig%luminosity, &
        !maxL1*usr_physunit%myconfig%luminosity
        !write(limiter_config,*)w(ixI^S,phys_ind%Lcool1_)*usr_physunit%myconfig%luminosity
        !close(limiter_config)
        !PRINT*,'at time t(yrs) = ', qt*usr_physunit%myconfig%time/3.1557600e7
        !PRINT*,' reset_flux_scheme_thresholdL1_max unit = ',usr_physunit%myconfig%luminosity
        !PRINT*,'reset_flux_scheme_thresholdL1_max (without unit) = ',usrconfig%reset_flux_scheme_thresholdL1_max
        !PRINT*,' reset_flux_scheme_thresholdL1_max (with unit) = ',usrconfig%reset_flux_scheme_thresholdL1_max*usr_physunit%myconfig%luminosity
        !PRINT*,'reset_flux_scheme_thresholdL1_min (without unit) = ',usrconfig%reset_flux_scheme_thresholdL1_min
        !PRINT*,' reset_flux_scheme_thresholdL1_min (with unit) = ',usrconfig%reset_flux_scheme_thresholdL1_min*usr_physunit%myconfig%luminosity
        !PRINT*,'L1<reset_flux_scheme_thresholL1_max, L1>reset_flux_scheme_thresholL1_min : ', &
        !any(w(ixI^S,phys_ind%Lcool1_)<usrconfig%reset_flux_scheme_thresholdL1_max), &
        !any(w(ixI^S,phys_ind%Lcool1_)>usrconfig%reset_flux_scheme_thresholdL1_min)
        !PRINT*,' min, max : ', &
        !MINVAL(w(ixI^S,phys_ind%Lcool1_))*usr_physunit%myconfig%luminosity, &
        !maxL1*usr_physunit%myconfig%luminosity
        !PRINT*,w(ixI^S,phys_ind%Lcool1_)
      !end if
    !end if
    if(usrconfig%reset_flux_scheme_on.or.usrconfig%reset_limiter_on)then
      if(usrconfig%reset_flux_scheme_on)then
        select case(idims)
        case(1)
          !bigtest = any(w(ixI^S,phys_ind%Lcool1_)>usrconfig%reset_flux_scheme_thresholdL1_min)
          !midtest=any((w(ixI^S,phys_ind%Lcool1_)>0.0_dp).and.(w(ixI^S,phys_ind%Lcool1_)<usrconfig%reset_flux_scheme_thresholdL1_max))
          !bigtest=bigtest.and.midtest
          if(phys_config%cool_saveL)then
            if(old_method/=usr_method)then
              usrconfig%reset_flux_scheme_old_method = trim(old_method)
            end if
            !if(bigtest)then
              usr_method = trim(usrconfig%reset_flux_scheme_diffuse)
              usr_limiter= limiter_minmod
            !else
            !  usr_method = trim(usrconfig%reset_flux_scheme_old_method)
            !end if
          else
            ! dummy
          end if
        case(2)
          !bigtest = any(w(ixI^S,phys_ind%Lcool1_)>usrconfig%reset_flux_scheme_thresholdL1_min)
          !midtest=any((w(ixI^S,phys_ind%Lcool1_)>0.0_dp).and.(w(ixI^S,phys_ind%Lcool1_)<usrconfig%reset_flux_scheme_thresholdL1_max))
          !bigtest=bigtest.and.midtest
          if(phys_config%cool_saveL)then
            if(old_method/=usr_method)then
              usrconfig%reset_flux_scheme_old_method = trim(old_method)
            end if
            !if(bigtest)then
              usr_method = trim(usrconfig%reset_flux_scheme_diffuse)
              usr_limiter= limiter_minmod
            !else
            !  usr_method = trim(usrconfig%reset_flux_scheme_old_method)
            !end if
          else
            ! dummy
          end if
        case(3)
          !dummy
        case default
        end select
      end if
    end if
  end if
  return
    theta_min=5.0_dp*dpi/180.0_dp
    usr_method  = old_method
    usr_limiter = old_limiter
    cond_pw : if(any(w(ixO^S,lfac_)>2.0_dp))then
     usr_method  = 'tvdlf'
     usr_limiter = limiter_minmod
    end if cond_pw
   {^NOONED
    where(x(ixO^S,zjet_)>smalldouble)
      theta(ixO^S) = atan(x(ixO^S,r_)/x(ixO^S,zjet_))
    elsewhere
      theta(ixO^S) = 0.0_dp
    end where
    cond_axis : if(any(theta(ixO^S)<20))then
      usr_method  = 'tvdlf'
      usr_limiter = old_limiter
    end if cond_axis
   }
  end subroutine special_reset_solver
  !---------------------------------------------------------------------
  !> subroutine to write simulation configuration
  subroutine usr_write_setting
    implicit none
    integer,parameter   :: unit_config =12
    character(len=75)   :: filename_config
    integer             :: i_cloud,i_ism,i_jet_yso
    real(kind=dp)       :: restst
    !-------------------------------------
    filename_config=trim(base_filename)//'.config'

    open(unit_config,file=trim(filename_config), status='replace')
    write(unit_config,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    write(unit_config,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    write(unit_config,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    write(unit_config,*)'%%%%%%%%%%% Simulation configuration %%%%%%%%%%%%'
    write(unit_config,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    write(unit_config,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    write(unit_config,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    if(usrconfig%physunit_on)call usr_physunit%write_setting(unit_config)
    write(unit_config,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    write(unit_config,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    write(unit_config,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    if(usrconfig%ism_on)then
      Loop_isms : do i_ism=0,usrconfig%ism_number-1
       call ism_surround(i_ism)%write_setting(unit_config)
      end do Loop_isms
    end if
    write(unit_config,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    write(unit_config,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    write(unit_config,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    if(usrconfig%cloud_on)then
      Loop_clouds : do i_cloud=0,usrconfig%cloud_number-1
       call cloud_medium(i_cloud)%write_setting(unit_config)
      end do Loop_clouds
    end if
    write(unit_config,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    write(unit_config,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    write(unit_config,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'


    write(unit_config,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    write(unit_config,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    write(unit_config,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    if(usrconfig%jet_yso_on)then
      Loop_jet_yso : do i_jet_yso=0,usrconfig%jet_yso_number-1
       call jet_yso(i_jet_yso)%write_setting(unit_config)
      end do Loop_jet_yso
    end if
    write(unit_config,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    write(unit_config,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    write(unit_config,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'

    call grackle_object%write_setting(unit_config)


    write(unit_config,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    write(unit_config,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    write(unit_config,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    restst = 7. / 6. / 5.
    #write(unit_config,*)'%%%%%7/6/5 =',restst
    write(unit_config,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    close(unit_config)

  end subroutine usr_write_setting

!----------------------------------------------------------------------
!> compute the total mass and volume in the cloud
subroutine usr_global_var
  use mod_global_parameters

end subroutine usr_global_var
!-------------------------------------------------------------------------------------------
subroutine usr_set_profile_cloud(ixI^L,ixO^L,x,cloud_isuse,cloud_profile)

  use mod_global_parameters
  implicit none
  integer, intent(in)          :: ixI^L,ixO^L
  real(kind=dp), intent(in)    :: x(ixI^S,1:ndim)
  type(cloud), intent(in)      :: cloud_isuse
  real(kind=dp), intent(inout) :: cloud_profile(ixI^S,1:nw)
  ! .. local ..
  real(kind=dp)                :: fprofile(ixI^S),distance(ixI^S)
  real(kind=dp)                :: standart_deviation,angle_theta
  character(len=20)            :: profile
  integer                      :: idir
  !------------------------------------------------------------
  cloud_profile=1.0_dp

  angle_theta = cloud_isuse%myconfig%extend(3)*&
                          usr_physunit%myconfig%length*dpi/180.0_dp

  distance(ixO^S) =   (x(ixO^S,zjet_)-&
                      ( dtan(angle_theta)*(x(ixO^S,r_)-xprobmin1))+&
                        cloud_isuse%myconfig%extend(2))*dsin(angle_theta)

  select case(usrconfig%cloud_structure)
  case(1)
   profile='tanh'
   standart_deviation = cloud_isuse%myconfig%extend(2)
  case(2)
   profile='linear'
   standart_deviation = 2.0_dp*cloud_isuse%myconfig%extend(2)
  case default
   profile='none'
   standart_deviation = cloud_isuse%myconfig%extend(2)
  end select


  call usr_mat_profile_dist(ixI^L,ixO^L,profile, distance,&
                            standart_deviation,fprofile)

  select case(usrconfig%cloud_structure)
  case(1)
    where(distance(ixO^S)>cloud_isuse%myconfig%extend(2))
      fprofile(ixO^S) = 10.0_dp * fprofile(ixO^S)
    end where
  case(2)
    where(distance(ixO^S)>cloud_isuse%myconfig%extend(2))
      fprofile(ixO^S) = 10.0_dp !* fprofile(ixO^S)
    end where

  case default
    fprofile = 1.0_dp
  end select

  if(usrconfig%cloud_profile_density_on)then
    where(cloud_isuse%patch(ixO^S))
     cloud_profile(ixO^S,phys_ind%rho_) = fprofile(ixO^S)
    end where
  else
    where(cloud_isuse%patch(ixO^S))
      cloud_profile(ixO^S,phys_ind%rho_) = 1.0_dp
    end where
  end if

  if(usrconfig%cloud_profile_pressure_on)then
    where(cloud_isuse%patch(ixO^S))
     cloud_profile(ixO^S,phys_ind%pressure_) = fprofile(ixO^S)
    end where
  else
    where(cloud_isuse%patch(ixO^S))
     cloud_profile(ixO^S,phys_ind%pressure_) = 1.0_dp
    end where
  end if
  if(usrconfig%cloud_profile_velocity_on)then
    Loop_idir_0 : do idir =1,ndir
     where(cloud_isuse%patch(ixO^S))
      cloud_profile(ixO^S,phys_ind%mom(idir)) = (1.0_dp-fprofile(ixO^S))/2.0_dp
     end where
   end do  Loop_idir_0
  else
    Loop_idir_1 : do idir =1,ndir
    where(cloud_isuse%patch(ixO^S))
      cloud_profile(ixO^S,phys_ind%mom(idir)) = 1.0_dp
    end where
  end do  Loop_idir_1
  end if
end subroutine usr_set_profile_cloud
end module mod_usr
