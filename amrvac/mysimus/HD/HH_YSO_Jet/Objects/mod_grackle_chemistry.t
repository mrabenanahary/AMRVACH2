module mod_grackle_chemistry
  use mod_constants
  use mod_global_parameters
  use mod_physics
  use mod_obj_dust, only : dust
  use mod_obj_global_parameters
  use mod_obj_usr_unit
  use mod_obj_mat
  use grackle_header
  use mod_grackle_parameters


  type gr_config
    character(len=78)    :: obj_name(max_num_parameters)       !> Obj name that call it
    character(len=20)    :: unit(max_num_parameters)           !> physical unit at parameter file
    integer              :: myindice(max_num_parameters)       !> ism indices associated with ism in use

    ! Parameters as taken from RAMSES implimentation with Grackle
    ! Chemistry with Grackle : on=1 ; off=0
    INTEGER :: use_grackle(max_num_parameters)

    ! Cooling with Grackle : on=1 ; off=0
    INTEGER :: gr_with_radiative_cooling(max_num_parameters)

    ! Molecular/atomic network solved by Grackle :
    ! 0: no chemistry network. Radiative cooling for primordial species is solved by interpolating from lookup tables calculated with Cloudy.
    ! 1: 6-species atomic H and He.
    ! Active species: H, H+, He, He+, He++, e-.
    ! 2: 9-species network including atomic species above and species
    ! for molecular hydrogen formation.
    ! This network includes formation from the H- and H2+ channels,
    ! three-body formation (H+H+H and H+H+H2),
    ! H2 rotational transitions, chemical heating,
    ! and collision-induced emission (optional).
    ! Active species: above + H-, H2, H2+.
    ! 3: 12-species network include all above plus HD rotation cooling.
    ! Active species: above + D, D+, HD.
    ! Reactions listed in Table 3 of Smith et al. 2017
    INTEGER :: gr_primordial_chemistry(max_num_parameters)

    ! Include metal cooling (Smith et al. 2017) :
    ! on = 1 ; off = 0
    ! If enabled, the cooling table to be
    ! used must be specified with the grackle_data_file
    INTEGER :: gr_metal_cooling(max_num_parameters)

    ! Include UV background :
    ! on = 1 ; off = 0
    ! If enabled, the cooling table
    ! to be used must be specified with the grackle_data_file parameter
    INTEGER :: gr_UVbackground(max_num_parameters)

    ! Flag to enable an effective CMB temperature floor.
    ! on = 1 ; off = 0
    ! This is implemented by subtracting the value
    ! of the cooling rate at TCMB from the total cooling rate.
    ! If enabled, the cooling table
    ! to be used must be specified with the grackle_data_file parameter
    INTEGER :: gr_cmb_temperature_floor(max_num_parameters)

    ! Flag to enable H2 formation on dust grains,
    ! dust cooling, and dust-gas heat transfer
    ! follow Omukai (2000). This assumes that the dust
    ! to gas ratio scales with the metallicity.
    INTEGER :: gr_h2_on_dust(max_num_parameters)


    ! Flag to control additional dust cooling
    ! and chemistry processes.
    ! 0: no dust-related processes included.
    ! 1: adds the following processes:
    !   1. photo-electric heating (sets photoelectric_heating to 2).
    !   2. cooling from electron recombination onto dust (equation 9
    !      from Wolfire et al. 1995). Both the photo-electric heating
    !      and recombination cooling are scaled by the value of the
    !      interstellar_radiation_field.
    !   3. H2 formation on dust (sets h2_on_dust to 1 if primordial_chemistry > 1).
    ! Setting dust_chemistry greater than 0 requires metal_cooling to be enabled.
    INTEGER :: gr_dust_chemistry(max_num_parameters)

    ! Flag to provide the dust density as a field using the dust_density
    ! pointer in the grackle_field_data struct. If set to 0,
    ! the dust density takes the value of local_dust_to_gas_ratio
    ! multiplied by the metallicity
    INTEGER :: gr_use_dust_density_field(max_num_parameters)

    ! Flag to enable photo-electric heating from irradiated dust grains. Default: 0.
    ! 0: no photo-electric heating.
    ! 1: a spatially uniform heating term
    ! from Tasker & Bryan (2008). The exact
    ! heating rate used must be specified
    ! with the photoelectric_heating_rate
    ! parameter. For temperatures above
    ! 20,000 K, the photo-electric heating
    ! rate is set to 0.
    ! 2: similar to option 1, except the
    ! heating rate is calculated using
    ! equation 1 of Wolfire et al. (1995)
    ! and the user must supply the intensity
    ! of the interstellar radiation field with the
    ! interstellar_radiation_field parameter.
    ! The value of epsilon is taken as a constant equal to 0.05
    ! for gas below 20,000 K and 0 otherwise.
    ! 3: similar to option 1, except the value of
    ! epsilon is calculated directly
    ! from equation 2 of Wolfire et al. (1995).
    INTEGER :: gr_photoelectric_heating(max_num_parameters)

    ! Flag to signal that an array of volumetric
    ! heating rates is being provided in the
    ! volumetric_heating_rate field of the
    ! grackle_field_data struct.
    ! on = 1 ; off = 0
    INTEGER :: gr_use_volumetric_heating_rate(max_num_parameters)

    ! Flag to signal that an array of
    ! specific heating rates is being
    ! provided in the specific_heating_rate
    ! field of the grackle_field_data struct.
    ! on = 1 ; off = 0
    INTEGER :: gr_use_specific_heating_rate(max_num_parameters)

    ! Flag to control which three-body H2 formation rate is used.
    !0: Abel, Bryan & Norman (2002)
    !1: Palla, Salpeter & Stahler (1983)
    !2: Cohen & Westberg (1983)
    !3: Flower & Harris (2007)
    !4: Glover (2008)
    !5: Forrey (2013).
    !The first five options are discussed in Turk et. al. (2011).
    INTEGER :: gr_three_body_rate(max_num_parameters)

    ! Flag to enable H2 collision-induced
    ! emission cooling from Ripamonti & Abel (2004).
    ! on = 1 ; off = 0
    INTEGER :: gr_cie_cooling(max_num_parameters)

    ! Flag to enable H2 cooling attenuation
    ! from Ripamonti & Abel (2004).
    ! on = 1 ; off = 0
    INTEGER :: gr_h2_optical_depth_approximation(max_num_parameters)

    INTEGER :: gr_ih2co(max_num_parameters)
    INTEGER :: gr_ipiht(max_num_parameters)
    INTEGER :: gr_NumberOfTemperatureBins(max_num_parameters)

    ! 0 : The recombination of H + , He + and He ++
    ! is modelled using the case A recombination rate coefficients
    ! (the optically-thin approximation in which recombination
    ! photons above 1 Ryd escape).
    ! 1 : case B rate coefficients (in which recombination photons above 1 Ryd
    ! are locally re-absorbed, Osterbrock 1989) can instead be se-
    ! lected by setting CaseBRecombination = 1.
    INTEGER :: gr_CaseBRecombination(max_num_parameters)

    ! Flag to enable Compton heating
    ! from an X-ray background following
    ! Madau & Efstathiou (1999).
    ! on = 1 ; off = 0
    INTEGER :: gr_Compton_xray_heating(max_num_parameters)

    ! Flag to enable suppression of Lyman-Werner flux due to Lyman-series
    ! absorption (giving a sawtooth pattern), taken from Haiman & Abel, & Rees (2000)
    ! on = 1 ; off = 0
    INTEGER :: gr_LWbackground_sawtooth_suppression(max_num_parameters)

    INTEGER :: gr_NumberOfDustTemperatureBins(max_num_parameters)


    ! Flag to signal that arrays of ionization and heating rates from
    ! radiative transfer solutions are being provided.
    ! Only available if primordial_chemistry is greater than 0.
    ! HI, HeI, and HeII ionization arrays are provided in
    ! RT_HI_ionization_rate, RT_HeI_ionization_rate, and RT_HeII_ionization_rate fields, respectively,
    ! of the grackle_field_data struct. Associated heating rate is
    ! provided in the RT_heating_rate field, and H2photodissociation
    ! rate can also be provided in the RT_H2_dissociation_rate field
    ! when primordial_chemistry is set to either 2 or 3.
    ! on = 1 ; off = 0
    INTEGER :: gr_use_radiative_transfer(max_num_parameters)


    ! When used with use_radiative_transfer set to 1,
    ! this flag makes it possible to solve the chemistry and cooling
    ! of the computational elements for which the radiation field is non-zero
    ! separately from those with no incident radiation. This allows radiation transfer calculations to
    ! be performed on a smaller timestep than the global timestep. The parameter, radiative_transfer_intermediate_step,
    ! is then used to toggle between updating the cells/particles receiving radiative input and those that
    ! are not.
    ! on = 1 ; off = 0
    INTEGER :: gr_radiative_transfer_coupled_rate_solver(max_num_parameters)

    ! Used in conjunction with radiative_transfer_coupled_rate_solver set to 1, setting this parameter to 1 tells the solver
    ! to only update cells/particles where the radiation field is non-zero. Setting this to 0 updates only those elements with
    ! no incident radiation. When radiative_transfer_coupled_rate_solver is set to 0, changing this parameter
    ! will have no effect.
    INTEGER :: gr_radiative_transfer_intermediate_step(max_num_parameters)

    ! Flag to only use hydrogen ionization and heating rates from the radiative transfer solutions.
    ! on = 1 ; off = 0
    INTEGER :: gr_radiative_transfer_hydrogen_only(max_num_parameters)

    ! Switch to enable approximate self-shielding from the UV background. All three of the below methods incorporate Eq. 13 and
    ! 14 from Rahmati et. al. 2013. These equations involve using the spectrum averaged photoabsorption cross for the given species (HI or HeI).
    ! These redshift dependent values are pre-computed for the HM2012 and FG2011 UV backgrounds and included in their respective cooling data tables.
    ! Care is advised in using any of these methods. The default behavior is to apply no self-shielding, but this is not necessarily the proper assumption,
    ! depending on the use case. If the user desires to turn on self-shielding, we strongly advise using option 3. All options include HI self-shielding,
    ! and vary only in treatment of HeI and HeII. In options 2 and 3, we approximately account for HeI self-shielding by applying the Rahmati et. al. 2013 relations,
    ! which are only strictly valid for HI, !to HeI under the assumption that it behaves similarly to HI. None of these options are completely correct in practice,
    ! but option 3 has produced the most reasonable results in test simulations. Repeating the analysis of Rahmati et. al. 2013 to directly parameterize HeI and HeII
    ! self-shielding behavior would be a valuable avenue of future research in developing a more complete self-shielding model. Each self-shielding option is described below.
    ! 0: No self shielding. Elements are optically thin to the UV background.
    ! 1: Not Recommended. Approximate self-shielding in HI only.
    ! HeI and HeII are left as optically thin.
    ! 2: Approximate self-shielding in both HI and HeI. HeII remains
    ! optically thin.
    ! 3: Approximate self-shielding in both HI and HeI, but ignoring
    ! HeII ionization and heating from the UV background entirely (HeII ionization and heating rates are set to zero).
    ! These methods only work in conjunction with using updated Cloudy cooling tables, denoted with “_shielding”. These tables properly account for the decrease
    ! in metal line cooling rates in self-shielded regions, which can be significant.
    !For consistency, when primordial_chemistry > 2, the self-shielding attenutation factors calculated for HI and HeI are applied to the H2ionization (15.4 eV) and H2+
    ! dissociation rates (30 eV) respectively. These reaction rates are distinct from the H2self-shielding computed using the H2_self_shielding flag.
    ! on = 1 ; off = 0
    INTEGER :: gr_self_shielding_method(max_num_parameters)
    INTEGER :: gr_use_isrf_field(max_num_parameters)


    ! The ratio of specific heats for an ideal gas.
    ! A direct calculation for the molecular component
    ! is used if primordial_chemistry > 1.
    REAL(kind=gr_rpknd)    :: gr_Gamma(max_num_parameters)

    ! If photoelectric_heating is enabled, the heating rate in units of (erg cm-3/s)
    ! n-1, where n is the total hydrogen number density. In other words, this is the
    ! volumetric heating rate at a hydrogen number density of n = 1 cm-3.
    REAL(kind=gr_rpknd)    :: gr_photoelectric_heating_rate(max_num_parameters)



    ! Temperature limits
    REAL(kind=gr_rpknd)    :: TemperatureStart(max_num_parameters)
    REAL(kind=gr_rpknd)    :: TemperatureEnd(max_num_parameters)

    ! Dust temperature limits
    REAL(kind=gr_rpknd)    :: DustTemperatureStart(max_num_parameters)
    REAL(kind=gr_rpknd)    :: DustTemperatureEnd(max_num_parameters)

    ! Intensity of a constant Lyman-Werner H2 photo-dissociating
    ! radiation field in units of 10-21 erg /s /cm2 Hz-1 sr-1. Default: 0.
    REAL(kind=gr_rpknd)    :: gr_LWbackground_intensity(max_num_parameters)

    ! Used in combination with UVbackground_redshift_fullon, UVbackground_redshift_drop,
    ! and UVbackground_redshift_off to set an attenuation factor
    ! for the photo-heating and photo-ionization rates of the UV background model.
    ! See the figure below for an illustration its behavior. If not set, this parameter will
    ! be set to the highest redshift of the UV background data being used.
    REAL(kind=gr_rpknd)    :: gr_UVbackground_redshift_on(max_num_parameters)

    ! Used in combination with UVbackground_redshift_on, UVbackground_redshift_fullon, and
    ! UVbackground_redshift_drop to set an attenuation factor for the photo-heating and
    ! photo-ionization rates of the UV background model. See the figure below for an illustration its behavior.
    ! If not set, this parameter will be set to the lowest redshift of the UV background data being used.
    REAL(kind=gr_rpknd)    :: gr_UVbackground_redshift_off(max_num_parameters)

    ! Used in combination with UVbackground_redshift_on, UVbackground_redshift_drop, and UVbackground_redshift_off
    ! to set an attenuation factor for the photo-heating and photo-ionization rates of the UV background model.
    ! See the figure below for an illustration its behavior. If not set, this parameter will be set to the highest
    ! redshift of the UV background data being used.
    REAL(kind=gr_rpknd)    :: gr_UVbackground_redshift_fullon(max_num_parameters)

    ! Used in combination with UVbackground_redshift_on, UVbackground_redshift_fullon, and UVbackground_redshift_off
    ! to set an attenuation factor for the photo-heating and photo-ionization rates of the UV background model.
    ! See the figure below for an illustration its behavior. If not set, this parameter will be set to the lowest
    ! redshift of the UV background data being used.
    REAL(kind=gr_rpknd)    :: gr_UVbackground_redshift_drop(max_num_parameters)




    ! Cloudy 07.02 abundances :
    ! A float value to account for additional electrons contributed by metals. This is only used with Cloudy datasets
    ! with dimension greater than or equal to 4. The value of this factor is calculated as the sum of (Ai * i) over all
    ! elements i heavier than He, where Ai is the solar number abundance relative to H. For the solar abundance pattern
    ! from the latest version of Cloudy, using all metals through Zn, this value is 9.153959e-3. Default: 9.153959e-3.
    REAL(kind=gr_rpknd)    :: cloudy_electron_fraction_factor(max_num_parameters)
    CHARACTER(LEN=128) :: data_filename(max_num_parameters)
    CHARACTER(LEN=128) :: data_dir(max_num_parameters)

    logical :: normalize_done(max_num_parameters)
    INTEGER :: gr_comoving_coordinates(max_num_parameters)
    REAL(kind=gr_rpknd)    :: gr_a_units(max_num_parameters)
    REAL(kind=gr_rpknd)    :: gr_density_units(max_num_parameters)
    REAL(kind=gr_rpknd)    :: gr_length_units(max_num_parameters)
    REAL(kind=gr_rpknd)    :: gr_time_units(max_num_parameters)
    REAL(kind=gr_rpknd)    :: gr_current_redshift(max_num_parameters)
      real*8               :: a_value(max_num_parameters)





    CHARACTER(LEN=257)   :: data_file(max_num_parameters)
    ! Gas to dust ratios
    real(kind=gr_rpknd)  :: chi_dust(max_num_parameters)
    real(kind=gr_rpknd)  :: xi_dust(max_num_parameters)

    ! Fraction by masses X,Y,Z in the baryonic gas (i.e. without electrons and dust)
    ! the values inferred from Asplund, Grevesse & Sauval (2005) by Asplund et al. (2009)
    ! for the protosolar mass fractions: X=0.7166 Y=0.2704 Z=0.0130
    !
    ! The fraction by mass of Hydrogen in the baryonic gas (i.e. without electrons and dust)
    ! This is the famous X parameter
    REAL(kind=gr_rpknd)    :: HydrogenFractionByMass(max_num_parameters)

    ! The fraction by mass of Helium in the baryonic gas (i.e. without electrons and dust)
    ! This is the famous Y parameter
    REAL(kind=gr_rpknd)    :: HeliumFractionByMass(max_num_parameters)

    ! The fraction of total gas mass in metals for a solar composition.
    ! Default: 0.01295 (consistent with the default abundances in the Cloudy code).
    REAL(kind=gr_rpknd)    :: SolarMetalFractionByMass(max_num_parameters) !Z_solar

    ! The fraction by mass of Metals in the baryonic gas (i.e. without electrons and dust)
    ! This is the famous Z parameter
    real(kind=gr_rpknd)  :: MetalFractionByMass(max_num_parameters)

    ! The fraction by mass of Hydrogen in the metal-free portion of the gas (i.e., just the H and He).
    ! In the non-equilibrium solver, this is used to ensure consistency in the densities of the individual species.
    ! In tabulated mode, this is used to calculate the H number density from the total gas density,
    ! which is a parameter of the heating/cooling tables. When using the non-equilibrium solver,
    ! a sensible default is 0.76. However, the tables for tabulated mode were created assuming nHe/nH = 0.1,
    ! which corresponds to an H mass fraction of about 0.716.
    ! When running in tabulated mode, this parameter will automatically be changed to this value.


    ! The ratio by mass of Deuterium to Hydrogen.
    ! Default: 6.8e-5 (the value from Burles & Tytler (1998)
    ! multiplied by 2 for the mass of Deuterium).
    REAL(kind=gr_rpknd)    :: DeuteriumToHydrogenRatio(max_num_parameters) !chi_D


    ! Ionization fraction by number density : x_ion = n(e-)/n_H
    REAL(kind=gr_rpknd)    :: IonizationFraction(max_num_parameters) !x_ion

    ! He to H abundance : x(He) = (n(He)+n(He+)+n(He++))/n_H
    REAL(kind=gr_rpknd)   :: He_abundance(max_num_parameters)  !x(He)

    ! Zeta constant : zeta =(rho-rhodust)/mH
    REAL(kind=gr_rpknd)   :: Zeta_nH_to_rho(max_num_parameters)  !zeta

    ! Zeta prime constant : zeta prime  = Zeta - x(He)w_He/(1-Z)
    REAL(kind=gr_rpknd)   :: ZetaPrime_nH_to_rho(max_num_parameters)  !zeta prime

    ! Varsigma constant : Varsigma =(rho-rhodust)/mH
    REAL(kind=gr_rpknd)   :: Varsigma_nH_to_rho(max_num_parameters)  !Varsigma

    ! Varsigma prime constant : Varsigma prime  = Varsigma - x(He)w_He/(1-Z)
    REAL(kind=gr_rpknd)   :: VarsigmaPrime_nH_to_rho(max_num_parameters)  !Varsigma prime

    logical :: correctdensities(max_num_parameters)
    REAL(kind=gr_rpknd) :: factor_to_density(max_num_parameters)
    REAL(kind=gr_rpknd) :: deviation_to_density(max_num_parameters)
    REAL(kind=gr_rpknd) :: deviation_to_density_limit(max_num_parameters)


    !Species abundances
    integer             :: number_of_solved_species(max_num_parameters)
    real(kind=gr_rpknd) :: x_HI(max_num_parameters)
    real(kind=gr_rpknd) :: x_HII(max_num_parameters)
    real(kind=gr_rpknd) :: x_HM(max_num_parameters)
    real(kind=gr_rpknd) :: x_HeI(max_num_parameters)
    real(kind=gr_rpknd) :: x_HeII(max_num_parameters)
    real(kind=gr_rpknd) :: x_HeIII(max_num_parameters)
    real(kind=gr_rpknd) :: x_H2I(max_num_parameters)
    real(kind=gr_rpknd) :: x_H2II(max_num_parameters)
    real(kind=gr_rpknd) :: x_DI(max_num_parameters)
    real(kind=gr_rpknd) :: x_DII(max_num_parameters)
    real(kind=gr_rpknd) :: x_HDI(max_num_parameters)
    real(kind=gr_rpknd) :: x_e(max_num_parameters)


    ! Traditional MPI-AMRVAC parameters
    real(kind=dp)   :: mean_nall_to_nH(max_num_parameters)
    real(kind=dp)   :: mean_mass(max_num_parameters)
    real(kind=dp)   :: mean_mup(max_num_parameters)
    real(kind=dp)   :: mean_ne_to_nH(max_num_parameters)

    real(dp)        :: density(max_num_parameters)        !rho=rho_total
    real(dp)        :: number_H_density(max_num_parameters) !> n_H
    real(dp)        :: density_dust(max_num_parameters)        !rhodust
    real(dp)        :: density_tot(max_num_parameters)        !rhogas
    real(dp)        :: density_bar(max_num_parameters)        !rhobaryonic
    real(dp)        :: density_X(max_num_parameters)        !rhoX
    real(dp)        :: density_Y(max_num_parameters)        !rhoY
    real(dp)        :: density_Z(max_num_parameters)        !rhoZ
    real(dp)        :: density_deut(max_num_parameters)        !rhodeut
    real(dp)        :: density_not_deut(max_num_parameters)        !rhonotdeut
    real(dp)        :: density_metal_free(max_num_parameters)        !rhometalfree (no deuterium)
    real(dp)        :: densityD(max_num_parameters)
    real(dp)        :: densityHD(max_num_parameters)
    real(dp)        :: densityDplusHD(max_num_parameters)
    real(dp)        :: densityH(max_num_parameters)
    real(dp)        :: densityHtwo(max_num_parameters)
    real(dp)        :: densityHplusH2(max_num_parameters)
    real(dp)        :: densityHe(max_num_parameters)
    real(dp)        :: densityElectrons(max_num_parameters)
    ! to autoscale rhoDI + rhoDII + rhoHDI + rhoHI + rhoHII + rhoHM + rhoH2I + rhoH2II
    real(dp)        :: densityDI(max_num_parameters)
    real(dp)        :: densityDII(max_num_parameters)
    real(dp)        :: densityHDI(max_num_parameters)
    real(dp)        :: densityHI(max_num_parameters)
    real(dp)        :: densityHII(max_num_parameters)
    real(dp)        :: densityHM(max_num_parameters)
    real(dp)        :: densityH2I(max_num_parameters)
    real(dp)        :: densityH2II(max_num_parameters)
    ! to autoscale rhoHeI + rhoHeII + rhoHeIII
    real(dp)        :: densityHeI(max_num_parameters)
    real(dp)        :: densityHeII(max_num_parameters)
    real(dp)        :: densityHeIII(max_num_parameters)

    logical :: boundary_on
    character(len=30)    :: boundary_cond(3,2)

  end type gr_config


  type gr_solver
    integer         :: number_of_objects
    type(gr_config) :: myconfig
    type(usrboundary_type)          :: myboundaries
    type(usrphysical_unit), pointer :: myphysunit
    contains

    !=BEGIN PREINITIALIZATION PROCESS SUBROUTINES=====================================================================================

    !=END PREINITIALIZATION PROCESS SUBROUTINES=====================================================================================
    PROCEDURE, PASS(self) :: set_default_config => grackle_set_default
    PROCEDURE, PASS(self) :: read_parameters => grackle_config_read
    PROCEDURE, PASS(self) :: write_setting        => grackle_write_setting
    PROCEDURE, PASS(self) :: set_complet        => grackle_set_complet
    PROCEDURE, PASS(self) :: normalize       => grackle_normalize
    PROCEDURE, PASS(self) :: denormalize       => grackle_denormalize

    !PROCEDURE, PASS(self) :: check_nH_abundances   => check_nH_abundances_consistency
    !PROCEDURE, PASS(self) :: check_xHe_abundances   => check_xhe_abundances_consistency
    !PROCEDURE, PASS(self) :: check_Xfraction_abundances   => check_Xfraction_abundances_consistency
    !PROCEDURE, PASS(self) :: check_densities   => check_densities_consistency
    !PROCEDURE, PASS(self) :: check_deviation   => check_deviation_by_density
    !PROCEDURE, PASS(self) :: set_w_init   => usr_gr_set_w_init
    !PROCEDURE, PASS(self)        :: clean_grid               => grackle_deallocate_wgr
    !PROCEDURE, PASS(self)        :: set_grid                 => grackle_allocate_wgr
    PROCEDURE, PASS(self)        :: link_par_to_gr           => grackle_solver_associate
    PROCEDURE, PASS(self)        :: set_parameters           => grackle_chemistry_set_parameters
    !PROCEDURE, PASS(self) :: associate_units          => grackle_solver_associate_units
    PROCEDURE, PASS(self)        :: grackle_source           => grackle_chemistry_static_source

  end type gr_solver



  character(len=257), TARGET :: general_grackle_filename





CONTAINS

!-------------------------------------------------------------------------
!> subroutine default setting for Grackle
subroutine grackle_set_default(self)
  implicit none
  class(gr_solver)            :: self
  integer :: istr
  !----------------------------------
  !TODO: check consistency with set_default_chemistry_parameters.c in grackle /src/clib
  self%myconfig%obj_name(1:max_num_parameters)              = 'grackle_chemistry_config'
  self%myconfig%unit(1:max_num_parameters)                  = 'cgs'
  self%myconfig%myindice(1:max_num_parameters)              = 0

  self%myconfig%use_grackle(1:max_num_parameters) = 1 !don t forget to change to 0 after coding
  self%myconfig%gr_with_radiative_cooling(1:max_num_parameters) = 1
  self%myconfig%gr_primordial_chemistry(1:max_num_parameters) = 0
  self%myconfig%gr_metal_cooling(1:max_num_parameters) = 1 !don t forget to change to 0 after coding
  self%myconfig%gr_UVbackground(1:max_num_parameters) = 0
  self%myconfig%gr_cmb_temperature_floor(1:max_num_parameters) = 1
  self%myconfig%gr_h2_on_dust(1:max_num_parameters) = 0
  self%myconfig%gr_dust_chemistry(1:max_num_parameters) = 0
  self%myconfig%gr_use_dust_density_field(1:max_num_parameters) = 0
  self%myconfig%gr_photoelectric_heating(1:max_num_parameters) = -1
  self%myconfig%gr_photoelectric_heating_rate(1:max_num_parameters) = 8.5D-26
  self%myconfig%gr_use_volumetric_heating_rate(1:max_num_parameters) = 0
  self%myconfig%gr_use_specific_heating_rate(1:max_num_parameters) = 0
  self%myconfig%gr_three_body_rate(1:max_num_parameters) = 0
  self%myconfig%gr_cie_cooling(1:max_num_parameters) = 0
  self%myconfig%gr_h2_optical_depth_approximation(1:max_num_parameters) = 0
  self%myconfig%gr_ih2co(1:max_num_parameters) = 1
  self%myconfig%gr_ipiht(1:max_num_parameters) = 1
  self%myconfig%gr_NumberOfTemperatureBins(1:max_num_parameters) = 600
  self%myconfig%gr_CaseBRecombination(1:max_num_parameters) = 0
  self%myconfig%gr_Compton_xray_heating (1:max_num_parameters)= 0
  self%myconfig%gr_LWbackground_sawtooth_suppression(1:max_num_parameters) = 0
  self%myconfig%gr_NumberOfDustTemperatureBins(1:max_num_parameters) = 250
  self%myconfig%gr_use_radiative_transfer(1:max_num_parameters) = 0
  self%myconfig%gr_radiative_transfer_coupled_rate_solver(1:max_num_parameters) = 0
  self%myconfig%gr_radiative_transfer_intermediate_step(1:max_num_parameters) = 0
  self%myconfig%gr_radiative_transfer_hydrogen_only(1:max_num_parameters) = 0
  self%myconfig%gr_self_shielding_method(1:max_num_parameters) = 0
  self%myconfig%gr_use_isrf_field(1:max_num_parameters)=0
  self%myconfig%gr_Gamma(1:max_num_parameters) = 5.d0/3.d0

  self%myconfig%TemperatureStart(1:max_num_parameters) = 1.0d0
  self%myconfig%TemperatureEnd(1:max_num_parameters) = 1.0D9
  self%myconfig%DustTemperatureStart(1:max_num_parameters) = 1.0d0
  self%myconfig%DustTemperatureEnd(1:max_num_parameters) = 1500.0d0
  self%myconfig%gr_LWbackground_intensity(1:max_num_parameters) = 0.0d0
  self%myconfig%gr_UVbackground_redshift_on(1:max_num_parameters) = -99999.0 !=FLOAT_UNDEFINED
  self%myconfig%gr_UVbackground_redshift_off(1:max_num_parameters) = -99999.0
  self%myconfig%gr_UVbackground_redshift_fullon(1:max_num_parameters) = -99999.0
  self%myconfig%gr_UVbackground_redshift_drop(1:max_num_parameters) = -99999.0
  self%myconfig%cloudy_electron_fraction_factor(1:max_num_parameters) = 9.153959D-3
  self%myconfig%data_dir(1:max_num_parameters) = "../../../src/grackle/input/"
  self%myconfig%data_filename(1:max_num_parameters) = "CloudyData_UVB=HM2012_shielded.h5"

  self%myconfig%data_file(1:max_num_parameters) = "../../../src/grackle/input/"//&
  "CloudyData_UVB=HM2012_shielded.h5"//C_NULL_CHAR

  self%myconfig%normalize_done(1:max_num_parameters) = .false.
  self%myconfig%gr_comoving_coordinates(1:max_num_parameters) = 0
  self%myconfig%gr_a_units(1:max_num_parameters) = 1.0d0
  self%myconfig%gr_density_units(1:max_num_parameters) = 1.0d0
  self%myconfig%gr_length_units(1:max_num_parameters) = 1.0d0
  self%myconfig%gr_time_units(1:max_num_parameters) = 1.0d0
  self%myconfig%gr_current_redshift(1:max_num_parameters) = 0.
  self%myconfig%a_value(1:max_num_parameters) = 1.0d0

  ! Parameters default values
  ! Gas to dust ratios
  self%myconfig%myindice(1:max_num_parameters) = 0

  ! Z = Z _solar :
  self%myconfig%HydrogenFractionByMass(1:max_num_parameters) = 0.76d0
  self%myconfig%HeliumFractionByMass(1:max_num_parameters) = 1.0d0-0.76d0 ! 0.22705d0
  self%myconfig%SolarMetalFractionByMass(1:max_num_parameters) = 0.01295d0 !Solar metallicity
  self%myconfig%MetalFractionByMass(1:max_num_parameters) = self%myconfig%SolarMetalFractionByMass(1:max_num_parameters)

  self%myconfig%chi_dust(1:max_num_parameters) = 0.009387d0
  self%myconfig%xi_dust(1:max_num_parameters) = 0.009387d0*&
  (self%myconfig%MetalFractionByMass(1:max_num_parameters)/&
  self%myconfig%SolarMetalFractionByMass(1:max_num_parameters))

  self%myconfig%DeuteriumToHydrogenRatio(1:max_num_parameters) = 2.0d0*3.4d-5
  self%myconfig%IonizationFraction(1:max_num_parameters) = 0.0d0
  self%myconfig%He_abundance(1:max_num_parameters) = 0.0d0
  self%myconfig%Zeta_nH_to_rho(1:max_num_parameters) = 0.0d0
  self%myconfig%ZetaPrime_nH_to_rho(1:max_num_parameters) = 0.0d0
  self%myconfig%Varsigma_nH_to_rho(1:max_num_parameters) = 0.0d0
  self%myconfig%VarsigmaPrime_nH_to_rho(1:max_num_parameters) = 0.0d0

  self%myconfig%correctdensities(1:max_num_parameters) = .false.
  self%myconfig%factor_to_density(1:max_num_parameters) = 0.0d0
  self%myconfig%deviation_to_density(1:max_num_parameters) = 0.0d0
  self%myconfig%deviation_to_density_limit(1:max_num_parameters) = 1.0d0



  !Species abundances
  self%myconfig%number_of_solved_species(1:max_num_parameters) = 0
  self%myconfig%x_HI(1:max_num_parameters) = 1.0d0
  self%myconfig%x_HII(1:max_num_parameters) = 0.0d0
  self%myconfig%x_HM(1:max_num_parameters) = 0.0d0
  self%myconfig%x_H2I(1:max_num_parameters) = 0.0d0
  self%myconfig%x_H2II(1:max_num_parameters) = 0.0d0
  self%myconfig%x_HeI(1:max_num_parameters) = 1.0d0
  self%myconfig%x_HeII(1:max_num_parameters) = 0.0d0
  self%myconfig%x_HeIII(1:max_num_parameters) = 0.0d0
  self%myconfig%x_DI(1:max_num_parameters) = 1.0d0
  self%myconfig%x_DII(1:max_num_parameters) = 0.0d0
  self%myconfig%x_HDI(1:max_num_parameters) = 0.0d0
  self%myconfig%x_e(1:max_num_parameters) = 1.0d0
  self%myconfig%densityDI(1:max_num_parameters) = 0.0d0
  self%myconfig%densityDII(1:max_num_parameters) = 0.0d0
  self%myconfig%densityHDI(1:max_num_parameters) = 0.0d0
  self%myconfig%densityHI(1:max_num_parameters) = 0.0d0
  self%myconfig%densityHII(1:max_num_parameters) = 0.0d0
  self%myconfig%densityHM(1:max_num_parameters) = 0.0d0
  self%myconfig%densityH2I(1:max_num_parameters) = 0.0d0
  self%myconfig%densityH2II(1:max_num_parameters) = 0.0d0
  self%myconfig%densityHeI(1:max_num_parameters) = 0.0d0
  self%myconfig%densityHeII(1:max_num_parameters) = 0.0d0
  self%myconfig%densityHeIII(1:max_num_parameters) = 0.0d0

  self%myconfig%density(1:max_num_parameters) = 0.0d0
  self%myconfig%number_H_density(1:max_num_parameters) = 0.0d0
  self%myconfig%density_dust(1:max_num_parameters) = 0.0d0
  self%myconfig%density_tot(1:max_num_parameters) = 0.0d0
  self%myconfig%density_bar(1:max_num_parameters) = 0.0d0
  self%myconfig%density_X(1:max_num_parameters) = 0.0d0
  self%myconfig%density_Y(1:max_num_parameters) = 0.0d0
  self%myconfig%density_Z(1:max_num_parameters) = 0.0d0
  self%myconfig%density_deut(1:max_num_parameters) = 0.0d0
  self%myconfig%density_not_deut(1:max_num_parameters) = 0.0d0
  self%myconfig%density_metal_free(1:max_num_parameters) = 0.0d0
  self%myconfig%densityD(1:max_num_parameters) = 0.0d0
  self%myconfig%densityHD(1:max_num_parameters) = 0.0d0
  self%myconfig%densityH(1:max_num_parameters) = 0.0d0
  self%myconfig%densityHtwo(1:max_num_parameters) = 0.0d0
  self%myconfig%densityHe(1:max_num_parameters) = 0.0d0
  self%myconfig%densityElectrons(1:max_num_parameters) = 0.0d0

  self%myconfig%boundary_on           = .false.
  self%myconfig%boundary_cond         = 'fix'

  call self%myboundaries%set_default

  !write(*,*) 'Grackle configuration defaulting successfully done !'
end subroutine grackle_set_default


subroutine grackle_normalize(self,physunit_inuse)
  use mod_obj_usr_unit
  implicit none
  class(gr_solver)                      :: self
  type(usrphysical_unit), target, intent(in),optional     :: physunit_inuse
  integer :: iobject
  !----------------------------------


  do iobject = 1, self%number_of_objects
    if(present(physunit_inuse))then
      self%myphysunit =>physunit_inuse
      if(trim(self%myconfig%unit(iobject))=='code'.or.self%myconfig%normalize_done(iobject))then
         if(self%myconfig%normalize_done(iobject))then
          write(*,*) 'WARNING: Second call for Grackle normalisation', &
                       'no new normalisation will be done'
          cycle
         end if
      end if
    else
      if(self%myconfig%normalize_done(iobject))then
       write(*,*) 'WARNING: Second call for Grackle normalisation of object #', iobject, &
                    ' no new normalisation will be done'
       cycle
      end if
    end if

    !write(*,*) ' density =', self%myconfig%density(iobject)
    self%myconfig%density(iobject)=self%myconfig%density(iobject)/self%myphysunit%myconfig%density
    self%myconfig%densityDI(iobject)=self%myconfig%densityDI(iobject)/self%myphysunit%myconfig%density
    self%myconfig%densityDII(iobject)=self%myconfig%densityDII(iobject)/self%myphysunit%myconfig%density
    self%myconfig%densityHDI(iobject)=self%myconfig%densityHDI(iobject)/self%myphysunit%myconfig%density
    self%myconfig%densityHI(iobject)=self%myconfig%densityHI(iobject)/self%myphysunit%myconfig%density
    self%myconfig%densityHII(iobject)=self%myconfig%densityHII(iobject)/self%myphysunit%myconfig%density
    self%myconfig%densityHM(iobject)=self%myconfig%densityHM(iobject)/self%myphysunit%myconfig%density
    self%myconfig%densityH2I(iobject)=self%myconfig%densityH2I(iobject)/self%myphysunit%myconfig%density
    self%myconfig%densityH2II(iobject)=self%myconfig%densityH2II(iobject)/self%myphysunit%myconfig%density
    self%myconfig%densityHeI(iobject)=self%myconfig%densityHeI(iobject)/self%myphysunit%myconfig%density
    self%myconfig%densityHeII(iobject)=self%myconfig%densityHeII(iobject)/self%myphysunit%myconfig%density
    self%myconfig%densityHeIII(iobject)=self%myconfig%densityHeIII(iobject)/self%myphysunit%myconfig%density
    self%myconfig%densityElectrons(iobject)=self%myconfig%densityElectrons(iobject)/self%myphysunit%myconfig%density
    self%myconfig%density_Z(iobject)=self%myconfig%density_Z(iobject)/self%myphysunit%myconfig%density
    self%myconfig%density_dust(iobject)=self%myconfig%density_dust(iobject)/self%myphysunit%myconfig%density

    self%myconfig%normalize_done(iobject)=.true.

  end do
end subroutine grackle_normalize

subroutine grackle_denormalize(self,physunit_inuse)
  use mod_obj_usr_unit
  implicit none
  class(gr_solver)                      :: self
  type(usrphysical_unit), target, intent(in),optional     :: physunit_inuse
  integer :: iobject
  !----------------------------------


  do iobject = 1, self%number_of_objects
    if(present(physunit_inuse))then
      self%myphysunit =>physunit_inuse
      if(trim(self%myconfig%unit(iobject))=='code'.or.(.not.self%myconfig%normalize_done(iobject)))then
         if(.not.self%myconfig%normalize_done(iobject))then
          write(*,*) 'WARNING: Second call for Grackle denormalisation', &
                       'no new denormalisation will be done'
         end if
         cycle
      end if
    else
      if(.not.self%myconfig%normalize_done(iobject))then
       write(*,*) 'WARNING: Second call for Grackle denormalisation of object #', iobject, &
                    ' no new denormalisation will be done'
      cycle
      end if

    end if
    self%myconfig%density(iobject)=self%myconfig%density(iobject)*self%myphysunit%myconfig%density
    self%myconfig%densityDI(iobject)=self%myconfig%densityDI(iobject)*self%myphysunit%myconfig%density
    self%myconfig%densityDII(iobject)=self%myconfig%densityDII(iobject)*self%myphysunit%myconfig%density
    self%myconfig%densityHDI(iobject)=self%myconfig%densityHDI(iobject)*self%myphysunit%myconfig%density
    self%myconfig%densityHI(iobject)=self%myconfig%densityHI(iobject)*self%myphysunit%myconfig%density
    self%myconfig%densityHII(iobject)=self%myconfig%densityHII(iobject)*self%myphysunit%myconfig%density
    self%myconfig%densityHM(iobject)=self%myconfig%densityHM(iobject)*self%myphysunit%myconfig%density
    self%myconfig%densityH2I(iobject)=self%myconfig%densityH2I(iobject)*self%myphysunit%myconfig%density
    self%myconfig%densityH2II(iobject)=self%myconfig%densityH2II(iobject)*self%myphysunit%myconfig%density
    self%myconfig%densityHeI(iobject)=self%myconfig%densityHeI(iobject)*self%myphysunit%myconfig%density
    self%myconfig%densityHeII(iobject)=self%myconfig%densityHeII(iobject)*self%myphysunit%myconfig%density
    self%myconfig%densityHeIII(iobject)=self%myconfig%densityHeIII(iobject)*self%myphysunit%myconfig%density
    self%myconfig%densityElectrons(iobject)=self%myconfig%densityElectrons(iobject)*self%myphysunit%myconfig%density
    self%myconfig%density_Z(iobject)=self%myconfig%density_Z(iobject)*self%myphysunit%myconfig%density
    self%myconfig%density_dust(iobject)=self%myconfig%density_dust(iobject)*self%myphysunit%myconfig%density

    self%myconfig%normalize_done(iobject)=.false.

  end do
end subroutine grackle_denormalize




!-------------------------------------------------------------------------
!> Read grackle parameters  from a parfile
subroutine grackle_config_read(grackle_config,files,self)
  implicit none
  class(gr_solver)            :: self
  character(len=*),intent(in)        :: files(:)
  type(gr_config), intent(out)  :: grackle_config

  ! .. local ..
  integer                            :: i_file,i_error_read
  integer                  :: idim,iside
  character(len=70)            :: error_message
  logical                  :: grackle_and_usr_boundary
  !-------------------------------------------------------------------------
  namelist /grackle_conf_list/  grackle_config
  namelist /grackle_conf1_list/ grackle_config
  namelist /grackle_conf2_list/ grackle_config
  namelist /grackle_conf3_list/ grackle_config


  error_message = 'In the procedure : grackle_config_read'

  if(mype==0)write(*,*)'Reading grackle_conf_list'
  do i_file = 1, size(files)
     open(unitpar, file=trim(files(i_file)), status="old")
     select case(grackle_config%myindice(1))
     case(1)
       read(unitpar, grackle_conf1_list, iostat=i_error_read)
     case(2)
       read(unitpar, grackle_conf2_list, iostat=i_error_read)
     case(3)
       read(unitpar, grackle_conf3_list, iostat=i_error_read)
     case default
       read(unitpar, grackle_conf_list, iostat=i_error_read)
       !write(*,*) 'use_grackle : ', grackle_config%use_grackle(1)
       !write(*,*) 'gr_primordial_chemistry : ', grackle_config%gr_primordial_chemistry(1)
       !write(*,*) 'gr_metal_cooling : ', grackle_config%gr_metal_cooling(1)
       !write(*,*) 'gr_dust_chemistry : ', grackle_config%gr_dust_chemistry(1)
     end select
     call usr_mat_read_error_message(i_error_read,grackle_config%myindice(1),&
                                     'grackle_chemistry_config')
     close(unitpar)

     !add a routine to check that the sum of each density is equal to one times
     ! the total density in ISM AND JET
  end do

  if(grackle_config%boundary_on)then
    self%myboundaries%myconfig%myindice =1
    call self%myboundaries%read_parameters(self%myboundaries%myconfig,files)
  end if

  if(grackle_config%boundary_on)then
      do idim=1,ndim
        do iside=1,2
          grackle_and_usr_boundary = (trim(grackle_config%boundary_cond(idim,iside))/=&
          trim(self%myboundaries%myconfig%boundary_type(idim,iside)))
          if(grackle_and_usr_boundary) then
             self%myboundaries%myconfig%boundary_type(idim,iside)=&
             grackle_config%boundary_cond(idim,iside)
          end if
        end do
      end do
    end if

  !write(*,*) 'Grackle usr configuration reading successfully done !'
end subroutine grackle_config_read

subroutine grackle_fields_config_allocate(number_of_isms,number_of_jets,number_of_clouds,&
gr_patches_name,gr_patches_indices_global,gr_patches_indices_local,gr_profiles,&
gr_epsilon_tol,gr_density_method)
 implicit none

 integer,intent(inout)        :: number_of_isms
 integer,intent(inout)        :: number_of_jets
 integer,intent(inout)        :: number_of_clouds
 !for namelist:
 CHARACTER(LEN=30) :: gr_patches_name(max_num_parameters)
 ! indices of the objects to be initialized with chemical species
 INTEGER :: gr_patches_indices_global(max_num_parameters)
 INTEGER :: gr_patches_indices_local(max_num_parameters)
 CHARACTER(LEN=64) :: gr_profiles(max_num_parameters)
 ! fields config in the fraction of ISM/JET n_H
 real(kind=gr_rpknd) :: gr_epsilon_tol(max_num_parameters) !density fraction tolerance
 CHARACTER(LEN=30) :: gr_density_method(max_num_parameters)

 ! Fields parameters
 ! parts to be initialized with chemical species
 ! 'ism' : ISM
 ! 'jet' : jet (on top of the ISM)
 integer :: i_zero_ism,i_zero_jet,i_zero_cloud
 integer :: i_end_ism,i_end_jet,i_end_cloud

 ! .. local ..
 integer                  :: n_total_objects
 integer                  :: i_ism,i_jet,i_cloud,i_all,counter,i_patch
 !-------------------------------------------------------------------------

 ! Defaults
 !default_number_of_isms = 1
 !default_number_of_jets = 0
 !default_number_of_clouds = 0
 !if(.not.present(number_of_isms))number_of_isms = default_number_of_isms
 !if(.not.present(number_of_jets))number_of_jets = default_number_of_jets
 !if(.not.present(number_of_clouds))number_of_clouds = default_number_of_clouds


 n_total_objects = 0
 i_zero_ism = 1
 i_zero_jet = 1
 i_zero_cloud = 1
 i_end_ism = 1
 i_end_jet = 1
 i_end_cloud = 1

 ! First, count the number of objects
 ! ex of slicing :
 ! 1 -- ISM --> 4; 5-- JET-->7; 8-- CLOUDS--> 9
 if(number_of_isms>0)then
   i_zero_ism = i_zero_ism ! ex: i_zero_ism = 1
   n_total_objects = n_total_objects + number_of_isms ! ex : n_tot = 0 + 4 = 4
   i_end_ism = n_total_objects
 end if
 if(number_of_jets>0)then
   i_zero_jet = i_zero_ism + number_of_isms ! ex: i_zero_jet = 1 + 4 = 5
   n_total_objects = n_total_objects + number_of_jets ! ex : n_tot = 0 + 4 + 3 = 7
   i_end_jet = n_total_objects
 end if
 if(number_of_clouds>0)then
   i_zero_cloud = i_zero_ism + number_of_isms + number_of_jets
   ! ex: i_zero_cloud = 1 + 4 + 3 = 8
   n_total_objects = n_total_objects + number_of_clouds ! ex : n_tot = 0 + 4 + 3 + 2 = 9
   i_end_cloud = n_total_objects
 end if
 if(n_total_objects==0)then
  WRITE(*,*) 'Error, the code found 0 objects ism+jets+clouds, n_total_objects = ', n_total_objects
  call mpistop('Incoherence in the number of objects. The code stops.')
 end if


 ! Default reference value of density is 0.0d0
 !grackle_ref_density(1:n_total_objects) = 0.0d0

 ! Fourth, fill in specifically for ism

 counter = 0
 ism_is_present :if(number_of_isms>0)then
   Loop_fill_ism : do i_ism=i_zero_ism,i_end_ism

     counter=counter+1
     gr_patches_name(i_ism) = 'ism'
     gr_patches_indices_global(i_ism) = i_ism
     gr_patches_indices_local(i_ism) = counter


   end do Loop_fill_ism
 end if ism_is_present

 ! do the same  specifically for jet
 counter = 0
 jet_is_present :if(number_of_jets>0)then
   Loop_fill_jet : do i_jet=i_zero_jet,i_end_jet

     counter=counter+1
     gr_patches_name(i_jet) = 'jet'
     gr_patches_indices_global(i_jet) = i_jet
     gr_patches_indices_local(i_jet) = counter


   end do Loop_fill_jet
 end if jet_is_present


 ! do the same  specifically for cloud
 counter = 0
 cloud_is_present :if(number_of_clouds>0)then
   Loop_fill_cloud : do i_cloud=i_zero_cloud,i_end_cloud

     counter=counter+1
     gr_patches_name(i_cloud) = 'cloud'
     gr_patches_indices_global(i_cloud) = i_cloud
     gr_patches_indices_local(i_cloud) = counter


   end do Loop_fill_cloud
 end if cloud_is_present

end subroutine grackle_fields_config_allocate


subroutine grackle_write_setting(self,unit_config,gr_patches_name,&
gr_patches_indices_global,gr_patches_indices_local,gr_profiles,&
gr_epsilon_tol,gr_density_method)
  implicit none
  class(gr_solver)            :: self
  integer,intent(in)                  :: unit_config
  CHARACTER(LEN=30) :: gr_patches_name(max_num_parameters)
  ! indices of the objects to be initialized with chemical species
  INTEGER :: gr_patches_indices_global(max_num_parameters)
  INTEGER :: gr_patches_indices_local(max_num_parameters)
  CHARACTER(LEN=64) :: gr_profiles(max_num_parameters)
  ! fields config in the fraction of ISM/JET n_H
  real(kind=gr_rpknd) :: gr_epsilon_tol(max_num_parameters) !density fraction tolerance
  CHARACTER(LEN=30) :: gr_density_method(max_num_parameters)
  integer                             :: idims2,iside2,iB2
  real(kind=dp)                       :: rto_print
  character(len=64)                   :: sto_print
  character(len=128)                   :: wto_print
  integer                             :: idim,iside,idims,iw2,i_all
  ! .. local ..

  !-----------------------------------

  if(self%myconfig%use_grackle(1)==1)then
    write(unit_config,*)'      ****** Grackle Parameters ************* '
    write(unit_config,*)'- In Physical Unit'

    write(unit_config,*)'* Grackle parameters'
    write(unit_config,*) 'use_grackle =',  self%myconfig%use_grackle(1)

    write(unit_config,*)'> Grackle configuration'
    write(unit_config,*) 'Unit :',  self%myconfig%unit(1)
    write(unit_config,*) 'gr_with_radiative_cooling = ',  self%myconfig%gr_with_radiative_cooling(1)
    write(unit_config,*) 'gr_primordial_chemistry = ',  self%myconfig%gr_primordial_chemistry(1)
    write(unit_config,*) 'gr_metal_cooling = ',  self%myconfig%gr_metal_cooling(1)
    write(unit_config,*) 'gr_UVbackground = ',  self%myconfig%gr_UVbackground(1)
    write(unit_config,*) 'gr_cmb_temperature_floor = ',  self%myconfig%gr_cmb_temperature_floor(1)
    write(unit_config,*) 'gr_h2_on_dust = ',  self%myconfig%gr_h2_on_dust(1)
    write(unit_config,*) 'gr_dust_chemistry = ',  self%myconfig%gr_dust_chemistry(1)
    write(unit_config,*) 'gr_use_dust_density_field = ',  self%myconfig%gr_use_dust_density_field(1)
    write(unit_config,*) 'gr_photoelectric_heating = ',  self%myconfig%gr_photoelectric_heating(1)
    write(unit_config,*) 'gr_use_volumetric_heating_rate = ',  self%myconfig%gr_use_volumetric_heating_rate(1)
    write(unit_config,*) 'gr_use_specific_heating_rate = ',  self%myconfig%gr_use_specific_heating_rate(1)
    write(unit_config,*) 'gr_three_body_rate = ',  self%myconfig%gr_three_body_rate(1)
    write(unit_config,*) 'gr_cie_cooling = ',  self%myconfig%gr_cie_cooling(1)
    write(unit_config,*) 'gr_h2_optical_depth_approximation = ',  self%myconfig%gr_h2_optical_depth_approximation(1)
    write(unit_config,*) 'gr_ih2co = ',  self%myconfig%gr_ih2co(1)
    write(unit_config,*) 'gr_ipiht = ',  self%myconfig%gr_ipiht(1)
    write(unit_config,*) 'gr_NumberOfTemperatureBins = ',  self%myconfig%gr_NumberOfTemperatureBins(1)
    write(unit_config,*) 'gr_CaseBRecombination = ',  self%myconfig%gr_CaseBRecombination(1)
    write(unit_config,*) 'gr_Compton_xray_heating = ',  self%myconfig%gr_Compton_xray_heating(1)
    write(unit_config,*) 'gr_LWbackground_sawtooth_suppression = ',  self%myconfig%gr_LWbackground_sawtooth_suppression(1)
    write(unit_config,*) 'gr_NumberOfDustTemperatureBins = ',  self%myconfig%gr_NumberOfDustTemperatureBins(1)
    write(unit_config,*) 'gr_use_radiative_transfer = ',  self%myconfig%gr_use_radiative_transfer(1)
    write(unit_config,*) 'gr_radiative_transfer_coupled_rate_solver = ',  self%myconfig%gr_radiative_transfer_coupled_rate_solver(1)
    write(unit_config,*) 'gr_radiative_transfer_intermediate_step = ',  self%myconfig%gr_radiative_transfer_intermediate_step(1)
    write(unit_config,*) 'gr_radiative_transfer_hydrogen_only = ',  self%myconfig%gr_radiative_transfer_hydrogen_only(1)
    write(unit_config,*) 'gr_self_shielding_method = ',  self%myconfig%gr_self_shielding_method(1)
    write(unit_config,*) 'gr_use_isrf_field = ',  self%myconfig%gr_use_isrf_field(1)
    write(unit_config,*) 'gr_Gamma = ',  self%myconfig%gr_Gamma(1)
    write(unit_config,*) 'gr_photoelectric_heating_rate = ',  self%myconfig%gr_photoelectric_heating_rate(1)

    write(unit_config,*) 'TemperatureStart = ',  self%myconfig%TemperatureStart(1)
    write(unit_config,*) 'TemperatureEnd = ',  self%myconfig%TemperatureEnd(1)
    write(unit_config,*) 'DustTemperatureStart = ',  self%myconfig%DustTemperatureStart(1)
    write(unit_config,*) 'DustTemperatureEnd = ',  self%myconfig%DustTemperatureEnd(1)
    write(unit_config,*) 'gr_LWbackground_intensity = ',  self%myconfig%gr_LWbackground_intensity(1)
    write(unit_config,*) 'gr_UVbackground_redshift_on = ',  self%myconfig%gr_UVbackground_redshift_on(1)
    write(unit_config,*) 'gr_UVbackground_redshift_off = ',  self%myconfig%gr_UVbackground_redshift_off(1)
    write(unit_config,*) 'gr_UVbackground_redshift_fullon = ',  self%myconfig%gr_UVbackground_redshift_fullon(1)
    write(unit_config,*) 'gr_UVbackground_redshift_drop = ',  self%myconfig%gr_UVbackground_redshift_drop(1)
    write(unit_config,*) 'cloudy_electron_fraction_factor = ',  self%myconfig%cloudy_electron_fraction_factor(1)
    write(unit_config,*) 'data_dir = ',  self%myconfig%data_dir(1)
    write(unit_config,*) 'data_filename = ',  self%myconfig%data_filename(1)
    write(unit_config,*) 'data_file = ',  trim(self%myconfig%data_dir(1)//self%myconfig%data_filename(1))
    !write(unit_config,*) 'normalize_done = ',  self%myconfig%normalize_done(1)
    write(unit_config,*) 'gr_comoving_coordinates = ',  self%myconfig%gr_comoving_coordinates(1)
    write(unit_config,*) 'gr_a_units = ',  self%myconfig%gr_a_units(1)
    write(unit_config,*) 'gr_density_units = ',  self%myconfig%gr_density_units(1)
    write(unit_config,*) 'gr_length_units = ',  self%myconfig%gr_length_units(1)
    write(unit_config,*) 'gr_time_units = ',  self%myconfig%gr_time_units(1)
    write(unit_config,*) 'gr_current_redshift = ',  self%myconfig%gr_current_redshift(1)

    write(unit_config,*)'> Physical parameters'
    write(unit_config,*) 'chi_dust = ',  self%myconfig%chi_dust(1)
    write(unit_config,*) 'xi_dust = ',  self%myconfig%xi_dust(1)
    write(unit_config,*) 'HydrogenFractionByMass = ',  self%myconfig%HydrogenFractionByMass(1)
    write(unit_config,*) 'HeliumFractionByMass = ',  self%myconfig%HeliumFractionByMass(1)
    write(unit_config,*) 'MetalFractionByMass = ',  self%myconfig%MetalFractionByMass(1)

    write(unit_config,*) 'SolarMetalFractionByMass = ',  self%myconfig%SolarMetalFractionByMass(1)

    write(unit_config,*) 'DeuteriumToHydrogenRatio = ',  self%myconfig%DeuteriumToHydrogenRatio(1)
    write(unit_config,*) 'IonizationFraction = ',  self%myconfig%IonizationFraction(1)
    write(unit_config,*) 'He_abundance = ',  self%myconfig%He_abundance(1)
    write(unit_config,*) 'Zeta_nH_to_rho = ',  self%myconfig%Zeta_nH_to_rho(1)
    write(unit_config,*) 'ZetaPrime_nH_to_rho = ',  self%myconfig%ZetaPrime_nH_to_rho(1)
    write(unit_config,*) 'Varsigma_nH_to_rho = ',  self%myconfig%Varsigma_nH_to_rho(1)
    write(unit_config,*) 'VarsigmaPrime_nH_to_rho = ',  self%myconfig%VarsigmaPrime_nH_to_rho(1)

    write(unit_config,*) 'correctdensities = ',  self%myconfig%correctdensities(1)
    write(unit_config,*) 'deviation_to_density = ',  self%myconfig%deviation_to_density(1), ' %'
    write(unit_config,*) 'maximum deviation_to_density allowed = ',  self%myconfig%deviation_to_density_limit(1), ' %'


    write(unit_config,*)'>> ',self%myconfig%number_of_solved_species(1),' solved species abundances : '


    do i_all=1,self%number_of_objects
    write(unit_config,*)'      ****** Parameters for ',gr_patches_name(gr_patches_indices_global(i_all)),&
    '#',gr_patches_indices_local(gr_patches_indices_global(i_all)),'******      '

    !Species abundances
    write(unit_config,*) 'Abundances : '

        write(unit_config,*)' + [HI/H] : ' , self%myconfig%x_HI(i_all)
        write(unit_config,*)' + [HII/H] : ' , self%myconfig%x_HII(i_all)
        write(unit_config,*)' + [HeI/H] : ' , self%myconfig%x_HeI(i_all)
        write(unit_config,*)' + [HeII/H] : ' , self%myconfig%x_HeII(i_all)
        write(unit_config,*)' + [HeIII/H] : ' , self%myconfig%x_HeIII(i_all)
        write(unit_config,*)' + [e-/H] : ' , self%myconfig%x_e(i_all)

        write(unit_config,*)' + [HM/H] : ' , self%myconfig%x_HM(i_all)
        write(unit_config,*)' + [H2I/H] : ' , self%myconfig%x_H2I(i_all)
        write(unit_config,*)' + [H2II/H] : ' , self%myconfig%x_H2II(i_all)


        write(unit_config,*)' + [DI/H] : ' , self%myconfig%x_DI(i_all)
        write(unit_config,*)' + [DII/H] : ' , self%myconfig%x_DII(i_all)
        write(unit_config,*)' + [HDI/H] : ' , self%myconfig%x_HDI(i_all)

        !Species densities
        write(unit_config,*) 'Densities : '

        write(unit_config,*)' + rho[HI] : ' , self%myconfig%densityHI(i_all)
        write(unit_config,*)' + rho[HII] : ' , self%myconfig%densityHII(i_all)
        write(unit_config,*)' + rho[HeI] : ' , self%myconfig%densityHeI(i_all)
        write(unit_config,*)' + rho[HeII] : ' , self%myconfig%densityHeII(i_all)
        write(unit_config,*)' + rho[HeIII] : ' , self%myconfig%densityHeIII(i_all)
        write(unit_config,*)' + rho[e-] : ' , self%myconfig%densityelectrons(i_all)
        write(unit_config,*)' + rho[HM] : ' , self%myconfig%densityHM(i_all)
        write(unit_config,*)' + rho[H2I] : ' , self%myconfig%densityH2I(i_all)
        write(unit_config,*)' + rho[H2II] : ' , self%myconfig%densityH2II(i_all)


        write(unit_config,*)' + rho[DI] : ' , self%myconfig%densityDI(i_all)
        write(unit_config,*)' + rho[DII] : ' , self%myconfig%densityDII(i_all)
        write(unit_config,*)' + rho[HDI] : ' , self%myconfig%densityHDI(i_all)


        write(unit_config,*) 'density = ',  self%myconfig%density(i_all)
        write(unit_config,*) 'number_H_density = ',  self%myconfig%number_H_density(i_all)
        write(unit_config,*) 'density_dust = ',  self%myconfig%density_dust(i_all)
        write(unit_config,*) 'density_tot = ',  self%myconfig%density_tot(i_all)
        write(unit_config,*) 'density_bar = ',  self%myconfig%density_bar(i_all)
        write(unit_config,*) 'density_X = ',  self%myconfig%density_X(i_all)
        write(unit_config,*) 'density_Y = ',  self%myconfig%density_Y(i_all)
        write(unit_config,*) 'density_Z = ',  self%myconfig%density_Z(i_all)
        write(unit_config,*) 'density_deut = ',  self%myconfig%density_deut(i_all)
        write(unit_config,*) 'density_not_deut = ',  self%myconfig%density_not_deut(i_all)
        write(unit_config,*) 'density_metal_free = ',  self%myconfig%density_metal_free(i_all)
        write(unit_config,*) 'densityD = ',  self%myconfig%densityD(i_all)
        write(unit_config,*) 'densityHD = ',  self%myconfig%densityHD(i_all)
        write(unit_config,*) 'densityH = ',  self%myconfig%densityH(i_all)
        write(unit_config,*) 'densityHtwo = ',  self%myconfig%densityHtwo(i_all)
        write(unit_config,*) 'densityHe = ',  self%myconfig%densityHe(i_all)
        write(unit_config,*) 'densityElectrons = ',  self%myconfig%densityElectrons(i_all)
    ! fields config
    write(unit_config,*)'> Objects parameters'
    write(unit_config,*) 'gr_patches_indices_global = ',gr_patches_indices_global
    write(unit_config,*) 'gr_patches_indices_local = ',gr_patches_indices_local
    write(unit_config,*) 'gr_patches_name = ',gr_patches_name
    write(unit_config,*) 'gr_epsilon_tol = ',gr_epsilon_tol
    write(unit_config,*) 'gr_profiles = ',gr_profiles
    write(unit_config,*) 'gr_density_method = ',gr_density_method

    end do

    write(unit_config,*) '======================================================='
    write(unit_config,*)'************************************'
    write(unit_config,*)'******** END of Grackle setting **********'
    write(unit_config,*)'************************************'
  else
    write(unit_config,*)'      ****** Grackle disabled************************'
  end if


end    subroutine grackle_write_setting






!I> COMPLETE EACH OBJECTS GRACKLE PARAMETERS

!--------------------------------------------------------------------
!> subroutine check the parfile setting for ism
subroutine grackle_set_complet(self,ref_density,iobject,normalized,mydensityunit)
  implicit none
  class(gr_solver)                                :: self
  real(kind=dp),intent(in) :: ref_density
  integer, intent(in)             :: iobject
  logical       :: normalized
  real(kind=dp),optional :: mydensityunit
  ! .. local ..
  !integer :: iobject
  real(kind=dp) :: physical_ref_density
  !-----------------------------------
  if(normalized.and.(.not.present(mydensityunit)))then
    call mpistop('normalized==true but no density unit given')
  end if

  physical_ref_density = ref_density !rho

  if(normalized)then
    physical_ref_density = ref_density*&
    mydensityunit
  end if

  write(*,*) ' set_complet iobject begining of set_complet =', iobject


  ! == rho = rhoH + rhoHe + rhoZ
  ! Z
  !if(self%myconfig%MetalFractionByMass(iobject)<smalldouble)then
    !self%myconfig%MetalFractionByMass(iobject)=self%myconfig%SolarMetalFractionByMass(iobject)
  !end if

  ! rhoZ = Z*rho
  self%myconfig%density_Z(iobject)=self%myconfig%MetalFractionByMass(iobject)*&
  physical_ref_density

  ! chi_dust and xi_dust
  if((self%myconfig%chi_dust(iobject)<smalldouble).and.&
  (self%myconfig%xi_dust(iobject)>=smalldouble))then
      self%myconfig%chi_dust(iobject)=DABS(self%myconfig%xi_dust(iobject))/&
      (self%myconfig%MetalFractionByMass(iobject)/self%myconfig%SolarMetalFractionByMass(iobject))
  elseif(self%myconfig%chi_dust(iobject)>=smalldouble)then
      self%myconfig%xi_dust(iobject)=DABS(self%myconfig%chi_dust(iobject))*&
      (self%myconfig%MetalFractionByMass(iobject)/self%myconfig%SolarMetalFractionByMass(iobject))
  elseif((self%myconfig%chi_dust(iobject)<smalldouble).and.&
  (self%myconfig%xi_dust(iobject)<smalldouble))then
    self%myconfig%xi_dust(iobject)=0.0d0
    self%myconfig%chi_dust(iobject)=0.0d0
  end if

  !rhodust = chiDust * rho
  self%myconfig%density_dust(iobject)=physical_ref_density*&
  self%myconfig%chi_dust(iobject)


  !rhoHplusH2 = rhoH
  self%myconfig%densityHplusH2(iobject)=self%myconfig%HydrogenFractionByMass(iobject)*&
  physical_ref_density

  !Y
  self%myconfig%HeliumFractionByMass(iobject) = 1.0_dp - &
  self%myconfig%HydrogenFractionByMass(iobject)


  !rhoY
  self%myconfig%densityHe(iobject) = self%myconfig%HeliumFractionByMass(iobject)*&
  physical_ref_density

  self%myconfig%density_Y(iobject)=self%myconfig%densityHe(iobject)

  !rhodeut
  self%myconfig%density_deut(iobject) = self%myconfig%DeuteriumToHydrogenRatio(iobject)*&
  self%myconfig%densityHplusH2(iobject)

  !rhoX = rhoD + rhoHD + rhoH + rhoH2
  self%myconfig%density_X(iobject) = self%myconfig%density_deut(iobject)+&
  self%myconfig%densityHplusH2(iobject)


  self%myconfig%density_bar(iobject) = self%myconfig%density_X(iobject)+&
  self%myconfig%density_Y(iobject)+ self%myconfig%density_Z(iobject)

    !Finalyy compute values from abundances
    !x_DI + x_DII +x_HDI = 1
    self%myconfig%densityDI(iobject) = self%myconfig%x_DI(iobject)*&
    self%myconfig%density_deut(iobject)
    self%myconfig%densityDII(iobject) = self%myconfig%x_DII(iobject)*&
    self%myconfig%density_deut(iobject)
    self%myconfig%densityHDI(iobject) = self%myconfig%x_HDI(iobject)*&
    self%myconfig%density_deut(iobject)
    ! x_HI + x_HII + x_HM + x_H2I + x_H2II = 1
    self%myconfig%densityHI(iobject) = self%myconfig%x_HI(iobject)*&
    self%myconfig%densityHplusH2(iobject)
    self%myconfig%densityHII(iobject) = self%myconfig%x_HII(iobject)*&
    self%myconfig%densityHplusH2(iobject)
    self%myconfig%densityHM(iobject) = self%myconfig%x_HM(iobject)*&
    self%myconfig%densityHplusH2(iobject)
    self%myconfig%densityH2I(iobject) = self%myconfig%x_H2I(iobject)*&
    self%myconfig%densityHplusH2(iobject)
    self%myconfig%densityH2II(iobject) = self%myconfig%x_H2II(iobject)*&
    self%myconfig%densityHplusH2(iobject)
    ! x_HeI + x_HeII + x_HeIII = 1
    self%myconfig%densityHeI(iobject) = self%myconfig%x_HeI(iobject)*&
    self%myconfig%density_Y(iobject)
    self%myconfig%densityHeII(iobject) = self%myconfig%x_HeII(iobject)*&
    self%myconfig%density_Y(iobject)
    self%myconfig%densityHeIII(iobject) = self%myconfig%x_HeIII(iobject)*&
    self%myconfig%density_Y(iobject)
    ! rhoD = rhoDI + rhoDII
    self%myconfig%densityD(iobject)=self%myconfig%densityDI(iobject)+&
    self%myconfig%densityDII(iobject)
    ! rhoHD
    self%myconfig%densityHD(iobject)=self%myconfig%densityHDI(iobject)
    !rhonotdeut
    self%myconfig%density_not_deut(iobject) = self%myconfig%densityHI(iobject)+&
    self%myconfig%densityHII(iobject) + self%myconfig%densityHM(iobject)+&
    self%myconfig%densityH2I(iobject) + self%myconfig%densityH2II(iobject)+&
    self%myconfig%densityHeI(iobject) + self%myconfig%densityHeII(iobject)+&
    self%myconfig%densityHeIII(iobject) + &
    self%myconfig%density_Z(iobject)
    !same for rhoH and rhoH2
    self%myconfig%densityH(iobject)=self%myconfig%densityHI(iobject) +&
    self%myconfig%densityHII(iobject) + self%myconfig%densityHM(iobject)
    self%myconfig%densityHtwo(iobject)=  self%myconfig%densityH2I(iobject) +&
      self%myconfig%densityH2II(iobject)

    !nH
    self%myconfig%number_H_density(iobject) = &
    self%myconfig%densityHI(iobject)/(w_HI*mh_gr)+&
    self%myconfig%densityHII(iobject)/(w_HII*mh_gr)+&
    self%myconfig%densityHM(iobject)/(w_HM*mh_gr)+&
    2.0_dp * self%myconfig%densityDI(iobject)/(w_DI*mh_gr)+&
    2.0_dp * self%myconfig%densityDII(iobject)/(w_DII*mh_gr)+&
    2.0_dp * self%myconfig%densityH2I(iobject)/(w_H2I*mh_gr)+&
    2.0_dp * self%myconfig%densityH2II(iobject)/(w_H2II*mh_gr)+&
    3.0_dp * self%myconfig%densityHDI(iobject)/(w_HDI*mh_gr)

    !rho_Electrons
    self%myconfig%densityElectrons(iobject) = w_e *&
    self%myconfig%IonizationFraction(iobject)*& !checked
    self%myconfig%number_H_density(iobject)*mh_gr

    !fix rhotot
    self%myconfig%density_tot(iobject) = &
    self%myconfig%densityHI(iobject)+&
    self%myconfig%densityHII(iobject)+&
    self%myconfig%densityHM(iobject)+&
    self%myconfig%densityDI(iobject)+&
    self%myconfig%densityDII(iobject)+&
    self%myconfig%densityH2I(iobject)+&
    self%myconfig%densityH2II(iobject)+&
    self%myconfig%densityHDI(iobject)+&
    self%myconfig%densityHeI(iobject)+&
    self%myconfig%densityHeII(iobject)+&
    self%myconfig%densityHeIII(iobject)+&
    self%myconfig%densityElectrons(iobject)+&
    self%myconfig%density_Z(iobject)+&
    self%myconfig%density_dust(iobject)

end subroutine grackle_set_complet


subroutine grackle_chemistry_static_source(ixI^L,ixO^L,x,qdt,qtC,wCT,qt,w,grid_dx22,&
solver_activated,pressure_field,gamma_field,temperature_field,self)
  implicit none
  integer, intent(in)                     :: ixI^L,ixO^L
  real(kind=dp), intent(in)               :: qdt,qtC,qt
  real(kind=dp), intent(in)               :: x(ixI^S,1:ndim)
  real(kind=dp), intent(in)               :: wCT(ixI^S,1:nw)
  real(kind=dp), intent(in)               :: grid_dx22(1:ndim)
  logical, intent(in) :: solver_activated
  real(kind=dp), intent(inout),optional   :: pressure_field(ixI^S)
  real(kind=dp), intent(inout),optional   :: temperature_field(ixI^S)
  real(kind=dp), intent(inout),optional   :: gamma_field(ixI^S)
  !type(usrphysical_unit), target     :: physunit_inuse
  !real(kind=dp), intent(in)               :: gamma
  real(kind=dp),intent(inout)             :: w(ixI^S,1:nw)
  class(gr_solver),TARGET                        :: self
  !type(gr_objects), intent(inout)  :: gr_obj
  ! .. local ..
  integer                                 :: idim,iresult,ifield^D,imesh^D,level,iw,idir,i
  character(len=150), TARGET :: filename
  real(kind=dp) :: temperature_units, pressure_units, dtchem, velocity_units
  type(grackle_field_data)                :: my_gr_fields
  real(kind=dp), target                   :: w_gr(ixI^S,1:nw)
  real(kind=dp)                   :: temperature(ixI^S)
  real(kind=dp)                           :: dt_old
  integer         :: field_size(1:ndim)
  real :: initial_redshift
  real(kind=dp),parameter :: fH = 0.76
  real(kind=dp), TARGET :: density({1:ixOmax^D-ixOmin^D+1|,}), energy({1:ixOmax^D-ixOmin^D+1|,}), &
  x_velocity({1:ixOmax^D-ixOmin^D+1|,}), y_velocity({1:ixOmax^D-ixOmin^D+1|,}), &
  z_velocity({1:ixOmax^D-ixOmin^D+1|,}), &
  HI_density({1:ixOmax^D-ixOmin^D+1|,}), HII_density({1:ixOmax^D-ixOmin^D+1|,}), &
  HM_density({1:ixOmax^D-ixOmin^D+1|,}), &
  HeI_density({1:ixOmax^D-ixOmin^D+1|,}), HeII_density({1:ixOmax^D-ixOmin^D+1|,}), &
  HeIII_density({1:ixOmax^D-ixOmin^D+1|,}), &
  H2I_density({1:ixOmax^D-ixOmin^D+1|,}), H2II_density({1:ixOmax^D-ixOmin^D+1|,}), &
  DI_density({1:ixOmax^D-ixOmin^D+1|,}), DII_density({1:ixOmax^D-ixOmin^D+1|,}), &
  HDI_density({1:ixOmax^D-ixOmin^D+1|,}), &
  e_density({1:ixOmax^D-ixOmin^D+1|,}), metal_density({1:ixOmax^D-ixOmin^D+1|,}), &
  dust_density({1:ixOmax^D-ixOmin^D+1|,}), &
  volumetric_heating_rate({1:ixOmax^D-ixOmin^D+1|,}), &
  specific_heating_rate({1:ixOmax^D-ixOmin^D+1|,}), &
  RT_HI_ionization_rate({1:ixOmax^D-ixOmin^D+1|,}), &
  RT_HeI_ionization_rate({1:ixOmax^D-ixOmin^D+1|,}), &
  RT_HeII_ionization_rate({1:ixOmax^D-ixOmin^D+1|,}), &
  RT_H2_dissociation_rate({1:ixOmax^D-ixOmin^D+1|,}), &
  RT_heating_rate({1:ixOmax^D-ixOmin^D+1|,})
  real(kind=dp), TARGET  :: pressure({1:ixOmax^D-ixOmin^D+1|,}),&
  gamma({1:ixOmax^D-ixOmin^D+1|,}),temperature_gr({1:ixOmax^D-ixOmin^D+1|,})
  INTEGER, TARGET :: grid_rank, grid_dimension(3), grid_start(3), grid_end(3)
  logical :: logicndir,logicenergy
  !real(kind=dp)           :: gamma_amrvac(ixI^S)
  real(kind=dp) :: grid_dx
  real(kind=dp) :: my_temperature(1),&
  my_gamma(1),my_pressure(1),density_gr(1),&
           dust_density_gr(1),&
           metal_density_gr(1),&
           e_density_gr(1),&
           DI_density_gr(1),&
           DII_density_gr(1),&
           HDI_density_gr(1),&
           HI_density_gr(1),&
           HII_density_gr(1),&
           HM_density_gr(1),&
           H2I_density_gr(1),&
           H2II_density_gr(1),&
           HeI_density_gr(1),&
           HeII_density_gr(1),&
           HeIII_density_gr(1)
  real(kind=dp) :: current_redshift
  TYPE (grackle_units) :: my_units
  TYPE (grackle_field_data) :: my_fields
  TYPE (grackle_chemistry_data) :: grackle_data
  !-----------------------------------------------------------------------------------
  dt_old = dt

  iresult = set_default_chemistry_parameters(grackle_data)
  call self%link_par_to_gr(grackle_data,my_units)
  iresult = initialize_chemistry_data(my_units)


  grackle_data%TemperatureStart=100.;
  !write(*,*) '===================================='
  !write(*,*) 'grackle_chemistry_static_source : '
  !write(*,*) 'grackle_data%TemperatureStart == ', grackle_data%TemperatureStart
  !write(*,*) '===================================='


  !Set field arrays
  {field_size(^D) = ixOmax^D-ixOmin^D+1|\}
  !If grid rank is less than 3, set the other dimensions,
  !start indices, and end indices to 0.
  grid_rank = 3
  do i = 1, grid_rank
     grid_dimension(i) = 1
     grid_start(i) = 0
     grid_end(i) = 0
  enddo
  grid_dx = grid_dx22(1)*unit_length/my_units%length_units

  {grid_dimension(^D) = field_size(^D)
  !0-based is necessary for grid_start and grid_end
  grid_end(^D) = field_size(^D)-1|\}
  temperature_units = get_temperature_units(my_units)
  velocity_units = get_velocity_units(my_units)

  w_gr=w
  call phys_get_temperature(ixI^L,ixI^L,w_gr,x,temperature)


  !do iw=1,nw
    !logicndir=.false.
    !do idir=1,ndir
      !logicndir = logicndir .or. (iw==phys_ind%mom(idir))
    !end do
    !if(logicndir)then
      ! conservative variable
      !w_gr(ixI^S,iw)=w_gr(ixI^S,iw)*unit_density*unit_velocity
    !end if
  !end do

  ! conservative variable
  w_gr(ixI^S,phys_ind%e_)=temperature(ixI^S)*unit_temperature/temperature_units

  !w_gr(ixI^S,phys_ind%e_) = w_gr(ixI^S,phys_ind%e_)*velocity_units**2.0_dp

  !write(*,*) ' temperature_units == ', velocity_units**2.0_dp
  !w(ixO^S,phys_ind%gamma_) = w_gr(ixI^S,phys_ind%e_)*velocity_units**2.0_dp


  w_gr(ixO^S,phys_ind%rho_) = w_gr(ixO^S,phys_ind%rho_)*w_convert_factor(phys_ind%rho_)/&
  my_units%density_units
  w_gr(ixO^S,phys_ind%HI_density_) = w_gr(ixO^S,phys_ind%HI_density_)*w_convert_factor(phys_ind%HI_density_)/&
  my_units%density_units
  w_gr(ixO^S,phys_ind%HII_density_) = w_gr(ixO^S,phys_ind%HII_density_)*w_convert_factor(phys_ind%HII_density_)/&
  my_units%density_units
  w_gr(ixO^S,phys_ind%HM_density_) = w_gr(ixO^S,phys_ind%HM_density_)*w_convert_factor(phys_ind%HM_density_)/&
  my_units%density_units
  w_gr(ixO^S,phys_ind%HeI_density_) = w_gr(ixO^S,phys_ind%HeI_density_)*w_convert_factor(phys_ind%HeI_density_)/&
  my_units%density_units
  w_gr(ixO^S,phys_ind%HeII_density_) = w_gr(ixO^S,phys_ind%HeII_density_)*w_convert_factor(phys_ind%HeII_density_)/&
  my_units%density_units
  w_gr(ixO^S,phys_ind%HeIII_density_) = w_gr(ixO^S,phys_ind%HeIII_density_)*w_convert_factor(phys_ind%HeIII_density_)/&
  my_units%density_units
  w_gr(ixO^S,phys_ind%H2I_density_) = w_gr(ixO^S,phys_ind%H2I_density_)*w_convert_factor(phys_ind%H2I_density_)/&
  my_units%density_units
  w_gr(ixO^S,phys_ind%H2II_density_) = w_gr(ixO^S,phys_ind%H2II_density_)*w_convert_factor(phys_ind%H2II_density_)/&
  my_units%density_units
  w_gr(ixO^S,phys_ind%DI_density_) = w_gr(ixO^S,phys_ind%DI_density_)*w_convert_factor(phys_ind%DI_density_)/&
  my_units%density_units
  w_gr(ixO^S,phys_ind%DII_density_) = w_gr(ixO^S,phys_ind%DII_density_)*w_convert_factor(phys_ind%DII_density_)/&
  my_units%density_units
  w_gr(ixO^S,phys_ind%HDI_density_) =  w_gr(ixO^S,phys_ind%HDI_density_)*w_convert_factor(phys_ind%HDI_density_)/&
  my_units%density_units
  w_gr(ixO^S,phys_ind%e_density_) = w_gr(ixO^S,phys_ind%e_density_)*(mh_gr/me_gr)*w_convert_factor(phys_ind%e_density_)/&
  my_units%density_units
  w_gr(ixO^S,phys_ind%metal_density_) =  w_gr(ixO^S,phys_ind%metal_density_)*w_convert_factor(phys_ind%metal_density_)/&
  my_units%density_units
  w_gr(ixO^S,phys_ind%dust_density_) =  w_gr(ixO^S,phys_ind%dust_density_)*w_convert_factor(phys_ind%dust_density_)/&
  my_units%density_units


  {do ifield^D = 1,field_size(^D)
    imesh^D=ixOmin^D+ifield^D-1|\}
    density(ifield^D) = w_gr(imesh^D,phys_ind%rho_)
    HI_density(ifield^D) = w_gr(imesh^D,phys_ind%HI_density_)
    HII_density(ifield^D) = w_gr(imesh^D,phys_ind%HII_density_)
    HM_density(ifield^D) = w_gr(imesh^D,phys_ind%HM_density_)
    HeI_density(ifield^D) = w_gr(imesh^D,phys_ind%HeI_density_)
    HeII_density(ifield^D) = w_gr(imesh^D,phys_ind%HeII_density_)
    HeIII_density(ifield^D) = w_gr(imesh^D,phys_ind%HeIII_density_)
    H2I_density(ifield^D) = w_gr(imesh^D,phys_ind%H2I_density_)
    H2II_density(ifield^D) = w_gr(imesh^D,phys_ind%H2II_density_)
    DI_density(ifield^D) =  w_gr(imesh^D,phys_ind%DI_density_)
    DII_density(ifield^D) = w_gr(imesh^D,phys_ind%DII_density_)
    HDI_density(ifield^D) = w_gr(imesh^D,phys_ind%HDI_density_)
    e_density(ifield^D) = w_gr(imesh^D,phys_ind%e_density_)
    metal_density(ifield^D) = w_gr(imesh^D,phys_ind%metal_density_)
    if(grackle_data%use_dust_density_field==1)then
      dust_density(ifield^D) = w_gr(imesh^D,phys_ind%dust_density_)
    elseif(grackle_data%use_dust_density_field==0)then
      dust_density(ifield^D) = grackle_data%local_dust_to_gas_ratio*density(ifield^D)
    end if
    energy(ifield^D) = w_gr(imesh^D,phys_ind%e_)
     x_velocity(ifield^D) = 0.0
     y_velocity(ifield^D) = 0.0
     z_velocity(ifield^D) = 0.0
     volumetric_heating_rate(ifield^D) = 0.0
     specific_heating_rate(ifield^D) = 0.0
     RT_HI_ionization_rate(ifield^D) = 0.0
     RT_HeI_ionization_rate(ifield^D) = 0.0
     RT_HeII_ionization_rate(ifield^D) = 0.0
     RT_H2_dissociation_rate(ifield^D) = 0.0
     RT_heating_rate(ifield^D) = 0.0
  {enddo^D&|\}

  !write(*,*) 'Energy field = ', energy

  !Fill in structure to be passed to Grackle

  my_fields%grid_rank = ndim
  my_fields%grid_dimension = C_LOC(grid_dimension)
  my_fields%grid_start = C_LOC(grid_start)
  my_fields%grid_end = C_LOC(grid_end)
  my_fields%grid_dx  = grid_dx

  my_fields%density = C_LOC(density)
  my_fields%HI_density = C_LOC(HI_density)
  my_fields%HII_density = C_LOC(HII_density)
  my_fields%HM_density = C_LOC(HM_density)
  my_fields%HeI_density = C_LOC(HeI_density)
  my_fields%HeII_density = C_LOC(HeII_density)
  my_fields%HeIII_density = C_LOC(HeIII_density)
  my_fields%H2I_density = C_LOC(H2I_density)
  my_fields%H2II_density = C_LOC(H2II_density)
  my_fields%DI_density = C_LOC(DI_density)
  my_fields%DII_density = C_LOC(DII_density)
  my_fields%HDI_density = C_LOC(HDI_density)
  my_fields%e_density = C_LOC(e_density)
  my_fields%metal_density = C_LOC(metal_density)
  my_fields%dust_density = C_LOC(dust_density)
  my_fields%internal_energy = C_LOC(energy)
  my_fields%x_velocity = C_LOC(x_velocity)
  my_fields%y_velocity = C_LOC(y_velocity)
  my_fields%z_velocity = C_LOC(z_velocity)
  my_fields%volumetric_heating_rate = &
                                  C_LOC(volumetric_heating_rate)
  my_fields%specific_heating_rate = C_LOC(specific_heating_rate)
  my_fields%RT_HI_ionization_rate = C_LOC(RT_HI_ionization_rate)
  my_fields%RT_HeI_ionization_rate = C_LOC(RT_HeI_ionization_rate)
  my_fields%RT_HeII_ionization_rate = &
                                    C_LOC(RT_HeII_ionization_rate)
  my_fields%RT_H2_dissociation_rate = C_LOC(RT_H2_dissociation_rate)
  my_fields%RT_heating_rate = C_LOC(RT_heating_rate)

  !Calling the chemistry solver
  !These routines can now be called during the simulation.
  !Evolving the chemistry.


  !write(*,*) " Before chemistry solving, HeI = ", HeI_density*my_units%density_units
  !write(*,*) " Before chemistry solving, e- =", e_density*my_units%density_units
  dtchem = qdt*time_convert_factor / my_units%time_units
  !write(*,*) " dtchem = ", dtchem
  if(solver_activated)then
    iresult = solve_chemistry(my_units, my_fields, dtchem)
  end if
  pressure_units = my_units%density_units * my_units%velocity_units**2.0_dp
  iresult = calculate_pressure(my_units, my_fields, pressure)
  iresult = calculate_gamma(my_units, my_fields, gamma)
  iresult =  calculate_temperature(my_units, my_fields, temperature_gr)

  !write(*,*) 'New temperature = ', temperature_gr


  !write(*,*) " gamma = ", gamma
  !write(*,*) " After chemistry solving, HeI = ", HeI_density*my_units%density_units
  !write(*,*) " After chemistry solving, e- =", e_density*my_units%density_units

  {do ifield^D = 1,field_size(^D)|\}
     {imesh^D=ixOmin^D+ifield^D-1|\}
     w_gr(imesh^D,phys_ind%HI_density_) = HI_density(ifield^D)
     w_gr(imesh^D,phys_ind%HII_density_) = HII_density(ifield^D)
     w_gr(imesh^D,phys_ind%HM_density_) = HM_density(ifield^D)
     w_gr(imesh^D,phys_ind%HeI_density_) = HeI_density(ifield^D)
     w_gr(imesh^D,phys_ind%HeII_density_) = HeII_density(ifield^D)
     w_gr(imesh^D,phys_ind%HeIII_density_) = HeIII_density(ifield^D)
     w_gr(imesh^D,phys_ind%H2I_density_) = H2I_density(ifield^D)
     w_gr(imesh^D,phys_ind%H2II_density_) = H2II_density(ifield^D)
     w_gr(imesh^D,phys_ind%DI_density_) = DI_density(ifield^D)
     w_gr(imesh^D,phys_ind%DII_density_)  = DII_density(ifield^D)
     w_gr(imesh^D,phys_ind%HDI_density_) = HDI_density(ifield^D)
     w_gr(imesh^D,phys_ind%e_density_) = e_density(ifield^D)
     w_gr(imesh^D,phys_ind%metal_density_) = metal_density(ifield^D)
     w_gr(imesh^D,phys_ind%dust_density_) = dust_density(ifield^D)
  {enddo^D&|\}


      !do not forget to scale back the electron field with (me/mH)

      w(ixO^S,phys_ind%HI_density_) = w_gr(ixO^S,phys_ind%HI_density_)*my_units%density_units/&
      w_convert_factor(phys_ind%HI_density_)
      w(ixO^S,phys_ind%HII_density_) = w_gr(ixO^S,phys_ind%HII_density_)*my_units%density_units/&
      w_convert_factor(phys_ind%HII_density_)
      w(ixO^S,phys_ind%HM_density_) = w_gr(ixO^S,phys_ind%HM_density_)*my_units%density_units/&
      w_convert_factor(phys_ind%HM_density_)
      w(ixO^S,phys_ind%HeI_density_) = w_gr(ixO^S,phys_ind%HeI_density_)*my_units%density_units/&
      w_convert_factor(phys_ind%HeI_density_)
      w(ixO^S,phys_ind%HeII_density_) = w_gr(ixO^S,phys_ind%HeII_density_)*my_units%density_units/&
      w_convert_factor(phys_ind%HeII_density_)
      w(ixO^S,phys_ind%HeIII_density_) = w_gr(ixO^S,phys_ind%HeIII_density_)*my_units%density_units/&
      w_convert_factor(phys_ind%HeIII_density_)
      w(ixO^S,phys_ind%H2I_density_) = w_gr(ixO^S,phys_ind%H2I_density_)*my_units%density_units/&
      w_convert_factor(phys_ind%H2I_density_)
      w(ixO^S,phys_ind%H2II_density_) = w_gr(ixO^S,phys_ind%H2II_density_)*my_units%density_units/&
      w_convert_factor(phys_ind%H2II_density_)
      w(ixO^S,phys_ind%DI_density_) = w_gr(ixO^S,phys_ind%DI_density_)*my_units%density_units/&
      w_convert_factor(phys_ind%DI_density_)
      w(ixO^S,phys_ind%DII_density_) = w_gr(ixO^S,phys_ind%DII_density_)*my_units%density_units/&
      w_convert_factor(phys_ind%DII_density_)
      w(ixO^S,phys_ind%HDI_density_) =  w_gr(ixO^S,phys_ind%HDI_density_)*my_units%density_units/&
      w_convert_factor(phys_ind%HDI_density_)
      w(ixO^S,phys_ind%e_density_) = w_gr(ixO^S,phys_ind%e_density_)*(me_gr/mh_gr)*&
      my_units%density_units/w_convert_factor(phys_ind%e_density_)
      w(ixO^S,phys_ind%metal_density_) =  w_gr(ixO^S,phys_ind%metal_density_)*my_units%density_units/&
      w_convert_factor(phys_ind%metal_density_)
      w(ixO^S,phys_ind%dust_density_) =  w_gr(ixO^S,phys_ind%dust_density_)*my_units%density_units/&
      w_convert_factor(phys_ind%dust_density_)
      w(ixO^S,phys_ind%rhoX_)=w(ixO^S,phys_ind%HI_density_)+&
      w(ixO^S,phys_ind%HII_density_)+w(ixO^S,phys_ind%HM_density_)+&
      w(ixO^S,phys_ind%H2I_density_)+w(ixO^S,phys_ind%H2II_density_)+&
      w(ixO^S,phys_ind%DI_density_)+w(ixO^S,phys_ind%DII_density_)+&
      w(ixO^S,phys_ind%HDI_density_)
      w(ixO^S,phys_ind%rhoY_)=w(ixO^S,phys_ind%HeI_density_)+&
      w(ixO^S,phys_ind%HeII_density_)+w(ixO^S,phys_ind%HeIII_density_)


    !do iw=1,nw
      !logicndir=.false.
      !do idir=1,ndim+1
        !logicndir = logicndir .or. (iw==phys_ind%mom(idir))
      !end do
      !if(logicndir)then
        ! conservative variable
        !w_gr(ixI^S,iw)=w_gr(ixI^S,iw)/(unit_density*unit_velocity)
      !end if
    !end do

    if(present(gamma_field))then
    {do ifield^D = 1,field_size(^D)|\}
       {imesh^D=ixOmin^D+ifield^D-1|\}
       ! Beware : in Grackle, rho*internal_energy_gr = pressure/(Gamma-1)
       ! whereas in MPI-AMRVAC : total_energy_amrvac = pressure/(Gamma-1) + 0.5rho*v^2
       ! so that internal_energy_amrvac =  rho*internal_energy_gr
       ! e = p/(gamma-1) = nk_b T/((gamma-1)*mu*mH)
       gamma_field(imesh^D) = gamma(ifield^D)
    {enddo^D&|\}
    end if


    if(present(temperature_field))then
    {do ifield^D = 1,field_size(^D)|\}
       {imesh^D=ixOmin^D+ifield^D-1|\}
       ! Beware : in Grackle, rho*internal_energy_gr = pressure/(Gamma-1)
       ! whereas in MPI-AMRVAC : total_energy_amrvac = pressure/(Gamma-1) + 0.5rho*v^2
       ! so that internal_energy_amrvac =  rho*internal_energy_gr
       ! e = p/(gamma-1) = nk_b T/((gamma-1)*mu*mH)
       temperature_field(imesh^D) = temperature_gr(ifield^D)/unit_temperature
    {enddo^D&|\}
    end if


    !write(*,*) " gamma = ", w(ixO^S,phys_ind%gamma_)


    !denormalize pressure

    ! pth = gamma_1 * (e-kin_en) ==> e = pth/gamma_1+kin_en


    if(present(pressure_field))then
      {do ifield^D = 1,field_size(^D)|\}
         {imesh^D=ixOmin^D+ifield^D-1|\}
          pressure_field(imesh^D) = (w_gr(imesh^D,phys_ind%rho_)*my_units%density_units/&
          w_convert_factor(phys_ind%rho_))*&
          (temperature_field(imesh^D)/(temperature_units*unit_pressure))
          !pressure(ifield^D)*pressure_units/&
          !w_convert_factor(phys_ind%pressure_)
      {enddo^D&|\}
    end if
    if(solver_activated)then
        !{do ifield^D = 1,field_size(^D)|\}
           !{imesh^D=ixOmin^D+ifield^D-1|\}
           !w_gr(imesh^D,phys_ind%e_) = pressure(ifield^D)*pressure_units/&
           !w_convert_factor(phys_ind%pressure_)
        !{enddo^D&|\}

        !w(ixO^S,phys_ind%e_) = (w_gr(ixO^S,phys_ind%e_)/&
        !(grackle_data%gamma-1.0_dp))+&
        !0.5_dp * sum(w_gr(ixO^S,phys_ind%mom(:))**2.0_dp,dim=ndim+1) / &
        !(w_gr(ixO^S,phys_ind%rho_)*my_units%density_units/&
        !w_convert_factor(phys_ind%rho_)) ! AMRVAC cgs denormalized specific energy
    end if

  !call self%normalize()

  !iresult = gr_free_memory(my_units, my_fields)


  dt = dt_old
end subroutine grackle_chemistry_static_source


subroutine grackle_solver_associate(gr_data,myunits,self)
  implicit none
  class(gr_solver)                        :: self
  type(grackle_chemistry_data),intent(inout)        :: gr_data
  TYPE (grackle_units),intent(inout) :: myunits
  integer :: iresult
  !------------------------------
  !     Create a grackle chemistry object for parameters and set defaults




    !     Set parameters

  gr_data%use_grackle = self%myconfig%use_grackle(1)
  gr_data%with_radiative_cooling = self%myconfig%gr_with_radiative_cooling(1)
  gr_data%primordial_chemistry = self%myconfig%gr_primordial_chemistry(1)
  gr_data%dust_chemistry = self%myconfig%gr_dust_chemistry(1)
  gr_data%metal_cooling = self%myconfig%gr_metal_cooling(1)
  gr_data%UVbackground                   = self%myconfig%gr_UVbackground(1)
  general_grackle_filename = "/home/mrabenanahary/MPI-AMRVAC/amrvac/"//&
  "src/grackle/input/CloudyData_UVB=HM2012_shielded.h5"//C_NULL_CHAR
  !CALL GETCWD(filename)
  !write(*,*) 'CWDDDDD = ', filename
  gr_data%grackle_data_file = C_LOC(general_grackle_filename(1:1))
  gr_data%h2_on_dust = self%myconfig%gr_h2_on_dust(1)
  gr_data%cmb_temperature_floor = self%myconfig%gr_cmb_temperature_floor(1)
  gr_data%Gamma = self%myconfig%gr_Gamma(1)
  gr_data%use_dust_density_field = self%myconfig%gr_use_dust_density_field(1)
  gr_data%three_body_rate = self%myconfig%gr_three_body_rate(1)
  gr_data%cie_cooling                    = self%myconfig%gr_cie_cooling(1)
  gr_data%h2_optical_depth_approximation = self%myconfig%gr_h2_optical_depth_approximation(1)
  gr_data%photoelectric_heating = self%myconfig%gr_photoelectric_heating(1)
  ! epsilon=0.05, G_0=1.7 (in erg s -1 cm-3)
  gr_data%photoelectric_heating_rate     = self%myconfig%gr_photoelectric_heating_rate(1)
  gr_data%use_volumetric_heating_rate    = self%myconfig%gr_use_volumetric_heating_rate(1)
  gr_data%use_specific_heating_rate      = self%myconfig%gr_use_specific_heating_rate(1)
  gr_data%UVbackground_redshift_on      = self%myconfig%gr_UVbackground_redshift_on(1)
  gr_data%UVbackground_redshift_off     = self%myconfig%gr_UVbackground_redshift_off(1)
  gr_data%UVbackground_redshift_fullon  = self%myconfig%gr_UVbackground_redshift_fullon(1)
  gr_data%UVbackground_redshift_drop    = self%myconfig%gr_UVbackground_redshift_drop(1)
  gr_data%Compton_xray_heating   = self%myconfig%gr_Compton_xray_heating(1)
  gr_data%LWbackground_intensity = self%myconfig%gr_LWbackground_intensity(1)
  gr_data%LWbackground_sawtooth_suppression = self%myconfig%gr_LWbackground_sawtooth_suppression(1)
  gr_data%NumberOfTemperatureBins      = self%myconfig%gr_NumberOfTemperatureBins(1)
  gr_data%ih2co                        = self%myconfig%gr_ih2co(1)
  gr_data%ipiht                        = self%myconfig%gr_ipiht(1)
  gr_data%TemperatureStart             = self%myconfig%TemperatureStart(1)
  gr_data%TemperatureEnd               = self%myconfig%TemperatureEnd(1)
  gr_data%CaseBRecombination           = self%myconfig%gr_CaseBRecombination(1)
  gr_data%NumberOfDustTemperatureBins  = self%myconfig%gr_NumberOfDustTemperatureBins(1)
  gr_data%DustTemperatureStart         = self%myconfig%DustTemperatureStart(1)
  gr_data%DustTemperatureEnd           = self%myconfig%DustTemperatureEnd(1)
  gr_data%cloudy_electron_fraction_factor = self%myconfig%cloudy_electron_fraction_factor(1)
  ! radiative transfer parameters
  gr_data%use_radiative_transfer = self%myconfig%gr_use_radiative_transfer(1)
  gr_data%radiative_transfer_coupled_rate_solver = self%myconfig%gr_radiative_transfer_coupled_rate_solver(1)
  gr_data%radiative_transfer_intermediate_step   = self%myconfig%gr_radiative_transfer_intermediate_step(1)
  gr_data%radiative_transfer_hydrogen_only       = self%myconfig%gr_radiative_transfer_hydrogen_only(1)
  ! approximate self-shielding
  gr_data%self_shielding_method                  = self%myconfig%gr_self_shielding_method(1)
  gr_data%use_isrf_field                         = self%myconfig%gr_use_isrf_field(1)

  myunits%comoving_coordinates = self%myconfig%gr_comoving_coordinates(1)
  myunits%density_units = self%myconfig%gr_density_units(1)
  myunits%length_units = self%myconfig%gr_length_units(1)
  myunits%time_units = self%myconfig%gr_time_units(1)
  myunits%a_units = self%myconfig%gr_a_units(1)
  !     Set initial expansion factor (for internal units).
  !     Set expansion factor to 1 for non-cosmological simulation.
  myunits%a_value = self%myconfig%a_value(1)
  call set_velocity_units(myunits)


end subroutine grackle_solver_associate

subroutine grackle_chemistry_set_parameters(temgas, &
out_temperature,out_pressure,out_gamma,out_meanmw,i_obj,self)
  implicit none
  class(gr_solver),TARGET                        :: self
  real(kind=dp),intent(in)                              :: temgas
  real(kind=dp),intent(inout)                           :: out_temperature
  real(kind=dp),intent(inout)                           :: out_pressure
  real(kind=dp),intent(inout)                           :: out_gamma
  real(kind=dp),intent(inout)                           :: out_meanmw
  integer                                 :: i_obj
  ! .. local ..
  integer                                 :: iresult,i
  real(kind=dp) :: temperature_units, pressure_units
  real(kind=dp)                           :: dt_old
  integer ,parameter        :: field_size=1
  real(kind=dp), TARGET :: density(field_size), energy(field_size), &
  x_velocity(field_size), y_velocity(field_size), &
  z_velocity(field_size), &
  HI_density(field_size), HII_density(field_size), &
  HM_density(field_size), &
  HeI_density(field_size), HeII_density(field_size), &
  HeIII_density(field_size), &
  H2I_density(field_size), H2II_density(field_size), &
  DI_density(field_size), DII_density(field_size), &
  HDI_density(field_size), &
  e_density(field_size), metal_density(field_size), &
  dust_density(field_size), &
  volumetric_heating_rate(field_size), &
  specific_heating_rate(field_size), &
  RT_HI_ionization_rate(field_size), &
  RT_HeI_ionization_rate(field_size), &
  RT_HeII_ionization_rate(field_size), &
  RT_H2_dissociation_rate(field_size), &
  RT_heating_rate(field_size), &
  H2_self_shielding_length(field_size), &
  H2_custom_shielding_factor(field_size), &
  isrf_habing(field_size)
  real(kind=dp), TARGET  :: pressure(field_size),&
  gamma(field_size),temperature_gr(field_size)
  INTEGER, TARGET :: grid_rank, grid_dimension(3), grid_start(3), grid_end(3)
  real(kind=dp) :: grid_dx
  real(kind=dp)  :: nH2, nother,  gammaH2, gammaH2_1, gammaeff, gammaeff_1
  real(kind=dp)  :: gammaoth, gammaoth_1
  real(kind=dp)  :: inv_gammaH2_1, inv_gamma_eff_1, inv_gammaoth_1,tvar
  real(kind=dp)  :: meanmw
  TYPE (grackle_units) :: my_units
  TYPE (grackle_field_data) :: my_fields
  TYPE (grackle_chemistry_data) :: my_grackle_data
  TYPE (grackle_inout) :: my_config
  !-----------------------------------------------------------------------------------
  dt_old = dt

  iresult = set_default_chemistry_parameters(my_grackle_data)
  call self%link_par_to_gr(my_grackle_data,my_units)
  iresult = initialize_chemistry_data(my_units)

  grid_rank = 3
  do i = 1, grid_rank
     grid_dimension(i) = 1
     grid_start(i) = 0
     grid_end(i) = 0
  enddo
  grid_dx = 0.0d0
  grid_dimension(1) = field_size
  !0-based is necessary for grid_start and grid_end
  grid_end(1) = field_size-1
  temperature_units = get_temperature_units(my_units)
  my_units%velocity_units = get_velocity_units(my_units)

  do i = 1, field_size
    density(i) = self%myconfig%density(i_obj)/my_units%density_units
    HI_density(i) = self%myconfig%densityHI(i_obj)/my_units%density_units
    HII_density(i) = self%myconfig%densityHII(i_obj)/my_units%density_units
    HM_density(i) = self%myconfig%densityHM(i_obj)/my_units%density_units
    HeI_density(i) = self%myconfig%densityHeI(i_obj)/my_units%density_units
    HeII_density(i) = self%myconfig%densityHeII(i_obj)/my_units%density_units
    HeIII_density(i) = self%myconfig%densityHeIII(i_obj)/my_units%density_units
    H2I_density(i) = self%myconfig%densityH2I(i_obj)/my_units%density_units
    H2II_density(i) = self%myconfig%densityH2II(i_obj)/my_units%density_units
    DI_density(i) = self%myconfig%densityDI(i_obj)/my_units%density_units
    DII_density(i) = self%myconfig%densityDII(i_obj)/my_units%density_units
    HDI_density(i) = self%myconfig%densityHDI(i_obj)/my_units%density_units
    e_density(i) = self%myconfig%densityElectrons(i_obj)/my_units%density_units
    metal_density(i) = self%myconfig%density_Z(i_obj)/my_units%density_units
    dust_density(i) = self%myconfig%density_dust(i_obj)/my_units%density_units
    x_velocity(i) = 0.0
    y_velocity(i) = 0.0
    z_velocity(i) = 0.0


    ! temgas = temperature in Kelvin
    write(*,*) 'T = ', temgas, ' K'
    nH2 = 0.5d0*(self%myconfig%densityH2I(i_obj) + &
    self%myconfig%densityH2II(i_obj))
    write(*,*) 'nH2 = ',nH2
    nother = (self%myconfig%densityHeI(i_obj) + self%myconfig%densityHeII(i_obj) +&
                 self%myconfig%densityHeIII(i_obj))/4.0d0 +&
                 self%myconfig%densityHI(i_obj) + self%myconfig%densityHII(i_obj) +&
                 self%myconfig%densityElectrons(i_obj)
    write(*,*) 'nother = ',nother
    if(nH2/nother > 1.0d-3)then
      tvar = 6100.0d0/temgas
      if (tvar > 10.0)then
        inv_gammaH2_1 = 0.5d0*5.0d0
      else
        inv_gammaH2_1 = 0.5d0*(5.0d0 +&
                           2.0d0 * tvar**2.0d0 *&
                        DEXP(tvar)/(DEXP(tvar)-1.0d0)**2.0d0)
      end if
    else
      inv_gammaH2_1 = 2.5
    end if
    write(*,*) 'inv_gammaH2_1 =' , inv_gammaH2_1
    gammaoth = self%myconfig%gr_Gamma(i_obj)
    write(*,*) 'gammaoth = ',gammaoth
    inv_gammaoth_1 = 1.0d0/(gammaoth-1.0d0)
    write(*,*) 'inv_gammaoth_1 = ',inv_gammaoth_1


    gammaeff = 1.0d0 + (nH2 + nother)/&
                 (nH2*inv_gammaH2_1 +&
              nother*inv_gammaoth_1)

    gammaeff_1 = gammaeff - 1.0d0
    inv_gamma_eff_1 = 1.0d0/gammaeff_1


    write(*,*) 'gammaeff = ', gammaeff
    write(*,*) 'gammaeff - 1 =', gammaeff_1
    write(*,*) 'inv_gamma_eff_1 = ', inv_gamma_eff_1

    meanmw = self%myconfig%density(i_obj)/&
        (self%myconfig%densityHI(i_obj)+self%myconfig%densityHII(i_obj)+&
        self%myconfig%densityHM(i_obj)+&
        0.25*(self%myconfig%densityHeI(i_obj)+self%myconfig%densityHeII(i_obj)+&
        self%myconfig%densityHeIII(i_obj))+&
        0.5*(self%myconfig%densityH2I(i_obj)+self%myconfig%densityH2II(i_obj))+&
        self%myconfig%density_Z(i_obj)/16.0)

    write(*,*) 'meanmw = ',meanmw
    out_meanmw = meanmw

    energy(i)=temgas/(meanmw * temperature_units * gammaeff_1)
    write(*,*) 'internal energy', density(i)*my_units%density_units*&
    energy(i)*my_units%velocity_units * my_units%velocity_units
    volumetric_heating_rate(i) = 0.0
    specific_heating_rate(i) = 0.0
    RT_HI_ionization_rate(i) = 0.0
    RT_HeI_ionization_rate(i) = 0.0
    RT_HeII_ionization_rate(i) = 0.0
    RT_H2_dissociation_rate(i) = 0.0
    RT_heating_rate(i) = 0.0
    H2_self_shielding_length(i) = 0.0
    H2_custom_shielding_factor(i) = 0.0
    isrf_habing(i) = my_grackle_data%interstellar_radiation_field
    write(*,*) 'isrf_habing field = ', isrf_habing
  end do

  !
  !     Fill in structure to be passed to Grackle
  !

  my_fields%grid_rank = 3
  my_fields%grid_dimension = C_LOC(grid_dimension)
  my_fields%grid_start = C_LOC(grid_start)
  my_fields%grid_end = C_LOC(grid_end)
  my_fields%grid_dx  = grid_dx

  my_fields%density = C_LOC(density)
  my_fields%HI_density = C_LOC(HI_density)
  my_fields%HII_density = C_LOC(HII_density)
  my_fields%HM_density = C_LOC(HM_density)
  my_fields%HeI_density = C_LOC(HeI_density)
  my_fields%HeII_density = C_LOC(HeII_density)
  my_fields%HeIII_density = C_LOC(HeIII_density)
  my_fields%H2I_density = C_LOC(H2I_density)
  my_fields%H2II_density = C_LOC(H2II_density)
  my_fields%DI_density = C_LOC(DI_density)
  my_fields%DII_density = C_LOC(DII_density)
  my_fields%HDI_density = C_LOC(HDI_density)
  my_fields%e_density = C_LOC(e_density)
  my_fields%metal_density = C_LOC(metal_density)
  my_fields%dust_density = C_LOC(dust_density)
  my_fields%internal_energy = C_LOC(energy)
  my_fields%x_velocity = C_LOC(x_velocity)
  my_fields%y_velocity = C_LOC(y_velocity)
  my_fields%z_velocity = C_LOC(z_velocity)
  my_fields%volumetric_heating_rate = &
                                  C_LOC(volumetric_heating_rate)
  my_fields%specific_heating_rate = C_LOC(specific_heating_rate)
  my_fields%RT_HI_ionization_rate = C_LOC(RT_HI_ionization_rate)
  my_fields%RT_HeI_ionization_rate = C_LOC(RT_HeI_ionization_rate)
  my_fields%RT_HeII_ionization_rate = &
                                    C_LOC(RT_HeII_ionization_rate)
  my_fields%RT_H2_dissociation_rate = C_LOC(RT_H2_dissociation_rate)
  my_fields%RT_heating_rate = C_LOC(RT_heating_rate)
  my_fields%H2_self_shielding_length = C_LOC(H2_self_shielding_length)
  my_fields%H2_custom_shielding_factor = &
      C_LOC(H2_custom_shielding_factor)
  my_fields%isrf_habing = C_LOC(isrf_habing)
!     These routines can now be called during the simulations
!     Calculate temperature.
  iresult = calculate_temperature(my_units, my_fields, temperature_gr)
  write(*,*) "Temperature = ", temperature_gr(1), "K."
  out_temperature = temperature_gr(1)
!     Calcualte pressure.
  pressure_units = my_units%density_units *&
       my_units%velocity_units**2.0_dp
  iresult = calculate_pressure(my_units, my_fields, pressure)
  write(*,*) "Pressure = ", pressure(1)*pressure_units, "dyne/cm^2."
  out_pressure = pressure(1)*pressure_units
!     Calculate gamma.
  iresult = calculate_gamma(my_units, my_fields, gamma)
  write(*,*) "Gamma = ", gamma(1)
  out_gamma = gamma(1)
  dt = dt_old
end subroutine grackle_chemistry_set_parameters




end module mod_grackle_chemistry
