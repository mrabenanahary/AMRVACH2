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

character(len=*),parameter :: mod_grackle_chemistry_name='mod_grackle_chemistry'


! Grackle AMRVAC-defined parameters for this solver time step
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
  REAL(kind=gr_rpknd)    :: gr_current_redshift(max_num_parameters)
end type gr_config

type gr_params
  integer              :: myindice(max_num_parameters)
  real*8               :: initial_redshift
  real*8               :: temperature_units, pressure_units, dtchem
  real*8               :: grid_dx


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
  REAL(kind=gr_rpknd)    :: chi_Hydrogen(max_num_parameters) ! = chi_H

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
  real(dp)        :: density_gas(max_num_parameters)        !rhogas
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

end type gr_params

integer,allocatable              :: tstst(:)

! Fields parameters
! parts to be initialized with chemical species
! 'ism' : ISM
! 'jet' : jet (on top of the ISM)
integer :: i_zero_ism,i_zero_jet,i_zero_cloud
integer :: i_end_ism,i_end_jet,i_end_cloud

!for namelist:
CHARACTER(LEN=30), allocatable :: gr_patches_name(:)
! indices of the objects to be initialized with chemical species
INTEGER, allocatable :: gr_patches_indices_global(:)
INTEGER, allocatable :: gr_patches_indices_local(:)
CHARACTER(LEN=64), allocatable :: gr_profiles(:)
! fields config in the fraction of ISM/JET n_H
real(kind=gr_rpknd), allocatable :: gr_epsilon_tol(:) !density fraction tolerance
CHARACTER(LEN=30), allocatable :: gr_density_method(:)

type gr_objects
  logical, allocatable            :: patch(:,:^D&)           !> spatial patch
  character(len=78)               :: subname                 !> subroutine name that call it
  !> let s not define too much derived type to store
  !> gr_objects list parameters (number of objects, etc...)
  integer         :: number_of_objects
  type(gr_config) :: myconfig
  type(gr_params) :: myparams
  type(gr_fields) :: myfields
  type(grackle_type) :: mygrtype
  type(usrphysical_unit), pointer :: myphysunit
  contains
  !PRIVATE
!=BEGIN PREINITIALIZATION PROCESS SUBROUTINES=====================================================================================

!=END PREINITIALIZATION PROCESS SUBROUTINES=====================================================================================
 PROCEDURE, PASS(self) :: set_default_config => grackle_set_default
 PROCEDURE, PASS(self) :: write_setting        => grackle_write_setting
end type gr_objects

type(gr_objects), allocatable,TARGET :: gr_objects_list(:)

contains

!-------------------------------------------------------------------------
!> subroutine default setting for Grackle
subroutine grackle_set_default(self)
  implicit none
  class(gr_objects),TARGET            :: self
  !----------------------------------
  self%myconfig%obj_name(1:max_num_parameters)              = 'grackle_chemistry_config'
  self%myconfig%unit(1:max_num_parameters)                  = 'cgs'
  self%myconfig%myindice(1:max_num_parameters)              = 0

  self%myconfig%use_grackle(1:max_num_parameters) = 1
  self%myconfig%gr_with_radiative_cooling(1:max_num_parameters) = 1
  self%myconfig%gr_primordial_chemistry(1:max_num_parameters) = 0
  self%myconfig%gr_metal_cooling(1:max_num_parameters) = 1
  self%myconfig%gr_UVbackground(1:max_num_parameters) = 0
  self%myconfig%gr_cmb_temperature_floor(1:max_num_parameters) = 0
  self%myconfig%gr_h2_on_dust(1:max_num_parameters) = 0
  self%myconfig%gr_dust_chemistry(1:max_num_parameters) = 0
  self%myconfig%gr_use_dust_density_field(1:max_num_parameters) = 0
  self%myconfig%gr_photoelectric_heating(1:max_num_parameters) = 0
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
  self%myconfig%gr_Gamma(1:max_num_parameters) = 5.d0/3.d0
  self%myconfig%gr_photoelectric_heating_rate(1:max_num_parameters) = 8.5D-26

  self%myconfig%TemperatureStart(1:max_num_parameters) = 1.0d0
  self%myconfig%TemperatureEnd(1:max_num_parameters) = 1.0D9
  self%myconfig%DustTemperatureStart(1:max_num_parameters) = 1.0d0
  self%myconfig%DustTemperatureEnd(1:max_num_parameters) = 1500.0d0
  self%myconfig%gr_LWbackground_intensity(1:max_num_parameters) = 0.0d0
  self%myconfig%gr_UVbackground_redshift_on(1:max_num_parameters) = 7.0d0
  self%myconfig%gr_UVbackground_redshift_off(1:max_num_parameters) = 0.0d0
  self%myconfig%gr_UVbackground_redshift_fullon(1:max_num_parameters) = 6.0d0
  self%myconfig%gr_UVbackground_redshift_drop(1:max_num_parameters) = 0.0d0
  self%myconfig%cloudy_electron_fraction_factor(1:max_num_parameters) = 9.153959D-3
  self%myconfig%data_dir(1:max_num_parameters) = "../../../src/grackle/input/"
  self%myconfig%data_filename(1:max_num_parameters) = "CloudyData_UVB=HM2012.h5"
  self%myparams%data_file(1:max_num_parameters) = self%myconfig%data_dir(1:max_num_parameters)//&
  self%myconfig%data_filename(1:max_num_parameters)//C_NULL_CHAR
  self%myconfig%normalize_done(1:max_num_parameters) = .false.
  self%myconfig%gr_comoving_coordinates(1:max_num_parameters) = 0
  self%myconfig%gr_a_units(1:max_num_parameters) = 0.0d0
  self%myconfig%gr_current_redshift(1:max_num_parameters) = 0.

  ! Parameters default values
  ! Gas to dust ratios
  self%myparams%myindice(1:max_num_parameters) = 0
  self%myparams%chi_dust(1:max_num_parameters) = -1.0d0 ! 0.0d0
  self%myparams%xi_dust(1:max_num_parameters) = -1.0d0 ! 0.0d0
  ! Z = Z _solar :
  self%myparams%HydrogenFractionByMass(1:max_num_parameters) = -1.0d0 ! 0.76d0
  self%myparams%HeliumFractionByMass(1:max_num_parameters) = -1.0d0 ! 0.22705d0
  self%myparams%SolarMetalFractionByMass(1:max_num_parameters) = -1.0d0 ! 0.01295d0 !Solar metallicity
  self%myparams%MetalFractionByMass(1:max_num_parameters) = self%myparams%SolarMetalFractionByMass(1:max_num_parameters)
  self%myparams%chi_Hydrogen(1:max_num_parameters) = -1.0d0 ! 0.76d0 !<- still to compute relative to HydrogenFractionByMass
  self%myparams%DeuteriumToHydrogenRatio(1:max_num_parameters) = -1.0d0 ! 2.0d0*3.4d-5
  self%myparams%IonizationFraction(1:max_num_parameters) = -1.0d0 !  0.0d0
  self%myparams%He_abundance(1:max_num_parameters) = -1.0d0 ! 0.0d0
  self%myparams%Zeta_nH_to_rho(1:max_num_parameters) = -1.0d0 ! 0.0d0
  self%myparams%ZetaPrime_nH_to_rho(1:max_num_parameters) = -1.0d0

  !Species abundances
  self%myparams%number_of_solved_species(1:max_num_parameters) = 0
  self%myparams%x_HI(1:max_num_parameters) = 0.0d0
  self%myparams%x_HII(1:max_num_parameters) = 0.0d0
  self%myparams%x_HM(1:max_num_parameters) = 0.0d0
  self%myparams%x_HeI(1:max_num_parameters) = 0.0d0
  self%myparams%x_HeII(1:max_num_parameters) = 0.0d0
  self%myparams%x_HeIII(1:max_num_parameters) = 0.0d0
  self%myparams%x_H2I(1:max_num_parameters) = 0.0d0
  self%myparams%x_H2II(1:max_num_parameters) = 0.0d0
  self%myparams%x_DI(1:max_num_parameters) = 0.0d0
  self%myparams%x_DII(1:max_num_parameters) = 0.0d0
  self%myparams%x_HDI(1:max_num_parameters) = 0.0d0
  self%myparams%x_e(1:max_num_parameters) = 0.0d0
  self%myparams%densityDI(1:max_num_parameters) = 0.0d0
  self%myparams%densityDII(1:max_num_parameters) = 0.0d0
  self%myparams%densityHDI(1:max_num_parameters) = 0.0d0
  self%myparams%densityHI(1:max_num_parameters) = 0.0d0
  self%myparams%densityHII(1:max_num_parameters) = 0.0d0
  self%myparams%densityHM(1:max_num_parameters) = 0.0d0
  self%myparams%densityH2I(1:max_num_parameters) = 0.0d0
  self%myparams%densityH2II(1:max_num_parameters) = 0.0d0
  self%myparams%densityHeI(1:max_num_parameters) = 0.0d0
  self%myparams%densityHeII(1:max_num_parameters) = 0.0d0
  self%myparams%densityHeIII(1:max_num_parameters) = 0.0d0

  self%myparams%density(1:max_num_parameters) = -1.0d0 ! 0.0d0
  self%myparams%number_H_density(1:max_num_parameters) = -1.0d0 ! 0.0d0
  self%myparams%density_dust(1:max_num_parameters) = -1.0d0 ! 0.0d0
  self%myparams%density_gas(1:max_num_parameters) = -1.0d0 ! 0.0d0
  self%myparams%density_bar(1:max_num_parameters) = -1.0d0 ! 0.0d0
  self%myparams%density_X(1:max_num_parameters) = -1.0d0 ! 0.0d0
  self%myparams%density_Y(1:max_num_parameters) = -1.0d0 ! 0.0d0
  self%myparams%density_Z(1:max_num_parameters) = -1.0d0 ! 0.0d0
  self%myparams%density_deut(1:max_num_parameters) = -1.0d0 ! 0.0d0
  self%myparams%density_not_deut(1:max_num_parameters) = -1.0d0 ! 0.0d0
  self%myparams%density_metal_free(1:max_num_parameters) = -1.0d0 ! 0.0d0
  self%myparams%densityD(1:max_num_parameters) = -1.0d0 ! 0.0d0
  self%myparams%densityHD(1:max_num_parameters) = -1.0d0 ! 0.0d0
  self%myparams%densityH(1:max_num_parameters) = -1.0d0 ! 0.0d0
  self%myparams%densityHtwo(1:max_num_parameters) = -1.0d0 ! 0.0d0
  self%myparams%densityHe(1:max_num_parameters) = -1.0d0 ! 0.0d0
  self%myparams%densityElectrons(1:max_num_parameters) = -1.0d0 ! 0.0d0

  write(*,*) 'Grackle configuration defaulting successfully done !'



end subroutine grackle_set_default


!-------------------------------------------------------------------------
subroutine grackle_fields_config_allocate(number_of_isms,number_of_jets,number_of_clouds)
 implicit none

 integer,intent(inout),optional        :: number_of_isms
 integer,intent(inout),optional        :: number_of_jets
 integer,intent(inout),optional        :: number_of_clouds
 integer        :: default_number_of_isms
 integer        :: default_number_of_jets
 integer        :: default_number_of_clouds
 ! .. local ..
 integer                  :: n_total_objects
 integer                  :: i_ism,i_jet,i_cloud,i_all,counter
 !-------------------------------------------------------------------------

 ! Defaults
 default_number_of_isms = 1
 default_number_of_jets = 0
 default_number_of_clouds = 0
 if(.not.present(number_of_isms))number_of_isms = default_number_of_isms
 if(.not.present(number_of_jets))number_of_jets = default_number_of_jets
 if(.not.present(number_of_clouds))number_of_clouds = default_number_of_clouds


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

 ! Second, deallocate already allocated arrays


 if(allocated(gr_objects_list))deallocate(gr_objects_list)
 if(allocated(gr_patches_indices_global))deallocate(gr_patches_indices_global)
 if(allocated(gr_patches_indices_local))deallocate(gr_patches_indices_local)
 if(allocated(gr_patches_name))deallocate(gr_patches_name)
 if(allocated(gr_profiles))deallocate(gr_profiles)
 if(allocated(gr_epsilon_tol))deallocate(gr_epsilon_tol)
 if(allocated(gr_density_method))deallocate(gr_density_method)


 ! Third, allocate the size of the allocatable arrays

 if(.not.allocated(gr_objects_list))allocate(gr_objects_list(1:n_total_objects))
 if(.not.allocated(gr_patches_indices_global))allocate(gr_patches_indices_global(1:n_total_objects))
 if(.not.allocated(gr_patches_indices_local))allocate(gr_patches_indices_local(1:n_total_objects))
 if(.not.allocated(gr_patches_name))allocate(gr_patches_name(1:n_total_objects))
 if(.not.allocated(gr_profiles))allocate(gr_profiles(1:n_total_objects))
 if(.not.allocated(gr_epsilon_tol))allocate(gr_epsilon_tol(1:n_total_objects))
 if(.not.allocated(gr_density_method))allocate(gr_density_method(1:n_total_objects))

 ! Store the total number of objects
 do i_all=1,size(gr_patches_indices_global)
  gr_objects_list(i_all)%number_of_objects = n_total_objects
 end do

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

subroutine grackle_fields_config_set_default(i_all)
 implicit none
 integer,intent(in)                            :: i_all
 !-------------------------------------------------------------------------

  ! Fill in default values for everything, independently

  ! Fields parameters


  gr_profiles(gr_patches_indices_global(i_all)) = 'uniform'
  ! fields config
  gr_epsilon_tol(gr_patches_indices_global(i_all)) = 14.01*tiny_number !density fraction tolerance

  ! Densities
  gr_density_method(gr_patches_indices_global(i_all)) = 'chemical_coefficient' !method to estimate following densities

end subroutine grackle_fields_config_set_default

!-------------------------------------------------------------------------
!> Read grackle parameters  from a parfile
subroutine grackle_config_read(grackle_config,files)
  implicit none

  character(len=*),intent(in)        :: files(:)
  type(gr_config), intent(out)  :: grackle_config

  ! .. local ..
  integer                            :: i_file,i_error_read
  integer                  :: idim,iside
  character(len=70)            :: error_message
  !-------------------------------------------------------------------------
  namelist /grackle_conf_list/  grackle_config,tstst
  namelist /grackle_conf1_list/ grackle_config,tstst
  namelist /grackle_conf2_list/ grackle_config,tstst
  namelist /grackle_conf3_list/ grackle_config,tstst


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
       write(*,*) 'use_grackle : ', grackle_config%use_grackle(1)
       write(*,*) 'gr_primordial_chemistry : ', grackle_config%gr_primordial_chemistry(1)
       write(*,*) 'gr_metal_cooling : ', grackle_config%gr_metal_cooling(1)
       write(*,*) 'gr_dust_chemistry : ', grackle_config%gr_dust_chemistry(1)
     end select
     call usr_mat_read_error_message(i_error_read,grackle_config%myindice(1),&
                                     'grackle_chemistry_config')
     close(unitpar)

     !add a routine to check that the sum of each density is equal to one times
     ! the total density in ISM AND JET
  end do

  write(*,*) 'tstst = ', tstst
  write(*,*) 'Grackle usr configuration reading successfully done !'
end subroutine grackle_config_read

subroutine grackle_params_read(grackle_par,files)
  implicit none

  character(len=*),intent(in)        :: files(:)
  type(gr_params), intent(out)  :: grackle_par

  ! .. local ..
  integer                            :: i_file,i_error_read
  integer                  :: idim,iside
  character(len=70)            :: error_message
  !-------------------------------------------------------------------------
  namelist /grackle_par_list/  grackle_par
  namelist /grackle_par1_list/ grackle_par
  namelist /grackle_par2_list/ grackle_par
  namelist /grackle_par3_list/ grackle_par


  error_message = 'In the procedure : grackle_params_read'

  if(mype==0)write(*,*)'Reading grackle_par_list'
  do i_file = 1, size(files)
     open(unitpar, file=trim(files(i_file)), status="old")
     select case(grackle_par%myindice(1))
     case(1)
       read(unitpar, grackle_par1_list, iostat=i_error_read)
     case(2)
       read(unitpar, grackle_par2_list, iostat=i_error_read)
     case(3)
       read(unitpar, grackle_par3_list, iostat=i_error_read)
     case default
       read(unitpar, grackle_par_list, iostat=i_error_read)
     end select
     call usr_mat_read_error_message(i_error_read,grackle_par%myindice(1),&
                                     'grackle_chemistry_par')
     close(unitpar)

     !add a routine to check that the sum of each density is equal to one times
     ! the total density in ISM AND JET
  end do
  write(*,*) 'Grackle usr parameters reading successfully done !'
end subroutine grackle_params_read

subroutine grackle_associate_parameters(n_obj,gr_conf_lst,gr_par_lst,gr_obj_lst,i_all)
  implicit none
  integer,intent(in)                            :: n_obj
  type(gr_config),intent(in)                    :: gr_conf_lst
  type(gr_params),intent(in)                    :: gr_par_lst
  type(gr_objects), intent(inout),TARGET        :: gr_obj_lst(1:n_obj)
  integer,intent(in)                            :: i_all
  !--------------------------------------------------------

  gr_obj_lst(i_all)%myconfig%obj_name(1)              = gr_conf_lst%obj_name(i_all)
  gr_obj_lst(i_all)%myconfig%unit(1)                  = gr_conf_lst%unit(i_all)
  !gr_obj_lst(i_all)%myconfig%myindice(1)              = gr_conf_lst%myindice(i_all)

  gr_obj_lst(i_all)%myconfig%use_grackle(1) = gr_conf_lst%use_grackle(i_all)
  gr_obj_lst(i_all)%myconfig%gr_with_radiative_cooling(1) = gr_conf_lst%gr_with_radiative_cooling(i_all)
  gr_obj_lst(i_all)%myconfig%gr_primordial_chemistry(1) = gr_conf_lst%gr_primordial_chemistry(i_all)
  gr_obj_lst(i_all)%myconfig%gr_metal_cooling(1) = gr_conf_lst%gr_metal_cooling(i_all)
  gr_obj_lst(i_all)%myconfig%gr_UVbackground(1) = gr_conf_lst%gr_UVbackground(i_all)
  gr_obj_lst(i_all)%myconfig%gr_cmb_temperature_floor(1) = gr_conf_lst%gr_cmb_temperature_floor(i_all)
  gr_obj_lst(i_all)%myconfig%gr_h2_on_dust(1) = gr_conf_lst%gr_h2_on_dust(i_all)
  gr_obj_lst(i_all)%myconfig%gr_dust_chemistry(1) = gr_conf_lst%gr_dust_chemistry(i_all)
  gr_obj_lst(i_all)%myconfig%gr_use_dust_density_field(1) = gr_conf_lst%gr_use_dust_density_field(i_all)
  gr_obj_lst(i_all)%myconfig%gr_photoelectric_heating(1) = gr_conf_lst%gr_photoelectric_heating(i_all)
  gr_obj_lst(i_all)%myconfig%gr_use_volumetric_heating_rate(1) = gr_conf_lst%gr_use_volumetric_heating_rate(i_all)
  gr_obj_lst(i_all)%myconfig%gr_use_specific_heating_rate(1) = gr_conf_lst%gr_use_specific_heating_rate(i_all)
  gr_obj_lst(i_all)%myconfig%gr_three_body_rate(1) = gr_conf_lst%gr_three_body_rate(i_all)
  gr_obj_lst(i_all)%myconfig%gr_cie_cooling(1) = gr_conf_lst%gr_cie_cooling(i_all)
  gr_obj_lst(i_all)%myconfig%gr_h2_optical_depth_approximation(1) = gr_conf_lst%gr_h2_optical_depth_approximation(i_all)
  gr_obj_lst(i_all)%myconfig%gr_ih2co(1) = gr_conf_lst%gr_ih2co(i_all)
  gr_obj_lst(i_all)%myconfig%gr_ipiht(1) = gr_conf_lst%gr_ipiht(i_all)
  gr_obj_lst(i_all)%myconfig%gr_NumberOfTemperatureBins(1) = gr_conf_lst%gr_NumberOfTemperatureBins(i_all)
  gr_obj_lst(i_all)%myconfig%gr_CaseBRecombination(1) = gr_conf_lst%gr_CaseBRecombination(i_all)
  gr_obj_lst(i_all)%myconfig%gr_Compton_xray_heating(1)= gr_conf_lst%gr_Compton_xray_heating(i_all)
  gr_obj_lst(i_all)%myconfig%gr_LWbackground_sawtooth_suppression(1) = gr_conf_lst%gr_LWbackground_sawtooth_suppression(i_all)
  gr_obj_lst(i_all)%myconfig%gr_NumberOfDustTemperatureBins(1) = gr_conf_lst%gr_NumberOfDustTemperatureBins(i_all)
  gr_obj_lst(i_all)%myconfig%gr_use_radiative_transfer(1) = gr_conf_lst%gr_use_radiative_transfer(i_all)
  gr_obj_lst(i_all)%myconfig%gr_radiative_transfer_coupled_rate_solver(1) = gr_conf_lst%gr_radiative_transfer_coupled_rate_solver(i_all)
  gr_obj_lst(i_all)%myconfig%gr_radiative_transfer_intermediate_step(1) = gr_conf_lst%gr_radiative_transfer_intermediate_step(i_all)
  gr_obj_lst(i_all)%myconfig%gr_radiative_transfer_hydrogen_only(1) = gr_conf_lst%gr_radiative_transfer_hydrogen_only(i_all)
  gr_obj_lst(i_all)%myconfig%gr_self_shielding_method(1) = gr_conf_lst%gr_self_shielding_method(i_all)
  gr_obj_lst(i_all)%myconfig%gr_Gamma(1) = gr_conf_lst%gr_Gamma(i_all)
  gr_obj_lst(i_all)%myconfig%gr_photoelectric_heating_rate(1) = gr_conf_lst%gr_photoelectric_heating_rate(i_all)

  gr_obj_lst(i_all)%myconfig%TemperatureStart(1) = gr_conf_lst%TemperatureStart(i_all)
  gr_obj_lst(i_all)%myconfig%TemperatureEnd(1) = gr_conf_lst%TemperatureEnd(i_all)
  gr_obj_lst(i_all)%myconfig%DustTemperatureStart(1) = gr_conf_lst%DustTemperatureStart(i_all)
  gr_obj_lst(i_all)%myconfig%DustTemperatureEnd(1) = gr_conf_lst%DustTemperatureEnd(i_all)
  gr_obj_lst(i_all)%myconfig%gr_LWbackground_intensity(1) = gr_conf_lst%gr_LWbackground_intensity(i_all)
  gr_obj_lst(i_all)%myconfig%gr_UVbackground_redshift_on(1) = gr_conf_lst%gr_UVbackground_redshift_on(i_all)
  gr_obj_lst(i_all)%myconfig%gr_UVbackground_redshift_off(1) = gr_conf_lst%gr_UVbackground_redshift_off(i_all)
  gr_obj_lst(i_all)%myconfig%gr_UVbackground_redshift_fullon(1) = gr_conf_lst%gr_UVbackground_redshift_fullon(i_all)
  gr_obj_lst(i_all)%myconfig%gr_UVbackground_redshift_drop(1) = gr_conf_lst%gr_UVbackground_redshift_drop(i_all)
  gr_obj_lst(i_all)%myconfig%cloudy_electron_fraction_factor(1) = gr_conf_lst%cloudy_electron_fraction_factor(i_all)
  gr_obj_lst(i_all)%myconfig%data_dir(1) = gr_conf_lst%data_dir(i_all)
  gr_obj_lst(i_all)%myconfig%data_filename(1) = gr_conf_lst%data_filename(i_all)
  gr_obj_lst(i_all)%myparams%data_file(1) = gr_obj_lst(i_all)%myconfig%data_dir(1)//&
  gr_obj_lst(i_all)%myconfig%data_filename(1)//C_NULL_CHAR
  gr_obj_lst(i_all)%myconfig%normalize_done(1) = gr_conf_lst%normalize_done(i_all)
  gr_obj_lst(i_all)%myconfig%gr_comoving_coordinates(1) = gr_conf_lst%gr_comoving_coordinates(i_all)
  gr_obj_lst(i_all)%myconfig%gr_a_units(1) = gr_conf_lst%gr_a_units(i_all)
  gr_obj_lst(i_all)%myconfig%gr_current_redshift(1) = gr_conf_lst%gr_current_redshift(i_all)

  ! Parameters default values
  ! Gas to dust ratios
  gr_obj_lst(i_all)%myparams%myindice(1) = gr_par_lst%myindice(i_all)
  gr_obj_lst(i_all)%myparams%chi_dust(1) = gr_par_lst%chi_dust(i_all)
  gr_obj_lst(i_all)%myparams%xi_dust(1) = gr_par_lst%xi_dust(i_all)
  gr_obj_lst(i_all)%myparams%HydrogenFractionByMass(1) = gr_par_lst%HydrogenFractionByMass(i_all)
  gr_obj_lst(i_all)%myparams%HeliumFractionByMass(1) = gr_par_lst%HeliumFractionByMass(i_all)
  gr_obj_lst(i_all)%myparams%MetalFractionByMass(1) = gr_par_lst%MetalFractionByMass(i_all)
  gr_obj_lst(i_all)%myparams%SolarMetalFractionByMass(1) = gr_par_lst%SolarMetalFractionByMass(i_all)
  gr_obj_lst(i_all)%myparams%chi_Hydrogen(1) = gr_par_lst%chi_Hydrogen(i_all)
  gr_obj_lst(i_all)%myparams%DeuteriumToHydrogenRatio(1) = gr_par_lst%DeuteriumToHydrogenRatio(i_all)
  gr_obj_lst(i_all)%myparams%IonizationFraction(1) = gr_par_lst%IonizationFraction(i_all)
  gr_obj_lst(i_all)%myparams%He_abundance(1) = gr_par_lst%He_abundance(i_all)
  gr_obj_lst(i_all)%myparams%Zeta_nH_to_rho(1) = gr_par_lst%Zeta_nH_to_rho(i_all)
  gr_obj_lst(i_all)%myparams%ZetaPrime_nH_to_rho(1) = gr_par_lst%ZetaPrime_nH_to_rho(i_all)



  !Species abundances
  gr_obj_lst(i_all)%myparams%number_of_solved_species(1) = &
   convert_logical_to_integer(gr_obj_lst(i_all)%myconfig%gr_primordial_chemistry(1)>0)*6&
  +convert_logical_to_integer(gr_obj_lst(i_all)%myconfig%gr_primordial_chemistry(1)>1)*3&
  +convert_logical_to_integer(gr_obj_lst(i_all)%myconfig%gr_primordial_chemistry(1)>2)*3
  if(gr_obj_lst(i_all)%myconfig%gr_primordial_chemistry(1)>0)then
    gr_obj_lst(i_all)%myparams%x_HI(1) = gr_par_lst%x_HI(i_all)
    gr_obj_lst(i_all)%myparams%x_HII(1) = gr_par_lst%x_HII(i_all)
    gr_obj_lst(i_all)%myparams%x_HeI(1) = gr_par_lst%x_HeI(i_all)
    gr_obj_lst(i_all)%myparams%x_HeII(1) = gr_par_lst%x_HeII(i_all)
    gr_obj_lst(i_all)%myparams%x_HeIII(1) = gr_par_lst%x_HeIII(i_all)
    gr_obj_lst(i_all)%myparams%x_e(1) = gr_par_lst%x_e(i_all)
    if(gr_obj_lst(i_all)%myconfig%gr_primordial_chemistry(1)>1)then
      gr_obj_lst(i_all)%myparams%x_HM(1) = gr_par_lst%x_HM(i_all)
      gr_obj_lst(i_all)%myparams%x_H2I(1) = gr_par_lst%x_H2I(i_all)
      gr_obj_lst(i_all)%myparams%x_H2II(1) = gr_par_lst%x_H2II(i_all)
      if(gr_obj_lst(i_all)%myconfig%gr_primordial_chemistry(1)>2)then
        gr_obj_lst(i_all)%myparams%x_DI(1) = gr_par_lst%x_DI(i_all)
        gr_obj_lst(i_all)%myparams%x_DII(1) = gr_par_lst%x_DII(i_all)
        gr_obj_lst(i_all)%myparams%x_HDI(1) = gr_par_lst%x_HDI(i_all)
      end if
    end if
  end if

  gr_obj_lst(i_all)%myparams%density(1) = gr_par_lst%density(i_all)
  gr_obj_lst(i_all)%myparams%number_H_density(1) = gr_par_lst%number_H_density(i_all)
  gr_obj_lst(i_all)%myparams%density_dust(1) = gr_par_lst%density_dust(i_all)
  gr_obj_lst(i_all)%myparams%density_gas(1) = gr_par_lst%density_gas(i_all)
  gr_obj_lst(i_all)%myparams%density_bar(1) = gr_par_lst%density_bar(i_all)
  gr_obj_lst(i_all)%myparams%density_X(1) = gr_par_lst%density_X(i_all)
  gr_obj_lst(i_all)%myparams%density_Y(1) = gr_par_lst%density_Y(i_all)
  gr_obj_lst(i_all)%myparams%density_Z(1) = gr_par_lst%density_Z(i_all)
  gr_obj_lst(i_all)%myparams%density_deut(1) = gr_par_lst%density_deut(i_all)
  gr_obj_lst(i_all)%myparams%density_not_deut(1) = gr_par_lst%density_not_deut(i_all)
  gr_obj_lst(i_all)%myparams%density_metal_free(1) = gr_par_lst%density_metal_free(i_all)
  gr_obj_lst(i_all)%myparams%densityD(1) = gr_par_lst%densityD(i_all)
  gr_obj_lst(i_all)%myparams%densityHD(1) = gr_par_lst%densityHD(i_all)
  gr_obj_lst(i_all)%myparams%densityH(1) = gr_par_lst%densityH(i_all)
  gr_obj_lst(i_all)%myparams%densityHtwo(1) = gr_par_lst%densityHtwo(i_all)
  gr_obj_lst(i_all)%myparams%densityHe(1) = gr_par_lst%densityHe(i_all)
  gr_obj_lst(i_all)%myparams%densityElectrons(1) = gr_par_lst%densityElectrons(i_all)

end subroutine grackle_associate_parameters


subroutine grackle_write_setting(self,unit_config,iobj)
  implicit none
  class(gr_objects),TARGET            :: self
  integer,intent(in)                  :: unit_config
  integer,intent(in)                  :: iobj
  integer                             :: idims2,iside2,iB2
  real(kind=dp)                       :: rto_print
  character(len=64)                   :: sto_print
  character(len=128)                   :: wto_print
  integer                             :: idim,iside,idims,iw2,i_all
  ! .. local ..

  !-----------------------------------

  if(self%myconfig%use_grackle(1)==1)then
    write(unit_config,*)'      ****** Parameters for ',gr_patches_name(gr_patches_indices_global(iobj)),&
    '#',gr_patches_indices_local(gr_patches_indices_global(iobj)),'******      '
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
    write(unit_config,*) 'normalize_done = ',  self%myconfig%normalize_done(1)
    write(unit_config,*) 'gr_comoving_coordinates = ',  self%myconfig%gr_comoving_coordinates(1)
    write(unit_config,*) 'gr_a_units = ',  self%myconfig%gr_a_units(1)
    write(unit_config,*) 'gr_current_redshift = ',  self%myconfig%gr_current_redshift(1)

    write(unit_config,*)'> Physical parameters'
    write(unit_config,*) 'chi_dust = ',  self%myparams%chi_dust(1)
    write(unit_config,*) 'xi_dust = ',  self%myparams%xi_dust(1)
    write(unit_config,*) 'HydrogenFractionByMass = ',  self%myparams%HydrogenFractionByMass(1)
    write(unit_config,*) 'HeliumFractionByMass = ',  self%myparams%HeliumFractionByMass(1)
    write(unit_config,*) 'MetalFractionByMass = ',  self%myparams%MetalFractionByMass(1)

    write(unit_config,*) 'SolarMetalFractionByMass = ',  self%myparams%SolarMetalFractionByMass(1)
    write(unit_config,*) 'chi_Hydrogen = ',  self%myparams%chi_Hydrogen(1)
    write(unit_config,*) 'DeuteriumToHydrogenRatio = ',  self%myparams%DeuteriumToHydrogenRatio(1)
    write(unit_config,*) 'IonizationFraction = ',  self%myparams%IonizationFraction(1)
    write(unit_config,*) 'He_abundance = ',  self%myparams%He_abundance(1)
    write(unit_config,*) 'Zeta_nH_to_rho = ',  self%myparams%Zeta_nH_to_rho(1)
    write(unit_config,*) 'ZetaPrime_nH_to_rho = ',  self%myparams%ZetaPrime_nH_to_rho(1)

    write(unit_config,*)'>> ',self%myparams%number_of_solved_species(1),' solved species abundances : '
    !Species abundances
    if(self%myconfig%gr_primordial_chemistry(1)>0)then
      write(unit_config,*)' + [HI/H] : ' , self%myparams%x_HI(1)
      write(unit_config,*)' + [HII/H] : ' , self%myparams%x_HII(1)
      write(unit_config,*)' + [HeI/H] : ' , self%myparams%x_HeI(1)
      write(unit_config,*)' + [HeII/H] : ' , self%myparams%x_HeII(1)
      write(unit_config,*)' + [HeIII/H] : ' , self%myparams%x_HeIII(1)
      write(unit_config,*)' + [e-/H] : ' , self%myparams%x_e(1)
      if(self%myconfig%gr_primordial_chemistry(1)>1)then
        write(unit_config,*)' + [HM/H] : ' , self%myparams%x_HM(1)
        write(unit_config,*)' + [H2I/H] : ' , self%myparams%x_H2I(1)
        write(unit_config,*)' + [H2II/H] : ' , self%myparams%x_H2II(1)
        if(self%myconfig%gr_primordial_chemistry(1)>2)then
          write(unit_config,*)' + [DI/H] : ' , self%myparams%x_DI(1)
          write(unit_config,*)' + [DII/H] : ' , self%myparams%x_DII(1)
          write(unit_config,*)' + [HDI/H] : ' , self%myparams%x_HDI(1)
        end if
      end if
    end if

    write(unit_config,*) 'density = ',  self%myparams%density(1)
    write(unit_config,*) 'number_H_density = ',  self%myparams%number_H_density(1)
    write(unit_config,*) 'density_dust = ',  self%myparams%density_dust(1)
    write(unit_config,*) 'density_gas = ',  self%myparams%density_gas(1)
    write(unit_config,*) 'density_bar = ',  self%myparams%density_bar(1)
    write(unit_config,*) 'density_X = ',  self%myparams%density_X(1)
    write(unit_config,*) 'density_Y = ',  self%myparams%density_Y(1)
    write(unit_config,*) 'density_Z = ',  self%myparams%density_Z(1)
    write(unit_config,*) 'density_deut = ',  self%myparams%density_deut(1)
    write(unit_config,*) 'density_not_deut = ',  self%myparams%density_not_deut(1)
    write(unit_config,*) 'density_metal_free = ',  self%myparams%density_metal_free(1)
    write(unit_config,*) 'densityD = ',  self%myparams%densityD(1)
    write(unit_config,*) 'densityHD = ',  self%myparams%densityHD(1)
    write(unit_config,*) 'densityH = ',  self%myparams%densityH(1)
    write(unit_config,*) 'densityHtwo = ',  self%myparams%densityHtwo(1)
    write(unit_config,*) 'densityHe = ',  self%myparams%densityHe(1)
    write(unit_config,*) 'densityElectrons = ',  self%myparams%densityElectrons(1)
    ! fields config
    write(unit_config,*)'> Objects parameters'
    write(unit_config,*) 'gr_patches_indices_global = ',gr_patches_indices_global(gr_patches_indices_global(iobj))
    write(unit_config,*) 'gr_epsilon_tol = ',gr_epsilon_tol(gr_patches_indices_global(iobj))
    write(unit_config,*) 'gr_profiles = ',gr_profiles(gr_patches_indices_global(iobj))
    write(unit_config,*) 'gr_density_method = ',gr_density_method(gr_patches_indices_global(iobj))

    write(unit_config,*) '======================================================='
    write(unit_config,*)'************************************'
    write(unit_config,*)'******** END of this object s Grackle setting **********'
    write(unit_config,*)'************************************'
  else
    write(unit_config,*)'      ****** Grackle disabled for ',gr_patches_name(gr_patches_indices_global(iobj)),&
    '(',gr_patches_indices_local(gr_patches_indices_global(iobj)),') ******      '
  end if


end    subroutine grackle_write_setting


!--------------------------------------------------------------------
!> subroutine check the parfile setting for ism
subroutine grackle_set_complet(self,ref_density,normalized,mydensityunit)
  implicit none
  class(gr_objects), TARGET                                :: self
  real(kind=dp),intent(in) :: ref_density
  logical,intent(in)       :: normalized
  real(kind=dp),optional :: mydensityunit
  ! .. local ..
  integer :: i_all,icorrect
  real(kind=dp) :: physical_ref_density
  logical :: HNotcorrected,HeNotcorrected
  !-----------------------------------
  if(normalized.and.(.not.present(mydensityunit)))then
    call mpistop('normalized==true but no density unit given')
  end if
  if(normalized)then
    physical_ref_density = ref_density*mydensityunit
  else
    physical_ref_density = ref_density
  end if
  HNotcorrected=.true.
  HeNotcorrected=.true.


  ! TO DO: negatify and ajust the parameters according to which components
  ! (composition,coolings,rates,etc...) is not enabled
  ! Default parameters values :
  ! * Z_solar = 0.01295
  ! * Z = Z_solar
  ! * chi_dust or xi_dust must be given
  ! * x_ion = 1e-3
  ! * chi_D = 2.0*3.4e-5 /* The DToHRatio is by mass in the code, so multiply by 2. */
  ! * x(He) = 10 %
  ! ==> fix Zeta and n_H

  Loop_through_all : do i_all = 1,self%number_of_objects


  do icorrect = 1, 2
  ! Z_solar
  if(self%myparams%SolarMetalFractionByMass(1)<-smalldouble)then
    self%myparams%SolarMetalFractionByMass(1) = 0.01295d0
  else
    self%myparams%SolarMetalFractionByMass(1)=DABS(self%myparams%SolarMetalFractionByMass(1))
  end if

  ! Z
  if(self%myparams%MetalFractionByMass(1)<-smalldouble)then
    self%myparams%MetalFractionByMass(1)=self%myparams%SolarMetalFractionByMass(1)
  end if


  ! chi_dust and xi_dust
  if(self%myparams%chi_dust(1)<-smalldouble.and.self%myparams%xi_dust(1)>=-smalldouble)then
      self%myparams%chi_dust(1)=DABS(self%myparams%xi_dust(1))
  elseif(self%myparams%chi_dust(1)>=-smalldouble.and.self%myparams%xi_dust(1)<-smalldouble)then
      self%myparams%xi_dust(1)=DABS(self%myparams%chi_dust(1))
  else
    call return_error(i_all,'chi_dust','xi_dust',&
                      self%myparams%chi_dust(1),&
                      self%myparams%xi_dust(1))
  end if

  !fix rhogas
  self%myparams%density_gas(1) = physical_ref_density/(1.0_dp+self%myparams%chi_dust(1))

  !fix rho_dust:
  self%myparams%density_dust(1) = physical_ref_density*&
  self%myparams%chi_dust(1)/(1.0_dp+self%myparams%chi_dust(1))

  ! x_ion
  if(self%myparams%IonizationFraction(1)<-smalldouble)then
     self%myparams%IonizationFraction(1) = 1.0d-3
  end if

  ! chi_D
  if(self%myparams%DeuteriumToHydrogenRatio(1)<-smalldouble)then
    self%myparams%DeuteriumToHydrogenRatio(1) = 2.0d0*3.4d-5
  end if

  ! x(He)
  if(self%myparams%He_abundance(1)<-smalldouble)then
    self%myparams%He_abundance(1)= 0.1_dp
  end if

  !Fix reference Zeta:
  self%myparams%Zeta_nH_to_rho(1) = &
  ((1.0_dp+self%myparams%DeuteriumToHydrogenRatio(1))/(1.0_dp-self%myparams%MetalFractionByMass(1)))*&
  (w_HI*self%myparams%x_HI(1) + w_HII*self%myparams%x_HII(1) + w_HM*self%myparams%x_HM(1)+&
  w_H2I*self%myparams%x_H2I(1) + w_H2II*self%myparams%x_H2II(1))+&
  (self%myparams%He_abundance(1)*w_HeI/(1.0_dp-self%myparams%MetalFractionByMass(1)))+&
  w_e*self%myparams%IonizationFraction(1)

  !Fix reference Zeta Prime:
  self%myparams%ZetaPrime_nH_to_rho(1) = &
  self%myparams%Zeta_nH_to_rho(1) -&
  self%myparams%He_abundance(1)*w_HeI/(1.0_dp -&
  self%myparams%MetalFractionByMass(1))

  !Fix reference n_H:
  self%myparams%number_H_density(1) = (physical_ref_density - self%myparams%density_dust(1))/&
  (self%myparams%Zeta_nH_to_rho(1)*mh_gr)

  !TODO: complete and check the species abundances

  !rho_Electrons
  self%myparams%densityElectrons(1) = w_e*self%myparams%IonizationFraction(1)*&
  self%myparams%number_H_density(1)*mh_gr

  !fix rhobar = rho - rhodust - rho_electrons
  self%myparams%density_bar(1) = self%myparams%density_gas(1) - &
  self%myparams%densityElectrons(1)

  !rhoZ
  self%myparams%density_Z(1) = self%myparams%MetalFractionByMass(1)*&
  self%myparams%density_bar(1)

  !Y
  self%myparams%HeliumFractionByMass(1) = &
  (self%myparams%He_abundance(1)/self%myparams%Zeta_nH_to_rho(1))*&
  w_HeI/(1.0_dp-(self%myparams%densityElectrons(1)/(physical_ref_density-&
  self%myparams%density_dust(1))))

  if(DABS(self%myparams%He_abundance(1)-&
  (self%myparams%x_HeI(1)+self%myparams%x_HeII(1)+&
  self%myparams%x_HeIII(1)))/DABS(self%myparams%He_abundance(1))>abundance_tolerance)then
    WRITE(*,*) 'Error: we have x(He) /= [He]+[He+]+[He++]'
    WRITE(*,*) 'Please correct this : '

    call return_error(i_all,'x(He)','sum of [He]',&
                      self%myparams%He_abundance(1),&
                      self%myparams%x_HeI(1)+self%myparams%x_HeII(1)+&
                      self%myparams%x_HeIII(1))
  end if

  !rho_He = rho_Y

  self%myparams%densityHe(1) = self%myparams%He_abundance(1)*w_HeI*&
  self%myparams%number_H_density(1)*mh_gr

  if(DABS(self%myparams%densityHe(1)-&
  (self%myparams%HeliumFractionByMass(1)*&
  self%myparams%density_bar(1)))/self%myparams%densityHe(1)>abundance_density)then
    WRITE(*,*) 'Error: we have rho(He) /= Y*rho'
    WRITE(*,*) 'Please debug this : '

    call return_error(i_all,'rho(He)','Y*rho',&
                      self%myparams%densityHe(1),&
                      self%myparams%HeliumFractionByMass(1)*&
                      self%myparams%density_bar(1))
  end if

  self%myparams%density_Y(1)=self%myparams%densityHe(1)

  if(HeNotcorrected)then
    !correct-scale helium-whole abundances so that
    ! rhoY = rhoHe
    ! = rhoHeI + rhoHeII + rhoHeIII

    !first compute uncorrected values from abundances
    self%myparams%densityHeI(1) = w_HeI*self%myparams%x_HeI(1)*self%myparams%number_H_density(1)*mh_gr
    self%myparams%densityHeII(1) = w_HeII*self%myparams%x_HeII(1)*self%myparams%number_H_density(1)*mh_gr
    self%myparams%densityHeIII(1) = w_HeIII*self%myparams%x_HeIII(1)*self%myparams%number_H_density(1)*mh_gr

    !second compute corrected-scaled densities
    self%myparams%densityHeI(1) = self%myparams%densityHeI(1)*self%myparams%density_Y(1)/&
    (self%myparams%densityHeI(1)+self%myparams%densityHeII(1)+self%myparams%densityHeIII(1))
    self%myparams%densityHeII(1) = self%myparams%densityHeII(1)*self%myparams%density_X(1)/&
    (self%myparams%densityHeI(1)+self%myparams%densityHeII(1)+self%myparams%densityHeIII(1))
    self%myparams%densityHeIII(1) = self%myparams%densityHeIII(1)*self%myparams%density_X(1)/&
    (self%myparams%densityHeI(1)+self%myparams%densityHeII(1)+self%myparams%densityHeIII(1))


    !Third compute corrected-scaled abundances
    self%myparams%x_HeI(1) = self%myparams%densityHeI(1)/&
    (w_DI*self%myparams%number_H_density(1)*mh_gr)
    self%myparams%x_HeII(1) = self%myparams%densityHeII(1)/&
    (w_DI*self%myparams%number_H_density(1)*mh_gr)
    self%myparams%x_HeIII(1) = self%myparams%densityHeIII(1)/&
    (w_DI*self%myparams%number_H_density(1)*mh_gr)

    HeNotcorrected = .false.
  end if



  ! X, chiH
  if(self%myparams%chi_Hydrogen(1)>=-smalldouble)then

    !X
    self%myparams%HydrogenFractionByMass(1)=&
    (DABS(self%myparams%chi_Hydrogen(1))*(1.0_dp+&
    self%myparams%DeuteriumToHydrogenRatio(1))*&
    (w_HI*self%myparams%x_HI(1) + w_HII*self%myparams%x_HII(1) + w_HM*self%myparams%x_HM(1)+&
    w_H2I*self%myparams%x_H2I(1) + w_H2II*self%myparams%x_H2II(1)+&
    (1.0_dp-self%myparams%MetalFractionByMass(1))/&
    (1.0_dp-(self%myparams%HeliumFractionByMass(1)+&
    self%myparams%MetalFractionByMass(1))/&
    (1.0_dp-self%myparams%densityElectrons(1)/&
    (physical_ref_density-self%myparams%density_dust(1)))))*&
    self%myparams%number_H_density(1)*mh_gr)/&
    (physical_ref_density-self%myparams%density_dust(1)-&
    self%myparams%densityElectrons(1))

  elseif(self%myparams%HydrogenFractionByMass(1)>=-smalldouble)then

    !chi_H
    self%myparams%chi_Hydrogen(1)=&
    DABS(self%myparams%HydrogenFractionByMass(1))*&
    (physical_ref_density-self%myparams%density_dust(1)-&
    self%myparams%densityElectrons(1))/&
    ((1.0_dp+&
    self%myparams%DeuteriumToHydrogenRatio(1))*&
    (w_HI*self%myparams%x_HI(1) + w_HII*self%myparams%x_HII(1) + w_HM*self%myparams%x_HM(1)+&
    w_H2I*self%myparams%x_H2I(1) + w_H2II*self%myparams%x_H2II(1)+&
    (1.0_dp-self%myparams%MetalFractionByMass(1))/&
    (1.0_dp-(self%myparams%HeliumFractionByMass(1)+&
    self%myparams%MetalFractionByMass(1))/&
    (1.0_dp-self%myparams%densityElectrons(1)/&
    (physical_ref_density-self%myparams%density_dust(1)))))*&
    self%myparams%number_H_density(1)*mh_gr)


  else
    call return_error(i_all,'chi_H','X',&
                      self%myparams%chi_Hydrogen(1),&
                      self%myparams%HydrogenFractionByMass(1))
  end if

  !rhoX = rhoD + rhoHD + rhoH + rhoH2

  self%myparams%density_X(1) = self%myparams%HydrogenFractionByMass(1)*&
  self%myparams%density_bar(1)

  !rhometalfree = rhoX/(1+chiD)+rhoY
  self%myparams%density_metal_free(1) = self%myparams%density_X(1)/&
  (1.0_dp + self%myparams%DeuteriumToHydrogenRatio(1))+&
  self%myparams%density_Y(1)

  !intermediate rhoHpH2=rhoH + rhoH2 = chiH * rhometalfree
  self%myparams%densityHplusH2(1) = self%myparams%chi_Hydrogen(1)*&
  self%myparams%density_metal_free(1)

  !intermediate rhoDpHD=rhoD + rhoHD = chiD * rhoHpH2
  self%myparams%densityDplusHD(1) = self%myparams%DeuteriumToHydrogenRatio(1)*&
  self%myparams%densityHplusH2(1)

  !rhodeut =rhoDpHD
  self%myparams%density_deut(1) = self%myparams%densityDplusHD(1)

  !rhonotdeut
  self%myparams%density_not_deut(1) = self%myparams%density_bar(1)-&
  self%myparams%density_deut(1)

  ! Parameters default values
  ! Gas to dust ratios
  self%myparams%myindice(1:max_num_parameters) = 0

  !rhoD,rhoHD
  self%myparams%densityD(1)=(w_DI*self%myparams%x_DI(1)+&
  w_DII*self%myparams%x_DII(1))*self%myparams%number_H_density(1)*mh_gr
  self%myparams%densityHD(1)=w_HDI*self%myparams%x_HDI(1)*&
  self%myparams%number_H_density(1)*mh_gr

  if(HNotcorrected)then
    !correct-scale densities of D+D^+ and HD so that (rhoD+rhoHD)=chiD*(rhoH+rhoH2)
    self%myparams%densityD(1)=self%myparams%densityD(1)*&
    (self%myparams%DeuteriumToHydrogenRatio(1)*&
    self%myparams%densityHplusH2(1))/self%myparams%densityDplusHD(1)
    self%myparams%densityHD(1)=self%myparams%densityHD(1)*&
    (self%myparams%DeuteriumToHydrogenRatio(1)*&
    self%myparams%densityHplusH2(1))/self%myparams%densityDplusHD(1)
  end if

  !same for rhoH and rhoH2
  self%myparams%densityH(1)=(w_HI*self%myparams%x_HI(1)+&
  w_HII*self%myparams%x_HII(1)+w_HM*self%myparams%x_HM(1))*&
  self%myparams%number_H_density(1)*mh_gr
  self%myparams%densityHtwo(1)=(w_H2I*self%myparams%x_H2I(1)+&
  w_H2II*self%myparams%x_H2II(1))*&
  self%myparams%number_H_density(1)*mh_gr

  if(HNotcorrected)then
    !correct-scale densities of H+H^++H^- and H2+H2^+ so that (rhoH+rhoH2)=chiH*rhometalfree
    self%myparams%densityH(1)=self%myparams%densityH(1)*&
    (self%myparams%chi_Hydrogen(1)*&
    self%myparams%density_metal_free(1))/self%myparams%densityHplusH2(1)
    self%myparams%densityHtwo(1)=self%myparams%densityHtwo(1)*&
    (self%myparams%chi_Hydrogen(1)*&
    self%myparams%density_metal_free(1))/self%myparams%densityHplusH2(1)
  end if

  !correct-scale hydrogen-whole abundances so that
  ! rhoX = (rhoD+rhoHD+rhoH+rhoH2)
  ! = rhoDI + rhoDII + rhoHDI + rhoHI + rhoHII + rhoHM + rhoH2I + rhoH2II

  !first compute uncorrected values from abundances
  self%myparams%densityDI(1) = w_DI*self%myparams%x_DI(1)*self%myparams%number_H_density(1)*mh_gr
  self%myparams%densityDII(1) = w_DII*self%myparams%x_DII(1)*self%myparams%number_H_density(1)*mh_gr
  self%myparams%densityHDI(1) = w_HDI*self%myparams%x_HDI(1)*self%myparams%number_H_density(1)*mh_gr
  self%myparams%densityHI(1) = w_HI*self%myparams%x_HI(1)*self%myparams%number_H_density(1)*mh_gr
  self%myparams%densityHII(1) = w_HII*self%myparams%x_HII(1)*self%myparams%number_H_density(1)*mh_gr
  self%myparams%densityHM(1) = w_HM*self%myparams%x_HM(1)*self%myparams%number_H_density(1)*mh_gr
  self%myparams%densityH2I(1) = w_H2I*self%myparams%x_H2I(1)*self%myparams%number_H_density(1)*mh_gr
  self%myparams%densityH2II(1) = w_H2II*self%myparams%x_H2II(1)*self%myparams%number_H_density(1)*mh_gr

  if(HNotcorrected)then
    !second compute corrected-scaled densities
    self%myparams%densityDI(1) = self%myparams%densityDI(1)*self%myparams%density_X(1)/&
    (self%myparams%densityDI(1)+self%myparams%densityDII(1)+self%myparams%densityHDI(1)+&
    self%myparams%densityHI(1)+self%myparams%densityHII(1)+self%myparams%densityHM(1)+&
    self%myparams%densityH2I(1)+self%myparams%densityH2II(1))
    self%myparams%densityDII(1) = self%myparams%densityDII(1)*self%myparams%density_X(1)/&
    (self%myparams%densityDI(1)+self%myparams%densityDII(1)+self%myparams%densityHDI(1)+&
    self%myparams%densityHI(1)+self%myparams%densityHII(1)+self%myparams%densityHM(1)+&
    self%myparams%densityH2I(1)+self%myparams%densityH2II(1))
    self%myparams%densityHDI(1) = self%myparams%densityHDI(1)*self%myparams%density_X(1)/&
    (self%myparams%densityDI(1)+self%myparams%densityDII(1)+self%myparams%densityHDI(1)+&
    self%myparams%densityHI(1)+self%myparams%densityHII(1)+self%myparams%densityHM(1)+&
    self%myparams%densityH2I(1)+self%myparams%densityH2II(1))
    self%myparams%densityHI(1) = self%myparams%densityHI(1)*self%myparams%density_X(1)/&
    (self%myparams%densityDI(1)+self%myparams%densityDII(1)+self%myparams%densityHDI(1)+&
    self%myparams%densityHI(1)+self%myparams%densityHII(1)+self%myparams%densityHM(1)+&
    self%myparams%densityH2I(1)+self%myparams%densityH2II(1))
    self%myparams%densityHII(1) = self%myparams%densityHII(1)*self%myparams%density_X(1)/&
    (self%myparams%densityDI(1)+self%myparams%densityDII(1)+self%myparams%densityHDI(1)+&
    self%myparams%densityHI(1)+self%myparams%densityHII(1)+self%myparams%densityHM(1)+&
    self%myparams%densityH2I(1)+self%myparams%densityH2II(1))
    self%myparams%densityHM(1) = self%myparams%densityHM(1)*self%myparams%density_X(1)/&
    (self%myparams%densityDI(1)+self%myparams%densityDII(1)+self%myparams%densityHDI(1)+&
    self%myparams%densityHI(1)+self%myparams%densityHII(1)+self%myparams%densityHM(1)+&
    self%myparams%densityH2I(1)+self%myparams%densityH2II(1))
    self%myparams%densityH2I(1) = self%myparams%densityH2I(1)*self%myparams%density_X(1)/&
    (self%myparams%densityDI(1)+self%myparams%densityDII(1)+self%myparams%densityHDI(1)+&
    self%myparams%densityHI(1)+self%myparams%densityHII(1)+self%myparams%densityHM(1)+&
    self%myparams%densityH2I(1)+self%myparams%densityH2II(1))
    self%myparams%densityH2II(1) = self%myparams%densityH2II(1)*self%myparams%density_X(1)/&
    (self%myparams%densityDI(1)+self%myparams%densityDII(1)+self%myparams%densityHDI(1)+&
    self%myparams%densityHI(1)+self%myparams%densityHII(1)+self%myparams%densityHM(1)+&
    self%myparams%densityH2I(1)+self%myparams%densityH2II(1))

    !Third compute corrected-scaled abundances
    self%myparams%x_DI(1) = self%myparams%densityDI(1)/&
    (w_DI*self%myparams%number_H_density(1)*mh_gr)
    self%myparams%x_DII(1) = self%myparams%densityDII(1)/&
    (w_DI*self%myparams%number_H_density(1)*mh_gr)
    self%myparams%x_HDI(1) = self%myparams%densityHDI(1)/&
    (w_DI*self%myparams%number_H_density(1)*mh_gr)
    self%myparams%x_HI(1) = self%myparams%densityHI(1)/&
    (w_DI*self%myparams%number_H_density(1)*mh_gr)
    self%myparams%x_HII(1) = self%myparams%densityHII(1)/&
    (w_DI*self%myparams%number_H_density(1)*mh_gr)
    self%myparams%x_HM(1) = self%myparams%densityHM(1)/&
    (w_DI*self%myparams%number_H_density(1)*mh_gr)
    self%myparams%x_H2I(1) = self%myparams%densityH2I(1)/&
    (w_DI*self%myparams%number_H_density(1)*mh_gr)
    self%myparams%x_H2II(1) = self%myparams%densityH2II(1)/&
    (w_DI*self%myparams%number_H_density(1)*mh_gr)
    !End of correction-scale for hydrogen-whole

    HNotcorrected = .false.
  end if

  end do






     select case(gr_density_method(gr_patches_indices_global(i_all)))
     !TO DO: other cases than chemical_coefficient

     case('chemical_coefficient')
       !DO NOTHING TO COMPLETE
     case default
       !DO NOTHING TO COMPLETE
     end select

  end do Loop_through_all

contains

  subroutine return_error(iobj,name1,name2,&
                           myvalue1,myvalue2)
    implicit none
    integer, intent(in)       :: iobj
    character(len=*),intent(in) :: name1,name2
    real(kind=dp), intent(in) :: myvalue1
    real(kind=dp), intent(in) :: myvalue2
    !---------------------------------------------
    WRITE(*,*) 'In object # ', iobj, ', we have : '
    WRITE(*,*)  name1, ' = ', myvalue1,' and ', name2, ' = ',myvalue2
    call mpistop('Inconsistent parameters for Grackle')
  end subroutine return_error

end subroutine grackle_set_complet





end module mod_grackle_chemistry
