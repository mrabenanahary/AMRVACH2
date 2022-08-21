!> module radiative cooling -- add optically thin radiative cooling for HD and MHD
!>
!> Assumptions: full ionize plasma dominated by H and He, ionization equilibrium
!> Formula: Q=-n_H*n_e*f(T), positive f(T) function is pre-computed and tabulated
!> Developed by Allard Jan van Marle, Rony Keppens, and Chun Xia
!> Cooling tables extend to 1O^9 K, (17.11.2009)  AJvM
!> Included table by Smith (18.11.2009)           AJvM
!> Included luminosity output option (18.11.2009) AJvM
!> Included two cooling curves from Cloudy code supplied by Wang Ye (12.11.2011) AJvM
!> Included a solar coronal cooling curve (JCcorona) supplied by J. Colgan (2008) ApJ
!> Modulized and simplified by Chun Xia (2017)
module mod_radiative_cooling
!
! these tables contain log_10 temperature values and corresponding
! log_10 luminosity values. The simulation-dependent temperature and luminosity
! scaling parameters are supposed to be provided in the user file.
! All tables have been extended to at least T=10^9 K using a pure Bremmstrahlung
! relationship of Lambda~sqrt(T). This to ensure that a purely explicit calculation
! without timestep check is only used for extremely high temperatures.
! (Except for the SPEX curve, which is more complicated and therefore simply stops
! at the official upper limit of log(T) = 8.16)
!
  use mod_constants
  use mod_global_parameters, only: std_len
  use mod_physics
  use mod_chemical
  implicit none
  ! parameters used for implicit cooling source calculations

  !> Coefficent of cooling time step
  real(kind=dp)   , private   :: cfrac

  !> Lower limit of temperature
  real(kind=dp)   , private   :: tlow

  !> Helium abundance over Hydrogen
  real(kind=dp)   , private    :: He_abundance

  !> Chemical_gas_type
  character(len=std_len), private :: rc_chemical_gas_type

  real(kind=dp)    , private               :: mean_mass
  real(kind=dp)    , private               :: mean_mup
  real(kind=dp)    , private               :: mean_ne_to_nH
  real(kind=dp)    , private               :: mean_nall_to_nH


  !> resolution of temperature in interpolated tables
  integer, private :: ncool

  !> Name of cooling curve
  character(len=std_len), private  :: coolcurve

  !> Name of cooling method
  character(len=std_len), private  :: coolmethod

  !> Fixed temperature not lower than tlow
  logical, private    :: Tfix

  !> Add cooling source in a split way (.true.) or un-split way (.false.)
  logical :: rc_split=.false.

  !> Index of the density (in the w array)
  integer, private, protected              :: rho_

  !> Indices of the momentum density
  integer, allocatable, private, protected :: mom(:)

  !> Index of the energy density
  integer, private, protected              :: e_

  !> Index of cooling erg/s
  integer, private, protected              :: Lcool1_

  !> Index of cooling time
  integer, private, protected              :: dtcool1_

  !!> Index of cooling time
  !integer, private, protected              :: mucool1_

  logical, private, protected              :: coolsaveL
  logical, private, protected              :: coolsavedT
  logical, private, protected              :: coolsavemu
  !> inonisation rate
  !real(kind=dp), private, protected        :: coolfn=1.0d-3
  type cooling_config
    logical                  :: saveL
    logical                  :: savedT
    logical                  :: savemu
      !> Coefficent of cooling time step
      real(kind=dp)           :: cfrac

      !> Lower limit of temperature
      real(kind=dp)           :: tlow

      !> Helium abundance over Hydrogen
      real(kind=dp)           :: He_abundance

      !> Chemical_gas_type
      character(len=std_len) :: rc_chemical_gas_type

      real(kind=dp)                   :: mean_mass
      real(kind=dp)                   :: mean_mup
      real(kind=dp)                   :: mean_ne_to_nH
      real(kind=dp)                   :: mean_nall_to_nH

      !> resolution of temperature in interpolated tables
      integer                 :: npoint

      !> Name of cooling curve
      character(len=std_len)  :: curve

      !> Name of cooling method
      character(len=std_len)  :: method

      !> Fixed temperature not lower than tlow
      logical, private        :: Tfix

      !> Add cooling source in a split way (.true.) or un-split way (.false.)
      logical                 :: rc_split


  end type cooling_config
  type(cooling_config), public, protected :: coolconfig
  !> The adiabatic index
  real(kind=dp)   , private :: rc_gamma
  real(kind=dp)   , private :: unit_luminosity
  real(kind=dp)   , allocatable :: tcool(:), Lcool(:), dLdtcool(:)
  real(kind=dp)   , allocatable :: Yc(:), invYc(:)
  real(kind=dp)     :: Tref, Lref, tcoolmin,tcoolmax
  real(kind=dp)     :: lgtcoolmin, lgtcoolmax, lgstep

  integer          :: n_DM      , n_MB      , n_MLcosmol &
                    , n_MLwc    , n_MLsolar1, n_SPEX     &
                    , n_JCcorona, n_cl_ism  , n_cl_solar &
                    , n_DM_2

  real(kind=dp)    :: t_DM(1:71),       t_MB(1:51),       t_MLcosmol(1:71) &
                    , t_MLwc(1:71),     t_MLsolar1(1:71), t_SPEX(1:110)    &
                    , t_JCcorona(1:45), t_cl_ism(1:151),  t_cl_solar(1:151)&
                    , t_DM_2(1:76)

  real(kind=dp)    :: l_DM(1:71),       l_MB(1:51),       l_MLcosmol(1:71) &
                    , l_MLwc(1:71),     l_MLsolar1(1:71), l_SPEX(1:110)    &
                    , l_JCcorona(1:45), l_cl_ism(1:151),  l_cl_solar(1:151)&
                    , l_DM_2(1:76)

  real(kind=dp)    :: nenh_SPEX(1:110)

  data    n_JCcorona / 45 /

  data    t_JCcorona / 4.00000, 4.14230, 4.21995, 4.29761, 4.37528, &
                       4.45294, 4.53061, 4.60827, 4.68593, 4.76359, &
                       4.79705, 4.83049, 4.86394, 4.89739, 4.93084, &
                       4.96428, 4.99773, 5.03117, 5.06461, 5.17574, &
                       5.28684, 5.39796, 5.50907, 5.62018, 5.73129, &
                       5.84240, 5.95351, 6.06461, 6.17574, 6.28684, &
                       6.39796, 6.50907, 6.62018, 6.73129, 6.84240, &
                       6.95351, 7.06461, 7.17574, 7.28684, 7.39796, &
                       7.50907, 7.62018, 7.73129, 7.84240, 7.95351  /

  data    l_JCcorona / -200.18883, -100.78630, -30.60384, -22.68481, -21.76445, &
                       -21.67936, -21.54218, -21.37958, -21.25172, -21.17584, &
                       -21.15783, -21.14491, -21.13527, -21.12837, -21.12485, &
                       -21.12439, -21.12642, -21.12802, -21.12548, -21.08965, &
                       -21.08812, -21.19542, -21.34582, -21.34839, -21.31701, &
                       -21.29072, -21.28900, -21.34104, -21.43122, -21.62448, &
                       -21.86694, -22.02897, -22.08051, -22.06057, -22.01973, &
                       -22.00000, -22.05161, -22.22175, -22.41452, -22.52581, &
                       -22.56914, -22.57486, -22.56151, -22.53969, -22.51490  /

  data    n_DM / 71 /

  data    t_DM / 2.0, 2.1, 2.2, 2.3, 2.4, &
                 2.5, 2.6, 2.7, 2.8, 2.9, &
                 3.0, 3.1, 3.2, 3.3, 3.4, &
                 3.5, 3.6, 3.7, 3.8, 3.9, &
                 4.0, 4.1, 4.2, 4.3, 4.4, &
                 4.5, 4.6, 4.7, 4.8, 4.9, &
                 5.0, 5.1, 5.2, 5.3, 5.4, &
                 5.5, 5.6, 5.7, 5.8, 5.9, &
                 6.0, 6.1, 6.2, 6.3, 6.4, &
                 6.5, 6.6, 6.7, 6.8, 6.9, &
                 7.0, 7.1, 7.2, 7.3, 7.4, &
                 7.5, 7.6, 7.7, 7.8, 7.9, &
                 8.0, 8.1, 8.2, 8.3, 8.4, &
                 8.5, 8.6, 8.7, 8.8, 8.9, &
                 9.0 /


  data    l_DM / -26.523, -26.398, -26.301, -26.222, -26.097 &
               , -26.011, -25.936, -25.866, -25.807, -25.754 &
               , -25.708, -25.667, -25.630, -25.595, -25.564 &
               , -25.534, -25.506, -25.479, -25.453, -25.429 &
               , -25.407, -23.019, -21.762, -21.742, -21.754 &
               , -21.730, -21.523, -21.455, -21.314, -21.229 &
               , -21.163, -21.126, -21.092, -21.060, -21.175 &
               , -21.280, -21.390, -21.547, -21.762, -22.050 &
               , -22.271, -22.521, -22.646, -22.660, -22.676 &
               , -22.688, -22.690, -22.662, -22.635, -22.609 &
               , -22.616, -22.646, -22.697, -22.740, -22.788 &
               , -22.815, -22.785, -22.754, -22.728, -22.703 &
               , -22.680, -22.630, -22.580, -22.530, -22.480 &
               , -22.430, -22.380, -22.330, -22.280, -22.230 &
               , -22.180 /

  data    n_MB / 51 /

  data    t_MB / 4.0, 4.1, 4.2, 4.3, 4.4, &
                 4.5, 4.6, 4.7, 4.8, 4.9, &
                 5.0, 5.1, 5.2, 5.3, 5.4, &
                 5.5, 5.6, 5.7, 5.8, 5.9, &
                 6.0, 6.1, 6.2, 6.3, 6.4, &
                 6.5, 6.6, 6.7, 6.8, 6.9, &
                 7.0, 7.1, 7.2, 7.3, 7.4, &
                 7.5, 7.6, 7.7, 7.8, 7.9, &
                 8.0, 8.1, 8.2, 8.3, 8.4, &
                 8.5, 8.6, 8.7, 8.8, 8.9, &
                 9.0 /

  data    l_MB  / -23.133, -22.895, -22.548, -22.285, -22.099 &
                , -21.970, -21.918, -21.826, -21.743, -21.638 &
                , -21.552, -21.447, -21.431, -21.418, -21.461 &
                , -21.570, -21.743, -21.832, -21.908, -21.981 &
                , -22.000, -21.998, -21.992, -22.013, -22.095 &
                , -22.262, -22.397, -22.445, -22.448, -22.446 &
                , -22.448, -22.465, -22.575, -22.725, -22.749 &
                , -22.768, -22.753, -22.717, -22.678, -22.637 &
                , -22.603, -22.553, -22.503, -22.453, -22.403 &
                , -22.353, -22.303, -22.253, -22.203, -22.153 &
                , -22.103 /

  data    n_MLcosmol / 71 /

  data    t_MLcosmol / &
                 2.0, 2.1, 2.2, 2.3, 2.4, &
                 2.5, 2.6, 2.7, 2.8, 2.9, &
                 3.0, 3.1, 3.2, 3.3, 3.4, &
                 3.5, 3.6, 3.7, 3.8, 3.9, &
                 4.0, 4.1, 4.2, 4.3, 4.4, &
                 4.5, 4.6, 4.7, 4.8, 4.9, &
                 5.0, 5.1, 5.2, 5.3, 5.4, &
                 5.5, 5.6, 5.7, 5.8, 5.9, &
                 6.0, 6.1, 6.2, 6.3, 6.4, &
                 6.5, 6.6, 6.7, 6.8, 6.9, &
                 7.0, 7.1, 7.2, 7.3, 7.4, &
                 7.5, 7.6, 7.7, 7.8, 7.9, &
                 8.0, 8.1, 8.2, 8.3, 8.4, &
                 8.5, 8.6, 8.7, 8.8, 8.9, &
                 9.0 /

  data    l_MLcosmol / &
                 -99.000, -99.000, -99.000, -99.000, -99.000 &
               , -99.000, -99.000, -99.000, -99.000, -99.000 &
               , -99.000, -99.000, -99.000, -99.000, -99.000 &
               , -99.000, -44.649, -38.362, -33.324, -29.292 &
               , -26.063, -23.532, -22.192, -22.195, -22.454 &
               , -22.676, -22.909, -22.925, -22.499, -22.276 &
               , -22.440, -22.688, -22.917, -23.116, -23.274 &
               , -23.394, -23.472, -23.516, -23.530, -23.525 &
               , -23.506, -23.478, -23.444, -23.408, -23.368 &
               , -23.328, -23.286, -23.244, -23.201, -23.157 &
               , -23.114, -23.070, -23.026, -22.981, -22.937 &
               , -22.893, -22.848, -22.803, -22.759, -22.714 &
               , -22.669, -22.619, -22.569, -22.519, -22.469 &
               , -22.419, -22.369, -22.319, -22.269, -22.190 &
               , -22.169 /

  data    n_MLwc / 71 /

  data    t_MLwc / &
                 2.0, 2.1, 2.2, 2.3, 2.4, &
                 2.5, 2.6, 2.7, 2.8, 2.9, &
                 3.0, 3.1, 3.2, 3.3, 3.4, &
                 3.5, 3.6, 3.7, 3.8, 3.9, &
                 4.0, 4.1, 4.2, 4.3, 4.4, &
                 4.5, 4.6, 4.7, 4.8, 4.9, &
                 5.0, 5.1, 5.2, 5.3, 5.4, &
                 5.5, 5.6, 5.7, 5.8, 5.9, &
                 6.0, 6.1, 6.2, 6.3, 6.4, &
                 6.5, 6.6, 6.7, 6.8, 6.9, &
                 7.0, 7.1, 7.2, 7.3, 7.4, &
                 7.5, 7.6, 7.7, 7.8, 7.9, &
                 8.0, 8.1, 8.2, 8.3, 8.4, &
                 8.5, 8.6, 8.7, 8.8, 8.9, &
                 9.0 /

  data    l_MLwc / &
                 -21.192, -21.160, -21.150, -21.150, -21.166 &
               , -21.191, -21.222, -21.264, -21.308, -21.357 &
               , -21.408, -21.449, -21.494, -21.544, -21.587 &
               , -21.638, -21.686, -21.736, -21.780, -21.800 &
               , -21.744, -21.547, -21.208, -20.849, -20.345 &
               , -19.771, -19.409, -19.105, -18.827, -18.555 &
               , -18.460, -18.763, -19.168, -19.334, -19.400 &
               , -19.701, -20.090, -20.288, -20.337, -20.301 &
               , -20.233, -20.275, -20.363, -20.508, -20.675 &
               , -20.856, -21.025, -21.159, -21.256, -21.320 &
               , -21.354, -21.366, -21.361, -21.343, -21.317 &
               , -21.285, -21.250, -21.212, -21.172, -21.131 &
               , -21.089, -21.039, -20.989, -20.939, -20.889 &
               , -20.839, -20.789, -20.739, -20.689, -20.639 &
               , -20.589 /

  data n_MLsolar1 / 71 /

  data    t_MLsolar1 / &
                       2.0, 2.1, 2.2, 2.3, 2.4, &
                       2.5, 2.6, 2.7, 2.8, 2.9, &
                        3.0, 3.1, 3.2, 3.3, 3.4, &
                       3.5, 3.6, 3.7, 3.8, 3.9, &
                       4.0, 4.1, 4.2, 4.3, 4.4, &
                       4.5, 4.6, 4.7, 4.8, 4.9, &
                       5.0, 5.1, 5.2, 5.3, 5.4, &
                       5.5, 5.6, 5.7, 5.8, 5.9, &
                       6.0, 6.1, 6.2, 6.3, 6.4, &
                       6.5, 6.6, 6.7, 6.8, 6.9, &
                       7.0, 7.1, 7.2, 7.3, 7.4, &
                       7.5, 7.6, 7.7, 7.8, 7.9, &
                       8.0, 8.1, 8.2, 8.3, 8.4, &
                       8.5, 8.6, 8.7, 8.8, 8.9, &
                       9.0 /

  data    l_MLsolar1 / &
                       -26.983, -26.951, -26.941, -26.940, -26.956 &
                     , -26.980, -27.011, -27.052, -27.097, -27.145 &
                     , -27.195, -27.235, -27.279, -27.327, -27.368 &
                     , -27.415, -27.456, -27.485, -27.468, -27.223 &
                     , -25.823, -23.501, -22.162, -22.084, -22.157 &
                     , -22.101, -21.974, -21.782, -21.542, -21.335 &
                     , -21.251, -21.275, -21.236, -21.173, -21.167 &
                     , -21.407, -21.670, -21.788, -21.879, -22.008 &
                     , -22.192, -22.912, -22.918, -22.887, -22.929 &
                     , -23.023, -23.094, -23.117, -23.108, -23.083 &
                     , -23.049, -23.011, -22.970, -22.928, -22.885 &
                     , -22.842, -22.798, -22.754, -22.709, -22.665 &
                     , -22.620, -22.570, -22.520, -22.470, -22.420 &
                     , -22.370, -22.320, -22.270, -22.220, -22.170 &
                     , -22.120 /


  data n_SPEX / 110 /

  data    t_SPEX / &
                       3.80, 3.84, 3.88, 3.92, 3.96, &
                       4.00, 4.04, 4.08, 4.12, 4.16, &
                       4.20, 4.24, 4.28, 4.32, 4.36, &
                       4.40, 4.44, 4.48, 4.52, 4.56, &
                       4.60, 4.64, 4.68, 4.72, 4.76, &
                       4.80, 4.84, 4.88, 4.92, 4.96, &
                       5.00, 5.04, 5.08, 5.12, 5.16, &
                       5.20, 5.24, 5.28, 5.32, 5.36, &
                       5.40, 5.44, 5.48, 5.52, 5.56, &
                       5.60, 5.64, 5.68, 5.72, 5.76, &
                       5.80, 5.84, 5.88, 5.92, 5.96, &
                       6.00, 6.04, 6.08, 6.12, 6.16, &
                       6.20, 6.24, 6.28, 6.32, 6.36, &
                       6.40, 6.44, 6.48, 6.52, 6.56, &
                       6.60, 6.64, 6.68, 6.72, 6.76, &
                       6.80, 6.84, 6.88, 6.92, 6.96, &
                       7.00, 7.04, 7.08, 7.12, 7.16, &
                       7.20, 7.24, 7.28, 7.32, 7.36, &
                       7.40, 7.44, 7.48, 7.52, 7.56, &
                       7.60, 7.64, 7.68, 7.72, 7.76, &
                       7.80, 7.84, 7.88, 7.92, 7.96, &
                       8.00, 8.04, 8.08, 8.12, 8.16  /

  data    l_SPEX / &
                      -25.7331, -25.0383, -24.4059, -23.8288, -23.3027 &
                    , -22.8242, -22.3917, -22.0067, -21.6818, -21.4529 &
                    , -21.3246, -21.3459, -21.4305, -21.5293, -21.6138 &
                    , -21.6615, -21.6551, -21.5919, -21.5092, -21.4124 &
                    , -21.3085, -21.2047, -21.1067, -21.0194, -20.9413 &
                    , -20.8735, -20.8205, -20.7805, -20.7547, -20.7455 &
                    , -20.7565, -20.7820, -20.8008, -20.7994, -20.7847 &
                    , -20.7687, -20.7590, -20.7544, -20.7505, -20.7545 &
                    , -20.7888, -20.8832, -21.0450, -21.2286, -21.3737 &
                    , -21.4573, -21.4935, -21.5098, -21.5345, -21.5863 &
                    , -21.6548, -21.7108, -21.7424, -21.7576, -21.7696 &
                    , -21.7883, -21.8115, -21.8303, -21.8419, -21.8514 &
                    , -21.8690, -21.9057, -21.9690, -22.0554, -22.1488 &
                    , -22.2355, -22.3084, -22.3641, -22.4033, -22.4282 &
                    , -22.4408, -22.4443, -22.4411, -22.4334, -22.4242 &
                    , -22.4164, -22.4134, -22.4168, -22.4267, -22.4418 &
                    , -22.4603, -22.4830, -22.5112, -22.5449, -22.5819 &
                    , -22.6177, -22.6483, -22.6719, -22.6883, -22.6985 &
                    , -22.7032, -22.7037, -22.7008, -22.6950, -22.6869 &
                    , -22.6769, -22.6655, -22.6531, -22.6397, -22.6258 &
                    , -22.6111, -22.5964, -22.5816, -22.5668, -22.5519 &
                    , -22.5367, -22.5216, -22.5062, -22.4912, -22.4753 /


  data    nenh_SPEX / &
            0.000013264, 0.000042428, 0.000088276, 0.00017967 &
          , 0.00084362,  0.0034295,   0.013283,    0.042008   &
          , 0.12138,     0.30481,     0.53386,     0.76622    &
          , 0.89459,     0.95414,     0.98342 &
          , 1.0046,    1.0291,    1.0547    &
          , 1.0767,    1.0888    &
          , 1.0945,    1.0972,    1.0988    &
          , 1.1004,    1.1034    &
          , 1.1102,    1.1233,    1.1433    &
          , 1.1638,    1.1791    &
          , 1.1885,    1.1937,    1.1966    &
          , 1.1983,    1.1993    &
          , 1.1999,    1.2004,    1.2008    &
          , 1.2012,    1.2015    &
          , 1.2020,    1.2025,    1.2030    &
          , 1.2035,    1.2037    &
          , 1.2039,    1.2040,    1.2041    &
          , 1.2042,    1.2044    &
          , 1.2045,    1.2046,    1.2047    &
          , 1.2049,    1.2050    &
          , 1.2051,    1.2053,    1.2055    &
          , 1.2056,    1.2058    &
          , 1.2060,    1.2062,    1.2065    &
          , 1.2067,    1.2070    &
          , 1.2072,    1.2075,    1.2077    &
          , 1.2078,    1.2079    &
          , 1.2080,    1.2081,    1.2082    &
          , 1.2083,    1.2083    &
          , 1.2084,    1.2084,    1.2085    &
          , 1.2085,    1.2086    &
          , 1.2086,    1.2087,    1.2087    &
          , 1.2088,    1.2088    &
          , 1.2089,    1.2089,    1.2089    &
          , 1.2089,    1.2089    &
          , 1.2090,    1.2090,    1.2090    &
          , 1.2090,    1.2090    &
          , 1.2090,    1.2090,    1.2090    &
          , 1.2090,    1.2090    &
          , 1.2090,    1.2090,    1.2090    &
          , 1.2090,    1.2090    &
          , 1.2090,    1.2090,    1.2090    &
          , 1.2090,    1.2090    /
  !
  ! To be used together with the SPEX table for the SPEX_DM option
  ! Assuming an ionization fraction of 10^-3
  !
  data    n_DM_2 / 76 /

  data    t_DM_2 / 1.00, 1.04, 1.08, 1.12, 1.16, 1.20 &
                 , 1.24, 1.28, 1.32, 1.36, 1.40 &
                 , 1.44, 1.48, 1.52, 1.56, 1.60 &
                 , 1.64, 1.68, 1.72, 1.76, 1.80 &
                 , 1.84, 1.88, 1.92, 1.96, 2.00 &
                 , 2.04, 2.08, 2.12, 2.16, 2.20 &
                 , 2.24, 2.28, 2.32, 2.36, 2.40 &
                 , 2.44, 2.48, 2.52, 2.56, 2.60 &
                 , 2.64, 2.68, 2.72, 2.76, 2.80 &
                 , 2.84, 2.88, 2.92, 2.96, 3.00 &
                 , 3.04, 3.08, 3.12, 3.16, 3.20 &
                 , 3.24, 3.28, 3.32, 3.36, 3.40 &
                 , 3.44, 3.48, 3.52, 3.56, 3.60 &
                 , 3.64, 3.68, 3.72, 3.76, 3.80 &
                 , 3.84, 3.88, 3.92, 3.96, 4.00 /


  data    l_DM_2 / -30.0377, -29.7062, -29.4055, -29.1331, -28.8864, -28.6631 &
                 , -28.4614, -28.2791, -28.1146, -27.9662, -27.8330 &
                 , -27.7129, -27.6052, -27.5088, -27.4225, -27.3454 &
                 , -27.2767, -27.2153, -27.1605, -27.1111, -27.0664 &
                 , -27.0251, -26.9863, -26.9488, -26.9119, -26.8742 &
                 , -26.8353, -26.7948, -26.7523, -26.7080, -26.6619 &
                 , -26.6146, -26.5666, -26.5183, -26.4702, -26.4229 &
                 , -26.3765, -26.3317, -26.2886, -26.2473, -26.2078 &
                 , -26.1704, -26.1348, -26.1012, -26.0692, -26.0389 &
                 , -26.0101, -25.9825, -25.9566, -25.9318, -25.9083 &
                 , -25.8857, -25.8645, -25.8447, -25.8259, -25.8085 &
                 , -25.7926, -25.7778, -25.7642, -25.7520, -25.7409 &
                 , -25.7310, -25.7222, -25.7142, -25.7071, -25.7005 &
                 , -25.6942, -25.6878, -25.6811, -25.6733, -25.6641 &
                 , -25.6525, -25.6325, -25.6080, -25.5367, -25.4806  /

  data n_cl_solar / 151 /

  data    t_cl_solar /    &
     1.000   ,   1.050   ,   1.100   ,   1.150   , &
     1.200   ,   1.250   ,   1.300   ,   1.350   , &
     1.400   ,   1.450   ,   1.500   ,   1.550   , &
     1.600   ,   1.650   ,   1.700   ,   1.750   , &
     1.800   ,   1.850   ,   1.900   ,   1.950   , &
     2.000   ,   2.050   ,   2.100   ,   2.150   , &
     2.200   ,   2.250   ,   2.300   ,   2.350   , &
     2.400   ,   2.450   ,   2.500   ,   2.550   , &
     2.600   ,   2.650   ,   2.700   ,   2.750   , &
     2.800   ,   2.850   ,   2.900   ,   2.950   , &
     3.000   ,   3.050   ,   3.100   ,   3.150   , &
     3.200   ,   3.250   ,   3.300   ,   3.350   , &
     3.400   ,   3.450   ,   3.500   ,   3.550   , &
     3.600   ,   3.650   ,   3.700   ,   3.750   , &
     3.800   ,   3.850   ,   3.900   ,   3.950   , &
     4.000   ,   4.050   ,   4.100   ,   4.150   , &
     4.200   ,   4.250   ,   4.300   ,   4.350   , &
     4.400   ,   4.450   ,   4.500   ,   4.550   , &
     4.600   ,   4.650   ,   4.700   ,   4.750   , &
     4.800   ,   4.850   ,   4.900   ,   4.950   , &
     5.000   ,   5.050   ,   5.100   ,   5.150   , &
     5.200   ,   5.250   ,   5.300   ,   5.350   , &
     5.400   ,   5.450   ,   5.500   ,   5.550   , &
     5.600   ,   5.650   ,   5.700   ,   5.750   , &
     5.800   ,   5.850   ,   5.900   ,   5.950   , &
     6.000   ,   6.050   ,   6.100   ,   6.150   , &
     6.200   ,   6.250   ,   6.300   ,   6.350   , &
     6.400   ,   6.450   ,   6.500   ,   6.550   , &
     6.600   ,   6.650   ,   6.700   ,   6.750   , &
     6.800   ,   6.850   ,   6.900   ,   6.950   , &
     7.000   ,   7.050   ,   7.100   ,   7.150   , &
     7.200   ,   7.250   ,   7.300   ,   7.350   , &
     7.400   ,   7.450   ,   7.500   ,   7.550   , &
     7.600   ,   7.650   ,   7.700   ,   7.750   , &
     7.800   ,   7.850   ,   7.900   ,   7.950   , &
     8.000   ,   8.100   ,   8.200   ,   8.300   , &
     8.400   ,   8.500   ,   8.600   ,   8.700   , &
     8.800   ,   8.900   ,   9.00  /


  data    l_cl_solar / &
   -28.375   , -28.251   , -28.137   , -28.029   , &
   -27.929   , -27.834   , -27.745   , -27.662   , &
   -27.584   , -27.512   , -27.445   , -27.383   , &
   -27.326   , -27.273   , -27.223   , -27.175   , &
   -27.128   , -27.079   , -27.027   , -26.972   , &
   -26.911   , -26.846   , -26.777   , -26.705   , &
   -26.632   , -26.554   , -26.479   , -26.407   , &
   -26.338   , -26.274   , -26.213   , -26.156   , &
   -26.101   , -26.049   , -25.999   , -25.949   , &
   -25.901   , -25.852   , -25.803   , -25.754   , &
   -25.707   , -25.662   , -25.621   , -25.588   , &
   -25.561   , -25.538   , -25.518   , -25.497   , &
   -25.475   , -25.452   , -25.426   , -25.400   , &
   -25.374   , -25.333   , -25.295   , -25.261   , &
   -25.228   , -25.189   , -25.136   , -25.053   , &
   -24.888   , -24.454   , -23.480   , -22.562   , &
   -22.009   , -21.826   , -21.840   , -21.905   , &
   -21.956   , -21.971   , -21.958   , -21.928   , &
   -21.879   , -21.810   , -21.724   , -21.623   , &
   -21.512   , -21.404   , -21.321   , -21.273   , &
   -21.250   , -21.253   , -21.275   , -21.287   , &
   -21.282   , -21.275   , -21.272   , -21.267   , &
   -21.281   , -21.357   , -21.496   , -21.616   , &
   -21.677   , -21.698   , -21.708   , -21.730   , &
   -21.767   , -21.793   , -21.794   , -21.787   , &
   -21.787   , -21.802   , -21.826   , -21.859   , &
   -21.911   , -21.987   , -22.082   , -22.173   , &
   -22.253   , -22.325   , -22.392   , -22.448   , &
   -22.487   , -22.512   , -22.524   , -22.528   , &
   -22.524   , -22.516   , -22.507   , -22.501   , &
   -22.502   , -22.511   , -22.533   , -22.565   , &
   -22.600   , -22.630   , -22.648   , -22.656   , &
   -22.658   , -22.654   , -22.647   , -22.634   , &
   -22.619   , -22.602   , -22.585   , -22.566   , &
   -22.546   , -22.525   , -22.505   , -22.480   , &
   -22.465   , -22.415   , -22.365   , -22.315   , &
   -22.265   , -22.215   , -22.165   , -22.115   , &
   -22.065   , -22.015   , -21.965   /


  data n_cl_ism / 151 /

  data t_cl_ism / &
     1.000   ,   1.050   ,   1.100   ,   1.150   , &
     1.200   ,   1.250   ,   1.300   ,   1.350   , &
     1.400   ,   1.450   ,   1.500   ,   1.550   , &
     1.600   ,   1.650   ,   1.700   ,   1.750   , &
     1.800   ,   1.850   ,   1.900   ,   1.950   , &
     2.000   ,   2.050   ,   2.100   ,   2.150   , &
     2.200   ,   2.250   ,   2.300   ,   2.350   , &
     2.400   ,   2.450   ,   2.500   ,   2.550   , &
     2.600   ,   2.650   ,   2.700   ,   2.750   , &
     2.800   ,   2.850   ,   2.900   ,   2.950   , &
     3.000   ,   3.050   ,   3.100   ,   3.150   , &
     3.200   ,   3.250   ,   3.300   ,   3.350   , &
     3.400   ,   3.450   ,   3.500   ,   3.550   , &
     3.600   ,   3.650   ,   3.700   ,   3.750   , &
     3.800   ,   3.850   ,   3.900   ,   3.950   , &
     4.000   ,   4.050   ,   4.100   ,   4.150   , &
     4.200   ,   4.250   ,   4.300   ,   4.350   , &
     4.400   ,   4.450   ,   4.500   ,   4.550   , &
     4.600   ,   4.650   ,   4.700   ,   4.750   , &
     4.800   ,   4.850   ,   4.900   ,   4.950   , &
     5.000   ,   5.050   ,   5.100   ,   5.150   , &
     5.200   ,   5.250   ,   5.300   ,   5.350   , &
     5.400   ,   5.450   ,   5.500   ,   5.550   , &
     5.600   ,   5.650   ,   5.700   ,   5.750   , &
     5.800   ,   5.850   ,   5.900   ,   5.950   , &
     6.000   ,   6.050   ,   6.100   ,   6.150   , &
     6.200   ,   6.250   ,   6.300   ,   6.350   , &
     6.400   ,   6.450   ,   6.500   ,   6.550   , &
     6.600   ,   6.650   ,   6.700   ,   6.750   , &
     6.800   ,   6.850   ,   6.900   ,   6.950   , &
     7.000   ,   7.050   ,   7.100   ,   7.150   , &
     7.200   ,   7.250   ,   7.300   ,   7.350   , &
     7.400   ,   7.450   ,   7.500   ,   7.550   , &
     7.600   ,   7.650   ,   7.700   ,   7.750   , &
     7.800   ,   7.850   ,   7.900   ,   7.950   , &
     8.000   ,   8.100   ,   8.200   ,   8.300   , &
     8.400   ,   8.500   ,   8.600   ,   8.700   , &
     8.800   ,   8.900   ,   9.00    /

  data l_cl_ism / &
   -28.365   , -28.242   , -28.127   , -28.020   , &
   -27.919   , -27.825   , -27.736   , -27.653   , &
   -27.575   , -27.504   , -27.437   , -27.376   , &
   -27.319   , -27.267   , -27.220   , -27.176   , &
   -27.134   , -27.095   , -27.058   , -27.021   , &
   -26.985   , -26.948   , -26.910   , -26.870   , &
   -26.827   , -26.775   , -26.721   , -26.664   , &
   -26.608   , -26.552   , -26.495   , -26.437   , &
   -26.378   , -26.317   , -26.255   , -26.190   , &
   -26.123   , -26.053   , -25.984   , -25.913   , &
   -25.847   , -25.786   , -25.736   , -25.702   , &
   -25.678   , -25.662   , -25.649   , -25.636   , &
   -25.621   , -25.604   , -25.587   , -25.571   , &
   -25.562   , -25.526   , -25.505   , -25.499   , &
   -25.499   , -25.491   , -25.468   , -25.410   , &
   -25.268   , -24.888   , -23.702   , -22.624   , &
   -22.036   , -21.843   , -21.854   , -21.924   , &
   -21.986   , -22.017   , -22.021   , -22.005   , &
   -21.964   , -21.896   , -21.806   , -21.699   , &
   -21.580   , -21.463   , -21.370   , -21.312   , &
   -21.284   , -21.290   , -21.322   , -21.345   , &
   -21.354   , -21.366   , -21.385   , -21.396   , &
   -21.414   , -21.483   , -21.600   , -21.696   , &
   -21.742   , -21.759   , -21.776   , -21.816   , &
   -21.885   , -21.939   , -21.946   , -21.918   , &
   -21.873   , -21.818   , -21.756   , -21.689   , &
   -21.618   , -21.547   , -21.475   , -21.403   , &
   -21.331   , -21.260   , -21.188   , -21.114   , &
   -21.039   , -20.963   , -20.887   , -20.810   , &
   -20.734   , -20.657   , -20.581   , -20.505   , &
   -20.429   , -20.352   , -20.276   , -20.200   , &
   -20.125   , -20.049   , -19.973   , -19.898   , &
   -19.822   , -19.747   , -19.671   , -19.596   , &
   -19.520   , -19.445   , -19.370   , -19.295   , &
   -19.220   , -19.144   , -19.069   , -18.994   , &
   -18.919   , -18.869   , -18.819   , -18.769   , &
   -18.719   , -18.669   , -18.619   , -18.569   , &
   -18.519   , -18.469   , -18.419   /
   type rad_cooling_type
     real(kind=dp)     :: gamma
     real(kind=dp)     :: tfloor
     real(kind=dp)     :: cfrac
     logical           :: Tfix
     character(len=30) :: coolcurve
     character(len=30) :: coolmethod
   end type rad_cooling_type

   real(kind=dp),private   :: rc_gammamin1
  contains

    !> Read this module's parameters from a file
    subroutine rc_params_read(files)
      use mod_global_parameters, only: unitpar
      character(len=*), intent(in) :: files(:)
      integer                      :: i_file,i_reason
      character(len=70)            :: error_message

      namelist /rc_list/ coolcurve, coolmethod, ncool, cfrac, tlow, Tfix,&
                         coolsaveL,coolsavedT,coolsavemu,&
                         He_abundance,rc_chemical_gas_type,mean_mass,&
                         mean_mup,mean_ne_to_nH,mean_nall_to_nH
    error_message = 'At '//' mod_radiative_cooling.t'//'  in the procedure : rc_read_params'
    Loop_iparfile : do i_file = 1, size(files)
      open(unitpar, file=trim(files(i_file)), status="old")
      read(unitpar, rc_list, iostat=i_reason)
      cond_ierror : if(i_reason>0)then
       write(*,*)' Error in reading the parameters file : ',trim(files(i_file))
       write(*,*)' Error at namelist: ', 'rc_list'
       write(*,*)' Error number = ',i_reason
       write(*,*)' The code stops now '
       call mpistop(trim(error_message))
      elseif(i_reason<0)then cond_ierror
       write(*,*)' Reache the end of the file  : ',trim(files(i_file))
       write(*,*)' Error at namelist: rc_list'
       write(*,*)' Error number = ',i_reason
       write(*,*)' The code stops now '
       call mpistop(trim(error_message))
      else cond_ierror
       write(*,*)' End of reading of the rc_list'
      end if cond_ierror
      close(unitpar)

    end do Loop_iparfile
    if(tlow>0.9*bigdouble)tlow=-1.0_dp

      ! reset the setting according to user configuration
      coolconfig%He_abundance = He_abundance
      coolconfig%npoint = ncool
      coolconfig%curve  = coolcurve
      coolconfig%method = coolmethod
      coolconfig%cfrac  = cfrac
      coolconfig%tlow   = tlow
      coolconfig%Tfix   = Tfix
      coolconfig%saveL  = coolsaveL
      coolconfig%savedT = coolsavedT
      coolconfig%savemu = coolsavemu
      coolconfig%mean_mass     = mean_mass
      coolconfig%mean_mup     = mean_mup
      coolconfig%mean_ne_to_nH     = mean_ne_to_nH
      coolconfig%mean_nall_to_nH     = mean_nall_to_nH
      coolconfig%He_abundance     = He_abundance
      coolconfig%rc_chemical_gas_type     = rc_chemical_gas_type



    end subroutine rc_params_read

    !> Radiative cooling set default parameters
    subroutine radiative_cooling_config_default(phys_config_inuse)
      type(physconfig)                :: phys_config_inuse
      !-------------------------------------
      rc_gamma=phys_config_inuse%gamma
      He_abundance=phys_config_inuse%He_abundance

      ! old setting to be deleted
      ncool=4000
      coolcurve='JCcorona'
      coolmethod='exact'
      cfrac=0.1d0
      tlow       = -1.0_dp
      Tfix       = .false.
      coolsaveL  = .false.
      coolsavedT = .false.
      coolsavemu  = .false.
      !mean_ne_to_nH = 1.0d-3
      mean_mass   = (1.0_dp+4.0_dp*He_abundance)
      mean_mup     = mean_mass/(2.0_dp+2.0_dp*He_abundance)
      mean_ne_to_nH     = (1.0_dp+He_abundance)
      mean_nall_to_nH     = (1.0_dp+2.0_dp*He_abundance)
      rc_chemical_gas_type = 'fullyionised'


      rc_gammamin1 = rc_gamma -1.0_dp

      ! new setting
      coolconfig%He_abundance=He_abundance
      coolconfig%npoint = 4000
      coolconfig%curve  = 'JCcorona'
      coolconfig%method = 'exact'
      coolconfig%cfrac  = 0.1_dp
      coolconfig%tlow   = -1.0_dp
      coolconfig%Tfix       = .false.
      coolconfig%saveL  = .false.
      coolconfig%savedT = .false.
      coolconfig%savemu  = .false.
      !coolconfig%mean_ne_to_nH     = 1.0d-3
      coolconfig%mean_mass   = (1.0_dp+4.0_dp*coolconfig%He_abundance)
      coolconfig%mean_mup     = coolconfig%mean_mass/(2.0_dp+2.0_dp*coolconfig%He_abundance)
      coolconfig%mean_ne_to_nH     = (1.0_dp+coolconfig%He_abundance)
      coolconfig%mean_nall_to_nH     = (1.0_dp+2.0_dp*coolconfig%He_abundance)
      coolconfig%rc_chemical_gas_type = 'fullyionised'
    end subroutine radiative_cooling_config_default


    !> Radiative cooling initialization
  subroutine radiative_cooling_init(phys_indices_inuse,phys_config_inuse,&
                                      phys_gamma_insuse,phys_He_abund_insuse)
    use mod_global_parameters

    type(phys_variables_indices)    :: phys_indices_inuse
    type(physconfig)                :: phys_config_inuse
    real(kind=dp)   , intent(in)    :: phys_gamma_insuse,phys_He_abund_insuse
      ! .. local ..
    real(kind=dp)    :: tstep, Lstep
    integer ::  i,j
    !---------------------------------------------------------




    call radiative_cooling_config_default(phys_config_inuse)
    call rc_params_read(par_files)
    phys_config_inuse%cool_saveL   = coolconfig%saveL
    phys_config_inuse%cool_savedT  = coolconfig%savedT

    ! set the cooling variables indices
    call cooling_fill_phys_indices(phys_indices_inuse)
    !call rad_cooling_physical_units()
    ! Scale the lower temperature imposed by the user
    if (tlow>0.0_dp)tlow = tlow / unit_temperature


    cond_analytical : if(trim(coolcurve)/='analytical')then
      ! fill Lambda_HD and associate temperature array
      call cooling_fill_curve()


      ! Scale both T and Lambda
      tcool(1:ncool) = tcool(1:ncool) / unit_temperature
      Lcool(1:ncool) = Lcool(1:ncool) / unit_luminosity
    !(1.d0+2.d0*He_abundance)
      tcoolmin       = tcool(1)+smalldouble  ! avoid pointless interpolation
      tcoolmax       = tcool(ncool)
      ! smaller value for lowest temperatures from cooling table and user's choice
      if (tlow<0.0_dp)tlow = tcoolmin

      lgtcoolmin = dlog10(tcoolmin)
      lgtcoolmax = dlog10(tcoolmax)
      lgstep = (lgtcoolmax-lgtcoolmin)/ dble(ncool-1)

      dLdtcool(1)     = (Lcool(2)-Lcool(1))/(tcool(2)-tcool(1))
      dLdtcool(ncool) = (Lcool(ncool)-Lcool(ncool-1))/(tcool(ncool)-tcool(ncool-1))

      !do i=2,ncool-1
      !  dLdtcool(i) = (Lcool(i+1)-Lcool(i-1))/(tcool(i+1)-tcool(i-1))
      !enddo
      dLdtcool(2:ncool-1) = (Lcool(3:ncool)-Lcool(1:ncool-2))/(tcool(3:ncool)-tcool(1:ncool-2))


      cond_exact: if( coolmethod == 'exact' ) then
        !> this method is set according to van marle & Keppens 2010; equation 17
            Tref = tcoolmax
            Lref = Lcool(ncool)
            Yc(ncool) = zero
            Loop_i: do i=ncool-1, 1, -1
               Yc(i) = Yc(i+1)
               Loop_j: do j=1,100
                  tstep = 1.0d-2*(tcool(i+1)-tcool(i))
                  call findL(tcool(i+1)-j*tstep, 1.0_dp,Lstep)
                  Yc(i) = Yc(i) + Lref/Tref*tstep/Lstep
               enddo Loop_j
            enddo Loop_i

      endif cond_exact
    else cond_analytical
      call mpistop ('is not implimented for the moment')
      ! dummy to be added
    end if   cond_analytical

    phys_config%cool_tlow=tlow

  end subroutine radiative_cooling_init


  subroutine rc_fill_chemical_ionisation(He_abundance_sub,rc_chemical_gas_type,mean_nall_to_nH &
                                        ,mean_mass,mean_mup,mean_ne_to_nH)
    use mod_global_parameters
    implicit none
    real(kind=dp), intent(in)    :: He_abundance_sub
    character(len=*), intent(in) :: rc_chemical_gas_type
    real(kind=dp), intent(out)   :: mean_nall_to_nH,mean_mass,mean_mup,mean_ne_to_nH
    !----------------------------------------------------------
    if(He_abundance_sub>0.0_dp)then
      mean_nall_to_nH = (1.0_dp+2.0_dp*He_abundance_sub)
      mean_mass       = (1.0_dp+4.0_dp*He_abundance_sub)
      select case(trim(rc_chemical_gas_type))
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
         write(*,*) 'The chemical gas type : ', trim(rc_chemical_gas_type)
          write (*,*) "Undefined gas chemical type entered in mod_hd_phys.t "
          call mpistop('The stops at rc_fill_chemical_ionisation in src/hd/mod_hd_phys.t')
      end select
      !hd_config%mean_mup          = (2.0_dp+3.0_dp*He_abundance_sub)

    else if(dabs(He_abundance_sub)<=smalldouble)then
      mean_mass         = 1.0_dp
      mean_mup          = 1.0_dp
      mean_ne_to_nH     = 1.0_dp
    end if
  end   subroutine rc_fill_chemical_ionisation

  subroutine rad_cooling_physical_units

    use mod_global_parameters
    use mod_physics
    implicit none
    !.......................
     unit_luminosity=  unit_pressure/( unit_numberdensity**2.0_dp * unit_time &
                                 *  phys_config%mean_mass**2.0_dp)



     if(coolsavedT)w_convert_factor(dtcool1_) = unit_time
     if(coolsaveL)w_convert_factor(Lcool1_)   = unit_luminosity

  end subroutine rad_cooling_physical_units

  !> set the cooling variables indices
  subroutine cooling_fill_phys_indices(phys_indices_inuse)
    use mod_global_parameters
    implicit none
    type(phys_variables_indices)    :: phys_indices_inuse
    integer                         :: idir
    ! Determine flux variables
    rho_= phys_indices_inuse%rho_
    allocate(mom(ndir))
    Loop_idir : do idir = 1, ndir
       mom(idir) = phys_indices_inuse%mom(idir)!nwx       ! momentum density
    end do Loop_idir
    e_     = phys_indices_inuse%e_!nwx          ! energy density


    if(coolconfig%saveL)then
      Lcool1_   = var_set_extravar('L', 'L', 1)
      phys_indices_inuse%Lcool1_   = Lcool1_
    end if
    if(coolconfig%savedT)then
      dtcool1_ = var_set_extravar('dtcool', 'dtcool', 1)
      phys_indices_inuse%dtcool1_ = dtcool1_
    end if

  end   subroutine cooling_fill_phys_indices

  !> fill the cooling curve Lambda_HD and T (Temperature)
  subroutine cooling_fill_curve
   use mod_global_parameters
   implicit none
   ! .. local ..
   real(kind=dp)   , dimension(:), allocatable :: t_table
   real(kind=dp)   , dimension(:), allocatable :: L_table
   real(kind=dp)    :: ratt, Lerror
   real(kind=dp)    ::fact1, fact2, fact3, dL1, dL2
   real(kind=dp)    :: tstep, Lstep
   integer :: ntable, i, j, ic, idir
   logical :: jump
   !------------------------------------------------------------

      allocate(tcool(1:ncool), Lcool(1:ncool), dLdtcool(1:ncool))
      allocate(Yc(1:ncool), invYc(1:ncool))

      tcool(1:ncool)    = zero
      Lcool(1:ncool)    = zero
      dLdtcool(1:ncool) = zero

      ! Read in the selected cooling curve
      select case(coolcurve)

      case('JCcorona')
         if(mype ==0) &
         print *,'Use Colgan & Feldman (2008) cooling curve'
         if(mype ==0) &
         print *,'This version only till 10000 K, beware for floor T treatment'

         ntable = n_JCcorona

         allocate(t_table(1:ntable))
         allocate(L_table(1:ntable))
         t_table(1:ntable) = t_JCcorona(1:n_JCcorona)
         L_table(1:ntable) = l_JCcorona(1:n_JCcorona)

      case('DM')
         if(mype ==0) &
         print *,'Use Delgano & McCray (1972) cooling curve'

         ntable = n_DM

         allocate(t_table(1:ntable))
         allocate(L_table(1:ntable))
         t_table(1:ntable) = t_DM(1:n_DM)
         L_table(1:ntable) = l_DM(1:n_DM)

      case('MB')
         if(mype ==0) &
         write(*,'(3a)') 'Use MacDonald & Bailey (1981) cooling curve '&
              ,'as implemented in ZEUS-3D for log(T)>=4.0 K, with the values '&
              ,'from Dalgarno & McCRay (1972) for low temperatures as'&
              ,'tabulated from Schure et al. (2006) at f_i=1e-3 between'&
              ,'log(T)=1.0 and 4.0 K (4.0 K excluded).'

         !ntable = n_MB + 20
         ntable = n_DM_2 + n_MB - 1

         allocate(t_table(1:ntable))
         allocate(L_table(1:ntable))

         !t_table(1:ntable) = t_DM(1:21)
         !L_table(1:ntable) = l_DM(1:21)
         !t_table(22:ntable) = t_MB(2:n_MB)
         !L_table(22:ntable) = l_MB(2:n_MB)
         !
         ! To be used together with the SPEX table for the SPEX_DM option
         ! Assuming an ionization fraction of 10^-3 with log(T_break) = 4.0
         !
         t_table(1:ntable) = t_DM_2(1:n_DM_2-1)
         L_table(1:ntable) = l_DM_2(1:n_DM_2-1)
         t_table(n_DM_2:ntable) = t_MB(1:n_MB)
         L_table(n_DM_2:ntable) = l_MB(1:n_MB)

         case('MB_original')
            if(mype ==0) &
            write(*,'(3a)') 'Use MacDonald & Bailey (1981) cooling curve '&
                 ,'as implemented in ZEUS-3D, with the values '&
                 ,'from Delgano & McCRay (1972) for low temperatures.'

            ntable = n_MB + 20


            allocate(t_table(1:ntable))
            allocate(L_table(1:ntable))

            t_table(1:ntable) = t_DM(1:21)
            L_table(1:ntable) = l_DM(1:21)
            t_table(22:ntable) = t_MB(2:n_MB)
            L_table(22:ntable) = l_MB(2:n_MB)

      case('MLcosmol')
         if(mype ==0) &
         print *,'Use Mellema & Lundqvist (2002) cooling curve '&
              ,'for zero metallicity '

         ntable = n_MLcosmol

         allocate(t_table(1:ntable))
         allocate(L_table(1:ntable))
         t_table(1:ntable) = t_MLcosmol(1:n_MLcosmol)
         L_table(1:ntable) = l_MLcosmol(1:n_MLcosmol)

      case('MLwc')
         if(mype ==0) &
         print *,'Use Mellema & Lundqvist (2002) cooling curve '&
              ,'for WC-star metallicity '

         ntable = n_MLwc

         allocate(t_table(1:ntable))
         allocate(L_table(1:ntable))
         t_table(1:ntable) = t_MLwc(1:n_MLwc)
         L_table(1:ntable) = l_MLwc(1:n_MLwc)

      case('MLsolar1')
         if(mype ==0) &
         print *,'Use Mellema & Lundqvist (2002) cooling curve '&
              ,'for solar metallicity '

         ntable = n_MLsolar1

         allocate(t_table(1:ntable))
         allocate(L_table(1:ntable))
         t_table(1:ntable) = t_MLsolar1(1:n_MLsolar1)
         L_table(1:ntable) = l_MLsolar1(1:n_MLsolar1)

      case('cloudy_ism')
         if(mype ==0) &
         print *,'Use Cloudy based cooling curve '&
              ,'for ism metallicity '

         ntable = n_cl_ism

         allocate(t_table(1:ntable))
         allocate(L_table(1:ntable))
         t_table(1:ntable) = t_cl_ism(1:n_cl_ism)
         L_table(1:ntable) = l_cl_ism(1:n_cl_ism)

      case('cloudy_solar')
         if(mype ==0) &
         print *,'Use Cloudy based cooling curve '&
              ,'for solar metallicity '

         ntable = n_cl_solar

         allocate(t_table(1:ntable))
         allocate(L_table(1:ntable))
         t_table(1:ntable) = t_cl_solar(1:n_cl_solar)
         L_table(1:ntable) = l_cl_solar(1:n_cl_solar)

      case('SPEX')
         if(mype ==0) &
         print *,'Use SPEX cooling curve (Schure et al. 2009) '&
              ,'for solar metallicity '

         ntable = n_SPEX

         allocate(t_table(1:ntable))
         allocate(L_table(1:ntable))
         t_table(1:ntable) = t_SPEX(1:n_SPEX)
         L_table(1:ntable) = l_SPEX(1:n_SPEX) + log10(nenh_SPEX(1:n_SPEX))

      case('SPEX_DM')
         if(mype ==0) then
            print *, 'Use SPEX cooling curve for solar metallicity above 10^4 K. '
            print *, 'At lower temperatures,use Dalgarno & McCray (1972), '
            print *, 'with a pre-set ionization fraction of 10^-3. '
            print *, 'as described by Schure et al. (2009). '
         endif

         ntable = n_SPEX + n_DM_2 - 6

         allocate(t_table(1:ntable))
         allocate(L_table(1:ntable))
         t_table(1:n_DM_2-1) = t_DM_2(1:n_DM_2-1)
         L_table(1:n_DM_2-1) = L_DM_2(1:n_DM_2-1)
         t_table(n_DM_2:ntable) = t_SPEX(6:n_SPEX)
         L_table(n_DM_2:ntable) = l_SPEX(6:n_SPEX) + log10(nenh_SPEX(6:n_SPEX))

      !  if(phys_config%chemical_on)L_table(1:n_DM_2-1)=L_table(1:n_DM_2-1)/coolconfig%fn
      case default
         call mpistop("This coolingcurve is unknown")
      end select

      ! create cooling table(s) for use in amrvac
      tcoolmax = t_table(ntable)
      tcoolmin = t_table(1)
      ratt     = (tcoolmax-tcoolmin)/( dble(ncool-1) + smalldouble)

      tcool(1) = tcoolmin
      Lcool(1) = L_table(1)

      tcool(ncool) = tcoolmax
      Lcool(ncool) = L_table(ntable)

      do i=2,ncool        ! loop to create one table
        tcool(i) = tcool(i-1)+ratt
        do j=1,ntable-1   ! loop to create one spot on a table
        ! Second order polynomial interpolation, except at the outer edge,
        ! or in case of a large jump.
          if(tcool(i) < t_table(j+1)) then
             if(j.eq. ntable-1 )then
               fact1 = (tcool(i)-t_table(j+1))     &
                     /(t_table(j)-t_table(j+1))

               fact2 = (tcool(i)-t_table(j))       &
                     /(t_table(j+1)-t_table(j))

               Lcool(i) = L_table(j)*fact1 + L_table(j+1)*fact2
               exit
             else
               dL1 = L_table(j+1)-L_table(j)
               dL2 = L_table(j+2)-L_table(j+1)
               jump =(max(dabs(dL1),dabs(dL2)) > 2*min(dabs(dL1),dabs(dL2)))
             endif

             if( jump ) then
               fact1 = (tcool(i)-t_table(j+1))     &
                     /(t_table(j)-t_table(j+1))

               fact2 = (tcool(i)-t_table(j))       &
                     /(t_table(j+1)-t_table(j))

               Lcool(i) = L_table(j)*fact1 + L_table(j+1)*fact2
               exit
             else
               fact1 = ((tcool(i)-t_table(j+1))     &
                     * (tcool(i)-t_table(j+2)))   &
                     / ((t_table(j)-t_table(j+1)) &
                     * (t_table(j)-t_table(j+2)))

               fact2 = ((tcool(i)-t_table(j))       &
                     * (tcool(i)-t_table(j+2)))   &
                     / ((t_table(j+1)-t_table(j)) &
                     * (t_table(j+1)-t_table(j+2)))

               fact3 = ((tcool(i)-t_table(j))       &
                     * (tcool(i)-t_table(j+1)))   &
                     / ((t_table(j+2)-t_table(j)) &
                     * (t_table(j+2)-t_table(j+1)))

               Lcool(i) = L_table(j)*fact1 + L_table(j+1)*fact2 &
                        + L_table(j+2)*fact3
               exit
             endif
          endif
        enddo  ! end loop to find create one spot on a table
      enddo    ! end loop to create one table

      ! Go from logarithmic to actual values.
      tcool(1:ncool) = 10.0_dp**tcool(1:ncool)
      Lcool(1:ncool) = 10.0_dp**Lcool(1:ncool)


      deallocate(t_table)
      deallocate(L_table)
    end subroutine cooling_fill_curve



    subroutine cooling_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
      use mod_global_parameters

      integer, intent(in)             :: ixI^L, ixO^L
      real(kind=dp)   , intent(in)    :: dx^D, x(ixI^S,1:ndim)
      real(kind=dp)   , intent(in)    :: w(ixI^S,1:nw)
      real(kind=dp)   , intent(inout) :: dtnew


      real(kind=dp)                   :: etherm(ixI^S),ndensity_electron(ixO^S)
      real(kind=dp)                   :: L1,Tlocal1, ptherm(ixI^S), lum(ixI^S)
      real(kind=dp)                   :: plocal, rholocal
      real(kind=dp)                   :: Lmax
      integer                         :: ix^D
      !
      ! Limit timestep to avoid cooling problems when using explicit cooling
      !
      if(phys_config%chemical_on) then
        call chemical_get_electron_density(ixI^L, ixO^L,w,ndensity_electron)
      else
        ndensity_electron(ixO^S) =   w(ixO^S,rho_)
      end if
      if(coolmethod == 'explicit1') then
       call phys_get_pthermal(w,x,ixI^L,ixO^L,ptherm)
       {do ix^DB = ixO^LIM^DB\}
         plocal   = ptherm(ix^D)
         rholocal = w(ix^D,rho_)
         !  Tlocal = P/rho
         Tlocal1       = max(plocal/(rholocal),smalldouble)
         !  Determine explicit cooling
         !  If temperature is below floor level, no cooling.
         !  Stop wasting time and go to next gridpoint.
         !  If the temperature is higher than the maximum,
         !  assume Bremmstrahlung
         if( Tlocal1<=max(tcoolmin,tlow) ) then
            L1 = zero
         else if( Tlocal1>=tcoolmax )then
            L1         = Lcool(ncool)*sqrt(Tlocal1/tcoolmax)
            L1         = L1*(rholocal*ndensity_electron(ix^D))
         else
            call findL(Tlocal1,ndensity_electron(ix^D),L1)
            L1         = L1*(rholocal*ndensity_electron(ix^D))
         endif
         lum(ix^D) = L1
        {enddo^D&\}
        etherm(ixO^S)=ptherm(ixO^S)/(rc_gammamin1)
        !if(coolsavedT)w(ixO^S,dtcool1_) = cfrac*(etherm(ixO^S)/max(lum(ixO^S),smalldouble))
        dtnew = cfrac*minval(etherm(ixO^S)/max(lum(ixO^S),smalldouble))
      endif

    end subroutine cooling_get_dt

    subroutine getvar_cooling(ixI^L,ixO^L,w,x,coolrate)
    !
    ! Create extra variable to show cooling rate in the output
    ! Uses a simple explicit scheme.
    ! N.B. Since there is no knowledge of the timestep size,
    ! there is no upper limit for the cooling rate.
      use mod_global_parameters

      integer, intent(in)          :: ixI^L,ixO^L
      real(kind=dp)   , intent(in) :: x(ixI^S,1:ndim)
      real(kind=dp)                :: w(ixI^S,1:nw)
      real(kind=dp)   , intent(out):: coolrate(ixI^S)

      real(kind=dp)    :: etherm(ixI^S),ndensity_electron(ixO^S)
      real(kind=dp)    :: L1,Tlocal1, ptherm(ixI^S)
      real(kind=dp)    :: plocal, rholocal
      real(kind=dp)    :: emin

      integer :: ix^D

      call phys_get_pthermal(w,x,ixI^L,ixO^L,ptherm)

      if(phys_config%chemical_on) then
        call chemical_get_electron_density(ixI^L, ixO^L,w,ndensity_electron)
      else
        ndensity_electron(ixO^S) =   w(ixO^S,rho_)
      end if

      {do ix^DB = ixO^LIM^DB\}
         plocal   = ptherm(ix^D)
         rholocal = w(ix^D,rho_)
         !  Tlocal = P/rho
         Tlocal1       = max(plocal/(rholocal),smalldouble)
         !
         !  Determine explicit cooling
         !
         if( Tlocal1<=max(tcoolmin,tlow) ) then
            L1         = zero
         else if( Tlocal1>=tcoolmax )then
            L1         = Lcool(ncool)*sqrt(Tlocal1/tcoolmax)
            L1         = L1*(rholocal*ndensity_electron(ix^D))
         else
            call findL(Tlocal1,ndensity_electron(ix^D),L1)
            L1         = L1*(rholocal*ndensity_electron(ix^D))
         endif
         coolrate(ix^D)= L1
      {enddo^D&\}

    end subroutine getvar_cooling

    subroutine radiative_cooling_add_source(qdt,ixI^L,ixO^L,wCT,w,x,&
         qsourcesplit,active)

    ! w[iw]=w[iw]+qdt*S[wCT,x] where S is the source based on wCT within ixO
      use mod_global_parameters

      integer, intent(in) :: ixI^L, ixO^L
      real(kind=dp)   , intent(in) :: qdt, x(ixI^S,1:ndim), wCT(ixI^S,1:nw)
      real(kind=dp)   , intent(inout) :: w(ixI^S,1:nw)
      logical, intent(in) :: qsourcesplit
      logical, intent(inout) :: active

      call rc_fill_chemical_ionisation(coolconfig%He_abundance,&
                    coolconfig%rc_chemical_gas_type,coolconfig%mean_nall_to_nH, &
                                            coolconfig%mean_mass,coolconfig%mean_mup,&
                                            coolconfig%mean_ne_to_nH)
      !write(*,*) coolconfig%He_abundance,&
      !              coolconfig%rc_chemical_gas_type,coolconfig%mean_nall_to_nH, &
      !                                      coolconfig%mean_mass,coolconfig%mean_mup,&
      !                                      coolconfig%mean_ne_to_nH

      if(qsourcesplit .eqv. rc_split) then
        active = .true.
        select case(coolmethod)
        case ('explicit1')
          if(mype==0)then
            if(it==1) then
              write(*,*)'Fully explicit cooling is not completely safe in this version'
              write(*,*)'PROCEED WITH CAUTION!'
            endif
          endif
          call cool_explicit1(qdt,ixI^L,ixO^L,wCT,w,x)
        case ('explicit2')
          call cool_explicit2(qdt,ixI^L,ixO^L,wCT,w,x)
        case ('semiimplicit')
          call cool_semiimplicit(qdt,ixI^L,ixO^L,wCT,w,x)
        case ('implicit')
          call cool_implicit(qdt,ixI^L,ixO^L,wCT,w,x)
        case ('exact')
          call cool_exact(qdt,ixI^L,ixO^L,wCT,w,x)
        case default
          call mpistop("This cooling method is unknown")
        end select
        if( Tfix ) call floortemperature(qdt,ixI^L,ixO^L,tlow,wCT,w,x)
      end if

    end subroutine radiative_cooling_add_source

    subroutine floortemperature(qdt,ixI^L,ixO^L,tfloor,wCT,w,x)
    !  Force minimum temperature to a fixed temperature
      use mod_global_parameters

      integer, intent(in)             :: ixI^L, ixO^L
      real(kind=dp)   , intent(in)    :: tfloor
      real(kind=dp)   , intent(in)    :: qdt, x(ixI^S,1:ndim), wCT(ixI^S,1:nw)
      real(kind=dp)   , intent(inout) :: w(ixI^S,1:nw)

      real(kind=dp)    :: etherm_emin(ixI^S)
      !--------------------------------------------
      !tfloor = tlow
      ! use 'etherm' to save the pressure
      call phys_get_pthermal(w,x,ixI^L,ixO^L,etherm_emin)
      ! compute the thermal energy min the minimum thermal energy
      etherm_emin(ixO^S) = etherm_emin(ixO^S)/(rc_gammamin1)-( w(ixO^S,rho_)*tfloor/(rc_gammamin1))
      where(etherm_emin(ixO^S)<0.0_dp )
          w(ixO^S,e_)=w(ixO^S,e_)-etherm_emin(ixO^S)
      end where

    end subroutine floortemperature

    subroutine cool_explicit1(qdt,ixI^L,ixO^L,wCT,w,x)
    ! explicit cooling routine that depends on getdt to
    ! adjust the timestep. Accurate but incredibly slow
      use mod_global_parameters

      integer, intent(in)             :: ixI^L, ixO^L
      real(kind=dp)   , intent(in)    :: qdt, x(ixI^S,1:ndim)
      real(kind=dp)   , intent(inout) :: w(ixI^S,1:nw)
      real(kind=dp)                   :: wCT(ixI^S,1:nw)

      real(kind=dp)    :: L_hd,Tlocal1
      real(kind=dp), dimension(ixI^S)    ::  pnew,temperature!,ndensity_electron
      real(kind=dp)    ::  nlocal !rholocal
      real(kind=dp)    :: emin, Lmax

      integer :: ix^D
      integer :: icool
      !-----------------------------------------------

      call phys_get_pthermal(w,x,ixI^L,ixO^L,pnew)
      call phys_get_temperature( ixI^L, ixO^L,w, x, temperature)
      ! if(phys_config%chemical_on) then
      !   call chemical_get_electron_density(ixI^L, ixO^L,wCT,ndensity_electron)
      ! else
      !   ndensity_electron(ixO^S) =   wCT(ixO^S,rho_)
      ! end if
      {do ix^DB = ixO^LIM^DB\}
         !plocal   = ptherm(ix^D)
         nlocal    = wCT(ix^D,rho_)
         emin     = nlocal*tlow/rc_gammamin1
         Lmax     = max(zero,pnew(ix^D)/rc_gammamin1-emin)/qdt
         !  Tlocal = P/rho
         Tlocal1  = max(temperature(ix^D),smalldouble)!max(plocal/nlocal,smalldouble)
         !
         !  Determine explicit cooling
         !
         !  If temperature is below floor level, no cooling.
         !  Stop wasting time and go to next gridpoint.
         !  If the temperature is higher than the maximum,
         !  assume Bremmstrahlung
         if( Tlocal1<=max(tcoolmin,tlow) ) then
            L_hd = zero
         else if( Tlocal1>=tcoolmax )then
            L_hd         = Lcool(ncool)*sqrt(Tlocal1/tcoolmax)
            L_hd         = L_hd*(nlocal*nlocal)
            L_hd         = min(L_hd,Lmax)
            w(ix^D,e_) = w(ix^D,e_)-L_hd*qdt
         else
            call findL(Tlocal1,nlocal,L_hd)
            L_hd         = L_hd*(nlocal*nlocal)
            L_hd         = min(L_hd,Lmax)
            w(ix^D,e_) = w(ix^D,e_)-L_hd*qdt
         endif
         if(coolsavedT)then
           w(ix^D,dtcool1_) = cfrac*(pnew(ix^D)/rc_gammamin1/max(L_hd,smalldouble))
         end if
         if(coolsaveL)w(ix^D,Lcool1_)   = L_hd
      {enddo^D&\}

    end subroutine cool_explicit1

    subroutine cool_explicit2(qdt,ixI^L,ixO^L,wCT,w,x)
    !
    ! explicit cooling routine that does a series
    ! of small forward integration steps, to make
    ! sure the amount of cooling remains correct
    !
    ! Not as accurate as 'explicit1', but a lot faster
    ! tends to overestimate cooling

      use mod_global_parameters

      integer, intent(in)             :: ixI^L, ixO^L
      real(kind=dp)   , intent(in)    :: qdt, x(ixI^S,1:ndim)
      real(kind=dp)   , intent(inout) :: w(ixI^S,1:nw)
      real(kind=dp)                   :: wCT(ixI^S,1:nw)
      real(kind=dp)    :: Ltest, etherm, de
      real(kind=dp)    :: dtmax, dtstep

      integer :: idt,ndtstep

      real(kind=dp)    :: L1,Tlocal1,ptherm(ixI^S),pnew(ixI^S),ndensity_electron(ixO^S)
      real(kind=dp)    :: plocal, rholocal
      real(kind=dp)    :: emin, Lmax

      integer :: ix^D
      integer :: icool

      call phys_get_pthermal(wCT,x,ixI^L,ixO^L,ptherm)
      call phys_get_pthermal(w,x,ixI^L,ixO^L,pnew )
      if(phys_config%chemical_on) then
        call chemical_get_electron_density(ixI^L, ixO^L,wCT,ndensity_electron)
      else
        ndensity_electron(ixO^S) =   wCT(ixO^S,rho_)
      end if

      {do ix^DB = ixO^LIM^DB\}
         !  Calculate explicit cooling value
         dtmax    = qdt
         plocal   = ptherm(ix^D)
         etherm   = plocal/(rc_gammamin1)

         rholocal = wCT(ix^D,rho_)
         emin     = rholocal*tlow/(rc_gammamin1)
         Lmax            = max(zero,(pnew(ix^D)/(rc_gammamin1))-emin)/qdt
         !  Tlocal = P/rho
         Tlocal1       = plocal/(rholocal)
         !
         !  Determine explicit cooling
         !
         !  If temperature is below floor level, no cooling.
         !  Stop wasting time and go to next gridpoint.
         !  If the temperature is higher than the maximum,
         !  assume Bremmstrahlung
         if( Tlocal1<=max(tcoolmin,tlow) ) then
            Ltest = zero
         else if( Tlocal1>=tcoolmax )then
            Ltest = Lcool(ncool)*sqrt(Tlocal1/tcoolmax)
            Ltest = L1*(rholocal*ndensity_electron(ix^D))
            Ltest = min(L1,Lmax)
            if( dtmax>cfrac*etherm/Ltest) dtmax = cfrac*etherm/Ltest
         else
            call findL(Tlocal1,ndensity_electron(ix^D),Ltest)
            Ltest = Ltest*(rholocal*ndensity_electron(ix^D))
            Ltest = min(Ltest,Lmax)
            if( dtmax>cfrac*etherm/Ltest) dtmax = cfrac*etherm/Ltest
         endif
         !
         !  Calculate number of steps for cooling
         !
         ndtstep = max(nint(qdt/dtmax),1)+1
         dtstep = qdt/ndtstep
         !
         !  Use explicit cooling value for first step
         !
         de     = Ltest*dtstep
         etherm = etherm - de

         do idt=2,ndtstep
            plocal = etherm*(rc_gammamin1)
            Lmax   = max(zero,etherm-emin)/dtstep
            !  Tlocal = P/rho
            Tlocal1 = plocal/(rholocal)
            if( Tlocal1<=max(tcoolmin,tlow) ) then
               L1 = zero
               exit
            else if( Tlocal1>=tcoolmax )then
               L1 = Lcool(ncool)*sqrt(Tlocal1/tcoolmax)
               L1 = L1*(rholocal*ndensity_electron(ix^D))
               L1 = min(L1,Lmax)
            else
               call findL(Tlocal1,ndensity_electron(ix^D),L1)
               L1 = L1*(rholocal*ndensity_electron(ix^D))
               L1 = min(L1,Lmax)
            endif

            de     = de + L1*dtstep
            etherm = etherm - L1*dtstep
         enddo
         w(ix^D,e_) = w(ix^D,e_) -de
         if(coolsavedT)w(ix^D,dtcool1_) = cfrac*(plocal/rc_gammamin1/max(L1,smalldouble))
         if(coolsaveL)w(ix^D,Lcool1_) = L1
      {enddo^D&\}

    end subroutine cool_explicit2

    subroutine cool_semiimplicit(qdt,ixI^L,ixO^L,wCT,w,x)
    !
    ! Semi-implicit cooling method based on a two point
    ! average
    !
    ! Fast, but tends to underestimate cooling
      use mod_global_parameters

      integer, intent(in)             :: ixI^L, ixO^L
      real(kind=dp)   , intent(in)    :: qdt, x(ixI^S,1:ndim)
      real(kind=dp)   , intent(inout) :: w(ixI^S,1:nw)
      real(kind=dp)                   ::  wCT(ixI^S,1:nw)
      real(kind=dp)    :: L1,L2,Tlocal1, Tlocal2
      real(kind=dp)    :: etemp

      real(kind=dp)    :: plocal, rholocal
      real(kind=dp)    :: emin, Lmax
      real(kind=dp)    :: mp,kB,miu0

      real(kind=dp)    :: ptherm(ixI^S),pnew(ixI^S),ndensity_electron(ixO^S)
      integer :: ix^D



      if(SI_unit) then
        mp=mp_SI
        kB=kB_SI
        miu0=miu0_SI
      else
        mp=mp_cgs
        kB=kB_cgs
        miu0=4.0_dp*dpi
      end if

      call phys_get_pthermal(wCT,x,ixI^L,ixO^L,ptherm)
      call phys_get_pthermal(w,x,ixI^L,ixO^L,pnew )
      if(phys_config%chemical_on) then
        call chemical_get_electron_density(ixI^L, ixO^L,wCT,ndensity_electron)
      else
        ndensity_electron(ixO^S) =   wCT(ixO^S,rho_)
      end if
      {do ix^DB = ixO^LIM^DB\}
         plocal   = ptherm(ix^D)
         rholocal = wCT(ix^D,rho_)
         emin     = kB*rholocal*tlow/(rc_gammamin1)
         Lmax            = max(zero,pnew(ix^D)/(rc_gammamin1)-emin)/qdt
         !  Tlocal = P/rho
         Tlocal1       = max(plocal/(kB*rholocal),smalldouble)
         !
         !  Determine explicit cooling at present temperature
         !
         !  If temperature is below floor level, no cooling.
         !  Stop wasting time and go to next gridpoint.
         !  If the temperature is higher than the maximum,
         !  assume Bremmstrahlung
         if( Tlocal1<=max(tcoolmin,tlow) ) then
            L1 = zero
            L2 = zero
         else
           if( Tlocal1>=tcoolmax )then
              L1 = Lcool(ncool)*sqrt(Tlocal1/tcoolmax)
           else
              call findL(Tlocal1,ndensity_electron(ix^D),L1)
           end if
           L1      = L1*(rholocal*ndensity_electron(ix^D))
           etemp   = plocal/(rc_gammamin1) - L1*qdt
           Tlocal2 = etemp*(rc_gammamin1)/(kB*rholocal)
           !
           !  Determine explicit cooling at new temperature
           !
           if( Tlocal2<=tcoolmin ) then
              L2 = zero
           else if( Tlocal2>=tcoolmax )then
              L2 = Lcool(ncool)*sqrt(Tlocal2/tcoolmax)
           else
              call findL(Tlocal2,ndensity_electron(ix^D),L2)
           end if
           L2  = L2*(rholocal*ndensity_electron(ix^D))
           w(ix^D,e_) = w(ix^D,e_) - min(half*(L1+L2),Lmax)*qdt
         endif
         if(coolsavedT)w(ix^D,dtcool1_) = cfrac*(plocal/rc_gammamin1/max(L1,smalldouble))
         if(coolsaveL)w(ix^D,Lcool1_) = L1
      {enddo^D&\}

    end subroutine cool_semiimplicit

    subroutine cool_implicit(qdt,ixI^L,ixO^L,wCT,w,x)
    !
    ! Implicit cooling method based on a half-step
    ! refinement algorithm
    !
      use mod_global_parameters

      integer, intent(in)             :: ixI^L, ixO^L
      real(kind=dp)   , intent(in)    :: qdt, x(ixI^S,1:ndim)
      real(kind=dp)   , intent(inout) :: w(ixI^S,1:nw)
      real(kind=dp)                   :: wCT(ixI^S,1:nw)
      real(kind=dp)    :: Ltemp,Tlocal1,Tnew,f1,f2
      real(kind=dp)    :: ptherm(ixI^S), pnew(ixI^S),ndensity_electron(ixO^S)

      real(kind=dp)    :: plocal, rholocal, elocal
      real(kind=dp)    :: emin, Lmax, eold, enew, estep
      real(kind=dp)    :: mp,kB,miu0
      integer, parameter :: maxiter = 100
      real(kind=dp)   , parameter :: e_error = 1.0D-6

      integer :: ix^D, j

      if(SI_unit) then
        mp=mp_SI
        kB=kB_SI
        miu0=miu0_SI
      else
        mp=mp_cgs
        kB=kB_cgs
        miu0=4.0_dp*dpi
      end if

      call phys_get_pthermal(wCT,x,ixI^L,ixO^L,ptherm)
      call phys_get_pthermal(w,x,ixI^L,ixO^L,pnew )
      if(phys_config%chemical_on) then
        call chemical_get_electron_density(ixI^L, ixO^L,wCT,ndensity_electron)
      else
        ndensity_electron(ixO^S) =   wCT(ixO^S,rho_)
      end if

      {do ix^DB = ixO^LIM^DB\}
         plocal   = ptherm(ix^D)
         elocal   = plocal/(rc_gammamin1)
         rholocal = wCT(ix^D,rho_)
         emin     = kB*rholocal*tlow/(rc_gammamin1)
         Lmax            = max(zero,pnew(ix^D)/(rc_gammamin1)-emin)/qdt
         !  Tlocal = P/rho
         Tlocal1       = max(plocal/(kB*rholocal),smalldouble)
         !
         !  Determine explicit cooling at present temperature
         !
         !  If temperature is below floor level, no cooling.
         !  Stop wasting time and go to next gridpoint.
         !  If the temperature is higher than the maximum,
         !  assume Bremmstrahlung
         if( Tlocal1<=max(tcoolmin,tlow) ) then
            Ltemp = zero
         else
           eold  = elocal
           enew  = elocal
           estep = -(smalldouble)
           f2    = 1.d0
           do j=1,maxiter+1
             if( j>maxiter ) call mpistop("Implicit cooling exceeds maximum iterations")
             Tnew  = enew*(rc_gammamin1)/(kB*rholocal)
             if( Tnew<=tcoolmin ) then
               Ltemp = Lmax
               exit
             else if( Tnew>=tcoolmax )then
               Ltemp = Lcool(ncool)*sqrt(Tnew/tcoolmax)
             else
               call findL(Tnew,ndensity_electron(ix^D),Ltemp)
             end if
             Ltemp = Ltemp*(rholocal*ndensity_electron(ix^D))
             eold  = enew + Ltemp*qdt

             f1 = elocal -eold
             if( abs(half*f1/(elocal+eold)) < e_error  ) exit

             if(j==1) estep = max((elocal-emin)*half,smalldouble)
             if(f1*f2 < zero ) estep = -half*estep
             f2 = f1
             enew = enew +estep
           enddo
         endif
         w(ix^D,e_)     = w(ix^D,e_) - min(Ltemp,Lmax)*qdt
         if(coolsavedT)w(ix^D,dtcool1_) = cfrac*(plocal/rc_gammamin1/max(min(Ltemp,Lmax),smalldouble))
         if(coolsaveL)w(ix^D,Lcool1_) = min(Ltemp,Lmax)
      {enddo^D&\}

    end subroutine cool_implicit

    subroutine cool_exact(qdt,ixI^L,ixO^L,wCT,w,x)
    !
    !  Cooling routine using exact integration method from Townsend 2009
    !
      use mod_global_parameters

      integer, intent(in)             :: ixI^L, ixO^L
      real(kind=dp)   , intent(in)    :: qdt, x(ixI^S,1:ndim)
      real(kind=dp)   , intent(inout) :: w(ixI^S,1:nw)
      real(kind=dp)   , intent(in)    :: wCT(ixI^S,1:nw)
      ! .. local ..
      real(kind=dp)    :: Y1, tc, Y2
      real(kind=dp)    :: L1,Tlocal1,  Tlocal2
      real(kind=dp), dimension(ixI^S) :: ptherm,pnew,temperature
      real(kind=dp)    :: nlocal,nHlocal
      real(kind=dp)    :: emin, Lmax, fact
      real(kind=dp)    :: mp,kB,miu0
      real(kind=dp)    :: ndensity_electron(ixO^S)

      integer :: ix^D
      integer :: icool
      !------------------------------------------


      if(SI_unit) then
        mp=mp_SI
        kB=kB_SI
        miu0=miu0_SI
      else
        mp=mp_cgs
        kB=kB_cgs
        miu0=4.0_dp*dpi
      end if

    !  call phys_get_pthermal(wCT,x,ixI^L,ixO^L,ptherm)
      call phys_get_pthermal(w,x,ixI^L,ixO^L,pnew )
      call phys_get_temperature( ixI^L, ixO^L,w, x, temperature)
      fact = Lref*qdt/Tref

      !write(*,*) coolconfig%He_abundance,&
      !              coolconfig%rc_chemical_gas_type,coolconfig%mean_nall_to_nH, &
      !                                      coolconfig%mean_mass,coolconfig%mean_mup,&
      !                                      coolconfig%mean_ne_to_nH


      !invgam=1.0_dp/(rc_gammamin1)
      ! if(phys_config%chemical_on) then
      !   call chemical_get_electron_density(ixI^L, ixO^L,wCT,ndensity_electron)
      ! else
      !   ndensity_electron(ixO^S) =   wCT(ixO^S,rho_)
      ! end if
      {do ix^DB = ixO^LIM^DB\}
        ! plocal   = ptherm(ix^D)
         !v<== it's noted rho_ but in fact the passed array wCT is equalt to nH
         nHlocal = wCT(ix^D,rho_)
         ! ntot = <m>/mu*nH
         nlocal = (coolconfig%mean_mass/coolconfig%mean_mup)*wCT(ix^D,rho_)
         ! input from mod_source.t that will call this subroutine
         ! is not the total mass rho but the hydrogen number density nH



         emin     = kB*nlocal*tlow/(rc_gammamin1) ! =P/(gamma-1)
         Lmax     = max(zero,(pnew(ix^D)/(rc_gammamin1)-emin)/qdt)!dimension of Lambda_Hd*nH*nH

         !  Tlocal = P/rho
         Tlocal1   = max(temperature(ix^D),smalldouble)
         !
         !  Determine explicit cooling
         !
         !  If temperature is below floor level, no cooling.
         !  Stop wasting time and go to next gridpoint.
         !  If the temperature is higher than the maximum,
         !  assume Bremmstrahlung
         if( Tlocal1<=max(tcoolmin,tlow) ) then
            L1 = zero
         else if( Tlocal1>=tcoolmax )then
            L1         = Lcool(ncool)*sqrt(Tlocal1/tcoolmax)
            L1         = L1*(nlocal*nlocal)
            L1         = min(L1,Lmax)
            w(ix^D,e_) = w(ix^D,e_)-L1*qdt
         else
            call findL(Tlocal1,nlocal,L1) !get Lambda_HD
            call findY(Tlocal1,nlocal,Y1) !get Y(T)
            !tc         = Tlocal1/(rc_gammamin1)/(nlocal*L1) !=t_cool in Townsend 2009
            !Y2         = Y1+(Tlocal1*fact)/(L1*tc)
            ! tc=t_cool in Townsend 2009
            ! tc = n_all * kB* T / ((gamma-1)*nH*nH*Lambda_HD)
            tc         = nlocal*kB*Tlocal1/(rc_gammamin1)/(nHlocal*nHlocal*L1)
            Y2         = Y1+(Tlocal1*fact)/(L1*tc) ! must be dimensionless
            call findT(Tlocal2,Y2)

            if(Tlocal2<=tcoolmin) then
              L1 = Lmax/(nHlocal*nHlocal) !since L1 will be multipled by (nHlocal*nHlocal)
            else
              L1 = ((Tlocal1-Tlocal2)*kB*nlocal)/(rc_gammamin1)/(nHlocal*nHlocal*qdt) !missed some terms
            endif

            L1          = L1*(nHlocal*nHlocal)!transform L1 from Lambda_HD to S=Lambda_HD*nH*nH
            L1          = min(L1,Lmax)
            !if(it==1) write(*,*) w(ix^D,e_), w_convert_factor(e_), w(ix^D,e_)*w_convert_factor(e_)
            w(ix^D,e_)  = w(ix^D,e_)-L1*qdt
         endif
         if(coolsavedT)w(ix^D,dtcool1_) = cfrac*(pnew(ix^D)/rc_gammamin1/max(L1,smalldouble))
         if(coolsaveL)w(ix^D,Lcool1_) = L1

      {enddo^D&\}

    end subroutine cool_exact

    subroutine findL (tpoint,rhopoint,Lpoint)
    !  Fast search option to find correct point
    !  in cooling curve
      use mod_global_parameters

      real(kind=dp)   ,intent(IN)   :: tpoint,rhopoint
      real(kind=dp)   , intent(OUT) :: Lpoint
      integer :: ipoint
      integer :: jl,jc,jh
      real(kind=dp)    :: lgtp

      lgtp = dlog10(tpoint)
      jl = int((lgtp - lgtcoolmin) / lgstep) + 1
      Lpoint = Lcool(jl)+ (tpoint-tcool(jl)) &
                * (Lcool(jl+1)-Lcool(jl)) &
                / (tcool(jl+1)-tcool(jl))

!      if (tpoint == tcoolmin) then
!        Lpoint = Lcool(1)
!      else if (tpoint == tcoolmax) then
!        Lpoint = Lcool(ncool)
!      else
!        jl=0
!        jh=ncool+1
!        do
!          if (jh-jl <= 1) exit
!          jc=(jh+jl)/2
!          if (tpoint >= tcool(jc)) then
!              jl=jc
!          else
!              jh=jc
!          end if
!        end do
!        ! Linear interpolation to obtain correct cooling
!        Lpoint = Lcool(jl)+ (tpoint-tcool(jl)) &
!                  * (Lcool(jl+1)-Lcool(jl)) &
!                  / (tcool(jl+1)-tcool(jl))
!      end if

    end subroutine findL

    subroutine findY (tpoint,rhopoint,Ypoint)

    !  Fast search option to find correct point in cooling time

      use mod_global_parameters

      real(kind=dp)   ,intent(IN)   :: tpoint,rhopoint
      real(kind=dp)   , intent(OUT) :: Ypoint
      integer :: ipoint
      integer :: jl,jc,jh
      real(kind=dp)    :: lgtp

      lgtp = dlog10(tpoint)
      jl = int((lgtp - lgtcoolmin) / lgstep) + 1
      Ypoint = Yc(jl)+ (tpoint-tcool(jl)) &
                * (Yc(jl+1)-Yc(jl)) &
                / (tcool(jl+1)-tcool(jl))

  !    integer i
  !
  !    if (tpoint == tcoolmin) then
  !      Ypoint = Yc(1)
  !    else if (tpoint == tcoolmax) then
  !      Ypoint = Yc(ncool)
  !    else
  !      jl=0
  !      jh=ncool+1
  !      do
  !        if (jh-jl <= 1) exit
  !        jc=(jh+jl)/2
  !        if (tpoint >= tcool(jc)) then
  !           jl=jc
  !        else
  !           jh=jc
  !        end if
  !      end do
  !      ! Linear interpolation to obtain correct value
  !      Ypoint = Yc(jl)+ (tpoint-tcool(jl)) &
  !                * (Yc(jl+1)-Yc(jl)) &
  !                / (tcool(jl+1)-tcool(jl))
  !    end if

    end subroutine findY

    subroutine findT (tpoint,Ypoint)

    !  Fast search option to find correct temperature
    !  from cooling time. Only possible this way because T is a monotonously
    !  decreasing function
      use mod_global_parameters

      real(kind=dp)   ,intent(OUT)   :: tpoint
      real(kind=dp)   , intent(IN) :: Ypoint
      integer :: ipoint
      integer :: jl,jc,jh
      integer i

      if (Ypoint >= Yc(1)) then
        tpoint = tcoolmin
      else if (Ypoint == Yc(ncool)) then
        tpoint = tcoolmax
      else
        jl=0
        jh=ncool+1
        do
          if (jh-jl <= 1) exit
          jc=(jh+jl)/2
          if (Ypoint <= Yc(jc)) then
              jl=jc
          else
              jh=jc
          end if
        end do
        ! Linear interpolation to obtain correct temperature
        tpoint = tcool(jl)+ (Ypoint-Yc(jl)) &
                  * (tcool(jl+1)-tcool(jl)) &
                  / (Yc(jl+1)-Yc(jl))
      end if
    end subroutine findT

    subroutine finddLdt (tpoint,dLpoint)

    !  Fast search option to find correct point
    !  in derivative of cooling curve
      use mod_global_parameters

      real(kind=dp)   ,intent(IN)   :: tpoint
      real(kind=dp)   , intent(OUT) :: dLpoint
      integer :: ipoint
      integer :: jl,jc,jh
      real(kind=dp)    :: lgtp

      lgtp = dlog10(tpoint)
      jl = int((lgtp - lgtcoolmin) / lgstep) + 1
      dLpoint = dLdtcool(jl)+ (tpoint-tcool(jl)) &
                * (dLdtcool(jl+1)-dLdtcool(jl)) &
                / (tcool(jl+1)-tcool(jl))

!      if (tpoint == tcoolmin) then
!        dLpoint = dLdtcool(1)
!      else if (tpoint == tcoolmax) then
!        dLpoint = dLdtcool(ncool)
!      else
!        jl=0
!        jh=ncool+1
!        do
!          if (jh-jl <= 1) exit
!          jc=(jh+jl)/2
!          if (tpoint >= tcool(jc)) then
!              jl=jc
!          else
!              jh=jc
!          end if
!        end do
!        ! Linear interpolation to obtain correct cooling derivative
!        dLpoint = dLdtcool(jl)+ (tpoint-tcool(jl)) &
!                  * (dLdtcool(jl+1)-dLdtcool(jl)) &
!                  / (tcool(jl+1)-tcool(jl))
!      end if
    end subroutine finddLdt


  subroutine chemical_get_electron_density(ixI^L, ixO^L,w,ndensity_electron)
    use mod_global_parameters
    implicit none
    integer, intent(in)            :: ixI^L, ixO^L
    real(kind=dp), intent(in)      :: w(ixI^S,1:nw)
    real(dp), intent(out)          :: ndensity_electron(ixO^S)
    !-------------------------------------
    ndensity_electron = w(ixO^S,rho_)
  end subroutine chemical_get_electron_density
end module mod_radiative_cooling
