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
  implicit none
  save
  real(dp) :: theta, kx, ly, vc

  type usr_config
    logical           :: physunit_on
    logical           :: ism_on
    logical           :: cloud_on
    logical           :: jet_yso_on
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


  type (ISM),allocatable,target      :: ism_surround(:)
  type (cloud),allocatable,target    :: cloud_medium(:)
  type (cla_jet),allocatable,target  :: jet_yso(:)
  type (ISM),target                  :: ism_default
  type (cloud),target                :: cloud_default
  type (cla_jet),target              :: jet_yso_default
  type (dust),target                 :: dust_mialy
  type (dust),allocatable,target     :: the_dust_inuse(:)

  !type(star) :: star_ms
  !type(star) :: sun

  type(usrphysical_unit) :: usr_physunit





contains

  !-----------------------------------------------------------
  !> Calculate gravitational acceleration in each dimension
  subroutine test_vacpp(ixI^L,ixO^L,wCT,x,gravity_field)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    real(dp), intent(in)    :: x(ixI^S,1:ndim)
    real(dp), intent(in)    :: wCT(ixI^S,1:nw)
    real(dp), intent(out)   :: gravity_field(ixI^S,ndim)
    real(dp)                        :: Ggrav, Mpoint
    real(kind=dp), dimension(ixI^S) :: r_distance
    real(kind=dp), dimension(ixI^S) :: theta_profile
    real(kind=dp), dimension(1:ndim) :: zero_dim
    integer                      :: i_ism,idim
    !-----------------------------------

    ! surbroutine to compute the matrix to transform the coordinates
    ! from an avalaible coordinate system to another
    ! 20/02/22 : avalaible coordinate systems :
    ! - Cartesian (or 'slab'): indexes x_, y_, z_ = 1,2,3
    ! - Cylindrical: indexes r_, z_, phi_ = 1,2,3
    ! - Polar: indexes r_, phi_, z_ = 1,2,3
    ! - Spherical: indexes r_, theta_ , phi_ = 1,2,3
    !
    ! Then there are 4x3 = 12 possible transformations ('X'=Done):
    ! 1) Cartesian -> Cylindrical
    ! 2) Cartesian -> Polar
    ! 3) Cartesian -> Spherical
    ! 4) Cylindrical -> Polar
    ! 5) Cylindrical -> Spherical
    ! 6) Cylindrical -> Cartesian
    ! 7) Polar -> Cylindrical
    ! 8) Polar -> Spherical
    ! 9) Polar -> Cartesian
    ! 10) Spherical -> Cylindrical
    ! 11) Spherical -> Polar
    ! 12) Spherical -> Cartesian
    ! Along with 4 trivial transformations
    ! 13) Cartesian -> Cartesian
    ! 14) Cylindrical -> Cylindrical
    ! 15) Polar -> Polar
    ! 16) Spherical -> Spherical
    subroutine transform_matrix(ixI^L,ixO^L,x,csI,csO,coord,tmatrix,xplus)
      use mod_global_parameters
      integer, intent(in)                 :: ixI^L, ixO^L
      real(dp), intent(in)                :: x(ixI^S,1:ndim)
      character(len=std_len), intent(in)  :: csI,csO,coord
      real(dp), intent(out),allocatable              :: tmatrix(:,:,:)
      real(dp), intent(in), optional      :: xplus(1:3) ! if angle, must be in radians
      ! .. local ..
      real(kind=dp), dimension(ixI^S) :: xr,xtheta,xphi
      real(kind=dp), dimension(1:ndim) :: zero_dim
      integer                      :: idim,nvecdim,i1,i2,i3
      !-----------------------------------


      !Determine nvecdim
      select case(coord)
        case('x','y','z','R','r','phi','theta')
          nvecdim=1
        case('xy','xz','yz','Rz','Rphi','zphi','phiz','rtheta','rphi','thetaphi')
          nvecdim=2
        case('xyz','Rzphi','Rphiz','rthetaphi')
          nvecdim=3
        case default
          write(*,*) 'Coordinate system in input is', csI
          call mpistop("Unknown coord specified in input")
      end select

      if(allocated(tmatrix))deallocate(tmatrix)
      allocate(tmatrix(ixO^S,nvecdim,nvecdim))

      ! Identity matrix
      tmatrix(ixO^S,1:nvecdim,1:nvecdim) = zero
      do idim=1,nvecdim
        tmatrix(ixO^S,idim,idim)=one
      enddo

      if(nvecdim>ndim)then

      elseif(nvecdim<ndim)then

      else !nvecdim==ndim
      end if




  end subroutine test_vacpp
end module mod_usr
