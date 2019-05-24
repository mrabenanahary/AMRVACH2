module mod_obj_sn_remnant
use mod_constants
use mod_global_parameters
use mod_srmhd_parameters
use mod_obj_dust, only : dust
use mod_obj_global_parameters
use mod_obj_mat
use mod_obj_usr_unit
implicit none
  type supernovae_remnant_parameters
    character(len=20)    :: unit !> physical unit at parameter file
    logical              :: normalize_done !> supernovae remnant is the normalisation is already done
    real(dp)             :: center(3)           !> supernovae remnant center
    real(dp)             :: r_in !> supernovae remnant inner radius
    real(dp)             :: r_out !> supernovae remnant outer radius
    real(dp)             :: velocity_init(3) !> supernovae remnant initial speed (spherical outward from center)
    real(dp)             :: lfac_init !> supernovae remnant initial Lorentz factor
    real(dp)             :: velocity_proper(3) !> supernovae remnant initial speed  (cartesian from center)
    real(dp)             :: density_init !> supernovae remnant initial density
    real(dp)             :: number_density_init !> supernovae remnant initial density
    real(dp)             :: pressure_init !> supernovae remnant initial pressure
    real(dp)             :: temperature_init !> supernovae remnant initial temperature
    real(dp)             :: magnetic_init(3) !> supernovae remnant initial magnetic field
    real(dp)             :: xisigma0_init !> supernovae remnant initial sigma
    real(dp)             :: mass                !> supernovae remnant mass
    real(dp)             :: energy              !> supernovae remnant energy
    real(dp)             :: energy_kinetic !> supernovae remnant kinetic energy
    real(dp)             :: energy_thermal !> supernovae remnant thermal energy
    real(dp)             :: power_speed_increase !> supernovae remnant power velocity decrease
    real(dp)             :: volume_init !> supernovae remnant initial  volume
    real(dp)             :: time_set !> supernovae remnant starting time
    logical              :: dust_on             !> cloud with dust in is true
    real(dp)             :: dust_frac           !> dust fraction
    character(len=20)    :: dust_profile        !> could dust inside profile
    logical              :: tracer_on           !> supernovae remnant tracer
    integer              :: itr                 !> supernovae indice
  end type supernovae_remnant_parameters

  type supernovae_remnant
    logical, allocatable                :: patch(:) !> supernovae remnant region
    logical, allocatable                :: patch_escape(:) !> supernovae remnant  not on cell
    type(supernovae_remnant_parameters) :: myconfig !>  supernovae remnant parameters to be read
    type(usrphysical_unit), pointer     :: myphysunit !> supernovae remnant physics unit in use
    type (dust)                         :: mydust !> supernovae remnant dust

    character(len=78)                   :: subname !> subroutine name that call it
  contains
   PROCEDURE, PASS(self)        :: set_default     => &
      usr_supernovae_remnant_set_default
   PROCEDURE, PASS(self)        :: set_complet     => &
      usr_supernovae_remnant_set_complet
   PROCEDURE, PASS(self)        :: normalize       => &
      usr_supernovae_remnant_normalize
   PROCEDURE, PASS(self)        :: set_w           => &
      usr_supernovae_remnant_set_w
   PROCEDURE, PASS(self)        :: read_parameters => &
      usr_supernovae_remnant_read_p
   PROCEDURE, PASS(self)        :: write_setting   => &
      usr_supernovae_remnant_write_setting
   PROCEDURE, PASS(self)        :: spd_rad_to_cart => &
      usr_supernovae_remnant_spd_rad_to_cart
   PROCEDURE, PASS(self)        :: clean_memory    => &
      usr_supernovae_remnant_clean_memory
   PROCEDURE, PASS(self)        :: get_dt          => &
      usr_supernovae_remnant_get_dt
   PROCEDURE, PASS(self)        :: get_patch       => &
      usr_supernovae_remnant_patch
  end type supernovae_remnant
contains

!-------------------------------------------------------------------------
 !> Read the ism parameters  from a parfile
  subroutine usr_supernovae_remnant_read_p(self,supernovae_remnant_config,&
     files)
    class(supernovae_remnant)              :: self
    character(len=*), intent(in)           :: files(:)
    type(supernovae_remnant_parameters)    :: supernovae_remnant_config
    ! .. local ..
    integer  :: i_file, i_error_read
    namelist /usr_supernovae_remnant_list/ supernovae_remnant_config

    if(mype==0)write(*,*)'Reading usr_supernovae_remnant_list'
    do i_file = 1, size(files)
       open(unitpar, file=trim(files(i_file)), status="old")
       read(unitpar, usr_supernovae_remnant_list, iostat=i_error_read)

       if(i_error_read>0)then
        write(*,*&
           )'At user side in mod_obj_sn_remnant: error at reading parfile'
        write(*,*)'Check input.  Something was wrong, it will stop'
        call mpistop('It stops at reading the parfile')
       elseif(i_error_read<0)then
        write(*,*&
           )'At user side in mod_obj_sn_remnant: error at reading parfile'
        write(*,*)'Reach the end of the file it will stop'
        call mpistop('It stops at reading the parfile')
       else
         write(*,*) 'En of Reading usr_supernovae_remnant_list'
       end if
    end do

    if(supernovae_remnant_config%dust_on)then
      self%mydust%myconfig%associated_medium = 'sn'
      call self%mydust%read_parameters(self%mydust%myconfig,files)
    end if
  end subroutine usr_supernovae_remnant_read_p

!------------------------------------------------------------------------
!> write the cloud setting
subroutine usr_supernovae_remnant_write_setting(self,unit_config)
  implicit none
  class(supernovae_remnant)                :: self
  integer,intent(in)                  :: unit_config
  ! .. local ..

  !-----------------------------------

  write(unit_config,*)'*********************************************'
  write(unit_config,*)'********** supenovae remnant setting ********'
  write(unit_config,*)'*********************************************'
  write(unit_config,*)'      ****** Code Unit *******      '
  write(unit_config,*) 'Density     = ', self%myconfig%density_init
  write(unit_config,*) 'Pressure    = ', self%myconfig%pressure_init
  write(unit_config,*) 'Energy      = ', self%myconfig%energy
  write(unit_config,*) 'Mass        = ', self%myconfig%mass
  write(unit_config,*) 'Speed       = ', self%myconfig%velocity_init
  write(unit_config,*)'      ****** Physical Unit *******   '
  write(unit_config,*) 'Density     = ', &
     self%myconfig%density_init*self%myphysunit%myconfig%density,'  ',&
     self%myphysunit%myunit%density
  write(unit_config,*) 'Pressure    = ', &
     self%myconfig%pressure_init*self%myphysunit%myconfig%pressure,'  ',&
     self%myphysunit%myunit%pressure
  write(unit_config,*) 'Energy      = ', &
     self%myconfig%energy*self%myphysunit%myconfig%energy,'  ',&
     self%myphysunit%myunit%energy
  write(unit_config,*) 'Mass        = ', &
     self%myconfig%mass*self%myphysunit%myconfig%mass,'  ',&
     self%myphysunit%myunit%mass
  write(unit_config,*) 'Speed       = ', self%myconfig%velocity_init  &
     *self%myphysunit%myconfig%velocity,'  ',self%myphysunit%myunit%velocity
  write(unit_config,*)'*********************************************'
  write(unit_config,*)'******* END supenovae remnant setting *******'
  write(unit_config,*)'*********************************************'
end    subroutine usr_supernovae_remnant_write_setting
!--------------------------------------------------------------------

!> subroutine default setting for cloud
 subroutine usr_supernovae_remnant_set_default(self)
  implicit none
  class(supernovae_remnant)          :: self
  !----------------------------------
  self%myconfig%unit                   = 'code'
  self%myconfig%density_init           = 0.0_dp
  self%myconfig%number_density_init    = 0.0_dp
  self%myconfig%temperature_init       = 0.0_dp
  self%myconfig%pressure_init          = 0.0_dp
  self%myconfig%center(:)              = 0.0_dp !(box_limit(2,:)-box_limit(1,:))/2.0_dp
  self%myconfig%r_in                   = 0.0_dp !(box_limit(2,:)+box_limit(1,:))/2.0_dp
  self%myconfig%r_out                  = 0.0_dp
  self%myconfig%velocity_init(:)       = 0.0_dp
  self%myconfig%lfac_init              = 0.0_dp
  self%myconfig%velocity_proper(:)     = 0.0_dp
  self%myconfig%magnetic_init(:)       = 0.0_dp
  self%myconfig%xisigma0_init          = 0.0_dp
  self%myconfig%power_speed_increase   = 0.0_dp
  self%myconfig%time_set               = 0.0_dp
  self%myconfig%tracer_on              = .false.
  self%myconfig%dust_on                = .false.
  self%myconfig%normalize_done         = .false.
 end subroutine usr_supernovae_remnant_set_default
 !--------------------------------------------------------------------
 !> subroutine check the parfile setting for cloud
 subroutine usr_supernovae_remnant_set_complet(self)
   implicit none
   class(supernovae_remnant)          :: self
   ! .. local ..
   real(dp)                 :: mp,kb
   logical                  :: dust_is_frac
   !-----------------------------------
   if(SI_unit) then
     mp=mp_SI
     kB=kB_SI
   else
     mp=mp_cgs
     kB=kB_cgs
   end if



  cond_lfac_set : if (self%myconfig%lfac_init>0.0) then
     self%myconfig%velocity_init(r_)     = &
        dsqrt(1.0_dp-1.0_dp/self%myconfig%lfac_init**2.0_dp)
     self%myconfig%velocity_init(theta_) = 0.0_dp
     self%myconfig%velocity_init(phi_)   = 0.0_dp
    select case(trim(self%myconfig%unit))
     case('code')
      ! do no thing
     case default
       self%myconfig%velocity_init=unit_velocity*self%myconfig%velocity_init
     end select
  else cond_lfac_set
     select case(trim(self%myconfig%unit))
      case('code')
       if(phys_config%isrel)then
         self%myconfig%lfac_init=1.0_dp/dsqrt(1.0_dp-&
            SUM(self%myconfig%velocity_init**2.0_dp))
       else
         self%myconfig%lfac_init=1.0_dp/dsqrt(1.0_dp-&
            SUM(self%myconfig%velocity_init**2.0_dp)* unit_velocity/const_c)
       end if
      case default
       self%myconfig%lfac_init=1.0_dp/dsqrt(1.0_dp-&
          SUM((self%myconfig%velocity_init/const_c)**2.0_dp))
     end select
  end if cond_lfac_set



  self%myconfig%volume_init = 4.0_dp/3.0_dp* &
     dpi*(self%myconfig%r_out**3.0_dp-self%myconfig%r_in**3.0_dp)



  cond_mass_set : if(self%myconfig%mass>0.0_dp) then
   self%myconfig%density_init = self%myconfig%mass*               &
      (1.0_dp-self%myconfig%dust_frac)  /self%myconfig%volume_init
   self%myconfig%number_density_init =self%myconfig%density_init/mp
  else cond_mass_set
   if(dabs(self%myconfig%density_init)<=0) then
    self%myconfig%density_init =self%myconfig%number_density_init*mp
   end if

   if(dabs(self%myconfig%number_density_init)<=0) then
    self%myconfig%number_density_init =self%myconfig%density_init/mp
   end if
   self%myconfig%mass = self%myconfig%density_init*self%myconfig%volume_init
  end if cond_mass_set



  cond_energy_set : if(dabs(self%myconfig%energy)>0.0_dp) then

    self%myconfig%energy_kinetic = 0.5_dp* self%myconfig%mass*                 &
          (1.0_dp-self%myconfig%dust_frac)*               &
       (1.0_dp-1.0_dp/self%myconfig%lfac_init**2.0_dp)*unit_velocity**2.0_dp


    self%myconfig%energy_kinetic = 3.0/5.0*  self%myconfig%energy_kinetic

    self%myconfig%energy_thermal = self%myconfig%energy - &
       self%myconfig%energy_kinetic

    check_thermal_energy : if(self%myconfig%energy_thermal<0.0_dp ) then
      write(*,*) ' The kinetic energy si large ', self%myconfig%energy_thermal
      write(*,*) ' You should decrease the initial SN speed '
      call mpistop&
         ('stop at the mod_usr.t in usr_supernovae_remnant_set_complet')
    end if check_thermal_energy


    self%myconfig%pressure_init  = (srmhd_gamma-1.0_dp) * &
       self%myconfig%energy_thermal/self%myconfig%volume_init

    self%myconfig%temperature_init=   self%myconfig%pressure_init / &
       (kb*self%myconfig%number_density_init)


  else   cond_energy_set
   if(dabs(self%myconfig%pressure_init)<=0.and.&
      dabs(self%myconfig%temperature_init)<=0) then
      write(*,*)' the temperature and pressure of envelope are not set'
      call mpistop('the code stop at envelope object in the user file')
   elseif(dabs(self%myconfig%pressure_init)<smalldouble) then
    self%myconfig%pressure_init&
        =self%myconfig%number_density_init*kB*self%myconfig%temperature_init
   end if
 end if cond_energy_set


  if(.not.self%myconfig%tracer_on)self%myconfig%itr=0

   cond_duston : if(self%myconfig%dust_on)then

     dust_is_frac=.true.
     call self%mydust%set_complet(dust_is_frac,self%myconfig%dust_frac,&
        self%myconfig%density_init,self%myconfig%velocity_init)
   end if cond_duston

  if(self%myconfig%tracer_on)then
       prim_wnames(self%myconfig%itr+(phys_ind%tracer(1)-1)) = 'tracer_sn'
       cons_wnames(self%myconfig%itr+(phys_ind%tracer(1)-1)) = 'tracer_sn'
  end if
 end subroutine usr_supernovae_remnant_set_complet
!--------------------------------------------------------------------
 subroutine usr_supernovae_remnant_normalize(self,physunit_inuse)
  use mod_obj_usr_unit
  implicit none
  class(supernovae_remnant)                       :: self
  type(usrphysical_unit),target, intent(in)       :: physunit_inuse
  !----------------------------------
  !----------------------------------
  self%myphysunit =>physunit_inuse
  if(trim(self%myconfig%unit)=='code'.or.self%myconfig%normalize_done)then
   if(self%myconfig%normalize_done)then
    write(*,*) 'WARNING: Second call for supernovae remnant normalisation ',&
        'no new normalisation will be done'
   end if
   return
  end if
  self%myconfig%density_init     = self%myconfig%density_init      &
     /physunit_inuse%myconfig%density
  self%myconfig%temperature_init = self%myconfig%temperature_init  &
     /physunit_inuse%myconfig%temperature
  self%myconfig%pressure_init    = self%myconfig%pressure_init     &
     /physunit_inuse%myconfig%pressure
  self%myconfig%velocity_init    = self%myconfig%velocity_init     &
     /physunit_inuse%myconfig%velocity
  self%myconfig%velocity_proper  = self%myconfig%velocity_proper   &
     /physunit_inuse%myconfig%velocity
  self%myconfig%magnetic_init    = self%myconfig%magnetic_init     &
     /physunit_inuse%myconfig%magnetic
  self%myconfig%center           = self%myconfig%center            &
     /physunit_inuse%myconfig%length
  self%myconfig%r_in             = self%myconfig%r_in              &
     /physunit_inuse%myconfig%length
  self%myconfig%r_out            = self%myconfig%r_out             &
     /physunit_inuse%myconfig%length
  self%myconfig%mass             = self%myconfig%mass              &
     /physunit_inuse%myconfig%mass
  self%myconfig%energy           = self%myconfig%energy            &
     /physunit_inuse%myconfig%energy
  self%myconfig%energy_kinetic   = self%myconfig%energy_kinetic    &
     /physunit_inuse%myconfig%energy
  self%myconfig%energy_thermal   = self%myconfig%energy_thermal    &
     /physunit_inuse%myconfig%energy
  self%myconfig%volume_init      = self%myconfig%volume_init       &
     /physunit_inuse%myconfig%volume

  self%myconfig%time_set         = self%myconfig%time_set          &
     /physunit_inuse%myconfig%time
  cond_duston : if(self%myconfig%dust_on)then
    call self%mydust%normalize(physunit_inuse)
    call self%mydust%to_phys()
  end if cond_duston
 end subroutine usr_supernovae_remnant_normalize
!--------------------------------------------------------------------
!--------------------------------------------------------------------
 !> subroutine patch for the cloud
 subroutine usr_supernovae_remnant_patch(ixImin1,ixImax1,ixOmin1,ixOmax1,qt,x,&
    self)
  implicit none
  integer, intent(in)        :: ixImin1,ixImax1,ixOmin1,ixOmax1
  real(kind=dp), intent(in)  :: qt
  real(kind=dp), intent(in)  :: x(ixImin1:ixImax1,1:ndir)
  class(supernovae_remnant)  :: self
  real(dp), dimension(ixImin1:ixImax1) :: dist
  !----------------------------------

  allocate(self%patch(ixGlo1:ixGhi1))
  ! E = 0.5 rho*v2 int (R^4/r_out^2 4 pi *dR)=0.5/5 *4 pi *rho*v^2*r_out^3 = 0.5  rho*v2*Volume *3/5=Ek*3/5
  call usr_distance(ixImin1,ixImax1,ixOmin1,ixOmax1,typeaxial,&
     self%myconfig%center,x,dist)
  self%patch(ixOmin1:ixOmax1) = Dist(ixOmin1:ixOmax1) >self%myconfig%r_in &
     .and. Dist(ixOmin1:ixOmax1) <self%myconfig%r_out
  if(allocated(self%patch_escape))then
    self%patch(ixOmin1:ixOmax1) = self%patch(ixOmin1:ixOmax1).and.&
       (.not.self%patch_escape(ixOmin1:ixOmax1))
  end if
 end subroutine usr_supernovae_remnant_patch
!--------------------------------------------------------------------
 !> subroutine setting for cloud
 subroutine usr_supernovae_remnant_set_w(ixImin1,ixImax1,ixOmin1,ixOmax1,qt,x,&
    w,self)
  implicit none
  integer, intent(in)          :: ixImin1,ixImax1,ixOmin1,ixOmax1
  real(kind=dp), intent(in)    :: qt
  real(kind=dp)                :: x(ixImin1:ixImax1,1:ndir)
  real(kind=dp)                :: w(ixImin1:ixImax1,1:nw)
  class(supernovae_remnant)         :: self
  ! .. local..
  integer                      :: idir
  real(kind=dp)                :: Dist(ixImin1:ixImax1)
  logical                      :: dust_is_frac
  real(kind=dp)              :: fprofile(ixImin1:ixImax1)
  !----------------------------------

  call usr_supernovae_remnant_patch(ixImin1,ixImax1,ixOmin1,ixOmax1,qt,x,self)
  call usr_distance(ixImin1,ixImax1,ixOmin1,ixOmax1,typeaxial,&
     self%myconfig%center,x,dist)

  where(self%patch(ixOmin1:ixOmax1))
   w(ixOmin1:ixOmax1,phys_ind%rho_)           = self%myconfig%density_init
   w(ixOmin1:ixOmax1,phys_ind%pressure_)             = &
      self%myconfig%pressure_init

   w(ixOmin1:ixOmax1,phys_ind%mom(r_))        = &
      self%myconfig%velocity_init(r_) * Dist(&
      ixOmin1:ixOmax1)/self%myconfig%r_out

   w(ixOmin1:ixOmax1,phys_ind%mom(theta_))    = &
      self%myconfig%velocity_init(theta_)
   w(ixOmin1:ixOmax1,phys_ind%mom(phi_))      = &
      self%myconfig%velocity_init(phi_)
  endwhere
  if(phys_config%isrel) then
   where(self%patch(ixOmin1:ixOmax1))w(ixOmin1:ixOmax1,&
      phys_ind%lfac_)          = self%myconfig%lfac_init
  end if
  if(phys_config%ismhd)then
   where(self%patch(ixOmin1:ixOmax1))
    w(ixOmin1:ixOmax1,phys_ind%mag(r_))        = &
       self%myconfig%magnetic_init(1)
    w(ixOmin1:ixOmax1,phys_ind%mag(theta_))    = &
       self%myconfig%magnetic_init(2)
    w(ixOmin1:ixOmax1,phys_ind%mag(phi_))      = &
       self%myconfig%magnetic_init(3)
   end where
  end if





  if(trim(typeaxial)/='spherical')call self%spd_rad_to_cart(ixImin1,ixImax1,&
     ixOmin1,ixOmax1,x,w)


  call usr_lorentz_transrmation_add_proper_speed(ixImin1,ixImax1,ixOmin1,&
     ixOmax1,self%myconfig%velocity_proper,w(ixImin1:ixImax1,&
     phys_ind%mom(1):phys_ind%mom(ndir)))

  cond_tracer_on :if(self%myconfig%tracer_on.and.phys_config%n_tracer>0.and.&
     self%myconfig%itr<=phys_config%n_tracer)then
   where(self%patch(ixOmin1:ixOmax1))
    w(ixOmin1:ixOmax1,phys_ind%tracer(self%myconfig%itr)) = 1.0d2 !w(ixOmin1:ixOmax1,rho_)
   elsewhere
    w(ixOmin1:ixOmax1,phys_ind%tracer(self%myconfig%itr)) = 0.0_dp
   end where
  end if cond_tracer_on


  cond_dust_on : if(self%myconfig%dust_on)then
      if(self%myconfig%dust_profile/='none') then
       call usr_mat_profile(ixImin1,ixImax1,ixOmin1,ixOmax1,typeaxial,&
          self%myconfig%dust_profile,self%myconfig%center,&
          (/self%myconfig%r_out,0.0_dp,0.0_dp/),x,fprofile)
      else
        fprofile = 1.0_dp
      end if
     call self%mydust%set_patch(ixImin1,ixImax1,ixOmin1,ixOmax1,self%patch)
     self%mydust%myconfig%velocity=self%myconfig%velocity_init
      dust_is_frac=.true.
    call self%mydust%set_w(ixImin1,ixImax1,ixOmin1,ixOmax1,qt,dust_is_frac,&
       self%myconfig%dust_frac,fprofile,x,w)
  end if cond_dust_on
 end subroutine usr_supernovae_remnant_set_w
!===========================================================
!> Subroutine to set time set before supernovae_remnant starts
subroutine usr_supernovae_remnant_get_dt(self,ixImin1,ixImax1,ixOmin1,ixOmax1,&
   dx1,x,w,qt,dtnew)
  use mod_global_parameters
  class(supernovae_remnant)       :: self
  integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
  double precision, intent(in)    :: dx1,qt, x(ixImin1:ixImax1,1:ndim)
  double precision, intent(in)    :: w(ixImin1:ixImax1,1:nw)
  double precision, intent(inout) :: dtnew
  !--------------------------------------------------------------
  if(qt<self%myconfig%time_set)then
   dtnew=min(self%myconfig%time_set-qt,dtnew)
  end if
end subroutine usr_supernovae_remnant_get_dt
!===========================================================

!===========================================================
!> Subroutine to convert radial wind from spherical to cartesian coordinates
subroutine usr_supernovae_remnant_spd_rad_to_cart(ixImin1,ixImax1,ixOmin1,&
   ixOmax1,x,w,self)
  use mod_global_parameters
  implicit none
  integer, intent(in)             :: ixImin1,ixImax1,ixOmin1,ixOmax1
  real(kind=dp), intent(in)       :: x(ixImin1:ixImax1,1:ndir)
  real(kind=dp), intent(inout)    :: w(ixImin1:ixImax1,1:nw)
  class(supernovae_remnant)            :: self
  ! .. local ..
  real(kind=dp), dimension(ixImin1:ixImax1,1:ndir) :: v_spherical
  real(kind=dp), dimension(ixImin1:ixImax1)        :: sin_theta,cos_theta
  real(kind=dp), dimension(ixImin1:ixImax1)        :: sin_phi,cos_phi
  real(kind=dp), dimension(ixImin1:ixImax1)        :: Dist
  !--------------------------------------------------

 call usr_distance(ixImin1,ixImax1,ixOmin1,ixOmax1,typeaxial,&
    self%myconfig%center,x,dist)

 select case(typeaxial)
  case('slab')
   v_spherical(ixOmin1:ixOmax1,1:ndir) = w(ixOmin1:ixOmax1,&
      phys_ind%mom(1):phys_ind%mom(ndir))
   if(z_in) then
     sin_theta(ixOmin1:ixOmax1)=dsqrt(x(ixOmin1:ixOmax1,&
        x_)**2.0_dp+x(ixOmin1:ixOmax1,y_)**2.0_dp)/Dist(ixOmin1:ixOmax1)
     cos_theta(ixOmin1:ixOmax1)=x(ixOmin1:ixOmax1,z_)/Dist(ixOmin1:ixOmax1)
   else
    sin_theta(ixOmin1:ixOmax1)=1.0_DP
    cos_theta(ixOmin1:ixOmax1)=0.0_dp
   end if
   if(y_in) then
     sin_phi(ixOmin1:ixOmax1)=dsin(x(ixOmin1:ixOmax1,phi_))
     cos_phi(ixOmin1:ixOmax1)=dcos(x(ixOmin1:ixOmax1,phi_))
   else
    sin_phi(ixOmin1:ixOmax1)=0.0_DP
    cos_phi(ixOmin1:ixOmax1)=1.0_dp
   end if
   where(self%patch(ixOmin1:ixOmax1))
    w(ixOmin1:ixOmax1,phys_ind%mom(x_)) = v_spherical(ixOmin1:ixOmax1,&
       r_)*cos_phi(ixOmin1:ixOmax1)*sin_theta(ixOmin1:ixOmax1)

    where(sin_phi(ixOmin1:ixOmax1)>0.0_dp)
     w(ixOmin1:ixOmax1,phys_ind%mom(y_)) = v_spherical(ixOmin1:ixOmax1,&
        r_)*sin_phi(ixOmin1:ixOmax1)*sin_theta(ixOmin1:ixOmax1)
    else where
     w(ixOmin1:ixOmax1,phys_ind%mom(y_)) = zero
    end where
    where(cos_phi(ixOmin1:ixOmax1)>0.0_dp)
     w(ixOmin1:ixOmax1,phys_ind%mom(z_)) = v_spherical(ixOmin1:ixOmax1,&
        r_)*cos_theta(ixOmin1:ixOmax1)
    else where
     w(ixOmin1:ixOmax1,phys_ind%mom(z_)) = zero
    end where

   end where
  case('cylindrical')
    v_spherical(ixOmin1:ixOmax1,1:ndir) = w(ixOmin1:ixOmax1,&
       phys_ind%mom(1):phys_ind%mom(ndir))
   if(z_in) then
     where(Dist(ixOmin1:ixOmax1)>0.0_dp)
      sin_theta(ixOmin1:ixOmax1)=x(ixOmin1:ixOmax1,r_)/Dist(ixOmin1:ixOmax1)
      cos_theta(ixOmin1:ixOmax1)=x(ixOmin1:ixOmax1,z_)/Dist(ixOmin1:ixOmax1)
     elsewhere
      sin_theta(ixOmin1:ixOmax1)=1.0_DP
      cos_theta(ixOmin1:ixOmax1)=0.0_dp
     endwhere
   else
    sin_theta(ixOmin1:ixOmax1)=1.0_DP
    cos_theta(ixOmin1:ixOmax1)=0.0_dp
   end if
   if(phi_in) then
     sin_phi(ixOmin1:ixOmax1)=dsin(x(ixOmin1:ixOmax1,phi_))
     cos_phi(ixOmin1:ixOmax1)=dcos(x(ixOmin1:ixOmax1,phi_))
   else
    sin_phi(ixOmin1:ixOmax1)=0.0_DP
    cos_phi(ixOmin1:ixOmax1)=1.0_dp
   end if

   where(self%patch(ixOmin1:ixOmax1))
    w(ixOmin1:ixOmax1,phys_ind%mom(r_)) = v_spherical(ixOmin1:ixOmax1,&
       r_)*cos_phi(ixOmin1:ixOmax1)*sin_theta(ixOmin1:ixOmax1)

    where(cos_phi(ixOmin1:ixOmax1)>0.0_dp.and.self%patch(ixOmin1:ixOmax1))
     w(ixOmin1:ixOmax1,phys_ind%mom(z_)) = v_spherical(ixOmin1:ixOmax1,&
        r_)*cos_theta(ixOmin1:ixOmax1)*cos_phi(ixOmin1:ixOmax1)
    else where(cos_phi(ixOmin1:ixOmax1)<=0.0_dp.and.&
       self%patch(ixOmin1:ixOmax1))
     w(ixOmin1:ixOmax1,phys_ind%mom(z_)) = zero
    end where
   end where
  case('spherical')
  ! Dummy
  case default
  write(*,*) ' is not implimented '
  call mpistop(' stop at mod_usr.t  et usr_supernovae_remnant_spd_rad_to_cart')

 end select
end subroutine usr_supernovae_remnant_spd_rad_to_cart


!--------------------------------------------------------------------
!> Subroutine to clean array memory of associated with cloud object
subroutine usr_supernovae_remnant_clean_memory(self)
  class(supernovae_remnant)    :: self

  if(allocated(self%patch_escape))deallocate(self%patch_escape)
  if(allocated(self%patch))deallocate(self%patch)
  if(self%myconfig%dust_on)call self%mydust%clean_memory
end subroutine usr_supernovae_remnant_clean_memory

end module mod_obj_sn_remnant
