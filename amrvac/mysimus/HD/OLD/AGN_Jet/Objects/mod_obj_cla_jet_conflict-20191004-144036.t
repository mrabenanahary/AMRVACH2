module mod_obj_cla_jet
  use mod_constants
  use mod_global_parameters
  use mod_physics
  use mod_obj_dust, only : dust
  use mod_obj_global_parameters
  use mod_obj_mat
  use mod_obj_usr_unit
  implicit none
  type cla_jet_parameters
    character(len=20)    :: unit                   !> physical unit at parameter file
    logical              :: normalize_done         !> cla_jet is the normalisation is already done
    real(dp)             :: time_cla_jet_on        !> initial time the cla_jet is set in the simulation box

    integer              :: myindice               !> cla_jet associated indice
    real(dp)             :: density                !> cla_jet density  (g/cm^3)
    real(dp)             :: number_density         !> cla_jet number density (1/cm^3)
    real(dp)             :: temperature            !> cla_jet temperature  (K)
    real(dp)             :: pressure               !> cla_jet pressure
    real(dp)             :: pressure_toism         !> cla_jet pressure relative to ism pressure
    real(dp)             :: pressure_associate_ism !> cla_jet pressure of associated ism pressure
    real(dp)             :: magnetic(1:3)          !> cla_jet magnetic field components
    real(dp)             :: xisigma                !> cla_jet magnetisation
    real(dp)             :: magn_anglePHItoPol     !> cla_jet magnetisation field incl


    real(dp)             :: mach_number            !> cla_jet  mach number
    real(dp)             :: c_sound                !> cla_jet  sound speed (cm/s)
    real(dp)             :: velocity(1:3)          !> cla_jet  velocity (cm/s)
    real(dp)             :: velocity_toroidal      !> cla_jet  toroidal velocity (cm/s)
    real(dp)             :: velocity_poloidal      !> cla_jet  poloidal velocity (cm/s)
    real(dp)             :: power                  !> wind power flux
    real(dp)             :: mass_flux              !> wind mass flux


    real(dp)             :: open_angle             !> cla_jet initial open angle  (degre)
    real(dp)             :: z_in                   !> cla_jet inner boundary position
    real(dp)             :: z_impos                !> cla_jet  impose r
    real(dp)             :: z_out_init             !> cla_jetinitial wind region
    real(dp)             :: r_out_init             !> cla_jet inner boundary wind position
    real(dp)             :: r_out_impos                !> cla_jet  impose radius
    real(dp)             :: r_in_init              !> initial wind region
    real(dp)             :: center(1:3)            !> cla_jet center position (cm)
    real(dp)             :: extend(1:3)            !> cla_jet region in space (cm)

    logical              :: tracer_on              !> cla_jet logical to set tracer
    integer              :: itr                    !> cla_jet tracer indice
    real(dp)             :: tracer_init_density    !> cla_jet tracer initial density
    real(dp)             :: tracer_small_density   !> cla_jet tracer small density cut
    character(len=20)    :: shape                  !> cla_jet shape
    character(len=20)    :: profile                !> cla_jet profile

    integer              :: refine_min_level   !> jet minimum refinent level
    integer              :: refine_max_level   !> jet maximum refinent level
    real(dp)             :: coarsen_distance    !> jet  distance to start coarsen
    real(dp)             :: coarsen_var_distance!> jetscaling distance for coarsen

    character(len=40)    :: variation_type     !> jet type of the power variation

    logical              :: dust_on                !> cla_jet with dust in is true
    real(dp)             :: dust_frac              !> cla_jet dust fraction
    character(len=20)    :: dust_profile           !> cla_jet dust inside profile
  end type cla_jet_parameters
  ! cla_jet features
  type cla_jet
    type(cla_jet_parameters)          :: myconfig           !> cla_jet configuration parameters
    type (dust)                       :: mydust             !> cla_jet dust
    type(usrphysical_unit), pointer   :: myphysunit         !> cla_jet physics unit in use
    logical, allocatable              :: patch(:^D&)        !> spatial patch
    logical, allocatable              :: escape_patch(:^D&) !> spatial patch
    character(len=78)                 :: subname            !> subroutine name that call it

    contains
     !PRIVATE
     PROCEDURE, PASS(self) :: set_default     => usr_cla_jet_set_default
     PROCEDURE, PASS(self) :: set_complet     => usr_cla_jet_set_complet
     PROCEDURE, PASS(self) :: normalize       => usr_cla_jet_normalize
     PROCEDURE, PASS(self) :: set_w           => usr_cla_jet_set_w
     PROCEDURE, PASS(self) :: read_parameters => usr_cla_jet_read_p
     PROCEDURE, PASS(self) :: write_setting   => usr_cla_jet_write_setting
     PROCEDURE, PASS(self) :: alloc_set_patch => usr_cla_jet_alloc_set_patch
     PROCEDURE, PASS(self) :: set_patch       => usr_cla_jet_set_patch
     PROCEDURE, PASS(self) :: clean_memory    => usr_cla_jet_clean_memory
  end type
  integer, save            ::  zjet_,thetajet_
contains

  !-------------------------------------------------------------------------
   !> Read the ism parameters  from a parfile
    subroutine usr_cla_jet_read_p(files,cla_jet_config,self)
      class(cla_jet)                          :: self
      character(len=*), intent(in)            :: files(:)
      type(cla_jet_parameters), intent(inout) :: cla_jet_config
      ! .. local ..
      integer                                 :: i_file
      !----------------------------------------------------------------
      namelist /usr_cla_jet_list/ cla_jet_config
      namelist /usr_cla_jet_1_list/ cla_jet_config
      namelist /usr_cla_jet_2_list/ cla_jet_config
      namelist /usr_cla_jet_3_list/ cla_jet_config



      select case(cla_jet_config%myindice)
      case(1)
       do i_file = 1, size(files)
          open(unitpar, file=trim(files(i_file)), status="old")
          read(unitpar, usr_cla_jet_1_list)
          close(unitpar)
       end do
      case(2)
       do i_file = 1, size(files)
          open(unitpar, file=trim(files(i_file)), status="old")
          read(unitpar, usr_cla_jet_2_list)
          close(unitpar)
       end do
      case(3)
       do i_file = 1, size(files)
          open(unitpar, file=trim(files(i_file)), status="old")
          read(unitpar, usr_cla_jet_3_list)
          close(unitpar)
       end do
      case default
       if(mype==0)write(*,*)'Reading usr_cla_jet_list'
       do i_file = 1, size(files)
          open(unitpar, file=trim(files(i_file)), status="old")
          read(unitpar, usr_cla_jet_list)
          close(unitpar)
       end do


     end select
     if(self%myconfig%dust_on)then
        self%mydust%myconfig%associated_medium = 'cloud'
        call self%mydust%read_parameters(self%mydust%myconfig,files)
    end if
    end subroutine usr_cla_jet_read_p

  !------------------------------------------------------------------------
  !> write the cla_jet setting
  subroutine usr_cla_jet_write_setting(self,unit_config)
    implicit none
    class(cla_jet)                        :: self
    integer,intent(in)                  :: unit_config
    ! .. local ..

    !-----------------------------------

    write(unit_config,*)'************************************'
    write(unit_config,*)'**********cla_jet setting **********'
    write(unit_config,*)'************************************'
    write(unit_config,*) 'power       = ', self%myconfig%power
    write(unit_config,*) 'mass flux   = ', self%myconfig%mass_flux
    write(unit_config,*) 'Density     = ', self%myconfig%density
    write(unit_config,*) 'Pressure    = ', self%myconfig%pressure
    write(unit_config,*) 'Temperature = ', self%myconfig%temperature
    write(unit_config,*) 'Speed       = ', self%myconfig%velocity
    write(unit_config,*) 'magnetic    = ', self%myconfig%magnetic
    write(unit_config,*) 'sound speed = ', self%myconfig%c_sound
    write(unit_config,*) 'open angle  = ', self%myconfig%open_angle

    if(self%myconfig%dust_on)  call self%mydust%write_setting(unit_config)
    write(unit_config,*)'************************************'
    write(unit_config,*)'******* END cla_jet setting ********'
    write(unit_config,*)'************************************'
  end    subroutine usr_cla_jet_write_setting
  !--------------------------------------------------------------------
!> subroutine default setting for cla_jet
 subroutine usr_cla_jet_set_default(self)
  implicit none
  class(cla_jet)          :: self
  !----------------------------------
  self%myconfig%unit                   = 'code'
  self%myconfig%myindice               = 0

  self%myconfig%density                = 0.0_dp
  self%myconfig%number_density         = 0.0_dp
  self%myconfig%temperature            = 0.0_dp
  self%myconfig%pressure               = 0.0_dp
  self%myconfig%pressure_toism         = 0.0_dp
  self%myconfig%pressure_associate_ism = 0.0_dp
  self%myconfig%velocity(:)            = 0.0_dp
  self%myconfig%velocity_poloidal      = 0.0_dp
  self%myconfig%velocity_toroidal      = 0.0_dp
  self%myconfig%magnetic(:)            = 0.0_dp
  self%myconfig%mach_number            = 0.0_dp
  self%myconfig%c_sound                = 0.0_dp

  self%myconfig%power                  = 0.0_dp
  self%myconfig%mass_flux              = 0.0_dp



  self%myconfig%z_in                   = 0.0_dp
  self%myconfig%z_out_init             = 0.0_dp
  self%myconfig%z_impos                = 0.0_dp
  self%myconfig%r_out_init             = 0.0_dp
  self%myconfig%center(:)              = 0.0_dp
  self%myconfig%extend(:)              = 0.0_dp


  self%myconfig%shape                  = 'cylinder'
  self%myconfig%profile                = 'uniform'
  self%myconfig%open_angle             =  0.0_dp
  self%myconfig%tracer_on              = .false.
  self%myconfig%itr                    = 0
  self%myconfig%tracer_init_density    = 0.0_dp
  self%myconfig%tracer_small_density   = 0.0_dp

  self%myconfig%refine_min_level       = 0
  self%myconfig%refine_max_level       = 0
  self%myconfig%coarsen_var_distance   = 0.0_dp
  self%myconfig%coarsen_distance       = 0.0_dp
  self%myconfig%variation_type         = 'none'

  self%myconfig%normalize_done         = .false.
  zjet_                                = 2
  thetajet_                            = 2
  call self%mydust%set_default()
 end subroutine usr_cla_jet_set_default
 !--------------------------------------------------------------------
 !> subroutine check the parfile setting for cla_jet
 subroutine usr_cla_jet_set_complet(self)
   implicit none
   class(cla_jet)             :: self
   ! .. local ..
   logical                    :: dust_is_frac
   real(dp)                   :: mp,kb, jet_surface_init,&
                                 Magnetic_poloidal
   !-----------------------------------
   zjet_     = merge(min(zjet_,ndir),ndir+1,ndir>1)
   thetajet_ = zjet_
   if(SI_unit) then
     mp=mp_SI
     kB=kB_SI
   else
     mp=mp_cgs
     kB=kB_cgs
   end if
   cond_vtor_set : if (self%myconfig%velocity_toroidal>0.0) then
     self%myconfig%velocity(phi_)   = self%myconfig%velocity_toroidal
   else  cond_vtor_set
     self%myconfig%velocity_toroidal     = self%myconfig%velocity(phi_)
   end if cond_vtor_set

   cond_vpol_set : if (self%myconfig%velocity_poloidal>0.0) then
      self%myconfig%velocity(r_)     = self%myconfig%velocity_poloidal*dsin(self%myconfig%open_angle)
      self%myconfig%velocity(theta_) = self%myconfig%velocity_poloidal*dcos(self%myconfig%open_angle)
      self%myconfig%velocity(phi_)   = 0.0_dp

   else cond_vpol_set
      self%myconfig%velocity_poloidal=dsqrt(self%myconfig%velocity(r_)**2.0_dp&
                                           +self%myconfig%velocity(theta_)**2.0_dp)
   end if cond_vpol_set

   if(self%myconfig%r_out_init>0.0_dp) then
     if(self%myconfig%open_angle>0.0_dp)then
      self%myconfig%z_in  = self%myconfig%r_out_init/dtan(self%myconfig%open_angle)
     else
      self%myconfig%z_in=min(xprobmin2/2.0_dp,self%myconfig%z_in)
     end if
   end if

   select case(typeaxial)
   case('cylindrical','spherical')
    jet_surface_init = dpi *(self%myconfig%r_out_init**2.0_dp&
                             -max(self%myconfig%r_in_init,0.0_dp)**2.0_dp)
   case('slab')
      if(self%myconfig%r_out_init*self%myconfig%r_in_init>0.0_dp&
       .or.self%myconfig%r_in_init>xprobmin1) then
        jet_surface_init = dpi *(self%myconfig%r_out_init**2.0_dp&
                             -self%myconfig%r_in_init**2.0_dp)
      else
        jet_surface_init = dpi *(self%myconfig%r_out_init-self%myconfig%r_in_init)**2.0_dp/2.0_dp
      end if
   end select








   if(dabs(self%myconfig%power)>smalldouble) then
     self%myconfig%mass_flux           = self%myconfig%power/unit_velocity**2.0_dp
   end if


   if (dabs(self%myconfig%mass_flux)>smalldouble)then
    self%myconfig%density        = self%myconfig%mass_flux/(jet_surface_init*  &
                                   dabs(self%myconfig%velocity(zjet_)))

    self%myconfig%number_density = self%myconfig%density/mp
    self%myconfig%power               = self%myconfig%mass_flux*unit_velocity
  else
   if (dabs(self%myconfig%density)<smalldouble*mp)then
    self%myconfig%density        = self%myconfig%number_density*mp

   else if (dabs(self%myconfig%number_density)<smalldouble*mp)then
      self%myconfig%number_density        = self%myconfig%density/mp
   end if
   self%myconfig%mass_flux      = self%myconfig%density*jet_surface_init *  &
                                    dabs(self%myconfig%velocity(zjet_))
   self%myconfig%power          = self%myconfig%mass_flux*unit_velocity**2.0_dp
   end if




   if(self%myconfig%pressure_toism>0.0_dp.and.self%myconfig%pressure_associate_ism>0.0_dp) then
     self%myconfig%pressure    = self%myconfig%pressure_toism*self%myconfig%pressure_associate_ism
   end if
   cond_Mach_set : if(self%myconfig%mach_number>0.0_dp) then
     self%myconfig%c_sound = dsqrt(sum(self%myconfig%velocity**2.0_dp))/&
       self%myconfig%mach_number
   else if(self%myconfig%c_sound>0.0_dp) then
     self%myconfig%mach_number = dsqrt(sum(self%myconfig%velocity**2.0_dp))/self%myconfig%c_sound
   end if cond_Mach_set
   cond_csound_set : if(self%myconfig%c_sound>0.0_dp) then
      self%myconfig%pressure = self%myconfig%c_sound**2.0_dp * self%myconfig%density /&
                               phys_config%gamma
   else  cond_csound_set

    if(dabs(self%myconfig%pressure)<=0.0_dp) then
    self%myconfig%pressure =self%myconfig%number_density*&
                                 kB*self%myconfig%temperature
    else
    self%myconfig%temperature = self%myconfig%pressure/(kB*self%myconfig%number_density)
    end if
    self%myconfig%c_sound = sqrt(phys_config%gamma*self%myconfig%pressure/self%myconfig%density)
   end if cond_csound_set




  if(self%myconfig%xisigma>0.0_dp)then
    if(dabs(dabs(dsin(self%myconfig%magn_anglePHItoPol*dpi/180.0d0))-1)>smalldouble)then
     Magnetic_poloidal =  dsqrt( self%myconfig%density*&
                          self%myconfig%xisigma&
                         /(1.0_dp+(dtan(self%myconfig%magn_anglePHItoPol*dpi/180.0_dp))**2.0_dp))


     self%myconfig%magnetic(r_)   = Magnetic_poloidal * dsin(self%myconfig%open_angle*dpi/180.0d0)
     self%myconfig%magnetic(zjet_)   = Magnetic_poloidal * dcos(self%myconfig%open_angle*dpi/180.0d0)
     self%myconfig%magnetic(phi_) = Magnetic_poloidal * dtan(self%myconfig%magn_anglePHItoPol*dpi/180.0d0)
    else
     self%myconfig%magnetic(r_)    =  0.0_dp
     self%myconfig%magnetic(zjet_)    =  0.0_dp
     self%myconfig%magnetic(phi_)  = sign(dsqrt( self%myconfig%density*&
                          self%myconfig%xisigma),self%myconfig%magn_anglePHItoPol)
    end if
  end if

   if(any(self%myconfig%extend(:)<=0.0_dp))then
     self%myconfig%extend(r_)                                 = self%myconfig%r_out_init
     if(z_in.or.zjet_<=ndim) self%myconfig%extend(thetajet_)  = self%myconfig%open_angle
     if(phi_in) self%myconfig%extend(phi_)                    = 360.0_dp
   end if



 if(self%myconfig%refine_min_level==0)then
    self%myconfig%refine_min_level=1
  end if
  if(self%myconfig%refine_max_level==0)then
    self%myconfig%refine_max_level=refine_max_level
  end if
  if(dabs(self%myconfig%coarsen_var_distance )<= 0.0_dp)then
   self%myconfig%coarsen_var_distance=HUGE(0.0_dp)
  end if
  if(dabs(self%myconfig%coarsen_distance )<= 0.0_dp)then
   self%myconfig%coarsen_distance=HUGE(0.0_dp)
  end if




   cond_duston : if(self%myconfig%dust_on)then
     dust_is_frac = .false.
     call self%mydust%set_complet(dust_is_frac,self%myconfig%dust_frac,self%myconfig%density&
               ,self%myconfig%velocity)
   end if cond_duston

   cond_traceron : if(self%myconfig%tracer_on)then
     prim_wnames(self%myconfig%itr+(phys_ind%tracer(1)-1)) = 'tracer_jet'
     cons_wnames(self%myconfig%itr+(phys_ind%tracer(1)-1)) = 'tracer_jet'
   else cond_traceron ! OFF NOW
     self%myconfig%itr=0
   end if cond_traceron


 end subroutine usr_cla_jet_set_complet
!--------------------------------------------------------------------
 subroutine usr_cla_jet_normalize(self,physunit_inuse)
  use mod_obj_usr_unit
  implicit none
  class(cla_jet)                                 :: self
  type(usrphysical_unit), target, intent(in)     :: physunit_inuse
  !----------------------------------

  self%myphysunit =>physunit_inuse
  if(trim(self%myconfig%unit)=='code'.or.self%myconfig%normalize_done)then
     if(self%myconfig%normalize_done)then
      write(*,*) 'WARNING: Second call for cloud normalisation', &
                   'no new normalisation will be done'
     end if
     return
  end if

  self%myconfig%density          = self%myconfig%density       /physunit_inuse%myconfig%density
  self%myconfig%temperature      = self%myconfig%temperature   /physunit_inuse%myconfig%temperature
  self%myconfig%pressure         = self%myconfig%pressure      /physunit_inuse%myconfig%pressure
  self%myconfig%pressure_associate_ism  =self%myconfig%pressure_associate_ism/physunit_inuse%myconfig%pressure
  self%myconfig%velocity    = self%myconfig%velocity /physunit_inuse%myconfig%velocity
  self%myconfig%magnetic         = self%myconfig%magnetic      /physunit_inuse%myconfig%magnetic
  self%myconfig%c_sound          =  self%myconfig%c_sound      /physunit_inuse%myconfig%velocity

  self%myconfig%open_angle       = self%myconfig%open_angle   *(dpi/180._dp)
  self%myconfig%magn_anglePHItoPol = self%myconfig%magn_anglePHItoPol *(dpi/180._dp)
  self%myconfig%center           = self%myconfig%center        /physunit_inuse%myconfig%length


  self%myconfig%z_in             = self%myconfig%z_in          /physunit_inuse%myconfig%length
  self%myconfig%z_impos          = self%myconfig%z_impos       /physunit_inuse%myconfig%length
  self%myconfig%z_out_init       = self%myconfig%z_out_init    /physunit_inuse%myconfig%length

  self%myconfig%r_out_init       = self%myconfig%r_out_init     /physunit_inuse%myconfig%length
  self%myconfig%r_out_impos      = self%myconfig%r_out_impos        /physunit_inuse%myconfig%length
  self%myconfig%r_in_init        = self%myconfig%r_in_init      /physunit_inuse%myconfig%length

  self%myconfig%extend(r_)       = self%myconfig%extend(r_)     /physunit_inuse%myconfig%length

  if(z_in)then
    select case(typeaxial)
    case('spherical')
     self%myconfig%extend(theta_) = self%myconfig%extend(theta_) *(dpi/180._dp)
    case('slab','cylindrical')
      self%myconfig%extend(zjet_) = self%myconfig%extend(zjet_)/physunit_inuse%myconfig%length
    end select
  end if
  if(phi_in)then
    select case(typeaxial)
    case('spherical','cylindrical')
       self%myconfig%extend(phi_) = self%myconfig%extend(phi_) *(dpi/180._dp)
    case('slab')
      self%myconfig%extend(phi_) = self%myconfig%extend(phi_)/physunit_inuse%myconfig%length
    end select
  end if

  self%myconfig%power            = self%myconfig%power            /unit_user%luminosity
  self%myconfig%mass_flux        = self%myconfig%mass_flux        /unit_user%mass_flux
  self%myconfig%time_cla_jet_on  = self%myconfig%time_cla_jet_on  / unit_time


self%myconfig%coarsen_distance     = self%myconfig%coarsen_distance /physunit_inuse%myconfig%length
  self%myconfig%coarsen_var_distance = self%myconfig%coarsen_var_distance /physunit_inuse%myconfig%length


  if(self%myconfig%dust_on)then
    call self%mydust%normalize(physunit_inuse)
    call self%mydust%to_phys()
  end if
  self%myconfig%normalize_done=.true.

 end subroutine usr_cla_jet_normalize


!--------------------------------------------------------------------

 !--------------------------------------------------------------------
 !> Subroutine to clean array memory of associated with cla_jet object
 subroutine usr_cla_jet_alloc_set_patch(ixI^L,ixO^L,qt,x,use_tracer,w,&
                                        escape_patch,self)
   implicit none
   integer, intent(in)                     :: ixI^L,ixO^L
   real(kind=dp), intent(in)               :: qt
   real(kind=dp), intent(in)               :: x(ixI^S,1:ndim)
   real(kind=dp), intent(in), optional     :: w(ixI^S,1:nw)
   logical, intent(in), optional           :: use_tracer
   logical, intent(in),optional            :: escape_patch(ixI^S)
   class(cla_jet)                          :: self
   ! .. local ..
   logical                                 :: usr_tracer_loc
   !---------------------------------------------------------
     if(present(use_tracer).and. .not.present(w)) then
       usr_tracer_loc = use_tracer
     else
       usr_tracer_loc = .false.
     end if
     cond_tracer : if(usr_tracer_loc) then
       if(allocated(self%patch))deallocate(self%patch)
       allocate(self%patch(ixI^S))
       self%patch(ixO^S) =.true.

     else cond_tracer
       cond_jettracer_on : if(self%myconfig%tracer_on)then
        if(allocated(self%patch))deallocate(self%patch)
        allocate(self%patch(ixI^S))
        where(w(ixO^S,phys_ind%tracer(self%myconfig%itr))>small_density)
          self%patch(ixO^S)=.true.
        else where
          self%patch(ixO^S)=.false.
        end where
      end if cond_jettracer_on

     end if cond_tracer

     if(allocated(self%escape_patch))deallocate(self%escape_patch)
     allocate(self%escape_patch(ixI^S))
     if(present(escape_patch))then
      self%escape_patch(ixO^S) = escape_patch(ixO^S)
     else
      self%escape_patch(ixO^S) =.false.
     end if
     self%patch(ixO^S)        = self%patch(ixO^S) .and. (.not.self%escape_patch(ixO^S))
 end subroutine usr_cla_jet_alloc_set_patch
!---------------------------------------------------------------------
!--------------------------------------------------------------------
 !> subroutine patch for the cla_jet
 subroutine usr_cla_jet_set_patch(ixI^L,ixO^L,qt,x,self,force_refine,dx_loc,&
                                  r_limit)
  implicit none
  integer, intent(in)       :: ixI^L,ixO^L
  real(kind=dp), intent(in) :: x(ixI^S,1:ndim)
  real(kind=dp), intent(in) :: qt
  class(cla_jet)            :: self
  integer, optional          :: force_refine
  real(kind=dp),optional     :: dx_loc(1:ndim)
  real(kind=dp),optional     :: r_limit(2)
  ! .. local ..
  integer                    :: idims
  real(dp)                   :: x_edge(ixI^S,1:ndim)
  real(dp), dimension(ixI^S) :: dist_edge
  real(dp)                   :: r_in,r_out,z_in,z_out,min_dist
  !----------------------------------
  ! check distance
  r_in =self%myconfig%r_in_init
  z_in =self%myconfig%z_in
  if(dabs(qt-self%myconfig%time_cla_jet_on)<=smalldouble)then
    r_out=self%myconfig%r_out_init
    z_out=self%myconfig%z_out_init
  else
    r_out=self%myconfig%r_out_impos
    z_out=self%myconfig%z_impos
  end if

  if(present(r_limit))then
   r_limit(1) = r_in
   r_limit(2) = r_out
  end if

  cond_from_refine : if(present(force_refine)) then
    !if(force_refine==-1)r_in=r_out
   Loop_idim : do idims=1,ndim
    x_edge(ixO^S,idims) = x(ixO^S,idims) - dx_loc(idims)/2.0_dp
   end do Loop_idim
    Dist_edge(ixO^S) = x_edge(ixO^S,1)

   need_refine: if(minval(Dist_edge(ixO^S)-r_out)&
                   *maxval(x(ixO^S,1)-r_out)<=smalldouble.or. &
                   minval(Dist_edge(ixO^S)-r_in)&
                   *maxval(x(ixO^S,1)-r_in)<=smalldouble) then
     min_dist=minval(Dist_edge(ixO^S),Dist_edge(ixO^S)>smalldouble)
     r_out=r_out+min_dist
     r_in=r_in-min_dist
   end if  need_refine
  end if  cond_from_refine

  if(.not.allocated(self%patch))allocate(self%patch(ixI^S))

  self%patch              = .false.
  select case(trim(self%myconfig%shape))
  case('conical')
    self%patch(ixO^S)  = dabs(x(ixO^S,r_))<=r_out&
                   +dabs(x(ixO^S,zjet_)-xprobmin2)*dtan(self%myconfig%open_angle) .and. &
                   dabs(x(ixO^S,r_))>=r_in&
                   +dabs(x(ixO^S,zjet_)-xprobmin2)*dtan(self%myconfig%open_angle)

    if(zjet_<=ndim) then
       self%patch(ixO^S)  = self%patch(ixO^S).and. x(ixO^S,zjet_)<=z_out .and. &
                         x(ixO^S,zjet_)>=z_in
    end if

  case('cylindrical')
    self%patch(ixO^S)  = dabs(x(ixO^S,r_))<=r_out .and. &
                         dabs(x(ixO^S,r_))>=r_in
    if(zjet_<=ndim) then
       self%patch(ixO^S)  = self%patch(ixO^S).and. x(ixO^S,zjet_)<=z_out .and. &
                         x(ixO^S,zjet_)>=z_in
    end if

  case('cartesian')
    self%patch(ixO^S)  = x(ixO^S,x_)<=r_out .and. &
                         x(ixO^S,x_)>=r_in
    if(zjet_<=ndim) then
       self%patch(ixO^S)  = self%patch(ixO^S).and. x(ixO^S,zjet_)<=z_out .and. &
                         x(ixO^S,zjet_)>=z_in
    end if

  case default
     write(*,*)'this cla_jet shape ',trim(self%myconfig%shape),' is not implimented'
     call mpistop('This cla_jet shape is not implimented in usr_cla_jet_patch at mod_usr.t')
  end select
end subroutine usr_cla_jet_set_patch
!--------------------------------------------------------------------
 !> subroutine setting for cla_jet
 subroutine usr_cla_jet_set_w(ixI^L,ixO^L,qt,x,w,self)
  implicit none
  integer, intent(in)        :: ixI^L,ixO^L
  real(kind=dp), intent(in)  :: qt
  real(kind=dp)              :: x(ixI^S,1:ndim)
  real(kind=dp)              :: w(ixI^S,1:nw)
  class(cla_jet)               :: self
  ! .. local..
  integer                        :: idir
  logical                        :: dust_is_frac
  real(kind=dp), dimension(ixI^S):: fprofile,angle_theta,&
                                    jet_surface,&
                                    jet_radius_in,jet_radius_out
  !----------------------------------

  if(.not.allocated(self%patch))then
    call self%set_patch(ixI^L,ixO^L,qt,x)
  end if


  cond_inside_cla_jet: if(any(self%patch(ixO^S))) then

   jet_radius_out(ixO^S) =  self%myconfig%r_out_init&
                   +dabs(x(ixO^S,zjet_))*dtan(self%myconfig%extend(theta_))

   jet_radius_in(ixO^S) =  self%myconfig%r_in_init&
                   +dabs(x(ixO^S,zjet_))*dtan(self%myconfig%extend(theta_))
    select case(typeaxial)
    case('cylindrical','spherical')
     jet_surface(ixO^S)  =  dpi* (jet_radius_out(ixO^S)**2.0_dp-max(jet_radius_in(ixO^S)**2.0_dp,0.0_dp))
    case('slab')
     jet_surface(ixO^S)  =  dpi* ((jet_radius_out(ixO^S)-jet_radius_in(ixO^S))/2.0_dp)**2.0_dp
    end select
    if(dabs(self%myconfig%open_angle)>0.0_dp) then
       where(self%patch(ixO^S))
        angle_theta(ixO^S)  =  datan(x(ixO^S,r_)/(dabs(self%myconfig%z_in+(x(ixO^S,zjet_)-xprobmin2))))
       end where
    else
       where(self%patch(ixO^S))
        angle_theta(ixO^S)  = 0.0_dp
       end where
    end if

   where(self%patch(ixO^S))
    w(ixO^S,phys_ind%mom(r_))    = self%myconfig%velocity(r_)!_poloidal * dsin(angle_theta(ixO^S))
    w(ixO^S,phys_ind%mom(zjet_)) = self%myconfig%velocity(zjet_)!_poloidal * dcos(angle_theta(ixO^S))
    w(ixO^S,phys_ind%mom(phi_))  = self%myconfig%velocity(phi_)
   end where
   where(self%patch(ixO^S))
    w(ixO^S,phys_ind%rho_)        = self%myconfig%mass_flux/(jet_surface(ixO^S)*  &
                                        dabs(w(ixO^S,phys_ind%mom(zjet_))))

    w(ixO^S,phys_ind%pressure_)   = self%myconfig%pressure
   end where

   cond_mhd0 : if(phys_config%ismhd)then
    where(self%patch(ixO^S))
     w(ixO^S,phys_ind%mag(r_))    = self%myconfig%magnetic(r_)
     w(ixO^S,phys_ind%mag(zjet_)) = self%myconfig%magnetic(zjet_)
     w(ixO^S,phys_ind%mag(phi_))  = self%myconfig%magnetic(phi_)
    end where
   end if cond_mhd0



  cond_profile : if(self%myconfig%profile/='none') then

      call usr_mat_profile(ixI^L,ixO^L,typeaxial,self%myconfig%profile,&
                           self%myconfig%center,self%myconfig%extend,&
                            x,fprofile)

      where(self%patch(ixO^S))
       w(ixO^S,phys_ind%rho_) = w(ixO^S,phys_ind%rho_) * fprofile(ixO^S)
       w(ixO^S,phys_ind%pressure_) = w(ixO^S,phys_ind%pressure_) * fprofile(ixO^S)
       w(ixO^S,phys_ind%mom(zjet_)) = w(ixO^S,phys_ind%mom(zjet_)) * fprofile(ixO^S)
      end where
  end if cond_profile

  cond_mhd : if(phys_config%ismhd)then
    where(self%patch(ixO^S))
     w(ixO^S,phys_ind%pressure_)   = w(ixO^S,phys_ind%pressure_) -&
                   ((1.0_dp-w(ixO^S,phys_ind%mom(zjet_))**2.0_dp)*&
                    0.5_dp*w(ixO^S,phys_ind%mag(phi_))**2.0_dp&
                    + w(ixO^S,phys_ind%mag(zjet_))**2.0_dp/2.0_dp)
    end where
   end if cond_mhd


   cond_tracer_on :if(self%myconfig%tracer_on.and.phys_config%n_tracer>0&
                      .and.self%myconfig%itr<=phys_config%n_tracer)then
     cond_jet_on : if(qt< self%myconfig%time_cla_jet_on)then
      where(.not.self%patch(ixO^S))
       w(ixO^S,phys_ind%tracer(self%myconfig%itr)) = 0.0_dp
      end where
     else cond_jet_on
      if(self%myconfig%tracer_init_density>0.0_dp) then
       where(self%patch(ixO^S))
        w(ixO^S,phys_ind%tracer(self%myconfig%itr)) = self%myconfig%tracer_init_density
       end where
      else
       where(self%patch(ixO^S))
        w(ixO^S,phys_ind%tracer(self%myconfig%itr)) =  w(ixO^S,phys_ind%rho_)
       end where
      end if
     end if cond_jet_on
     itr=itr+1
   end if cond_tracer_on

   cond_dust_on : if(self%myconfig%dust_on)then
       if(self%myconfig%dust_profile/='none') then
         call usr_mat_profile(ixI^L,ixO^L,typeaxial,self%myconfig%dust_profile,&
                              self%myconfig%center,self%myconfig%extend,&
                               x,fprofile)
       else
         fprofile = 1.0_dp
       end if
       call self%mydust%set_patch(ixI^L,ixO^L,self%patch)
       self%mydust%myconfig%velocity=self%myconfig%velocity
       dust_is_frac=.false.
       call self%mydust%set_w(ixI^L,ixO^L,qt,dust_is_frac,self%myconfig%dust_frac,fprofile,x,w)
    end if cond_dust_on
   end if cond_inside_cla_jet


 end subroutine usr_cla_jet_set_w

!--------------------------------------------------------------------
!> Subroutine to clean array memory of associated with cla_jet object
subroutine usr_cla_jet_clean_memory(self)
  class(cla_jet)    :: self
  if(allocated(self%patch))deallocate(self%patch)
  if(self%myconfig%dust_on)call self%mydust%clean_memory()
end subroutine usr_cla_jet_clean_memory
end module mod_obj_cla_jet
