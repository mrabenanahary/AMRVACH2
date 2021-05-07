module mod_obj_rel_jet
  use mod_constants
  use mod_global_parameters
  use mod_srmhd_parameters
  use mod_obj_dust, only : dust
  use mod_obj_global_parameters
  use mod_obj_mat
  implicit none
type rel_jet_parameters
    character(len=20)    :: unit            !> physical unit at parameter file
    real(dp)             :: time_rel_jet_on !> initial time the rel_jet is set in the simulation box

    real(dp)             :: density         !> rel_jet density  (g/cm^3)
    real(dp)             :: number_density  !> rel_jet number density (1/cm^3)
    real(dp)             :: temperature     !> rel_jet temperature  (K)
    real(dp)             :: pressure        !> rel_jet pressure
    real(dp)             :: pressure_toism  !> rel_jet pressure relative to ism pressure
    real(dp)             :: pressure_associate_ism !> rel_jet pressure of associated ism pressure
    real(dp)             :: magnetic(1:3)   !> rel_jet magnetic field components
    real(dp)             :: xisigma         !> rel_jet magnetisation
    real(dp)             :: magn_anglePHItoPol!> rel_jet magnetisation field incl

    real(dp)             :: c_sound         !> rel_jet  sound speed (cm/s)
    real(dp)             :: velocity(1:3)   !> rel_jet  velocity (cm/s)
    real(dp)             :: Lorentz_factor  !> rel_jet Lorentz factor
    real(dp)             :: power           !> wind power flux
    real(dp)             :: mass_flux       !> wind mass flux


    real(dp)             :: open_angle      !> rel_jet initial open angle  (degre)
    real(dp)             :: z_in            !> rel_jet inner boundary position
    real(dp)             :: z_impos         !> rel_jet  impose r
    real(dp)             :: z_out_init      !> rel_jetinitial wind region
    real(dp)             :: r_out_init      !> rel_jet inner boundary wind position
    real(dp)             :: r_impos         !> rel_jet  impose radius

    real(dp)             :: r_in_init      !> initial wind region
    real(dp)             :: center(1:3)     !> rel_jet center position (cm)
    real(dp)             :: extend(1:3)     !> rel_jet region in space (cm)
    logical              :: radial_structure!> rel_jet s
    logical              :: tracer_on       !> rel_jet logical to set tracer
    integer              :: itr             !> rel_jet tracer indice
    character(len=20)    :: shape           !> rel_jet shape
    character(len=20)    :: profile         !> rel_jet profile
end type rel_jet_parameters
! rel_jet features
type rel_jet
  type(rel_jet_parameters)   :: myconfig           !> rel_jet configuration parameters

  logical, allocatable     :: patch(:^D&)        !> spatial patch
  logical, allocatable     :: escape_patch(:^D&) !> spatial patch
  character(len=78)        :: subname            !> subroutine name that call it
  contains
   !PRIVATE
   PROCEDURE, PASS(self) :: set_default     => usr_rel_jet_set_default
   PROCEDURE, PASS(self) :: set_complet     => usr_rel_jet_set_complet
   PROCEDURE, PASS(self) :: normalize       => usr_rel_jet_normalize
   PROCEDURE, PASS(self) :: set_w           => usr_rel_jet_set_w
   PROCEDURE, PASS(self) :: read_parameters => usr_rel_jet_read_p
   PROCEDURE, PASS(self) :: write_setting   => usr_rel_jet_write_setting
   PROCEDURE, PASS(self) :: alloc_set_patch => usr_rel_jet_alloc_set_patch
   PROCEDURE, PASS(self) :: the_patch       => usr_rel_jet_the_patch
   PROCEDURE, PASS(self) :: clean_memory    => usr_rel_jet_clean_memory
end type

contains

  !-------------------------------------------------------------------------
   !> Read the ism parameters  from a parfile
    subroutine usr_rel_jet_read_p(rel_jet_config,files,i_jet,self)
      class(rel_jet)                          :: self
      character(len=*), intent(in)            :: files(:)
      type(rel_jet_parameters), intent(inout) :: rel_jet_config
      integer, optional                       :: i_jet
      ! .. local ..
      integer                                 :: i_file,i_jet_loc
      !----------------------------------------------------------------

      namelist /usr_rel_jet_list/ rel_jet_config
      namelist /usr_rel_jet_1_list/ rel_jet_config
      namelist /usr_rel_jet_2_list/ rel_jet_config
      namelist /usr_rel_jet_3_list/ rel_jet_config


      if(present(i_jet))then
        i_jet_loc = i_jet
      else
        i_jet_loc = 0
      end if
      select case(i_jet_loc)
      case(1)
       do i_file = 1, size(files)
          open(unitpar, file=trim(files(i_file)), status="old")
          read(unitpar, usr_rel_jet_1_list)
          close(unitpar)
       end do
      case(2)
       do i_file = 1, size(files)
          open(unitpar, file=trim(files(i_file)), status="old")
          read(unitpar, usr_rel_jet_2_list)
          close(unitpar)
       end do
      case(3)
       do i_file = 1, size(files)
          open(unitpar, file=trim(files(i_file)), status="old")
          read(unitpar, usr_rel_jet_3_list)
          close(unitpar)
       end do
      case default
       if(mype==0)write(*,*)'Reading usr_rel_jet_list'
       do i_file = 1, size(files)
          open(unitpar, file=trim(files(i_file)), status="old")
          read(unitpar, usr_rel_jet_list)
          close(unitpar)
       end do
     end select

    end subroutine usr_rel_jet_read_p

  !------------------------------------------------------------------------
  !> write the rel_jet setting
  subroutine usr_rel_jet_write_setting(self,unit_config)
    implicit none
    class(rel_jet)                        :: self
    integer,intent(in)                  :: unit_config
    ! .. local ..

    !-----------------------------------

    write(unit_config,*)'************************************'
    write(unit_config,*)'*********rel_jet setting ***********'
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
    write(unit_config,*)'************************************'
    write(unit_config,*)'******* END rel_jet setting ********'
    write(unit_config,*)'************************************'
  end    subroutine usr_rel_jet_write_setting
  !--------------------------------------------------------------------
!> subroutine default setting for rel_jet
 subroutine usr_rel_jet_set_default(self)
  implicit none
  class(rel_jet)          :: self
  !----------------------------------
  self%myconfig%unit                   = 'code'

  self%myconfig%density                = 0.0_dp
  self%myconfig%number_density         = 0.0_dp
  self%myconfig%temperature            = 0.0_dp
  self%myconfig%pressure               = 0.0_dp
  self%myconfig%pressure_toism         = 0.0_dp
  self%myconfig%pressure_associate_ism = 0.0_dp
  self%myconfig%velocity(:)            = 0.0_dp
  self%myconfig%magnetic(:)            = 0.0_dp
  self%myconfig%Lorentz_factor         = 1.0_dp
  self%myconfig%c_sound                = 0.0_dp

  self%myconfig%power                  = 0.0_dp
  self%myconfig%mass_flux              = 0.0_dp



  self%myconfig%z_in                   = 0.0_dp
  self%myconfig%z_out_init             = 0.0_dp
  self%myconfig%z_impos            = 0.0_dp
  self%myconfig%r_out_init                 = 0.0_dp
  self%myconfig%center(:)              = 0.0_dp
  self%myconfig%extend(:)              = 0.0_dp


  self%myconfig%shape                  = 'cylinder'
  self%myconfig%profile                = 'uniform'
  self%myconfig%open_angle             =  0.0_dp
  self%myconfig%radial_structure       = .false.
  self%myconfig%tracer_on              = .false.
  self%myconfig%itr                    = 0
 end subroutine usr_rel_jet_set_default
 !--------------------------------------------------------------------
 !> subroutine check the parfile setting for rel_jet
 subroutine usr_rel_jet_set_complet(self)
   implicit none
   class(rel_jet)             :: self
   ! .. local ..
   real(dp)                   :: mp,kb,beta_speed, jet_surface_init,&
                                 Magnetic_poloidal
   !-----------------------------------
   if(SI_unit) then
     mp=mp_SI
     kB=kB_SI
   else
     mp=mp_cgs
     kB=kB_cgs
   end if


   if(self%myconfig%r_out_init>0.0_dp) then
     if(self%myconfig%open_angle>0.0_dp)then
      self%myconfig%z_in  = self%myconfig%r_out_init/dtan(self%myconfig%open_angle)
     else
      self%myconfig%z_in=xprobmin2/2.0_dp
     end if
   end if


   jet_surface_init = dpi *(self%myconfig%r_out_init**2.0_dp-self%myconfig%r_in_init**2.0_dp)


   if (self%myconfig%Lorentz_factor>0.0_dp) then
     beta_speed = dsqrt(1.0_dp-1.0_dp                             &
                      /self%myconfig%Lorentz_factor**2.0_dp)
     self%myconfig%velocity(r_)     = beta_speed* dsin(self%myconfig%open_angle)


     self%myconfig%velocity(z_)     = beta_speed*dcos(self%myconfig%open_angle)
     self%myconfig%velocity         = unit_velocity*self%myconfig%velocity
   else
    self%myconfig%Lorentz_factor=1.0_dp / &
             dsqrt(1.0_dp-(sum(self%myconfig%velocity)/unit_velocity)**2.0_dp)
   end if

   if(self%myconfig%pressure_toism>0.0_dp.and.self%myconfig%pressure_associate_ism>0.0_dp) then
     self%myconfig%pressure    = self%myconfig%pressure_toism*self%myconfig%pressure_associate_ism
   end if



   if(dabs(self%myconfig%power)>smalldouble) then
     self%myconfig%mass_flux           = self%myconfig%power/unit_velocity**2.0_dp
   end if


   if (dabs(self%myconfig%mass_flux)>smalldouble)then
    self%myconfig%density        = self%myconfig%mass_flux/(jet_surface_init*  &
                                   dabs(self%myconfig%velocity(z_))*self%myconfig%Lorentz_factor**2.0_dp)

    self%myconfig%number_density = self%myconfig%density/mp
    self%myconfig%power               = self%myconfig%mass_flux*unit_velocity
  else if (dabs(self%myconfig%density)<smalldouble*mp)then
    self%myconfig%density        = self%myconfig%number_density*mp
    self%myconfig%mass_flux      = self%myconfig%density*jet_surface_init *  &
                                    dabs(self%myconfig%velocity(z_))*self%myconfig%Lorentz_factor**2.0_dp
    self%myconfig%power          = self%myconfig%mass_flux*unit_velocity**2.0_dp
   else

    if (dabs(self%myconfig%number_density)<smalldouble*mp)then
      self%myconfig%number_density        = self%myconfig%density/mp
    end if

   end if




   if(dabs(self%myconfig%pressure)<=0.0_dp) then
    self%myconfig%pressure =self%myconfig%number_density*&
                                 kB*self%myconfig%temperature
   else
    self%myconfig%temperature = self%myconfig%pressure/(kB*self%myconfig%number_density)
   end if


  self%myconfig%c_sound = sqrt(srmhd_gamma*self%myconfig%pressure/self%myconfig%density)


  if(self%myconfig%xisigma>0.0_dp)then
    if(dabs(dabs(dsin(self%myconfig%magn_anglePHItoPol*dpi/180.0d0))-1)>smalldouble)then
     Magnetic_poloidal =  dsqrt(self%myconfig%Lorentz_factor**2.0_dp* self%myconfig%density*&
                          self%myconfig%xisigma&
                         /(1.0_dp+(dtan(self%myconfig%magn_anglePHItoPol*dpi/180.0_dp))**2.0_dp))


     self%myconfig%magnetic(r_)   = Magnetic_poloidal * dsin(self%myconfig%open_angle*dpi/180.0d0)
     self%myconfig%magnetic(z_)   = Magnetic_poloidal * dcos(self%myconfig%open_angle*dpi/180.0d0)
     self%myconfig%magnetic(phi_) = Magnetic_poloidal * dtan(self%myconfig%magn_anglePHItoPol*dpi/180.0d0)
    else
     self%myconfig%magnetic(r_)    =  0.0_dp
     self%myconfig%magnetic(z_)    =  0.0_dp
     self%myconfig%magnetic(phi_)  = sign(dsqrt(self%myconfig%Lorentz_factor**2.0_dp* self%myconfig%density*&
                          self%myconfig%xisigma),self%myconfig%magn_anglePHItoPol)
    end if
  end if

   if(any(self%myconfig%extend(:)<0.0_dp))then
     self%myconfig%extend(r_)               = self%myconfig%r_out_init
     if(z_in) self%myconfig%extend(theta_)  = self%myconfig%open_angle
     if(phi_in) self%myconfig%extend(phi_)  = 360.0_dp
   end if

   if(.not.self%myconfig%tracer_on)self%myconfig%itr=0

 end subroutine usr_rel_jet_set_complet
!--------------------------------------------------------------------
 subroutine usr_rel_jet_normalize(self,physunit_inuse)
  use mod_obj_usr_unit
  implicit none
  class(rel_jet)                                 :: self
  type(usr_physical_unit_values), intent(in)     :: physunit_inuse
  !----------------------------------
  if(trim(self%myconfig%unit)=='code')return
  self%myconfig%density          = self%myconfig%density       /unit_density
  self%myconfig%temperature      = self%myconfig%temperature   /unit_temperature
  self%myconfig%pressure         = self%myconfig%pressure      /unit_pressure
  self%myconfig%pressure_associate_ism  =self%myconfig%pressure_associate_ism/unit_pressure
  self%myconfig%velocity         = self%myconfig%velocity      /unit_velocity
  self%myconfig%magnetic         = self%myconfig%magnetic      /unit_magneticfield
  self%myconfig%c_sound          =  self%myconfig%c_sound      /unit_velocity

  self%myconfig%open_angle       = self%myconfig%open_angle   *(dpi/180._dp)
  self%myconfig%magn_anglePHItoPol = self%myconfig%magn_anglePHItoPol *(dpi/180._dp)
  self%myconfig%center           = self%myconfig%center        /unit_length
  self%myconfig%extend           = self%myconfig%extend        /unit_length

  self%myconfig%z_in             = self%myconfig%z_in          /unit_length
  self%myconfig%z_impos          = self%myconfig%z_impos       /unit_length
  self%myconfig%z_out_init       = self%myconfig%z_out_init    /unit_length

  self%myconfig%r_out_init       = self%myconfig%r_out_init     /unit_length
  self%myconfig%r_impos          = self%myconfig%r_impos        /unit_length
  self%myconfig%r_in_init        = self%myconfig%r_in_init      /unit_length

  self%myconfig%extend(r_)       = self%myconfig%extend(r_)     /unit_length

  if(z_in) self%myconfig%extend(theta_) = self%myconfig%extend(theta_) *(dpi/180._dp)
  if(phi_in) self%myconfig%extend(phi_) = self%myconfig%extend(phi_) *(dpi/180._dp)


  self%myconfig%power            = self%myconfig%power            /unit_user%luminosity
  self%myconfig%mass_flux        = self%myconfig%mass_flux        /unit_user%mass_flux
  self%myconfig%time_rel_jet_on  = self%myconfig%time_rel_jet_on  / unit_time

  print*, 'pressure in obj_rel_jet', self%myconfig%pressure, self%myconfig%pressure_associate_ism
 end subroutine usr_rel_jet_normalize


!--------------------------------------------------------------------

 !--------------------------------------------------------------------
 !> Subroutine to clean array memory of associated with rel_jet object
 subroutine usr_rel_jet_alloc_set_patch(ixI^L,ixO^L,qt,x,use_tracer,w,&
                                        escape_patch,self)
   implicit none
   integer, intent(in)                     :: ixI^L,ixO^L
   real(kind=dp), intent(in)               :: qt
   real(kind=dp), intent(in)               :: x(ixI^S,1:ndim)
   real(kind=dp), intent(in), optional     :: w(ixI^S,1:nw)
   logical, intent(in), optional           :: use_tracer
   logical, intent(in),optional            :: escape_patch(ixI^S)
   class(rel_jet)                          :: self
   !---------------------------------------------------------
     cond_tracer : if(.not.present(use_tracer).and. .not.present(w)) then
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

     if(present(escape_patch))then
      self%escape_patch(ixO^S) = escape_patch(ixO^S)
     else
      if(allocated(self%escape_patch))deallocate(self%escape_patch)
      allocate(self%escape_patch(ixI^S))
      self%escape_patch(ixO^S) =.false.
     end if
     where(self%patch(ixO^S))self%patch(ixO^S)        = .not.self%escape_patch(ixO^S)
 end subroutine usr_rel_jet_alloc_set_patch
!---------------------------------------------------------------------
!--------------------------------------------------------------------
 !> subroutine patch for the rel_jet
 subroutine usr_rel_jet_the_patch(ixI^L,ixO^L,x,self)
  implicit none
  integer, intent(in)  :: ixI^L,ixO^L
  real(kind=dp)        :: x(ixI^S,1:ndir)
  class(rel_jet)         :: self
  !----------------------------------

  if(.not.allocated(self%patch))allocate(self%patch(ixG^T))

  self%patch              = .false.
  select case(trim(self%myconfig%shape))
  case('conical')
    self%patch(ixO^S)  = dabs(x(ixO^S,r_))<=self%myconfig%r_out_init&
                   +dabs(x(ixO^S,z_)-xprobmin2)*dtan(self%myconfig%open_angle) .and. &
                   dabs(x(ixO^S,r_))>=self%myconfig%r_in_init&
                   +dabs(x(ixO^S,z_)-xprobmin2)*dtan(self%myconfig%open_angle)


  case('cylindrical')
    self%patch(ixO^S)  = dabs(x(ixO^S,r_))<=self%myconfig%r_out_init .and. &
                 dabs(x(ixO^S,r_))>=self%myconfig%r_in_init

  case default
     write(*,*)'this rel_jet shape ',trim(self%myconfig%shape),' is not implimented'
     call mpistop('This rel_jet shape is not implimented in usr_rel_jet_patch at mod_usr.t')
  end select
end subroutine usr_rel_jet_the_patch
!--------------------------------------------------------------------
 !> subroutine setting for rel_jet
 subroutine usr_rel_jet_set_w(ixI^L,ixO^L,qt,x,w,self)
  implicit none
  integer, intent(in)        :: ixI^L,ixO^L
  real(kind=dp), intent(in)  :: qt
  real(kind=dp)              :: x(ixI^S,1:ndir)
  real(kind=dp)              :: w(ixI^S,1:nw)
  class(rel_jet)               :: self
  ! .. local..
  integer                        :: idir
  real(kind=dp)                  :: beta_speed
  real(kind=dp), dimension(ixI^S):: fprofile,angle_theta,&
                                    jet_surface,&
                                    jet_radius_in,jet_radius_out
  !----------------------------------

  if(.not.allocated(self%patch))then
    call self%the_patch(ixI^L,ixO^L,x)
  end if


  cond_inside_rel_jet: if(any(self%patch(ixO^S))) then
   beta_speed = dsqrt(1.0_dp-1.0_dp                             &
                      /self%myconfig%Lorentz_factor**2.0_dp)


   jet_radius_out(ixO^S) =  self%myconfig%r_out_init&
                   +dabs(x(ixO^S,z_))*dtan(self%myconfig%extend(theta_))

   jet_radius_in(ixO^S) =  self%myconfig%r_in_init&
                   +dabs(x(ixO^S,z_))*dtan(self%myconfig%extend(theta_))
   jet_surface(ixO^S)  =  dpi* (jet_radius_out(ixO^S)**2.0_dp-jet_radius_in(ixO^S)**2.0_dp)

    if(dabs(self%myconfig%open_angle)>0.0_dp) then
       where(self%patch(ixO^S))
        angle_theta(ixO^S)  =  datan(x(ixO^S,r_)/(dabs(self%myconfig%z_in+(x(ixO^S,z_)-xprobmin2))))
       end where
    else
       where(self%patch(ixO^S))
        angle_theta(ixO^S)  = 0.0_dp
       end where
    end if

   where(self%patch(ixO^S))
    w(ixO^S,mom(r_))    = beta_speed * dsin(angle_theta(ixO^S))
    w(ixO^S,mom(z_))    = beta_speed * dcos(angle_theta(ixO^S))
    w(ixO^S,mom(phi_))  = self%myconfig%velocity(phi_)

    w(ixO^S,mag(r_))    = self%myconfig%magnetic(r_)
    w(ixO^S,mag(z_))    = self%myconfig%magnetic(z_)
    w(ixO^S,mag(phi_))  = self%myconfig%magnetic(phi_)
   end where

   where(self%patch(ixO^S))
    w(ixO^S,rho_) = self%myconfig%mass_flux/(jet_surface(ixO^S)*  &
                    dabs(w(ixO^S,mom(z_)))*self%myconfig%Lorentz_factor)



    w(ixO^S,p_)   = self%myconfig%pressure -((1.0d0-w(ixO^S,mom(z_))**2.0D0)*&
                    0.5d0*w(ixO^S,mag(phi_))**2.0D0&
                    + w(ixO^S,mag(z_))**2.0D0/2.0D0)
   end where




   cond_tracer_on :if(self%myconfig%tracer_on.and.phys_config%n_tracer>0&
                     .and.self%myconfig%itr<=phys_config%n_tracer)then

    where(self%patch(ixO^S))
     w(ixO^S,phys_ind%tracer(self%myconfig%itr)) = w(ixO^S,phys_ind%rho_)
    elsewhere
     w(ixO^S,phys_ind%tracer(self%myconfig%itr)) = 0.0_dp
    end where
   end if cond_tracer_on


   end if cond_inside_rel_jet



 end subroutine usr_rel_jet_set_w

!--------------------------------------------------------------------
!> Subroutine to clean array memory of associated with rel_jet object
subroutine usr_rel_jet_clean_memory(self)
  class(rel_jet)    :: self
  if(allocated(self%patch))deallocate(self%patch)
end subroutine usr_rel_jet_clean_memory
end module mod_obj_rel_jet
