module mod_obj_usr_unit
 use mod_constants
 use mod_global_parameters
 use mod_physics

 type usr_physical_unit_values
   character(len=20) :: unit
   real(kind=dp)     :: length
   real(kind=dp)     :: volume
   real(kind=dp)     :: time

   real(kind=dp)     :: number_density
   real(kind=dp)     :: density
   real(kind=dp)     :: pressure
   real(kind=dp)     :: temperature


   real(kind=dp)     :: velocity
   real(kind=dp)     :: angular_velocity
   real(kind=dp)     :: momentum
   real(kind=dp)     :: magnetic

   real(kind=dp)     :: mass
   real(kind=dp)     :: mass_flux

   real(kind=dp)     :: energy_density
   real(kind=dp)     :: energy
   real(kind=dp)     :: energy_flux
   real(kind=dp)     :: luminosity
 end type usr_physical_unit_values
 type usr_physical_unit_dim
   character(len=10)     :: length
   character(len=10)     :: volume
   character(len=10)     :: time

   character(len=10)     :: number_density
   character(len=10)     :: density
   character(len=10)     :: pressure
   character(len=10)     :: temperature


   character(len=10)     :: velocity
   character(len=10)     :: angular_velocity
   character(len=10)     :: momentum
   character(len=10)     :: magnetic

   character(len=10)     :: mass
   character(len=10)     :: mass_flux

   character(len=10)     :: energy_density
   character(len=10)     :: energy
   character(len=10)     :: energy_flux
   character(len=10)     :: luminosity
 end type usr_physical_unit_dim
 type usrphysical_unit
   logical                          :: is_on
   type(usr_physical_unit_dim)      :: myunit
   type(usr_physical_unit_values)   :: myconfig !> physical unit values
 contains
   PROCEDURE, PASS(self) :: set_default        => usr_physical_unit_set_default
   PROCEDURE, PASS(self) :: set_unit           => usr_physical_unit_set_unit
   PROCEDURE, PASS(self) :: set_complet        => usr_physical_unit_complet
   PROCEDURE, PASS(self) :: read_parameters    => usr_physical_unit_read_p
   PROCEDURE, PASS(self) :: write_setting      => usr_physical_unit_write_setting
   PROCEDURE, PASS(self) :: fillphysunit       => usr_physical_unit_fillphysunit
   PROCEDURE, PASS(self) :: set_to_one         => usr_physical_unit_set_to_one
 end type usrphysical_unit

 contains


!-------------------------------------------------------------------------
 !> Read the ism parameters  from a parfile
  subroutine usr_physical_unit_read_p(self,phys_unit,&
     files)
    class(usrphysical_unit)               :: self
    character(len=*), intent(in)          :: files(:)
    type(usr_physical_unit_values)        :: phys_unit

    integer  :: i_file
    namelist /usr_physical_unit_list/ phys_unit

    if(mype==0)write(*,*)'Reading usr_physical_values_list'
    Loop_ifile: do i_file = 1, size(files)
       open(unitpar, file=trim(files(i_file)), status="old")
       read(unitpar, usr_physical_unit_list)
       close(unitpar)
    end do Loop_ifile
    call self%set_unit
  end subroutine usr_physical_unit_read_p

  !------------------------------------------------------------------------
  !> write the cloud setting
  subroutine usr_physical_unit_write_setting(self,unit_config)
    implicit none
    class(usrphysical_unit)                :: self
    integer,intent(in)                     :: unit_config
    ! .. local ..

    !-----------------------------------

    write(unit_config,*)'************************************'
    write(unit_config,*)'********** Normalisation ********'
    write(unit_config,*)'************************************'
    write(unit_config,*) 'The unit      = ', self%myconfig%unit
    write(unit_config,*) 'Length        = ', self%myconfig%length, '   ', trim(self%myunit%length)
    write(unit_config,*) 'Volume        = ', self%myconfig%volume, '   ',self%myunit%volume
    write(unit_config,*) 'Time          = ', self%myconfig%time, '   ',self%myunit%time
    write(unit_config,*) 'Density       = ', self%myconfig%density,'   ',self%myunit%density
    write(unit_config,*) 'Number Density= ', self%myconfig%number_density,'   ',  self%myunit%number_density
    write(unit_config,*) 'Pressure      = ', self%myconfig%pressure, '   ', self%myunit%pressure
    write(unit_config,*) 'Temperature   = ', self%myconfig%temperature, '   ',  self%myunit%temperature
    write(unit_config,*) 'Energy        = ', self%myconfig%energy, '   ',  self%myunit%energy
    write(unit_config,*) 'Mass          = ', self%myconfig%mass, '   ', self%myunit%mass
    write(unit_config,*) 'Speed         = ', self%myconfig%velocity, '   ', self%myunit%velocity
    write(unit_config,*) 'Angular Speed = ', self%myconfig%angular_velocity,'   ', self%myunit%angular_velocity
    write(unit_config,*) 'Moment        = ', self%myconfig%momentum,'   ', self%myunit%momentum
    write(unit_config,*) 'Magnetic      = ', self%myconfig%magnetic,'   ', self%myunit%magnetic
    write(unit_config,*) 'Mass          = ', self%myconfig%mass_flux,'   ',self%myunit%mass_flux
    write(unit_config,*) 'Luminosity    = ', self%myconfig%luminosity,'   ',self%myunit%luminosity
    write(unit_config,*)'************************************'
    write(unit_config,*)'******* END Normalisation *******'
    write(unit_config,*)'************************************'
  end    subroutine usr_physical_unit_write_setting
  !--------------------------------------------------------------------
!> subroutine default setting for cloud
 subroutine usr_physical_unit_set_default(self)
  implicit none
  class(usrphysical_unit)          :: self
  !----------------------------------
  self%myconfig%unit                    =  'code'
  self%myconfig%density                 = 0.0_DP
  self%myconfig%number_density          = 0.0_DP
  self%myconfig%pressure                = 0.0_dp
  self%myconfig%temperature             = 0.0_dp
  self%myconfig%velocity                = 0.0_dp
  self%myconfig%angular_velocity        = 0.0_dp
  self%myconfig%momentum                = 0.0_dp

  self%myconfig%magnetic                = 0.0_dp

  self%myconfig%mass                    = 0.0_dp
  self%myconfig%length                  = 0.0_dp
  self%myconfig%volume                  = 0.0_dp
  self%myconfig%time                    = 0.0_dp
  self%myconfig%energy                  = 0.0_dp
  self%myconfig%energy_density          = 0.0_dp
  self%myconfig%mass_flux               = 0.0_dp
  self%myconfig%luminosity              = 0.0_dp
 end subroutine usr_physical_unit_set_default
!> subroutine unit to code unit
 subroutine usr_physical_unit_set_to_one(self)
  implicit none
  class(usrphysical_unit)          :: self
  !----------------------------------
  self%myconfig%unit                    =  'code'
  self%myconfig%density                 = 1.0_DP
  self%myconfig%number_density          = 1.0_DP
  self%myconfig%pressure                = 1.0_dp
  self%myconfig%temperature             = 1.0_dp
  self%myconfig%velocity                = 1.0_dp
  self%myconfig%angular_velocity        = 1.0_dp
  self%myconfig%momentum                = 1.0_dp
  self%myconfig%magnetic                = 1.0_dp
  self%myconfig%mass                    = 1.0_dp
  self%myconfig%length                  = 1.0_dp
  self%myconfig%volume                  = 1.0_dp
  self%myconfig%time                    = 1.0_dp
  self%myconfig%energy                  = 1.0_dp
  self%myconfig%energy_density          = 1.0_dp
  self%myconfig%mass_flux               = 1.0_dp
  self%myconfig%luminosity              = 1.0_dp
end subroutine usr_physical_unit_set_to_one
!> subroutine default setting for cloud
 subroutine usr_physical_unit_set_unit(self)
  implicit none
  class(usrphysical_unit)          :: self
  !----------------------------------
  select case(trim(self%myconfig%unit))
  case('cgs')
  self%myunit%density                 = 'g/cm3'
  self%myunit%number_density          = '1/cm^3'
  self%myunit%pressure                = 'g/(s cm)'
  self%myunit%temperature             = 'K'
  self%myunit%velocity                = 'cm/s'
  self%myunit%angular_velocity        = '1/s'
  self%myunit%momentum                = 'g/s'

  self%myunit%magnetic                = 'Gauss'

  self%myunit%mass                    = 'g'
  self%myunit%length                  = 'cm'
  self%myunit%volume                  = 'cm^3'
  self%myunit%time                    = 's'
  self%myunit%energy                  = 'g cm^2/s^2'
  self%myunit%energy_density          = 'g/(cm s^2)'
  self%myunit%mass_flux               = 'g/s'
  self%myunit%luminosity              = 'g cm^2/s^3'
case('mks')
  self%myunit%density                 = 'kg/m3'
  self%myunit%number_density          = '1/m^3'
  self%myunit%pressure                = 'kg/(s m)'
  self%myunit%temperature             = 'K'
  self%myunit%velocity                = 'm/s'
  self%myunit%angular_velocity        = '1/s'
  self%myunit%momentum                = 'kg/s'

  self%myunit%magnetic                = 'Tesla'

  self%myunit%mass                    = 'kg'
  self%myunit%length                  = 'm'
  self%myunit%volume                  = 'm^3'
  self%myunit%time                    = 's'
  self%myunit%energy                  = 'kg m^2/s^2'
  self%myunit%energy_density          = 'kg/(m s^2)'
  self%myunit%mass_flux               = 'kg/s'
  self%myunit%luminosity              = 'kg m^2/s^3'
end select
end subroutine usr_physical_unit_set_unit
 !--------------------------------------------------------------------
 !> subroutine check the parfile setting for cloud
 subroutine usr_physical_unit_complet(the_physics_type,self)
   implicit none
   character(len=*)                  :: the_physics_type
   class(usrphysical_unit)           :: self
   ! .. local ..
   real(dp)                          :: mp,kb,miu0,mean_mass
   integer                           :: i_phys,n_phys
   character(len=20)                 :: phys_array(4)
   !-----------------------------------
   if(trim(self%myconfig%unit) == 'code')return
   if(SI_unit) then
     mp=mp_SI
     kB=kB_SI
     miu0=miu0_SI
     mean_mass = (2.d0)
   else
     mp=mp_cgs
     kB=kB_cgs
     miu0=4.d0*dpi
     mean_mass = (2.d0)
   end if
   phys_array = [character(len=20):: 'srmhd','srhd','grmhd','grhd']
   n_phys=size(phys_array)
   Loop_phys : do i_phys = 1,n_phys
    if(index(trim(the_physics_type),trim(phys_array(i_phys)))>0)then
      self%myconfig%velocity = const_c
    end if
   end do Loop_phys

   if(self%myconfig%length>smalldouble) then
     if(self%myconfig%velocity>smalldouble) then
      self%myconfig%time = self%myconfig%length/dabs(self%myconfig%velocity)
     else if(self%myconfig%time>smalldouble)then
      self%myconfig%velocity = self%myconfig%length/self%myconfig%time
     end if
   elseif(self%myconfig%time>smalldouble)then
      self%myconfig%length     = self%myconfig%time*dabs(self%myconfig%velocity)
   end if

   if(self%myconfig%number_density>smalldouble)then
     self%myconfig%density        = self%myconfig%number_density*mp
   else if(self%myconfig%density>smalldouble)  then
     self%myconfig%number_density = self%myconfig%density/mp
   end if

   if(self%myconfig%temperature>smalldouble)then
    self%myconfig%pressure      = mean_mass*self%myconfig%number_density&
                                  *kB*self%myconfig%temperature
   else
    self%myconfig%pressure      = self%myconfig%density*self%myconfig%velocity**2.0_dp
    self%myconfig%temperature   = self%myconfig%pressure&
                                 / (mean_mass*Kb*self%myconfig%number_density)
   end if

   self%myconfig%momentum         = self%myconfig%density*self%myconfig%velocity
   self%myconfig%angular_velocity = self%myconfig%velocity/self%myconfig%length

   self%myconfig%volume           = self%myconfig%length**3.0_dp
   self%myconfig%mass             = self%myconfig%volume*self%myconfig%density

   self%myconfig%energy           = self%myconfig%mass*self%myconfig%velocity**2.0_dp
   self%myconfig%energy_density   = self%myconfig%density*self%myconfig%velocity**2.0_dp

   self%myconfig%magnetic         = dsqrt(self%myconfig%energy_density)

   self%myconfig%mass_flux        = self%myconfig%mass/self%myconfig%time
   self%myconfig%energy_flux      = self%myconfig%energy/self%myconfig%time
   self%myconfig%luminosity       = self%myconfig%energy/self%myconfig%time


     unit_length        = self%myconfig%length
     unit_time          = self%myconfig%time

     unit_velocity      = self%myconfig%velocity
     unit_density       = self%myconfig%density
     unit_numberdensity = self%myconfig%number_density
     unit_pressure      = self%myconfig%pressure
     unit_temperature   = self%myconfig%temperature
 end subroutine usr_physical_unit_complet

 subroutine usr_physical_unit_fillphysunit(self)
   class(usrphysical_unit)    :: self
   ! .. local ..
   integer :: idust
   !---------------------------------

   w_convert_factor =1.0_dp


   w_convert_factor(phys_ind%rho_)                  = self%myconfig%density
   if(phys_config%energy)then
     w_convert_factor(phys_ind%e_)                  = self%myconfig%energy_density
   end if
   if(saveprim)then
    w_convert_factor(phys_ind%mom(:))               = self%myconfig%velocity
   else
    w_convert_factor(phys_ind%mom(:))               = self%myconfig%momentum
   end if


   time_convert_factor                              = self%myconfig%time
   length_convert_factor                            = self%myconfig%length


   if(phys_config%dust_on)then
      w_convert_factor(phys_ind%dust_rho(:))        = self%myconfig%density
      if(saveprim)then
       do idust = 1, phys_config%dust_n_species
        w_convert_factor(phys_ind%dust_mom(:,idust))= self%myconfig%velocity
       end do
      else
       do idust = 1, phys_config%dust_n_species
        w_convert_factor(phys_ind%dust_mom(:,idust)) = self%myconfig%momentum
       end do
      end if
   end if
 end subroutine usr_physical_unit_fillphysunit
end module mod_obj_usr_unit
