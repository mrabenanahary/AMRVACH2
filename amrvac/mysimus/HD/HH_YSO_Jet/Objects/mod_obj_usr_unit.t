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


   character(len=30) :: chemical_gas_type
   real(kind=dp)     :: chemical_He_abundance
   real(kind=dp)     :: mean_mass
   real(kind=dp)     :: mean_mup
   real(kind=dp)     :: mean_ne_to_nH
   real(kind=dp)     :: mean_nall_to_nH
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

    write(unit_config,*) 'He abundance  = ', self%myconfig%chemical_He_abundance
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

  self%myconfig%chemical_gas_type       = 'fullyionised'
  self%myconfig%chemical_He_abundance   = 0.1_dp
  call unit_chemical_ionisation(self%myconfig%chemical_He_abundance&
                               ,self%myconfig%chemical_gas_type&
                               ,self%myconfig%mean_nall_to_nH &
                               ,self%myconfig%mean_mass,self%myconfig%mean_mup&
                               , self%myconfig%mean_ne_to_nH, &
                               'usr_physical_unit_set_default')

  !write(*,*) ' mod_obj_usr_unit.t--> usr_physical_unit_set_default'
  !write(*,*) 'self: chemical_composition | xHe | mean_nall_to_nH |',&
  !' mean_mass | mean_ne_to_nH | mean_mup '
  !write(*,*) self%myconfig%chemical_gas_type,&
  !self%myconfig%chemical_He_abundance,&
  !self%myconfig%mean_nall_to_nH,&
  !self%myconfig%mean_mass,&
  !self%myconfig%mean_ne_to_nH,&
  !self%myconfig%mean_mup
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
  self%myunit%energy                  = 'g cm^2/s^2'  ! 'erg'
  self%myunit%energy_density          = 'g/(cm s^2)'
  self%myunit%mass_flux               = 'g/s'
  self%myunit%luminosity              = 'g cm^5/s^3'
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
  self%myunit%luminosity              = 'kg m^5/s^3'
end select
end subroutine usr_physical_unit_set_unit
 !--------------------------------------------------------------------
 !> subroutine check the parfile setting for cloud
 subroutine usr_physical_unit_complet(the_physics_type,self)
   use mod_physics, only: phys_config
   implicit none

   character(len=*)                  :: the_physics_type
   class(usrphysical_unit)           :: self
   ! .. local ..
   real(dp)                          :: mp,kb,miu0
   integer                           :: i_phys,n_phys
   character(len=20)                 :: phys_array(4)
   !-----------------------------------

   if(trim(self%myconfig%unit) == 'code')return
   if(SI_unit) then
     mp=mp_SI
     kB=kB_SI
     miu0=miu0_SI
   else
     mp=mp_cgs
     kB=kB_cgs
     miu0=4.0_dp*dpi
   end if



   if(self%myconfig%chemical_He_abundance>0.0_dp)then
     call unit_chemical_ionisation(self%myconfig%chemical_He_abundance,self%myconfig%chemical_gas_type&
                                  ,self%myconfig%mean_nall_to_nH &
                                  ,self%myconfig%mean_mass,self%myconfig%mean_mup&
                                  , self%myconfig%mean_ne_to_nH,&
                                  'usr_physical_unit_complet')

      else if(dabs(self%myconfig%chemical_He_abundance)<=smalldouble)then
        self%myconfig%mean_nall_to_nH   = 1.0_dp
        self%myconfig%mean_mass         = 1.0_dp
        self%myconfig%mean_ne_to_nH     = 1.0_dp
        self%myconfig%mean_mup          = 1.0_dp
      end if

   !write(*,*) ' mod_obj_usr_unit.t-->usr_physical_unit_complet'
   !write(*,*) 'self: chemical_composition | xHe | mean_nall_to_nH |',&
   !' mean_mass | mean_ne_to_nH | mean_mup '
   !write(*,*) self%myconfig%chemical_gas_type,&
   !self%myconfig%chemical_He_abundance,&
   !self%myconfig%mean_nall_to_nH,&
   !self%myconfig%mean_mass,&
   !self%myconfig%mean_ne_to_nH,&
   !self%myconfig%mean_mup

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
     self%myconfig%density        = self%myconfig%number_density*mp*self%myconfig%mean_mass
   else if(self%myconfig%density>smalldouble)  then
     self%myconfig%number_density = self%myconfig%density/(mp*self%myconfig%mean_mass)
   end if

   if(self%myconfig%temperature>smalldouble)then
    self%myconfig%pressure      = self%myconfig%density&
            *kB*self%myconfig%temperature / (self%myconfig%mean_mup*mp)
   else if(self%myconfig%velocity>smalldouble)then
    self%myconfig%pressure      = self%myconfig%density*self%myconfig%velocity**2.0_dp
    self%myconfig%temperature   = self%myconfig%mean_mup*mp*self%myconfig%pressure&
              / ( Kb*self%myconfig%density)
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
   self%myconfig%luminosity       = self%myconfig%pressure &
                          /(self%myconfig%number_density**2.0_dp*unit_time * self%myconfig%mean_mass**2.0_dp)

    ! in mod_radiative_cooling : unit_luminosity=  unit_pressure/( unit_numberdensity**2.0_dp * unit_time &
    !                            *  phys_config%mean_mass**2.0_dp)

    ! set the physical unit
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

   !HI
   w_convert_factor(phys_ind%HI_density_) = self%myconfig%density
   !HII
   w_convert_factor(phys_ind%HII_density_) = self%myconfig%density
   !HeI
   w_convert_factor(phys_ind%HeI_density_) = self%myconfig%density
   !HeII
   w_convert_factor(phys_ind%HeII_density_) = self%myconfig%density
   !HeIII
   w_convert_factor(phys_ind%HeIII_density_) = self%myconfig%density
   !electrons
   w_convert_factor(phys_ind%e_density_) = self%myconfig%density
   !HM
   w_convert_factor(phys_ind%HM_density_) = self%myconfig%density
   !H2I
   w_convert_factor(phys_ind%H2I_density_) = self%myconfig%density
   !H2II
   w_convert_factor(phys_ind%H2II_density_) = self%myconfig%density
   !DI
   w_convert_factor(phys_ind%DI_density_) = self%myconfig%density
   !DII
   w_convert_factor(phys_ind%DII_density_) = self%myconfig%density
   !HDI
   w_convert_factor(phys_ind%HDI_density_) = self%myconfig%density
   !metal
   w_convert_factor(phys_ind%metal_density_) = self%myconfig%density
   !dust
   w_convert_factor(phys_ind%dust_density_) = self%myconfig%density

   w_convert_factor(phys_ind%rhoX_) = self%myconfig%density
   w_convert_factor(phys_ind%rhoY_) = self%myconfig%density


   !write(*,*) ' mod_obj_usr_unit.t--> usr_physical_unit_fillphysunit'
   !write(*,*) 'phys_config (before change): chemical_composition | ',&
   !'xHe | mean_nall_to_nH |',&
   !' mean_mass | mean_ne_to_nH | mean_mup '
   !write(*,*) phys_config%chemical_gas_type,&
   !phys_config%He_abundance,&
   !phys_config%mean_nall_to_nH,&
   !phys_config%mean_mass,&
   !phys_config%mean_ne_to_nH,&
   !phys_config%mean_mup

   phys_config%He_abundance = self%myconfig%chemical_He_abundance
   phys_config%mean_mass    = self%myconfig%mean_mass
   phys_config%mean_mup     = self%myconfig%mean_mup

   !write(*,*) ' mod_obj_usr_unit.t--> usr_physical_unit_fillphysunit'
   !write(*,*) 'phys_config (after change): chemical_composition | ',&
   !'xHe | mean_nall_to_nH |',&
   !' mean_mass | mean_ne_to_nH | mean_mup '
   !write(*,*) phys_config%chemical_gas_type,&
   !phys_config%He_abundance,&
   !phys_config%mean_nall_to_nH,&
   !phys_config%mean_mass,&
   !phys_config%mean_ne_to_nH,&
   !phys_config%mean_mup

   if (phys_config%radiative_cooling) then

     if(phys_config%cool_savedT)w_convert_factor(phys_ind%dtcool1_) = self%myconfig%time
     if(phys_config%cool_saveL)w_convert_factor(phys_ind%Lcool1_)   =    &
           1.0_dp/(self%myconfig%number_density**2.0_dp * self%myconfig%time &
                 * self%myconfig%mean_mass**2.0_dp /self%myconfig%pressure)
     if(phys_config%mean_mup_on)w_convert_factor(phys_ind%mup_) = self%myconfig%mean_mup
   end if
  if(mype==0)then
    print*, 'phys_config%He_abundance', phys_config%He_abundance
    print*, 'in mod_obj_usr_unit.t usr_physical_unit_fillphysunit subroutine: '
    print*, 'ism density', self%myconfig%density
    print*, 'ism ndensity', self%myconfig%number_density
    print*, 'ism pressure', self%myconfig%pressure
    print*, 'ism temperature', self%myconfig%temperature
  end if
 end subroutine usr_physical_unit_fillphysunit


   subroutine unit_chemical_ionisation(He_abundance_sub,chemical_gas_type,mean_nall_to_nH &
                                        ,mean_mass,mean_mup,mean_ne_to_nH,subroutine_name)
    use mod_global_parameters
    implicit none
    real(kind=dp), intent(in)    :: He_abundance_sub
    character(len=*), intent(in) :: chemical_gas_type
    real(kind=dp), intent(out)   :: mean_nall_to_nH,mean_mass,mean_mup,mean_ne_to_nH
    character(len=*),  optional, intent(in) :: subroutine_name
    character(len=72)             :: subroutine_name_='undefined'
    !----------------------------------------------------------
    if(present(subroutine_name)) then
      subroutine_name_ = subroutine_name
    end if
    if(He_abundance_sub>0.0_dp)then
      mean_nall_to_nH = 1.0_dp+2.0_dp*He_abundance_sub
      mean_mass       = 1.0_dp+4.0_dp*He_abundance_sub
      select case(trim(chemical_gas_type))
        case('fullymolecular')
          !M> meaning that n(H+H^+)\simeq n(e^-)\simeq0 then n_H=0+2n(H_2) ==> n(H_2)=0.5n_H
          !M> and n(He)=He_abundance_sub*n_H, n(He+)=n(He++)=0
          !M> ==> total particle density n_tot = n(H+H^++H_2+He+e^-+He^++He^++)
          !M> = 0+0.5n_H+He_abundance_sub*n_H+0+0+0= (0.5+He_abundance_sub)*n_H
          !M> ==> total to H-formed species particle density ratio mean_mass/mean_mup = mean_mass/mu = n_tot*m_H/(n_H*m_H) = 0.5+He_abundance_sub
          !M> Fully molecular fluid ensures n(e^-)\simeq 0 so that mean_ne_to_nH = 0
          mean_mup = mean_mass/(0.5_dp+He_abundance_sub)
          mean_ne_to_nH =0.0  ! is the value used in implimentation of low temperarure cooling table DM2
        case('fullyatomic')
          !M> meaning that n(H^+)\simeq n(e^-)\simeq 0 and n(H_2)=0 then n_H=n(H)+0 ==> n(H)=1.0n_H
          !M> and n(He)=He_abundance_sub*n_H and n(He^+)=n(He^++)=0
          !M> ==> total particle density n_tot = n(H+H^++H_2+He+e^-+He^++He^++)
          !M> = 1.0n_H+0+0+He_abundance_sub*n_H+0+0+0= (1.0+He_abundance_sub)*n_H
          !M> ==> total to H-formed species particle density ration mean_mass/mean_mup = mean_mass/mu = n_tot*m_H/(n_H*m_H) =  1.0+He_abundance_sub
          !M> Fully atomi fluid ensures n(e^-)\simeq 0 so that mean_ne_to_nH = 0
          mean_mup = mean_mass/(1.0_dp+He_abundance_sub)
          mean_ne_to_nH =0.0 ! is the value used in implimentation of low temperarure cooling table DM2
        case('ionised')
          !M> equivalent to a fluid where He is once ionized,
          !M> meaning that n(H)\simeq n(H_2)\simeq 0 then n_H=n(H^+) ==> n(H^+)=1.0n_H
          !M> and n(He^+)=He_abundance_sub*n_H and n(He)=n(He^++)=0 and n(e^-)=n(H^+)+n(He^+)=1.0n_H+He_abundance_sub*nH
          !M> ==> total particle density n_tot = n(H+H^++H_2+He+e^-+He^++He^++)
          !M> = 0+1.0n_H+0+0+(1.0n_H+He_abundance_sub*n_H)+He_abundance_sub*n_H+0= (2.0+2*He_abundance_sub)*n_H
          !M> ==> total to H-formed species particle density ration mean_mass/mean_mup = mean_mass/mu = n_tot*m_H/(n_H*m_H) =  2.0+2*He_abundance_sub
          !M> once ionized H and He fluid ensures n(e^-)=n(H^+)+n(He+)=1.0n_H+He_abundance_sub*n_H so that mean_ne_to_nH = 1.0+He_abundance_sub
          mean_mup = mean_mass/(2.0_dp+2.0_dp*He_abundance_sub)
          mean_ne_to_nH = 1.0_dp+He_abundance_sub
        case('fullyionised')
          !M> equivalent to a fluid where He is twice ionized,
          !M> meaning that n(H)\simeq n(H_2)\simeq 0 then n_H=n(H^+) ==> n(H^+)=1.0n_H
          !M> and n(He^++)=2*He_abundance_sub*n_H and n(He)=n(He^+)=0 and n(e^-)=n(H^+)+2*n(He^++)=1.0n_H+2*He_abundance_sub*nH
          !M> ==> total particle density n_tot = n(H+H^++H_2+He+e^-+He^++He^++)
          !M> = 0+1.0n_H+0+0+(1.0n_H+2*He_abundance_sub*n_H)+He_abundance_sub*n_H+0= (2.0+3*He_abundance_sub)*n_H
          !M> ==> total to H-formed species particle density ration mean_mass/mean_mup = mean_mass/mu = n_tot*m_H/(n_H*m_H) =  2.0+3*He_abundance_sub
          !M> once ionized H and twice ionized He fluid ensures n(e^-)=n(H^+)+n(He^++)=1.0n_H+2*He_abundance_sub*n_H so that mean_ne_to_nH = 1.0+2*He_abundance_sub
          mean_mup = mean_mass/(2.0_dp+3.0_dp*He_abundance_sub)
          mean_ne_to_nH = 1.0_dp+2.0_dp*He_abundance_sub
       case default
          write(*,*) 'The chemical gas type : ', trim(chemical_gas_type)
          write (*,*) "Undefined gas chemical type entered in mod_obj_usr_unit.t "
          call mpistop('The code stops at mod_obj_usr_unit.t')
      end select
      !hd_config%mean_mup          = (2.0_dp+3.0_dp*He_abundance_sub)

    else if(dabs(He_abundance_sub)<=smalldouble)then
      mean_mass         = 1.0_dp
      mean_mup          = 1.0_dp
      mean_ne_to_nH     = 1.0_dp
    end if

    !write(*,*) ' mod_obj_usr_unit.t--> unit_chemical_ionisation'
    !write(*,*) ' function private: chemical_composition |',&
    !' xHe | mean_nall_to_nH |',&
    !' mean_mass | mean_ne_to_nH | mean_mup | sub_name'
    !write(*,*) chemical_gas_type,&
    !He_abundance_sub,&
    !mean_nall_to_nH,&
    !mean_mass,&
    !mean_ne_to_nH,&
    !mean_mup,&
    !subroutine_name_

  end   subroutine unit_chemical_ionisation
end module mod_obj_usr_unit
