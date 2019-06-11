module mod_obj_rel_jet
  use mod_constants
  use mod_global_parameters
  use mod_hd
  use mod_srmhd_parameters, only: mag,lfac_,psi_,xi_
  use mod_obj_dust, only : dust
  use mod_obj_global_parameters
  use mod_obj_mat
  implicit none
type rel_jet_parameters
    character(len=20)    :: unit            !> physical unit at parameter file
    real(dp)             :: time_rel_jet_on !> initial time the rel_jet is set in the simulation box
    real(dp)             :: density         !> rel_jet density  (g/cm^3)
    real(dp)             :: number_density  !> rel_jet number density (1/cm^3)
    real(dp)             :: mass            !> rel_jet mass (g)
    real(dp)             :: temperature     !> rel_jet temperature  (K)
    real(dp)             :: pressure        !> rel_jet pressure
    real(dp)             :: center(3)       !> rel_jet center position (cm)
    real(dp)             :: extend(3)       !> rel_jet region in space (cm)
    real(dp)             :: velocity(3)     !> rel_jet  velocity (cm/s)
    real(dp)             :: Lorentz_factor  !> rel_jet Lorentz factor
    real(dp)             :: open_angle      !> rel_jet initial open angle  (degre)
    logical              :: tracer_on       !> logical to set tracer
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
   PROCEDURE, PASS(self) :: the_patch       => usr_rel_jet_patch
   PROCEDURE, PASS(self) :: clean_memory    => usr_rel_jet_clean_memory
end type

contains

  !-------------------------------------------------------------------------
   !> Read the ism parameters  from a parfile
    subroutine usr_rel_jet_read_p(rel_jet_config,files,self)
      class(rel_jet)                          :: self
      character(len=*), intent(in)          :: files(:)
      type(rel_jet_parameters), intent(inout) :: rel_jet_config
      ! .. local ..
      integer                               :: i_file
      !----------------------------------------------------------------

      namelist /usr_rel_jet_list/ rel_jet_config

      if(mype==0)write(*,*)'Reading usr_rel_jet_list'
      do i_file = 1, size(files)
         open(unitpar, file=trim(files(i_file)), status="old")
         read(unitpar, usr_rel_jet_list, end=112)
  112    close(unitpar)
      end do


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
    write(unit_config,*)'************rel_jet setting ************'
    write(unit_config,*)'************************************'
    write(unit_config,*) 'Density     = ', self%myconfig%density
    write(unit_config,*) 'Pressure    = ', self%myconfig%pressure
    write(unit_config,*) 'Temperature = ', self%myconfig%temperature
    write(unit_config,*) 'Speed       = ', self%myconfig%velocity
    write(unit_config,*)'************************************'
    write(unit_config,*)'******** END rel_jet setting **********'
    write(unit_config,*)'************************************'
  end    subroutine usr_rel_jet_write_setting
  !--------------------------------------------------------------------
!> subroutine default setting for rel_jet
 subroutine usr_rel_jet_set_default(self)
  implicit none
  class(rel_jet)          :: self
  !----------------------------------
  self%myconfig%unit             = 'code'
  self%myconfig%density          = 0.0_dp
  self%myconfig%number_density   = 0.0_dp
  self%myconfig%mass             = 0.0_dp
  self%myconfig%temperature      = 0.0_dp
  self%myconfig%pressure         = 0.0_dp
  self%myconfig%center(:)        = 0.0_dp!(box_limit(2,:)-box_limit(1,:))/2.0_dp
  self%myconfig%extend(:)        = 0.0_dp!(box_limit(2,:)+box_limit(1,:))/2.0_dp
  self%myconfig%velocity(:)      = 0.0_dp
  self%myconfig%shape            = 'cylinder'
  self%myconfig%profile          = 'uniform'
  self%myconfig%open_angle      =  0.0_dp
  self%myconfig%tracer_on        = .false.
  call self%mydust%set_default()
 end subroutine usr_rel_jet_set_default
 !--------------------------------------------------------------------
 !> subroutine check the parfile setting for rel_jet
 subroutine usr_rel_jet_set_complet(self)
   implicit none
   class(rel_jet)             :: self
   ! .. local ..
   real(dp)                 :: mp,kb,rel_jet_volume
   !-----------------------------------
   if(SI_unit) then
     mp=mp_SI
     kB=kB_SI
   else
     mp=mp_cgs
     kB=kB_cgs
   end if
   call usr_get_volume(self%myconfig%extend,self%myconfig%shape, rel_jet_volume)

   if (dabs(self%myconfig%mass)>smalldouble)then
    self%myconfig%density        = self%myconfig%mass*(1.0_dp-self%myconfig%dust_frac)/rel_jet_volume
    self%myconfig%number_density = self%myconfig%density/mp
  else if (dabs(self%myconfig%density)<smalldouble*mp)then
    self%myconfig%density        = self%myconfig%number_density*mp
    if(self%myconfig%dust_frac<1.0_dp)then
       self%myconfig%mass        = self%myconfig%density*rel_jet_volume/(1.0_dp-self%myconfig%dust_frac)
     else
       self%myconfig%mass        = 0.0_dp
     end if
   else
    if (dabs(self%myconfig%number_density)<smalldouble*mp)then
      self%myconfig%number_density        = self%myconfig%density/mp
    end if
    if(self%myconfig%dust_frac<1.0_dp)then
      self%myconfig%mass        = self%myconfig%density*rel_jet_volume/(1.0_dp-self%myconfig%dust_frac)
    else
      self%myconfig%mass        = 0.0_dp
    end if
   end if

   if(dabs(self%myconfig%pressure)<smalldouble) then
    self%myconfig%pressure =(2.d0+3.d0*He_abundance)*self%myconfig%number_density*kB*self%myconfig%temperature
   end if





   if(self%myconfig%dust_on)then
    if(self%myconfig%mass>smalldouble) then
      dust_is_frac=.true.
    else
      dust_is_frac=.false.
    end if
     call self%mydust%set_complet(dust_is_frac,self%myconfig%dust_frac,self%myconfig%density,self%myconfig%velocity)
   end if

 end subroutine usr_rel_jet_set_complet
!--------------------------------------------------------------------
 subroutine usr_rel_jet_normalize(self)
  implicit none
  class(rel_jet)          :: self
  !----------------------------------
  if(trim(self%myconfig%unit)=='code')return
  self%myconfig%density          = self%myconfig%density       /unit_density
  self%myconfig%temperature      = self%myconfig%temperature   /unit_temperature
  self%myconfig%pressure         = self%myconfig%pressure      /unit_pressure
  self%myconfig%velocity         = self%myconfig%velocity      /unit_velocity
  self%myconfig%open_angle       = self%myconfig%open_angle   *(dpi/180._dp)
  self%myconfig%center           = self%myconfig%center        /unit_length
  self%myconfig%extend           = self%myconfig%extend        /unit_length


 end subroutine usr_rel_jet_normalize
!--------------------------------------------------------------------

 !--------------------------------------------------------------------
 !> Subroutine to clean array memory of associated with rel_jet object
 subroutine usr_rel_jet_alloc_set_patch(ixI^L,ixO^L,qt,x,escape_patch,self)
   implicit none
   integer, intent(in)           :: ixI^L,ixO^L
   real(kind=dp), intent(in)     :: qt
   real(kind=dp), intent(in)     :: x(ixI^S,1:ndir)
   logical, intent(in),optional  :: escape_patch(ixI^S)
   class(rel_jet)                  :: self
   !---------------------------------------------------------
   if(allocated(self%patch))deallocate(self%patch)
   allocate(self%patch(ixI^S))
   if(allocated(self%escape_patch))deallocate(self%escape_patch)
   allocate(self%escape_patch(ixI^S))
   if(present(escape_patch))then
     self%escape_patch(ixO^S) = escape_patch(ixO^S)
   else
     self%escape_patch(ixO^S) =.false.
   end if


   call self%the_patch(ixI^L,ixO^L,x)
   self%patch(ixO^S)        = self%patch(ixO^S).and. (.not.self%escape_patch(ixO^S))
 end subroutine usr_rel_jet_alloc_set_patch
!---------------------------------------------------------------------
!--------------------------------------------------------------------
 !> subroutine patch for the rel_jet
 subroutine usr_rel_jet_patch(ixI^L,ixO^L,x,self)
  implicit none
  integer, intent(in)  :: ixI^L,ixO^L
  real(kind=dp)        :: x(ixI^S,1:ndir)
  class(rel_jet)         :: self
  !----------------------------------

  if(.not.allocated(self%patch))allocate(self%patch(ixG^T))

  self%patch              = .false.
  select case(trim(self%myconfig%shape))
  case('sphere')
    call usr_set_patch_sphere(ixI^L,ixO^L,typeaxial,self%myconfig%center,    &
                              self%myconfig%extend,x,self%patch)

  case('cylinder')
    call usr_set_patch_cylinder(ixI^L,ixO^L,typeaxial,self%myconfig%center,  &
                                self%myconfig%extend,x,self%patch)
  case('cube')
    call usr_set_patch_cube(ixI^L,ixO^L,typeaxial,self%myconfig%center,      &
                            self%myconfig%extend,x,self%patch)
  case default
     write(*,*)'this rel_jet shape ',trim(self%myconfig%shape),' is not implimented'
     call mpistop('This rel_jet shape is not implimented in usr_rel_jet_patch at mod_usr.t')
  end select
 end subroutine usr_rel_jet_patch
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
  integer                    :: idir
  real(kind=dp)              :: fprofile(ixI^S)
  !----------------------------------
!  call usr_rel_jet_patch(ixI^L,ixO^L,x,self)
  if(.not.allocated(self%patch))then
    call usr_rel_jet_patch(ixI^L,ixO^L,x,self)
  end if
  cond_inside_rel_jet: if(any(self%patch(ixO^S))) then
   where(self%patch(ixO^S))
    w(ixO^S,rho_) = self%myconfig%density
    w(ixO^S,p_)   = self%myconfig%pressure
   end where
   Loop_idir : do idir=1,ndir
       where(self%patch(ixO^S))
         w(ixO^S,mom(idir)) = self%myconfig%velocity(idir)
       end where
   end do Loop_idir
   if(self%myconfig%profile/='none') then
    call usr_mat_profile(ixI^L,ixO^L,typeaxial,self%myconfig%profile,self%myconfig%center,self%myconfig%extend,&
                         x,fprofile)
    where(self%patch(ixO^S))
      w(ixO^S,rho_)=w(ixO^S,rho_)*fprofile(ixO^S)
      !w(ixO^S,p_)=w(ixO^S,p_)*fprofile(ixO^S)
    end where
   end if


   cond_tracer_on :if(self%myconfig%tracer_on.and.phys_n_tracer>0&
                     .and.itr<=phys_n_tracer)then
    where(self%patch(ixO^S))
     w(ixO^S,tracer(itr)) = w(ixO^S,rho_)
    elsewhere
     w(ixO^S,tracer(itr)) = 0.0_dp
    end where
    itr=itr+1
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
