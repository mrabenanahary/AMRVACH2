module mod_obj_pulsar
use mod_constants
use mod_global_parameters
use mod_obj_global_parameters
use mod_obj_mat
use mod_physics
use mod_srmhd_parameters
use mod_obj_star_envelope
use mod_obj_sn_remnant
use mod_obj_relativistic_wind
use mod_obj_star
implicit none

  type pulsar_parametres
    character(len=20)    :: unit                 !> physical unit at parameter file
    real(dp)             :: velocity(1:3)         !> moving speed of the pulsar
    real(dp)             :: center(1:3)           !> moving speed of the pulsar
    real(dp)             :: radius                !> pulsar radius
    real(dp)             :: omega                 !> pulsar rotation speed
    real(dp)             :: period                !> pulsar period
    real(dp)             :: mass                  !> pulsar mass
    logical              :: wind_on               !> pulsar to set wind
    logical              :: envelope_on           !> pulsar to set envelope
    logical              :: star_on               !> pulsar to set star
    logical              :: supernovae_remnant_on !> pulsar to set star
    real(dp)             :: t_start_pulsar_wind   !> pulsar starting wind time
    real(dp)             :: t_end_pulsar_wind     !> pulsar starting wind time
    real(dp)             :: t_start_star_wind     !> pulsar starting stellar wind time
    real(dp)             :: t_end_star_wind       !> pulsar starting stellar wind time
    real(dp)             :: t_sn                  !> supernovae time
  end  type

  ! the pulsar type
  type pulsar
    logical, allocatable    :: patch(:^D&)           !> pulsar is on cell
    logical, allocatable    :: patch_escape(:^D&)    !> pulsar not on cell
    type(pulsar_parametres) :: myconfig              !> pulsar parameters to be read
    character(len=78)       :: subname               !> subroutine name that call it

    type(rel_wind)          :: mywind                !> pulsar associate wind
    type(star_envelope)     :: myenvelope            !> pulsar associate envelope
    type(star)              :: mystar                !> pulsar associate progeniture star
    type(supernovae_remnant):: mysupernovae_remnant  !> pulsar associate supernovae
    contains
     PROCEDURE, PASS(self) :: set_default          => usr_pulsar_set_default
     PROCEDURE, PASS(self) :: set_complet          => usr_pulsar_set_complet
     PROCEDURE, PASS(self) :: normalize            => usr_pulsar_normalize
     PROCEDURE, PASS(self) :: set_w                => usr_pulsar_set_w
     PROCEDURE, PASS(self) :: read_parameters      => usr_pulsar_read_p
     PROCEDURE, PASS(self) :: write_setting        => usr_pulsar_write_setting
     PROCEDURE, PASS(self) :: clean_memory         => usr_pulsar_clean_memory
  end type

contains
!---------------------------------------------------------------------

!> Read the ism parameters  from a parfile
 subroutine usr_pulsar_read_p(self,pulsar_config,files)
   class(pulsar)                          :: self
   character(len=*),intent(in)            :: files(:)
   type(pulsar_parametres)                :: pulsar_config
   integer                                :: i_file

   namelist /usr_pulsar_list/ pulsar_config

   if(mype==0)write(*,*)'Reading usr_pulsar_list'
   Loop_read_file : do i_file = 1, size(files)
      open(unitpar, file=trim(files(i_file)), status="old")
      read(unitpar, usr_pulsar_list, end=113)
113    close(unitpar)
   end do Loop_read_file


   if(pulsar_config%wind_on)call self%mywind%read_parameters(&
                            self%mywind%myconfig,files)
   if(pulsar_config%envelope_on)call self%myenvelope%read_parameters(&
                              self%myenvelope%myconfig,files)
   if(pulsar_config%supernovae_remnant_on)call self%mysupernovae_remnant%read_parameters(&
                                self%mysupernovae_remnant%myconfig,files)
 end subroutine usr_pulsar_read_p
!------------------------------------------------------------------------
subroutine usr_pulsar_write_setting(self,unit_config)
  implicit none
  class(pulsar)                       :: self
  integer,intent(in)                  :: unit_config
  ! .. local ..

  !-----------------------------------

  write(unit_config,*)'************************************'
  write(unit_config,*)'**********Pulsar setting ***********'
  write(unit_config,*)'************************************'
  write(unit_config,*) 'position     = ', self%myconfig%center
  write(unit_config,*) 'radius       = ', self%myconfig%radius
  write(unit_config,*) 'Omega        = ', self%myconfig%omega
  write(unit_config,*) 'Speed       = ', self%myconfig%velocity

  if(self%myconfig%wind_on)    call self%mywind%write_setting(unit_config)
  if(self%myconfig%envelope_on)call self%myenvelope%write_setting(unit_config)
  if(self%myconfig%supernovae_remnant_on)call self%myenvelope%write_setting(unit_config)
!  if(self%myconfig%star_on)    call self%mystar%write_setting(unit_config)
  write(unit_config,*)'************************************'
  write(unit_config,*)'******* END Pulsar setting *********'
  write(unit_config,*)'************************************'
end    subroutine usr_pulsar_write_setting
!-------------------------------------------------------------------------
!> subroutine default setting for ISM
subroutine usr_pulsar_set_default(self)
 implicit none
 class(pulsar)            :: self
 !----------------------------------
 self%myconfig%unit                     = 'code'
 self%myconfig%velocity(1:3)             = 0.0_dp
 self%myconfig%omega                     = 0.0_dp
 self%myconfig%radius                    = 0.0_dp
 self%myconfig%center                    = 0.0_dp

 self%myconfig%wind_on                   = .false.
 self%myconfig%envelope_on               = .false.
 self%myconfig%star_on                   = .false.
 self%myconfig%supernovae_remnant_on     = .false.

 self%myconfig%t_sn                      = -1.0_dp
 self%myconfig%t_start_pulsar_wind       = -1.0_dp
 self%myconfig%t_end_pulsar_wind         = 1d20

 self%myconfig%t_start_star_wind          = -1.0_dp
 self%myconfig%t_end_star_wind            = -1.0_dp

 !if(self%myconfig%wind_on)
 call self%mywind%set_default
 !if(self%myconfig%envelope_on)
 call self%myenvelope%set_default
 !if(self%myconfig%supernovae_remnant_on)
 call self%mysupernovae_remnant%set_default
 !if(self%myconfig%star_on)call self%mystar%set_default
 !call self%mystar%set_default
end subroutine usr_pulsar_set_default
!--------------------------------------------------------------------
!> subroutine check the parfile setting for ism
subroutine usr_pulsar_set_complet(self)
 implicit none
 class(pulsar)            :: self
 ! .. local ..
 real(dp)                 :: mp,kb
 !-----------------------------------
 if(SI_unit) then
   mp=mp_SI
   kB=kB_SI
 else
   mp=mp_cgs
   kB=kB_cgs
 end if

 if(dabs(self%myconfig%omega)<smalldouble.and.self%myconfig%period>smalldouble) then
   self%myconfig%omega=1.0_dp/self%myconfig%period
 end if

 if(dabs(self%myconfig%omega)>smalldouble.and.self%myconfig%period<smalldouble) then
   self%myconfig%period=1.0_dp/self%myconfig%omega
 end if


 if(self%myconfig%wind_on)then
   self%mywind%myconfig%time_start = self%myconfig%t_start_pulsar_wind
   self%mywind%myconfig%time_end   = self%myconfig%t_end_pulsar_wind
   call self%mywind%set_complet(self%myconfig%velocity)
 end if

 if(self%myconfig%supernovae_remnant_on)then
    self%mysupernovae_remnant%myconfig%time_set =  self%myconfig%t_sn
    call self%mysupernovae_remnant%set_complet
 end if

 if(self%myconfig%envelope_on)then
    call self%myenvelope%set_complet
 end if


 if(self%myconfig%star_on)then
  !  call self%mystar%set_complet(self%velocity)
 end if
end subroutine usr_pulsar_set_complet
!--------------------------------------------------------------------
!> subroutine normalize setting for ISM
subroutine usr_pulsar_normalize(self,physunit_inuse)
 use mod_obj_usr_unit
 implicit none
 class(pulsar)                      :: self
 type(usr_physical_unit_values)     :: physunit_inuse
 !----------------------------------
 if(trim(self%myconfig%unit)=='code')return
 self%myconfig%center           = self%myconfig%center    / unit_length
 self%myconfig%radius           = self%myconfig%radius    / unit_length
 self%myconfig%velocity         = self%myconfig%velocity  / unit_velocity
 self%myconfig%omega            = self%myconfig%omega     / (unit_velocity/unit_length)
 if(self%myconfig%wind_on)call self%mywind%normalize()

 if(self%myconfig%envelope_on)call self%myenvelope%normalize()

 self%myconfig%t_start_star_wind   = self%myconfig%t_start_star_wind    / unit_time
 self%myconfig%t_end_star_wind     = self%myconfig%t_end_star_wind      / unit_time
 self%myconfig%t_sn                = self%myconfig%t_sn                 / unit_time
 self%myconfig%t_start_pulsar_wind = self%myconfig%t_start_pulsar_wind  / unit_time
 self%myconfig%t_end_pulsar_wind   = self%myconfig%t_end_pulsar_wind    / unit_time


 if(self%myconfig%wind_on)call self%mywind%normalize(physunit_inuse)
 if(self%myconfig%envelope_on)call self%myenvelope%normalize(physunit_inuse)
 if(self%myconfig%supernovae_remnant_on)call self%mysupernovae_remnant%normalize(physunit_inuse)
  !if(self%myconfig%star_on)call self%mystar%normalize

 if(mype==0)then
  WRITE(*,*)'********************************************************'
  WRITE(*,*)' is in pulsar'
  WRITE(*,*)'center : ',  self%myconfig%center
    WRITE(*,*)'the velocity : ', self%myconfig%velocity
  WRITE(*,*)'********************************************************'
 end if
end subroutine usr_pulsar_normalize
!--------------------------------------------------------------------
!> subroutine setting for ISM
subroutine usr_pulsar_set_w(ixI^L,ixO^L,qt,x,w,self)
 implicit none
 integer, intent(in)      :: ixI^L,ixO^L
 real(kind=dp), intent(in):: qt
 real(kind=dp)            :: x(ixI^S,1:ndir)
 real(kind=dp)            :: w(ixI^S,1:nw)
 class(pulsar)            :: self
 ! .. local..
 integer                  :: idir

 !----------------------------------
 allocate(self%patch(ixI^S))
 self%patch(ixO^S) = .false.


 cond_star_on : if(self%myconfig%star_on  &
   .and.qt>=self%myconfig%t_start_star_wind.and. &
         qt<=self%myconfig%t_start_pulsar_wind)then

!  call self%mystar%set_w(ixI^L,ixO^L,.false.,qt,x,w)

 end if cond_star_on


 cond_sn_on : if(self%myconfig%supernovae_remnant_on &
   .and.qt==self%myconfig%t_sn)then
  self%mysupernovae_remnant%subname = self%subname
  if(allocated(self%patch_escape))then
    allocate(self%mysupernovae_remnant%patch_escape(ixI^S))
    self%mysupernovae_remnant%patch_escape(ixO^S) = self%patch_escape(ixO^S)
  end if
  call self%mysupernovae_remnant%set_w(ixI^L,ixO^L,qt,x,w)
  self%patch(ixO^S) = self%patch(ixO^S).or.self%mysupernovae_remnant%patch(ixO^S)

 end if cond_sn_on


 cond_wind_on : if(self%myconfig%wind_on &
   .and.qt>=self%myconfig%t_start_pulsar_wind)then
  self%mywind%subname = self%subname
  if(allocated(self%patch_escape))then
    allocate(self%mywind%patch_escape(ixI^S))
    self%mywind%patch_escape(ixO^S) = self%patch_escape(ixO^S)
  end if
  call self%mywind%set_w(ixI^L,ixO^L,qt,x,w)
  self%patch(ixO^S) = self%patch(ixO^S).or.self%mywind%patch(ixO^S)
 end if cond_wind_on



 cond_envelope_on : if(self%myconfig%envelope_on)then
  self%myenvelope%subname = self%subname
  if(allocated(self%patch_escape))then
    allocate(self%myenvelope%patch_escape(ixI^S))
    self%myenvelope%patch_escape(ixO^S) = self%patch_escape(ixO^S)
  end if
  call self%myenvelope%set_w(ixI^L,ixO^L,qt,x,w)
  self%patch(ixO^S) = self%patch(ixO^S).or.self%myenvelope%patch(ixO^S)
 end if cond_envelope_on


end subroutine usr_pulsar_set_w


!--------------------------------------------------------------------
!> Subroutine to clean array memory of associated with cloud object
subroutine usr_pulsar_clean_memory(self)
  class(pulsar)    :: self
  !-------------------------------------

  if(allocated(self%patch))deallocate(self%patch)
  if(allocated(self%patch_escape))deallocate(self%patch_escape)
  if(self%myconfig%wind_on)call self%mywind%clean_memory
  if(self%myconfig%envelope_on)call self%myenvelope%clean_memory
  if(self%myconfig%supernovae_remnant_on)call self%mysupernovae_remnant%clean_memory
!  if(self%myconfig%star_on)call self%mystar%clean_memory
end subroutine usr_pulsar_clean_memory
!---------------------------------------------------------------------


end module mod_obj_pulsar
