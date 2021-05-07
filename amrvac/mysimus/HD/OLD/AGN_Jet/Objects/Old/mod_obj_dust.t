module mod_obj_dust
  use mod_constants
  use mod_global_parameters
  use mod_physics
  use mod_hd, only: hd_dust
  use mod_srmhd_parameters !, only: mag,lfac_,psi_,xi_
  use mod_dust
  use mod_obj_global_parameters
  use mod_obj_mat
  use mod_obj_usr_unit
  implicit none

  integer, parameter  :: n_dust_max = 20

  type dust_parameters
     character(len=20)    :: unit               !> physical unit at parameter file
     integer              :: myindice           !> dust associated indice
     character(len=78)    :: obj_name           !> dust obj name
     logical              :: normalize_done     !> dust is the normalisation is already done
     character(len=20)    :: distrub_func       !> dust distribution function
     integer              :: n_species          !> number of dust species
     real(dp)             :: sizes(1:n_dust_max)!> dust size
     real(dp)             :: min_radius         !> dust smaller radius
     real(dp)             :: max_radius         !> dust bigger radius
     real(dp)             :: power_a            !> dust power indice  distrbution func
     real(dp)             :: grain_density(n_dust_max) !> dust grain density (g/cm)
     real(dp)             :: density(n_dust_max)!> dust density
     real(dp)             :: number_density(n_dust_max)!> dust density
     real(dp)             :: temperature        !> dust temperature  (K)
     real(dp)             :: min_limit_rel      !> dust floor relative to gas
     real(dp)             :: min_limit_abs      !> dust floor abs (g/cm^3)
     real(dp)             :: max_limit_rel      !> dust high relative to gas
     real(dp)             :: extend(1:2,1:3)    !> region in space (cm)
     real(dp)             :: velocity(1:3)      !> dust velocity (cm/s)
     logical              :: tracer_on          !> logical to set tracer
     character(len=78)    :: associated_medium  !> dust associated medium
     logical              :: associate_med      !> logical to associate dust to medium
     integer              :: idust_first        !> indice of first dust species
     integer              :: idust_last         !> indice of first dust species
  end type dust_parameters

  type dust_ispecies
    logical, allocatable :: patch(:^D&)
  end type dust_ispecies
  ! the dust type
  type dust
   logical, allocatable              :: patch(:^D&)        !> dust is on cell
   character(len=78)                 :: subname            !> dust subroutine name that call it
   type(dust_ispecies), allocatable  :: the_ispecies(:)    !> dust configuration for each dust especies
   type(dust_parameters)             :: myconfig           !> dust parameters
   type(usrphysical_unit), pointer   :: myphysunit         !> dust physicq unit in use
   contains
   !PRIVATE
   !PUBLIC
   PROCEDURE, PASS(self) :: set_default     => usr_dust_set_default
   PROCEDURE, PASS(self) :: set_complet     => usr_dust_set_complet
   PROCEDURE, PASS(self) :: to_phys         => usr_dust_to_phys
   PROCEDURE, PASS(self) :: normalize       => usr_dust_normalize
   PROCEDURE, PASS(self) :: set_w           => usr_dust_set_w
   PROCEDURE, PASS(self) :: set_w_zero      => usr_dust_set_w_zero
   PROCEDURE, PASS(self) :: handel_small_val=> usr_dust_handel_small_val
   PROCEDURE, PASS(self) :: read_parameters => usr_dust_read_p
   PROCEDURE, PASS(self) :: write_setting   => usr_dust_write_setting
   PROCEDURE, PASS(self) :: set_patch       => usr_dust_set_patch
   PROCEDURE, PASS(self) :: set_allpatch    => usr_dust_set_allpatch
   PROCEDURE, PASS(self) :: clean_memory    => usr_dust_clean_memory
  end type

contains

  !---------------------------------------------------------------------

   !> Read the ism parameters  from a parfile
    subroutine usr_dust_read_p(self,dust_config,files)
      class(dust)                     :: self
      type(dust_parameters)           :: dust_config
      character(len=*), intent(in)    :: files(:)
      ! .. local ..
      integer                         :: i_file,i_error_read
      !--------------------------------------------------
      namelist /usr_dust_ism_list/ dust_config
      namelist /usr_dust_cloud_list/ dust_config
      namelist /usr_dust_sn_list/ dust_config
      namelist /usr_dust_list/ dust_config

      if(mype==0)write(*,*)'Reading usr_list'
      select case(self%myconfig%associated_medium)
      case('ism')
       Loop_read_files_dustism : do i_file = 1, size(files)
         open(unitpar, file=trim(files(i_file)), status="old")
         read(unitpar, usr_dust_ism_list, iostat=i_error_read)
         call usr_mat_read_error_message(i_error_read,self%myconfig%myindice,&
                                      self%myconfig%obj_name)
         close(unitpar)
       end do Loop_read_files_dustism
     case('cloud')
       Loop_read_files_dustcloud : do i_file = 1, size(files)
         open(unitpar, file=trim(files(i_file)), status="old")
         read(unitpar, usr_dust_cloud_list, iostat=i_error_read)
         call usr_mat_read_error_message(i_error_read,self%myconfig%myindice,&
                                      self%myconfig%obj_name)
         close(unitpar)
       end do Loop_read_files_dustcloud
     case('sn')
        Loop_read_files_dustsn : do i_file = 1, size(files)
          open(unitpar, file=trim(files(i_file)), status="old")
          read(unitpar, usr_dust_sn_list, iostat=i_error_read)
          close(unitpar)
        end do Loop_read_files_dustsn
     case default
       Loop_read_files_dust : do i_file = 1, size(files)
         open(unitpar, file=trim(files(i_file)), status="old")
         read(unitpar, usr_dust_list, iostat=i_error_read)
         call usr_mat_read_error_message(i_error_read,self%myconfig%myindice,&
                                      self%myconfig%obj_name)
         close(unitpar)
       end do Loop_read_files_dust
     end select
    end subroutine usr_dust_read_p
  !------------------------------------------------------------------------
  !> write the cloud setting
  subroutine usr_dust_write_setting(self,unit_config)
    implicit none
    class(dust)                         :: self
    integer,intent(in)                  :: unit_config
    ! .. local ..

    !-----------------------------------

    write(unit_config,*)'**************************'
    write(unit_config,*)'**** Dust  setting *******'
    write(unit_config,*)'**************************'
    write(unit_config,*) 'Density       = ', self%myconfig%density(1:self%myconfig%idust_last-self%myconfig%idust_first+1)
    write(unit_config,*) 'Grain density = ', self%myconfig%grain_density(1:self%myconfig%idust_last-self%myconfig%idust_first+1)
    write(unit_config,*) 'Temperature   = ', self%myconfig%temperature
    write(unit_config,*) 'Speed         = ', self%myconfig%velocity
    write(unit_config,*) 'n_species     = ', self%myconfig%n_species
    write(unit_config,*)'**************************'
    write(unit_config,*)'**** End dust setting ***'
    write(unit_config,*)'**************************'
  end    subroutine usr_dust_write_setting
  !--------------------------------------------------------------------
  !> subroutine default setting for ISM
   subroutine usr_dust_set_default(self)
    implicit none
    class(dust)          :: self
    !----------------------------------

    self%myconfig%unit               = 'code'
    self%myconfig%obj_name           = 'dust'
    self%myconfig%myindice           = 1
    self%myconfig%distrub_func       = 'uniform'
    self%myconfig%associated_medium  = 'default'
    self%myconfig%n_species          = 0
    self%myconfig%sizes              = 0.0_dp
    self%myconfig%density            = 0.0_dp
    self%myconfig%number_density     = 0.0_dp
    self%myconfig%temperature        = 0.0_dp
    self%myconfig%velocity           = 0.0_dp
    self%myconfig%min_limit_rel      = 0.0_dp
    self%myconfig%min_limit_abs      = 0.0_dp
    self%myconfig%max_limit_rel      = 1.0d9
    self%myconfig%extend(1:2,1:3)    = 0.0_dp!box_limit(1:2,1:ndim) ! it zero at this time
    self%myconfig%associate_med      = .true.
    self%myconfig%grain_Density      = 0.0_dp
    self%myconfig%min_radius         =-1.0_dp
    self%myconfig%max_radius         =-1.0_dp
    self%myconfig%idust_first        = 1
    self%myconfig%idust_last         = 1
    self%myconfig%tracer_on          = .false.
    self%myconfig%normalize_done     = .false.


   end subroutine usr_dust_set_default
  !--------------------------------------------------------------------
  !> subroutine setting for dust
  subroutine usr_dust_set_w_zero(ixI^L,ixO^L,x,w,self)
    use mod_dust
   implicit none
   integer, intent(in)          :: ixI^L,ixO^L
   real(kind=dp), intent(in)    :: x(ixI^S,1:ndim)
   real(kind=dp), intent(inout) :: w(ixI^S,1:nw)
   class(dust)                  :: self
   ! .. local ..
   integer                      :: idust,idir
   !----------------------------------------------------

    Loop_idust01 : do idust = self%myconfig%idust_first,self%myconfig%idust_last
      if(all(self%the_ispecies(idust)%patch(ixO^S)))cycle Loop_idust01
      where(.not.self%patch(ixO^S))
        w(ixO^S,phys_ind%dust_rho(idust)) = 0.0_dp!1.1d-2*self%myconfig%min_limit_abs
      end where
      Loop_idir : do idir=1,ndir
          where(.not.self%patch(ixO^S))
           w(ixO^S,phys_ind%dust_mom(idir,idust)) = 0.0_dp!self%myconfig%velocity(idir)
          end where
      end do Loop_idir
    end do Loop_idust01

  end  subroutine usr_dust_set_w_zero
  !--------------------------------------------------------------------

   !> subroutine setting W for dust
   subroutine usr_dust_set_w(ixI^L,ixO^L,qt,is_frac,dust_frac,fprofile,x,w,self)
    implicit none
    integer, intent(in)          :: ixI^L,ixO^L
    real(kind=dp), intent(in)    :: qt
    logical                      :: is_frac
    real(kind=dp), intent(in)    :: dust_frac
    real(kind=dp), intent(in)    :: fprofile(ixI^S)
    real(kind=dp), intent(in)    :: x(ixI^S,1:ndim)
    real(kind=dp), intent(inout) :: w(ixI^S,1:nw)
    class(dust)                  :: self
    ! .. local..
    real(kind=dp)                :: slop,dr_dust,dust_frac_loc
    integer                      :: idir,idust,istart_dust,istop_dust,i_dust_loc
    !----------------------------------

    if(.not.any(self%patch(ixO^S)))return
    istart_dust=1
    istop_dust = (istart_dust-1)+self%myconfig%n_species
    if(is_frac) then
      dust_frac_loc=dust_frac/(1.0_dp-dust_frac)
    else
      dust_frac_loc=dust_frac
    end if


     select case(self%myconfig%distrub_func)
     case('uniform')
      cond_fracdust : if(dust_frac_loc<0.0_dp) then

       Loop_idust2 : do idust = self%myconfig%idust_first, self%myconfig%idust_last
        i_dust_loc = idust-self%myconfig%idust_first+1
        where(self%patch(ixO^S))
         w(ixO^S,phys_ind%dust_rho(idust)) = self%myconfig%density(i_dust_loc)&
                                    / real(self%myconfig%n_species,kind=dp)
        end where
      end do Loop_idust2
      else cond_fracdust

       Loop_idust3 : do idust = self%myconfig%idust_first, self%myconfig%idust_last
        where(self%patch(ixO^S))
         w(ixO^S,phys_ind%dust_rho(idust)) = dust_frac_loc*w(ixO^S,phys_ind%rho_)&
                                    / real(self%myconfig%n_species,kind=dp)
        end where
      end do Loop_idust3


      end if cond_fracdust

     case('powerlaw')



      Loop_idust4 : do idust = self%myconfig%idust_first, self%myconfig%idust_last
       dr_dust=0.5_dp*(dust_size(max(idust-1,self%myconfig%idust_first))-&
                       dust_size(min(idust+1,self%myconfig%idust_last)))
       where(self%patch(ixO^S))
        w(ixO^S,phys_ind%dust_rho(idust)) = w(ixO^S,phys_ind%rho_)*dust_frac_loc*&
           0.5_dp*(dr_dust/dsqrt(dust_size(idust))&
           / (dsqrt(dust_size(self%myconfig%idust_first))-dsqrt(dust_size(self%myconfig%idust_last))))
        w(ixO^S,phys_ind%dust_rho(idust)) = w(ixO^S,phys_ind%dust_rho(idust)) *fprofile(ixO^S)
       end where
      end do Loop_idust4

     case default
      write(*,*)'The ', self%myconfig%distrub_func, ' is not implimented'
      call mpistop('is not implimented in usr_dust_set_w at mod_usr.t')
     end select



    cond_isassociated : if(self%myconfig%associate_med)then
     Loop_idust5 : do idust = self%myconfig%idust_first, self%myconfig%idust_last
        Loop_idir0 : do idir=1,ndir
          where(self%patch(ixO^S))
           w(ixO^S,phys_ind%dust_mom(idir,idust)) = w(ixO^S,phys_ind%mom(idir))
          end where
        end do Loop_idir0
      end do Loop_idust5
    else cond_isassociated
     Loop_idust6 : do idust = self%myconfig%idust_first, self%myconfig%idust_last
      Loop_idir : do idir=1,ndir
       where(self%patch(ixO^S))
       w(ixO^S,phys_ind%dust_mom(idir,idust)) = self%myconfig%velocity(idir)
       end where
      end do Loop_idir
     end do Loop_idust6
    end if cond_isassociated

   end subroutine usr_dust_set_w
   !--------------------------------------------------------------------
   !> subroutine check the parfile setting for cloud
   subroutine usr_dust_handel_small_val(ixI^L,ixO^L,qt,x,w,self)
     implicit none
     integer, intent(in)          :: ixI^L,ixO^L
     real(kind=dp), intent(in)    :: qt
     real(kind=dp), intent(in)    :: x(ixI^S,1:ndim)
     real(kind=dp), intent(inout) :: w(ixI^S,1:nw)
     class(dust)                  :: self
     ! .. local ..
     integer                      :: idust,idir
     real(kind=dp)                :: small_dust_rho
     logical, dimension(ixI^S)    :: patch_correct,patch_slow
     !-----------------------------------------------------

     small_dust_rho = self%myconfig%min_limit_rel


     ! handel small density dust
     Loop_idust : do idust =self%myconfig%idust_first, self%myconfig%idust_last
      where(w(ixI^S, phys_ind%dust_rho(idust))<max(small_dust_rho*w(ixI^S,phys_ind%rho_),&
         self%myconfig%min_limit_abs))
        w(ixI^S, phys_ind%dust_rho(idust))= 0.8* min(small_dust_rho*w(ixI^S,phys_ind%rho_),&
             self%myconfig%min_limit_abs)
        patch_correct(ixI^S) = .true.
      elsewhere
        patch_correct(ixI^S) = .false.
      end where
      ! handel large density dust
      where(w(ixI^S, phys_ind%dust_rho(idust))>self%myconfig%max_limit_rel*w(ixI^S,phys_ind%rho_))
        w(ixI^S, phys_ind%dust_rho(idust))=0.8*self%myconfig%max_limit_rel*w(ixI^S,phys_ind%rho_)
        patch_slow(ixI^S) = .true.
      elsewhere
        patch_slow(ixI^S) =.false.
      end where



      Loop_idir1 : do idir = 1,ndir
       where(patch_correct(ixI^S))
        w(ixI^S, dust_mom(idir,idust))=0.0_dp
       elsewhere(patch_slow(ixI^S))
        w(ixI^S, phys_ind%dust_mom(idir,idust))=w(ixI^S,phys_ind%mom(idir))
       end where
      end do Loop_idir1
     end do  Loop_idust

   end  subroutine usr_dust_handel_small_val
   !--------------------------------------------------------------------
   !> subroutine check the parfile setting for dust
   subroutine usr_dust_set_complet(self,is_frac,dust_frac,med_density,med_velocity)
     implicit none
     class(dust)              :: self
     logical                  :: is_frac
     real(dp), intent(in)     :: dust_frac,med_density,med_velocity(1:ndir)
     ! .. local ..
     integer                  :: idust,i_dust_loc
     real(dp)                 :: mp,kb,slop,dr_dust,dust_frac_loc
     !-----------------------------------
     if(SI_unit) then
       mp=mp_SI
       kB=kB_SI
     else
       mp=mp_cgs
       kB=kB_cgs
     end if
     self%myconfig%idust_last=(self%myconfig%idust_first-1) +self%myconfig%n_species
     if(is_frac) then
       dust_frac_loc=dust_frac/(1.0_dp-dust_frac)
     else
       dust_frac_loc=dust_frac
     end if
     cond_size_off : if(any(self%myconfig%sizes(1:self%myconfig%n_species)<=smalldouble))then

      select case(self%myconfig%distrub_func)
      case('uniform')

       slop=(self%myconfig%max_radius-self%myconfig%min_radius)/dble(self%myconfig%n_species-1)
       Loop_idust1 : do idust = self%myconfig%idust_first, self%myconfig%idust_last
         i_dust_loc = idust-self%myconfig%idust_first+1
         self%myconfig%sizes(i_dust_loc)   = self%myconfig%min_radius+slop*real(idust-1,kind=dp)
         self%myconfig%density(i_dust_loc) = &
              dust_frac_loc* med_density/real(self%myconfig%n_species,kind=dp)
       end do Loop_idust1

      case('powerlaw')
       call usr_mat_powerlaw_withX(self%myconfig%n_species,self%myconfig%power_a,self%myconfig%min_radius,&
                                   self%myconfig%max_radius,   &
                                   self%myconfig%sizes)


       Loop_idust2  : do idust = self%myconfig%idust_first, self%myconfig%idust_last
        i_dust_loc = idust-self%myconfig%idust_first+1
        dr_dust=0.5_dp*(self%myconfig%sizes(max(i_dust_loc-1,1))-&
                        self%myconfig%sizes(min(i_dust_loc+1,&
                      self%myconfig%idust_last-self%myconfig%idust_first+1)))

        self%myconfig%density(i_dust_loc) =med_density*dust_frac_loc*   &
               half*(dr_dust/dsqrt(self%myconfig%sizes(i_dust_loc))     &
          / (dsqrt(self%myconfig%sizes(1))                              &
          -dsqrt(self%myconfig%sizes(self%myconfig%idust_last           &
                                    -self%myconfig%idust_first+1))))

       end do  Loop_idust2
      case default
       write(*,*)'The ', self%myconfig%distrub_func, ' is not implimented'
       call mpistop('is not implimented in usr_dust_set_w at mod_usr.t')
      end select
     end if cond_size_off
     if(self%myconfig%associate_med)then
       self%myconfig%velocity = med_velocity
     end if
     allocate(self%the_ispecies(self%myconfig%idust_first:self%myconfig%idust_last))

   end subroutine usr_dust_set_complet
   !--------------------------------------------------------------------
   !> subroutine normalize setting for ISM
    subroutine usr_dust_to_phys(self)
      implicit none
      class(dust)                        :: self
      ! .. local ..
      integer                            :: idust,i_dust_loc
      !----------------------------------

      Loop_idust1 : do idust = self%myconfig%idust_first, self%myconfig%idust_last
        i_dust_loc            = idust-self%myconfig%idust_first+1
        dust_size(idust)      = self%myconfig%sizes(i_dust_loc)
        dust_density(idust)   = self%myconfig%grain_density(i_dust_loc)
      end do Loop_idust1



    end subroutine usr_dust_to_phys
  !--------------------------------------------------------------------
  !> subroutine normalize setting for ISM
   subroutine usr_dust_normalize(self,physunit_inuse)
    use mod_obj_usr_unit
    implicit none
    class(dust)                                   :: self
    type(usrphysical_unit),target, intent(in)     :: physunit_inuse
    !----------------------------------
    self%myphysunit =>physunit_inuse

    if(trim(self%myconfig%unit)=='code'.or.self%myconfig%normalize_done)then
     if(self%myconfig%normalize_done)then
      write(*,*) 'WARNING: Second call for dust normalisation:: ', &
                   '  no new normalisation will be done'
      write(*,*)'the reason : ' , 'code unit: ', self%myconfig%unit,&
                'normalisation is done ',   self%myconfig%normalize_done
      write(*,*)' called from :: ',self%myconfig%associated_medium
     end if
     return
    end if
    self%myconfig%sizes            = self%myconfig%sizes         / physunit_inuse%myconfig%length
    self%myconfig%density          = self%myconfig%density       / physunit_inuse%myconfig%density
    self%myconfig%number_density   = self%myconfig%number_density/ physunit_inuse%myconfig%number_density
    self%myconfig%grain_density    = self%myconfig%grain_density / physunit_inuse%myconfig%density
    self%myconfig%temperature      = self%myconfig%temperature   / physunit_inuse%myconfig%temperature
    self%myconfig%velocity         = self%myconfig%velocity      / physunit_inuse%myconfig%velocity
    self%myconfig%min_limit_abs    = self%myconfig%min_limit_abs / physunit_inuse%myconfig%density
    self%myconfig%extend           = self%myconfig%extend        / physunit_inuse%myconfig%length

    self%myconfig%normalize_done =.true.
   end subroutine usr_dust_normalize
  !--------------------------------------------------------------------
  !> Subroutine set the patch array memory of associated with cloud object

   subroutine usr_dust_set_patch(ixI^L,ixO^L,imposed_patch,self)
    implicit none
    integer, intent(in)          :: ixI^L,ixO^L
    logical, intent(in)          :: imposed_patch(ixI^S)
    class(dust)                  :: self
    !.. local ..
    integer                      :: idust
    !-------------------------------------------------------------

    if(allocated(self%patch))deallocate(self%patch)
    allocate(self%patch(ixI^S))
    self%patch(ixO^S) = imposed_patch(ixO^S)
    Loop_idust :  do idust = self%myconfig%idust_first, self%myconfig%idust_last
     if(allocated(self%the_ispecies(idust)%patch))then
       deallocate(self%the_ispecies(idust)%patch)
     end if
     allocate(self%the_ispecies(idust)%patch(ixO^S))
     self%the_ispecies(idust)%patch(ixO^S)=self%patch(ixO^S)
    end do Loop_idust

   end subroutine usr_dust_set_patch
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !> Subroutine set the patch array memory of associated with cloud object
    subroutine usr_dust_set_allpatch(ixI^L,ixO^L,all_dust_obj,self)
     implicit none
     integer, intent(in)          :: ixI^L,ixO^L
     type(dust), intent(in)       :: all_dust_obj(:)
     class(dust)                  :: self
     !.. local ..
     integer                      :: idust,i_obj,n_obj,i_dust_loc,i_obj_start
     !-------------------------------------------------------------

     n_obj=size(all_dust_obj)


     if(allocated(self%patch))deallocate(self%patch)
     allocate(self%patch(ixI^S))
     self%patch(ixO^S) = .false.

     Loop_allobj_1: do i_obj=1,n_obj
       if(allocated(all_dust_obj(i_obj)%patch))then
        self%patch(ixO^S)=self%patch(ixO^S).or.all_dust_obj(i_obj)%patch(ixO^S)
       end if
     end do Loop_allobj_1


     self%myconfig%idust_first=minval(all_dust_obj(1:n_obj)%myconfig%idust_first)
     self%myconfig%idust_last=maxval(all_dust_obj(1:n_obj)%myconfig%idust_last)
     if(allocated(self%the_ispecies)) then
       deallocate(self%the_ispecies)
     end if

     allocate(self%the_ispecies(self%myconfig%idust_first:self%myconfig%idust_last))


     Loop_idust :  do idust = self%myconfig%idust_first, self%myconfig%idust_last
      if(allocated(self%the_ispecies(idust)%patch))then
        deallocate(self%the_ispecies(idust)%patch)
      end if
      allocate(self%the_ispecies(idust)%patch(ixO^S))
      self%the_ispecies(idust)%patch(ixO^S)=.false.
      if(n_obj>1)then
       Loop_allobj_2: do i_obj=1,n_obj
        if(idust>=all_dust_obj(i_obj)%myconfig%idust_first.and.&
           idust<=all_dust_obj(i_obj)%myconfig%idust_last) then
           !i_dust_loc = idust-all_dust_obj(i_obj)%myconfig%idust_first+1
           cond_obj2_on : if(allocated(all_dust_obj(i_obj)%the_ispecies(idust)%patch))then
            self%the_ispecies(idust)%patch(ixO^S)=self%the_ispecies(idust)%patch(ixO^S).or.&
                all_dust_obj(i_obj)%the_ispecies(idust)%patch(ixO^S)
           end if  cond_obj2_on

        end if
       end do  Loop_allobj_2
      end if
     end do Loop_idust

   end subroutine usr_dust_set_allpatch
    !--------------------------------------------------------------------
   !> Subroutine to clean array memory of associated with cloud object
   subroutine usr_dust_clean_memory(self)
     class(dust)    :: self
     !.. local ..
     integer        :: idust
     !-------------------------------------
     if(allocated(self%patch))deallocate(self%patch)

     if(allocated(self%the_ispecies))then
      Loop_idust : do idust = self%myconfig%idust_first, self%myconfig%idust_last
       if(allocated(self%the_ispecies(idust)%patch))then
         deallocate(self%the_ispecies(idust)%patch)
       end if
      end do Loop_idust
     end if
   end subroutine usr_dust_clean_memory
  !-------------------------------------------------------------------------
end module mod_obj_dust
