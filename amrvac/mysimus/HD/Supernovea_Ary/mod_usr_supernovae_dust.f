module mod_usr
  use mod_hd!, only: mhd_n_tracer,mhd_dust
  use mod_physics
  use mod_global_parameters
  use mod_obj_global_parameters
  use mod_obj_mat
  use mod_obj_cloud
  use mod_obj_ism
  use mod_obj_sn_remnant
  use mod_obj_usr_unit
  implicit none
  save
  real(dp) :: theta, kx, ly, vc

  type usr_config
    logical           :: physunit_on
    logical           :: sn_on
    logical           :: ism_on
    logical           :: cloud_on
    logical           :: ism_list_diff
    logical           :: cloud_list_diff
    integer           :: cloud_number,ism_number
    character(len=30) :: coordinate_system
  end type usr_config
  type(usr_config) :: usrconfig
  integer, parameter  :: n_dust_max = 20
  real(dp) :: SUM_MASS   = 0.0_dp
  real(dp) :: SUM_VOLUME = 0.0_dp


  type (ISM),allocatable,target      :: ism_surround(:)
  type (cloud),allocatable,target    :: cloud_medium(:)
  type (ISM),target                  :: ism_default
  type (cloud),target                :: cloud_default
  type (dust),target                 :: dust_ary
  type (dust),allocatable,target     :: the_dust_inuse(:)
  type (supernovae_remnant), target  :: sn_wdust

  !type(star) :: star_ms
  !type(star) :: sun

  type(usrphysical_unit) :: usr_physunit





contains
  subroutine usr_init
    ! .. local ..
    integer :: i_cloud,i_ism
    !-------------------------------------------
    ! configuration of procedures to be used in this project
    usr_set_parameters  => initglobaldata_usr
    usr_init_one_grid   => initonegrid_usr
    usr_special_bc      => specialbound_usr
    usr_aux_output      => specialvar_output
    usr_add_aux_names   => specialvarnames_output
    usr_source          => specialsource_usr
    usr_refine_grid     => specialrefine_usr
    usr_special_global  => usr_global_var
    usr_process_grid    => process_grid_usr
    usr_get_dt          => special_get_dt

    call usr_set_default_parameters



    call usr_physunit%set_default

    ! set default values for supernovae remnant configuration
    call sn_wdust%set_default

    ! set default values for ISMs configuration
    call ism_default%set_default


    ! set default values for clouds configuration
    call cloud_default%set_default



    call usr_params_read(par_files)



    ! complet all physical unit in use
    if(usrconfig%physunit_on) then
     call usr_physunit%set_complet
    end if
    call usr_physical_unit
    call set_coordinate_system(trim(usrconfig%coordinate_system))
    call hd_activate


    call usr_check_conflict


  end subroutine usr_init
  !------------------------------------------------------------------
  !> default usr parameters from a file
  subroutine usr_set_default_parameters
    !-------------------------------------
    usrconfig%physunit_on         = .false.
    usrconfig%sn_on               = .false.
    usrconfig%ism_on              = .false.
    usrconfig%cloud_on            = .false.
    usrconfig%cloud_number        = 1
    usrconfig%ism_number          = 1
    usrconfig%ism_list_diff       = .false.
    usrconfig%cloud_list_diff     = .false.
    usrconfig%coordinate_system   = 'slab'
  end subroutine usr_set_default_parameters
  !------------------------------------------------------------------
  !> Read this module s parameters from a file
  subroutine usr_params_read(files)
    character(len=*), intent(in) :: files(:)
    integer                      :: i_file,i_cloud,i_ism
    !-------------------------------------
    namelist /usr_list/ usrconfig


    if(mype==0)write(*,*)'Reading usr_list'
    Loop_ifile : do i_file = 1, size(files)
       open(unitpar, file=trim(files(i_file)), status="old")
       read(unitpar, usr_list)
       close(unitpar)
    end do Loop_ifile





    if(usrconfig%physunit_on)then
      call usr_physunit%read_parameters(usr_physunit%myconfig,files)
    else
      call usr_unit_read(files)
      call usr_physunit%set_to_one
    end if

    if(usrconfig%sn_on)call sn_wdust%read_parameters(sn_wdust%myconfig,files)

    if(usrconfig%ism_on)then
      allocate(ism_surround(0:usrconfig%ism_number-1))
      Loop_allism : do i_ism =0,usrconfig%ism_number-1
       ism_surround(i_ism)%myconfig = ism_default%myconfig
       ism_surround(i_ism)%myconfig%myindice=i_ism
       call ism_surround(i_ism)%read_parameters(ism_surround(i_ism)%myconfig,&
          files)
      end do Loop_allism
    end if

    if(usrconfig%cloud_on)then
      allocate(cloud_medium(0:usrconfig%cloud_number-1))
      Loop_allcloud : do i_cloud =0,usrconfig%ism_number
       cloud_medium(i_cloud)%myconfig          = cloud_default%myconfig
       cloud_medium(i_cloud)%myconfig%myindice = i_cloud

       call cloud_medium(i_cloud)%read_parameters(files,&
          cloud_medium(i_cloud)%myconfig)
      end do Loop_allcloud
    end if

  end subroutine usr_params_read

  !> subroutine to clean memory at the end
  subroutine usr_clean_memory_final
    if(usrconfig%ism_on)then
      if(allocated(ism_surround))deallocate(ism_surround)
    end if
    if(usrconfig%cloud_on)then
      if(allocated(cloud_medium))deallocate(cloud_medium)
    end if
    if(allocated(the_dust_inuse))deallocate(the_dust_inuse)
  end subroutine usr_clean_memory_final

!-------------------------------------------------------------------
!> subroutine read unit used in the code
  subroutine usr_unit_read(files)
   implicit none
   character(len=*), intent(in) :: files(:)
   integer                      :: i_file

   namelist /usr_unit_list/ unit_length , unit_time,unit_velocity,&
                unit_density, unit_numberdensity,&
                           unit_pressure,unit_temperature


  if(mype==0)write(*,*)'Reading usr_unit_list'
  Loop_read_usrfile : do i_file = 1, size(files)
         open(unitpar, file=trim(files(i_file)), status="old")
         read(unitpar, usr_unit_list, end=109)
  109    close(unitpar)
  end do Loop_read_usrfile
 end subroutine usr_unit_read
  !-----------------------------------------------------------
  !> subroutine to check configuration conflits
  subroutine usr_check_conflict
    implicit none
    ! .. local ..
    integer  :: i_ism,i_cloud
    !------------------------------
    cond_dust_on : if(.not.phys_config%dust_on)then
      if(usrconfig%ism_on)then
       Loop_isms : do i_ism=0,usrconfig%ism_number-1
        ism_surround(i_ism)%myconfig%dust_on =.false.
       end do Loop_isms
      end if
      if(usrconfig%cloud_on)then
       Loop_clouds : do i_cloud=0,usrconfig%cloud_number-1
        cloud_medium(i_cloud)%myconfig%dust_on =.false.
       end do Loop_clouds
      end if
    end if  cond_dust_on
  end   subroutine usr_check_conflict
  !-----------------------------------------------------------
  !> subroutine to normalize parameters in the code
  subroutine usr_normalise_parameters
   implicit none
   integer            :: idust



   constusr%G         = constusr%G*(unit_density*(unit_length/unit_velocity)**&
      (2.0_dp))


   constusr%clight                      = constusr%clight/unit_velocity

   ! complet all physical unit in use
   if(usrconfig%physunit_on) then
      call usr_physunit%fillphysunit
    end if


    w_convert_factor(phys_ind%rho_)              = unit_density
    if(srmhd_energy)w_convert_factor(phys_ind%e_)= &
       unit_density*unit_velocity**2.0
    if(saveprim)then
     w_convert_factor(phys_ind%mom(:))           = unit_velocity
    else
     w_convert_factor(phys_ind%mom(:))           = unit_density*unit_velocity
    end if
    time_convert_factor                 = unit_time
    length_convert_factor               = unit_length




    if(phys_config%dust_on)then
          w_convert_factor(phys_ind%dust_rho(:))              = unit_density
          if(saveprim)then
            do idust = 1, phys_config%dust_n_species
             w_convert_factor(phys_ind%dust_mom(:,idust))     = unit_velocity
            end do
          else
            do idust = 1, phys_config%dust_n_species
             w_convert_factor(phys_ind%dust_mom(:,&
                idust))     = unit_density*unit_velocity
           end do
          end if
    end if

  end subroutine usr_normalise_parameters


!-------------------------------------------------------------------------
  subroutine initglobaldata_usr
   use mod_variables
   implicit none
   ! .. local ..
   integer   :: i_cloud,i_ism,n_objects
   !------------------------------------
    n_objects =0
    itr=1
   ! complet ism parameters
   if(usrconfig%ism_on)then
    Loop_isms : do i_ism=0,usrconfig%ism_number-1
     ism_surround(i_ism)%myconfig%itr=itr
     call ism_surround(i_ism)%set_complet
     call ism_surround(i_ism)%normalize(usr_physunit)
    end do Loop_isms
    itr=ism_surround(usrconfig%ism_number-1)%myconfig%itr+1
    n_objects = n_objects + usrconfig%ism_number
   end if



   ! complet cloud parameters
   if(usrconfig%cloud_on) then
     Loop_clouds : do i_cloud=0,usrconfig%cloud_number-1
      cloud_medium(i_cloud)%myconfig%itr=itr
      call cloud_medium(i_cloud)%set_complet
      call cloud_medium(i_cloud)%normalize(usr_physunit)
     end do Loop_clouds
     itr=cloud_medium(usrconfig%cloud_number-1)%myconfig%itr+1
     n_objects = n_objects + usrconfig%cloud_number
   end if

   if(usrconfig%sn_on)then
     sn_wdust%myconfig%itr=itr
     call sn_wdust%set_complet
     call sn_wdust%normalize(usr_physunit)
     itr=sn_wdust%myconfig%itr+1
     n_objects = n_objects + 1
   end if

   allocate(the_dust_inuse(n_objects))
   call usr_normalise_parameters
   if(mype==0)call usr_write_setting

  end subroutine initglobaldata_usr
  !> The initial conditions
  subroutine initonegrid_usr(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x)
    ! initialize one grid

    implicit none

    integer, intent(in)     :: ixImin1,ixImax1,ixOmin1,ixOmax1
    real(dp), intent(in)    :: x(ixImin1:ixImax1,1:ndim)
    real(dp), intent(inout) :: w(ixImin1:ixImax1,1:nw)
    !.. local ..
    real(dp)      :: res
    integer       :: ix1,na,flag(ixImin1:ixImax1)
    integer       :: i_cloud,i_ism
    logical, save :: first=.true.
    logical       :: patch_all(ixImin1:ixImax1)
    type(dust)    :: dust_dummy
    integer       :: i_object
    ! .. only test ..
    real(dp)      ::old_w(ixOmin1:ixOmax1,1:nw)
    !-----------------------------------------
    patch_all(ixOmin1:ixOmax1) = .true.
    if(first)then
      if(mype==0) then
        write(*,*)'supernovae start :-)'
      endif
      first=.false.
    endif

    i_object=1
    ! set the ism
    if(usrconfig%ism_on) then
      Loop_isms : do i_ism=0,usrconfig%ism_number-1
       call ism_surround(i_ism)%set_w(ixImin1,ixImax1,ixOmin1,ixOmax1,&
          global_time,x,w)
       patch_all(ixOmin1:ixOmax1) =  patch_all(ixOmin1:ixOmax1) .and. &
          .not.ism_surround(i_ism)%patch(ixOmin1:ixOmax1)
       the_dust_inuse(i_object)=ism_surround(i_ism)%mydust
       i_object = i_object +1
      end do Loop_isms
    end if
    ! set one cloud
    if(usrconfig%cloud_on)then
      Loop_clouds : do i_cloud=0,usrconfig%cloud_number-1
       call cloud_medium(i_cloud)%set_w(ixImin1,ixImax1,ixOmin1,ixOmax1,&
          global_time,x,w)

       if(usrconfig%ism_on)then
         if(ism_surround(0)%myconfig%tracer_on)then
           where(cloud_medium(i_cloud)%patch(ixOmin1:ixOmax1))
             w(ixOmin1:ixOmax1,phys_ind%tracer(ism_surround(&
                0)%myconfig%itr))=0.0_dp
           end where
         end if
       end if

       patch_all(ixOmin1:ixOmax1) =  patch_all(ixOmin1:ixOmax1) .and. &
          .not.cloud_medium(i_cloud)%patch(ixOmin1:ixOmax1)
       the_dust_inuse(i_object)=cloud_medium(i_cloud)%mydust
       i_object = i_object +1
      end do Loop_clouds
    end if

    ! set the pulsar and associated wind + envelope if they are on
    if(usrconfig%sn_on)then
      sn_wdust%subname='initonegrid_usr'
      call sn_wdust%set_w(ixImin1,ixImax1,ixOmin1,ixOmax1,global_time,x,w)
      if(usrconfig%ism_on)then
        if(ism_surround(0)%myconfig%tracer_on)then
          where(sn_wdust%patch(ixOmin1:ixOmax1))
            w(ixOmin1:ixOmax1,phys_ind%tracer(ism_surround(&
               0)%myconfig%itr))=0.0_dp
          endwhere
        end if
      end if
      patch_all(ixOmin1:ixOmax1) =  patch_all(ixOmin1:ixOmax1) .and. &
         .not.sn_wdust%patch(ixOmin1:ixOmax1)
      the_dust_inuse(i_object)=sn_wdust%mydust
      i_object = i_object +1
    end if
    if(any(patch_all(ixOmin1:ixOmax1)))then
      call usr_fill_empty_region(ixImin1,ixImax1,ixOmin1,ixOmax1,0.0_dp,&
         patch_all,x,w)
    end if


  ! put dust to zero in all other zones
    cond_dust_on : if(phys_config%dust_on) then
      call dust_dummy%set_allpatch(ixImin1,ixImax1,ixOmin1,ixOmax1,&
         the_dust_inuse)
      call dust_dummy%set_w_zero(ixImin1,ixImax1,ixOmin1,ixOmax1,x,w)
      call dust_dummy%clean_memory
    end if cond_dust_on

    ! check is if initial setting is correct
    call  phys_check_w(.true., ixImin1,ixImax1, ixOmin1,ixOmax1, w, flag)

    if(any(flag(ixOmin1:ixOmax1)>0)) PRINT*,' is error',&
       maxval(flag(ixOmin1:ixOmax1)),minval(w(ixOmin1:ixOmax1,&
       phys_ind%pressure_))


    ! get conserved variables to be used in the code
    call phys_to_conserved(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x)


   call usr_clean_memory
  end subroutine initonegrid_usr
!----------------------------------------------------------------
  subroutine usr_clean_memory
    implicit none
    ! .. local ..
    integer   :: i_cloud,i_ism
    !------------------------------
        if(usrconfig%ism_on)then
          Loop_isms : do i_ism=0,usrconfig%ism_number-1
           call ism_surround(i_ism)%clean_memory
          end do Loop_isms
        end if
        if(usrconfig%cloud_on)then
          Loop_clouds : do i_cloud=0,usrconfig%cloud_number-1
           call cloud_medium(i_cloud)%clean_memory
          end do Loop_clouds
        end if
        if(usrconfig%sn_on)then
           call sn_wdust%clean_memory
        end if
  end subroutine usr_clean_memory

  subroutine specialsource_usr(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,iwmin,iwmax,&
     qtC,wCT,qt,w,x)
    use mod_dust
    implicit none

    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1, iwmin,&
       iwmax
    real(dp), intent(in)            :: qdt, qtC, qt
    real(dp), intent(in)            :: wCT(ixImin1:ixImax1,1:nw),&
        x(ixImin1:ixImax1,1:ndim)
    real(dp), intent(inout)         :: w(ixImin1:ixImax1,1:nw)
    ! .. local ..

    !----------------------------------------------------------

  end subroutine specialsource_usr
  !-------------------------------------------------------------------------
  subroutine specialbound_usr(qt,ixImin1,ixImax1,ixOmin1,ixOmax1,iB,w,x)
    ! special boundary types, user defined
    integer, intent(in)     :: ixOmin1,ixOmax1, iB, ixImin1,ixImax1
    real(dp), intent(in)    :: qt, x(ixImin1:ixImax1,1:ndim)
    real(dp), intent(inout) :: w(ixImin1:ixImax1,1:nw)
    ! .. local ..
    integer                 :: flag(ixImin1:ixImax1)
    integer                 :: i_cloud,i_ism
    logical                 :: patch_all(ixImin1:ixImax1)
    !-------------------------------------

    patch_all(ixOmin1:ixOmax1) = .true.

  ! set the ism
    if(usrconfig%ism_on)then
     Loop_isms : do i_ism=0,usrconfig%ism_number-1
      call ism_surround(i_ism)%set_w(ixImin1,ixImax1,ixOmin1,ixOmax1,qt,x,w)
      patch_all(ixOmin1:ixOmax1) =  patch_all(ixOmin1:ixOmax1) &
         .and.(.not.ism_surround(i_ism)%patch(ixOmin1:ixOmax1))
     end do Loop_isms
    end if
  ! set one cloud
    if(usrconfig%cloud_on)then
     Loop_clouds : do i_cloud=0,usrconfig%cloud_number-1
      call cloud_medium(i_cloud)%set_w(ixImin1,ixImax1,ixOmin1,ixOmax1,qt,x,w)
      patch_all(ixOmin1:ixOmax1) =  patch_all(ixOmin1:ixOmax1) &
         .and.(.not.cloud_medium(i_cloud)%patch(ixOmin1:ixOmax1))
     end do Loop_clouds
    end if


    if(any(patch_all(ixOmin1:ixOmax1)))then
     call usr_fill_empty_region(ixImin1,ixImax1,ixOmin1,ixOmax1,qt,patch_all,x,&
        w)
    end if


  ! get conserved variables to be used in the code

  call phys_to_conserved(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x)
  call usr_clean_memory


  end subroutine specialbound_usr




     !> Enforce additional refinement or coarsening
     !> One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.
     !> you must set consistent values for integers refine/coarsen:
     !> refine = -1 enforce to not refine
     !> refine =  0 doesn't enforce anything
     !> refine =  1 enforce refinement
     !> coarsen = -1 enforce to not coarsen
     !> coarsen =  0 doesn't enforce anything
     !> coarsen =  1 enforce coarsen
     !> e.g. refine for negative first coordinate x < 0 as
     !> if (any(x(ix^S,1) < zero)) refine=1
     subroutine specialrefine_usr(igrid,level,ixImin1,ixImax1,ixOmin1,ixOmax1,&
        qt,w,x,refine,coarsen)
       use mod_global_parameters
       integer, intent(in)          :: igrid, level, ixImin1,ixImax1, ixOmin1,&
          ixOmax1
       real(dp), intent(in)         :: qt, w(ixImin1:ixImax1,1:nw),&
           x(ixImin1:ixImax1,1:ndim)
       integer, intent(inout)       :: refine, coarsen
       integer                      :: level_min,level_max
       logical                      :: patch_cond
      !----------------------------------------

      ! supernovae_remnant

      cond_init_t: if(qt==0.0_dp) then
        call sn_wdust%get_patch(ixImin1,ixImax1,ixOmin1,ixOmax1,qt,x)
        if(any(sn_wdust%patch(ixOmin1:ixOmax1)))then
        refine  =  1
        coarsen = - 1

        level_min = refine_max_level/2
        level_max = refine_max_level
        patch_cond=.true.
        call user_fixrefineregion(level,level_min,level_max,patch_cond,refine,&
           coarsen)
        end if
        call sn_wdust%clean_memory
      else
        if(any(w(ixOmin1:ixOmax1,phys_ind%rho_)&
           >sn_wdust%myconfig%density_init/2.0_dp))then
         refine  =  1
         coarsen = - 1
        end if
      end if cond_init_t

      ! if(qt<wn_pulsar%myconfig%t_end_pulsar_wind.and.&
      !    qt<wn_pulsar%myconfig%t_start_pulsar_wind)return
      ! cond_pulsar_on : if(usrconfig%pulsar_on)then
      !
      ! end if cond_pulsar_on
     end subroutine specialrefine_usr
  !=====================================================================
     subroutine user_fixrefineregion(level,level_min,level_max,patch_cond,&
        refine,coarsen)
     integer, intent(in)    :: level,level_min,level_max
     logical,intent(in)     :: patch_cond
     integer, intent(inout) :: refine, coarsen
     ! .. local ..
     !-------------------------------------------
     if(patch_cond)then
      if(level>level_max)then
        coarsen = 1
        refine  = -1
      else if(level==level_max)then
        coarsen = 0
        refine  = -1
      end if
      if(level<level_min)then
        coarsen = -1
        refine  =  1
      else if(level==level_min)then
        coarsen = -1
        refine  =  0
      end if
     end if
     end subroutine user_fixrefineregion

     subroutine special_get_dt(w,ixImin1,ixImax1,ixOmin1,ixOmax1,qt,dtnew,dx1,&
        x)
       use mod_global_parameters
       integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
       double precision, intent(in)    :: dx1,qt, x(ixImin1:ixImax1,1:ndim)
       double precision, intent(in)    :: w(ixImin1:ixImax1,1:nw)
       double precision, intent(inout) :: dtnew
       !--------------------------------------------------------------

     end subroutine special_get_dt

  !> special output
  subroutine specialvar_output(ixImin1,ixImax1,ixOmin1,ixOmax1,win,x,normconv)
  ! this subroutine can be used in convert, to add auxiliary variables to the
  ! converted output file, for further analysis using tecplot, paraview, ....
  ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
  ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
  ! corresponding normalization values (default value 1)
    use mod_physics
    use mod_dust
    implicit none
    integer, intent(in)        :: ixImin1,ixImax1,ixOmin1,ixOmax1
    real(dp), intent(in)       :: x(ixImin1:ixImax1,1:ndim)
    real(dp)                   :: win(ixImin1:ixImax1,nw+nwauxio)
    real(dp)                   :: normconv(0:nw+nwauxio)
    ! .. local ..
    real(dp)                   :: w(ixImin1:ixImax1,nw)
    !----------------------------------------------------

   !w(ixI^S,1:nw)=win(ixI^S,1:nw)


  end subroutine specialvar_output


  !> this subroutine is ONLY to be used for computing auxiliary variables
  !> which happen to be non-local (like div v), and are in no way used for
  !> flux computations. As auxiliaries, they are also not advanced
  subroutine process_grid_usr(igrid,level,ixImin1,ixImax1,ixOmin1,ixOmax1,qt,w,&
     x)
    use mod_global_parameters
    implicit none
    integer, intent(in)             :: igrid,level,ixImin1,ixImax1,ixOmin1,&
       ixOmax1
    real(kind=dp)   , intent(in)    :: qt,x(ixImin1:ixImax1,1:ndim)
    real(kind=dp)   , intent(inout) :: w(ixImin1:ixImax1,1:nw)

    ! .. local ..

    !---------------------------------------------------




  end subroutine process_grid_usr
  !---------------------------------------------------------------------
  subroutine specialvarnames_output(varnames)
  ! newly added variables need to be concatenated with the w_names/primnames string
    character(len=*), intent(inout) :: varnames(:)
  !  varnames(1)  = 'deltav1'

  end subroutine specialvarnames_output
  !---------------------------------------------------------------------
  !> subroutine to fill the space regions that are not filled by the model
  subroutine usr_fill_empty_region(ixImin1,ixImax1,ixOmin1,ixOmax1,qt,&
     patchw_empty,x,w)
    use mod_global_parameters
    implicit none
    integer, intent(in)         :: ixImin1,ixImax1,ixOmin1,ixOmax1
    real(kind=dp), intent(in)   :: qt
    logical, intent(in)         :: patchw_empty(ixImin1:ixImax1)
    real(kind=dp),intent(in)    :: x(ixImin1:ixImax1,1:ndir)
    real(kind=dp),intent(inout) :: w(ixImin1:ixImax1,1:nw)
    ! .. local ..
    integer                     :: idir
    !------------------------------------------------
    where(patchw_empty(ixOmin1:ixOmax1))
      w(ixOmin1:ixOmax1,phys_ind%rho_)      = 1.0_DP
      w(ixOmin1:ixOmax1,phys_ind%pressure_)        = 1.0d-2
    end where
    Loop_idir_v : do idir=1,ndir
     where(patchw_empty(ixOmin1:ixOmax1))
      w(ixOmin1:ixOmax1,phys_ind%mom(idir)) = 0.0_dp
     end where
     if(phys_config%ismhd) then
        where(patchw_empty(ixOmin1:ixOmax1))
          w(ixOmin1:ixOmax1,phys_ind%mag(idir)) = 0.0_dp
        end where
      end if
    end do Loop_idir_v
  end subroutine usr_fill_empty_region
  !---------------------------------------------------------------------
  !> subroutine to write simulation configuration
  subroutine usr_write_setting
    implicit none
    integer,parameter   :: unit_config =12
    character(len=75)   :: filename_config
    integer             :: i_cloud,i_ism
    !-------------------------------------
    filename_config=trim(base_filename)//'.config'

    open(unit_config,file=trim(filename_config), status='replace')
    write(unit_config,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    write(unit_config,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    write(unit_config,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    write(unit_config,*)'%%%%%%%%%%% Simulation configuration %%%%%%%%%%%%'
    write(unit_config,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    write(unit_config,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    write(unit_config,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    if(usrconfig%physunit_on)call usr_physunit%write_setting(unit_config)

    write(unit_config,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    write(unit_config,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    write(unit_config,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    if(usrconfig%sn_on)call sn_wdust%write_setting(unit_config)
    if(usrconfig%ism_on)then
      Loop_isms : do i_ism=0,usrconfig%ism_number-1
       call ism_surround(i_ism)%write_setting(unit_config)
      end do Loop_isms
    end if
    write(unit_config,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    write(unit_config,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    write(unit_config,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    if(usrconfig%cloud_on)then
      Loop_clouds : do i_cloud=0,usrconfig%cloud_number-1
       call cloud_medium(i_cloud)%write_setting(unit_config)
      end do Loop_clouds
    end if
    write(unit_config,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    write(unit_config,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    write(unit_config,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    close(unit_config)

  end subroutine usr_write_setting

!----------------------------------------------------------------------
!> compute the total mass and volume in the cloud
subroutine usr_global_var
  use mod_global_parameters

end subroutine usr_global_var



end module mod_usr
