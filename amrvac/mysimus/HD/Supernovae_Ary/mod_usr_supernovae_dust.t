module mod_usr
  use mod_hd, only : hd_activate
  use mod_dust
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
      physics_type='hd'
      call usr_physunit%set_complet(physics_type)
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

       ism_surround(i_ism)%myconfig        = ism_default%myconfig
       ism_surround(i_ism)%mydust%myconfig = ism_default%mydust%myconfig

       ism_surround(i_ism)%myconfig%myindice=i_ism
       call ism_surround(i_ism)%read_parameters(ism_surround(i_ism)%myconfig,files)
      end do Loop_allism
    end if

    if(usrconfig%cloud_on)then
      allocate(cloud_medium(0:usrconfig%cloud_number-1))
      Loop_allcloud : do i_cloud =0,usrconfig%ism_number
       cloud_medium(i_cloud)%myconfig          = cloud_default%myconfig
       cloud_medium(i_cloud)%mydust%myconfig   = cloud_default%mydust%myconfig
       cloud_medium(i_cloud)%myconfig%myindice = i_cloud

       call cloud_medium(i_cloud)%read_parameters(files,cloud_medium(i_cloud)%myconfig)
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

   namelist /usr_unit_list/ unit_length , unit_time,unit_velocity,          &
                      unit_density, unit_numberdensity,                     &
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
      if(usrconfig%sn_on)then
        sn_wdust%myconfig%dust_on =.false.
      end if
    end if  cond_dust_on
  end   subroutine usr_check_conflict
  !-----------------------------------------------------------
  !> subroutine to normalize parameters in the code
  subroutine usr_normalise_parameters
   implicit none
   integer            :: idust



   constusr%G         = constusr%G*&
                      (unit_density*(unit_length/unit_velocity)**(2.0_dp))


   constusr%clight                      = constusr%clight/unit_velocity

   ! complet all physical unit in use
   if(usrconfig%physunit_on) then
      call usr_physunit%fillphysunit
    end if


    w_convert_factor(phys_ind%rho_)              = unit_density
    if(phys_config%energy)w_convert_factor(phys_ind%e_)= unit_density*unit_velocity**2.0
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
             w_convert_factor(phys_ind%dust_mom(:,idust))     = unit_density*unit_velocity
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

   if(phys_config%dust_on)allocate(the_dust_inuse(n_objects))
   call usr_normalise_parameters
   if(mype==0)call usr_write_setting

  end subroutine initglobaldata_usr
  !> The initial conditions
  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    ! initialize one grid

    implicit none

    integer, intent(in)     :: ixI^L,ixO^L
    real(dp), intent(in)    :: x(ixI^S,1:ndim)
    real(dp), intent(inout) :: w(ixI^S,1:nw)
    !.. local ..
    real(dp)      :: res
    integer       :: ix^D,na,flag(ixI^S)
    integer       :: i_cloud,i_ism,i_dust,i_start,i_end
    logical, save :: first=.true.
    logical       :: patch_all(ixI^S)
    type(dust)    :: dust_dummy
    integer       :: i_object
    ! .. only test ..
    real(dp)      ::old_w(ixO^S,1:nw)
    !-----------------------------------------
    patch_all(ixO^S) = .true.
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
       call ism_surround(i_ism)%set_w(ixI^L,ixO^L,global_time,x,w)
       patch_all(ixO^S) =  patch_all(ixO^S) .and. .not.ism_surround(i_ism)%patch(ixO^S)
       if(ism_surround(i_ism)%myconfig%dust_on)the_dust_inuse(i_object)=ism_surround(i_ism)%mydust
       i_object = i_object +1
      end do Loop_isms
    end if


    ! set one cloud
    if(usrconfig%cloud_on)then
      Loop_clouds : do i_cloud=0,usrconfig%cloud_number-1
       call cloud_medium(i_cloud)%set_w(ixI^L,ixO^L,global_time,x,w)

       if(usrconfig%ism_on)then
         if(ism_surround(0)%myconfig%tracer_on)then
           where(cloud_medium(i_cloud)%patch(ixO^S))
             w(ixO^S,phys_ind%tracer(ism_surround(0)%myconfig%itr))=0.0_dp
           end where

         end if
       end if

       patch_all(ixO^S) =  patch_all(ixO^S) .and. .not.cloud_medium(i_cloud)%patch(ixO^S)
       if(cloud_medium(i_cloud)%myconfig%dust_on)the_dust_inuse(i_object)=cloud_medium(i_cloud)%mydust
       i_object = i_object +1
      end do Loop_clouds
    end if

    ! set the pulsar and associated wind + envelope if they are on
    if(usrconfig%sn_on)then
      sn_wdust%subname='initonegrid_usr'
      call sn_wdust%set_w(ixI^L,ixO^L,global_time,x,w)

      if(usrconfig%ism_on)then

        if(ism_surround(0)%myconfig%tracer_on)then
          where(sn_wdust%patch(ixO^S))
            w(ixO^S,phys_ind%tracer(ism_surround(0)%myconfig%itr))=0.0_dp
          endwhere
        end if
        Loop_isms_sn : do i_ism=0,usrconfig%ism_number-1
          if(ism_surround(i_ism)%myconfig%dust_on) then
            ism_surround(i_ism)%mydust%patch(ixO^S)=.not.sn_wdust%patch(ixO^S)
            i_start= ism_surround(i_ism)%mydust%myconfig%idust_first
            i_end  = ism_surround(i_ism)%mydust%myconfig%idust_last
            Loop_idust_1:  do i_dust=i_start,i_end
              ism_surround(i_ism)%mydust%the_ispecies(i_dust)%patch(ixO^S)=.not.sn_wdust%patch(ixO^S)
            end do Loop_idust_1
            call ism_surround(i_ism)%mydust%set_w_zero(ixI^L,ixO^L,x,w)
          end if
        end do   Loop_isms_sn
      end if

      patch_all(ixO^S) =  patch_all(ixO^S) .and. .not.sn_wdust%patch(ixO^S)
      if(sn_wdust%myconfig%dust_on)the_dust_inuse(i_object)=sn_wdust%mydust
      i_object = i_object +1
    end if



    if(any(patch_all(ixO^S)))then
      call usr_fill_empty_region(ixI^L,ixO^L,0.0_dp,patch_all,x,w)
    end if


  ! put dust to zero in all other zones
    cond_dust_on : if(phys_config%dust_on) then
      dust_dummy%myconfig%idust_first = 1
      dust_dummy%myconfig%idust_last  = dust_n_species
      call dust_dummy%set_allpatch(ixI^L,ixO^L,the_dust_inuse)
      call dust_dummy%set_w_zero(ixI^L,ixO^L,x,w)
      call dust_dummy%clean_memory
    end if cond_dust_on

    ! check is if initial setting is correct
    call  phys_check_w(.true., ixI^L, ixO^L, w, flag)

    if(any(flag(ixO^S)>0)) PRINT*,' is error',maxval(flag(ixO^S)),minval(w(ixO^S,phys_ind%pressure_))


    ! get conserved variables to be used in the code
    call phys_to_conserved(ixI^L,ixO^L,w,x)


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

  subroutine specialsource_usr(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
    use mod_dust
    implicit none

    integer, intent(in)             :: ixI^L, ixO^L, iw^LIM
    real(dp), intent(in)            :: qdt, qtC, qt
    real(dp), intent(in)            :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    real(dp), intent(inout)         :: w(ixI^S,1:nw)
    ! .. local ..

    !----------------------------------------------------------

  end subroutine specialsource_usr
  !-------------------------------------------------------------------------
  subroutine specialbound_usr(qt,ixI^L,ixO^L,iB,w,x)
    ! special boundary types, user defined
    integer, intent(in)     :: ixO^L, iB, ixI^L
    real(dp), intent(in)    :: qt, x(ixI^S,1:ndim)
    real(dp), intent(inout) :: w(ixI^S,1:nw)
    ! .. local ..
    integer                 :: flag(ixI^S)
    integer                 :: i_cloud,i_ism
    logical                 :: patch_all(ixI^S)
    !-------------------------------------

    patch_all(ixO^S) = .true.

  ! set the ism
    if(usrconfig%ism_on)then
     Loop_isms : do i_ism=0,usrconfig%ism_number-1
      call ism_surround(i_ism)%set_w(ixI^L,ixO^L,qt,x,w)
      patch_all(ixO^S) =  patch_all(ixO^S) .and.(.not.ism_surround(i_ism)%patch(ixO^S))
     end do Loop_isms
    end if
  ! set one cloud
    if(usrconfig%cloud_on)then
     Loop_clouds : do i_cloud=0,usrconfig%cloud_number-1
      call cloud_medium(i_cloud)%set_w(ixI^L,ixO^L,qt,x,w)
      patch_all(ixO^S) =  patch_all(ixO^S) .and.(.not.cloud_medium(i_cloud)%patch(ixO^S))
     end do Loop_clouds
    end if


    if(any(patch_all(ixO^S)))then
     call usr_fill_empty_region(ixI^L,ixO^L,qt,patch_all,x,w)
    end if


  ! get conserved variables to be used in the code

  call phys_to_conserved(ixI^L,ixO^L,w,x)
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
     subroutine specialrefine_usr(igrid,level,ixI^L,ixO^L,qt,w,x,refine,coarsen)
       use mod_global_parameters
       integer, intent(in)          :: igrid, level, ixI^L, ixO^L
       real(dp), intent(in)         :: qt, w(ixI^S,1:nw), x(ixI^S,1:ndim)
       integer, intent(inout)       :: refine, coarsen
       ! .. local ..
       integer                      :: level_min,level_max,level_need
       logical                      :: patch_cond
       real(dp)                     :: dx_loc(1:ndim)
      !----------------------------------------

      ! supernovae_remnant
      cond_init_t: if(dabs(qt-0.0_dp)<smalldouble) then
        ^D&dx_loc(^D)=rnode(rpdx^D_,igrid);
        call sn_wdust%get_patch(ixI^L,ixO^L,qt,x,force_refine=1,dx_loc=dx_loc)
        if(any(sn_wdust%patch(ixO^S)))then
         level_need = nint(dlog((xprobmax1/sn_wdust%myconfig%r_out)/domain_nx1)/dlog(2.0_dp))
         level_min = level_need+1
         level_max = refine_max_level
         patch_cond=.true.
         call user_fixrefineregion(level,level_min,level_max,patch_cond,refine,coarsen)
        else
         refine  = -1
         coarsen = 1
        end if
        call sn_wdust%clean_memory
      end if cond_init_t

      ! if(qt<wn_pulsar%myconfig%t_end_pulsar_wind.and.&
      !    qt<wn_pulsar%myconfig%t_start_pulsar_wind)return
      ! cond_pulsar_on : if(usrconfig%pulsar_on)then
      !
      ! end if cond_pulsar_on
     end subroutine specialrefine_usr
  !=====================================================================
     subroutine user_fixrefineregion(level,level_min,level_max,patch_cond,refine,coarsen)
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

     subroutine special_get_dt(w,ixI^L,ixO^L,qt,dtnew,dx^D,x)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L
       double precision, intent(in)    :: dx^D,qt, x(ixI^S,1:ndim)
       double precision, intent(in)    :: w(ixI^S,1:nw)
       double precision, intent(inout) :: dtnew
       !--------------------------------------------------------------

     end subroutine special_get_dt

  !> special output
  subroutine specialvar_output(ixI^L,ixO^L,win,x,normconv)
  ! this subroutine can be used in convert, to add auxiliary variables to the
  ! converted output file, for further analysis using tecplot, paraview, ....
  ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
  ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
  ! corresponding normalization values (default value 1)
    use mod_physics
    use mod_dust
    implicit none
    integer, intent(in)        :: ixI^L,ixO^L
    real(dp), intent(in)       :: x(ixI^S,1:ndim)
    real(dp)                   :: win(ixI^S,nw+nwauxio)
    real(dp)                   :: normconv(0:nw+nwauxio)
    ! .. local ..
    real(dp)                   :: w(ixI^S,nw)
    real(dp)                   :: error_var(ixM^T)
    integer                    :: iw, level,idir,idust
    integer, parameter         :: iz_                = 1
    integer, parameter         :: ilevel_            = 2
    integer, parameter         :: ierror_lohner_rho_ = 3
    integer, parameter         :: ierror_lohner_p_   = 4
    integer, parameter         :: idragfroce11_      = 5
    integer, parameter         :: idragfroce12_      = 6
    integer, parameter         :: idragfroce13_      = 7
    integer, parameter         :: idragfroce14_      = 8
    integer, parameter         :: idv11_             = 9
    integer, parameter         :: idv12_             = 10
    integer, parameter         :: idv13_             = 11
    integer, parameter         :: idv14_             = 12
    integer, parameter         :: idtgas1_           = 13
    integer, parameter         :: idtdust11_         = 14
    integer, parameter         :: idtdust12_         = 15
    integer, parameter         :: idtdust13_         = 16
    integer, parameter         :: idtdust14_         = 17
    integer, parameter         :: istincking1_       = 18
    integer, parameter         :: istincking2_       = 19
    integer, parameter         :: istincking3_       = 20
    integer, parameter         :: istincking4_       = 21
    integer, parameter         :: idtgas2_           = 22
    integer, parameter         :: idtdust2_          = 23
    real(dp)                   :: z_ism,z_sn
    real(dp), dimension(ixO^S) :: fraction_sn,rho_sn,rho_ism
    real(dp)                   :: epsilon_ism,epsilon_sn,dx^D
    logical,dimension(ixI^S)          :: fd_flag
    real(dp),dimension(ixI^S)         :: ptherm,cmax
    real(dp),dimension(ixI^S,1:ndir)  :: vgas
    real(dp),dimension(ixI^S,1:phys_config%dust_n_species) ::alpha_T
    real(dp),dimension(ixI^S,1:ndir,1:phys_config%dust_n_species) :: vdust,fdrag
    !----------------------------------------------------
    w(ixI^S,1:nw) = win(ixI^S,1:nw)
    level = node(plevel_,saveigrid)
    call phys_get_pthermal(w, x, ixI^L, ixI^L, ptherm)
    call dust_get_sticking(w, x, ixI^L, ixO^L, alpha_T, ptherm)
    Loop_idir1 : do idir=1,ndir
      vgas(ixI^S,idir)=w(ixI^S,phys_ind%mom(idir))/w(ixI^S,phys_ind%rho_)
      Loop_idust1 :  do idust = 1,phys_config%dust_n_species
        where(w(ixI^S,phys_ind%dust_rho(idust))>phys_config%dust_small_density)
         vdust(ixI^S,idir,idust)=w(ixI^S,phys_ind%dust_mom(idir, idust))&
                              /w(ixI^S,phys_ind%dust_rho(idust))
        elsewhere
          vdust(ixI^S,idir,idust)=0.0_dp
        endwhere
      end do Loop_idust1
    end do Loop_idir1
    dx^D=rnode(rpdx^D_,saveigrid);
    Loop_iw :  do iw = 1,nwauxio
    select case(iw)
    case(iz_)
      epsilon_ism = 0.3_dp
      epsilon_sn  = 6.0_dp
      z_ism = constusr%z_solar*epsilon_ism
      z_sn  = constusr%z_solar*epsilon_sn

      if(phys_config%n_tracer>0)then
       fraction_sn = w(ixO^S,phys_ind%tracer(sn_wdust%myconfig%itr))&
          /SUM(w(ixO^S,phys_ind%tracer(1):phys_ind%tracer(phys_config%n_tracer)))
        rho_sn =   fraction_sn*w(ixO^S,phys_ind%rho_)
        rho_ism=(1.0_dp-fraction_sn)*w(ixO^S,phys_ind%rho_)
        win(ixO^S,nw+iz_) = (z_ism*rho_ism+z_sn*rho_sn)/(rho_sn+rho_ism)
      end if
    case(ilevel_)
      normconv(nw+ilevel_)     = 1.0_dp
      win(ixO^S,nw+ilevel_) = node(plevel_,saveigrid)
     case(ierror_lohner_rho_)
       normconv(nw+iw)     = 1.0_dp
       win(ixG^T,nw+ierror_lohner_rho_) = 0.0
       call usr_mat_get_Lohner_error(ixI^L, ixM^LL,level,phys_ind%rho_,w,error_var)
       win(ixM^T,nw+ierror_lohner_rho_) = error_var(ixM^T)
     case(ierror_lohner_p_)
       normconv(nw+iw)     = 1.0_dp
       win(ixG^T,nw+ierror_lohner_p_) = 0.0_dp
       call usr_mat_get_Lohner_error(ixI^L, ixM^LL,level,phys_ind%pressure_,w,error_var)
       win(ixM^T,nw+ierror_lohner_p_) = error_var(ixM^T)
    case(idragfroce11_)
       normconv(nw+idragfroce11_)     = 1.0_dp
       call get_3d_dragforce(ixI^L, ixO^L, w, x, fdrag, ptherm, vgas,fd_flag)
       win(ixO^S,nw+idragfroce11_) = fdrag(ixO^S,1,1)
    case(idragfroce12_)
       normconv(nw+idragfroce12_)     = 1.0_dp
       call get_3d_dragforce(ixI^L, ixO^L, w, x, fdrag, ptherm, vgas,fd_flag)
       win(ixO^S,nw+idragfroce12_) = fdrag(ixO^S,1,2)
    case(idragfroce13_)
       normconv(nw+idragfroce13_)     = 1.0_dp
       call get_3d_dragforce(ixI^L, ixO^L, w, x, fdrag, ptherm, vgas,fd_flag)
       win(ixO^S,nw+idragfroce13_) = fdrag(ixO^S,1,3)
    case(idragfroce14_)
       normconv(nw+idragfroce14_)     = 1.0_dp
       call get_3d_dragforce(ixI^L, ixO^L, w, x, fdrag, ptherm, vgas,fd_flag)
       win(ixO^S,nw+idragfroce14_) = fdrag(ixO^S,1,4)
    case(idv11_)
      normconv(nw+idv11_)            = 1.0_dp
      where(w(ixO^S,phys_ind%dust_rho(1))>phys_config%dust_small_density)
       win(ixO^S,nw+idv11_)          =  vgas(ixO^S,1)-vdust(ixO^S,1,1)
      elsewhere
       win(ixO^S,nw+idv11_)          = 0.0_dp
      end where
    case(idv12_)
      normconv(nw+idv12_)            = 1.0_dp
      where(w(ixO^S,phys_ind%dust_rho(1))>phys_config%dust_small_density)
       win(ixO^S,nw+idv12_)          =  vgas(ixO^S,1)-vdust(ixO^S,1,2)
      elsewhere
       win(ixO^S,nw+idv12_)          = 0.0_dp
      end where
    case(idv13_)
      normconv(nw+idv13_)            = 1.0_dp
      where(w(ixO^S,phys_ind%dust_rho(1))>phys_config%dust_small_density)
       win(ixO^S,nw+idv13_)          =  vgas(ixO^S,1)-vdust(ixO^S,1,3)
      elsewhere
       win(ixO^S,nw+idv13_)          = 0.0_dp
      end where
    case(idv14_)
      normconv(nw+idv14_)            = 1.0_dp
      where(w(ixO^S,phys_ind%dust_rho(1))>phys_config%dust_small_density)
       win(ixO^S,nw+idv14_)          =  vgas(ixO^S,1)-vdust(ixO^S,1,4)
      elsewhere
       win(ixO^S,nw+idv14_)          = 0.0_dp
      end where
    case(idtgas1_)
      normconv(nw+idtgas1_)            = unit_time
      call phys_get_cmax(w,x,ixI^L,ixO^L,1,cmax)
      win(ixO^S,nw+idtgas1_) = dx1/cmax(ixO^S)
    case(idtdust11_)
      normconv(nw+idtdust11_)            = unit_time
      call dust_get_tstop(ixI^L, ixO^L, 1,1,dx^D, x,w,cmax)
      win(ixO^S,nw+idtdust11_) = cmax(ixO^S)
    case(idtdust12_)
      normconv(nw+idtdust12_)            = unit_time
      call dust_get_tstop(ixI^L, ixO^L, 2,1,dx^D, x,w,cmax)
      win(ixO^S,nw+idtdust12_) = cmax(ixO^S)
    case(idtdust13_)
      normconv(nw+idtdust13_)            = unit_time
      call dust_get_tstop(ixI^L, ixO^L, 3,1,dx^D, x,w,cmax)
      win(ixO^S,nw+idtdust13_) = cmax(ixO^S)
    case(idtdust14_)
      normconv(nw+idtdust14_)            = unit_time
      call dust_get_tstop(ixI^L, ixO^L, 4,1,dx^D, x,w,cmax)
      win(ixO^S,nw+idtdust14_) = cmax(ixO^S)
    case(istincking1_)
      normconv(nw+istincking1_)            = 1.0_dp
      win(ixO^S,nw+istincking1_) = alpha_T(ixO^S,1)
    case(istincking2_)
      normconv(nw+istincking2_)            = 1.0_dp
      win(ixO^S,nw+istincking2_) = alpha_T(ixO^S,2)
    case(istincking3_)
      normconv(nw+istincking3_)            = 1.0_dp
      win(ixO^S,nw+istincking3_) = alpha_T(ixO^S,3)
    case(istincking4_)
      normconv(nw+istincking4_)            = 1.0_dp
      win(ixO^S,nw+istincking4_) = alpha_T(ixO^S,4)
    case(idtgas2_)
      normconv(nw+idtgas2_)            = unit_time
      call phys_get_cmax(w,x,ixI^L,ixO^L,2,cmax)
      win(ixO^S,nw+idtgas2_) = cmax(ixO^S)
    case(idtdust2_)
      normconv(nw+idtdust2_)            = unit_time
      call phys_get_cmax(w,x,ixI^L,ixO^L,2,cmax)
      win(ixO^S,nw+idtdust2_) = cmax(ixO^S)
    case default
       write(*,*)'is not implimented at specialvar_output in mod_user'
    end select
  end do Loop_iw
    !----------------------------------------------------

   !w(ixI^S,1:nw)=win(ixI^S,1:nw)


  end subroutine specialvar_output


  !> this subroutine is ONLY to be used for computing auxiliary variables
  !> which happen to be non-local (like div v), and are in no way used for
  !> flux computations. As auxiliaries, they are also not advanced
  subroutine process_grid_usr(igrid,level,ixI^L,ixO^L,qt,w,x)
    use mod_global_parameters
    implicit none
    integer, intent(in)             :: igrid,level,ixI^L,ixO^L
    real(kind=dp)   , intent(in)    :: qt,x(ixI^S,1:ndim)
    real(kind=dp)   , intent(inout) :: w(ixI^S,1:nw)

    ! .. local ..
    integer                         :: idust,idir
    real(dp)                        :: small_dust_rho,coef
    logical, dimension(ixI^S)       :: patch_correct,patch_slow
    !---------------------------------------------------


    if(phys_config%dust_on) then
     small_dust_rho = sn_wdust%mydust%myconfig%min_limit_rel

     call phys_to_primitive(ixI^L,ixI^L,w,x)
     ! handel small density dust
     Loop_idust : do idust =1, dust_n_species
      where(w(ixI^S, phys_ind%dust_rho(idust))<max(small_dust_rho*w(ixI^S,phys_ind%rho_),&
         ism_surround(0)%mydust%myconfig%min_limit_abs))
        w(ixI^S, phys_ind%dust_rho(idust))= 0.8* min(small_dust_rho*w(ixI^S,phys_ind%rho_),&
             ism_surround(0)%mydust%myconfig%min_limit_abs)
        patch_correct(ixI^S) = .true.
      elsewhere
        patch_correct(ixI^S) = .false.
      end where
      ! handel large density dust
      where(w(ixI^S,phys_ind%rho_)<0.9*ism_surround(0)%myconfig%density)
       where(w(ixI^S, phys_ind%dust_rho(idust))>ism_surround(0)%mydust%myconfig%max_limit_rel*w(ixI^S,phys_ind%rho_))
        w(ixI^S, phys_ind%dust_rho(idust))=0.8*ism_surround(0)%mydust%myconfig%max_limit_rel*w(ixI^S,phys_ind%rho_)
        patch_slow(ixI^S) = .true.
       elsewhere
        patch_slow(ixI^S) =.false.
       end where
      end where
      ! handel large cmax in rarefied region
      ! do idir = 1, ndim
      !   call phys_get_cmax(w,x,ixI^L,ixO^L,idir,cmax)
      !
      ! end do


    !  new_dvflag(ixI^S)=.true.
    !  new_dfflag(ixI^S)=.true.

    !  vt2(ixI^S) = 3.0d0*w(ixI^S,e_)/w(ixI^S, rho_)
      Loop_idir1 : do idir = 1,ndim
       where(patch_correct(ixI^S))
        w(ixI^S, phys_ind%dust_mom(idir,idust))=0.0_dp
       end where
       where(patch_slow(ixI^S))
               w(ixI^S, phys_ind%dust_mom(idir,idust))=w(ixI^S,phys_ind%mom(idir))
       end where
      end do   Loop_idir1
     end do Loop_idust
     call phys_to_conserved(ixI^L,ixI^L,w,x)
    end if

  end subroutine process_grid_usr
  !---------------------------------------------------------------------
  subroutine specialvarnames_output(varnames)
  ! newly added variables need to be concatenated with the w_names/primnames string
    character(len=*), intent(inout) :: varnames(:)
    ! .. local ..
    integer                         :: iw
    integer, parameter              :: iz_                = 1
    integer, parameter              :: ilevel_            = 2
    integer, parameter              :: ierror_lohner_rho_ = 3
    integer, parameter              :: ierror_lohner_p_   = 4
    integer, parameter              :: idragfroce11_      = 5
    integer, parameter              :: idragfroce12_      = 6
    integer, parameter              :: idragfroce13_      = 7
    integer, parameter              :: idragfroce14_      = 8
    integer, parameter              :: idv11_             = 9
    integer, parameter              :: idv12_             = 10
    integer, parameter              :: idv13_             = 11
    integer, parameter              :: idv14_             = 12
    integer, parameter              :: idtgas1_           = 13
    integer, parameter              :: idtdust11_         = 14
    integer, parameter              :: idtdust12_         = 15
    integer, parameter              :: idtdust13_         = 16
    integer, parameter              :: idtdust14_         = 17
    integer, parameter              :: istincking1_       = 18
    integer, parameter              :: istincking2_       = 19
    integer, parameter              :: istincking3_       = 20
    integer, parameter              :: istincking4_       = 21
    integer, parameter              :: idtgas2_           = 22
    integer, parameter              :: idtdust2_          = 23
    !----------------------------------------------------
    Loop_iw : do  iw = 1,nwauxio
    select case(iw)
    case(iz_)
      varnames(iz_)               = 'zmetalicity'
    case(ilevel_)
      varnames(ilevel_)           = 'level'
    case(ierror_lohner_rho_)
      varnames(ierror_lohner_rho_)='erroramrrho'
    case(ierror_lohner_p_)
      varnames(ierror_lohner_p_)  ='erroamrp'
    case(idragfroce11_)
      varnames(idragfroce11_)     ='dragfroce11'
    case(idragfroce12_)
      varnames(idragfroce12_)     ='dragfroce12'
    case(idragfroce13_)
      varnames(idragfroce13_)     ='dragfroce13'
    case(idragfroce14_)
      varnames(idragfroce14_)     ='dragfroce14'
    case(idv11_)
      varnames(idv11_)            ='dvdust11'
    case(idv12_)
      varnames(idv12_)            ='dvdust12'
    case(idv13_)
      varnames(idv13_)            ='dvdust13'
    case(idv14_)
      varnames(idv14_)            ='dvdust14'
    case(idtgas1_)
      varnames(idtgas1_)          ='dtgas1'
    case(idtgas2_)
      varnames(idtgas2_)          ='dtgas2'
    case(idtdust11_)
      varnames(idtdust11_)        ='dtdust11'
    case(idtdust12_)
      varnames(idtdust12_)        ='dtdust12'
    case(idtdust13_)
      varnames(idtdust13_)        ='dtdust13'
    case(idtdust14_)
      varnames(idtdust14_)        ='dtdust14'
    case(idtdust2_)
      varnames(idtdust2_)         ='dtgas2'
    case(istincking1_)
      varnames(istincking1_)      ='stincking1'
    case(istincking2_)
      varnames(istincking2_)      ='stincking2'
    case(istincking3_)
      varnames(istincking3_)      ='stincking3'
    case(istincking4_)
      varnames(istincking4_)      ='stincking4'
    end select
    end do Loop_iw



  end subroutine specialvarnames_output
  !---------------------------------------------------------------------
  !> subroutine to fill the space regions that are not filled by the model
  subroutine usr_fill_empty_region(ixI^L,ixO^L,qt,patchw_empty,x,w)
    use mod_global_parameters
    implicit none
    integer, intent(in)         :: ixI^L,ixO^L
    real(kind=dp), intent(in)   :: qt
    logical, intent(in)         :: patchw_empty(ixI^S)
    real(kind=dp),intent(in)    :: x(ixI^S,1:ndir)
    real(kind=dp),intent(inout) :: w(ixI^S,1:nw)
    ! .. local ..
    integer                     :: idir
    !------------------------------------------------
    where(patchw_empty(ixO^S))
      w(ixO^S,phys_ind%rho_)      = 1.0_DP
      w(ixO^S,phys_ind%pressure_)        = 1.0d-2
    end where
    Loop_idir_v : do idir=1,ndir
     where(patchw_empty(ixO^S))
      w(ixO^S,phys_ind%mom(idir)) = 0.0_dp
     end where
     if(phys_config%ismhd) then
        where(patchw_empty(ixO^S))
          w(ixO^S,phys_ind%mag(idir)) = 0.0_dp
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
