module mod_usr
  use mod_srmhd_parameters
  use mod_srmhd
  use mod_physics
  use mod_dust
  use mod_global_parameters
  use mod_obj_global_parameters
  use mod_obj_mat
  use mod_obj_cloud
  use mod_obj_ism
  use mod_obj_sn_remnant
  use mod_obj_relativistic_wind
  use mod_obj_star_envelope
  use mod_obj_pulsar
  use mod_obj_star
  use mod_obj_usr_unit
  implicit none
  save
  real(dp) :: theta, kx, ly, vc

  type usr_config
    logical           :: physunit_on
    logical           :: ism_on
    logical           :: cloud_on
    logical           :: pulsar_on
    logical           :: pulsar_wind_on
    logical           :: pulsar_envelope_on
    logical           :: ism_list_diff
    logical           :: cloud_list_diff
    integer           :: cloud_number,ism_number
    real(kind=dp)     :: coarsen_power_slop
    real(kind=dp)     :: coarsen_distance
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

  type (pulsar), target  :: wn_pulsar
  type (rel_wind), target:: wn_wind

  type(star) :: star_ms
  type(star) :: sun

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
    usr_internal_bc     => usr_special_internal_bc
    usr_get_dt          => special_get_dt

    call usr_set_default_parameters



    call usr_physunit%set_default

    ! set default values for pulsar configuration
    call wn_pulsar%set_default

    ! set default values for ISMs configuration
    call ism_default%set_default


    ! set default values for clouds configuration
    call cloud_default%set_default



    call usr_params_read(par_files)



    ! complet all physical unit in use
    if(usrconfig%physunit_on) then
     physics_type='srmhd'
     call usr_physunit%set_complet(physics_type)
    end if
    call usr_physical_unit
    call set_coordinate_system(trim(usrconfig%coordinate_system))
    call srmhd_activate

    phys_n_tracer=srmhd_n_tracer
    call usr_check_conflict


  end subroutine usr_init
  !------------------------------------------------------------------
  !> default usr parameters from a file
  subroutine usr_set_default_parameters
    !-------------------------------------
    usrconfig%physunit_on         = .false.
    usrconfig%ism_on              = .false.
    usrconfig%cloud_on            = .false.
    usrconfig%pulsar_on           = .false.
    usrconfig%pulsar_wind_on      = .false.
    usrconfig%pulsar_envelope_on  = .false.
    usrconfig%cloud_number        = 1
    usrconfig%ism_number          = 1
    usrconfig%ism_list_diff       = .false.
    usrconfig%cloud_list_diff     = .false.
    usrconfig%coarsen_distance    = 0
    usrconfig%coarsen_power_slop  = 0
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

    if(usrconfig%pulsar_on)call wn_pulsar%read_parameters(wn_pulsar%myconfig,files)

    if(usrconfig%ism_on)then
      allocate(ism_surround(0:usrconfig%ism_number-1))
      Loop_allism : do i_ism =0,usrconfig%ism_number-1
       ism_surround(i_ism)%myconfig = ism_default%myconfig
       ism_surround(i_ism)%myconfig%myindice=i_ism
       call ism_surround(i_ism)%read_parameters(ism_surround(i_ism)%myconfig,files)
      end do Loop_allism
    end if

    if(usrconfig%cloud_on)then
      allocate(cloud_medium(0:usrconfig%cloud_number-1))
      Loop_allcloud : do i_cloud =0,usrconfig%ism_number
       cloud_medium(i_cloud)%myconfig          = cloud_default%myconfig
       cloud_medium(i_cloud)%myconfig%myindice = i_cloud

       call cloud_medium(i_cloud)%read_parameters(files,cloud_medium(i_cloud)%myconfig)
      end do Loop_allcloud
    end if

  end subroutine usr_params_read


  subroutine usr_clean_memory_final
    if(usrconfig%ism_on)then
      if(allocated(ism_surround))deallocate(ism_surround)
    end if
    if(usrconfig%cloud_on)then
      if(allocated(cloud_medium))deallocate(cloud_medium)
    end if
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
    if(srmhd_energy)w_convert_factor(phys_ind%e_)   = unit_density*unit_velocity**2.0
    if(saveprim)then
     w_convert_factor(phys_ind%mom(:))           = unit_velocity
    else
     w_convert_factor(phys_ind%mom(:))           = unit_density*unit_velocity
    end if
    time_convert_factor                 = unit_time
    length_convert_factor               = unit_length




    if(srmhd_dust)then
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
   integer   :: i_cloud,i_ism
   !------------------------------------

    itr=1
   ! complet ism parameters
   if(usrconfig%ism_on)then
    Loop_isms : do i_ism=0,usrconfig%ism_number-1
     ism_surround(i_ism)%myconfig%itr=itr
     call ism_surround(i_ism)%set_complet
     call ism_surround(i_ism)%normalize(usr_physunit)
    end do Loop_isms
    itr=ism_surround(usrconfig%ism_number-1)%myconfig%itr+1
   end if



   ! complet cloud parameters
   if(usrconfig%cloud_on) then
     Loop_clouds : do i_cloud=0,usrconfig%cloud_number-1
      cloud_medium(i_cloud)%myconfig%itr=itr
      call cloud_medium(i_cloud)%set_complet
      call cloud_medium(i_cloud)%normalize(usr_physunit)
     end do Loop_clouds
     itr=cloud_medium(usrconfig%cloud_number-1)%myconfig%itr+1
   end if


   if(usrconfig%pulsar_on)then
     wn_pulsar%myconfig%itr=itr
     call wn_pulsar%set_complet
     call wn_pulsar%normalize(usr_physunit)
   end if

   call usr_normalise_parameters
   if(mype==0)call usr_write_setting

  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    ! initialize one grid

    implicit none

    integer, intent(in)     :: ixI^L,ixO^L
    real(dp), intent(in)    :: x(ixI^S,1:ndim)
    real(dp), intent(inout) :: w(ixI^S,1:nw)
    !.. local ..
    real(dp)      :: res
    integer       :: ix^D,na,flag(ixI^S)
    integer       :: i_cloud,i_ism
    logical, save :: first=.true.
    logical       :: patch_all(ixI^S)
    real(dp)::old_w(ixO^S,1:nw)

    !-----------------------------------------
    patch_all(ixO^S) = .true.
    if(first)then
      if(mype==0) then
        write(*,*)"----------------------------------------------------------"
        write(*,*)'pulsar wind start :-)'
        write(*,*)"----------------------------------------------------------"
      endif
      first=.false.
    endif


    ! set the ism
    if(usrconfig%ism_on) then
      Loop_isms : do i_ism=0,usrconfig%ism_number-1
       call ism_surround(i_ism)%set_w(ixI^L,ixO^L,global_time,x,w)
       patch_all(ixO^S) =  patch_all(ixO^S) .and. .not.ism_surround(i_ism)%patch(ixO^S)
      end do Loop_isms
    end if
    ! set one cloud
    if(usrconfig%cloud_on)then
      Loop_clouds : do i_cloud=0,usrconfig%cloud_number-1
       call cloud_medium(i_cloud)%set_w(ixI^L,ixO^L,global_time,x,w)

       if(usrconfig%ism_on)then
         if(ism_surround(0)%myconfig%tracer_on)then
           where(cloud_medium(i_cloud)%patch(ixO^S))w(ixO^S,tracer(ism_surround(0)%myconfig%itr))=0.0_dp
         end if
       end if

       patch_all(ixO^S) =  patch_all(ixO^S) .and. .not.cloud_medium(i_cloud)%patch(ixO^S)
      end do Loop_clouds
    end if
    ! set the pulsar and associated wind + envelope if they are on
    if(usrconfig%pulsar_on)then
      wn_pulsar%subname='initonegrid_usr'
      call wn_pulsar%set_w(ixI^L,ixO^L,global_time,x,w)
      if(usrconfig%ism_on)then
        if(ism_surround(0)%myconfig%tracer_on)then
          where(wn_pulsar%patch(ixO^S))w(ixO^S,tracer(ism_surround(0)%myconfig%itr))=0.0_dp
        end if
      end if
      patch_all(ixO^S) =  patch_all(ixO^S) .and. .not.wn_pulsar%patch(ixO^S)
    end if

    if(any(patch_all(ixO^S)))then
      call usr_fill_empty_region(ixI^L,ixO^L,0.0_dp,patch_all,x,w)
    end if
  !  w(ixO^S,p_) = wn_pulsar%myenvelope%myconfig%pressure_init

    call srmhd_get_4u_from_3v(ixI^L,ixO^L,w(ixI^S,mom(1):mom(ndir)),&
                              w(ixI^S,lfac_))

    ! check is if initial setting is correct
    call  phys_check_w(.true., ixI^L, ixO^L, w, flag)

    if(any(flag(ixO^S)/=0)) PRINT*,' is error',maxval(abs(flag(ixO^S))),minval(w(ixO^S,p_))


    ! get conserved variables to be used in the code
    call phys_to_conserved(ixI^L,ixO^L,w,x)
    old_w(ixO^S,1:nw)=w(ixO^S,1:nw)
    call phys_to_primitive(ixI^L,ixO^L,w,x)
    if(any(dabs(old_w(ixO^S,1:nw)-w(ixO^S,1:nw))>0.0_dp))PRINT*,' difference at consreve primitive',maxval(dabs(old_w(ixO^S,1:nw)-w(ixO^S,1:nw))),&
   maxloc(dabs(old_w(ixO^S,1:nw)-w(ixO^S,1:nw)))
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
        if(usrconfig%pulsar_on)call wn_pulsar%clean_memory
  end subroutine usr_clean_memory

  subroutine specialsource_usr(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
    use mod_dust
    implicit none

    integer, intent(in)             :: ixI^L, ixO^L, iw^LIM
    real(dp), intent(in)            :: qdt, qtC, qt
    real(dp), intent(in)            :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    real(dp), intent(inout)         :: w(ixI^S,1:nw)
    ! .. local ..
    integer                         :: i_tracer
    integer                         :: iw,i_ism
    integer                         :: flag(ixI^S)
    logical, dimension(ixI^S)       :: patch_all,patch
    real(kind=dp)                   :: r_inner_boundary
    real(kind=dp)                   :: Dist(ixI^S)
    real(kind=dp)                   :: w_inboundary(ixI^S,1:nw)
    logical                         :: set_variation
    !----------------------------------------------------------
    ! better controle of
    cond_tracer_on : if(phys_config%n_tracer>0) then
     Loop_tracer : do i_tracer =1,phys_config%n_tracer
      where(w(ixO^S,phys_ind%tracer(i_tracer))>2d2)
         w(ixO^S,phys_ind%tracer(i_tracer)) = 1d2
      elsewhere(w(ixO^S,phys_ind%tracer(i_tracer))<1e-2)
         w(ixO^S,phys_ind%tracer(i_tracer)) = 0.0_dp
      end where

     end do Loop_tracer
     patch(ixO^S) = .true.
     if(usrconfig%ism_on)then
      Loop_isms : do i_ism=0,usrconfig%ism_number-1
        where(w(ixO^S,phys_ind%tracer(ism_surround(i_ism)%myconfig%itr))>smalldouble)
          patch(ixO^S) = .false.
        end where
      end do Loop_isms
     end if
     where(w(ixO^S,phys_ind%tracer(wn_pulsar%mywind%myconfig%itr))>smalldouble)
      patch(ixO^S) = .false.
     endwhere
     where(patch(ixO^S))
      w(ixO^S,phys_ind%tracer(wn_pulsar%mysupernovae_remnant%myconfig%itr)) = 1d2
     end where
    end if cond_tracer_on


    set_variation = .false.
    cond_wind : if(dabs(qtC-wn_pulsar%myconfig%t_start_pulsar_wind)<smalldouble) then
      call usr_distance(ixI^L,ixO^L,typeaxial,wn_pulsar%mywind%myconfig%center&
                        ,x,dist)
      patch(ixO^S) = Dist(ixO^S) <=   wn_pulsar%mywind%myconfig%r_out_init.and. &
                     Dist(ixO^S) <   wn_pulsar%mywind%myconfig%r_in
      set_variation = .true.
    end if cond_wind
    cond_rebuild : if(dabs(qtC-wn_pulsar%myconfig%t_rebuild)<smalldouble) then
      call usr_distance(ixI^L,ixO^L,typeaxial,wn_pulsar%myreadrebuild%myconfig%center&
                        ,x,dist)
      patch(ixO^S) = Dist(ixO^S) <=   wn_pulsar%myreadrebuild%myconfig%extend(r_).and. &
                     Dist(ixO^S) <   wn_pulsar%myreadrebuild%myconfig%center(r_)
      set_variation = .true.
    end if cond_rebuild
    cond_variation : if(set_variation) then
      cond_inside_wind : if(any(patch(ixO^S))) then

        w_inboundary = w

        if(any(patch(ixO^S))) then
         call usr_fill_empty_region(ixI^L,ixO^L,qt,patch(ixO^S),x,w_inboundary)
         call phys_to_primitive(ixI^L,ixO^L,w_inboundary,x)
         call srmhd_get_3v_from_4u(ixI^L,ixO^L,w_inboundary(ixI^S,lfac_),&
                                    w_inboundary(ixI^S,mom(1):mom(ndir)))
        end if

        wn_pulsar%subname='specialsource_usr'
        call wn_pulsar%set_w(ixI^L,ixO^L,qtC,x,w_inboundary)
        patch_all(ixO^S) =  patch_all(ixO^S)&
                            .and.(.not.wn_pulsar%patch(ixO^S).and.patch(ixO^S))



        if(any(patch_all(ixO^S)))then
          call usr_fill_empty_region(ixI^L,ixO^L,qt,patch_all,x,w_inboundary)
        end if
        call srmhd_get_4u_from_3v(ixI^L,ixO^L,&
                                    w_inboundary(ixI^S,mom(1):mom(ndir)),&
                                    w_inboundary(ixI^S,lfac_))

        ! check is if initial setting is correct
        call  phys_check_w(.true., ixI^L, ixO^L, w_inboundary, flag)

        ! get conserved variables to be used in the code

        call phys_to_conserved(ixI^L,ixO^L,w_inboundary,x)

        Loop_iw_inbound :  do iw = 1,nw
         where(patch(ixO^S))
          w(ixO^S,iw) = w_inboundary(ixO^S,iw)
         endwhere
        end do Loop_iw_inbound

        call wn_pulsar%clean_memory

      end if   cond_inside_wind
    end if cond_variation
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
    if(usrconfig%pulsar_on)then
      wn_pulsar%subname='specialbound_usr'
      call wn_pulsar%set_w(ixI^L,ixO^L,qt,x,w)
      patch_all(ixO^S) =  patch_all(ixO^S) .and.(.not.wn_pulsar%patch(ixO^S))
    end if

    if(any(patch_all(ixO^S)))then
     call usr_fill_empty_region(ixI^L,ixO^L,qt,patch_all,x,w)
    end if

  call srmhd_get_4u_from_3v(ixI^L,ixO^L,w(ixI^S,mom(1):mom(ndir)),&
                             w(ixI^S,lfac_))
  ! get conserved variables to be used in the code

  call phys_to_conserved(ixI^L,ixO^L,w,x)
  call usr_clean_memory


  end subroutine specialbound_usr


   !> This subroutine can be used to artificially overwrite ALL conservative
   !> variables in a user-selected region of the mesh, and thereby act as
   !> an internal boundary region. It is called just before external (ghost cell)
   !> boundary regions will be set by the BC selection. Here, you could e.g.
   !> want to introduce an extra variable (nwextra, to be distinguished from nwaux)
   !> which can be used to identify the internal boundary region location.
   !> Its effect should always be local as it acts on the mesh.
   subroutine usr_special_internal_bc(level,qt,ixI^L,ixO^L,w,x)
     use mod_global_parameters
     integer, intent(in)             :: ixI^L,ixO^L,level
     real(kind=dp)   , intent(in)    :: qt
     real(kind=dp)   , intent(inout) :: w(ixI^S,1:nw)
     real(kind=dp)   , intent(in)    :: x(ixI^S,1:ndim)
     ! .. local ..
     integer                  :: iw
     integer                  :: flag(ixI^S)
     logical, dimension(ixI^S):: patch_all,patch
     real(kind=dp)            :: r_inner_boundary
     real(kind=dp)            :: Dist(ixI^S)
     real(kind=dp)            :: w_inboundary(ixI^S,1:nw)
     logical,save             :: sn_start          = .false.
     logical,save             :: pulsar_wind_start = .false.
    !-------------------------------------
    if(qt<wn_pulsar%myconfig%t_end_pulsar_wind.and.&
       qt<wn_pulsar%myconfig%t_start_pulsar_wind)return
    cond_pulsar_on : if(usrconfig%pulsar_on)then

      cond_wind_on : if(wn_pulsar%myconfig%wind_on.or.&
                        wn_pulsar%myconfig%supernovae_remnant_on) then
        patch_all(ixO^S) = .true.
        cond_wind : if(qt>=wn_pulsar%myconfig%t_start_pulsar_wind .and. &
           qt<wn_pulsar%myconfig%t_end_pulsar_wind) then
          if(.not.pulsar_wind_start)then
            write(*,*) 'the pulsar wind is on'
            write(*,*) 'physical time : ',qt*usr_physunit%myconfig%time
          end if
          pulsar_wind_start=.true.
          wn_pulsar%mywind%subname='specialbound_usr'
          call wn_pulsar%mywind%get_patch(ixI^L,ixO^L,qt,x,r_out_limit=r_inner_boundary)
          !r_inner_boundary = r_wind
          call wn_pulsar%clean_memory
        end if cond_wind

        cond_sn : if(qt>=wn_pulsar%myconfig%t_sn .and. &
                     qt<wn_pulsar%myconfig%t_start_pulsar_wind ) then

          if(.not.sn_start)then
            write(*,*) 'the supernovae is on'
            write(*,*) 'physical time : ',qt*usr_physunit%myconfig%time
          end if
          sn_start=.true.
          r_inner_boundary   =  0.8_dp*wn_pulsar%mysupernovae_remnant%myconfig%r_in
        end if cond_sn

        call usr_distance(ixI^L,ixO^L,typeaxial,wn_pulsar%mywind%myconfig%center&
                          ,x,dist)
        patch(ixO^S) = Dist(ixO^S) <   r_inner_boundary

        cond_inside_wind : if(any(patch(ixO^S))) then

          w_inboundary = w

          if(any(patch(ixO^S))) then
           call usr_fill_empty_region(ixI^L,ixO^L,qt,patch(ixO^S),x,w_inboundary)
           call phys_to_primitive(ixI^L,ixO^L,w_inboundary,x)
           call srmhd_get_3v_from_4u(ixI^L,ixO^L,w_inboundary(ixI^S,lfac_),&
                                      w_inboundary(ixI^S,mom(1):mom(ndir)))
          end if

          wn_pulsar%subname='usr_special_internal_bc'
          call wn_pulsar%set_w(ixI^L,ixO^L,qt,x,w_inboundary)
          patch_all(ixO^S) =  patch_all(ixO^S)&
                              .and.(.not.wn_pulsar%patch(ixO^S).and.patch(ixO^S))



          if(any(patch_all(ixO^S)))then
            call usr_fill_empty_region(ixI^L,ixO^L,qt,patch_all,x,w_inboundary)
          end if
          call srmhd_get_4u_from_3v(ixI^L,ixO^L,&
                                      w_inboundary(ixI^S,mom(1):mom(ndir)),&
                                      w_inboundary(ixI^S,lfac_))

          ! check is if initial setting is correct
          call  phys_check_w(.true., ixI^L, ixO^L, w_inboundary, flag)

          ! get conserved variables to be used in the code

          call phys_to_conserved(ixI^L,ixO^L,w_inboundary,x)

          Loop_iw_inbound :  do iw = 1,nw
           where(patch(ixO^S))
            w(ixO^S,iw) = w_inboundary(ixO^S,iw)
           endwhere
          end do Loop_iw_inbound

          call wn_pulsar%clean_memory

        end if   cond_inside_wind
      end if cond_wind_on
    end if cond_pulsar_on


   end subroutine usr_special_internal_bc

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
       integer                      :: level_min,level_max,level_need,idir
       logical                      :: patch_cond
       real(kind=dp)                :: coef_coarsen
       real(kind=dp)                :: dist(ixI^S)
       real(kind=dp)                :: dx_loc(1:ndim)
       real(kind=dp)                :: dv(ixI^S),v_idim(ixI^S)
      !----------------------------------------

      ! at the initial condition
      cond_init_t: if(dabs(qt-wn_pulsar%myconfig%t_sn)<smalldouble) then
        ^D&dx_loc(^D)=rnode(rpdx^D_,igrid);
        wn_pulsar%mysupernovae_remnant%subname='specialrefine_usr first'
        call wn_pulsar%mysupernovae_remnant%get_patch(ixI^L,ixO^L,qt,x,force_refine=1,dx_loc=dx_loc)

        if(any(wn_pulsar%mysupernovae_remnant%patch(ixO^S)))then
         level_need= nint(dlog((xprobmax1/wn_pulsar%mysupernovae_remnant%myconfig%r_out)/domain_nx1)/dlog(2.0_dp))
         level_min = max(wn_pulsar%mysupernovae_remnant%myconfig%refine_max_level*2/3,level_need)
         level_max = wn_pulsar%mysupernovae_remnant%myconfig%refine_max_level
         patch_cond=.true.
         call user_fixrefineregion(level,level_min,level_max,patch_cond,refine,coarsen)
         else
         refine =-1
         coarsen= 1
        end if
        call wn_pulsar%mysupernovae_remnant%clean_memory
      else cond_init_t
      ! supernovae_remnant

       call usr_distance(ixI^L,ixO^L,typeaxial,wn_pulsar%mywind%myconfig%center&
                           ,x,dist)


       cond_pulsar : if(any(w(ixO^S,lfac_)>1.5.and.&
                        w(ixO^S,phys_ind%tracer(wn_pulsar%mywind%myconfig%itr))&
                        >0.01.and.&
                        Dist(ixO^S)>=wn_pulsar%mywind%myconfig%r_out_impos*0.9_dp)) then
         call usr_mat_profile_tanh_scalar_maxdist(dist(ixOmax^D),&
                     wn_pulsar%mywind%myconfig%coarsen_distance&
                     ,wn_pulsar%mywind%myconfig%coarsen_var_distance,coef_coarsen)


         level_min = nint(wn_pulsar%mywind%myconfig%refine_min_level * coef_coarsen)
         level_max = nint(wn_pulsar%mywind%myconfig%refine_max_level * coef_coarsen)

         patch_cond=.true.
         call user_fixrefineregion(level,level_min,level_max,patch_cond,&
                                   refine,coarsen)
       elseif(all(Dist(ixO^S)<wn_pulsar%mywind%myconfig%r_out_impos*0.9_dp .and.&
                  w(ixO^S,lfac_)>wn_pulsar%mywind%myconfig%lfac/2.and.&
                  w(ixO^S,phys_ind%tracer(wn_pulsar%mywind%myconfig%itr))&
                  >0.1)) then cond_pulsar
         call usr_mat_profile_tanh_scalar_maxdist(maxval(dist(ixO^S),&
                         dist(ixO^S)<wn_pulsar%mywind%myconfig%r_out_impos),&
                     wn_pulsar%mywind%myconfig%coarsen_distance&
                     ,wn_pulsar%mywind%myconfig%coarsen_var_distance,coef_coarsen)

         level_min = 1
         level_max = nint( wn_pulsar%mywind%myconfig%refine_max_level*coef_coarsen*2.0_dp/3.0_dp)
         patch_cond=.true.
         call user_fixrefineregion(level,level_min,level_max,patch_cond,refine,coarsen)
       elseif(all(w(ixO^S,phys_ind%tracer(wn_pulsar%mywind%myconfig%itr))&
                  <smalldouble.and.&
                  w(ixO^S,phys_ind%tracer(wn_pulsar%mysupernovae_remnant%myconfig%itr))&
                  <smalldouble.and.w(ixO^S,lfac_)-1.0_dp<smalldouble))then
        refine  =-1
        coarsen = 1
       else
         call phys_get_v_idim(w,x,ixI^L,ixO^L,1,v_idim)
         dv(ixO^S)=(v_idim(ixO^S)&
         -wn_pulsar%mysupernovae_remnant%myconfig%velocity_proper(1))**2.0_dp
         if(ndir>1) then
           Loop_idir : do idir=2,ndir
             call phys_get_v_idim(w,x,ixI^L,ixO^L,idir,v_idim)
             dv(ixO^S)=dv(ixO^S)+(v_idim(ixO^S)&
       -wn_pulsar%mysupernovae_remnant%myconfig%velocity_proper(idir))**2.0_dp
           end do Loop_idir
         end if

         cond_free_ism : if(all(dv(ixO^S)<smalldouble))then
           refine  = -1!
           coarsen = 1
         else cond_free_ism ! not
         call usr_mat_profile_tanh_scalar_maxdist(dist(ixOmax^D),&
              wn_pulsar%mysupernovae_remnant%myconfig%coarsen_distance,&
              wn_pulsar%mysupernovae_remnant%myconfig%coarsen_var_distance,&
              coef_coarsen)


          level_min = nint(wn_pulsar%mysupernovae_remnant%myconfig%refine_min_level&
                     *coef_coarsen)
          level_max = max(nint(wn_pulsar%mysupernovae_remnant%myconfig%refine_max_level&
                     *coef_coarsen),wn_pulsar%mysupernovae_remnant%myconfig%refine_max_level/2)
          patch_cond=.true.
          call user_fixrefineregion(level,level_min,level_max,patch_cond,refine,coarsen)
         end if cond_free_ism
       end if cond_pulsar

    cond_pulsar_start :if(qt>wn_pulsar%myconfig%t_start_pulsar_wind-10.0_dp*dt.and.&
                              qt<=wn_pulsar%myconfig%t_start_pulsar_wind) then
        ^D&dx_loc(^D)=rnode(rpdx^D_,igrid);
        wn_pulsar%mywind%subname='specialrefine_usr first'
        call wn_pulsar%mywind%get_patch(ixI^L,ixO^L,qt,x,force_refine=1,dx_loc=dx_loc)
        cond_pulsar_first :if(any(wn_pulsar%mywind%patch(ixO^S)))then
         level_min = wn_pulsar%mywind%myconfig%refine_max_level
         level_max = wn_pulsar%mywind%myconfig%refine_max_level
         patch_cond=.true.
         call user_fixrefineregion(level,level_min,level_max,patch_cond,&
                                   refine,coarsen)
        end if cond_pulsar_first
        call wn_pulsar%mywind%clean_memory
       elseif(qt>wn_pulsar%myconfig%t_start_pulsar_wind)then
        ^D&dx_loc(^D)=rnode(rpdx^D_,igrid);
        wn_pulsar%mywind%subname='specialrefine_usr second'
        call wn_pulsar%mywind%get_patch(ixI^L,ixO^L,qt,x,force_refine=-1,dx_loc=dx_loc)
        cond_pulsar_second :if(any(wn_pulsar%mywind%patch(ixO^S)))then
         level_min = wn_pulsar%mywind%myconfig%refine_max_level
         level_max = wn_pulsar%mywind%myconfig%refine_max_level
         patch_cond=.true.
         call user_fixrefineregion(level,level_min,level_max,patch_cond,&
                                   refine,coarsen)
        end if cond_pulsar_second
        call wn_pulsar%mywind%clean_memory
       end if cond_pulsar_start

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
       cond_pulsar: if(usrconfig%pulsar_on)then
         call wn_pulsar%get_dt(ixI^L,ixO^L,dx^D,x,w,qt,dtnew)
       end if cond_pulsar
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
    real(dp)                   :: w(ixI^S,1:nw)
    real(dp)                   :: error_var(ixM^T)
    integer                    :: iw,level
    integer, parameter         :: ilevel_            = 1
    integer, parameter         :: ierror_lohner_rho_ = 2
    integer, parameter         :: ierror_lohner_p_   = 3
    integer, parameter         :: inertie_           = 4
    !----------------------------------------------------
    w(ixI^S,1:nw) = win(ixI^S,1:nw)
    level = node(plevel_,saveigrid)
    Loop_iw :  do iw = 1,nwauxio
    select case(iw)
     case(ilevel_)
       win(ixO^S,nw+ilevel_) = node(plevel_,saveigrid)
     case(ierror_lohner_rho_)
       win(ixG^T,nw+ierror_lohner_rho_) = 0.0
       call usr_mat_get_Lohner_error(ixI^L, ixM^LL,level,rho_,w,error_var)
       win(ixM^T,nw+ierror_lohner_rho_) = error_var(ixM^T)
     case(ierror_lohner_p_)
       win(ixG^T,nw+ierror_lohner_p_) = 0.0_dp
       call usr_mat_get_Lohner_error(ixI^L, ixM^LL,level,p_,w,error_var)
       win(ixM^T,nw+ierror_lohner_p_) = error_var(ixM^T)
     case(inertie_)
       win(ixO^S,nw+inertie_) = w(ixO^S,lfac_)*w(ixO^S,rho_)
     case default
       write(*,*)'is not implimented at specialvar_output in mod_user'
     end select
    end do Loop_iw

  end subroutine specialvar_output



  !---------------------------------------------------------------------
  subroutine specialvarnames_output(varnames)
  ! newly added variables need to be concatenated with the w_names/primnames string
    character(len=*), intent(inout) :: varnames(:)
    integer                         :: iw
    integer, parameter              :: ilevel_            = 1
    integer, parameter              :: ierror_lohner_rho_ = 2
    integer, parameter              :: ierror_lohner_p_   = 3
    integer, parameter              :: inertie_           = 4
    !----------------------------------------------------
    Loop_iw : do  iw = 1,nwauxio
    select case(iw)
    case(ilevel_)
      varnames(ilevel_) = 'level'
    case(ierror_lohner_rho_)
      varnames(ierror_lohner_rho_) ='erroramrrho'
    case(ierror_lohner_p_)
      varnames(ierror_lohner_p_) ='erroamrp'
    case(inertie_)
      varnames(inertie_) ='inertie'
    end select
    end do Loop_iw

  end subroutine specialvarnames_output

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

      !---------------------------------------------------




    end subroutine process_grid_usr
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
      w(ixO^S,rho_)      = 1.0_DP
      w(ixO^S,p_)        = 1.0d-2
    end where
    Loop_idir_v : do idir=1,ndir
     where(patchw_empty(ixO^S))
      w(ixO^S,mom(idir)) = 0.0_dp
      w(ixO^S,mag(idir)) = 0.0_dp
     end where
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
    if(usrconfig%pulsar_on)call wn_pulsar%write_setting(unit_config)
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
