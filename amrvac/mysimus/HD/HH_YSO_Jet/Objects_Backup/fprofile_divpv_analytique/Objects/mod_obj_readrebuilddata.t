module mod_obj_readrebuilddata
  use mod_constants
  use mod_obj_global_parameters
  use mod_obj_usr_unit
  use mod_obj_mat
  use mod_physics
  implicit none


  type readrebuild_parameters
    character(len=20)     :: unit            !> rebuild physical unit at parameter file
    character(len=78)     :: obj_name        !> rebuild object name
    logical               :: normalize_done  !> rebuild is the normalisation is already done
    logical               :: normalizevariable_done
    character(len=78)     :: filename         !> rebuild starting file name
    integer               :: myindice         !> rebuild indice
    integer               :: read_npoint
    integer               :: read_ndir
    integer               :: read_ndim
    integer               :: read_nw
    character(len=78)     :: read_geo
    logical               :: read_ascii

    integer               :: build_ndir
    integer               :: build_ndim
    integer               :: i_typerebluid
    logical               :: read_prim
    logical               :: build_prim
    logical               :: read_volume_ison
    logical               :: read_level_multi
    character(len=20)     :: shape           !> rebuild zone shape
    real(kind=dp)         :: center(3)       !> rebuild zone center position (cm)
    real(kind=dp)         :: extend(3)       !> rebuild zone extentionin space (cm)
    real(kind=dp)         :: time_start
    logical               :: tracer_on
    integer               :: itr
  end type readrebuild_parameters
 type readrebuild
   type(readrebuild_parameters)   :: myconfig
   type(usrphysical_unit), pointer:: myphysunit         !> rebuild physics unit in use
   character(len=20),allocatable  :: read_var_names(:)
   character(len=20),allocatable  :: read_axis_names(:)
   integer, allocatable           :: read_npoint_dim(:)
   real(kind=dp), allocatable     :: read_w(:,:)
   real(kind=dp), allocatable     :: read_x(:,:)
   real(kind=dp), allocatable     :: read_volume(:)
   real(kind=dp), allocatable     :: read_level(:)
   real(kind=dp), allocatable     :: read_physunit(:)
   integer, allocatable           :: read_iw(:)

   logical, allocatable           :: patch(:^D&)
   logical, allocatable           :: patch_escape(:^D&)
   character(len=78)              :: subname            !> subroutine name that call it
   contains
    !PRIVATE
    PROCEDURE, PASS(self) :: set_default             => usr_readrebuild_default
    PROCEDURE, PASS(self) :: set_complet             => usr_readrebuild_set_complet
    PROCEDURE, PASS(self) :: normalize               => usr_readrebuild_normalize
    PROCEDURE, PASS(self) :: normalize_w             => usr_readrebuild_normalize_w
    PROCEDURE, PASS(self) :: read_initdatafile_ascii => usr_readrebuild_initdatafile_ascii
    PROCEDURE, PASS(self) :: read_initdatafile_binary=> usr_readrebuild_initdatafile_binary
    PROCEDURE, PASS(self) :: set_w                   => usr_readrebuild_set_w
    PROCEDURE, PASS(self) :: read_parameters         => usr_readrebuild_read_p
    PROCEDURE, PASS(self) :: write_setting           => usr_readrebuild_write_setting
    PROCEDURE, PASS(self) :: alloc_set_patch         => usr_readrebuild_alloc_set_patch
    PROCEDURE, PASS(self) :: get_patch_escape        => usr_readrebuild_get_patch_escape
    PROCEDURE, PASS(self) :: clean_memory            => usr_readrebuild_clean_memory
    PROCEDURE, PASS(self) :: clean_memory_all        => usr_readrebuild_clean_memory_all
    PROCEDURE, PASS(self) :: get_dt                  => usr_readrebuild_get_dt
    PROCEDURE, PASS(self) :: set_patch               => usr_readrebuild_set_patch
 end type readrebuild
 contains
  !--------------------------------------------------------------------
   !> Read the readrebuild parameters  from a parfile
  subroutine usr_readrebuild_read_p(self,readrebuild_config,files)
    use mod_obj_mat
    implicit none
    class(readrebuild)                         :: self
    character(len=*),intent(in)                :: files(:)
    type(readrebuild_parameters)               :: readrebuild_config

    ! .. local ..
    integer                                    :: i_file,i_error_read
    !-------------------------------------------------
    namelist /usr_readrebuild_list / readrebuild_config
    namelist /usr_readrebuild1_list/ readrebuild_config
    namelist /usr_readrebuild2_list/ readrebuild_config
    namelist /usr_readrebuild3_list/ readrebuild_config

    if(mype==0)write(*,*)'Reading usr_readrebuild_list'
    i_error_read = 0
    Loop_ifile : do i_file = 1, size(files)
       open(unitpar, file=trim(files(i_file)), status="old", iostat=i_error_read)
       select case(self%myconfig%myindice)
       case(1)
         read(unitpar, usr_readrebuild1_list, iostat=i_error_read)
       case(2)
         read(unitpar, usr_readrebuild2_list, iostat=i_error_read)
       case(3)
         read(unitpar, usr_readrebuild3_list, iostat=i_error_read)
       case default
         read(unitpar, usr_readrebuild_list, iostat=i_error_read)
       end select

       call usr_mat_read_error_message(i_error_read,&
                       self%myconfig%myindice,self%myconfig%obj_name)
       close(unitpar)
    end do Loop_ifile
  end subroutine usr_readrebuild_read_p
  !------------------------------------------------------------------------
   !------------------------------------------------------------------------
   subroutine usr_readrebuild_write_setting(self,unit_config)
     implicit none
    class(readrebuild)                          :: self
     integer,intent(in)                         :: unit_config
     ! .. local ..

     !-----------------------------------

     write(unit_config,*)'********************************************'
     write(unit_config,*)'************readrebuild setting ************'
     write(unit_config,*)'********************************************'
     write(unit_config,*) '** The configuation of the initial file **'
     write(unit_config,*) 'unit in use                      = ',  self%myconfig%unit
     write(unit_config,*) 'starting filename                = ',  self%myconfig%filename
     write(unit_config,*) 'starting type of file is ascii   = ',  self%myconfig%read_ascii
     write(unit_config,*) 'starting number dimension        = ',  self%myconfig%read_ndim
     write(unit_config,*) 'starting number direction        = ',  self%myconfig%read_ndir
     write(unit_config,*) 'starting number of variables     = ',  self%myconfig%read_nw
     write(unit_config,*) 'starting geometry                = ',  self%myconfig%read_geo
     write(unit_config,*) 'starting variables are primitive = ',  self%myconfig%read_prim
     write(unit_config,*) 'starting with volume             = ',  self%myconfig%read_volume_ison
     write(unit_config,*) 'starting with level              = ',  self%myconfig%read_level_multi

     write(unit_config,*) '** The configuation of the rebluid data **'
     write(unit_config,*) 'new number dimension             = ',  self%myconfig%build_ndim
     write(unit_config,*) 'new number direction             = ',  self%myconfig%build_ndir
     write(unit_config,*) 'rebuild zone shape               = ',  self%myconfig%shape
     write(unit_config,*) 'rebuild zone center              = ',  self%myconfig%center
     write(unit_config,*) 'rebuild zone extend              = ',  self%myconfig%extend
     write(unit_config,*) 'rebuild starting time            = ',  self%myconfig%time_start

     write(unit_config,*)'*******************************************'
     write(unit_config,*)'******** END readrebuild setting **********'
     write(unit_config,*)'*******************************************'
   end    subroutine usr_readrebuild_write_setting
   !-------------------------------------------------------------------------
  subroutine usr_readrebuild_default(self)
    implicit none
    class(readrebuild)     ::  self
    !-----------------------------------------
     self%myconfig%obj_name   = 'readrebuild'
     self%myconfig%myindice   = 0
     self%myconfig%filename   = 'none'
     self%myconfig%read_ndir  = 0
     self%myconfig%read_ndim  = 0
     self%myconfig%read_nw  = 0
     self%myconfig%read_geo   = 'slab'
     self%myconfig%read_ascii = .true.

     self%myconfig%build_ndir = 0
     self%myconfig%build_ndim = 0
     self%myconfig%read_prim  = .false.
     self%myconfig%build_prim = .false.

     self%myconfig%read_volume_ison  = .false.
     self%myconfig%read_level_multi  = .false.

     self%myconfig%center(:)         = 0.0_dp
     self%myconfig%extend(:)         = 0.0_dp
     self%myconfig%shape             = 'sphere'
     self%myconfig%time_start        = 0.0_dp
     self%myconfig%tracer_on         = .false.
     self%myconfig%itr               = 0
     self%myconfig%normalizevariable_done   = .false.
     self%myconfig%normalize_done   = .false.
  end subroutine usr_readrebuild_default

  subroutine usr_readrebuild_set_complet(self)
  implicit none
  class(readrebuild)     ::  self
  !-----------------------------------------
  ! no need for the rebuild for the moment
  end subroutine usr_readrebuild_set_complet

  subroutine usr_readrebuild_normalize(self,physunit_inuse)
    use mod_obj_usr_unit
    implicit none
    class(readrebuild)                             :: self
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
    self%myconfig%center           = self%myconfig%center        /physunit_inuse%myconfig%length
    self%myconfig%extend           = self%myconfig%extend        /physunit_inuse%myconfig%length
    self%myconfig%time_start       = self%myconfig%time_start    /physunit_inuse%myconfig%time
    self%myconfig%normalize_done   = .true.
  end subroutine usr_readrebuild_normalize
  !----------------------------------------------------------------------
  subroutine usr_readrebuild_normalize_w(self)
    use mod_obj_usr_unit
    implicit none
    class(readrebuild)             :: self
    ! .. local ..
    integer                        :: idims,iw

    !----------------------------------

    Loop_idims: do idims=1,self%myconfig%read_ndim
      self%read_x(idims,:) = self%read_x(idims,:)/self%read_physunit(idims)
    end do Loop_idims

    Loop_iw: do iw=1,self%myconfig%read_nw
      self%read_w(idims,:) = self%read_w(iw,:)/self%read_physunit(iw)
    end do Loop_iw
    self%myconfig%normalizevariable_done   = .true.

  end subroutine usr_readrebuild_normalize_w
  !----------------------------------------------------------------------
  subroutine usr_readrebuild_initdatafile_ascii(self)
    implicit none
    class(readrebuild)            :: self
    integer                       :: unit_file_read
    ! .. local ..
    integer                       :: configdata_sizes(4)
    logical                       :: config_logicals(3)
    integer, dimension(2)         :: sizes, subsizes, start
    integer                       :: MPI_type_readX_block,MPI_type_readw_block
    integer                       :: iw, read_ix
    character(len=78)             :: the_header
    integer                       :: size_character_name
    !------------------------------------------------------
    unit_file_read=12
    size_character_name = 20
    cond_open_config : if(mype==0)then
      open(unit=unit_file_read,file=self%myconfig%filename)
      !read(unit=unit_file_read,fmt="(A)")the_header
      read(unit=unit_file_read,fmt="(2(i2.1))")self%myconfig%read_ndim,self%myconfig%read_ndir
      read(unit=unit_file_read,fmt="(i4.4)")self%myconfig%read_nw
      read(unit=unit_file_read,fmt="(A)")self%myconfig%read_geo

      if(allocated(self%read_axis_names))deallocate(self%read_axis_names)
      allocate(self%read_axis_names(1:self%myconfig%read_ndim))
      if(allocated(self%read_var_names))deallocate(self%read_var_names)
      allocate(self%read_var_names(1:self%myconfig%read_nw))
      if(allocated(self%read_npoint_dim))deallocate(self%read_npoint_dim)
      allocate(self%read_npoint_dim(self%myconfig%read_ndim))
      if(allocated(self%read_physunit))deallocate(self%read_physunit)
      allocate(self%read_physunit(1:self%myconfig%read_ndim+self%myconfig%read_nw))

      read(unit=unit_file_read,fmt="(100(A16))")self%read_axis_names,self%read_var_names

      read(unit=unit_file_read,fmt="(3(L1))")self%myconfig%read_prim,self%myconfig%read_volume_ison,&
                                self%myconfig%read_level_multi
      read(unit=unit_file_read,fmt="(4(i7.7))")self%myconfig%read_npoint,self%read_npoint_dim
      read(unit=unit_file_read,fmt="(100(e14.6))")self%read_physunit
    end if cond_open_config


    cond_multicpu_0 : if (npe>1)then
      cond_master_0 : if(mype==0) then
       configdata_sizes(1)       = self%myconfig%read_ndim
       configdata_sizes(2)       = self%myconfig%read_ndir
       configdata_sizes(3)       = self%myconfig%read_nw
       configdata_sizes(4)       = self%myconfig%read_npoint
      end if cond_master_0
      call MPI_BCAST(configdata_sizes,4,MPI_INTEGER,0,icomm,ierrmpi)
      cond_slave_0 : if(mype/=0) then
       self%myconfig%read_ndim   = configdata_sizes(1)
       self%myconfig%read_ndir   = configdata_sizes(2)
       self%myconfig%read_nw     = configdata_sizes(3)
       self%myconfig%read_npoint = configdata_sizes(4)
      end if cond_slave_0
      cond_master_1 : if(mype==0) then
       config_logicals(1)        = self%myconfig%read_prim
       config_logicals(2)        = self%myconfig%read_volume_ison
       config_logicals(3)        = self%myconfig%read_level_multi
      end if cond_master_1

      call MPI_BCAST(config_logicals,3,MPI_LOGICAL,0,icomm,ierrmpi)

      cond_slave_1 : if(mype/=0) then
        self%myconfig%read_prim        = config_logicals(1)
        self%myconfig%read_volume_ison = config_logicals(2)
        self%myconfig%read_level_multi = config_logicals(3)

        if(allocated(self%read_axis_names))deallocate(self%read_axis_names)
        allocate(self%read_axis_names(self%myconfig%read_ndim))

        if(allocated(self%read_var_names))deallocate(self%read_var_names)
        allocate(self%read_var_names(self%myconfig%read_nw))

        if(allocated(self%read_physunit))deallocate(self%read_physunit)
        allocate(self%read_physunit(self%myconfig%read_nw+self%myconfig%read_ndim))
      end if cond_slave_1

      call MPI_BCAST(self%read_axis_names,&
          size_character_name*self%myconfig%read_ndim,&
          MPI_CHARACTER,0,icomm,ierrmpi)

      call MPI_BCAST(self%read_var_names,&
          size_character_name*self%myconfig%read_nw,&
          MPI_CHARACTER,0,icomm,ierrmpi)

      call MPI_BCAST(self%myconfig%read_geo,&
          len(trim(self%myconfig%read_geo)),&
          MPI_CHARACTER,0,icomm,ierrmpi)

      call MPI_BCAST(self%read_physunit,&
       self%myconfig%read_nw+self%myconfig%read_ndim,MPI_DOUBLE_PRECISION,&
       0,icomm,ierrmpi)
    end if cond_multicpu_0

    allocate(self%read_x(self%myconfig%read_npoint,self%myconfig%read_ndim),&
          self%read_w(self%myconfig%read_npoint,self%myconfig%read_nw))
    allocate(self%read_iw(1:self%myconfig%read_nw))
    cond_read_data : if(mype==0)then
      cond_read_var: if(.not.self%myconfig%read_volume_ison .and. .not.self%myconfig%read_level_multi)then
        Loop_ix_read_0 :  do read_ix = 1,self%myconfig%read_npoint
          read(unit=unit_file_read,fmt="(100(e14.6))")self%read_x(read_ix,1:self%myconfig%read_ndim),&
                                                     self%read_w(read_ix,1:self%myconfig%read_nw)
        end do  Loop_ix_read_0
      elseif(self%myconfig%read_volume_ison .and. .not.self%myconfig%read_level_multi)then
        allocate(self%read_volume(self%myconfig%read_npoint))
        Loop_ix_read_1 :  do read_ix = 1,self%myconfig%read_npoint
          read(unit=unit_file_read,fmt="(100(e14.6))")self%read_x(read_ix,1:self%myconfig%read_ndim),&
                                                     self%read_w(read_ix,1:self%myconfig%read_nw),&
                                                     self%read_volume(read_ix)
        end do  Loop_ix_read_1

      elseif(.not.self%myconfig%read_volume_ison .and. self%myconfig%read_level_multi)then
        allocate(self%read_level(self%myconfig%read_npoint))
        Loop_ix_read_2 :  do read_ix = 1,self%myconfig%read_npoint
          read(unit=unit_file_read,fmt="(100(e14.6))")self%read_x(read_ix,1:self%myconfig%read_ndim),&
                                                     self%read_w(read_ix,1:self%myconfig%read_nw),&
                                                     self%read_level(read_ix)
        end do  Loop_ix_read_2
      elseif(self%myconfig%read_volume_ison .and.self%myconfig%read_level_multi)then
        allocate(self%read_volume(self%myconfig%read_npoint),&
                self%read_level(self%myconfig%read_npoint))
        Loop_ix_read_3 :  do read_ix = 1,self%myconfig%read_npoint
          read(unit=unit_file_read,fmt="(100(e14.6))")self%read_x(read_ix,1:self%myconfig%read_ndim),&
                                                     self%read_w(read_ix,1:self%myconfig%read_nw),&
                                                     self%read_volume(read_ix),&
                                                     self%read_level(read_ix)
        end do  Loop_ix_read_3
      end if cond_read_var
    end if cond_read_data
    cond_multicpu_1 : if(npe>1) then
      sizes(1)=self%myconfig%read_npoint;
      sizes(2)=self%myconfig%read_ndim
      subsizes(1)=self%myconfig%read_npoint;
      subsizes(2)=self%myconfig%read_ndim
      start(1)=0;
      start(2)=0
      call MPI_TYPE_CREATE_SUBARRAY(2,sizes,subsizes,start, &
                                    MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION, &
                                    MPI_type_readX_block,ierrmpi)
      call MPI_TYPE_COMMIT(MPI_type_readX_block,ierrmpi)

      sizes(1)=self%myconfig%read_npoint;
      sizes(2)=self%myconfig%read_nw
      subsizes(1)=self%myconfig%read_npoint;
      subsizes(2)=self%myconfig%read_nw
      start(1)=0;
      start(2)=0
      call MPI_TYPE_CREATE_SUBARRAY(2,sizes,subsizes,start, &
                                  MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION, &
                                  MPI_type_readw_block,ierrmpi)
      call MPI_TYPE_COMMIT(MPI_type_readw_block,ierrmpi)

      call MPI_BCAST(self%read_x,1,MPI_type_readx_block,0,icomm,ierrmpi)
      call MPI_BCAST(self%read_w,1,MPI_type_readw_block,0,icomm,ierrmpi)
    end if cond_multicpu_1

    self%read_iw(phys_ind%rho_)   = phys_ind%rho_
    Loop_imom : do iw = 1,self%myconfig%read_ndir
      self%read_iw(1+iw) = phys_ind%mom(1)+iw-1
    end do Loop_imom
    if(self%myconfig%read_nw>self%myconfig%read_ndir+1.and.phys_config%nw>ndir+1) then
     Loop_ienergy : do iw = self%myconfig%read_ndir+2,self%myconfig%read_nw
      self%read_iw(iw) = self%read_iw(self%myconfig%read_ndir+1)+iw-(self%myconfig%read_ndir+1)
    end do Loop_ienergy
    end if
    call self%normalize_w()
  end subroutine usr_readrebuild_initdatafile_ascii

!----------------------------------------------------------------------
subroutine usr_readrebuild_initdatafile_binary(self)
  implicit none
  class(readrebuild)            :: self
  integer                       :: unit_file_read
  ! .. local ..
  integer                       :: configdata_sizes(4)
  logical                       :: config_logicals(3)
  integer, dimension(2)         :: sizes, subsizes, start
  integer                       :: MPI_type_readX_block,MPI_type_readw_block
  integer                       :: iw,read_ix
  !------------------------------------------------------
  unit_file_read=12
  cond_open_config : if(mype==0)then
    open(unit=unit_file_read,file=self%myconfig%filename,form='unformatted')
    read(unit=unit_file_read)self%myconfig%read_ndim,self%myconfig%read_ndir
    read(unit=unit_file_read)self%myconfig%read_nw
    read(unit=unit_file_read)self%myconfig%read_geo

    if(allocated(self%read_axis_names))deallocate(self%read_axis_names)
    allocate(self%read_axis_names(1:self%myconfig%read_ndim))
    if(allocated(self%read_var_names))deallocate(self%read_var_names)
    allocate(self%read_var_names(1:self%myconfig%read_nw))
    if(allocated(self%read_npoint_dim))deallocate(self%read_npoint_dim)
    allocate(self%read_npoint_dim(self%myconfig%read_ndim))

    read(unit=unit_file_read)self%read_axis_names,self%read_var_names

    read(unit=unit_file_read)self%myconfig%read_prim,self%myconfig%read_volume_ison,&
                              self%myconfig%read_level_multi
    read(unit=unit_file_read)self%myconfig%read_npoint,self%read_npoint_dim

  end if cond_open_config

  cond_multicpu_0 : if (npe>1)then
    configdata_sizes(1)       = self%myconfig%read_ndim
    configdata_sizes(2)       = self%myconfig%read_ndir
    configdata_sizes(3)       = self%myconfig%read_nw
    configdata_sizes(4)       = self%myconfig%read_npoint

    call MPI_BCAST(configdata_sizes,4,MPI_INTEGER,0,icomm,ierrmpi)

    self%myconfig%read_ndim   = configdata_sizes(1)
    self%myconfig%read_ndir   = configdata_sizes(2)
    self%myconfig%read_nw     = configdata_sizes(3)
    self%myconfig%read_npoint = configdata_sizes(4)

    config_logicals(1)        = self%myconfig%read_prim
    config_logicals(2)        = self%myconfig%read_volume_ison
    config_logicals(3)        = self%myconfig%read_level_multi
    call MPI_BCAST(config_logicals,3,MPI_LOGICAL,0,icomm,ierrmpi)
    self%myconfig%read_prim        = config_logicals(1)
    self%myconfig%read_volume_ison = config_logicals(2)
    self%myconfig%read_level_multi = config_logicals(3)
  end if cond_multicpu_0

  allocate(self%read_x(self%myconfig%read_npoint,self%myconfig%read_ndim),&
        self%read_w(self%myconfig%read_npoint,self%myconfig%read_nw))
  allocate(self%read_iw(1:self%myconfig%read_nw))

  cond_read_data : if(mype==0)then
    cond_read_var: if(.not.self%myconfig%read_volume_ison .and. .not.self%myconfig%read_level_multi)then
      Loop_ix_read_0 :  do read_ix = 1,self%myconfig%read_npoint
        read(unit=unit_file_read)self%read_x(read_ix,1:self%myconfig%read_ndim),&
                                                   self%read_w(read_ix,1:self%myconfig%read_nw)
      end do  Loop_ix_read_0
    elseif(self%myconfig%read_volume_ison .and. .not.self%myconfig%read_level_multi)then
      allocate(self%read_volume(self%myconfig%read_npoint))
      Loop_ix_read_1 :  do read_ix = 1,self%myconfig%read_npoint
        read(unit=unit_file_read)self%read_x(read_ix,1:self%myconfig%read_ndim),&
                                                   self%read_w(read_ix,1:self%myconfig%read_nw),&
                                                   self%read_volume(read_ix)
      end do  Loop_ix_read_1

    elseif(.not.self%myconfig%read_volume_ison .and. self%myconfig%read_level_multi)then
      allocate(self%read_level(self%myconfig%read_npoint))
      Loop_ix_read_2 :  do read_ix = 1,self%myconfig%read_npoint
        read(unit=unit_file_read)self%read_x(read_ix,1:self%myconfig%read_ndim),&
                                                   self%read_w(read_ix,1:self%myconfig%read_nw),&
                                                   self%read_level(read_ix)
      end do  Loop_ix_read_2
    elseif(self%myconfig%read_volume_ison .and.self%myconfig%read_level_multi)then
      allocate(self%read_volume(self%myconfig%read_npoint),&
              self%read_level(self%myconfig%read_npoint))
      Loop_ix_read_3 :  do read_ix = 1,self%myconfig%read_npoint
        read(unit=unit_file_read)self%read_x(read_ix,1:self%myconfig%read_ndim),&
                                                   self%read_w(read_ix,1:self%myconfig%read_nw),&
                                                   self%read_volume(read_ix),&
                                                   self%read_level(read_ix)
      end do  Loop_ix_read_3
    end if cond_read_var
  end if cond_read_data
  cond_multicpu_1 : if(npe>1) then
    sizes(1)=self%myconfig%read_npoint;
    sizes(2)=self%myconfig%read_ndim
    subsizes(1)=self%myconfig%read_npoint;
    subsizes(2)=self%myconfig%read_ndim
    start(1)=0;
    start(2)=0
    call MPI_TYPE_CREATE_SUBARRAY(2,sizes,subsizes,start, &
                                  MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION, &
                                  MPI_type_readX_block,ierrmpi)
    call MPI_TYPE_COMMIT(MPI_type_readX_block,ierrmpi)

    sizes(1)=self%myconfig%read_npoint;
    sizes(2)=self%myconfig%read_nw
    subsizes(1)=self%myconfig%read_npoint;
    subsizes(2)=self%myconfig%read_nw
    start(1)=0;
    start(2)=0
    call MPI_TYPE_CREATE_SUBARRAY(2,sizes,subsizes,start, &
                                MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION, &
                                MPI_type_readw_block,ierrmpi)
    call MPI_TYPE_COMMIT(MPI_type_readw_block,ierrmpi)

    call MPI_BCAST(self%read_x,1,MPI_type_readx_block,0,icomm,ierrmpi)
    call MPI_BCAST(self%read_w,1,MPI_type_readw_block,0,icomm,ierrmpi)
  end if cond_multicpu_1

  self%read_iw(phys_ind%rho_)   = phys_ind%rho_
  Loop_imom : do iw = 1,self%myconfig%read_ndir
    self%read_iw(1+iw) = phys_ind%mom(1)+iw
  end do Loop_imom
  cond_energy_in : if(self%myconfig%read_nw>self%myconfig%read_ndir+1&
                     .and.phys_config%nw>ndir+1) then
   Loop_ienergy : do iw = self%myconfig%read_ndir+2,self%myconfig%read_nw
    self%read_iw(iw) = self%read_iw(self%myconfig%read_ndir+1)+iw
   end do Loop_ienergy
 end if cond_energy_in
 call self%normalize_w()

end subroutine usr_readrebuild_initdatafile_binary

!--------------------------------------------------------------------
 !> subroutine patch for the cloud
 subroutine usr_readrebuild_set_patch(ixI^L,ixO^L,qt,x,self,force_refine,dx_loc)
  implicit none
  integer, intent(in)        :: ixI^L,ixO^L
  real(kind=dp), intent(in)  :: qt
  real(kind=dp), intent(in)  :: x(ixI^S,1:ndim)
  class(readrebuild)         :: self
  integer, optional          :: force_refine
  real(kind=dp),optional     :: dx_loc(1:ndim)
  ! .. local ..
  integer                    :: idims
  real(dp)                   :: x_edge(ixI^S,1:ndim)
  real(dp)                   :: min_dist,r_in,r_out
  real(dp), dimension(ixI^S) :: dist,dist_edge
  !----------------------------------

  if(allocated(self%patch))deallocate(self%patch)
  allocate(self%patch(ixI^S))
  ! E = 0.5 rho*v2 int (R^4/r_out^2 4 pi *dR)=0.5/5 *4 pi *rho*v^2*r_out^3 = 0.5  rho*v2*Volume *3/5=Ek*3/5
  call usr_distance(ixI^L,ixO^L,typeaxial,self%myconfig%center,x,dist)

  ! check distance
  r_in =self%myconfig%center(r_)
  r_out=self%myconfig%extend(r_)



  self%patch(ixO^S) = Dist(ixO^S) >r_in &
                      .and. Dist(ixO^S) <r_out!+min_dist

  if(allocated(self%patch_escape))then
    self%patch(ixO^S) = self%patch(ixO^S).and.(.not.self%patch_escape(ixO^S))
  end if
end subroutine usr_readrebuild_set_patch
!--------------------------------------------------------------------
subroutine usr_readrebuild_set_w(ixI^L,ixO^L,qt,x,w,self,interpw)
    implicit none
    integer, intent(in)           :: ixI^L,ixO^L
    real(kind=dp), intent(in)     :: qt
    real(kind=dp), intent(in)     :: x(ixI^S,1:ndim)
    real(kind=dp), intent(inout)  :: w(ixI^S,1:nw)
    class(readrebuild)            :: self
    logical, optional             :: interpw(1:nw)
    ! .. local ..
    integer                       :: ix^D,iw
    logical                       :: interpw_loc(1:nw)
    integer                       :: rightnei(1),leftnei(1)
    !--------------------------------------------------------------------
    if(present(interpw))then
      interpw_loc = interpw
    else
      interpw_loc = .true.
    end if
    call self%set_patch(ixI^L,ixO^L,qt,x)
    cond_read_oned : if(self%myconfig%read_ndim==1)then
      Loop_1D : do ix1=ixOmin1,ixOmax1
        cond_patch : if(self%patch(ix1,ixOmin^DE))then
          rightnei=minloc(dabs(self%read_x(1:self%myconfig%read_npoint,1)-x(ix1,ixOmin^DE,1))&
                           ,self%read_x(1:self%myconfig%read_npoint,1)>=x(ix1,ixOmin^DE,1))
          leftnei=minloc(dabs(self%read_x(1:self%myconfig%read_npoint,1)-x(ix1,ixOmin^DE,1))&
                             ,self%read_x(1:self%myconfig%read_npoint,1)<=x(ix1,ixOmin^DE,1))

          if (leftnei(1)==self%myconfig%read_npoint) rightnei(1)= leftnei(1)
          if (rightnei(1)==1)                        leftnei(1) = rightnei(1)

          cond_i_same: if(leftnei(1)==rightnei(1))then
            Loop_iw_1D0 : do iw =1,phys_config%nw
              if (interpw_loc(iw).and.self%read_iw(iw)>0) then
                 if(self%patch(ix1,ixOmin^DE))w(ix1,ixO^SE,iw)=self%read_w(leftnei(1),self%read_iw(iw))
              end if
            end do Loop_iw_1D0
          else cond_i_same
            cond_x_same :if(dabs(self%read_x(rightnei(1),1)-self%read_x(leftnei(1),1) )&
                 >smalldouble) then
              if(self%patch(ix1,ixOmin^DE)) then
               Loop_iw_1D : do iw =1,phys_config%nw
                if (interpw_loc(iw).and.self%read_iw(iw)>0) then
                 

                 w(ix1,ixO^SE,iw)=(self%read_w(rightnei(1),self%read_iw(iw))&
                               -self%read_w(leftnei(1),self%read_iw(iw)))&
                               /(self%read_x(rightnei(1),1)-self%read_x(leftnei(1),1))&
                               *(x(ix1,ixO^SE,1)-self%read_x(leftnei(1),1))&
                               +self%read_w(leftnei(1),self%read_iw(iw))
                 end if
               end do Loop_iw_1D
              end if
              else if(dabs(self%read_x(rightnei(1),1)-self%read_x(leftnei(1),1))&
                    <smalldouble) then cond_x_same
              Loop_iw_1D2 : do iw =1,phys_config%nw
               if (interpw_loc(iw).and.self%read_iw(iw)>0) then
                 if(self%patch(ix1,ixOmin^DE))w(ix1,ixO^SE,iw)=self%read_w(leftnei(1),self%read_iw(iw))
               end if
             end do Loop_iw_1D2
            end if cond_x_same
          end if cond_i_same
        end if cond_patch
      end do Loop_1D
    end if cond_read_oned
  end subroutine usr_readrebuild_set_w

 !--------------------------------------------------------------------
 !> Subroutine to clean array memory of associated with cloud object
 subroutine usr_readrebuild_alloc_set_patch(ixI^L,ixO^L,qt,x,self,&
                                            use_tracer,w,patch_escape,patch_usr)
   implicit none
   integer, intent(in)                     :: ixI^L,ixO^L
   real(kind=dp), intent(in)               :: qt
   real(kind=dp), intent(in)               :: x(ixI^S,1:ndim)
   real(kind=dp), intent(in), optional     :: w(ixI^S,1:nw)
   logical, intent(in), optional           :: use_tracer
   logical, intent(in),optional            :: patch_escape(ixI^S)
   logical, intent(in),optional            :: patch_usr(ixI^S)
  class(readrebuild)                              :: self
   !---------------------------------------------------------
   cond_tracer : if(.not.present(use_tracer).and. .not.present(w)) then
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
     case('usr')
       if(present(patch_usr))then
        self%patch(ixO^S) = patch_usr(ixO^S)
       else
        write(*,*)'this rebuild zone shape ','usr',' is not implimented in mod_obj_readrebuilddata.t',&
                   ' at subroutine  usr_readrebuild_alloc_set_patch'
        call mpistop('This cloud shape is not implimented in mod_obj_readrebuilddata.t')
       end if
     case default
        write(*,*)'this rebuild zone shape ',trim(self%myconfig%shape),' is not implimented'
        call mpistop('This cloud shape is not implimented in usr_readrebuild_alloc_set_patch')
     end select

   else cond_tracer

     cond_readrebuildtracer_on : if(self%myconfig%tracer_on)then
      if(allocated(self%patch))deallocate(self%patch)
      allocate(self%patch(ixI^S))
      where(w(ixO^S,phys_ind%tracer(self%myconfig%itr))>small_density)
        self%patch(ixO^S)=.true.
      else where
        self%patch(ixO^S)=.false.
      end where
     end if cond_readrebuildtracer_on

   end if cond_tracer


   if(present(patch_escape))then
    call self%get_patch_escape(ixI^L,ixO^L,.true.,patch_escape)
   else
    if(allocated(self%patch_escape))deallocate(self%patch_escape)
    allocate(self%patch_escape(ixI^S))
    self%patch_escape(ixO^S) =.false.
   end if

   where(self%patch(ixO^S))self%patch(ixO^S)         = .not.self%patch_escape(ixO^S)
 end subroutine usr_readrebuild_alloc_set_patch
!---------------------------------------------------------------------
 subroutine  usr_readrebuild_get_patch_escape(ixI^L,ixO^L,need_dealloc,patch_escape,self)
   implicit none
   integer, intent(in)           :: ixI^L,ixO^L
   logical, intent(in)           :: need_dealloc
   logical, intent(in)           :: patch_escape(ixI^S)
  class(readrebuild)                    :: self
   !----------------------------------------------------------
   if(allocated(self%patch_escape))deallocate(self%patch_escape)
    allocate(self%patch_escape(ixI^S))
    self%patch_escape(ixO^S) = .false.

   self%patch_escape(ixO^S)=self%patch_escape(ixO^S).or.patch_escape(ixO^S)
 end subroutine  usr_readrebuild_get_patch_escape


 !===========================================================
 !> Subroutine to set time set before wind starts and stops
 subroutine usr_readrebuild_get_dt(self,ixI^L,ixO^L,dx^D,x,w,qt,dtnew)
   use mod_global_parameters
   class(readrebuild)              :: self
   integer, intent(in)             :: ixI^L, ixO^L
   double precision, intent(in)    :: dx^D,qt, x(ixI^S,1:ndim)
   double precision, intent(in)    :: w(ixI^S,1:nw)
   double precision, intent(inout) :: dtnew
   !--------------------------------------------------------------
   if(qt<self%myconfig%time_start)then
    dtnew=min(self%myconfig%time_start-qt,dtnew)
   end if
 end subroutine usr_readrebuild_get_dt
 !===========================================================
  subroutine usr_readrebuild_clean_memory(self)
    implicit none
    class(readrebuild) ::   self
    !-----------------------------------------------------------
    if(allocated(self%patch))deallocate(self%patch)
    if(allocated(self%patch_escape))deallocate(self%patch_escape)
  end subroutine usr_readrebuild_clean_memory

  subroutine usr_readrebuild_clean_memory_all(self)
    implicit none
    class(readrebuild) ::   self
    !-----------------------------------------------------------
    if(allocated(self%patch))deallocate(self%patch)
    if(allocated(self%patch_escape))deallocate(self%patch_escape)
    if(allocated(self%read_w))deallocate(self%read_w)
    if(allocated(self%read_x))deallocate(self%read_x)
    if(allocated(self%read_volume))deallocate(self%read_volume)
    if(allocated(self%read_level))deallocate(self%read_level)
    if(allocated(self%read_var_names))deallocate(self%read_var_names)
    if(allocated(self%read_axis_names))deallocate(self%read_axis_names)
  end subroutine usr_readrebuild_clean_memory_all
end module mod_obj_readrebuilddata
