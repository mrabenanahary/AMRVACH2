!> Module for including dust species, which interact with the gas through a drag
!> force
module mod_chemical
  use mod_global_parameters, only: std_len,dp,nw,ndim,unit_temperature,unit_density,unit_time
  use mod_constants
  use mod_physics
  implicit none
  private


  integer, protected              :: gas_rho_ = -1
  integer, allocatable, protected :: gas_mom(:)
  integer, protected              :: gas_e_   = -1


  type chemical_parameters
      !> The number of dust species
      integer             :: n_species

      !> Mean molecular weight of gas molecules
      real(kind=dp)       :: gas_mu

      !> Set small chemical densities to zero to avoid numerical problems
      logical             :: small_to_zero

      !> Minimum chemical density
      real(kind=dp)       :: small_density

      !> TODO: 1. Introduce this generically in advance, 2: document
      logical             :: source_split
      character(len=20)   :: method
      real(kind=dp)       :: neutral_fraction
      real(kind=dp)       :: relative_fractionmin
      integer             :: H_
      integer             :: H2_
      integer             :: He_
  end type chemical_parameters



  type chemical_coef_reaction
    real(kind=dp) :: alpha_a0
    real(kind=dp) :: colf_a0
    real(kind=dp) :: colf_T0
  end type chemical_coef_reaction

  type(chemical_coef_reaction)           :: chemical_const
  type theion
    integer                              :: myind_
    type(chemical_coef_reaction)         :: myreccoef
  end type theion

  type associate_element
    logical                              :: react_on
    logical                              :: allowed
    type(chemical_coef_reaction)         :: myreccoef
  end type associate_element


  type chemical_element
    integer                              :: ionmin
    integer                              :: ionmax
    type(theion),allocatable             :: ion(:)
    character(len=5)                     :: elem_name
    type(associate_element),allocatable  :: asso_elem(:)
  end type chemical_element

  type thechemicalices
    type(chemical_element),allocatable       :: rho(:)
  end type thechemicalices

  type chemical_phys
    integer                                  :: nvar
    type(chemical_parameters)                :: myconfig
    type(chemical_element),allocatable       :: element(:)
  end type chemical_phys
  type(chemical_phys),target :: Thechemical
  ! Public methods
  public :: chemical_init
  public :: chemical_get_dt
  public :: chemical_get_flux
  public :: chemical_get_flux_prim
  public :: chemical_add_source
  public :: chemical_to_conserved
  public :: chemical_to_primitive
  public :: chemical_check_params
  public :: chemical_check_w
  public :: chemical_set_floor
contains

  subroutine chemical_init(phys_indices_inuse,phys_config_inuse,g_rho, g_mom, g_energy)
    use mod_global_parameters
    use mod_physics
    type(phys_variables_indices)    :: phys_indices_inuse
    type(physconfig)                :: phys_config_inuse
    integer, intent(in)             :: g_rho
    integer, intent(in)             :: g_mom(ndir)
    integer, intent(in)             :: g_energy ! Negative value if not present
    integer                         :: ichemical, idir,i_ion
    character(len=3)                :: elem

    call chemical_set_default(thechemical%myconfig)
    call chemical_read_params(par_files,thechemical%myconfig)
    allocate(gas_mom(ndir))
    gas_rho_ = g_rho
    gas_mom  = g_mom
    gas_e_   = g_energy



    allocate(thechemical%element(thechemical%myconfig%n_species))
    call chemical_elements_init(thechemical)
    allocate(phys_indices_inuse%chemical_element(thechemical%nvar))
      ! Set index of chemical densities
    call chemical_set_physvariable(thechemical,phys_indices_inuse)

    phys_config_inuse%chemical_n_species     =  thechemical%myconfig%n_species
    phys_config_inuse%chemical_small_density = thechemical%myconfig%small_density



  !  chemical_const%colf_a0 = 5.83d-11/dsqrt(unit_temperature)
  !  chemical_const%colf_T0 = 157828.0_dp/unit_temperature
  !  chemical_const%alpha_a0 =(2.55d-13*(1.0d4*unit_temperature)**0.79_dp)
  end subroutine chemical_init

  subroutine chemical_set_default(self)
    implicit none
    type(chemical_parameters) :: self

    self%n_species            = 1
    self%small_to_zero        = .true.
    self%small_density        = 0.0_dp
    self%gas_mu               = 1.0_dp
    self%source_split         = .true.
    self%neutral_fraction     = 0.0_dp
    self%relative_fractionmin = 1.0d-6
    self%H_                   = 1
    self%H2_                  = 2
    self%He_                  = 3
  end subroutine chemical_set_default

  subroutine chemical_elements_init(thechemical_inuse)
    implicit none
    type(chemical_phys), intent(inout)    :: thechemical_inuse
    ! .. local ..
    integer                               :: I_ionmax,i_element,i_ion
    !-------------------------------------------------------------------
    Loop_i_chemical0 : do i_element=1,thechemical_inuse%myconfig%n_species
      if(i_element==thechemical_inuse%myconfig%H_)then
        thechemical_inuse%element(i_element)%ionmin=0
        thechemical_inuse%element(i_element)%ionmax=1
        thechemical_inuse%element(i_element)%elem_name='H'
      elseif(i_element==thechemical_inuse%myconfig%H2_)then
        thechemical_inuse%element(i_element)%ionmin=0
        thechemical_inuse%element(i_element)%ionmax=4
        thechemical_inuse%element(i_element)%elem_name='H2'
      elseif(i_element==thechemical_inuse%myconfig%He_)then
        thechemical_inuse%element(i_element)%ionmin=0
        thechemical_inuse%element(i_element)%ionmax=2
        thechemical_inuse%element(i_element)%elem_name='He'
      else
        write(*,*)' This element : ',i_element,' is not implimented'
        call mpistop('Stop at chemical_elements_init in mod_chemical.t')
      end if
      thechemical_inuse%nvar = thechemical_inuse%nvar+&
                               (thechemical_inuse%element(i_element)%ionmax-&
                                thechemical_inuse%element(i_element)%ionmin+1)
      I_ionmax = thechemical_inuse%element(i_element)%ionmax
      thechemical_inuse%element(i_element)%ionmax=merge(thechemical_inuse%element(i_element)%ionmax&
                                                      ,thechemical_inuse%element(i_element)%ionmax-1,i_element>1)
      allocate(thechemical_inuse%element(i_element)%ion(thechemical_inuse%element(i_element)%ionmin&
                :thechemical_inuse%element(i_element)%ionmax))

      Loop_ionisation0 :     do i_ion =thechemical_inuse%element(i_element)%ionmin,I_ionmax-1
       if(i_element==thechemical_inuse%myconfig%H_)then
        thechemical_inuse%element(i_element)%ion(i_ion)%myreccoef%alpha_a0= &
                      (2.55d-13*(1.0d4*unit_temperature)**0.79_dp)
        thechemical_inuse%element(i_element)%ion(i_ion)%myreccoef%colf_a0 = &
                      5.83d-11/dsqrt(unit_temperature)
        thechemical_inuse%element(i_element)%ion(i_ion)%myreccoef%colf_T0 = &
                      157828.0_dp/unit_temperature
       elseif(i_element==thechemical_inuse%myconfig%H2_)then
        thechemical_inuse%element(i_element)%ion(i_ion)%myreccoef%alpha_a0= 0.0_dp
        thechemical_inuse%element(i_element)%ion(i_ion)%myreccoef%colf_a0 = 0.0_dp
        thechemical_inuse%element(i_element)%ion(i_ion)%myreccoef%colf_T0 = 0.0_dp
       elseif(i_element==thechemical_inuse%myconfig%He_)then
        thechemical_inuse%element(i_element)%ion(i_ion)%myreccoef%alpha_a0= 0.0_dp
        thechemical_inuse%element(i_element)%ion(i_ion)%myreccoef%colf_a0 = 0.0_dp
        thechemical_inuse%element(i_element)%ion(i_ion)%myreccoef%colf_T0 = 0.0_dp
       else
        write(*,*)' This element : ',i_element,' is not implimented'
        call mpistop('Stop at chemical_elements_init in mod_chemical.t')
      end if
      end do Loop_ionisation0
    end do Loop_i_chemical0
  end subroutine chemical_elements_init

  subroutine chemical_set_physvariable(self_ind,phys_indices_inuse)
    use mod_global_parameters
    use mod_physics
    implicit none
    type(chemical_phys), intent(inout)           :: self_ind
    type(phys_variables_indices), intent(inout)  :: phys_indices_inuse
    ! .. local ..
    integer                                      :: i_elem,i_chemical,i_ion
    character(len=20)                            :: elem
  !--------------------------------------------
      i_elem = 1
      Loop_i_chemical1 : do i_chemical = 1, self_ind%myconfig%n_species
        write(elem, "(A,A)") trim(self_ind%element(i_chemical)%elem_name), "_"
        Loop_ionisation0 : do i_ion =self_ind%element(i_chemical)%ionmin,self_ind%element(i_chemical)%ionmax
         self_ind%element(i_chemical)%ion(i_ion)%myind_ = var_set_fluxvar(trim(elem), trim(elem), i_ion)
         phys_indices_inuse%chemical_element(i_elem) =  self_ind%element(i_chemical)%ion(i_ion)%myind_
         i_elem = i_elem+1
        end do  Loop_ionisation0
      end do Loop_i_chemical1
  end   subroutine chemical_set_physvariable

  !> Read this module s parameters from a file
  subroutine chemical_read_params(files,chemical_config)
    use mod_global_parameters, only: unitpar
    implicit none

    character(len=*), intent(in)             :: files(:)
    type(chemical_parameters), intent(inout) :: chemical_config
    ! .. local ..
    integer                      :: i_file, i_reason

    namelist /chemical_list/ chemical_config
    !----------------------------------------------------------

    Loop_files : do i_file = 1, size(files)
      open(unitpar, file=trim(files(i_file)), status="old")
      read(unitpar, chemical_list, iostat=i_reason)
      cond_ierror : if(i_reason>0)then
        write(*,*)' Error in reading the parameters file : ',trim(files(i_file))
        write(*,*)' Error at namelist: chemical_list'
        write(*,*)' The code stops now '
        call mpistop('At mod_chemical.t in the procedure : chemical_read_params')
      elseif(i_reason<0)then cond_ierror
        write(*,*)' Reache the end of the file  : ',trim(files(i_file))
        write(*,*)' Error at namelist: chemical_list'
        write(*,*)' The code stops now '
        call mpistop('At mod_chemical.t in the procedure : chemical_read_params')
      else cond_ierror
        write(*,*)' End of reading of the  chemical_list'
      end if  cond_ierror
      close(unitpar)
    end do Loop_files
  end subroutine chemical_read_params

  subroutine chemical_check_params()
  !  if (gas_mu <= 0.0d0) call mpistop ("Chemical error: gas_mu (molecular weight)"//&
  !       "negative or not set")

  end subroutine chemical_check_params

  subroutine chemical_to_conserved(ixI^L, ixO^L, w, x)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    real(kind=dp), intent(inout)    :: w(ixI^S, nw)
    real(kind=dp), intent(in)       :: x(ixI^S, 1:ndim)
    integer                         :: n, idir

    ! dummy
  end subroutine chemical_to_conserved

  subroutine chemical_to_primitive(ixI^L, ixO^L, w, x)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    real(kind=dp), intent(inout)    :: w(ixI^S, nw)
    real(kind=dp), intent(in)       :: x(ixI^S, 1:ndim)
    integer                         :: n, idir

    ! dummy
  end subroutine chemical_to_primitive

  subroutine chemical_get_flux(w, x, ixI^L, ixO^L, idim, f)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L, idim
    real(kind=dp), intent(in)       :: w(ixI^S, 1:nw), x(ixI^S, 1:^ND)
    real(kind=dp), intent(inout)    :: f(ixI^S, nwflux)
    ! .. local ..
    integer                         :: i_chemical,i_ion,ind_elem
    !----------------------------------------------
    Loop_chemical_element : do i_chemical = 1, thechemical%myconfig%n_species
      Loop_chemical_ion   : do i_ion = thechemical%element(i_chemical)%ionmin,&
                                       thechemical%element(i_chemical)%ionmax
       ind_elem = thechemical%element(i_chemical)%ion(i_ion)%myind_
       where (w(ixO^S,ind_elem) > thechemical%myconfig%small_density)
        f(ixO^S, ind_elem ) = w(ixO^S, gas_mom(idim))
       elsewhere
        f(ixO^S, ind_elem ) = 0.0d0
       end where


     end do  Loop_chemical_ion
    end do Loop_chemical_element

  end subroutine chemical_get_flux

  subroutine chemical_get_flux_prim(w, x, ixI^L, ixO^L, idim, f)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L, idim
    real(kind=dp), intent(in)       :: w(ixI^S, 1:nw), x(ixI^S, 1:^ND)
    real(kind=dp), intent(inout)    :: f(ixI^S, nwflux)
    integer                         :: i_chemical,i_ion,ind_elem
    !---------------------------------------------------------------
    Loop_chemical_element : do i_chemical = 1, thechemical%myconfig%n_species
     Loop_chemical_ion   : do i_ion = thechemical%element(i_chemical)%ionmin,&
                                     thechemical%element(i_chemical)%ionmax
      ind_elem = thechemical%element(i_chemical)%ion(i_ion)%myind_
      where (w(ixO^S, ind_elem) > thechemical%myconfig%small_density)
        f(ixO^S, ind_elem) = w(ixO^S, gas_mom(idim))*w(ixO^S, ind_elem)
      elsewhere             ! TODO: remove?
        f(ixO^S, ind_elem) = 0.0d0
      end where
     end do  Loop_chemical_ion
    end do Loop_chemical_element
  end subroutine chemical_get_flux_prim



  ! Force dust density to zero if dust_rho <= dust_min_rho
  subroutine set_chemicaltozero(qdt, ixI^L, ixO^L,  wCT,  w, x)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    real(kind=dp), intent(in)       :: qdt
    real(kind=dp), intent(in)       :: wCT(ixI^S, 1:nw), x(ixI^S, 1:ndim)
    real(kind=dp), intent(inout)    :: w(ixI^S, 1:nw)
    integer                         :: n, idir
    integer                         :: flag(ixI^S)
    flag(ixO^S) = 0
    call chemical_check_w(.false., ixI^L, ixO^L, flag, w)
    call chemical_set_floor(.false.,ixI^L,ixO^L,flag,x,w)
  end subroutine set_chemicaltozero



  !> w[iw]= w[iw]+qdt*S[wCT,  x] where S is the source based on wCT within ixO
  subroutine chemical_add_source(qdt, ixI^L, ixO^L, wCT,w, x, qsourcesplit, active)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    real(kind=dp), intent(in)       :: qdt
    real(kind=dp), intent(in)       :: wCT(ixI^S, 1:nw), x(ixI^S, 1:ndim)
    real(kind=dp), intent(inout)    :: w(ixI^S, 1:nw)
    logical, intent(in)             :: qsourcesplit
    logical, intent(inout)          :: active

    ! .. local ..
    logical          :: fd_flag(ixI^S)
    real(kind=dp)    :: ptherm(ixI^S), vgas(ixI^S, ndir),vdust(ixI^S),dv(ixI^S)
    real(kind=dp)    :: fdrag(ixI^S, ndir, thechemical%nvar)
    integer          :: i_chemical,idir, it_diff

    select case(TRIM(thechemical%myconfig%method) )
    case('molecular')
      call chemical_atomic(ixI^L,ixO^L,qdt,x,wCT,w)
    case( 'none' )
      !do nothing here
    case default !all regular chemical methods here
        call chemical_atomic(ixI^L,ixO^L,qdt,x,wCT,w)
        cond_floor: if (thechemical%myconfig%small_to_zero) then
          call set_chemicaltozero(qdt, ixI^L, ixO^L,  wCT,  w, x)
        end if cond_floor
    end select

  end subroutine chemical_add_source

  !> Get dt related to dust and gas stopping time (Laibe 2011)
  subroutine chemical_get_dt(w, ixI^L, ixO^L, dtnew, dx^D, x)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    real(kind=dp), intent(in)       :: dx^D, x(ixI^S, 1:^ND)
    real(kind=dp), intent(in)       :: w(ixI^S, 1:nw)
    real(kind=dp), intent(inout)    :: dtnew

    real(kind=dp)                                      :: vgas(ixI^S, ndir)
    real(kind=dp), dimension(thechemical%nvar)         :: dtdust
    real(kind=dp), dimension(ixI^S)                    :: vt2, deltav, tstop, vdust,ptherm
    real(kind=dp), dimension(ixI^S,1:thechemical%nvar) :: alpha_T
    real(kind=dp)                                      :: K
    integer                                            :: n, idir
    ! dummy for the moment

  end subroutine chemical_get_dt


  !--------------------------------------------------------------------
  !> subroutine check error in dust
  subroutine chemical_check_w(primitive, ixI^L, ixO^L, flag, w)
    use mod_global_parameters
    implicit none
    logical, intent(in)       :: primitive
    integer, intent(in)       :: ixI^L, ixO^L
    integer, intent(inout)    :: flag(ixI^S)
    real(dp), intent(in)      :: w(ixI^S,1:nw)
    ! .. local ..
    integer                   :: i_chemical,i_ion,ind_elem
    !------------------------------------
    Loop_chemical :  do i_chemical=1,thechemical%myconfig%n_species
      Loop_chemical_ion   : do i_ion = thechemical%element(i_chemical)%ionmin,&
                                      thechemical%element(i_chemical)%ionmax
       ind_elem = thechemical%element(i_chemical)%ion(i_ion)%myind_
       where(w(ixO^S,ind_elem)<=thechemical%myconfig%small_density)
        flag(ixO^S) = - ind_elem
       end where
     end do Loop_chemical_ion
    end do Loop_chemical
  end subroutine chemical_check_w
  !--------------------------------------------------------------------
  !> subroutine to handel small values
  subroutine chemical_set_floor(primitive,ixI^L,ixO^L,flag,x,w)
    use mod_global_parameters
    implicit none
    logical, intent(in)       :: primitive
    integer, intent(in)       :: ixI^L, ixO^L
    integer, intent(in)       :: flag(ixI^S)
    real(dp), intent(in)      :: x(ixI^S,1:ndim)
    real(dp), intent(inout)   :: w(ixI^S,1:nw)
    ! .. local ..
    integer                   :: i_chemical,idir,i_ion,ind_elem
    logical, dimension(ixI^S) :: patch_correct
    !------------------------------------
    Loop_chemical :  do i_chemical=1,thechemical%myconfig%n_species
      Loop_chemical_ion   : do i_ion = thechemical%element(i_chemical)%ionmin,&
                                            thechemical%element(i_chemical)%ionmax
       ind_elem = thechemical%element(i_chemical)%ion(i_ion)%myind_
       patch_correct(ixO^S) = (flag(ixO^S)==-ind_elem).or.(flag(ixO^S)>0)
       where(patch_correct(ixO^S))
         w(ixO^S,ind_elem) = 0.0_dp
       end where
      end do  Loop_chemical_ion
    end do Loop_chemical
  end subroutine chemical_set_floor
  !------------------------------------------------------------------------------------------
  subroutine chemical_cooling_set_alpha(ixI^L,ixO^L,i_element,i_ion,Temperature,alpha)
    use mod_global_parameters
    implicit none
    integer, intent(in)       :: ixI^L, ixO^L
    integer, intent(in)       :: i_element,i_ion
    real(dp), intent(in)      :: Temperature(ixI^S)
    real(dp), intent(inout)   :: alpha(ixO^S)
    !........................................


    alpha(ixO^S)=thechemical%element(i_element)%ion(i_ion)%myreccoef%alpha_a0/Temperature(ixO^S)**0.79_dp
  end   subroutine chemical_cooling_set_alpha

!------------------------------------------------------------------------------------------
  subroutine chemical_cooling_set_colf(ixI^L,ixO^L,i_element,i_ion,Temperature,colf)
    use mod_global_parameters
    implicit none
    integer, intent(in)       :: ixI^L, ixO^L
    integer, intent(in)       :: i_element,i_ion
    real(dp), intent(in)      :: Temperature(ixI^S)
    real(dp), intent(inout)   :: colf(ixO^S)
    !........................................

    colf(ixO^S)=thechemical%element(i_element)%ion(i_ion)%myreccoef%colf_a0&
                *dexp(-thechemical%element(i_element)%ion(i_ion)%myreccoef%colf_T0/Temperature(ixO^S))
  end   subroutine chemical_cooling_set_colf
!------------------------------------------------------------------------------------------
  subroutine chemical_coronal_ionisationequilibrium_H(ixI^L,ixO^L,qdt,&
                  temperature,density_H,&
                  fractionHI_CT,fractionHI_new)
    use mod_global_parameters
    implicit none
    integer, intent(in)       :: ixI^L, ixO^L
    real(dp), intent(in)      :: qdt
    real(dp), intent(in)      :: Temperature(ixI^S)
    real(dp), intent(in)      :: density_H(ixO^S)
    real(dp), intent(in)      :: fractionHI_CT(ixO^S)
    real(dp), intent(inout)   :: fractionHI_new(ixO^S)
    ! .. local ..
    real(dp), dimension(ixO^S):: colf,alpha
    real(dp), dimension(ixO^S):: equ2_A,equ2_B,equ2_C,equ2_Detr
    real(dp), dimension(ixO^S):: ode_e,ode_g0
    !--------------------------------------------------

    call chemical_cooling_set_alpha(ixI^L,ixO^L,thechemical%myconfig%H_,0,Temperature,alpha)
    call chemical_cooling_set_colf(ixI^L,ixO^L,thechemical%myconfig%H_,0,Temperature,colf)
    ! solve second order ODE
    equ2_A = alpha+colf
    equ2_B = -((2.0_dp+thechemical%myconfig%neutral_fraction)*alpha+&
                      (1.0_dp+thechemical%myconfig%neutral_fraction)*colf)
    equ2_C = (1.0_dp+thechemical%myconfig%neutral_fraction)*alpha
    equ2_Detr = dsqrt(equ2_B**2.0_dp-4.0_dp*equ2_A*equ2_C)

    ode_g0=(2.*equ2_A*fractionHI_CT+equ2_B+equ2_Detr)&
     /(2.*equ2_A**fractionHI_CT+equ2_B-equ2_Detr)
    ode_e=dexp( -equ2_Detr*density_H*qdt )
  !
    fractionHI_new=(-equ2_B-equ2_Detr*(1.0_dp+ode_g0*ode_e) / &
                   (1.0_dp-ode_g0*ode_e))/(2.0_dp*equ2_A) !# the new neutral fraction
    fractionHI_new=max(min(fractionHI_new,1.0_dp),thechemical%myconfig%relative_fractionmin)

  end subroutine chemical_coronal_ionisationequilibrium_H


  subroutine chemical_atomic(ixI^L,ixO^L,qdt,x,wCT,w)
    use mod_global_parameters
    implicit none
    integer, intent(in)            :: ixI^L, ixO^L
    real(kind=dp), intent(in)      :: x(ixI^S,1:ndim)
    real(kind=dp), intent(in)      :: wCT(ixI^S,1:nw)
    real(kind=dp), intent(inout)   :: w(ixI^S,1:nw)
    real(dp), intent(in)           :: qdt
    ! .. local ..
    real(dp),dimension(ixO^S)      :: density_H
    real(dp),dimension(ixI^S)      :: temperature
    real(dp),dimension(ixO^S)      :: fractionHI_CT,fractionHI_new
    !----------------------------------------------------

    call phys_get_temperature(ixI^L,ixO^L,wCT,x,temperature)
    call chemical_get_density_H(ixI^L,ixO^L,wCT,density_H)
    call chemical_get_fractionHI(ixI^L,ixO^L,wCT,density_H,fractionHI_CT)


    call chemical_coronal_ionisationequilibrium_H(ixI^L,ixO^L,qdt,&
                    temperature,density_H,&
                    fractionHI_CT,fractionHI_new)
    call chemical_get_HI_density(ixI^L, ixO^L,fractionHI_new,w)
  end subroutine chemical_atomic

  subroutine  chemical_get_density_H(ixI^L,ixO^L,w,density_H)
    use mod_global_parameters
    implicit none
    integer, intent(in)         :: ixI^L, ixO^L
    real(kind=dp), intent(in)   :: w(ixI^S,1:nw)
    real(kind=dp),intent(inout) :: density_H(ixO^S)
    !....................................
    density_H = w(ixO^S,gas_rho_)
  end subroutine  chemical_get_density_H

  subroutine chemical_get_HI_density(ixI^L, ixO^L,fractionHI,w)
    use mod_global_parameters
    implicit none
    integer, intent(in)            :: ixI^L, ixO^L
    real(kind=dp), intent(in)      :: fractionHI(ixO^S)
    real(kind=dp), intent(inout)   :: w(ixI^S,1:nw)
    real(kind=dp)                  :: density_H(ixO^S)
    ! .. local ..
    integer                        :: ind_elem,H_
    !-----------------------------------------

    call   chemical_get_density_H(ixI^L,ixO^L,w,density_H)
    H_ = thechemical%myconfig%H_
    ind_elem = thechemical%element(H_)%ion(0)%myind_
    w(ixO^S,ind_elem) = fractionHI(ixO^S)*density_H(ixO^S)
  end subroutine chemical_get_HI_density

  subroutine chemical_get_fractionHI(ixI^L,ixO^L,w,density_H,fractionHI)
    use mod_global_parameters
    implicit none
    integer, intent(in)            :: ixI^L, ixO^L
    real(kind=dp), intent(in)      :: w(ixI^S,1:nw)
    real(dp), intent(in)           :: density_H(ixO^S)
    real(dp), intent(out)          :: fractionHI(ixO^S)
    real(dp)                       :: neutral_density(ixO^S)
    !---------------------------------------------------------
    call chemical_get_neutral_density(ixI^L, ixO^L,thechemical%myconfig%H_,w,neutral_density)
    fractionHI = neutral_density(ixO^S)/density_H(ixO^S)

    where(fractionHI<0.0_dp)fractionHI=0.0_dp
  end subroutine chemical_get_fractionHI

  subroutine chemical_get_neutral_density(ixI^L, ixO^L,i_chemical_inuse, w,neutral_density)
    use mod_global_parameters
    implicit none
    integer, intent(in)            :: ixI^L, ixO^L
    integer, intent(in)            :: i_chemical_inuse
    real(kind=dp), intent(in)      :: w(ixI^S,1:nw)
    real(dp), intent(out)          :: neutral_density(ixO^S)
    ! .. local ..
    integer                        :: i_chemical,i_ion,i_min_ion,i_max_ion
    integer                        :: ind_elem
    !--------------------------------------------------------------


    if(thechemical%element(i_chemical_inuse)%ionmin==0) then
        ind_elem = thechemical%element(i_chemical_inuse)%ion(0)%myind_
        neutral_density(ixO^S) = w(ixO^S,ind_elem)
    else


      neutral_density(ixO^S) = 0.0_dp
      Loop_chemical :  do i_chemical=1,thechemical%myconfig%n_species
        if(i_chemical_inuse==i_chemical)then
          i_min_ion =  max(thechemical%element(i_chemical_inuse)%ionmin,1)
          i_max_ion = thechemical%element(i_chemical)%ionmax
        else
         i_min_ion = 0
         i_max_ion = thechemical%element(i_chemical)%ionmax
        end if
        Loop_chemical_ion   : do i_ion = i_min_ion,i_max_ion
         ind_elem = thechemical%element(i_chemical_inuse)%ion(i_ion)%myind_
         neutral_density(ixO^S) = neutral_density(ixO^S)+w(ixO^S,ind_elem)
        end do Loop_chemical_ion
      end do Loop_chemical
      neutral_density(ixO^S) = w(ixO^S,phys_ind%rho_) - neutral_density(ixO^S)
     end if


   end subroutine chemical_get_neutral_density

  subroutine chemical_get_electron_density(ixI^L, ixO^L,w,rho_electron)

    implicit none
    integer, intent(in)            :: ixI^L, ixO^L
    real(kind=dp), intent(in)      :: w(ixI^S,1:nw)
    real(dp), intent(out)          :: rho_electron(ixO^S)
    real(dp)                       :: neutral_density(ixO^S)
    ! .. local ..
    integer                        :: i_chemical,i_ion,ind_elem
    !-------------------------------------

    if(thechemical%myconfig%neutral_fraction<smalldouble) then
      rho_electron = 0.0_dp
    else
      call chemical_get_neutral_density(ixI^L, ixO^L,thechemical%myconfig%H_,w,neutral_density)
      rho_electron = thechemical%myconfig%neutral_fraction*neutral_density(ixO^S)
    end if
    Loop_chemical :  do i_chemical=1,thechemical%myconfig%n_species
          Loop_chemical_ion   : do i_ion = max(thechemical%element(i_chemical)%ionmin,1),&
                                           thechemical%element(i_chemical)%ionmax
            ind_elem = thechemical%element(i_chemical)%ion(i_ion)%myind_
            rho_electron = i_ion*w(ixO^S,ind_elem)+&
                   rho_electron
          end do Loop_chemical_ion
    end do  Loop_chemical
    end subroutine chemical_get_electron_density
end module mod_chemical
