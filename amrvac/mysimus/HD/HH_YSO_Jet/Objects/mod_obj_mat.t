module mod_obj_mat
use mod_constants
use mod_global_parameters
use mod_obj_global_parameters
implicit none

type usrboundary_parameters
    character(len=78)             :: obj_name               !> Obj name that call it
    integer                       :: myindice       !> boundary indices associated with ism in use
    character(len=30)             :: boundary_type(3,2)
    logical                       :: useprimitive
    integer                       :: nghostcells(3)
    real(kind=dp)                 :: flux_frac
    logical                       :: special_origin_theta
    logical                       :: mixed_fixed_bound(3,2,100) !> ism 'fix' flux variables :
    !TO DO    : f95 does not allow to namelist an allocatable array, but fortran 2003+ does
    !So, for now, take a big enough useful array and the same size as in mod_obj_ism.t
    integer                       :: LHS
    integer                       :: RHS
    integer                       :: ngc
end type usrboundary_parameters
type usrboundary_type
  character(len=30), allocatable    :: variable_typebound(:,:,:)
  logical,allocatable               :: mixed_fixed_bound(:,:,:) !> ism 'fix' flux variables
  type(usrboundary_parameters)      :: myconfig
  !character(len=30)                 :: boundary_type(3,2)
  !logical                           :: useprimitive
  !integer                           :: nghostcells(3)
  !real(kind=dp)                     :: flux_frac
   contains
   !PRIVATE
   PROCEDURE, PASS(self) :: read_parameters    => usr_boundaries_read_p
   PROCEDURE, PASS(self) :: set_default        => usr_boundaries_set_default
   PROCEDURE, PASS(self) :: set_complet        => usr_boundaries_set_complet
   PROCEDURE, PASS(self) :: set_w              => usr_boundaries_set_w
end type usrboundary_type
contains

!--------------------------------------------------------------------

subroutine usr_boundaries_set_default(self)
  implicit none
  class(usrboundary_type)              :: self
   !------------------------------------
   self%myconfig%useprimitive       = .false.
   self%myconfig%boundary_type      = 'open'
   self%myconfig%nghostcells        = 2
   self%myconfig%flux_frac          = 1.0d-2
   self%myconfig%obj_name           = 'boundary'
   self%myconfig%myindice           = 1
   self%myconfig%LHS           = 0
   self%myconfig%RHS           = 0
   self%myconfig%ngc           = 0
   self%myconfig%special_origin_theta = .false.
   if(.not.allocated(self%mixed_fixed_bound))&
    allocate(self%mixed_fixed_bound(1:ndim,1:2,1:nwfluxbc))
   self%mixed_fixed_bound(1:ndim,1:2,1:nwfluxbc) = .false.
   self%myconfig%mixed_fixed_bound(1:ndim,1:2,:) = .false.
end subroutine usr_boundaries_set_default


!--------------------------------------------------------------------
   !> Read the ism parameters  from a parfile
    subroutine usr_boundaries_read_p(self,usrboundary_config,files)

      implicit none
      class(usrboundary_type)                    :: self
      character(len=*),intent(in)                :: files(:)
      type(usrboundary_parameters), intent(out)  :: usrboundary_config
      ! .. local ..
      integer                            :: i_file,i_error_read
      namelist /usr_usrboundary_list/  usrboundary_config
      namelist /usr_usrboundary_ism_list/  usrboundary_config
      namelist /usr_usrboundary_cloud_list/ usrboundary_config
      namelist /usr_usrboundary_clajet_list/ usrboundary_config
      namelist /usr_usrboundary_reljet_list/ usrboundary_config
      namelist /usr_usrboundary_relwind_list/ usrboundary_config

      if(mype==0)write(*,*)'Reading usr_boundary_list',usrboundary_config%myindice,self%myconfig%obj_name
      do i_file = 1, size(files)
         open(unitpar, file=trim(files(i_file)), status="old")
         select case(usrboundary_config%myindice)
         case(1)
           read(unitpar, usr_usrboundary_ism_list, iostat=i_error_read)
         case(2)
           read(unitpar, usr_usrboundary_cloud_list, iostat=i_error_read)
         case(3)
           read(unitpar, usr_usrboundary_clajet_list, iostat=i_error_read)
         case(4)
           read(unitpar, usr_usrboundary_reljet_list, iostat=i_error_read)
         case(5)
             read(unitpar, usr_usrboundary_relwind_list, iostat=i_error_read)
         case default
           read(unitpar, usr_usrboundary_list, iostat=i_error_read)
         end select
         call usr_mat_read_error_message(i_error_read,usrboundary_config%myindice,&
                                         self%myconfig%obj_name)
         close(unitpar)
      end do


    end subroutine usr_boundaries_read_p
   !------------------------------------------------------------------------
!> Subtroutine initialise the boundary conditions
subroutine usr_boundaries_set_complet(self)
  use mod_physics
  implicit none
  class(usrboundary_type)               :: self
  ! .. local ..
  integer                               :: idims,iside,idir,iw
   !------------------------------------
  if(allocated(self%variable_typebound))deallocate(self%variable_typebound)
  allocate(self%variable_typebound(ndim,2,nwfluxbc))


  self%myconfig%nghostcells=nghostcells
  Loop_idims : do idims=1,ndim
    Loop_iside : do iside=1,2
      !initialize mixed_fixed_bound
      do iw=1,nwfluxbc
        !self%myconfig%mixed_fixed_bound may be fixed and transmitted by ISMs
        !so pass it to self%mixed_fixed_bound
        self%myconfig%mixed_fixed_bound(idims,iside,iw)=self%mixed_fixed_bound(idims,iside,iw)
      end do

     select case(self%myconfig%boundary_type(idims,iside))
     case('wall')
      self%variable_typebound(idims,iside,:)='cont'
      Loop_iw_wall: do idir = 1,ndir
        if(idir==idims) then
          self%variable_typebound(idims,iside,phys_ind%mom(idir))='asymm'
        else
          self%variable_typebound(idims,iside,phys_ind%mom(idir))='symm'
        end if
      end do Loop_iw_wall
    case('noinflow')
      self%variable_typebound(idims,iside,:)='cont'
      self%variable_typebound(idims,iside,phys_ind%mom(idims))='noinflow'
     case('limitinflow')
      self%variable_typebound(idims,iside,:)='cont'
      self%variable_typebound(idims,iside,phys_ind%mom(idims))='limitinflow'
     case('limitinflowdisc')
      self%variable_typebound(idims,iside,:)='symm'
      self%variable_typebound(idims,iside,phys_ind%mom(r_))='symm'
      self%variable_typebound(idims,iside,phys_ind%mom(phi_))='symm'
      self%variable_typebound(idims,iside,phys_ind%mom(z_))='limitinflow'
     case('nooutflow')
      self%variable_typebound(idims,iside,:)='cont'
      self%variable_typebound(idims,iside,phys_ind%mom(idims))='nooutflow'
     case('limitoutflow')
      self%variable_typebound(idims,iside,:)='cont'
      self%variable_typebound(idims,iside,phys_ind%mom(idims))='limitoutflow'
     case('limitoutflowdisc')
      self%variable_typebound(idims,iside,:)='symm'
      self%variable_typebound(idims,iside,phys_ind%mom(r_))='symm'
      self%variable_typebound(idims,iside,phys_ind%mom(phi_))='symm'
      self%variable_typebound(idims,iside,phys_ind%mom(z_))='limitoutflow'
    case('grad')
       self%variable_typebound(idims,iside,:)='cont'
       self%variable_typebound(idims,iside,phys_ind%mom(idims))='grad'
     case('open')
       self%variable_typebound(idims,iside,:)='cont'
     case('disc')
       self%variable_typebound(idims,iside,:)='symm'
       self%variable_typebound(idims,iside,2)='symm'
       self%variable_typebound(idims,iside,3)='asymm'
       self%variable_typebound(idims,iside,4)='symm'
     case('axis')
      self%variable_typebound(idims,iside,:)='symm'
      self%variable_typebound(idims,iside,phys_ind%mom(r_))='asymm'
      self%variable_typebound(idims,iside,phys_ind%mom(phi_))='asymm'
      self%variable_typebound(idims,iside,phys_ind%mom(z_))='symm'

     end select

     !set the variables to fix
     do iw=1,nwfluxbc
      if(self%mixed_fixed_bound(idims,iside,iw))then
        self%variable_typebound(idims,iside,iw)='fix'
      end if
     end do

    end do Loop_iside
  end do Loop_idims
end subroutine usr_boundaries_set_complet
!--------------------------------------------------------------------
!> Subtroutine set the boundary conditions
subroutine usr_boundaries_set_w(ixI^L,ixO^L,iB,idims,iside,&
                              patchw,x,w,self)

  use mod_physics
  implicit none

  integer, intent(in)                    :: ixI^L,ixO^L,iB,idims,iside
  logical,intent(in)                     :: patchw(ixI^S)
  real(kind=dp), intent(in)              :: x(ixI^S,1:ndim)
  real(kind=dp), intent(inout)           :: w(ixI^S,1:nw)
  class(usrboundary_type)                :: self
  ! ... local ...
  integer                                :: ixG^L,ix^D,iw,isymm,ixBC^L
  real(kind=dp)                          :: wp(ixI^S,1:nw)


  !-----------------------------------------------------------------------------
  ! * According to subroutine usr_ism_set_w in which it is called twice in mod_obj_ism.t :
  ! the input integers are ixI^L=ixG^L=ixGmin1,ixGmin2,ixGmax1,ixGmax2 (or ixCoGmin1,ixCoGmin2,ixCoGmax1,ixCoGmax2
  ! when working with a coarsened grid)
  ! = the range delimiting the whole domain (including boundary ghost cells), i.e. between integers
  ! 1 and the highest possible indices for the coordinates for the grid for each dimension
  ! and ixO^L=ixI^L defined from ixB^L as :
  ! > idims = 1, iside==2 (maximal boundary)
  ! ==> ixImin1=ixBmax1+1-nghostcells;ixImin2=ixBmin2;ixImax1=ixBmax1;ixImax2=ixBmax2;
  ! > idims = 1, iside==1 (minimal boundary)
  ! ==> ixImin1=ixBmin1;ixImin2=ixBmin2;ixImax1=ixBmin1-1+nghostcells;ixImax2=ixBmax2;
  ! > idims = 2, iside==2 (maximal boundary)
  ! ==> ixImin1=ixBmin1;ixImin2=ixBmax2+1-nghostcells;ixImax1=ixBmax1;ixImax2=ixBmax2;
  ! > idims = 2, iside==1 (minimal boundary)
  ! ==> ixImin1=ixBmin1;ixImin2=ixBmin2;ixImax1=ixBmax1;ixImax2=ixBmin2-1+nghostcells;
  ! ==> no need to redefine ixO^L (or ixI^L in bc_phys) in this subroutine usr_boundaries_set_w,
  ! manipulate w directly !
  ! to sum up, the range delimiting each boundary individually
  ! idims = the direction along which the boundary conditions must be applied
  ! iside = side (min or max) along which the bc must be applied


  ! 04-03-22 : must be corrected to be identical to bc_phys in boundary_conditions.t
  ! * make sure that the correct operations are done on the correct part(s) of the whole grid




  select case (idims)
  {case (^D)
   if (iside==2) then
      ! maximal boundary
      ! 05-03-22 : the way this subroutine is called from usr_special_bc      => specialbound_usr
      ! in mod_obj_usr_yso_jet.t, there is no need to redefine the input/output integers:

      if(self%myconfig%useprimitive)then
       call phys_to_primitive(ixI^L,ixO^L,w,x)
      end if
      ! cont/symm/asymm types
      do iw=1,nwfluxbc
         select case (self%variable_typebound(idims,iside,iw))
         case ("symm")
            w(ixO^S,iw) = w(ixOmin^D-1:ixOmin^D-self%myconfig%nghostcells(^D):-1^D%ixO^S,iw)
         case ("asymm")
            w(ixO^S,iw) = -w(ixOmin^D-1:ixOmin^D-self%myconfig%nghostcells(^D):-1^D%ixO^S,iw)
         case ("cont")
           do ix^D=ixOmin^D,ixOmax^D
              w(ix^D^D%ixO^S,iw) = w(ixOmin^D-1^D%ixO^S,iw)
           end do
         case("noinflow")
           if (iw==1+^D)then
             do ix^D=ixOmin^D,ixOmax^D
                 w(ix^D^D%ixO^S,iw) = max(w(ixOmin^D-1^D%ixO^S,iw),zero)
             end do
           else
             do ix^D=ixOmin^D,ixOmax^D
                 w(ix^D^D%ixO^S,iw) = w(ixOmin^D-1^D%ixO^S,iw)
             end do
           end if
         case("nooutflow")
           if (iw==1+^D)then
             do ix^D=ixOmin^D,ixOmax^D
                 w(ix^D^D%ixO^S,iw) = min(w(ixOmin^D-1^D%ixO^S,iw),zero)
             end do
           else
             do ix^D=ixOmin^D,ixOmax^D
                 w(ix^D^D%ixO^S,iw) = w(ixOmin^D-1^D%ixO^S,iw)
             end do
           end if
         case("limitinflow")
           if (iw==1+^D)then
             do ix^D=ixOmin^D,ixOmax^D
                 w(ix^D^D%ixO^S,iw) = max(w(ixOmin^D-1^D%ixO^S,iw),&
                                         self%myconfig%flux_frac*w(ixOmin^D-1^D%ixO^S,iw))
             end do
           else
             do ix^D=ixOmin^D,ixOmax^D
                 w(ix^D^D%ixO^S,iw) = w(ixOmin^D-1^D%ixO^S,iw)
             end do
           end if
         case("limitoutflow")
           if (iw==1+^D)then
             do ix^D=ixOmin^D,ixOmax^D
                 w(ix^D^D%ixO^S,iw) = min(w(ixOmin^D-1^D%ixO^S,iw),&
                                         self%myconfig%flux_frac*w(ixOmin^D-1^D%ixO^S,iw))
             end do
           else
             do ix^D=ixOmin^D,ixOmax^D
                 w(ix^D^D%ixO^S,iw) = w(ixOmin^D-1^D%ixO^S,iw)
             end do
           end if
          case('grad')
              do ix^D=ixOmin^D,ixOmax^D
                  w(ix^D^D%ixO^S,iw) = 2.0_dp*w(ixOmin^D-1^D%ixO^S,iw) &
                                        -w(ixOmin^D-2^D%ixO^S,iw)
              end do
          case('fix')
              !skip it here, do after the whole domain has been set before
              !the call of this subroutine in mod_obj_usr_yso_jet.t:specialboundary_usr
         case default
            write (unitterm,*) "Undefined boundarytype ",&
               self%variable_typebound(idims,iside,iw), &
               "for variable iw=",iw," and side iB=",iB
         end select
      end do
   else
      ! minimal boundary
      ! 05-03-22 : the way this subroutine is called from usr_special_bc      => specialbound_usr
      ! in mod_obj_usr_yso_jet.t, there is no need to redefine the input/output integers:
      ! Legend:
      ! * ixO^L = ixI^L of bc_phys = defined in boundary_conditions.t:getbc:bc_phys:usr_special_bc as the limits
      ! of only the boundaries (i.e. the ghost cells only)
      ! * ixI^L = general ixG^L = the whole domain (including the boundaries)
      ! if(iB/=ismin^D)call mpistop("iB is broken with ismin^D !!!")


      if(self%myconfig%useprimitive)then
       call phys_to_primitive(ixI^L,ixO^L,w,x)
      end if
      ! cont/symm/asymm types
      !print*, 'nwfluxbc = ',nwfluxbc
      do iw=1,nwfluxbc
         select case (self%variable_typebound(idims,iside,iw))
         case ("symm")
               w(ixO^S,iw) = w(ixOmax^D+self%myconfig%nghostcells(^D):ixOmax^D+1:-1^D%ixO^S,iw)

               !if(idims==2.and.iside==1)then
               !print*,'ixImin^D=',ixImin^DD
               !print*,'ixImax^D=',ixImax^DD
               !print*,'ixOmin^D=',ixOmin^DD
               !print*,'ixOmax^D=',ixOmax^DD
               !print*,'ixGmin^D=',ixGmin^DD
               !print*,'ixGmax^D=',ixGmax^DD
               !print*,'ixBCmin^D=',ixBCmin^DD
               !print*,'ixBCmax^D=',ixBCmax^DD
               !end if
         case ("asymm")
             w(ixO^S,iw) = -w(ixOmax^D+self%myconfig%nghostcells(^D):ixOmax^D+1:-1^D%ixO^S,iw)
             !print*,cons_wnames(iw),' = ',w(ixO^S,iw)
             !print*,'min(',cons_wnames(iw),') = ',MINVAL(w(ixO^S,iw))
             !print*,'max(',cons_wnames(iw),') = ',MAXVAL(w(ixO^S,iw))
         case ("cont")
             do ix^D=ixOmin^D,ixOmax^D
                w(ix^D^D%ixO^S,iw) = w(ixOmax^D+1^D%ixO^S,iw)
             end do
            !print*,cons_wnames(iw),' = ',w(ixO^S,iw)
            !print*,'min(',cons_wnames(iw),') = ',MINVAL(w(ixO^S,iw))
            !print*,'max(',cons_wnames(iw),') = ',MAXVAL(w(ixO^S,iw))
         case("noinflow")
           if (iw==1+^D)then
             do ix^D=ixOmin^D,ixOmax^D
                 w(ix^D^D%ixO^S,iw) = min(w(ixOmax^D+1^D%ixO^S,iw),zero)
             end do
           else
             do ix^D=ixOmin^D,ixOmax^D
                 w(ix^D^D%ixO^S,iw) = w(ixOmax^D+1^D%ixO^S,iw)
             end do
           end if
          case("nooutflow")
            if (iw==1+^D)then
              do ix^D=ixOmin^D,ixOmax^D
                  w(ix^D^D%ixO^S,iw) = max(w(ixOmax^D+1^D%ixO^S,iw),zero)
              end do
            else
              do ix^D=ixOmin^D,ixOmax^D
                  w(ix^D^D%ixO^S,iw) = w(ixOmax^D+1^D%ixO^S,iw)
              end do
            end if
         case("limitinflow")
             if (iw==1+^D)then
               do ix^D=ixOmin^D,ixOmax^D
                   w(ix^D^D%ixO^S,iw) = min(w(ixOmax^D+1^D%ixO^S,iw),&
                                           self%myconfig%flux_frac*w(ixOmax^D+1^D%ixO^S,iw))
               end do
             else
               do ix^D=ixOmin^D,ixOmax^D
                   w(ix^D^D%ixO^S,iw) = w(ixOmax^D+1^D%ixO^S,iw)
               end do
             end if
            !print*,cons_wnames(iw),' = ',w(ixO^S,iw)
            !print*,'min(',cons_wnames(iw),') = ',MINVAL(w(ixO^S,iw))
            !print*,'max(',cons_wnames(iw),') = ',MAXVAL(w(ixO^S,iw))
         case("limitoutflow")
           if (iw==1+^D)then
             do ix^D=ixOmin^D,ixOmax^D
                 w(ix^D^D%ixO^S,iw) = max(w(ixOmax^D+1^D%ixO^S,iw),&
                                         self%myconfig%flux_frac*w(ixOmax^D+1^D%ixO^S,iw))
             end do
           else
             do ix^D=ixOmin^D,ixOmax^D
                 w(ix^D^D%ixO^S,iw) = w(ixOmax^D+1^D%ixO^S,iw)
             end do
           end if
          case('grad')
             do ix^D=ixOmin^D,ixOmax^D
                  w(ix^D^D%ixO^S,iw) = 2.0_dp*w(ixOmax^D+1^D%ixO^S,iw) &
                                        -w(ixOmax^D+2^D%ixO^S,iw)
             end do
          case('fix')
                 !skip it here, do after the whole domain has been set before
                 !the call of this subroutine in mod_obj_usr_yso_jet.t:specialboundary_usr
         case default
            write (unitterm,*) "Undefined boundarytype ",&
               self%variable_typebound(idims,iside,iw), &
               "for variable iw=",iw," and side iB=",iB
         end select
      end do
   end if \}
  end select

  if(self%myconfig%useprimitive)then
   call phys_to_conserved(ixI^L,ixO^L,w,x)
  end if
end  subroutine usr_boundaries_set_w
!--------------------------------------------------------------------
!> Subtroutine for power law distribution,here the dust sizes are defined. Ndust bins, with all bins having equal total mass.

subroutine usr_mat_powerlaw_withX(n_point,power_a,min_var,max_var,   &
                                  var_r)
 implicit none
 integer, intent(in)    :: n_point
 real(dp), intent(in)   :: power_a
 real(dp), intent(in)   :: min_var,max_var
 real(dp), intent(inout):: var_r(n_point)
 ! .. local ..
 integer                :: i
 real(dp)               :: r(0:n_point)
 !------------------------------------
 if(dabs(power_a-1.0_dp)<smalldouble)then
  call mpistop('power_a ==1 at usr_mat_powerlaw_withX in mod_usr.t')
 end if
 r(0) = min_var

 Loop_point : do i = 1,n_point
  r(i)      = (dsqrt(r(i-1))+(dsqrt(max_var)-dsqrt(min_var))/n_point)**2.0_dp
  !dvar_r(i) = r(i)-r(i-1)
  var_r(i)  = -1.0_dp/(power_a-1.0_dp)*&
              (r(i)**(power_a+2.0_dp)-r(i-1)**(power_a+2.0_dp))/&
              (r(i)**(power_a+1)-r(i-1)**(power_a+1))
 end do Loop_point
end subroutine usr_mat_powerlaw_withX

!--------------------------------------------------------------------
!> subroutien to set profile de distance r
subroutine usr_mat_profile(ixI^L,ixO^L,typeaxial_loc,profile, &
                          center,extend,x,fprofile)
 implicit none
 integer,intent(in)           :: ixI^L,ixO^L
 character(len=*), intent(in) :: typeaxial_loc !< Name of the coordinate system
 character(len=*), intent(in) :: profile
 real(kind=dp), intent(in)    :: center(1:3),extend(1:3)
 real(kind=dp), intent(in)    :: x(ixI^S,1:ndim)
 real(kind=dp), intent(inout) :: fprofile(ixI^S)
 !.. local ..
 real(dp)                     :: Dist(ixI^S)
 !--------------------------------
 call usr_distance(ixI^L,ixO^L,typeaxial_loc,center,x,dist)
 select case(trim(profile))
   case('gaussian')
     ! where(Dist(ixO^S)<extend(1)*0.2)
     !   fprofile(ixO^S) = 1.0_dp
     ! else where
     !21.04d5 pour extend(1)
     !21.025d5
    fprofile(ixO^S) = 21.025d5*dexp(-Dist(ixO^S)**2.0/(2.0*(extend(1))**2.0_dp))&
                   / (2.0*dpi*(extend(1))**2.0_dp)
     ! end where
  case('tanh')
    fprofile(ixO^S) = tanh(x(ixO^S,r_)/(3*extend(1)))
  case('tanh_radial')
    fprofile(ixO^S) = (1.0_dp-tanh(dabs(x(ixO^S,r_)**2.0_dp)/(5.0*extend(1)**2.0_dp)))/2.0_dp
  case default
    fprofile(ixO^S) =1.0_dp
 end select
end subroutine usr_mat_profile
!--------------------------------------------------------------------

!> procedure will be used to set scalar profile
subroutine usr_mat_profile_scalar(pos_t,standart_deviation,variation_type,Vprofile)
   implicit none
   real(kind=dp) , intent(in)   :: pos_t,standart_deviation
   character(len=*), intent(in) :: variation_type
   real(kind=dp), intent(out)   :: Vprofile
   ! .. local ..
   integer                      :: nstep
   real(kind=dp)                :: local_pos_t
   !----------------------------------------------------------------
      select case(trim(variation_type))
       case('sin')
        Vprofile = dsin(2.0_dp*dpi*pos_t/standart_deviation)
       case('cosin')
        Vprofile = 1.0_dp/dsin(2.0_dp*dpi*pos_t/standart_deviation)
       case('door')
        nstep = floor(pos_t/standart_deviation)
        local_pos_t = pos_t-nstep*standart_deviation
        if(local_pos_t/standart_deviation<0.5_dp) then
          Vprofile = 1.0_dp
        else
          Vprofile = -1.0_dp
        end if
       case('sawtooth')
        nstep = floor(pos_t/standart_deviation)
        local_pos_t = pos_t-nstep*standart_deviation
        Vprofile = -local_pos_t/standart_deviation  ! ma version :
        ! celle qui va de Vmax=150 à Vmin=100=Vmax-DeltaV sur [n*standart_deviation;(n+1)*standart_deviation]
        ! Vprofile = 1.0_dp -2.0_dp*local_pos_t/standart_deviation ! version de Zakaria
       case('sawtooth2')
        local_pos_t = pos_t-0.5_dp*standart_deviation
        Vprofile = 1.0_dp-2.0_dp*modulo(local_pos_t/standart_deviation,1.0_dp)
       case('sawtooth3')
         local_pos_t = pos_t-0.5_dp*standart_deviation
         Vprofile = -1.0_dp+2.0_dp*modulo(local_pos_t/standart_deviation,1.0_dp)
       case default
        Vprofile = 1.0_dp
      end select
end subroutine usr_mat_profile_scalar
!> subroutien to set profile de distance r
subroutine usr_mat_profile_dist(ixI^L,ixO^L,profile, dist,&
                                standart_deviation,fprofile)
 implicit none
 integer,intent(in)           :: ixI^L,ixO^L

 character(len=*), intent(in) :: profile
 real(dp), intent(in)         :: Dist(ixI^S)
 real(dp), intent(in)         :: standart_deviation
 real(kind=dp), intent(inout) :: fprofile(ixI^S)
 !.. local ..

 !--------------------------------

 fprofile(ixO^S) =    1.0_dp
 if(standart_deviation==0.0_dp .and.trim(profile)/='none')then
   call mpistop('standard deviation ==0 at usr_mat_profile_dist in mod_obj_mat.t')
 end if

 select case(trim(profile))
   case('gaussian')
    where(dist(ixO^S)>0.0_dp)
     fprofile(ixO^S) =   dexp(-Dist(ixO^S)**2.0/(2.0*(standart_deviation)**2.0_dp))&
                        / (2.0*dpi*standart_deviation**2.0_dp)
    end where
  case('tanh')
    where(dist(ixO^S)>0.0_dp)
     fprofile(ixO^S) = tanh(Dist(ixO^S)/standart_deviation)
    end where
  case('tanh_inv')
    where(dist(ixO^S)>0.0_dp)
     fprofile(ixO^S) = (1.0_dp-tanh(dabs(Dist(ixO^S)**2.0_dp)/(standart_deviation**2.0_dp)))/2.0_dp
    end where

  case('power4')

    where(dist(ixO^S)>0.0_dp)
     fprofile(ixO^S) = (Dist(ixO^S)/standart_deviation)**4.0_dp
    elsewhere
      fprofile(ixO^S) = 1.0_dp
    end where
  case('linear')

    where(dist(ixO^S)>0.0_dp)
     fprofile(ixO^S) = (Dist(ixO^S)/standart_deviation)
    elsewhere
      fprofile(ixO^S) = 1.0_dp
    end where
  case('sin')
    where(dist(ixO^S)>0.0_dp)
     fprofile(ixO^S) = dsin(Dist(ixO^S)/standart_deviation)
    elsewhere
      fprofile(ixO^S) = 1.0_dp
    end where
  case default
    fprofile(ixO^S) =1.0_dp
 end select
end subroutine usr_mat_profile_dist
!--------------------------------------------------------------------

  !-----------------------------------------------------------------------



  subroutine usr_set_patch_sphere(ixI^L,ixO^L,typeaxial_loc,center,extend,x,patch)
    implicit none
    integer, intent(in)          :: ixI^L,ixO^L
    character(len=*), intent(in) :: typeaxial_loc !< Name of the coordinate system
    real(dp), intent(in)         :: x(ixI^S,1:ndim)
    real(dp), intent(in)         :: center(1:ndim),extend(1:ndim)
    logical, intent(inout)       :: patch(ixI^S)
    !.. local ..
    real(dp)                     :: Dist(ixI^S)
    !-----------------------------------------------
    call usr_distance(ixI^L,ixO^L,typeaxial_loc,center,x,dist)
    patch(ixO^S)=(Dist(ixO^S)<=extend(r_))
  end  subroutine usr_set_patch_sphere
  !-------------------------------------------------------------------------

  subroutine usr_set_patch_conical(ixI^L,ixO^L,typeaxial_loc,center,extend,x,patch)
    implicit none
    integer, intent(in)          :: ixI^L,ixO^L
    character(len=*), intent(in) :: typeaxial_loc !< Name of the coordinate system
    real(dp), intent(in)         :: x(ixI^S,1:ndim)
    real(dp), intent(in)         :: center(:),extend(:)
    logical, intent(inout)       :: patch(ixI^S)

    !-----------------------------------------------



    if(z_in)then
      patch(ixO^S)=dabs(x(ixO^S,r_))<=extend(r_)&
                   +dabs(x(ixO^S,z_))*dtan(extend(theta_))
    else
      patch(ixO^S)=dabs(x(ixO^S,r_))<=extend(r_)
    end if
  end  subroutine usr_set_patch_conical
  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------

  subroutine usr_set_patch_cylinder(ixI^L,ixO^L,typeaxial_loc,center,extend,x,patch)
    use mod_global_parameters, only : z_in,phi_in
    implicit none
    integer, intent(in)          :: ixI^L,ixO^L
    character(len=*), intent(in) :: typeaxial_loc !< Name of the coordinate system
    real(dp), intent(in)         :: x(ixI^S,1:ndim)
    real(dp), intent(in)         :: center(:),extend(:)
    logical, intent(inout)       :: patch(ixI^S)
    !.. local ..
    real(dp)                     :: Dist(ixI^S)
    !-----------------------------------------------
    call usr_distance(ixI^L,ixO^L,typeaxial_loc,center,x,dist)

    patch(ixO^S)=(Dist(ixI^S)<=extend(r_))
    if(z_in)then
      patch(ixO^S)=patch(ixO^S).and.dabs(x(ixO^S,z_)-center(z_))<extend(z_)
    end if
  end  subroutine usr_set_patch_cylinder
  !-------------------------------------------------------------------------

  subroutine usr_set_patch_cube(ixI^L,ixO^L,typeaxial_loc,center,extend,x,patch)
    use mod_global_parameters, only : z_in,phi_in
    implicit none
    integer, intent(in)          :: ixI^L,ixO^L
    character(len=*), intent(in) :: typeaxial_loc !< Name of the coordinate system
    real(dp), intent(in)         :: x(ixI^S,1:ndim)
    real(dp), intent(in)         :: center(:),extend(:)
    logical, intent(inout)       :: patch(ixI^S)
    !.. local ..
    integer                      :: idims
    real(dp)                     :: Dist(ixI^S)
    !-----------------------------------------------
    select case(trim(typeaxial_loc))
    case('slab')
      patch(ixO^S)=dabs(x(ixO^S,1)-center(1))<extend(1)
      if(ndim>=2)then
        Loop_idim : do idims=2,ndim
          patch(ixO^S)=patch(ixO^S).and.dabs(x(ixO^S,idims)-center(idims))<extend(idims)
        end do Loop_idim
      end if
    case('cylindrical')
      patch(ixO^S)=dabs(x(ixO^S,r_)-center(r_))<extend(r_)
      if(z_in)then
        patch(ixO^S)=patch(ixO^S).and.dabs(x(ixO^S,z_)-center(z_))<extend(z_)
      end if
      if(phi_in)then
        patch(ixO^S)=patch(ixO^S).and.dabs(x(ixO^S,r_)*dsin(x(ixO^S,phi_)))&
        -center(r_)*dsin(center(phi_))<dabs(extend(r_)*dsin(extend(phi_)))
      end if

    case('spherical')
      call mpistop('not implimented')
    end select
  end  subroutine usr_set_patch_cube
  !-------------------------------------------------------------------------
    !> Distance between a cells and point 'center'
    subroutine usr_distance(ixI^L,ixO^L,typeaxial_loc,center,x,dist)
      implicit none
      integer, intent(in)          :: ixI^L,ixO^L
      character(len=*), intent(in) :: typeaxial_loc !< Name of the coordinate system
      real(dp), intent(in)         :: x(ixI^S,1:ndim)
      real(dp), intent(in)         :: center(1:ndim)
      real(dp), intent(inout)      :: Dist(ixI^S)
      !.. local ..
      integer                      :: idims
      real(dp)                     :: XX(ixI^S,1:ndim)
      !------------------------------------------------
    select case(trim(typeaxial_loc))
    case('slab')
      FORALL (idims=1:ndim)XX(ixO^S,idims) = x(ixO^S,idims)-center(idims)
      Dist(ixO^S) = dsqrt(SUM(xx(ixO^S,1:ndim)**2.0_DP,dim=ndim+1))
    case('cylindrical')
      if(phi_<=ndim)then
        Dist(ixO^S) =   x(ixO^S,r_)**2.0_dp+center(r_)**2.0_dp - &
                        2.0_dp*x(ixO^S,r_)*center(r_)*dcos(x(ixO^S,phi_)*center(phi_))
      else
        Dist(ixO^S) =   (x(ixO^S,r_)-center(r_))**2.0_dp
      end if
      if(z_<=ndim) then
        Dist(ixO^S) =Dist(ixO^S) +(x(ixO^S,z_)-center(z_))**2.0_dp
      end if

      Dist(ixO^S) = dsqrt(Dist(ixO^S))
    case('spherical')
      if(ndim == 1)    then
        Dist(ixO^S) =   (x(ixO^S,r_)-center(r_))**2.0_dp
      else if(ndim==2) then
        if(phi_<=ndim)then
          Dist(ixO^S) =   x(ixO^S,r_)**2.0_dp+center(r_)**2.0_dp - &
                          2.0_dp*x(ixO^S,r_)*center(r_)&
                          *dcos(x(ixO^S,phi_)*center(phi_))
        else if(z_<=ndim) then
          Dist(ixO^S) =   x(ixO^S,r_)**2.0_dp+center(r_)**2.0_dp - &
                          2.0_dp*x(ixO^S,r_)*center(r_)&
                          *dcos(x(ixO^S,z_)*center(z_))
        end if
      else if(ndim==3) then
        Dist(ixO^S) = x(ixO^S,r_)**2.0_DP+center(r_)**2.0_DP&
                     -2.0_DP*x(ixO^S,r_)*center(r_)*&
                     ( dsin(x(ixO^S,z_))*dsin(center(z_))&
                     *dcos(x(ixO^S,phi_)-center(phi_))&
                      +dcos(x(ixO^S,z_))*dcos(center(z_)) )
      end if
      Dist(ixO^S)=dsqrt(Dist(ixO^S))
    case default
    end select
  end subroutine usr_distance





  !--------------------------------------------------------------------
  subroutine usr_ulrich1976_costheta_zero(ixI^L, ixO^L,x,r_normalized,theta_profile,cos_theta_zero)
   implicit none
   integer, intent(in)             :: ixI^L,ixO^L
   real(kind=dp), intent(in)       :: x(ixI^S,1:ndim),r_normalized(ixI^S)
   real(kind=dp), intent(in)       :: theta_profile(ixI^S)
   real(kind=dp), intent(inout)    :: cos_theta_zero(ixI^S)
   real(kind=dp), dimension(1:ndim) :: zero_dim
   !----------------------------------------------------

   {zero_dim(^D)=0.0_dp\}!-dx(1,1)

   cos_theta_zero(ixO^S) = 0.0_dp

     where(dabs(r_normalized(ixO^S)-1.0_dp) < smalldouble)

       cos_theta_zero(ixO^S) = (DCOS(theta_profile(ixO^S)))**(1.0_dp/3.0_dp)

     elsewhere (r_normalized(ixO^S) > 1.0_dp)

       cos_theta_zero(ixO^S) = 2.0_dp * ( ( ( r_normalized(ixO^S) - 1.0_dp ) / &
       3.0_dp ) ** ( 0.5_dp ) ) * DSINH( ( 1.0_dp / 3.0_dp ) * &
       DASINH( ( r_normalized(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / &
       ( 2.0_dp * ( ( ( r_normalized(ixO^S) - 1.0_dp ) / &
       3.0_dp ) ** ( 1.5_dp ) ) ) ) )

     elsewhere (r_normalized(ixO^S) < 1.0_dp)
       where ((((r_normalized(ixO^S)*DCOS(theta_profile(ixO^S)))/2.0_dp)**2.0_dp)-&
       (((1.0_dp-r_normalized(ixO^S))/3.0_dp)**3.0_dp)>0.0_dp)

         cos_theta_zero(ixO^S) = 2.0_dp * ( ( ( 1.0_dp - r_normalized(ixO^S) ) / &
         3.0_dp ) ** ( 0.5_dp ) ) * DCOSH( ( 1.0_dp / 3.0_dp ) * &
         DACOSH( ( r_normalized(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / &
         ( 2.0_dp * ( ( ( 1.0_dp - r_normalized(ixO^S) ) / &
         3.0_dp ) ** ( 1.5_dp ) ) ) ) )

       elsewhere ((((r_normalized(ixO^S)*DCOS(theta_profile(ixO^S)))/2.0_dp)**2.0_dp)-&
       (((1.0_dp-r_normalized(ixO^S))/3.0_dp)**3.0_dp)<0.0_dp)

         cos_theta_zero(ixO^S) = 2.0_dp * ( ( ( 1.0_dp - r_normalized(ixO^S) ) / &
         3.0_dp ) ** ( 0.5_dp ) ) * DCOS( ( 1.0_dp / 3.0_dp ) * &
         DACOS( ( r_normalized(ixO^S) * DCOS( theta_profile(ixO^S) ) ) / &
         ( 2.0_dp * ( ( ( 1.0_dp - r_normalized(ixO^S) ) / &
         3.0_dp ) ** ( 1.5_dp ) ) ) ) )
       elsewhere(dabs(r_normalized(ixO^S)-1.0_dp) < smalldouble)

         cos_theta_zero(ixO^S) = (DCOS(theta_profile(ixO^S)))**(1.0_dp/3.0_dp)
       end where
     end where

  end subroutine usr_ulrich1976_costheta_zero

  !--------------------------------------------------------------------
  subroutine usr_get_theta(ixI^L, ixO^L,x,theta_profile,specialorigintheta)
   implicit none
   integer, intent(in)             :: ixI^L,ixO^L
   real(kind=dp), intent(in)       :: x(ixI^S,1:ndim)
   real(kind=dp), intent(inout)    :: theta_profile(ixI^S)
   real(kind=dp)                   :: dstce(ixI^S)
   logical,optional                :: specialorigintheta
   !----------------------------------------------------



   select case(typeaxial)
   case('spherical')
        theta_profile(ixO^S)          = x(ixO^S,theta_)
   case('cylindrical')
      !Do not use arctan approach which is tricky with z<0
      !but use arccos instead !!
      if(ndim>1)then
            where(DABS(x(ixO^S,r_))<=smalldouble.and.DABS(x(ixO^S,z_))<=smalldouble)
                  theta_profile(ixO^S)        = 0.0_dp
            elsewhere
              dstce(ixO^S)=dsqrt((x(ixO^S,r_)**2.0_dp)+(x(ixO^S,z_)**2.0_dp))
              where(x(ixO^S,z_)>smalldouble)
                theta_profile(ixO^S)        = DACOS(x(ixO^S,z_)/dstce(ixO^S))
              elsewhere ! negative case
                theta_profile(ixO^S)        = dpi-DACOS(x(ixO^S,z_)/dstce(ixO^S))
              end where
            end where

      else
          theta_profile(ixO^S)        = 0.0_dp
      end if
   case('slab','slabstretch')
     where(DABS(x(ixO^S,x_))<=smalldouble.and.DABS(x(ixO^S,y_))<=smalldouble)
       theta_profile(ixO^S)        = 0.0_dp
     elsewhere
       where((x(ixO^S,x_)/x(ixO^S,y_))>smalldouble)
         theta_profile(ixO^S)        = DATAN(DABS(x(ixO^S,x_)/x(ixO^S,y_)))
       elsewhere ! negative case
         theta_profile(ixO^S)        = -DATAN(DABS(x(ixO^S,x_)/x(ixO^S,y_)))
       end where
     end where
   case default
      call mpistop('Unknown typeaxial')
   end select

  end subroutine usr_get_theta





  !------------------------------------------------------------------------
  subroutine usr_get_spherical(ixI^L,ixO^L,typeaxial_loc,center,x,x_sphere)
    use mod_global_parameters, only : z_in,phi_in
    implicit none
    integer, intent(in)          :: ixI^L,ixO^L
    character(len=*), intent(in) :: typeaxial_loc !< Name of the coordinate system
    real(dp), intent(in)         :: x(ixI^S,1:ndim)
    real(dp), intent(in)         :: center(1:ndim)
    real(dp), intent(inout)      :: x_sphere(ixI^S,1:ndim)
    !.. local ..
    integer                      :: idims
    real(dp)                     :: XX(ixI^S,1:ndim),Dist(ixI^S)
      !------------------------------------------------

    call usr_distance(ixI^L,ixO^L,typeaxial_loc,center,x,x_sphere(ixI^S,r_))

    select case(trim(typeaxial_loc))
    case('slab')
      FORALL (idims=1:ndim)XX(ixO^S,idims) = x(ixO^S,idims)-center(idims)
      if(z_in)then
       where(x_sphere(ixO^S,r_)>smalldouble)
         x_sphere(ixO^S,theta_) = dacos((x(ixO^S,z_)-center(z_))/x_sphere(ixO^S,r_))
       elsewhere
         x_sphere(ixO^S,theta_) = 0.0_dp
       end where
      end if
    case('cylindrical')
      if(phi_<=ndim)then
       where(x_sphere(ixO^S,r_)>smalldouble)
        x_sphere(ixO^S,phi_) = dacos((x(ixO^S,x_)-center(x_))/x_sphere(ixO^S,r_))
       elsewhere
         x_sphere(ixO^S,phi_) = 0.0_dp
       end where
      end if

      if(z_in) then
       where(x_sphere(ixO^S,r_)>smalldouble)
        x_sphere(ixO^S,theta_) = dacos((x(ixO^S,z_)-center(z_))/x_sphere(ixO^S,r_))
       elsewhere
         x_sphere(ixO^S,theta_) = 0.0_dp
       end where
      end if


    case('spherical')
      if(ndim == 1)    then
       ! dummy
      else if(ndim==2) then
        if(phi_in)then
          x_sphere(ixO^S,phi_)   = x(ixO^S,phi_)
        else if(z_in) then
          x_sphere(ixO^S,theta_)   = x(ixO^S,theta_)
        end if
      else if(ndim==3) then
       x_sphere(ixO^S,theta_)   = x(ixO^S,theta_)
       x_sphere(ixO^S,phi_)     = x(ixO^S,phi_)
      end if

    case default
    end select
  end subroutine usr_get_spherical
  !> subroutine get the volume
  subroutine usr_get_volume(extend,shape, volume)
   implicit none
   real(dp), intent(in)         :: extend(:)
   character(len=*), intent(in) :: shape
   real(dp), intent(inout)      :: volume
   !---------------------------------------------------
    select case(shape)
    case('sphere')
     volume = 4.0_dp/3.0_dp*dpi*extend(1)**3.0_dp
   case ('cube')
      select case(ndim)
       case(1)
        volume =extend(1)
       case(2)
        volume =(extend(1)*extend(2))
       case(3)
        volume =extend(1)*extend(2)*extend(3)
     end select
   case('cylinder')
     select case(ndim)
       case(1)
        volume = dpi*extend(r_)**2.0_dp
      case(2)
        if(phi_<=ndim)then
          volume = dpi*extend(r_)**2.0_dp
        elseif(z_<=ndim) then
          volume = dpi*extend(r_)**2.0_dp *extend(z_)
        end if
      case(3)
        volume = dpi*extend(r_)**2.0_dp *extend(z_)
     end select
   case default
     volume = 1.0d0
   end select
  end subroutine usr_get_volume
  !----------------------------------------------------------------------
!------------------------------------------------------------------------------------------
  !> subroutine to compute large variation in of an array
  subroutine usr_meanvalue_of_array(ixI^L,ixO^L,array,mean_array,mean_flag)
    use mod_global_parameters
    use mod_small_values
    implicit none
    integer, intent(in)             :: ixI^L, ixO^L
    real(kind=dp), intent(in)       :: array(ixI^S)
    real(kind=dp), intent(out)      :: mean_array(ixI^S)
    logical, intent(out)            :: mean_flag(ixI^S)
    ! .. local ..
    integer                         :: ix^D,kxO^L,n_cell_daverage,n_cell_true
    real(kind=dp)                   :: sum_array,usr_max_deviation
    !----------------------------------------------------------------
    usr_max_deviation = 100.0d0
    n_cell_daverage=2
    mean_flag(ixI^S)  = .true.

    {do ix^DB= ixO^LIM^DB\}

      {kxOmin^D= max(ix^D-n_cell_daverage, ixImin^D);
      kxOmax^D= min(ix^D+n_cell_daverage, ixImax^D);\}

      SUM_array=sum(array(kxO^S),mean_flag(kxO^S))

      n_cell_true = count(mean_flag(kxO^S))

      if(n_cell_true>1) then
        mean_array(ix^D)=(sum_array-array(ix^D))&
           /(count(mean_flag(kxO^S))-1.0_dp)
      else
        mean_array(ix^D)= array(ix^D)
      end if

      if(dabs(array(ix^D))>smalldouble&
         .and.dabs(mean_array(ix^D))>smalldouble)then

          mean_flag(ix^D)=dabs((array(ix^D)-mean_array(ix^D))/mean_array(ix^D))&
                           <usr_max_deviation
      else
          mean_flag(ix^D)=.true.
      end if

    {end do^D&\}

  end subroutine usr_meanvalue_of_array
  !------------------------------------------------------------------------------------------


  !------------------------------------------------------------------------------------------
  !> subroutine to compute large variation in of an array
  subroutine usr_medianvalue_of_array(ixI^L,ixO^L,array,median_array,indices_median,&
                                       mean_flag)
    use mod_global_parameters
    use mod_small_values
    implicit none
    integer, intent(in)             :: ixI^L, ixO^L
    real(kind=dp), intent(in)       :: array(ixI^S)
    real(kind=dp), intent(out)      :: median_array(ixI^S)
    logical, intent(out)            :: mean_flag(ixI^S)
    integer, intent(out)            :: indices_median(ixI^S,1:ndim,2)
    ! .. local ..
    integer                         :: ix^D,kxO^L,n_cell_daverage,n_cell_true
    integer                         :: loc_indices(1:3,1:2)
    real(kind=dp)                   :: median_value,usr_max_deviation
    real(kind=dp)                   :: usr_min_deviation_abs,usr_max_variation_abs
    logical                         :: logic_sort
    !----------------------------------------------------------------
    usr_max_deviation     = 2.0d0
    usr_min_deviation_abs = 1.0e-4
    usr_max_variation_abs = 2.0d-2
    n_cell_daverage       = 1
    mean_flag(ixI^S)      = .true.




    {do ix^DB= ixO^LIM^DB\}

      {kxOmin^D= max(ix^D-n_cell_daverage, ixImin^D);
      kxOmax^D= min(ix^D+n_cell_daverage, ixImax^D);\}
      if(dabs(array(ix^D))<usr_min_deviation_abs)then
        median_array(ix^D)=array(ix^D)
        mean_flag(ix^D)   =.true.
        cycle
      end if
      if(maxval(array(kxO^S))-minval(array(kxO^S))<usr_max_variation_abs)then
        median_array(ix^D)=array(ix^D)
        mean_flag(ix^D)   =.true.
        cycle
      end if

      call usr_median_from_w(ixI^L,kxO^L,array,median_value,&
                             loc_indices,logic_sort)

      if(logic_sort)then
        indices_median(ix^D,1:ndim,1:2) =  loc_indices(1:ndim,1:2)
         median_array(ix^D) = median_value

       if(dabs(array(ix^D))>dabs(median_array(ix^D))&
         .and.dabs(median_value)>smalldouble)then

          mean_flag(ix^D)=dabs((array(ix^D)-median_value)/median_value)&
                           <usr_max_deviation

       else if(dabs(median_value)<smalldouble) then
          mean_flag(ix^D)=dabs(array(ix^D))<100.0_dp
       else
          mean_flag(ix^D)=.true.
       end if
      else
        median_array(ix^D) =array(ix^D)
      end if
    {end do^D&\}

end subroutine usr_medianvalue_of_array
  !------------------------------------------------------------------------------------------
  subroutine usr_median_from_w(ixI^L,ixO^L,w_array,median_value,indice_median,logic_sort)
    use mod_global_parameters
    use mod_small_values
    implicit none
    integer, intent(in)             :: ixI^L, ixO^L
    real(kind=dp), intent(in)       :: w_array(ixI^S)
    real(kind=dp), intent(out)      :: median_value
    integer, intent(out)            :: indice_median(3,2)
    logical, intent(out)            :: logic_sort
    ! .. local ..
    integer                         :: n_point, n_half,n_R,n_colonne
    integer                         :: n_cells(ndim),ix_half(ndim)
    real(kind=dp)                   :: array({^D&(ixOmax^D-ixOmin^D+1)*})
    integer                         :: array_ind({^D&(ixOmax^D-ixOmin^D+1)*})
    !----------------------------------------------------
     n_point={^D&(ixOmax^D-ixOmin^D+1)*}
     array=reshape(w_array(ixO^S),(/n_point/))
     call usr_median_array(n_point,array,array_ind,median_value,logic_sort)
     if(.not.logic_sort)return

     {^D&n_cells(^D)=ixOmax^D-ixOmin^D+1;}
     n_half    = n_point/2

     n_colonne = array_ind(n_half)/n_cells(1)
     if(mod(array_ind(n_half),n_cells(1))==0) then
      indice_median(1,1) = n_cells(1)
      indice_median(2,1) = n_colonne
     else
      indice_median(1,1) = array_ind(n_half)-n_colonne*n_cells(1)
      indice_median(2,1) = n_colonne +1
     end if

     if(mod(n_point,2)==0) then
      n_R       = n_half+1
      n_colonne = array_ind(n_R)/n_cells(1)
      if(mod(array_ind(n_R),n_cells(1))==0) then
       indice_median(1,2) = n_cells(1)
       indice_median(2,2) = n_colonne
      else
       indice_median(1,1) = array_ind(n_R)-n_colonne*n_cells(1)
       indice_median(2,1) = n_colonne +1
      end if
     else
      indice_median(:,2) = indice_median(:,1)
     end if
  end subroutine usr_median_from_w

  !------------------------------------------------------------------------------------------
  subroutine usr_median_array(n_point,array,array_ind,median_value,logic_sort)
    implicit none
    integer, intent(in)             :: n_point
    real(kind=dp), intent(inout)    :: array(n_point)
    real(kind=dp), intent(out)      :: median_value
    integer, intent(out)            :: array_ind(n_point)
    logical, intent(out)            :: logic_sort
    ! .. local ..

    integer                         :: i
    !----------------------------------------------------
    array_ind=(/(i,i=1,n_point)/)
    call  usr_sort_array(n_point,array,array_ind,logic_sort)
    if(mod(n_point,2)==0)then
      median_value  =(array(n_point/2+1)+array(n_point/2-1))/2.0_dp
    else
      median_value=array(n_point/2+1)
    end if

  end subroutine usr_median_array
  !-------------------------------------------------
  subroutine usr_sort_array(n_point,array,array_ind,logic_sort)

    implicit none
    integer, intent(in)             :: n_point
    real(kind=dp), intent(inout)    :: array(n_point)
    integer, intent(out)            :: array_ind(n_point)
    logical, intent(out)            :: logic_sort

    ! .. local ..
    real(kind=dp)                   :: temp_data
    integer                         :: i, iter,temp_indice
    !----------------------------------------------------
    if(maxval(array(:))-minval(array(:))<smalldouble)then
     logic_sort=.false.
    else
     logic_sort=.true.
     iter=0
    Loop_while : do while(iter<n_point-2)
      Loop_sort : do i=1,n_point-1
       if(array(i)>array(i+1))then
        temp_data = array(i)
        array(i)  = array(i+1)
        array(i+1)=temp_data
        temp_indice    = array_ind(i)
        array_ind(i)   = array_ind(i+1)
        array_ind(i+1) = temp_indice
        iter = 0
       else
        iter=iter+1
       end if
      end do Loop_sort
     end do Loop_while
    end if
  end subroutine usr_sort_array



  !----------------------------------------------------------------------
  !> subroutine for lorentz transrmation to add proper motion
  subroutine usr_lorentz_transformation_add_proper_speed(ixI^L,ixO^L,patch,&
    velocity_proper,v)
    ! Eqaution in use https://en.wikipedia.org/wiki/Lorentz_transformation in Transformation of velocities
    implicit none
    integer, intent(in)           :: ixI^L,ixO^L
    logical, intent(in)           :: patch(ixI^S)
    real(kind=dp), intent(in)     :: velocity_proper(:)
    real(kind=dp), intent(inout)  :: v(ixI^S,1:ndir)
    !.. local ..
    integer                    :: idir
    real(kind=dp)              :: uv(ixO^S)
    real(kind=dp)              :: lfac_proper
    !-------------------------------------------------------------
    if(all(velocity_proper==0.0_dp))return
    where(patch(ixO^S))
     uv(ixO^S)  = v(ixO^S,1)*velocity_proper(1)
    endwhere
    cond_idir2: if(idir>1) then
     Loop_idir : do idir = 2,ndir
      where(patch(ixO^S))
       uv(ixO^S) = uv(ixO^S) + v(ixO^S,idir)*velocity_proper(idir)
      endwhere
     end do Loop_idir
    end if cond_idir2
    lfac_proper   = 1.0_dp/dsqrt(1.0_dp - SUM(velocity_proper**2.0_dp))

    Loop_idir_newv : do idir = 2,ndir
     where(patch(ixO^S))
      v(ixO^S,idir) = 1.0_dp/(1.0_dp+uv(ixO^S))*(v(ixO^S,idir)/lfac_proper&
                    +velocity_proper(idir)&
                    +lfac_proper/(lfac_proper+1.0_dp)*uv(ixO^S)*velocity_proper(idir) )
     endwhere
    end do Loop_idir_newv
  end subroutine usr_lorentz_transformation_add_proper_speed

  subroutine usr_mat_get_Lohner_error(ixI^L,ixO^L,level,ivar,w,error_lohner)
    use mod_global_parameters
    implicit none
    integer, intent(in)             :: ixI^L,ixO^L,level,ivar
    real(kind=dp), intent(in)       :: w(ixI^S,1:nw)
    real(kind=dp), intent(inout)    :: error_lohner(ixM^T)
    ! .. local ..
    integer                         :: ix^L,hx^L, jx^L, h2x^L, j2x^L, ix^D
    integer                         :: idims,idims2
    real(kind=dp)                   :: epsilon
    real(kind=dp), dimension(ixG^T) :: tmp,tmp1,tmp2
    real(kind=dp), dimension(ixM^T) :: numerator,denominator
    !-------------------------------------------------------
     epsilon = 1.0d-6
     ix^L=ixO^L^LADD1;

     error_lohner=zero
     numerator=zero
     denominator=zero

     Loop_idims_num : do idims=1,ndim
        hx^L=ix^L-kr(^D,idims);
        jx^L=ix^L+kr(^D,idims);
        if (ivar<=nw) then
          if (logflag(ivar)) then
            tmp(ix^S)=dlog10(w(jx^S,ivar))-dlog10(w(hx^S,ivar))
          else
            tmp(ix^S)=w(jx^S,ivar)-w(hx^S,ivar)
          end if
        else
          if (logflag(ivar)) then
            tmp(ix^S)=dlog10(tmp1(jx^S))-dlog10(tmp1(hx^S))
          else
            tmp(ix^S)=tmp1(jx^S)-tmp1(hx^S)
          end if
        end if
        Loop_idims1: do idims2=1,ndim
           h2x^L=ixM^LL-kr(^D,idims2);
           j2x^L=ixM^LL+kr(^D,idims2);
           numerator=numerator+(tmp(j2x^S)-tmp(h2x^S))**2.0d0
        end do Loop_idims1
     end do Loop_idims_num

     Loop_idims_dem : do idims=1,ndim
        if (ivar<=nw) then
           if (logflag(ivar)) then
            tmp=dabs(dlog10(w(ixG^T,ivar)))
           else
            tmp=dabs(w(ixG^T,ivar))
           end if
        else
           if (logflag(ivar)) then
            tmp=dabs(dlog10(tmp1(ixG^T)))
           else
            tmp=dabs(tmp1(ixG^T))
           end if
        end if
        hx^L=ix^L-kr(^D,idims);
        jx^L=ix^L+kr(^D,idims);
        tmp2(ix^S)=tmp(jx^S)+tmp(hx^S)
        hx^L=ixM^LL-2*kr(^D,idims);
        jx^L=ixM^LL+2*kr(^D,idims);
        if (ivar<=nw) then
          if (logflag(ivar)) then
            tmp(ixM^T)=dabs(dlog10(w(jx^S,ivar))&
                       -dlog10(w(ixM^T,ivar))) &
                       +dabs(dlog10(w(ixM^T,ivar))&
                       -dlog10(w(hx^S,ivar)))
          else
             tmp(ixM^T)=dabs(w(jx^S,ivar)-w(ixM^T,ivar)) &
                        +dabs(w(ixM^T,ivar)-w(hx^S,ivar))
          end if
        else
          if (logflag(ivar)) then
            tmp(ixM^T)=dabs(dlog10(tmp1(jx^S))-dlog10(tmp1(ixM^T))) &
                      +dabs(dlog10(tmp1(ixM^T))-dlog10(tmp1(hx^S)))
          else
             tmp(ixM^T)=dabs(tmp1(jx^S)-tmp1(ixM^T)) &
                        +dabs(tmp1(ixM^T)-tmp1(hx^S))
          end if
        end if
        Loop_idims2: do idims2=1,ndim
           h2x^L=ixM^LL-kr(^D,idims2);
           j2x^L=ixM^LL+kr(^D,idims2);
           denominator=denominator &
                      +(tmp(ixM^T)+amr_wavefilter(level)*(tmp2(j2x^S)+tmp2(h2x^S)))**2
        end do Loop_idims2
     end do Loop_idims_dem
     error_lohner=error_lohner+dsqrt(numerator/max(denominator,epsilon))

  end subroutine usr_mat_get_Lohner_error

  subroutine usr_mat_profile_tanh_scalar_maxdist(dist_impos,&
                                  dist_experance,dist_ecart,coef)
    implicit none

    real(kind=dp), intent(in)          :: dist_impos
    real(kind=dp), intent(in)          :: dist_experance,dist_ecart
    real(kind=dp), intent(out)         :: coef
    coef = (1.0_dp-tanh((dist_impos-dist_experance)/ dist_ecart))/2.0_dp

  end subroutine usr_mat_profile_tanh_scalar_maxdist


  subroutine usr_mat_read_error_message(i_error_read,myindice,objname)

    implicit none
    integer, intent(in)             :: i_error_read,myindice
    character(len=*), intent(in)    ::  objname
    ! .. local ..
    character(len=78)               :: modname
    character(len=78)               :: thelistname
    !---------------------------------------------------
    modname='mod_obj_'//trim(objname)//'.t'
    thelistname='usr_'//trim(objname)//'list'
    write(thelistname,'(A,I1)')trim(thelistname),myindice

    if(i_error_read>0)then
      write(*,*)'At user side in ', modname,': error at reading parfile'
      write(*,*)'Check input.  Something was wrong, it will stop'
      call mpistop('It stops at reading the parfile')
      elseif(i_error_read<0)then
      write(*,*)'At user side in ', modname,': error at reading parfile'
      write(*,*)'Reach the end of the file it will stop'
      call mpistop('It stops at reading the parfile')
    else
      if(mype==0)write(*,*) 'End of Reading ' ,thelistname
    end if
  end subroutine usr_mat_read_error_message
end module mod_obj_mat
