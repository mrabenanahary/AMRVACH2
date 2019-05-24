module mod_obj_mat
use mod_constants
use mod_obj_global_parameters
implicit none
contains


!--------------------------------------------------------------------
!> Subtroutine for power law distribution,here the dust sizes are defined. Ndust bins, with all bins having equal total mass.

subroutine usr_mat_powerlaw_withX(n_point,power_a,min_var,max_var,   var_r)
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
  var_r(i)  = -1.0_dp/(power_a-1.0_dp)*(r(i)**(power_a+&
     2.0_dp)-r(i-1)**(power_a+2.0_dp))/(r(i)**(power_a+1)-r(i-1)**(power_a+1))
 end do Loop_point
end subroutine usr_mat_powerlaw_withX

!--------------------------------------------------------------------
!> subroutien to set profile de distance r
subroutine usr_mat_profile(ixImin1,ixImax1,ixOmin1,ixOmax1,typeaxial_loc,&
   profile, center,extend,x,fprofile)
 implicit none
 integer,intent(in)           :: ixImin1,ixImax1,ixOmin1,ixOmax1
 character(len=*), intent(in) :: typeaxial_loc !< Name of the coordinate system
 character(len=*), intent(in) :: profile
 real(kind=dp), intent(in)    :: center(1:ndim),extend(1:ndim)
 real(kind=dp), intent(in)    :: x(ixImin1:ixImax1,1:ndim)
 real(kind=dp), intent(inout) :: fprofile(ixImin1:ixImax1)
 !.. local ..
 real(dp)                     :: Dist(ixImin1:ixImax1)
 !--------------------------------
 call usr_distance(ixImin1,ixImax1,ixOmin1,ixOmax1,typeaxial_loc,center,x,&
    dist)
 select case(trim(profile))
   case('gaussian')
     ! where(Dist(ixO^S)<extend(1)*0.2)
     !   fprofile(ixO^S) = 1.0_dp
     ! else where
     !21.04d5 pour extend(1)
     !21.025d5
    fprofile(ixOmin1:ixOmax1) = 21.025d5*dexp(-&
       Dist(ixOmin1:ixOmax1)**2.0/(2.0*(extend(1))**2.0_dp))/ &
       (2.0*dpi*(extend(1))**2.0_dp)
     ! end where
  case('tanh')
    fprofile(ixOmin1:ixOmax1) = tanh(x(ixOmin1:ixOmax1,r_)/(3*extend(1)))
   case default
    fprofile(ixOmin1:ixOmax1) =1.0_dp
 end select
end subroutine usr_mat_profile
!--------------------------------------------------------------------
  !-----------------------------------------------------------------------



  subroutine usr_set_patch_sphere(ixImin1,ixImax1,ixOmin1,ixOmax1,&
     typeaxial_loc,center,extend,x,patch)
    implicit none
    integer, intent(in)          :: ixImin1,ixImax1,ixOmin1,ixOmax1
    character(len=*), intent(in) :: typeaxial_loc !< Name of the coordinate system
    real(dp), intent(in)         :: x(ixImin1:ixImax1,1:ndim)
    real(dp), intent(in)         :: center(1:ndim),extend(1:ndim)
    logical, intent(inout)       :: patch(ixImin1:ixImax1)
    !.. local ..
    real(dp)                     :: Dist(ixImin1:ixImax1)
    !-----------------------------------------------
    call usr_distance(ixImin1,ixImax1,ixOmin1,ixOmax1,typeaxial_loc,center,x,&
       dist)
    patch(ixOmin1:ixOmax1)=(Dist(ixOmin1:ixOmax1)<=extend(r_))
  end  subroutine usr_set_patch_sphere
  !-------------------------------------------------------------------------

  subroutine usr_set_patch_conical(ixImin1,ixImax1,ixOmin1,ixOmax1,&
     typeaxial_loc,center,extend,x,patch)
    implicit none
    integer, intent(in)          :: ixImin1,ixImax1,ixOmin1,ixOmax1
    character(len=*), intent(in) :: typeaxial_loc !< Name of the coordinate system
    real(dp), intent(in)         :: x(ixImin1:ixImax1,1:ndim)
    real(dp), intent(in)         :: center(:),extend(:)
    logical, intent(inout)       :: patch(ixImin1:ixImax1)

    !-----------------------------------------------



    if(z_in)then
      patch(ixOmin1:ixOmax1)=dabs(x(ixOmin1:ixOmax1,&
         r_))<=extend(r_)+dabs(x(ixOmin1:ixOmax1,z_))*dtan(extend(theta_))
    else
      patch(ixOmin1:ixOmax1)=dabs(x(ixOmin1:ixOmax1,r_))<=extend(r_)
    end if
  end  subroutine usr_set_patch_conical
  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------

  subroutine usr_set_patch_cylinder(ixImin1,ixImax1,ixOmin1,ixOmax1,&
     typeaxial_loc,center,extend,x,patch)
    implicit none
    integer, intent(in)          :: ixImin1,ixImax1,ixOmin1,ixOmax1
    character(len=*), intent(in) :: typeaxial_loc !< Name of the coordinate system
    real(dp), intent(in)         :: x(ixImin1:ixImax1,1:ndim)
    real(dp), intent(in)         :: center(:),extend(:)
    logical, intent(inout)       :: patch(ixImin1:ixImax1)
    !.. local ..
    real(dp)                     :: Dist(ixImin1:ixImax1)
    !-----------------------------------------------
    call usr_distance(ixImin1,ixImax1,ixOmin1,ixOmax1,typeaxial_loc,center,x,&
       dist)

    patch(ixOmin1:ixOmax1)=(Dist(ixImin1:ixImax1)<=extend(r_))
    if(z_in)then
      patch(ixOmin1:ixOmax1)=patch(ixOmin1:ixOmax1).and.dabs(x(ixOmin1:ixOmax1,&
         z_)-center(z_))<extend(z_)
    end if
  end  subroutine usr_set_patch_cylinder
  !-------------------------------------------------------------------------
  subroutine usr_set_patch_cube(ixImin1,ixImax1,ixOmin1,ixOmax1,typeaxial_loc,&
     center,extend,x,patch)
    implicit none
    integer, intent(in)          :: ixImin1,ixImax1,ixOmin1,ixOmax1
    character(len=*), intent(in) :: typeaxial_loc !< Name of the coordinate system
    real(dp), intent(in)         :: x(ixImin1:ixImax1,1:ndim)
    real(dp), intent(in)         :: center(:),extend(:)
    logical, intent(inout)       :: patch(ixImin1:ixImax1)
    !.. local ..
    integer                      :: idims
    real(dp)                     :: Dist(ixImin1:ixImax1)
    !-----------------------------------------------
    select case(trim(typeaxial_loc))
    case('slab')
      patch(ixOmin1:ixOmax1)=dabs(x(ixOmin1:ixOmax1,1)-center(1))<extend(1)
      if(ndim>=2)then
        Loop_idim : do idims=2,ndim
          patch(ixOmin1:ixOmax1)=patch(ixOmin1:ixOmax1).and.&
             dabs(x(ixOmin1:ixOmax1,idims)-center(idims))<extend(idims)
        end do Loop_idim
      end if
    case('cylindrical')
      patch(ixOmin1:ixOmax1)=dabs(x(ixOmin1:ixOmax1,r_)-center(r_))<extend(r_)
      if(z_in)then
        patch(ixOmin1:ixOmax1)=patch(ixOmin1:ixOmax1).and.&
           dabs(x(ixOmin1:ixOmax1,z_)-center(z_))<extend(z_)
      end if
      if(phi_in)then
        patch(ixOmin1:ixOmax1)=patch(ixOmin1:ixOmax1).and.&
           dabs(x(ixOmin1:ixOmax1,r_)*dsin(x(ixOmin1:ixOmax1,&
           phi_)))-center(r_)*dsin(center(phi_))<dabs(extend(r_)*dsin(extend(&
           phi_)))
      end if

    case('spherical')
      call mpistop('not implimented')
    end select
  end  subroutine usr_set_patch_cube
  !-------------------------------------------------------------------------
    !> Distance between a cells and point 'center'
    subroutine usr_distance(ixImin1,ixImax1,ixOmin1,ixOmax1,typeaxial_loc,&
       center,x,dist)
      implicit none
      integer, intent(in)          :: ixImin1,ixImax1,ixOmin1,ixOmax1
      character(len=*), intent(in) :: typeaxial_loc !< Name of the coordinate system
      real(dp), intent(in)         :: x(ixImin1:ixImax1,1:ndim)
      real(dp), intent(in)         :: center(1:ndim)
      real(dp), intent(inout)      :: Dist(ixImin1:ixImax1)
      !.. local ..
      integer                      :: idims
      real(dp)                     :: XX(ixImin1:ixImax1,1:ndim)
      !------------------------------------------------
    select case(trim(typeaxial_loc))
    case('slab')
      FORALL (idims=1:ndim)XX(ixOmin1:ixOmax1,idims) = x(ixOmin1:ixOmax1,&
         idims)-center(idims)
      Dist(ixOmin1:ixOmax1) = dsqrt(SUM(xx(ixOmin1:ixOmax1,1:ndim)**2.0_DP,&
         dim=ndim+1))
    case('cylindrical')
      if(phi_<=ndim)then
        Dist(ixOmin1:ixOmax1) =   x(ixOmin1:ixOmax1,&
           r_)**2.0_dp+center(r_)**2.0_dp - 2.0_dp*x(ixOmin1:ixOmax1,&
           r_)*center(r_)*dcos(x(ixOmin1:ixOmax1,phi_)*center(phi_))
      else
        Dist(ixOmin1:ixOmax1) =   (x(ixOmin1:ixOmax1,r_)-center(r_))**2.0_dp
      end if
      if(z_<=ndim) then
        Dist(ixOmin1:ixOmax1) =Dist(ixOmin1:ixOmax1) +(x(ixOmin1:ixOmax1,&
           z_)-center(z_))**2.0_dp
      end if

      Dist(ixOmin1:ixOmax1) = dsqrt(Dist(ixOmin1:ixOmax1))
    case('spherical')
      if(ndim == 1)    then
        Dist(ixOmin1:ixOmax1) =   (x(ixOmin1:ixOmax1,r_)-center(r_))**2.0_dp
      else if(ndim==2) then
        if(phi_<=ndim)then
          Dist(ixOmin1:ixOmax1) =   x(ixOmin1:ixOmax1,&
             r_)**2.0_dp+center(r_)**2.0_dp - 2.0_dp*x(ixOmin1:ixOmax1,&
             r_)*center(r_)*dcos(x(ixOmin1:ixOmax1,phi_)*center(phi_))
        else if(z_<=ndim) then
          Dist(ixOmin1:ixOmax1) =   x(ixOmin1:ixOmax1,&
             r_)**2.0_dp+center(r_)**2.0_dp - 2.0_dp*x(ixOmin1:ixOmax1,&
             r_)*center(r_)*dcos(x(ixOmin1:ixOmax1,z_)*center(z_))
        end if
      else if(ndim==3) then
        Dist(ixOmin1:ixOmax1) = x(ixOmin1:ixOmax1,&
           r_)**2.0_DP+center(r_)**2.0_DP&
                     -2.0_DP*x(ixOmin1:ixOmax1,&
                        r_)*center(r_)*( dsin(x(ixOmin1:ixOmax1,&
                        z_))*dsin(center(z_))*dcos(x(ixOmin1:ixOmax1,&
                        phi_)-center(phi_))+dcos(x(ixOmin1:ixOmax1,&
                        z_))*dcos(center(z_)) )
      end if
      Dist(ixOmin1:ixOmax1)=dsqrt(Dist(ixOmin1:ixOmax1))
    case default
    end select
  end subroutine usr_distance
  !------------------------------------------------------------------------
  subroutine usr_get_spherical(ixImin1,ixImax1,ixOmin1,ixOmax1,typeaxial_loc,&
     center,x,x_sphere)
    implicit none
    integer, intent(in)          :: ixImin1,ixImax1,ixOmin1,ixOmax1
    character(len=*), intent(in) :: typeaxial_loc !< Name of the coordinate system
    real(dp), intent(in)         :: x(ixImin1:ixImax1,1:ndim)
    real(dp), intent(in)         :: center(1:ndim)
    real(dp), intent(inout)      :: x_sphere(ixImin1:ixImax1,1:ndim)
    !.. local ..
    integer                      :: idims
    real(dp)                     :: XX(ixImin1:ixImax1,1:ndim),&
       Dist(ixImin1:ixImax1)
      !------------------------------------------------

    call usr_distance(ixImin1,ixImax1,ixOmin1,ixOmax1,typeaxial_loc,center,x,&
       x_sphere(ixImin1:ixImax1,r_))

    select case(trim(typeaxial_loc))
    case('slab')
      FORALL (idims=1:ndim)XX(ixOmin1:ixOmax1,idims) = x(ixOmin1:ixOmax1,&
         idims)-center(idims)
      if(z_in)then
       where(x_sphere(ixOmin1:ixOmax1,r_)>smalldouble)
         x_sphere(ixOmin1:ixOmax1,theta_) = dacos((x(ixOmin1:ixOmax1,&
            z_)-center(z_))/x_sphere(ixOmin1:ixOmax1,r_))
       elsewhere
         x_sphere(ixOmin1:ixOmax1,theta_) = 0.0_dp
       end where
      end if
    case('cylindrical')
      if(phi_<=ndim)then
       where(x_sphere(ixOmin1:ixOmax1,r_)>smalldouble)
        x_sphere(ixOmin1:ixOmax1,phi_) = dacos((x(ixOmin1:ixOmax1,&
           x_)-center(x_))/x_sphere(ixOmin1:ixOmax1,r_))
       elsewhere
         x_sphere(ixOmin1:ixOmax1,phi_) = 0.0_dp
       end where
      end if

      if(z_in) then
       where(x_sphere(ixOmin1:ixOmax1,r_)>smalldouble)
        x_sphere(ixOmin1:ixOmax1,theta_) = dacos((x(ixOmin1:ixOmax1,&
           z_)-center(z_))/x_sphere(ixOmin1:ixOmax1,r_))
       elsewhere
         x_sphere(ixOmin1:ixOmax1,theta_) = 0.0_dp
       end where
      end if


    case('spherical')
      if(ndim == 1)    then
       ! dummy
      else if(ndim==2) then
        if(phi_in)then
          x_sphere(ixOmin1:ixOmax1,phi_)   = x(ixOmin1:ixOmax1,phi_)
        else if(z_in) then
          x_sphere(ixOmin1:ixOmax1,theta_)   = x(ixOmin1:ixOmax1,theta_)
        end if
      else if(ndim==3) then
       x_sphere(ixOmin1:ixOmax1,theta_)   = x(ixOmin1:ixOmax1,theta_)
       x_sphere(ixOmin1:ixOmax1,phi_)     = x(ixOmin1:ixOmax1,phi_)
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
   end select
  end subroutine usr_get_volume
  !----------------------------------------------------------------------
!------------------------------------------------------------------------------------------
  !> subroutine to compute large variation in of an array
  subroutine usr_meanvalue_of_array(ixImin1,ixImax1,ixOmin1,ixOmax1,array,&
     mean_array,mean_flag)
    use mod_global_parameters
    use mod_small_values
    implicit none
    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    real(kind=dp), intent(in)       :: array(ixImin1:ixImax1)
    real(kind=dp), intent(out)      :: mean_array(ixImin1:ixImax1)
    logical, intent(out)            :: mean_flag(ixImin1:ixImax1)
    ! .. local ..
    integer                         :: ix1,kxOmin1,kxOmax1,n_cell_daverage,&
       n_cell_true
    real(kind=dp)                   :: sum_array,usr_max_deviation
    !----------------------------------------------------------------
    usr_max_deviation = 100.0d0
    n_cell_daverage=2
    mean_flag(ixImin1:ixImax1)  = .true.

    do ix1= ixOmin1,ixOmax1

      kxOmin1= max(ix1-n_cell_daverage, ixImin1);
      kxOmax1= min(ix1+n_cell_daverage, ixImax1);

      SUM_array=sum(array(kxOmin1:kxOmax1),mean_flag(kxOmin1:kxOmax1))

      n_cell_true = count(mean_flag(kxOmin1:kxOmax1))

      if(n_cell_true>1) then
        mean_array(ix1)=(sum_array-array(ix1))/(count(mean_flag(&
           kxOmin1:kxOmax1))-1.0_dp)
      else
        mean_array(ix1)= array(ix1)
      end if

      if(dabs(array(ix1))>smalldouble.and.dabs(mean_array(&
         ix1))>smalldouble)then

          mean_flag(ix1)=dabs((array(ix1)-&
             mean_array(ix1))/mean_array(ix1))<usr_max_deviation
      else
          mean_flag(ix1)=.true.
      end if

    end do

  end subroutine usr_meanvalue_of_array
  !------------------------------------------------------------------------------------------


  !------------------------------------------------------------------------------------------
  !> subroutine to compute large variation in of an array
  subroutine usr_medianvalue_of_array(ixImin1,ixImax1,ixOmin1,ixOmax1,array,&
     median_array,indices_median,mean_flag)
    use mod_global_parameters
    use mod_small_values
    implicit none
    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    real(kind=dp), intent(in)       :: array(ixImin1:ixImax1)
    real(kind=dp), intent(out)      :: median_array(ixImin1:ixImax1)
    logical, intent(out)            :: mean_flag(ixImin1:ixImax1)
    integer, intent(out)            :: indices_median(ixImin1:ixImax1,1:ndim,&
       2)
    ! .. local ..
    integer                         :: ix1,kxOmin1,kxOmax1,n_cell_daverage,&
       n_cell_true
    integer                         :: loc_indices(1:3,1:2)
    real(kind=dp)                   :: median_value,usr_max_deviation
    real(kind=dp)                   :: usr_min_deviation_abs,&
       usr_max_variation_abs
    logical                         :: logic_sort
    !----------------------------------------------------------------
    usr_max_deviation     = 2.0d0
    usr_min_deviation_abs = 1.0e-4
    usr_max_variation_abs = 2.0d-2
    n_cell_daverage       = 1
    mean_flag(ixImin1:ixImax1)      = .true.




    do ix1= ixOmin1,ixOmax1

      kxOmin1= max(ix1-n_cell_daverage, ixImin1);
      kxOmax1= min(ix1+n_cell_daverage, ixImax1);
      if(dabs(array(ix1))<usr_min_deviation_abs)then
        median_array(ix1)=array(ix1)
        mean_flag(ix1)   =.true.
        cycle
      end if
      if(maxval(array(kxOmin1:kxOmax1))-minval(array(&
         kxOmin1:kxOmax1))<usr_max_variation_abs)then
        median_array(ix1)=array(ix1)
        mean_flag(ix1)   =.true.
        cycle
      end if

      call usr_median_from_w(ixImin1,ixImax1,kxOmin1,kxOmax1,array,&
         median_value,loc_indices,logic_sort)

      if(logic_sort)then
        indices_median(ix1,1:ndim,1:2) =  loc_indices(1:ndim,1:2)
         median_array(ix1) = median_value

       if(dabs(array(ix1))>dabs(median_array(ix1)).and.&
          dabs(median_value)>smalldouble)then

          mean_flag(ix1)=dabs((array(ix1)-&
             median_value)/median_value)<usr_max_deviation

       else if(dabs(median_value)<smalldouble) then
          mean_flag(ix1)=dabs(array(ix1))<100.0_dp
       else
          mean_flag(ix1)=.true.
       end if
      else
        median_array(ix1) =array(ix1)
      end if
    end do

end subroutine usr_medianvalue_of_array
  !------------------------------------------------------------------------------------------
  subroutine usr_median_from_w(ixImin1,ixImax1,ixOmin1,ixOmax1,w_array,&
     median_value,indice_median,logic_sort)
    use mod_global_parameters
    use mod_small_values
    implicit none
    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    real(kind=dp), intent(in)       :: w_array(ixImin1:ixImax1)
    real(kind=dp), intent(out)      :: median_value
    integer, intent(out)            :: indice_median(3,2)
    logical, intent(out)            :: logic_sort
    ! .. local ..
    integer                         :: n_point, n_half,n_R,n_colonne
    integer                         :: n_cells(ndim),ix_half(ndim)
    real(kind=dp)                   :: array((ixOmax1-ixOmin1+1))
    integer                         :: array_ind((ixOmax1-ixOmin1+1))
    !----------------------------------------------------
     n_point=(ixOmax1-ixOmin1+1)
     array=reshape(w_array(ixOmin1:ixOmax1),(/n_point/))
     call usr_median_array(n_point,array,array_ind,median_value,logic_sort)
     if(.not.logic_sort)return

     n_cells(1)=ixOmax1-ixOmin1+1;
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
  subroutine usr_lorentz_transrmation_add_proper_speed(ixImin1,ixImax1,ixOmin1,&
     ixOmax1,velocity_proper,v)
    ! Eqaution in use https://en.wikipedia.org/wiki/Lorentz_transformation in Transformation of velocities
    implicit none
    integer, intent(in)           :: ixImin1,ixImax1,ixOmin1,ixOmax1
    real(kind=dp), intent(in)     :: velocity_proper(1:ndir)
    real(kind=dp), intent(inout)  :: v(ixImin1:ixImax1,1:ndir)
    !.. local ..
    integer                    :: idir
    real(kind=dp)              :: uv(ixOmin1:ixOmax1)
    real(kind=dp)              :: lfac_proper
    !-------------------------------------------------------------
    if(all(velocity_proper==0.0_dp))return
    uv(ixOmin1:ixOmax1)  = v(ixImin1:ixImax1,1)*velocity_proper(1)
    cond_idir2: if(idir>1) then
     Loop_idir : do idir = 2,ndir
      uv(ixOmin1:ixOmax1) = uv(ixOmin1:ixOmax1) + v(ixImin1:ixImax1,&
         idir)*velocity_proper(idir)
     end do Loop_idir
    end if cond_idir2
    lfac_proper   = 1.0_dp/dsqrt(1.0_dp - SUM(velocity_proper**2.0_dp))

    Loop_idir_newv : do idir = 2,ndir
     v(ixOmin1:ixOmax1,idir) = 1.0_dp/(1.0_dp+&
        uv(ixOmin1:ixOmax1))*(v(ixOmin1:ixOmax1,&
        idir)/lfac_proper+velocity_proper(idir)+lfac_proper/(lfac_proper+&
        1)*uv(ixOmin1:ixOmax1)*velocity_proper(idir) )
    end do Loop_idir_newv
  end subroutine usr_lorentz_transrmation_add_proper_speed
end module mod_obj_mat
