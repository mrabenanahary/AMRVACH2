module mod_usr
  use mod_mhd
  implicit none
  real(dp) :: theta, kx, ly, vc

  ! Star features
  real(dp)  :: star_radius     ! star radius (cm)
  real(dp)  :: star_mass       ! star mass  (g)
  real(dp)  :: star_luminosity ! star luminosity in erg/s 
  real(dp)  :: star_Eddington  ! Eddington limit erg/s
  real(dp)  :: star_temperature! star temperature
  real(dp)  :: star_vrotation  ! rotation parameter cm/s
  real(dp)  :: star_magnetic   ! Magnetic field strength at star surface (gauss)
  real(dp)  :: star_eta

  ! cak constants:
  real(dp)  :: cak_Q ! energy exchange radiation-thermal energy
  real(dp)  :: cak_alpha ! power low distribution  of clumbs
  real(dp)  :: cak_Qbar
  real(dp)  :: cak_kappa

  
  ! wind
  real(dp)  :: wind_massflux 
  real(dp)  :: wind_surface_star_density ! density
  real(dp)  :: wind_density_factor
  real(dp)  :: wind_sdspeed ! sound speed


  ! local normalisation
  real(dp)  :: unit_inv_distance
  real(dp)  :: unit_inv_time
  real(dp)  :: unit_inv_speed
  real(dp)  :: unit_inv_volum
  real(dp)  :: unit_inv_mass
  real(dp)  :: unit_inv_energy 
  ! scaling
  real(dp)  :: unit_luminosity


  logical   :: use1Dfile
  ! physical constantes:
  real(dp),save  :: const_kappae = 0.345D0
  real(dp),save  :: const_newt   = 6.67259D-8
  real(dp),save  :: const_clight = 2.99792458d10   ! cm s^-1
  ! geometry
  character(len=30):: coordinate_system

  ! reconstruction from old 1D file
  logical               :: reuse_old1Dfile
  integer               :: type_oldfile
  integer               :: old_nw,old_ndim,old_npoint
  integer ,allocatable  :: old_iw(:)
  real(dp),allocatable  :: old_w(:,:)
  real(dp),allocatable  :: old_x(:,:)
  character(len=30)     :: old_file_filenameini
contains
  subroutine usr_init()
    call set_coordinate_system("spherical_3D")

    usr_set_parameters  => initglobaldata_usr
    usr_init_one_grid   => initonegrid_usr
    usr_special_bc      => specialbound_usr
    usr_aux_output      => specialvar_output
    usr_add_aux_names   => specialvarnames_output 
    usr_source          => specialsource
    usr_get_dt          => getdt_special

    call usr_params_read(par_files)
    call set_coordinate_system(trim(coordinate_system))
    call mhd_activate()
    call usr_physical_unit()
  end subroutine usr_init

  !> Read this module s parameters from a file
  subroutine usr_params_read(files)
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /usr_list/ star_radius, &
                        star_mass, &
                        star_luminosity, &
                        star_vrotation,&
                        star_magnetic,&
                        wind_density_factor,&
                        cak_alpha, cak_Qbar,&
                        coordinate_system,&
                        reuse_old1Dfile,&
                        type_oldfile,&
                        old_file_filenameini


    if(mype==0)write(*,*)'Reading usr_list'
    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, usr_list, end=111)
111    close(unitpar)
    end do

  end subroutine usr_params_read
!--------------------------------------------------------------------
  subroutine usr_physical_unit()

   unit_inv_distance   = 1.0_dp/unit_length
   unit_inv_time       = 1.0_dp/unit_time
   unit_inv_speed      = 1.0_dp/unit_velocity 
   unit_inv_volum      = unit_inv_distance**3.0_dp
   unit_inv_mass       = unit_inv_volum /unit_density
   unit_inv_energy     = unit_inv_mass * unit_inv_speed**2.0_dp

   unit_luminosity     = unit_inv_distance**5.0_dp &
                         /(unit_inv_mass*unit_inv_time**3.0_dp*const_mp) 
  end subroutine usr_physical_unit


!------------------------------------------------------------------------
  subroutine usr_normalise_parameters()

   const_newt         = const_newt*&
                        (unit_density*(unit_length/unit_velocity)**2.0_dp)

   const_kappae       = const_kappae*(unit_length*unit_density) 

   const_clight       = const_clight/unit_velocity

   star_luminosity    = star_luminosity/(unit_density*unit_length**2.0_dp&
                                         *unit_velocity**3.0_dp)

   star_radius        = star_radius/unit_length

   star_mass          = star_mass*unit_inv_mass

   star_temperature   = star_temperature/unit_temperature

  
   star_eta           = (star_magnetic*star_radius)**2.0&
                        /(400.0_dp*wind_massflux*wind_sdspeed)

   wind_sdspeed   = wind_sdspeed*unit_inv_speed      

   star_vrotation = star_vrotation*unit_inv_speed  

   wind_surface_star_density = wind_surface_star_density/unit_density


  end subroutine usr_normalise_parameters


!-------------------------------------------------------------------------
  subroutine initglobaldata_usr()
   use mod_variables

   star_Eddington = const_kappae*star_luminosity&
                    /(4.0_dp*dpi*const_newt*star_mass*const_c)

   star_temperature=(star_luminosity&
                   /(0.0000567*4.0_dp*dpi*star_radius**2.0_dp))**0.25_dp

   wind_sdspeed  = dsqrt(mhd_gamma*star_temperature*(const_kB/const_mp))

   wind_massflux = star_luminosity/const_c**2.0&
     *cak_alpha/(1.0_dp-cak_alpha)&
     *(star_Eddington*cak_Qbar/(1.0-star_Eddington))**(1.0_dp/cak_alpha&
                                                       -1.0_dp)

   wind_surface_star_density = wind_massflux*wind_density_factor/&
                               (4.0_dp*dpi*star_radius**2.0_dp*wind_sdspeed)

   call usr_normalise_parameters()
   !-----------------------------
   if(reuse_old1Dfile)call read_oldfile_oneblock(22)
   !-----------------------------
  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    ! initialize one grid
    integer, intent(in)     :: ixI^L,ixO^L
    real(dp), intent(in)    :: x(ixI^S,1:ndim)
    real(dp), intent(inout) :: w(ixI^S,1:nw)

    real(dp) :: res
    integer :: ix^D,na
    logical, save :: first=.true.

    if(first)then
      if(mype==0) then
        write(*,*)'Simulating 3D rotating dipole'
      endif
      first=.false.
    endif
    call usr_set_stellarwind(ixI^L,ixO^L,x,w)
    call mhd_to_conserved(ixI^L,ixO^L,w,x)
  end subroutine initonegrid_usr
!----------------------------------------------------------------
  subroutine specialbound_usr(qt,ixI^L,ixO^L,iB,w,x)
    ! special boundary types, user defined
    integer, intent(in) :: ixO^L, iB, ixI^L
    real(dp), intent(in) :: qt, x(ixI^S,1:ndim)
    real(dp), intent(inout) :: w(ixI^S,1:nw)

    real(dp), dimension(ixI^S) :: pth,tmp,ggrid,invT,dv
    real(dp) :: delydelx
    integer :: ix^D,idir,ix1plus1,ix1plus2

    select case(iB)
    case(1) 
      ix1plus1=ixOmax1+1
      ix1plus2=ix1plus1+1        
      w(ixO^S,rho_)   = wind_surface_star_density 
      w(ixO^S,p_)     = w(ixO^S,rho_)*star_temperature
      Loop_ix1 : do ix1=ixOmin1,ixOmax1
       dv(ix1,ixO^SE) = SUM(block%dx(ix1:ixOmax1,ixO^SE,1),dim=1)&
                      /(block%dx(ix1plus1,ixO^SE,1)*(ixOmax1-ix1+1))&
                     *(w(ix1plus1,ixO^SE,mom(1))/w(ix1plus1,ixO^SE,rho_)&
                      -w(ix1plus2,ixO^SE,mom(1))/w(ix1plus2,ixO^SE,rho_))

       {^IFONED
       if(dv(ix1)>0)then
        w(ix1,ixO^SE,mom(1)) = min( w(ix1plus1,ixO^SE,mom(1))&
                          /w(ix1plus1,ixO^SE,rho_)+dv(ix1)& 
                          ,w(ix1plus1,ixO^SE,mom(1))/w(ix1plus1,ixO^SE,rho_))
       else
        w(ix1,ixO^SE,mom(1)) = sign(min(dabs(w(ix1plus1,ixO^SE,mom(1))&
                          /w(ix1plus1,ixO^SE,rho_)+dv(ix1))&
                          ,wind_sdspeed)&
                          ,w(ix1plus1,ixO^SE,mom(1))/w(ix1plus1,ixO^SE,rho_))
       end if
       }
       {^NOONED
       where(dv(ix1,ixO^SE)>0.0_dp)
        w(ix1,ixO^SE,mom(1)) = min( w(ix1plus1,ixO^SE,mom(1))&
            /w(ix1plus1,ixO^SE,rho_)&
            +dv(ix1,ixO^SE),w(ix1plus1,ixO^SE,mom(1))/w(ix1plus1,ixO^SE,rho_))

       else where
        w(ix1,ixO^SE,mom(1)) = sign(min(dabs(w(ix1plus1,ixO^SE,mom(1))&
              /w(ix1plus1,ixO^SE,rho_)+dv(ix1,ixO^SE)),wind_sdspeed)&
              ,w(ix1plus1,ixO^SE,mom(1))/w(ix1plus1,ixO^SE,rho_))

       end where
       }

       {^IFONED
       if(w(ix1,mom(1))<0.0_dp.and.&
             w(ix1,mom(1))>-wind_sdspeed)&
        w(ix1,mom(1))=w(ix1plus1,mom(1))/w(ix1plus1,rho_)&
               *(x(ixOmax1+1,1)/x(ix1,1))**2.0_dp
       if(w(ix1,mom(1))>wind_sdspeed)w(ix1,ixO^SE,mom(1)) = wind_sdspeed
       if(w(ix1,mom(1))<-wind_sdspeed)w(ix1,mom(1))=-wind_sdspeed
       if(dabs(w(ix1,mom(1))) > dabs(w(ix1plus1,mom(1))))&
        w(ix1,mom(1)) = w(ix1plus1,mom(1))
       if(w(ix1,mom(1))*w(ix1plus1,mom(1))<0.0_dp)w(ix1,mom(1)) = 0.0_DP
       }

       {^NOONED
       where(w(ix1,ixO^SE,mom(1))<0.0_dp.and.&
             w(ix1,ixO^SE,mom(1))>-wind_sdspeed)
        w(ix1,ixO^SE,mom(1))=w(ix1plus1,ixO^SE,mom(1))/w(ix1plus1,ixO^SE,rho_)&
               *(x(ixOmax1+1,ixO^SE,1)/x(ix1,ixO^SE,1))**2.0_dp
       end where
       where(w(ix1,ixO^SE,mom(1))>wind_sdspeed)
         w(ix1,ixO^SE,mom(1)) = wind_sdspeed
       end where
       where(w(ix1,ixO^SE,mom(1))<-wind_sdspeed)
         w(ix1,ixO^SE,mom(1)) = -wind_sdspeed
       end where
       where(dabs(w(ix1,ixO^SE,mom(1))) > dabs(w(ix1plus1,ixO^SE,mom(1))))
         w(ix1,ixO^SE,mom(1)) = w(ix1plus1,ixO^SE,mom(1))
       endwhere
       where(w(ix1,ixO^SE,mom(1))*w(ix1plus1,ixO^SE,mom(1))<0.0_dp)
         w(ix1,ixO^SE,mom(1)) = 0.0_DP
       endwhere
       }
      end do Loop_ix1

      if(star_eta<0) then
        w(ixO^S,mom(2)) = 0.0_DP
      else
       if(ndir>1)then
        Loop_ix1_mom2: do ix1=ixOmax1,ixOmin1,-1
         tmp(1,ixO^SE ) = sum(block%dx(ix1:ixOmax1,ixO^SE,1),dim=1)&
                            /( block%dx(ix1plus1,ixO^SE,1)*(ixOmax1-ix1+1))& 
                        * (w(ix1plus1,ixO^SE,mom(2))/w(ix1plus1,ixO^SE,rho_)&
                      -w(ix1plus2,ixO^SE,mom(2))/w(ix1plus2,ixO^SE,rho_))
         w(ix1,ixO^SE,mom(2))=w(ix1plus1,ixO^SE,mom(2))/w(ix1plus1,ixO^SE,rho_)&
                              +tmp(1,ixO^SE)
         {^IFONED
         if(w(ix1,mom(2))*w(ix1plus1,mom(2))>0._dp)w(ix1,mom(2)) = 0.0_dp
         }
         {^NOONED
         where(w(ix1,ixO^SE,mom(2))*w(ix1plus1,ixO^SE,mom(2))>0)
          w(ix1,ixO^SE,mom(2)) = 0.0_dp
         end where
         }
        end do Loop_ix1_mom2
       end if
      end if
      if(ndir>2)w(ixO^S,mom(3))=star_vrotation*x(ixO^S,1)/star_radius{^NOONED &
                      * dsin(x(ixO^S,2))}
      w(ixO^S,p_) = w(ixO^S,rho_)*star_temperature
      w(ixO^S,mag(1))=0.0_dp
   
      cond_no1dirB2 : if(ndir>1)then
       Loop_b2: do ix1=ixOmax1,ixOmin1,-1
         tmp(1,ixO^SE) = sum(block%dx(ix1:ixOmax1,ixO^SE,1),dim=1)&
                            / block%dx(ix1plus1,ixO^SE,1)&
                        * (w(ix1plus1,ixO^SE,mag(2))&
                           -w(ix1plus2,ixO^SE,mag(2)))

         w(ix1,ixO^SE,mag(2)) = w(ix1plus1,ixO^SE,mag(2))+tmp(1,ixO^SE)
         {^IFONED
         if(w(ix1,mag(2)) * w(ix1plus1,mag(2))<0.0_dp)&
          w(ix1,mag(2)) = 0.0_dp
         }
         {^NOONED
         where(w(ix1,ixO^SE,mag(2)) * w(ix1plus1,ixO^SE,mag(2))<0.0_dp)
          w(ix1,ixO^SE,mag(2)) = 0.0_dp
         end where
         }
       end do Loop_b2
      end if cond_no1dirB2
      cond_3dirB3 : if(ndir>2)then
       Loop_b3: do ix1=ixOmax1,ixOmin1,-1
         tmp(1,ixO^SE) = sum(block%dx(ix1:ixOmax1,ixO^SE,1),dim=1)&
                            / block%dx(ix1plus1,ixO^SE,1)&
                         * (w(ix1plus1,ixO^SE,mag(3))&
                            -w(ix1plus2,ixO^SE,mag(3)))

         w(ix1,ixO^SE,mag(3)) = w(ix1plus1,ixO^SE,mag(3))+tmp(1,ixO^SE)
         {^IFONED
         if(w(ix1,mag(3)) * w(ix1plus1,mag(3))<0.0_dp)&
           w(ix1,mag(3)) = 0.0_dp
         }
         {^NOONED
         where(w(ix1,ixO^SE,mag(3)) * w(ix1plus1,ixO^SE,mag(3))<0.0_dp)
          w(ix1,ixO^SE,mag(3)) = 0.0_dp
         end where
         }
       end do Loop_b3
      end if cond_3dirB3
      call mhd_to_conserved(ixI^L,ixO^L,w,x)
    case default
       write(*,*) ' The boundary IB = ',IB, '  is not set at user side'
       call mpistop("Special boundary is not defined for this region")
    end select
    
  end subroutine specialbound_usr

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
  ! this subroutine can be used in convert, to add auxiliary variables to the
  ! converted output file, for further analysis using tecplot, paraview, ....
  ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
  ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
  ! corresponding normalization values (default value 1)
    integer, intent(in)                :: ixI^L,ixO^L
    real(dp), intent(in)       :: x(ixI^S,1:ndim)
    real(dp)                   :: w(ixI^S,nw+nwauxio)
    real(dp)                   :: normconv(0:nw+nwauxio)

    real(dp) :: pth(ixI^S),B2(ixI^S),divb(ixI^S)
    real(dp) :: Btotal(ixI^S,1:ndir)
    integer :: idir

    ! output temperature
    call mhd_get_pthermal(w,x,ixI^L,ixO^L,pth)
    w(ixO^S,nw+1)=pth(ixO^S)/w(ixO^S,rho_)

    do idir=1,ndir
      if(B0field) then
        Btotal(ixI^S,idir)=w(ixI^S,mag(idir))+block%B0(ixI^S,idir,0)
      else
        Btotal(ixI^S,idir)=w(ixI^S,mag(idir))
      endif
    end do
    ! B^2
    B2(ixO^S)=sum((Btotal(ixO^S,:))**2,dim=ndim+1)

    ! output Alfven wave speed B/sqrt(rho)
    w(ixO^S,nw+2)=dsqrt(B2(ixO^S)/w(ixO^S,rho_))

    ! output divB1
    call divvector(Btotal,ixI^L,ixO^L,divb)
    w(ixO^S,nw+3)=0.5_dp*divb(ixO^S)/dsqrt(B2(ixO^S))/(^D&1.0_dp/dxlevel(^D)+)
    ! output the plasma beta p*2/B**2
    w(ixO^S,nw+4)=pth(ixO^S)*two/B2(ixO^S)

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
  ! newly added variables need to be concatenated with the w_names/primnames string
    character(len=*) :: varnames
    varnames='Te Alfv divB beta'

  end subroutine specialvarnames_output


!-----------------------------------------------------------------------
  subroutine usr_set_stellarwind(ixI^L,ixO^L,x,w)
   implicit none
   integer, intent(in)       :: ixI^L,ixO^L
   real(dp), intent(in)      :: x(ixI^S,1:ndim)
   real(dp), intent(inout)   :: w(ixI^S,1:nw)
   
   logical                   :: interpw(1:nw)
   logical, dimension(ixG^T) :: patchw

    is_reuse1D : if(.not.use1Dfile) then
     w(ixO^S,mom(1)) = 2.4d8*(1.0_DP-star_radius/x(ixO^S,1))/unit_velocity
     w(ixO^S, rho_)  = wind_surface_star_density&
                     *(wind_sdspeed/max(w(ixO^S,mom(1)),smalldouble)&
                     *(star_radius/x(ixO^S,1))**2.0_dp)

     w(ixO^S,mag(:)) = 0.0
     if(ndir>1)then
      w(ixO^S,mom(2))=0.0_dp
      if(ndir==3)then
       w(ixO^S,mom(3))=0.0_dp
      end if
     end if

    else is_reuse1D
     interpw = .true.
     patchw  = .true.
     call interpolation_olddata(ixI^L,ixO^L,.false.,interpw,'other',x,w,patchw)     
    end if is_reuse1D
    w(ixO^S,p_)     = w(ixO^S,rho_)*star_temperature 

    if(mhd_glm)w(ixO^S,psi_)=0.0_dp
  end subroutine usr_set_stellarwind





!#############################################################################
! module amrvacusr - cakdriving
!=============================================================================
! Project : line driven winds
! pre-compile: include this routine in the amrvacusr.t.****** file
! based on the version  : 15/07/2009, Allard Jan
! state  : testing
! parameters :
! main variables
!   rho          : density
!   m^C=rho v^C  : momentum
!   e            : total energy
!============================================================================

subroutine cak_pointgrav(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

!
!  Calculate velocity change due to point source gravity
!  and CAK line driving force.
!  As it stands now, the points source has to be in the same
!  'dimensionality' as the grid. E.g. It is not possible to
!  have an xy grid, with the point source localized in the
!  z-direction.
!  On the other hand, it IS possible to have the point source
!  in a location that is not part of the grid. This will in
!  fact often be the case in spherical symmetries with the
!  gravitational source at r=0
!
!  the scaled value for the central mass: eqpar(Mstar_) has to be set in
!  the parfile, together with the coordinates x1ptms,x2ptms,x3ptms
!

use mod_physics


integer, intent(in)             :: ixI^L, ixO^L, iw^LIM
double precision, intent(in)    :: qdt, qtC, qt, x(ixI^S,1:ndim)
double precision, intent(in)    :: wCT(ixI^S,1:nw)
double precision, intent(inout) :: w(ixI^S,1:nw)

double precision                :: gforce(ixG^T)


!-----------------------------------------------------------------------------
if (star_mass == zero) return

call cak_accel(ixI^L,ixO^L,wCT,x,gforce)
!
!  update momentum and energy
!

w(ixO^S,mom(1)) =  w(ixO^S,mom(1)) + qdt*gforce(ixO^S)*wCT(ixO^S,rho_)
w(ixO^S,e_)  =  w(ixO^S,e_)  + qdt*gforce(ixO^S)*wCT(ixO^S,mom(1))
call phys_handle_small_values(.false.,w,x,ixI^L,ixO^L,"pointgravcak")
if(it<1) PRINT*,' the accelerations ',gforce(ixO^S)
end subroutine cak_pointgrav

!==============================================================================
subroutine cak_pointgrav_getdt(w,ixG^L,ixO^L,dtnew,dx^D,x)

!
! Limits timestep for gravitational pointsource cak force
!

integer, intent(in)             :: ixG^L, ixO^L
double precision, intent(in)    :: dx^D, x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw), dtnew
double precision                :: gforce(ixG^T)

! .. local ..
double precision                :: r2i(ixG^T,1:ndim)
double precision                :: dtgravcak
integer                         :: idims
!-----------------------------------------------------------------------------

if( star_Mass == zero ) return
!
! call cak_accel with w instead of wCT, because you need values for the next time step
!
call cak_accel(ixG^LL,ixO^L,w,x,gforce)
dtnew = min(bigdouble,{^D&minval(sqrt(dabs(block%dx(ixO^S,^D)/max(gforce(ixO^S),smalldouble))))|,})

!dtnew = dtgravcak

end subroutine cak_pointgrav_getdt
!===========================================================================
subroutine cak_accel(ixI^L,ixO^L,wCT,x,gforce)

use mod_physics
integer, intent(in)             :: ixI^L, ixO^L
double precision, intent(in)    :: x(ixI^S,1:ndim)
double precision, intent(in)    :: wCT(ixI^S,1:nw)

double precision, intent(out)   :: gforce(ixG^T)

! ... local ...
double precision  :: cak_opa, cak_oma,x^Dptms
double precision,dimension(ixG^T)  :: v1, dvdr
double precision,dimension(ixG^T)  :: Tlocal,tt, gthin, sigma, &
                                      fdisk, gcak, kqbareff, ptherm
double precision  :: r2i(ixG^T,1:ndim)
integer           :: idims
!-----------------------------------------------------------------------------
x^Dptms=0.0_dp;
call getgvector(ixI^L,ixO^L,x,x^Dptms,r2i)
idims=1

call phys_get_v_idim(wCT,x,ixI^L,ixI^L,idims,v1)

select case(typegrad)
case("central")
     call gradient(v1,ixI^L,ixO^L,idims,dvdr)
case("limited")
     call gradientS(v1,ixI^L,ixO^L,idims,dvdr)
end select

dvdr(ixO^S) = max(dvdr(ixO^S), smalldouble)


call phys_get_pthermal(wCT,x,ixI^L,ixO^L,ptherm)
cak_opa=1.0_dp + cak_alpha
cak_oma=1.0_dp - cak_alpha

Tlocal(ixO^S)   = ptherm(ixO^S)/wCT(ixO^S,rho_)
tt(ixO^S)= Tlocal(ixO^S)/star_temperature
! exp(-100)==0 here for numerical stability, compiler dependent, fairly arbitrary
! cutoff.
where(tt(ixO^S)<100.0)
 kqbareff(ixO^S)=merge(dexp(-tt(ixO^S)+1.0_dp),1.0_dp,tt(ixO^S)>1.0_dp)
else where
 kqbareff(ixO^S)=zero
end where
! Eq. 12 Asif
gthin(ixO^S) = (cak_Qbar*const_kappae)**cak_oma &
              *star_luminosity*r2i(ixO^S,1) &
              /(4.0_dp*dpi*cak_oma*const_c**cak_opa)

sigma(ixO^S) = (1.0_dp-v1(ixO^S)*dsqrt(r2i(ixO^S,1))/dvdr(ixO^S)) &
               * (star_radius**2.0_dp) *r2i(ixO^S,1)

!
!  Finite disk correction factor
!


where(sigma(ixO^S) >= 1.0_dp)
      fdisk(ixO^S)= 1.0_dp/cak_opa
else where (sigma(ixO^S) < -1.0d10)
      fdisk(ixO^S) = ((-sigma(ixO^S))**cak_alpha/cak_opa)
else where (dabs(sigma(ixO^S)) > 1.0d-3)
      fdisk(ixO^S) = (1.0_dp-(1.0_dp-sigma(ixO^S))**cak_opa)&
                     /(sigma(ixO^S)*cak_opa)
else where
      fdisk(ixO^S) = 1.0_dp-0.5_dp*cak_alpha*sigma(ixO^S)   &
            * (1.0_dp+(1.0_dp-cak_alpha)*sigma(ixO^S)/3.0_dp)
end where

!
!  Calculates cak and gravity at the same time.
!

gforce(ixO^S) = gthin(ixO^S)*fdisk(ixO^S)&
                *(dvdr(ixO^S)/wCT(ixO^S,rho_))**cak_alpha

gforce(ixO^S) = gforce(ixO^S)*kqbareff(ixO^S) &
                + (const_newt*star_mass*r2i(ixO^S,1))*(star_Eddington-1.0_dp)

if(it<1)PRINT*,' is at cak accel',gthin(ixO^S)*fdisk(ixO^S)&
                *(dvdr(ixO^S)/wCT(ixO^S,rho_))**cak_alpha*kqbareff(ixO^S),' is you',(const_newt*star_mass*r2i(ixO^S,1))*(star_Eddington-1.0_dp),'r2iii',r2i(ixO^S,1), 'grad vitesse',dvdr(ixO^S), 'eeee',star_Eddington,'cons',const_newt,' mass',star_mass,'unit density',unit_density,'unit_volum',unit_inv_volum,'unit_mass',unit_inv_mass,'unit_length',unit_length

return
end subroutine cak_accel

!============================================================================

!===========================================================================
subroutine getgvector(ixI^L,ixO^L,x,x^Dptms,r2i)

!
!  returns the vector components of the 1/r^2 vector:
!  dx1/r^3,
!  dx2/r^3,
!  dx3/r^3,
!  between a pre-determined pointsource and a point in the grid.
!



integer, intent(in)           :: ixO^L,ixI^L
double precision, intent(in)  :: x(ixI^S,1:ndim),x^Dptms
double precision, intent(out) :: r2i(ixG^T,1:ndim)

double precision              :: rstptms, rctptms

double precision, allocatable :: r2(:^D&)
{^NOONED
double precision, allocatable :: rst(:^D&), rct(:^D&)
}
!-----------------------------------------------------------------------------
allocate(r2(ixO^S))


select case (typeaxial)
case('slab')
{^IFONED
   r2(ixO^S) = (x(ixO^S,1)-x1ptms)**2
}
{^IFTWOD
   r2(ixO^S) = (x(ixO^S,1)-x1ptms)**2  &
             + (x(ixO^S,2)-x2ptms)**2
}
{^IFTHREED
   r2(ixO^S) = (x(ixO^S,1)-x1ptms)**2  &
             + (x(ixO^S,2)-x2ptms)**2  &
             + (x(ixO^S,3)-x3ptms)**2
}


   r2i(ixO^S,1) = (x(ixO^S,1)-x1ptms) &
                / (r2(ixO^S) * sqrt(r2(ixO^S)) + smalldouble)
{^NOONED
   r2i(ixO^S,2) = (x(ixO^S,2)-x2ptms) &
                / (r2(ixO^S) * sqrt(r2(ixO^S)) + smalldouble)
}
{^IFTHREED
   r2i(ixO^S,3) = (x(ixO^S,3)-x3ptms) &
                / (r2(ixO^S) * sqrt(r2(ixO^S)) + smalldouble)
}

case('cylindrical')
{^IFONED
   r2(ixO^S)    = (x(ixO^S,r_)-x1ptms)**2
   r2i(ixO^S,1) = one / (r2(ixO^S) + smalldouble)
}
{^IFTWOD
   if(^PHI == 2) then
      r2(ixO^S)    = (x(ixO^S,r_)-x1ptms)**2 &
                   + two*x(ixO^S,r_)*x1ptms  &
                   * (one -cos(x(ixO^S,2)-x2ptms))

      r2i(ixO^S,1) = (x(ixO^S,r_)-x1ptms*cos(x(ixO^S,2)-x2ptms)) &
                   / (r2(ixO^S) * sqrt(r2(ixO^S)) + smalldouble)

      r2i(ixO^S,2) = x1ptms * sin(x(ixO^S,2)-x2ptms)             &
                   / (r2(ixO^S) * sqrt(r2(ixO^S)) + smalldouble)
   else
      r2(ixO^S)    = (x(ixO^S,r_)-x1ptms)**2 &
                   + (x(ixO^S,2)-x2ptms)**2

      r2i(ixO^S,1) = (x(ixO^S,r_)-x1ptms)    &
                   / (r2(ixO^S) * sqrt(r2(ixO^S)) + smalldouble)
      r2i(ixO^S,2) = (x(ixO^S,2)-x2ptms)    &
                   / (r2(ixO^S) * sqrt(r2(ixO^S)) + smalldouble)
   endif
}
{^IFTHREED
   r2(ixO^S) = (x(ixO^S,r_)-x1ptms)**2  &
             + (x(ixO^S,z_)-x2ptms)**2  &
             + two*x(ixO^S,r_)*x1ptms   &
             * (one -cos(x(ixO^S,phi_)-x3ptms))

   r2i(ixO^S,1) = (x(ixO^S,r_)-x1ptms*cos(x(ixO^S,phi_)-x3ptms)) &
                / (r2(ixO^S) * sqrt(r2(ixO^S)) + smalldouble)
   r2i(ixO^S,2) = (x(ixO^S,z_)-x2ptms)                           &
                / (r2(ixO^S) * sqrt(r2(ixO^S)) + smalldouble)
   r2i(ixO^S,3) = x1ptms * sin(x(ixO^S,phi_)-x3ptms)             &
                / (r2(ixO^S) * sqrt(r2(ixO^S)) + smalldouble)
}

case('spherical')

{^IFONED
   r2(ixO^S)   = (x(ixO^S,r_)-x1ptms)**2
   r2i(ixO^S,1) = one / (r2(ixO^S) + smalldouble)
}

{^NOONED
   rstptms = x1ptms * sin(x2ptms)
   rctptms = x1ptms * cos(x2ptms)
}

{^IFTWOD
   allocate(rct(ixO^S))
   allocate(rst(ixO^S))

   rct(ixO^S) = x(ixO^S,r_)*cos(x(ixO^S,2))
   rst(ixO^S) = x(ixO^S,r_)*sin(x(ixO^S,2))

   r2(ixO^S) = (rct(ixO^S)-rctptms)**2  &
             + (rst(ixO^S)-rstptms)**2


   r2i(ixO^S,1) = ((rct(ixO^S)-rctptms) * cos(x(ixO^S,2))     &
                + (rst(ixO^S)-rstptms)  * sin(x(ixO^S,2)))    &
                / (r2(ixO^S) * sqrt(r2(ixO^S)) + smalldouble)
   r2i(ixO^S,2) = (-(rct(ixO^S)-rctptms) * sin(x(ixO^S,2))    &
                + (rst(ixO^S) - rstptms) * cos(x(ixO^S,2)))   &
                / (r2(ixO^S) * sqrt(r2(ixO^S)) + smalldouble)

   deallocate(rst)
   deallocate(rct)
}
{^IFTHREED
   allocate(rct(ixO^S))
   allocate(rst(ixO^S))

   rct(ixO^S) = x(ixO^S,r_)*cos(x(ixO^S,2))
   rst(ixO^S) = x(ixO^S,r_)*sin(x(ixO^S,2))

   r2(ixO^S) = (rct(ixO^S)-rctptms)**2    &
             + (rst(ixO^S)-rstptms)**2    &
             + two * rst(ixO^S) * rstptms &
             * (one-cos(x(ixO^S,phi_)-x3ptms))

   r2i(ixO^S,1) = ((rct(ixO^S)-rctptms) * cos(x(ixO^S,2))        &
                + (rst(ixO^S)-rstptms*cos(x(ixO^S,phi_)-x3ptms)) &
                * sin(x(ixO^S,2)))    &
                / (r2(ixO^S) * sqrt(r2(ixO^S)) + smalldouble)
   r2i(ixO^S,2) = (-(rct(ixO^S)-rctptms) * sin(x(ixO^S,2))       &
                + (rst(ixO^S)-rstptms*cos(x(ixO^S,phi_)-x3ptms)) &
                * cos(x(ixO^S,2)))  &
                / (r2(ixO^S) * sqrt(r2(ixO^S)) + smalldouble)
   r2i(ixO^S,3) = rstptms * sin(x(ixO^S,phi_)-x3ptms)            &
                / (r2(ixO^S) * sqrt(r2(ixO^S)) + smalldouble)


   deallocate(rct)
   deallocate(rst)
}
   end select

deallocate(r2)

return

end subroutine getgvector

!=============================================================================
subroutine specialsource(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

! Calculate w(iw)=w(iw)+qdt*SOURCE[wCT,qtC,x] within ixO for all indices
! iw=iwmin...iwmax.  wCT is at time qCT


integer, intent(in) :: ixI^L, ixO^L, iw^LIM
double precision, intent(in) :: qdt, qtC, qt, x(ixI^S,1:ndim)
double precision, intent(inout) :: wCT(ixI^S,1:nw), w(ixI^S,1:nw)

integer :: ix
!-----------------------------------------------------------------------------
!if(it<1)print*,'is your in ',w(ixO^S,:),saveigrid,it

call cak_pointgrav(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

!if(it<1)print*,'is your end ',w(ixO^S,:)
!if(it>3)call MPISTOP('your lasttest')
end subroutine specialsource

!=============================================================================
subroutine getdt_special(w,ixG^L,ix^L,dtnew,dx^D,x)

! Limit "dt" further if necessary, e.g. due to the special source terms.
! The getdt_courant (CFL condition) and the getdt subroutine in the AMRVACPHYS
! module have already been called.


integer, intent(in) :: ixG^L, ix^L
double precision, intent(in) :: dx^D, x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw), dtnew
!-----------------------------------------------------------------------------


call cak_pointgrav_getdt(w,ixG^L,ix^L,dtnew,dx^D,x)

end subroutine getdt_special



subroutine read_oldfile_oneblock(unitoldfile)
! type_oldfile > 0 hydrodynamics case
  use mod_physics
  integer, intent(in)  :: unitoldfile
! .. local variables ..
  character(len=1024) :: outfilehead
  integer, dimension(2) :: sizes, subsizes, start
  integer                    :: type_readoldX_block,type_readoldw_block
  integer                    :: ix1
!-----------------------------------------------------------------------------
useold : if (type_oldfile/=0) then
  allocate(old_iw(1:nw))
  old_iw(1:nw)=-1
  if (type_oldfile==11) then
     old_ndim= 1
     old_nw  = 3
     old_iw(rho_) = rho_;old_iw(mom(1)) = mom(1);old_iw(e_) = e_;
  end if

  readfile :if(mype==0) then
     open(unit=unitoldfile,file=old_file_filenameini)
     read(unitoldfile,"(a)")outfilehead
     read(unitoldfile,*)old_npoint
     allocate(old_x(1:old_npoint,1:ndim),old_w(1:old_npoint,1:nw))
     do ix1=1,old_npoint
       read(unitoldfile,*)old_x(ix1,1:old_ndim),old_w(ix1,1:old_nw)
     end do
     close(unitoldfile)
  end if readfile

multicpu :  if (npe>1) then
  sizes(1)=old_npoint;
  sizes(2)=ndim
  subsizes(1)=old_npoint;
  subsizes(2)=old_ndim
  start(1)=0;
  start(2)=0
  call MPI_TYPE_CREATE_SUBARRAY(2,sizes,subsizes,start, &
                                MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION, &
                                type_readoldX_block,ierrmpi)
  call MPI_TYPE_COMMIT(type_readoldX_block,ierrmpi)

  sizes(1)=old_npoint;
  sizes(2)=old_nw
  subsizes(1)=old_npoint;
  subsizes(2)=old_nw
  start(1)=0;
  start(2)=0
  call MPI_TYPE_CREATE_SUBARRAY(2,sizes,subsizes,start, &
                              MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION, &
                              type_readoldw_block,ierrmpi)
  call MPI_TYPE_COMMIT(type_readoldw_block,ierrmpi)

  call MPI_BCAST(old_npoint,1,MPI_INTEGER,0,icomm,ierrmpi)
  call MPI_BCAST(old_x,type_readoldw_block,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
  call MPI_BCAST(old_w,type_readoldX_block,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
 end if multicpu

end if useold
end subroutine read_oldfile_oneblock
!=============================================================================
subroutine interpolation_olddata(ixI^L,ix^L,conservar,interpw,typeinter,x,w,patchw)

  use mod_physics
  integer, intent(in)               :: ixI^L,ix^L
  logical, intent(in)               :: conservar
  logical, intent(in)               :: interpw(1:nw)
  character(len=*), intent(in)      :: typeinter
  double precision, intent(in)      :: x(ixI^S,1:ndim)
  double precision, intent(inout)   :: w(ixI^S,1:nw)
  logical, intent(inout)            :: patchw(ixG^T)


! .. local variables ..
  logical                           :: patchwint(ixG^T)
  integer                           :: rightnei(1),leftnei(1),ix^D,iw
  double precision                  :: R1,R2,R3,R4, normR
  double precision, allocatable     :: Dist(:)
  integer :: nmmi(1), nmpi(1), npmi(1), nppi(1) ,ni(1)
  double precision                  :: tmpw(ixG^T,1:nw)
!----------------------------------------------------------------------------
if (conservar) then
  patchwint(ix^S)=patchw(ix^S)
else
  patchwint(ix^S)=.true.
end if
if (abs(type_oldfile)<20 ) then
 Loop_1D : do ix1=ixmin1,ixmax1
   if (patchwint(ix1,ixmin^DE)) then
       rightnei=minloc(dabs(old_x(1:old_npoint,1)-x(ix1,ixmin^DE,1))&
                       ,old_x(1:old_npoint,1)>=x(ix1,ixmin^DE,1))
       leftnei=minloc(dabs(old_x(1:old_npoint,1)-x(ix1,ixmin^DE,1))&
                      ,old_x(1:old_npoint,1)<=x(ix1,ixmin^DE,1))
       Loop_iw_1D : do iw =1,nw
        if (interpw(iw).and.old_iw(iw)>0) then
         if(dabs(old_x(rightnei(1),1)-old_x(leftnei(1),1) )&
               >smalldouble) then
            w(ix1,ix^SE,iw)=(old_w(rightnei(1),iw)&
                            -old_w(leftnei(1),iw))&
                            /(old_x(rightnei(1),1)-old_x(leftnei(1),1))&
                            *(x(ix1,ix^SE,1)-old_x(leftnei(1),1))&
                            +old_w(leftnei(1),iw)
         else if(dabs(old_x(rightnei(1),1)-old_x(leftnei(1),1))&
                  <smalldouble) then
             w(ix1,ix^SE,iw)=old_w(leftnei(1),old_iw(iw))
          end if
        end if
       end do Loop_iw_1D
   end if
  end do Loop_1D
  if (any(.not.interpw(1:nw).or.old_iw(1:nw)<0)) then
    tmpw(ix^S,1:nw)=0
    do iw =1,nw
     if (interpw(iw).and.old_iw(iw)>0) tmpw(ix^S,iw)=w(ix^S,iw)
    end do
    call phys_to_primitive(ixI^L,ix^L,tmpw,x)
    do iw =1,nw
     if (interpw(iw).and.old_iw(iw)>0) w(ix^S,iw)=tmpw(ix^S,iw)
    end do
  end if
 end if
end subroutine interpolation_olddata
end module mod_usr
