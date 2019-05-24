!> Module with geometry-related routines (e.g., divergence, curl)
module mod_geometry

  use mod_constants
  implicit none
  public

contains

  !> Set the coordinate system to be used
  subroutine set_coordinate_system(geom)
    use mod_global_parameters
    character(len=*), intent(in) :: geom !< Name of the coordinate system

    select case (geom)
    case ("Cartesian","Cartesian_1D","Cartesian_2D","Cartesian_3D")
      ndir = ndim
      typeaxial='slab'
    case ("Cartesian_1.5D")
      if (ndim /= 1) call mpistop("Geometry Cartesian_1.5D but ndim /= 1")
      typeaxial='slab'
      ndir = 2
      x_=1; r_ = 1
      y_=2
    case ("Cartesian_1.75D")
      if (ndim /= 1) call mpistop("Geometry Cartesian_1.75D but ndim /= 1")
      typeaxial='slab'
      ndir = 3
      x_=1; r_      = 1
      y_=2; phi_    = 2
      z_=3; theta_  = 3
    case ("Cartesian_2.5D")
      if (ndim /= 2) call mpistop("Geometry Cartesian_2.5D but ndim /= 2")
      typeaxial='slab'
      ndir = 3
      x_=1; r_      = 1
      y_=2; phi_    = 2
      z_=3; theta_  = 3
    case ("cylindrical","cylindrical_2D","cylindrical_3D")
      ndir = ndim
      r_   = 1; x_     = 1
      z_   = 2; theta_ = 2
      if(ndir==3)then
         phi_ = 3; y_  = 3
      end if
      typeaxial='cylindrical'
    case ("cylindrical_2.5D")
      if (ndim /= 2) call mpistop("Geometry cylindrical_2.5D but ndim /= 2")
      ndir = 3
      r_   = 1; x_    = 1
      z_   = 2;theta_ = 2
      phi_ = 3; y_    = 3
      typeaxial='cylindrical'
    case ("polar","polar_2D","polar_3D")
      ndir = ndim
      r_   = 1;  x_  = 1
      phi_ = 2;  y_  = 2
      if(ndir==3)then
         z_     = 3
         theta_ = 3
      end if
      typeaxial='cylindrical'
    case ("polar_2.5D")
      if (ndim /= 2) call mpistop("Geometry polar_2.5D but ndim /= 2")
      ndir = 3
      r_   = 1; x_     = 1
      phi_ = 2; y_     = 2
      z_   = 3; theta_ = 3
      typeaxial='cylindrical'
    case ("spherical","spherical_2D","spherical_3D")
      ndir = ndim
      r_   = 1; x_     = 1
      if(ndir==3)then
         phi_ = 3; y_  = 3
      end if
      z_   = 3; theta_ = 2
      typeaxial='spherical'
    case ("spherical_2.5D")
      if (ndim /= 2) call mpistop("Geometry spherical_2.5D requires ndim == 2")

      ndir = 3
      r_   = 1; x_    = 1
      phi_ = 3; y_    = 3
      z_   = 2; theta_= 2
      typeaxial='spherical'
    case ("spherical_1.75D")
      if (ndim /= 1) call mpistop&
         ("Geometry spherical_1.75D requires ndim == 2")
      ndir = 3
      r_   = 1; x_     = 1
      phi_ = 3; y_     = 3
      z_   = 2; theta_ = 2
      typeaxial='spherical'
    case default
      call mpistop("Unknown geometry specified")
    end select
    r_in     = r_>0.and.r_<ndim
    phi_in   = phi_>0.and.phi_<=ndim
    theta_in = theta_>0.and.theta_<=ndim

    x_in     = x_>0.and.x_<=ndim
    y_in     = y_>0.and.y_<=ndim
    z_in     = z_>0.and.z_<=ndim
  end subroutine set_coordinate_system

  subroutine set_pole
    use mod_global_parameters

    select case (typeaxial)
    case ("spherical") 
    case ("cylindrical")
      
      if (1 == phi_ .and. periodB(1)) then
        if(mod(ng1(1),2)/=0) then
          call mpistop("Number of meshes in phi-direction should be even!")
        end if

        if(abs(xprobmin1)<smalldouble) then
          if (mype==0) then
            write(unitterm,*) "Will apply pi-periodic conditions at r=0"
          end if
          poleB(1,1)=.true.
        else
          if (mype==0) then
            write(unitterm,*) "There is no cylindrical axis!"
          end if
        end if
      end if
      
      if (2 == phi_ .and. periodB(2)) then
        if(mod(ng2(1),2)/=0) then
          call mpistop("Number of meshes in phi-direction should be even!")
        end if

        if(abs(xprobmin1)<smalldouble) then
          if (mype==0) then
            write(unitterm,*) "Will apply pi-periodic conditions at r=0"
          end if
          poleB(1,1)=.true.
        else
          if (mype==0) then
            write(unitterm,*) "There is no cylindrical axis!"
          end if
        end if
      end if
    end select

  end subroutine set_pole

  !> Deallocate geometry-related variables
  subroutine putgridgeo(igrid)
    use mod_global_parameters
    integer, intent(in) :: igrid

    deallocate(pw(igrid)%surfaceC,pw(igrid)%surface,pw(igrid)%dvolume,&
       pw(igrid)%dx,pw(igrid)%dxcoarse,pw(igrid)%ds,pw(igrid)%dvolumecoarse)

  end subroutine putgridgeo

  subroutine fillgeo(igrid,ixGmin1,ixGmin2,ixGmax1,ixGmax2)
    use mod_global_parameters

    integer, intent(in) :: igrid, ixGmin1,ixGmin2,ixGmax1,ixGmax2

    integer :: ixmin1,ixmin2,ixmax1,ixmax2, ixCmin1,ixCmin2,ixCmax1,ixCmax2
    double precision :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ndim),&
        drs(ixGmin1:ixGmax1,ixGmin2:ixGmax2), dx2(ixGmin1:ixGmax1,&
       ixGmin2:ixGmax2), dx3(ixGmin1:ixGmax1,ixGmin2:ixGmax2)
    !-----------------------------------------------------------------------------
    ixmin1=ixGmin1+1;ixmin2=ixGmin2+1;ixmax1=ixGmax1-1;ixmax2=ixGmax2-1;

    select case (typeaxial)
    case ("slabstretch")
      drs(ixGmin1:ixGmax1,ixGmin2:ixGmax2)=pw(igrid)%dx(ixGmin1:ixGmax1,&
         ixGmin2:ixGmax2,1)
      
      dx2(ixGmin1:ixGmax1,ixGmin2:ixGmax2)=pw(igrid)%dx(ixGmin1:ixGmax1,&
         ixGmin2:ixGmax2,2)
      

      
      
      ixCmin1=ixmin1-kr(1,1);ixCmin2=ixmin2-kr(2,1); ixCmax1=ixmax1
      ixCmax2=ixmax2;
      pw(igrid)%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         1)=dx2(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
      pw(igrid)%surface(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         1) =dx2(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
      ixCmin1=ixmin1-kr(1,2);ixCmin2=ixmin2-kr(2,2); ixCmax1=ixmax1
      ixCmax2=ixmax2;
      pw(igrid)%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         2)=drs(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
      pw(igrid)%surface(ixCmin1:ixCmax1,ixCmin2:ixCmax2,2)=drs(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)
     
      

    case ("spherical")
      x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1)=pw(igrid)%x(ixGmin1:ixGmax1,&
         ixGmin2:ixGmax2,1)
      
      x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,2)=pw(igrid)%x(ixGmin1:ixGmax1,&
         ixGmin2:ixGmax2,2)

      drs(ixGmin1:ixGmax1,ixGmin2:ixGmax2)=pw(igrid)%dx(ixGmin1:ixGmax1,&
         ixGmin2:ixGmax2,1)
      
      dx2(ixGmin1:ixGmax1,ixGmin2:ixGmax2)=pw(igrid)%dx(ixGmin1:ixGmax1,&
         ixGmin2:ixGmax2,2)
      

      ixCmin1=ixmin1-kr(1,1);ixCmin2=ixmin2-kr(2,1); ixCmax1=ixmax1
      ixCmax2=ixmax2;

      pw(igrid)%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,1)=(x(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,1)+half*drs(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2))**2  *two*dsin(x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         2))*dsin(half*dx2(ixCmin1:ixCmax1,ixCmin2:ixCmax2))

      
      ixCmin1=ixmin1-kr(1,2);ixCmin2=ixmin2-kr(2,2); ixCmax1=ixmax1
      ixCmax2=ixmax2;
      pw(igrid)%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,2)=x(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,1)*drs(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)*dsin(x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         2)+half*dx2(ixCmin1:ixCmax1,ixCmin2:ixCmax2))

      

      ixCmin1=ixmin1-kr(1,1);ixCmin2=ixmin2-kr(2,1); ixCmax1=ixmax1
      ixCmax2=ixmax2;
      pw(igrid)%surface(ixCmin1:ixCmax1,ixCmin2:ixCmax2,1)=x(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,1)**2  *two*dsin(x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         2))*dsin(half*dx2(ixCmin1:ixCmax1,ixCmin2:ixCmax2))
      
      ixCmin1=ixmin1-kr(1,2);ixCmin2=ixmin2-kr(2,2); ixCmax1=ixmax1
      ixCmax2=ixmax2;
      pw(igrid)%surface(ixCmin1:ixCmax1,ixCmin2:ixCmax2,2)=x(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,1)*drs(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)*dsin(x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,2))

      

    case ("cylindrical")
      x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1)=pw(igrid)%x(ixGmin1:ixGmax1,&
         ixGmin2:ixGmax2,1)
      drs(ixGmin1:ixGmax1,ixGmin2:ixGmax2)=pw(igrid)%dx(ixGmin1:ixGmax1,&
         ixGmin2:ixGmax2,1)
      
      dx2(ixGmin1:ixGmax1,ixGmin2:ixGmax2)=pw(igrid)%dx(ixGmin1:ixGmax1,&
         ixGmin2:ixGmax2,2)
      

      ixCmin1=ixmin1-kr(1,1);ixCmin2=ixmin2-kr(2,1); ixCmax1=ixmax1
      ixCmax2=ixmax2;
      pw(igrid)%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         1)=dabs(x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,1)+half*drs(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2))*dx2(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
      
      ixCmin1=ixmin1-kr(1,2);ixCmin2=ixmin2-kr(2,2); ixCmax1=ixmax1
      ixCmax2=ixmax2;
      if (z_==2) pw(igrid)%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         2)=x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,1)*drs(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)
      if (phi_==2) pw(igrid)%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         2)=drs(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
      

      ixCmin1=ixmin1-kr(1,1);ixCmin2=ixmin2-kr(2,1); ixCmax1=ixmax1
      ixCmax2=ixmax2;
      pw(igrid)%surface(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         1)=dabs(x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,1))*dx2(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)
      
      ixCmin1=ixmin1-kr(1,2);ixCmin2=ixmin2-kr(2,2); ixCmax1=ixmax1
      ixCmax2=ixmax2;
      if (z_==2) pw(igrid)%surface(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         2)=x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,1)*drs(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)
      if (phi_==2) pw(igrid)%surface(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         2)=drs(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
      

    case default
      call mpistop("Sorry, typeaxial unknown")
    end select

  end subroutine fillgeo

  !> Calculate gradient of a scalar q within ixL in direction idir
  subroutine gradient(q,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,idir,gradq)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2, idir
    double precision, intent(in)    :: q(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision, intent(inout) :: gradq(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision                :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    integer                         :: jxOmin1,jxOmin2,jxOmax1,jxOmax2,&
        hxOmin1,hxOmin2,hxOmax1,hxOmax2

    x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)=block%x(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:ndim)

    hxOmin1=ixOmin1-kr(idir,1);hxOmin2=ixOmin2-kr(idir,2)
    hxOmax1=ixOmax1-kr(idir,1);hxOmax2=ixOmax2-kr(idir,2);
    jxOmin1=ixOmin1+kr(idir,1);jxOmin2=ixOmin2+kr(idir,2)
    jxOmax1=ixOmax1+kr(idir,1);jxOmax2=ixOmax2+kr(idir,2);
    if(slab) then
      gradq(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=half*(q(jxOmin1:jxOmax1,&
         jxOmin2:jxOmax2)-q(hxOmin1:hxOmax1,hxOmin2:hxOmax2))/dxlevel(idir)
    else
      select case(typeaxial)
      case('slabstretch')
        gradq(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(q(jxOmin1:jxOmax1,&
           jxOmin2:jxOmax2)-q(hxOmin1:hxOmax1,&
           hxOmin2:hxOmax2))/(x(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
           idir)-x(hxOmin1:hxOmax1,hxOmin2:hxOmax2,idir))
      case('spherical')
        select case(idir)
        case(1)
          gradq(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(q(jxOmin1:jxOmax1,&
             jxOmin2:jxOmax2)-q(hxOmin1:hxOmax1,&
             hxOmin2:hxOmax2))/((x(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
             1)-x(hxOmin1:hxOmax1,hxOmin2:hxOmax2,1)))
          
        case(2)
          gradq(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(q(jxOmin1:jxOmax1,&
             jxOmin2:jxOmax2)-q(hxOmin1:hxOmax1,&
             hxOmin2:hxOmax2))/((x(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
             2)-x(hxOmin1:hxOmax1,hxOmin2:hxOmax2,2))*x(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,1))
         
          
        end select
      case('cylindrical')
        if(idir==phi_) then
          gradq(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(q(jxOmin1:jxOmax1,&
             jxOmin2:jxOmax2)-q(hxOmin1:hxOmax1,&
             hxOmin2:hxOmax2))/((x(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
             phi_)-x(hxOmin1:hxOmax1,hxOmin2:hxOmax2,phi_))*x(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,r_))
        else
          gradq(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(q(jxOmin1:jxOmax1,&
             jxOmin2:jxOmax2)-q(hxOmin1:hxOmax1,&
             hxOmin2:hxOmax2))/(x(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
             idir)-x(hxOmin1:hxOmax1,hxOmin2:hxOmax2,idir))
        end if
      case default
        call mpistop('Unknown geometry')
      end select
    end if

  end subroutine gradient

  !> Calculate gradient of a scalar q within ixL in direction idir
  !> first use limiter to go from cell center to edge
  subroutine gradientS(q,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,idir,gradq)
    use mod_global_parameters
    use mod_limiter
    use mod_ppm

    integer, intent(in)                :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2, idir
    double precision, intent(in)       :: q(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision, intent(inout)    :: gradq(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    double precision ,dimension(ixImin1:ixImax1,ixImin2:ixImax2) :: qC,qL,qR,&
       dqC,ldq,rdq

    double precision :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision :: invdx
    integer          :: hxOmin1,hxOmin2,hxOmax1,hxOmax2,ixCmin1,ixCmin2,&
       ixCmax1,ixCmax2,jxCmin1,jxCmin2,jxCmax1,jxCmax2,gxCmin1,gxCmin2,gxCmax1,&
       gxCmax2,hxCmin1,hxCmin2,hxCmax1,hxCmax2

    x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)=block%x(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:ndim)

    invdx=1.d0/dxlevel(idir)
    hxOmin1=ixOmin1-kr(idir,1);hxOmin2=ixOmin2-kr(idir,2)
    hxOmax1=ixOmax1-kr(idir,1);hxOmax2=ixOmax2-kr(idir,2);
    ixCmin1=hxOmin1;ixCmin2=hxOmin2;ixCmax1=ixOmax1;ixCmax2=ixOmax2;
    jxCmin1=ixCmin1+kr(idir,1);jxCmin2=ixCmin2+kr(idir,2)
    jxCmax1=ixCmax1+kr(idir,1);jxCmax2=ixCmax2+kr(idir,2);
    gxCmin1=ixCmin1-kr(idir,1);gxCmin2=ixCmin2-kr(idir,2);gxCmax1=jxCmax1
    gxCmax2=jxCmax2;
    hxCmin1=gxCmin1+kr(idir,1);hxCmin2=gxCmin2+kr(idir,2)
    hxCmax1=gxCmax1+kr(idir,1);hxCmax2=gxCmax2+kr(idir,2);

    ! set the gradient limiter here
    qR(gxCmin1:gxCmax1,gxCmin2:gxCmax2) = q(hxCmin1:hxCmax1,hxCmin2:hxCmax2)
    qL(gxCmin1:gxCmax1,gxCmin2:gxCmax2) = q(gxCmin1:gxCmax1,gxCmin2:gxCmax2)
    if (typegradlimiter/=limiter_ppm) then
      dqC(gxCmin1:gxCmax1,gxCmin2:gxCmax2)= qR(gxCmin1:gxCmax1,&
         gxCmin2:gxCmax2)-qL(gxCmin1:gxCmax1,gxCmin2:gxCmax2)
      call dwlimiter2(dqC,ixImin1,ixImin2,ixImax1,ixImax2,gxCmin1,gxCmin2,&
         gxCmax1,gxCmax2,idir,typegradlimiter,ldw=ldq,rdw=rdq)
      qL(ixCmin1:ixCmax1,ixCmin2:ixCmax2) = qL(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2) + half*ldq(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
      qR(ixCmin1:ixCmax1,ixCmin2:ixCmax2) = qR(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2) - half*rdq(jxCmin1:jxCmax1,jxCmin2:jxCmax2)
    else
      call PPMlimitervar(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,idir,q,q,qL,qR)
    endif

    if(slab) then
      gradq(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=half*(qR(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)-qL(hxOmin1:hxOmax1,hxOmin2:hxOmax2))*invdx
    else
      gradq(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(qR(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)-qL(hxOmin1:hxOmax1,&
         hxOmin2:hxOmax2))/block%dx(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir)
      select case(typeaxial)
      case('spherical')
        select case(idir)
        case(2)
          gradq(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=gradq(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2)/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)
          
        end select
      case('cylindrical')
        if(idir==phi_) gradq(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2)=gradq(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2)/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)
      end select
    end if

  end subroutine gradientS

  !> Calculate divergence of a vector qvec within ixL
  subroutine divvector(qvec,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,divq, fourthorder)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: qvec(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndir)
    double precision, intent(inout) :: divq(ixImin1:ixImax1,ixImin2:ixImax2)
    logical, intent(in), optional   :: fourthorder !< Default: false
    logical                         :: use_4th_order
    double precision                :: qC(ixImin1:ixImax1,ixImin2:ixImax2),&
        invdx(1:ndim)
    integer                         :: jxOmin1,jxOmin2,jxOmax1,jxOmax2,&
        hxOmin1,hxOmin2,hxOmax1,hxOmax2, ixCmin1,ixCmin2,ixCmax1,ixCmax2,&
        jxCmin1,jxCmin2,jxCmax1,jxCmax2
    integer                         :: idims, ixmin1,ixmin2,ixmax1,ixmax2,&
        gxOmin1,gxOmin2,gxOmax1,gxOmax2, kxOmin1,kxOmin2,kxOmax1,kxOmax2

    use_4th_order = .false.
    if (present(fourthorder)) use_4th_order = fourthorder

    if (use_4th_order) then
      if (.not. slab) call mpistop&
         ("divvector: 4th order only supported for slab geometry")
      ! Fourth order, stencil width is two
      ixmin1=ixOmin1-2;ixmin2=ixOmin2-2;ixmax1=ixOmax1+2;ixmax2=ixOmax2+2;
    else
      ! Second order, stencil width is one
      ixmin1=ixOmin1-1;ixmin2=ixOmin2-1;ixmax1=ixOmax1+1;ixmax2=ixOmax2+1;
    end if

    if (ixImin1>ixmin1.or.ixImax1<ixmax1.or.ixImin2>ixmin2.or.ixImax2<ixmax2) &
       call mpistop("Error in divvector: Non-conforming input limits")

    invdx=1.d0/dxlevel
    divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=0.0d0

    if (slab) then
      do idims=1,ndim
        if (.not. use_4th_order) then
          ! Use second order scheme
          jxOmin1=ixOmin1+kr(idims,1);jxOmin2=ixOmin2+kr(idims,2)
          jxOmax1=ixOmax1+kr(idims,1);jxOmax2=ixOmax2+kr(idims,2);
          hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
          hxOmax1=ixOmax1-kr(idims,1);hxOmax2=ixOmax2-kr(idims,2);
          divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=divq(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2)+half*(qvec(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
             idims) - qvec(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
             idims))*invdx(idims)
        else
          ! Use fourth order scheme
          kxOmin1=ixOmin1+2*kr(idims,1);kxOmin2=ixOmin2+2*kr(idims,2)
          kxOmax1=ixOmax1+2*kr(idims,1);kxOmax2=ixOmax2+2*kr(idims,2);
          jxOmin1=ixOmin1+kr(idims,1);jxOmin2=ixOmin2+kr(idims,2)
          jxOmax1=ixOmax1+kr(idims,1);jxOmax2=ixOmax2+kr(idims,2);
          hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
          hxOmax1=ixOmax1-kr(idims,1);hxOmax2=ixOmax2-kr(idims,2);
          gxOmin1=ixOmin1-2*kr(idims,1);gxOmin2=ixOmin2-2*kr(idims,2)
          gxOmax1=ixOmax1-2*kr(idims,1);gxOmax2=ixOmax2-2*kr(idims,2);
          divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=divq(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2)+(-qvec(kxOmin1:kxOmax1,kxOmin2:kxOmax2,&
             idims) + 8.0d0 * qvec(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
             idims) - 8.0d0 * qvec(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
             idims) + qvec(gxOmin1:gxOmax1,gxOmin2:gxOmax2,&
             idims))/(12.0d0 * dxlevel(idims))
        end if
      end do
    else
      do idims=1,ndim
        hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
        hxOmax1=ixOmax1-kr(idims,1);hxOmax2=ixOmax2-kr(idims,2);
        ixCmin1=hxOmin1;ixCmin2=hxOmin2;ixCmax1=ixOmax1;ixCmax2=ixOmax2;
        jxCmin1=ixCmin1+kr(idims,1);jxCmin2=ixCmin2+kr(idims,2)
        jxCmax1=ixCmax1+kr(idims,1);jxCmax2=ixCmax2+kr(idims,2);
        if(stretched_dim(idims).or.stretched_symm_dim(idims)) then
          ! linear interpolation at cell interface along stretched dimension
          qC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=block%surfaceC(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,idims)*(qvec(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             idims)+0.5d0*block%dx(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             idims)*(qvec(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
             idims)-qvec(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             idims))/(block%x(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
             idims)-block%x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idims)))
        else
          qC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=block%surfaceC(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,idims)*half*(qvec(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             idims)+qvec(jxCmin1:jxCmax1,jxCmin2:jxCmax2,idims))
        end if
        divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=divq(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2)+qC(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2)-qC(hxOmin1:hxOmax1,hxOmin2:hxOmax2)
      end do
      divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=divq(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)/block%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    end if


  end subroutine divvector

  !> Calculate divergence of a vector qvec within ixL
  !> using limited extrapolation to cell edges
  subroutine divvectorS(qvec,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,divq)
    use mod_global_parameters
    use mod_limiter
    use mod_ppm

    integer, intent(in)                :: ixImin1,ixImin2,ixImax1,ixImax2,&
       ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)       :: qvec(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndir)
    double precision, intent(inout)    :: divq(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2) :: qL,qR,dqC,&
       ldq,rdq

    double precision :: invdx(1:ndim)
    integer          :: hxOmin1,hxOmin2,hxOmax1,hxOmax2,ixCmin1,ixCmin2,&
       ixCmax1,ixCmax2,jxCmin1,jxCmin2,jxCmax1,jxCmax2,idims,ixmin1,ixmin2,&
       ixmax1,ixmax2,gxCmin1,gxCmin2,gxCmax1,gxCmax2,hxCmin1,hxCmin2,hxCmax1,&
       hxCmax2

    ixmin1=ixOmin1-2;ixmin2=ixOmin2-2;ixmax1=ixOmax1+2;ixmax2=ixOmax2+2;

    if (ixImin1>ixmin1.or.ixImax1<ixmax1.or.ixImin2>ixmin2.or.ixImax2<ixmax2) &
       call mpistop("Error in divvectorS: Non-conforming input limits")

    invdx=1.d0/dxlevel
    divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=zero
    do idims=1,ndim
      hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
      hxOmax1=ixOmax1-kr(idims,1);hxOmax2=ixOmax2-kr(idims,2);
      ixCmin1=hxOmin1;ixCmin2=hxOmin2;ixCmax1=ixOmax1;ixCmax2=ixOmax2;
      jxCmin1=ixCmin1+kr(idims,1);jxCmin2=ixCmin2+kr(idims,2)
      jxCmax1=ixCmax1+kr(idims,1);jxCmax2=ixCmax2+kr(idims,2);
      gxCmin1=ixCmin1-kr(idims,1);gxCmin2=ixCmin2-kr(idims,2);gxCmax1=jxCmax1
      gxCmax2=jxCmax2;
      hxCmin1=gxCmin1+kr(idims,1);hxCmin2=gxCmin2+kr(idims,2)
      hxCmax1=gxCmax1+kr(idims,1);hxCmax2=gxCmax2+kr(idims,2);
      qR(gxCmin1:gxCmax1,gxCmin2:gxCmax2) = qvec(hxCmin1:hxCmax1,&
         hxCmin2:hxCmax2,idims)
      qL(gxCmin1:gxCmax1,gxCmin2:gxCmax2) = qvec(gxCmin1:gxCmax1,&
         gxCmin2:gxCmax2,idims)
      if(typegradlimiter/=limiter_ppm) then
        dqC(gxCmin1:gxCmax1,gxCmin2:gxCmax2)= qR(gxCmin1:gxCmax1,&
           gxCmin2:gxCmax2)-qL(gxCmin1:gxCmax1,gxCmin2:gxCmax2)
        call dwlimiter2(dqC,ixImin1,ixImin2,ixImax1,ixImax2,gxCmin1,gxCmin2,&
           gxCmax1,gxCmax2,idims,typegradlimiter,ldw=ldq,rdw=rdq)
        qL(ixCmin1:ixCmax1,ixCmin2:ixCmax2) = qL(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2) + half*ldq(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
        qR(ixCmin1:ixCmax1,ixCmin2:ixCmax2) = qR(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2) - half*rdq(jxCmin1:jxCmax1,jxCmin2:jxCmax2)
      else
        dqC(ixImin1:ixImax1,ixImin2:ixImax2)=qvec(ixImin1:ixImax1,&
           ixImin2:ixImax2,idims)
        call PPMlimitervar(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
           ixOmax1,ixOmax2,idims,dqC,dqC,qL,qR)
      endif

      if (slab) then
        divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=divq(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2)+half*(qR(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2)-qL(hxOmin1:hxOmax1,hxOmin2:hxOmax2))*invdx(idims)
      else
        qR(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=block%surfaceC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,idims)*qR(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
        qL(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=block%surfaceC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,idims)*qL(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
        divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=divq(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2)+qR(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2)-qL(hxOmin1:hxOmax1,hxOmin2:hxOmax2)
      end if
    end do
    if(.not.slab) divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=divq(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)/block%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

  end subroutine divvectorS

  !> Calculate curl of a vector qvec within ixL
  !> Options to
  !>        employ standard second order CD evaluations
  !>        use Gauss theorem for non-Cartesian grids
  !>        use Stokes theorem for non-Cartesian grids
  subroutine curlvector(qvec,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,curlvec,idirmin,idirmin0,ndir0)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    integer, intent(in)             :: ndir0, idirmin0
    integer, intent(inout)          :: idirmin
    double precision, intent(in)    :: qvec(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndir0)
    double precision, intent(inout) :: curlvec(ixImin1:ixImax1,ixImin2:ixImax2,&
       idirmin0:3)

    integer          :: ixCmin1,ixCmin2,ixCmax1,ixCmax2,jxCmin1,jxCmin2,&
       jxCmax1,jxCmax2,idir,jdir,kdir,hxOmin1,hxOmin2,hxOmax1,hxOmax2,jxOmin1,&
       jxOmin2,jxOmax1,jxOmax2
    double precision :: invdx(1:ndim)
    double precision :: tmp(ixImin1:ixImax1,ixImin2:ixImax2),&
       tmp2(ixImin1:ixImax1,ixImin2:ixImax2),xC(ixImin1:ixImax1,&
       ixImin2:ixImax2),surface(ixImin1:ixImax1,ixImin2:ixImax2)
    !-----------------------------------------------------------------------------

    ! Calculate curl within ixL: CurlV_i=eps_ijk*d_j V_k
    ! Curl can have components (idirmin:3)
    ! Determine exact value of idirmin while doing the loop.

    idirmin=4
    curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idirmin0:3)=zero

    if(slab) then ! Cartesian case
      invdx=1.d0/dxlevel
      do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir0
        if(lvc(idir,jdir,kdir)/=0)then
          tmp(ixImin1:ixImax1,ixImin2:ixImax2)=qvec(ixImin1:ixImax1,&
             ixImin2:ixImax2,kdir)
          hxOmin1=ixOmin1-kr(jdir,1);hxOmin2=ixOmin2-kr(jdir,2)
          hxOmax1=ixOmax1-kr(jdir,1);hxOmax2=ixOmax2-kr(jdir,2);
          jxOmin1=ixOmin1+kr(jdir,1);jxOmin2=ixOmin2+kr(jdir,2)
          jxOmax1=ixOmax1+kr(jdir,1);jxOmax2=ixOmax2+kr(jdir,2);
          ! second order centered differencing
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=half*(tmp(jxOmin1:jxOmax1,&
             jxOmin2:jxOmax2)-tmp(hxOmin1:hxOmax1,&
             hxOmin2:hxOmax2))*invdx(jdir)
          !> \todo allow for 4th order CD evaluation here as well
          if(lvc(idir,jdir,kdir)==1)then
            curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               idir)=curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               idir)+tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
          else
            curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               idir)=curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               idir)-tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
          endif
          if(idir<idirmin)idirmin=idir
        endif
      enddo; enddo; enddo;
      return
    endif

    ! all non-Cartesian cases
    select case(typeaxial)
      case('slabstretch') ! stretched Cartesian grids
        do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir0
          if(lvc(idir,jdir,kdir)/=0)then
            select case(typecurl)
              case('central')
                tmp(ixImin1:ixImax1,ixImin2:ixImax2)=qvec(ixImin1:ixImax1,&
                   ixImin2:ixImax2,kdir)
                hxOmin1=ixOmin1-kr(jdir,1);hxOmin2=ixOmin2-kr(jdir,2)
                hxOmax1=ixOmax1-kr(jdir,1);hxOmax2=ixOmax2-kr(jdir,2);
                jxOmin1=ixOmin1+kr(jdir,1);jxOmin2=ixOmin2+kr(jdir,2)
                jxOmax1=ixOmax1+kr(jdir,1);jxOmax2=ixOmax2+kr(jdir,2);
                ! second order centered differencing
                tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(tmp(jxOmin1:jxOmax1,&
                   jxOmin2:jxOmax2)-tmp(hxOmin1:hxOmax1,&
                   hxOmin2:hxOmax2))/(block%x(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
                   jdir)-block%x(hxOmin1:hxOmax1,hxOmin2:hxOmax2,jdir))
              case('Gaussbased')
                hxOmin1=ixOmin1-kr(jdir,1);hxOmin2=ixOmin2-kr(jdir,2)
                hxOmax1=ixOmax1-kr(jdir,1);hxOmax2=ixOmax2-kr(jdir,2);
                ixCmin1=hxOmin1;ixCmin2=hxOmin2;ixCmax1=ixOmax1
                ixCmax2=ixOmax2;
                jxCmin1=ixCmin1+kr(jdir,1);jxCmin2=ixCmin2+kr(jdir,2)
                jxCmax1=ixCmax1+kr(jdir,1);jxCmax2=ixCmax2+kr(jdir,2);
                tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=block%surfaceC(&
                   ixCmin1:ixCmax1,ixCmin2:ixCmax2,jdir)*(qvec(ixCmin1:ixCmax1,&
                   ixCmin2:ixCmax2,kdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                   ixCmin2:ixCmax2,jdir)*(qvec(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                   kdir)-qvec(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                   kdir))/(block%x(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                   jdir)-block%x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,jdir)))
                tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(tmp(ixOmin1:ixOmax1,&
                   ixOmin2:ixOmax2)-tmp(hxOmin1:hxOmax1,&
                   hxOmin2:hxOmax2))/block%dvolume(ixOmin1:ixOmax1,&
                   ixOmin2:ixOmax2)
              case('Stokesbased')
                hxOmin1=ixOmin1-kr(jdir,1);hxOmin2=ixOmin2-kr(jdir,2)
                hxOmax1=ixOmax1-kr(jdir,1);hxOmax2=ixOmax2-kr(jdir,2);
                ixCmin1=hxOmin1;ixCmin2=hxOmin2;ixCmax1=ixOmax1
                ixCmax2=ixOmax2;
                jxCmin1=ixCmin1+kr(jdir,1);jxCmin2=ixCmin2+kr(jdir,2)
                jxCmax1=ixCmax1+kr(jdir,1);jxCmax2=ixCmax2+kr(jdir,2);
                if(kdir<=ndim)then
                  tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=block%ds(&
                     ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     kdir)*(qvec(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                     kdir)+0.5d0*block%dx(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                     jdir)*(qvec(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                     kdir)-qvec(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                     kdir))/(block%x(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                     jdir)-block%x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,jdir)))
                else
                  tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=(qvec(ixCmin1:ixCmax1,&
                     ixCmin2:ixCmax2,kdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                     ixCmin2:ixCmax2,jdir)*(qvec(jxCmin1:jxCmax1,&
                     jxCmin2:jxCmax2,kdir)-qvec(ixCmin1:ixCmax1,&
                     ixCmin2:ixCmax2,kdir))/(block%x(jxCmin1:jxCmax1,&
                     jxCmin2:jxCmax2,jdir)-block%x(ixCmin1:ixCmax1,&
                     ixCmin2:ixCmax2,jdir)))
                endif
                if(idir<=ndim)then
                  tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(tmp(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2)-tmp(hxOmin1:hxOmax1,&
                     hxOmin2:hxOmax2))/block%surface(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2,idir)
                else ! essentially for 2.5D case, idir=3 and jdir,kdir<=2
                  tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(tmp(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2)-tmp(hxOmin1:hxOmax1,&
                     hxOmin2:hxOmax2))/(block%ds(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2,jdir)*block%ds(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2,kdir))
                endif
              case default
                call mpistop('no such curl evaluator')
            end select
            if(lvc(idir,jdir,kdir)==1)then
              curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 idir)=curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 idir)+tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
            else
              curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 idir)=curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 idir)-tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
            endif
            if(idir<idirmin)idirmin=idir
          endif
        enddo; enddo; enddo;
      case('spherical') ! possibly stretched spherical grids
        select case(typecurl)
          case('central') ! ok for any dimensionality
            do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir0
              if(lvc(idir,jdir,kdir)/=0)then
                tmp(ixImin1:ixImax1,ixImin2:ixImax2)=qvec(ixImin1:ixImax1,&
                   ixImin2:ixImax2,kdir)
                hxOmin1=ixOmin1-kr(jdir,1);hxOmin2=ixOmin2-kr(jdir,2)
                hxOmax1=ixOmax1-kr(jdir,1);hxOmax2=ixOmax2-kr(jdir,2);
                jxOmin1=ixOmin1+kr(jdir,1);jxOmin2=ixOmin2+kr(jdir,2)
                jxOmax1=ixOmax1+kr(jdir,1);jxOmax2=ixOmax2+kr(jdir,2);
                select case(jdir)
                case(1)
                tmp(ixImin1:ixImax1,ixImin2:ixImax2)=tmp(ixImin1:ixImax1,&
                   ixImin2:ixImax2)*block%x(ixImin1:ixImax1,ixImin2:ixImax2,1)
                tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(tmp(jxOmin1:jxOmax1,&
                   jxOmin2:jxOmax2)-tmp(hxOmin1:hxOmax1,&
                   hxOmin2:hxOmax2))/((block%x(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
                   1)-block%x(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
                   1))*block%x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1))
                    case(2)
                if(idir==1) tmp(ixImin1:ixImax1,&
                   ixImin2:ixImax2)=tmp(ixImin1:ixImax1,&
                   ixImin2:ixImax2)*dsin(block%x(ixImin1:ixImax1,&
                   ixImin2:ixImax2,2))
                tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(tmp(jxOmin1:jxOmax1,&
                   jxOmin2:jxOmax2)-tmp(hxOmin1:hxOmax1,&
                   hxOmin2:hxOmax2))/((block%x(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
                   2)-block%x(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
                   2))*block%x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1))
                if(idir==1) tmp2(ixOmin1:ixOmax1,&
                   ixOmin2:ixOmax2)=tmp2(ixOmin1:ixOmax1,&
                   ixOmin2:ixOmax2)/dsin(block%x(ixOmin1:ixOmax1,&
                   ixOmin2:ixOmax2,2))
               
                
                end select
                if(lvc(idir,jdir,kdir)==1)then
                  curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     idir)=curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     idir)+tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
                else
                  curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     idir)=curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     idir)-tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
                endif
                if(idir<idirmin)idirmin=idir
              endif
            enddo; enddo; enddo;
          case('Gaussbased')
            do idir=idirmin0,3;
            do jdir=1,ndim; do kdir=1,ndir0
              if(lvc(idir,jdir,kdir)/=0)then
                hxOmin1=ixOmin1-kr(jdir,1);hxOmin2=ixOmin2-kr(jdir,2)
                hxOmax1=ixOmax1-kr(jdir,1);hxOmax2=ixOmax2-kr(jdir,2);
                ixCmin1=hxOmin1;ixCmin2=hxOmin2;ixCmax1=ixOmax1
                ixCmax2=ixOmax2;
                jxCmin1=ixCmin1+kr(jdir,1);jxCmin2=ixCmin2+kr(jdir,2)
                jxCmax1=ixCmax1+kr(jdir,1);jxCmax2=ixCmax2+kr(jdir,2);
                tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=block%surfaceC(&
                   ixCmin1:ixCmax1,ixCmin2:ixCmax2,jdir)*(qvec(ixCmin1:ixCmax1,&
                   ixCmin2:ixCmax2,kdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                   ixCmin2:ixCmax2,jdir)*(qvec(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                   kdir)-qvec(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                   kdir))/(block%x(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                   jdir)-block%x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,jdir)))
                tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(tmp(ixOmin1:ixOmax1,&
                   ixOmin2:ixOmax2)-tmp(hxOmin1:hxOmax1,&
                   hxOmin2:hxOmax2))/block%dvolume(ixOmin1:ixOmax1,&
                   ixOmin2:ixOmax2)
                if(lvc(idir,jdir,kdir)==1)then
                  curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     idir)=curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     idir)+tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
                else
                  curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     idir)=curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     idir)-tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
                 endif
                 if(idir<idirmin)idirmin=idir
              endif
            enddo; enddo;
            ! geometric terms
            if(idir==2.and.phi_>0) curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               2)=curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               2)+qvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               phi_)/block%x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,r_)
            
            if(idir==phi_) curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               phi_)=curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               phi_)-qvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               2)/block%x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               r_) +qvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               r_)*dcos(block%x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               2))/(block%x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               r_)*dsin(block%x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)))
           
            enddo;
          case('Stokesbased')
            !if(ndim<3) call mpistop("Stokesbased for 3D spherical only")
            do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir0
             if(lvc(idir,jdir,kdir)/=0)then
              select case(idir)
              case(1)
                if(jdir<kdir) then
                  ! idir=1,jdir=2,kdir=3
                  !! integral along 3rd dimension
                  hxOmin1=ixOmin1-kr(jdir,1);hxOmin2=ixOmin2-kr(jdir,2)
                  hxOmax1=ixOmax1-kr(jdir,1);hxOmax2=ixOmax2-kr(jdir,2);
                  ixCmin1=hxOmin1;ixCmin2=hxOmin2;ixCmax1=ixOmax1
                  ixCmax2=ixOmax2;
                  jxCmin1=ixCmin1+kr(jdir,1);jxCmin2=ixCmin2+kr(jdir,2)
                  jxCmax1=ixCmax1+kr(jdir,1);jxCmax2=ixCmax2+kr(jdir,2);
                  ! qvec(3) at cell interface along 2nd dimension
                  if(stretched_dim(jdir).or.stretched_symm_dim(jdir)) then
                    tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=qvec(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,kdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,jdir)*(qvec(jxCmin1:jxCmax1,&
                       jxCmin2:jxCmax2,kdir)-qvec(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,kdir))/(block%x(jxCmin1:jxCmax1,&
                       jxCmin2:jxCmax2,jdir)-block%x(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,jdir))
                  else
                    tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=0.5d0*(qvec(&
                       ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                       kdir)+qvec(jxCmin1:jxCmax1,jxCmin2:jxCmax2,kdir))
                  end if
                  ! 2nd coordinate at cell interface along 2nd dimension
                  xC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=block%x(ixCmin1:ixCmax1,&
                     ixCmin2:ixCmax2,jdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                     ixCmin2:ixCmax2,jdir)
                  curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     idir)=(dsin(xC(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2))*tmp(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2)-dsin(xC(hxOmin1:hxOmax1,&
                     hxOmin2:hxOmax2))*tmp(hxOmin1:hxOmax1,&
                     hxOmin2:hxOmax2))*block%dx(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2,kdir)
                  !! integral along 2nd dimension
                  hxOmin1=ixOmin1-kr(kdir,1);hxOmin2=ixOmin2-kr(kdir,2)
                  hxOmax1=ixOmax1-kr(kdir,1);hxOmax2=ixOmax2-kr(kdir,2);
                  ixCmin1=hxOmin1;ixCmin2=hxOmin2;ixCmax1=ixOmax1
                  ixCmax2=ixOmax2;
                  jxCmin1=ixCmin1+kr(kdir,1);jxCmin2=ixCmin2+kr(kdir,2)
                  jxCmax1=ixCmax1+kr(kdir,1);jxCmax2=ixCmax2+kr(kdir,2);
                  ! qvec(2) at cell interface along 3rd dimension
                  if(stretched_dim(kdir).or.stretched_symm_dim(kdir)) then
                    tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=qvec(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,jdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,kdir)*(qvec(jxCmin1:jxCmax1,&
                       jxCmin2:jxCmax2,jdir)-qvec(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,jdir))/(block%x(jxCmin1:jxCmax1,&
                       jxCmin2:jxCmax2,kdir)-block%x(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,kdir))
                  else
                    tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=0.5d0*(qvec(&
                       ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                       jdir)+qvec(jxCmin1:jxCmax1,jxCmin2:jxCmax2,jdir))
                  end if
                  curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     idir)=(curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     idir)+(tmp(hxOmin1:hxOmax1,&
                     hxOmin2:hxOmax2)-tmp(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2))*block%dx(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2,jdir))/block%surface(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2,idir)*block%x(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2,idir)
                end if
              case(2)
                if(jdir<kdir) then
                  ! idir=2,jdir=1,kdir=3
                  !! integral along 1st dimension
                  hxOmin1=ixOmin1-kr(kdir,1);hxOmin2=ixOmin2-kr(kdir,2)
                  hxOmax1=ixOmax1-kr(kdir,1);hxOmax2=ixOmax2-kr(kdir,2);
                  ixCmin1=hxOmin1;ixCmin2=hxOmin2;ixCmax1=ixOmax1
                  ixCmax2=ixOmax2;
                  jxCmin1=ixCmin1+kr(kdir,1);jxCmin2=ixCmin2+kr(kdir,2)
                  jxCmax1=ixCmax1+kr(kdir,1);jxCmax2=ixCmax2+kr(kdir,2);
                  ! qvec(1) at cell interface along 3rd dimension
                  if(stretched_dim(kdir).or.stretched_symm_dim(kdir)) then
                    tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=qvec(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,jdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,kdir)*(qvec(jxCmin1:jxCmax1,&
                       jxCmin2:jxCmax2,jdir)-qvec(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,jdir))/(block%x(jxCmin1:jxCmax1,&
                       jxCmin2:jxCmax2,kdir)-block%x(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,kdir))
                  else
                    tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=0.5d0*(qvec(&
                       ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                       jdir)+qvec(jxCmin1:jxCmax1,jxCmin2:jxCmax2,jdir))
                  end if
                  curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     idir)=(tmp(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2)-tmp(hxOmin1:hxOmax1,&
                     hxOmin2:hxOmax2))*block%dx(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2,1)
                  !! integral along 3rd dimension
                  hxOmin1=ixOmin1-kr(jdir,1);hxOmin2=ixOmin2-kr(jdir,2)
                  hxOmax1=ixOmax1-kr(jdir,1);hxOmax2=ixOmax2-kr(jdir,2);
                  ixCmin1=hxOmin1;ixCmin2=hxOmin2;ixCmax1=ixOmax1
                  ixCmax2=ixOmax2;
                  jxCmin1=ixCmin1+kr(jdir,1);jxCmin2=ixCmin2+kr(jdir,2)
                  jxCmax1=ixCmax1+kr(jdir,1);jxCmax2=ixCmax2+kr(jdir,2);
                  ! qvec(3) at cell interface along 1st dimension
                  if(stretched_dim(jdir).or.stretched_symm_dim(jdir)) then
                    tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=qvec(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,kdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,jdir)*(qvec(jxCmin1:jxCmax1,&
                       jxCmin2:jxCmax2,kdir)-qvec(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,kdir))/(block%x(jxCmin1:jxCmax1,&
                       jxCmin2:jxCmax2,jdir)-block%x(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,jdir))
                  else
                    tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=0.5d0*(qvec(&
                       ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                       kdir)+qvec(jxCmin1:jxCmax1,jxCmin2:jxCmax2,kdir))
                  end if
                  ! 1st coordinate at cell interface along 1st dimension
                  xC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=block%x(ixCmin1:ixCmax1,&
                     ixCmin2:ixCmax2,jdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                     ixCmin2:ixCmax2,jdir)
                  curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     idir)=(curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     idir)+(xC(hxOmin1:hxOmax1,&
                     hxOmin2:hxOmax2)*tmp(hxOmin1:hxOmax1,&
                     hxOmin2:hxOmax2)-xC(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2)*tmp(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2))*dsin(block%x(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2,idir))*block%dx(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2,kdir))/block%surface(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2,idir)
                end if
              case(3)
                if(jdir<kdir) then
                  ! idir=3,jdir=1,kdir=2
                  !! integral along 1st dimension
                  hxOmin1=ixOmin1-kr(kdir,1);hxOmin2=ixOmin2-kr(kdir,2)
                  hxOmax1=ixOmax1-kr(kdir,1);hxOmax2=ixOmax2-kr(kdir,2);
                  ixCmin1=hxOmin1;ixCmin2=hxOmin2;ixCmax1=ixOmax1
                  ixCmax2=ixOmax2;
                  jxCmin1=ixCmin1+kr(kdir,1);jxCmin2=ixCmin2+kr(kdir,2)
                  jxCmax1=ixCmax1+kr(kdir,1);jxCmax2=ixCmax2+kr(kdir,2);
                  ! qvec(1) at cell interface along 2nd dimension
                  if(stretched_dim(kdir).or.stretched_symm_dim(kdir)) then
                    tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=qvec(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,jdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,kdir)*(qvec(jxCmin1:jxCmax1,&
                       jxCmin2:jxCmax2,jdir)-qvec(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,jdir))/(block%x(jxCmin1:jxCmax1,&
                       jxCmin2:jxCmax2,kdir)-block%x(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,kdir))
                  else
                    tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=0.5d0*(qvec(&
                       ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                       jdir)+qvec(jxCmin1:jxCmax1,jxCmin2:jxCmax2,jdir))
                  end if
                  curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     idir)=(tmp(hxOmin1:hxOmax1,&
                     hxOmin2:hxOmax2)-tmp(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2))*block%dx(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2,jdir)
                  !! integral along 2nd dimension
                  hxOmin1=ixOmin1-kr(jdir,1);hxOmin2=ixOmin2-kr(jdir,2)
                  hxOmax1=ixOmax1-kr(jdir,1);hxOmax2=ixOmax2-kr(jdir,2);
                  ixCmin1=hxOmin1;ixCmin2=hxOmin2;ixCmax1=ixOmax1
                  ixCmax2=ixOmax2;
                  jxCmin1=ixCmin1+kr(jdir,1);jxCmin2=ixCmin2+kr(jdir,2)
                  jxCmax1=ixCmax1+kr(jdir,1);jxCmax2=ixCmax2+kr(jdir,2);
                  ! qvec(2) at cell interface along 1st dimension
                  if(stretched_dim(jdir).or.stretched_symm_dim(jdir)) then
                    tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=qvec(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,kdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,jdir)*(qvec(jxCmin1:jxCmax1,&
                       jxCmin2:jxCmax2,kdir)-qvec(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,kdir))/(block%x(jxCmin1:jxCmax1,&
                       jxCmin2:jxCmax2,jdir)-block%x(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,jdir))
                  else
                    tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=0.5d0*(qvec(&
                       ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                       kdir)+qvec(jxCmin1:jxCmax1,jxCmin2:jxCmax2,kdir))
                  end if
                  ! 1st coordinate at cell interface along 1st dimension
                  xC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=block%x(ixCmin1:ixCmax1,&
                     ixCmin2:ixCmax2,jdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                     ixCmin2:ixCmax2,jdir)
                  if(ndim==3) then
                    surface(ixOmin1:ixOmax1,&
                       ixOmin2:ixOmax2)=block%surface(ixOmin1:ixOmax1,&
                       ixOmin2:ixOmax2,idir)
                  else
                    surface(ixOmin1:ixOmax1,&
                       ixOmin2:ixOmax2)=block%x(ixOmin1:ixOmax1,&
                       ixOmin2:ixOmax2,jdir)*block%dx(ixOmin1:ixOmax1,&
                       ixOmin2:ixOmax2,kdir)*block%dx(ixOmin1:ixOmax1,&
                       ixOmin2:ixOmax2,jdir)
                  end if
                  curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     idir)=(curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     idir)+(xC(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2)*tmp(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2)-xC(hxOmin1:hxOmax1,&
                     hxOmin2:hxOmax2)*tmp(hxOmin1:hxOmax1,&
                     hxOmin2:hxOmax2))*block%dx(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2,kdir))/surface(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2)
                end if
              end select
              if(idir<idirmin)idirmin=idir
             endif
            enddo; enddo; enddo;
          case default
            call mpistop('no such curl evaluator')
        end select
      case('cylindrical') ! possibly stretched cylindrical grids
        select case(typecurl)
          case('central')  ! works for any dimensionality, polar/cylindrical
            do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir0
              if(lvc(idir,jdir,kdir)/=0)then
                tmp(ixImin1:ixImax1,ixImin2:ixImax2)=qvec(ixImin1:ixImax1,&
                   ixImin2:ixImax2,kdir)
                hxOmin1=ixOmin1-kr(jdir,1);hxOmin2=ixOmin2-kr(jdir,2)
                hxOmax1=ixOmax1-kr(jdir,1);hxOmax2=ixOmax2-kr(jdir,2);
                jxOmin1=ixOmin1+kr(jdir,1);jxOmin2=ixOmin2+kr(jdir,2)
                jxOmax1=ixOmax1+kr(jdir,1);jxOmax2=ixOmax2+kr(jdir,2);
                if(z_==3.or.z_==-1) then
                  ! Case Polar_2D, Polar_2.5D or Polar_3D, i.e. R,phi,Z
                  select case(jdir)
                  case(1)
                  if(idir==z_) tmp(ixImin1:ixImax1,&
                     ixImin2:ixImax2)=tmp(ixImin1:ixImax1,&
                     ixImin2:ixImax2)*block%x(ixImin1:ixImax1,ixImin2:ixImax2,&
                     1) !R V_phi
                  ! computes d(R V_phi)/dR or d V_Z/dR
                  tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(tmp(jxOmin1:jxOmax1,&
                     jxOmin2:jxOmax2)-tmp(hxOmin1:hxOmax1,&
                     hxOmin2:hxOmax2))/(block%x(jxOmin1:jxOmax1,&
                     jxOmin2:jxOmax2,1)-block%x(hxOmin1:hxOmax1,&
                     hxOmin2:hxOmax2,1))
                  if(idir==z_) tmp2(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2)=tmp2(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2)/block%x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     1) !(1/R)*d(R V_phi)/dR
                        case(2)
                  ! handles (1/R)d V_Z/dphi or (1/R)d V_R/dphi
                  tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(tmp(jxOmin1:jxOmax1,&
                     jxOmin2:jxOmax2)-tmp(hxOmin1:hxOmax1,&
                     hxOmin2:hxOmax2))/((block%x(jxOmin1:jxOmax1,&
                     jxOmin2:jxOmax2,2)-block%x(hxOmin1:hxOmax1,&
                     hxOmin2:hxOmax2,2))*block%x(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2,1))
                 
                  
                  end select
                end if
                if(phi_==3.or.phi_==-1) then
                  ! Case Cylindrical_2D, Cylindrical_2.5D or Cylindrical_3D, i.e. R,Z,phi
                  select case(jdir)
                  case(1)
                  if(idir==z_) tmp(ixImin1:ixImax1,&
                     ixImin2:ixImax2)=tmp(ixImin1:ixImax1,&
                     ixImin2:ixImax2)*block%x(ixImin1:ixImax1,ixImin2:ixImax2,&
                     1) !R V_phi
                  ! computes d(R V_phi)/dR or d V_Z/dR
                  tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(tmp(jxOmin1:jxOmax1,&
                     jxOmin2:jxOmax2)-tmp(hxOmin1:hxOmax1,&
                     hxOmin2:hxOmax2))/(block%x(jxOmin1:jxOmax1,&
                     jxOmin2:jxOmax2,1)-block%x(hxOmin1:hxOmax1,&
                     hxOmin2:hxOmax2,1))
                  if(idir==z_) tmp2(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2)=tmp2(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2)/block%x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     1) !(1/R)*d(R V_phi)/dR
                        case(2)
                  ! handles d V_phi/dZ or d V_R/dZ
                  tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(tmp(jxOmin1:jxOmax1,&
                     jxOmin2:jxOmax2)-tmp(hxOmin1:hxOmax1,&
                     hxOmin2:hxOmax2))/(block%x(jxOmin1:jxOmax1,&
                     jxOmin2:jxOmax2,2)-block%x(hxOmin1:hxOmax1,&
                     hxOmin2:hxOmax2,2))
                 
                  
                  end select
                end if
                if(lvc(idir,jdir,kdir)==1)then
                  curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     idir)=curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     idir)+tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
                else
                  curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     idir)=curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     idir)-tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
                 endif
                 if(idir<idirmin)idirmin=idir
              endif
            enddo; enddo; enddo;
          case('Gaussbased') ! works for any dimensionality, polar/cylindrical
            if(ndim<2) call mpistop&
               ("Gaussbased for 2D, 2.5D or 3D polar or cylindrical only")
            do idir=idirmin0,3;
            do jdir=1,ndim; do kdir=1,ndir0
              if(lvc(idir,jdir,kdir)/=0)then
                hxOmin1=ixOmin1-kr(jdir,1);hxOmin2=ixOmin2-kr(jdir,2)
                hxOmax1=ixOmax1-kr(jdir,1);hxOmax2=ixOmax2-kr(jdir,2);
                ixCmin1=hxOmin1;ixCmin2=hxOmin2;ixCmax1=ixOmax1
                ixCmax2=ixOmax2;
                jxCmin1=ixCmin1+kr(jdir,1);jxCmin2=ixCmin2+kr(jdir,2)
                jxCmax1=ixCmax1+kr(jdir,1);jxCmax2=ixCmax2+kr(jdir,2);
                tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=block%surfaceC(&
                   ixCmin1:ixCmax1,ixCmin2:ixCmax2,jdir)*(qvec(ixCmin1:ixCmax1,&
                   ixCmin2:ixCmax2,kdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                   ixCmin2:ixCmax2,jdir)*(qvec(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                   kdir)-qvec(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                   kdir))/(block%x(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                   jdir)-block%x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,jdir)))
                tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(tmp(ixOmin1:ixOmax1,&
                   ixOmin2:ixOmax2)-tmp(hxOmin1:hxOmax1,&
                   hxOmin2:hxOmax2))/block%dvolume(ixOmin1:ixOmax1,&
                   ixOmin2:ixOmax2)
                if(lvc(idir,jdir,kdir)==1)then
                  curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     idir)=curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     idir)+tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
                else
                  curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     idir)=curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     idir)-tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
                 endif
                 if(idir<idirmin)idirmin=idir
              endif
            enddo; enddo;
            ! geometric term from d e_R/d phi= e_phi for unit vectors e_R, e_phi
            !       but minus sign appears due to R,Z,phi ordering (?)
            ! note that in cylindrical 2D (RZ), phi_ is -1
            ! note that in polar 2D     (Rphi), z_ is -1
            if((idir==phi_.or.(phi_==-1.and.idir==3)).and.z_>0) then
              ! cylindrical
              if(  z_==2) curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 idir)=curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 idir)-qvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 z_)/block%x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,r_)
              ! polar
              if(phi_==2) curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 idir)=curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 idir)+qvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 z_)/block%x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,r_)
            endif
            enddo;
          case('Stokesbased')
            !if(ndim<3) call mpistop("Stokesbased for 3D cylindrical only")
            do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir0
              if(lvc(idir,jdir,kdir)/=0)then
               if(idir==r_) then
                if(jdir==phi_) then
                  ! idir=r,jdir=phi,kdir=z
                  !! integral along z dimension
                  hxOmin1=ixOmin1-kr(jdir,1);hxOmin2=ixOmin2-kr(jdir,2)
                  hxOmax1=ixOmax1-kr(jdir,1);hxOmax2=ixOmax2-kr(jdir,2);
                  ixCmin1=hxOmin1;ixCmin2=hxOmin2;ixCmax1=ixOmax1
                  ixCmax2=ixOmax2;
                  jxCmin1=ixCmin1+kr(jdir,1);jxCmin2=ixCmin2+kr(jdir,2)
                  jxCmax1=ixCmax1+kr(jdir,1);jxCmax2=ixCmax2+kr(jdir,2);
                  ! qvec(z) at cell interface along phi dimension
                  if(stretched_dim(jdir).or.stretched_symm_dim(jdir)) then
                    tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=qvec(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,kdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,jdir)*(qvec(jxCmin1:jxCmax1,&
                       jxCmin2:jxCmax2,kdir)-qvec(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,kdir))/(block%x(jxCmin1:jxCmax1,&
                       jxCmin2:jxCmax2,jdir)-block%x(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,jdir))
                  else
                    tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=0.5d0*(qvec(&
                       ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                       kdir)+qvec(jxCmin1:jxCmax1,jxCmin2:jxCmax2,kdir))
                  end if
                  curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     idir)=(tmp(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2)-tmp(hxOmin1:hxOmax1,&
                     hxOmin2:hxOmax2))*block%dx(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2,kdir)
                  !! integral along phi dimension
                  hxOmin1=ixOmin1-kr(kdir,1);hxOmin2=ixOmin2-kr(kdir,2)
                  hxOmax1=ixOmax1-kr(kdir,1);hxOmax2=ixOmax2-kr(kdir,2);
                  ixCmin1=hxOmin1;ixCmin2=hxOmin2;ixCmax1=ixOmax1
                  ixCmax2=ixOmax2;
                  jxCmin1=ixCmin1+kr(kdir,1);jxCmin2=ixCmin2+kr(kdir,2)
                  jxCmax1=ixCmax1+kr(kdir,1);jxCmax2=ixCmax2+kr(kdir,2);
                  ! qvec(phi) at cell interface along z dimension
                  if(stretched_dim(kdir).or.stretched_symm_dim(kdir)) then
                    tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=qvec(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,jdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,kdir)*(qvec(jxCmin1:jxCmax1,&
                       jxCmin2:jxCmax2,jdir)-qvec(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,jdir))/(block%x(jxCmin1:jxCmax1,&
                       jxCmin2:jxCmax2,kdir)-block%x(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,kdir))
                  else
                    tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=0.5d0*(qvec(&
                       ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                       jdir)+qvec(jxCmin1:jxCmax1,jxCmin2:jxCmax2,jdir))
                  end if
                  curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     idir)=(curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     idir)+(tmp(hxOmin1:hxOmax1,&
                     hxOmin2:hxOmax2)-tmp(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2))*block%x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     idir)*block%dx(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     jdir))/block%surface(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     idir)
                end if
               else if(idir==phi_) then
                if(jdir<kdir) then
                  ! idir=phi,jdir=r,kdir=z
                  !! integral along r dimension
                  hxOmin1=ixOmin1-kr(kdir,1);hxOmin2=ixOmin2-kr(kdir,2)
                  hxOmax1=ixOmax1-kr(kdir,1);hxOmax2=ixOmax2-kr(kdir,2);
                  ixCmin1=hxOmin1;ixCmin2=hxOmin2;ixCmax1=ixOmax1
                  ixCmax2=ixOmax2;
                  jxCmin1=ixCmin1+kr(kdir,1);jxCmin2=ixCmin2+kr(kdir,2)
                  jxCmax1=ixCmax1+kr(kdir,1);jxCmax2=ixCmax2+kr(kdir,2);
                  ! qvec(r) at cell interface along z dimension
                  if(stretched_dim(kdir).or.stretched_symm_dim(kdir)) then
                    tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=qvec(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,jdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,kdir)*(qvec(jxCmin1:jxCmax1,&
                       jxCmin2:jxCmax2,jdir)-qvec(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,jdir))/(block%x(jxCmin1:jxCmax1,&
                       jxCmin2:jxCmax2,kdir)-block%x(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,kdir))
                  else
                    tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=0.5d0*(qvec(&
                       ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                       jdir)+qvec(jxCmin1:jxCmax1,jxCmin2:jxCmax2,jdir))
                  end if
                  curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     idir)=(tmp(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2)-tmp(hxOmin1:hxOmax1,&
                     hxOmin2:hxOmax2))*block%dx(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2,jdir)
                  !! integral along z dimension
                  hxOmin1=ixOmin1-kr(jdir,1);hxOmin2=ixOmin2-kr(jdir,2)
                  hxOmax1=ixOmax1-kr(jdir,1);hxOmax2=ixOmax2-kr(jdir,2);
                  ixCmin1=hxOmin1;ixCmin2=hxOmin2;ixCmax1=ixOmax1
                  ixCmax2=ixOmax2;
                  jxCmin1=ixCmin1+kr(jdir,1);jxCmin2=ixCmin2+kr(jdir,2)
                  jxCmax1=ixCmax1+kr(jdir,1);jxCmax2=ixCmax2+kr(jdir,2);
                  ! qvec(z) at cell interface along r dimension
                  if(stretched_dim(jdir).or.stretched_symm_dim(jdir)) then
                    tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=qvec(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,kdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,jdir)*(qvec(jxCmin1:jxCmax1,&
                       jxCmin2:jxCmax2,kdir)-qvec(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,kdir))/(block%x(jxCmin1:jxCmax1,&
                       jxCmin2:jxCmax2,jdir)-block%x(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,jdir))
                  else
                    tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=0.5d0*(qvec(&
                       ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                       kdir)+qvec(jxCmin1:jxCmax1,jxCmin2:jxCmax2,kdir))
                  end if
                  if(ndim==3) then
                    surface(ixOmin1:ixOmax1,&
                       ixOmin2:ixOmax2)=block%surface(ixOmin1:ixOmax1,&
                       ixOmin2:ixOmax2,idir)
                  else
                    surface(ixOmin1:ixOmax1,&
                       ixOmin2:ixOmax2)=block%dx(ixOmin1:ixOmax1,&
                       ixOmin2:ixOmax2,jdir)*block%dx(ixOmin1:ixOmax1,&
                       ixOmin2:ixOmax2,kdir)
                  end if
                  curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     idir)=(curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     idir)+(tmp(hxOmin1:hxOmax1,&
                     hxOmin2:hxOmax2)-tmp(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2))*block%dx(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2,kdir))/surface(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2)
                end if
               else ! idir==z_
                if(jdir<kdir) then
                  ! idir=z,jdir=r,kdir=phi
                  !! integral along r dimension
                  hxOmin1=ixOmin1-kr(kdir,1);hxOmin2=ixOmin2-kr(kdir,2)
                  hxOmax1=ixOmax1-kr(kdir,1);hxOmax2=ixOmax2-kr(kdir,2);
                  ixCmin1=hxOmin1;ixCmin2=hxOmin2;ixCmax1=ixOmax1
                  ixCmax2=ixOmax2;
                  jxCmin1=ixCmin1+kr(kdir,1);jxCmin2=ixCmin2+kr(kdir,2)
                  jxCmax1=ixCmax1+kr(kdir,1);jxCmax2=ixCmax2+kr(kdir,2);
                  ! qvec(r) at cell interface along phi dimension
                  if(stretched_dim(kdir).or.stretched_symm_dim(kdir)) then
                    tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=qvec(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,jdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,kdir)*(qvec(jxCmin1:jxCmax1,&
                       jxCmin2:jxCmax2,jdir)-qvec(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,jdir))/(block%x(jxCmin1:jxCmax1,&
                       jxCmin2:jxCmax2,kdir)-block%x(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,kdir))
                  else
                    tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=0.5d0*(qvec(&
                       ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                       jdir)+qvec(jxCmin1:jxCmax1,jxCmin2:jxCmax2,jdir))
                  end if
                  curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     idir)=(tmp(hxOmin1:hxOmax1,&
                     hxOmin2:hxOmax2)-tmp(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2))*block%dx(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2,jdir)
                  !! integral along phi dimension
                  hxOmin1=ixOmin1-kr(jdir,1);hxOmin2=ixOmin2-kr(jdir,2)
                  hxOmax1=ixOmax1-kr(jdir,1);hxOmax2=ixOmax2-kr(jdir,2);
                  ixCmin1=hxOmin1;ixCmin2=hxOmin2;ixCmax1=ixOmax1
                  ixCmax2=ixOmax2;
                  jxCmin1=ixCmin1+kr(jdir,1);jxCmin2=ixCmin2+kr(jdir,2)
                  jxCmax1=ixCmax1+kr(jdir,1);jxCmax2=ixCmax2+kr(jdir,2);
                  ! qvec(phi) at cell interface along r dimension
                  if(stretched_dim(jdir).or.stretched_symm_dim(jdir)) then
                    tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=qvec(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,kdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,jdir)*(qvec(jxCmin1:jxCmax1,&
                       jxCmin2:jxCmax2,kdir)-qvec(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,kdir))/(block%x(jxCmin1:jxCmax1,&
                       jxCmin2:jxCmax2,jdir)-block%x(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,jdir))
                  else
                    tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=0.5d0*(qvec(&
                       ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                       kdir)+qvec(jxCmin1:jxCmax1,jxCmin2:jxCmax2,kdir))
                  end if
                  ! r coordinate at cell interface along r dimension
                  xC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=block%x(ixCmin1:ixCmax1,&
                     ixCmin2:ixCmax2,jdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                     ixCmin2:ixCmax2,jdir)
                  if(ndim==3) then
                    surface(ixOmin1:ixOmax1,&
                       ixOmin2:ixOmax2)=block%surface(ixOmin1:ixOmax1,&
                       ixOmin2:ixOmax2,idir)
                  else
                    surface(ixOmin1:ixOmax1,&
                       ixOmin2:ixOmax2)=block%x(ixOmin1:ixOmax1,&
                       ixOmin2:ixOmax2,jdir)*block%dx(ixOmin1:ixOmax1,&
                       ixOmin2:ixOmax2,kdir)*block%dx(ixOmin1:ixOmax1,&
                       ixOmin2:ixOmax2,jdir)
                  end if
                  curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     idir)=(curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     idir)+(xC(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2)*tmp(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2)-xC(hxOmin1:hxOmax1,&
                     hxOmin2:hxOmax2)*tmp(hxOmin1:hxOmax1,&
                     hxOmin2:hxOmax2))*block%dx(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2,kdir))/surface(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2)
                end if
               end if
               if(idir<idirmin)idirmin=idir
              endif
            enddo; enddo; enddo;
          case default
            call mpistop('no such curl evaluator')
        end select
    end select

  end subroutine curlvector





end module mod_geometry
