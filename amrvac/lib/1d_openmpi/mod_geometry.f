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
    end select

  end subroutine set_pole

  !> Deallocate geometry-related variables
  subroutine putgridgeo(igrid)
    use mod_global_parameters
    integer, intent(in) :: igrid

    deallocate(pw(igrid)%surfaceC,pw(igrid)%surface,pw(igrid)%dvolume,&
       pw(igrid)%dx,pw(igrid)%dxcoarse,pw(igrid)%ds,pw(igrid)%dvolumecoarse)

  end subroutine putgridgeo

  subroutine fillgeo(igrid,ixGmin1,ixGmax1)
    use mod_global_parameters

    integer, intent(in) :: igrid, ixGmin1,ixGmax1

    integer :: ixmin1,ixmax1, ixCmin1,ixCmax1
    double precision :: x(ixGmin1:ixGmax1,ndim), drs(ixGmin1:ixGmax1),&
        dx2(ixGmin1:ixGmax1), dx3(ixGmin1:ixGmax1)
    !-----------------------------------------------------------------------------
    ixmin1=ixGmin1+1;ixmax1=ixGmax1-1;

    select case (typeaxial)
    case ("slabstretch")
      drs(ixGmin1:ixGmax1)=pw(igrid)%dx(ixGmin1:ixGmax1,1)
      
      

      
      ixCmin1=ixmin1-kr(1,1); ixCmax1=ixmax1;
      pw(igrid)%surfaceC(ixCmin1:ixCmax1,1)=1.d0
      pw(igrid)%surface(ixCmin1:ixCmax1,1) =1.d0
     
      
      

    case ("spherical")
      x(ixGmin1:ixGmax1,1)=pw(igrid)%x(ixGmin1:ixGmax1,1)
      

      drs(ixGmin1:ixGmax1)=pw(igrid)%dx(ixGmin1:ixGmax1,1)
      
      

      ixCmin1=ixmin1-kr(1,1); ixCmax1=ixmax1;

      pw(igrid)%surfaceC(ixCmin1:ixCmax1,1)=(x(ixCmin1:ixCmax1,&
         1)+half*drs(ixCmin1:ixCmax1))**2 

      

      

      ixCmin1=ixmin1-kr(1,1); ixCmax1=ixmax1;
      pw(igrid)%surface(ixCmin1:ixCmax1,1)=x(ixCmin1:ixCmax1,1)**2 
      

      

    case ("cylindrical")
      x(ixGmin1:ixGmax1,1)=pw(igrid)%x(ixGmin1:ixGmax1,1)
      drs(ixGmin1:ixGmax1)=pw(igrid)%dx(ixGmin1:ixGmax1,1)
      
      

      ixCmin1=ixmin1-kr(1,1); ixCmax1=ixmax1;
      pw(igrid)%surfaceC(ixCmin1:ixCmax1,1)=dabs(x(ixCmin1:ixCmax1,&
         1)+half*drs(ixCmin1:ixCmax1))
      
      

      ixCmin1=ixmin1-kr(1,1); ixCmax1=ixmax1;
      pw(igrid)%surface(ixCmin1:ixCmax1,1)=dabs(x(ixCmin1:ixCmax1,1))
      
      

    case default
      call mpistop("Sorry, typeaxial unknown")
    end select

  end subroutine fillgeo

  !> Calculate gradient of a scalar q within ixL in direction idir
  subroutine gradient(q,ixImin1,ixImax1,ixOmin1,ixOmax1,idir,gradq)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1, idir
    double precision, intent(in)    :: q(ixImin1:ixImax1)
    double precision, intent(inout) :: gradq(ixImin1:ixImax1)
    double precision                :: x(ixImin1:ixImax1,1:ndim)
    integer                         :: jxOmin1,jxOmax1, hxOmin1,hxOmax1

    x(ixImin1:ixImax1,1:ndim)=block%x(ixImin1:ixImax1,1:ndim)

    hxOmin1=ixOmin1-kr(idir,1);hxOmax1=ixOmax1-kr(idir,1);
    jxOmin1=ixOmin1+kr(idir,1);jxOmax1=ixOmax1+kr(idir,1);
    if(slab) then
      gradq(ixOmin1:ixOmax1)=half*(q(jxOmin1:jxOmax1)-&
         q(hxOmin1:hxOmax1))/dxlevel(idir)
    else
      select case(typeaxial)
      case('slabstretch')
        gradq(ixOmin1:ixOmax1)=(q(jxOmin1:jxOmax1)-&
           q(hxOmin1:hxOmax1))/(x(jxOmin1:jxOmax1,idir)-x(hxOmin1:hxOmax1,&
           idir))
      case('spherical')
        select case(idir)
        case(1)
          gradq(ixOmin1:ixOmax1)=(q(jxOmin1:jxOmax1)-&
             q(hxOmin1:hxOmax1))/((x(jxOmin1:jxOmax1,1)-x(hxOmin1:hxOmax1,1)))
          
          
        end select
      case('cylindrical')
        if(idir==phi_) then
          gradq(ixOmin1:ixOmax1)=(q(jxOmin1:jxOmax1)-&
             q(hxOmin1:hxOmax1))/((x(jxOmin1:jxOmax1,phi_)-x(hxOmin1:hxOmax1,&
             phi_))*x(ixOmin1:ixOmax1,r_))
        else
          gradq(ixOmin1:ixOmax1)=(q(jxOmin1:jxOmax1)-&
             q(hxOmin1:hxOmax1))/(x(jxOmin1:jxOmax1,idir)-x(hxOmin1:hxOmax1,&
             idir))
        end if
      case default
        call mpistop('Unknown geometry')
      end select
    end if

  end subroutine gradient

  !> Calculate gradient of a scalar q within ixL in direction idir
  !> first use limiter to go from cell center to edge
  subroutine gradientS(q,ixImin1,ixImax1,ixOmin1,ixOmax1,idir,gradq)
    use mod_global_parameters
    use mod_limiter
    use mod_ppm

    integer, intent(in)                :: ixImin1,ixImax1, ixOmin1,ixOmax1,&
        idir
    double precision, intent(in)       :: q(ixImin1:ixImax1)
    double precision, intent(inout)    :: gradq(ixImin1:ixImax1)
    double precision ,dimension(ixImin1:ixImax1) :: qC,qL,qR,dqC,ldq,rdq

    double precision :: x(ixImin1:ixImax1,1:ndim)
    double precision :: invdx
    integer          :: hxOmin1,hxOmax1,ixCmin1,ixCmax1,jxCmin1,jxCmax1,&
       gxCmin1,gxCmax1,hxCmin1,hxCmax1

    x(ixImin1:ixImax1,1:ndim)=block%x(ixImin1:ixImax1,1:ndim)

    invdx=1.d0/dxlevel(idir)
    hxOmin1=ixOmin1-kr(idir,1);hxOmax1=ixOmax1-kr(idir,1);
    ixCmin1=hxOmin1;ixCmax1=ixOmax1;
    jxCmin1=ixCmin1+kr(idir,1);jxCmax1=ixCmax1+kr(idir,1);
    gxCmin1=ixCmin1-kr(idir,1);gxCmax1=jxCmax1;
    hxCmin1=gxCmin1+kr(idir,1);hxCmax1=gxCmax1+kr(idir,1);

    ! set the gradient limiter here
    qR(gxCmin1:gxCmax1) = q(hxCmin1:hxCmax1)
    qL(gxCmin1:gxCmax1) = q(gxCmin1:gxCmax1)
    if (typegradlimiter/=limiter_ppm) then
      dqC(gxCmin1:gxCmax1)= qR(gxCmin1:gxCmax1)-qL(gxCmin1:gxCmax1)
      call dwlimiter2(dqC,ixImin1,ixImax1,gxCmin1,gxCmax1,idir,typegradlimiter,&
         ldw=ldq,rdw=rdq)
      qL(ixCmin1:ixCmax1) = qL(ixCmin1:ixCmax1) + half*ldq(ixCmin1:ixCmax1)
      qR(ixCmin1:ixCmax1) = qR(ixCmin1:ixCmax1) - half*rdq(jxCmin1:jxCmax1)
    else
      call PPMlimitervar(ixImin1,ixImax1,ixOmin1,ixOmax1,idir,q,q,qL,qR)
    endif

    if(slab) then
      gradq(ixOmin1:ixOmax1)=half*(qR(ixOmin1:ixOmax1)-&
         qL(hxOmin1:hxOmax1))*invdx
    else
      gradq(ixOmin1:ixOmax1)=(qR(ixOmin1:ixOmax1)-&
         qL(hxOmin1:hxOmax1))/block%dx(ixOmin1:ixOmax1,idir)
      select case(typeaxial)
      case('spherical')
        select case(idir)
        case(2)
          gradq(ixOmin1:ixOmax1)=gradq(ixOmin1:ixOmax1)/x(ixOmin1:ixOmax1,1)
          
        end select
      case('cylindrical')
        if(idir==phi_) gradq(ixOmin1:ixOmax1)=gradq(ixOmin1:ixOmax1)/x(&
           ixOmin1:ixOmax1,1)
      end select
    end if

  end subroutine gradientS

  !> Calculate divergence of a vector qvec within ixL
  subroutine divvector(qvec,ixImin1,ixImax1,ixOmin1,ixOmax1,divq, fourthorder)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImax1,ixOmin1,ixOmax1
    double precision, intent(in)    :: qvec(ixImin1:ixImax1,1:ndir)
    double precision, intent(inout) :: divq(ixImin1:ixImax1)
    logical, intent(in), optional   :: fourthorder !< Default: false
    logical                         :: use_4th_order
    double precision                :: qC(ixImin1:ixImax1), invdx(1:ndim)
    integer                         :: jxOmin1,jxOmax1, hxOmin1,hxOmax1,&
        ixCmin1,ixCmax1, jxCmin1,jxCmax1
    integer                         :: idims, ixmin1,ixmax1, gxOmin1,gxOmax1,&
        kxOmin1,kxOmax1

    use_4th_order = .false.
    if (present(fourthorder)) use_4th_order = fourthorder

    if (use_4th_order) then
      if (.not. slab) call mpistop&
         ("divvector: 4th order only supported for slab geometry")
      ! Fourth order, stencil width is two
      ixmin1=ixOmin1-2;ixmax1=ixOmax1+2;
    else
      ! Second order, stencil width is one
      ixmin1=ixOmin1-1;ixmax1=ixOmax1+1;
    end if

    if (ixImin1>ixmin1.or.ixImax1<ixmax1) call &
       mpistop("Error in divvector: Non-conforming input limits")

    invdx=1.d0/dxlevel
    divq(ixOmin1:ixOmax1)=0.0d0

    if (slab) then
      do idims=1,ndim
        if (.not. use_4th_order) then
          ! Use second order scheme
          jxOmin1=ixOmin1+kr(idims,1);jxOmax1=ixOmax1+kr(idims,1);
          hxOmin1=ixOmin1-kr(idims,1);hxOmax1=ixOmax1-kr(idims,1);
          divq(ixOmin1:ixOmax1)=divq(ixOmin1:ixOmax1)+&
             half*(qvec(jxOmin1:jxOmax1,idims) - qvec(hxOmin1:hxOmax1,&
             idims))*invdx(idims)
        else
          ! Use fourth order scheme
          kxOmin1=ixOmin1+2*kr(idims,1);kxOmax1=ixOmax1+2*kr(idims,1);
          jxOmin1=ixOmin1+kr(idims,1);jxOmax1=ixOmax1+kr(idims,1);
          hxOmin1=ixOmin1-kr(idims,1);hxOmax1=ixOmax1-kr(idims,1);
          gxOmin1=ixOmin1-2*kr(idims,1);gxOmax1=ixOmax1-2*kr(idims,1);
          divq(ixOmin1:ixOmax1)=divq(ixOmin1:ixOmax1)+(-qvec(kxOmin1:kxOmax1,&
             idims) + 8.0d0 * qvec(jxOmin1:jxOmax1,&
             idims) - 8.0d0 * qvec(hxOmin1:hxOmax1,&
             idims) + qvec(gxOmin1:gxOmax1,idims))/(12.0d0 * dxlevel(idims))
        end if
      end do
    else
      do idims=1,ndim
        hxOmin1=ixOmin1-kr(idims,1);hxOmax1=ixOmax1-kr(idims,1);
        ixCmin1=hxOmin1;ixCmax1=ixOmax1;
        jxCmin1=ixCmin1+kr(idims,1);jxCmax1=ixCmax1+kr(idims,1);
        if(stretched_dim(idims).or.stretched_symm_dim(idims)) then
          ! linear interpolation at cell interface along stretched dimension
          qC(ixCmin1:ixCmax1)=block%surfaceC(ixCmin1:ixCmax1,&
             idims)*(qvec(ixCmin1:ixCmax1,&
             idims)+0.5d0*block%dx(ixCmin1:ixCmax1,&
             idims)*(qvec(jxCmin1:jxCmax1,idims)-qvec(ixCmin1:ixCmax1,&
             idims))/(block%x(jxCmin1:jxCmax1,idims)-block%x(ixCmin1:ixCmax1,&
             idims)))
        else
          qC(ixCmin1:ixCmax1)=block%surfaceC(ixCmin1:ixCmax1,&
             idims)*half*(qvec(ixCmin1:ixCmax1,idims)+qvec(jxCmin1:jxCmax1,&
             idims))
        end if
        divq(ixOmin1:ixOmax1)=divq(ixOmin1:ixOmax1)+&
           qC(ixOmin1:ixOmax1)-qC(hxOmin1:hxOmax1)
      end do
      divq(ixOmin1:ixOmax1)=divq(ixOmin1:ixOmax1)/block%dvolume(&
         ixOmin1:ixOmax1)
    end if


  end subroutine divvector

  !> Calculate divergence of a vector qvec within ixL
  !> using limited extrapolation to cell edges
  subroutine divvectorS(qvec,ixImin1,ixImax1,ixOmin1,ixOmax1,divq)
    use mod_global_parameters
    use mod_limiter
    use mod_ppm

    integer, intent(in)                :: ixImin1,ixImax1,ixOmin1,ixOmax1
    double precision, intent(in)       :: qvec(ixImin1:ixImax1,1:ndir)
    double precision, intent(inout)    :: divq(ixImin1:ixImax1)
    double precision, dimension(ixImin1:ixImax1) :: qL,qR,dqC,ldq,rdq

    double precision :: invdx(1:ndim)
    integer          :: hxOmin1,hxOmax1,ixCmin1,ixCmax1,jxCmin1,jxCmax1,idims,&
       ixmin1,ixmax1,gxCmin1,gxCmax1,hxCmin1,hxCmax1

    ixmin1=ixOmin1-2;ixmax1=ixOmax1+2;

    if (ixImin1>ixmin1.or.ixImax1<ixmax1) call &
       mpistop("Error in divvectorS: Non-conforming input limits")

    invdx=1.d0/dxlevel
    divq(ixOmin1:ixOmax1)=zero
    do idims=1,ndim
      hxOmin1=ixOmin1-kr(idims,1);hxOmax1=ixOmax1-kr(idims,1);
      ixCmin1=hxOmin1;ixCmax1=ixOmax1;
      jxCmin1=ixCmin1+kr(idims,1);jxCmax1=ixCmax1+kr(idims,1);
      gxCmin1=ixCmin1-kr(idims,1);gxCmax1=jxCmax1;
      hxCmin1=gxCmin1+kr(idims,1);hxCmax1=gxCmax1+kr(idims,1);
      qR(gxCmin1:gxCmax1) = qvec(hxCmin1:hxCmax1,idims)
      qL(gxCmin1:gxCmax1) = qvec(gxCmin1:gxCmax1,idims)
      if(typegradlimiter/=limiter_ppm) then
        dqC(gxCmin1:gxCmax1)= qR(gxCmin1:gxCmax1)-qL(gxCmin1:gxCmax1)
        call dwlimiter2(dqC,ixImin1,ixImax1,gxCmin1,gxCmax1,idims,&
           typegradlimiter,ldw=ldq,rdw=rdq)
        qL(ixCmin1:ixCmax1) = qL(ixCmin1:ixCmax1) + half*ldq(ixCmin1:ixCmax1)
        qR(ixCmin1:ixCmax1) = qR(ixCmin1:ixCmax1) - half*rdq(jxCmin1:jxCmax1)
      else
        dqC(ixImin1:ixImax1)=qvec(ixImin1:ixImax1,idims)
        call PPMlimitervar(ixImin1,ixImax1,ixOmin1,ixOmax1,idims,dqC,dqC,qL,&
           qR)
      endif

      if (slab) then
        divq(ixOmin1:ixOmax1)=divq(ixOmin1:ixOmax1)+&
           half*(qR(ixOmin1:ixOmax1)-qL(hxOmin1:hxOmax1))*invdx(idims)
      else
        qR(ixCmin1:ixCmax1)=block%surfaceC(ixCmin1:ixCmax1,&
           idims)*qR(ixCmin1:ixCmax1)
        qL(ixCmin1:ixCmax1)=block%surfaceC(ixCmin1:ixCmax1,&
           idims)*qL(ixCmin1:ixCmax1)
        divq(ixOmin1:ixOmax1)=divq(ixOmin1:ixOmax1)+&
           qR(ixOmin1:ixOmax1)-qL(hxOmin1:hxOmax1)
      end if
    end do
    if(.not.slab) divq(ixOmin1:ixOmax1)=divq(ixOmin1:ixOmax1)/block%dvolume(&
       ixOmin1:ixOmax1)

  end subroutine divvectorS

  !> Calculate curl of a vector qvec within ixL
  !> Options to
  !>        employ standard second order CD evaluations
  !>        use Gauss theorem for non-Cartesian grids
  !>        use Stokes theorem for non-Cartesian grids
  subroutine curlvector(qvec,ixImin1,ixImax1,ixOmin1,ixOmax1,curlvec,idirmin,&
     idirmin0,ndir0)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImax1,ixOmin1,ixOmax1
    integer, intent(in)             :: ndir0, idirmin0
    integer, intent(inout)          :: idirmin
    double precision, intent(in)    :: qvec(ixImin1:ixImax1,1:ndir0)
    double precision, intent(inout) :: curlvec(ixImin1:ixImax1,idirmin0:3)

    integer          :: ixCmin1,ixCmax1,jxCmin1,jxCmax1,idir,jdir,kdir,hxOmin1,&
       hxOmax1,jxOmin1,jxOmax1
    double precision :: invdx(1:ndim)
    double precision :: tmp(ixImin1:ixImax1),tmp2(ixImin1:ixImax1),&
       xC(ixImin1:ixImax1),surface(ixImin1:ixImax1)
    !-----------------------------------------------------------------------------

    ! Calculate curl within ixL: CurlV_i=eps_ijk*d_j V_k
    ! Curl can have components (idirmin:3)
    ! Determine exact value of idirmin while doing the loop.

    idirmin=4
    curlvec(ixOmin1:ixOmax1,idirmin0:3)=zero

    if(slab) then ! Cartesian case
      invdx=1.d0/dxlevel
      do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir0
        if(lvc(idir,jdir,kdir)/=0)then
          tmp(ixImin1:ixImax1)=qvec(ixImin1:ixImax1,kdir)
          hxOmin1=ixOmin1-kr(jdir,1);hxOmax1=ixOmax1-kr(jdir,1);
          jxOmin1=ixOmin1+kr(jdir,1);jxOmax1=ixOmax1+kr(jdir,1);
          ! second order centered differencing
          tmp(ixOmin1:ixOmax1)=half*(tmp(jxOmin1:jxOmax1)-&
             tmp(hxOmin1:hxOmax1))*invdx(jdir)
          !> \todo allow for 4th order CD evaluation here as well
          if(lvc(idir,jdir,kdir)==1)then
            curlvec(ixOmin1:ixOmax1,idir)=curlvec(ixOmin1:ixOmax1,&
               idir)+tmp(ixOmin1:ixOmax1)
          else
            curlvec(ixOmin1:ixOmax1,idir)=curlvec(ixOmin1:ixOmax1,&
               idir)-tmp(ixOmin1:ixOmax1)
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
                tmp(ixImin1:ixImax1)=qvec(ixImin1:ixImax1,kdir)
                hxOmin1=ixOmin1-kr(jdir,1);hxOmax1=ixOmax1-kr(jdir,1);
                jxOmin1=ixOmin1+kr(jdir,1);jxOmax1=ixOmax1+kr(jdir,1);
                ! second order centered differencing
                tmp2(ixOmin1:ixOmax1)=(tmp(jxOmin1:jxOmax1)-&
                   tmp(hxOmin1:hxOmax1))/(block%x(jxOmin1:jxOmax1,&
                   jdir)-block%x(hxOmin1:hxOmax1,jdir))
              case('Gaussbased')
                hxOmin1=ixOmin1-kr(jdir,1);hxOmax1=ixOmax1-kr(jdir,1);
                ixCmin1=hxOmin1;ixCmax1=ixOmax1;
                jxCmin1=ixCmin1+kr(jdir,1);jxCmax1=ixCmax1+kr(jdir,1);
                tmp(ixCmin1:ixCmax1)=block%surfaceC(ixCmin1:ixCmax1,&
                   jdir)*(qvec(ixCmin1:ixCmax1,&
                   kdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                   jdir)*(qvec(jxCmin1:jxCmax1,kdir)-qvec(ixCmin1:ixCmax1,&
                   kdir))/(block%x(jxCmin1:jxCmax1,&
                   jdir)-block%x(ixCmin1:ixCmax1,jdir)))
                tmp2(ixOmin1:ixOmax1)=(tmp(ixOmin1:ixOmax1)-&
                   tmp(hxOmin1:hxOmax1))/block%dvolume(ixOmin1:ixOmax1)
              case('Stokesbased')
                hxOmin1=ixOmin1-kr(jdir,1);hxOmax1=ixOmax1-kr(jdir,1);
                ixCmin1=hxOmin1;ixCmax1=ixOmax1;
                jxCmin1=ixCmin1+kr(jdir,1);jxCmax1=ixCmax1+kr(jdir,1);
                if(kdir<=ndim)then
                  tmp(ixCmin1:ixCmax1)=block%ds(ixOmin1:ixOmax1,&
                     kdir)*(qvec(ixCmin1:ixCmax1,&
                     kdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                     jdir)*(qvec(jxCmin1:jxCmax1,kdir)-qvec(ixCmin1:ixCmax1,&
                     kdir))/(block%x(jxCmin1:jxCmax1,&
                     jdir)-block%x(ixCmin1:ixCmax1,jdir)))
                else
                  tmp(ixCmin1:ixCmax1)=(qvec(ixCmin1:ixCmax1,&
                     kdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                     jdir)*(qvec(jxCmin1:jxCmax1,kdir)-qvec(ixCmin1:ixCmax1,&
                     kdir))/(block%x(jxCmin1:jxCmax1,&
                     jdir)-block%x(ixCmin1:ixCmax1,jdir)))
                endif
                if(idir<=ndim)then
                  tmp2(ixOmin1:ixOmax1)=(tmp(ixOmin1:ixOmax1)-&
                     tmp(hxOmin1:hxOmax1))/block%surface(ixOmin1:ixOmax1,idir)
                else ! essentially for 2.5D case, idir=3 and jdir,kdir<=2
                  tmp2(ixOmin1:ixOmax1)=(tmp(ixOmin1:ixOmax1)-&
                     tmp(hxOmin1:hxOmax1))/(block%ds(ixOmin1:ixOmax1,&
                     jdir)*block%ds(ixOmin1:ixOmax1,kdir))
                endif
              case default
                call mpistop('no such curl evaluator')
            end select
            if(lvc(idir,jdir,kdir)==1)then
              curlvec(ixOmin1:ixOmax1,idir)=curlvec(ixOmin1:ixOmax1,&
                 idir)+tmp2(ixOmin1:ixOmax1)
            else
              curlvec(ixOmin1:ixOmax1,idir)=curlvec(ixOmin1:ixOmax1,&
                 idir)-tmp2(ixOmin1:ixOmax1)
            endif
            if(idir<idirmin)idirmin=idir
          endif
        enddo; enddo; enddo;
      case('spherical') ! possibly stretched spherical grids
        select case(typecurl)
          case('central') ! ok for any dimensionality
            do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir0
              if(lvc(idir,jdir,kdir)/=0)then
                tmp(ixImin1:ixImax1)=qvec(ixImin1:ixImax1,kdir)
                hxOmin1=ixOmin1-kr(jdir,1);hxOmax1=ixOmax1-kr(jdir,1);
                jxOmin1=ixOmin1+kr(jdir,1);jxOmax1=ixOmax1+kr(jdir,1);
                select case(jdir)
                case(1)
                tmp(ixImin1:ixImax1)=tmp(ixImin1:ixImax1)*block%x(&
                   ixImin1:ixImax1,1)
                tmp2(ixOmin1:ixOmax1)=(tmp(jxOmin1:jxOmax1)-&
                   tmp(hxOmin1:hxOmax1))/((block%x(jxOmin1:jxOmax1,&
                   1)-block%x(hxOmin1:hxOmax1,1))*block%x(ixOmin1:ixOmax1,1))
                
                
                end select
                if(lvc(idir,jdir,kdir)==1)then
                  curlvec(ixOmin1:ixOmax1,idir)=curlvec(ixOmin1:ixOmax1,&
                     idir)+tmp2(ixOmin1:ixOmax1)
                else
                  curlvec(ixOmin1:ixOmax1,idir)=curlvec(ixOmin1:ixOmax1,&
                     idir)-tmp2(ixOmin1:ixOmax1)
                endif
                if(idir<idirmin)idirmin=idir
              endif
            enddo; enddo; enddo;
          case('Gaussbased')
            do idir=idirmin0,3;
            do jdir=1,ndim; do kdir=1,ndir0
              if(lvc(idir,jdir,kdir)/=0)then
                hxOmin1=ixOmin1-kr(jdir,1);hxOmax1=ixOmax1-kr(jdir,1);
                ixCmin1=hxOmin1;ixCmax1=ixOmax1;
                jxCmin1=ixCmin1+kr(jdir,1);jxCmax1=ixCmax1+kr(jdir,1);
                tmp(ixCmin1:ixCmax1)=block%surfaceC(ixCmin1:ixCmax1,&
                   jdir)*(qvec(ixCmin1:ixCmax1,&
                   kdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                   jdir)*(qvec(jxCmin1:jxCmax1,kdir)-qvec(ixCmin1:ixCmax1,&
                   kdir))/(block%x(jxCmin1:jxCmax1,&
                   jdir)-block%x(ixCmin1:ixCmax1,jdir)))
                tmp2(ixOmin1:ixOmax1)=(tmp(ixOmin1:ixOmax1)-&
                   tmp(hxOmin1:hxOmax1))/block%dvolume(ixOmin1:ixOmax1)
                if(lvc(idir,jdir,kdir)==1)then
                  curlvec(ixOmin1:ixOmax1,idir)=curlvec(ixOmin1:ixOmax1,&
                     idir)+tmp2(ixOmin1:ixOmax1)
                else
                  curlvec(ixOmin1:ixOmax1,idir)=curlvec(ixOmin1:ixOmax1,&
                     idir)-tmp2(ixOmin1:ixOmax1)
                 endif
                 if(idir<idirmin)idirmin=idir
              endif
            enddo; enddo;
            ! geometric terms
            if(idir==2.and.phi_>0) curlvec(ixOmin1:ixOmax1,&
               2)=curlvec(ixOmin1:ixOmax1,2)+qvec(ixOmin1:ixOmax1,&
               phi_)/block%x(ixOmin1:ixOmax1,r_)
            
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
                  hxOmin1=ixOmin1-kr(jdir,1);hxOmax1=ixOmax1-kr(jdir,1);
                  ixCmin1=hxOmin1;ixCmax1=ixOmax1;
                  jxCmin1=ixCmin1+kr(jdir,1);jxCmax1=ixCmax1+kr(jdir,1);
                  ! qvec(3) at cell interface along 2nd dimension
                  if(stretched_dim(jdir).or.stretched_symm_dim(jdir)) then
                    tmp(ixCmin1:ixCmax1)=qvec(ixCmin1:ixCmax1,&
                       kdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                       jdir)*(qvec(jxCmin1:jxCmax1,kdir)-qvec(ixCmin1:ixCmax1,&
                       kdir))/(block%x(jxCmin1:jxCmax1,&
                       jdir)-block%x(ixCmin1:ixCmax1,jdir))
                  else
                    tmp(ixCmin1:ixCmax1)=0.5d0*(qvec(ixCmin1:ixCmax1,&
                       kdir)+qvec(jxCmin1:jxCmax1,kdir))
                  end if
                  ! 2nd coordinate at cell interface along 2nd dimension
                  xC(ixCmin1:ixCmax1)=block%x(ixCmin1:ixCmax1,&
                     jdir)+0.5d0*block%dx(ixCmin1:ixCmax1,jdir)
                  curlvec(ixOmin1:ixOmax1,&
                     idir)=(dsin(xC(ixOmin1:ixOmax1))*tmp(ixOmin1:ixOmax1)-&
                     dsin(xC(hxOmin1:hxOmax1))*tmp(hxOmin1:hxOmax1))*block%dx(&
                     ixOmin1:ixOmax1,kdir)
                  !! integral along 2nd dimension
                  hxOmin1=ixOmin1-kr(kdir,1);hxOmax1=ixOmax1-kr(kdir,1);
                  ixCmin1=hxOmin1;ixCmax1=ixOmax1;
                  jxCmin1=ixCmin1+kr(kdir,1);jxCmax1=ixCmax1+kr(kdir,1);
                  ! qvec(2) at cell interface along 3rd dimension
                  if(stretched_dim(kdir).or.stretched_symm_dim(kdir)) then
                    tmp(ixCmin1:ixCmax1)=qvec(ixCmin1:ixCmax1,&
                       jdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                       kdir)*(qvec(jxCmin1:jxCmax1,jdir)-qvec(ixCmin1:ixCmax1,&
                       jdir))/(block%x(jxCmin1:jxCmax1,&
                       kdir)-block%x(ixCmin1:ixCmax1,kdir))
                  else
                    tmp(ixCmin1:ixCmax1)=0.5d0*(qvec(ixCmin1:ixCmax1,&
                       jdir)+qvec(jxCmin1:jxCmax1,jdir))
                  end if
                  curlvec(ixOmin1:ixOmax1,idir)=(curlvec(ixOmin1:ixOmax1,&
                     idir)+(tmp(hxOmin1:hxOmax1)-&
                     tmp(ixOmin1:ixOmax1))*block%dx(ixOmin1:ixOmax1,&
                     jdir))/block%surface(ixOmin1:ixOmax1,&
                     idir)*block%x(ixOmin1:ixOmax1,idir)
                end if
              case(2)
                if(jdir<kdir) then
                  ! idir=2,jdir=1,kdir=3
                  !! integral along 1st dimension
                  hxOmin1=ixOmin1-kr(kdir,1);hxOmax1=ixOmax1-kr(kdir,1);
                  ixCmin1=hxOmin1;ixCmax1=ixOmax1;
                  jxCmin1=ixCmin1+kr(kdir,1);jxCmax1=ixCmax1+kr(kdir,1);
                  ! qvec(1) at cell interface along 3rd dimension
                  if(stretched_dim(kdir).or.stretched_symm_dim(kdir)) then
                    tmp(ixCmin1:ixCmax1)=qvec(ixCmin1:ixCmax1,&
                       jdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                       kdir)*(qvec(jxCmin1:jxCmax1,jdir)-qvec(ixCmin1:ixCmax1,&
                       jdir))/(block%x(jxCmin1:jxCmax1,&
                       kdir)-block%x(ixCmin1:ixCmax1,kdir))
                  else
                    tmp(ixCmin1:ixCmax1)=0.5d0*(qvec(ixCmin1:ixCmax1,&
                       jdir)+qvec(jxCmin1:jxCmax1,jdir))
                  end if
                  curlvec(ixOmin1:ixOmax1,&
                     idir)=(tmp(ixOmin1:ixOmax1)-&
                     tmp(hxOmin1:hxOmax1))*block%dx(ixOmin1:ixOmax1,1)
                  !! integral along 3rd dimension
                  hxOmin1=ixOmin1-kr(jdir,1);hxOmax1=ixOmax1-kr(jdir,1);
                  ixCmin1=hxOmin1;ixCmax1=ixOmax1;
                  jxCmin1=ixCmin1+kr(jdir,1);jxCmax1=ixCmax1+kr(jdir,1);
                  ! qvec(3) at cell interface along 1st dimension
                  if(stretched_dim(jdir).or.stretched_symm_dim(jdir)) then
                    tmp(ixCmin1:ixCmax1)=qvec(ixCmin1:ixCmax1,&
                       kdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                       jdir)*(qvec(jxCmin1:jxCmax1,kdir)-qvec(ixCmin1:ixCmax1,&
                       kdir))/(block%x(jxCmin1:jxCmax1,&
                       jdir)-block%x(ixCmin1:ixCmax1,jdir))
                  else
                    tmp(ixCmin1:ixCmax1)=0.5d0*(qvec(ixCmin1:ixCmax1,&
                       kdir)+qvec(jxCmin1:jxCmax1,kdir))
                  end if
                  ! 1st coordinate at cell interface along 1st dimension
                  xC(ixCmin1:ixCmax1)=block%x(ixCmin1:ixCmax1,&
                     jdir)+0.5d0*block%dx(ixCmin1:ixCmax1,jdir)
                  curlvec(ixOmin1:ixOmax1,idir)=(curlvec(ixOmin1:ixOmax1,&
                     idir)+(xC(hxOmin1:hxOmax1)*tmp(hxOmin1:hxOmax1)-&
                     xC(ixOmin1:ixOmax1)*tmp(ixOmin1:ixOmax1))*dsin(block%x(&
                     ixOmin1:ixOmax1,idir))*block%dx(ixOmin1:ixOmax1,&
                     kdir))/block%surface(ixOmin1:ixOmax1,idir)
                end if
              case(3)
                if(jdir<kdir) then
                  ! idir=3,jdir=1,kdir=2
                  !! integral along 1st dimension
                  hxOmin1=ixOmin1-kr(kdir,1);hxOmax1=ixOmax1-kr(kdir,1);
                  ixCmin1=hxOmin1;ixCmax1=ixOmax1;
                  jxCmin1=ixCmin1+kr(kdir,1);jxCmax1=ixCmax1+kr(kdir,1);
                  ! qvec(1) at cell interface along 2nd dimension
                  if(stretched_dim(kdir).or.stretched_symm_dim(kdir)) then
                    tmp(ixCmin1:ixCmax1)=qvec(ixCmin1:ixCmax1,&
                       jdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                       kdir)*(qvec(jxCmin1:jxCmax1,jdir)-qvec(ixCmin1:ixCmax1,&
                       jdir))/(block%x(jxCmin1:jxCmax1,&
                       kdir)-block%x(ixCmin1:ixCmax1,kdir))
                  else
                    tmp(ixCmin1:ixCmax1)=0.5d0*(qvec(ixCmin1:ixCmax1,&
                       jdir)+qvec(jxCmin1:jxCmax1,jdir))
                  end if
                  curlvec(ixOmin1:ixOmax1,&
                     idir)=(tmp(hxOmin1:hxOmax1)-&
                     tmp(ixOmin1:ixOmax1))*block%dx(ixOmin1:ixOmax1,jdir)
                  !! integral along 2nd dimension
                  hxOmin1=ixOmin1-kr(jdir,1);hxOmax1=ixOmax1-kr(jdir,1);
                  ixCmin1=hxOmin1;ixCmax1=ixOmax1;
                  jxCmin1=ixCmin1+kr(jdir,1);jxCmax1=ixCmax1+kr(jdir,1);
                  ! qvec(2) at cell interface along 1st dimension
                  if(stretched_dim(jdir).or.stretched_symm_dim(jdir)) then
                    tmp(ixCmin1:ixCmax1)=qvec(ixCmin1:ixCmax1,&
                       kdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                       jdir)*(qvec(jxCmin1:jxCmax1,kdir)-qvec(ixCmin1:ixCmax1,&
                       kdir))/(block%x(jxCmin1:jxCmax1,&
                       jdir)-block%x(ixCmin1:ixCmax1,jdir))
                  else
                    tmp(ixCmin1:ixCmax1)=0.5d0*(qvec(ixCmin1:ixCmax1,&
                       kdir)+qvec(jxCmin1:jxCmax1,kdir))
                  end if
                  ! 1st coordinate at cell interface along 1st dimension
                  xC(ixCmin1:ixCmax1)=block%x(ixCmin1:ixCmax1,&
                     jdir)+0.5d0*block%dx(ixCmin1:ixCmax1,jdir)
                  if(ndim==3) then
                    surface(ixOmin1:ixOmax1)=block%surface(ixOmin1:ixOmax1,&
                       idir)
                  else
                    surface(ixOmin1:ixOmax1)=block%x(ixOmin1:ixOmax1,&
                       jdir)*block%dx(ixOmin1:ixOmax1,&
                       kdir)*block%dx(ixOmin1:ixOmax1,jdir)
                  end if
                  curlvec(ixOmin1:ixOmax1,idir)=(curlvec(ixOmin1:ixOmax1,&
                     idir)+(xC(ixOmin1:ixOmax1)*tmp(ixOmin1:ixOmax1)-&
                     xC(hxOmin1:hxOmax1)*tmp(hxOmin1:hxOmax1))*block%dx(&
                     ixOmin1:ixOmax1,kdir))/surface(ixOmin1:ixOmax1)
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
                tmp(ixImin1:ixImax1)=qvec(ixImin1:ixImax1,kdir)
                hxOmin1=ixOmin1-kr(jdir,1);hxOmax1=ixOmax1-kr(jdir,1);
                jxOmin1=ixOmin1+kr(jdir,1);jxOmax1=ixOmax1+kr(jdir,1);
                if(z_==3.or.z_==-1) then
                  ! Case Polar_2D, Polar_2.5D or Polar_3D, i.e. R,phi,Z
                  select case(jdir)
                  case(1)
                  if(idir==z_) tmp(ixImin1:ixImax1)=tmp(&
                     ixImin1:ixImax1)*block%x(ixImin1:ixImax1,1) !R V_phi
                  ! computes d(R V_phi)/dR or d V_Z/dR
                  tmp2(ixOmin1:ixOmax1)=(tmp(jxOmin1:jxOmax1)-&
                     tmp(hxOmin1:hxOmax1))/(block%x(jxOmin1:jxOmax1,&
                     1)-block%x(hxOmin1:hxOmax1,1))
                  if(idir==z_) tmp2(ixOmin1:ixOmax1)=tmp2(&
                     ixOmin1:ixOmax1)/block%x(ixOmin1:ixOmax1,1) !(1/R)*d(R V_phi)/dR
                  
                  
                  end select
                end if
                if(phi_==3.or.phi_==-1) then
                  ! Case Cylindrical_2D, Cylindrical_2.5D or Cylindrical_3D, i.e. R,Z,phi
                  select case(jdir)
                  case(1)
                  if(idir==z_) tmp(ixImin1:ixImax1)=tmp(&
                     ixImin1:ixImax1)*block%x(ixImin1:ixImax1,1) !R V_phi
                  ! computes d(R V_phi)/dR or d V_Z/dR
                  tmp2(ixOmin1:ixOmax1)=(tmp(jxOmin1:jxOmax1)-&
                     tmp(hxOmin1:hxOmax1))/(block%x(jxOmin1:jxOmax1,&
                     1)-block%x(hxOmin1:hxOmax1,1))
                  if(idir==z_) tmp2(ixOmin1:ixOmax1)=tmp2(&
                     ixOmin1:ixOmax1)/block%x(ixOmin1:ixOmax1,1) !(1/R)*d(R V_phi)/dR
                  
                  
                  end select
                end if
                if(lvc(idir,jdir,kdir)==1)then
                  curlvec(ixOmin1:ixOmax1,idir)=curlvec(ixOmin1:ixOmax1,&
                     idir)+tmp2(ixOmin1:ixOmax1)
                else
                  curlvec(ixOmin1:ixOmax1,idir)=curlvec(ixOmin1:ixOmax1,&
                     idir)-tmp2(ixOmin1:ixOmax1)
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
                hxOmin1=ixOmin1-kr(jdir,1);hxOmax1=ixOmax1-kr(jdir,1);
                ixCmin1=hxOmin1;ixCmax1=ixOmax1;
                jxCmin1=ixCmin1+kr(jdir,1);jxCmax1=ixCmax1+kr(jdir,1);
                tmp(ixCmin1:ixCmax1)=block%surfaceC(ixCmin1:ixCmax1,&
                   jdir)*(qvec(ixCmin1:ixCmax1,&
                   kdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                   jdir)*(qvec(jxCmin1:jxCmax1,kdir)-qvec(ixCmin1:ixCmax1,&
                   kdir))/(block%x(jxCmin1:jxCmax1,&
                   jdir)-block%x(ixCmin1:ixCmax1,jdir)))
                tmp2(ixOmin1:ixOmax1)=(tmp(ixOmin1:ixOmax1)-&
                   tmp(hxOmin1:hxOmax1))/block%dvolume(ixOmin1:ixOmax1)
                if(lvc(idir,jdir,kdir)==1)then
                  curlvec(ixOmin1:ixOmax1,idir)=curlvec(ixOmin1:ixOmax1,&
                     idir)+tmp2(ixOmin1:ixOmax1)
                else
                  curlvec(ixOmin1:ixOmax1,idir)=curlvec(ixOmin1:ixOmax1,&
                     idir)-tmp2(ixOmin1:ixOmax1)
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
              if(  z_==2) curlvec(ixOmin1:ixOmax1,&
                 idir)=curlvec(ixOmin1:ixOmax1,idir)-qvec(ixOmin1:ixOmax1,&
                 z_)/block%x(ixOmin1:ixOmax1,r_)
              ! polar
              if(phi_==2) curlvec(ixOmin1:ixOmax1,&
                 idir)=curlvec(ixOmin1:ixOmax1,idir)+qvec(ixOmin1:ixOmax1,&
                 z_)/block%x(ixOmin1:ixOmax1,r_)
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
                  hxOmin1=ixOmin1-kr(jdir,1);hxOmax1=ixOmax1-kr(jdir,1);
                  ixCmin1=hxOmin1;ixCmax1=ixOmax1;
                  jxCmin1=ixCmin1+kr(jdir,1);jxCmax1=ixCmax1+kr(jdir,1);
                  ! qvec(z) at cell interface along phi dimension
                  if(stretched_dim(jdir).or.stretched_symm_dim(jdir)) then
                    tmp(ixCmin1:ixCmax1)=qvec(ixCmin1:ixCmax1,&
                       kdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                       jdir)*(qvec(jxCmin1:jxCmax1,kdir)-qvec(ixCmin1:ixCmax1,&
                       kdir))/(block%x(jxCmin1:jxCmax1,&
                       jdir)-block%x(ixCmin1:ixCmax1,jdir))
                  else
                    tmp(ixCmin1:ixCmax1)=0.5d0*(qvec(ixCmin1:ixCmax1,&
                       kdir)+qvec(jxCmin1:jxCmax1,kdir))
                  end if
                  curlvec(ixOmin1:ixOmax1,&
                     idir)=(tmp(ixOmin1:ixOmax1)-&
                     tmp(hxOmin1:hxOmax1))*block%dx(ixOmin1:ixOmax1,kdir)
                  !! integral along phi dimension
                  hxOmin1=ixOmin1-kr(kdir,1);hxOmax1=ixOmax1-kr(kdir,1);
                  ixCmin1=hxOmin1;ixCmax1=ixOmax1;
                  jxCmin1=ixCmin1+kr(kdir,1);jxCmax1=ixCmax1+kr(kdir,1);
                  ! qvec(phi) at cell interface along z dimension
                  if(stretched_dim(kdir).or.stretched_symm_dim(kdir)) then
                    tmp(ixCmin1:ixCmax1)=qvec(ixCmin1:ixCmax1,&
                       jdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                       kdir)*(qvec(jxCmin1:jxCmax1,jdir)-qvec(ixCmin1:ixCmax1,&
                       jdir))/(block%x(jxCmin1:jxCmax1,&
                       kdir)-block%x(ixCmin1:ixCmax1,kdir))
                  else
                    tmp(ixCmin1:ixCmax1)=0.5d0*(qvec(ixCmin1:ixCmax1,&
                       jdir)+qvec(jxCmin1:jxCmax1,jdir))
                  end if
                  curlvec(ixOmin1:ixOmax1,idir)=(curlvec(ixOmin1:ixOmax1,&
                     idir)+(tmp(hxOmin1:hxOmax1)-&
                     tmp(ixOmin1:ixOmax1))*block%x(ixOmin1:ixOmax1,&
                     idir)*block%dx(ixOmin1:ixOmax1,&
                     jdir))/block%surface(ixOmin1:ixOmax1,idir)
                end if
               else if(idir==phi_) then
                if(jdir<kdir) then
                  ! idir=phi,jdir=r,kdir=z
                  !! integral along r dimension
                  hxOmin1=ixOmin1-kr(kdir,1);hxOmax1=ixOmax1-kr(kdir,1);
                  ixCmin1=hxOmin1;ixCmax1=ixOmax1;
                  jxCmin1=ixCmin1+kr(kdir,1);jxCmax1=ixCmax1+kr(kdir,1);
                  ! qvec(r) at cell interface along z dimension
                  if(stretched_dim(kdir).or.stretched_symm_dim(kdir)) then
                    tmp(ixCmin1:ixCmax1)=qvec(ixCmin1:ixCmax1,&
                       jdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                       kdir)*(qvec(jxCmin1:jxCmax1,jdir)-qvec(ixCmin1:ixCmax1,&
                       jdir))/(block%x(jxCmin1:jxCmax1,&
                       kdir)-block%x(ixCmin1:ixCmax1,kdir))
                  else
                    tmp(ixCmin1:ixCmax1)=0.5d0*(qvec(ixCmin1:ixCmax1,&
                       jdir)+qvec(jxCmin1:jxCmax1,jdir))
                  end if
                  curlvec(ixOmin1:ixOmax1,&
                     idir)=(tmp(ixOmin1:ixOmax1)-&
                     tmp(hxOmin1:hxOmax1))*block%dx(ixOmin1:ixOmax1,jdir)
                  !! integral along z dimension
                  hxOmin1=ixOmin1-kr(jdir,1);hxOmax1=ixOmax1-kr(jdir,1);
                  ixCmin1=hxOmin1;ixCmax1=ixOmax1;
                  jxCmin1=ixCmin1+kr(jdir,1);jxCmax1=ixCmax1+kr(jdir,1);
                  ! qvec(z) at cell interface along r dimension
                  if(stretched_dim(jdir).or.stretched_symm_dim(jdir)) then
                    tmp(ixCmin1:ixCmax1)=qvec(ixCmin1:ixCmax1,&
                       kdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                       jdir)*(qvec(jxCmin1:jxCmax1,kdir)-qvec(ixCmin1:ixCmax1,&
                       kdir))/(block%x(jxCmin1:jxCmax1,&
                       jdir)-block%x(ixCmin1:ixCmax1,jdir))
                  else
                    tmp(ixCmin1:ixCmax1)=0.5d0*(qvec(ixCmin1:ixCmax1,&
                       kdir)+qvec(jxCmin1:jxCmax1,kdir))
                  end if
                  if(ndim==3) then
                    surface(ixOmin1:ixOmax1)=block%surface(ixOmin1:ixOmax1,&
                       idir)
                  else
                    surface(ixOmin1:ixOmax1)=block%dx(ixOmin1:ixOmax1,&
                       jdir)*block%dx(ixOmin1:ixOmax1,kdir)
                  end if
                  curlvec(ixOmin1:ixOmax1,idir)=(curlvec(ixOmin1:ixOmax1,&
                     idir)+(tmp(hxOmin1:hxOmax1)-&
                     tmp(ixOmin1:ixOmax1))*block%dx(ixOmin1:ixOmax1,&
                     kdir))/surface(ixOmin1:ixOmax1)
                end if
               else ! idir==z_
                if(jdir<kdir) then
                  ! idir=z,jdir=r,kdir=phi
                  !! integral along r dimension
                  hxOmin1=ixOmin1-kr(kdir,1);hxOmax1=ixOmax1-kr(kdir,1);
                  ixCmin1=hxOmin1;ixCmax1=ixOmax1;
                  jxCmin1=ixCmin1+kr(kdir,1);jxCmax1=ixCmax1+kr(kdir,1);
                  ! qvec(r) at cell interface along phi dimension
                  if(stretched_dim(kdir).or.stretched_symm_dim(kdir)) then
                    tmp(ixCmin1:ixCmax1)=qvec(ixCmin1:ixCmax1,&
                       jdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                       kdir)*(qvec(jxCmin1:jxCmax1,jdir)-qvec(ixCmin1:ixCmax1,&
                       jdir))/(block%x(jxCmin1:jxCmax1,&
                       kdir)-block%x(ixCmin1:ixCmax1,kdir))
                  else
                    tmp(ixCmin1:ixCmax1)=0.5d0*(qvec(ixCmin1:ixCmax1,&
                       jdir)+qvec(jxCmin1:jxCmax1,jdir))
                  end if
                  curlvec(ixOmin1:ixOmax1,&
                     idir)=(tmp(hxOmin1:hxOmax1)-&
                     tmp(ixOmin1:ixOmax1))*block%dx(ixOmin1:ixOmax1,jdir)
                  !! integral along phi dimension
                  hxOmin1=ixOmin1-kr(jdir,1);hxOmax1=ixOmax1-kr(jdir,1);
                  ixCmin1=hxOmin1;ixCmax1=ixOmax1;
                  jxCmin1=ixCmin1+kr(jdir,1);jxCmax1=ixCmax1+kr(jdir,1);
                  ! qvec(phi) at cell interface along r dimension
                  if(stretched_dim(jdir).or.stretched_symm_dim(jdir)) then
                    tmp(ixCmin1:ixCmax1)=qvec(ixCmin1:ixCmax1,&
                       kdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                       jdir)*(qvec(jxCmin1:jxCmax1,kdir)-qvec(ixCmin1:ixCmax1,&
                       kdir))/(block%x(jxCmin1:jxCmax1,&
                       jdir)-block%x(ixCmin1:ixCmax1,jdir))
                  else
                    tmp(ixCmin1:ixCmax1)=0.5d0*(qvec(ixCmin1:ixCmax1,&
                       kdir)+qvec(jxCmin1:jxCmax1,kdir))
                  end if
                  ! r coordinate at cell interface along r dimension
                  xC(ixCmin1:ixCmax1)=block%x(ixCmin1:ixCmax1,&
                     jdir)+0.5d0*block%dx(ixCmin1:ixCmax1,jdir)
                  if(ndim==3) then
                    surface(ixOmin1:ixOmax1)=block%surface(ixOmin1:ixOmax1,&
                       idir)
                  else
                    surface(ixOmin1:ixOmax1)=block%x(ixOmin1:ixOmax1,&
                       jdir)*block%dx(ixOmin1:ixOmax1,&
                       kdir)*block%dx(ixOmin1:ixOmax1,jdir)
                  end if
                  curlvec(ixOmin1:ixOmax1,idir)=(curlvec(ixOmin1:ixOmax1,&
                     idir)+(xC(ixOmin1:ixOmax1)*tmp(ixOmin1:ixOmax1)-&
                     xC(hxOmin1:hxOmax1)*tmp(hxOmin1:hxOmax1))*block%dx(&
                     ixOmin1:ixOmax1,kdir))/surface(ixOmin1:ixOmax1)
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
