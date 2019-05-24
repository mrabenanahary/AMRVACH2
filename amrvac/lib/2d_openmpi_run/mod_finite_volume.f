!> Module with finite volume methods for fluxes
module mod_finite_volume
  implicit none
  private

  public :: finite_volume
  public :: hancock
  public :: reconstruct_LR

contains

  !> The non-conservative Hancock predictor for TVDLF
  !>
  !> on entry:
  !> input available on ixI^L=ixG^L asks for output on ixO^L=ixG^L^LSUBnghostcells
  !> one entry: (predictor): wCT -- w_n        wnew -- w_n   qdt=dt/2
  !> on exit :  (predictor): wCT -- w_n        wnew -- w_n+1/2
  subroutine hancock(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,idimmin,idimmax,qtC,wCT,qt,wnew,dx1,dx2,x)
    use mod_physics
    use mod_global_parameters
    use mod_source, only: addsource2

    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2, idimmin,idimmax
    double precision, intent(in) :: qdt, qtC, qt, dx1,dx2, x(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:ndim)
    double precision, intent(inout) :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw), wnew(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nw) :: wprim,&
        wLC, wRC
    ! left and right constructed status in primitive form, needed for better performance
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nw) :: wLp,&
        wRp
    double precision :: fLC(ixImin1:ixImax1,ixImin2:ixImax2, nwflux),&
        fRC(ixImin1:ixImax1,ixImin2:ixImax2, nwflux)
    double precision :: dxinv(1:ndim)
    integer :: idim, iw, ixmin1,ixmin2,ixmax1,ixmax2, hxOmin1,hxOmin2,hxOmax1,&
       hxOmax2

    ! Expand limits in each idim direction in which fluxes are added
    ixmin1=ixOmin1;ixmin2=ixOmin2;ixmax1=ixOmax1;ixmax2=ixOmax2;
    do idim= idimmin,idimmax
       ixmin1=ixmin1-kr(idim,1);ixmin2=ixmin2-kr(idim,2)
       ixmax1=ixmax1+kr(idim,1);ixmax2=ixmax2+kr(idim,2);
    end do
    if (ixImin1>ixmin1.or.ixImin2>ixmin2.or.ixImax1<ixmax1.or.ixImax2<ixmax2) &
       call mpistop("Error in Hancock: Nonconforming input limits")

    wprim=wCT
    call phys_to_primitive(ixImin1,ixImin2,ixImax1,ixImax2,ixImin1,ixImin2,&
       ixImax1,ixImax2,wprim,x)

    dxinv(1)=-qdt/dx1;dxinv(2)=-qdt/dx2;
    do idim= idimmin,idimmax
       block%iw0=idim
       ! Calculate w_j+g_j/2 and w_j-g_j/2
       ! First copy all variables, then upwind wLC and wRC.
       ! wLC is to the left of ixO, wRC is to the right of wCT.
       hxOmin1=ixOmin1-kr(idim,1);hxOmin2=ixOmin2-kr(idim,2)
       hxOmax1=ixOmax1-kr(idim,1);hxOmax2=ixOmax2-kr(idim,2);

       wRp(hxOmin1:hxOmax1,hxOmin2:hxOmax2,1:nwflux)=wprim(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,1:nwflux)
       wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwflux)=wprim(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,1:nwflux)

       call reconstruct_LR(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
          ixOmax1,ixOmax2,hxOmin1,hxOmin2,hxOmax1,hxOmax2,idim,wprim,wprim,wLC,&
          wRC,wLp,wRp,x,.false.)

       ! Calculate the fLC and fRC fluxes
       call phys_get_flux(wRC,wRp,x,ixImin1,ixImin2,ixImax1,ixImax2,hxOmin1,&
          hxOmin2,hxOmax1,hxOmax2,idim,fRC)
       call phys_get_flux(wLC,wLp,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
          ixOmin2,ixOmax1,ixOmax2,idim,fLC)

       ! Advect w(iw)
       do iw=1,nwflux
          if (slab) then
             wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=wnew(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,iw)+dxinv(idim)* (fLC(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2, iw)-fRC(hxOmin1:hxOmax1,hxOmin2:hxOmax2, iw))
          else
             wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=wnew(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,iw)-qdt/block%dvolume(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2) *(block%surfaceC(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,idim)*fLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 iw) -block%surfaceC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
                idim)*fRC(hxOmin1:hxOmax1,hxOmin2:hxOmax2, iw))
          end if
       end do
    end do ! next idim
    block%iw0=0

    if (.not.slab.and.idimmin==1) call phys_add_source_geom(qdt,ixImin1,&
       ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,wCT,wnew,x)
    call addsource2(qdt*dble(idimmax-idimmin+1)/dble(ndim), ixImin1,ixImin2,&
       ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,1,nw,qtC,wCT,qt,wnew,x,&
       .false.)

    ! check and optionally correct unphysical values
    call phys_handle_small_values(.false.,wnew,x,ixImin1,ixImin2,ixImax1,&
       ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,'finite_volume')

  end subroutine hancock

  !> finite volume method
  subroutine finite_volume(method,qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,idimmin,idimmax, qtC,wCT,qt,wnew,wold,fC,dx1,dx2,&
     x)

    use mod_physics
    use mod_global_parameters
    use mod_tvd, only:tvdlimit2
    use mod_source, only: addsource2

    character(len=*), intent(in)                         :: method
    double precision, intent(in)                         :: qdt, qtC, qt, dx1,&
       dx2
    integer, intent(in)                                  :: ixImin1,ixImin2,&
       ixImax1,ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, idimmin,idimmax
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim),&
        intent(in) ::  x
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw)               :: wCT, wnew, wold
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nwflux,&
       1:ndim)  :: fC

    ! primitive w at cell center
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nw) :: wprim
    ! left and right constructed status in conservative form
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nw) :: wLC,&
        wRC
    ! left and right constructed status in primitive form, needed for better performance
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nw) :: wLp,&
        wRp
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
        nwflux) :: fLC, fRC
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2)      :: cmaxC
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2)      :: cminC
    double precision, dimension(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)      :: inv_volume
    double precision, dimension(1:ndim)     :: dxinv
    integer, dimension(ixImin1:ixImax1,ixImin2:ixImax2)               :: &
       patchf
    integer :: idim, iw, ixmin1,ixmin2,ixmax1,ixmax2, hxOmin1,hxOmin2,hxOmax1,&
       hxOmax2, ixCmin1,ixCmin2,ixCmax1,ixCmax2, ixCRmin1,ixCRmin2,ixCRmax1,&
       ixCRmax2, kxCmin1,kxCmin2,kxCmax1,kxCmax2, kxRmin1,kxRmin2,kxRmax1,&
       kxRmax2

    !----------------------------------------------------------------
    if (idimmax>idimmin .and. typelimited=='original')call &
       mpistop("Error in fv: Unsplit dim. and original is limited")

    fC=0.d0

    ! The flux calculation contracts by one in the idim direction it is applied.
    ! The limiter contracts the same directions by one more, so expand ixO by 2.
    ixmin1=ixOmin1;ixmin2=ixOmin2;ixmax1=ixOmax1;ixmax2=ixOmax2;
    do idim= idimmin,idimmax
       ixmin1=ixmin1-2*kr(idim,1);ixmin2=ixmin2-2*kr(idim,2)
       ixmax1=ixmax1+2*kr(idim,1);ixmax2=ixmax2+2*kr(idim,2);
    end do
    if (ixImin1>ixmin1.or.ixImin2>ixmin2.or.ixImax1<ixmax1.or.ixImax2<ixmax2) &
       call mpistop("Error in fv : Nonconforming input limits")


    wprim=wCT
    call phys_to_primitive(ixImin1,ixImin2,ixImax1,ixImax2,ixImin1,ixImin2,&
       ixImax1,ixImax2,wprim,x)


    dxinv(1)=-qdt/dx1;dxinv(2)=-qdt/dx2;
    do idim= idimmin,idimmax
       ! use interface value of w0 at idim
       block%iw0=idim

       hxOmin1=ixOmin1-kr(idim,1);hxOmin2=ixOmin2-kr(idim,2)
       hxOmax1=ixOmax1-kr(idim,1);hxOmax2=ixOmax2-kr(idim,2);
       ! ixC is centered index in the idim direction from ixOmin-1/2 to ixOmax+1/2
       ixCmax1=ixOmax1;ixCmax2=ixOmax2; ixCmin1=hxOmin1;ixCmin2=hxOmin2;

       kxCmin1=ixImin1;kxCmin2=ixImin2; kxCmax1=ixImax1-kr(idim,1)
       kxCmax2=ixImax2-kr(idim,2);
       kxRmin1=kxCmin1+kr(idim,1);kxRmin2=kxCmin2+kr(idim,2)
       kxRmax1=kxCmax1+kr(idim,1);kxRmax2=kxCmax2+kr(idim,2);

       ! wRp and wLp are defined at the same locations, and will correspond to
       ! the left and right reconstructed values at a cell face. Their indexing
       ! is similar to cell-centered values, but in direction idim they are
       ! shifted half a cell towards the 'lower' direction.
       wRp(kxCmin1:kxCmax1,kxCmin2:kxCmax2,1:nw)=wprim(kxRmin1:kxRmax1,&
          kxRmin2:kxRmax2,1:nw)
       wLp(kxCmin1:kxCmax1,kxCmin2:kxCmax2,1:nw)=wprim(kxCmin1:kxCmax1,&
          kxCmin2:kxCmax2,1:nw)

       ! Determine stencil size
       ixCRmin1 = ixCmin1 - phys_wider_stencil
       ixCRmin2 = ixCmin2 - phys_wider_stencil
       ixCRmax1 = ixCmax1 + phys_wider_stencil
       ixCRmax2 = ixCmax2 + phys_wider_stencil

       ! apply limited reconstruction for left and right status at cell interfaces
       select case (typelimited)
       case ('previous')
          call reconstruct_LR(ixImin1,ixImin2,ixImax1,ixImax2,ixCRmin1,&
             ixCRmin2,ixCRmax1,ixCRmax2,ixCRmin1,ixCRmin2,ixCRmax1,ixCRmax2,&
             idim,wold,wprim,wLC,wRC,wLp,wRp,x,.true.)
       case ('predictor')
          call reconstruct_LR(ixImin1,ixImin2,ixImax1,ixImax2,ixCRmin1,&
             ixCRmin2,ixCRmax1,ixCRmax2,ixCRmin1,ixCRmin2,ixCRmax1,ixCRmax2,&
             idim,wprim,wprim,wLC,wRC,wLp,wRp,x,.false.)
       case default
          call mpistop("Error in reconstruction: no such base for limiter")
       end select

       ! special modification of left and right status before flux evaluation
       call phys_modify_wLR(wLp, wRp, ixImin1,ixImin2,ixImax1,ixImax2, ixCmin1,&
          ixCmin2,ixCmax1,ixCmax2, idim)

       ! evaluate physical fluxes according to reconstructed status
       call phys_get_flux(wLC,wLp,x,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,&
          ixCmin2,ixCmax1,ixCmax2,idim,fLC)
       call phys_get_flux(wRC,wRp,x,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,&
          ixCmin2,ixCmax1,ixCmax2,idim,fRC)
       ! estimating bounds for the minimum and maximum signal velocities
       if(method=='tvdlf'.or.method=='tvdmu') then
         call phys_get_cbounds(wLC,wRC,wLp,wRp,x,ixImin1,ixImin2,ixImax1,&
            ixImax2,ixCmin1,ixCmin2,ixCmax1,ixCmax2,idim,cmaxC)
       else
         call phys_get_cbounds(wLC,wRC,wLp,wRp,x,ixImin1,ixImin2,ixImax1,&
            ixImax2,ixCmin1,ixCmin2,ixCmax1,ixCmax2,idim,cmaxC,cminC)
       end if

       ! use approximate Riemann solver to get flux at interfaces
       select case(method)
       case('tvdmu')
         call get_Riemann_flux_tvdmu()
       case('tvdlf')
         call get_Riemann_flux_tvdlf()
       case('hll')
         call get_Riemann_flux_hll()
       case('hllc','hllcd')
         call get_Riemann_flux_hllc()
       case('hlld')
         call get_Riemann_flux_hlld()
       case default
         call mpistop('unkown Riemann flux')
       end select
    end do ! Next idim
    block%iw0=0

    do idim= idimmin,idimmax
       hxOmin1=ixOmin1-kr(idim,1);hxOmin2=ixOmin2-kr(idim,2)
       hxOmax1=ixOmax1-kr(idim,1);hxOmax2=ixOmax2-kr(idim,2);

       ! Multiply the fluxes by -dt/dx since Flux fixing expects this
       if (slab) then
          fC(ixImin1:ixImax1,ixImin2:ixImax2,1:nwflux,&
             idim)=dxinv(idim)*fC(ixImin1:ixImax1,ixImin2:ixImax2,1:nwflux,&
             idim)
          wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwflux)=wnew(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,1:nwflux) + (fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             1:nwflux,idim)-fC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,1:nwflux,idim))
       else
          fC(ixImin1:ixImax1,ixImin2:ixImax2,1:nwflux,&
             idim)=-qdt*fC(ixImin1:ixImax1,ixImin2:ixImax2,1:nwflux,idim)
          if (.not. angmomfix) then ! default case
            inv_volume = 1.0d0/block%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
            do iw=1,nwflux
              wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=wnew(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,iw) + (fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw,&
                 idim)-fC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,iw,&
                 idim)) * inv_volume
            enddo
          else
            ! If angular momentum conserving way to solve the equations,
            ! some fluxes additions need to be treated specifically
            call phys_angmomfix(fC,x,wnew,ixImin1,ixImin2,ixImax1,ixImax2,&
               ixOmin1,ixOmin2,ixOmax1,ixOmax2,idim)
          endif
       end if

       ! For the MUSCL scheme apply the characteristic based limiter
       if (method=='tvdmu') call tvdlimit2(method,qdt,ixImin1,ixImin2,ixImax1,&
          ixImax2,ixCmin1,ixCmin2,ixCmax1,ixCmax2,ixOmin1,ixOmin2,ixOmax1,&
          ixOmax2,idim,wLC,wRC,wnew,x,fC,dx1,dx2)

    end do ! Next idim

    if (.not.slab.and.idimmin==1) call phys_add_source_geom(qdt,ixImin1,&
       ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,wCT,wnew,x)

    ! add source terms
    call addsource2(qdt*dble(idimmax-idimmin+1)/dble(ndim), ixImin1,ixImin2,&
       ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,1,nw,qtC,wCT,qt,wnew,x,&
       .false.)

    ! check and optionally correct unphysical values
    call phys_handle_small_values(.false.,wnew,x,ixImin1,ixImin2,ixImax1,&
       ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,'finite_volume')

  contains

    subroutine get_Riemann_flux_tvdmu()

      do iw=1,nwflux
         ! To save memory we use fLC to store (F_L+F_R)/2=half*(fLC+fRC)
         fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2, iw)=half*(fLC(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2, iw)+fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2, iw))
         if (slab) then
           fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,idim)=fLC(ixCmin1:ixCmax1,&
              ixCmin2:ixCmax2, iw)
         else
           fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,&
              idim)=block%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
              idim)*fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2, iw)
         end if
      end do
    end subroutine get_Riemann_flux_tvdmu

    subroutine get_Riemann_flux_tvdlf()
      double precision :: fac(ixCmin1:ixCmax1,ixCmin2:ixCmax2)

      fac = -0.5d0*tvdlfeps*cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)

      ! Calculate fLC=f(uL_j+1/2) and fRC=f(uR_j+1/2) for each iw
      do iw=1,nwflux

         ! To save memory we use fLC to store (F_L+F_R)/2=half*(fLC+fRC)
         fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2, iw)=0.5d0*(fLC(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2, iw)+fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2, iw))

         ! Add TVDLF dissipation to the flux
         if (flux_type(idim, iw) /= flux_no_dissipation) then
            fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2, iw)=fLC(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2, iw) + fac*(wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               iw)-wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw))
         end if

         if (slab) then
           fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,idim)=fLC(ixCmin1:ixCmax1,&
              ixCmin2:ixCmax2, iw)
         else
           fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,&
              idim)=block%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
              idim)*fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2, iw)
         end if

      end do ! Next iw
    end subroutine get_Riemann_flux_tvdlf

    subroutine get_Riemann_flux_hll()

      double precision :: fac(ixCmin1:ixCmax1,ixCmin2:ixCmax2),&
          div(ixCmin1:ixCmax1,ixCmin2:ixCmax2)

      where(cminC(ixCmin1:ixCmax1,ixCmin2:ixCmax2) >= zero)
        patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2) = -2
      elsewhere(cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2) <= zero)
        patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2) =  2
      elsewhere
        patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2) =  1
      endwhere

      fac = tvdlfeps*cminC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)*cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
      div = 1/(cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)-cminC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2))

      ! Calculate fLC=f(uL_j+1/2) and fRC=f(uR_j+1/2) for each iw
      do iw=1,nwflux
         if (flux_type(idim, iw) == flux_tvdlf) then
            fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                iw) = half*(fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                iw) + fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                iw) -tvdlfeps*max(cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2),&
                dabs(cminC(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2))) * (wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               iw)-wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)))
         else
            where(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2)==1)
               ! Add hll dissipation to the flux
               fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                   iw) = (cmaxC(ixCmin1:ixCmax1,&
                  ixCmin2:ixCmax2)*fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                   iw)-cminC(ixCmin1:ixCmax1,&
                  ixCmin2:ixCmax2) * fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                   iw) +fac*(wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                  iw)-wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw))) * div
            elsewhere(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2)== 2)
               fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2, iw)=fRC(ixCmin1:ixCmax1,&
                  ixCmin2:ixCmax2, iw)
            elsewhere(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2)==-2)
               fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2, iw)=fLC(ixCmin1:ixCmax1,&
                  ixCmin2:ixCmax2, iw)
            endwhere
         endif

         if (slab) then
           fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,idim)=fLC(ixCmin1:ixCmax1,&
              ixCmin2:ixCmax2, iw)
         else
           fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,&
              idim)=block%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
              idim)*fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2, iw)
         end if

      end do ! Next iw
    end subroutine get_Riemann_flux_hll

    subroutine get_Riemann_flux_hllc()
      implicit none
      double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
         1:nwflux)     :: whll, Fhll, fCD
      double precision, dimension(ixImin1:ixImax1,&
         ixImin2:ixImax2)              :: lambdaCD

      patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2) =  1
      where(cminC(ixCmin1:ixCmax1,ixCmin2:ixCmax2) >= zero)
         patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2) = -2
      elsewhere(cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2) <= zero)
         patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2) =  2
      endwhere
      ! Use more diffusive scheme, is actually TVDLF and selected by patchf=4
      if(method=='hllcd') call phys_diffuse_hllcd(ixImin1,ixImin2,ixImax1,&
         ixImax2,ixCmin1,ixCmin2,ixCmax1,ixCmax2,idim,wLC,wRC,fLC,fRC,patchf)


      ! now patchf may be -1 or 1 due to phys_get_lCD
      if(any(abs(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2))== 1))then
         !---- calculate speed lambda at CD ----!
         call phys_get_lCD(wLC,wRC,fLC,fRC,cminC,cmaxC,idim,ixImin1,ixImin2,&
            ixImax1,ixImax2,ixCmin1,ixCmin2,ixCmax1,ixCmax2, whll,Fhll,&
            lambdaCD,patchf)
         !======== flux at intermediate state ========!
         call phys_get_wCD(x,wLC,wRC,whll,fRC,fLC,Fhll,patchf,lambdaCD,cminC,&
            cmaxC,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,ixCmax1,&
            ixCmax2,idim,fCD)
      endif ! Calculate the CD flux

      do iw=1,nwflux
         if (flux_type(idim, iw) == flux_tvdlf) then
            fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               iw) = 0.5d0 * (fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               iw) + fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               iw) - tvdlfeps *  max(cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2),&
                abs(cminC(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2))) * (wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               iw) - wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)))
         else
            where(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2)==-2)
               fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=fLC(ixCmin1:ixCmax1,&
                  ixCmin2:ixCmax2,iw)
            elsewhere(abs(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2))==1)
               fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=fCD(ixCmin1:ixCmax1,&
                  ixCmin2:ixCmax2,iw)
            elsewhere(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2)==2)
               fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=fRC(ixCmin1:ixCmax1,&
                  ixCmin2:ixCmax2,iw)
            elsewhere(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2)==3)
               ! fallback option, reducing to HLL flux
               fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=Fhll(ixCmin1:ixCmax1,&
                  ixCmin2:ixCmax2,iw)
            elsewhere(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2)==4)
               ! fallback option, reducing to TVDLF flux
               fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                  iw) = half*((fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                  iw)+fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                  iw)) -tvdlfeps * max(cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2),&
                   dabs(cminC(ixCmin1:ixCmax1,&
                  ixCmin2:ixCmax2))) * (wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                  iw)-wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)))
            endwhere
         end if

         if (slab) then
           fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,idim)=fLC(ixCmin1:ixCmax1,&
              ixCmin2:ixCmax2,iw)
         else
           fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,&
              idim)=block%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
              idim)*fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)
         end if

      end do ! Next iw
    end subroutine get_Riemann_flux_hllc

    !> HLLD Riemann flux from Miyoshi 2005 JCP, 208, 315
    subroutine get_Riemann_flux_hlld()
      use mod_mhd_phys
      implicit none
      double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
         1:nwflux) :: w1R,w1L,f1R,f1L
      double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
         1:nwflux) :: w2R,w2L
      double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2) :: sm,s1R,&
         s1L,suR,suL,Bx
      double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2) :: pts,ptR,&
         ptL,signBx,r1L,r1R,tmp
      double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ndir) :: vRC,&
          vLC
      integer :: ip1,ip2,ip3,idir

      f1R=0.d0
      f1L=0.d0
      ip1=idim
      ip3=3
      vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,:)=wRp(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,mom(:))
      vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,:)=wLp(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,mom(:))
      ! estimate normal magnetic field at cell interfaces
      Bx(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=0.5d0*(wRC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,mag(ip1))+wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         mag(ip1)))
      suR(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=(cmaxC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)-vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ip1))*wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,rho_)
      suL(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=(cminC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)-vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ip1))*wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,rho_)
      ptR(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=wRp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         e_)+0.5d0*sum(wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(:))**2,&
         dim=ndim+1)
      ptL(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=wLp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         e_)+0.5d0*sum(wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(:))**2,&
         dim=ndim+1)
      ! equation (38)
      sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=(suR(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)*vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ip1)-suL(ixCmin1:ixCmax1,ixCmin2:ixCmax2)*vLC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ip1)-ptR(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)+ptL(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2))/(suR(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)-suL(ixCmin1:ixCmax1,ixCmin2:ixCmax2))
      ! equation (39)
      w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(ip1))=sm(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)
      w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(ip1))=sm(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)
      w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(ip1))=sm(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)
      w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(ip1))=sm(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)
      w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip1))=Bx(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)
      w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip1))=Bx(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)
      w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip1))=Bx(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)
      w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip1))=Bx(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)
      ! equation (41)
      pts(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=(suR(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)*ptL(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)-suL(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)*ptR(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)+suR(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)*suL(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)*(vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ip1)-vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ip1)))/(suR(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)-suL(ixCmin1:ixCmax1,ixCmin2:ixCmax2))
      ! equation (43)
      w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,rho_)=suR(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)/(cmaxC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)-sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2))
      w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,rho_)=suL(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)/(cminC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)-sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2))
      ! equation (44) ~ (47)
      ip2=mod(ip1+1,ndir)
      if(ip2==0) ip2=ndir
      r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=suR(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)*(cmaxC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)-sm(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2))-Bx(ixCmin1:ixCmax1,ixCmin2:ixCmax2)**2
      where(r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2)/=0.d0)
        r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=1.d0/r1R(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)
      endwhere
      r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=suL(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)*(cminC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)-sm(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2))-Bx(ixCmin1:ixCmax1,ixCmin2:ixCmax2)**2
      where(r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2)/=0.d0)
        r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=1.d0/r1L(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)
      endwhere
      w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(ip2))=vRC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ip2)-Bx(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)*wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         mag(ip2))*(sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2)-vRC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ip1))*r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
      w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(ip2))=vLC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ip2)-Bx(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)*wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         mag(ip2))*(sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2)-vLC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ip1))*r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
      w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip2))=(suR(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)*(cmaxC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)-vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ip1))-Bx(ixCmin1:ixCmax1,ixCmin2:ixCmax2)**2)*r1R(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)
      w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip2))=(suL(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)*(cminC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)-vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ip1))-Bx(ixCmin1:ixCmax1,ixCmin2:ixCmax2)**2)*r1L(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)
      if(ndir==3) then
        ip3=mod(ip1+2,ndir)
        if(ip3==0) ip3=ndir
        w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(ip3))=vRC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,ip3)-Bx(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           mag(ip3))*(sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2)-vRC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,ip1))*r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
        w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(ip3))=vLC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,ip3)-Bx(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           mag(ip3))*(sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2)-vLC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,ip1))*r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
        w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip3))=wRC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,mag(ip3))*w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           mag(ip2))
        w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip3))=wLC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,mag(ip3))*w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           mag(ip2))
      end if
      w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip2))=wRC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,mag(ip2))*w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         mag(ip2))
      w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip2))=wLC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,mag(ip2))*w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         mag(ip2))
      ! equation (48)
      if(mhd_energy) then
        w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,e_)=((cmaxC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)-vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ip1))*wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,e_)-ptR(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ip1)+pts(ixCmin1:ixCmax1,ixCmin2:ixCmax2)*sm(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)+Bx(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*(sum(vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           :)*wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(:)),&
           dim=ndim+1)-sum(w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           mom(:))*w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(:)),&
           dim=ndim+1)))/(cmaxC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)-sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2))
        w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,e_)=((cminC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)-vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ip1))*wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,e_)-ptL(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ip1)+pts(ixCmin1:ixCmax1,ixCmin2:ixCmax2)*sm(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)+Bx(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*(sum(vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           :)*wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(:)),&
           dim=ndim+1)-sum(w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           mom(:))*w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(:)),&
           dim=ndim+1)))/(cminC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)-sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2))
      end if
      ! equation (49)
      w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,rho_)=w1R(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,rho_)
      w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,rho_)=w1L(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,rho_)
      r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=sqrt(w1R(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,rho_))
      r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=sqrt(w1L(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,rho_))
      tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=1.d0/(r1R(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)+r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2))
      signBx(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=sign(1.d0,Bx(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2))
      ! equation (51)
      s1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=sm(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)+abs(Bx(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2))/r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
      s1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=sm(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)-abs(Bx(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2))/r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
      ! equation (59)
      w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(ip2))=(r1L(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)*w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         mom(ip2))+r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2)*w1R(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,mom(ip2))+(w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         mag(ip2))-w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         mag(ip2)))*signBx(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2))*tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
      w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(ip2))=w2R(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,mom(ip2))
      ! equation (61)
      w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip2))=(r1L(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)*w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         mag(ip2))+r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2)*w1L(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,mag(ip2))+r1L(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)*r1R(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)*(w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         mom(ip2))-w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         mom(ip2)))*signBx(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2))*tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
      w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip2))=w2R(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,mag(ip2))
      if(ndir==3) then
        ! equation (60)
        w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(ip3))=(r1L(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           mom(ip3))+r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2)*w1R(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,mom(ip3))+(w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           mag(ip3))-w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           mag(ip3)))*signBx(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2))*tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
        w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(ip3))=w2R(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,mom(ip3))
        ! equation (62)
        w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip3))=(r1L(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           mag(ip3))+r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2)*w1L(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,mag(ip3))+r1L(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*r1R(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*(w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           mom(ip3))-w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           mom(ip3)))*signBx(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2))*tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
        w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip3))=w2R(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,mag(ip3))
      end if
      ! equation (63)
      if(mhd_energy) then
        w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,e_)=w1R(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,e_)+r1R(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*(sum(w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           mom(:))*w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(:)),&
           dim=ndim+1)-sum(w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           mom(:))*w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(:)),&
           dim=ndim+1))*signBx(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
        w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,e_)=w1L(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,e_)-r1L(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*(sum(w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           mom(:))*w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(:)),&
           dim=ndim+1)-sum(w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           mom(:))*w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(:)),&
           dim=ndim+1))*signBx(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
      end if
      ! convert velocity to momentum
      do idir=1,ndir
        w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(idir))=w1R(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,mom(idir))*w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           rho_)
        w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(idir))=w1L(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,mom(idir))*w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           rho_)
        w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(idir))=w2R(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,mom(idir))*w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           rho_)
        w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(idir))=w2L(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,mom(idir))*w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           rho_)
      end do
      ! get fluxes of intermedate states
      do iw=1,nwflux
        if (flux_type(idim, iw) == flux_tvdlf) then
          fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,&
             ip1)=0.5d0*(fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             iw) + fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw))
          !fC(ixC^S,iw,ip1) = 0.5d0 * (fLC(ixC^S,iw) + fRC(ixC^S,iw) - tvdlfeps * &
          !     max(cmaxC(ixC^S), abs(cminC(ixC^S))) * &
          !     (wRC(ixC^S,iw)-wLC(ixC^S,iw)))
          cycle
        end if
        f1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=fLC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,iw)+cminC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*(w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           iw)-wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw))
        f1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=fRC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,iw)+cmaxC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*(w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           iw)-wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw))
        where(cminC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)>0.d0)
          fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,ip1)=fLC(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,iw)
        else where(s1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2)>=0.d0)
          fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,ip1)=f1L(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,iw)
        else where(sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2)>=0.d0)
          fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,ip1)=f1L(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,iw)+s1L(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2)*(w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             iw)-w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw))
        else where(s1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2)>=0.d0)
          fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,ip1)=f1R(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,iw)+s1R(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2)*(w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             iw)-w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw))
        else where(cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)>=0.d0)
          fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,ip1)=f1R(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,iw)
        else where(cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)<0.d0)
          fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,ip1)=fRC(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,iw)
        end where
        if(.not.slab) then
          fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,&
             ip1)=block%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ip1)*fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,ip1)
        end if
      end do

    end subroutine get_Riemann_flux_hlld

  end subroutine finite_volume

  !> Determine the upwinded wLC(ixL) and wRC(ixR) from w.
  !> the wCT is only used when PPM is exploited.
  subroutine reconstruct_LR(ixImin1,ixImin2,ixImax1,ixImax2,ixLmin1,ixLmin2,&
     ixLmax1,ixLmax2,ixRmin1,ixRmin2,ixRmax1,ixRmax2,idim,w,wCT,wLC,wRC,wLp,&
     wRp,x,needprim)
    use mod_physics
    use mod_global_parameters
    use mod_limiter

    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixLmin1,ixLmin2,&
       ixLmax1,ixLmax2, ixRmin1,ixRmin2,ixRmax1,ixRmax2, idim
    logical, intent(in) :: needprim
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nw) :: w,&
        wCT
    ! left and right constructed status in conservative form
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nw) :: wLC,&
        wRC
    ! left and right constructed status in primitive form
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nw) :: wLp,&
        wRp
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim) :: x

    integer            :: jxRmin1,jxRmin2,jxRmax1,jxRmax2, ixCmin1,ixCmin2,&
       ixCmax1,ixCmax2, jxCmin1,jxCmin2,jxCmax1,jxCmax2, iw
    double precision   :: ldw(ixImin1:ixImax1,ixImin2:ixImax2),&
        rdw(ixImin1:ixImax1,ixImin2:ixImax2), dwC(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    ! integer            :: flagL(ixI^S), flagR(ixI^S)

    ! Transform w,wL,wR to primitive variables
    if (needprim) then
       call phys_to_primitive(ixImin1,ixImin2,ixImax1,ixImax2,ixImin1,ixImin2,&
          ixImax1,ixImax2,w,x)
    end if

    if (typelimiter == limiter_mp5) then
       call MP5limiter(ixImin1,ixImin2,ixImax1,ixImax2,ixLmin1,ixLmin2,ixLmax1,&
          ixLmax2,idim,w,wLp,wRp)
    else if (typelimiter == limiter_ppm) then
       call PPMlimiter(ixImin1,ixImin2,ixImax1,ixImax2,ixMlo1,ixMlo2,ixMhi1,&
          ixMhi2,idim,w,wCT,wLp,wRp)
    else
       jxRmin1=ixRmin1+kr(idim,1);jxRmin2=ixRmin2+kr(idim,2)
       jxRmax1=ixRmax1+kr(idim,1);jxRmax2=ixRmax2+kr(idim,2);
       ixCmax1=jxRmax1;ixCmax2=jxRmax2; ixCmin1=ixLmin1-kr(idim,1)
       ixCmin2=ixLmin2-kr(idim,2);
       jxCmin1=ixCmin1+kr(idim,1);jxCmin2=ixCmin2+kr(idim,2)
       jxCmax1=ixCmax1+kr(idim,1);jxCmax2=ixCmax2+kr(idim,2);

       do iw=1,nwflux
          if (loglimit(iw)) then
             w(ixCmin1:jxCmax1,ixCmin2:jxCmax2,iw)=dlog10(w(ixCmin1:jxCmax1,&
                ixCmin2:jxCmax2,iw))
             wLp(ixLmin1:ixLmax1,ixLmin2:ixLmax2,&
                iw)=dlog10(wLp(ixLmin1:ixLmax1,ixLmin2:ixLmax2,iw))
             wRp(ixRmin1:ixRmax1,ixRmin2:ixRmax2,&
                iw)=dlog10(wRp(ixRmin1:ixRmax1,ixRmin2:ixRmax2,iw))
          end if

          dwC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=w(jxCmin1:jxCmax1,&
             jxCmin2:jxCmax2,iw)-w(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)

          ! limit flux from left and/or right
          call dwlimiter2(dwC,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,&
             ixCmax1,ixCmax2,idim,typelimiter,ldw,rdw)
          wLp(ixLmin1:ixLmax1,ixLmin2:ixLmax2,iw)=wLp(ixLmin1:ixLmax1,&
             ixLmin2:ixLmax2,iw)+half*ldw(ixLmin1:ixLmax1,ixLmin2:ixLmax2)
          wRp(ixRmin1:ixRmax1,ixRmin2:ixRmax2,iw)=wRp(ixRmin1:ixRmax1,&
             ixRmin2:ixRmax2,iw)-half*rdw(jxRmin1:jxRmax1,jxRmin2:jxRmax2)

          if (loglimit(iw)) then
             w(ixCmin1:jxCmax1,ixCmin2:jxCmax2,iw)=10.0d0**w(ixCmin1:jxCmax1,&
                ixCmin2:jxCmax2,iw)
             wLp(ixLmin1:ixLmax1,ixLmin2:ixLmax2,&
                iw)=10.0d0**wLp(ixLmin1:ixLmax1,ixLmin2:ixLmax2,iw)
             wRp(ixRmin1:ixRmax1,ixRmin2:ixRmax2,&
                iw)=10.0d0**wRp(ixRmin1:ixRmax1,ixRmin2:ixRmax2,iw)
          end if
       end do

       !! TODO: does this actually help? if not, remove
       !call phys_check_w(.true., ixI^L, ixL^L, wLtmp, flagL)
       !call phys_check_w(.true., ixI^L, ixR^L, wRtmp, flagR)

       !do iw=1,nwflux
       !   where (flagL(ixL^S) == 0 .and. flagR(ixR^S) == 0)
       !      wLC(ixL^S,iw)=wLtmp(ixL^S,iw)
       !      wRC(ixR^S,iw)=wRtmp(ixR^S,iw)
       !   end where

       !   ! Elsewhere, we still need to convert back when using loglimit
       !   if (loglimit(iw)) then
       !      where (flagL(ixL^S) /= 0 .or. flagR(ixR^S) /= 0)
       !         wLC(ixL^S,iw)=10.0d0**wLC(ixL^S,iw)
       !         wRC(ixR^S,iw)=10.0d0**wRC(ixR^S,iw)
       !      end where
       !   end if
       !enddo
    endif

    ! Transform w,wL,wR back to conservative variables
    if(needprim)then
       call phys_to_conserved(ixImin1,ixImin2,ixImax1,ixImax2,ixImin1,ixImin2,&
          ixImax1,ixImax2,w,x)
    endif
    wLC(ixLmin1:ixLmax1,ixLmin2:ixLmax2,1:nw)=wLp(ixLmin1:ixLmax1,&
       ixLmin2:ixLmax2,1:nw)
    wRC(ixRmin1:ixRmax1,ixRmin2:ixRmax2,1:nw)=wRp(ixRmin1:ixRmax1,&
       ixRmin2:ixRmax2,1:nw)
    call phys_to_conserved(ixImin1,ixImin2,ixImax1,ixImax2,ixLmin1,ixLmin2,&
       ixLmax1,ixLmax2,wLC,x)
    call phys_to_conserved(ixImin1,ixImin2,ixImax1,ixImax2,ixRmin1,ixRmin2,&
       ixRmax1,ixRmax2,wRC,x)

  end subroutine reconstruct_LR

end module mod_finite_volume
