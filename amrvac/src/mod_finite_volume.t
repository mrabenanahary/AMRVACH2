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
  subroutine hancock(qdt,ixI^L,ixO^L,idim^LIM,qtC,wCT,qt,wnew,dx^D,x)
    use mod_physics
    use mod_global_parameters
    use mod_source, only: addsource2

    integer, intent(in) :: ixI^L, ixO^L, idim^LIM
    real(kind=dp)   , intent(in) :: qdt, qtC, qt, dx^D, x(ixI^S,1:ndim)
    real(kind=dp)   , intent(inout) :: wCT(ixI^S,1:nw), wnew(ixI^S,1:nw)

    real(kind=dp)   , dimension(ixI^S,1:nw) :: wprim, wLC, wRC
    ! left and right constructed status in primitive form, needed for better performance
    real(kind=dp)   , dimension(ixI^S,1:nw) :: wLp, wRp
    real(kind=dp)    :: fLC(ixI^S, nwflux), fRC(ixI^S, nwflux)
    real(kind=dp)    :: dxinv(1:ndim)
    integer :: idim, iw, ix^L, hxO^L

    ! Expand limits in each idim direction in which fluxes are added
    ix^L=ixO^L;
    do idim= idim^LIM
       ix^L=ix^L^LADDkr(idim,^D);
    end do
    if (ixI^L^LTix^L|.or.|.or.) &
         call mpistop("Error in Hancock: Nonconforming input limits")

    wprim=wCT
    call phys_to_primitive(ixI^L,ixI^L,wprim,x)

    ^D&dxinv(^D)=-qdt/dx^D;
    do idim= idim^LIM
       block%iw0=idim
       ! Calculate w_j+g_j/2 and w_j-g_j/2
       ! First copy all variables, then upwind wLC and wRC.
       ! wLC is to the left of ixO, wRC is to the right of wCT.
       hxO^L=ixO^L-kr(idim,^D);

       wRp(hxO^S,1:nwflux)=wprim(ixO^S,1:nwflux)
       wLp(ixO^S,1:nwflux)=wprim(ixO^S,1:nwflux)

       call reconstruct_LR(ixI^L,ixO^L,hxO^L,idim,wprim,wprim,wLC,wRC,wLp,wRp,x,.false.)

       ! Calculate the fLC and fRC fluxes
       call phys_get_flux(wRC,wRp,x,ixI^L,hxO^L,idim,fRC)
       call phys_get_flux(wLC,wLp,x,ixI^L,ixO^L,idim,fLC)

       ! Advect w(iw)
       do iw=1,nwflux
          if (slab) then
             wnew(ixO^S,iw)=wnew(ixO^S,iw)+dxinv(idim)* &
                  (fLC(ixO^S, iw)-fRC(hxO^S, iw))
          else
             wnew(ixO^S,iw)=wnew(ixO^S,iw)-qdt/block%dvolume(ixO^S) &
                  *(block%surfaceC(ixO^S,idim)*fLC(ixO^S, iw) &
                  -block%surfaceC(hxO^S,idim)*fRC(hxO^S, iw))
          end if
       end do
    end do ! next idim
    block%iw0=0

    if (.not.slab.and.idimmin==1) call phys_add_source_geom(qdt,ixI^L,ixO^L,wCT,wnew,x)
    call addsource2(qdt*dble(idimmax-idimmin+1)/dble(ndim), &
         ixI^L,ixO^L,1,nw,qtC,wCT,qt,wnew,x,.false.)

    ! check and optionally correct unphysical values
    call phys_handle_small_values(.false.,wnew,x,ixI^L,ixO^L,'finite_volume')

  end subroutine hancock

  !> finite volume method
  subroutine finite_volume(method,qdt,ixI^L,ixO^L,idim^LIM, &
       qtC,wCT,qt,wnew,wold,fC,dx^D,x)

    use mod_physics
    use mod_global_parameters
    use mod_tvd, only:tvdlimit2
    use mod_source, only: addsource2
    use mod_usr_methods, only: usr_reset_solver
    character(len=*), intent(in)                         :: method
    real(kind=dp)   , intent(in)                         :: qdt, qtC, qt, dx^D
    integer, intent(in)                                  :: ixI^L, ixO^L, idim^LIM
    real(kind=dp)   , dimension(ixI^S,1:ndim), intent(in) ::  x
    real(kind=dp)   , dimension(ixI^S,1:nw)               :: wCT, wnew, wold
    real(kind=dp)   , dimension(ixI^S,1:nwflux,1:ndim)  :: fC

    ! primitive w at cell center
    real(kind=dp)   , dimension(ixI^S,1:nw) :: wprim
    ! left and right constructed status in conservative form
    real(kind=dp)   , dimension(ixI^S,1:nw) :: wLC, wRC
    ! left and right constructed status in primitive form, needed for better performance
    real(kind=dp)   , dimension(ixI^S,1:nw) :: wLp, wRp
    real(kind=dp)   , dimension(ixI^S, nwflux) :: fLC, fRC
    real(kind=dp)   , dimension(ixI^S)      :: cmaxC
    real(kind=dp)   , dimension(ixI^S)      :: cminC
    real(kind=dp)   , dimension(ixO^S)      :: inv_volume
    real(kind=dp)   , dimension(1:ndim)     :: dxinv
    integer, dimension(ixI^S)               :: patchf
    integer :: idim, iw, ix^L, hxO^L, ixC^L, ixCR^L, kxC^L, kxR^L
    character(len=len(method))              :: method_loc
    !----------------------------------------------------------------
    if (idimmax>idimmin .and. typelimited=='original')&
         call mpistop("Error in fv: Unsplit dim. and original is limited")

    fC=0.d0

    ! The flux calculation contracts by one in the idim direction it is applied.
    ! The limiter contracts the same directions by one more, so expand ixO by 2.
    ix^L=ixO^L;
    do idim= idim^LIM
       ix^L=ix^L^LADD2*kr(idim,^D);
    end do
    if (ixI^L^LTix^L|.or.|.or.) &
         call mpistop("Error in fv : Nonconforming input limits")


    wprim=wCT
    call phys_to_primitive(ixI^L,ixI^L,wprim,x)


    ^D&dxinv(^D)=-qdt/dx^D;
    ! initialise the flux scheme
    method_loc=trim(method)

    Loop_idims_flux: do idim= idim^LIM
       !allow the user to change localy the numerical method
       !Mialy : this should be at the end of the finite volume treatment
       !after cooling source term and then L1 is being properly computed
       !so that solver is changed for the next time step
       call usr_reset_solver(ixI^L,ixI^L,idim,qt,wprim,x,method,&
                type_limiter(node(plevel_,saveigrid)),method_loc,typelimiter)
       ! use interface value of w0 at idim
       block%iw0=idim

       hxO^L=ixO^L-kr(idim,^D);
       ! ixC is centered index in the idim direction from ixOmin-1/2 to ixOmax+1/2
       ixCmax^D=ixOmax^D; ixCmin^D=hxOmin^D;

       kxCmin^D=ixImin^D; kxCmax^D=ixImax^D-kr(idim,^D);
       kxR^L=kxC^L+kr(idim,^D);

       ! wRp and wLp are defined at the same locations, and will correspond to
       ! the left and right reconstructed values at a cell face. Their indexing
       ! is similar to cell-centered values, but in direction idim they are
       ! shifted half a cell towards the 'lower' direction.
       wRp(kxC^S,1:nw)=wprim(kxR^S,1:nw)
       wLp(kxC^S,1:nw)=wprim(kxC^S,1:nw)

       ! Determine stencil size
       {ixCRmin^D = ixCmin^D - phys_wider_stencil\}
       {ixCRmax^D = ixCmax^D + phys_wider_stencil\}

       ! apply limited reconstruction for left and right status at cell interfaces
       select case (typelimited)
       case ('previous')
          call reconstruct_LR(ixI^L,ixCR^L,ixCR^L,idim,wold,wprim,wLC,wRC,wLp,wRp,x,.true.)
       case ('predictor')
          call reconstruct_LR(ixI^L,ixCR^L,ixCR^L,idim,wprim,wprim,wLC,wRC,wLp,wRp,x,.false.)
       case default
          call mpistop("Error in reconstruction: no such base for limiter")
       end select

       ! special modification of left and right status before flux evaluation
       call phys_modify_wLR(wLp, wRp, ixI^L, ixC^L, idim)

       ! evaluate physical fluxes according to reconstructed status
       call phys_get_flux(wLC,wLp,x,ixI^L,ixC^L,idim,fLC)
       call phys_get_flux(wRC,wRp,x,ixI^L,ixC^L,idim,fRC)
       ! estimating bounds for the minimum and maximum signal velocities
       if(method_loc=='tvdlf'.or.method_loc=='tvdmu') then
         call phys_get_cbounds(wLC,wRC,wLp,wRp,x,ixI^L,ixC^L,idim,cmaxC)
       else
         call phys_get_cbounds(wLC,wRC,wLp,wRp,x,ixI^L,ixC^L,idim,cmaxC,&
                               cminC)
       end if

       ! use approximate Riemann solver to get flux at interfaces
       select case(method_loc)
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
    end do Loop_idims_flux ! Next idim
    block%iw0=0

    Loop_idims_wnew : do idim= idim^LIM
       hxO^L=ixO^L-kr(idim,^D);

       ! Multiply the fluxes by -dt/dx since Flux fixing expects this
       if (slab) then
          fC(ixI^S,1:nwflux,idim)=dxinv(idim)*fC(ixI^S,1:nwflux,idim)
          wnew(ixO^S,1:nwflux)=wnew(ixO^S,1:nwflux) &
               + (fC(ixO^S,1:nwflux,idim)-fC(hxO^S,1:nwflux,idim))
       else
          fC(ixI^S,1:nwflux,idim)=-qdt*fC(ixI^S,1:nwflux,idim)
          if (.not. angmomfix) then ! default case
            inv_volume = 1.0d0/block%dvolume(ixO^S)
            Loop_iw_new : do iw=1,nwflux
              wnew(ixO^S,iw)=wnew(ixO^S,iw) + (fC(ixO^S,iw,idim)-fC(hxO^S,iw,idim)) * &
                  inv_volume
            end do Loop_iw_new
          else
            ! If angular momentum conserving way to solve the equations,
            ! some fluxes additions need to be treated specifically
            call phys_angmomfix(fC,x,wnew,ixI^L,ixO^L,idim)
          endif
       end if

       ! For the MUSCL scheme apply the characteristic based limiter
       if (method_loc=='tvdmu') &
            call tvdlimit2(method,qdt,ixI^L,ixC^L,ixO^L,idim,wLC,wRC,wnew,x,fC,dx^D)

    end do Loop_idims_wnew ! Next idim

    if (.not.slab.and.idimmin==1) &
         call phys_add_source_geom(qdt,ixI^L,ixO^L,wCT,wnew,x)

    ! add source terms
    call addsource2(qdt*dble(idimmax-idimmin+1)/dble(ndim), &
         ixI^L,ixO^L,1,nw,qtC,wCT,qt,wnew,x,.false.)

    ! check and optionally correct unphysical values
    call phys_handle_small_values(.false.,wnew,x,ixI^L,ixO^L,&
                                  'finite_volume')

    !allow the user to change localy the numerical method
    !call usr_reset_solver(ixI^L,ixI^L,idim,qt,wnew,x,method,&
    !         type_limiter(node(plevel_,saveigrid)),method_loc,typelimiter)

  contains

    subroutine get_Riemann_flux_tvdmu()

      do iw=1,nwflux
         ! To save memory we use fLC to store (F_L+F_R)/2=half*(fLC+fRC)
         fLC(ixC^S, iw)=half*(fLC(ixC^S, iw)+fRC(ixC^S, iw))
         if (slab) then
           fC(ixC^S,iw,idim)=fLC(ixC^S, iw)
         else
           fC(ixC^S,iw,idim)=block%surfaceC(ixC^S,idim)*fLC(ixC^S, iw)
         end if
      end do
    end subroutine get_Riemann_flux_tvdmu

    subroutine get_Riemann_flux_tvdlf()
      real(kind=dp)    :: fac(ixC^S)

      fac = -0.5d0*tvdlfeps*cmaxC(ixC^S)

      ! Calculate fLC=f(uL_j+1/2) and fRC=f(uR_j+1/2) for each iw
      do iw=1,nwflux

         ! To save memory we use fLC to store (F_L+F_R)/2=half*(fLC+fRC)
         fLC(ixC^S, iw)=0.5_dp*(fLC(ixC^S, iw)+fRC(ixC^S, iw))

         ! Add TVDLF dissipation to the flux
         if (flux_type(idim, iw) /= flux_no_dissipation) then
            fLC(ixC^S, iw)=fLC(ixC^S, iw) + fac*(wRC(ixC^S,iw)-wLC(ixC^S,iw))
         end if

         if (slab) then
           fC(ixC^S,iw,idim)=fLC(ixC^S, iw)
         else
           fC(ixC^S,iw,idim)=block%surfaceC(ixC^S,idim)*fLC(ixC^S, iw)
         end if

      end do ! Next iw
    end subroutine get_Riemann_flux_tvdlf

    subroutine get_Riemann_flux_hll()

      real(kind=dp)    :: fac(ixC^S), div(ixC^S)

      where(cminC(ixC^S) >= zero)
        patchf(ixC^S) = -2
      elsewhere(cmaxC(ixC^S) <= zero)
        patchf(ixC^S) =  2
      elsewhere
        patchf(ixC^S) =  1
      endwhere

      fac = tvdlfeps*cminC(ixC^S)*cmaxC(ixC^S)
      div = 1/(cmaxC(ixC^S)-cminC(ixC^S))

      ! Calculate fLC=f(uL_j+1/2) and fRC=f(uR_j+1/2) for each iw
      do iw=1,nwflux
         if (flux_type(idim, iw) == flux_tvdlf) then
            fLC(ixC^S, iw) = half*(fLC(ixC^S, iw) + fRC(ixC^S, iw) &
                 -tvdlfeps*max(cmaxC(ixC^S), dabs(cminC(ixC^S))) * &
                 (wRC(ixC^S,iw)-wLC(ixC^S,iw)))
         else
            where(patchf(ixC^S)==1)
               ! Add hll dissipation to the flux
               fLC(ixC^S, iw) = (cmaxC(ixC^S)*fLC(ixC^S, iw)-cminC(ixC^S) * fRC(ixC^S, iw) &
                    +fac*(wRC(ixC^S,iw)-wLC(ixC^S,iw))) * div
            elsewhere(patchf(ixC^S)== 2)
               fLC(ixC^S, iw)=fRC(ixC^S, iw)
            elsewhere(patchf(ixC^S)==-2)
               fLC(ixC^S, iw)=fLC(ixC^S, iw)
            endwhere
         endif

         if (slab) then
           fC(ixC^S,iw,idim)=fLC(ixC^S, iw)
         else
           fC(ixC^S,iw,idim)=block%surfaceC(ixC^S,idim)*fLC(ixC^S, iw)
         end if

      end do ! Next iw
    end subroutine get_Riemann_flux_hll

    subroutine get_Riemann_flux_hllc()
      implicit none
      real(kind=dp)   , dimension(ixI^S,1:nwflux)     :: whll, Fhll, fCD
      real(kind=dp)   , dimension(ixI^S)              :: lambdaCD

      patchf(ixC^S) =  1
      where(cminC(ixC^S) >= zero)
         patchf(ixC^S) = -2
      elsewhere(cmaxC(ixC^S) <= zero)
         patchf(ixC^S) =  2
      endwhere
      ! Use more diffusive scheme, is actually TVDLF and selected by patchf=4
      if(method=='hllcd') &
           call phys_diffuse_hllcd(ixI^L,ixC^L,idim,wLC,wRC,fLC,fRC,patchf)


      ! now patchf may be -1 or 1 due to phys_get_lCD
      if(any(abs(patchf(ixC^S))== 1))then
         !---- calculate speed lambda at CD ----!
         call phys_get_lCD(wLC,wRC,fLC,fRC,cminC,cmaxC,idim,ixI^L,ixC^L, &
                           whll,Fhll,lambdaCD,patchf)
         !======== flux at intermediate state ========!
         call phys_get_wCD(x,wLC,wRC,whll,fRC,fLC,Fhll,patchf,lambdaCD,&
              cminC,cmaxC,ixI^L,ixC^L,idim,fCD)
      endif ! Calculate the CD flux

      do iw=1,nwflux
         if (flux_type(idim, iw) == flux_tvdlf) then
            fLC(ixC^S,iw) = 0.5_dp * (fLC(ixC^S,iw) + fRC(ixC^S,iw) &
                          - tvdlfeps *  max(cmaxC(ixC^S), abs(cminC(ixC^S))) &
                     * (wRC(ixC^S,iw) - wLC(ixC^S,iw)))
         else
            where(patchf(ixC^S)==-2)
               fLC(ixC^S,iw)=fLC(ixC^S,iw)
            elsewhere(abs(patchf(ixC^S))==1)
               fLC(ixC^S,iw)=fCD(ixC^S,iw)
            elsewhere(patchf(ixC^S)==2)
               fLC(ixC^S,iw)=fRC(ixC^S,iw)
            elsewhere(patchf(ixC^S)==3)
               ! fallback option, reducing to HLL flux
               fLC(ixC^S,iw)=Fhll(ixC^S,iw)
            elsewhere(patchf(ixC^S)==4)
               ! fallback option, reducing to TVDLF flux
               fLC(ixC^S,iw) = half*((fLC(ixC^S,iw)+fRC(ixC^S,iw)) &
                    -tvdlfeps * max(cmaxC(ixC^S), dabs(cminC(ixC^S))) * &
                    (wRC(ixC^S,iw)-wLC(ixC^S,iw)))
            endwhere
         end if

         if (slab) then
           fC(ixC^S,iw,idim)=fLC(ixC^S,iw)
         else
           fC(ixC^S,iw,idim)=block%surfaceC(ixC^S,idim)*fLC(ixC^S,iw)
         end if

      end do ! Next iw
    end subroutine get_Riemann_flux_hllc

    !> HLLD Riemann flux from Miyoshi 2005 JCP, 208, 315
    subroutine get_Riemann_flux_hlld()
      use mod_mhd_phys
      implicit none
      real(kind=dp)   , dimension(ixI^S,1:nwflux) :: w1R,w1L,f1R,f1L
      real(kind=dp)   , dimension(ixI^S,1:nwflux) :: w2R,w2L
      real(kind=dp)   , dimension(ixI^S) :: sm,s1R,s1L,suR,suL,Bx
      real(kind=dp)   , dimension(ixI^S) :: pts,ptR,ptL,signBx,r1L,r1R,tmp
      real(kind=dp)   , dimension(ixI^S,ndir) :: vRC, vLC
      integer :: ip1,ip2,ip3,idir

      f1R=0.d0
      f1L=0.d0
      ip1=idim
      ip3=3
      vRC(ixC^S,:)=wRp(ixC^S,mom(:))
      vLC(ixC^S,:)=wLp(ixC^S,mom(:))
      ! estimate normal magnetic field at cell interfaces
      Bx(ixC^S)=0.5d0*(wRC(ixC^S,mag(ip1))+wLC(ixC^S,mag(ip1)))
      suR(ixC^S)=(cmaxC(ixC^S)-vRC(ixC^S,ip1))*wRC(ixC^S,rho_)
      suL(ixC^S)=(cminC(ixC^S)-vLC(ixC^S,ip1))*wLC(ixC^S,rho_)
      ptR(ixC^S)=wRp(ixC^S,e_)+0.5d0*sum(wRC(ixC^S,mag(:))**2,dim=ndim+1)
      ptL(ixC^S)=wLp(ixC^S,e_)+0.5d0*sum(wLC(ixC^S,mag(:))**2,dim=ndim+1)
      ! equation (38)
      sm(ixC^S)=(suR(ixC^S)*vRC(ixC^S,ip1)-suL(ixC^S)*vLC(ixC^S,ip1)-&
                 ptR(ixC^S)+ptL(ixC^S))/(suR(ixC^S)-suL(ixC^S))
      ! equation (39)
      w1R(ixC^S,mom(ip1))=sm(ixC^S)
      w1L(ixC^S,mom(ip1))=sm(ixC^S)
      w2R(ixC^S,mom(ip1))=sm(ixC^S)
      w2L(ixC^S,mom(ip1))=sm(ixC^S)
      w1R(ixC^S,mag(ip1))=Bx(ixC^S)
      w1L(ixC^S,mag(ip1))=Bx(ixC^S)
      w2R(ixC^S,mag(ip1))=Bx(ixC^S)
      w2L(ixC^S,mag(ip1))=Bx(ixC^S)
      ! equation (41)
      pts(ixC^S)=(suR(ixC^S)*ptL(ixC^S)-suL(ixC^S)*ptR(ixC^S)+suR(ixC^S)*suL(ixC^S)*&
                 (vRC(ixC^S,ip1)-vLC(ixC^S,ip1)))/(suR(ixC^S)-suL(ixC^S))
      ! equation (43)
      w1R(ixC^S,rho_)=suR(ixC^S)/(cmaxC(ixC^S)-sm(ixC^S))
      w1L(ixC^S,rho_)=suL(ixC^S)/(cminC(ixC^S)-sm(ixC^S))
      ! equation (44) ~ (47)
      ip2=mod(ip1+1,ndir)
      if(ip2==0) ip2=ndir
      r1R(ixC^S)=suR(ixC^S)*(cmaxC(ixC^S)-sm(ixC^S))-Bx(ixC^S)**2
      where(r1R(ixC^S)/=0.d0)
        r1R(ixC^S)=1.d0/r1R(ixC^S)
      endwhere
      r1L(ixC^S)=suL(ixC^S)*(cminC(ixC^S)-sm(ixC^S))-Bx(ixC^S)**2
      where(r1L(ixC^S)/=0.d0)
        r1L(ixC^S)=1.d0/r1L(ixC^S)
      endwhere
      w1R(ixC^S,mom(ip2))=vRC(ixC^S,ip2)-Bx(ixC^S)*wRC(ixC^S,mag(ip2))*&
        (sm(ixC^S)-vRC(ixC^S,ip1))*r1R(ixC^S)
      w1L(ixC^S,mom(ip2))=vLC(ixC^S,ip2)-Bx(ixC^S)*wLC(ixC^S,mag(ip2))*&
        (sm(ixC^S)-vLC(ixC^S,ip1))*r1L(ixC^S)
      w1R(ixC^S,mag(ip2))=(suR(ixC^S)*(cmaxC(ixC^S)-vRC(ixC^S,ip1))-Bx(ixC^S)**2)*r1R(ixC^S)
      w1L(ixC^S,mag(ip2))=(suL(ixC^S)*(cminC(ixC^S)-vLC(ixC^S,ip1))-Bx(ixC^S)**2)*r1L(ixC^S)
      if(ndir==3) then
        ip3=mod(ip1+2,ndir)
        if(ip3==0) ip3=ndir
        w1R(ixC^S,mom(ip3))=vRC(ixC^S,ip3)-Bx(ixC^S)*wRC(ixC^S,mag(ip3))*&
          (sm(ixC^S)-vRC(ixC^S,ip1))*r1R(ixC^S)
        w1L(ixC^S,mom(ip3))=vLC(ixC^S,ip3)-Bx(ixC^S)*wLC(ixC^S,mag(ip3))*&
          (sm(ixC^S)-vLC(ixC^S,ip1))*r1L(ixC^S)
        w1R(ixC^S,mag(ip3))=wRC(ixC^S,mag(ip3))*w1R(ixC^S,mag(ip2))
        w1L(ixC^S,mag(ip3))=wLC(ixC^S,mag(ip3))*w1L(ixC^S,mag(ip2))
      end if
      w1R(ixC^S,mag(ip2))=wRC(ixC^S,mag(ip2))*w1R(ixC^S,mag(ip2))
      w1L(ixC^S,mag(ip2))=wLC(ixC^S,mag(ip2))*w1L(ixC^S,mag(ip2))
      ! equation (48)
      if(mhd_energy) then
        w1R(ixC^S,e_)=((cmaxC(ixC^S)-vRC(ixC^S,ip1))*wRC(ixC^S,e_)-ptR(ixC^S)*vRC(ixC^S,ip1)+&
          pts(ixC^S)*sm(ixC^S)+Bx(ixC^S)*(sum(vRC(ixC^S,:)*wRC(ixC^S,mag(:)),dim=ndim+1)-&
          sum(w1R(ixC^S,mom(:))*w1R(ixC^S,mag(:)),dim=ndim+1)))/(cmaxC(ixC^S)-sm(ixC^S))
        w1L(ixC^S,e_)=((cminC(ixC^S)-vLC(ixC^S,ip1))*wLC(ixC^S,e_)-ptL(ixC^S)*vLC(ixC^S,ip1)+&
          pts(ixC^S)*sm(ixC^S)+Bx(ixC^S)*(sum(vLC(ixC^S,:)*wLC(ixC^S,mag(:)),dim=ndim+1)-&
          sum(w1L(ixC^S,mom(:))*w1L(ixC^S,mag(:)),dim=ndim+1)))/(cminC(ixC^S)-sm(ixC^S))
      end if
      ! equation (49)
      w2R(ixC^S,rho_)=w1R(ixC^S,rho_)
      w2L(ixC^S,rho_)=w1L(ixC^S,rho_)
      r1R(ixC^S)=sqrt(w1R(ixC^S,rho_))
      r1L(ixC^S)=sqrt(w1L(ixC^S,rho_))
      tmp(ixC^S)=1.d0/(r1R(ixC^S)+r1L(ixC^S))
      signBx(ixC^S)=sign(1.d0,Bx(ixC^S))
      ! equation (51)
      s1R(ixC^S)=sm(ixC^S)+abs(Bx(ixC^S))/r1R(ixC^S)
      s1L(ixC^S)=sm(ixC^S)-abs(Bx(ixC^S))/r1L(ixC^S)
      ! equation (59)
      w2R(ixC^S,mom(ip2))=(r1L(ixC^S)*w1L(ixC^S,mom(ip2))+r1R(ixC^S)*w1R(ixC^S,mom(ip2))+&
          (w1R(ixC^S,mag(ip2))-w1L(ixC^S,mag(ip2)))*signBx(ixC^S))*tmp(ixC^S)
      w2L(ixC^S,mom(ip2))=w2R(ixC^S,mom(ip2))
      ! equation (61)
      w2R(ixC^S,mag(ip2))=(r1L(ixC^S)*w1R(ixC^S,mag(ip2))+r1R(ixC^S)*w1L(ixC^S,mag(ip2))+&
          r1L(ixC^S)*r1R(ixC^S)*(w1R(ixC^S,mom(ip2))-w1L(ixC^S,mom(ip2)))*signBx(ixC^S))*tmp(ixC^S)
      w2L(ixC^S,mag(ip2))=w2R(ixC^S,mag(ip2))
      if(ndir==3) then
        ! equation (60)
        w2R(ixC^S,mom(ip3))=(r1L(ixC^S)*w1L(ixC^S,mom(ip3))+r1R(ixC^S)*w1R(ixC^S,mom(ip3))+&
            (w1R(ixC^S,mag(ip3))-w1L(ixC^S,mag(ip3)))*signBx(ixC^S))*tmp(ixC^S)
        w2L(ixC^S,mom(ip3))=w2R(ixC^S,mom(ip3))
        ! equation (62)
        w2R(ixC^S,mag(ip3))=(r1L(ixC^S)*w1R(ixC^S,mag(ip3))+r1R(ixC^S)*w1L(ixC^S,mag(ip3))+&
            r1L(ixC^S)*r1R(ixC^S)*(w1R(ixC^S,mom(ip3))-w1L(ixC^S,mom(ip3)))*signBx(ixC^S))*tmp(ixC^S)
        w2L(ixC^S,mag(ip3))=w2R(ixC^S,mag(ip3))
      end if
      ! equation (63)
      if(mhd_energy) then
        w2R(ixC^S,e_)=w1R(ixC^S,e_)+r1R(ixC^S)*(sum(w1R(ixC^S,mom(:))*w1R(ixC^S,mag(:)),dim=ndim+1)-&
          sum(w2R(ixC^S,mom(:))*w2R(ixC^S,mag(:)),dim=ndim+1))*signBx(ixC^S)
        w2L(ixC^S,e_)=w1L(ixC^S,e_)-r1L(ixC^S)*(sum(w1L(ixC^S,mom(:))*w1L(ixC^S,mag(:)),dim=ndim+1)-&
          sum(w2L(ixC^S,mom(:))*w2L(ixC^S,mag(:)),dim=ndim+1))*signBx(ixC^S)
      end if
      ! convert velocity to momentum
      do idir=1,ndir
        w1R(ixC^S,mom(idir))=w1R(ixC^S,mom(idir))*w1R(ixC^S,rho_)
        w1L(ixC^S,mom(idir))=w1L(ixC^S,mom(idir))*w1L(ixC^S,rho_)
        w2R(ixC^S,mom(idir))=w2R(ixC^S,mom(idir))*w2R(ixC^S,rho_)
        w2L(ixC^S,mom(idir))=w2L(ixC^S,mom(idir))*w2L(ixC^S,rho_)
      end do
      ! get fluxes of intermedate states
      do iw=1,nwflux
        if (flux_type(idim, iw) == flux_tvdlf) then
          fC(ixC^S,iw,ip1)=0.5d0*(fLC(ixC^S,iw) + fRC(ixC^S,iw))
          !fC(ixC^S,iw,ip1) = 0.5d0 * (fLC(ixC^S,iw) + fRC(ixC^S,iw) - tvdlfeps * &
          !     max(cmaxC(ixC^S), abs(cminC(ixC^S))) * &
          !     (wRC(ixC^S,iw)-wLC(ixC^S,iw)))
          cycle
        end if
        f1L(ixC^S,iw)=fLC(ixC^S,iw)+cminC(ixC^S)*(w1L(ixC^S,iw)-wLC(ixC^S,iw))
        f1R(ixC^S,iw)=fRC(ixC^S,iw)+cmaxC(ixC^S)*(w1R(ixC^S,iw)-wRC(ixC^S,iw))
        where(cminC(ixC^S)>0.d0)
          fC(ixC^S,iw,ip1)=fLC(ixC^S,iw)
        else where(s1L(ixC^S)>=0.d0)
          fC(ixC^S,iw,ip1)=f1L(ixC^S,iw)
        else where(sm(ixC^S)>=0.d0)
          fC(ixC^S,iw,ip1)=f1L(ixC^S,iw)+s1L(ixC^S)*(w2L(ixC^S,iw)-w1L(ixC^S,iw))
        else where(s1R(ixC^S)>=0.d0)
          fC(ixC^S,iw,ip1)=f1R(ixC^S,iw)+s1R(ixC^S)*(w2R(ixC^S,iw)-w1R(ixC^S,iw))
        else where(cmaxC(ixC^S)>=0.d0)
          fC(ixC^S,iw,ip1)=f1R(ixC^S,iw)
        else where(cmaxC(ixC^S)<0.d0)
          fC(ixC^S,iw,ip1)=fRC(ixC^S,iw)
        end where
        if(.not.slab) then
          fC(ixC^S,iw,ip1)=block%surfaceC(ixC^S,ip1)*fC(ixC^S,iw,ip1)
        end if
      end do

    end subroutine get_Riemann_flux_hlld

  end subroutine finite_volume

  !> Determine the upwinded wLC(ixL) and wRC(ixR) from w.
  !> the wCT is only used when PPM is exploited.
  subroutine reconstruct_LR(ixI^L,ixL^L,ixR^L,idim,w,wCT,wLC,&
                            wRC,wLp,wRp,x,needprim)
    use mod_physics
    use mod_global_parameters
    use mod_limiter

    integer, intent(in) :: ixI^L, ixL^L, ixR^L, idim
    logical, intent(in) :: needprim
    real(kind=dp)   , dimension(ixI^S,1:nw) :: w, wCT
    ! left and right constructed status in conservative form
    real(kind=dp)   , dimension(ixI^S,1:nw) :: wLC, wRC
    ! left and right constructed status in primitive form
    real(kind=dp)   , dimension(ixI^S,1:nw) :: wLp, wRp
    real(kind=dp)   , dimension(ixI^S,1:ndim) :: x

    integer            :: jxR^L, ixC^L, jxC^L, iw
    real(kind=dp)      :: ldw(ixI^S), rdw(ixI^S), dwC(ixI^S)
    logical            :: limiter_need_log
    ! integer            :: flagL(ixI^S), flagR(ixI^S)
     limiter_need_log=.false.
    ! Transform w,wL,wR to primitive variables
    if (needprim) then
       call phys_to_primitive(ixI^L,ixI^L,w,x)
    end if

    if (typelimiter == limiter_mp5) then
       call MP5limiter(ixI^L,ixL^L,idim,w,wLp,wRp)
    else if (typelimiter == limiter_ppm) then
       call PPMlimiter(ixI^L,ixM^LL,idim,w,wCT,wLp,wRp)
    else
       jxR^L=ixR^L+kr(idim,^D);
       ixCmax^D=jxRmax^D; ixCmin^D=ixLmin^D-kr(idim,^D);
       jxC^L=ixC^L+kr(idim,^D);

       do iw=1,nwflux

          if (loglimit(iw)) then
             limiter_need_log=maxval(dabs(wRp(ixR^S,iw)-wLp(ixL^S,iw))&
                                    /min(wRp(ixR^S,iw),wLp(ixL^S,iw)))>limiter_log_threshold
             if(limiter_need_log)then
              w(ixCmin^D:jxCmax^D,iw)=dlog10(w(ixCmin^D:jxCmax^D,iw))
              wLp(ixL^S,iw)=dlog10(wLp(ixL^S,iw))
              wRp(ixR^S,iw)=dlog10(wRp(ixR^S,iw))
             end if
          end if

          dwC(ixC^S)=w(jxC^S,iw)-w(ixC^S,iw)

          ! limit flux from left and/or right
          call dwlimiter2(dwC,ixI^L,ixC^L,idim,typelimiter,ldw,rdw)
          wLp(ixL^S,iw)=wLp(ixL^S,iw)+half*ldw(ixL^S)
          wRp(ixR^S,iw)=wRp(ixR^S,iw)-half*rdw(jxR^S)

          if (loglimit(iw).and.limiter_need_log) then
             w(ixCmin^D:jxCmax^D,iw)=10.0d0**w(ixCmin^D:jxCmax^D,iw)
             wLp(ixL^S,iw)=10.0d0**wLp(ixL^S,iw)
             wRp(ixR^S,iw)=10.0d0**wRp(ixR^S,iw)
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
       call phys_to_conserved(ixI^L,ixI^L,w,x)
    endif
    wLC(ixL^S,1:nw)=wLp(ixL^S,1:nw)
    wRC(ixR^S,1:nw)=wRp(ixR^S,1:nw)
    call phys_to_conserved(ixI^L,ixL^L,wLC,x)
    call phys_to_conserved(ixI^L,ixR^L,wRC,x)

  end subroutine reconstruct_LR

end module mod_finite_volume
