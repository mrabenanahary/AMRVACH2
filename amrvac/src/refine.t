!=============================================================================
subroutine refine_grids(child_igrid,child_ipe,igrid,ipe,active)

use mod_global_parameters

integer, dimension(2^D&), intent(in) :: child_igrid, child_ipe
integer, intent(in) :: igrid, ipe
logical, intent(in) :: active

integer :: ic^D
!-----------------------------------------------------------------------------

! allocate solution space for new children
{do ic^DB=1,2\}
   call alloc_node(child_igrid(ic^D))
{end do\}

if ((time_advance .and. active).or.convert.or.firstprocess) then
   ! prolong igrid to new children
   call prolong_grid(child_igrid,child_ipe,igrid,ipe)
else
   ! Fill new created children with initial condition
   {do ic^DB=1,2\}
      call initial_condition(child_igrid(ic^D))
   {end do\}
end if

! remove solution space of igrid
!call dealloc_node(igrid)

end subroutine refine_grids
!=============================================================================
subroutine prolong_grid(child_igrid,child_ipe,igrid,ipe)
  use mod_physics, only: phys_to_primitive, phys_to_conserved
  use mod_global_parameters

  integer, dimension(2^D&), intent(in) :: child_igrid, child_ipe
  integer, intent(in) :: igrid, ipe

  integer :: ix^L, ichild, ixCo^L, ic^D
  double precision :: dxCo^D, xComin^D, dxFi^D, xFimin^D

  if (prolongation_method=="linear") then
     ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
     {#IFDEF EVOLVINGBOUNDARY
     if (phyboundblock(igrid)) then
        ix^L=ixG^LL;
     else
        ix^L=ixM^LL^LADD1;
     end if
     }{#IFNDEF EVOLVINGBOUNDARY
     ix^L=ixM^LL^LADD1;
     }

     if(prolongprimitive) call phys_to_primitive(ixG^LL,ix^L,pw(igrid)%w,pw(igrid)%x)

     xComin^D=rnode(rpxmin^D_,igrid)\
     dxCo^D=rnode(rpdx^D_,igrid)\
  end if

  {do ic^DB=1,2\}
  ichild=child_igrid(ic^D)

  ixComin^D=ixMlo^D+(ic^D-1)*block_nx^D/2\
  ixComax^D=ixMhi^D+(ic^D-2)*block_nx^D/2\

  if (prolongation_method=="linear") then
     xFimin^D=rnode(rpxmin^D_,ichild)\
     dxFi^D=rnode(rpdx^D_,ichild)\
     call prolong_2nd(pw(igrid)%w,pw(igrid)%x,ixCo^L,pw(ichild)%w,pw(ichild)%x, &
          dxCo^D,xComin^D,dxFi^D,xFimin^D,igrid,ichild)
  else
     call prolong_1st(pw(igrid)%w,ixCo^L,pw(ichild)%w,pw(ichild)%x)
  end if
  {end do\}

  if (prolongation_method=="linear" .and. prolongprimitive) then
     call phys_to_conserved(ixG^LL,ix^L,pw(igrid)%w,pw(igrid)%x)
  end if

end subroutine prolong_grid
!=============================================================================
subroutine prolong_2nd(wCo,xCo,ixCo^L,wFi,xFi,dxCo^D,xComin^D,dxFi^D,xFimin^D,igridCo,igridFi)

use mod_physics, only: phys_to_conserved
use mod_global_parameters

integer, intent(in) :: ixCo^L, igridFi, igridCo
double precision, intent(in) :: dxCo^D, xComin^D, dxFi^D, xFimin^D
double precision, intent(in) :: wCo(ixG^T,nw), xCo(ixG^T,1:ndim), xFi(ixG^T,1:ndim)
double precision, intent(inout) :: wFi(ixG^T,nw)

integer :: ixCo^D, jxCo^D, hxCo^D, ixFi^D, ix^D, idim, iw, ixCg^L, el
double precision :: slopeL, slopeR, slopeC, signC, signR
double precision :: slope(nw,ndim)
double precision :: eta^D

!-----------------------------------------------------------------------------
ixCg^L=ixCo^L;
{#IFDEF EVOLVINGBOUNDARY
if (phyboundblock(ichild)) then
  el=ceiling(real(nghostcells)/2.)
  ixCgmin^D=ixComin^D-el\
  ixCgmax^D=ixComax^D+el\
end if
}
{do ixCo^DB = ixCg^LIM^DB
   ! lower left grid index in finer child block
   ixFi^DB=2*(ixCo^DB-ixComin^DB)+ixMlo^DB\}

   do idim=1,ndim
      hxCo^D=ixCo^D-kr(^D,idim)\
      jxCo^D=ixCo^D+kr(^D,idim)\

      do iw=1,nw
         slopeL=wCo(ixCo^D,iw)-wCo(hxCo^D,iw)
         slopeR=wCo(jxCo^D,iw)-wCo(ixCo^D,iw)
         slopeC=half*(slopeR+slopeL)

         ! get limited slope
         signR=sign(one,slopeR)
         signC=sign(one,slopeC)
         select case(typeprolonglimit)
         case('unlimit')
           slope(iw,idim)=slopeC
         case('minmod')
           slope(iw,idim)=signR*max(zero,min(dabs(slopeR), &
                                             signR*slopeL))
         case('woodward')
           slope(iw,idim)=two*signR*max(zero,min(dabs(slopeR), &
                              signR*slopeL,signR*half*slopeC))
         case('koren')
           slope(iw,idim)=signR*max(zero,min(two*signR*slopeL, &
            (dabs(slopeR)+two*slopeL*signR)*third,two*dabs(slopeR)))
         case default
           slope(iw,idim)=signC*max(zero,min(dabs(slopeC), &
                             signC*slopeL,signC*slopeR))
         end select
      end do
   end do
   ! cell-centered coordinates of coarse grid point
   !^D&xCo^D=xCo({ixCo^DD},^D)\
   !if(.not.slab) then
   !  ^D&invdxCo^D=1.d0/pw(igridCo)%dx(ixCo^DD,^D)\
   !endif
   {do ix^DB=ixFi^DB,ixFi^DB+1 \}
      ! cell-centered coordinates of fine grid point
      !^D&xFi^D=xFi({ix^DD},^D)\
      if(slab) then
        ! normalized distance between fine/coarse cell center
        ! in coarse cell: ranges from -0.5 to 0.5 in each direction
        ! (origin is coarse cell center)
        ! hence this is +1/4 or -1/4 on cartesian mesh
        !eta^D=(xFi^D-xCo^D)*invdxCo^D;
        eta^D=0.5d0*(dble(ix^D-ixFi^D)-0.5d0);
      else
        {! forefactor is -0.5d0 when ix=ixFi and +0.5d0 for ixFi+1
        eta^D=(dble(ix^D-ixFi^D)-0.5d0)*(one-pw(igridFi)%dvolume(ix^DD) &
              /sum(pw(igridFi)%dvolume(ixFi^D:ixFi^D+1^D%ix^DD)))  \}
      end if
      wFi(ix^D,1:nw) = wCo(ixCo^D,1:nw) &
                            + {(slope(1:nw,^D)*eta^D)+}
   {end do\}
{end do\}


if(prolongprimitive) call phys_to_conserved(ixG^LL,ixM^LL,wFi,xFi)

end subroutine prolong_2nd
!=============================================================================
subroutine prolong_1st(wCo,ixCo^L,wFi,xFi)

use mod_global_parameters

integer, intent(in) :: ixCo^L
double precision, intent(in) :: wCo(ixG^T,nw), xFi(ixG^T,1:ndim)
double precision, intent(out) :: wFi(ixG^T,nw)

integer :: ixCo^D, ixFi^D, iw
integer :: ixFi^L
!-----------------------------------------------------------------------------
{do ixCo^DB = ixCo^LIM^DB
   ixFi^DB=2*(ixCo^DB-ixComin^DB)+ixMlo^DB\}
   forall(iw=1:nw) wFi(ixFi^D:ixFi^D+1,iw)=wCo(ixCo^D,iw)
{end do\}

end subroutine prolong_1st
!=============================================================================
