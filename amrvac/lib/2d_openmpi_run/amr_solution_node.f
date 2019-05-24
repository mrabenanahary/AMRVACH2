!> Get first available igrid on processor ipe
integer function getnode(ipe)
use mod_forest, only: igrid_inuse
use mod_global_parameters

integer, intent(in) :: ipe
integer :: igrid, igrid_available

igrid_available=0

do igrid=1,max_blocks
   if (igrid_inuse(igrid,ipe)) cycle

   igrid_available=igrid
   exit
end do

if (igrid_available == 0) then
   getnode = -1
   print *, "Current maximum number of grid blocks:", max_blocks
   call mpistop("Insufficient grid blocks; increase max_blocks in meshlist")
else
   getnode=igrid_available
   igrid_inuse(igrid,ipe)=.true.
end if

if (ipe==mype) then
   ! initialize nodal block
   node(1:nodehi,getnode) = 0
   rnode(1:rnodehi,getnode) = zero
end if

end function getnode
!=============================================================================
subroutine putnode(igrid,ipe)
use mod_forest
implicit none

! putnode = return igrid node on processor ipe
 
integer, intent(in) :: igrid, ipe
!----------------------------------------------------------------------------
igrid_inuse(igrid,ipe)=.false.

end subroutine putnode
!=============================================================================
subroutine alloc_node(igrid)
use mod_forest
use mod_global_parameters
use mod_geometry

integer, intent(in) :: igrid

integer :: level, ig1,ig2, ign1,ign2, ixCoGmin1,ixCoGmin2,ixCoGmax1,ixCoGmax2,&
    ix, i1,i2
integer :: imin, imax, index, igCo1,igCo2, ixshift, offset, ifirst
integer:: icase, ixGextmin1,ixGextmin2,ixGextmax1,ixGextmax2
double precision :: rXmin1,rXmin2, dx1,dx2, summeddx, sizeuniformpart1,&
   sizeuniformpart2
double precision :: xext(ixGlo1-1:ixGhi1+1,ixGlo2-1:ixGhi2+1,1:ndim)
!-----------------------------------------------------------------------------
ixCoGmin1=1;ixCoGmin2=1;
!ixCoGmax^D=ixGhi^D/2+nghostcells;
ixCoGmax1=(ixGhi1-2*nghostcells)/2+2*nghostcells
ixCoGmax2=(ixGhi2-2*nghostcells)/2+2*nghostcells;

icase=mod(nghostcells,2)
select case(icase)
   case(0)
    ixGextmin1=ixGlo1;ixGextmin2=ixGlo2;ixGextmax1=ixGhi1;ixGextmax2=ixGhi2;
   case(1)
    ! for ghost cell related prolongations, we need
    ! an extra layer with known volumes and dx-intervals
    ! in case the number of ghost cells is odd
    ixGextmin1=ixGlo1-1;ixGextmin2=ixGlo2-1;ixGextmax1=ixGhi1+1
    ixGextmax2=ixGhi2+1;
   case default
     call mpistop("no such case")
end select

if(.not. allocated(pw(igrid)%w)) then
  
  ! initialize solution space
  allocate(pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nw),&
      pw(igrid)%wold(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nw),&
      pw(igrid)%w1(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nw),&
      pw(igrid)%wcoarse(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,1:nw))
  
  ! wio for visualization data
  allocate(pw(igrid)%wio(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nw+nwauxio))
  
  ! allocate temperary solution space
  select case (time_integrator)
  case("threestep","fourstep","jameson","twostep_trapezoidal")
    allocate(pw(igrid)%w2(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nw))
  case("rk4","ssprk43")
    allocate(pw(igrid)%w2(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nw))
    allocate(pw(igrid)%w3(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nw))
  case("ssprk54")
    allocate(pw(igrid)%w2(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nw))
    allocate(pw(igrid)%w3(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nw))
    allocate(pw(igrid)%w4(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nw))
  end select
  
  ! allocate coordinates
  allocate(pw(igrid)%x(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:ndim),&
      pw(igrid)%xcoarse(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,1:ndim))
  if(.not.slab)then
      allocate(pw(igrid)%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
         1:ndim), pw(igrid)%dxcoarse(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,&
         1:ndim),pw(igrid)%ds(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
         1:ndim))
      allocate(pw(igrid)%dvolume(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2),&
          pw(igrid)%dvolumecoarse(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2))
      allocate(pw(igrid)%surfaceC(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:ndim),&
          pw(igrid)%surface(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:ndim))
  endif
end if

pw(igrid)%w=0.d0
pw(igrid)%wold=0.d0
pw(igrid)%w1=0.d0
pw(igrid)%wcoarse=0.d0
pw(igrid)%igrid=igrid
pw(igrid)%wio=0.d0

! wb is w by default
pw(igrid)%wb=>pw(igrid)%w


! block pointer to current block
block=>pw(igrid)

if(phys_energy .and. solve_internal_e) then
  block%e_is_internal=.true.
endif

! set level information
level=igrid_to_node(igrid,mype)%node%level
ig1=igrid_to_node(igrid,mype)%node%ig1;ig2=igrid_to_node(igrid,mype)%node%ig2;

node(plevel_,igrid)=level
node(pig1_,igrid)=ig1
node(pig2_,igrid)=ig2

! set dx information
rnode(rpdx1_,igrid)=dx(1,level)
rnode(rpdx2_,igrid)=dx(2,level)
dxlevel(:)=dx(:,level)

! uniform cartesian case as well as all unstretched coordinates
! determine the minimal and maximal corners
rnode(rpxmin1_,igrid)=xprobmin1+dble(ig1-1)*dg1(level)
rnode(rpxmin2_,igrid)=xprobmin2+dble(ig2-1)*dg2(level)
rnode(rpxmax1_,igrid)=xprobmax1-dble(ng1(level)-ig1)*dg1(level)
rnode(rpxmax2_,igrid)=xprobmax2-dble(ng2(level)-ig2)*dg2(level)

dx1=rnode(rpdx1_,igrid)
dx2=rnode(rpdx2_,igrid)
rXmin1=rnode(rpxmin1_,igrid)-nghostcells*dx1
rXmin2=rnode(rpxmin2_,igrid)-nghostcells*dx2
do ix=ixGlo1,ixGhi1
   pw(igrid)%x(ix,ixGlo2:ixGhi2,1)=rXmin1+(dble(ix)-half)*dx1
end do
do ix=ixGlo2,ixGhi2
   pw(igrid)%x(ixGlo1:ixGhi1,ix,2)=rXmin2+(dble(ix)-half)*dx2
end do
dx1=2.0d0*rnode(rpdx1_,igrid)
dx2=2.0d0*rnode(rpdx2_,igrid)
rXmin1=rnode(rpxmin1_,igrid)-nghostcells*dx1
rXmin2=rnode(rpxmin2_,igrid)-nghostcells*dx2
do ix=ixCoGmin1,ixCoGmax1
   pw(igrid)%xcoarse(ix,ixCoGmin2:ixCoGmax2,1)=rXmin1+(dble(ix)-half)*dx1
end do
do ix=ixCoGmin2,ixCoGmax2
   pw(igrid)%xcoarse(ixCoGmin1:ixCoGmax1,ix,2)=rXmin2+(dble(ix)-half)*dx2
end do

if (.not.slab) then
  pw(igrid)%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,1)=rnode(rpdx1_,&
     igrid)
  pw(igrid)%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,2)=rnode(rpdx2_,&
     igrid);
  pw(igrid)%dxcoarse(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,&
     1)=2.0d0*rnode(rpdx1_,igrid)
  pw(igrid)%dxcoarse(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,&
     2)=2.0d0*rnode(rpdx2_,igrid);
  dx1=rnode(rpdx1_,igrid)
  dx2=rnode(rpdx2_,igrid)
  rXmin1=rnode(rpxmin1_,igrid)-nghostcells*dx1
  rXmin2=rnode(rpxmin2_,igrid)-nghostcells*dx2
  do ix=ixGextmin1,ixGextmax1
      xext(ix,ixGextmin2:ixGextmax2,1)=rXmin1+(dble(ix)-half)*dx1
  end do
  do ix=ixGextmin2,ixGextmax2
      xext(ixGextmin1:ixGextmax1,ix,2)=rXmin2+(dble(ix)-half)*dx2
  end do
endif

if(stretched_grid) then
  if(stretched_dim(1))then
    imin=(ig1-1)*block_nx1
    imax=ig1*block_nx1
    rnode(rpxmin1_,igrid)=xprobmin1+dxfirst_1mq(level,&
       1) *(1.0d0-qstretch(level,1)**imin)
    rnode(rpxmax1_,igrid)=xprobmin1+dxfirst_1mq(level,&
       1) *(1.0d0-qstretch(level,1)**imax)
    ! fix possible out of bound due to precision
    if(rnode(rpxmax1_,igrid)>xprobmax1) rnode(rpxmax1_,igrid)=xprobmax1
    ixshift=(ig1-1)*block_nx1-nghostcells
    do ix=ixGextmin1,ixGextmax1
      index=ixshift+ix
      pw(igrid)%dx(ix,ixGlo2:ixGhi2,1)=dxfirst(level,1)*qstretch(level,&
         1)**(index-1)
    enddo
    igCo1=(ig1-1)/2
    ixshift=igCo1*block_nx1+(1-modulo(ig1,2))*block_nx1/2-nghostcells
    do ix=ixCoGmin1,ixCoGmax1
      index=ixshift+ix
      pw(igrid)%dxcoarse(ix,ixCoGmin2:ixCoGmax2,1)=dxfirst(level-1,&
         1)*qstretch(level-1,1)**(index-1)
      pw(igrid)%xcoarse(ix,ixCoGmin2:ixCoGmax2,&
         1)=xprobmin1+dxfirst_1mq(level-1,1)*(1.0d0-qstretch(level-1,&
         1)**(index-1)) + 0.5d0*dxfirst(level-1,1)*qstretch(level-1,&
         1)**(index-1)
    end do
    ! now that dx and grid boundaries are known: fill cell centers
    ifirst=nghostcells+1
    ! first fill the mesh
    summeddx=0.0d0
    do ix=ixMlo1,ixMhi1
       pw(igrid)%x(ix,ixGlo2:ixGhi2,1)=rnode(rpxmin1_,&
          igrid)+summeddx+0.5d0*pw(igrid)%dx(ix,ixGlo2:ixGhi2,1)
       summeddx=summeddx+pw(igrid)%dx(ix,ifirst,1)
    enddo
    ! then ghost cells to left
    summeddx=0.0d0
    do ix=nghostcells,1,-1
       pw(igrid)%x(ix,ixGlo2:ixGhi2,1)=rnode(rpxmin1_,&
          igrid)-summeddx-0.5d0*pw(igrid)%dx(ix,ixGlo2:ixGhi2,1)
       summeddx=summeddx+pw(igrid)%dx(ix,ifirst,1)
    enddo
    ! then ghost cells to right
    summeddx=0.0d0
    do ix=ixGhi1-nghostcells+1,ixGhi1
       pw(igrid)%x(ix,ixGlo2:ixGhi2,1)=rnode(rpxmax1_,&
          igrid)+summeddx+0.5d0*pw(igrid)%dx(ix,ixGlo2:ixGhi2,1)
       summeddx=summeddx+pw(igrid)%dx(ix,ifirst,1)
    enddo
    select case(icase)
      case(0)
        ! if even number of ghost cells: xext is just copy of local x
        xext(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
           1)=pw(igrid)%x(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,1)
      case(1)
        ! if uneven number of ghost cells: extra layer left/right
        summeddx=0.0d0
        do ix=ixMlo1,ixMhi1
          xext(ix,ixGextmin2:ixGextmax2,1)=rnode(rpxmin1_,&
             igrid)+summeddx+0.5d0*pw(igrid)%dx(ix,ixGextmin2:ixGextmax2,1)
         summeddx=summeddx+pw(igrid)%dx(ix,ifirst,1)
        enddo
        ! then ghost cells to left
        summeddx=0.0d0
        do ix=nghostcells,ixGextmin1,-1
          xext(ix,ixGextmin2:ixGextmax2,1)=rnode(rpxmin1_,&
             igrid)-summeddx-0.5d0*pw(igrid)%dx(ix,ixGextmin2:ixGextmax2,1)
          summeddx=summeddx+pw(igrid)%dx(ix,ifirst,1)
        enddo
       ! then ghost cells to right
       summeddx=0.0d0
       do ix=ixGhi1-nghostcells+1,ixGextmax1
          xext(ix,ixGextmin2:ixGextmax2,1)=rnode(rpxmax1_,&
             igrid)+summeddx+0.5d0*pw(igrid)%dx(ix,ixGextmin2:ixGextmax2,1)
          summeddx=summeddx+pw(igrid)%dx(ix,ifirst,1)
       enddo
      case default
        call mpistop("no such case")
    end select
   endif
  if(stretched_dim(2))then
    imin=(ig2-1)*block_nx2
    imax=ig2*block_nx2
    rnode(rpxmin2_,igrid)=xprobmin2+dxfirst_1mq(level,&
       2) *(1.0d0-qstretch(level,2)**imin)
    rnode(rpxmax2_,igrid)=xprobmin2+dxfirst_1mq(level,&
       2) *(1.0d0-qstretch(level,2)**imax)
    ! fix possible out of bound due to precision
    if(rnode(rpxmax2_,igrid)>xprobmax2) rnode(rpxmax2_,igrid)=xprobmax2
    ixshift=(ig2-1)*block_nx2-nghostcells
    do ix=ixGextmin2,ixGextmax2
      index=ixshift+ix
      pw(igrid)%dx(ixGlo1:ixGhi1,ix,2)=dxfirst(level,2)*qstretch(level,&
         2)**(index-1)
    enddo
    igCo2=(ig2-1)/2
    ixshift=igCo2*block_nx2+(1-modulo(ig2,2))*block_nx2/2-nghostcells
    do ix=ixCoGmin2,ixCoGmax2
      index=ixshift+ix
      pw(igrid)%dxcoarse(ixCoGmin1:ixCoGmax1,ix,2)=dxfirst(level-1,&
         2)*qstretch(level-1,2)**(index-1)
      pw(igrid)%xcoarse(ixCoGmin1:ixCoGmax1,ix,&
         2)=xprobmin2+dxfirst_1mq(level-1,2)*(1.0d0-qstretch(level-1,&
         2)**(index-1)) + 0.5d0*dxfirst(level-1,2)*qstretch(level-1,&
         2)**(index-1)
    end do
    ! now that dx and grid boundaries are known: fill cell centers
    ifirst=nghostcells+1
    ! first fill the mesh
    summeddx=0.0d0
    do ix=ixMlo2,ixMhi2
       pw(igrid)%x(ixGlo1:ixGhi1,ix,2)=rnode(rpxmin2_,&
          igrid)+summeddx+0.5d0*pw(igrid)%dx(ixGlo1:ixGhi1,ix,2)
       summeddx=summeddx+pw(igrid)%dx(ifirst,ix,2)
    enddo
    ! then ghost cells to left
    summeddx=0.0d0
    do ix=nghostcells,1,-1
       pw(igrid)%x(ixGlo1:ixGhi1,ix,2)=rnode(rpxmin2_,&
          igrid)-summeddx-0.5d0*pw(igrid)%dx(ixGlo1:ixGhi1,ix,2)
       summeddx=summeddx+pw(igrid)%dx(ifirst,ix,2)
    enddo
    ! then ghost cells to right
    summeddx=0.0d0
    do ix=ixGhi2-nghostcells+1,ixGhi2
       pw(igrid)%x(ixGlo1:ixGhi1,ix,2)=rnode(rpxmax2_,&
          igrid)+summeddx+0.5d0*pw(igrid)%dx(ixGlo1:ixGhi1,ix,2)
       summeddx=summeddx+pw(igrid)%dx(ifirst,ix,2)
    enddo
    select case(icase)
      case(0)
        ! if even number of ghost cells: xext is just copy of local x
        xext(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
           2)=pw(igrid)%x(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,2)
      case(1)
        ! if uneven number of ghost cells: extra layer left/right
        summeddx=0.0d0
        do ix=ixMlo2,ixMhi2
          xext(ixGextmin1:ixGextmax1,ix,2)=rnode(rpxmin2_,&
             igrid)+summeddx+0.5d0*pw(igrid)%dx(ixGextmin1:ixGextmax1,ix,2)
         summeddx=summeddx+pw(igrid)%dx(ifirst,ix,2)
        enddo
        ! then ghost cells to left
        summeddx=0.0d0
        do ix=nghostcells,ixGextmin2,-1
          xext(ixGextmin1:ixGextmax1,ix,2)=rnode(rpxmin2_,&
             igrid)-summeddx-0.5d0*pw(igrid)%dx(ixGextmin1:ixGextmax1,ix,2)
          summeddx=summeddx+pw(igrid)%dx(ifirst,ix,2)
        enddo
       ! then ghost cells to right
       summeddx=0.0d0
       do ix=ixGhi2-nghostcells+1,ixGextmax2
          xext(ixGextmin1:ixGextmax1,ix,2)=rnode(rpxmax2_,&
             igrid)+summeddx+0.5d0*pw(igrid)%dx(ixGextmin1:ixGextmax1,ix,2)
          summeddx=summeddx+pw(igrid)%dx(ifirst,ix,2)
       enddo
      case default
        call mpistop("no such case")
    end select
   endif
  if(stretched_symm_dim(1))then
    ! here we distinguish three kinds of grid blocks
    ! depending on their ig-index, set per level 
    !      the first n_stretchedblocks/2  will stretch to the left
    !      the middle ntotal-n_stretchedblocks will be uniform
    !      the last  n_stretchedblocks/2  will stretch to the right
    if(ig1<=nstretchedblocks(level,1)/2)then
      ! stretch to the left
      offset=block_nx1*nstretchedblocks(level,1)/2
      imin=(ig1-1)*block_nx1
      imax=ig1*block_nx1
      rnode(rpxmin1_,igrid)=xprobmin1+xstretch1-dxfirst_1mq(level,&
         1) *(1.0d0-qstretch(level,1)**(offset-imin))
      rnode(rpxmax1_,igrid)=xprobmin1+xstretch1-dxfirst_1mq(level,&
         1) *(1.0d0-qstretch(level,1)**(offset-imax))
      ! fix possible out of bound due to precision
      if(rnode(rpxmin1_,igrid)<xprobmin1) rnode(rpxmin1_,igrid)=xprobmin1
      ixshift=(ig1-1)*block_nx1-nghostcells
      do ix=ixGextmin1,ixGextmax1
         index=ixshift+ix
         pw(igrid)%dx(ix,ixGlo2:ixGhi2,1)=dxfirst(level,1)*qstretch(level,&
            1)**(offset-index)
      enddo
      ixshift=(nstretchedblocks(level,1)/2-ig1)*(block_nx1/2)+block_nx1/2+&
         nghostcells
      do ix=ixCoGmin1,ixCoGmax1
         index=ixshift-ix
         pw(igrid)%dxcoarse(ix,ixCoGmin2:ixCoGmax2,1)=dxfirst(level-1,&
            1)*qstretch(level-1,1)**index
      enddo
      ! last block: to modify ghost cells!!!
      if(ig1==nstretchedblocks(level,1)/2)then
        if(ng1(level)==nstretchedblocks(level,1))then
           ! if middle blocks do not exist then use symmetry
           do ix=ixGhi1-nghostcells+1,ixGextmax1
              pw(igrid)%dx(ix,ixGlo2:ixGhi2,&
                 1)= pw(igrid)%dx(2*(ixGhi1-nghostcells)+1-ix,ixGlo2:ixGhi2,1)
           enddo
           do ix=ixCoGmax1-nghostcells+1,ixCoGmax1
              pw(igrid)%dxcoarse(ix,ixCoGmin2:ixCoGmax2,&
                 1)= pw(igrid)%dxcoarse(2*(ixCoGmax1-nghostcells)+1-ix,&
                 ixCoGmin2:ixCoGmax2,1)
           enddo
        else
           ! if middle blocks exist then use same as middle blocks: 
           do ix=ixGhi1-nghostcells+1,ixGextmax1
              pw(igrid)%dx(ix,ixGlo2:ixGhi2,1)=dxmid(level,1)
           enddo
           do ix=ixCoGmax1-nghostcells+1,ixCoGmax1
              pw(igrid)%dxcoarse(ix,ixCoGmin2:ixCoGmax2,1)=dxmid(level-1,1)
           enddo
        endif
      endif
      ! first block: make ghost cells symmetric (to allow periodicity)
      if(ig1==1)then
         do ix=ixGextmin1,nghostcells
            pw(igrid)%dx(ix,ixGextmin2:ixGextmax2,&
               1)=pw(igrid)%dx(2*nghostcells+1-ix,ixGextmin2:ixGextmax2,1)
         enddo
         do ix=1,nghostcells
            pw(igrid)%dxcoarse(ix,ixCoGmin2:ixCoGmax2,&
               1)=pw(igrid)%dxcoarse(2*nghostcells+1-ix,ixCoGmin2:ixCoGmax2,1)
         enddo
      endif
    else 
      if (ig1<=ng1(level)-nstretchedblocks(level,1)/2) then
         ! keep uniform
         pw(igrid)%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
            1)=dxmid(level,1)
         pw(igrid)%dxcoarse(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,&
            1)=dxmid(level-1,1)
         rnode(rpxmin1_,igrid)=xprobmin1+xstretch1+(ig1-nstretchedblocks(level,&
            1)/2-1)*block_nx1*dxmid(level,1)
         rnode(rpxmax1_,igrid)=xprobmin1+xstretch1+(ig1-nstretchedblocks(level,&
            1)/2)  *block_nx1*dxmid(level,1)
         ! first and last block: to modify the ghost cells!!!
         if(ig1==nstretchedblocks(level,1)/2+1)then
            do ix=ixGextmin1,nghostcells
               pw(igrid)%dx(ix,ixGextmin2:ixGextmax2,1)=dxfirst(level,&
                  1)*qstretch(level,1)**(nghostcells-ix)
            enddo
            do ix=1,nghostcells
               pw(igrid)%dxcoarse(ix,ixCoGmin2:ixCoGmax2,1)=dxfirst(level-1,&
                  1)*qstretch(level-1,1)**(nghostcells-ix)
            enddo
         endif
         if(ig1==ng1(level)-nstretchedblocks(level,1))then
            do ix=ixGhi1-nghostcells+1,ixGextmax1
              pw(igrid)%dx(ix,ixGlo2:ixGhi2,1)=dxfirst(level,1)*qstretch(level,&
                 1)**(ix-block_nx1-nghostcells-1)
            enddo
            do ix=ixCoGmax1-nghostcells+1,ixCoGmax1
              pw(igrid)%dxcoarse(ix,ixCoGmin2:ixCoGmax2,1)=dxfirst(level-1,&
                 1)*qstretch(level-1,1)**(ix-ixCoGmax1+nghostcells-1)
            enddo
         endif
      else
         ! stretch to the right
         offset=block_nx1*(ng1(level)-nstretchedblocks(level,1)/2)
         sizeuniformpart1=dxmid(1,1)*(domain_nx1-&
            nstretchedblocks_baselevel(1)*block_nx1)
         imin=(ig1-1)*block_nx1-offset
         imax=ig1*block_nx1-offset
         rnode(rpxmin1_,igrid)=xprobmin1+xstretch1+sizeuniformpart1+&
            dxfirst_1mq(level,1) *(1.0d0-qstretch(level,1)**imin)
         rnode(rpxmax1_,igrid)=xprobmin1+xstretch1+sizeuniformpart1+&
            dxfirst_1mq(level,1) *(1.0d0-qstretch(level,1)**imax)
         ! fix possible out of bound due to precision
         if(rnode(rpxmax1_,igrid)>xprobmax1) rnode(rpxmax1_,igrid)=xprobmax1
         ixshift=(ig1-1)*block_nx1-nghostcells-offset
         do ix=ixGextmin1,ixGextmax1
            index=ixshift+ix
            pw(igrid)%dx(ix,ixGextmin2:ixGextmax2,1)=dxfirst(level,&
               1)*qstretch(level,1)**(index-1)
         enddo
         ixshift=(ig1+nstretchedblocks(level,&
            1)/2-ng1(level)-1)*(block_nx1/2)-nghostcells
         do ix=ixCoGmin1,ixCoGmax1
            index=ixshift+ix
            pw(igrid)%dxcoarse(ix,ixCoGmin2:ixCoGmax2,1)=dxfirst(level-1,&
               1)*qstretch(level-1,1)**(index-1)
         enddo
         ! first block: modify ghost cells!!!
         if(ig1==ng1(level)-nstretchedblocks(level,1)+1)then
            if(ng1(level)==nstretchedblocks(level,1))then
               ! if middle blocks do not exist then use symmetry
               do ix=ixGextmin1,nghostcells
                  pw(igrid)%dx(ix,ixGextmin2:ixGextmax2,&
                     1)=pw(igrid)%dx(2*nghostcells+1-ix,ixGextmin2:ixGextmax2,&
                     1)
               enddo
               do ix=1,nghostcells
                  pw(igrid)%dxcoarse(ix,ixCoGmin2:ixCoGmax2,&
                     1)=pw(igrid)%dxcoarse(2*nghostcells+1-ix,&
                     ixCoGmin2:ixCoGmax2,1)
               enddo
            else
               ! if middle blocks exist then use same as middle blocks: 
               do ix=ixGextmin1,nghostcells
                  pw(igrid)%dx(ix,ixGextmin2:ixGextmax2,1)=dxmid(level,1)
               enddo
               do ix=1,nghostcells
                  pw(igrid)%dxcoarse(ix,ixCoGmin2:ixCoGmax2,1)=dxmid(level-1,&
                     1)
               enddo
            endif
         endif
         ! last block: make ghost cells symmetric (to allow periodicity)
         if(ig1==ng1(level))then
            do ix=ixGhi1-nghostcells+1,ixGextmax1
               pw(igrid)%dx(ix,ixGextmin2:ixGextmax2,&
                  1)=pw(igrid)%dx(2*(ixGhi1-nghostcells)+1-ix,&
                  ixGextmin2:ixGextmax2,1)
            enddo
           do ix=ixCoGmax1-nghostcells+1,ixCoGmax1
              pw(igrid)%dxcoarse(ix,ixCoGmin2:ixCoGmax2,&
                 1)=pw(igrid)%dxcoarse(2*(ixCoGmax1-nghostcells)+1-ix,&
                 ixCoGmin2:ixCoGmax2,1)
           enddo
         endif
      endif
    endif
    ! now that dx and grid boundaries are known: fill cell centers
    ifirst=nghostcells+1
    ! first fill the mesh
    summeddx=0.0d0
    do ix=ixMlo1,ixMhi1
       pw(igrid)%x(ix,ixGlo2:ixGhi2,1)=rnode(rpxmin1_,&
          igrid)+summeddx+0.5d0*pw(igrid)%dx(ix,ixGlo2:ixGhi2,1)
       summeddx=summeddx+pw(igrid)%dx(ix,ifirst,1)
    enddo
    ! then ghost cells to left
    summeddx=0.0d0
    do ix=nghostcells,1,-1
       pw(igrid)%x(ix,ixGlo2:ixGhi2,1)=rnode(rpxmin1_,&
          igrid)-summeddx-0.5d0*pw(igrid)%dx(ix,ixGlo2:ixGhi2,1)
       summeddx=summeddx+pw(igrid)%dx(ix,ifirst,1)
    enddo
    ! then ghost cells to right
    summeddx=0.0d0
    do ix=ixGhi1-nghostcells+1,ixGhi1
       pw(igrid)%x(ix,ixGlo2:ixGhi2,1)=rnode(rpxmax1_,&
          igrid)+summeddx+0.5d0*pw(igrid)%dx(ix,ixGlo2:ixGhi2,1)
       summeddx=summeddx+pw(igrid)%dx(ix,ifirst,1)
    enddo
    ! and next for the coarse representation
    ! first fill the mesh
    summeddx=0.0d0
    do ix=nghostcells+1,ixCoGmax1-nghostcells
       pw(igrid)%xcoarse(ix,ixCoGmin2:ixCoGmax2,1)=rnode(rpxmin1_,&
          igrid)+summeddx+0.5d0*pw(igrid)%dxcoarse(ix,ixCoGmin2:ixCoGmax2,1)
       summeddx=summeddx+pw(igrid)%dxcoarse(ix,ifirst,1)
    enddo
    ! then ghost cells to left
    summeddx=0.0d0
    do ix=nghostcells,1,-1
       pw(igrid)%xcoarse(ix,ixCoGmin2:ixCoGmax2,1)=rnode(rpxmin1_,&
          igrid)-summeddx-0.5d0*pw(igrid)%dxcoarse(ix,ixCoGmin2:ixCoGmax2,1)
       summeddx=summeddx+pw(igrid)%dxcoarse(ix,ifirst,1)
    enddo
    ! then ghost cells to right
    summeddx=0.0d0
    do ix=ixCoGmax1-nghostcells+1,ixCoGmax1
       pw(igrid)%xcoarse(ix,ixCoGmin2:ixCoGmax2,1)=rnode(rpxmax1_,&
          igrid)+summeddx+0.5d0*pw(igrid)%dxcoarse(ix,ixCoGmin2:ixCoGmax2,1)
       summeddx=summeddx+pw(igrid)%dxcoarse(ix,ifirst,1)
    enddo
    select case(icase)
      case(0)
        ! if even number of ghost cells: xext is just copy of local x
        xext(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
           1)=pw(igrid)%x(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,1)
      case(1)
        ! if uneven number of ghost cells: extra layer left/right
        summeddx=0.0d0
        do ix=ixMlo1,ixMhi1
          xext(ix,ixGextmin2:ixGextmax2,1)=rnode(rpxmin1_,&
             igrid)+summeddx+0.5d0*pw(igrid)%dx(ix,ixGextmin2:ixGextmax2,1)
         summeddx=summeddx+pw(igrid)%dx(ix,ifirst,1)
        enddo
        ! then ghost cells to left
        summeddx=0.0d0
        do ix=nghostcells,ixGextmin1,-1
          xext(ix,ixGextmin2:ixGextmax2,1)=rnode(rpxmin1_,&
             igrid)-summeddx-0.5d0*pw(igrid)%dx(ix,ixGextmin2:ixGextmax2,1)
          summeddx=summeddx+pw(igrid)%dx(ix,ifirst,1)
        enddo
       ! then ghost cells to right
       summeddx=0.0d0
       do ix=ixGhi1-nghostcells+1,ixGextmax1
          xext(ix,ixGextmin2:ixGextmax2,1)=rnode(rpxmax1_,&
             igrid)+summeddx+0.5d0*pw(igrid)%dx(ix,ixGextmin2:ixGextmax2,1)
          summeddx=summeddx+pw(igrid)%dx(ix,ifirst,1)
       enddo
      case default
        call mpistop("no such case")
    end select
   endif
  if(stretched_symm_dim(2))then
    ! here we distinguish three kinds of grid blocks
    ! depending on their ig-index, set per level 
    !      the first n_stretchedblocks/2  will stretch to the left
    !      the middle ntotal-n_stretchedblocks will be uniform
    !      the last  n_stretchedblocks/2  will stretch to the right
    if(ig2<=nstretchedblocks(level,2)/2)then
      ! stretch to the left
      offset=block_nx2*nstretchedblocks(level,2)/2
      imin=(ig2-1)*block_nx2
      imax=ig2*block_nx2
      rnode(rpxmin2_,igrid)=xprobmin2+xstretch2-dxfirst_1mq(level,&
         2) *(1.0d0-qstretch(level,2)**(offset-imin))
      rnode(rpxmax2_,igrid)=xprobmin2+xstretch2-dxfirst_1mq(level,&
         2) *(1.0d0-qstretch(level,2)**(offset-imax))
      ! fix possible out of bound due to precision
      if(rnode(rpxmin2_,igrid)<xprobmin2) rnode(rpxmin2_,igrid)=xprobmin2
      ixshift=(ig2-1)*block_nx2-nghostcells
      do ix=ixGextmin2,ixGextmax2
         index=ixshift+ix
         pw(igrid)%dx(ixGlo1:ixGhi1,ix,2)=dxfirst(level,2)*qstretch(level,&
            2)**(offset-index)
      enddo
      ixshift=(nstretchedblocks(level,2)/2-ig2)*(block_nx2/2)+block_nx2/2+&
         nghostcells
      do ix=ixCoGmin2,ixCoGmax2
         index=ixshift-ix
         pw(igrid)%dxcoarse(ixCoGmin1:ixCoGmax1,ix,2)=dxfirst(level-1,&
            2)*qstretch(level-1,2)**index
      enddo
      ! last block: to modify ghost cells!!!
      if(ig2==nstretchedblocks(level,2)/2)then
        if(ng2(level)==nstretchedblocks(level,2))then
           ! if middle blocks do not exist then use symmetry
           do ix=ixGhi2-nghostcells+1,ixGextmax2
              pw(igrid)%dx(ixGlo1:ixGhi1,ix,2)= pw(igrid)%dx(ixGlo1:ixGhi1,&
                 2*(ixGhi2-nghostcells)+1-ix,2)
           enddo
           do ix=ixCoGmax2-nghostcells+1,ixCoGmax2
              pw(igrid)%dxcoarse(ixCoGmin1:ixCoGmax1,ix,&
                 2)= pw(igrid)%dxcoarse(ixCoGmin1:ixCoGmax1,&
                 2*(ixCoGmax2-nghostcells)+1-ix,2)
           enddo
        else
           ! if middle blocks exist then use same as middle blocks: 
           do ix=ixGhi2-nghostcells+1,ixGextmax2
              pw(igrid)%dx(ixGlo1:ixGhi1,ix,2)=dxmid(level,2)
           enddo
           do ix=ixCoGmax2-nghostcells+1,ixCoGmax2
              pw(igrid)%dxcoarse(ixCoGmin1:ixCoGmax1,ix,2)=dxmid(level-1,2)
           enddo
        endif
      endif
      ! first block: make ghost cells symmetric (to allow periodicity)
      if(ig2==1)then
         do ix=ixGextmin2,nghostcells
            pw(igrid)%dx(ixGextmin1:ixGextmax1,ix,&
               2)=pw(igrid)%dx(ixGextmin1:ixGextmax1,2*nghostcells+1-ix,2)
         enddo
         do ix=1,nghostcells
            pw(igrid)%dxcoarse(ixCoGmin1:ixCoGmax1,ix,&
               2)=pw(igrid)%dxcoarse(ixCoGmin1:ixCoGmax1,2*nghostcells+1-ix,2)
         enddo
      endif
    else 
      if (ig2<=ng2(level)-nstretchedblocks(level,2)/2) then
         ! keep uniform
         pw(igrid)%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
            2)=dxmid(level,2)
         pw(igrid)%dxcoarse(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,&
            2)=dxmid(level-1,2)
         rnode(rpxmin2_,igrid)=xprobmin2+xstretch2+(ig2-nstretchedblocks(level,&
            2)/2-1)*block_nx2*dxmid(level,2)
         rnode(rpxmax2_,igrid)=xprobmin2+xstretch2+(ig2-nstretchedblocks(level,&
            2)/2)  *block_nx2*dxmid(level,2)
         ! first and last block: to modify the ghost cells!!!
         if(ig2==nstretchedblocks(level,2)/2+1)then
            do ix=ixGextmin2,nghostcells
               pw(igrid)%dx(ixGextmin1:ixGextmax1,ix,2)=dxfirst(level,&
                  2)*qstretch(level,2)**(nghostcells-ix)
            enddo
            do ix=1,nghostcells
               pw(igrid)%dxcoarse(ixCoGmin1:ixCoGmax1,ix,2)=dxfirst(level-1,&
                  2)*qstretch(level-1,2)**(nghostcells-ix)
            enddo
         endif
         if(ig2==ng2(level)-nstretchedblocks(level,2))then
            do ix=ixGhi2-nghostcells+1,ixGextmax2
              pw(igrid)%dx(ixGlo1:ixGhi1,ix,2)=dxfirst(level,2)*qstretch(level,&
                 2)**(ix-block_nx2-nghostcells-1)
            enddo
            do ix=ixCoGmax2-nghostcells+1,ixCoGmax2
              pw(igrid)%dxcoarse(ixCoGmin1:ixCoGmax1,ix,2)=dxfirst(level-1,&
                 2)*qstretch(level-1,2)**(ix-ixCoGmax2+nghostcells-1)
            enddo
         endif
      else
         ! stretch to the right
         offset=block_nx2*(ng2(level)-nstretchedblocks(level,2)/2)
         sizeuniformpart2=dxmid(1,2)*(domain_nx2-&
            nstretchedblocks_baselevel(2)*block_nx2)
         imin=(ig2-1)*block_nx2-offset
         imax=ig2*block_nx2-offset
         rnode(rpxmin2_,igrid)=xprobmin2+xstretch2+sizeuniformpart2+&
            dxfirst_1mq(level,2) *(1.0d0-qstretch(level,2)**imin)
         rnode(rpxmax2_,igrid)=xprobmin2+xstretch2+sizeuniformpart2+&
            dxfirst_1mq(level,2) *(1.0d0-qstretch(level,2)**imax)
         ! fix possible out of bound due to precision
         if(rnode(rpxmax2_,igrid)>xprobmax2) rnode(rpxmax2_,igrid)=xprobmax2
         ixshift=(ig2-1)*block_nx2-nghostcells-offset
         do ix=ixGextmin2,ixGextmax2
            index=ixshift+ix
            pw(igrid)%dx(ixGextmin1:ixGextmax1,ix,2)=dxfirst(level,&
               2)*qstretch(level,2)**(index-1)
         enddo
         ixshift=(ig2+nstretchedblocks(level,&
            2)/2-ng2(level)-1)*(block_nx2/2)-nghostcells
         do ix=ixCoGmin2,ixCoGmax2
            index=ixshift+ix
            pw(igrid)%dxcoarse(ixCoGmin1:ixCoGmax1,ix,2)=dxfirst(level-1,&
               2)*qstretch(level-1,2)**(index-1)
         enddo
         ! first block: modify ghost cells!!!
         if(ig2==ng2(level)-nstretchedblocks(level,2)+1)then
            if(ng2(level)==nstretchedblocks(level,2))then
               ! if middle blocks do not exist then use symmetry
               do ix=ixGextmin2,nghostcells
                  pw(igrid)%dx(ixGextmin1:ixGextmax1,ix,&
                     2)=pw(igrid)%dx(ixGextmin1:ixGextmax1,2*nghostcells+1-ix,&
                     2)
               enddo
               do ix=1,nghostcells
                  pw(igrid)%dxcoarse(ixCoGmin1:ixCoGmax1,ix,&
                     2)=pw(igrid)%dxcoarse(ixCoGmin1:ixCoGmax1,&
                     2*nghostcells+1-ix,2)
               enddo
            else
               ! if middle blocks exist then use same as middle blocks: 
               do ix=ixGextmin2,nghostcells
                  pw(igrid)%dx(ixGextmin1:ixGextmax1,ix,2)=dxmid(level,2)
               enddo
               do ix=1,nghostcells
                  pw(igrid)%dxcoarse(ixCoGmin1:ixCoGmax1,ix,2)=dxmid(level-1,&
                     2)
               enddo
            endif
         endif
         ! last block: make ghost cells symmetric (to allow periodicity)
         if(ig2==ng2(level))then
            do ix=ixGhi2-nghostcells+1,ixGextmax2
               pw(igrid)%dx(ixGextmin1:ixGextmax1,ix,&
                  2)=pw(igrid)%dx(ixGextmin1:ixGextmax1,&
                  2*(ixGhi2-nghostcells)+1-ix,2)
            enddo
           do ix=ixCoGmax2-nghostcells+1,ixCoGmax2
              pw(igrid)%dxcoarse(ixCoGmin1:ixCoGmax1,ix,&
                 2)=pw(igrid)%dxcoarse(ixCoGmin1:ixCoGmax1,&
                 2*(ixCoGmax2-nghostcells)+1-ix,2)
           enddo
         endif
      endif
    endif
    ! now that dx and grid boundaries are known: fill cell centers
    ifirst=nghostcells+1
    ! first fill the mesh
    summeddx=0.0d0
    do ix=ixMlo2,ixMhi2
       pw(igrid)%x(ixGlo1:ixGhi1,ix,2)=rnode(rpxmin2_,&
          igrid)+summeddx+0.5d0*pw(igrid)%dx(ixGlo1:ixGhi1,ix,2)
       summeddx=summeddx+pw(igrid)%dx(ifirst,ix,2)
    enddo
    ! then ghost cells to left
    summeddx=0.0d0
    do ix=nghostcells,1,-1
       pw(igrid)%x(ixGlo1:ixGhi1,ix,2)=rnode(rpxmin2_,&
          igrid)-summeddx-0.5d0*pw(igrid)%dx(ixGlo1:ixGhi1,ix,2)
       summeddx=summeddx+pw(igrid)%dx(ifirst,ix,2)
    enddo
    ! then ghost cells to right
    summeddx=0.0d0
    do ix=ixGhi2-nghostcells+1,ixGhi2
       pw(igrid)%x(ixGlo1:ixGhi1,ix,2)=rnode(rpxmax2_,&
          igrid)+summeddx+0.5d0*pw(igrid)%dx(ixGlo1:ixGhi1,ix,2)
       summeddx=summeddx+pw(igrid)%dx(ifirst,ix,2)
    enddo
    ! and next for the coarse representation
    ! first fill the mesh
    summeddx=0.0d0
    do ix=nghostcells+1,ixCoGmax2-nghostcells
       pw(igrid)%xcoarse(ixCoGmin1:ixCoGmax1,ix,2)=rnode(rpxmin2_,&
          igrid)+summeddx+0.5d0*pw(igrid)%dxcoarse(ixCoGmin1:ixCoGmax1,ix,2)
       summeddx=summeddx+pw(igrid)%dxcoarse(ifirst,ix,2)
    enddo
    ! then ghost cells to left
    summeddx=0.0d0
    do ix=nghostcells,1,-1
       pw(igrid)%xcoarse(ixCoGmin1:ixCoGmax1,ix,2)=rnode(rpxmin2_,&
          igrid)-summeddx-0.5d0*pw(igrid)%dxcoarse(ixCoGmin1:ixCoGmax1,ix,2)
       summeddx=summeddx+pw(igrid)%dxcoarse(ifirst,ix,2)
    enddo
    ! then ghost cells to right
    summeddx=0.0d0
    do ix=ixCoGmax2-nghostcells+1,ixCoGmax2
       pw(igrid)%xcoarse(ixCoGmin1:ixCoGmax1,ix,2)=rnode(rpxmax2_,&
          igrid)+summeddx+0.5d0*pw(igrid)%dxcoarse(ixCoGmin1:ixCoGmax1,ix,2)
       summeddx=summeddx+pw(igrid)%dxcoarse(ifirst,ix,2)
    enddo
    select case(icase)
      case(0)
        ! if even number of ghost cells: xext is just copy of local x
        xext(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
           2)=pw(igrid)%x(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,2)
      case(1)
        ! if uneven number of ghost cells: extra layer left/right
        summeddx=0.0d0
        do ix=ixMlo2,ixMhi2
          xext(ixGextmin1:ixGextmax1,ix,2)=rnode(rpxmin2_,&
             igrid)+summeddx+0.5d0*pw(igrid)%dx(ixGextmin1:ixGextmax1,ix,2)
         summeddx=summeddx+pw(igrid)%dx(ifirst,ix,2)
        enddo
        ! then ghost cells to left
        summeddx=0.0d0
        do ix=nghostcells,ixGextmin2,-1
          xext(ixGextmin1:ixGextmax1,ix,2)=rnode(rpxmin2_,&
             igrid)-summeddx-0.5d0*pw(igrid)%dx(ixGextmin1:ixGextmax1,ix,2)
          summeddx=summeddx+pw(igrid)%dx(ifirst,ix,2)
        enddo
       ! then ghost cells to right
       summeddx=0.0d0
       do ix=ixGhi2-nghostcells+1,ixGextmax2
          xext(ixGextmin1:ixGextmax1,ix,2)=rnode(rpxmax2_,&
             igrid)+summeddx+0.5d0*pw(igrid)%dx(ixGextmin1:ixGextmax1,ix,2)
          summeddx=summeddx+pw(igrid)%dx(ifirst,ix,2)
       enddo
      case default
        call mpistop("no such case")
    end select
   endif
endif

if (.not.slab) then
   call fillgeo(igrid,ixGlo1,ixGlo2,ixGhi1,ixGhi2)
   select case (typeaxial)
      case ("slabstretch")
         pw(igrid)%dvolume(ixGextmin1:ixGextmax1,&
            ixGextmin2:ixGextmax2)= pw(igrid)%dx(ixGextmin1:ixGextmax1,&
            ixGextmin2:ixGextmax2,1)*pw(igrid)%dx(ixGextmin1:ixGextmax1,&
            ixGextmin2:ixGextmax2,2)
         pw(igrid)%dvolumecoarse(ixCoGmin1:ixCoGmax1,&
            ixCoGmin2:ixCoGmax2)= pw(igrid)%dxcoarse(ixCoGmin1:ixCoGmax1,&
            ixCoGmin2:ixCoGmax2,1)*pw(igrid)%dxcoarse(ixCoGmin1:ixCoGmax1,&
            ixCoGmin2:ixCoGmax2,2)
         pw(igrid)%ds(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
            1:ndim)=pw(igrid)%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
            1:ndim)
      case ("spherical")
         pw(igrid)%dvolume(ixGextmin1:ixGextmax1,&
            ixGextmin2:ixGextmax2)=(xext(ixGextmin1:ixGextmax1,&
            ixGextmin2:ixGextmax2,1)**2 +pw(igrid)%dx(ixGextmin1:ixGextmax1,&
            ixGextmin2:ixGextmax2,1)**2/12.0d0)*pw(igrid)%dx(&
            ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
            1) *two*dabs(dsin(xext(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
            2))) *dsin(half*pw(igrid)%dx(ixGextmin1:ixGextmax1,&
            ixGextmin2:ixGextmax2,2))
         pw(igrid)%dvolumecoarse(ixCoGmin1:ixCoGmax1,&
            ixCoGmin2:ixCoGmax2)=(pw(igrid)%xcoarse(ixCoGmin1:ixCoGmax1,&
            ixCoGmin2:ixCoGmax2,1)**2 +pw(igrid)%dxcoarse(ixCoGmin1:ixCoGmax1,&
            ixCoGmin2:ixCoGmax2,1)**2/12.0d0)*pw(igrid)%dxcoarse(&
            ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,&
            1) *two*dabs(dsin(pw(igrid)%xcoarse(ixCoGmin1:ixCoGmax1,&
            ixCoGmin2:ixCoGmax2,2))) *dsin(half*pw(igrid)%dxcoarse(&
            ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,2))
         pw(igrid)%ds(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
            1)=pw(igrid)%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,1)
          pw(igrid)%ds(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
             2)=xext(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
             1)*pw(igrid)%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,2)
         
      case ("cylindrical")
         pw(igrid)%dvolume(ixGextmin1:ixGextmax1,&
            ixGextmin2:ixGextmax2)=dabs(xext(ixGextmin1:ixGextmax1,&
            ixGextmin2:ixGextmax2,1)) *pw(igrid)%dx(ixGextmin1:ixGextmax1,&
            ixGextmin2:ixGextmax2,1)*pw(igrid)%dx(ixGextmin1:ixGextmax1,&
            ixGextmin2:ixGextmax2,2)
         pw(igrid)%dvolumecoarse(ixCoGmin1:ixCoGmax1,&
            ixCoGmin2:ixCoGmax2)=dabs(pw(igrid)%xcoarse(ixCoGmin1:ixCoGmax1,&
            ixCoGmin2:ixCoGmax2,1)) *pw(igrid)%dxcoarse(ixCoGmin1:ixCoGmax1,&
            ixCoGmin2:ixCoGmax2,1)*pw(igrid)%dxcoarse(ixCoGmin1:ixCoGmax1,&
            ixCoGmin2:ixCoGmax2,2)
         pw(igrid)%ds(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
            r_)=pw(igrid)%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,r_)
         if(z_>0.and.z_<=ndim) pw(igrid)%ds(ixGextmin1:ixGextmax1,&
            ixGextmin2:ixGextmax2,z_)=pw(igrid)%dx(ixGextmin1:ixGextmax1,&
            ixGextmin2:ixGextmax2,z_)
         if (phi_ > 0) then
           if (2==phi_) pw(igrid)%ds(ixGextmin1:ixGextmax1,&
              ixGextmin2:ixGextmax2,2)= xext(ixGextmin1:ixGextmax1,&
              ixGextmin2:ixGextmax2,1)*pw(igrid)%dx(ixGextmin1:ixGextmax1,&
              ixGextmin2:ixGextmax2,2)
         end if
      case default
         call mpistop("Sorry, typeaxial unknown")
      end select
endif

if (B0field) then
   ! initialize background non-evolving solution
   call alloc_B0_grid(igrid)
   call set_B0_grid(igrid,global_time)
end if

! find the blocks on the boundaries
pw(igrid)%is_physical_boundary=.false.

do i1=-1,1
  if(i1==0) cycle
  ign1=ig1+i1
  ! blocks at periodic boundary have neighbors in the physical domain
  ! thus threated at internal blocks with no physical boundary
  if (periodB(1)) ign1=1+modulo(ign1-1,ng1(level))
  if (ign1 > ng1(level)) then
     if(phi_ > 0 .and. poleB(2,1)) then
       ! if at a pole, the boundary is not physical boundary
       pw(igrid)%is_physical_boundary(2*1)=.false.
     else
       pw(igrid)%is_physical_boundary(2*1)=.true.
     end if
  else if (ign1 < 1) then
     if(phi_ > 0 .and. poleB(1,1)) then
       ! if at a pole, the boundary is not physical boundary
       pw(igrid)%is_physical_boundary(2*1-1)=.false.
     else
       pw(igrid)%is_physical_boundary(2*1-1)=.true.
     end if
  end if
end do


do i2=-1,1
  if(i2==0) cycle
  ign2=ig2+i2
  ! blocks at periodic boundary have neighbors in the physical domain
  ! thus threated at internal blocks with no physical boundary
  if (periodB(2)) ign2=1+modulo(ign2-1,ng2(level))
  if (ign2 > ng2(level)) then
     if(phi_ > 0 .and. poleB(2,2)) then
       ! if at a pole, the boundary is not physical boundary
       pw(igrid)%is_physical_boundary(2*2)=.false.
     else
       pw(igrid)%is_physical_boundary(2*2)=.true.
     end if
  else if (ign2 < 1) then
     if(phi_ > 0 .and. poleB(1,2)) then
       ! if at a pole, the boundary is not physical boundary
       pw(igrid)%is_physical_boundary(2*2-1)=.false.
     else
       pw(igrid)%is_physical_boundary(2*2-1)=.true.
     end if
  end if
end do

if(any(pw(igrid)%is_physical_boundary)) then
  phyboundblock(igrid)=.true.
else
  phyboundblock(igrid)=.false.
end if

end subroutine alloc_node
!=============================================================================
subroutine dealloc_node(igrid)

use mod_global_parameters
use mod_geometry

integer, intent(in) :: igrid
!-----------------------------------------------------------------------------
if (igrid==0) then
   call mpistop("trying to delete a non-existing grid in dealloc_node")
end if

deallocate(pw(igrid)%w,pw(igrid)%w1,pw(igrid)%x)
deallocate(pw(igrid)%wcoarse,pw(igrid)%xcoarse)
deallocate(pw(igrid)%wold,pw(igrid)%wio)
! deallocate temperary solution space
select case (time_integrator)
case("threestep","fourstep","jameson","twostep_trapezoidal")
  deallocate(pw(igrid)%w2)
case("rk4","ssprk43")
  deallocate(pw(igrid)%w2)
  deallocate(pw(igrid)%w3)
case("ssprk54")
  deallocate(pw(igrid)%w2)
  deallocate(pw(igrid)%w3)
  deallocate(pw(igrid)%w4)
end select
if(allocated(pw(igrid)%w2)) deallocate(pw(igrid)%w2)
if(allocated(pw(igrid)%w3)) deallocate(pw(igrid)%w3)

if (.not.slab) call putgridgeo(igrid)

if (B0field) call dealloc_B0_grid(igrid)

end subroutine dealloc_node
!=============================================================================
