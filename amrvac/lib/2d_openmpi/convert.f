!=============================================================================
subroutine generate_plotfile
use mod_usr_methods, only: usr_special_convert,usr_special_global
use mod_global_parameters
use mod_ghostcells_update
use mod_physics, only: phys_req_diagonal
!-----------------------------------------------------------------------------

if(mype==0.and.level_io>0) write(unitterm,*)'reset tree to fixed level=',&
   level_io
if(level_io>0 .or. level_io_min.ne.1 .or. level_io_max.ne.nlevelshi) then
   call resettree_convert
else if(.not. phys_req_diagonal) then
   call getbc(global_time,0.d0,0,nwflux+nwaux)
end if
 call usr_special_global

select case(convert_type)
  case('tecplot','tecplotCC','tecline')
   call tecplot(unitconvert)
  case('tecplotmpi','tecplotCCmpi','teclinempi')
   call tecplot_mpi(unitconvert)
  case('vtu','vtuCC')
   call unstructuredvtk(unitconvert)
  case('vtumpi','vtuCCmpi')
   call unstructuredvtk_mpi(unitconvert)
  case('vtuB','vtuBCC','vtuBmpi','vtuBCCmpi')
   call unstructuredvtkB(unitconvert)
  case('pvtumpi','pvtuCCmpi')
   call punstructuredvtk_mpi(unitconvert)
  case('pvtuBmpi','pvtuBCCmpi')
   call punstructuredvtkB_mpi(unitconvert)
  case('vtimpi','vtiCCmpi')
   call ImageDataVtk_mpi(unitconvert)
  case('onegrid','onegridmpi')
   call onegrid(unitconvert)
  case('oneblock','oneblockB')
   call oneblock(unitconvert)
  case('user','usermpi')
     if (.not. associated(usr_special_convert)) then
        call mpistop("usr_special_convert not defined")
     else
        call usr_special_convert(unitconvert)
     end if
  case default
   call mpistop("Error in generate_plotfile: Unknown convert_type")
end select

end subroutine generate_plotfile
!=============================================================================
subroutine calc_grid(qunit,igrid,xC,xCC,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,normconv,&
   ixCmin1,ixCmin2,ixCmax1,ixCmax2,ixCCmin1,ixCCmin2,ixCCmax1,ixCCmax2,first)

! this subroutine computes both corner as well as cell-centered values
! it handles how we do the center to corner averaging, as well as
! whether we switch to cartesian or want primitive or conservative output,
! handling the addition of B0 in B0+B1 cases, ...
!
! the normconv is passed on to usr_aux_output for extending with
! possible normalization values for the nw+1:nw+nwauxio entries
use mod_usr_methods, only: usr_aux_output
use mod_global_parameters
use mod_limiter
use mod_physics, only: physics_type, phys_to_primitive

integer, intent(in) :: qunit, igrid
double precision, intent(in), dimension(ixMlo1-1:ixMhi1,ixMlo2-1:ixMhi2,&
   ndim) :: xC
double precision, intent(in), dimension(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
   ndim)   :: xCC
integer :: ixCmin1,ixCmin2,ixCmax1,ixCmax2,ixCCmin1,ixCCmin2,ixCCmax1,ixCCmax2
logical, intent(in) :: first

double precision, dimension(ixMlo1-1:ixMhi1,ixMlo2-1:ixMhi2,ndim) :: xC_TMP
double precision, dimension(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ndim)   :: xCC_TMP
double precision, dimension(ixMlo1-1:ixMhi1,ixMlo2-1:ixMhi2,&
   nw+nwauxio)   :: wC_TMP
double precision, dimension(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
   nw+nwauxio)     :: wCC_TMP
double precision,dimension(0:nw+nwauxio),intent(out)       :: normconv

double precision :: ldw(ixGlo1:ixGhi1,ixGlo2:ixGhi2), dwC(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2)
double precision, dimension(ixMlo1-1:ixMhi1,ixMlo2-1:ixMhi2,&
   nw+nwauxio)   :: wC
double precision, dimension(ixMlo1:ixMhi1,ixMlo2:ixMhi2,nw+nwauxio)     :: wCC
double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nw+nwauxio)   :: w
double precision :: dx1,dx2
integer :: nxCC1,nxCC2,idims,jxCmin1,jxCmin2,jxCmax1,jxCmax2,iwe
integer :: nx1,nx2, nxC1,nxC2, ix1,ix2, ix, iw, level, idir
logical, save :: subfirst=.true.
!-----------------------------------------------------------------------------


ixCmin1=ixMlo1-1;ixCmin2=ixMlo2-1; ixCmax1=ixMhi1;ixCmax2=ixMhi2; !Corner indices
ixCCmin1=ixMlo1;ixCCmin2=ixMlo2; ixCCmax1=ixMhi1;ixCCmax2=ixMhi2; !Center indices

saveigrid=igrid
nx1=ixMhi1-ixMlo1+1;nx2=ixMhi2-ixMlo2+1;
level=node(plevel_,igrid)
dx1=dx(1,level);dx2=dx(2,level);

normconv(0) = length_convert_factor
normconv(1:nw) = w_convert_factor

w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nw)=pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
   1:nw)

if (nwextra>0) then
 ! here we actually fill the ghost layers for the nwextra variables using
 ! continuous extrapolation (as these values do not exist normally in ghost cells)
 do idims=1,ndim
  select case(idims)
   case(1)
     jxCmin1=ixGhi1+1-nghostcells;jxCmin2=ixGlo2;
     jxCmax1=ixGhi1;jxCmax2=ixGhi2;
     do ix1=jxCmin1,jxCmax1
         w(ix1,jxCmin2:jxCmax2,nw-nwextra+1:nw) = w(jxCmin1-1,jxCmin2:jxCmax2,&
            nw-nwextra+1:nw)
     end do
     jxCmin1=ixGlo1;jxCmin2=ixGlo2;
     jxCmax1=ixGlo1-1+nghostcells;jxCmax2=ixGhi2;
     do ix1=jxCmin1,jxCmax1
         w(ix1,jxCmin2:jxCmax2,nw-nwextra+1:nw) = w(jxCmax1+1,jxCmin2:jxCmax2,&
            nw-nwextra+1:nw)
     end do 
   case(2)
     jxCmin1=ixGlo1;jxCmin2=ixGhi2+1-nghostcells;
     jxCmax1=ixGhi1;jxCmax2=ixGhi2;
     do ix2=jxCmin2,jxCmax2
         w(jxCmin1:jxCmax1,ix2,nw-nwextra+1:nw) = w(jxCmin1:jxCmax1,jxCmin2-1,&
            nw-nwextra+1:nw)
     end do
     jxCmin1=ixGlo1;jxCmin2=ixGlo2;
     jxCmax1=ixGhi1;jxCmax2=ixGlo2-1+nghostcells;
     do ix2=jxCmin2,jxCmax2
         w(jxCmin1:jxCmax1,ix2,nw-nwextra+1:nw) = w(jxCmin1:jxCmax1,jxCmax2+1,&
            nw-nwextra+1:nw)
     end do 
  end select
 end do
end if

! next lines needed when usr_aux_output uses gradients
! and later on when dwlimiter2 is used
typelimiter=type_limiter(node(plevel_,igrid))
typegradlimiter=type_gradient_limiter(node(plevel_,igrid))
dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid);
block=>pw(igrid)
if(nwauxio>0)then
  ! auxiliary io variables can be computed and added by user
  ! next few lines ensure correct usage of routines like divvector etc
  ! default (no) normalization for auxiliary variables
  normconv(nw+1:nw+nwauxio)=one

  if (.not. associated(usr_aux_output)) then
     call mpistop("usr_aux_output not defined")
  else
     call usr_aux_output(ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixMlo1-1,ixMlo2-1,&
        ixMhi1+1,ixMhi2+1,w,pw(igrid)%x,normconv)
  end if
endif

! In case primitives to be saved: use primitive subroutine
!  extra layer around mesh only needed when storing corner values and averaging
if(saveprim.and.first) call phys_to_primitive(ixGlo1,ixGlo2,ixGhi1,ixGhi2,&
   ixMlo1-1,ixMlo2-1,ixMhi1+1,ixMhi2+1,w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nw),&
   pw(igrid)%x)


if(allocated(pw(igrid)%B0)) then
! B0+B1 split handled here
  if(.not.saveprim.and.phys_energy) then
    w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,iw_e)=w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       iw_e)+0.5d0*sum(pw(igrid)%B0(ixGlo1:ixGhi1,ixGlo2:ixGhi2,:,0)**2,&
       dim=ndim+1) + sum(w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       iw_mag(:))*pw(igrid)%B0(ixGlo1:ixGhi1,ixGlo2:ixGhi2,:,0),dim=ndim+1)
  end if
  w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,iw_mag(:))=w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
     iw_mag(:))+pw(igrid)%B0(ixGlo1:ixGhi1,ixGlo2:ixGhi2,:,0)
end if
! compute the cell-center values for w first
! cell center values obtained from mere copy
wCC(ixCCmin1:ixCCmax1,ixCCmin2:ixCCmax2,:)=w(ixCCmin1:ixCCmax1,&
   ixCCmin2:ixCCmax2,:)

! compute the corner values for w now by averaging

if(slab) then
   ! for slab symmetry: no geometrical info required
   do iw=1,nw+nwauxio
     do ix2=ixCmin2,ixCmax2
     do ix1=ixCmin1,ixCmax1
        wC(ix1,ix2,iw)=sum(w(ix1:ix1+1,ix2:ix2+1,iw))/dble(2**ndim)
     end do
     end do
   end do
else
   do iw=1,nw+nwauxio
     do ix2=ixCmin2,ixCmax2
     do ix1=ixCmin1,ixCmax1
       wC(ix1,ix2,iw)=sum(w(ix1:ix1+1,ix2:ix2+1,&
          iw)*pw(igrid)%dvolume(ix1:ix1+1,&
          ix2:ix2+1)) /sum(pw(igrid)%dvolume(ix1:ix1+1,ix2:ix2+1))
     end do
     end do
   end do
endif

if(nocartesian) then
  ! keep the coordinate and vector components
  xC_TMP(ixCmin1:ixCmax1,ixCmin2:ixCmax2,1:ndim)          = xC(ixCmin1:ixCmax1,&
     ixCmin2:ixCmax2,1:ndim)
  wC_TMP(ixCmin1:ixCmax1,ixCmin2:ixCmax2,1:nw+nwauxio)    = wC(ixCmin1:ixCmax1,&
     ixCmin2:ixCmax2,1:nw+nwauxio)
  xCC_TMP(ixCCmin1:ixCCmax1,ixCCmin2:ixCCmax2,&
     1:ndim)        = xCC(ixCCmin1:ixCCmax1,ixCCmin2:ixCCmax2,1:ndim)
  wCC_TMP(ixCCmin1:ixCCmax1,ixCCmin2:ixCCmax2,&
     1:nw+nwauxio)  = wCC(ixCCmin1:ixCCmax1,ixCCmin2:ixCCmax2,1:nw+nwauxio)
else
  ! do all conversions to cartesian coordinates and vector components
  ! start for the corner values
  call cartesian(xC_TMP,wC_TMP,ixCmin1,ixCmin2,ixCmax1,ixCmax2,xC,wC)
  ! then cell center values
  call cartesian(xCC_TMP,wCC_TMP,ixCCmin1,ixCCmin2,ixCCmax1,ixCCmax2,xCC,wCC)
endif



! Warning: differentiate between idl/idlCC/tecplot/tecplotCC/vtu(B)/vtu(B)CC
if(nwaux>0 .and. mype==0 .and. first.and.subfirst) then
  ! when corner values are computed and auxiliaries present: warn!
  if(convert_type=='idl'.or.convert_type=='tecplot' &
     .or.convert_type=='vtu'.or.convert_type=='vtuB') write(*,&
     *) 'Warning: also averaged auxiliaries within calc_grid'
  subfirst=.false.
endif

end subroutine calc_grid
!=============================================================================
subroutine cartesian(x_TMP,w_TMP,ixmin1,ixmin2,ixmax1,ixmax2,xC,wC)

! conversion of coordinate and vector components from cylindrical/spherical
! to cartesian coordinates and components done here
! Also: nullifying values lower than smalldouble

use mod_global_parameters

integer :: ixmin1,ixmin2,ixmax1,ixmax2, ix1,ix2, idim, iw, ivector, iw0
integer, dimension(nw) :: vectoriw
double precision :: x_TEC(ndim), w_TEC(nw+nwauxio)
double precision, dimension(ndim,ndim) :: normal

double precision, dimension(ixmin1:ixmax1,ixmin2:ixmax2,ndim) :: xC
double precision, dimension(ixmin1:ixmax1,ixmin2:ixmax2,nw+nwauxio)   :: wC

double precision, dimension(ixmin1:ixmax1,ixmin2:ixmax2,ndim) :: x_TMP
double precision, dimension(ixmin1:ixmax1,ixmin2:ixmax2,nw+nwauxio)   :: w_TMP
!-----------------------------------------------------------------------------

iw0=0
vectoriw=-1
if(nvector>0) then
  do ivector=1,nvector
     do idim=1,ndim
        vectoriw(iw_vector(ivector)+idim)=iw_vector(ivector)
     end do
  end do
endif
x_TEC=0.d0
do ix2=ixmin2,ixmax2
do ix1=ixmin1,ixmax1
   select case (typeaxial)
   case ("slab","slabstretch")
      x_TEC(1:ndim)=xC(ix1,ix2,1:ndim)
      w_TEC(1:nw+nwauxio)=wC(ix1,ix2,1:nw+nwauxio)
   case ("cylindrical")
      
      
      select case (phi_)
      case (2)
         x_TEC(1)=xC(ix1,ix2,1)*cos(xC(ix1,ix2,2))
         x_TEC(2)=xC(ix1,ix2,1)*sin(xC(ix1,ix2,2))
      case default
         x_TEC(1)=xC(ix1,ix2,1)
         x_TEC(2)=xC(ix1,ix2,2)
      end select
      

      if (nvector>0) then
         

         
         select case (phi_)
         case (2)
            normal(1,1)=cos(xC(ix1,ix2,2))
            normal(1,2)=-sin(xC(ix1,ix2,2))
            normal(2,1)=sin(xC(ix1,ix2,2))
            normal(2,2)=cos(xC(ix1,ix2,2))
         case default
            normal(1,1)=one
            normal(1,2)=zero
            normal(2,1)=zero
            normal(2,2)=one
         end select

         
      end if
      do iw=1,nw+nwauxio
         if (iw<=nw) iw0=vectoriw(iw)
         if (iw0>=0.and.iw<=iw0+ndim.and.iw<=nw) then
            idim=iw-iw0
            w_TEC(iw0+idim)=wC(ix1,ix2,iw0+1)*normal(idim,1)+wC(ix1,ix2,&
               iw0+2)*normal(idim,2)
         else
            w_TEC(iw)=wC(ix1,ix2,iw)
         end if
      end do
   case ("spherical")
      x_TEC(1)=xC(ix1,ix2,1)*sin(xC(ix1,ix2,2))
      
      x_TEC(2)=xC(ix1,ix2,1)*cos(xC(ix1,ix2,2))
      

      if (nvector>0) then
         
         
         normal(1,1)=sin(xC(ix1,ix2,2))
         normal(1,2)=cos(xC(ix1,ix2,2))

         
         normal(2,1)=cos(xC(ix1,ix2,2))
         normal(2,2)=-sin(xC(ix1,ix2,2))
         
      end if
      do iw=1,nw+nwauxio
         if (iw<=nw) iw0=vectoriw(iw)
         if (iw0>=0.and.iw<=iw0+ndim.and.iw<=nw) then
            idim=iw-iw0
            w_TEC(iw0+idim)=wC(ix1,ix2,iw0+1)*normal(idim,1)+wC(ix1,ix2,&
               iw0+2)*normal(idim,2)
         else
            w_TEC(iw)=wC(ix1,ix2,iw)
         end if
      end do
   case default
      write(*,*) "No converter for typeaxial=",typeaxial
   end select
   x_TMP(ix1,ix2,1:ndim)=x_TEC(1:ndim)
   w_TMP(ix1,ix2,1:nw+nwauxio)=w_TEC(1:nw+nwauxio)
   ! Be aware that small values are nullified here!!!
   where(dabs(w_TMP(ix1,ix2,1:nw+nwauxio))<smalldouble)
         w_TMP(ix1,ix2,1:nw+nwauxio)=zero
   endwhere
end do
end do

end subroutine cartesian
!=============================================================================
subroutine getheadernames(wnamei,xandwnamei,outfilehead)

! this collects all variables names in the wnamei character array, getting the info from
! the prim_wnames/cons_wnames strings (depending on saveprim). It combines this info with names
! for the dimensional directions in the xandwnamei array. In the outfilehead, it collects
! the dimensional names, and only those names from the nw variables for output (through w_write)
! together with the added names for nwauxio variables

  use mod_usr_methods, only: usr_add_aux_names
  use mod_global_parameters

character(len=name_len)   :: wnamei(1:nw+nwauxio),&
   xandwnamei(1:ndim+nw+nwauxio)
character(len=1024) :: outfilehead

integer                           ::  space_position,iw
character(len=name_len)           ::  wname
character(len=std_len),allocatable:: aux_variable_names(:)
!character(len=std_len)::  scanstring

logical, save:: first=.true.
!-----------------------------------------------------------------------------

! in case additional variables are computed and stored for output
if (nwauxio>0) then
   if (.not. associated(usr_add_aux_names)) then
      call mpistop("usr_add_aux_names not defined")
   else
      if(nwauxio>0)allocate(aux_variable_names(1:nwauxio))
      call usr_add_aux_names(aux_variable_names)
   end if
end if

! --- part to provide variable names from prim_wnames/cons_wnames strings
if(saveprim) then
  ! scanstring=TRIM(aux_variable_names)
   wnamei(1:nw)=prim_wnames(1:nw)
else
  ! scanstring=TRIM(aux_variable_names)
   wnamei(1:nw)=cons_wnames(1:nw)
endif

! space_position=index(scanstring,' ')
! do iw=nw+1,nw+nwauxio
!    do while (space_position==1)
!      scanstring=scanstring(2:)
!      space_position=index(scanstring,' ')
!    enddo
!    wname=scanstring(:space_position-1)
!    scanstring=scanstring(space_position+1:)
!    space_position=index(scanstring,' ')
!
!    ! fill all names, even those that we will not write away from the first nw
!    wnamei(iw)=TRIM(wname)
! enddo
if(nwauxio>0) then
  do iw=nw+1,nw+nwauxio
    wnamei(iw)=TRIM(aux_variable_names(iw-nw))
  end do
end if
! --- end of part to provide variable names

select case (typeaxial)
   case( "spherical" )
      xandwnamei(1)="r"; xandwnamei(2)="Theta";
   case( "cylindrical" )
      xandwnamei(1)="R";
      
      if( phi_ == 2 )then
         xandwnamei(2)="Phi"
      else
         xandwnamei(2)="Z"
      endif
      
   case default
      xandwnamei(1)="X"; xandwnamei(2)="Y";
end select

xandwnamei(ndim+1:ndim+nw+nwauxio)=wnamei(1:nw+nwauxio)

! in outfilehead, collect the dimensional names, and all output variable names
! first all dimensions
write(outfilehead,'(a)') TRIM(xandwnamei(1))

do iw=2,ndim
   wname=xandwnamei(iw)
write(outfilehead,'(a)')outfilehead(1:len_trim(outfilehead))//" "//TRIM(wname)
enddo

! then all nw variables, with w_write control for inclusion
do iw=ndim+1,ndim+nw
   wname=xandwnamei(iw)
   if(w_write(iw-ndim)) then
write(outfilehead,'(a)')outfilehead(1:len_trim(outfilehead))//" "//TRIM(wname)
   endif
enddo
! then all nwauxio variables
if(nwauxio>0) then
  do iw=ndim+nw+1,ndim+nw+nwauxio
     wname=xandwnamei(iw)
write(outfilehead,'(a)')outfilehead(1:len_trim(outfilehead))//" "//TRIM(wname)
  enddo
endif

if(first.and.mype==0)then
  print*,'-------------------------------------------------------------------------------'
  write(unitterm,*&
     )'Saving visual data. Coordinate directions and variable names are:'
  do iw=1,ndim
    print *,iw,xandwnamei(iw)
  enddo
  do iw=ndim+1,ndim+nw+nwauxio
    print *,iw,wnamei(iw-ndim),xandwnamei(iw)
  enddo
  write(unitterm,*)'time =', global_time
  print*,'-------------------------------------------------------------------------------'
  first=.false.
endif

end subroutine getheadernames
!=============================================================================
subroutine oneblock(qunit)
! this is for turning an AMR run into a single block
! the data will be all on selected level level_io

! this version should work for any dimension
! only writes w_write selected 1:nw variables, also nwauxio
! may use saveprim to switch to primitives
! this version can not work on multiple CPUs
! does not renormalize variables

! header info differs from onegrid below

! ASCII or binary output

use mod_forest, only: Morton_start, Morton_stop, sfc_to_igrid, igrid_to_node
use mod_global_parameters
use mod_usr_methods, only: usr_aux_output
use mod_physics
integer, intent(in) :: qunit

integer             :: Morton_no,igrid,ix1,ix2,ig1,ig2,level
integer, pointer    :: ig_to_igrid(:,:,:)
logical             :: fileopen,writeblk(max_blocks)
character(len=80)   :: filename
integer             :: filenr,ncells,ncells1,ncells2,ncellg,ncellx1,ncellx2,&
   jg1,jg2,jig1,jig2

character(len=name_len) :: wnamei(1:nw+nwauxio),xandwnamei(1:ndim+nw+nwauxio)
character(len=1024) :: outfilehead

double precision :: wval1,xval1
double precision, dimension(1:1,1:1,1:nw+nwauxio)   :: wval
double precision, dimension(1:1,1:1,1:ndim)         :: xval
double precision:: normconv(0:nw+nwauxio)

integer           :: iw,iiw,writenw,iwrite(1:nw+nwauxio),iigrid,idim
logical :: patchw(ixGlo1:ixGhi1,ixGlo2:ixGhi2)
!-----------------------------------------------------------------------------

if(level_io<1)then
 call mpistop('please specify level_io>0 for usage with oneblock')
end if

if(npe>1)then
 if(mype==0) PRINT *,'ONEBLOCK as yet to be parallelized'
 call mpistop('npe>1, oneblock')
end if

! only variables selected by w_write will be written out
normconv(0:nw+nwauxio)=one
normconv(0) = length_convert_factor
normconv(1:nw) = w_convert_factor
writenw=count(w_write(1:nw))+nwauxio
iiw=0
do iw =1,nw
 if (.not.w_write(iw))cycle
 iiw=iiw+1
 iwrite(iiw)=iw
end do
if(nwauxio>0)then
  do iw =nw+1,nw+nwauxio
   iiw=iiw+1
   iwrite(iiw)=iw
  end do
endif

allocate(ig_to_igrid(ng1(level_io),ng2(level_io),0:npe-1))
ig_to_igrid=-1
writeblk=.false.
do Morton_no=Morton_start(mype),Morton_stop(mype)
  igrid=sfc_to_igrid(Morton_no)
  level=node(plevel_,igrid)
  ig1=igrid_to_node(igrid,mype)%node%ig1
  ig2=igrid_to_node(igrid,mype)%node%ig2;
  ig_to_igrid(ig1,ig2,mype)=igrid
  if((rnode(rpxmin1_,igrid)>=xprobmin1+(xprobmax1-xprobmin1)*writespshift(1,&
     1).and.rnode(rpxmin2_,igrid)>=xprobmin2+&
     (xprobmax2-xprobmin2)*writespshift(2,1)).and.(rnode(rpxmax1_,&
     igrid)<=xprobmax1-(xprobmax1-xprobmin1)*writespshift(1,&
     2).and.rnode(rpxmax2_,igrid)<=xprobmax2-(xprobmax2-&
     xprobmin2)*writespshift(2,2))) then
    writeblk(igrid)=.true.
  end if
end do

call getheadernames(wnamei,xandwnamei,outfilehead)
ncells=0
ncells1=0;ncells2=0;
ncellg=((ixMhi1-ixMlo1+1)*(ixMhi2-ixMlo2+1))
ncellx1=ixMhi1-ixMlo1+1
ncellx2=ixMhi2-ixMlo2+1
do ig1=1,ng1(level_io)
do ig2=1,ng2(level_io)
  igrid=ig_to_igrid(ig1,ig2,mype)
  if(writeblk(igrid)) go to 20
end do
end do
20 continue
jg1=ig1;jg2=ig2;

jig1=jg1;jig2=jg2;
do ig1=1,ng1(level_io)
  jig1=ig1
  igrid=ig_to_igrid(jig1,jig2,mype)
  if(writeblk(igrid)) ncells1=ncells1+ncellx1
end do


jig1=jg1;jig2=jg2;
do ig2=1,ng2(level_io)
  jig2=ig2
  igrid=ig_to_igrid(jig1,jig2,mype)
  if(writeblk(igrid)) ncells2=ncells2+ncellx2
end do


do iigrid=1,igridstail; igrid=igrids(iigrid)
   if(.not.writeblk(igrid)) cycle
   ncells=ncells+ncellg
   pw(igrid)%wio(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nw)=pw(igrid)%w(ixGlo1:ixGhi1,&
      ixGlo2:ixGhi2,1:nw)

   if (nwauxio > 0) then
      if (.not. associated(usr_aux_output)) then
         call mpistop("usr_aux_output not defined")
      else
         call usr_aux_output(ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixMlo1-1,ixMlo2-1,&
            ixMhi1+1,ixMhi2+1, pw(igrid)%wio,pw(igrid)%x,normconv)
      end if
   end if
end do

if (saveprim) then
  do iigrid=1,igridstail; igrid=igrids(iigrid)
    if (.not.writeblk(igrid)) cycle
    call phys_to_primitive(ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixGlo1+1,ixGlo2+1,&
       ixGhi1-1,ixGhi2-1,pw(igrid)%wio,pw(igrid)%x)
    if (allocated(pw(igrid)%B0)) then
      ! add background magnetic field B0 to B
      pw(igrid)%wio(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
         iw_mag(:))=pw(igrid)%wio(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
         iw_mag(:))+pw(igrid)%B0(ixGlo1:ixGhi1,ixGlo2:ixGhi2,:,0)
    end if
  end do
else
  do iigrid=1,igridstail; igrid=igrids(iigrid)
    if (.not.writeblk(igrid)) cycle
    if (allocated(pw(igrid)%B0)) then
      ! add background magnetic field B0 to B
      if(phys_energy) pw(igrid)%wio(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
         iw_e)=pw(igrid)%wio(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
         iw_e)+0.5d0*sum(pw(igrid)%B0(ixGlo1:ixGhi1,ixGlo2:ixGhi2,:,0)**2,&
         dim=ndim+1) + sum(pw(igrid)%wio(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
         iw_mag(:))*pw(igrid)%B0(ixGlo1:ixGhi1,ixGlo2:ixGhi2,:,0),dim=ndim+1)
      pw(igrid)%wio(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
         iw_mag(:))=pw(igrid)%wio(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
         iw_mag(:))+pw(igrid)%B0(ixGlo1:ixGhi1,ixGlo2:ixGhi2,:,0)
    end if
  end do
end if

Master_cpu_open : if (mype == 0) then
 inquire(qunit,opened=fileopen)
 if (.not.fileopen) then
   ! generate filename
    filenr=snapshotini
    if (autoconvert) filenr=snapshotnext
   write(filename,'(a,i4.4,a)') TRIM(base_filename),filenr,".blk"
   select case(convert_type)
    case("oneblock")
     open(qunit,file=filename,status='unknown')
     write(qunit,*) TRIM(outfilehead)
     write(qunit,*) ncells,ncells1,ncells2
     write(qunit,*) global_time*time_convert_factor
    case("oneblockB")
     open(qunit,file=filename,form='unformatted',status='unknown')
     write(qunit) outfilehead
     write(qunit) ncells,ncells1,ncells2
     write(qunit) global_time*time_convert_factor
   end select
 end if
end if Master_cpu_open



   
   do ig2=1,ng2(level_io)
     do ix2=ixMlo2,ixMhi2

       do ig1=1,ng1(level_io)
         igrid=ig_to_igrid(ig1,ig2,mype)
         if(.not.writeblk(igrid)) cycle
         do ix1=ixMlo1,ixMhi1
           Master_write : if(mype==0) then
             select case(convert_type)
               case("oneblock")
                 write(qunit,fmt="(100(e14.6))") pw(igrid)%x(ix1,ix2,&
                    1:ndim)*normconv(0),(pw(igrid)%wio(ix1,ix2,&
                    iwrite(iw))*normconv(iwrite(iw)),iw=1,writenw)
               case("oneblockB")
                 write(qunit) real(pw(igrid)%x(ix1,ix2,1:ndim)*normconv(0)),&
                    (real(pw(igrid)%wio(ix1,ix2,&
                    iwrite(iw))*normconv(iwrite(iw))),iw=1,writenw)
             end select
           end if Master_write
         end do
       end do
    
     end do
   end do
 

close(qunit)

end subroutine oneblock
!=============================================================================
subroutine onegrid(qunit)

! this is for turning an AMR run into a single grid
! this version should work for any dimension, can be in parallel
! in 1D, should behave much like oneblock, except for header info

! only writes all 1:nw variables, no nwauxio
! may use saveprim to switch to primitives
! this version can work on multiple CPUs
! does not renormalize variables
! ASCII output

use mod_forest, only: Morton_start, Morton_stop, sfc_to_igrid
use mod_global_parameters
use mod_physics

integer, intent(in) :: qunit

integer             :: Morton_no,igrid,ix1,ix2,iw
logical             :: fileopen
character(len=80)   :: filename
integer             :: filenr

character(len=name_len) :: wnamei(1:nw+nwauxio),xandwnamei(1:ndim+nw+nwauxio)
character(len=1024)     :: outfilehead

!.. MPI variables ..
integer           :: igrid_recv,ipe
double precision  :: w_recv(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nw),&
   x_recv(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:ndim)
integer, allocatable :: intstatus(:,:)

!-----------------------------------------------------------------------------

if(nwauxio>0)then
 if(mype==0) PRINT *,'ONEGRID to be used without nwauxio'
 call mpistop('nwauxio>0, onegrid')
end if

if(saveprim)then
 if(mype==0.and.nwaux>0) PRINT *&
    ,'warning: ONEGRID used with saveprim, check auxiliaries'
end if



Master_cpu_open : if (mype == 0) then
 call getheadernames(wnamei,xandwnamei,outfilehead)
 write(outfilehead,'(a)') "#"//" "//TRIM(outfilehead)
 inquire(qunit,opened=fileopen)
 if (.not.fileopen) then
   ! generate filename
    filenr=snapshotini
    if (autoconvert) filenr=snapshotnext
   write(filename,'(a,i4.4,a)') TRIM(base_filename),filenr,".blk"
   open(qunit,file=filename,status='unknown')
 end if
 write(qunit,"(a)")outfilehead
 write(qunit,"(i7)") ( (ixMhi1-ixMlo1+1)*(ixMhi2-ixMlo2+1) &
    )*(Morton_stop(npe-1)-Morton_start(0)+1)
end if Master_cpu_open

do Morton_no=Morton_start(mype),Morton_stop(mype)
  igrid=sfc_to_igrid(Morton_no)
  if(saveprim) call phys_to_primitive(ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixMlo1,&
     ixMlo2,ixMhi1,ixMhi2,pw(igrid)%w,pw(igrid)%x)
  if (mype/=0)then
      itag=Morton_no
      call MPI_SEND(igrid,1,MPI_INTEGER, 0,itag,icomm,ierrmpi)
      call MPI_SEND(pw(igrid)%x,1,type_block_xcc_io, 0,itag,icomm,ierrmpi)
      itag=igrid
      call MPI_SEND(pw(igrid)%w,1,type_block_io, 0,itag,icomm,ierrmpi)
  else
   do ix2=ixMlo2,ixMhi2
   do ix1=ixMlo1,ixMhi1
      do iw=1,nw
        if( dabs(pw(igrid)%w(ix1,ix2,iw)) < 1.0d-32 ) pw(igrid)%w(ix1,ix2,&
           iw) = zero
      enddo
       write(qunit,fmt="(100(e14.6))") pw(igrid)%x(ix1,ix2,1:ndim),&
          pw(igrid)%w(ix1,ix2,1:nw)
   end do
   end do
  end if
end do

if(mype==0.and.npe>1) allocate(intstatus(MPI_STATUS_SIZE,1))

Manycpu : if (npe>1) then
 if (mype==0) then
  loop_cpu : do ipe =1, npe-1
   loop_Morton : do Morton_no=Morton_start(ipe),Morton_stop(ipe)
         itag=Morton_no
         call MPI_RECV(igrid_recv,1,MPI_INTEGER, ipe,itag,icomm,intstatus(:,1),&
            ierrmpi)
         call MPI_RECV(x_recv,1,type_block_xcc_io, ipe,itag,icomm,intstatus(:,&
            1),ierrmpi)
         itag=igrid_recv
         call MPI_RECV(w_recv,1,type_block_io, ipe,itag,icomm,intstatus(:,1),&
            ierrmpi)
         do ix2=ixMlo2,ixMhi2
         do ix1=ixMlo1,ixMhi1
            do iw=1,nw
              if( dabs(pw(igrid)%w(ix1,ix2,&
                 iw)) < smalldouble ) pw(igrid)%w(ix1,ix2,iw) = zero
            enddo
            write(qunit,fmt="(100(e14.6))") x_recv(ix1,ix2,1:ndim),w_recv(ix1,&
               ix2,1:nw)
         end do
         end do
   end do loop_Morton
  end do loop_cpu
 end if
end if Manycpu

if (npe>1) then
  call MPI_BARRIER(icomm,ierrmpi)
  if(mype==0)deallocate(intstatus)
endif

if(mype==0) close(qunit)
end subroutine onegrid
!============================================================================
subroutine tecplot(qunit)

! output for tecplot (ASCII format)
! not parallel, uses calc_grid to compute nwauxio variables
! allows renormalizing using convert factors

use mod_global_parameters

integer, intent(in) :: qunit

integer::               igrid,iigrid,level,igonlevel,iw,idim,ix1,ix2
integer::               NumGridsOnLevel(1:nlevelshi)
integer :: nx1,nx2,nxC1,nxC2,nodesonlevel,elemsonlevel,ixCmin1,ixCmin2,ixCmax1,&
   ixCmax2,ixCCmin1,ixCCmin2,ixCCmax1,ixCCmax2

integer ::              nodes, elems


double precision :: x_TEC(ndim), w_TEC(nw+nwauxio)

double precision, dimension(ixMlo1-1:ixMhi1,ixMlo2-1:ixMhi2,ndim) :: xC_TMP
double precision, dimension(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ndim)   :: xCC_TMP
double precision, dimension(ixMlo1-1:ixMhi1,ixMlo2-1:ixMhi2,ndim) :: xC
double precision, dimension(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ndim)   :: xCC

double precision, dimension(ixMlo1-1:ixMhi1,ixMlo2-1:ixMhi2,&
   nw+nwauxio)   :: wC_TMP
double precision, dimension(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
   nw+nwauxio)     :: wCC_TMP
double precision, dimension(0:nw+nwauxio)                   :: normconv
logical :: fileopen,first
character(len=80) :: filename
integer  :: filenr

!!! possible length conflict
character(len=1024) :: tecplothead

character(len=name_len) :: wnamei(1:nw+nwauxio),xandwnamei(1:ndim+nw+nwauxio)
character(len=1024) :: outfilehead
!-----------------------------------------------------------------------------

if(npe>1)then
 if(mype==0) PRINT *,'tecplot not parallel, use tecplotmpi'
 call mpistop('npe>1, tecplot')
end if

if(nw/=count(w_write(1:nw)))then
 if(mype==0) PRINT *,'tecplot does not use w_write=F'
 call mpistop('w_write, tecplot')
end if

if(nocartesian)then
 if(mype==0) PRINT *,'tecplot with nocartesian and typeaxial=',typeaxial
endif

inquire(qunit,opened=fileopen)
if (.not.fileopen) then
   ! generate filename
   filenr=snapshotini
   if (autoconvert) filenr=snapshotnext
   write(filename,'(a,i4.4,a)') TRIM(base_filename),filenr,".plt"
   open(qunit,file=filename,status='unknown')
end if

call getheadernames(wnamei,xandwnamei,outfilehead)

write(tecplothead,'(a)') "VARIABLES = "//TRIM(outfilehead)

write(qunit,'(a)') tecplothead(1:len_trim(tecplothead))

NumGridsOnLevel(1:nlevelshi)=0
do level=levmin,levmax
   NumGridsOnLevel(level)=0
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      if (node(plevel_,igrid)/=level) cycle
      NumGridsOnLevel(level)=NumGridsOnLevel(level)+1
   end do
end do

nx1=ixMhi1-ixMlo1+1;nx2=ixMhi2-ixMlo2+1;
nxC1=nx1+1;nxC2=nx2+1;


do level=levmin,levmax
   nodesonlevel=NumGridsOnLevel(level)*nxC1*nxC2
   elemsonlevel=NumGridsOnLevel(level)*nx1*nx2
   ! for all tecplot variants coded up here, we let the TECPLOT ZONES coincide
   ! with the AMR grid LEVEL. Other options would be
   !    let each grid define a zone: inefficient for TECPLOT internal workings
   !       hence not implemented
   !    let entire octree define 1 zone: no difference in interpolation
   !       properties across TECPLOT zones detected as yet, hence not done
   select case(convert_type)
     case('tecplot')
       ! in this option, we store the corner coordinates, as well as the corner
       ! values of all variables (obtained by averaging). This allows POINT packaging,
       ! and thus we can save full grid info by using one call to calc_grid
       write(qunit,"(a,i7,a,a,i7,a,i7,a,f25.16,a,a)") 'ZONE T="',level,'"',&
          ', N=',nodesonlevel,', E=',elemsonlevel, ', SOLUTIONTIME=',&
          global_time*time_convert_factor,', DATAPACKING=POINT, ZONETYPE=',&
            'FEQUADRILATERAL'
       do iigrid=1,igridstail; igrid=igrids(iigrid);
         if (node(plevel_,igrid)/=level) cycle
         call calc_x(igrid,xC,xCC)
         call calc_grid(qunit,igrid,xC,xCC,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,&
            normconv,ixCmin1,ixCmin2,ixCmax1,ixCmax2,ixCCmin1,ixCCmin2,&
            ixCCmax1,ixCCmax2,.true.)
         do ix2=ixCmin2,ixCmax2
         do ix1=ixCmin1,ixCmax1
            x_TEC(1:ndim)=xC_TMP(ix1,ix2,1:ndim)*normconv(0)
            w_TEC(1:nw+nwauxio)=wC_TMP(ix1,ix2,&
               1:nw+nwauxio)*normconv(1:nw+nwauxio)
            write(qunit,fmt="(100(e14.6))") x_TEC, w_TEC
         end do
         end do
       enddo
     case('tecplotCC')
       ! in this option, we store the corner coordinates, and the cell center
       ! values of all variables. Due to this mix of corner/cell center, we must
       ! use BLOCK packaging, and thus we have enormous overhead by using
       ! calc_grid repeatedly to merely fill values of cell corner coordinates
       ! and cell center values per dimension, per variable
       if(ndim+nw+nwauxio>99) call mpistop(&
          "adjust format specification in writeout")
       if(nw+nwauxio==1)then
         ! to make tecplot happy: avoid [ndim+1-ndim+1] in varlocation varset
         ! and just set [ndim+1]
         write(qunit,"(a,i7,a,a,i7,a,i7,a,f25.16,a,i1,a,a)") 'ZONE T="',level,&
            '"',', N=',nodesonlevel,', E=',elemsonlevel, ', SOLUTIONTIME=',&
            global_time*time_convert_factor,&
            ', DATAPACKING=BLOCK, VARLOCATION=([', ndim+1,&
            ']=CELLCENTERED), ZONETYPE=',  'FEQUADRILATERAL'
       else
        if(ndim+nw+nwauxio<10) then
         ! difference only in length of integer format specification for ndim+nw+nwauxio
         write(qunit,"(a,i7,a,a,i7,a,i7,a,f25.16,a,i1,a,i1,a,a)") 'ZONE T="',&
            level,'"',', N=',nodesonlevel,', E=',elemsonlevel,&
             ', SOLUTIONTIME=',global_time*time_convert_factor,&
            ', DATAPACKING=BLOCK, VARLOCATION=([', ndim+1,'-',ndim+nw+nwauxio,&
            ']=CELLCENTERED), ZONETYPE=',  'FEQUADRILATERAL'
        else
         write(qunit,"(a,i7,a,a,i7,a,i7,a,f25.16,a,i1,a,i2,a,a)") 'ZONE T="',&
            level,'"',', N=',nodesonlevel,', E=',elemsonlevel,&
             ', SOLUTIONTIME=',global_time*time_convert_factor,&
            ', DATAPACKING=BLOCK, VARLOCATION=([', ndim+1,'-',ndim+nw+nwauxio,&
            ']=CELLCENTERED), ZONETYPE=',  'FEQUADRILATERAL'
        endif
       endif
       do idim=1,ndim
         first=(idim==1)
         do iigrid=1,igridstail; igrid=igrids(iigrid);
            if (node(plevel_,igrid)/=level) cycle
            call calc_x(igrid,xC,xCC)
            call calc_grid(qunit,igrid,xC,xCC,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,&
               normconv,ixCmin1,ixCmin2,ixCmax1,ixCmax2,ixCCmin1,ixCCmin2,&
               ixCCmax1,ixCCmax2,first)
            write(qunit,fmt="(100(e14.6))") xC_TMP(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2,idim)*normconv(0)
         enddo
       enddo
       do iw=1,nw+nwauxio
         do iigrid=1,igridstail; igrid=igrids(iigrid);
            if (node(plevel_,igrid)/=level) cycle
            call calc_x(igrid,xC,xCC)
            call calc_grid(qunit,igrid,xC,xCC,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,&
               normconv,ixCmin1,ixCmin2,ixCmax1,ixCmax2,ixCCmin1,ixCCmin2,&
               ixCCmax1,ixCCmax2,.true.)
            write(qunit,fmt="(100(e14.6))") wCC_TMP(ixCCmin1:ixCCmax1,&
               ixCCmin2:ixCCmax2,iw)*normconv(iw)
         enddo
       enddo
     case default
       call mpistop('no such tecplot type')
   end select
   igonlevel=0
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      if (node(plevel_,igrid)/=level) cycle
      igonlevel=igonlevel+1
      call save_conntec(qunit,igrid,igonlevel)
   end do
end do


close(qunit)

end subroutine tecplot
!=============================================================================
subroutine save_conntec(qunit,igrid,igonlevel)

! this saves the basic line, quad and brick connectivity,
! as used by TECPLOT file outputs for unstructured grid

use mod_global_parameters

integer, intent(in) :: qunit, igrid, igonlevel

integer :: nx1,nx2, nxC1,nxC2, ix1,ix2

   integer, external:: nodenumbertec2D 

!-----------------------------------------------------------------------------
nx1=ixMhi1-ixMlo1+1;nx2=ixMhi2-ixMlo2+1;
nxC1=nx1+1;nxC2=nx2+1;

! connectivity list
do ix2=1,nx2
do ix1=1,nx1
   
   
   ! basic quadrilateral connectivity
   write(qunit,'(4(i7,1x))') nodenumbertec2D(ix1,  ix2-1,nxC1,nxC2,igonlevel,&
      igrid),nodenumbertec2D(ix1+1,ix2-1,nxC1,nxC2,igonlevel,igrid),&
      nodenumbertec2D(ix1+1,ix2  ,nxC1,nxC2,igonlevel,igrid),&
      nodenumbertec2D(ix1  ,ix2  ,nxC1,nxC2,igonlevel,igrid)
  
   
end do
end do

end subroutine save_conntec
!=============================================================================
integer function nodenumbertec1D(i1,nx1,ig,igrid)

integer, intent(in):: i1,nx1,ig,igrid
!-----------------------------------------------------------------------------
nodenumbertec1D=i1+(ig-1)*nx1

if(nodenumbertec1D>9999999)call mpistop("too large nodenumber")
end function nodenumbertec1D
!====================================================================================
integer function nodenumbertec2D(i1,i2,nx1,nx2,ig,igrid)

integer, intent(in):: i1,i2,nx1,nx2,ig,igrid
!-----------------------------------------------------------------------------
nodenumbertec2D=i1+i2*nx1+(ig-1)*nx1*nx2

if(nodenumbertec2D>9999999)call mpistop("too large nodenumber")
end function nodenumbertec2D
!====================================================================================
integer function nodenumbertec3D(i1,i2,i3,nx1,nx2,nx3,ig,igrid)

integer, intent(in):: i1,i2,i3,nx1,nx2,nx3,ig,igrid
!-----------------------------------------------------------------------------
nodenumbertec3D=i1+i2*nx1+i3*nx1*nx2+(ig-1)*nx1*nx2*nx3

if(nodenumbertec3D>9999999)call mpistop("too large nodenumber")
end function nodenumbertec3D
!=============================================================================
subroutine unstructuredvtk(qunit)

! output for vtu format to paraview
! not parallel, uses calc_grid to compute nwauxio variables
! allows renormalizing using convert factors

use mod_global_parameters

integer, intent(in) ::    qunit

double precision ::  x_VTK(1:3)

double precision, dimension(ixMlo1-1:ixMhi1,ixMlo2-1:ixMhi2,ndim) :: xC_TMP
double precision, dimension(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ndim)   :: xCC_TMP
double precision, dimension(ixMlo1-1:ixMhi1,ixMlo2-1:ixMhi2,ndim) :: xC
double precision, dimension(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ndim)   :: xCC

double precision, dimension(ixMlo1-1:ixMhi1,ixMlo2-1:ixMhi2,&
   nw+nwauxio)   :: wC_TMP
double precision, dimension(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
   nw+nwauxio)     :: wCC_TMP
double precision, dimension(0:nw+nwauxio)                   :: normconv
integer::               igrid,iigrid,level,igonlevel,icel,ixCmin1,ixCmin2,&
   ixCmax1,ixCmax2,ixCCmin1,ixCCmin2,ixCCmax1,ixCCmax2,iw
integer::               NumGridsOnLevel(1:nlevelshi)
integer :: nx1,nx2,nxC1,nxC2,nodesonlevel,elemsonlevel,nc,np,VTK_type,ix1,ix2

character(len=80)::  filename
integer          :: filenr

character(len=name_len) :: wnamei(1:nw+nwauxio),xandwnamei(1:ndim+nw+nwauxio)
character(len=1024) :: outfilehead

logical :: fileopen
!-----------------------------------------------------------------------------

if(npe>1)then
 if(mype==0) PRINT *,'unstructuredvtk not parallel, use vtumpi'
 call mpistop('npe>1, unstructuredvtk')
end if

inquire(qunit,opened=fileopen)
if(.not.fileopen)then
  ! generate filename
   filenr=snapshotini
   if (autoconvert) filenr=snapshotnext
  write(filename,'(a,i4.4,a)') TRIM(base_filename),filenr,".vtu"
  ! Open the file for the header part
  open(qunit,file=filename,status='unknown')
endif

call getheadernames(wnamei,xandwnamei,outfilehead)

! generate xml header
write(qunit,'(a)')'<?xml version="1.0"?>'
write(qunit,'(a)',advance='no') '<VTKFile type="UnstructuredGrid"'
write(qunit,'(a)')' version="0.1" byte_order="LittleEndian">'
write(qunit,'(a)')'<UnstructuredGrid>'
write(qunit,'(a)')'<FieldData>'
write(qunit,'(2a)')'<DataArray type="Float32" Name="TIME" ',&
   'NumberOfTuples="1" format="ascii">'
write(qunit,*) dble(dble(global_time)*time_convert_factor)
write(qunit,'(a)')'</DataArray>'
write(qunit,'(a)')'</FieldData>'

! number of cells, number of corner points, per grid.
nx1=ixMhi1-ixMlo1+1;nx2=ixMhi2-ixMlo2+1;
nxC1=nx1+1;nxC2=nx2+1;
nc=nx1*nx2
np=nxC1*nxC2

! Note: using the w_write, writelevel, writespshift
! we can clip parts of the grid away, select variables, levels etc.
do level=levmin,levmax
 if (writelevel(level)) then
   do iigrid=1,igridstail; igrid=igrids(iigrid);
    if (node(plevel_,igrid)/=level) cycle
    ! only output a grid when fully within clipped region selected
    ! by writespshift array
    if ((rnode(rpxmin1_,igrid)>=xprobmin1+(xprobmax1-xprobmin1)*writespshift(1,&
       1).and.rnode(rpxmin2_,igrid)>=xprobmin2+&
       (xprobmax2-xprobmin2)*writespshift(2,1)).and.(rnode(rpxmax1_,&
       igrid)<=xprobmax1-(xprobmax1-xprobmin1)*writespshift(1,&
       2).and.rnode(rpxmax2_,igrid)<=xprobmax2-(xprobmax2-&
       xprobmin2)*writespshift(2,2))) then
      call calc_x(igrid,xC,xCC)
      call calc_grid(qunit,igrid,xC,xCC,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,normconv,&
         ixCmin1,ixCmin2,ixCmax1,ixCmax2,ixCCmin1,ixCCmin2,ixCCmax1,ixCCmax2,&
         .true.)
      select case(convert_type)
       case('vtu')
         ! we write out every grid as one VTK PIECE
         write(qunit,'(a,i7,a,i7,a)') '<Piece NumberOfPoints="',np,&
            '" NumberOfCells="',nc,'">'
         write(qunit,'(a)')'<PointData>'
         do iw=1,nw+nwauxio
         if(iw<=nw) then
            if(.not.w_write(iw)) cycle
         endif

            write(qunit,'(a,a,a)')'<DataArray type="Float64" Name="',&
               TRIM(wnamei(iw)),'" format="ascii">'
            write(qunit,'(200(1pe14.6))') ((wC_TMP(ix1,ix2,iw)*normconv(iw),&
               ix1=ixCmin1,ixCmax1),ix2=ixCmin2,ixCmax2)
            write(qunit,'(a)')'</DataArray>'
         enddo
         write(qunit,'(a)')'</PointData>'

         write(qunit,'(a)')'<Points>'
         write(qunit,'(a)'&
)'<DataArray type="Float32" NumberOfComponents="3" format="ascii">'
         ! write cell corner coordinates in a backward dimensional loop, always 3D output
         do ix2=ixCmin2,ixCmax2 
         do ix1=ixCmin1,ixCmax1 
            x_VTK(1:3)=zero;
            x_VTK(1:ndim)=xC_TMP(ix1,ix2,1:ndim)*normconv(0);
            write(qunit,'(3(1pe14.6))') x_VTK
         end do 
         end do 
         write(qunit,'(a)')'</DataArray>'
         write(qunit,'(a)')'</Points>'
       case('vtuCC')
         ! we write out every grid as one VTK PIECE
         write(qunit,'(a,i7,a,i7,a)') '<Piece NumberOfPoints="',np,&
            '" NumberOfCells="',nc,'">'
         write(qunit,'(a)')'<CellData>'
         do iw=1,nw+nwauxio
         if(iw<=nw) then
            if(.not.w_write(iw)) cycle
         endif

            write(qunit,'(a,a,a)')'<DataArray type="Float64" Name="',&
               TRIM(wnamei(iw)),'" format="ascii">'
            write(qunit,'(200(1pe14.6))') ((wCC_TMP(ix1,ix2,iw)*normconv(iw),&
               ix1=ixCCmin1,ixCCmax1),ix2=ixCCmin2,ixCCmax2)
            write(qunit,'(a)')'</DataArray>'
         enddo
         write(qunit,'(a)')'</CellData>'

         write(qunit,'(a)')'<Points>'
         write(qunit,'(a)'&
)'<DataArray type="Float32" NumberOfComponents="3" format="ascii">'
         ! write cell corner coordinates in a backward dimensional loop, always 3D output
         do ix2=ixCmin2,ixCmax2 
         do ix1=ixCmin1,ixCmax1 
            x_VTK(1:3)=zero;
            x_VTK(1:ndim)=xC_TMP(ix1,ix2,1:ndim)*normconv(0);
            write(qunit,'(3(1pe14.6))') x_VTK
         end do 
         end do 
         write(qunit,'(a)')'</DataArray>'
         write(qunit,'(a)')'</Points>'
      end select


      write(qunit,'(a)')'<Cells>'

      ! connectivity part
      write(qunit,'(a)'&
         )'<DataArray type="Int32" Name="connectivity" format="ascii">'
      call save_connvtk(qunit,igrid)
      write(qunit,'(a)')'</DataArray>'

      ! offsets data array
      write(qunit,'(a)'&
         )'<DataArray type="Int32" Name="offsets" format="ascii">'
      do icel=1,nc
         write(qunit,'(i7)') icel*(2**2)
      end do
      write(qunit,'(a)')'</DataArray>'

      ! VTK cell type data array
      write(qunit,'(a)')'<DataArray type="Int32" Name="types" format="ascii">'
      ! VTK_LINE=3; VTK_PIXEL=8; VTK_VOXEL=11 -> vtk-syntax
      
       VTK_type=8 
      
      do icel=1,nc
         write(qunit,'(i2)') VTK_type
      enddo
      write(qunit,'(a)')'</DataArray>'

      write(qunit,'(a)')'</Cells>'

      write(qunit,'(a)')'</Piece>'
    endif
   enddo
 endif
enddo

write(qunit,'(a)')'</UnstructuredGrid>'
write(qunit,'(a)')'</VTKFile>'
close(qunit)

end subroutine unstructuredvtk
!====================================================================================
subroutine unstructuredvtkB(qunit)

! output for vtu format to paraview, binary version output
! not parallel, uses calc_grid to compute nwauxio variables
! allows renormalizing using convert factors
use mod_forest, only: Morton_start, Morton_stop, sfc_to_igrid
use mod_global_parameters

integer, intent(in) ::    qunit

double precision ::  x_VTK(1:3)

double precision, dimension(ixMlo1-1:ixMhi1,ixMlo2-1:ixMhi2,ndim) :: xC_TMP
double precision, dimension(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ndim)   :: xCC_TMP
double precision, dimension(ixMlo1-1:ixMhi1,ixMlo2-1:ixMhi2,ndim) :: xC
double precision, dimension(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ndim)   :: xCC

double precision, dimension(ixMlo1-1:ixMhi1,ixMlo2-1:ixMhi2,&
   nw+nwauxio):: wC_TMP
double precision, dimension(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
   nw+nwauxio)  :: wCC_TMP
double precision :: normconv(0:nw+nwauxio)

integer, allocatable :: intstatus(:,:)
integer :: ipe,igrid,level,icel,ixCmin1,ixCmin2,ixCmax1,ixCmax2,ixCCmin1,&
   ixCCmin2,ixCCmax1,ixCCmax2,Morton_no,Morton_length
integer :: nx1,nx2,nxC1,nxC2,nc,np,VTK_type,ix1,ix2,filenr
integer*8 :: offset

integer::  k,iw
integer::  length,lengthcc,length_coords,length_conn,length_offsets
character::  buf
character(len=80)::  filename
character(len=name_len) :: wnamei(1:nw+nwauxio),xandwnamei(1:ndim+nw+nwauxio)
character(len=1024) :: outfilehead

logical ::   fileopen,cell_corner=.true.
logical, allocatable :: Morton_aim(:),Morton_aim_p(:)
!-----------------------------------------------------------------------------

normconv=one
Morton_length=Morton_stop(npe-1)-Morton_start(0)+1
allocate(Morton_aim(Morton_start(0):Morton_stop(npe-1)))
allocate(Morton_aim_p(Morton_start(0):Morton_stop(npe-1)))
Morton_aim=.false.
Morton_aim_p=.false.
do Morton_no=Morton_start(mype),Morton_stop(mype)
  igrid=sfc_to_igrid(Morton_no)
  level=node(plevel_,igrid)
  ! we can clip parts of the grid away, select variables, levels etc.
  if(writelevel(level)) then
   ! only output a grid when fully within clipped region selected
   ! by writespshift array
   if((rnode(rpxmin1_,igrid)>=xprobmin1+(xprobmax1-xprobmin1)*writespshift(1,&
      1).and.rnode(rpxmin2_,igrid)>=xprobmin2+&
      (xprobmax2-xprobmin2)*writespshift(2,1)).and.(rnode(rpxmax1_,&
      igrid)<=xprobmax1-(xprobmax1-xprobmin1)*writespshift(1,&
      2).and.rnode(rpxmax2_,igrid)<=xprobmax2-(xprobmax2-&
      xprobmin2)*writespshift(2,2))) then
     Morton_aim_p(Morton_no)=.true.
   end if
  end if
end do
call MPI_ALLREDUCE(Morton_aim_p,Morton_aim,Morton_length,MPI_LOGICAL,MPI_LOR,&
   icomm,ierrmpi)
select case(convert_type)
 case('vtuB','vtuBmpi')
   cell_corner=.true.
 case('vtuBCC','vtuBCCmpi')
   cell_corner=.false.
end select
if (mype /= 0) then
 do Morton_no=Morton_start(mype),Morton_stop(mype)
   if(.not. Morton_aim(Morton_no)) cycle
   igrid=sfc_to_igrid(Morton_no)
   call calc_x(igrid,xC,xCC)
   call calc_grid(qunit,igrid,xC,xCC,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,normconv,&
      ixCmin1,ixCmin2,ixCmax1,ixCmax2,ixCCmin1,ixCCmin2,ixCCmax1,ixCCmax2,&
      .true.)
   itag=Morton_no
   call MPI_SEND(xC_TMP,1,type_block_xc_io, 0,itag,icomm,ierrmpi)
   if(cell_corner) then
     call MPI_SEND(wC_TMP,1,type_block_wc_io, 0,itag,icomm,ierrmpi)
   else
     call MPI_SEND(wCC_TMP,1,type_block_wcc_io, 0,itag,icomm,ierrmpi)
   endif
 end do

else
 ! mype==0
 offset=0

 inquire(qunit,opened=fileopen)
 if(.not.fileopen)then
   ! generate filename
    filenr=snapshotini
    if (autoconvert) filenr=snapshotnext
   write(filename,'(a,i4.4,a)') TRIM(base_filename),filenr,".vtu"
   ! Open the file for the header part
   open(qunit,file=filename,status='replace')
 endif

 call getheadernames(wnamei,xandwnamei,outfilehead)

 ! generate xml header
 write(qunit,'(a)')'<?xml version="1.0"?>'
 write(qunit,'(a)',advance='no') '<VTKFile type="UnstructuredGrid"'
 write(qunit,'(a)')' version="0.1" byte_order="LittleEndian">'
 write(qunit,'(a)')'<UnstructuredGrid>'
 write(qunit,'(a)')'<FieldData>'
 write(qunit,'(2a)')'<DataArray type="Float32" Name="TIME" ',&
    'NumberOfTuples="1" format="ascii">'
 write(qunit,*) real(global_time*time_convert_factor)
 write(qunit,'(a)')'</DataArray>'
 write(qunit,'(a)')'</FieldData>'

 ! number of cells, number of corner points, per grid.
 nx1=ixMhi1-ixMlo1+1;nx2=ixMhi2-ixMlo2+1;
 nxC1=nx1+1;nxC2=nx2+1;
 nc=nx1*nx2
 np=nxC1*nxC2

 length=np*size_real
 lengthcc=nc*size_real

 length_coords=3*length
 length_conn=2**2*size_int*nc
 length_offsets=nc*size_int

 ! Note: using the w_write, writelevel, writespshift
 do Morton_no=Morton_start(0),Morton_stop(0)
    if(.not. Morton_aim(Morton_no)) cycle
    if(cell_corner) then
      ! we write out every grid as one VTK PIECE
      write(qunit,'(a,i7,a,i7,a)') '<Piece NumberOfPoints="',np,&
         '" NumberOfCells="',nc,'">'
      write(qunit,'(a)')'<PointData>'
      do iw=1,nw+nwauxio
         if(iw<=nw) then
            if(.not.w_write(iw)) cycle
         endif

         write(qunit,'(a,a,a,i16,a)')'<DataArray type="Float32" Name="',&
            TRIM(wnamei(iw)), '" format="appended" offset="',offset,'">'
         write(qunit,'(a)')'</DataArray>'
         offset=offset+length+size_int
      enddo
      write(qunit,'(a)')'</PointData>'

      write(qunit,'(a)')'<Points>'
      write(qunit,'(a,i16,a)')&
'<DataArray type="Float32" NumberOfComponents="3" format="appended" offset="',&
         offset,'"/>'
      ! write cell corner coordinates in a backward dimensional loop, always 3D output
      offset=offset+length_coords+size_int
      write(qunit,'(a)')'</Points>'
    else
      ! we write out every grid as one VTK PIECE
      write(qunit,'(a,i7,a,i7,a)') '<Piece NumberOfPoints="',np,&
         '" NumberOfCells="',nc,'">'
      write(qunit,'(a)')'<CellData>'
      do iw=1,nw+nwauxio
         if(iw<=nw) then
            if(.not.w_write(iw)) cycle
         endif

         write(qunit,'(a,a,a,i16,a)')'<DataArray type="Float32" Name="',&
            TRIM(wnamei(iw)), '" format="appended" offset="',offset,'">'
         write(qunit,'(a)')'</DataArray>'
         offset=offset+lengthcc+size_int
      enddo
      write(qunit,'(a)')'</CellData>'

      write(qunit,'(a)')'<Points>'
      write(qunit,'(a,i16,a)')&
'<DataArray type="Float32" NumberOfComponents="3" format="appended" offset="',&
         offset,'"/>'
      ! write cell corner coordinates in a backward dimensional loop, always 3D output
      offset=offset+length_coords+size_int
      write(qunit,'(a)')'</Points>'
    end if

    write(qunit,'(a)')'<Cells>'

    ! connectivity part
    write(qunit,'(a,i16,a)'&
)'<DataArray type="Int32" Name="connectivity" format="appended" offset="',&
       offset,'"/>'
    offset=offset+length_conn+size_int

    ! offsets data array
    write(qunit,'(a,i16,a)')&
        '<DataArray type="Int32" Name="offsets" format="appended" offset="',&
       offset,'"/>'
    offset=offset+length_offsets+size_int

    ! VTK cell type data array
    write(qunit,'(a,i16,a)')&
        '<DataArray type="Int32" Name="types" format="appended" offset="',&
       offset,'"/>'
    offset=offset+size_int+nc*size_int

    write(qunit,'(a)')'</Cells>'

    write(qunit,'(a)')'</Piece>'
 end do
 ! write metadata communicated from other processors
 if(npe>1)then
  do ipe=1, npe-1
    do Morton_no=Morton_start(ipe),Morton_stop(ipe)
      if(.not. Morton_aim(Morton_no)) cycle
      if(cell_corner) then
        ! we write out every grid as one VTK PIECE
        write(qunit,'(a,i7,a,i7,a)') '<Piece NumberOfPoints="',np,&
           '" NumberOfCells="',nc,'">'
        write(qunit,'(a)')'<PointData>'
        do iw=1,nw+nwauxio
         if(iw<=nw) then
            if(.not.w_write(iw)) cycle
         endif

           write(qunit,'(a,a,a,i16,a)')'<DataArray type="Float32" Name="',&
              TRIM(wnamei(iw)), '" format="appended" offset="',offset,'">'
           write(qunit,'(a)')'</DataArray>'
           offset=offset+length+size_int
        enddo
        write(qunit,'(a)')'</PointData>'

        write(qunit,'(a)')'<Points>'
        write(qunit,'(a,i16,a)')&
'<DataArray type="Float32" NumberOfComponents="3" format="appended" offset="',&
           offset,'"/>'
        ! write cell corner coordinates in a backward dimensional loop, always 3D output
        offset=offset+length_coords+size_int
        write(qunit,'(a)')'</Points>'
      else
        ! we write out every grid as one VTK PIECE
        write(qunit,'(a,i7,a,i7,a)') '<Piece NumberOfPoints="',np,&
           '" NumberOfCells="',nc,'">'
        write(qunit,'(a)')'<CellData>'
        do iw=1,nw+nwauxio
         if(iw<=nw) then
            if(.not.w_write(iw)) cycle
         endif

           write(qunit,'(a,a,a,i16,a)')'<DataArray type="Float32" Name="',&
              TRIM(wnamei(iw)), '" format="appended" offset="',offset,'">'
           write(qunit,'(a)')'</DataArray>'
           offset=offset+lengthcc+size_int
        enddo
        write(qunit,'(a)')'</CellData>'

        write(qunit,'(a)')'<Points>'
        write(qunit,'(a,i16,a)')&
'<DataArray type="Float32" NumberOfComponents="3" format="appended" offset="',&
           offset,'"/>'
        ! write cell corner coordinates in a backward dimensional loop, always 3D output
        offset=offset+length_coords+size_int
        write(qunit,'(a)')'</Points>'
      end if

      write(qunit,'(a)')'<Cells>'

      ! connectivity part
      write(qunit,'(a,i16,a)'&
)'<DataArray type="Int32" Name="connectivity" format="appended" offset="',&
         offset,'"/>'
      offset=offset+length_conn+size_int

      ! offsets data array
      write(qunit,'(a,i16,a)')&
          '<DataArray type="Int32" Name="offsets" format="appended" offset="',&
         offset,'"/>'
      offset=offset+length_offsets+size_int

      ! VTK cell type data array
      write(qunit,'(a,i16,a)')&
          '<DataArray type="Int32" Name="types" format="appended" offset="',&
         offset,'"/>'
      offset=offset+size_int+nc*size_int

      write(qunit,'(a)')'</Cells>'

      write(qunit,'(a)')'</Piece>'
    end do
  end do
 end if

 write(qunit,'(a)')'</UnstructuredGrid>'
 write(qunit,'(a)')'<AppendedData encoding="raw">'
 close(qunit)
 open(qunit,file=filename,access='stream',form='unformatted',&
    position='append')
 buf='_'
 write(qunit) TRIM(buf)

 do Morton_no=Morton_start(0),Morton_stop(0)
   if(.not. Morton_aim(Morton_no)) cycle
   igrid=sfc_to_igrid(Morton_no)

   call calc_x(igrid,xC,xCC)
   call calc_grid(qunit,igrid,xC,xCC,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,normconv,&
      ixCmin1,ixCmin2,ixCmax1,ixCmax2,ixCCmin1,ixCCmin2,ixCCmax1,ixCCmax2,&
      .true.)

   do iw=1,nw+nwauxio
         if(iw<=nw) then
            if(.not.w_write(iw)) cycle
         endif
     if(cell_corner) then

       write(qunit) length
       write(qunit) ((real(wC_TMP(ix1,ix2,iw)*normconv(iw)),ix1=ixCmin1,&
          ixCmax1),ix2=ixCmin2,ixCmax2)
     else
       write(qunit) lengthcc
       write(qunit) ((real(wCC_TMP(ix1,ix2,iw)*normconv(iw)),ix1=ixCCmin1,&
          ixCCmax1),ix2=ixCCmin2,ixCCmax2)
     end if
   enddo

   write(qunit) length_coords
   do ix2=ixCmin2,ixCmax2 
   do ix1=ixCmin1,ixCmax1 
     x_VTK(1:3)=zero;
     x_VTK(1:ndim)=xC_TMP(ix1,ix2,1:ndim)*normconv(0);
     do k=1,3
      write(qunit) real(x_VTK(k))
     end do
   end do 
   end do 

   write(qunit) length_conn
   do ix2=1,nx2
   do ix1=1,nx1
   
   
   write(qunit)(ix2-1)*nxC1+ix1-1, (ix2-1)*nxC1+ix1,ix2*nxC1+ix1-1,&
      ix2*nxC1+ix1
    
   
   end do
   end do

   write(qunit) length_offsets
   do icel=1,nc
     write(qunit) icel*(2**2)
   end do


  
   VTK_type=8 
  
   write(qunit) size_int*nc
   do icel=1,nc
     write(qunit) VTK_type
   end do
 end do
 allocate(intstatus(MPI_STATUS_SIZE,1))
 if(npe>1)then
  ixCCmin1=ixMlo1;ixCCmin2=ixMlo2; ixCCmax1=ixMhi1;ixCCmax2=ixMhi2;
  ixCmin1=ixMlo1-1;ixCmin2=ixMlo2-1; ixCmax1=ixMhi1;ixCmax2=ixMhi2;
  do ipe=1, npe-1
    do Morton_no=Morton_start(ipe),Morton_stop(ipe)
      if(.not. Morton_aim(Morton_no)) cycle
      itag=Morton_no
      call MPI_RECV(xC_TMP,1,type_block_xc_io, ipe,itag,icomm,intstatus(:,1),&
         ierrmpi)
      if(cell_corner) then
        call MPI_RECV(wC_TMP,1,type_block_wc_io, ipe,itag,icomm,intstatus(:,1),&
           ierrmpi)
      else
        call MPI_RECV(wCC_TMP,1,type_block_wcc_io, ipe,itag,icomm,intstatus(:,&
           1),ierrmpi)
      end if
      do iw=1,nw+nwauxio
         if(iw<=nw) then
            if(.not.w_write(iw)) cycle
         endif
        if(cell_corner) then
          write(qunit) length
          write(qunit) ((real(wC_TMP(ix1,ix2,iw)*normconv(iw)),ix1=ixCmin1,&
             ixCmax1),ix2=ixCmin2,ixCmax2)
        else
          write(qunit) lengthcc
          write(qunit) ((real(wCC_TMP(ix1,ix2,iw)*normconv(iw)),ix1=ixCCmin1,&
             ixCCmax1),ix2=ixCCmin2,ixCCmax2)
        end if
      enddo

      write(qunit) length_coords
      do ix2=ixCmin2,ixCmax2 
      do ix1=ixCmin1,ixCmax1 
        x_VTK(1:3)=zero;
        x_VTK(1:ndim)=xC_TMP(ix1,ix2,1:ndim)*normconv(0);
        do k=1,3
         write(qunit) real(x_VTK(k))
        end do
      end do 
      end do 

      write(qunit) length_conn
      do ix2=1,nx2
      do ix1=1,nx1
      
      
      write(qunit)(ix2-1)*nxC1+ix1-1, (ix2-1)*nxC1+ix1,ix2*nxC1+ix1-1,&
         ix2*nxC1+ix1
       
      
      end do
      end do

      write(qunit) length_offsets
      do icel=1,nc
        write(qunit) icel*(2**2)
      end do
      
       VTK_type=8 
      
      write(qunit) size_int*nc
      do icel=1,nc
        write(qunit) VTK_type
      end do
    end do
  end do
 end if
 close(qunit)
 open(qunit,file=filename,status='unknown',form='formatted',position='append')
 write(qunit,'(a)')'</AppendedData>'
 write(qunit,'(a)')'</VTKFile>'
 close(qunit)
 deallocate(intstatus)
end if

deallocate(Morton_aim,Morton_aim_p)
if (npe>1) then
  call MPI_BARRIER(icomm,ierrmpi)
endif

end subroutine unstructuredvtkB
!====================================================================================
subroutine save_connvtk(qunit,igrid)

! this saves the basic line, pixel and voxel connectivity,
! as used by VTK file outputs for unstructured grid

use mod_global_parameters

integer, intent(in) :: qunit, igrid

integer :: nx1,nx2, nxC1,nxC2, ix1,ix2
!-----------------------------------------------------------------------------
nx1=ixMhi1-ixMlo1+1;nx2=ixMhi2-ixMlo2+1;
nxC1=nx1+1;nxC2=nx2+1;

do ix2=1,nx2
do ix1=1,nx1
        
        
        write(qunit,'(4(i7,1x))')(ix2-1)*nxC1+ix1-1, (ix2-1)*nxC1+ix1,&
           ix2*nxC1+ix1-1,ix2*nxC1+ix1
        
        
end do
end do

end subroutine save_connvtk
!=============================================================================
subroutine ImageDataVtk_mpi(qunit)

! output for vti format to paraview, non-binary version output
! parallel, uses calc_grid to compute nwauxio variables
! allows renormalizing using convert factors
! allows skipping of w_write selected variables

! implementation such that length of ASCII output is identical when
! run on 1 versus multiple CPUs (however, the order of the vtu pieces can differ)

use mod_forest, only: Morton_start, Morton_stop, tree_node_ptr, igrid_to_node,&
    sfc_to_igrid
use mod_global_parameters

integer, intent(in) ::    qunit

double precision, dimension(ixMlo1-1:ixMhi1,ixMlo2-1:ixMhi2,ndim) :: xC_TMP,&
   xC_TMP_recv
double precision, dimension(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ndim)   :: xCC_TMP,&
   xCC_TMP_recv
double precision, dimension(ixMlo1-1:ixMhi1,ixMlo2-1:ixMhi2,ndim) :: xC
double precision, dimension(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ndim)   :: xCC

double precision, dimension(ixMlo1-1:ixMhi1,ixMlo2-1:ixMhi2,&
   nw+nwauxio)   :: wC_TMP,wC_TMP_recv
double precision, dimension(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
   nw+nwauxio)     :: wCC_TMP,wCC_TMP_recv
double precision, dimension(0:nw+nwauxio)                   :: normconv
integer::               igrid,iigrid,level,ixCmin1,ixCmin2,ixCmax1,ixCmax2,&
   ixCCmin1,ixCCmin2,ixCCmax1,ixCCmax2
integer::               NumGridsOnLevel(1:nlevelshi)
integer :: nx1,nx2

character(len=80)::  filename
integer ::           filenr

integer, allocatable :: intstatus(:,:)
logical, allocatable :: Morton_aim(:),Morton_aim_p(:)

character(len=name_len) :: wnamei(1:nw+nwauxio),xandwnamei(1:ndim+nw+nwauxio)
character(len=1024) :: outfilehead

logical :: fileopen
integer :: ipe,Morton_no,Morton_length
integer :: ixrvCmin1,ixrvCmin2,ixrvCmax1,ixrvCmax2, ixrvCCmin1,ixrvCCmin2,&
   ixrvCCmax1,ixrvCCmax2, siz_ind, ind_send(5*2), ind_recv(5*2)
double precision    :: origin(1:3), spacing(1:3)
integer :: wholeExtent(1:6), ig1,ig2
type(tree_node_ptr) :: tree
!-----------------------------------------------------------------------------
if(levmin/=levmax) call mpistop&
   ('ImageData can only be used when levmin=levmax')

normconv(0) = length_convert_factor
normconv(1:nw) = w_convert_factor
siz_ind=5*2
Morton_length=Morton_stop(npe-1)-Morton_start(0)+1
allocate(Morton_aim(Morton_start(0):Morton_stop(npe-1)))
allocate(Morton_aim_p(Morton_start(0):Morton_stop(npe-1)))
Morton_aim=.false.
Morton_aim_p=.false.
do Morton_no=Morton_start(mype),Morton_stop(mype)
  igrid=sfc_to_igrid(Morton_no)
  level=node(plevel_,igrid)
  ! we can clip parts of the grid away, select variables, levels etc.
  if(writelevel(level)) then
   ! only output a grid when fully within clipped region selected
   ! by writespshift array
   if((rnode(rpxmin1_,igrid)>=xprobmin1+(xprobmax1-xprobmin1)*writespshift(1,&
      1).and.rnode(rpxmin2_,igrid)>=xprobmin2+&
      (xprobmax2-xprobmin2)*writespshift(2,1)).and.(rnode(rpxmax1_,&
      igrid)<=xprobmax1-(xprobmax1-xprobmin1)*writespshift(1,&
      2).and.rnode(rpxmax2_,igrid)<=xprobmax2-(xprobmax2-&
      xprobmin2)*writespshift(2,2))) then
     Morton_aim_p(Morton_no)=.true.
   end if
  end if
end do
call MPI_ALLREDUCE(Morton_aim_p,Morton_aim,Morton_length,MPI_LOGICAL,MPI_LOR,&
   icomm,ierrmpi)


if (mype /= 0) then
 do Morton_no=Morton_start(mype),Morton_stop(mype)
   if(.not. Morton_aim(Morton_no)) cycle
   igrid=sfc_to_igrid(Morton_no)
   call calc_x(igrid,xC,xCC)
   call calc_grid(qunit,igrid,xC,xCC,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,normconv,&
      ixCmin1,ixCmin2,ixCmax1,ixCmax2,ixCCmin1,ixCCmin2,ixCCmax1,ixCCmax2,&
      .true.)
   tree%node => igrid_to_node(igrid, mype)%node
    ig1 = tree%node%ig1;  ig2 = tree%node%ig2;
   itag=Morton_no
   ind_send=(/ ixCmin1,ixCmin2,ixCmax1,ixCmax2,ixCCmin1,ixCCmin2,ixCCmax1,&
      ixCCmax2, ig1,ig2 /)
   call MPI_SEND(ind_send,siz_ind,MPI_INTEGER, 0,itag,icomm,ierrmpi)
   call MPI_SEND(wC_TMP,1,type_block_wc_io, 0,itag,icomm,ierrmpi)
   call MPI_SEND(wCC_TMP,1,type_block_wcc_io, 0,itag,icomm,ierrmpi)
 end do

else

 inquire(qunit,opened=fileopen)
 if(.not.fileopen)then
    ! generate filename
    filenr=snapshotini
    if (autoconvert) filenr=snapshotnext
    write(filename,'(a,i4.4,a)') TRIM(base_filename),filenr,".vti"
   ! Open the file for the header part
   open(qunit,file=filename,status='unknown',form='formatted')
 endif

call getheadernames(wnamei,xandwnamei,outfilehead)

! number of cells per grid.
nx1=ixMhi1-ixMlo1+1;nx2=ixMhi2-ixMlo2+1;

origin      = 0
 origin(1) = xprobmin1*normconv(0);  origin(2) = xprobmin2*normconv(0);
spacing     = zero
spacing(1) = dxlevel(1)*normconv(0); spacing(2) = dxlevel(2)*normconv(0);

wholeExtent = 0
! if we use writespshift, the whole extent has to be calculated:
wholeExtent(1*2-1) = nx1 * ceiling(((xprobmax1-xprobmin1)*writespshift(1,&
   1)) /(nx1*dxlevel(1))) 
wholeExtent(2*2-1) = nx2 * ceiling(((xprobmax2-xprobmin2)*writespshift(2,&
   1)) /(nx2*dxlevel(2))) 
wholeExtent(1*2)   = nx1 * floor(((xprobmax1-xprobmin1)*(1.0d0-writespshift(1,&
   2))) /(nx1*dxlevel(1))) 
wholeExtent(2*2)   = nx2 * floor(((xprobmax2-xprobmin2)*(1.0d0-writespshift(2,&
   2))) /(nx2*dxlevel(2))) 

! generate xml header
write(qunit,'(a)')'<?xml version="1.0"?>'
write(qunit,'(a)',advance='no') '<VTKFile type="ImageData"'
write(qunit,'(a)')' version="0.1" byte_order="LittleEndian">'
write(qunit,'(a,3(1pe14.6),a,6(i10),a,3(1pe14.6),a)')'  <ImageData Origin="',&
   origin,'" WholeExtent="',wholeExtent,'" Spacing="',spacing,'">'
 write(qunit,'(a)')'<FieldData>'
 write(qunit,'(2a)')'<DataArray type="Float32" Name="TIME" ',&
    'NumberOfTuples="1" format="ascii">'
 write(qunit,*) real(global_time*time_convert_factor)
 write(qunit,'(a)')'</DataArray>'
 write(qunit,'(a)')'</FieldData>'

! write the data from proc 0
do Morton_no=Morton_start(0),Morton_stop(0)
   if(.not. Morton_aim(Morton_no)) cycle
   igrid=sfc_to_igrid(Morton_no)
   tree%node => igrid_to_node(igrid, 0)%node
    ig1 = tree%node%ig1;  ig2 = tree%node%ig2;
   call calc_x(igrid,xC,xCC)
   call calc_grid(qunit,igrid,xC,xCC,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,normconv,&
      ixCmin1,ixCmin2,ixCmax1,ixCmax2,ixCCmin1,ixCCmin2,ixCCmax1,ixCCmax2,&
      .true.)
   call write_vti(qunit,ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixCmin1,ixCmin2,ixCmax1,&
      ixCmax2,ixCCmin1,ixCCmin2,ixCCmax1,ixCCmax2,ig1,ig2,nx1,nx2,normconv,&
      wnamei,wC_TMP,wCC_TMP)
end do

if(npe>1)then
   allocate(intstatus(MPI_STATUS_SIZE,1))
   do ipe=1, npe-1
      do Morton_no=Morton_start(ipe),Morton_stop(ipe)
         if(.not. Morton_aim(Morton_no)) cycle
         itag=Morton_no
         call MPI_RECV(ind_recv,siz_ind, MPI_INTEGER, ipe,itag,icomm,&
            intstatus(:,1),ierrmpi)
         ixrvCmin1=ind_recv(1);ixrvCmin2=ind_recv(2);ixrvCmax1=ind_recv(2+1)
         ixrvCmax2=ind_recv(2+2);
         ixrvCCmin1=ind_recv(2*2+1);ixrvCCmin2=ind_recv(2*2+2)
         ixrvCCmax1=ind_recv(3*2+1);ixrvCCmax2=ind_recv(3*2+2);
         ig1=ind_recv(4*2+1);ig2=ind_recv(4*2+2);
         call MPI_RECV(wC_TMP,1,type_block_wc_io, ipe,itag,icomm,intstatus(:,&
            1),ierrmpi)
         call MPI_RECV(wCC_TMP,1,type_block_wcc_io, ipe,itag,icomm,intstatus(:,&
            1),ierrmpi)
         call write_vti(qunit,ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixrvCmin1,ixrvCmin2,&
            ixrvCmax1,ixrvCmax2,ixrvCCmin1,ixrvCCmin2,ixrvCCmax1,ixrvCCmax2,&
            ig1,ig2,nx1,nx2,normconv,wnamei,wC_TMP,wCC_TMP)
      end do
   end do
end if

write(qunit,'(a)')'</ImageData>'
write(qunit,'(a)')'</VTKFile>'
close(qunit)
if(npe>1) deallocate(intstatus)
endif

deallocate(Morton_aim,Morton_aim_p)
if (npe>1) then
  call MPI_BARRIER(icomm,ierrmpi)
endif

end subroutine ImageDataVtk_mpi
!============================================================================
subroutine punstructuredvtk_mpi(qunit)

! Write one pvtu and vtu files for each processor
! Otherwise like unstructuredvtk_mpi

use mod_forest, only: Morton_start, Morton_stop, sfc_to_igrid
use mod_global_parameters

integer, intent(in) ::    qunit
!
double precision, dimension(0:nw+nwauxio)                   :: normconv
double precision, dimension(ixMlo1-1:ixMhi1,ixMlo2-1:ixMhi2,&
   ndim)         :: xC_TMP
double precision, dimension(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
   ndim)           :: xCC_TMP
double precision, dimension(ixMlo1-1:ixMhi1,ixMlo2-1:ixMhi2,&
   ndim)         :: xC
double precision, dimension(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ndim)           :: xCC
double precision, dimension(ixMlo1-1:ixMhi1,ixMlo2-1:ixMhi2,&
   nw+nwauxio)   :: wC_TMP
double precision, dimension(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
   nw+nwauxio)     :: wCC_TMP
character(len=name_len)   :: wnamei(1:nw+nwauxio),&
   xandwnamei(1:ndim+nw+nwauxio)
character(len=1024) :: outfilehead
integer             :: nx1,nx2,nxC1,nxC2,nc,np, igrid,ixCmin1,ixCmin2,ixCmax1,&
   ixCmax2,ixCCmin1,ixCCmin2,ixCCmax1,ixCCmax2,level,Morton_no
character(len=80)   :: pfilename
integer             :: filenr
logical             :: fileopen,conv_grid
!----------------------------------------------------------------------------

! Write pvtu-file:
if (mype==0) then
   call write_pvtu(qunit)
endif
! Now write the Source files:

inquire(qunit,opened=fileopen)
if(.not.fileopen)then
   ! generate filename
   filenr=snapshotini
   if (autoconvert) filenr=snapshotnext
   ! Open the file for the header part
   write(pfilename,'(a,i4.4,a,i4.4,a)') TRIM(base_filename),filenr,"p",mype,&
      ".vtu"
   open(qunit,file=pfilename,status='unknown',form='formatted')
endif
! generate xml header
write(qunit,'(a)')'<?xml version="1.0"?>'
write(qunit,'(a)',advance='no') '<VTKFile type="UnstructuredGrid"'
write(qunit,'(a)')' version="0.1" byte_order="LittleEndian">'
write(qunit,'(a)')'  <UnstructuredGrid>'
write(qunit,'(a)')'<FieldData>'
write(qunit,'(2a)')'<DataArray type="Float32" Name="TIME" ',&
   'NumberOfTuples="1" format="ascii">'
write(qunit,*) real(global_time*time_convert_factor)
write(qunit,'(a)')'</DataArray>'
write(qunit,'(a)')'</FieldData>'

call getheadernames(wnamei,xandwnamei,outfilehead)

! number of cells, number of corner points, per grid.
nx1=ixMhi1-ixMlo1+1;nx2=ixMhi2-ixMlo2+1;
nxC1=nx1+1;nxC2=nx2+1;
nc=nx1*nx2
np=nxC1*nxC2

! Note: using the w_write, writelevel, writespshift
! we can clip parts of the grid away, select variables, levels etc.
do level=levmin,levmax
   if (.not.writelevel(level)) cycle
   do Morton_no=Morton_start(mype),Morton_stop(mype)
    igrid=sfc_to_igrid(Morton_no)
    if (node(plevel_,igrid)/=level) cycle

    ! only output a grid when fully within clipped region selected
    ! by writespshift array
    conv_grid=(rnode(rpxmin1_,igrid)>=xprobmin1+&
       (xprobmax1-xprobmin1)*writespshift(1,1).and.rnode(rpxmin2_,&
       igrid)>=xprobmin2+(xprobmax2-xprobmin2)*writespshift(2,&
       1)).and.(rnode(rpxmax1_,igrid)<=xprobmax1-(xprobmax1-&
       xprobmin1)*writespshift(1,2).and.rnode(rpxmax2_,&
       igrid)<=xprobmax2-(xprobmax2-xprobmin2)*writespshift(2,2))
    if (.not.conv_grid) cycle

    call calc_x(igrid,xC,xCC)
    call calc_grid(qunit,igrid,xC,xCC,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,normconv,&
       ixCmin1,ixCmin2,ixCmax1,ixCmax2,ixCCmin1,ixCCmin2,ixCCmax1,ixCCmax2,&
       .true.)

    call write_vtk(qunit,ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixCmin1,ixCmin2,ixCmax1,&
       ixCmax2,ixCCmin1,ixCCmin2,ixCCmax1,ixCCmax2,igrid,nc,np,nx1,nx2,nxC1,&
       nxC2,normconv,wnamei,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP)
   end do ! Morton_no loop
end do ! level loop

 write(qunit,'(a)')'  </UnstructuredGrid>'
 write(qunit,'(a)')'</VTKFile>'
 close(qunit)

if (npe>1) then
  call MPI_BARRIER(icomm,ierrmpi)
endif
end subroutine punstructuredvtk_mpi
!============================================================================
subroutine unstructuredvtk_mpi(qunit)

! output for vtu format to paraview, non-binary version output
! parallel, uses calc_grid to compute nwauxio variables
! allows renormalizing using convert factors
! allows skipping of w_write selected variables

! implementation such that length of ASCII output is identical when
! run on 1 versus multiple CPUs (however, the order of the vtu pieces can differ)

use mod_forest, only: Morton_start, Morton_stop, sfc_to_igrid
use mod_global_parameters

integer, intent(in) ::    qunit

double precision ::  x_VTK(1:3)

double precision, dimension(ixMlo1-1:ixMhi1,ixMlo2-1:ixMhi2,ndim) :: xC_TMP,&
   xC_TMP_recv
double precision, dimension(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ndim)   :: xCC_TMP,&
   xCC_TMP_recv
double precision, dimension(ixMlo1-1:ixMhi1,ixMlo2-1:ixMhi2,ndim) :: xC
double precision, dimension(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ndim)   :: xCC

double precision, dimension(ixMlo1-1:ixMhi1,ixMlo2-1:ixMhi2,&
   nw+nwauxio)   :: wC_TMP,wC_TMP_recv
double precision, dimension(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
   nw+nwauxio)     :: wCC_TMP,wCC_TMP_recv
double precision, dimension(0:nw+nwauxio)                   :: normconv
integer::               igrid,iigrid,level,ixCmin1,ixCmin2,ixCmax1,ixCmax2,&
   ixCCmin1,ixCCmin2,ixCCmax1,ixCCmax2
integer::               NumGridsOnLevel(1:nlevelshi)
integer :: nx1,nx2,nxC1,nxC2,nodesonlevel,elemsonlevel,nc,np,ix1,ix2

character(len=80)::  filename
integer ::           filenr

integer, allocatable :: intstatus(:,:)

character(len=name_len) :: wnamei(1:nw+nwauxio),xandwnamei(1:ndim+nw+nwauxio)
character(len=1024) :: outfilehead

logical :: fileopen,conv_grid,cond_grid_recv
integer :: ipe,Morton_no,siz_ind
integer :: ind_send(4*2),ind_recv(4*2)
integer :: levmin_recv,levmax_recv,level_recv,igrid_recv,ixrvCmin1,ixrvCmin2,&
   ixrvCmax1,ixrvCmax2,ixrvCCmin1,ixrvCCmin2,ixrvCCmax1,ixrvCCmax2
!-----------------------------------------------------------------------------
if (mype==0) then
 inquire(qunit,opened=fileopen)
 if(.not.fileopen)then
    ! generate filename
    filenr=snapshotini
    if (autoconvert) filenr=snapshotnext
    write(filename,'(a,i4.4,a)') TRIM(base_filename),filenr,".vtu"
   ! Open the file for the header part
   open(qunit,file=filename,status='unknown',form='formatted')
 endif
 ! generate xml header
 write(qunit,'(a)')'<?xml version="1.0"?>'
 write(qunit,'(a)',advance='no') '<VTKFile type="UnstructuredGrid"'
 write(qunit,'(a)')' version="0.1" byte_order="LittleEndian">'
 write(qunit,'(a)')'<UnstructuredGrid>'
 write(qunit,'(a)')'<FieldData>'
 write(qunit,'(2a)')'<DataArray type="Float32" Name="TIME" ',&
    'NumberOfTuples="1" format="ascii">'
 write(qunit,*) real(global_time*time_convert_factor)
 write(qunit,'(a)')'</DataArray>'
 write(qunit,'(a)')'</FieldData>'
end if

call getheadernames(wnamei,xandwnamei,outfilehead)
! number of cells, number of corner points, per grid.
nx1=ixMhi1-ixMlo1+1;nx2=ixMhi2-ixMlo2+1;
nxC1=nx1+1;nxC2=nx2+1;
nc=nx1*nx2
np=nxC1*nxC2

! all slave processors send their minmal/maximal levels
if  (mype/=0) then
 if (Morton_stop(mype)==0) call mpistop("nultag")
 itag=1000*Morton_stop(mype)
 !print *,'ype,itag for levmin=',mype,itag,levmin
 call MPI_SEND(levmin,1,MPI_INTEGER, 0,itag,icomm,ierrmpi)
 itag=2000*Morton_stop(mype)
 !print *,'mype,itag for levmax=',mype,itag,levmax
 call MPI_SEND(levmax,1,MPI_INTEGER, 0,itag,icomm,ierrmpi)
end if


! Note: using the w_write, writelevel, writespshift
! we can clip parts of the grid away, select variables, levels etc.
do level=levmin,levmax
   if (.not.writelevel(level)) cycle
   do Morton_no=Morton_start(mype),Morton_stop(mype)
    igrid=sfc_to_igrid(Morton_no)
    if (mype/=0)then
      itag=Morton_no
      call MPI_SEND(igrid,1,MPI_INTEGER, 0,itag,icomm,ierrmpi)
      itag=igrid
      call MPI_SEND(node(plevel_,igrid),1,MPI_INTEGER, 0,itag,icomm,ierrmpi)
    end if
    if (node(plevel_,igrid)/=level) cycle

    ! only output a grid when fully within clipped region selected
    ! by writespshift array
    conv_grid=(rnode(rpxmin1_,igrid)>=xprobmin1+&
       (xprobmax1-xprobmin1)*writespshift(1,1).and.rnode(rpxmin2_,&
       igrid)>=xprobmin2+(xprobmax2-xprobmin2)*writespshift(2,&
       1)).and.(rnode(rpxmax1_,igrid)<=xprobmax1-(xprobmax1-&
       xprobmin1)*writespshift(1,2).and.rnode(rpxmax2_,&
       igrid)<=xprobmax2-(xprobmax2-xprobmin2)*writespshift(2,2))
    if (mype/=0)then
      call MPI_SEND(conv_grid,1,MPI_LOGICAL,0,itag,icomm,ierrmpi)
    end if
    if (.not.conv_grid) cycle

    call calc_x(igrid,xC,xCC)
    call calc_grid(qunit,igrid,xC,xCC,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,normconv,&
       ixCmin1,ixCmin2,ixCmax1,ixCmax2,ixCCmin1,ixCCmin2,ixCCmax1,ixCCmax2,&
       .true.)

    if (mype/=0) then
       itag=Morton_no
       ind_send=(/ ixCmin1,ixCmin2,ixCmax1,ixCmax2,ixCCmin1,ixCCmin2,ixCCmax1,&
          ixCCmax2 /)
       siz_ind=4*2
       call MPI_SEND(ind_send,siz_ind,MPI_INTEGER, 0,itag,icomm,ierrmpi)
       call MPI_SEND(normconv,nw+nwauxio+1,MPI_DOUBLE_PRECISION, 0,itag,icomm,&
          ierrmpi)

       call MPI_SEND(wC_TMP,1,type_block_wc_io, 0,itag,icomm,ierrmpi)
       call MPI_SEND(xC_TMP,1,type_block_xc_io, 0,itag,icomm,ierrmpi)
       itag=igrid
       call MPI_SEND(wCC_TMP,1,type_block_wcc_io, 0,itag,icomm,ierrmpi)
       call MPI_SEND(xCC_TMP,1,type_block_xcc_io, 0,itag,icomm,ierrmpi)
    else
       call write_vtk(qunit,ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixCmin1,ixCmin2,&
          ixCmax1,ixCmax2,ixCCmin1,ixCCmin2,ixCCmax1,ixCCmax2,igrid,nc,np,nx1,&
          nx2,nxC1,nxC2,normconv,wnamei,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP)
    end if
   end do ! Morton_no loop
end do ! level loop


if (mype==0) then
 allocate(intstatus(MPI_STATUS_SIZE,1))
 if(npe>1)then
  do ipe=1,npe-1
   itag=1000*Morton_stop(ipe)
   call MPI_RECV(levmin_recv,1,MPI_INTEGER, ipe,itag,icomm,intstatus(:,1),&
      ierrmpi)

   itag=2000*Morton_stop(ipe)
   call MPI_RECV(levmax_recv,1,MPI_INTEGER, ipe,itag,icomm,intstatus(:,1),&
      ierrmpi)

   do level=levmin_recv,levmax_recv
    if (.not.writelevel(level)) cycle
    do  Morton_no=Morton_start(ipe),Morton_stop(ipe)
     itag=Morton_no
     call MPI_RECV(igrid_recv,1,MPI_INTEGER, ipe,itag,icomm,intstatus(:,1),&
        ierrmpi)
     itag=igrid_recv
     call MPI_RECV(level_recv,1,MPI_INTEGER, ipe,itag,icomm,intstatus(:,1),&
        ierrmpi)
     if (level_recv/=level) cycle

     call MPI_RECV(cond_grid_recv,1,MPI_LOGICAL, ipe,itag,icomm,intstatus(:,1),&
        ierrmpi)
     if(.not.cond_grid_recv)cycle

     itag=Morton_no
     siz_ind=4*2
     call MPI_RECV(ind_recv,siz_ind, MPI_INTEGER, ipe,itag,icomm,intstatus(:,&
        1),ierrmpi)
     ixrvCmin1=ind_recv(1);ixrvCmin2=ind_recv(2);ixrvCmax1=ind_recv(2+1)
     ixrvCmax2=ind_recv(2+2);
     ixrvCCmin1=ind_recv(2*2+1);ixrvCCmin2=ind_recv(2*2+2)
     ixrvCCmax1=ind_recv(3*2+1);ixrvCCmax2=ind_recv(3*2+2);
     call MPI_RECV(normconv,nw+nwauxio+1, MPI_DOUBLE_PRECISION,ipe,itag,icomm,&
        intstatus(:,1),ierrmpi)

     call MPI_RECV(wC_TMP_recv,1,type_block_wc_io, ipe,itag,icomm,intstatus(:,&
        1),ierrmpi)
     call MPI_RECV(xC_TMP_recv,1,type_block_xc_io, ipe,itag,icomm,intstatus(:,&
        1),ierrmpi)

     itag=igrid_recv
     call MPI_RECV(wCC_TMP_recv,1,type_block_wcc_io, ipe,itag,icomm,&
        intstatus(:,1),ierrmpi)
     call MPI_RECV(xCC_TMP_recv,1,type_block_xcc_io, ipe,itag,icomm,&
        intstatus(:,1),ierrmpi)
     call write_vtk(qunit,ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixrvCmin1,ixrvCmin2,&
        ixrvCmax1,ixrvCmax2,ixrvCCmin1,ixrvCCmin2,ixrvCCmax1,ixrvCCmax2,&
        igrid_recv,nc,np,nx1,nx2,nxC1,nxC2,normconv,wnamei,xC_TMP_recv,&
        xCC_TMP_recv,wC_TMP_recv,wCC_TMP_recv)
    enddo ! Morton_no loop
   enddo ! level loop
  enddo ! processor loop
 endif ! multiple processors
 write(qunit,'(a)')'</UnstructuredGrid>'
 write(qunit,'(a)')'</VTKFile>'
 close(qunit)
endif

if (npe>1) then
  call MPI_BARRIER(icomm,ierrmpi)
  if(mype==0)deallocate(intstatus)
endif

end subroutine unstructuredvtk_mpi
!============================================================================
subroutine write_vtk(qunit,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,&
   ixCmax1,ixCmax2,ixCCmin1,ixCCmin2,ixCCmax1,ixCCmax2,igrid,nc,np,nx1,nx2,&
   nxC1,nxC2,normconv,wnamei,xC,xCC,wC,wCC)

use mod_global_parameters

integer, intent(in) :: qunit
integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,ixCmax1,&
   ixCmax2,ixCCmin1,ixCCmin2,ixCCmax1,ixCCmax2
integer, intent(in) :: igrid,nc,np,nx1,nx2,nxC1,nxC2
double precision, intent(in) :: normconv(0:nw+nwauxio)
character(len=name_len), intent(in)::  wnamei(1:nw+nwauxio)

double precision, dimension(ixMlo1-1:ixMhi1,ixMlo2-1:ixMhi2,ndim) :: xC
double precision, dimension(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ndim)   :: xCC

double precision, dimension(ixMlo1-1:ixMhi1,ixMlo2-1:ixMhi2,&
   nw+nwauxio)   :: wC
double precision, dimension(ixMlo1:ixMhi1,ixMlo2:ixMhi2,nw+nwauxio)     :: wCC

double precision ::  x_VTK(1:3)
integer :: iw,ix1,ix2,icel,VTK_type
!----------------------------------------------------------------------------

select case(convert_type)
    case('vtumpi','pvtumpi')
         ! we write out every grid as one VTK PIECE
      write(qunit,'(a,i7,a,i7,a)') '<Piece NumberOfPoints="',np,&
         '" NumberOfCells="',nc,'">'
      write(qunit,'(a)')'<PointData>'
      do iw=1,nw+nwauxio
         if(iw<=nw) then
            if(.not.w_write(iw)) cycle
         endif

            write(qunit,'(a,a,a)')'<DataArray type="Float64" Name="',&
               TRIM(wnamei(iw)),'" format="ascii">'
            write(qunit,'(200(1pe14.6))') ((wC(ix1,ix2,iw)*normconv(iw),&
               ix1=ixCmin1,ixCmax1),ix2=ixCmin2,ixCmax2)
            write(qunit,'(a)')'</DataArray>'
      enddo
      write(qunit,'(a)')'</PointData>'

      write(qunit,'(a)')'<Points>'
      write(qunit,'(a)'&
         )'<DataArray type="Float32" NumberOfComponents="3" format="ascii">'
         ! write cell corner coordinates in a backward dimensional loop, always 3D output
      do ix2=ixCmin2,ixCmax2 
      do ix1=ixCmin1,ixCmax1 
            x_VTK(1:3)=zero;
            x_VTK(1:ndim)=xC(ix1,ix2,1:ndim)*normconv(0);
            write(qunit,'(3(1pe14.6))') x_VTK
      end do 
      end do 
      write(qunit,'(a)')'</DataArray>'
      write(qunit,'(a)')'</Points>'

    case('vtuCCmpi','pvtuCCmpi')
      ! we write out every grid as one VTK PIECE
      write(qunit,'(a,i7,a,i7,a)') '<Piece NumberOfPoints="',np,&
         '" NumberOfCells="',nc,'">'
      write(qunit,'(a)')'<CellData>'
      do iw=1,nw+nwauxio
         if(iw<=nw) then
            if(.not.w_write(iw)) cycle
         endif
            write(qunit,'(a,a,a)')'<DataArray type="Float64" Name="',&
               TRIM(wnamei(iw)),'" format="ascii">'
            write(qunit,'(200(1pe14.6))') ((wCC(ix1,ix2,iw)*normconv(iw),&
               ix1=ixCCmin1,ixCCmax1),ix2=ixCCmin2,ixCCmax2)
            write(qunit,'(a)')'</DataArray>'
      enddo
      write(qunit,'(a)')'</CellData>'

      write(qunit,'(a)')'<Points>'
      write(qunit,'(a)'&
         )'<DataArray type="Float32" NumberOfComponents="3" format="ascii">'
      ! write cell corner coordinates in a backward dimensional loop, always 3D output
      do ix2=ixCmin2,ixCmax2 
      do ix1=ixCmin1,ixCmax1 
            x_VTK(1:3)=zero;
            x_VTK(1:ndim)=xC(ix1,ix2,1:ndim)*normconv(0);
            write(qunit,'(3(1pe14.6))') x_VTK
      end do 
      end do 
      write(qunit,'(a)')'</DataArray>'
      write(qunit,'(a)')'</Points>'
end select

write(qunit,'(a)')'<Cells>'

! connectivity part
write(qunit,'(a)'&
   )'<DataArray type="Int32" Name="connectivity" format="ascii">'
call save_connvtk(qunit,igrid)
write(qunit,'(a)')'</DataArray>'

! offsets data array
write(qunit,'(a)')'<DataArray type="Int32" Name="offsets" format="ascii">'
do icel=1,nc
    write(qunit,'(i7)') icel*(2**2)
end do
write(qunit,'(a)')'</DataArray>'

! VTK cell type data array
write(qunit,'(a)')'<DataArray type="Int32" Name="types" format="ascii">'
! VTK_LINE=3; VTK_PIXEL=8; VTK_VOXEL=11 -> vtk-syntax

 VTK_type=8 

do icel=1,nc
   write(qunit,'(i2)') VTK_type
enddo
write(qunit,'(a)')'</DataArray>'

write(qunit,'(a)')'</Cells>'

write(qunit,'(a)')'</Piece>'

end subroutine write_vtk
!============================================================================
subroutine write_vti(qunit,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,&
   ixCmax1,ixCmax2,ixCCmin1,ixCCmin2,ixCCmax1,ixCCmax2,ig1,ig2,nx1,nx2,&
   normconv,wnamei,wC,wCC)
use mod_global_parameters

integer, intent(in) :: qunit
integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,ixCmax1,&
   ixCmax2,ixCCmin1,ixCCmin2,ixCCmax1,ixCCmax2
integer, intent(in) :: ig1,ig2,nx1,nx2
double precision, intent(in) :: normconv(0:nw+nwauxio)
character(len=name_len), intent(in)::  wnamei(1:nw+nwauxio)

double precision, dimension(ixMlo1-1:ixMhi1,ixMlo2-1:ixMhi2,&
   nw+nwauxio)   :: wC
double precision, dimension(ixMlo1:ixMhi1,ixMlo2:ixMhi2,nw+nwauxio)     :: wCC

integer :: iw,ix1,ix2
integer :: extent(1:6)
!----------------------------------------------------------------------------

extent = 0
 extent(1*2-1) = (ig1-1) * nx1;  extent(2*2-1) = (ig2-1) * nx2;
 extent(1*2)   = (ig1)   * nx1;  extent(2*2)   = (ig2)   * nx2;


select case(convert_type)
    case('vtimpi','pvtimpi')
         ! we write out every grid as one VTK PIECE
      write(qunit,'(a,6(i10),a)') '<Piece Extent="',extent,'">'
      write(qunit,'(a)')'<PointData>'
      do iw=1,nw+nwauxio
         if(iw<=nw) then
            if(.not.w_write(iw)) cycle
         endif

            write(qunit,'(a,a,a)')'<DataArray type="Float64" Name="',&
               TRIM(wnamei(iw)),'" format="ascii">'
            write(qunit,'(200(1pe20.12))') ((wC(ix1,ix2,iw)*normconv(iw),&
               ix1=ixCmin1,ixCmax1),ix2=ixCmin2,ixCmax2)
            write(qunit,'(a)')'</DataArray>'
      enddo
      write(qunit,'(a)')'</PointData>'

    case('vtiCCmpi','pvtiCCmpi')
      ! we write out every grid as one VTK PIECE
      write(qunit,'(a,6(i10),a)') '<Piece Extent="',extent,'">'
      write(qunit,'(a)')'<CellData>'
      do iw=1,nw+nwauxio
         if(iw<=nw) then
            if(.not.w_write(iw)) cycle
         endif
            write(qunit,'(a,a,a)')'<DataArray type="Float64" Name="',&
               TRIM(wnamei(iw)),'" format="ascii">'
            write(qunit,'(200(1pe20.12))') ((wCC(ix1,ix2,iw)*normconv(iw),&
               ix1=ixCCmin1,ixCCmax1),ix2=ixCCmin2,ixCCmax2)
            write(qunit,'(a)')'</DataArray>'
      enddo
      write(qunit,'(a)')'</CellData>'
end select

write(qunit,'(a)')'</Piece>'

end subroutine write_vti
!=============================================================================
subroutine write_pvtu(qunit)

use mod_global_parameters

integer, intent(in) :: qunit

character(len=name_len)   :: wnamei(1:nw+nwauxio),&
   xandwnamei(1:ndim+nw+nwauxio),outtype
character(len=1024) :: outfilehead
character(len=80)   :: filename,pfilename
integer             :: filenr,iw,ipe,iscalars
logical             :: fileopen

select case(convert_type)
case('pvtumpi','pvtuBmpi')
   outtype="PPointData"
case('pvtuCCmpi','pvtuBCCmpi')
   outtype="PCellData"
end select
inquire(qunit,opened=fileopen)
if(.not.fileopen)then
   ! generate filename
   filenr=snapshotini
   if (autoconvert) filenr=snapshotnext
   write(filename,'(a,i4.4,a)') TRIM(base_filename),filenr,".pvtu"
   ! Open the file
   open(qunit,file=filename,status='unknown',form='formatted')
endif

call getheadernames(wnamei,xandwnamei,outfilehead)

! Get the default selection:
iscalars=1
do iw=nw,1, -1
   if (w_write(iw)) iscalars=iw
end do


! generate xml header
write(qunit,'(a)')'<?xml version="1.0"?>'
write(qunit,'(a)',advance='no') '<VTKFile type="PUnstructuredGrid"'
write(qunit,'(a)')' version="0.1" byte_order="LittleEndian">'
write(qunit,'(a)')'  <PUnstructuredGrid GhostLevel="0">'
! Either celldata or pointdata:
write(qunit,'(a,a,a,a,a)')'    <',TRIM(outtype),' Scalars="',&
   TRIM(wnamei(iscalars))//'">'
do iw=1,nw
   if(.not.w_write(iw))cycle
   write(qunit,'(a,a,a)')'      <PDataArray type="Float32" Name="',&
      TRIM(wnamei(iw)),'"/>'
end do
do iw=nw+1,nw+nwauxio
   write(qunit,'(a,a,a)')'      <PDataArray type="Float32" Name="',&
      TRIM(wnamei(iw)),'"/>'
end do
write(qunit,'(a,a,a)')'    </',TRIM(outtype),'>'

write(qunit,'(a)')'    <PPoints>'
write(qunit,'(a)')'      <PDataArray type="Float32" NumberOfComponents="3"/>'
write(qunit,'(a)')'    </PPoints>'

do ipe=0,npe-1
   write(pfilename,'(a,i4.4,a,i4.4,a)') TRIM(base_filename(INDEX &
      (base_filename, '/', BACK = .TRUE.)+1:LEN(base_filename))),filenr,"p",&
      ipe,".vtu"
   write(qunit,'(a,a,a)')'    <Piece Source="',TRIM(pfilename),'"/>'
end do
write(qunit,'(a)')'  </PUnstructuredGrid>'
write(qunit,'(a)')'</VTKFile>'

close(qunit)

end subroutine write_pvtu
!=============================================================================
subroutine tecplot_mpi(qunit)

! output for tecplot (ASCII format)
! parallel, uses calc_grid to compute nwauxio variables
! allows renormalizing using convert factors

! the current implementation is such that tecplotmpi and tecplotCCmpi will
! create different length output ASCII files when used on 1 versus multiple CPUs
! in fact, on 1 CPU, there will be as many zones as there are levels
! on multiple CPUs, there will be a number of zones up to the number of
! levels times the number of CPUs (can be less, when some level not on a CPU)

use mod_forest, only: Morton_start, Morton_stop, sfc_to_igrid
use mod_global_parameters

integer, intent(in) :: qunit

integer::               igrid,iigrid,level,igonlevel,iw,idim,ix1,ix2
integer::               NumGridsOnLevel(1:nlevelshi)
integer :: nx1,nx2,nxC1,nxC2,nodesonlevel,elemsonlevel,ixCmin1,ixCmin2,ixCmax1,&
   ixCmax2,ixCCmin1,ixCCmin2,ixCCmax1,ixCCmax2
integer :: nodesonlevelmype,elemsonlevelmype

integer ::              nodes, elems

integer, allocatable :: intstatus(:,:)

double precision :: x_TEC(ndim), w_TEC(nw+nwauxio)

double precision, dimension(ixMlo1-1:ixMhi1,ixMlo2-1:ixMhi2,ndim) :: xC_TMP,&
   xC_TMP_recv
double precision, dimension(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ndim)   :: xCC_TMP,&
   xCC_TMP_recv
double precision, dimension(ixMlo1-1:ixMhi1,ixMlo2-1:ixMhi2,ndim) :: xC
double precision, dimension(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ndim)   :: xCC

double precision, dimension(ixMlo1-1:ixMhi1,ixMlo2-1:ixMhi2,&
   nw+nwauxio)   :: wC_TMP,wC_TMP_recv
double precision, dimension(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
   nw+nwauxio)     :: wCC_TMP,wCC_TMP_recv
double precision, dimension(0:nw+nwauxio)                   :: normconv
logical :: fileopen,first
integer :: Morton_no,ipe,levmin_recv,levmax_recv,igrid_recv,level_recv
integer :: ixrvCmin1,ixrvCmin2,ixrvCmax1,ixrvCmax2,ixrvCCmin1,ixrvCCmin2,&
   ixrvCCmax1,ixrvCCmax2
integer :: ind_send(2*2),ind_recv(2*2),siz_ind,igonlevel_recv
integer :: NumGridsOnLevel_mype(1:nlevelshi,0:npe-1)
character(len=80) :: filename
integer ::           filenr
character(len=1024) :: tecplothead

character(len=name_len) :: wnamei(1:nw+nwauxio),xandwnamei(1:ndim+nw+nwauxio)
character(len=1024) :: outfilehead
!-----------------------------------------------------------------------------
if(nw/=count(w_write(1:nw)))then
 if(mype==0) PRINT *,'tecplot_mpi does not use w_write=F'
 call mpistop('w_write, tecplot')
end if

if(nocartesian)then
 if(mype==0) PRINT *,'tecplot_mpi with nocartesian and typeaxial=',typeaxial
endif

Master_cpu_open : if (mype == 0) then
 inquire(qunit,opened=fileopen)
 if (.not.fileopen) then
   ! generate filename
    filenr=snapshotini
    if (autoconvert) filenr=snapshotnext
   write(filename,'(a,i4.4,a)') TRIM(base_filename),filenr,".plt"
   open(qunit,file=filename,status='unknown')
 end if

 call getheadernames(wnamei,xandwnamei,outfilehead)

 write(tecplothead,'(a)') "VARIABLES = "//TRIM(outfilehead)
 write(qunit,'(a)') tecplothead(1:len_trim(tecplothead))
end if  Master_cpu_open


! determine overall number of grids per level, and the same info per CPU
NumGridsOnLevel(1:nlevelshi)=0
do level=levmin,levmax
   NumGridsOnLevel(level)=0
   do Morton_no=Morton_start(mype),Morton_stop(mype)
      igrid = sfc_to_igrid(Morton_no)
      if (node(plevel_,igrid)/=level) cycle
      NumGridsOnLevel(level)=NumGridsOnLevel(level)+1
   end do
   NumGridsOnLevel_mype(level,0:npe-1)=0
   NumGridsOnLevel_mype(level,mype) = NumGridsOnLevel(level)
   call MPI_ALLREDUCE(MPI_IN_PLACE,NumGridsOnLevel_mype(level,0:npe-1),npe,&
      MPI_INTEGER,MPI_MAX,icomm,ierrmpi)
   call MPI_ALLREDUCE(MPI_IN_PLACE,NumGridsOnLevel(level),1,MPI_INTEGER,&
      MPI_SUM, icomm,ierrmpi)
end do





nx1=ixMhi1-ixMlo1+1;nx2=ixMhi2-ixMlo2+1;
nxC1=nx1+1;nxC2=nx2+1;

if(mype==0.and.npe>1) allocate(intstatus(MPI_STATUS_SIZE,1))



if  (mype/=0) then
 itag=1000*Morton_stop(mype)
 call MPI_SEND(levmin,1,MPI_INTEGER, 0,itag,icomm,ierrmpi)
 itag=2000*Morton_stop(mype)
 call MPI_SEND(levmax,1,MPI_INTEGER, 0,itag,icomm,ierrmpi)
end if

do level=levmin,levmax
   nodesonlevelmype=NumGridsOnLevel_mype(level,mype)*nxC1*nxC2
   elemsonlevelmype=NumGridsOnLevel_mype(level,mype)*nx1*nx2
   nodesonlevel=NumGridsOnLevel(level)*nxC1*nxC2
   elemsonlevel=NumGridsOnLevel(level)*nx1*nx2
   ! for all tecplot variants coded up here, we let the TECPLOT ZONES coincide
   ! with the AMR grid LEVEL. Other options would be
   !    let each grid define a zone: inefficient for TECPLOT internal workings
   !       hence not implemented
   !    let entire octree define 1 zone: no difference in interpolation
   !       properties across TECPLOT zones detected as yet, hence not done
   select case(convert_type)
     case('tecplotmpi')
       ! in this option, we store the corner coordinates, as well as the corner
       ! values of all variables (obtained by averaging). This allows POINT packaging,
       ! and thus we can save full grid info by using one call to calc_grid
       if (mype==0.and.(nodesonlevelmype>0.and.elemsonlevelmype>0))write(qunit,&
          "(a,i7,a,a,i7,a,i7,a,f25.16,a,a)") 'ZONE T="',level,'"',', N=',&
          nodesonlevelmype,', E=',elemsonlevelmype, ', SOLUTIONTIME=',&
          global_time*time_convert_factor,', DATAPACKING=POINT, ZONETYPE=',&
            'FEQUADRILATERAL'
      do Morton_no=Morton_start(mype),Morton_stop(mype)
         igrid = sfc_to_igrid(Morton_no)
         if (mype/=0)then
           itag=Morton_no
           call MPI_SEND(igrid,1,MPI_INTEGER, 0,itag,icomm,ierrmpi)
           itag=igrid
           call MPI_SEND(node(plevel_,igrid),1,MPI_INTEGER, 0,itag,icomm,&
              ierrmpi)
         end if
         if (node(plevel_,igrid)/=level) cycle
         call calc_x(igrid,xC,xCC)
         call calc_grid(qunit,igrid,xC,xCC,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,&
            normconv,ixCmin1,ixCmin2,ixCmax1,ixCmax2,ixCCmin1,ixCCmin2,&
            ixCCmax1,ixCCmax2,.true.)
         if (mype/=0) then
            itag=Morton_no
            ind_send=(/ ixCmin1,ixCmin2,ixCmax1,ixCmax2 /)
            siz_ind=2*2
            call MPI_SEND(ind_send,siz_ind, MPI_INTEGER, 0,itag,icomm,ierrmpi)
            call MPI_SEND(normconv,nw+nwauxio+1,MPI_DOUBLE_PRECISION, 0,itag,&
               icomm,ierrmpi)

            call MPI_SEND(wC_TMP,1,type_block_wc_io, 0,itag,icomm,ierrmpi)
            call MPI_SEND(xC_TMP,1,type_block_xc_io, 0,itag,icomm,ierrmpi)
         else
           do ix2=ixCmin2,ixCmax2
           do ix1=ixCmin1,ixCmax1
              x_TEC(1:ndim)=xC_TMP(ix1,ix2,1:ndim)*normconv(0)
              w_TEC(1:nw+nwauxio)=wC_TMP(ix1,ix2,&
                 1:nw+nwauxio)*normconv(1:nw+nwauxio)
              write(qunit,fmt="(100(e14.6))") x_TEC, w_TEC
           end do
           end do
         end if
       enddo

     case('tecplotCCmpi')
       ! in this option, we store the corner coordinates, and the cell center
       ! values of all variables. Due to this mix of corner/cell center, we must
       ! use BLOCK packaging, and thus we have enormous overhead by using
       ! calc_grid repeatedly to merely fill values of cell corner coordinates
       ! and cell center values per dimension, per variable
       if(ndim+nw+nwauxio>99) call mpistop(&
          "adjust format specification in writeout")
       if(nw+nwauxio==1)then
         ! to make tecplot happy: avoid [ndim+1-ndim+1] in varlocation varset
         ! and just set [ndim+1]
         if (mype==0.and.(nodesonlevelmype>0.and.&
            elemsonlevelmype>0))write(qunit,&
            "(a,i7,a,a,i7,a,i7,a,f25.16,a,i1,a,a)") 'ZONE T="',level,'"',&
            ', N=',nodesonlevelmype,', E=',elemsonlevelmype, ', SOLUTIONTIME=',&
            global_time*time_convert_factor,&
            ', DATAPACKING=BLOCK, VARLOCATION=([', ndim+1,&
            ']=CELLCENTERED), ZONETYPE=',  'FEQUADRILATERAL'
       else
        if(ndim+nw+nwauxio<10) then
         ! difference only in length of integer format specification for ndim+nw+nwauxio
         if (mype==0.and.(nodesonlevelmype>0.and.&
            elemsonlevelmype>0))write(qunit,&
            "(a,i7,a,a,i7,a,i7,a,f25.16,a,i1,a,i1,a,a)") 'ZONE T="',level,'"',&
            ', N=',nodesonlevelmype,', E=',elemsonlevelmype, ', SOLUTIONTIME=',&
            global_time*time_convert_factor,&
            ', DATAPACKING=BLOCK, VARLOCATION=([', ndim+1,'-',ndim+nw+nwauxio,&
            ']=CELLCENTERED), ZONETYPE=',  'FEQUADRILATERAL'
        else
         if (mype==0.and.(nodesonlevelmype>0.and.&
            elemsonlevelmype>0))write(qunit,&
            "(a,i7,a,a,i7,a,i7,a,f25.16,a,i1,a,i2,a,a)") 'ZONE T="',level,'"',&
            ', N=',nodesonlevelmype,', E=',elemsonlevelmype, ', SOLUTIONTIME=',&
            global_time*time_convert_factor,&
            ', DATAPACKING=BLOCK, VARLOCATION=([', ndim+1,'-',ndim+nw+nwauxio,&
            ']=CELLCENTERED), ZONETYPE=',  'FEQUADRILATERAL'
        endif
       endif

       do idim=1,ndim
         first=(idim==1)
         do Morton_no=Morton_start(mype),Morton_stop(mype)
          igrid = sfc_to_igrid(Morton_no)
          if (mype/=0)then
           itag=Morton_no*idim
           call MPI_SEND(igrid,1,MPI_INTEGER, 0,itag,icomm,ierrmpi)
           itag=igrid*idim
           call MPI_SEND(node(plevel_,igrid),1,MPI_INTEGER, 0,itag,icomm,&
              ierrmpi)
          end if
          if (node(plevel_,igrid)/=level) cycle

          call calc_x(igrid,xC,xCC)
          call calc_grid(qunit,igrid,xC,xCC,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,&
             normconv,ixCmin1,ixCmin2,ixCmax1,ixCmax2,ixCCmin1,ixCCmin2,&
             ixCCmax1,ixCCmax2,first)
          if (mype/=0)then
            ind_send=(/ ixCmin1,ixCmin2,ixCmax1,ixCmax2 /)
            siz_ind=2*2
            itag=igrid*idim
            call MPI_SEND(ind_send,siz_ind, MPI_INTEGER, 0,itag,icomm,ierrmpi)
            call MPI_SEND(normconv,nw+nwauxio+1,MPI_DOUBLE_PRECISION, 0,itag,&
               icomm,ierrmpi)
            call MPI_SEND(xC_TMP,1,type_block_xc_io, 0,itag,icomm,ierrmpi)
          else
            write(qunit,fmt="(100(e14.6))") xC_TMP(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2,idim)*normconv(0)
          end if
         enddo
       enddo

       do iw=1,nw+nwauxio
        do Morton_no=Morton_start(mype),Morton_stop(mype)
         igrid = sfc_to_igrid(Morton_no)
         if (mype/=0)then
           itag=Morton_no*(ndim+iw)
           call MPI_SEND(igrid,1,MPI_INTEGER, 0,itag,icomm,ierrmpi)
           itag=igrid*(ndim+iw)
           call MPI_SEND(node(plevel_,igrid),1,MPI_INTEGER, 0,itag,icomm,&
              ierrmpi)
         end if
         if (node(plevel_,igrid)/=level) cycle

         call calc_x(igrid,xC,xCC)
         call calc_grid(qunit,igrid,xC,xCC,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,&
            normconv,ixCmin1,ixCmin2,ixCmax1,ixCmax2,ixCCmin1,ixCCmin2,&
            ixCCmax1,ixCCmax2,.true.)

         if (mype/=0)then
            ind_send=(/ ixCCmin1,ixCCmin2,ixCCmax1,ixCCmax2 /)
            siz_ind=2*2
            itag=igrid*(ndim+iw)
            call MPI_SEND(ind_send,siz_ind, MPI_INTEGER, 0,itag,icomm,ierrmpi)
            call MPI_SEND(normconv,nw+nwauxio+1,MPI_DOUBLE_PRECISION, 0,itag,&
               icomm,ierrmpi)
            call MPI_SEND(wCC_TMP,1,type_block_wcc_io, 0,itag,icomm,ierrmpi)
         else
            write(qunit,fmt="(100(e14.6))") wCC_TMP(ixCCmin1:ixCCmax1,&
               ixCCmin2:ixCCmax2,iw)*normconv(iw)
         endif
        enddo
       enddo
     case default
       call mpistop('no such tecplot type')
   end select


   igonlevel=0
   do Morton_no=Morton_start(mype),Morton_stop(mype)
      igrid = sfc_to_igrid(Morton_no)
      if (mype/=0)then
          itag=Morton_no
          call MPI_SEND(igrid,1,MPI_INTEGER, 0,itag,icomm,ierrmpi)
          itag=igrid
          call MPI_SEND(node(plevel_,igrid),1,MPI_INTEGER, 0,itag,icomm,&
             ierrmpi)
      end if
      if (node(plevel_,igrid)/=level) cycle

      igonlevel=igonlevel+1
      if (mype/=0)then
          itag=igrid
          call MPI_SEND(igonlevel,1,MPI_INTEGER, 0,itag,icomm,ierrmpi)
      end if
      if(mype==0)then
        call save_conntec(qunit,igrid,igonlevel)
      endif
   end do
end do

if (mype==0) then
 if (npe>1) then
  do ipe=1,npe-1
   itag=1000*Morton_stop(ipe)
   call MPI_RECV(levmin_recv,1,MPI_INTEGER, ipe,itag,icomm,intstatus(:,1),&
      ierrmpi)
   itag=2000*Morton_stop(ipe)
   call MPI_RECV(levmax_recv,1,MPI_INTEGER, ipe,itag,icomm,intstatus(:,1),&
      ierrmpi)
   do level=levmin_recv,levmax_recv
    nodesonlevelmype=NumGridsOnLevel_mype(level,ipe)*nxC1*nxC2
    elemsonlevelmype=NumGridsOnLevel_mype(level,ipe)*nx1*nx2
    nodesonlevel=NumGridsOnLevel(level)*nxC1*nxC2
    elemsonlevel=NumGridsOnLevel(level)*nx1*nx2
    select case(convert_type)
     case('tecplotmpi')
        ! in this option, we store the corner coordinates, as well as the corner
        ! values of all variables (obtained by averaging). This allows POINT packaging,
        ! and thus we can save full grid info by using one call to calc_grid
        if(nodesonlevelmype>0.and.elemsonlevelmype>0) write(qunit,&
           "(a,i7,a,a,i7,a,i7,a,f25.16,a,a)") 'ZONE T="',level,'"',', N=',&
           nodesonlevelmype,', E=',elemsonlevelmype, ', SOLUTIONTIME=',&
           global_time*time_convert_factor,', DATAPACKING=POINT, ZONETYPE=',&
             'FEQUADRILATERAL'
        do  Morton_no=Morton_start(ipe),Morton_stop(ipe)
         itag=Morton_no
         call MPI_RECV(igrid_recv,1,MPI_INTEGER, ipe,itag,icomm,intstatus(:,1),&
            ierrmpi)
         itag=igrid_recv
         call MPI_RECV(level_recv,1,MPI_INTEGER, ipe,itag,icomm,intstatus(:,1),&
            ierrmpi)
         if (level_recv/=level) cycle

         itag=Morton_no
         siz_ind=2*2
         call MPI_RECV(ind_recv,siz_ind, MPI_INTEGER, ipe,itag,icomm,&
            intstatus(:,1),ierrmpi)
         ixrvCmin1=ind_recv(1);ixrvCmin2=ind_recv(2);ixrvCmax1=ind_recv(2+1)
         ixrvCmax2=ind_recv(2+2);
         call MPI_RECV(normconv,nw+nwauxio+1, MPI_DOUBLE_PRECISION,ipe,itag,&
            icomm,intstatus(:,1),ierrmpi)

         call MPI_RECV(wC_TMP_recv,1,type_block_wc_io, ipe,itag,icomm,&
            intstatus(:,1),ierrmpi)
         call MPI_RECV(xC_TMP_recv,1,type_block_xc_io, ipe,itag,icomm,&
            intstatus(:,1),ierrmpi)

         do ix2=ixrvCmin2,ixrvCmax2
         do ix1=ixrvCmin1,ixrvCmax1
              x_TEC(1:ndim)=xC_TMP_RECV(ix1,ix2,1:ndim)*normconv(0)
              w_TEC(1:nw+nwauxio)=wC_TMP_RECV(ix1,ix2,&
                 1:nw+nwauxio)*normconv(1:nw+nwauxio)
              write(qunit,fmt="(100(e14.6))") x_TEC, w_TEC
         end do
         end do
        end do
     case('tecplotCCmpi')
       ! in this option, we store the corner coordinates, and the cell center
       ! values of all variables. Due to this mix of corner/cell center, we must
       ! use BLOCK packaging, and thus we have enormous overhead by using
       ! calc_grid repeatedly to merely fill values of cell corner coordinates
       ! and cell center values per dimension, per variable
       if(ndim+nw+nwauxio>99) call mpistop(&
          "adjust format specification in writeout")
       if(nw+nwauxio==1)then
         ! to make tecplot happy: avoid [ndim+1-ndim+1] in varlocation varset
         ! and just set [ndim+1]
         if(nodesonlevelmype>0.and.elemsonlevelmype>0) write(qunit,&
            "(a,i7,a,a,i7,a,i7,a,f25.16,a,i1,a,a)") 'ZONE T="',level,'"',&
            ', N=',nodesonlevelmype,', E=',elemsonlevelmype, ', SOLUTIONTIME=',&
            global_time*time_convert_factor,&
            ', DATAPACKING=BLOCK, VARLOCATION=([', ndim+1,&
            ']=CELLCENTERED), ZONETYPE=',  'FEQUADRILATERAL'
       else
        if(ndim+nw+nwauxio<10) then
         ! difference only in length of integer format specification for ndim+nw+nwauxio
         if(nodesonlevelmype>0.and.elemsonlevelmype>0) write(qunit,&
            "(a,i7,a,a,i7,a,i7,a,f25.16,a,i1,a,i1,a,a)") 'ZONE T="',level,'"',&
            ', N=',nodesonlevelmype,', E=',elemsonlevelmype, ', SOLUTIONTIME=',&
            global_time*time_convert_factor,&
            ', DATAPACKING=BLOCK, VARLOCATION=([', ndim+1,'-',ndim+nw+nwauxio,&
            ']=CELLCENTERED), ZONETYPE=',  'FEQUADRILATERAL'
        else
         if(nodesonlevelmype>0.and.elemsonlevelmype>0) write(qunit,&
            "(a,i7,a,a,i7,a,i7,a,f25.16,a,i1,a,i2,a,a)") 'ZONE T="',level,'"',&
            ', N=',nodesonlevelmype,', E=',elemsonlevelmype, ', SOLUTIONTIME=',&
            global_time*time_convert_factor,&
            ', DATAPACKING=BLOCK, VARLOCATION=([', ndim+1,'-',ndim+nw+nwauxio,&
            ']=CELLCENTERED), ZONETYPE=',  'FEQUADRILATERAL'
        endif
       endif

       do idim=1,ndim
         do  Morton_no=Morton_start(ipe),Morton_stop(ipe)
           itag=Morton_no*idim
           call MPI_RECV(igrid_recv,1,MPI_INTEGER, ipe,itag,icomm,intstatus(:,&
              1),ierrmpi)
           itag=igrid_recv*idim
           call MPI_RECV(level_recv,1,MPI_INTEGER, ipe,itag,icomm,intstatus(:,&
              1),ierrmpi)
           if (level_recv/=level) cycle

           siz_ind=2*2
           itag=igrid_recv*idim
           call MPI_RECV(ind_recv,siz_ind, MPI_INTEGER, ipe,itag,icomm,&
              intstatus(:,1),ierrmpi)
           ixrvCmin1=ind_recv(1);ixrvCmin2=ind_recv(2)
           ixrvCmax1=ind_recv(2+1);ixrvCmax2=ind_recv(2+2);
           call MPI_RECV(normconv,nw+nwauxio+1, MPI_DOUBLE_PRECISION,ipe,itag,&
              icomm,intstatus(:,1),ierrmpi)
           call MPI_RECV(xC_TMP_recv,1,type_block_xc_io, ipe,itag,icomm,&
              intstatus(:,1),ierrmpi)
           write(qunit,fmt="(100(e14.6))") xC_TMP_recv(ixrvCmin1:ixrvCmax1,&
              ixrvCmin2:ixrvCmax2,idim)*normconv(0)
         end do
       end do

       do iw=1,nw+nwauxio
        do Morton_no=Morton_start(ipe),Morton_stop(ipe)
           itag=Morton_no*(ndim+iw)
           call MPI_RECV(igrid_recv,1,MPI_INTEGER, ipe,itag,icomm,intstatus(:,&
              1),ierrmpi)
           itag=igrid_recv*(ndim+iw)
           call MPI_RECV(level_recv,1,MPI_INTEGER, ipe,itag,icomm,intstatus(:,&
              1),ierrmpi)
           if (level_recv/=level) cycle

           siz_ind=2*2
           itag=igrid_recv*(ndim+iw)
           call MPI_RECV(ind_recv,siz_ind, MPI_INTEGER, ipe,itag,icomm,&
              intstatus(:,1),ierrmpi)
           ixrvCCmin1=ind_recv(1);ixrvCCmin2=ind_recv(2)
           ixrvCCmax1=ind_recv(2+1);ixrvCCmax2=ind_recv(2+2);
           call MPI_RECV(normconv,nw+nwauxio+1, MPI_DOUBLE_PRECISION,ipe,itag,&
              icomm,intstatus(:,1),ierrmpi)
           call MPI_RECV(wCC_TMP_recv,1,type_block_wcc_io, ipe,itag,icomm,&
              intstatus(:,1),ierrmpi)
           write(qunit,fmt="(100(e14.6))") wCC_TMP_recv(ixrvCCmin1:ixrvCCmax1,&
              ixrvCCmin2:ixrvCCmax2,iw)*normconv(iw)
        enddo
       enddo
     case default
       call mpistop('no such tecplot type')
    end select

    do Morton_no=Morton_start(ipe),Morton_stop(ipe)
      itag=Morton_no
      call MPI_RECV(igrid_recv,1,MPI_INTEGER, ipe,itag,icomm,intstatus(:,1),&
         ierrmpi)
      itag=igrid_recv
      call MPI_RECV(level_recv,1,MPI_INTEGER, ipe,itag,icomm,intstatus(:,1),&
         ierrmpi)
      if (level_recv/=level) cycle

      itag=igrid_recv
      call MPI_RECV(igonlevel_recv,1,MPI_INTEGER, ipe,itag,icomm,intstatus(:,&
         1),ierrmpi)
      call save_conntec(qunit,igrid_recv,igonlevel_recv)
    end do ! morton loop
   end do ! level loop
  end do ! ipe loop
 end if ! npe>1 if
end if ! mype=0 if


if (npe>1) then
  call MPI_BARRIER(icomm,ierrmpi)
  if(mype==0)deallocate(intstatus)
endif

end subroutine tecplot_mpi
!=============================================================================
subroutine punstructuredvtkB_mpi(qunit)

! Write one pvtu and vtu files for each processor
! Otherwise like unstructuredvtk_mpi
! output for vtu format to paraview, binary version output
! uses calc_grid to compute nwauxio variables
! allows renormalizing using convert factors

use mod_forest, only: Morton_start, Morton_stop, sfc_to_igrid
use mod_global_parameters

integer, intent(in) ::    qunit

double precision ::  x_VTK(1:3)

double precision, dimension(ixMlo1-1:ixMhi1,ixMlo2-1:ixMhi2,ndim) :: xC_TMP
double precision, dimension(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ndim)   :: xCC_TMP
double precision, dimension(ixMlo1-1:ixMhi1,ixMlo2-1:ixMhi2,ndim) :: xC
double precision, dimension(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ndim)   :: xCC

double precision, dimension(ixMlo1-1:ixMhi1,ixMlo2-1:ixMhi2,&
   nw+nwauxio)   :: wC_TMP
double precision, dimension(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
   nw+nwauxio)     :: wCC_TMP

integer :: igrid,iigrid,level,igonlevel,icel,ixCmin1,ixCmin2,ixCmax1,ixCmax2,&
   ixCCmin1,ixCCmin2,ixCCmax1,ixCCmax2,Morton_no
integer ::               NumGridsOnLevel(1:nlevelshi)
integer :: nx1,nx2,nxC1,nxC2,nodesonlevel,elemsonlevel,nc,np,VTK_type,ix1,ix2
double precision :: normconv(0:nw+nwauxio)
character(len=80) :: pfilename
character(len=name_len) :: wnamei(1:nw+nwauxio),xandwnamei(1:ndim+nw+nwauxio)
character(len=1024) :: outfilehead

integer*8 :: offset
integer::  recsep,k,iw,filenr
integer::  length,lengthcc,offset_points,offset_cells, length_coords,&
   length_conn,length_offsets
character::  buf
character(len=6)::  bufform

logical ::   fileopen
!-----------------------------------------------------------------------------

! Write pvtu-file:
if (mype==0) then
   call write_pvtu(qunit)
endif
! Now write the Source files:
inquire(qunit,opened=fileopen)
if(.not.fileopen)then
   ! generate filename
   filenr=snapshotini
   if (autoconvert) filenr=snapshotnext
   ! Open the file for the header part
   write(pfilename,'(a,i4.4,a,i4.4,a)') TRIM(base_filename),filenr,"p",mype,&
      ".vtu"
   open(qunit,file=pfilename,status='unknown',form='formatted')
endif
! generate xml header
write(qunit,'(a)')'<?xml version="1.0"?>'
write(qunit,'(a)',advance='no') '<VTKFile type="UnstructuredGrid"'
write(qunit,'(a)')' version="0.1" byte_order="LittleEndian">'
write(qunit,'(a)')'  <UnstructuredGrid>'
write(qunit,'(a)')'<FieldData>'
write(qunit,'(2a)')'<DataArray type="Float32" Name="TIME" ',&
   'NumberOfTuples="1" format="ascii">'
write(qunit,*) real(global_time*time_convert_factor)
write(qunit,'(a)')'</DataArray>'
write(qunit,'(a)')'</FieldData>'


offset=0
recsep=4

call getheadernames(wnamei,xandwnamei,outfilehead)

! number of cells, number of corner points, per grid.
nx1=ixMhi1-ixMlo1+1;nx2=ixMhi2-ixMlo2+1;
nxC1=nx1+1;nxC2=nx2+1;
nc=nx1*nx2
np=nxC1*nxC2

length=np*size_real
lengthcc=nc*size_real

length_coords=3*length
length_conn=2**2*size_int*nc
length_offsets=nc*size_int

! Note: using the w_write, writelevel, writespshift
! we can clip parts of the grid away, select variables, levels etc.
do level=levmin,levmax
 if (writelevel(level)) then
   do Morton_no=Morton_start(mype),Morton_stop(mype)
    igrid=sfc_to_igrid(Morton_no)
    if (node(plevel_,igrid)/=level) cycle
    ! only output a grid when fully within clipped region selected
    ! by writespshift array
    if ((rnode(rpxmin1_,igrid)>=xprobmin1+(xprobmax1-xprobmin1)*writespshift(1,&
       1).and.rnode(rpxmin2_,igrid)>=xprobmin2+&
       (xprobmax2-xprobmin2)*writespshift(2,1)).and.(rnode(rpxmax1_,&
       igrid)<=xprobmax1-(xprobmax1-xprobmin1)*writespshift(1,&
       2).and.rnode(rpxmax2_,igrid)<=xprobmax2-(xprobmax2-&
       xprobmin2)*writespshift(2,2))) then
      select case(convert_type)
       case('pvtuBmpi')
         ! we write out every grid as one VTK PIECE
         write(qunit,'(a,i7,a,i7,a)') '<Piece NumberOfPoints="',np,&
            '" NumberOfCells="',nc,'">'
         write(qunit,'(a)')'<PointData>'
         do iw=1,nw
            if(.not.w_write(iw))cycle

            write(qunit,'(a,a,a,i16,a)')'<DataArray type="Float32" Name="',&
               TRIM(wnamei(iw)), '" format="appended" offset="',offset,'">'
            write(qunit,'(a)')'</DataArray>'
            offset=offset+length+size_int
         enddo
         do iw=nw+1,nw+nwauxio

            write(qunit,'(a,a,a,i16,a)')'<DataArray type="Float32" Name="',&
               TRIM(wnamei(iw)), '" format="appended" offset="',offset,'">'
            write(qunit,'(a)')'</DataArray>'
            offset=offset+length+size_int
         enddo
         write(qunit,'(a)')'</PointData>'

         write(qunit,'(a)')'<Points>'
         write(qunit,'(a,i16,a)')&
'<DataArray type="Float32" NumberOfComponents="3" format="appended" offset="',&
            offset,'"/>'
         ! write cell corner coordinates in a backward dimensional loop, always 3D output
         offset=offset+length_coords+size_int
         write(qunit,'(a)')'</Points>'
       case('pvtuBCCmpi')
         ! we write out every grid as one VTK PIECE
         write(qunit,'(a,i7,a,i7,a)') '<Piece NumberOfPoints="',np,&
            '" NumberOfCells="',nc,'">'
         write(qunit,'(a)')'<CellData>'
         do iw=1,nw
            if(.not.w_write(iw))cycle

            write(qunit,'(a,a,a,i16,a)')'<DataArray type="Float32" Name="',&
               TRIM(wnamei(iw)), '" format="appended" offset="',offset,'">'
            write(qunit,'(a)')'</DataArray>'
            offset=offset+lengthcc+size_int
         enddo
         do iw=nw+1,nw+nwauxio

            write(qunit,'(a,a,a,i16,a)')'<DataArray type="Float32" Name="',&
               TRIM(wnamei(iw)), '" format="appended" offset="',offset,'">'
            write(qunit,'(a)')'</DataArray>'
            offset=offset+lengthcc+size_int
         enddo
         write(qunit,'(a)')'</CellData>'

         write(qunit,'(a)')'<Points>'
         write(qunit,'(a,i16,a)')&
'<DataArray type="Float32" NumberOfComponents="3" format="appended" offset="',&
            offset,'"/>'
         ! write cell corner coordinates in a backward dimensional loop, always 3D output
         offset=offset+length_coords+size_int
         write(qunit,'(a)')'</Points>'
      end select


      write(qunit,'(a)')'<Cells>'

      ! connectivity part
      write(qunit,'(a,i16,a)'&
)'<DataArray type="Int32" Name="connectivity" format="appended" offset="',&
         offset,'"/>'
      offset=offset+length_conn+size_int

      ! offsets data array
      write(qunit,'(a,i16,a)')&
          '<DataArray type="Int32" Name="offsets" format="appended" offset="',&
         offset,'"/>'
      offset=offset+length_offsets+size_int

      ! VTK cell type data array
      write(qunit,'(a,i16,a)')&
          '<DataArray type="Int32" Name="types" format="appended" offset="',&
         offset,'"/>'
      offset=offset+size_int+nc*size_int

      write(qunit,'(a)')'</Cells>'

      write(qunit,'(a)')'</Piece>'
    endif
   enddo
 endif
enddo

write(qunit,'(a)')'</UnstructuredGrid>'
write(qunit,'(a)')'<AppendedData encoding="raw">'

close(qunit)
! next to make gfortran compiler happy, as it does not know
! form='binary' and produces error on compilation
!bufform='binary'
!open(qunit,file=pfilename,form=bufform,position='append')
!This should in principle do also for gfortran (tested with gfortran 4.6.0 and Intel 11.1):
open(qunit,file=pfilename,access='stream',form='unformatted',&
   position='append')
buf='_'
write(qunit) TRIM(buf)

do level=levmin,levmax
 if (writelevel(level)) then
   do Morton_no=Morton_start(mype),Morton_stop(mype)
    igrid=sfc_to_igrid(Morton_no)
    if (node(plevel_,igrid)/=level) cycle
      ! only output a grid when fully within clipped region selected
      ! by writespshift array
      if ((rnode(rpxmin1_,igrid)>=xprobmin1+&
         (xprobmax1-xprobmin1)*writespshift(1,1).and.rnode(rpxmin2_,&
         igrid)>=xprobmin2+(xprobmax2-xprobmin2)*writespshift(2,&
         1)).and.(rnode(rpxmax1_,igrid)<=xprobmax1-(xprobmax1-&
         xprobmin1)*writespshift(1,2).and.rnode(rpxmax2_,&
         igrid)<=xprobmax2-(xprobmax2-xprobmin2)*writespshift(2,2))) then
        call calc_x(igrid,xC,xCC)
        call calc_grid(qunit,igrid,xC,xCC,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,&
           normconv,ixCmin1,ixCmin2,ixCmax1,ixCmax2,ixCCmin1,ixCCmin2,ixCCmax1,&
           ixCCmax2,.true.)
        do iw=1,nw
          if(.not.w_write(iw))cycle
          select case(convert_type)
            case('pvtuBmpi')
              write(qunit) length
              write(qunit) ((real(wC_TMP(ix1,ix2,iw)*normconv(iw)),ix1=ixCmin1,&
                 ixCmax1),ix2=ixCmin2,ixCmax2)
            case('pvtuBCCmpi')
              write(qunit) lengthcc
              write(qunit) ((real(wCC_TMP(ix1,ix2,iw)*normconv(iw)),&
                 ix1=ixCCmin1,ixCCmax1),ix2=ixCCmin2,ixCCmax2)
          end select
        enddo
        do iw=nw+1,nw+nwauxio
          select case(convert_type)
            case('pvtuBmpi')
              write(qunit) length
              write(qunit) ((real(wC_TMP(ix1,ix2,iw)*normconv(iw)),ix1=ixCmin1,&
                 ixCmax1),ix2=ixCmin2,ixCmax2)
            case('pvtuBCCmpi')
              write(qunit) lengthcc
              write(qunit) ((real(wCC_TMP(ix1,ix2,iw)*normconv(iw)),&
                 ix1=ixCCmin1,ixCCmax1),ix2=ixCCmin2,ixCCmax2)
          end select
        enddo

        write(qunit) length_coords
        do ix2=ixCmin2,ixCmax2 
        do ix1=ixCmin1,ixCmax1 
          x_VTK(1:3)=zero;
          x_VTK(1:ndim)=xC_TMP(ix1,ix2,1:ndim)*normconv(0);
          do k=1,3
           write(qunit) real(x_VTK(k))
          end do
        end do 
        end do 

        write(qunit) length_conn
        do ix2=1,nx2
        do ix1=1,nx1
        
        
        write(qunit)(ix2-1)*nxC1+ix1-1, (ix2-1)*nxC1+ix1,ix2*nxC1+ix1-1,&
           ix2*nxC1+ix1
         
        
        end do
        end do

        write(qunit) length_offsets
        do icel=1,nc
           write(qunit) icel*(2**2)
        end do


       
        VTK_type=8 
       
        write(qunit) size_int*nc
        do icel=1,nc
         write(qunit) VTK_type
       enddo
    endif
  end do
 endif
end do

close(qunit)
open(qunit,file=pfilename,status='unknown',form='formatted',position='append')

write(qunit,'(a)')'</AppendedData>'
write(qunit,'(a)')'</VTKFile>'
close(qunit)

end subroutine punstructuredvtkB_mpi
!=============================================================================
subroutine calc_x(igrid,xC,xCC)
  use mod_global_parameters

  integer, intent(in)               :: igrid
  double precision, intent(out)     :: xC(ixMlo1-1:ixMhi1,ixMlo2-1:ixMhi2,&
     ndim)
  double precision, intent(out)     :: xCC(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ndim)
  ! .. local ..
  integer                           :: ixCmin1,ixCmin2,ixCmax1,ixCmax2, idims,&
      level, ix

  level=node(plevel_,igrid)

  ! coordinates of cell centers
  xCC(ixMlo1:ixMhi1,ixMlo2:ixMhi2,:)=pw(igrid)%x(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
     :)

  ! coordinates of cell corners
  ixCmin1=ixMlo1-1;ixCmin2=ixMlo2-1; ixCmax1=ixMhi1;ixCmax2=ixMhi2;
  if(slab)then
     do idims=1,ndim
       xC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idims)=pw(igrid)%x(ixCmin1:ixCmax1,&
          ixCmin2:ixCmax2,idims)+0.5d0*dx(idims,level)
     end do
  else
     ! for any non-cartesian or stretched coordinate (allow multiple stretched directions)
     do ix=ixCmin1,ixCmax1
       xC(ix,ixCmin2:ixCmax2,1)=pw(igrid)%x(ix,ixCmin2:ixCmax2,&
          1)+0.5d0*pw(igrid)%dx(ix,ixCmin2:ixCmax2,1)
     end do
     do ix=ixCmin2,ixCmax2
       xC(ixCmin1:ixCmax1,ix,2)=pw(igrid)%x(ixCmin1:ixCmax1,ix,&
          2)+0.5d0*pw(igrid)%dx(ixCmin1:ixCmax1,ix,2)
     end do
  endif

end subroutine calc_x
!=============================================================================
