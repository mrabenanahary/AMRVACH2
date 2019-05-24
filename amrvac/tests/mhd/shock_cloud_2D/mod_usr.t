!=============================================================================
! amrvacusr.t.shockcloud
!=============================================================================
module mod_usr
  use mod_mhd
  implicit none
  double precision :: beta,chi,nval,rcore,rbound,xshock,machs

contains

  subroutine usr_init()
    usr_set_parameters=> initglobaldata_usr
    usr_init_one_grid => initonegrid_usr
    usr_special_bc    => specialbound_usr
    usr_aux_output    => specialvar_output
    usr_add_aux_names => specialvarnames_output 

    call set_coordinate_system("Cartesian")

    call mhd_activate()

  end subroutine usr_init

  subroutine initglobaldata_usr
    mhd_gamma=1.66666667d0
    mhd_eta=zero
    mhd_etah=zero

    chi=10.0d0
    nval=8.0d0
    rcore=0.62d0
    rbound=1.77d0
    xshock=-2.66d0
    machs=10.d0

    select case(iprob)
      case(1,11)
        beta=1.0d0
      case(2,12)
        beta=10.0d0
      case(3,13)
        beta=0.5d0
      case default
        call mpistop('iprob to implement')
    endselect
  
  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixG^L,ix^L,w,x)
  ! initialize one grid
    integer, intent(in) :: ixG^L,ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    double precision:: rval(ixG^S),pleft,rholeft,vleft,Prat,alfa,cleft
    double precision:: pright,rhoright,vright,deltax
    logical, save :: first=.true.

    ! uniform cloud surroundings
    pleft=one
    rholeft=one+(chi-one)/(one+(rbound/rcore)**nval)
    vleft=zero

    ! set uniform parallel B field
    w(ix^S,mag(1))=dsqrt(two*pleft/beta)
    w(ix^S,mag(2))=zero
    {^IFTHREED
    w(ix^S,mag(3))=zero
    }
    
    ! set trivial velocity components
    w(ix^S,mom(2))=zero
    {^IFTHREED
    w(ix^S,mom(3))=zero
    }
    
    ! compute the RH related states
    Prat=one/(one+(machs**2-one)*two*mhd_gamma/(mhd_gamma+one))
    alfa=(mhd_gamma+1)/(mhd_gamma-one)
    cleft=dsqrt(mhd_gamma*pleft/rholeft)
    rhoright=rholeft*(alfa+Prat)/(alfa*Prat+one)
    pright=pleft/Prat
    vright=cleft*machs*(one-(alfa*Prat+one)/(alfa+Prat))
    
    deltax=2.0d-2
    select case(iprob)
      case(11,12,13)
        w(ix^S,rho_)=rhoright+(rholeft-rhoright) &
                *half*(dtanh((x(ix^S,1)-xshock)/deltax)+one)
        w(ix^S,mom(1)) =vright  +(vleft-vright)     &
                *half*(dtanh((x(ix^S,1)-xshock)/deltax)+one)
        w(ix^S,p_)  =pright  +(pleft-pright)     &
                *half*(dtanh((x(ix^S,1)-xshock)/deltax)+one)
      case(1,2,3)
        where(x(ix^S,1)<xshock)
          w(ix^S,rho_)=rhoright
          w(ix^S,mom(1)) =vright
          w(ix^S,p_)  =pright
        elsewhere
          w(ix^S,rho_)=rholeft
          w(ix^S,mom(1)) =vleft
          w(ix^S,p_)  =pleft
        endwhere
      case default
        call mpistop('iprob to implement')
    end select
    
    ! overwrite cloud region density, set tracer values
    rval(ix^S)=dsqrt(^D&x(ix^S,^D)**2+)
    where(rval(ix^S)<rbound)
      w(ix^S,rho_)=one+(chi-one)/(one+(rval(ix^S)/rcore)**nval)
    endwhere
    if(mhd_n_tracer>0) then
      where(rval(ix^S)<rbound)
        w(ix^S,tracer(1))=w(ix^S,rho_)
      elsewhere
        w(ix^S,tracer(1))=zero
      endwhere
    end if
    
    if(mype==0.and.first)then
       write(*,*)'Doing shock-cloud challenge, ideal MHD'
       write(*,*)'density is=',rholeft,' pressure is=',pleft
       write(*,*)'plasma beta is set to=',beta
       write(*,*)'shock Mach is=',machs
       first=.false.
    endif

    call mhd_to_conserved(ixG^L,ix^L,w,x)

  end subroutine initonegrid_usr

  subroutine specialbound_usr(qt,ixG^L,ixO^L,iB,w,x)
    ! special boundary types, user defined
    ! user must assign conservative variables in bounderies
    integer, intent(in) :: ixG^L, ixO^L, iB
    double precision, intent(in) :: qt, x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    double precision:: pleft,rholeft,Prat,alfa,cleft

    pleft=one
    rholeft=one+(chi-one)/(one+(rbound/rcore)**nval)
    Prat=one/(one+(machs**2-one)*two*mhd_gamma/(mhd_gamma+one))
    alfa=(mhd_gamma+1)/(mhd_gamma-one)
    cleft=dsqrt(mhd_gamma*pleft/rholeft)
    
    select case(iB)
     case(1)
      ! fix the postshock values
      w(ixO^S,rho_)=rholeft*(alfa+Prat)/(alfa*Prat+one)
      w(ixO^S,p_)=pleft/Prat
      w(ixO^S,mag(1))=dsqrt(two*pleft/beta)
      w(ixO^S,mag(2))=zero
      {^IFTHREED
      w(ixO^S,mag(3))=zero
      }
      w(ixO^S,mom(1))=cleft*machs*(one-(alfa*Prat+one)/(alfa+Prat))
      w(ixO^S,mom(2))=zero
      {^IFTHREED
      w(ixO^S,mom(3))=zero
      }
      if(mhd_n_tracer>0) w(ixO^S,tracer(1))=zero
      call mhd_to_conserved(ixG^L,ixO^L,w,x)
     case default
       call mpistop('boundary not defined')
    end select

  end subroutine specialbound_usr

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
  ! this subroutine can be used in convert, to add auxiliary variables to the
  ! converted output file, for further analysis using tecplot, paraview, ....
  ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
  !
  ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
  ! corresponding normalization values (default value 1)
    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    double precision                   :: tmp(ixI^S) 

    call mhd_get_pthermal(w,x,ixI^L,ixO^L,tmp)
    ! output the temperature p/rho
    w(ixO^S,nw+1)=tmp(ixO^S)/w(ixO^S,rho_)
    !! output the plasma beta p*2/B**2
    if(B0field)then
      w(ixO^S,nw+2)=tmp(ixO^S)*two/sum((w(ixO^S,mag(:))+&
                    block%B0(ixO^S,:,0))**2,dim=ndim+1)
    else
      w(ixO^S,nw+2)=tmp(ixO^S)*two/sum(w(ixO^S,mag(:))**2,dim=ndim+1)
    endif
    ! output divB1
    call get_divb(w,ixI^L,ixO^L,tmp)
    w(ixO^S,nw+3)=tmp(ixO^S)
    
  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
  ! newly added variables need to be concatenated with the w_names/primnames string
    character(len=*) :: varnames
    varnames='Te beta divb'

  end subroutine specialvarnames_output

end module mod_usr
