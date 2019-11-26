!> Module for physical and numeric constants
!> Created: 01.09.2012 Oliver Porth (Physical constants)
!> improved by Z. Meliani  add iso_fortran_env
module mod_constants
  use, intrinsic :: iso_fortran_env
  implicit none
  public

  ! Real variables setting
  integer, parameter :: sp = REAL32
  integer, parameter :: dp = REAL64
  integer, parameter :: qp = REAL128

  ! A very small real number (but above machine precision)
  real(kind=dp), parameter :: smalldouble = 1.D-12

  !> A very large real number
  real(kind=dp), parameter :: bigdouble = 1.D+99

  !> A very large integer
  integer, parameter :: biginteger = 10000000

  !> \todo Remove these
  real(kind=dp), parameter :: zero    = 0.0d0
  real(kind=dp), parameter :: one     = 1.0d0
  real(kind=dp), parameter :: two     = 2.0d0
  real(kind=dp), parameter :: half    = 0.5d0
  real(kind=dp), parameter :: quarter = 0.25d0
  real(kind=dp), parameter :: third   = 1/3.0d0

  !> Pi
  real(kind=dp), parameter :: dpi= 4.0_dp*datan(1.0_dp)
  !3.141592653589793238462643383279502884197169399375105d0

  !> Proton mass in cgs
  real(kind=dp), parameter :: mp_cgs = 1.672621777d-24 ! g

  !> Hydrogen mass in cgs
  real(kind=dp), parameter :: mH_cgs = 1.6733D-24 ! g

  !> Boltzmann constant in cgs
  real(kind=dp), parameter :: kB_cgs = 1.3806488d-16 ! erg K^-1

  !> Proton mass in SI
  real(kind=dp), parameter :: mp_SI = 1.672621777d-27 ! kg

  !> Boltzmann constant in SI
  real(kind=dp), parameter :: kB_SI = 1.3806488d-23 ! J K^-1

  !> Permeability in SI
  real(kind=dp), parameter :: miu0_SI = 1.2566370614d-6 ! H m^-1

  real(kind=dp), PARAMETER :: const_c     = 2.99792458d10   ! cm s^-1           ; Speed of light
  real(kind=dp), PARAMETER :: const_me    = 9.1093897d-28   ! g                 ; Electron mass
  real(kind=dp), PARAMETER :: const_mp    = 1.672621777d-24 ! g                 ; Proton mass
  real(kind=dp), PARAMETER :: const_e     = 4.8032068d-10   ! g^1/2 cm^3/2 s^-1 ; Electron charge
  real(kind=dp), PARAMETER :: const_MSun  = 1.98892d33      ! g                 ; Solar mass
  real(kind=dp), PARAMETER :: const_kB    = 1.3806488d-16   ! erg K^-1          ; Boltzmann constant
  real(kind=dp), PARAMETER :: const_h     = 6.6260755d-27   ! erg s             ; Planck constant
  ! Conversion factors:
  real(kind=dp), PARAMETER :: const_eV    = 1.6021772d-12   ! erg/eV            ; Electron volt
  real(kind=dp), PARAMETER :: const_Tera  = 1.d12           ! -                 ; Tera
  real(kind=dp), PARAMETER :: const_Peta  = 1.d15           ! -                 ; Peta
  real(kind=dp), PARAMETER :: const_years = 3.1536d7        ! s year^-1         ; seconds in a year
  real(kind=dp), PARAMETER :: const_sigmaT= 6.6524d-25      ! cm^2              ; thomson cross section
end module mod_constants
