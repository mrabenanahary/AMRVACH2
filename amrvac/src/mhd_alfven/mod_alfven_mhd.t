module mod_mhd
  use mod_mhd_alfven_phys
  use mod_mhd_alfven_hllc
  use mod_mhd_alfven_roe
  use mod_mhd_alfven_ppm

  use mod_amrvac

  implicit none
  public

contains

  subroutine mhd_alfven_activate()
    call mhd_alfven_phys_init()
    call mhd_alfven_hllc_init()
    call mhd_alfven_roe_init()
    call mhd_alfven_ppm_init()
  end subroutine mhd_alfven_activate

end module mod_mhd
