      module grackle_header

      USE ISO_C_BINDING
      implicit None
#include "grackle.def"

#ifdef CONFIG_BFLOAT_4
#define R_PREC_KIND 4
#endif
#ifdef CONFIG_BFLOAT_8
#define R_PREC_KIND 8
#endif

c     AMRVAC-DEFINED environment type for grackle use

      integer, parameter :: gr_RKIND=RKIND
      integer, parameter :: gr_DKIND=DKIND
      integer, parameter :: gr_DIKIND=DIKIND
      integer, parameter :: gr_rpknd=R_PREC_KIND
      R_PREC, parameter :: gr_tiny = tiny
      R_PREC, parameter :: gr_huge = huge

c     Define storage for grackle units, fields and parameter data

      contains

      subroutine test_grackle_header()



      PRINT *, 'grackle_header:gr_RKIND='
      PRINT *, gr_RKIND
      PRINT *, 'grackle_header:gr_DKIND='
      PRINT *, gr_DKIND
      PRINT *, 'grackle_header:gr_DIKIND='
      PRINT *, gr_DIKIND
      PRINT *, 'grackle_header:R_PREC_KIND='
      PRINT *, gr_rpknd


      end subroutine test_grackle_header


      end module grackle_header
