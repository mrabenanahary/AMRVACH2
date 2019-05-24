!============================================================================
Module mod_interface_metric
!============================================================================
type first_derive
  logical                                        :: elm_on
  logical                                        :: elm_isone
  logical                                        :: elm_needed
  double precision, dimension(:^D&), allocatable :: dwg
end type first_derive

type lapse_metric
 double precision, dimension(:^D&), allocatable  :: wg
 type(first_derive), dimension(:),  pointer      :: drv
 logical                                         :: elm_on
 logical                                         :: elm_isone
 logical                                         :: elm_needed
end type lapse_metric

type conformal_decomposition_metric
 double precision, dimension(:^D&), allocatable  :: wg
 type(first_derive), dimension(:),  pointer      :: drv
 logical                                         :: elm_on
 logical                                         :: elm_isone
 logical                                         :: elm_needed
end type conformal_decomposition_metric

type plapse_metric
 type (lapse_metric), pointer                    :: lpspt
end type plapse_metric
type shift_metric
 logical                                         :: elm_on
 logical                                         :: elm_isone
 logical                                         :: elm_needed
 double precision, dimension(:^D&), allocatable  :: wg
 type(first_derive), dimension(:),  pointer      :: drv
end type shift_metric

type pshift_metric
 type(shift_metric), pointer                  :: sftpt
end type pshift_metric

type elements_metric
 logical                                         :: elm_on
 logical                                         :: elm_isone
 logical                                         :: elm_needed
 double precision, dimension(:^D&), allocatable  :: wg 
 type(first_derive), dimension(:),  pointer      :: drv
end type elements_metric

type pelements_metric
 type(elements_metric),pointer                   :: elpt
 logical                                         :: isset
end type pelements_metric


type determinant_metric
 double precision, dimension(:^D&), allocatable  :: Gama,g
 logical                                         :: Gama_isone,G_isone
 logical                                         :: Gama_needed,G_needed
end type determinant_metric


type themetric
 type(pelements_metric), dimension(:,:),pointer          :: elem
 type(pelements_metric), dimension(:,:),pointer          :: elem_inv
 type(pelements_metric), dimension(:,:),pointer          :: elemsp
 type(pelements_metric), dimension(:,:),pointer          :: elemsp_inv
 type(pelements_metric), dimension(:,:),pointer          :: elemsp_updwn
 type(pelements_metric), dimension(:,:),pointer          :: Qsp
 type(pelements_metric), dimension(:,:),pointer          :: Qsp_inv
 type(pshift_metric), dimension(:), pointer              :: bt
 type(pshift_metric), dimension(:), pointer              :: bt_cont
 type (lapse_metric), pointer                            :: alfa
 type(determinant_metric), pointer                       :: sqrtdetr
 type(conformal_decomposition_metric)                    :: scalCD
 logical                                                 :: cell_center
 character(len=20)				         :: names
end type themetric 

type(themetric),allocatable,target,save  :: pm_cell(:),pm_face^D(:),&
                                            pmCoarse(:),pmCoarse_face^D(:)

type(themetric), pointer,save            :: mypm,mypm_cell,mypm_face^D


type tensor_wF
 logical                    :: elem_needed
 double precision, pointer  :: w(:^D&)
end type tensor_wF


type tensor_BLtoKS
 logical                                         :: elm_on
 logical                                         :: elm_isone
 double precision, dimension(:^D&), allocatable  :: wg
end type tensor_BLtoKS
type(tensor_BLtoKS), save, allocatable           :: MblTOks(:,:)
!============================================================================
End module mod_interface_metric
!============================================================================
