      module histo_data
      implicit none
      integer                     :: maxbin
      integer,dimension(30)       :: freq
      end module histo_data


      module input_data
      implicit none
      integer,parameter           :: maxFrag = 99
      integer                     :: nFragments
      integer,dimension(maxFrag)  :: nFrozen,nVec
      real(kind=8)                :: threshold
      logical                     :: extra_info,all_epsilons
      logical                     :: debug,noAvgCore
      character (len=20)          :: project
      end module input_data

      module input_data2
      implicit none
      integer   :: nbase,nstates,extraDim,nBlocks
      integer   :: reduced_extraDim
      integer, allocatable  :: ncombv(:,:),dimens(:),blocks(:,:)
      integer,allocatable  :: ovlp_mebfs(:)
      real (kind=8), allocatable   :: ecorr(:)
      real (kind=8), allocatable   :: h(:,:),s(:,:)
      character(len=80)     :: project
      logical :: extra,dressed_coupling
      end module input_data2
      
      module fileInfo
      implicit none
      character(len=80)     :: title
      character(len=25)     :: vecfile1,vecfile2
      end module

      module genInfo
      implicit none
      integer,allocatable    :: nBas(:)
      integer,allocatable    :: nOcc(:)
      integer                :: nIrrep,mxBas,totOrb
      end module

      module decisions
      implicit none
      logical                :: add,symmetry_is_lowered
      end module

      module coefficients
      implicit none
      real(kind=8),allocatable  :: coeff(:,:,:)
      end module
