      program common_MO
      implicit none

      integer, parameter              :: maxSym = 8
      integer, parameter              :: maxFrag = 99
      integer                         :: iFrag,nFragments
      integer                         :: iSym,nSym
      integer                         :: j,k
      integer                         :: nMaxBas
      integer, dimension(maxSym)      :: nBas,nOcc,nBasFrag            ! nOcc, number of MOs in fragment common basis
      integer, dimension(maxFrag)     :: nVec                          ! number of vector sets in a fragment
      integer, allocatable            :: startBas(:)
      integer, allocatable            :: startOrb(:)

      real (kind=8)                   :: threshold                     ! threshold for linear dependency in common MO basis
      real (kind=8), allocatable      :: superBasis(:,:,:)             ! final common MO basis to be used in MOTRA
      real (kind=8), allocatable      :: commonMOs(:,:,:)              ! common orbital basis of a fragment
      real (kind=8), allocatable      :: occNu(:,:)                    ! eigenvalues of the MO overlap matrix, used as dummy occupations
      real (kind=8), allocatable      :: sdiag(:,:)                    ! eigenvalues of the MO overlap matrix of the fragments

      character (len=20)              :: project

      logical                         :: gamess

      nVec = 0
      call readin(threshold,nFragments,nVec,gamess,project)
      write(*,'(A,E10.2)')'common MO basis with threshold:',threshold
      write(*,*)'R. K. Kathir, C. de Graaf, R. Broer, R. W. A. Havenith'
      write(*,*)'submitted, 2019'
      write(*,*)
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*                                                                                               *
*      IN                                                                                       *
*                                                                                               *
* RUNFILE      : Number of basis functions for the supermolecule                                *
* ONEINTX      : AO overlap matrix of fragment X                                                *
* RUNFILX      : Number of basis functions, etc. for each fragment                              *
* INPORB.X_Y   : Y different orbital sets of fragment X (one for each state o the fragment)     *
*                                                                                               *
*      OUT                                                                                      *
*                                                                                               *
* COMMONORB    : vectors of the common MO basis of the supermolecule (for motra)                *
* COMMON.X_Y   : vectors of the different states expressed in the common MO basis (for GronOR)  *
*                                                                                               *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

      nBas = 0
      call NameRun('RUNFILE')
      call Get_iScalar('nSym',nSym)
      Call Get_iArray('nBas',nBas,nSym)
      nMaxBas = maxval(nBas)
      allocate( superBasis(nSym,nMaxBas,nMaxBas) )
      allocate( commonMOs(nSym,nMaxBas,nMaxBas) )
      allocate( startOrb(nSym) )
      allocate( startBas(nSym) )
      allocate( sdiag(nSym,nMaxBas) )
      allocate( occNu(nSym,nMaxBas) )

      if ( gamess ) then
        open(12,file = 'common.basis')
      else
        open(12,file = 'COMMONORB')
        write(12,'(A11)') '#INPORB 2.2'
        write(12,'(A5)')  '#INFO'
        write(12,'(A48,E9.3)') 
     &    '* Common molecular orbital basis with tau_MO = ',threshold
        write(12,'(3I8)')0,nSym,0
        write(12,'(3I8)')(nBas(iSym),iSym = 1, nSym)
        write(12,'(3I8)')(nBas(iSym),iSym = 1, nSym)
        write(12,'(A4)') '#ORB'
      end if

      superBasis = 0.0
      occNu = 0.0
      startOrb = 0
      startBas = 0
      do iFrag = 1, nFragments
        commonMOs = 0.0
        sdiag = 0.0
        call common_basis(iFrag,nVec,commonMOs,nOcc,nBasFrag,          ! get the common basis for fragment iFrag, the eigen-
     &                       sdiag,threshold,gamess,project)           ! values of sMO will be used as occupation numbers
        do iSym = 1, nSym
          do j = 1, nOcc(iSym)
            do k = 1, nBasFrag(iSym)
              superBasis(iSym, startOrb(iSym) + j, startBas(iSym) + k) ! dump the common basis in the final set of 
     &                               = commonMOs(iSym,j,k)             ! MOs of the supermolecule
            end do
            occNu(iSym,startOrb + j) = sdiag(iSym,j)
          end do
          startOrb(iSym) = startOrb(iSym) + nOcc(iSym)
          startBas(iSym) = startBas(iSym) + nBasFrag(iSym)
        end do
      end do
      write(*,'(A,8I4)') 'Final number of orbitals in MO basis    : ',
     &                       (startOrb(iSym), iSym = 1, nSym)
      write(*,'(A,8I4)') 'Sum of the number of AO basis functions : ',
     &                       (startBas(iSym), iSym = 1, nSym)
      write(*,'(A,8I4)') 'Deleted orbitals in MOTRA               : ',
     &                 (startBas(iSym) - startOrb(iSym), iSym = 1, nSym)
      if ( sum(nBas) .ne. sum(startBas) ) then
        write(*,*) 'Inconsistency detected'
        write(*,*) 'The number of AO basis functions of the',
     &             '  supermolecule is not equal to the sum'
        write(*,*) 'of the the number of AO basis functions',
     &             ' of the fragments'
        write(*,*) sum(nBAs),' not equal to ',sum(startBas)
      end if

* write the final MO basis of the supermolecule to a file
      if ( gamess ) then
        do iSym = 1, nSym
          write(12,'(A11,2I5)')'read input ',nBas(iSym),nbas(iSym)
          do j = 1, nBas(iSym)
            write(12,'(6F20.15)')(superBasis(iSym,j,k),k=1,nBas(iSym))
          end do
        end do
      else
        do iSym = 1, nSym
          do j = 1, nBas(iSym)
            write(12,'(A9,2I5)')'* ORBITAL',iSym,j
            write(12,'(5E22.14)')(superBasis(iSym,j,k),k=1,nBas(iSym))
          end do
        end do
        write(12,'(A4)') '#OCC'
        write(12,'(A20)') '* OCCUPATION NUMBERS'
        do iSym = 1, nSym
          write(12,'(5E22.14)') (occNu(iSym,j),j=1,nBas(iSym))
        end do
      end if
      close(12)
      deallocate(superBasis)
      deallocate(commonMOs)
      deallocate(startOrb)
      deallocate(startBas)
      deallocate(sdiag)
      deallocate(occNu)

      contains    ! all subroutines

      subroutine common_basis(iFrag,nVec,basisMotra,lDim,nBasFrag,sdiag,
     &                           threshold,gamess,project) 
* generation o fthe common MO basis for fragment iFrag using nVec vectors sets
      implicit none 
      integer, parameter              :: maxSym = 8
      integer                         :: iFrag,nSym,iSym,j,k,l,iVec
      integer                         :: iRC,iOpt,iComponent,iSymLbl
      integer                         :: iCounter,lwork,startVec
      integer                         :: nBasFragTot,nBasFragMax
      integer, intent(in)             :: nVec(:)
      integer, intent(out)            :: nBasFrag(:),lDim(:)           ! lDim counts the number of linear independent MOs in the basis
      integer, dimension(maxSym)      :: fragOrbTot,nOcc               ! nOcc, the non-empty orbitals of a vector set

      real (kind=8)                   :: threshold
      real (kind=8), intent(out)      :: basisMotra(:,:,:)             ! common MO basis expressed in AOs
      real (kind=8), intent(inout)    :: sdiag(:,:)
      real (kind=8), allocatable      :: linDep(:,:,:)                 ! adding all non-empty orbitals of the different vectos sets
      real (kind=8), allocatable      :: basis(:,:,:)                  ! common MO basis, eigenvectors of sMO
      real (kind=8), allocatable      :: vec(:,:,:)
      real (kind=8), allocatable      :: sMO(:,:,:),aMatrix(:,:)
      real (kind=8), allocatable      :: sAO(:,:,:)
      real (kind=8), allocatable      :: sU(:,:)                       ! Intermediate matrix in the transformation of the vectors to 
                                                                       ! the common MO basis
      real (kind=8), allocatable      :: VsU(:,:)                      ! The vectors of the different states expressed in the common MO basis
      real (kind=8), allocatable      :: aVector(:)
      real (kind=8), allocatable      :: work(:)

      character (len=12)              :: runfileName,oneintName
      character (len=25)              :: vecFileName
      character (len=20)              :: base,project

      logical                         :: gamess

      if ( iFrag .le. 9 ) then
        write(runfileName,'(A6,I1)')'RUNFIL',iFrag
        write(oneintName,'(A6,I1)') 'ONEINT',iFrag
      else
        write(runfileName,'(A6,I2)')'RUNFIL',iFrag
        write(oneintName,'(A6,I2)') 'ONEINT',iFrag
      end if
      call NameRun(runfileName)
      call Get_iScalar('nSym',nSym)
      nBasFrag = 0
      Call Get_iArray('nBas',nBasFrag,nSym)
      nBasFragTot = sum(nBasFrag)
      nBasFragMax = nVec(iFrag) * nBasFragTot
      allocate( linDep(nSym,nBasFragMax,nBasFragMax) )
      allocate( vec(nSym,nBasFragTot,nBasFragTot) )
      allocate( sMO(nSym,nBasFragMax,nBasFragMax) )
      allocate( basis(nSym,nBasFragMax,nBasFragMax) )
      linDep = 0.0
      fragOrbTot = 0
      do iVec = 1, nVec(iFrag)
        nOcc = 0
        vec = 0.0
        call read_vec(iFrag,iVec,vec,nOcc,nBasFrag,nSym)
* add the vectors to the (linear dependent) common MO basis
        do iSym = 1, nSym
          do j = 1, nOcc(isym)
            do k = 1, nBasFrag(iSym)
              linDep(iSym,j+fragOrbTot(iSym),k) = vec(iSym,j,k) 
            end do
          end do
          fragOrbTot(iSym) = fragOrbTot(iSym) + nOcc(iSym)
        end do
      end do
* calculate the overlap matrix of the MOs
      call calculate_sMO(sMO,oneintName,nSym,fragOrbTot,nBasFrag,linDep)
      do iSym = 1, nSym
        lwork = 4*fragOrbTot(iSym)
        allocate( aMatrix(fragOrbTot(iSym),fragOrbTot(iSym)) )
        allocate( aVector(fragOrbTot(iSym)) )
        allocate ( work(lwork) )
        aMatrix = 0.0
        aVector = 0.0
        lDim = 0
        do j = 1, fragOrbTot(iSym)
          do k = 1, fragOrbTot(iSym)
            aMatrix(j,k) = sMO(iSym,j,k)
          end do
        end do
* diagonalize sMO (per symmetry, stored in aMatrix)
        iRc = 1 
        call dsyev('V','L',fragOrbTot(iSym),aMatrix,fragOrbTot(iSym),
     &                                   aVector,work,lwork,iRc)
        aMatrix = transpose(aMatrix)
        if ( iRc .ne. 0 ) then
          write(*,*) 'Something went wrong in dsyev'
          write(*,*) 'iRc = ',iRc
        end if
        call reverse_order(aMatrix,aVector,fragOrbTot(iSym))
        do j = 1, fragOrbTot(iSym)
          if ( abs(aVector(j)) .gt. threshold ) then
            lDim(iSym) = lDim(iSym) + 1
            sdiag(iSym,lDim(iSym)) = aVector(j)
            do k = 1, fragOrbTot(iSym)
              basis(iSym,lDim(iSym),k) = aMatrix(j,k)
            end do
          end if
        end do
        deallocate(work)
        deallocate(aMatrix)
        deallocate(aVector)
      end do
      if ( iFrag .eq. 1 ) then
        write(*,'(2A,8I4)')'Dimension of the common MO basis of ',
     &      'fragment  1:',(lDim(iSym),iSym = 1, nSym)
      else
        write(*,'(44x,I3,A,8I4)') iFrag,':',(lDim(iSym),iSym = 1, nSym)
      end if
* express the common MO basis in the AO basis functions
      basisMotra = 0.0
      do iSym = 1, nSym
        do j = 1, lDim(iSym)                    ! loop over the MOs in the common basis after removing lin. dep.
          do k = 1, FragOrbtot(iSym)
            do l = 1, nBasFrag(iSym)
              basisMotra (iSym,j,l) = basisMotra(iSym,j,l) + 
     &          basis(iSym,j,k)*linDep(iSym,k,l)/sqrt(sdiag(iSym,j))
            end do
          end do
        end do
      end do
*      call calculate_sMO(sMO,oneintName,nSym,lDim,nBasFrag,basisMotra)      ! why should we do this
* express all the states of the frgment in the common basis
      allocate( sU(nBasFragMax,nBasFragMax) )
      allocate( VsU(nBasFragMax,nBasFragMax) )
      allocate( sAO(nSym,nBasFragMax,nBasFragMax) )
      call  getAtomicOverlap(oneintName,sAO,nBasFrag,nSym)
      startVec = 0
      base=project
      do j =1, iFrag - 1
        startVec = startVec + nVec(j)
      end do
      do iVec = 1, nVec(iFrag)
        call getVecFilename(iVec+startVec,vecFilename,base)
        open(36,file=vecFilename)
        vec = 0.0
        call read_vec(iFrag,iVec,vec,nOcc,nBasFrag,nSym)
        do iSym = 1, nSym
          sU = 0.0
          VsU = 0.0
          do j = 1, nOcc(iSym)
            do k = 1, nBasFrag(iSym)
              do l = 1, nBasFrag(iSym)
                sU(j,k) = sU(j,k) + sAO(iSym,k,l) * vec(iSym,j,l)
              end do
            end do
          end do
          do j = 1, lDim(iSym)
            do k = 1, nOcc(iSym)
              do l = 1, nBasFrag(iSym)
                VsU(j,k) = VsU(j,k) + basisMotra(iSym,j,l) * sU(k,l)
              end do
            end do
          end do
          VsU = transpose(VsU)
          write(36,'(I4)') lDim(iSym)
          do j = 1, lDim(iSym)          ! actually, it should be nOcc, but we need to add some null vectors
            write(36,'(4F18.14)')(VsU(j,k),k=1,lDim(iSym))
          end do
        end do
        close(36)
      end do 

*don't forget to deallocate
      deallocate( linDep )
      deallocate( vec )
      deallocate( sMO )
      deallocate( basis )
      deallocate( sAO )
      deallocate( sU )
      deallocate( VsU )
      return
      end subroutine common_basis


      subroutine read_vec(iFrag,iVec,vec,nOcc,nBas,nSym) 
* read the different vector files of fragment iFrag 
* get the number of inactive and active orbitals from the vector file (#INDEX)
      implicit none

      integer                               :: iFrag,iVec,nSym
      integer                               :: iSym,j,k,nMaxBas
      integer, intent(in)                   :: nBas(:)
      integer, intent(out)                  :: nOcc(:)

      real (kind=8), intent(out)            :: vec(:,:,:)

      character (len=6)                     :: mark
      character (len=12)                    :: filename,base
      character (len=132)                   :: line
      character (len = 1 ), allocatable     :: orbLabel(:)

      nMaxBas = maxval(nBas) 
      allocate( orblabel(nMaxBas) )
      base = 'INPORB.'
      call getFileName(iVec,iFrag,filename,base)
      open( 35, file = filename, status = 'old' )

      nOcc = 0
      mark = '#INDEX'
 46   read(35,'(A132)') line
      if (line(1:6).ne.mark) goto 46
      read(35,'(A132)') line
      do iSym = 1, nSym
        orbLabel = ' '
        if ( nBas(iSym) .ne. 0 ) then
          read(35,'(2x,10A)')(orbLabel(j),j=1,nBas(iSym))
        end if
        do j = 1, nBas(iSym)
          if (orbLabel(j) .eq. 'f') nOcc(iSym) = nOcc(iSym) +1     ! frozen in RASSCF
          if (orbLabel(j) .eq. 'i') nOcc(iSym) = nOcc(iSym) +1     ! inactive
          if (orbLabel(j) .eq. '1') nOcc(iSym) = nOcc(iSym) +1     ! active (ras1)
          if (orbLabel(j) .eq. '2') nOcc(iSym) = nOcc(iSym) +1     ! active (ras2)
          if (orbLabel(j) .eq. '3') nOcc(iSym) = nOcc(iSym) +1     ! active (ras3)
        end do
      end do
      rewind(35)
      mark = '#ORB'
 47   read(35,'(A132)') line
      if (line(1:4).ne.mark) goto 47
      do iSym = 1, nSym
        if (nBas(iSym) .ne. 0) then
          do j = 1, nOcc(iSym)
            read(35,'(A132)') line
            read(35,'(5E22.14)') (vec(iSym,j,k),k=1,nBas(iSym))
          end do
        end if
      end do
      close(35)
      deallocate( orblabel )
      return
      end subroutine read_vec



      subroutine calculate_sMO(sMO,oneintName,nSym,nOrb,nBas,MOs)
      implicit none

      integer                         :: nSym,nBasMax
      integer                         :: iSym,j,k,l,m
      integer, intent(in)             :: nOrb(:),nBas(:)

      real (kind=8), intent(out)      :: sMO(:,:,:)
      real (kind=8), intent(in)       :: MOs(:,:,:)
      real (kind=8), allocatable      :: sAO(:,:,:)

      character (len=12)              :: oneintName

      nBasMax = maxval(nBas)
      allocate( sAO(nSym,nMaxBas,nMaxBas) )
      call getAtomicOverlap(oneintname,sAO,nBas,nSym)
      sMO = 0.0
      do iSym = 1, nSym
        do j = 1, nOrb(iSym)
          do k = 1, j
            do l = 1, nBas(iSym)
              do m = 1, nBas(iSym)
                sMO(iSym,j,k) = sMo(iSym,j,k) + MOs(iSym,j,l) 
     &                        * MOs(iSym,k,m) * sAO(iSym,l,m)
              end do
            end do
            sMO(iSym,k,j) = sMo(iSym,j,k)
          end do
        end do
      end do
      deallocate(sAO)
      return
      end subroutine calculate_sMO

      subroutine reverse_order(a,b,length)
      implicit none

      integer                     :: length,i,j,halfway
      real (kind=8), intent(inout):: a(:,:)
      real (kind=8), intent(inout):: b(:)
      real (kind=8)               :: aux

      halfway = int(length/2)
      do i = 1, halfway
        aux = b(i)
        b(i) = b(length+1-i)
        b(length+1-i) = aux
        do j = 1, length
          aux = a(i,j)
          a(i,j) = a(length+1-i,j)
          a(length+1-i,j) = aux  
        end do
      end do
      return
      end subroutine reverse_order

      subroutine readin(threshold,nFragments,nVec,gamess,project)
      implicit none

      integer, parameter           :: nKeys = 4
      integer, parameter           :: maxFrag = 99                ! 99 fragments seems ridiculously large
      integer                      :: nFragments
      integer                      :: jj,iKey,iFrag
      integer, intent(out)         :: nVec(:)  

      real (kind=8)                :: threshold

      character (len=4)            :: key
      character (len=4), dimension(nKeys)  :: keyword
      character (len=20)           :: project
      character (len=132)          :: line

      logical                      :: all_ok = .true.
      logical                      :: gamess 
      logical, dimension(nkeys)    :: hit = .false.

      data keyword /'THRE','FRAG','GAME','PROJ'/

      threshold = 1.0e-6
      nFragments = 1
      gamess = .false.
      project = 'unknown.'
      do while (all_ok)
        read(5,*,iostat=jj) line
        key = adjustl(line)
        call capitalize(key)
        do iKey = 1, nKeys
          if ( key .eq. keyword(iKey) ) hit(iKey) = .true.
        end do
        if (  jj .lt. 0 ) all_ok = .false.
      end do
      do iKey = 1, nKeys
        if ( hit(iKey) ) then
          select case(iKey)
            case(1)
              call locate('THRE')
              read(*,*) threshold
            case(2)
              call locate('FRAG')
              read(*,*) nFragments
              if ( nFragments .gt. maxFrag ) then
                write(*,*) 'Error: number of fragments is too large'
                write(*,'(A,I4)') 'Present value: ',nFragments
                write(*,'(A,I4)') 'Maximum value: ',maxFrag
                write(*,*) 'Change maxFrag and recompile'
                stop
              end if
              read(*,*) (nVec(iFrag), iFrag = 1, nFragments)
            case(3)
              gamess = .true.  
            case(4)
              call locate('PROJ')
              read(*,*) project
              project = trim(project)//'_'
          end select 
        end if
      end do
      return
      end subroutine readin

*  INPUT EXAMPLE, only the first four characters of the keyword are relevant (case insensitive)
*
*  THREshold            ! threshold for linear dependencies in the common MO basis
*     1.0e-3
*  FRAGments            ! number of different fragments, followed by the number of vector sets of each fragment
*     2
*     3   4
*  PROJect              ! root of the vector files that will be generated 
*     tetracene
*

      subroutine capitalize(string)
      implicit none
      integer      :: i
      character(*) string

      do i = 1, len(string)
        if (ichar(string(i:i)).gt.96) then
          string(i:i) = char(ichar(string(i:i))-32)
        endif
      end do
      return
      end subroutine capitalize

      subroutine locate(string)
      implicit none
      character(4)   ::  string,string2
      character(132) ::  line
      rewind(5)
 40   read(5,*) line
      string2=adjustl(line)
      call capitalize(string2)
      if (string2.ne.string) goto 40
      return
      end subroutine locate


      subroutine getAtomicOverlap(filename,sAO,nBas,nSym)
      implicit none

      integer                     :: nBasTot,lTriangle,nSym
      integer                     :: iCounter,iComponent
      integer                     :: iRC,LuOne,iOpt,iSymLbl,iSym,j,k
      integer, intent(in)         :: nBas(:)

      real (kind=8), intent(inout):: sAO(:,:,:)
      real (kind=8), allocatable  :: s(:)

      character (len=12)          :: filename

      nBasTot = sum(nBas)
      lTriangle = ( nBasTot * (nBasTot + 1) ) / 2
      allocate( s(lTriangle) )
      s = 0.0
      sAO = 0.0
      LuOne = 77
      iRc=-1
      iOpt=0
      Call OpnOne(iRC,iOpt,filename,LuOne)
      if (iRC.ne.0) write(6,*)'Something went wrong opening ',filename
      iRC =  0
      iOpt = 2
      iComponent = 1
      iSymLbl = 1
      Call RdOne(iRC,iOpt,'Mltpl  0',iComponent,s,iSymLbl)
      iCounter = 1
      do iSym = 1, nSym
        do j = 1, nBas(iSym)
          do k = 1, j
            sAO(iSym,j,k) = s(iCounter)
            sAO(iSym,k,j) = s(iCounter)
            iCounter = iCounter + 1
          end do
        end do
      end do
      deallocate( s )
      Call ClsOne(iRc,iOpt)
      return
      end subroutine getAtomicOverlap


      subroutine getFilename(iVec,iFrag,filename,base)
      implicit none
      integer              :: iVec,iFrag
      character (len=7)    :: base
      character (len=12)   :: filename

      if ( iFrag .le. 9 .and. iVec .le. 9 ) then
        write(filename,'(A7,I1,A1,I1)') base,iFrag,'_',iVec
      end if
      if ( iFrag .gt. 9 .and. iVec .le. 9 ) then
        write(filename,'(A7,I2,A1,I1)') base,iFrag,'_',iVec
      end if
      if ( iFrag .le. 9 .and. iVec .gt. 9 ) then
        write(filename,'(A7,I1,A1,I2)') base,iFrag,'_',iVec
      end if
      if ( iFrag .gt. 9 .and. iVec .gt. 9 ) then
        write(filename,'(A7,I2,A1,I2)') base,iFrag,'_',iVec
      end if
      return
      end subroutine getFilename

      subroutine getVecFilename(iVec,filename,base)
      implicit none
      integer              :: iVec
      character (len=20)   :: base
      character (len=25)   :: filename

      write(filename,'(A,I0.3,A)') trim(base),iVec,'.vec'
      filename = trim(filename)
      return
      end subroutine getVecFilename


      end program common_MO


* If we want to eliminate (freeze) the core orbitals in the NOCI calculations, they should
* already be eliminated from the construction of the common MO basis. This is in itself not
* so difficult, but we have to find a way to let the program know how many orbitals should be
* skipped when reading the different vector sets: nOCC = (nDouble - nFrozen) + nActive
* nFrozen is not on the RunFile, unless a MOTRA has been done, but this is not what we want.
* The common basis has to be constructed after the CASSCF calculations. There is of course a frozen
* option in CASSCF (core orbitals are not be optimized) but that is not the same. Normal RASSCF should
* be done (full orbital optimization) and then some orbitals should be eliminated from the integral
* transformation. I see two straightforward possibilities:
* 1. Use the standard number of frozen orbitals as internally defined in Molcas
* 2. User input
*
* The first one may eliminate quite some orbitals: for example, the core for TM-3d ions includes 
* the 3s orbital. The second requires some extra work for the user: define the frozen orbitals
* for each fragment. Maybe we can use the following strategy:
*
* 1. When there is no extra input, no orbitals will be eliminated
* 2. Introduces two keywords: 'STANdard core'     --> default core of Molcas
*                             'NONStandard core'  --> followed by explicit definition of the
*                                                     frozen orbitals per fragment 
*

