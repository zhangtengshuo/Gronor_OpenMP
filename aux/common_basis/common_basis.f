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

* ========= E N D   O F   M O D U L E S =========================================

      program common_MO
      use histo_data
      use input_data
      implicit none

      integer      :: nBas,nSym,lDim,totalFrozen
      integer      :: iFrag,nBasFrag,nBasFragMax
      integer      :: startOrb,startBas
      integer      :: i,j,k

      real(kind=8),allocatable  :: superBasis(:,:)                 ! collection of common MO basis sets of the different fragments
      real(kind=8),allocatable  :: commonMOs (:,:)                 ! common MO basis of one fragment
      real(kind=8),allocatable  :: frozenOrbs(:,:)                 ! all frozen orbitals
      real(kind=8),allocatable  :: frzFragOrb(:,:)                 ! average frozen orbitals of one fragment
      real(kind=8),allocatable  :: sDiag     (:)
      real(kind=8),allocatable  :: occNu     (:)
      real                      :: thrs

      character (len=12)         :: oneintName,runfileName

      call readin
      write(*,'(A,E10.2)')'common MO basis with threshold:',threshold
      write(*,*)'R. K. Kathir, C. de Graaf, R. Broer, R. W. A. Havenith'
      write(*,*)'J. Chem. Theory Comput. 16, 2941-2951 (2020)'
      write(*,*)
      write(*,*)' program written by Coen de Graaf, URV (2019)'
      write(*,*)' frozen orbitals by Aitor Sanchez-Mansilla, URV (2020)'
      write(*,*)

      call NameRun('RUNFILE')
      call Get_iScalar('nSym',nSym)
      if ( nSym .ne. 1 ) then
        write(*,*) '  Symmetry is not implemented in GronOR'
        write(*,*) 'Remove the symmetry elements and come back'
        stop
      end if
      Call Get_iArray('nBas',nBas,nSym)              ! get the number of basis functions of the full system
      allocate ( superBasis(nBas,nBas) )
      allocate ( commonMOs (nBas,nBas) )
      allocate ( frozenOrbs(nBas,nBas) )
      allocate ( frzFragOrb(nBas,nBas) )
      allocate ( sDiag     (nBas)      )
      allocate ( occNu     (nBas)      )

      open(12,file = 'COMMONORB')
      write(12,'(A11)') '#INPORB 2.2'
      write(12,'(A5)')  '#INFO'
      write(12,'(A47,E9.3)')
     &  '* Common molecular orbital basis with tau_MO = ',threshold
      write(12,'(3I8)')0,nSym,0
      write(12,'(I8)') nBas
      write(12,'(I8)') nBas
      write(12,'(A4)') '#ORB'

      superBasis = 0.0
      frozenOrbs = 0.0
      occNu = 0.0
      startOrb = 0
      startBas = 0
      totalFrozen = 0
      lDim = 0
      do iFrag = 1, nFragments
        if ( iFrag .le. 9 ) then
          write(runfileName,'(A6,I1)')'RUNFIL',iFrag
          write(oneintName,'(A6,I1)') 'ONEINT',iFrag
        else
          write(runfileName,'(A6,I2)')'RUNFIL',iFrag
          write(oneintName,'(A6,I2)') 'ONEINT',iFrag
        end if
        call NameRun(runfileName)
        call Get_iScalar('nSym',nSym)
        if ( nSym .ne. 1 ) then
          write(*,*) '  Symmetry is not implemented in GronOR'
          write(*,*) 'Remove the symmetry elements and come back'
          stop
        end if
        Call Get_iArray('nBas',nBasFrag,nSym)
        nBasFragMax = nBasFrag * nVec(iFrag)

        commonMOs = 0.0
        sDiag = 0.0
        call common_basis(iFrag,frzFragOrb,commonMOs,sDiag,nBas,
     &                         oneintName,nBasFrag,nBasFragMax,lDim)
        do j = 1, lDim
          do k = 1, nBasFrag
            superBasis(startOrb + j, startBas + k)                     ! dump the common basis in the final set of 
     &                               = commonMOs(j,k)                  ! MOs of the supermolecule
          end do
          occNu(startOrb + j) = sDiag(j)
        end do
        do j = 1, nFrozen(iFrag)
          do k = 1, nBasFrag
             frozenOrbs(totalFrozen + j, startBas + k)
     &                              = frzFragorb(j,k)
          end do
        end do
        totalFrozen = totalFrozen + nFrozen(iFrag)
        startOrb = startOrb + lDim
        startBas = startBas + nBasFrag
      end do
      if ( .not. noAvgCore .and. ( totalFrozen .ne. 0 ) ) then
        call ortho_frozen(frozenOrbs,totalFrozen,nBas)
      end if

      write(*,'(A,I4)') 'Final number of orbitals in MO basis    : ',
     &                       startOrb
      write(*,'(A,I4)') 'Sum of the number of AO basis functions : ',
     &                       startBas
      write(*,'(A,I4)') 'Frozen orbitals in MOTRA                : ',
     &                       totalFrozen
      write(*,'(A,I4)') 'Deleted orbitals in MOTRA               : ',
     &                       startBas - startOrb - totalFrozen
      if ( nBas .ne. startBas ) then
        write(*,*) 'Inconsistency detected'
        write(*,*) 'The number of AO basis functions of the',
     &             '  supermolecule is not equal to the sum'
        write(*,*) 'of the the number of AO basis functions',
     &             ' of the fragments'
        write(*,*) nBas,' not equal to ',startBas
      end if

      if (extra_info) then
        write(*,*)
        write(*,'(2A)') 'Dimension of the common MO ',
     &                      'basis as function of tau_MO'
        write(*,'(A)') ' tau_MO       # basis functions'
        thrs = 0.1
        do i = 1, maxbin - 1
          write (*,'(E10.2,3x,I5)') thrs,freq(i)
          thrs = thrs / 10
        end do
      end if
      
* write the final MO basis of the supermolecule to a file
      do j = 1, totalFrozen
        write(12,'(A9,2I5)')'* ORBITAL',j
        write(12,'(5E22.14)')(frozenOrbs(j,k),k=1,nBas)
      end do
      do j = 1, nBas - totalFrozen
        write(12,'(A9,2I5)')'* ORBITAL',nSym,j + totalFrozen
        write(12,'(5E22.14)')(superBasis(j,k),k=1,nBas)
      end do
      do j = nBas, 1 + totalFrozen, -1                              ! insert occNu = 2 for the frozen orbitals
        occNu(j) = occNu(j - totalFrozen)
      end do
      do j = 1, totalFrozen
        occNu(j) = 2.0
      end do
      write(12,'(A4)') '#OCC'
      write(12,'(A20)') '* OCCUPATION NUMBERS'
      write(12,'(5E22.14)') (occNu(j),j=1,nBas)
      close(12)
* Instead of occupation numbers, we write eigenvalues of the sMO matrix
* For the frozen orbitals we use 2 as occupation number.


      deallocate( superBasis )
      deallocate( commonMOs  )
      deallocate( frozenOrbs )
      deallocate( frzFragOrb )
      deallocate( sDiag      )

      end program common_MO

* ===============================================================================
      
      subroutine common_basis(iFrag,frzFragOrb,commonMOs,sDiag,nBas,
     &                            oneintName,nBasFrag,nBasFragMax,lDim)
      use input_data
      implicit none

      integer,intent(in)           :: iFrag,nBas
      integer,intent(out)          :: lDim
      integer                      :: iVec
      integer                      :: nBasFrag                    ! number of basis functions of the whole system and the fragment
      integer                      :: fragOrbTot                  ! number of linear dependent orbitals in the common MO basis
      integer                      :: lwork,iRc
      integer                      :: startVec,nBasFragMax,nOcc
      integer                      :: j,k,l
      integer                      :: luOne

      real (kind=8),intent(out)    :: commonMOs  (nBas,nBas)
      real (kind=8),intent(out)    :: frzFragOrb (nBas,nBas)
      real (kind=8),intent(out)    :: sDiag      (nBas)

      real (kind=8),allocatable  :: frzVec(:,:)       ! MO coefficients of the frozen orbitals of the electronic states
      real (kind=8),allocatable  :: froVec_1st(:,:)   ! MO coefficients of the frozen orbitals of the first state (NOAV)
      real (kind=8),allocatable  :: frzDensity(:,:)   ! Accumulative density for average frozen orbitals
      real (kind=8),allocatable  :: frzFragAvg(:,:)   ! Average frozen fragment orbitals
      real (kind=8),allocatable  :: linDep(:,:)       ! linear dependent common MO basis
      real (kind=8),allocatable  :: vec(:,:)          ! MO coefficients of the different electronic states of the fragment
      real (kind=8),allocatable  :: sAO(:,:)          ! Atomic basis overlap matrix of the fragment


      real (kind=8),allocatable   :: sMO        (:,:)              ! Overlap matrix of a set of MOs
      real (kind=8),allocatable   :: linDep2    (:,:)              ! linear dependent common MO basis, without the zero vectors
      real (kind=8),allocatable   :: work       (:)
      real (kind=8),allocatable   :: eigenValues(:)
      real (kind=8),allocatable   :: sU         (:,:)              ! Intermediate matrix in the transformation of the vectors to the common MOs basis
      real (kind=8),allocatable   :: VsU        (:,:)              ! The vectors of the different states expressed in the common MO basis
      real (kind=8),allocatable   :: tVsU       (:,:)              ! transpose(VsU)
      real (kind=8),allocatable   :: basis      (:,:)              ! Linear independent basis before transformation to AO basis

      real (kind=8),allocatable       :: commonMO_debug (:,:)      ! only for debugging

      real(kind=8) :: dummy

      character (len=25)              :: vecFileName
      character (len=20)              :: base
      character (len=12)              :: oneintName

      allocate(     linDep(nBasFragMax   ,nBasFrag)       )
      allocate(        vec(nBasFrag      ,nBasFrag)       )
      allocate(        sAO(nBasFrag      ,nBasFrag)       )
      allocate(     frzVec(nFrozen(iFrag),nBasFrag)       )
      allocate( froVec_1st(nFrozen(iFrag),nBasFrag)       )
      allocate( frzDensity(nFrozen(iFrag),nFrozen(iFrag)) )
      allocate( frzFragAvg(nFrozen(iFrag),nBasFrag)       )
      luOne = 87
      call getAtomicOverlap(oneintName,luOne,nBasFrag,sAO)
      fragOrbTot = 0
      linDep     = 0.0
      frzDensity = 0.0
      froVec_1st = 0.0
      frzFragOrb = 0.0
      do iVec = 1, nVec(iFrag)
        vec    = 0.0
        frzVec = 0.0
        call read_vec(iFrag,iVec,nBasFrag,frzVec,vec,nOcc)
* Save the frozen orbitals of the first state, used to express the frozen density of the
* other states. In case of NOAV, these orbitals are directly used in the common basis.
        if ( nFrozen(iFrag) .ne. 0 ) then
          if ( iVec .eq. 1 ) then
            do j = 1, nFrozen(iFrag)
              do k = 1, nBasFrag
                froVec_1st(j,k) = frzVec(j,k)
              end do
            end do
          end if
          if (.not. noAvgCore ) then                             ! accumulate the density of the frozen orbitals
            call calc_frz_density(nVec(iFrag),nBasFrag,nFrozen(iFrag),
     &                                 sAO,froVec_1st,frzVec,frzDensity)
          end if
        end if
* add the other occupied vectors to the (linear dependent) common MO basis
        do j = 1, nOcc
          do k = 1, nBasFrag
            linDep(j+fragOrbTot,k) = vec(j,k)
          end do
        end do
        fragOrbTot = fragOrbTot + nOcc
      end do
      allocate( linDep2(fragOrbTot,nBasFrag) )
      do j = 1, fragOrbTot              ! remove the zero vectors from lindep
        lindep2(j,:) = lindep(j,:)
      end do
*
* After collecting all the vectors of fragment iFrag: generate average set of frozen orbitals
* by diagonalizing the frozen density matrix. Next, the common basis set of MO's is constructed
* using the other occupied (inactive and active) orbitals of all the states of this fragment
*
      if ( nFrozen(iFrag) .ne. 0 ) then
        if ( .not. noAvgCore ) then
          call average_frozen(nFrozen(iFrag),nBasFrag,frzDensity,
     &                                           froVec_1st,frzFragAvg)
          do j = 1, nFrozen(iFrag)
            do k = 1, nBasFrag
              frzFragOrb(j,k) = frzFragAvg(j,k)
            end do
          end do
        else
          if ( iFrag .eq. 1 ) then
            write(*,*) 'No averaging of the core orbitals',
     &               ' --- Orbitals of the first root are used'
            write(*,*)
          end if
          do j = 1, nFrozen(iFrag)
            do k = 1, nBasFrag
              frzFragOrb(j,k) = froVec_1st(j,k)
            end do
          end do
        end if
      end if
* diagonalize sMO
      lwork = 4*fragOrbTot
      allocate ( sMO        ( fragOrbTot,fragOrbTot ) )
      allocate ( work       ( lwork                 ) )
      allocate ( eigenValues( fragOrbTot            ) )
      allocate ( basis      ( fragOrbTot,fragorbTot ) )
      call calculate_sMO(linDep2,sMO,sAO,fragOrbTot,nBasFrag)
      work        = 0.0
      eigenValues = 0.0
      iRc = 1
      call dsyev('V','L',fragOrbTot,sMO,fragOrbTot,
     &                      eigenValues,work,lwork,iRc)
      sMO = transpose(sMO)
      if ( iRc .ne. 0 ) then
        write(*,*) 'Something went wrong in dsyev'
        write(*,*) 'iRc = ',iRc
      end if
      call reverse_order(eigenValues,sMO,fragOrbTot)
      if ( extra_info ) call printDim(eigenValues,fragOrbTot)
      if ( all_epsilons ) call printepsilons(eigenValues,fragOrbTot,
     &                                                         iFrag)
* remove the linear dependencies
      lDim = 0
      do j = 1, fragOrbTot
        if ( abs(eigenValues(j)) .gt. threshold ) then
          lDim = lDim + 1
          sDiag(lDim) = eigenValues(j)
          do k = 1, fragOrbTot
            basis(lDim,k) = sMO(j,k)
          end do
        end if
      end do
      deallocate ( sMO )
      if ( iFrag .eq. 1 ) then
        write(*,'(2A,I4)')'Dimension of the common MO basis of ',
     &      'fragment  1:',lDim
      else
        write(*,'(44x,I3,A,I4)') iFrag,':',lDim
      end if
* express the non-linear dependent basis in the AO basis
      commonMOs = 0.0
      do j = 1, lDim                          ! loop over the MOs in the common basis after removing lin. dep.
        do k = 1, fragOrbTot
          do l = 1, nBasFrag
            commonMOs(j,l) = commonMOs(j,l) +
     &            basis(j,k)*linDep(k,l)/sqrt(sdiag(j))
          end do
        end do
      end do
      if (debug) then
        write(*,*) 'Linear independent MO basis of fragment',iFrag
        allocate( sMO           (lDim,lDim    ) )
        allocate( commonMO_debug(lDim,nBasFrag) )
        do j = 1, lDim
          do k = 1, nBasFrag
            commonMO_debug(j,k) = commonMOs(j,k)
          end do
        end do
        call calculate_sMO(commonMO_debug,sMO,sAO,lDim,nBasFrag)
        deallocate(commonMO_debug)
        deallocate(sMO)
      end if
* express all the states of the fragment in the common basis
      allocate( sU(nBasFragMax,nBasFragMax) )
      allocate( VsU(nBasFragMax,nBasFragMax) )
      allocate(tVsU(nBasFragMax,nBasFragMax) )
      startVec = 0
      base=project
      do j = 1, iFrag - 1
        startVec = startVec + nVec(j)
      end do
      do iVec = 1, nVec(iFrag)
        call getVecFilename(iVec+startVec,vecFilename,base)
        open(36,file=vecFilename)
        vec = 0.0
        call read_vec(iFrag,iVec,nBasFrag,frzVec,vec,nOcc)
        sU = 0.0
        VsU = 0.0
        do j = 1, nOcc
          do k = 1, nBasFrag
            do l = 1, nBasFrag
              sU(j,k) = sU(j,k) + sAO(k,l) * vec(j,l)
            end do
          end do
        end do
        do j = 1, lDim
          do k = 1, nOcc
            do l = 1, nBasFrag
              VsU(j,k) = VsU(j,k) + commonMOs(j,l) * sU(k,l)
            end do
          end do
        end do
        do j = 1, nOcc
          do k = 1, lDim
            tVsU(j,k) = VsU(k,j)
          end do
        end do
*        VsU = transpose(VsU)
        write(36,'(I4)') lDim
        do j = 1, lDim                ! actually, it should be nOcc, but we need to add some null vectors
          write(36,'(4F18.14)')(tVsU(j,k),k=1,lDim)
        end do
        close(36)
      end do

      deallocate( work       )
      deallocate( eigenValues)
      deallocate( basis      )
      deallocate( linDep     )
      deallocate( linDep2    )
      deallocate( sU         )
      deallocate( VsU        )


      end subroutine common_basis


* ===============================================================================

      subroutine readin
      use input_data
      implicit none

      integer, parameter                   :: nKeys = 8
      integer                              :: jj,iKey,iFrag

      character (len=4)                    :: key
      character (len=4), dimension(nKeys)  :: keyword
      character (len=132)                  :: line

      logical                              :: all_ok = .true.
      logical, dimension(nkeys)            :: hit = .false.

      data keyword /'THRE','FRAG','PROJ','EXTR','ALLE','FROZ',
     &              'NOAV','DEBU'/

* Defaults
      threshold    = 1.0e-6
      nFragments   = 1
      nFrozen      = 0
      project      = 'unknown.'
      extra_info   = .false.
      all_epsilons = .false.
      debug        = .false.
      noAvgCore    = .false.

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
              call locate('PROJ')
              read(*,*) project
              project = trim(project)//'_'
            case(4)
              extra_info = .true.
            case(5)
              all_epsilons = .true.
            case(6)
              call locate('FROZ')
              read(*,*) (nFrozen(iFrag), iFrag = 1, nFragments)
            case(7)
              noAvgCore = .true.
            case(8)
              debug = .true.
          end select
        end if
      end do
      return
      end subroutine readin

* -----------------------------------------------------------------------------------------------------------------
* -----------------------------------------------------------------------------------------------------------------
*  INPUT EXAMPLE, only the first four characters of the keyword are
*  relevant (case insensitive)
*
*  THREshold            ! threshold for linear dependencies in the common MO basis
*     1.0e-3
*  FRAGments            ! number of different fragments, followed by the number of vector sets of each fragment
*     2
*     3   4
*  PROJect              ! root of the vector files that will be generated 
*     tetracene
*  FROZen               ! excluding the C-1s orbitals from the calculation
*   18 18
*  EXTRa                ! print information on the size of the common MO basis for different thresholds
*  ALLEpsilons          ! print all the eigenvalues of the sMO of the different fragments
*  DEBUg                ! vast amount of useless info (unless you're debugging)
* -----------------------------------------------------------------------------------------------------------------
* -----------------------------------------------------------------------------------------------------------------


* ===============================================================================

      subroutine read_vec(iFrag,iVec,n,frzVec,vec,nOcc)
      use input_data, only : nFrozen
* read the different vector files of fragment iFrag 
* get the number of inactive and active orbitals from the vector file * (#INDEX)
*
* n       : length of the vectors (number of basis functions)
* frzVec  : frozen vectors (typically the core)
* vec     : the occupied (inactive + active) orbitals used to construct the common MO basis
* nOcc    : number of occupied orbitals (corrected for the number of frozen)
      implicit none

      integer,intent(in)                    :: iFrag,iVec,n
      integer,intent(out)                   :: nOcc
      integer                               :: j,k

      real(kind=8),intent(out)              :: frzVec(nFrozen(iFrag),n)
      real(kind=8),intent(out)              :: vec(n,n)

      character (len=6)                     :: mark
      character (len=12)                    :: filename,base
      character (len=132)                   :: line
      character (len = 1 )                  :: orbLabel(n)

      base = 'INPORB.'
      call getFileName(iVec,iFrag,filename,base)
      open( 35, file = filename, status = 'old' )

      nOcc = 0
      mark = '#INDEX'
 46   read(35,'(A132)') line
      if (line(1:6).ne.mark) goto 46
      read(35,'(A132)') line
      orbLabel = ' '
      read(35,'(2x,10A)')(orbLabel(j),j=1,n)
      do j = 1, n
        if (orbLabel(j) .eq. 'f') nOcc = nOcc + 1     ! frozen in RASSCF
        if (orbLabel(j) .eq. 'i') nOcc = nOcc + 1     ! inactive
        if (orbLabel(j) .eq. '1') nOcc = nOcc + 1     ! active (ras1)
        if (orbLabel(j) .eq. '2') nOcc = nOcc + 1     ! active (ras2)
        if (orbLabel(j) .eq. '3') nOcc = nOcc + 1     ! active (ras3)
      end do
      nOcc = nOcc - nFrozen(iFrag)
      rewind(35)
      mark = '#ORB'
 47   read(35,'(A132)') line
      if (line(1:4).ne.mark) goto 47
      if (nFrozen(iFrag) .ne. 0) then                             ! read the frozen (if any)
        do j = 1, nFrozen(iFrag)
          read(35,'(A132)') line
          read(35,'(5E22.14)') (frzVec(j,k),k=1,n)
        end do
      endif
      do j = 1, nOcc
        read(35,'(A132)') line
        read(35,'(5E22.14)') (vec(j,k),k=1,n)
      end do

      close(35)
      return
      end subroutine read_vec


* ===============================================================================

      subroutine getAtomicOverlap(filename,luOne,n,sAO)
      implicit none

      integer,intent(in)                :: n,luOne
      integer                           :: iCounter,iComponent
      integer                           :: iRC,iOpt,iSymLbl,j,k

      real (kind=8)                     :: s( n * (n + 1 ) / 2 )
      real (kind=8),intent(out)         :: sAO(n,n)

      character (len=12),intent(in)     :: filename

      s = 0.0
      sAO = 0.0
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
      do j = 1, n
        do k = 1, j
          sAO(j,k) = s(iCounter)
          sAO(k,j) = s(iCounter)
          iCounter = iCounter + 1
        end do
      end do
      iOpt = 0
      Call ClsOne(iRc,iOpt)
      return
      end subroutine getAtomicOverlap

* ===============================================================================

      subroutine calculate_sMO(a,sMO,sAO,n,m)
      use input_data, only : debug 
      implicit none

      real(kind=8),intent(in)  :: sAO(m,m)
      real(kind=8),intent(in)  :: a(n,m)
      real(kind=8),intent(out) :: sMO(n,n)
      integer,intent(in)       :: m,n
      integer                  :: i,j,k,l

      sMO = 0.0
      do i = 1, n
        do j = 1, i
          do k = 1, m
            do l = 1, m
              sMO(i,j) = sMO(i,j) + a(i,k) * a(j,l) * sAO(k,l)
            end do
          end do
          sMO(j,i) = sMO(i,j)
        end do
      end do
      if ( debug ) then
        write(*,*) 'MO overlap matrix'
        do i = 1, n
          write(*,'(20E18.8)') sMO(i,:)
        end do
      end if
      return
      end subroutine calculate_sMO

* ===============================================================================

      subroutine reverse_order(a,b,n)
      implicit none

      integer,intent(in)          :: n
      integer                     :: i,j,halfway
      real (kind=8)               :: aux
      real (kind=8),intent(inout) :: a(n),b(n,n)

      halfway = int(n/2)
      do i = 1, halfway
        aux = a(i)
        a(i) = a(n+1-i)
        a(n+1-i) = aux
        do j = 1, n
          aux = b(i,j)
          b(i,j) = b(n+1-i,j)
          b(n+1-i,j) = aux
        end do
      end do
      return
      end subroutine reverse_order

* ===============================================================================

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

* ===============================================================================

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

* ===============================================================================

      subroutine getFilename(iVec,iFrag,filename,base)
      implicit none
      integer,intent(in)              :: iVec,iFrag
      character (len=7),intent(in)    :: base
      character (len=12),intent(out)  :: filename

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

* ===============================================================================

      subroutine getVecFilename(iVec,filename,base)
      implicit none
      integer,intent(in)              :: iVec
      character (len=20),intent(in)   :: base
      character (len=25),intent(out)  :: filename

      write(filename,'(A,I0.3,A)') trim(base),iVec,'.vec'
      filename = trim(filename)
      return
      end subroutine getVecFilename

* ===============================================================================

      subroutine printDim(a,n)
      use histo_data
      implicit none

      integer,intent(in)       :: n
      integer                  :: i,bin
      real (kind=8),intent(in) :: a(n)
      real (kind=8)            :: thrs

      thrs = 0.1
      bin = 1
      do i = 1, n
        if ( a(i) .lt. thrs ) then
          freq(bin) = freq(bin) + (i-1)
          thrs = thrs / 10
          bin = bin + 1
        end if
      end do
      if (bin .gt. maxbin) maxbin = bin
      return
      end subroutine printDim

* ===============================================================================

      subroutine printepsilons(a,n,i)
      implicit none

      integer,intent(in)       :: n
      integer                  :: i,j
      real (kind=8),intent(in) :: a(n)

      write(*,*)
      write(*,'(2A,I4)')'List of all eigenvalues of the MO overlap',
     &           ' matrix of fragment ',i
      do j = 1, n
        write(*,'(I5,E20.8)') j,a(j)
      end do
      return
      end subroutine printepsilons

* ===============================================================================

      subroutine ortho_frozen(frozenOrbs,totalFrozen,nBas)
      use input_data, only : debug
      implicit none

      integer,intent(in)         :: nBas,totalFrozen
      integer                    :: j,k,l
      integer                    :: iRc,lwork
      integer                    :: luOne

      real(kind=8),intent(inout) :: frozenOrbs (nBas       ,nBas       )
      real(kind=8)               :: sAO        (nBas       ,nBas       )
      real(kind=8)               :: sMO        (totalFrozen,totalFrozen)
      real(kind=8)               :: aMatrix    (totalFrozen,nBas       )
      real(kind=8)               :: eigenValues(totalFrozen            )
      real(kind=8)               :: work       (4 * totalFrozen        )
      real(kind=8)               :: orbs_debug (totalFrozen,nBas       )

      character (len=12)         :: oneintName

      lwork = 4 * totalFrozen

      do j = 1, totalFrozen
        do k = 1, nBas
          aMatrix(j,k) = frozenOrbs(j,k)  ! save the non-orthogonal frozen orbitals for later use
        end do
      end do
      frozenOrbs = 0.0

      call NameRun('RUNFILE')
      write(oneintName,'(A6)') 'ONEINT'
      luOne = 88
      call getAtomicOverlap(oneintName,luOne,nBas,sAO)
      call calculate_sMO(aMatrix,sMO,sAO,totalFrozen,nBas)
      iRc = 1
      call dsyev('V','L',totalFrozen,sMO,totalFrozen,
     &   eigenValues,work,lwork,iRc)
      if ( iRc .ne. 0 ) then
        write(*,*) 'Something went wrong in dsyev (sMO frozen)'
        write(*,*) 'iRc = ',iRc
      end if
      sMO = transpose(sMO)
      if ( debug ) then
        write (*,*) 'eigenvalues and eigenvectors of the average',
     &              ' frozen SMO matrix'
        do j = 1, totalFrozen
          write(*,'(F15.8,20E18.8)')eigenValues(j),sMO(j,:)
        end do
      end if

*     Lowdin orthonormalization (expressing the eigenvectors of sMO in AO basis)
      do j = 1, totalFrozen
        do k = 1, totalFrozen
          do l = 1, nBas
            frozenOrbs(j,l) = frozenOrbs(j,l) +
     &               sMO(j,k)*aMatrix(k,l) / sqrt(eigenValues(j))
          enddo
        enddo
      enddo
      if (debug) then
        sMO = 0.0
        write (*,*) 'Lowdin orthgonalized average frozen orbitals'
        do j = 1, totalFrozen
          write(*,*) 'Average frozen orbital ',j
          write(*,'(5E22.14)') frozenOrbs(j,:)
          orbs_debug(j,:) = frozenOrbs(j,:)
        end do
        call calculate_sMO(orbs_debug,sMO,sAO,totalFrozen,nBas)
      end if

      return
      end subroutine ortho_frozen

* ===============================================================================
 
      subroutine calc_frz_density(n,nB,nF,sAO,froVec_1st,
     &                                             frzVec,frzDensity)
      use input_data, only : debug
      implicit none

      integer,intent(in)         :: n,nB,nF
      integer                    :: j
      real(kind=8),intent(in)    :: sAO(nB,nB)
      real(kind=8),intent(in)    :: frzVec(nF,nB),froVec_1st(nF,nB)
      real(kind=8),intent(inout) :: frzDensity(nF,nF)
      real(kind=8)               :: mprod(nB,nF)

      mprod = matmul(sAO,transpose(frovec_1st))
      frzDensity = frzDensity + matmul(frzVec,mprod) / n
      if ( debug ) then
        write(*,'(2A)') '*  Accumulating the density matrices of ',
     &      'the frozen orbitals'
        do j = 1,nF
          write(*,'(20E18.8)') frzDensity(j,:)
        end do
      end if

      return
      end subroutine calc_frz_density
 
* ===============================================================================

      subroutine average_frozen(nF,nB,frzDensity,froVec_1st,frzFragAvg)
      use input_data, only : debug
      implicit none

      integer,intent(in)    :: nF,nB
      integer               :: j
      integer               :: iRc

      real(kind=8),intent(inout)   :: frzDensity(nF,nF)
      real(kind=8),intent(in)      :: froVec_1st(nF,nB)
      real(kind=8),intent(out)     :: frzFragAvg(nF,nB)
      real(kind=8)                 :: work(4*nF)
      real(kind=8)                 :: eigenValues(nF)

      iRc = 1
      work = 0.0
      call dsyev('V','L',nF,frzDensity,nF,eigenValues,work,4*nF,iRc)
      if ( iRc .ne. 0 ) then
        write(*,*) 'Something went wrong in dsyev (frzDensity)'
        write(*,*) 'iRc = ',iRc
      end if
      frzDensity = transpose(frzDensity)            ! eigenvectors as rows
      if (debug ) then
        write(*,*)'eigenvectors of SMO(frozen)'
        do j = 1, nF
          write(*,'(A,F18.10)') 'eigenvalue :',eigenValues(j)
          write(*,'(20E18.8)') frzDensity(j,:)
        end do
      end if
      frzFragAvg = matmul(frzDensity,froVec_1st)
      if ( debug ) then
        write(*,*)'average frozen orbitals' 
        do j = 1, nF
          write(*,'(5E22.14)') frzFragAvg(j,:)
        end do
      end if

      return
      end subroutine average_frozen

* ===============================================================================

