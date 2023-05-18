module gcommon_histo_data
  implicit none
  integer   :: maxbin
  integer,dimension(30) :: freq
end module gcommon_histo_data


module gcommon_fragment_data
  implicit none
  integer,allocatable     :: nElectrons(:,:)
  real(kind=8),allocatable:: ener(:,:),enerpt(:,:)
end module gcommon_fragment_data

module gcommon_input_data
  implicit none
  integer,parameter :: maxFrag=99
  integer     :: nFragments
  integer,dimension(maxFrag)    :: nFrozen,nVec
  real(kind=8):: threshold
  logical     :: extra_info,all_epsilons
  logical     :: debug,noAvgCore,AOIntegrals
  logical     :: fragLabels,energy_on_INPORB
  character (len=4),allocatable :: fragName(:),fragState(:)
  character (len=255),allocatable :: fragLabel(:)
  character (len=32)      :: mebfLabel
  character (len=255)     :: project
end module gcommon_input_data


program gcommon
  use gcommon_histo_data
  use gcommon_input_data
  implicit none

  external :: gcommon_readin,gcommon_add_detinfo,gcommon_ortho_frozen
  external :: gcommon_generate_vecfiles,gcommon_common_basis
  external :: NameRun,Get_iScalar,Get_iArray
  
  integer:: nBas,nSym,lDim,totalFrozen
  integer:: iFrag,nBasFrag,nBasFragMax
  integer:: startOrb,startBas
  integer:: i,j,k

  real(kind=8),allocatable  :: superBasis(:,:)     ! collection of common MO basis sets of the different fragments
  real(kind=8),allocatable  :: commonMOs(:,:)     ! common MO basis of one fragment
  real(kind=8),allocatable  :: frozenOrbs(:,:)     ! all frozen orbitals
  real(kind=8),allocatable  :: frzFragOrb(:,:)     ! average frozen orbitals of one fragment
  real(kind=8),allocatable  :: sDiag    (:)
  real(kind=8),allocatable  :: occNu    (:)
  real(kind=8)  :: thrs

  character (len=255)  :: oneintName,runfileName

  call gcommon_readin
  write(*,'(A,E10.2)')'common MO basis with threshold:',threshold
  write(*,*)'R. K. Kathir, C. de Graaf, R. Broer, R. W. A. Havenith'
  write(*,*)'J. Chem. Theory Comput. 16, 2941-2951 (2020)'
  write(*,*)
  write(*,*)' program written by Coen de Graaf, URV (2019)'
  write(*,*)' frozen orbitals by Aitor Sanchez-Mansilla, URV (2020)'
  write(*,*)

  call NameRun('RUNFIL')
  call Get_iScalar('nSym',nSym)
  if(nSym.ne.1) then
    write(*,*) '  Symmetry is not implemented in GronOR'
    write(*,*) 'Remove the symmetry elements and come back'
    stop
  endif
  Call Get_iArray('nBas',nBas,nSym)  ! get the number of basis functions of the full system
  allocate(superBasis(nBas,nBas))
  allocate(commonMOs(nBas,nBas))
  allocate(frozenOrbs(nBas,nBas))
  allocate(frzFragOrb(nBas,nBas))
  allocate(sDiag    (nBas))
  allocate(occNu    (nBas))

  open(12,file='COMMONORB')
  write(12,'(A11)') '#INPORB 2.2'
  write(12,'(A5)')  '#INFO'
  write(12,'(A47,E9.3)') '* Common molecular orbital basis with tau_MO = ',threshold
  write(12,'(3I8)')0,nSym,0
  write(12,'(I8)') nBas
  write(12,'(I8)') nBas
  write(12,'(A4)') '#ORB'

  superBasis=0.0
  frozenOrbs=0.0
  occNu=0.0
  startOrb=0
  startBas=0
  totalFrozen=0
  lDim=0
  if(AOIntegrals) then
    write(*,*) 'A set of unit vectors will be written to COMMONORB'
    write(*,*)
    do i=1,nBas
      superBasis(i,i)=1.0
    enddo
    startBas=nBas
    startOrb=nBas
  endif
  do iFrag=1,nFragments
    if(iFrag .le. 9) then
      write(runfileName,'(a6,a)')'RUNFIL',trim(fragName(iFrag))
      write(oneintName,'(a6,a)') 'ONEINT',trim(fragName(iFrag))
    else
      write(runfileName,'(a6,a)')'RUNFIL',trim(fragName(iFrag))
      write(oneintName,'(a6,a)') 'ONEINT',trim(fragName(iFrag))
    endif
    call NameRun(runfileName)
    call Get_iScalar('nSym',nSym)
    if(nSym.ne.1) then
      write(*,*) '  Symmetry is not implemented in GronOR'
      write(*,*) 'Remove the symmetry elements and come back'
      stop
    endif
    Call Get_iArray('nBas',nBasFrag,nSym)
    nBasFragMax=nBasFrag*nVec(iFrag)

    commonMOs=0.0
    sDiag=0.0
    if(.not. AOIntegrals) then
      call gcommon_common_basis(iFrag,frzFragOrb,commonMOs,sDiag,nBas, &
         oneintName,nBasFrag,nBasFragMax,lDim)
      do j=1,lDim
        do k=1,nBasFrag
          ! dump the common basis in the final set of MOs of the supermolecule
          superBasis(startOrb+j,startBas+k)=commonMOs(j,k) 
        enddo
        occNu(startOrb+j)=sDiag(j)
      enddo
      do j=1,nFrozen(iFrag)
        do k=1,nBasFrag
          frozenOrbs(totalFrozen+j,startBas+k)=frzFragorb(j,k)
        enddo
      enddo
      totalFrozen=totalFrozen+nFrozen(iFrag)
      startOrb=startOrb+lDim
      startBas=startBas+nBasFrag
    else
      call gcommon_generate_vecfiles(iFrag,nBasfrag)
    endif
  enddo
  if(.not. noAvgCore .and.(totalFrozen.ne.0)) then
    call gcommon_ortho_frozen(frozenOrbs,totalFrozen,nBas)
  endif

! Adding number of dets, fragment labels and changing the number of inactives (if frozen.ne.0)
  call gcommon_add_detinfo()
  
  write(*,'(A,I4)') 'Final number of orbitals in MO basis    : ',startOrb
  write(*,'(A,I4)') 'Sum of the number of AO basis functions : ',startBas
  write(*,'(A,I4)') 'Frozen orbitals in MOTRA    : ',totalFrozen
  write(*,'(A,I4)') 'Deleted orbitals in MOTRA   : ',startBas-startOrb-totalFrozen
  if(nBas.ne.startBas) then
    write(*,*) 'Inconsistency detected'
    write(*,*) 'The number of AO basis functions of the supermolecule is not equal to the sum'
    write(*,*) 'of the the number of AO basis functions of the fragments'
    write(*,*) nBas,' not equal to ',startBas
  endif

  if(extra_info) then
    write(*,*)
    write(*,'(a)') 'Dimension of the common MO basis as function of tau_MO'
    write(*,'(a)') ' tau_MO # basis functions'
    thrs=0.1
    do i=1,maxbin-1
      write (*,'(E10.2,3x,I5)') thrs,freq(i)
      thrs=thrs/10
    enddo
  endif
  
! write the final MO basis of the supermolecule to a file
  
  do j=1,totalFrozen
    write(12,'(A9,2I5)')'* ORBITAL',j
    write(12,'(5E22.14)')(frozenOrbs(j,k),k=1,nBas)
  enddo
  do j=1,nBas-totalFrozen
    write(12,'(A9,2I5)')'* ORBITAL',nSym,j+totalFrozen
    write(12,'(5E22.14)')(superBasis(j,k),k=1,nBas)
  enddo
  do j=nBas,1+totalFrozen,-1! insert occNu = 2 for the frozen orbitals
    occNu(j)=occNu(j-totalFrozen)
  enddo
  do j=1,totalFrozen
    occNu(j)=2.0
  enddo
  write(12,'(A4)') '#OCC'
  write(12,'(A20)') '* OCCUPATION NUMBERS'
  write(12,'(5E22.14)') (occNu(j),j=1,nBas)
  close(12)
  
! Instead of occupation numbers, we write eigenvalues of the sMO matrix
! For the frozen orbitals we use 2 as occupation number.

  deallocate(superBasis)
  deallocate(commonMOs)
  deallocate(frozenOrbs)
  deallocate(frzFragOrb)
  deallocate(sDiag)

  if(.not.fragLabels) then
    write(*,*)
    write(*,'(a)')  'Note that fragment labels were not provided.'
    write(*,'(a)') 'Fragment labels are optional, but strongly recommended.'
  endif

end program gcommon


subroutine gcommon_common_basis &
    (iFrag,frzFragOrb,commonMOs,sDiag,nBas,oneintName,nBasFrag,nBasFragMax,lDim)
  use gcommon_input_data
  implicit none

  external :: dsyev
  external :: gcommon_getAtomicOverlap,gcommon_getFilename
  external :: gcommon_read_vec,gcommon_calculate_sMO
  external :: NameRun
  external :: gcommon_printepsilons,gcommon_printDim,gcommon_reverse_order
  external :: gcommon_calc_frz_density,gcommon_average_frozen
  
  integer,intent(in)     :: iFrag,nBas
  integer,intent(out)    :: lDim
  integer    :: iVec
  integer    :: nBasFrag  ! number of basis functions of the whole system and the fragment
  integer    :: fragOrbTot! number of linear dependent orbitals in the common MO basis
  integer    :: lwork,iRc
  integer    :: startVec,nBasFragMax,nOcc
  integer    :: j,k,l
  integer    :: luOne

  real (kind=8),intent(out)    :: commonMOs  (nBas,nBas)
  real (kind=8),intent(out)    :: frzFragOrb (nBas,nBas)
  real (kind=8),intent(out)    :: sDiag(nBas)

  real (kind=8),allocatable  :: frzVec(:,:) ! MO coefficients of the frozen orbitals of the electronic states
  real (kind=8),allocatable  :: froVec_1st(:,:)   ! MO coefficients of the frozen orbitals of the first state (NOAV)
  real (kind=8),allocatable  :: frzDensity(:,:)   ! Accumulative density for average frozen orbitals
  real (kind=8),allocatable  :: frzFragAvg(:,:)   ! Average frozen fragment orbitals
  real (kind=8),allocatable  :: linDep(:,:) ! linear dependent common MO basis
  real (kind=8),allocatable  :: vec(:,:)    ! MO coefficients of the different electronic states of the fragment
  real (kind=8),allocatable  :: sAO(:,:)    ! Atomic basis overlap matrix of the fragment


  real (kind=8),allocatable   :: sMO  (:,:)  ! Overlap matrix of a set of MOs
  real (kind=8),allocatable   :: linDep2    (:,:)  ! linear dependent common MO basis, without the zero vectors
  real (kind=8),allocatable   :: work (:)
  real (kind=8),allocatable   :: eigenValues(:)
  real (kind=8),allocatable   :: sU   (:,:)  ! Intermediate matrix in the transformation of the vectors to the common MOs basis
  real (kind=8),allocatable   :: VsU  (:,:)  ! The vectors of the different states expressed in the common MO basis
  real (kind=8),allocatable   :: tVsU (:,:)  ! transpose(VsU)
  real (kind=8),allocatable   :: basis(:,:)  ! Linear independent basis before transformation to AO basis

  real (kind=8),allocatable :: commonMO_debug (:,:)! only for debugging

  character (len=3)   :: suffix
  character (len=255) :: vecFileName
  character (len=255) :: base
  character (len=12)  :: oneintName

  allocate(linDep(nBasFragMax,nBasFrag))
  allocate(vec(nBasFrag,nBasFrag))
  allocate(sAO(nBasFrag,nBasFrag))
  allocate(frzVec(nFrozen(iFrag),nBasFrag))
  allocate(froVec_1st(nFrozen(iFrag),nBasFrag))
  allocate(frzDensity(nFrozen(iFrag),nFrozen(iFrag)))
  allocate(frzFragAvg(nFrozen(iFrag),nBasFrag))
  luOne=87
  call gcommon_getAtomicOverlap(oneintName,luOne,nBasFrag,sAO)
  fragOrbTot=0
  linDep    =0.0
  frzDensity=0.0
  froVec_1st=0.0
  frzFragOrb=0.0
  do iVec=1,nVec(iFrag)
    vec   =0.0
    frzVec=0.0
    call gcommon_read_vec(iFrag,iVec,nBasFrag,frzVec,vec,nOcc)
!  Save the frozen orbitals of the first state, used to express the frozen density of the
!  other states. In case of NOAV, these orbitals are directly used in the common basis.
    if(nFrozen(iFrag).ne.0) then
      if(iVec.eq.1) then
        do j=1,nFrozen(iFrag)
          do k=1,nBasFrag
            froVec_1st(j,k)=frzVec(j,k)
          enddo
        enddo
      endif
      if(.not. noAvgCore) then     ! accumulate the density of the frozen orbitals
        call gcommon_calc_frz_density(nVec(iFrag),nBasFrag,nFrozen(iFrag),sAO,froVec_1st,frzVec,frzDensity)
      endif
    endif
! add the other occupied vectors to the (linear dependent) common MO basis
    do j=1,nOcc
      do k=1,nBasFrag
        linDep(j+fragOrbTot,k)=vec(j,k)
      enddo
    enddo
    fragOrbTot=fragOrbTot+nOcc
  enddo
  allocate(linDep2(fragOrbTot,nBasFrag))
  do j=1,fragOrbTot  ! remove the zero vectors from lindep
    lindep2(j,:)=lindep(j,:)
  enddo
  
!  After collecting all the vectors of fragment iFrag: generate average set of frozen orbitals  
!  by diagonalizing the frozen density matrix. Next, the common basis set of MO's is constructed
!  using the other occupied (inactive and active) orbitals of all the states of this fragment

  if(nFrozen(iFrag).ne.0) then
    if(.not. noAvgCore) then
      call gcommon_average_frozen(nFrozen(iFrag),nBasFrag,frzDensity,froVec_1st,frzFragAvg)
      do j=1,nFrozen(iFrag)
        do k=1,nBasFrag
          frzFragOrb(j,k)=frzFragAvg(j,k)
        enddo
      enddo
    else
      if(iFrag.eq.1) then
        write(*,*) 'No averaging of the core orbitals --- Orbitals of the first root are used'
        write(*,*)
      endif
      do j=1,nFrozen(iFrag)
        do k=1,nBasFrag
          frzFragOrb(j,k)=froVec_1st(j,k)
        enddo
      enddo
    endif
  endif
 ! diagonalize sMO
  lwork=4*fragOrbTot
  allocate(sMO (fragOrbTot,fragOrbTot))
  allocate(work(lwork))
  allocate(eigenValues(fragOrbTot))
  allocate(basis(fragOrbTot,fragorbTot))
  call gcommon_calculate_sMO(linDep2,sMO,sAO,fragOrbTot,nBasFrag)
  work =0.0
  eigenValues=0.0
  iRc=1
  call dsyev('V','L',fragOrbTot,sMO,fragOrbTot,eigenValues,work,lwork,iRc)
  sMO=transpose(sMO)
  if(iRc.ne.0) then
    write(*,*) 'Something went wrong in dsyev'
    write(*,*) 'iRc=',iRc
  endif
  call gcommon_reverse_order(eigenValues,sMO,fragOrbTot)
  if(extra_info) call gcommon_printDim(eigenValues,fragOrbTot)
  if(all_epsilons) call gcommon_printepsilons(eigenValues,fragOrbTot,iFrag)
! remove the linear dependencies
  lDim=0
  do j=1,fragOrbTot
    if(abs(eigenValues(j)).gt.threshold) then
      lDim=lDim+1
      sDiag(lDim)=eigenValues(j)
      do k=1,fragOrbTot
        basis(lDim,k)=sMO(j,k)
      enddo
    endif
  enddo
  deallocate(sMO)
  if(iFrag.eq.1) then
    write(*,'(a,I4)')'Dimension of the common MO basis of fragment  1:',lDim
  else
    write(*,'(44x,I3,A,I4)') iFrag,':',lDim
  endif
! express the non-linear dependent basis in the AO basis
  commonMOs=0.0
  do j=1,lDim  ! loop over the MOs in the common basis after removing lin. dep.
    do k=1,fragOrbTot
      do l=1,nBasFrag
        commonMOs(j,l)=commonMOs(j,l)+basis(j,k)*linDep(k,l)/sqrt(sdiag(j))
      enddo
    enddo
  enddo
  if(debug) then
    write(*,*) 'Linear independent MO basis of fragment',iFrag
    allocate(sMO     (lDim,lDim))
    allocate(commonMO_debug(lDim,nBasFrag))
    do j=1,lDim
      do k=1,nBasFrag
        commonMO_debug(j,k)=commonMOs(j,k)
      enddo
    enddo
    call gcommon_calculate_sMO(commonMO_debug,sMO,sAO,lDim,nBasFrag)
    deallocate(commonMO_debug)
    deallocate(sMO)
  endif
! express all the states of the fragment in the common basis
  allocate(sU(nBasFragMax,nBasFragMax))
  allocate(VsU(nBasFragMax,nBasFragMax))
  allocate(tVsU(nBasFragMax,nBasFragMax))
  startVec=0
  base=project
  do j=1,iFrag-1
    startVec=startVec+nVec(j)
  enddo
  suffix='vec'
  do iVec=1,nVec(iFrag)
    call gcommon_getFilename(iVec+startVec,vecFilename,base,suffix)
    open(36,file=vecFilename)
    vec=0.0
    call gcommon_read_vec(iFrag,iVec,nBasFrag,frzVec,vec,nOcc)
    sU=0.0
    VsU=0.0
    tVsU=0.0
    do j=1,nOcc
      do k=1,nBasFrag
        do l=1,nBasFrag
          sU(j,k)=sU(j,k)+sAO(k,l)*vec(j,l)
        enddo
      enddo
    enddo
    do j=1,lDim
      do k=1,nOcc
        do l=1,nBasFrag
          VsU(j,k)=VsU(j,k)+commonMOs(j,l)*sU(k,l)
        enddo
      enddo
    enddo
    do j=1,nOcc
      do k=1,lDim
        tVsU(j,k)=VsU(k,j)
      enddo
    enddo
 !  VsU=transpose(VsU)
    write(36,'(I4)') lDim
    do j=1,lDim    ! actually, it should be nOcc, but we need to add some null vectors
      write(36,'(4F18.14)')(tVsU(j,k),k=1,lDim)
    enddo
    close(36)
  enddo
  deallocate(work)
  deallocate(eigenValues)
  deallocate(basis)
  deallocate(linDep)
  deallocate(linDep2)
  deallocate(sU)
  deallocate(VsU)

  return
end subroutine gcommon_common_basis


subroutine gcommon_readin
  use gcommon_input_data
  use gcommon_fragment_data
  implicit none

  external :: gcommon_capitalize,gcommon_locate
  
  integer, parameter :: nKeys=10
  integer:: j,jj,iKey,iFrag
  integer:: start

  character (len=4)  :: key
  character (len=4), dimension(nKeys)  :: keyword
  character (len=132):: line,string

  logical:: all_ok=.true.
  logical:: fragNameOK
  logical, dimension(nkeys):: hit=.false.

  data keyword /'THRE','FRAG','PROJ','EXTR','ALLE','FROZ','NOAV','DEBU','LABE','AOIN'/

! Defaults
  threshold   =1.0e-6
  nFragments  =1
  nFrozen= 0
  project= 'unknown.'
  extra_info  =.false.
  all_epsilons=.false.
  debug =.false.
  noAvgCore   =.false.
  fragLabels  =.false.
  AOIntegrals =.false.

  do while(all_ok)
    read(5,*,iostat=jj) line
    key=adjustl(line)
    call gcommon_capitalize(key)
    do iKey=1,nKeys
      if(key.eq.keyword(iKey)) hit(iKey)=.true.
    enddo
    if(jj.lt.0) all_ok=.false.
  enddo

  do iKey=1,nKeys
    if(hit(iKey)) then
      select case(iKey)
      case(1)
        call gcommon_locate('THRE')
        read(*,*) threshold
      case(2)
        call gcommon_locate('FRAG')
        read(*,*) nFragments
        if(nFragments.gt.maxFrag) then
          write(*,*) 'Error: number of fragments is too large'
          write(*,'(A,I4)') 'Present value: ',nFragments
          write(*,'(A,I4)') 'Maximum value: ',maxFrag
          write(*,*) 'Change maxFrag and recompile'
          stop
        endif
        read(*,*) (nVec(iFrag),iFrag=1,nFragments)
        allocate(nElectrons(nFragments,maxval(nVec)))
        allocate(ener(nFragments,maxval(nVec)))
        allocate(enerpt(nFragments,maxval(nVec)))
        ener=0.0
        enerpt=0.0
      case(3)
        call gcommon_locate('PROJ')
        read(*,*) project
 ! project=trim(project)//'_'
      case(4)
        extra_info=.true.
      case(5)
        all_epsilons=.true.
      case(6)
        call gcommon_locate('FROZ')
        read(*,*) (nFrozen(iFrag),iFrag=1,nFragments)
      case(7)
        noAvgCore=.true.
      case(8)
        debug=.true.
      case(9)
        call gcommon_locate('LABE')
        fragLabels=.true.
        allocate(fragName(nFragments))
        allocate(fragLabel(sum(nVec)))
        allocate(fragState(sum(nVec)))
        fragName(:) =''
        fragLabel(:)=''
        fragState(:)=''
        start=0
        mebfLabel=''
        do j=1,nFragments
          start=start+1
          read(*,'(A)') line
          string=adjustl(line)
          fragNameOK=.false.
          do jj=1,len(trim(string))
            if(string(jj:jj).ne.' ') then
              if(.not.fragNameOK) then
                write(fragName(j),'(a,a)') trim(adjustl(fragName(j))),string(jj:jj)
              else
                write(fragState(start),'(a,a)') trim(adjustl(fragState(start))),string(jj:jj)
              endif
            else
              if(fragNameOK) then
                if(string(jj-1:jj-1).ne.' ') start=start+1
              else
                fragnameOK=.true.
              endif
            endif
          enddo
          write(mebfLabel,'(a,a)') trim(mebfLabel),trim(fragName(j))
        enddo
        start=1
        do j=1,nFragments
          do jj=start,start+nVec(j)-1
            write(fragLabel(jj),'(a,a,a)') trim(fragName(j)),'_',trim(fragState(jj))
          enddo
          start=start+nVec(j)
        enddo
      case(10)
        AOIntegrals=.true.
      end select
    endif
  enddo
  return
end subroutine gcommon_readin


!  INPUT EXAMPLE, only the first four characters of the keyword are
!  relevant (case insensitive)
!
!  THREshold! threshold for linear dependencies in the common MO basis
!     1.0e-3
!  FRAGments! number of different fragments, followed by the number of vector sets of each fragment
!     2
!     3   4
!  PROJect  ! root of the vector files that will be generated 
!     tetracene
!  FROZen   ! excluding the C-1s orbitals from the calculation
!   18 18
!  Label    ! Labels of the fragment states, used for the definition of the MEBFs
!   A  S0 S1 T1   !    labels of fragment A
!   B  S0 S1 T1 T2!    labels of fragment B
!  EXTRa    ! print information on the size of the common MO basis for different thresholds
!  ALLEpsilons    ! print all the eigenvalues of the sMO of the different fragments
!  DEBUg    ! vast amount of useless info (unless you're debugging)


subroutine gcommon_read_vec(iFrag,iVec,n,frzVec,vec,nOcc)
  use gcommon_input_data, only : nFrozen,energy_on_INPORB,project,nVec
  use gcommon_fragment_data
! read the different vector files of fragment iFrag 
! get the number of inactive and active orbitals from the vector file*(#INDEX)

! n : length of the vectors (number of basis functions)
! frzVec  : frozen vectors (typically the core)
! vec     : the occupied (inactive+active) orbitals used to construct the common MO basis
! nOcc    : number of occupied orbitals (corrected for the number of frozen)
  
  implicit none

  external :: gcommon_getFilename
  
  integer,intent(in)  :: iFrag,iVec,n
  integer,intent(out) :: nOcc
  integer :: j,k,startVec

  real(kind=8),intent(out)  :: frzVec(nFrozen(iFrag),n)
  real(kind=8),intent(out)  :: vec(n,n)
  real(kind=8)  :: occ(n)

  character (len=6)   :: mark
  character (len=20)  :: base
  character (len=255) :: orbFilename
  character (len=132) :: line
  character (len=3)   :: suffix
  character (len=1):: orbLabel(n)

  base=project
  suffix='orb'
  startVec=0
  do j=1,iFrag-1
    startVec=startVec+nVec(j)
  enddo
  call gcommon_getFilename(iVec+startVec,orbFilename,project,suffix)
  open(35,file=orbFilename,status='old')

  nOcc=0
  mark='#INDEX'
46 read(35,'(A132)') line
  if(line(1:6).ne.mark) goto 46
  read(35,'(A132)') line
  orbLabel=' '
  read(35,'(2x,10A)')(orbLabel(j),j=1,n)
  do j=1,n
    if(orbLabel(j).eq.'f') nOcc=nOcc+1     ! frozen in RASSCF
    if(orbLabel(j).eq.'i') nOcc=nOcc+1     ! inactive
    if(orbLabel(j).eq.'1') nOcc=nOcc+1     ! active (ras1)
    if(orbLabel(j).eq.'2') nOcc=nOcc+1     ! active (ras2)
    if(orbLabel(j).eq.'3') nOcc=nOcc+1     ! active (ras3)
  enddo
  nOcc=nOcc-nFrozen(iFrag)
  rewind(35)
  mark='#ORB'
47 read(35,'(A132)') line
  if(line(1:4).ne.mark) goto 47
  if(nFrozen(iFrag).ne.0) then     ! read the frozen (if any)
    do j=1,nFrozen(iFrag)
      read(35,'(A132)') line
      read(35,'(5E22.14)') (frzVec(j,k),k=1,n)
    enddo
  endif
  do j=1,nOcc
    read(35,'(A132)') line
    read(35,'(5E22.14)') (vec(j,k),k=1,n)
  enddo
  rewind(35)
  energy_on_INPORB=.true.
  mark='#INFO'
48 read(35,'(A132)') line
  if(line(1:5).ne.mark) goto 48
  read(35,'(A132)') line
  if(line(10:34).ne.'natural orbitals for root') then
    energy_on_INPORB=.false.
    print*,"FALSE"
  else
    read(line,'(48x,f22.12)') ener(iFrag,ivec)
    write(*,'(2f22.12)') ener(iFrag,ivec)
  endif
  rewind(35)
  mark='#OCC'
49 read(35,'(A132)') line
  if(line(1:4).ne.mark) goto 49
  read(35,*) line
  read(35,'(5e22.14)')(occ(j),j=1,n)
  nElectrons(iFrag,ivec)=nint(sum(occ))
  close(35)
  return
end subroutine gcommon_read_vec


subroutine gcommon_getAtomicOverlap(filename,luOne,n,sAO)
  implicit none

  external :: OpnOne,RdOne,ClsOne
  
  integer,intent(in)    :: n,luOne
  integer   :: iCounter,iComponent
  integer   :: iRC,iOpt,iSymLbl,j,k

  real (kind=8)   :: s(n*(n+1)/2+4)
  real (kind=8),intent(out)   :: sAO(n,n)

  character (len=12),intent(in)     :: filename

  s=0.0
  sAO=0.0
  iRc=-1
  iOpt=0
  Call OpnOne(iRC,iOpt,filename,LuOne)
  if(iRC.ne.0) write(6,*)'Something went wrong opening ',filename
  iRC= 0
  iOpt=2
  iComponent=1
  iSymLbl=1
  Call RdOne(iRC,iOpt,'Mltpl  0',iComponent,s,iSymLbl)
  iCounter=1
  do j=1,n
    do k=1,j
      sAO(j,k)=s(iCounter)
      sAO(k,j)=s(iCounter)
      iCounter=iCounter+1
    enddo
  enddo
  iOpt=0
  Call ClsOne(iRc,iOpt)
  return
end subroutine gcommon_getAtomicOverlap


subroutine gcommon_calculate_sMO(a,sMO,sAO,n,m)
  use gcommon_input_data, only : debug 
  implicit none

  real(kind=8),intent(in)  :: sAO(m,m)
  real(kind=8),intent(in)  :: a(n,m)
  real(kind=8),intent(out) :: sMO(n,n)
  integer,intent(in) :: m,n
  integer:: i,j,k,l

  sMO=0.0
  do i=1,n
    do j=1,i
      do k=1,m
        do l=1,m
          sMO(i,j)=sMO(i,j)+a(i,k)*a(j,l)*sAO(k,l)
        enddo
      enddo
      sMO(j,i)=sMO(i,j)
    enddo
  enddo
  if(debug) then
    write(*,*) 'MO overlap matrix'
    do i=1,n
      write(*,'(20E18.8)') sMO(i,:)
    enddo
  endif
  return
end subroutine gcommon_calculate_sMO


subroutine gcommon_reverse_order(a,b,n)
  implicit none

  integer,intent(in)    :: n
  integer   :: i,j,halfway
  real (kind=8)   :: aux
  real (kind=8),intent(inout) :: a(n),b(n,n)

  halfway=int(n/2)
  do i=1,halfway
    aux=a(i)
    a(i)=a(n+1-i)
    a(n+1-i)=aux
    do j=1,n
      aux=b(i,j)
      b(i,j)=b(n+1-i,j)
      b(n+1-i,j)=aux
    enddo
  enddo
  return
end subroutine gcommon_reverse_order


subroutine gcommon_capitalize(string)
  implicit none
  integer:: i
  character(*) string

  do i=1,len(string)
    if(ichar(string(i:i)).gt.96) then
      string(i:i)=char(ichar(string(i:i))-32)
    endif
  enddo
  return
end subroutine gcommon_capitalize


subroutine gcommon_locate(string)
  implicit none
  external :: gcommon_capitalize
  character(4)   ::  string,string2
  character(132) ::  line
  rewind(5)
40 read(5,*) line
  string2=adjustl(line)
  call gcommon_capitalize(string2)
  if(string2.ne.string) goto 40
  return
end subroutine gcommon_locate


subroutine gcommon_getFilename(iVec,filename,base,suffix)
  use gcommon_input_data, only : fragLabel, mebfLabel
  implicit none
  integer,intent(in)  :: iVec
  character (len=3),intent(in)    :: suffix
  character (len=20),intent(in)   :: base
  character (len=255),intent(out) :: filename

  if(trim(suffix).eq.'vec') then
      write(filename,'(5A)') trim(base),trim(mebfLabel),trim(fragLabel(iVec)),'.',suffix
    else
      write(filename,'(4A)') trim(base),trim(fragLabel(iVec)),'.',suffix
    endif
  filename=trim(filename)
  return
end subroutine gcommon_getFilename


subroutine gcommon_printDim(a,n)
  use gcommon_histo_data
  implicit none

  integer,intent(in) :: n
  integer:: i,bin
  real (kind=8),intent(in) :: a(n)
  real (kind=8):: thrs

  thrs=0.1
  bin=1
  do i=1,n
    if(a(i).lt.thrs) then
      freq(bin)=freq(bin)+(i-1)
      thrs=thrs/10
      bin=bin+1
    endif
  enddo
  if(bin.gt.maxbin) maxbin=bin
  return
end subroutine gcommon_printDim


subroutine gcommon_printepsilons(a,n,i)
  implicit none

  integer,intent(in) :: n
  integer:: i,j
  real (kind=8),intent(in) :: a(n)

  write(*,*)
  write(*,'(a,I4)')'List of all eigenvalues of the MO overlap matrix of fragment ',i
  do j=1,n
    write(*,'(I5,E20.12)') j,a(j)
  enddo
  return
end subroutine gcommon_printepsilons


subroutine gcommon_ortho_frozen(frozenOrbs,totalFrozen,nBas)
  use gcommon_input_data, only : debug
  implicit none

  external :: NameRun,dsyev
  external :: gcommon_getAtomicOverlap,gcommon_calculate_sMO
  
  integer,intent(in)   :: nBas,totalFrozen
  integer  :: j,k,l
  integer  :: iRc,lwork
  integer  :: luOne

  real(kind=8),intent(inout) :: frozenOrbs (nBas,nBas)
  real(kind=8),allocatable   :: sAO  (:,:)
  real(kind=8),allocatable   :: sMO  (:,:)
  real(kind=8),allocatable   :: aMatrix    (:,:)
  real(kind=8),allocatable   :: eigenValues(:)
  real(kind=8),allocatable   :: work (:)
  real(kind=8),allocatable   :: orbs_debug (:,:)

  character (len=12)   :: oneintName

  lwork=4*totalFrozen


  allocate(sAO(nBas,nBas))
  allocate(sMO(totalFrozen,totalFrozen))
  allocate(aMatrix(totalFrozen,nBas))
  allocate(eigenValues(totalFrozen))
  allocate(work(lwork))
  allocate(orbs_debug(totalFrozen,nBas))

  do j=1,totalFrozen
    do k=1,nBas
      aMatrix(j,k)=frozenOrbs(j,k)  ! save the non-orthogonal frozen orbitals for later use
    enddo
  enddo
  frozenOrbs=0.0

  call NameRun('RUNFIL')
  write(oneintName,'(A6)') 'ONEINT'
  luOne=88
  call gcommon_getAtomicOverlap(oneintName,luOne,nBas,sAO)
  call gcommon_calculate_sMO(aMatrix,sMO,sAO,totalFrozen,nBas)
  iRc=1
  call dsyev('V','L',totalFrozen,sMO,totalFrozen,eigenValues,work,lwork,iRc)
  if(iRc.ne.0) then
    write(*,*) 'Something went wrong in dsyev (sMO frozen)'
    write(*,*) 'iRc=',iRc
  endif
  sMO=transpose(sMO)
  if(debug) then
    write (*,*) 'eigenvalues and eigenvectors of the average frozen SMO matrix'
    do j=1,totalFrozen
      write(*,'(F15.8,20E18.8)')eigenValues(j),sMO(j,:)
    enddo
  endif

! Lowdin orthonormalization (expressing the eigenvectors of sMO in AO basis)
  do j=1,totalFrozen
    do k=1,totalFrozen
      do l=1,nBas
        frozenOrbs(j,l)=frozenOrbs(j,l)+sMO(j,k)*aMatrix(k,l)/sqrt(eigenValues(j))
      enddo
    enddo
  enddo
  if(debug) then
    sMO=0.0
    write (*,*) 'Lowdin orthgonalized average frozen orbitals'
    do j=1,totalFrozen
      write(*,*) 'Average frozen orbital ',j
      write(*,'(5E22.14)') frozenOrbs(j,:)
      orbs_debug(j,:)=frozenOrbs(j,:)
    enddo
    call gcommon_calculate_sMO(orbs_debug,sMO,sAO,totalFrozen,nBas)
  endif
  deallocate(sAO)
  deallocate(sMO)
  deallocate(aMatrix)
  deallocate(eigenValues)
  deallocate(work)
  deallocate(orbs_debug)

  return
end subroutine gcommon_ortho_frozen


subroutine gcommon_calc_frz_density(n,nB,nF,sAO,froVec_1st,frzVec,frzDensity)
  use gcommon_input_data, only : debug
  implicit none

  integer,intent(in)   :: n,nB,nF
  integer  :: j
  real(kind=8),intent(in)    :: sAO(nB,nB)
  real(kind=8),intent(in)    :: frzVec(nF,nB),froVec_1st(nF,nB)
  real(kind=8),intent(inout) :: frzDensity(nF,nF)
  real(kind=8)   :: mprod(nB,nF)

  mprod=matmul(sAO,transpose(frovec_1st))
  frzDensity=frzDensity+matmul(frzVec,mprod)/n
  if(debug) then
    write(*,'(a)') '*  Accumulating the density matrices of the frozen orbitals'
    do j=1,nF
      write(*,'(20E18.8)') frzDensity(j,:)
    enddo
  endif

  return
end subroutine gcommon_calc_frz_density


subroutine gcommon_average_frozen(nF,nB,frzDensity,froVec_1st,frzFragAvg)
  use gcommon_input_data, only : debug
  implicit none

  external :: dsyev
  
  integer,intent(in)    :: nF,nB
  integer               :: j
  integer               :: iRc

  real(kind=8),intent(inout)   :: frzDensity(nF,nF)
  real(kind=8),intent(in)      :: froVec_1st(nF,nB)
  real(kind=8),intent(out)     :: frzFragAvg(nF,nB)
  real(kind=8)                 :: work(4*nF)
  real(kind=8)                 :: eigenValues(nF)

  iRc=1
  work=0.0
  call dsyev('V','L',nF,frzDensity,nF,eigenValues,work,4*nF,iRc)
  if(iRc.ne.0) then
    write(*,*) 'Something went wrong in dsyev (frzDensity)'
    write(*,*) 'iRc=',iRc
  endif
  frzDensity=transpose(frzDensity)! eigenvectors as rows
  if(debug) then
    write(*,*)'eigenvectors of SMO(frozen)'
    do j=1,nF
      write(*,'(A,F18.10)') 'eigenvalue :',eigenValues(j)
      write(*,'(20E18.8)') frzDensity(j,:)
    enddo
  endif
  frzFragAvg=matmul(frzDensity,froVec_1st)
  if(debug) then
    write(*,*)'average frozen orbitals' 
    do j=1,nF
      write(*,'(5E22.14)') frzFragAvg(j,:)
    enddo
  endif

  return
end subroutine gcommon_average_frozen


subroutine gcommon_add_detinfo()
  use gcommon_input_data
  use gcommon_fragment_data
  implicit none

  external :: gcommon_getFilename,gcommon_quicksort
  
  integer :: idet,ndet,inactm,i,j,iVec,istat,ndet_unique
  real(kind=8),allocatable  :: coeff(:),coeff_unique(:)
  character(len=3)      :: suffix
  character(len=255)    :: detFilename
  character (len=255), allocatable :: occ(:),occ_unique(:)

  suffix = 'det'
  iVec = 0
  do i = 1, nFragments
    do j = 1, nVec(i)
      ndet = 0
      iVec = iVec + 1
      call gcommon_getFilename(iVec,detFilename,project,suffix)
      open(37,file=detFilename,status='old')
      read(37,1500) inactm,ener(i,j),enerpt(i,j)
1500 format(i5,27x,f22.12,15x,f22.12)
      do
        read(37,*,iostat=istat)
        if (istat .ne. 0) then
          rewind(37)
          exit
        else
          ndet = ndet + 1
        endif
      enddo
      allocate(coeff(ndet))
      allocate(coeff_unique(ndet))
      allocate(occ(ndet))
      allocate(occ_unique(ndet))
      read(37,*)
      do idet = 1, ndet
        read(37,*) coeff(idet),occ(idet)
      end do
      rewind(37)
      call gcommon_quicksort(coeff,occ,ndet)
      ndet_unique = 1
      occ_unique(1) = occ(1)
      coeff_unique(1) = coeff(1)
      do idet = 2, ndet
        if ( occ(idet) .ne. occ(idet-1) ) then
          ndet_unique = ndet_unique + 1
          occ_unique(ndet_unique) = occ(idet)
          coeff_unique(ndet_unique) = coeff(idet)
        else
          coeff_unique(ndet_unique) = coeff_unique(ndet_unique) + &
              coeff(idet)
        end if
      end do
      write(*,'(2i4,2f22.12)') i,j,ener(i,j),enerpt(i,j)
      if (fragLabels.and.energy_on_INPORB) then
        write(37,1600) inactm,nFrozen(i),ndet_unique, &
            fragState(iVec),ener(i,j),nElectrons(i,j),threshold,enerpt(i,j)
      elseif (fragLabels) then
        write(37,1600) inactm,nFrozen(i),ndet_unique, &
            fragState(iVec),ener(i,j),nElectrons(i,j),threshold,enerpt(i,j)
      elseif (energy_on_INPORB) then
        write(37,1600) inactm,nFrozen(i),ndet_unique,'no_label', &
            ener(i,j),nElectrons(i,j),threshold,enerpt(i,j)
      else
        write(37,1600) inactm,nFrozen(i),ndet_unique,'no_label', &
            ener(i,j),nElectrons(i,j),threshold,enerpt(i,j)
      end if
      do idet = 1, ndet_unique
        write(37,'(e15.8,6x,A)') coeff_unique(idet),trim(occ_unique(idet))
      end do
      close(37)
      deallocate(coeff,occ)
      deallocate(coeff_unique,occ_unique)
    end do
1600 format(2i5,i12,4x,a2,4x,f22.12,i5,1pe10.3,0pf22.12)
  end do
  
end subroutine gcommon_add_detinfo

subroutine gcommon_generate_vecfiles(iFrag,nBas)
  use gcommon_input_data, only : project,nVec
  implicit none

  external :: gcommon_getFilename
  
  integer :: i,j,iFrag,nBas,iVec,start
  real(kind=8),allocatable  :: vec(:,:)
  character(len=255)  :: filename

  start=0
  if(iFrag.gt.1) then
    do i=2,iFrag
      start=start+nVec(i)
    enddo
  endif
  do iVec=1,nVec(iFrag)
    call gcommon_getFilename(start+iVec,filename,project,'vec')
    open(10,file=filename)
    write(10,'(i4)')nBas
    allocate(vec(nBas,nBas))
    do i=1,nBas
      vec(i,i)=1.0
    enddo
    do i=1,nBas
      write(10,'(4F18.14)')(vec(i,j),j=1,nBas)
    enddo
    deallocate(vec)
  enddo
end subroutine gcommon_generate_vecfiles


recursive subroutine gcommon_quicksort(coef,occu,n)
  implicit none

  real (kind=8)    :: coef(n)
  character(len=255)  :: occu(n),pivot,aux2
  real (kind=8)    :: random
  real (kind=8)    :: aux1
  integer:: n,left,right,marker

  if(n.gt.1) then
    call random_number(random)
    pivot=occu(int(random*real(n-1))+1)
    left=0
    right=n+1
    do while(left.lt.right)
      right=right-1
      do while(occu(right).lt.pivot)
        right=right-1
      enddo
      left=left+1
      do while(occu(left).gt.pivot)
        left=left+1
      enddo
      if(left.lt.right) then
        aux1=coef(left)
        coef(left)=coef(right)
        coef(right)=aux1
        aux2=occu(left)
        occu(left)=occu(right)
        occu(right)=aux2
      endif
    enddo
    if(left.eq.right) then
      marker=left+1
    else
      marker=left
    endif
    call gcommon_quicksort(coef(:marker-1),occu(:marker-1),marker-1)
    call gcommon_quicksort(coef(marker:),occu(marker:),n-marker+1)
  endif
end subroutine gcommon_quicksort
