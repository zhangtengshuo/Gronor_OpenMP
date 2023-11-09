******************************************************************************
*                                                                            *
*     program for overlapping fragments                                      *
*                                                                            *
* Corresponding orbitals of fragment A and B are determined and written      *
*   in SUPERORB. For corresponding orbitals with a singular value larger     *
*   than 'thrs' (default = 0.2), only the (A+B)/2 orbital is written. In     *
*   this way, the overlap between A and B is eliminated.                     *
*   The program can also eliminate basis functions for atoms that are        *
*   present in the fragments but do not appear in the supermolecule.         *
*   Typically, hydrogens that were used to saturate bond that were broken    *
*   in the definition of the fragments from the supermolecule.               *
*                                                                            *
* INPUT                                                                      *
*   - INPORB1     :   Orbitals of fragment A                                 *
*   - INPORB2     :   Orbitals of fragment B                                 *
*      overlapping atoms should be the last ones in fragment A and the       *
*      first ones in fragment B                                              *
*   - RUNFILE     :   RunFile from the supermolecule                         *
*   - ONEINT      :   One-electron integrals of the supermolecule            *
*   - GUESSORB    :   Any vector file of the supermolecule (e.g. GssOrb)     *
*      Only used to write the virtual orbitals in SUPERORB                   *
*                                                                            *
*  OUTPUT                                                                    *
*   - SUPERORB    : Vector file with the superimposed fragments              *
*   - CORRORB.A   : Corresponding orbitals of fragment A                     *
*   - CORRORB.B   : Corresponding orbitals of fragment B                     *
*                                                                            *
*      The atoms that are to be eliminated from the supermolecule            *
*      should be labeled with the letter 'Q'                                 *
*                                                                            *
* To be linked to molcas libraries for accessing ONEINT and RUNFILE, the     *
*    program also uses the 'dgesvd' routine from lapack                      *
*                                                                            *
*    Coen de Graaf, Sept. 2019, beta version                                 *
*        the program seems to work, but more testing is highly desirable     * 
*                                                                            *
******************************************************************************
      program main
      implicit none

      external :: dgesvd,ClsOne,RdOne,OpnOne,Get_cArray,Get_iArray
      external :: NameRun
      
      integer                       :: iSym,nSym
      integer                       :: lDim
      integer                       :: j,k,l,m1,m2,m3,jOrb
      integer                       :: LuOne,iRC,iOpt,iComponent,iSymLbl
      integer                       :: nBasTot,iCounter,to_be_averaged
      integer                       :: lTriangle
      integer                       :: dummy
      integer                       :: maxInact,nBasFinal
      integer, allocatable          :: nBas1(:),nBas2(:),nBas(:)
      integer, allocatable          :: nInact1(:),nInact2(:)
      integer, allocatable          :: nAct1(:),nAct2(:)
      integer, allocatable          :: n_elim(:)

      real (kind = 8)               :: thrs
      real (kind = 8), allocatable  :: c1(:,:,:),c2(:,:,:),c3(:,:,:)
      real (kind = 8), allocatable  :: c_1(:,:),c_2(:,:)
      real (kind = 8), allocatable  :: c1c(:,:),c2c(:,:)
      real (kind = 8), allocatable  :: occNu1(:,:),occNu2(:,:)
      real (kind = 8), allocatable  :: occNu3(:,:)
      real (kind = 8), allocatable  :: s(:),sAO(:,:,:),sMO(:,:)
      real (kind = 8), allocatable  :: U(:,:),VT(:,:),V(:,:)
      real (kind = 8), allocatable  :: aux(:,:),work(:)
      real (kind = 8), allocatable  :: sigma(:)

      character (len = 132)         :: line,title,command_line
      character (len = 6)           :: mark,dummy_string
      character (len = 14) , allocatable    :: basLabel(:)
      character (len = 1 ) , allocatable    :: orbLabel(:)

      logical , allocatable         :: eliminate(:,:)
      logical                       :: debug

      call get_command(command_line)
      read(command_line,*) dummy_string,thrs
      write(*,'(A,F8.3)') 'Running ovlf with threshold: ',thrs
      debug = .false.

*...Open the files where the corresponding orbitals will be stored
      open (12,file='SUPERORB',status='unknown')
      open (13,file='CORRORB.A',status='unknown')
      open (14,file='CORRORB.B',status='unknown')
*
*...Retrieving info from  INPORB1
      open (9,file='INPORB1',status='old')
      read(9,'(A132)') line
      write(12,'(A132)') line
      write(13,'(A132)') line
      read(9,'(A132)') line
      write(12,'(A132)') line
      write(13,'(A132)') line
      read(9,'(A132)') title
      write(6,*) "Title vector file 1 : ",title
      write(12,'(A32)')'* Superimposed fragment orbitals'
      write(13,'(A35)')'* Corresponding orbitals fragment A'
      read(9,*) dummy, nSym, dummy
      write(12,'(3I8)')0,nSym,0
      write(13,'(3I8)')0,nSym,0
      allocate ( nBas1(nSym) )
      nBas1 = 0.0
      read(9,*)(nBas1(iSym),iSym=1,nSym)
      write(6,'(A,8I5)')"Number of basis functions : ",
     &                         (nBas1(iSym),iSym=1,nSym)
      write(13,'(8I8)')(nBas1(iSym),iSym = 1, nSym)
      write(13,'(8I8)')(nBas1(iSym),iSym = 1, nSym)
      write(13,'(A4)')'#ORB'
      lDim = maxval(nBas1)
      allocate( c1(nSym,lDim,lDim) )
      allocate( occNu1(nSym,lDim) )
      c1 = 0.0
      occNu1 = 0.0
      mark = '#ORB'
 44   read(9,'(A132)') line
      if (line(1:4).ne.mark) goto 44
      do iSym = 1, nSym
        if (nBas1(iSym).ne.0) then
          do j = 1, nBas1(iSym)
            read(9,'(A132)') line
            read(9,'(5E22.14)') (c1(iSym,j,k),k=1,nBas1(iSym))
          end do
        endif
        if ( debug ) then
          write(6,*) 'MO coefficients, irrep ',iSym
          do j = 1, nBas1(isym)
            write(6,'(20F10.5)')(c1(iSym,j,k),k=1,nBas1(iSym))
          end do
        end if
      end do
      rewind(9)
      mark = '#OCC'
 45   read(9,'(A132)') line
      if (line(1:4).ne.mark) goto 45
      read(9,'(A132)') line
      do iSym = 1, nSym
        if (nBas1(iSym).ne.0) then
          read(9,'(5E22.14)') (occNu1(iSym,j),j=1,nBas1(iSym))
        endif
      end do
      rewind(9)
      mark = '#INDEX'
      allocate( orbLabel(lDim) )
      allocate(  nInact1(nSym) )
      allocate(    nAct1(nSym) )
      nInact1 = 0
      nAct1 = 0
 46   read(9,'(A132)') line
      if (line(1:6).ne.mark) goto 46
      read(9,'(A132)') line
      do iSym = 1, nSym
        orbLabel = ' '
        if ( nBas1(iSym) .ne. 0 ) then
          read(9,'(2x,10A)')(orbLabel(j),j=1,nBas1(iSym))
        end if
        do j = 1, nBas1(iSym)
          if (orbLabel(j) .eq. 'i') nInact1(iSym) = nInact1(iSym) + 1
          if (orbLabel(j) .eq. '2') nAct1(iSym) = nAct1(iSym) + 1
        end do
      end do
      write(6,'(A,8I5)') '        inactive orbitals : ',
     &                         (nInact1(iSym),iSym=1,nSym)
      write(6,'(A,8I5)') '          active orbitals : ',
     &                         (nAct1(iSym),iSym=1,nSym)
      write(6,*)
      deallocate( orbLabel )
      close(9)

*...Retrieving info from  INPORB2
      open (10,file='INPORB2',status='old')
      read(10,'(A132)') line
      write(14,'(A132)') line
      read(10,'(A132)') line
      write(14,'(A132)') line
      read(10,'(A132)') title
      write(6,*) "Title vector file 2 : ",title
      write(14,'(A35)')'* Corresponding orbitals fragment B'
      read(10,*) dummy,nSym,dummy
      write(14,'(3I8)')0,nSym,0
      allocate ( nBas2(nSym) )
      nBas2 = 0.0
      read(10,*) (nBas2(iSym),iSym=1,nSym)
      write(6,'(A,8I5)')"Number of basis functions : ",
     &                         (nBas2(iSym),iSym=1,nSym)
      write(14,'(8I8)')(nBas2(iSym),iSym = 1, nSym)
      write(14,'(8I8)')(nBas2(iSym),iSym = 1, nSym)
      write(14,'(A4)')'#ORB'
      lDim = maxval(nBas2)
      allocate( orbLabel(lDim) )
      allocate( c2(nSym,lDim,lDim) )
      allocate( occNu2(nSym,lDim) )
      c2 = 0.0
      occNu2 = 0.0
      mark = '#ORB'
 54   read(10,'(A132)') line
      if (line(1:4).ne.mark) goto 54
      do iSym = 1, nSym
        if (nBas2(iSym).ne.0) then
          do j = 1, nBas2(iSym)
            read(10,'(A132)') line
            read(10,'(5E22.14)') (c2(iSym,j,k),k=1,nBas2(iSym))
          end do
        endif
        if ( debug ) then
          write(6,*) 'MO coefficients, irrep ',iSym
          do j = 1, nBas2(isym)
            write(6,'(20F10.5)')(c2(iSym,j,k),k=1,nBas2(iSym))
          end do
        end if
      end do
      rewind(10)
      mark = '#OCC'
 55   read(10,'(A132)') line
      if (line(1:4).ne.mark) goto 55
      read(10,'(A132)') line
      do iSym = 1, nSym
        if (nBas2(iSym).ne.0) then
          read(10,'(5E22.14)') (occNu2(iSym,j),j=1,nBas2(iSym))
        endif
      end do
      rewind(10)
      mark = '#INDEX'
 56   read(10,'(A132)') line
      if (line(1:6).ne.mark) goto 56
      read(10,'(A132)') line
      allocate(  nInact2(nSym) )
      allocate(    nAct2(nSym) )
      nInact2 = 0
      nAct2 = 0
      do iSym = 1, nSym
        orbLabel = ' '
        if ( nBas1(iSym) .ne. 0 ) then
          read(10,'(2x,10A)')(orbLabel(j),j=1,nBas2(iSym))
        end if
        do j = 1, nBas2(iSym)
          if (orbLabel(j) .eq. 'i') nInact2(iSym) = nInact2(iSym) + 1
          if (orbLabel(j) .eq. '2') nAct2(iSym) = nAct2(iSym) + 1
        end do
      end do
      write(6,'(A,8I5)') '        inactive orbitals : ',
     &                         (nInact2(iSym),iSym=1,nSym)
      write(6,'(A,8I5)') '          active orbitals : ',
     &                         (nAct2(iSym),iSym=1,nSym)
      write(6,*)
      close(10)
*...Open the RunFile and read some info (number of basis functions and basis labels)
      allocate ( nBas(nSym) )
      nBas = 0.0
      Call NameRun('RUNFILE')
      Call Get_iArray('nBas',nBas,nSym)
      nBasTot = sum(nBas)
      lTriangle = ( nBasTot * ( nBasTot + 1 ) ) / 2
      lDim = maxval(nBas)
      allocate ( basLabel(nBasTot) )
      allocate ( eliminate(nSym,lDim) )
      allocate ( n_elim(nSym) )
      n_elim = 0
      Call Get_cArray('Unique Basis Names',basLabel,14*nBasTot)
      do iSym = 1, nSym
        do j = 1, nBas(iSym)
          if (basLabel(j)(1:1) .eq. 'Q') then
            eliminate(iSym,j) = .True.
            n_elim(iSym) = n_elim(iSym) + 1 
          else
            eliminate(iSym,j) = .False.
          end if
        end do
      end do
      write(12,'(8I8)')(nBas(iSym)-n_elim(nSym),iSym = 1, nSym)
      write(12,'(8I8)')(nBas(iSym)-n_elim(nSym),iSym = 1, nSym)
      write(12,'(A4)')'#ORB'
*...Open the OneInt file and read the overlap matrix of the ao-basis
      allocate ( s(lTriangle) )
      s = 0.0
      allocate ( sAO(nSym,lDim,lDim) )
      sAO = 0.0
      LuOne = 77
      iRc=-1
      iOpt=0
      Call OpnOne(iRC,iOpt,'ONEINT',LuOne)
      if (iRC.ne.0) write(6,*) 'Something went wrong opening ONEINT'
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
        if ( debug ) then
          write(6,*) 'AO overlap'
          do j = 1,nBas(iSym)
            write(6,'(20F10.5)')(sAO(iSym,j,k),k=1,nBas(iSym))
          end do
          write(6,*)
        end if
      end do
      deallocate( s )
      Call ClsOne(iRc,iOpt)
*...Open GUESSORB of the supermolecule for some reasonable virtual orbitals
      open(11,file='GUESSORB',status='old')
      allocate( c3(nSym,lDim,lDim) )
      allocate( occNu3(nSym,lDim) )
      c3 = 0.0
      occNu3 = 0.0
      mark = '#ORB'
 64   read(11,'(A132)') line
      if (line(1:4).ne.mark) goto 64
      do iSym = 1, nSym
        if (nBas(iSym).ne.0) then
          do j = 1, nBas(iSym)
            read(11,'(A132)') line
            read(11,'(5E22.14)') (c3(iSym,j,k),k=1,nBas(iSym))
          end do
        endif
      end do
      mark = '#OCC'
 65   read(11,'(A132)') line
      if (line(1:4).ne.mark) goto 65
      read(11,'(A132)') line
      do iSym = 1, nSym
        if (nBas(iSym).ne.0) then
          read(11,'(5E22.14)') (occNu3(iSym,j),j=1,nBas(iSym))
        endif
      end do
      close(11)
*...Replace the first vectors by the active orbitals of the fragments 
*     to save them for the final vector file SUPERORB
      do iSym = 1, nSym
        do j = 1, nAct1(iSym)
          do k = 1, nBas(iSym)
            c3(iSym,j,k) = 0.0
          end do
          do k = 1, nBas1(iSym)
            c3(iSym,j,k) = c1(iSym,j+nInact1(iSym),k)
          end do
        end do
        do j = nAct1(iSym) + 1, nAct1(iSym) + nAct2(iSym)
          do k = 1, nBas(iSym)
            c3(iSym,j,k) = 0.0
          end do
          do k = 1, nBas2(iSym)
            c3(iSym,j,k+nBas(iSym)-nBas2(iSym)) =
     &                 c2(iSym,j-nAct1(iSym)+nInact2(iSym),k)
          end do
        end do
      end do
*...End of data collection
*
*
*...Construct overlap matrix of the MOs, sMO = (c2)^T x sAO x c1 
      do iSym = 1, nSym
        m1 = nInact1(iSym)+nInact2(iSym)
        m2 = nBas(iSym)
        m3 = nBas(iSym)-nBas2(iSym)
        if (debug) write (6,*) 'm1 m2 m3 : ',m1,m2,m3
        allocate ( aux(m1,m2) )
        allocate ( sMO(m1,m1) )
        allocate (   U(m1,m1) )
        allocate (  VT(m1,m1) )
        allocate (   V(m1,m1) )
        allocate (  sigma(m1) )
        allocate ( work(5*m1) )
        allocate ( c_1(m1,m2) )
        allocate ( c_2(m1,m2) )
        allocate ( c1c(m1,m2) )
        allocate ( c2c(m1,m2) )
        c_1 = 0.0
        c_2 = 0.0
        c1c = 0.0
        c2c = 0.0
* Adding zeros to c1 and c2 to fit the dimension of the "super molecule"
        do j = 1, nInact1(iSym)
          do k = 1, nBas1(iSym)
            c_1(j,k) = c1(iSym,j,k)
          end do
        end do
        do j = 1, nInact2(iSym)
          do k = 1, nBas2(iSym)
            c_2(j+nInact1(iSym),k + m3) = c2(iSym,j,k)
          end do
        end do
        if ( debug ) then
          write(6,*)'c_1'
          do j = 1, m1
            write(6,'(20F10.5)') (c_1(j,k),k=1,m2)
          end do
          write(6,*)'c_2'
          do j = 1, m1
            write(6,'(20F10.5)') (c_2(j,k),k=1,m2)
          end do
        end if
        aux = 0.0
        sMO = 0.0
        do j = 1, m1
          do k = 1, m2
            do l = 1, m2
              aux(j,k) = aux(j,k) + c_2(j,l) * sAO(iSym,l,k) 
            end do
          end do
        end do
        do j = 1, m1
          do k = 1, m1
            do l = 1, m2
              sMO(j,k) = sMO(j,k) + aux(j,l) * c_1(k,l)
            end do
          end do
        end do
        if (debug) then
          write(*,*) 'MO overlap'
          do j = 1, m1
            write(6,'(20F10.5)')(sMO(j,k),k=1,m1)
          end do
          write(6,*)
        end if
* Singular value decomposition of sMO
        call dgesvd('A','A',m1,m1,sMO,m1,sigma,U,m1,VT,m1,work,5*m1,iRC)
        V = Transpose(VT)
        if ( debug ) then
          write(6,*)'Singular Value Decomposition'
          write(6,*) 'left sigular matrix, U'
          do j = 1, m1
            write(6,'(20F10.5)') (U(j,k),k = 1, m1)
          end do
          write(6,*) 'right sigular matrix, V'
          do j = 1, m1
            write(6,'(20F10.5)') (V(j,k),k = 1, m1)
          end do
        end if
        write(6,*) 'Singular values, sigma'
        write(6,'(10F10.5)') (sigma(j),j = 1, m1)
        write(6,*)
        do j = 1, m1
          do k = 1, m2
            do l = 1, m1
              c1c(j,k) = c1c(j,k) + c_1(l,k) * V(l,j)
              c2c(j,k) = c2c(j,k) + c_2(l,k) * U(l,j)
            end do
          end do
        end do
        if ( debug ) then
          write(6,*)
          write(6,*) 'Corresponding orbitals'
          write(6,*) 'transformed c1'
          do j = 1, m1
            write(6,'(20F10.5)') (c1c(j,k),k = 1, m2)
          end do
          write(6,*) 'transformed c2'
          do j = 1, m1
            write(*,'(20F10.5)') (c2c(j,k),k = 1, m2)
          end do
        end if
* Dump the corresponding orbitals on CORRORB.A and CORRORB.B
        do j = 1, nInact1(iSym)
          write(13,'(A,2I5)') '* ORBITAL',iSym,j
          write(13,'(5E22.14)')(c1c(j,k),k = 1, nBas1(iSym))
        end do
        do j = nInact1(iSym)+1,nBas1(iSym)
          write(13,'(A,2I5)') '* ORBITAL',iSym,j
          write(13,'(5E22.14)')(c1(iSym,j,k),k = 1, nBas1(iSym))
        end do
        do j = 1, nInact2(iSym)
          write(14,'(A,2I5)') '* ORBITAL',iSym,j
          write(14,'(5E22.14)')(c2c(j,k),k = m3+1, m2)
        end do
        do j = nInact2(iSym)+1, nBas2(iSym)
          write(14,'(A,2I5)') '* ORBITAL',iSym,j
          write(14,'(5E22.14)')(c2(iSym,j,k),k = 1, nBas2(iSym))
       end do
* Check that we have <phi_A,i | phi_B,j > =  sigma_ij (only non-zero for i = j)
       if ( debug ) then
          aux = 0.0
          do j = 1, m1 
            do k = 1, m2 
              do l = 1, m2 
                aux(j,k) = aux(j,k) + c2c(j,l) * sAO(iSym,l,k)
              end do
            end do
          end do
          sMO = 0.0
          do j = 1, m1 
            do k = 1, m1 
              do l = 1, m2
                sMO(j,k) = sMO(j,k) + aux(j,l) * c1c(k,l)
              end do
            end do
          end do
          write(6,*) 'Corresponding orbitals overlap'
          do j = 1, m1 
            write(6,'(20F10.5)')(sMO(j,k),k=1,m1)
          end do
          write(*,*)
        end if
* Before writing the corresponding orbitals to SUPERORB (unit 12), we first
* eliminate (if needed) the basis function(s) associated to atoms not present
* in the final supermolecule (typically saturating hydrogens)
        l = 0
        do k = 1, nBas(iSym)
          if ( .not.eliminate(iSym,k) ) then
            l = l + 1
            do j = 1, nInact1(iSym)
              c1c(j,l) = c1c(j,k)
            end do
            do j = 1, nInact2(iSym)
              c2c(j,l) = c2c(j,k)
            end do
            do j = 1, nBas(iSym)
              c3(iSym,j,l) = c3(iSym,j,k)
            end do
          end if
        end do
* dump corresponding orbitals in orbital file of the supermolecule 
* eliminating those orbitals for which the singular value is larger than 0.3
        to_be_averaged = 0
*        thrs = 0.10
        nBasFinal = nBas(iSym) - n_elim(iSym)
        do j = 1, m1
          if ( sigma(j) .ge. thrs ) to_be_averaged = to_be_averaged + 1
        end  do
        write (6,'(I4,A)') m1 - to_be_averaged, 
     &                  ' inactive orbitals will be dumped on SUPERORB'
        write(6,*) 'Please, check!!',
     &             ' If not correct, change --thrs-- and run again'
        maxInact = max(nInact1(iSym),nInact2(iSym))
        jOrb = 1
        do j = 1, maxInact
          if ( sigma(jOrb) .gt. thrs ) then
            write(12,'(A,2I5)') '* ORBITAL',iSym,jOrb
            write(12,'(5E22.14)')((c1c(j,k)+c2c(j,k))/2,k=1,nBasFinal)
            jOrb = jOrb + 1
          else
            if (j .le. nInact1(iSym)) then
              write(12,'(A,2I5)') '* ORBITAL',iSym,jOrb
              write(12,'(5E22.14)')(c1c(j,k),k=1,nBasFinal)
              jOrb = jOrb + 1
            endif
            if (j .le. nInact2(iSym)) then
              write(12,'(A,2I5)') '* ORBITAL',iSym,jOrb
              write(12,'(5E22.14)')(c2c(j,k),k=1,nBasFinal)
              jOrb = jOrb + 1
            endif
          endif
        end do
        do j = 1, nAct1(iSym) + nAct2(iSym)
          write(12,'(A,2I5)')'* ORBITAL',iSym,j+jOrb-1
          write(12,'(5E22.14)')(c3(iSym,j,k),k=1,nBasFinal)
        end do
        jOrb = jOrb + nAct1(iSym) + nAct2(iSym)
        do j = jOrb, nBasFinal
          write(12,'(A,2I5)')'* ORBITAL',iSym,j
          write(12,'(5E22.14)')(c3(iSym,j,k),k=1,nBasFinal)
        end do
* deallocate for next symmetry
        deallocate ( aux )
        deallocate ( sMO )
        deallocate (   U )
        deallocate (  VT )
        deallocate (   V )
        deallocate ( sigma )
        deallocate ( work )
        deallocate ( c_1 )
        deallocate ( c_2 )
        deallocate ( c1c )
        deallocate ( c2c )
      end do
* writing some dummy occupation numbers
      write(6,'(A,8I4)') 
     &   'Number of basis functions of the supermolecule: ',
     &   (nBas(iSym)-n_elim(iSym),iSym = 1, nSym)   
      write(12,'(A4)')'#OCC'
      write(12,'(A20)') '* OCCUPATION NUMBERS'
      write(13,'(A4)')'#OCC'
      write(13,'(A20)') '* OCCUPATION NUMBERS'
      write(14,'(A4)')'#OCC'
      write(14,'(A20)') '* OCCUPATION NUMBERS'
      do iSym = 1, nSym
       write(12,'(5E22.14)')(occNu3(iSym,j),j=1,nBas(iSym)-n_elim(iSym))
       write(13,'(5E22.14)')(occNu1(iSym,j),j=1,nBas1(iSym))
       write(14,'(5E22.14)')(occNu2(iSym,j),j=1,nBas2(iSym))
      end do
*
* clean up, deallocate all the rest
      close(12)
      close(13)
      close(14)
      deallocate (nBas)
      deallocate (nBas1)
      deallocate (nBas2)
      deallocate (c1)
      deallocate (c2)
      deallocate (c3)
      deallocate (occNu1)
      deallocate (occNu2)
      deallocate (occNu3)
      deallocate (nInact1)
      deallocate (nInact2)
      deallocate (nAct1)
      deallocate (nAct2)
      deallocate (sAO)
      deallocate (n_elim)
      deallocate (eliminate)
      deallocate (basLabel)
      deallocate (orbLabel)
      end



