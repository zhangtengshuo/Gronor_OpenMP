      program read_cholesky
*
* program to read the cholesky vectors transformed to the MO basis,
* reconstruct the 2-electron integrals and write them to TRAINT
*
* the following files are required as input
*    CHRED, CHVEC1, CHORST, CHOMAP and _CHMOT1
* the first four are present in the scratch dir of the molcas job as
*    $Project.ChRed, $Project.ChVec1, $Project.ChRst, $Project.ChMap
*
#ifdef _OPENMP
      use omp_lib
#endif
      use Cholesky

      implicit none

      external :: daclos,ddafile,Cho_X_nVecRS,GetMem,daname
      external :: Get_iArray,NameRun,Cho_X_Init
      
      integer                    :: iVrs,nVrs,iRc,luChVec,idisk
      integer                    :: luTra,iad50
      integer                    :: i,j,i1,i2
      integer                    :: ipq,irs,npq,lpqrs,totLength
      integer                    :: nBas,nDel,nOrb,nFrozen
      integer                    :: iBuf,nBuf,bufLength,lastBuf
      integer                    :: nTransCho,iOffset,iVec
      integer                    :: nodes,iCounter,inode,packages
      integer                    :: extra,iCounter2,ipack
      integer, parameter         :: nTraToc = 106
      integer, dimension(nTraToc):: iTraToc
      integer, allocatable       :: length(:)

      real (kind=8)              :: fracMem
      real (kind=8), allocatable :: coeff(:),integral(:),coeff2(:)

      integer :: numint,nthreads,mythread
      integer, allocatable :: ndx(:)

      integer  :: ip_Max,l_Max   ! to resolve a compilation issue

      character(len=10)          :: num_of_nodes

! #include "cholesky.fh"
! #include "choptr.fh"
#include "WrkSpc.fh"

      call get_command_argument(1,num_of_nodes)
      read(num_of_nodes,*) nodes
      call NameRun('RUNFILE')
      call Get_iArray('nBas',nBas,1)
      call Get_iArray('nDel',nDel,1)
      call Get_iArray('nFro',nFrozen,1)
* for the moment, no frozen orbitals, just eliminating the 'virtuals' (redundant orbitals in the common basis)
      nOrb = nBas - nDel -nFrozen
      npq = ( nOrb * (nOrb + 1) ) / 2
      lpqrs = (npq * (npq + 1) ) / 2
      bufLength = 9600
      nBuf = int(lpqrs / bufLength)
      lastBuf = lpqrs - nBuf * bufLength
      write(*,101) 'Number of nodes in MOTRA   : ',nodes
      write(*,101) 'Number of basis functions  : ',nBas
      write(*,101) 'Number of deleted orbitals : ',nDel
      write(*,101) 'Number of full buffers     : ',nBuf
      write(*,101) 'Length of last buffer      : ',lastBuf
      write(*,101) 'Number of integrals        : ',lpqrs
      write(*,*)
 101  format (a,i10)

      luTra = 50
* the code crashes if there is no TRAINT
      open(luTra,file='TRAINT')
      close(luTra)
      call daname(luTra,'TRAINT')
      iad50 = 0
      iTraToc = 0
      call ddafile(luTra,1,iTraToc,nTraToc,iad50)
      luChVec = 70
      call daname(luChVec,'_CHMOT1')
      
* Initialize the whole Cholesky business (for this we need CHRED, CHVEC1, CHORST and CHOMAP)
      iRc = 0
      fracMem = 0.0
* Due to compilation problems, this call was included
      Call GetMem('CXI_MX1','Max ','Real',ip_Max,l_Max)
* 
      call Cho_X_Init(iRc,fracMem)
      write(*,101) 'Number of Cholesky vectors : ',numCho(1)

* First call for finding out the number of transformed Cholesky vectors
      nTransCho = 0
      iVrs = 0
      nVrs = 0
      do while ( (iVrs + nVrs) .le. numCho(1) )
        nTransCho = nTransCho + 1
        call Cho_X_nVecRS(nTransCho,1,iVrs,nVrs)
      end do
      write(*,101) 'Transformed Chol. vectors  : ',nTransCho
      write(*,*)
* Second call for storing the length of the different vectors
      allocate( length(nTransCho) )
      do i = 1, nTransCho
        call Cho_X_nVecRS(i,1,iVrs,nVrs)
        length(i) = nVrs
      end do
      totLength = numCho(1) * npq 
      allocate( coeff (totLength) )
      allocate( coeff2(totLength) )
      coeff2 = 0.0
* Read the vectors (all of them)
      idisk = 0
      call ddafile(luChVec,2,coeff,totLength,idisk)
* If nodes larger than one, coeff needs to be ordered
      if ( nodes .ne. 1 ) then
        packages = int( totLength / (nodes * npq) )
        extra = ( totLength - packages * nodes * npq ) / npq
*        write(*,*) 'totLength npq ',totLength,npq
*        write(*,*) 'distribution over the nodes : ',packages,extra
        iCounter = 0
        do inode = 1, nodes
          iCounter2= inode
          do ipq = 1, npq
            do ipack = 1, packages
              iCounter = iCounter + 1
              coeff2(iCounter2) = coeff(iCounter)
              iCounter2= iCounter2 + nodes
            end do
            if (inode .le. extra ) then
              iCounter = iCounter + 1
              coeff2(iCounter2) = coeff(iCounter)
              iCounter2 = iCounter2 + extra
            else
              iCounter2 = iCounter2 + extra
            end if
          end do
        end do
*        write(*,*) 'raw data after ordering'     !  These lines can be activated for debugging
        do i = 1, totLength
*          write(*,'(I10,F18.8)')i,coeff2(i)      !  Here for motra in parallel
          coeff(i) = coeff2(i)
        end do
        deallocate(coeff2)
      end if
*      if ( nodes .eq. 1 ) then
*        write(*,*) 'raw data'                      ! This is the reference list of coefficients
*        do i = 1, totLength                        ! for motra in sequential
*          write(*,'(I10,F18.8)')i,coeff(i)
*        end do
*        write(*,*) 
*      end if
*     Reconstruction of the integrals

      numint=npq*(npq+1)/2
      allocate(ndx(npq))
      i=0
      do ipq=1,npq
        ndx(ipq)=i
        i=i+npq-ipq+1
      enddo
      allocate(integral(numint))
!$omp parallel private(nthreads,mythread,ipq,i,irs,iOffset,j,i1,i2,iVec)
!$omp& shared(integral,npq,nTransCho)
      nthreads=omp_get_num_threads()
      mythread=omp_get_thread_num()
      do ipq=1,npq
        if(mod(ipq,nthreads).eq.mythread) then
          i=ndx(ipq)
          do irs=ipq,npq
            i=i+1
            integral(i)=0
            iOffset=0
            do iVec=1,nTransCho
              do j=1,length(iVec)
                i1=(ipq-1)*numCho(1)+iOffset+j
                i2=(irs-1)*numCho(1)+iOffset+j
                integral(i)=integral(i)+coeff(i1)*coeff(i2)
              enddo
              iOffset=iOffset+length(iVec)
            enddo
          enddo
        endif
      enddo
!$omp end parallel
      deallocate(ndx)
      i=1
      iBuf=0
      do while(numint.gt.0)
        j=min(bufLength,numint)
        call ddafile(luTra,1,integral(i),j,iad50)
        i=i+j
        numint=numint-j
        iBuf=iBuf+1
        if(nBuf.ge.10) then
          if(mod(iBuf,int(nBuf/10)).eq.0 ) then
            write(*,'(A,I8,A,I8)') 'Writing buffer ',iBuf,
     &           ' of ',nBuf
          endif
        endif
      enddo
        write (*,'(A,I4,A)') 'Writing the last buffer with ',
     &                          j,' integrals'
      deallocate(integral)

      goto 9999
      allocate(integral(bufLength) )
      i = 0
      iBuf = 0
      do ipq = 1, npq
        do irs = ipq, npq
          i = i + 1
          integral(i) = 0
          iOffset = 0
          do iVec = 1, nTransCho
            do j = 1, length(iVec)
              i1 = (ipq - 1) * numCho(1) + iOffset + j
              i2 = (irs - 1) * numCho(1) + iOffset + j
              integral(i) = integral(i) + coeff(i1) * coeff(i2)
            end do
            iOffset = iOffset + length(iVec)
          end do
* Dump the integrals on TRAINT when the buffer is full and reset the counter (i)
          if ( i .eq. bufLength ) then
            i = 0
            iBuf = iBuf + 1
            if (nBuf.ge.10) then
              if ( mod(iBuf,int(nBuf/10)) .eq. 0 ) then
                write(*,'(A,I8,A,I8)') 'writing buffer ',iBuf,
     &                               ' of ',nBuf
              end if
            end if
*            if ( iBuf .eq. 1 ) then
*              write(*,*) '2-el integrals of the first buffer'
*              do j = 1, bufLength
*                write(*,'(i8,F25.15)')j,integral(j)
*              end do
*            end if
            call ddafile(luTra,1,integral,bufLength,iad50)
          end if
        end do
      end do
* Last (incomplete) buffer, skip it in the rare case that the number of integrals is a multiple of the buffer size
      if ( i .ne. 0) then
        write (*,'(A,I4,A)') 'Writing the last buffer with ',
     &                          i,' integrals'
*        do j = 1, i
*           write(*,'(i8,F25.15)')j,integral(j)
*        end do
       call ddafile(luTra,1,integral,bufLength,iad50)
      end if
      deallocate( integral )
 9999 continue

      call daclos(luTra)
      deallocate( length )
      deallocate( coeff )
c      deallocate( integral )
      end program read_cholesky
