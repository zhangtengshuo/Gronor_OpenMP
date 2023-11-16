! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !
!                                                               !
!  Program to generate natural transiton orbitals and prepare   !
!  the files for AIFD calculations with GronOR                  !
!                                                               !
!                                                               !
!  written by Aitor Sanchez-Mansilla, URV (2020)                !
!                                                               !
!                                                               !
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !

module nto_data
implicit none
integer                        :: nFrags
integer                        :: nras1,nras3
integer                        :: iRC,istat,lwork,ldiag
integer                        :: nnto,spin
integer,allocatable            :: nStates(:)
integer,allocatable            :: nBasFrag(:),nOrbFrag(:)
integer,allocatable            :: nInactFrag(:),nActFrag(:)
integer,allocatable            :: nDetFrag(:)
character(len=1),allocatable   :: orblabel(:)
character(len=4),allocatable   :: fragName(:),fragState(:)
character(len=25)              :: gsfile,excfile,detfile
character(len=25)              :: ntofile, ntodetfile, gsdetfile
character(len=25)              :: runfile,oneintfile
character(len=132)             :: Project
character(len=255),allocatable :: detocc(:),fragLabel(:)
real(kind=8)                   :: thresh_CI,thresh_NTO
real(kind=8)                   :: norm
real(kind=8),allocatable       :: gs(:,:),exc(:,:),occNu(:)
real(kind=8),allocatable       :: hole(:,:),part(:,:)
real(kind=8),allocatable       :: detcoef(:)
real(kind=8),allocatable       :: sAO(:,:),sHP(:,:)
real(kind=8),allocatable       :: u(:,:),vt(:,:),work(:),sval(:)
logical                        :: debug
end module

! ======================================================================= !

program nto_transformation
use nto_data
implicit none
integer              :: iFrag,iState

call nto_readin
allocate(nBasFrag(nFrags))
allocate(nOrbFrag(nFrags))
allocate(nInactFrag(nFrags))
allocate(nActFrag(nFrags))
allocate(nDetFrag(nFrags))
do iFrag = 1,nFrags
  write(*,*)
  write(*,*) 'Processing fragment ',fragName(iFrag)
  call nto_getnBas(iFrag)
  allocate(gs(nBasFrag(iFrag),nBasFrag(iFrag)))
  allocate(exc(nBasFrag(iFrag),nBasFrag(iFrag)))
  allocate(sAO(nBasFrag(iFrag),nBasFrag(iFrag)))
  allocate(occNu(nBasFrag(iFrag)))
  allocate(orblabel(nBasFrag(iFrag)))
  do iState = 2,nStates(iFrag)
  write(*,'(A,I4)') '   State ', iState
    gs = 0.0
    exc = 0.0
    sAO = 0.0
    occNu = 0
    call nto_getFilenames(iFrag,iState)
    if (iState.eq.2) call nto_getAtomicOverlap(iFrag)
    call nto_read_vecs(iFrag)
    call nto_read_dets(iFrag)
    if (debug) write(6,'(A)')'Vectors and determinants have been read and stored'
    call nto_getNTO(iFrag,iState)
    call nto_write_vecs(iFrag)
    call nto_write_dets(iFrag)
    write(*,'(A)') 'Natural transition orbitals have been calculated'
    write(*,'(A,I4,4X,I4)')'Original active space:',nRas1,nRas3
    write(*,'(A,I4)')'Number of significant NTO pairs:',nnto
    write(*,*)
    deallocate(hole,part)
  enddo
  deallocate(gs,exc,sAO,occNu,orblabel)
enddo
deallocate(nStates,nBasFrag,nOrbFrag)
deallocate(nInactFrag,nActFrag,nDetFrag)
deallocate(fragName,fragState,fragLabel)

end program nto_transformation

! ======================================================================= !

subroutine nto_getNTO(iFrag,iState)
use nto_data
implicit none
integer                  :: iFrag,iState
integer                  :: i,j,k,l,m,n
integer                  :: iOcc,iVir,ndiff
real(kind=8)             :: dum,overlap,mxOverlap
real(kind=8),allocatable :: tmatrix(:,:)
real(kind=8),allocatable :: aux(:,:)

character(len=2550) :: occStr

allocate(tmatrix(nRas1,nRas3))
allocate(hole(nRas1,nBasFrag(iFrag)),part(nRas3,nBasFrag(iFrag)) )
tmatrix = 0.0
if (debug ) then
  write(6,'(A,I2,A,I5)')'Number of basis fns of fragment ',iFrag,':',nBasFrag(iFrag)
  write(6,'(A,I3)')'RAS1 (occupied) orbitals:',nRas1
  write(6,'(A,I3)')'RAS3 (virtual)  orbitals:',nRas3
endif

occStr = detocc(1)
spin=1
do i=1,nActFrag(iFrag)
  if(occstr(i:i).eq.'a')spin=spin+1
  if(occstr(i:i).eq.'b')spin=spin-1
enddo
if (debug) then
  write(6,'(A,I1)')'State:',iState
  write(6,'(A,I1)')'Spin:',spin
endif

do i=1,nDetFrag(iFrag)
  if (debug) write(6,'(A,I5,A,I5)')'Determinant',i,' of ',nDetFrag(iFrag)
  iOcc=0
  iVir=0
  occStr = detocc(i)
  ndiff=0
  if (debug) write(6,'(2A)')'occstring: ',trim(occStr(1:nActFrag(iFrag)))
  if (spin.eq.1) then
    do j=1,nRas1
      if (occStr(j:j).eq.'b') then
        iOcc=j
        ndiff=ndiff+1
      endif
      if (occStr(j:j).eq.'a') then
        ndiff=ndiff+1
      endif
    enddo
    do j=nRas1+1,nRas1+nRas3
      if (occStr(j:j).eq.'a') then
        iVir=j
        ndiff=ndiff+1
      endif
      if (occStr(j:j).eq.'b') then
        ndiff=ndiff+1
      endif
    enddo
  else if (spin.eq.3) then
    do j=1,nRas1
      if (occStr(j:j).eq.'a') then
        iOcc=j
        ndiff=ndiff+1
      endif
    enddo
    do j=nRas1+1,nRas1+nRas3
      if (occStr(j:j).eq.'a') then
        iVir=j
        ndiff=ndiff+1
      endif
    enddo
    do j=1,nRas1+nRas3
      if (occStr(j:j).eq.'b') ndiff=ndiff+1
    enddo
  else
    write(6,'(A)')'spin configuration cannot be described with single excitations'
    write(6,'(A,I1)')'wrong spin:',spin
  endif
  if (iOcc.ne.0 .and. iVir.ne.0 .and. ndiff.eq.2) then
    if (debug) write(6,'(A,I3,A,I3)')'Excitation',iOcc,'  to',iVir
    if (debug) write(6,*)'with coefficient',detcoef(i)
    tmatrix(iOcc,iVir-nRas1) = tmatrix(iOcc,iVir-nRas1)+detcoef(i)
  else
    if (debug) then
      write(6,'(A)')'No single excitation'
    endif
  endif
  if (debug) write(6,*)'----------------------------------------------'
end do
if (debug) then
  write(6,'(A,I2,1x,I2)')'1-TDM FOR FRAGMENT:',iFrag,iState
  do i=1,nRas1
    write(6,'(200F7.3)')(tmatrix(i,j),j=1,nRas3)
  enddo
end if
if (nRas3.gt.nRas1) then
  ldiag = nRas1
  m = nRas3
  n = nRas1
  allocate(aux(nRas3,nRas1))
  do i=1,nRas3
    do j=1,nRas1
      aux(i,j)=tmatrix(j,i)
    enddo
  enddo
else
  ldiag = nRas3
  m = nRas1
  n = nRas3
  allocate(aux(nRas1,nRas3))
  do i=1,nRas1
    do j=1,nRas3
      aux(i,j)=tmatrix(i,j)
    enddo
  enddo
endif
allocate(u(m,m),vt(n,n))
allocate(work(5*max(m,n)),sval(ldiag))
u = 0.0
vt = 0.0
work = 0.0
sval = 0.0
call dgesvd('A','A',m,n,aux,max(m,n),sval,u,m,vt,n,work,5*max(m,n),iRC)
if (iRC.eq.0) then
  if (debug) write(6,*)'SVD successful exit',iRC
else if (iRC.lt.0) then
  write(6,*)'SVD error parameter:',abs(iRC)
else if (iRC.gt.0) then
  write(6,*)'SVD did not converge:',iRC
endif
deallocate(aux)
if (debug) then
  write(6,'(A)')'U MATRIX:'
  do i=1,m
    write(6,'(200F7.3)')u(i,:)
  enddo
  write(6,'(A)')'VT MATRIX:'
  do i=1,n
    write(6,'(200F7.3)')vt(i,:)
  enddo
  write(6,'(A)')'SINGULAR VALUES (EXCITATION AMPLITUDES):'
  write(6,'(12F14.8)')sval(:)
  dum = 0.0d0
  do i=1, ldiag
    dum = dum + sval(i)**2
  enddo
  write(6,'(A)')'SUM OF SQUARED AMPLITUDES:'
  write(6,'(F14.8)')dum
endif
nnto=0
do i=1,ldiag
  if (sval(i).ge.thresh_NTO) nnto=nnto+1
enddo

if (debug) write(6,'(A)')'CALCULATE NTOS'
hole=0.0
part=0.0
if (nRas1.ge.nRas3) then
  do i=1,m
    do j=1,nBasFrag(iFrag)
      do k=1,m
        hole(i,j)=hole(i,j)+u(k,i) * gs(k+nInactFrag(iFrag),j)
      enddo
    enddo
  enddo
else
  do i=1,m
    do j=1,nBasFrag(iFrag)
      do k=1,m
        part(i,j)=part(i,j)+u(k,i) * gs(k+nInactFrag(iFrag)+n,j)
      enddo
    enddo
  enddo
endif
if (nRas1.ge.nRas3) then
  do i=1,n
    do j=1,nBasFrag(iFrag)
      do k=1,n
        part(i,j)=part(i,j)+vt(i,k)*gs(k+nInactFrag(iFrag)+m,j)
      enddo
    enddo
  enddo
else
  do i=1,n
    do j=1,nBasFrag(iFrag)
      do k=1,n
        hole(i,j)=hole(i,j)+vt(i,k)*gs(k+nInactFrag(iFrag),j)
      enddo
    enddo
  enddo
endif
deallocate(u,vt,work)

if (debug) then
  do i=1,nRas1
    write(6,'(A,I3)')'Hole orbital',i
    write(6,'(10F15.8)')(hole(i,j),j=1,nBasFrag(iFrag))
  enddo
  do i=1,nRas3
    write(6,'(A,I3)')'particle orbital',i
    write(6,'(10F15.8)')(part(i,j),j=1,nBasFrag(iFrag))
  enddo
  allocate(sHP(nRas1,nRas3))
  sHP = 0.0
  do i=1,nRas1
    do j=1,nRas3
      do k=1,nBasFrag(iFrag)
        do l=1,nBasFrag(iFrag)
          sHP(i,j)=sHP(i,j)+hole(i,k)*part(j,l)*sAO(k,l)
        enddo
      enddo
    enddo
  enddo
  write(6,'(A)')'Hole-Particle overlap:'
  do i=1,nRas1
    write(6,'(12F14.8)')(sHP(i,j),j=1,nRas3)
  enddo
  deallocate(sHP)
  mxOverlap = 0.0d0
  do i = 1, nInactFrag(iFrag)
    do j = 1, nRas1
      overlap = 0.0d0
      do k = 1, nBasFrag(iFrag)
        do l = 1, nBasFrag(iFrag)
          overlap = overlap + gs(i,k) * hole(j,l) * sAO(k,l)
        end do
      end do
      if (abs(overlap) .gt. abs(mxOverlap)) mxOverlap = overlap
    end do
  end do
  do i = 1, nInactFrag(iFrag)
    do j = 1, nRas3
      overlap = 0.0d0
      do k = 1, nBasFrag(iFrag)
        do l = 1, nBasFrag(iFrag)
          overlap = overlap + gs(i,k) * part(j,l) * sAO(k,l)
        end do
      end do
      if (abs(overlap) .gt. abs(mxOverlap)) mxOverlap = overlap
    end do
  end do
  write(*,'(A,F20.10)') 'Max overlap between inactive and nto ',mxOverlap
endif

deallocate(tmatrix,detcoef,detocc)

end subroutine nto_getNTO

! ======================================================================= !

subroutine nto_write_vecs(iFrag)
use nto_data
implicit none

integer   :: iFrag,i,iorb,j,k

if (debug) write(6,'(2A)')'Write NTOs to file: ',ntofile
open(12,file=ntofile)
write(12,'(A11)')'#INPORB 2.2'
write(12,'(A5)')'#INFO'
write(12,'(A)')'* CIS excited state with active orbitals replaced by NTOs'
write(12,'(3I8)')0,1,0
write(12,'(I8)') nBasFrag(iFrag)
write(12,'(I8)') nBasFrag(iFrag)
write(12,'(A4)') '#ORB'
!     write occupied inactive orbitals
do i=1,nInactFrag(iFrag)
  write(12,'(A9,2I5)')'* ORBITAL',i
  write(12,'(5E22.14)')(gs(i,j),j=1,nBasFrag(iFrag))
enddo
!     write hole orbitals with sigma below the threshold
do i=1,nRas1-nnto
  write(12,'(A9,2I5)')'* ORBITAL',i+nInactFrag(iFrag)
  write(12,'(5E22.14)')(hole(i+nnto,j),j=1,nBasFrag(iFrag))
enddo
!     write hole orbitals significant for the excitation
do i=1,nnto
  write(12,'(A9,2I5)')'* ORBITAL',i+nInactFrag(iFrag)+nRas1-nnto
  write(12,'(5E22.14)')(hole(i,j),j=1,nBasFrag(iFrag))
enddo
!     write particle orbitals significant for the excitation
do i=1,nnto
  write(12,'(A9,2I5)')'* ORBITAL',i+nInactFrag(iFrag)+nRas1
  write(12,'(5E22.14)')(part(i,j),j=1,nBasFrag(iFrag))
enddo
!     write particle orbitals with sigma below the threshold
do i=1,nRas3-nnto
  write(12,'(A9,2I5)')'* ORBITAL',i+nInactFrag(iFrag)+nRas1+nnto
  write(12,'(5E22.14)')(part(i+nnto,j),j=1,nBasFrag(iFrag))
enddo
!     write virtual orbitals
do i=nInactFrag(iFrag)+nActFrag(iFrag)+1,nBasFrag(iFrag)
  write(12,'(A9,2I5)')'* ORBITAL',i
  write(12,'(5E22.14)')(gs(i,j),j=1,nBasFrag(iFrag))
enddo

write(12,'(A4)') '#OCC'
write(12,'(A20)') '* OCCUPATION NUMBERS'
write(12,'(5E22.14)') (occNu(j),j=1,nBasFrag(iFrag))

write(12,'(A)')'#INDEX'
write(12,'(A)')'* 1234567890'
iorb = 0
orblabel = ''
do i=1,nInactFrag(iFrag)
  iorb = iorb + 1
  orblabel(iorb) = 'i'
enddo
do i=1,nRas1-nnto
  iorb = iorb + 1
  orblabel(iorb) = 'i'
enddo
do i=1,nnto
  iorb = iorb + 1
  orblabel(iorb) = '1'
enddo
do i=1,nnto
  iorb = iorb + 1
  orblabel(iorb) = '3'
enddo
do i=1,nRas3-nnto
  iorb = iorb + 1
  orblabel(iorb) = 's'
enddo
do i=1,nOrbFrag(iFrag)-nActFrag(iFrag)-nInactFrag(iFrag)
  iorb = iorb + 1
  orblabel(iorb) = 's'
enddo
if (debug) write(*,'(10A1)')(orblabel(k),k=1,nBasFrag(iFrag))
do j = 1, nBasFrag(iFrag), 10
  if ( j + 9 .le. nBasFrag(iFrag)) then
    write(12,102)mod(int(j/10),10),(orbLabel(k),k=j,j+9)
  else
    write(12,102)mod(int(j/10),10),(orbLabel(k),k=j,nBasFrag(iFrag))
  end if
end do
102 format(I1,x,10A1)
close(unit=12)
if (debug) write(6,'(2A)')'CLOSE FILE:',ntofile

end subroutine nto_write_vecs

! ======================================================================= !

subroutine nto_write_dets(iFrag)
use nto_data
implicit none

integer             :: iFrag,i,j
integer             :: detheader(8)
character(len=2550) :: occStr

if (debug) write(*,'(2A)')'write NTO detlist to file: ',ntodetfile
open(13,file=ntodetfile)
detheader = 0
detheader(1) = nInactFrag(iFrag) + nRas1 - nnto
write(13,'(8(I4))') (detheader(i),i=1,8)

do i=1,nnto
  if (spin.eq.1) then
    norm=0.0d0
    do j=1,nnto
      norm=norm+2*sval(j)**2
    enddo
    norm=dsqrt(norm)
    occStr=''
    do j=1,nnto
      if(j.eq.i)then
        occStr = trim(occStr)//'b'
      else
        occStr = trim(occStr)//'2'
      endif
    enddo
    do j=1,nnto
      if(j.eq.i)then
        occStr = trim(occStr)//'a'
      else
        occStr = trim(occStr)//'0'
      endif
    enddo
    write(13,'(e15.8,6x,A)')sval(i)/norm,trim(occStr)
    occStr=''
    do j=1,nnto
      if(j.eq.i)then
        occStr = trim(occStr)//'a'
      else
        occStr = trim(occStr)//'2'
      endif
    enddo
    do j=1,nnto
      if(j.eq.i)then
        occStr = trim(occStr)//'b'
      else
        occStr = trim(occStr)//'0'
      endif
    enddo
    write(13,'(e15.8,6x,A)')sval(i)*(-1)/norm,trim(occStr)
  else if(spin.eq.3)then
    norm=0.0d0
    do j=1,nnto
      norm=norm+sval(j)**2
    enddo
    norm=dsqrt(norm)
    occStr=''
    do j=1,nnto
      if(j.eq.i)then
        occStr = trim(occStr)//'a'
      else
        occStr = trim(occStr)//'2'
      endif
    enddo
    do j=1,nnto
      if(j.eq.i)then
        occStr = trim(occStr)//'a'
      else
        occStr = trim(occStr)//'0'
      endif
    enddo
    write(13,'(e15.8,6x,A)')sval(i)/norm,trim(occStr)
  else
    write(6,'(A)')'spin configuration cannot be described with single excitations'
    write(6,'(A,I1)')'wrong spin:',spin
  endif
enddo
deallocate(sval)
close(13)

if (debug) write(6,'(2A)')'Ground state detfile : ',gsdetfile
open(14,file=gsdetfile)
write(14,'(8(I4))') (detheader(i),i=1,8)
occStr = ''
do i=1,nnto
  occStr = trim(occStr)//'2'
enddo
do i=1,nnto
  occStr = trim(occStr)//'0'
enddo
write(14,'(e15.8,6x,A)')1.0d0,trim(occStr)
close(unit=14)

end subroutine nto_write_dets

! ======================================================================= !

subroutine nto_readin
use nto_data
implicit none

integer, parameter                 :: nKeys = 5
integer                            :: iKey,iFrag,j,jj,start
logical                            :: all_ok,fragNameOK
logical, dimension(nKeys)          :: hit = .false.
character (len=132)                :: line,string
character (len=4)                  :: key
character (len=4),dimension(nKeys) :: keyword

data keyword /'PROJ','FRAG','LABE','THRE','DEBU'/

! some defaults
Project='unknown'
thresh_CI = 0.0d0
thresh_NTO = 0.0d0

all_ok=.true.
do while (all_ok)
  read(5,*,iostat=jj) line
  line = adjustl(line)
  key = line(1:4)
  call nto_capitalize(key)
  do iKey = 1, nKeys
    if ( key .eq. keyword(iKey) ) hit(iKey) = .true.
  end do
  if (  jj .lt. 0 ) all_ok = .false.
end do

do iKey = 1, nKeys
  if (hit(iKey)) then
    select case(iKey)
      case(1)
        call nto_locate('PROJ')
        read(*,*) Project
        Project = trim(Project)
      case(2)
        call nto_locate('FRAG')
        read(*,*) nFrags
        allocate(nStates(nFrags))
        read(*,*) (nStates(iFrag),iFrag=1,nFrags)
      case(3)
        call nto_locate('LABE')
        allocate(fragName(nFrags))
        allocate(fragLabel(sum(nStates)))
        allocate(fragState(sum(nStates)))
        fragName(:) =''
        fragLabel(:)=''
        fragState(:)=''
        start = 0
        do j = 1, nFrags
          start = start + 1
          read(*,'(A)') line
          string = adjustl(line)
          fragNameOK = .false.
          do jj = 1, len(trim(string))
            if (string(jj:jj) .ne. ' ') then
              if (.not. fragNameOK) then
                write(fragName(j),'(a,a)') trim(adjustl(fragName(j))),string(jj:jj)
              else
                write(fragState(start),'(a,a)') trim(adjustl(fragState(start))),string(jj:jj)
              endif
            else
              if (fragNameOK) then
                if (string(jj-1:jj-1) .ne. ' ') start = start + 1
              else
                fragNameOK = .true.
              endif
            endif
          enddo
        enddo
        start = 1
        do j = 1, nFrags
          do jj = start, start+nStates(j)-1
            write(fragLabel(jj),'(a,a,a)') trim(fragName(j)),'_',trim(fragState(jj))
          enddo
          start = start + nStates(j)
        enddo
      case(4)
        call nto_locate('THRE')
        read(*,*) thresh_NTO,thresh_CI
      case(5)
        debug = .true.
    end select
  endif
enddo

end subroutine nto_readin

! ======================================================================= !

subroutine nto_locate(string)
implicit none
external :: nto_capitalize
character(4)   ::  string,string2
character(132) ::  line
rewind(5)
40 read(5,*) line
line = adjustl(line)
string2=line(1:4)
call nto_capitalize(string2)
if (string2.ne.string) goto 40
return
end subroutine nto_locate

! ======================================================================= !

subroutine nto_capitalize(string)
implicit none
integer      :: i
character(*) string

do i = 1, len(string)
  if (ichar(string(i:i)).gt.96) then
    string(i:i) = char(ichar(string(i:i))-32)
  endif
end do
return
end subroutine nto_capitalize

! ======================================================================= !

subroutine nto_getnBas(iFrag)
use nto_data
implicit none
integer       :: iFrag,nSym
write(runfile,'(A6,A1)')'RUNFIL',fragName(iFrag)
call NameRun(runfile)
call Get_iScalar('nSym',nSym)
if ( nSym .ne. 1 ) then
  write(*,*) '  Symmetry is not implemented in GronOR'
  write(*,*) 'Remove the symmetry elements and come back'
  stop
end if
call Get_iArray('nBas',nBasFrag(iFrag),nSym)
end subroutine nto_getnBas

! ======================================================================= !

subroutine nto_getFilenames(iFrag,iState)
use nto_data
implicit none
integer              :: iFrag,iState,start
character (len=3)    :: suffix
start = (iFrag-1)*nStates(iFrag) 

suffix = 'orb'
write(gsfile,'(4A)')  trim(project),trim(fragLabel(start+1)),'.',suffix
write(excfile,'(4A)') trim(project),trim(fragLabel(start+iState)),'_org.',suffix
write(ntofile,'(4A)') trim(project),trim(fragLabel(start+iState)),'.',suffix

suffix = 'det'
write(detfile,'(4A)') trim(project),trim(fragLabel(start+iState)),'_org.',suffix
write(gsdetfile,'(4A)') trim(project),trim(fragLabel(start+1)),'.',suffix
write(ntodetfile,'(4A)') trim(project),trim(fragLabel(start+iState)),'.',suffix

write(runfile,'(2A)') 'RUNFIL',fragName(iFrag)
write(oneintfile,'(2A)') 'ONEINT',fragName(iFrag)

if (debug) then
  write(*,*) 'Filenames....'
  write(*,'(A)')trim(gsfile),trim(excfile),trim(ntofile),     &
                  trim(detfile),trim(gsdetfile),trim(ntodetfile)
endif

return
end subroutine nto_getFilenames

! ======================================================================= !

subroutine nto_getAtomicOverlap(iFrag)
use nto_data
implicit none
integer                           :: iFrag
integer                           :: iCounter,iComponent,LuOne
integer                           :: iOpt,iSymLbl,j,k
integer                           :: lTriangle
real (kind=8),allocatable         :: s(:)

call NameRun(runfile)
lTriangle = (nBasFrag(iFrag)*(nBasFrag(iFrag)+1))/2
allocate(s(lTriangle))
s = 0.0
sAO = 0.0
iRc=-1
iOpt=0
luOne = 77
Call OpnOne(iRC,iOpt,oneintfile,LuOne)
if (iRC.ne.0) write(6,*)'Something went wrong opening',oneintfile
iRC =  0
iOpt = 2
iComponent = 1
iSymLbl = 1
Call RdOne(iRC,iOpt,'Mltpl  0',iComponent,s,iSymLbl)
iCounter = 1
sAO = 0.0
do j = 1, nBasFrag(iFrag)
  do k = 1, j
    sAO(j,k) = s(iCounter)
    sAO(k,j) = s(iCounter)
    iCounter = iCounter + 1
  end do
end do
iOpt = 0
Call ClsOne(iRc,iOpt)
deallocate(s)
return
end subroutine nto_getAtomicOverlap

! ======================================================================= !

subroutine nto_read_vecs(iFrag)
use nto_data
implicit none
integer                         :: iFrag
integer                         :: i,j
character(len=6)                :: mark
character(len=132)              :: line

! ground state
open(9,file=gsfile, status='old')

gs = 0.0
exc = 0.0
occNu = 0.0

mark  = '#INDEX'
46 read(9,'(A132)') line
if (line(1:6).ne.mark) goto 46
read(9,'(A132)') line
orbLabel = ' '
read(9,'(2x,10A)')(orbLabel(i),i=1,nBasFrag(iFrag))
rewind(9)

mark = '#ORB'
47 read(9,'(A132)') line
if (line(1:4).ne.mark) goto 47
do i = 1, nBasFrag(iFrag)
  read(9,'(A132)') line
  read(9,'(5E22.14)') (gs(i,j),j=1,nBasFrag(iFrag))
end do

close(unit=9)

! excited state
open(10,file=excfile,status='old')

nOrbFrag(iFrag) = 0
nInactFrag(iFrag) = 0
nActFrag(iFrag) = 0
nRas1 = 0
nRas3 = 0

mark = '#INDEX'
48 read(10,'(A132)') line
if (line(1:6).ne.mark) goto 48
read(10,'(A132)') line
orbLabel = ' '
read(10,'(2x,10A)')(orbLabel(i),i=1,nBasFrag(iFrag))
do i = 1, nBasFrag(iFrag)
  if (orbLabel(i) .eq. 'f') then
    nInactFrag(iFrag) = nInactFrag(iFrag) + 1
    nOrbFrag(iFrag) = nOrbFrag(iFrag) + 1
  endif
  if (orbLabel(i) .eq. 'i') then
    nInactFrag(iFrag) = nInactFrag(iFrag) + 1
    nOrbFrag(iFrag) = nOrbFrag(iFrag) + 1
  endif
  if (orbLabel(i) .eq. '1') then
    nActFrag(iFrag) = nActFrag(iFrag) + 1
    nOrbFrag(iFrag) = nOrbFrag(iFrag) + 1
    nRas1 = nRas1 + 1
  endif
  if (orbLabel(i) .eq. '2') then
    nActFrag(iFrag) = nActFrag(iFrag) + 1
    nOrbFrag(iFrag) = nOrbFrag(iFrag) + 1
  endif
  if (orbLabel(i) .eq. '3') then
    nActFrag(iFrag) = nActFrag(iFrag) + 1
    nOrbFrag(iFrag) = nOrbFrag(iFrag) + 1
    nRas3 = nRas3 +1
  endif
  if (orbLabel(i) .eq. 's') then
    nOrbFrag(iFrag) = nOrbFrag(iFrag) + 1
  endif
end do

rewind(10)
mark = '#ORB'
49 read(10,'(A132)') line
if (line(1:4).ne.mark) goto 49
do i = 1, nBasFrag(iFrag)
  read(10,'(A132)') line
  read(10,'(5E22.14)') (exc(i,j),j=1,nBasFrag(iFrag))
end do

rewind(10)
mark = '#OCC'
50 read(10,'(A132)') line
if (line(1:4).ne.mark) goto 50
read(10,*) line
read(10,'(5e22.14)')(occNu(j),j=1,nBasFrag(iFrag))

close(unit=10)

end subroutine nto_read_vecs

! ======================================================================= !

subroutine nto_read_dets(iFrag)
use nto_data
implicit none

integer                         :: ndet_raw,i,j,iFrag
real(kind=8),allocatable        :: detcoef_raw(:)
character(len=255),allocatable  :: detocc_raw(:)

open(11,file=detfile, status='old')
ndet_raw = 0

read(11,*)
do
  read(11,*,iostat=istat)
  if(istat.ne.0)then
    rewind(11)
    exit
  else
    ndet_raw = ndet_raw + 1
  endif
enddo
allocate(detcoef_raw(ndet_raw) )
allocate(detocc_raw(ndet_raw) )
read(11,*)
do i = 1,ndet_raw
  read(11,*) detcoef_raw(i),detocc_raw(i)
enddo
nDetFrag(iFrag) = 0
do i=1,ndet_raw
  if(abs(detcoef_raw(i)).ge.thresh_CI) then
    nDetFrag(iFrag) = nDetFrag(iFrag) + 1
  endif
enddo
allocate(detcoef(nDetFrag(iFrag)) )
allocate(detocc(nDetFrag(iFrag)) )
j=0
do i=1,ndet_raw
  if(abs(detcoef_raw(i)).ge.thresh_CI)then
    j = j+1
    detcoef(j)=detcoef_raw(i)
    detocc(j)=detocc_raw(i)
  endif
enddo
deallocate(detcoef_raw,detocc_raw)
if (ndet_raw.ne.nDetFrag(iFrag)) then
  write(6,'(3(a,i4))')'number of determinants reduced from',ndet_raw,' to',nDetFrag(iFrag),' for fragment',iFrag
endif
if (debug) then
  write(6,'(a)')'*DETERMINANT LIST*'
  do i = 1,nDetFrag(iFrag)
    write(6,*) i,detcoef(i),trim(detocc(i))
  enddo
endif
close(unit=11)

end subroutine nto_read_dets

! ======================================================================= !
