program Effective_Hamiltonian
implicit real*8 (a-h,o-z)
parameter (nmax=40)
character*8 e_unit
dimension Etot(nmax),proj(nmax,nmax),E(nmax),smat(nmax,nmax),       &
          w(nmax),work(2*nmax),rwork(3*nmax),a(nmax,nmax),          &
          a2(nmax,nmax),rint(nmax,nmax),rdiag(nmax,nmax),           &
          r12inv(nmax,nmax),rs12inv(nmax,nmax),rs12int(nmax,nmax),  &
          pket(nmax,nmax),pbra(nmax,nmax),Heff(nmax,nmax),ipiv(nmax)
open(12,file='heff_info')
read(12,*)ndet,nstates
read(12,*)e_unit
iunit = 0
if (adjustl(e_unit(1:2)) .eq. 'cm')  iunit = 1
if (adjustl(e_unit(1:2)) .eq. 'eV')  iunit = 2
if (adjustl(e_unit(1:3)) .eq. 'meV') iunit = 3
if (iunit .eq. 0) then
  write(*,*) 'unit not implemented'
  write(*,*) 'please use cm, eV or meV'
  write(*,*) '   (case sensitive)     '
  stop
endif
do iss=1,nstates
  read(12,*) Etot(iss)
end do
do iss=1,nstates
  do jss=1,ndet
    read(12,*) proj(iss,jss)
  end do
end do
do iss=1,nstates
  write(*,*)
  write(*,'(A,I4)') 'State',iss
  if (iunit .eq. 1) then
    E(iss)=(etot(iss)-etot(1))*219474.631
    write(*,'(A,F10.4,A)') 'Energy:',E(iss),' cm-1'
  elseif (iunit .eq. 2) then
    E(iss)=(etot(iss)-etot(1))*27.2114
    write(*,'(A,F10.2,A)') 'Energy:',E(iss),' eV'
  elseif (iunit .eq. 3) then
    E(iss)=(etot(iss)-etot(1))*27211.4
    write(*,'(A,F12.1,A)') 'Energy:',E(iss),' meV'
  else
    write(*,*) 'not implemented yet'
    stop
  endif
  write(*,*) 'Projection on the model space'
  do jss=1,ndet
    write(*,100) '|',jss,'>',proj(iss,jss)
  end do
  do kss=1,nstates
    smat(iss,kss)=0.0
  end do
end do
100 format(A,I2,A,F12.6)

do iss=1,nstates
  do jss=1,nstates
    do kss=1,ndet
      smat(iss,jss)=smat(iss,jss)+proj(iss,kss)*proj(jss,kss)
    end do
  end do
end do
write(*,*)
write(*,*) 'OVERLAP MATRIX'
do iss=1,nstates
  write(*,101)(smat(iss,jss),jss=1,nstates)
  do jss=1,nstates
    a(iss,jss)=smat(iss,jss)
  end do
  w(iss)=0.0
  rwork(iss)=0.0
  rwork(iss+ndet)=0.0
  rwork(iss+2*ndet)=0.0
end do
101 format(24(F12.8))
 
call dsyev('V','U',nstates,a,nmax,w,rwork,3*nstates,info)
write (*,*) 
!write (*,*)'   Eigenvalues : Eigenvectors'
do i=1,nstates
!  write(*,102) i,w(i),(a(i,j),j=1,nstates)           ! DEBUG
  do j=1,nstates
    rint(i,j)=0.0
    rdiag(i,j)=0.0
    r12inv(i,j)=0.0
    rs12inv(i,j)=0.0
    rs12int(i,j)=0.0
  end do
  do k=1,ndet
    pket(i,k)=0.0
  end do
end do
!102 format(I4,2x,F9.6,' : ',24(F9.6,4x))
do i=1,nstates
  do j=1,nstates
    do k=1,nstates
      rint(i,j)=rint(i,j)+smat(i,k)*a(k,j)
    end do
    a2(i,j)=a(i,j)
  end do
  ipiv(i)=0
  work(i)=0.0
  work(i+nstates)=0.0
end do
call dgetrf(nstates,nstates,a2,nmax,ipiv,info)
call dgetri(nstates,a2,nmax,ipiv,work,2*nstates,info)
!do i=1,nstates
!  write(*,103) i,(a2(i,j),j=1,nstates)           ! DEBUG
!end do
!103 format(I4,3x,24(F9.6,12x))
do i=1,nstates
  do j=1,nstates
    do k=1,nstates
      rdiag(i,j)=rdiag(i,j)+a2(i,k)*rint(k,j)
    end do
  end do
end do
do i=1,nstates
!  write(*,103) i,(rdiag(i,j),j=1,nstates)          ! DEBUG
  r12inv(i,i)=1/sqrt(rdiag(i,i))
end do
do i=1,nstates
!  write(*,103) i,(r12inv(i,j),j=1,nstates)            ! DEBUG
end do

do i=1,nstates
  do j=1,nstates
    do k=1,nstates
      rs12int(i,j)=rs12int(i,j)+r12inv(i,k)*a2(k,j)
    end do
  end do
end do
do i=1,nstates
  do j=1,nstates
    do k=1,nstates
      rs12inv(i,j)=rs12inv(i,j)+a(i,k)*rs12int(k,j)
    end do
  end do
end do
!do i=1,nstates
!  write(*,103) i,(rs12inv(i,j),j=1,nstates)        ! DEBUG
!end do
do i=1,nstates
  do j=1,ndet
    do k=1,nstates
      pket(i,j)=pket(i,j)+rs12inv(k,i)*proj(k,j)
    end do
    pbra(i,j)=pket(i,j)
    HEff(i,j)=0.0
  end do
end do
write(*,*)
write(*,*)'|P> (orthonormalized projections)'
write(*,'(10x,24(A,I2,A,6x))')('|',i,'> ',i=1,ndet)
do i=1,nstates
  write(*,104) 'Psi',i,(pket(i,j),j=1,ndet)
end do
104 format(A,I2,24(F11.6))
do i=1,ndet
  do j=1,ndet
    do k=1,nstates
      Heff(i,j)=Heff(i,j)+pket(k,i)*E(k)*pbra(k,j)
    end do
  end do
end do
write(*,*)
if (iunit .eq. 1) then
  write(*,*)'Effective Hamiltonian (in cm-1)'
  write(*,'(14x,24(A,I2,A,7x))') ('|',i,'> ',i=1,ndet)
  do i=1,ndet
    write(*,105) '<',i,'|  ',(Heff(i,j),j=1,ndet)
  end do
elseif (iunit .eq. 2) then
  write(*,*)
  write(*,*)'Effective Hamiltonian (in eV)'
  write(*,'(14x,24(A,I2,A,7x))') ('|',i,'> ',i=1,ndet)
  do i=1,ndet
    write(*,106) '<',i,'|  ',(Heff(i,j),j=1,ndet)
  end do
else
  write(*,*)
  write(*,*)'Effective Hamiltonian (in meV)'
  write(*,'(14x,24(A,I2,A,7x))') ('|',i,'> ',i=1,ndet)
  do i=1,ndet
    write(*,107) '<',i,'|  ',(Heff(i,j),j=1,ndet)
  end do
endif
105 format(A,I2,A,24(F12.1))
106 format(A,I2,A,24(F12.3))
107 format(A,I2,A,24(F12.1))
end program Effective_Hamiltonian


