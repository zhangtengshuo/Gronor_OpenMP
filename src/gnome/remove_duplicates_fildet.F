!> @brief   Removes the duplicates from the raw list of determinants
!!          and eliminates the determinant with coefficient smaller than
!!          thresh_CI
!! @author  Aitor Sanchez, URV
!! @author  Coen de Graaf, URV
!! @date    October 2020
!! @param   ndetin number of raw determinants
!! @param   ndetout number of unique determinants
!! @param   nact number of active orbitals
!! @param   cicoef CI coefficients
!! @param   occ string with the orbital occupation
!! @param   ndetaux number of unique determinants before filtering
!!          against thresh_CI
!! @param   occaux auxiliary occupation string

      subroutine remove_duplicates_fildet
     &     (ndetin,ndetout,nact,cicoef,occ,thresh_CI)
      
      implicit none
      
      integer      :: ndetin,ndetout,nact
      integer      :: ndetaux
      integer      :: i,j
      real(kind=8) :: cicoef(ndetin)
      real(kind=8) :: ciaux(ndetin)
      real(kind=8) :: thresh_CI
      character(len=255)    :: occ(ndetin)
      character(len=255)    :: occaux(ndetin)      
      
      do i=1,ndetin
         ciaux(i)=0.0d0
         occaux(i)=''
      enddo
      ndetout=1
      ciaux(1)=cicoef(1)
      occaux(1)=occ(1)
      outer : do i=2,ndetin
         do j=1,ndetout
            if(occaux(j).eq.occ(i))then
               ciaux(j)=ciaux(j)+cicoef(i)               
               cycle outer
            endif
         enddo
         ndetout=ndetout+1
         ciaux(ndetout)=cicoef(i)
         occaux(ndetout)=occ(i)
      enddo outer
 
      do i=1,ndetout   
         cicoef(i)=0.0d0
         occ(i)=''
      enddo
      ndetaux=ndetout
      ndetout=0      
      do i=1,ndetaux
         if(dabs(ciaux(i)).ge.thresh_CI)then
            ndetout=ndetout+1
            cicoef(ndetout)=ciaux(i)
            occ(ndetout)=occaux(i)
         endif
      enddo            
      return
      end
