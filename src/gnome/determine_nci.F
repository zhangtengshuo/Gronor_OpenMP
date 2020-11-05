!> @brief   Update the number of determinants for a monomer function
!!          by the action of the S- operator 
!! @author  Coen de Graaf, URV
!! @date    November 2020

      subroutine determine_nci

      use cidef

      implicit none

      integer     :: i,j,k,pascal,ms
      integer     :: state,det,act,alpha

      integer,allocatable  :: spin(:)
      integer,allocatable  :: ndets(:)

      
      allocate( spin(mstates) )
      do state = 1, mstates
        spin(state) = 0                          ! spin of the monomer function
        do act = 1, nactm(state)
          if ( abs(ioccm(act,1,state)) .eq. 1 ) then
            spin(state) = spin(state) + ioccm(act,1,state)
          end if
        end do
      end do
      allocate( ndets(mstates) )

      nci = 0
      do state = 1, mstates
        ndets(state) = 0
        if ( spin(state) .eq. 0 ) then
          ndets(state) = idetm(state)
        else
          do det = 1, idetm(state)
            alpha = 0
            do act = 1, nactm(state)
              if ( ioccm(act,det,state) .eq. 1 ) alpha = alpha + 1
            end do
            do ms = 1, spin(state) + 1               ! from ms_max to -ms_max
              ndets(state) = ndets(state) + pascal(alpha+1,ms)
            end do
          end do
        end if
        nci = max(nci,ndets(state))
      end do
      
      deallocate(spin)
      deallocate(ndets)
      return
      end subroutine determine_nci
