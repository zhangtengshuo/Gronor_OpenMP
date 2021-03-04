!     This file is part of the GronOR software

!     GronOR is free software, and can be used, re-distributed and/or modified under
!     the Apache License version 2.0 (http://www.apache.org/licenses/LICENSE-2.0)
!     Any use of the software has to be in compliance with this license. Unless required
!     by applicable law or agreed to in writing, software distributed under the license
!     is distributed on an ‘as is’ basis, without warranties or conditions of any kind,
!     either express or implied.
!     See the license for the specific language governing permissions and limitations
!     under the license.

!     GronOR is copyright of the University of Groningen
      
!> @brief
!! Open CML file and dump basic info
!!     
!! @author  T. P. Straatsma, ORNL
!! @author  C. de Graaf, URV / ICREA
!! @date    2020
!!
  

      subroutine gronor_prelude_cml()
      use cidef

      implicit none
      character(len=132)   :: info_cml

      write(lfncml,'(a38)')'<?xml version="1.0" encoding="UTF-8"?>'
      info_cml = '<module xmlns="http://www.xml-cml.org/schema"'
      write(lfncml,'(a)') trim(info_cml)
      info_cml = 
     &  '  xmlns:cc="http://www.xml-cml.org/dictionary/compchem/"'
      write(lfncml,'(a)') trim(info_cml)
      info_cml = '  xmlns:cml="http://www.xml-cml.org/schema"'
      write(lfncml,'(a)') trim(info_cml)
      info_cml = '  xmlns:cmlx="http://www.xml-cml.org/schema/cmlx"'
      write(lfncml,'(a)') trim(info_cml)
      info_cml =
     &  '  xmlns:convention="http://www.xml-cml.org/convention/"'
      write(lfncml,'(a)') trim(info_cml)
      info_cml = 
     &  '  xmlns:gr="http://www.iochem-bd.org/dictionary/gronor/"'
      write(lfncml,'(a)') trim(info_cml)
      info_cml = '  xmlns:nonsi="http://www.xml-cml.org/unit/nonSi/"'
      write(lfncml,'(a)') trim(info_cml)
      info_cml = 
     &  '  xmlns:nonsi2="http://www.iochem-bd.org/unit/nonSi2/"'
      write(lfncml,'(a)') trim(info_cml)
      info_cml = '  xmlns:si="http://www.xml-cml.org/unit/si/"'
      write(lfncml,'(a)') trim(info_cml)
      info_cml = '  xmlns:xi="http://www.w3.org/2001/XInclude"'
      write(lfncml,'(a)') trim(info_cml)
      info_cml = '  xmlns:xsd="http://www.w3.org/2001/XMLSchema"'
      write(lfncml,'(a)') trim(info_cml)
      info_cml = '  convention="convention:compchem>"'
      write(lfncml,'(a)') trim(info_cml)

      return

      end subroutine gronor_prelude_cml



      subroutine gronor_env_cml(version,version_type)
      use cidef
      use gnome_parameters
      use cidist

      implicit none
      character(len=5)    :: version
      character(len=18)   :: runDate
      character(len=20)   :: label
      character(len=64)   :: version_type
      character(len=132)  :: info_cml
      integer             :: indent
      integer,external    :: getcpucount

* Open the calculation
      label = 'module'
      indent = 1
      info_cml = 'id="job" dictRef="cc:job"'
      call open_tag(lfncml,label,info_cml,indent)
* Block 1: general job info_cml
      indent = 2
      info_cml = 'dictRef="cc:environment" id="environment"'
      call open_tag(lfncml,label,info_cml,indent)
      label = 'parameterList'
      indent = 3
      info_cml = 'empty'
      call open_tag(lfncml,label,info_cml,indent)
      indent = 4
      label = 'parameter'
      info_cml = 'dictRef="cc:program"'
      call open_tag(lfncml,label,info_cml,indent)
      indent = 5
      info_cml = 'empty'
      call writetag_scalar_string(lfncml,info_cml,indent,
     &             version)
      call close_tag(lfncml,label,4)
      info_cml = 'dictRef="cc:programVersion"'
      call open_tag(lfncml,label,info_cml,4)
      info_cml='empty'
      call writetag_scalar_string(lfncml,info_cml,5,
     &            version_type)
      call close_tag(lfncml,label,4)
      info_cml = 'dictRef="cc:executable"'
      call open_tag(lfncml,label,info_cml,4)
      info_cml='empty'
      call writetag_scalar_string(lfncml,info_cml,5,command)
      call close_tag(lfncml,label,4)
      info_cml = 'dictRef="cc:hostName"'
      call open_tag(lfncml,label,info_cml,4)
      info_cml='empty'
      call writetag_scalar_string(lfncml,info_cml,5,host)
      call close_tag(lfncml,label,4)
      info_cml = 'dictRef="gr:numProcCPU"'
      call open_tag(lfncml,label,info_cml,4)
      info_cml='empty'
      call writetag_scalar_integer(lfncml,info_cml,5,getcpucount())
      call close_tag(lfncml,label,4)
      info_cml = 'dictRef="gr:numProcGPU"'
      call open_tag(lfncml,label,info_cml,4)
      info_cml='empty'
      call writetag_scalar_integer(lfncml,info_cml,5,numdev)
      call close_tag(lfncml,label,4)
      info_cml = 'dictRef="gr:numRanks"'
      call open_tag(lfncml,label,info_cml,4)
      info_cml='empty'
      call writetag_scalar_integer(lfncml,info_cml,5,
     &                                  nrsets-nummps*numdev)
      call close_tag(lfncml,label,4)
      info_cml = 'dictRef="gr:numAccRanks"'
      call open_tag(lfncml,label,info_cml,4)
      info_cml='empty'
      call writetag_scalar_integer(lfncml,info_cml,5,nummps*numdev)
      call close_tag(lfncml,label,4)
      info_cml = 'dictRef="gr:taskSize"'
      call open_tag(lfncml,label,info_cml,4)
      info_cml='empty'
      call writetag_scalar_integer(lfncml,info_cml,5,ntask)
      call close_tag(lfncml,label,4)
      call swatch(date,time)
      runDate = date(1:8)//'  '//time(1:8)
      info_cml = 'dictRef="cc:runDate"'
      call open_tag(lfncml,label,info_cml,4)
      info_cml='empty'
      call writetag_scalar_string(lfncml,info_cml,5,runDate)
      call close_tag(lfncml,label,4)
      label='parameterList'
      call close_tag(lfncml,label,3)
      label='module'
      call close_tag(lfncml,label,2)

      return
      end subroutine gronor_env_cml



      subroutine gronor_init_cml()
      use gnome_data
      use gnome_parameters
      use cidef

      implicit none
      character(len=1)    :: sep
      character(len=20)   :: label,method,pointgroup,wavefunctionType
      character(len=132)  :: info_cml
      integer             :: indent
      logical             :: scientific_format

* block 2, job specifics

      method = 'NOCI-Fragments'
      pointgroup = 'C1'
      wavefunctionType = 'NOCI'
* some day we might want to add 'Properties' as a new method
     

      label = 'module'
      indent = 2
      info_cml = 'dictRef="cc:initialization" id="jobInitialization"'
      call open_tag(lfncml,label,info_cml,indent)
      label = 'parameterList'
      indent = 3
      info_cml = 'empty'
      call open_tag(lfncml,label,info_cml,indent)
      indent = 4
      label = 'parameter'
      info_cml = 'dictRef="cc:method"'
      call open_tag(lfncml,label,info_cml,4)
      info_cml='empty'
      call writetag_scalar_string(lfncml,info_cml,5,method)
      call close_tag(lfncml,label,4)
      info_cml = 'dictRef="cc:pointGroup"'
      call open_tag(lfncml,label,info_cml,4)
      info_cml='empty'
      call writetag_scalar_string(lfncml,info_cml,5,pointgroup)
      call close_tag(lfncml,label,4)
      info_cml='dictRef="cc:wavefunctionType"'
      call open_tag(lfncml,label,info_cml,4)
      info_cml='empty'
      call writetag_scalar_string(lfncml,info_cml,5,wavefunctionType)
      call close_tag(lfncml,label,4)
      info_cml='dictRef="gr:vectorTitle"'
      call open_tag(lfncml,label,info_cml,4)
      info_cml='empty'
      sep='|'
      call writetag_array_string(lfncml,info_cml,5,sep,name,2)
      call close_tag(lfncml,label,4)
      info_cml='dictRef="gr:integralTitle"'
      call open_tag(lfncml,label,info_cml,4)
      info_cml='empty'
      sep='|'
      call writetag_array_string(lfncml,info_cml,5,sep,namint,2)
      call close_tag(lfncml,label,4)
* we need to calcuate the charge of the system: sum(znuc) - (nelec + 2*nFrozen)
*      info_cml='dictRef="cc:charge"'
*      call open_tag(lfncml,label,info_cml,4)
*      info_cml='units="nonsi:elementary_charge"'
*      call writetag_scalar_real(lfncml,info_cml,5,charge,.false.)
*      call close_tag(lfncml,label,4)
      info_cml='dictRef="cc:spinMultiplicity"'
      call open_tag(lfncml,label,info_cml,4)
      info_cml='empty'
      call writetag_scalar_integer(lfncml,info_cml,5,nspin)
      call close_tag(lfncml,label,4)
      info_cml='dictRef="gr:tauCI"'
      call open_tag(lfncml,label,info_cml,4)
      info_cml='empty'
      scientific_format=.true.
      call writetag_scalar_real(lfncml,info_cml,5,thresh,
     &       scientific_format)
      call close_tag(lfncml,label,4)
      info_cml='dictRef="gr:fragments"'
      call open_tag(lfncml,label,info_cml,4)
      info_cml='empty'
      call writetag_scalar_integer(lfncml,info_cml,5,nmol)
      call close_tag(lfncml,label,4)
      info_cml='dictRef="gr:nMEBFs"'
      call open_tag(lfncml,label,info_cml,4)
      info_cml='empty'
      call writetag_scalar_integer(lfncml,info_cml,5,nbase)
      call close_tag(lfncml,label,4)
      info_cml='dictRef="gr:MEBFs"'
      call open_tag(lfncml,label,info_cml,4)
      info_cml='empty'
      call writetag_matrix_integer(lfncml,info_cml,5,sep,ncombv,
     &   size(ncombv,1),size(ncombv,2))
      call close_tag(lfncml,label,4)
      sep =' '
      info_cml='dictRef="gr:frozen"'
      call open_tag(lfncml,label,info_cml,4)
      info_cml='empty'
      call writetag_array_integer(lfncml,info_cml,5,sep,nFrozen,mstates)
      call close_tag(lfncml,label,4)
      info_cml='dictRef="gr:inactive"'
      call open_tag(lfncml,label,info_cml,4)
      info_cml='empty'
      call writetag_array_integer(lfncml,info_cml,5,sep,inactm,mstates)
      call close_tag(lfncml,label,4)
      info_cml='dictRef="gr:active"'
      call open_tag(lfncml,label,info_cml,4)
      info_cml='empty'
      call writetag_array_integer(lfncml,info_cml,5,sep,nactm,mstates)
      call close_tag(lfncml,label,4)
      info_cml='dictRef="gr:determinants"'
      call open_tag(lfncml,label,info_cml,4)
      info_cml='empty'
      call writetag_array_integer(lfncml,info_cml,5,sep,idetm,mstates)
      call close_tag(lfncml,label,4)
      label='parameterList'
      call close_tag(lfncml,label,3)
* remains to dump the molecule on the cml file
      call gronor_molecule_cml()
      label='module'
      call close_tag(lfncml,label,2)

      return
      end subroutine gronor_init_cml


      subroutine gronor_results_header_cml
      use cidef

      implicit none
      character(len=132)   :: info_cml
      character(len=20)    :: label

      label = 'module'
      info_cml ='dictRef="cc:calculation" id="NOCI"'
      call open_tag(lfncml,label,info_cml,2)

      end subroutine gronor_results_header_cml
 
