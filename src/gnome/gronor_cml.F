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
      info_cml = '  convention="convention:compchem">'
      write(lfncml,'(a)') trim(info_cml)
      return

      end subroutine gronor_prelude_cml



      subroutine gronor_env_cml(version,version_type)
      use cidef
      use gnome_parameters
      use cidist

      implicit none
      character(len=5)    :: version
      character(len=6)    :: programName
      character(len=18)   :: runDate
      character(len=20)   :: label
      character(len=64)   :: version_type
      character(len=132)  :: info_cml
      integer             :: indent
      integer,external    :: getcpucount

* Open the calculation
      label = 'module'
      indent = 1
      info_cml = 'id="job" dictRef="cc:jobList"'
      call open_tag(lfncml,label,info_cml,indent)
      label = 'module'
      indent = 2
      info_cml = 'id="job" dictRef="cc:job"'
      call open_tag(lfncml,label,info_cml,indent)
* Block 1: general job info_cml
      indent = 3
      info_cml = 'dictRef="cc:environment" id="environment"'
      call open_tag(lfncml,label,info_cml,indent)
      label = 'parameterList'
      indent = 4
      info_cml = 'empty'
      call open_tag(lfncml,label,info_cml,indent)
      indent = 5
      label = 'parameter'
      info_cml = 'dictRef="cc:program"'
      call open_tag(lfncml,label,info_cml,indent)
      indent = 6
      info_cml = 'empty'
      programName = 'GronOR'
      call writetag_scalar_string(lfncml,info_cml,indent,
     &             programName)
      call close_tag(lfncml,label,5)
      indent = 5
      label = 'parameter'
      info_cml = 'dictRef="cc:programVersion"'
      call open_tag(lfncml,label,info_cml,indent)
      indent = 6
      info_cml = 'empty'
      call writetag_scalar_string(lfncml,info_cml,indent,
     &             version)
      call close_tag(lfncml,label,5)
      info_cml = 'dictRef="cc:programSubversion"'
      call open_tag(lfncml,label,info_cml,5)
      info_cml='empty'
      call writetag_scalar_string(lfncml,info_cml,6,
     &            version_type)
      call close_tag(lfncml,label,5)
      info_cml = 'dictRef="cc:executable"'
      call open_tag(lfncml,label,info_cml,5)
      info_cml='empty'
      call writetag_scalar_string(lfncml,info_cml,6,command)
      call close_tag(lfncml,label,5)
      info_cml = 'dictRef="cc:hostName"'
      call open_tag(lfncml,label,info_cml,5)
      info_cml='empty'
      call writetag_scalar_string(lfncml,info_cml,6,host)
      call close_tag(lfncml,label,5)
      info_cml = 'dictRef="cc:numProc"'
      call open_tag(lfncml,label,info_cml,5)
      info_cml='empty'
      call writetag_scalar_integer(lfncml,info_cml,6,getcpucount())
      call close_tag(lfncml,label,5)
      info_cml = 'dictRef="gr:numProcGPU"'
      call open_tag(lfncml,label,info_cml,5)
      info_cml='empty'
      call writetag_scalar_integer(lfncml,info_cml,6,numdev)
      call close_tag(lfncml,label,5)
      info_cml = 'dictRef="gr:numRanks"'
      call open_tag(lfncml,label,info_cml,5)
      info_cml='empty'
      call writetag_scalar_integer(lfncml,info_cml,6,
     &                                  nrsets-nummps*numdev)
      call close_tag(lfncml,label,5)
      info_cml = 'dictRef="gr:numAccRanks"'
      call open_tag(lfncml,label,info_cml,5)
      info_cml='empty'
      call writetag_scalar_integer(lfncml,info_cml,6,nummps*numdev)
      call close_tag(lfncml,label,5)
      info_cml = 'dictRef="gr:taskSize"'
      call open_tag(lfncml,label,info_cml,5)
      info_cml='empty'
      call writetag_scalar_integer(lfncml,info_cml,6,ntask)
      call close_tag(lfncml,label,5)
      call swatch(date,time)
      runDate = date(1:8)//'  '//time(1:8)
      info_cml = 'dictRef="cc:runDate"'
      call open_tag(lfncml,label,info_cml,5)
      info_cml='empty'
      call writetag_scalar_string(lfncml,info_cml,6,runDate)
      call close_tag(lfncml,label,5)
      label='parameterList'
      call close_tag(lfncml,label,4)
      label='module'
      call close_tag(lfncml,label,3)

      return
      end subroutine gronor_env_cml



      subroutine gronor_init_cml()
      use gnome_data
      use gnome_parameters
      use cidef

      implicit none
      character(len=1)    :: sep
      character(len=20)   :: label,method,pointgroup,wavefunctionType
      character(len=132)  :: info_cml,fmt_1
      integer             :: indent,charge,i

* block 2, job specifics

      method = 'NOCI-Fragments'
      pointgroup = 'C1'
      wavefunctionType = 'NOCI'
* some day we might want to add 'Properties' as a new method
     

      label = 'module'
      indent = 3
      info_cml = 'dictRef="cc:initialization" id="jobInitialization"'
      call open_tag(lfncml,label,info_cml,indent)
      label = 'parameterList'
      indent = 4
      info_cml = 'empty'
      call open_tag(lfncml,label,info_cml,indent)
      indent = 5
      label = 'parameter'
      info_cml = 'dictRef="cc:method"'
      call open_tag(lfncml,label,info_cml,indent)
      info_cml='empty'
      call writetag_scalar_string(lfncml,info_cml,6,method)
      call close_tag(lfncml,label,5)
      info_cml = 'dictRef="cc:pointGroup"'
      call open_tag(lfncml,label,info_cml,5)
      info_cml='empty'
      call writetag_scalar_string(lfncml,info_cml,6,pointgroup)
      call close_tag(lfncml,label,5)
      info_cml='dictRef="cc:wavefunctionType"'
      call open_tag(lfncml,label,info_cml,5)
      info_cml='empty'
      call writetag_scalar_string(lfncml,info_cml,6,wavefunctionType)
      call close_tag(lfncml,label,5)
      info_cml='dictRef="gr:vectorTitle"'
      call open_tag(lfncml,label,info_cml,5)
      info_cml='empty'
      sep='|'
      call writetag_array_string(lfncml,info_cml,6,sep,name,2)
      call close_tag(lfncml,label,5)
      info_cml='dictRef="gr:integralTitle"'
      call open_tag(lfncml,label,info_cml,5)
      info_cml='empty'
      sep='|'
      call writetag_array_string(lfncml,info_cml,6,sep,namint,2)
      call close_tag(lfncml,label,5)
      charge = zNucTot
      do i = 1, nmol
        charge = charge - nElectrons(ncombv(i,1))
      end do
      info_cml='dictRef="cc:charge"'
      call open_tag(lfncml,label,info_cml,5)
      info_cml='units="nonsi:elementary_charge"'
      fmt_1 = 'f10.4'
      call writetag_scalar_real(lfncml,info_cml,6,charge,fmt_1)
      call close_tag(lfncml,label,5)
      info_cml='dictRef="cc:spinMultiplicity"'
      call open_tag(lfncml,label,info_cml,5)
      info_cml='empty'
      call writetag_scalar_integer(lfncml,info_cml,6,nspin)
      call close_tag(lfncml,label,5)
      info_cml='dictRef="gr:tauCI"'
      call open_tag(lfncml,label,info_cml,5)
      info_cml='empty'
      fmt_1='e8.3'
      call writetag_scalar_real(lfncml,info_cml,6,thresh,fmt_1)
      call close_tag(lfncml,label,5)
      info_cml='dictRef="cc:nuclearRepulsionEnergy"'
      call open_tag(lfncml,label,info_cml,5)
      info_cml = 'units="nonsi:hartree"'
      fmt_1='f15.8'
      call writetag_scalar_real(lfncml,info_cml,6,potnuc,fmt_1)
      call close_tag(lfncml,label,5)
      info_cml='dictRef="gr:fragments"'
      call open_tag(lfncml,label,info_cml,5)
      info_cml='empty'
      call writetag_scalar_integer(lfncml,info_cml,6,nmol)
      call close_tag(lfncml,label,5)
      info_cml='dictRef="gr:nMEBFs"'
      call open_tag(lfncml,label,info_cml,5)
      info_cml='empty'
      call writetag_scalar_integer(lfncml,info_cml,6,nbase)
      call close_tag(lfncml,label,5)
      info_cml='dictRef="gr:MEBFs"'
      call open_tag(lfncml,label,info_cml,5)
      info_cml='empty'
      call writetag_matrix_integer(lfncml,info_cml,6,sep,ncombv,
     &   size(ncombv,1),size(ncombv,2))
      call close_tag(lfncml,label,5)
      if (mebfLabels) then
        info_cml='dictRef="gr:mebfLabels"'
        call open_tag(lfncml,label,info_cml,5)
        info_cml='empty'
        sep='|'
        call writetag_array_string(lfncml,info_cml,6,sep,mebfLabel,
     &     size(mebfLabel))
        call close_tag(lfncml,label,5)
      end if
      sep =' '
      info_cml='dictRef="gr:frozen"'
      call open_tag(lfncml,label,info_cml,5)
      info_cml='empty'
      call writetag_array_integer(lfncml,info_cml,6,sep,nFrozen,mstates)
      call close_tag(lfncml,label,5)
      info_cml='dictRef="gr:inactive"'
      call open_tag(lfncml,label,info_cml,5)
      info_cml='empty'
      call writetag_array_integer(lfncml,info_cml,6,sep,inactm,mstates)
      call close_tag(lfncml,label,5)
      info_cml='dictRef="gr:active"'
      call open_tag(lfncml,label,info_cml,5)
      info_cml='empty'
      call writetag_array_integer(lfncml,info_cml,6,sep,nactm,mstates)
      call close_tag(lfncml,label,5)
      info_cml='dictRef="gr:determinants"'
      call open_tag(lfncml,label,info_cml,5)
      info_cml='empty'
      call writetag_array_integer(lfncml,info_cml,6,sep,idetm,mstates)
      call close_tag(lfncml,label,5)
      label='parameterList'
      call close_tag(lfncml,label,4)
* remains to dump the molecule on the cml file
      indent = 4
      call gronor_molecule_cml(indent)
      label='module'
      call close_tag(lfncml,label,3)

      return
      end subroutine gronor_init_cml


      subroutine gronor_results_header_cml
      use cidef

      implicit none
      character(len=132)   :: info_cml
      character(len=20)    :: label

      label = 'module'
      info_cml ='dictRef="cc:calculation" id="NOCI"'
      call open_tag(lfncml,label,info_cml,3)

      end subroutine gronor_results_header_cml
 

      subroutine gronor_finalize_cml
      use cidef

      implicit none
      integer              :: indent
      character(len=132)   :: info_cml,fmt_1
      character(len=20)    :: label
      character(len=18)    :: stopDate
      real(kind=8), external :: timer_wall_total

      label='module'
      indent = 3
      call close_tag(lfncml,label,indent)
      indent = 3
      info_cml = 'dictRef="cc:finalization"'
      call open_tag(lfncml,label,info_cml,indent)
      label = 'parameterList'
      indent = 4
      info_cml = 'empty'
      call open_tag(lfncml,label,info_cml,indent)
      label = 'parameter'
      call swatch(date,time)
      stopDate = date(1:8)//'  '//time(1:8)
      info_cml = 'dictRef="cc:jobdatetime.end"'
      call open_tag(lfncml,label,info_cml,5)
      info_cml='empty'
      call writetag_scalar_string(lfncml,info_cml,6,stopDate)
      call close_tag(lfncml,label,5)
      info_cml = 'dictRef="cc.wallTime"'
      call open_tag(lfncml,label,info_cml,5)
      info_cml='units="si:s"'
      fmt_1='f12.3'
      call writetag_scalar_real(lfncml,info_cml,6,timer_wall_total(99),
     &                                fmt_1)
      call close_tag(lfncml,label,5)
      label = 'parameterList'
      call close_tag(lfncml,label,4)
      label='module'
      indent = 3
      call close_tag(lfncml,label,indent)
      indent = 2
      call close_tag(lfncml,label,indent)
      indent = 1
      call close_tag(lfncml,label,indent)
      indent = 0
      call close_tag(lfncml,label,indent)

      end subroutine gronor_finalize_cml
