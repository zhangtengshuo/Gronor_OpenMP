      program make_input

      implicit none

      character (len=128) :: ofile,ifile,card,xfile
      integer :: i,j,num,nt,nt0,ne,no,nelecs,nel(25),mol,numa
      character (len=2) :: atom(100),elem(25),nam(25)
      real (kind=8) :: x(100),y(100),z(100),shift(3)
      character (len=20) :: basis(25),project
      integer :: iunit,junit

      numa=10
      
      elem(1)='H '
      elem(2)='B '
      elem(3)='C '
      elem(4)='N '
      elem(5)='O '
      elem(6)='F '
      elem(7)='P '
      elem(8)='S '
      elem(9)='Cl'
      elem(10)='Br'
      elem(11)='I '

      nam(1)=' H'
      nam(2)=' B'
      nam(3)=' C'
      nam(4)=' N'
      nam(5)=' O'
      nam(6)=' F'
      nam(7)=' P'
      nam(8)=' S'
      nam(9)='Cl'
      nam(10)='Br'
      nam(11)=' I'

      nel(1)=1
      nel(2)=5
      nel(3)=6
      nel(4)=7
      nel(5)=8
      nel(6)=9
      nel(7)=15
      nel(8)=16
      nel(9)=17
      nel(10)=35
      nel(11)=53
      
      basis(1)="h.ano-s...3s2p.     "
      basis(2)="b.ano-s...4s3p2d.   "
      basis(3)="c.ano-s...4s3p2d.   "
      basis(4)="n.ano-s...4s3p2d.   "
      basis(5)="o.ano-s...4s3p2d.   "
      basis(6)="f.ano-s...4s3p2d.   "
      basis(7)="p.ano-s...5s4p3d.   "
      basis(8)="s.ano-s...5s4p3d.   "
      basis(9)="cl.ano-s...5s4p3d.  "
      basis(10)="br.ano-s...6s5p4d.  "
      basis(11)="i.ano-s...6s5p4d.   "
      
      print*,' Select project name:'
      read*,project
      print*,' NWChem optimization output file name'
      read*,ofile
      print*,' CAS, provide #electrons and #active orbitals'
      read*,ne,no
      print*,' Provide distance between molecules, dx, dy and dz'
      read*,shift

      if(index(ofile,'.nwout').gt.0) then
        iunit=10
        open(unit=iunit,file=ofile)
 1      continue
        read(iunit,100,end=99) card
 100    format(a)
        if(card(1:19).eq.' Output coordinates') then
          num=0
          read(iunit,100) card
          read(iunit,100) card
          read(iunit,100) card
 2        continue
          num=num+1
          read(iunit,101) i,atom(num),x(num),y(num),z(num)
 101      format(i5,1x,a2,25x,3f15.8)
          if(i.eq.0) goto 1
          goto 2
        else
          goto 1
        endif
 99     continue
        close(unit=iunit)
        num=num-1
        junit=12
        xfile=ofile(1:index(ofile,'.nwout'))//'xyz'
        print*,'xfile=',trim(xfile)
        open(unit=junit,file=xfile)
        write(junit,102) num
 102    format(i5,/,'Coordinates in Angstrom')
        do i=1,num
          write(junit,103) atom(i),x(i),y(i),z(i)
 103      format(a2,1x,3f15.8)
        enddo
        close(unit=junit)
      else
        if(index(ofile,'.xyz').eq.0) then
          xfile=ofile//'.xyz'
        else
          xfile=ofile
        endif
      endif

      print*,'xfile=',trim(xfile)
      iunit=10
      open(unit=iunit,file=xfile)
      read(iunit,*) num
      read(iunit,104) card
 104  format(a)
      do i=1,num
        read(iunit,*) atom(i),x(i),y(i),z(i)
        write(*,103) atom(i),x(i),y(i),z(i)
      enddo
      close(unit=iunit)

      nelecs=0
      do j=1,numa
        do i=1,num
          if(atom(i).eq.elem(j)) nelecs=nelecs+nel(j)
        enddo
      enddo

      print*,'nelecs=',nelecs

      do i=1,3
        shift(i)=0.5*shift(i)
      enddo

      mol=1

      iunit=11
      open(unit=iunit,file=trim(project)//'_A.input')
      
      write(iunit,200)
 200  format('&seward',/,'high cholesky')
      do j=1,numa
        nt=0
        do i=1,num
          if(atom(i).eq.elem(j)) nt=nt+1
        enddo
        if(nt.gt.0) then
          write(iunit,201) basis(j)
 201      format('basis set',/,a20)
          nt=0
          do i=1,num
            if(atom(i).eq.elem(j)) then
              nt=nt+1
              if(nt.le.9) then
                write(iunit,202) nam(j),nt,
     &               x(i)-shift(1),y(i)-shift(2),z(i)-shift(3)
 202            format(a2,i1,' ',3f15.8,' Angstrom')
              else
                write(iunit,203) nam(j),nt,
     &               x(i)-shift(1),y(i)-shift(2),z(i)-shift(3)
 203            format(a2,i2,3f15.8,' Angstrom')
              endif
            endif
          enddo
          write(iunit,204)
 204      format('end basis')
        endif
      enddo
      write(iunit,205)
 205  format(/,'>>> COPY $Project.OneInt  $CurrDir/ONEINT1',/,
     &     '>>> COPY $Project.RunFile $CurrDir/RUNFIL1',//,
     &     '&scf',/)
      
      write(iunit,206) ne,1,(nelecs-ne)/2,no,mol,1,1,trim(project),1
 206  format('&rasscf',/,'nactel',/,i3,/,'spin',/,i3,/,
     &     'inactive',/,i3,/,'ras2',/,i3,/,
     &     'prwf',/,'  0',/,'prsd',//,
     &     ">>> COPY $Project.RasOrb.1 $CurrDir/INPORB.",i1,'_',i1,/,
     &     '>>> COPY $Project.VecDet.',i1,
     &     ' $CurrDir/',a,'_',i3.3,'.det',/)
      
      write(iunit,208) ne,1,(nelecs-ne)/2,no,mol,2,2,trim(project),2
 208  format('&rasscf',/,'nactel',/,i3,/,'spin',/,i3,/,
     &     'inactive',/,i3,/,'ras2',/,i3,/,
     &     'CIRoot',/,'  1 2',/,'  2',/,
     &     'prwf',/,'  0',/,'prsd',//,
     &     '>>> COPY $Project.RasOrb.1 $CurrDir/INPORB.',i1,'_',i1,/,
     &     '>>> COPY $Project.VecDet.',i1,
     &     ' $CurrDir/',a,'_',i3.3,'.det',/)
      
      write(iunit,206) ne,3,(nelecs-ne)/2,no,mol,3,1,trim(project),3
      write(iunit,206) ne-1,2,(nelecs-ne)/2,no,mol,4,1,trim(project),4
      write(iunit,206) ne+1,2,(nelecs-ne)/2,no,mol,5,1,trim(project),5

      close(unit=iunit)

      mol=2
      open(unit=iunit,file=trim(project)//'_B.input')
      
      write(iunit,200)
      do j=1,numa
        nt=0
        do i=1,num
          if(atom(i).eq.elem(j)) nt=nt+1
        enddo
        nt0=nt
        if(nt.gt.0) then
          write(iunit,201) basis(j)
          nt=0
          do i=1,num
            if(atom(i).eq.elem(j)) then
              nt=nt+1
              if(nt+nt0.le.9) then
                write(iunit,202) nam(j),nt+nt0,
     &               x(i)+shift(1),y(i)+shift(2),z(i)+shift(3)
              else
                write(iunit,203) nam(j),nt+nt0,
     &               x(i)+shift(1),y(i)+shift(2),z(i)+shift(3)
              endif
            endif
          enddo
          write(iunit,204)
        endif
      enddo
      write(iunit,209)
 209  format('>>> COPY $Project.OneInt  $CurrDir/ONEINT2',/,
     &     '>>> COPY $Project.RunFile $CurrDir/RUNFIL2',//,
     &     '&scf',/)
      
      write(iunit,206) ne,1,(nelecs-ne)/2,no,mol,1,1,trim(project),6      
      write(iunit,208) ne,1,(nelecs-ne)/2,no,mol,2,2,trim(project),7
      write(iunit,206) ne,3,(nelecs-ne)/2,no,mol,3,1,trim(project),8
      write(iunit,206) ne-1,2,(nelecs-ne)/2,no,mol,4,1,trim(project),9
      write(iunit,206) ne+1,2,(nelecs-ne)/2,no,mol,5,1,trim(project),10

      close(unit=iunit)

      open(unit=iunit,file=trim(project)//'_D.input')
      open(unit=junit,file=trim(project)//'_dimer.xyz')
      write(junit,102) 2*num
      
      write(iunit,200)
      do j=1,numa
        nt=0
        do i=1,num
          if(atom(i).eq.elem(j)) nt=nt+1
        enddo
        if(nt.gt.0) then
          write(iunit,201) basis(j)
          nt=0
          do i=1,num
            if(atom(i).eq.elem(j)) then
              nt=nt+1
              if(nt.le.9) then
                write(iunit,202) nam(j),nt,
     &               x(i)-shift(1),y(i)-shift(2),z(i)-shift(3)
              else
                write(iunit,203) nam(j),nt,
     &               x(i)-shift(1),y(i)-shift(2),z(i)-shift(3)
              endif
              write(junit,103) atom(i),
     &               x(i)-shift(1),y(i)-shift(2),z(i)-shift(3)
            endif
          enddo
          write(iunit,204)
        endif
      enddo
      do j=1,numa
        nt=0
        do i=1,num
          if(atom(i).eq.elem(j)) nt=nt+1
        enddo
        nt0=nt
        if(nt.gt.0) then
          write(iunit,201) basis(j)
          nt=0
          do i=1,num
            if(atom(i).eq.elem(j)) then
              nt=nt+1
              if(nt+nt0.le.9) then
                write(iunit,202) nam(j),nt+nt0,
     &               x(i)+shift(1),y(i)+shift(2),z(i)+shift(3)
              else
                write(iunit,203) nam(j),nt+nt0,
     &               x(i)+shift(1),y(i)+shift(2),z(i)+shift(3)
              endif
              write(junit,103) atom(i),
     &               x(i)+shift(1),y(i)+shift(2),z(i)+shift(3)
            endif
          enddo
          write(iunit,204)
        endif
      enddo
      
      write(iunit,210)
 210  format(/,'oneonly',/,
     &     '>>> COPY $Project.RunFile $CurrDir/RUNFILE')

      close(unit=iunit)
      close(unit=junit)
      
      open(unit=iunit,file=trim(project)//'_M.input')
      
      write(iunit,200)
      do j=1,numa
        nt=0
        do i=1,num
          if(atom(i).eq.elem(j)) nt=nt+1
        enddo
        if(nt.gt.0) then
          write(iunit,201) basis(j)
          nt=0
          do i=1,num
            if(atom(i).eq.elem(j)) then
              nt=nt+1
              if(nt.le.9) then
                write(iunit,202) nam(j),nt,
     &               x(i)-shift(1),y(i)-shift(2),z(i)-shift(3)
              else
                write(iunit,203) nam(j),nt,
     &               x(i)-shift(1),y(i)-shift(2),z(i)-shift(3)
              endif
            endif
          enddo
          write(iunit,204)
        endif
      enddo
      do j=1,numa
        nt=0
        do i=1,num
          if(atom(i).eq.elem(j)) nt=nt+1
        enddo
        nt0=nt
        if(nt.gt.0) then
          write(iunit,201) basis(j)
          nt=0
          do i=1,num
            if(atom(i).eq.elem(j)) then
              nt=nt+1
              if(nt+nt0.le.9) then
                write(iunit,202) nam(j),nt+nt0,
     &               x(i)+shift(1),y(i)+shift(2),z(i)+shift(3)
              else
                write(iunit,203) nam(j),nt+nt0,
     &               x(i)+shift(1),y(i)+shift(2),z(i)+shift(3)
              endif
            endif
          enddo
          write(iunit,204)
        endif
      enddo
      
      write(iunit,211)
 211  format('>>> COPY $CurrDir/COMMONORB INPORB',//,
     &     '&motra',/,'LumOrb',/,'frozen',/,'  0',/,
     &     'deleted',/,' $DELETED',/,'ctonly',/,'kpq',//,
     &     '>>> COPY $Project.RunFile $CurrDir/RUNFILE',/,
     &     '>>> COPY $Project.OneInt  $CurrDir/ONEINT',/,
     &     '>>> COPY $Project.TraOne  $CurrDir/TRAONE',/,
     &     '>>> COPY $Project.ChVec1  $CurrDir/CHVEC1',/,
     &     '>>> COPY $Project.ChRed   $CurrDir/CHRED',/,
     &     '>>> COPY $Project.ChRst   $CurrDir/CHORST',/,
     &     '>>> COPY $Project.ChMap   $CurrDir/CHOMAP',/,
     &     '>>> COPY _CHMOT1 $CurrDir/_CHMOT1',/,
     &     '>>> eval NPROCS = $MOLCAS_NPROCS - 1',/,
     &     '>>> foreach L in (1 .. $NPROCS )',/,
     &     '>>> shell cat tmp_$L/_CHMOT1 >> $CurrDir/_CHMOT1',/,
     &     '>>> enddo')
      
      close(unit=iunit)
      
      open(unit=iunit,file=trim(project)//'_CB.input')

      write(iunit,212) project
 212  format('Project',/,1x,a,/,'Fragments',/,'  2',/,'  5 5',/,
     &     'Threshold',/,' 1e-5',/,'Labels',/,
     &     ' S0 S1 T1 D+ D- S0 S1 T1 D+ D-',/,'Energies')
      
      close(unit=iunit)
      
      open(unit=iunit,file=trim(project)//'.run')

      write(iunit,213) trim(project)
 213  format('#!/usr/bin/tcsh',/,
     &     'setenv MOLCAS_NPROCS 48',/,
     &     'setenv PROJECT "',a,'"',/,
     &     'pymolcas -clean $PROJECT"_A.input" > $PROJECT"_A.output"',/,
     &     'pymolcas -clean $PROJECT"_B.input" > $PROJECT"_B.output"',/,
     &     'pymolcas -clean $PROJECT"_D.input" > $PROJECT"_D.output"',/,
     &     'common_basis < $PROJECT"_CB.input" > $PROJECT"_CB.output"',/,
     &     'setenv DELETED ` grep "Deleted orbitals in MOTRA"',
     &     ' $PROJECT"_CB.output" | cut -b42- `',/,
     &     'touch TRAINT',/,
     &     'setenv MOLCAS_MEM 4096',/,
     &     'pymolcas -clean $PROJECT"_M.input" > $PROJECT"_M.output"',/,
     &     'setenv OMP_NUM_THREADS 12',/,
     &     'rdcho $MOLCAS_NPROCS',/,
     &     'rdtraint < $PROJECT"_CB.input" > $PROJECT"_RT.output"',/,
     &     'unsetenv DELETED',/,
     &     'mpirun -n 48 gronor $PROJECT"_test"')
      
      close(unit=iunit)
      
      open(unit=iunit,file=trim(project)//'_test.inp')

      write(iunit,214)
 214  format('MEBFs 7',/,
     &     '  1  2  1  2  3  4  5',/,
     &     '  6  6  7  7  8 10  9',/,
     &     'Batch 32',/,
     &     'Threshold 1e-4',/,
     &     'Timings 4',/,
     &     'Size 1',/,
     &     'Print medium')
      
      close(unit=iunit)
      
      end
      
