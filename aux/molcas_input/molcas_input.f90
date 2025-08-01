      program molcas_input

        implicit none

        external molcas_rotate

      integer (kind=4) :: iarg
      integer (kind=8) :: lfninp, lfnout, lfnxyz, lfnnwo
      character (len=255) :: filename, filinp, filxyz, filnwo
      character (len=255) :: card, options
      integer :: i,j,k,m,num,nt,nt0,ne,no,nelecs,nel(25),mol1,mol2,numa
      character (len=2) :: atom(100),elem(25),nam(25)
      real (kind=8) :: a,x1(100,3),x2(100,3),xc(3),v(3),w(3),angle,x(3),y(3)
      character (len=255) :: basis(25),project
      integer :: iunit,junit,kunit,lunit,nr,memory,nalter,ialter(3,25)
      character (len=1) :: operation
      logical :: done,aonly,caspt2

      numa=12
      nr=12
      memory=1024
      nalter=0
      aonly=.false.
      caspt2=.false.
      
      elem(1)='H '
      elem(2)='B '
      elem(3)='C '
      elem(4)='N '
      elem(5)='O '
      elem(6)='F '
      elem(7)='P '
      elem(8)='S '
      elem(9)='Cl'
      elem(10)='Se'
      elem(11)='Br'
      elem(12)='I '

      nam(1)=' H'
      nam(2)=' B'
      nam(3)=' C'
      nam(4)=' N'
      nam(5)=' O'
      nam(6)=' F'
      nam(7)=' P'
      nam(8)=' S'
      nam(9)='Cl'
      nam(10)='Se'
      nam(11)='Br'
      nam(12)=' I'

      nel(1)=1
      nel(2)=5
      nel(3)=6
      nel(4)=7
      nel(5)=8
      nel(6)=9
      nel(7)=15
      nel(8)=16
      nel(9)=17
      nel(10)=34
      nel(11)=35
      nel(12)=53
      
      basis(1)="h.ano-s...3s2p.     "
      basis(2)="b.ano-s...4s3p2d.   "
      basis(3)="c.ano-s...4s3p2d.   "
      basis(4)="n.ano-s...4s3p2d.   "
      basis(5)="o.ano-s...4s3p2d.   "
      basis(6)="f.ano-s...4s3p2d.   "
      basis(7)="p.ano-s...5s4p3d.   "
      basis(8)="s.ano-s...5s4p3d.   "
      basis(9)="cl.ano-s...5s4p3d.  "
      basis(10)="se.ano-s...6s5p4d.  "
      basis(11)="br.ano-s...6s5p4d.  "
      basis(12)="i.ano-s...6s5p4d.   "

      mol1=1
      mol2=2
      
      iarg=1
      call getarg(iarg,filename)

      lfninp=10
      lfnout=11
      lfnnwo=12
      lfnxyz=13
!      filinp=filename(1:index(filename,' ')-1)//'.min '
      project=trim(filename)
      filnwo=filename(1:index(filename,' ')-1)//'.nwout '
      filinp=filename(1:index(filename,' ')-1)//'.min '
!      filinp='molcas.min'
      
      open(unit=lfninp,file=filinp(1:index(filinp,' ')-1))

      done=.false.
      operation='n'
      
      do while(.not.done)
        
!     Read operation from input file

        if(operation.eq.'&') then
          read(lfninp,1000,end=1009) card
1000      format(a)

          if(card(2:2).eq.' ') then
            operation=card(1:1)
            options(1:254)=card(2:255)
          endif
          if(card(1:6).eq.'Rotate'.or.card(1:6).eq.'rotate') then
            operation='R'
            options(1:249)=card(7:255)
          endif
          if(card(1:6).eq.'Center'.or.card(1:6).eq.'center') then
            operation='C'
            options(1:249)=card(7:255)
          endif
          if(card(1:6).eq.'CASPT2'.or.card(1:6).eq.'caspt2') then
            operation='D'
            options(1:249)=card(7:255)
          endif
          if(card(1:6).eq.'NWChem'.or.card(1:6).eq.'nwchem') then
            operation='N'
            options(1:249)=card(7:255)
          endif
          if(card(1:7).eq.'Project'.or.card(1:7).eq.'project') then
            operation='P'
            options(1:248)=card(8:255)
          endif
          if(card(1:3).eq.'XYZ'.or.card(1:3).eq.'xyz') then
            operation='X'
            options(1:252)=card(4:255)

          endif
          if(card(1:5).eq.'Aonly'.or.card(1:5).eq.'aonly') then
            operation='A'
            options(1:250)=card(6:255)
          endif
          if(card(1:5).eq.'Ranks'.or.card(1:5).eq.'ranks') then
            operation='r'
            options(1:250)=card(6:255)
          endif
          if(card(1:6).eq.'Memory'.or.card(1:5).eq.'memory') then
            operation='M'
            options(1:249)=card(7:255)
          endif
          if(card(1:5).eq.'Alter'.or.card(1:5).eq.'alter') then
            operation='a'
            options(1:250)=card(6:255)
          endif
          if(card(1:6).eq.'Orient'.or.card(1:5).eq.'orient') then
            operation='O'
            options(1:249)=card(7:255)
          endif
          if(card(1:5).eq.'Space'.or.card(1:5).eq.'space') then
            operation='S'
            options(1:250)=card(6:255)
          endif
          if(card(1:9).eq.'Translate'.or.card(1:9).eq.'translate') then
            operation='T'
            options(1:246)=card(10:255)
          endif
          if(card(1:5).eq.'Write'.or.card(1:5).eq.'write') then
            operation='W'
            options(1:250)=card(6:255)
          endif
          if(card(1:4).eq.'Quit'.or.card(1:4).eq.'quit') then
            operation='Q'
            options(1:251)=card(5:255)
          endif
          if(card(1:4).eq.'Exit'.or.card(1:4).eq.'exit') then
            operation='Q'
            options(1:251)=card(5:255)
          endif
          if(card(1:1).eq.' '.or.card(1:1).eq.' ') then
            operation='Q'
            options(1:254)=card(2:255)
          endif
        endif
        if(operation.eq.'*') operation='&'
        if(operation.eq.'!') operation='&'
        if(operation.eq.'#') operation='&'

!     Set project name
        
        if(operation.eq.'P') then
          do while(options(1:1).eq.' ')
            options(1:254)=options(2:255)
          enddo
          project=trim(options)
          write(*,6000) trim(project)
6000      format(' Project name:',a)
          operation='&'
        endif
        
        if(operation.eq.'p') then
          write(*,6000) trim(project)
          operation='&'
        endif
        
        if(operation.eq.'D') then
          caspt2=.true.
          operation='&'
        endif

!     Read NWChem output file

        if(operation.eq.'N') then
          do while(options(1:1).eq.' ')
            options(1:254)=options(2:255)
          enddo
          filnwo=trim(options)
          write(*,6001) trim(filnwo)
 6001     format(' NWChem file:',a)
          open(unit=lfnnwo,file=trim(filnwo))
 1        continue
          read(lfnnwo,100,end=99) card
 100      format(a)
          if(card(1:19).eq.' Output coordinates') then
            num=0
            read(lfnnwo,100) card
            read(lfnnwo,100) card
            read(lfnnwo,100) card
 2          continue
            num=num+1
            read(lfnnwo,101) i,atom(num),x1(num,1),x1(num,2),x1(num,3)
 101        format(i5,1x,a2,25x,3f15.8)
            if(i.eq.0) goto 1
            goto 2
          else
            goto 1
          endif
 99       continue
          close(unit=lfnnwo)
          num=num-1
          filxyz=filnwo(1:index(filnwo,'.nwout'))//'.xyz'
          open(unit=lfnxyz,file=filxyz)
          write(lfnxyz,102) num
 102      format(i5,/,'Coordinates in Angstrom')
          do i=1,num
            write(lfnxyz,103) atom(i),x1(i,1),x1(i,2),x1(i,3)
 103        format(a2,1x,3f15.8)
          enddo
          close(unit=lfnxyz)
        endif
        
        if(operation.eq.'n') then
          write(*,6001) trim(filnwo)
          open(unit=lfnnwo,file=trim(filnwo))
 11        continue
          read(lfnnwo,100,end=991) card
          if(card(1:19).eq.' Output coordinates') then
            num=0
            read(lfnnwo,100) card
            read(lfnnwo,100) card
            read(lfnnwo,100) card
 21          continue
            num=num+1
            read(lfnnwo,101) i,atom(num),x1(num,1),x1(num,2),x1(num,3)
            if(i.eq.0) goto 11
            goto 21
          else
            goto 11
          endif
 991       continue
          close(unit=lfnnwo)
          num=num-1
          filxyz=filnwo(1:index(filnwo,'.nwout'))//'.xyz'
          open(unit=lfnxyz,file=filxyz)
          write(lfnxyz,102) num
          do i=1,num
            write(lfnxyz,103) atom(i),x1(i,1),x1(i,2),x1(i,3)
          enddo
          close(unit=lfnxyz)
        endif
        
!     Read XYZ coordinate file
        if(operation.eq.'X') then
          filxyz=options(2:index(options,'.xyz')+3)
!          write(*,6002) trim(filxyz)
 6002     format(' XYZ file:',a)
          open(unit=lfnxyz,file=trim(filxyz))
          rewind(lfnxyz)
          read(lfnxyz,*) num
          read(lfnxyz,104) card
 104      format(a)
          do i=1,num
            read(lfnxyz,*) atom(i),x1(i,1),x1(i,2),x1(i,3)
!            write(*,103) atom(i),x1(i,2),x1(i,2),x1(i,3)
          enddo
          close(unit=lfnxyz)
        endif
        
!     Set A only

        if(operation.eq.'A') then
          aonly=.true.
          operation='&'
        endif

!     Set number of ranks

        if(operation.eq.'r') then
          read(options,*) nr
          operation='&'
        endif
        
        if(operation.eq.'M') then
          read(options,*) memory
          operation='&'
        endif
        
        if(operation.eq.'a') then
          read(options,*) nalter,((ialter(i,j),i=1,3),j=1,nalter)
          operation='&'
        endif
        
        if(operation.eq.'O') then
          do k=1,3
            xc(k)=0.0d0
            do i=1,num
              xc(k)=xc(k)+x1(i,k)+x2(i,k)
            enddo
            xc(k)=xc(k)/dble(num)
          enddo
          do k=1,3
            do i=1,num
              x1(i,k)=x1(i,k)-xc(k)
            enddo
          enddo
          a=0.0d0
          do i=1,num
            if(abs(x1(i,1)).gt.a) then
              m=i
              a=x1(i,1)
            endif
            if(abs(x1(i,2)).gt.a) then
              m=i
              a=x1(i,2)
            endif
          enddo
          angle=-atan(x1(m,2)/x1(m,1))
          w(1)=0.0d0
          w(2)=0.0d0
          w(3)=1.0d0
          v(1)=0.0d0
          v(2)=0.0d0
          v(3)=0.0d0
          do i=1,num
            x(1)=x1(i,1)
            x(2)=x1(i,2)
            x(3)=x1(i,3)
            call molcas_rotate(v,w,angle,x,y)
            x1(i,1)=y(1)
            x1(i,2)=y(2)
            x1(i,3)=y(3)
          enddo
          do k=1,3
            do i=1,num
              x2(i,k)=x2(i,k)+xc(k)
            enddo
          enddo
        endif
        
!     Determine number of electrons
        
        if(operation.eq.'N'.or.operation.eq.'n'.or.operation.eq.'X'.or.operation.eq.'O') then
          nelecs=0
          do j=1,numa
            do i=1,num
              if(atom(i).eq.elem(j)) nelecs=nelecs+nel(j)
              x2(i,1)=x1(i,1)
              x2(i,2)=x1(i,2)
              x2(i,3)=x1(i,3)
            enddo
          enddo
          operation='&'
        endif

!     Define CAS space
        
        if(operation.eq.'S') then
          read(options,*) ne,no
          operation='&'
        endif

        
!     Center molecule(s)
        
        if(operation.eq.'C') then
          do k=1,3
            xc(k)=0.0d0
            do i=1,num
              xc(k)=xc(k)+x1(i,k)+x2(i,k)
            enddo
            xc(k)=xc(k)/dble(num)
          enddo
          do k=1,3
            do i=1,num
              x1(i,k)=x1(i,k)-xc(k)
              x2(i,k)=x2(i,k)-xc(k)
            enddo
          enddo
          operation='&'
        endif

!     Translate molecule 2

        if(operation.eq.'T') then
          read(options,*) xc
          do k=1,3
            do i=1,num
              x2(i,k)=x2(i,k)+xc(k)
            enddo
          enddo
          operation='&'
        endif

!     Rotate molecule 2

        if(operation.eq.'R') then
          read(options,*) w,angle
          if(abs(angle).gt.6.29) angle=angle*0.01745329252d0
        endif
        if(operation.eq.'R'.or.operation.eq.'r') then
          do k=1,3
            xc(k)=0.0d0
            do i=1,num
              xc(k)=xc(k)+x2(i,k)
            enddo
            xc(k)=xc(k)/dble(num)
            v(i)=0.0d0
          enddo
          do k=1,3
            do i=1,num
              x2(i,k)=x2(i,k)-xc(k)
            enddo
          enddo
          do i=1,num
            x(1)=x2(i,1)
            x(2)=x2(i,2)
            x(3)=x2(i,3)
            call molcas_rotate(v,w,angle,x,y)
            x2(i,1)=y(1)
            x2(i,2)=y(2)
            x2(i,3)=y(3)
          enddo
          do k=1,3
            do i=1,num
              x2(i,k)=x2(i,k)+xc(k)
            enddo
          enddo
          operation='&'
        endif

        if(operation.eq.' ') done=.true.
        if(operation.eq.'Q') done=.true.
        
!     Write results

        if(operation.eq.'W'.or.done) then
          iunit=11
          
          open(unit=iunit,file=trim(project)//'_SCF.input')
          write(iunit,200)
          do j=1,numa
            nt=0
            do i=1,num
              if(atom(i).eq.elem(j)) nt=nt+1
            enddo
            if(nt.gt.0) then
              write(iunit,201) trim(basis(j))
              nt=0
              do i=1,num
                if(atom(i).eq.elem(j)) then
                  nt=nt+1
                  if(nt.le.9) then
                    write(iunit,202) nam(j),nt,x1(i,1),x1(i,2),x1(i,3)
                  else
                    write(iunit,203) nam(j),nt,x1(i,1),x1(i,2),x1(i,3)
                  endif
                endif
              enddo
              write(iunit,204)
            endif
          enddo
          write(iunit,805)
 805      format(/,'&scf',/)

          close(unit=iunit)
          
          open(unit=iunit,file=trim(project)//'_A.input')
          write(iunit,200)
 200      format('&seward',/,'high cholesky')
          do j=1,numa
            nt=0
            do i=1,num
              if(atom(i).eq.elem(j)) nt=nt+1
            enddo
            if(nt.gt.0) then
              write(iunit,201) trim(basis(j))
 201          format('basis set',/,a)
              nt=0
              do i=1,num
                if(atom(i).eq.elem(j)) then
                  nt=nt+1
                  if(nt.le.9) then
                    write(iunit,202) nam(j),nt,x1(i,1),x1(i,2),x1(i,3)
 202                format(a2,i1,' ',3f15.8,' Angstrom')
                  else
                    write(iunit,203) nam(j),nt,x1(i,1),x1(i,2),x1(i,3)
 203                format(a2,i2,3f15.8,' Angstrom')
                  endif
                endif
              enddo
              write(iunit,204)
 204          format('end basis')
            endif
          enddo
          write(iunit,205)
 205      format(/,'>>> COPY $Project.OneInt $CurrDir/ONEINT1',/, &
               '>>> COPY $Project.RunFile $CurrDir/RUNFIL1')
          write(iunit,305)
 305      format(/,'>>> COPY $Project.OneInt $CurrDir/ONEINT2',/, &
               '>>> COPY $Project.RunFile $CurrDir/RUNFIL2',//,'&scf',/)

          write(iunit,606)
606       format('&rasscf')
          if(nalter.gt.0) then
            write(iunit,607) nalter
607         format('alter',/,i3)
            write(iunit,608) ((ialter(i,j),i=1,3),j=1,nalter)
608         format(3i3)
          endif
          write(iunit,206) ne,1,(nelecs-ne)/2,no,mol1,1,1,trim(project),1
206       format('nactel',/,i3,/,'spin',/,i3,/, &
              'inactive',/,i3,/,'ras2',/,i3,/, &
              'prwf',/,'  0',/,'prsd',//, &
              ">>> COPY $Project.RasOrb.1 $CurrDir/INPORB.",i1,'_',i1,/, &
              '>>> COPY $Project.VecDet.',i1, &
              ' $CurrDir/',a,'_',i3.3,'.det',/) 
          write(iunit,306) mol2,1,1,trim(project),6
306       format(">>> COPY $Project.RasOrb.1 $CurrDir/INPORB.",i1,'_',i1,/, &
              '>>> COPY $Project.VecDet.',i1,' $CurrDir/',a,'_',i3.3,'.det',/) 
          if(caspt2) write(iunit,219)
219       format('&caspt2',/)
229       format('&caspt2',/,'Multistate=  1  2',/)
          
          write(iunit,208) ne,1,(nelecs-ne)/2,no,mol1,2,2,trim(project),2
 208      format('&rasscf',/,'nactel',/,i3,/,'spin',/,i3,/, &
               'inactive',/,i3,/,'ras2',/,i3,/, &
               'CIRoot',/,'  1 2',/,'  2',/, &
               'prwf',/,'  0',/,'prsd',//, &
               '>>> COPY $Project.RasOrb.1 $CurrDir/INPORB.',i1,'_',i1,/, &
               '>>> COPY $Project.VecDet.',i1,' $CurrDir/',a,'_',i3.3,'.det',/)
          write(iunit,306) mol2,2,2,trim(project),7
          if(caspt2) write(iunit,229)
          
          write(iunit,606)
          write(iunit,206) ne,3,(nelecs-ne)/2,no,mol1,3,1,trim(project),3
          write(iunit,306) mol2,3,1,trim(project),8
          if(caspt2) write(iunit,219)
          write(iunit,606)
          write(iunit,206) ne-1,2,(nelecs-ne)/2,no,mol1,4,1,trim(project),4
          write(iunit,306) mol2,4,1,trim(project),9
          if(caspt2) write(iunit,219)
          write(iunit,606)
          write(iunit,206) ne+1,2,(nelecs-ne)/2,no,mol1,5,1,trim(project),5
          write(iunit,306) mol2,5,1,trim(project),10
          if(caspt2) write(iunit,219)

          close(unit=iunit)
          
          open(unit=iunit,file=trim(project)//'_SA.input')
          write(iunit,400)
 400      format('&seward',/,'high cholesky')
          do j=1,numa
            nt=0
            do i=1,num
              if(atom(i).eq.elem(j)) nt=nt+1
            enddo
            if(nt.gt.0) then
              write(iunit,401) trim(basis(j))
 401          format('basis set',/,a)
              nt=0
              do i=1,num
                if(atom(i).eq.elem(j)) then
                  nt=nt+1
                  if(nt.le.9) then
                    write(iunit,402) nam(j),nt,x1(i,1),x1(i,2),x1(i,3)
 402                format(a2,i1,' ',3f15.8,' Angstrom')
                  else
                    write(iunit,403) nam(j),nt,x1(i,1),x1(i,2),x1(i,3)
 403                format(a2,i2,3f15.8,' Angstrom')
                  endif
                endif
              enddo
              write(iunit,404)
 404          format('end basis')
            endif
          enddo
          write(iunit,405)
 405      format(/,'>>> COPY $Project.OneInt $CurrDir/ONEINT1',/, &
               '>>> COPY $Project.RunFile $CurrDir/RUNFIL1')
          write(iunit,505)
 505      format(/,'>>> COPY $Project.OneInt $CurrDir/ONEINT2',/, &
               '>>> COPY $Project.RunFile $CurrDir/RUNFIL2',//,'&scf',/)
          
          write(iunit,606)
          if(nalter.gt.0) then
            write(iunit,607) nalter
            write(iunit,608) ((ialter(i,j),i=1,3),j=1,nalter)
          endif
          write(iunit,406) ne,1,(nelecs-ne)/2,no,mol1,1,mol1,2,&
              1,trim(project),1,2,trim(project),2
406       format('nactel',/,i3,/,'spin',/,i3,/, &
              'inactive',/,i3,/,'ras2',/,i3,/, &
              'ciroot',/,' 2 2',/,' 1 2',/,' 1 1',/,&
              'prwf',/,'  0',/,'prsd',//, &
              ">>> COPY $Project.RasOrb $CurrDir/INPORB.",i1,'_',i1,/, &
              ">>> COPY $Project.RasOrb $CurrDir/INPORB.",i1,'_',i1,/, &
              '>>> COPY $Project.VecDet.',i1, &
              ' $CurrDir/',a,'_',i3.3,'.det',/,&
              '>>> COPY $Project.VecDet.',i1, &
              ' $CurrDir/',a,'_',i3.3,'.det',/) 
          write(iunit,506) mol2,1,1,trim(project),6
          write(iunit,506) mol2,2,2,trim(project),7
506       format(">>> COPY $Project.RasOrb $CurrDir/INPORB.",i1,'_',i1,/, &
              '>>> COPY $Project.VecDet.',i1,' $CurrDir/',a,'_',i3.3,'.det',/) 
          if(caspt2) write(iunit,419)
419       format('&caspt2',/,'Multistate=  1  2',/)
          
!          write(iunit,408) ne,1,(nelecs-ne)/2,no,mol1,2,2,trim(project),2
 408      format('&rasscf',/,'nactel',/,i3,/,'spin',/,i3,/, &
               'inactive',/,i3,/,'ras2',/,i3,/, &
               'CIRoot',/,'  1 2',/,'  2',/, &
               'prwf',/,'  0',/,'prsd',//, &
               '>>> COPY $Project.RasOrb.1 $CurrDir/INPORB.',i1,'_',i1,/, &
               '>>> COPY $Project.VecDet.',i1,' $CurrDir/',a,'_',i3.3,'.det',/)
!          write(iunit,506) mol2,2,2,trim(project),7
!          if(caspt2) write(iunit,229)
          
          write(iunit,606)
          write(iunit,206) ne,3,(nelecs-ne)/2,no,mol1,3,1,trim(project),3
          write(iunit,306) mol2,3,1,trim(project),8
          if(caspt2) write(iunit,219)
          write(iunit,606)
          write(iunit,206) ne-1,2,(nelecs-ne)/2,no,mol1,4,1,trim(project),4
          write(iunit,306) mol2,4,1,trim(project),9
          if(caspt2) write(iunit,219)
          write(iunit,606)
          write(iunit,206) ne+1,2,(nelecs-ne)/2,no,mol1,5,1,trim(project),5
          write(iunit,306) mol2,5,1,trim(project),10
          if(caspt2) write(iunit,219)

          close(unit=iunit)

          mol2=2
          
          open(unit=iunit,file=trim(project)//'_B.input')
          
          write(iunit,200)
          do j=1,numa
            nt=0
            do i=1,num
              if(atom(i).eq.elem(j)) nt=nt+1
            enddo
            nt0=nt
            if(nt.gt.0) then
              write(iunit,201) trim(basis(j))
              nt=0
              do i=1,num
                if(atom(i).eq.elem(j)) then
                  nt=nt+1
                  if(nt+nt0.le.9) then
                    write(iunit,202) nam(j),nt+nt0,x2(i,1),x2(i,2),x2(i,3)
                  else
                    write(iunit,203) nam(j),nt+nt0,x2(i,1),x2(i,2),x2(i,3)
                  endif
                endif
              enddo
              write(iunit,204)
            endif
          enddo
          write(iunit,209)
 209      format('>>> COPY $Project.OneInt  $CurrDir/ONEINT2',/, &
               '>>> COPY $Project.RunFile $CurrDir/RUNFIL2',//,'&scf',/)
          
          write(iunit,606)
          if(nalter.gt.0) then
            write(iunit,607) nalter
            write(iunit,608) ((ialter(i,j),i=1,3),j=1,nalter)
          endif
          write(iunit,206) ne,1,(nelecs-ne)/2,no,mol2,1,1,trim(project),6      
          write(iunit,208) ne,1,(nelecs-ne)/2,no,mol2,2,2,trim(project),7
          write(iunit,606)
          write(iunit,206) ne,3,(nelecs-ne)/2,no,mol2,3,1,trim(project),8
          write(iunit,606)
          write(iunit,206) ne-1,2,(nelecs-ne)/2,no,mol2,4,1,trim(project),9
          write(iunit,606)
          write(iunit,206) ne+1,2,(nelecs-ne)/2,no,mol2,5,1,trim(project),10

          close(unit=iunit)

          open(unit=iunit,file=trim(project)//'_SB.input')
          
          write(iunit,400)
          do j=1,numa
            nt=0
            do i=1,num
              if(atom(i).eq.elem(j)) nt=nt+1
            enddo
            nt0=nt
            if(nt.gt.0) then
              write(iunit,401) trim(basis(j))
              nt=0
              do i=1,num
                if(atom(i).eq.elem(j)) then
                  nt=nt+1
                  if(nt+nt0.le.9) then
                    write(iunit,402) nam(j),nt+nt0,x2(i,1),x2(i,2),x2(i,3)
                  else
                    write(iunit,403) nam(j),nt+nt0,x2(i,1),x2(i,2),x2(i,3)
                  endif
                endif
              enddo
              write(iunit,404)
            endif
          enddo
          write(iunit,409)
 409      format('>>> COPY $Project.OneInt  $CurrDir/ONEINT2',/, &
               '>>> COPY $Project.RunFile $CurrDir/RUNFIL2',//,'&scf',/)
          
          write(iunit,606)
          if(nalter.gt.0) then
            write(iunit,607) nalter
            write(iunit,608) ((ialter(i,j),i=1,3),j=1,nalter)
          endif
          write(iunit,406) ne,1,(nelecs-ne)/2,no,mol2,1,mol2,2,&
           1,trim(project),6,2,trim(project),7         
!          write(iunit,408) ne,1,(nelecs-ne)/2,no,mol2,2,2,trim(project),7
          write(iunit,606)
          write(iunit,206) ne,3,(nelecs-ne)/2,no,mol2,3,1,trim(project),8
          write(iunit,606)
          write(iunit,206) ne-1,2,(nelecs-ne)/2,no,mol2,4,1,trim(project),9
          write(iunit,606)
          write(iunit,206) ne+1,2,(nelecs-ne)/2,no,mol2,5,1,trim(project),10

          close(unit=iunit)
          
          open(unit=iunit,file=trim(project)//'_D.input')
          junit=12
          open(unit=junit,file=trim(project)//'_dimer.xyz')
          kunit=13
          open(unit=kunit,file=trim(project)//'_A.xyz')
          lunit=14
          open(unit=lunit,file=trim(project)//'_B.xyz')
          write(junit,102) 2*num
          write(kunit,102) num
          write(lunit,102) num
          
          write(iunit,200)
          do j=1,numa
            nt=0
            do i=1,num
              if(atom(i).eq.elem(j)) nt=nt+1
            enddo
            if(nt.gt.0) then
              write(iunit,201) trim(basis(j))
              nt=0
              do i=1,num
                if(atom(i).eq.elem(j)) then
                  nt=nt+1
                  if(nt.le.9) then
                    write(iunit,202) nam(j),nt,x1(i,1),x1(i,2),x1(i,3)
                  else
                    write(iunit,203) nam(j),nt,x1(i,1),x1(i,2),x1(i,3)
                  endif
                  write(junit,103) atom(i),x1(i,1),x1(i,2),x1(i,3)
                  write(kunit,103) atom(i),x1(i,1),x1(i,2),x1(i,3)
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
              write(iunit,201) trim(basis(j))
              nt=0
              do i=1,num
                if(atom(i).eq.elem(j)) then
                  nt=nt+1
                  if(nt+nt0.le.9) then
                    write(iunit,202) nam(j),nt+nt0,x2(i,1),x2(i,2),x2(i,3)
                  else
                    write(iunit,203) nam(j),nt+nt0,x2(i,1),x2(i,2),x2(i,3)
                  endif
                  write(junit,103) atom(i),x2(i,1),x2(i,2),x2(i,3)
                  write(lunit,103) atom(i),x2(i,1),x2(i,2),x2(i,3)
                endif
              enddo
              write(iunit,204)
            endif
          enddo
          
          write(iunit,210)
 210      format(/,'oneonly',/,'>>> COPY $Project.RunFile $CurrDir/RUNFILE')

          close(unit=iunit)
          close(unit=junit)
          close(unit=kunit)
          close(unit=lunit)
          
          open(unit=iunit,file=trim(project)//'_M.input')
          
          write(iunit,200)
          do j=1,numa
            nt=0
            do i=1,num
              if(atom(i).eq.elem(j)) nt=nt+1
            enddo
            if(nt.gt.0) then
              write(iunit,201) trim(basis(j))
              nt=0
              do i=1,num
                if(atom(i).eq.elem(j)) then
                  nt=nt+1
                  if(nt.le.9) then
                    write(iunit,202) nam(j),nt,x1(i,1),x1(i,2),x1(i,3)
                  else
                    write(iunit,203) nam(j),nt,x1(i,1),x1(i,2),x1(i,3)
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
              write(iunit,201) trim(basis(j))
              nt=0
              do i=1,num
                if(atom(i).eq.elem(j)) then
                  nt=nt+1
                  if(nt+nt0.le.9) then
                    write(iunit,202) nam(j),nt+nt0,x2(i,1),x2(i,2),x2(i,3)
                  else
                    write(iunit,203) nam(j),nt+nt0,x2(i,1),x2(i,2),x2(i,3)
                  endif
                endif
              enddo
              write(iunit,204)
            endif
          enddo
          
          write(iunit,211)
 211      format('>>> COPY $CurrDir/COMMONORB INPORB',//, &
               '&motra',/,'LumOrb',/,'frozen',/,'  0',/, &
               'deleted',/,' $DELETED',/,'ctonly',/,'kpq',//, &
               '>>> COPY $Project.RunFile $CurrDir/RUNFILE',/, &
               '>>> COPY $Project.OneInt  $CurrDir/ONEINT',/, &
               '>>> COPY $Project.TraOne  $CurrDir/TRAONE',/, &
               '>>> COPY $Project.ChVec1  $CurrDir/CHVEC1',/, &
               '>>> COPY $Project.ChRed   $CurrDir/CHRED',/, &
               '>>> COPY $Project.ChRst   $CurrDir/CHORST',/, &
               '>>> COPY $Project.ChMap   $CurrDir/CHOMAP',/, &
               '>>> COPY _CHMOT1 $CurrDir/_CHMOT1',/, &
               '>>> eval NPROCS = $MOLCAS_NPROCS - 1',/, &
               '>>> foreach L in (1 .. $NPROCS )',/, &
               '>>> shell cat tmp_$L/_CHMOT1 >> $CurrDir/_CHMOT1',/, &
               '>>> enddo')
          
          close(unit=iunit)
          
          open(unit=iunit,file=trim(project)//'_CB.input')

          write(iunit,212) trim(project)
 212      format('Project',/,1x,a,/,'Fragments',/,'  2',/,'  5 5',/, &
               'Threshold',/,' 1e-5',/,'Labels',/, &
               ' S0 S1 T1 D+ D- S0 S1 T1 D+ D-',/,'Energies')
          
          close(unit=iunit)
          
          open(unit=iunit,file=trim(project)//'_R.input')
          write(iunit,701)
701       format('xyzfiles',/,' fragA.xyz',/,' fragB.xyz',/, &
              ' 1 1',/,' 2 2',/,' 3 3')
          close(unit=iunit)
          
          open(unit=iunit,file=trim(project)//'_ss.run')

          if(aonly) then
            write(iunit,313) nr,memory,trim(project),nr,nr
          else
            write(iunit,213) nr,memory,trim(project),nr,nr
          endif
 213      format('#!/usr/bin/tcsh',/, &
              'setenv MOLCAS_NPROCS ',i3,/, &
              'setenv MOLCAS_MEM ',i5,/, &
              'setenv PROJECT "',a,'"',/, &
              'pymolcas -clean $PROJECT"_A.input" > $PROJECT".output"',/, &
              'pymolcas -clean $PROJECT"_B.input" >> $PROJECT".output"',/, &
              'pymolcas -clean $PROJECT"_D.input" >> $PROJECT".output"',/, &
              'common_basis < $PROJECT"_CB.input" >> $PROJECT".output"',/, &
              'setenv DELETED ` grep "Deleted orbitals in MOTRA"', &
              ' $PROJECT".output" | cut -b42- `',/, &
              'touch TRAINT',/, &
              'pymolcas -clean $PROJECT"_M.input" >> $PROJECT".output"',/, &
              'setenv OMP_NUM_THREADS ',i3,/, &
              'rdcho $MOLCAS_NPROCS >> $PROJECT".output"',/, &
              'rdtraint < $PROJECT"_CB.input" >> $PROJECT".output"',/, &
              'unsetenv DELETED',/, &
              '#mpirun -n',i3,' gronor $PROJECT"_dimer" >> $PROJECT".output"')
 313      format('#!/usr/bin/tcsh',/, &
              'setenv MOLCAS_NPROCS ',i3,/, &
              'setenv MOLCAS_MEM ',i5,/, &
              'setenv PROJECT "',a,'"',/, &
              'pymolcas -clean $PROJECT"_A.input" > $PROJECT".output"',/, &
              'pymolcas -clean $PROJECT"_B.input" >> $PROJECT".output"',/, &
              'pymolcas -clean $PROJECT"_D.input" >> $PROJECT".output"',/, &
              'common_basis < $PROJECT"_CB.input" >> $PROJECT".output"',/, &
              'setenv DELETED ` grep "Deleted orbitals in MOTRA"', &
              ' $PROJECT".output" | cut -b42- `',/, &
              'touch TRAINT',/, &
              'pymolcas -clean $PROJECT"_M.input" >> $PROJECT".output"',/, &
              'setenv OMP_NUM_THREADS ',i3,/, &
              'rdcho $MOLCAS_NPROCS >> $PROJECT".output"',/, &
              'rdtraint < $PROJECT"_CB.input" >> $PROJECT".output"',/, &
              'unsetenv DELETED',/, &
              '#mpirun -n',i3,' gronor $PROJECT"_dimer" >> $PROJECT".output"')
          
          close(unit=iunit)
          
          open(unit=iunit,file=trim(project)//'_sa.run')

          if(aonly) then
            write(iunit,513) nr,memory,trim(project),nr,nr
          else
            write(iunit,413) nr,memory,trim(project),nr,nr
          endif
 413      format('#!/usr/bin/tcsh',/, &
              'setenv MOLCAS_NPROCS ',i3,/, &
              'setenv MOLCAS_MEM ',i5,/, &
              'setenv PROJECT "',a,'"',/, &
              'pymolcas -clean $PROJECT"_SA.input" > $PROJECT".output"',/, &
              'pymolcas -clean $PROJECT"_SB.input" >> $PROJECT".output"',/, &
              'pymolcas -clean $PROJECT"_D.input" >> $PROJECT".output"',/, &
              'common_basis < $PROJECT"_CB.input" >> $PROJECT".output"',/, &
              'setenv DELETED ` grep "Deleted orbitals in MOTRA"', &
              ' $PROJECT".output" | cut -b42- `',/, &
              'touch TRAINT',/, &
              'pymolcas -clean $PROJECT"_M.input" >> $PROJECT".output"',/, &
              'setenv OMP_NUM_THREADS ',i3,/, &
              'rdcho $MOLCAS_NPROCS >> $PROJECT".output"',/, &
              'rdtraint < $PROJECT"_CB.input" >> $PROJECT".output"',/, &
              'unsetenv DELETED',/, &
              '#mpirun -n',i3,' gronor $PROJECT"_dimer" >> $PROJECT".output"')
 513      format('#!/usr/bin/tcsh',/, &
              'setenv MOLCAS_NPROCS ',i3,/, &
              'setenv MOLCAS_MEM ',i5,/, &
              'setenv PROJECT "',a,'"',/, &
              'pymolcas -clean $PROJECT"_SA.input" > $PROJECT".output"',/, &
              'pymolcas -clean $PROJECT"_SB.input" >> $PROJECT".output"',/, &
              'pymolcas -clean $PROJECT"_D.input" >> $PROJECT".output"',/, &
              'common_basis < $PROJECT"_CB.input" >> $PROJECT".output"',/, &
              'setenv DELETED ` grep "Deleted orbitals in MOTRA"', &
              ' $PROJECT".output" | cut -b42- `',/, &
              'touch TRAINT',/, &
              'pymolcas -clean $PROJECT"_M.input" >> $PROJECT".output"',/, &
              'setenv OMP_NUM_THREADS ',i3,/, &
              'rdcho $MOLCAS_NPROCS >> $PROJECT".output"',/, &
              'rdtraint < $PROJECT"_CB.input" >> $PROJECT".output"',/, &
              'unsetenv DELETED',/, &
              '#mpirun -n',i3,' gronor $PROJECT"_dimer" >> $PROJECT".output"')
          
          close(unit=iunit)
          
          open(unit=iunit,file=trim(project)//'_ssr.run')

          if(aonly) then
            write(iunit,1313) nr,memory,trim(project),nr,nr
          else
            write(iunit,1213) nr,memory,trim(project),nr,nr
          endif
1213      format('#!/usr/bin/tcsh',/, &
              'setenv MOLCAS_NPROCS ',i3,/, &
              'setenv MOLCAS_MEM ',i5,/, &
              'setenv PROJECT "',a,'"',/, &
              'pymolcas -clean $PROJECT"_A.input" > $PROJECT".output"',/, &
              'cp $PROJECT"_A.xyz" "fragA.xyz"',/, &
              'cp $PROJECT"_B.xyz" "fragB.xyz"',/, &
              'cp RUNFIL1 RUNFILE',/, &
              'foreach n ( 1 2 3 4 5 )',/, &
              ' cp "INPORB.1_"$n INPORB',/, &
              ' rotharm3 < $PROJECT"_R.input"',/, &
              ' cp ROTORB "INPORB.2_"$n',/, &
              ' cp transrot.xyz $PROJECT"_2_"$n".xyz"',/, &
              'end',/,&
              'pymolcas -clean $PROJECT"_D.input" >> $PROJECT".output"',/, &
              'common_basis < $PROJECT"_CB.input" >> $PROJECT".output"',/, &
              'setenv DELETED ` grep "Deleted orbitals in MOTRA"', &
              ' $PROJECT".output" | cut -b42- `',/, &
              'touch TRAINT',/, &
              'pymolcas -clean $PROJECT"_M.input" >> $PROJECT".output"',/, &
              'setenv OMP_NUM_THREADS ',i3,/, &
              'rdcho $MOLCAS_NPROCS >> $PROJECT".output"',/, &
              'rdtraint < $PROJECT"_CB.input" >> $PROJECT".output"',/, &
              'unsetenv DELETED',/, &
              '#mpirun -n',i3,' gronor $PROJECT"_dimer" >> $PROJECT".output"')
1313      format('#!/usr/bin/tcsh',/, &
              'setenv MOLCAS_NPROCS ',i3,/, &
              'setenv MOLCAS_MEM ',i5,/, &
              'setenv PROJECT "',a,'"',/, &
              'pymolcas -clean $PROJECT"_A.input" > $PROJECT".output"',/, &
              'cp $PROJECT"_A.xyz" "fragA.xyz"',/, &
              'cp $PROJECT"_B.xyz" "fragB.xyz"',/, &
              'cp RUNFIL1 RUNFILE',/, &
              'foreach n ( 1 2 3 4 5 )',/, &
              ' cp "INPORB.1_"$n INPORB',/, &
              ' rotharm3 < $PROJECT"_R.input"',/, &
              ' cp ROTORB "INPORB.2_"$n',/, &
              ' cp transrot.xyz $PROJECT"_2_"$n".xyz"',/, &
              'end',/,&
              'pymolcas -clean $PROJECT"_D.input" >> $PROJECT".output"',/, &
              'common_basis < $PROJECT"_CB.input" >> $PROJECT".output"',/, &
              'setenv DELETED ` grep "Deleted orbitals in MOTRA"', &
              ' $PROJECT".output" | cut -b42- `',/, &
              'touch TRAINT',/, &
              'pymolcas -clean $PROJECT"_M.input" >> $PROJECT".output"',/, &
              'setenv OMP_NUM_THREADS ',i3,/, &
              'rdcho $MOLCAS_NPROCS >> $PROJECT".output"',/, &
              'rdtraint < $PROJECT"_CB.input" >> $PROJECT".output"',/, &
              'unsetenv DELETED',/, &
              '#mpirun -n',i3,' gronor $PROJECT"_dimer" >> $PROJECT".output"')
          
          close(unit=iunit)
          
          open(unit=iunit,file=trim(project)//'_sar.run')

          if(aonly) then
            write(iunit,1513) nr,memory,trim(project),nr,nr
          else
            write(iunit,1413) nr,memory,trim(project),nr,nr
          endif
1413      format('#!/usr/bin/tcsh',/, &
              'setenv MOLCAS_NPROCS ',i3,/, &
              'setenv MOLCAS_MEM ',i5,/, &
              'setenv PROJECT "',a,'"',/, &
              'pymolcas -clean $PROJECT"_SA.input" > $PROJECT".output"',/, &
              'cp $PROJECT"_A.xyz" "fragA.xyz"',/, &
              'cp $PROJECT"_B.xyz" "fragB.xyz"',/, &
              'cp RUNFIL1 RUNFILE',/, &
              'foreach n ( 1 2 3 4 5 )',/, &
              ' cp "INPORB.1_"$n INPORB',/, &
              ' rotharm3 < $PROJECT"_R.input"',/, &
              ' cp ROTORB "INPORB.2_"$n',/, &
              ' cp transrot.xyz $PROJECT"_2_"$n".xyz"',/, &
              'end',/,&
              'pymolcas -clean $PROJECT"_D.input" >> $PROJECT".output"',/, &
              'common_basis < $PROJECT"_CB.input" >> $PROJECT".output"',/, &
              'setenv DELETED ` grep "Deleted orbitals in MOTRA"', &
              ' $PROJECT".output" | cut -b42- `',/, &
              'touch TRAINT',/, &
              'pymolcas -clean $PROJECT"_M.input" >> $PROJECT".output"',/, &
              'setenv OMP_NUM_THREADS ',i3,/, &
              'rdcho $MOLCAS_NPROCS >> $PROJECT".output"',/, &
              'rdtraint < $PROJECT"_CB.input" >> $PROJECT".output"',/, &
              'unsetenv DELETED',/, &
              '#mpirun -n',i3,' gronor $PROJECT"_dimer" >> $PROJECT".output"')
1513      format('#!/usr/bin/tcsh',/, &
              'setenv MOLCAS_NPROCS ',i3,/, &
              'setenv MOLCAS_MEM ',i5,/, &
              'setenv PROJECT "',a,'"',/, &
              'pymolcas -clean $PROJECT"_SA.input" > $PROJECT".output"',/, &
              'cp $PROJECT"_A.xyz" "fragA.xyz"',/, &
              'cp $PROJECT"_B.xyz" "fragB.xyz"',/, &
              'cp RUNFIL1 RUNFILE',/, &
              'foreach n ( 1 2 3 4 5 )',/, &
              ' cp "INPORB.1_"$n INPORB',/, &
              ' rotharm3 < $PROJECT"_R.input"',/, &
              ' cp ROTORB "INPORB.2_"$n',/, &
              ' cp transrot.xyz $PROJECT"_2_"$n".xyz"',/, &
              'end',/,&
              'pymolcas -clean $PROJECT"_D.input" >> $PROJECT".output"',/, &
              'common_basis < $PROJECT"_CB.input" >> $PROJECT".output"',/, &
              'setenv DELETED ` grep "Deleted orbitals in MOTRA"', &
              ' $PROJECT".output" | cut -b42- `',/, &
              'touch TRAINT',/, &
              'pymolcas -clean $PROJECT"_M.input" >> $PROJECT".output"',/, &
              'setenv OMP_NUM_THREADS ',i3,/, &
              'rdcho $MOLCAS_NPROCS >> $PROJECT".output"',/, &
              'rdtraint < $PROJECT"_CB.input" >> $PROJECT".output"',/, &
              'unsetenv DELETED',/, &
              '#mpirun -n',i3,' gronor $PROJECT"_dimer" >> $PROJECT".output"')
          
          close(unit=iunit)
          
          open(unit=iunit,file=trim(project)//'_scf.run')

            write(iunit,813) nr,memory,trim(project)

 813      format('#!/usr/bin/tcsh',/, &
              'setenv MOLCAS_NPROCS ',i3,/, &
              'setenv MOLCAS_MEM ',i5,/, &
              'setenv PROJECT "',a,'"',/, &
              'pymolcas -clean $PROJECT"_SCF.input" > $PROJECT"_SCF.output"') 
          
          close(unit=iunit)
          
          open(unit=iunit,file=trim(project)//'_dimer.inp')

          write(iunit,214)
 214      format('MEBFs 7',/, &
               '  1  2  1  2  3  4  5',/, &
               '  6  6  7  7  8 10  9',/, &
               'Threshold 1e-4',/, &
               'Timings 4',/, &
               'Size 1',/, &
               'Print medium')
          
          close(unit=iunit)
          operation='&'
        endif
        
      enddo
 1009 continue

      end
      
      subroutine molcas_rotate(v,w,g,x,y)

        implicit none

        real (kind=8) :: v(3),w(3),x(3),y(3),xx(3),t(3,3)
        real (kind=8) :: pi,r,a,b,ca,cb,cg,sa,sb,sg,g
        integer :: i
        real (kind=8) , parameter :: small=1.0d-24

        if(abs(w(2)-v(2)).lt.small) then
          if(abs(w(1)-v(1)).lt.small) then
            a=0.0d0
          else
            if(w(1)-v(1).gt.0.0d0) then
              a=2.0d0*datan(1.0d0)
            else
              a=(-2.0d0)*datan(1.0d0)
            endif
          endif
        else
          a=atan(abs(w(1)-v(1))/abs(w(2)-v(2)))
          pi=4.0d0*atan(1.0d0)
          if(w(1)-v(1).gt.0.0d0.and.w(2)-v(2).lt.0.0d0) a=pi-a
          if(w(1)-v(1).lt.0.0d0.and.w(2)-v(2).gt.0.0d0) a=-a
          if(w(1)-v(1).lt.0.0d0.and.w(2)-v(2).lt.0.0d0) a=pi+a
        endif
        r=0.0d0
        do i=1,3
          r=r+(w(i)-v(i))**2
          xx(i)=x(i)-v(i)
        enddo
        if(r.lt.small) then
          y(1)=x(1)
          y(2)=x(2)
          y(3)=x(3)
          return
        endif
        b=acos((w(3)-v(3))/sqrt(r))
        sa=sin(a)
        ca=cos(a)
        sb=sin(b)
        cb=cos(b)
        sg=sin(g)
        cg=cos(g)
        t(1,1)=ca*ca*cg-sa*ca*cb*sg+sa*ca*cb*sg+sa*sa*cb*cb*cg+sa*sa*sb*sb
        t(1,2)=(-sa)*ca*cg-ca*ca*cb*sg-sa*sa*cb*sg+sa*ca*cb*cb*cg+sa*ca*sb*sb
        t(1,3)=ca*sb*sg-sa*sb*cb*cg+sa*sb*cb
        t(2,1)=(-sa)*ca*cg+sa*sa*cb*sg+ca*ca*cb*sg+sa*ca*cb*cb*cg+sa*ca*sb*sb
        t(2,2)=sa*sa*cg+sa*ca*cb*sg-sa*ca*cb*sg+ca*ca*cb*cb*cg+ca*ca*sb*sb
        t(2,3)=(-sa)*sb*sg-ca*sb*cb*cg+ca*sb*cb
        t(3,1)=(-ca)*sb*sg-sa*sb*cb*cg+sa*sb*cb
        t(3,2)=sa*sb*sg-ca*sb*cb*cg+ca*sb*cb
        t(3,3)=sb*sb*cg+cb*cb
        do i=1,3
          y(i)=xx(1)*t(i,1)+xx(2)*t(i,2)+xx(3)*t(i,3)+v(i)
        enddo
        return
      end subroutine molcas_rotate
