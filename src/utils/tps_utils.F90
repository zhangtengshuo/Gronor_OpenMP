integer function atom_number(element)

  implicit none

  character*2 :: element

  character*2 :: aname(0:105),bname(0:105)

  data aname / '  ', &
      ' H','He','Li','Be',' B',' C',' N',' O',' F','Ne', &
      'Na','Mg','Al','Si',' P',' S','Cl','Ar',' K','Ca', &
      'Sc','Ti',' V','Cr','Mn','Fe','Co','Ni','Cu','Zn', &
      'Ga','Ge','As','Se','Br',' R','Rb','Sr',' Y','Zr', &
      'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn', &
      'Sb','Te',' I','Xe','Cs','Ba','La','Ce','Pr','Nd', &
      'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb', &
      'Lu','Hf','Ta',' W','Re','Os','Ir','Pt','Au','Hg', &
      'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th', &
      'Pa',' U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm', &
      'Md','No','Lr','Rf','Ha' /
  data bname / '  ', &
      ' H','HE','LI','BE',' B',' C',' N',' O',' F','NE', &
      'NA','MG','AL','SI',' P',' S','CL','AR',' K','CA', &
      'SC','TI',' V','CR','MN','FE','CO','NI','CU','ZN', &
      'GA','GE','AS','SE','BR',' R','RB','SR',' Y','ZR', &
      'NB','MO','TC','RU','RH','PD','AG','CD','IN','SN', &
      'SB','TE',' I','XE','CS','BA','LA','CE','PR','ND', &
      'PM','SM','EU','GD','TB','DY','HO','ER','TM','YB', &
      'LU','HF','TA',' W','RE','OS','IR','PT','AU','HG', &
      'TL','PB','BI','PO','AT','RN','FR','RA','AC','TH', &
      'PA',' U','Np','PU','AM','CM','BK','CF','ES','FM', &
      'MD','NO','LR','RF','HA' /

  integer i

  do i=1,105
    if(element.eq.aname(i).or.element.eq.bname(i)) then
      atom_number=i
      return
    endif
  enddo

  atom_number=0
  if(element.eq.'1H') atom_number=1
  if(element.eq.'2H') atom_number=1
  if(element.eq.'3H') atom_number=1
  if(element.eq.'4H') atom_number=1

  return
end function atom_number

real (kind=8) function atom_radius(number)

  implicit none

  integer :: number

  real (kind=8) :: radius(0:105)

  data radius / 99999.99, &
      0.35, 1.22, 1.23, 0.89, 0.88, 0.77, 0.70, 0.66, 0.58, 1.60, &
      1.40, 1.36, 1.25, 1.17, 1.10, 1.04, 0.99, 1.91, 2.03, 1.74, &
      1.44, 1.32, 1.22, 1.19, 1.17,1.165, 1.16, 1.15, 1.17, 1.25, &
      1.25, 1.22, 1.21, 1.17, 1.14, 1.98, 2.22, 1.92, 1.62, 1.45, &
      1.34, 1.29, 1.27, 1.24, 1.25, 1.28, 1.34, 1.41, 1.50, 1.40, &
      1.41, 1.37, 1.33, 2.09, 2.35, 1.98, 1.69, 1.65, 1.65, 1.64, &
      1.65, 1.66, 1.65, 1.61, 1.59, 1.59, 1.58, 1.57, 1.56, 1.56, &
      1.56, 1.44, 1.34, 1.30, 1.28, 1.26, 1.26, 1.29, 1.34, 1.44, &
      1.55, 1.54, 1.52, 1.53, 1.50, 2.20, 3.24, 2.68, 2.25, 2.16, &
      1.93, 1.66, 1.57, 1.81, 2.21, 1.43, 1.42, 1.40, 1.39, 1.38, &
      1.37, 1.36, 1.34, 1.30, 1.30 /

  atom_radius=1.1d-01*radius(number)

  return
end function atom_radius

character (len=255) function atom_color(number)

  implicit none
  integer :: number

  integer :: indexc(0:105)

  data indexc / 12, &
       7,  5, 14, 12, 13,  0,  1,  2,  6, 12, &
      13, 15,  9,  6,  8,  3, 13, 12, 16, 16, &
      12,  9, 12,  9,  9,  8, 12, 10, 10, 10, &
      12, 12, 12, 12, 10, 12, 12, 12, 12, 12, &
      12, 12, 12, 12, 12, 12, 16, 12, 12, 12, &
      12, 12, 11, 12, 12,  8, 12, 12, 12, 12, &
      12, 12, 12, 12, 12, 12, 12, 12, 12, 12, &
      12, 12, 12, 12, 12, 12, 12, 16, 17, 16, &
      12, 12, 12, 12, 12, 12, 12, 12, 12, 12, &
      12, 12, 12, 12, 12, 12, 12, 12, 12, 12, &
      12, 12, 12, 12, 12 /

  atom_color='white '
  if(indexc(number).eq.0) atom_color='LightGrey '
  if(indexc(number).eq.1) atom_color='SkyBlue '
  if(indexc(number).eq.2) atom_color='Red '
  if(indexc(number).eq.3) atom_color='Yellow '
  if(indexc(number).eq.4) atom_color='White '
  if(indexc(number).eq.5) atom_color='Pink '
  if(indexc(number).eq.6) atom_color='Goldenrod '
  if(indexc(number).eq.7) atom_color='Blue '
  if(indexc(number).eq.8) atom_color='Orange '
  if(indexc(number).eq.9) atom_color='Gray85 '
  if(indexc(number).eq.10) atom_color='Brown '
  if(indexc(number).eq.11) atom_color='Purple '
  if(indexc(number).eq.12) atom_color='SpicyPink '
  if(indexc(number).eq.13) atom_color='Green '
  if(indexc(number).eq.14) atom_color='Firebrick '
  if(indexc(number).eq.15) atom_color='DarkGreen '
  if(indexc(number).eq.16) atom_color='Silver '
  if(indexc(number).eq.17) atom_color='Gold '

  return
end function atom_color

character*2 function element(number)

  implicit none
  integer :: number

  character*2 aname(0:116)

  data aname / '  ', &
      ' H','He','Li','Be',' B',' C',' N',' O',' F','Ne', &
      'Na','Mg','Al','Si',' P',' S','Cl','Ar',' K','Ca', &
      'Sc','Ti',' V','Cr','Mn','Fe','Co','Ni','Cu','Zn', &
      'Ga','Ge','As','Se','Br',' R','Rb','Sr',' Y','Zr', &
      'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn', &
      'Sb','Te',' I','Xe','Cs','Ba','La','Ce','Pr','Nd', &
      'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb', &
      'Lu','Hf','Ta',' W','Re','Os','Ir','Pt','Au','Hg', &
      'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th', &
      'Pa',' U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm', &
      'Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds', &
      'Rg','Cn','  ','Fl','  ','Lv'/

  element=aname(number)

  return
end function element

subroutine tps_rotate(v,w,g,x,y)

  ! $Id: util_md.F 19707 2010-10-29 17:59:36Z d3y133 $

  implicit none

  real (kind=8) :: v(3),w(3),x(3),y(3),xx(3),t(3,3)
  real (kind=8) :: small,pi,r,a,b,ca,cb,cg,sa,sb,sg,g
  integer :: i
  parameter (small=1.0d-24)

  !     rotation with angle g around vector from v to w
  !     of point x giving result in y

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
end subroutine tps_rotate

real (kind=8) function angl(x,y,z)

  implicit none

  real (kind=8) :: x(3),y(3),z(3)

  real (kind=8) :: xy(3),zy(3),rxy,rzy,phi
  integer :: i

  rxy=0.0d0
  rzy=0.0d0
  do i=1,3
    xy(i)=x(i)-y(i)
    zy(i)=z(i)-y(i)
  enddo
  rxy=xy(1)*xy(1)+xy(2)*xy(2)+xy(3)*xy(3)
  rzy=zy(1)*zy(1)+zy(2)*zy(2)+zy(3)*zy(3)
  phi=(xy(1)*zy(1)+xy(2)*zy(2)+xy(3)*zy(3))/sqrt(rxy*rzy)
  if(phi.lt.-1.0d0) phi=-1.0d0
  if(phi.gt.1.0d0) phi=1.0d0
  angl=acos(phi)

  return
end function angl

subroutine povinc(lfn,xmin,xmax,ymin,ymax,zmin,zmax)

  implicit none

  integer :: lfn
  real (kind=8) :: xmin,xmax,ymin,ymax,zmin,zmax

  open(unit=lfn,file='camera.inc',form='formatted',status='new',err=9)
  write(lfn,1000) 0.0d0,0.0d0,-1.0d1*max(abs(xmax),abs(ymax),abs(zmax), &
      abs(xmin),abs(ymin),abs(zmin)),5.0d-1*(xmax+xmin),5.0d-1*(ymax+ymin),5.0d-1*(zmax+zmin)
1000 format('camera {',/,' location <',f12.6,',',f12.6,',',f12.6,'>',/, &
      ' look_at <',f12.6,',',f12.6,',',f12.6,'>',/,' angle 20.0',/,'}')
  write(lfn,1001)  0.0d1, 0.0d1,-5.0d1, 1.0d0,1.0d0,1.0d0
  write(lfn,1001) -1.0d1, 2.0d1,-1.0d1, 1.0d0,1.0d0,1.0d0
  write(lfn,1001)  1.0d1, 2.0d1,-1.0d1, 1.0d0,1.0d0,1.0d0
  write(lfn,1001)  0.0d1, 1.0d1,-2.0d1, 1.0d0,1.0d0,1.0d0
1001 format('light_source { <',f12.6,',',f12.6,',',f12.6, &
      '> color rgb <',f4.2,',',f4.2,',',f4.2,'> }')
  write(lfn,1002) 0.0d0,0.0d0,0.0d0
1002 format('background { color rgb <',f4.2,',',f4.2,',',f4.2,'> }')
  close(unit=lfn)

9 continue
  open(unit=lfn,file='colors.inc',form='formatted',status='new',err=99)
  write(lfn,2000)
2000 format('#declare Colors_Inc_Temp = version ;',/,'#version 2.0 ;')
  write(lfn,2001) 'Gray05',0.05,0.05,0.05
  write(lfn,2001) 'Gray05',0.05,0.05,0.05
  write(lfn,2001) 'Gray10',0.10,0.10,0.10
  write(lfn,2001) 'Gray15',0.15,0.15,0.15
  write(lfn,2001) 'Gray20',0.20,0.20,0.20
  write(lfn,2001) 'Gray25',0.25,0.25,0.25
  write(lfn,2001) 'Gray30',0.30,0.30,0.30
  write(lfn,2001) 'Gray35',0.35,0.35,0.35
  write(lfn,2001) 'Gray40',0.40,0.40,0.40
  write(lfn,2001) 'Gray45',0.45,0.45,0.45
  write(lfn,2001) 'Gray50',0.50,0.50,0.50
  write(lfn,2001) 'Gray55',0.55,0.55,0.55
  write(lfn,2001) 'Gray60',0.60,0.60,0.60
  write(lfn,2001) 'Gray65',0.65,0.65,0.65
  write(lfn,2001) 'Gray70',0.70,0.70,0.70
  write(lfn,2001) 'Gray75',0.75,0.75,0.75
  write(lfn,2001) 'Gray80',0.80,0.80,0.80
  write(lfn,2001) 'Gray85',0.85,0.85,0.85
  write(lfn,2001) 'Gray90',0.90,0.90,0.90
  write(lfn,2001) 'Gray95',0.95,0.95,0.95
  write(lfn,2001) 'DimGray',0.329412,0.329412,0.329412
  write(lfn,2001) 'DimGrey',0.329412,0.329412,0.329412
  write(lfn,2001) 'Gray',0.752941,0.752941,0.752941
  write(lfn,2001) 'Grey',0.752941,0.752941,0.752941
  write(lfn,2001) 'LightGray',0.658824,0.658824,0.658824
  write(lfn,2001) 'LightGrey',0.658824,0.658824,0.658824
  write(lfn,2001) 'VLightGrey',0.80,0.80,0.80
  write(lfn,2001) 'White',1.0,1.0,1.0
  write(lfn,2001) 'Red',1.0,0.0,0.0
  write(lfn,2001) 'Green',0.0,1.0,0.0
  write(lfn,2001) 'Blue',0.0,0.0,1.0
  write(lfn,2001) 'Yellow',1.0,1.0,0.0
  write(lfn,2001) 'Cyan',0.0,1.0,1.0
  write(lfn,2001) 'Magenta',1.0,0.0,1.0
  write(lfn,2001) 'Black',0.0,0.0,0.0
  write(lfn,2001) 'Aquamarine',0.439216,0.858824,0.576471
  write(lfn,2001) 'BlueViolet',0.62352,0.372549,0.623529
  write(lfn,2001) 'Brown',0.647059,0.164706,0.164706
  write(lfn,2001) 'CadetBlue',0.372549,0.623529,0.623529
  write(lfn,2001) 'Coral',1.0,0.498039,0.0
  write(lfn,2001) 'CornflowerBlue',0.258824,0.258824,0.435294
  write(lfn,2001) 'DarkGreen',0.184314,0.309804,0.184314
  write(lfn,2001) 'DarkOliveGreen',0.309804,0.309804,0.184314
  write(lfn,2001) 'DarkOrchid',0.6,0.196078,0.8
  write(lfn,2001) 'DarkSlateBlue',0.419608,0.137255,0.556863
  write(lfn,2001) 'DarkSlateGray',0.184314,0.309804,0.309804
  write(lfn,2001) 'DarkSlateGrey',0.184314,0.309804,0.309804
  write(lfn,2001) 'DarkTurquoise',0.439216,0.576471,0.858824
  write(lfn,2001) 'Firebrick',0.556863,0.137255,0.137255
  write(lfn,2001) 'ForestGreen',0.137255,0.556863,0.137255
  write(lfn,2001) 'Gold',0.8,0.498039,0.196078
  write(lfn,2001) 'Goldenrod',0.858824,0.858824,0.439216
  write(lfn,2001) 'GreenYellow',0.576471,0.858824,0.439216
  write(lfn,2001) 'IndianRed',0.309804,0.184314,0.184314
  write(lfn,2001) 'Khaki',0.623529,0.623529,0.372549
  write(lfn,2001) 'LightBlue',0.74902,0.847059,0.847059
  write(lfn,2001) 'LightSteelBlue',0.560784,0.560784,0.737255
  write(lfn,2001) 'LimeGreen',0.196078,0.8,0.196078
  write(lfn,2001) 'Maroon',0.556863,0.137255,0.419608
  write(lfn,2001) 'MediumAquamarine',0.196078,0.8,0.6
  write(lfn,2001) 'MediumBlue',0.196078,0.196078,0.8
  write(lfn,2001) 'MediumForestGreen',0.419608,0.556863,0.137255
  write(lfn,2001) 'MediumGoldenrod',0.917647,0.917647,0.678431
  write(lfn,2001) 'MediumOrchid',0.576471,0.439216,0.858824
  write(lfn,2001) 'MediumSeaGreen',0.258824,0.435294,0.258824
  write(lfn,2001) 'MediumSlateBlue',0.498039,1.0,0.0
  write(lfn,2001) 'MediumSpringGreen',0.498039,1.0,0.0
  write(lfn,2001) 'MediumTurquoise',0.439216,0.858824,0.858824
  write(lfn,2001) 'MediumVioletRed',0.858824,0.439216,0.576471
  write(lfn,2001) 'MidnightBlue',0.184314,0.184314,0.309804
  write(lfn,2001) 'Navy',0.137255,0.137255,0.556863
  write(lfn,2001) 'NavyBlue',0.137255,0.137255,0.556863
  write(lfn,2001) 'Orange',1,0.5,0.0
  write(lfn,2001) 'OrangeRed',1.0,0.498039,0.0
  write(lfn,2001) 'Orchid',0.858824,0.439216,0.858824
  write(lfn,2001) 'PaleGreen',0.560784,0.737255,0.560784
  write(lfn,2001) 'Pink',0.737255,0.560784,0.560784
  write(lfn,2001) 'Plum',0.917647,0.678431,0.917647
  write(lfn,2001) 'Salmon',0.435294,0.258824,0.258824
  write(lfn,2001) 'SeaGreen',0.137255,0.556863,0.419608
  write(lfn,2001) 'Sienna',0.556863,0.419608,0.137255
  write(lfn,2001) 'SkyBlue',0.196078,0.6,0.8
  write(lfn,2001) 'SlateBlue',0.0,0.498039,1.0
  write(lfn,2001) 'SpringGreen',0.0,1.0,0.498039
  write(lfn,2001) 'SteelBlue',0.137255,0.419608,0.556863
  write(lfn,2001) 'Tan',0.858824,0.576471,0.439216
  write(lfn,2001) 'Thistle',0.847059,0.74902,0.847059
  write(lfn,2001) 'Turquoise',0.678431,0.917647,0.917647
  write(lfn,2001) 'Violet',0.309804,0.184314,0.309804
  write(lfn,2001) 'VioletRed',0.8,0.196078,0.6
  write(lfn,2001) 'Wheat',0.847059,0.847059,0.74902
  write(lfn,2001) 'YellowGreen',0.6,0.8,0.196078
  write(lfn,2001) 'SummerSky',0.22,0.69,0.87
  write(lfn,2001) 'RichBlue',0.35,0.35,0.67
  write(lfn,2001) 'Brass',0.71,0.65,0.26
  write(lfn,2001) 'Copper',0.72,0.45,0.20
  write(lfn,2001) 'Bronze',0.55,0.47,0.14
  write(lfn,2001) 'Bronze2',0.65,0.49,0.24
  write(lfn,2001) 'Silver',0.90,0.91,0.98
  write(lfn,2001) 'BrightGold',0.85,0.85,0.10
  write(lfn,2001) 'OldGold',0.81,0.71,0.23
  write(lfn,2001) 'Feldspar',0.82,0.57,0.46
  write(lfn,2001) 'Quartz',0.85,0.85,0.95
  write(lfn,2001) 'Mica',0.0,0.0,0.0
  write(lfn,2001) 'NeonPink',1.00,0.43,0.78
  write(lfn,2001) 'DarkPurple',0.53,0.12,0.47
  write(lfn,2001) 'NeonBlue',0.30,0.30,1.00
  write(lfn,2001) 'CoolCopper',0.85,0.53,0.10
  write(lfn,2001) 'MandarinOrange',0.89,0.47,0.20
  write(lfn,2001) 'LightWood',0.91,0.76,0.65
  write(lfn,2001) 'MediumWood',0.65,0.50,0.39
  write(lfn,2001) 'DarkWood',0.52,0.37,0.26
  write(lfn,2001) 'SpicyPink',1.00,0.11,0.68
  write(lfn,2001) 'SemiSweetChoc',0.42,0.26,0.15
  write(lfn,2001) 'BakersChoc',0.36,0.20,0.09
  write(lfn,2001) 'Flesh',0.96,0.80,0.69
  write(lfn,2001) 'NewTan',0.92,0.78,0.62
  write(lfn,2001) 'NewMidnightBlue',0.00,0.00,0.61
  write(lfn,2001) 'VeryDarkBrown',0.35,0.16,0.14
  write(lfn,2001) 'DarkBrown',0.36,0.25,0.20
  write(lfn,2001) 'DarkTan',0.59,0.41,0.31
  write(lfn,2001) 'GreenCopper',0.32,0.49,0.46
  write(lfn,2001) 'DkGreenCopper',0.29,0.46,0.43
  write(lfn,2001) 'DustyRose',0.52,0.39,0.39
  write(lfn,2001) 'HuntersGreen',0.13,0.37,0.31
  write(lfn,2001) 'Scarlet',0.55,0.09,0.09
  write(lfn,2002) 'Clear',1.0,1.0,1.0,1.0
2001 format('#declare ',a,' = color red ',f8.6,' green ',f8.6,' blue ',f8.6,' ;')
2002 format('#declare ',a,' = color red ',f8.6,' green ',f8.6,' blue ',f8.6,' filter ',f8.6,' ;')
  write(lfn,2003)
2003 format('#declare Plane_Map = 0 ;',/,'#declare Sphere_Map = 1 ;',/, &
      '#declare Cylinder_Map = 2 ;',/,'#declare Torus_Map = 5 ;',/, &
      '#declare Bi   = 2 ;',/,'#declare Norm = 4 ;',/,'#version Colors_Inc_Temp ;')
  close(unit=lfn)

99 continue
  open(unit=lfn,file='plane.inc',form='formatted',status='new',err=999)
  write(lfn,3001)
3001 format('plane{ <0,1,0>,-1 pigment { White } }')
  close(unit=lfn)

999 continue

  return
end subroutine povinc

real (kind=8) function torsion(a,b,c,d)

  implicit none

  real (kind=8) :: dummy

  real (kind=8) :: a(3),b(3),c(3),d(3)

  dummy=a(1)+b(1)+c(1)+d(1)

  torsion=0.0d0

  return
end function torsion

subroutine argos_jacobi(a,n,na,d,v,nrot)

  !     compute eigenvectors and eigenvalues for real symmetric
  !     matrix using the Jacobi diagonalization

  implicit none

  external :: tps_abort

  integer :: nrmax
  parameter(nrmax=100)

  real (kind=8) :: zero,half,one,two
  parameter(zero=0.0d0)
  parameter(half=0.5d0)
  parameter(one=1.0d0)
  parameter(two=2.0d0)

  integer :: n,na,nrot
  real (kind=8) :: a(na,na),d(na),v(na,na)
  real (kind=8) :: at,b,dma,q

  integer :: i,j,k,l
  real (kind=8) :: c,s,t,sum,temp

  do i=1,n
    do j=1,n
      v(i,j)=zero
    enddo
    v(i,i)=one
    d(i)=a(i,i)
  enddo

  nrot=0
  do l=1,nrmax
    nrot=nrot+1
    sum=zero
    do i=1,n-1
      do j=i+1,n
        sum=sum+abs(a(i,j))
      enddo
    enddo
    if(sum.eq.zero) then
      do i=1,n-1
        do j=i+1,n
          if(d(i).gt.d(j)) then
            temp=d(i)
            d(i)=d(j)
            d(j)=temp
            do k=1,n
              temp=v(k,i)
              v(k,i)=v(k,j)
              v(k,j)=temp
            enddo
          endif
        enddo
      enddo

      return
    endif
    do j=2,n
      do i=1,j-1
        b=a(i,j)
        if(abs(b).gt.zero) then
          dma=d(j)-d(i)
          if(abs(dma)+abs(b).le.abs(dma)) then
            t=b/dma
          else
            q=half*dma/b
            t=sign(one/(abs(q)+sqrt(one+q*q)),q)
          endif
          c=one/sqrt(t*t+one)
          s=t*c
          a(i,j)=zero
          do k=1,i-1
            at=c*a(k,i)-s*a(k,j)
            a(k,j)=s*a(k,i)+c*a(k,j)
            a(k,i)=at
          enddo
          do k=i+1,j-1
            at=c*a(i,k)-s*a(k,j)
            a(k,j)=s*a(i,k)+c*a(k,j)
            a(i,k)=at
          enddo
          do k=j+1,n
            at=c*a(i,k)-s*a(j,k)
            a(j,k)=s*a(i,k)+c*a(j,k)
            a(i,k)=at
          enddo
          do k=1,n
            at=c*v(k,i)-s*v(k,j)
            v(k,j)=s*v(k,i)+c*v(k,j)
            v(k,i)=at
          enddo
          at=c*c*d(i)+s*s*d(j)-two*c*s*b
          d(j)=s*s*d(i)+c*c*d(j)+two*c*s*b
          d(i)=at
        endif
      enddo
    enddo
  enddo

  call tps_abort('argos_jacobi: maximum iterations reached',0)

  return
end subroutine argos_jacobi

subroutine suser(user)

  implicit none

#ifdef USE_POSIXF
  external :: pxfgetlogin
#else
  external :: getenv
#endif

  character*32 user

#ifdef USE_POSIXF
  integer (kind=4) :: ilen,ierror
#endif

#ifdef USE_POSIXF
  call pxfgetlogin(user,ilen,ierror)
#else
  !      call getlog(user)
  call getenv("USER",user)
#endif

  return
end subroutine suser

subroutine swatch(today,now)

  implicit none

#if defined(LINUX) && !defined(IBM)
  external :: linux_date,linux_time
#endif
  external :: tps_abort

  character*10 today,now

  !#if defined(LINUX)
  !      character*26 string
  !#endif
#if defined(IBM)
  character*26 dstring,tstring,tzone
  integer :: dtvalue(8)
#endif
#if defined(KSR)
  integer :: time
  character*24 ctime,string
#endif
#if defined(SP1) || defined(CRAY_T3D) || defined(CRAY_T3E) || defined(SOLARIS)
  character*26 string
#endif
#if defined(SGI)
  character*9 string
#endif

  today='00/00/00  '
  now='00:00:00   '

#if defined(IBM)
  !     call fdate(string)
  call date_and_time(dstring,tstring,tzone,dtvalue)
  today(1:2)=dstring(5:6)
  today(4:5)=dstring(7:8)
  today(7:8)=dstring(3:4)
  !      if(string(4:6).eq.'Jan') today(1:2)='01'
  !      if(string(4:6).eq.'Feb') today(1:2)='02'
  !      if(string(4:6).eq.'Mar') today(1:2)='03'
  !      if(string(4:6).eq.'Apr') today(1:2)='04'
  !      if(string(4:6).eq.'May') today(1:2)='05'
  !      if(string(4:6).eq.'Jun') today(1:2)='06'
  !      if(string(4:6).eq.'Jul') today(1:2)='07'
  !      if(string(4:6).eq.'Aug') today(1:2)='08'
  !      if(string(4:6).eq.'Sep') today(1:2)='09'
  !      if(string(4:6).eq.'Oct') today(1:2)='10'
  !      if(string(4:6).eq.'Nov') today(1:2)='11'
  !      if(string(4:6).eq.'Dec') today(1:2)='12'
  !      today(7:8)=string(8:9)
  !      today(4:5)=string(1:2)
  now(1:2)=tstring(1:2)
  now(4:5)=tstring(3:4)
  now(7:8)=tstring(5:6)
#endif
#if defined(KSR)
  string=ctime(time())
  if(string(5:7).eq.'Jan') today(1:2)='01'
  if(string(5:7).eq.'Feb') today(1:2)='02'
  if(string(5:7).eq.'Mar') today(1:2)='03'
  if(string(5:7).eq.'Apr') today(1:2)='04'
  if(string(5:7).eq.'May') today(1:2)='05'
  if(string(5:7).eq.'Jun') today(1:2)='06'
  if(string(5:7).eq.'Jul') today(1:2)='07'
  if(string(5:7).eq.'Aug') today(1:2)='08'
  if(string(5:7).eq.'Sep') today(1:2)='09'
  if(string(5:7).eq.'Oct') today(1:2)='10'
  if(string(5:7).eq.'Nov') today(1:2)='11'
  if(string(5:7).eq.'Dec') today(1:2)='12'
  today(7:8)=string(23:24)
  today(4:5)=string(9:10)
  now=string(11:20)
#endif
#if defined(CRAY_T3D) || defined(SP1) || defined(CRAY_T3E) || defined(SOLARIS)
  call date(string)
  if(string(5:7).eq.'Jan') today(1:2)='01'
  if(string(5:7).eq.'Feb') today(1:2)='02'
  if(string(5:7).eq.'Mar') today(1:2)='03'
  if(string(5:7).eq.'Apr') today(1:2)='04'
  if(string(5:7).eq.'May') today(1:2)='05'
  if(string(5:7).eq.'Jun') today(1:2)='06'
  if(string(5:7).eq.'Jul') today(1:2)='07'
  if(string(5:7).eq.'Aug') today(1:2)='08'
  if(string(5:7).eq.'Sep') today(1:2)='09'
  if(string(5:7).eq.'Oct') today(1:2)='10'
  if(string(5:7).eq.'Nov') today(1:2)='11'
  if(string(5:7).eq.'Dec') today(1:2)='12'
  today(7:8)=string(23:24)
  today(4:5)=string(9:10)
  now=string(11:20)
#endif
#if defined(LINUX) && !defined(IBM)
  call linux_date(today)
  call linux_time(now)
#endif
#if defined(SGI)
  call date(string)
  if(string(4:6).eq.'Jan') today(1:2)='01'
  if(string(4:6).eq.'Feb') today(1:2)='02'
  if(string(4:6).eq.'Mar') today(1:2)='03'
  if(string(4:6).eq.'Apr') today(1:2)='04'
  if(string(4:6).eq.'May') today(1:2)='05'
  if(string(4:6).eq.'Jun') today(1:2)='06'
  if(string(4:6).eq.'Jul') today(1:2)='07'
  if(string(4:6).eq.'Aug') today(1:2)='08'
  if(string(4:6).eq.'Sep') today(1:2)='09'
  if(string(4:6).eq.'Oct') today(1:2)='10'
  if(string(4:6).eq.'Nov') today(1:2)='11'
  if(string(4:6).eq.'Dec') today(1:2)='12'
  today(7:8)=string(8:9)
  today(4:5)=string(1:2)
  call time(now(1:8))
  now(9:10)='  '
#endif
  if(today(4:4).eq.' ') today(4:4)='0'
  return
end subroutine swatch


subroutine matinv(a,n,ndim)

  implicit none

  external :: tps_abort

  integer :: maxdim
  real (kind=8) :: zero,small,one
  parameter(maxdim=3)
  parameter(zero=0.0d0)
  parameter(small=1.0d-6)
  parameter(one=1.0d0)

  integer :: n,ndim
  real (kind=8) :: a(ndim,ndim)
  integer :: ia(2,maxdim),ib(maxdim),ic(maxdim)
  real (kind=8) :: d(maxdim)

  integer :: idim,i,j,k,l,m
  real (kind=8) :: b,e

  if(ndim.gt.maxdim) call tps_abort('matinv dimension error',0)

  do idim=1,n
    ia(1,idim)=0
    ia(2,idim)=0
  enddo

  do idim=1,n
    b=zero
    do l=1,n
      do m=1,n
        if(ia(1,l).ne.1.and.ia(2,m).ne.1) then
          e=dabs(a(l,m))
          if(e.ge.b) then
            i=l
            k=m
          endif
          b=dmax1(b,e)
        endif
      enddo
    enddo
    ia(1,i)=1
    ia(2,k)=1
    ib(k)=i
    ic(i)=k
    b=a(i,k)

    if(dabs(b).lt.small) call tps_abort('arg_matinv singular matrix',0)
    a(i,k)=one/b
    do l=1,n
      if(l.ne.k) a(i,l)=-a(i,l)/b
    enddo
    do l=1,n
      do m=1,n
        if(l.ne.i.and.m.ne.k) a(l,m)=a(l,m)+a(l,k)*a(i,m)
      enddo
    enddo
    do l=1,n
      if(l.ne.i) a(l,k)=a(l,k)/b
    enddo
  enddo

  do l=1,n
    do j=1,n
      k=ib(j)
      d(j)=a(k,l)
    enddo
    do j=1,n
      a(j,l)=d(j)
    enddo
  enddo

  do l=1,n
    do j=1,n
      k=ic(j)
      d(j)=a(l,k)
    enddo
    do j=1,n
      a(l,j)=d(j)
    enddo
  enddo

  return
end subroutine matinv



logical function frequency(istep,nstep)

  implicit none

  integer :: istep,nstep

  if(nstep.le.0) then
    frequency=.false.
  else
    frequency=mod(istep,nstep).eq.0
  endif

  return
end function frequency

subroutine rolex(elaps,cputim)

  implicit none

  real (kind=8) :: elaps,cputim

  real (kind=8) :: linux_walltime,linux_cputime
  external linux_walltime,linux_cputime

#if defined(IBM)
  call cpu_time(elaps)
  call cpu_time(cputim)
#else
  elaps=linux_walltime()
  cputim=linux_cputime()
#endif

  return
end subroutine rolex
    
      subroutine timer_init()

      implicit none

      external :: timer

      call timer(0,0)

      return
    end subroutine timer_init
    
      subroutine timer_reset(itime)

      implicit none

      external :: timer

      integer :: itime

      call timer(itime,-2)

      return
    end subroutine timer_reset
    
      subroutine timer_start(itime)

      implicit none

      external :: timer

      integer :: itime

      call timer(itime,0)

      return
    end subroutine timer_start
    
      subroutine timer_stop(itime)

      implicit none

      external :: timer

      integer :: itime

      call timer(itime,1)

      return
    end subroutine timer_stop
    
      subroutine timer(itime,iopt)

      implicit none

      external :: rolex,tps_abort

      integer :: itime,iopt

      real (kind=8) :: elaps,cputim

      integer :: mtime
      parameter(mtime=250)
      integer :: ncall(250)
      real (kind=8) :: ttime(250,3),ctime(250,3)
      common/tim/ncall,ttime,ctime

      integer :: i

      call rolex(elaps,cputim)

      if(itime.eq.0) then
      do i=1,mtime
      ncall(i)=0
      ctime(i,1)=0.0d0
      ttime(i,1)=0.0d0
      ctime(i,2)=0.0d0
      ttime(i,2)=0.0d0
      ctime(i,3)=1.0d9
      ttime(i,3)=1.0d9
      enddo
      elseif(itime.le.0.or.itime.gt.mtime) then
      call tps_abort('Timer index out of range',0)
      elseif(iopt.eq.-2) then
      ncall(itime)=0
      ctime(itime,1)=0.0d0
      ttime(itime,1)=0.0d0
      ctime(itime,2)=0.0d0
      ttime(itime,2)=0.0d0
      ctime(itime,3)=1.0d9
      ttime(itime,3)=1.0d9
      elseif(iopt.eq.-1) then
      ncall(itime)=0
      ctime(itime,1)=-cputim
      ttime(itime,1)=-elaps
      ctime(itime,2)=-cputim
      ttime(itime,2)=-elaps
      ctime(itime,3)=1.0d9
      ttime(itime,3)=1.0d9
      elseif(iopt.eq.0) then
      ctime(itime,1)=ctime(itime,1)-cputim
      ttime(itime,1)=ttime(itime,1)-elaps
      ctime(itime,2)=-cputim
      ttime(itime,2)=-elaps
      elseif(iopt.eq.1) then
      ncall(itime)=ncall(itime)+1
      ctime(itime,1)=ctime(itime,1)+cputim
      ttime(itime,1)=ttime(itime,1)+elaps
      ctime(itime,2)=ctime(itime,2)+cputim
      ttime(itime,2)=ttime(itime,2)+elaps
      ctime(itime,3)=min(ctime(itime,2),ctime(itime,3))
      ttime(itime,3)=min(ttime(itime,2),ttime(itime,3))
      else
      call tps_abort('Unimplemented timer option',0)
      endif

      return
    end subroutine timer
    
      real (kind=8) function timer_cpu(itime)

      implicit none

      external :: tps_abort

      integer :: itime

      integer :: mtime
      parameter(mtime=250)
      integer :: ncall(250)
      real (kind=8) :: ttime(250,3),ctime(250,3)
      common/tim/ncall,ttime,ctime

      if(itime.le.0.or.itime.gt.mtime) call tps_abort('Illegal timer index',0)

      if(ncall(itime).le.0) then
      timer_cpu=0.0d0
      else
      timer_cpu=ctime(itime,2)
      endif

      return
    end function timer_cpu
    
      real (kind=8) function timer_wall(itime)

      implicit none

      external :: tps_abort

      integer :: itime

      integer :: mtime
      parameter(mtime=250)
      integer :: ncall(250)
      real (kind=8) :: ttime(250,3),ctime(250,3)
      common/tim/ncall,ttime,ctime

      if(itime.le.0.or.itime.gt.mtime) call tps_abort('Illegal timer index',0)

      if(ncall(itime).le.0) then
      timer_wall=0.0d0
      else
      timer_wall=ttime(itime,2)
      endif

      return
      end
      integer function timer_calls(itime)

      implicit none

      external :: tps_abort
      
      integer :: itime
      integer :: mtime
      parameter(mtime=250)
      integer :: ncall(250)
      real (kind=8) :: ttime(250,3),ctime(250,3)
      common/tim/ncall,ttime,ctime

      if(itime.le.0.or.itime.gt.mtime) call tps_abort('Illegal timer index',0)

      timer_calls=ncall(itime)

      return
      end
      real (kind=8) function timer_cpu_minimum(itime)

      implicit none

      external :: tps_abort

      integer :: itime

      integer :: mtime
      parameter(mtime=250)
      integer :: ncall(250)
      real (kind=8) :: ttime(250,3),ctime(250,3)
      common/tim/ncall,ttime,ctime

      if(itime.le.0.or.itime.gt.mtime) call tps_abort('Illegal timer index',0)

      if(ncall(itime).le.0) then
      timer_cpu_minimum=0.0d0
      else
      timer_cpu_minimum=ctime(itime,3)
      endif

      return
      end
      real (kind=8) function timer_wall_minimum(itime)

      implicit none

      external :: tps_abort

      integer :: itime

      integer :: mtime
      parameter(mtime=250)
      integer :: ncall(250)
      real (kind=8) :: ttime(250,3),ctime(250,3)
      common/tim/ncall,ttime,ctime

      if(itime.le.0.or.itime.gt.mtime) call tps_abort('Illegal timer index',0)

      if(ncall(itime).le.0) then
      timer_wall_minimum=0.0d0
      else
      timer_wall_minimum=ttime(itime,3)
      endif

      return
      end
      real (kind=8) function timer_cpu_average(itime)

      implicit none

      external :: tps_abort

      integer :: itime

      integer :: mtime
      parameter(mtime=250)
      integer :: ncall(250)
      real (kind=8) :: ttime(250,3),ctime(250,3)
      common/tim/ncall,ttime,ctime

      if(itime.le.0.or.itime.gt.mtime) call tps_abort('Illegal timer index',0)

      if(ncall(itime).le.0) then
      timer_cpu_average=0.0d0
      else
      timer_cpu_average=ctime(itime,1)/dble(ncall(itime))
      endif

      return
      end
      real (kind=8) function timer_wall_average(itime)

      implicit none

      external :: tps_abort

      integer :: itime

      integer :: mtime
      parameter(mtime=250)
      integer :: ncall(250)
      real (kind=8) :: ttime(250,3),ctime(250,3)
      common/tim/ncall,ttime,ctime

      if(itime.le.0.or.itime.gt.mtime) call tps_abort('Illegal timer index',0)

      if(ncall(itime).le.0) then
      timer_wall_average=0.0d0
      else
      timer_wall_average=ttime(itime,1)/dble(ncall(itime))
      endif

      return
      end
      real (kind=8) function timer_cpu_total(itime)

      implicit none

      external :: tps_abort

      integer :: itime

      integer :: mtime
      parameter(mtime=250)
      integer :: ncall(250)
      real (kind=8) :: ttime(250,3),ctime(250,3)
      common/tim/ncall,ttime,ctime

      if(itime.le.0.or.itime.gt.mtime) call tps_abort('Illegal timer index',0)

      if(ncall(itime).le.0) then
      timer_cpu_total=0.0d0
      else
      timer_cpu_total=ctime(itime,1)
      endif

      return
      end
      real (kind=8) function timer_wall_total(itime)

      implicit none

      external :: tps_abort

      integer :: itime

      integer :: mtime
      parameter(mtime=250)
      integer :: ncall(250)
      real (kind=8) :: ttime(250,3),ctime(250,3)
      common/tim/ncall,ttime,ctime

      if(itime.le.0.or.itime.gt.mtime) call tps_abort('Illegal timer index',0)

      if(ncall(itime).le.0) then
      timer_wall_total=0.0d0
      else
      timer_wall_total=ttime(itime,1)
      endif

      return
      end
      subroutine tracer_init(itr,ltr,lfn)

      implicit none

      external :: rolex

      integer :: ntrace,itr,ltr,lfn,md0,md1,lfntrc
      real (kind=8) :: strace,ttrace(1000)
      integer :: itrace(1000)
      common/trc/strace,ttrace,ntrace,itrace,lfntrc,md0,md1

      real (kind=8) :: cputim

      ntrace=0
      call rolex(strace,cputim)
      md0=itr
      md1=ltr
      lfntrc=lfn

      return
      end
      subroutine tracer(md,id)

      implicit none

      external :: rolex

      integer :: ntrace
      real (kind=8) :: strace,ttrace(1000)
      integer :: md,id,lfntrc,md0,md1,itrace(1000)
      common/trc/strace,ttrace,ntrace,itrace,lfntrc,md0,md1

      integer :: i
      real (kind=8) :: elaps,cputim


      if(md.lt.md0.or.md.gt.md1) return

!      if(id.eq.0)
!     + write(lfntrc,100) (itrace(i),i=1,ntrace)
!  100 format(6i12)

      if(ntrace.lt.950) then
      ntrace=ntrace+1
      elseif(id.eq.0) then
      write(lfntrc,1000) (ttrace(i),i=1,ntrace)
 1000 format(100f12.6)
      ntrace=1
      else
      ntrace=ntrace+1
      endif

      itrace(ntrace)=100*md+id
      call rolex(elaps,cputim)
      ttrace(ntrace)=elaps-strace
      strace=elaps

      return
      end
      subroutine tracer_end(itr)

      implicit none

      integer :: itr

      integer :: ntrace
      real (kind=8) :: strace,ttrace(1000)
      integer :: lfntrc,md0,md1,itrace(1000)
      common/trc/strace,ttrace,ntrace,itrace,lfntrc,md0,md1

      integer :: i

      if(itr.gt.0) write(lfntrc,1000) (ttrace(i),i=1,ntrace)
 1000 format(100f12.6)

      return
      end
      subroutine auto_corr(data,ndata,aver,acf,lacf,ratio)

      implicit none

      real (kind=8) :: zero,half,one,two
      parameter(zero=0.0d0)
      parameter(half=5.0d-1)
      parameter(one=1.0d0)
      parameter(two=2.0d0)

      integer :: ndata,lacf,nacf
      real (kind=8) :: data(ndata),acf(lacf),aver,ratio

      integer :: i,j
      real (kind=8) :: dsum,dvar

      nacf=min(ndata,lacf)

      dsum=zero
      do 1 i=1,ndata
      dsum=dsum+(data(i)-aver)**2
    1 continue
      dvar=dble(ndata)/dsum
      do 2 i=1,nacf-1
      dsum=zero
      do 3 j=1,ndata-i
      dsum=dsum+(data(j)-aver)*(data(i+j)-aver)
    3 continue
      acf(i)=dvar*(dsum/dble(ndata-i))*(one-dble(i-1)/dble(ndata-2))
    2 continue

      ratio=half
      do 4 i=1,nacf-1
      ratio=ratio+(one-(dble(i)/dble(ndata))**2)*abs(acf(i))
    4 continue
      ratio=sqrt(two*abs(ratio))
      if(ratio.lt.one) ratio=one

      lacf=nacf

      return
      end
      subroutine acf_approx(acf,acfapp,lacf)

      implicit none

      integer :: lacf
      real (kind=8) :: acf(lacf),acfapp(lacf)
      integer :: i

      do i=1,lacf
        acfapp(i)=acf(i)
      enddo
      
!      integer :: kapprox(15)

!      rfact=two*sqrt(dble(klarge))/dble(nacfa-1)
!      wsum1=zero
!      wsum2=zero
!      warg=zero
!      do 5 i=1,nacfa-1
!      acf(i)=abs(acf(i))
!      if(acf(i).gt.zero) then
!      wterm=exp(weight*dble(nacfa-i)/dble(nacfa-1))
!      wsum1=wsum1+wterm
!      wsum2=wsum2+wterm*log(acf(i))/(rfact*dble(i))
!      endif
!    5 continue
!      if(abs(wsum1).gt.small) warg=wsum2/wsum1
!      do 55 i=1,nacfa-1
!      acf(i)=acf(i)-exp(warg*dble(i)*rfact)
!   55 continue
!      nfunc=iapprx
!      if(nfunc.le.0) nfunc=15
!      do 6 k=1,nfunc
!      do 7 l=1,nfunc
!      cdap(k,l)=zero
!      do 8 m=1,nacfa-1
!      xappm=dble(m)*rfact
!      cdawgt=exp(weight*dble(nacfa-m)/dble(nacfa-1))
!      cdap(k,l)=cdap(k,l)+cdawgt*exp((-xappm)*xappm)*
!     + approx(k)*approx(l)*(xappm**(kapprx(k)+kapprx(l)))
!    8 continue
!      cdaq(k,l)=cdap(k,l)
!    7 continue
!    6 continue
!      do 9 i=1,nfunc
!      cdad(i)=zero
!      do 10 j=1,nacfa-1
!      xappj=dble(j)*rfact
!      cdawgt=exp(weight*dble(nacfa-j)/dble(nacfa-1))
!      cdad(i)=cdad(i)+cdawgt*exp((-half)*xappj*xappj)*
!     + acf(j)*approx(i)*(xappj**kapprx(i))
!   10 continue
!    9 continue
!      do 11 k=1,nfunc
!      do 12 i=k,nfunc
!      cdaq(i,k)=cdap(i,k)
!      do 13 l=1,k-1
!      cdaq(i,k)=cdaq(i,k)-cdaq(i,l)*cdaq(l,k)
!   13 continue
!   12 continue
!      do 14 i=k+1,nfunc
!      cdaq(k,i)=cdap(k,i)
!      do 15 l=1,k-1
!      cdaq(k,i)=cdaq(k,i)-cdaq(k,l)*cdaq(l,i)
!   15 continue
!      cdaq(k,i)=cdaq(k,i)/cdaq(k,k)
!   14 continue
!   11 continue
!      do 16 j=1,nfunc
!      cdac(j)=cdad(j)
!      do 17 i=1,j-1
!      cdac(j)=cdac(j)-cdaq(j,i)*cdac(i)
!   17 continue
!      cdac(j)=cdac(j)/cdaq(j,j)
!   16 continue
!      do 18 k=1,nfunc
!      j=nfunc+1-k
!      do 19 i=j+1,nfunc
!      cdac(j)=cdac(j)-cdaq(j,i)*cdac(i)
!   19 continue
!   18 continue
!      do 20 i=1,nacfa-1
!      xappi=dble(i)*rfact
!      acf(i)=exp(warg*xappi)
!      do 21 j=1,nfunc
!      acf(i)=acf(i)+exp((-half)*xappi*xappi)*
!     + approx(j)*cdac(j)*(xappi**kapprx(j))
!   21 continue
!   20 continue
!
!     autocorrelation function upto lag min(nacf,ndata)
!
!      if(iopt.gt.0) then
!      warg=zero
!      nacfa=nacf
!      if(nacfa.gt.ndata) nacfa=ndata
!      dfsum=zero
!      do 2 i=1,ndata
!      dfsum=dfsum+(data(i)-aver)**2
!    2 continue
!      dvar=dble(ndata)/dfsum
!      do 3 i=1,nacfa-1
!      dfsum=zero
!      do 4 j=1,ndata-i
!      dfsum=dfsum+(data(j)-aver)*(data(i+j)-aver)
!    4 continue
!      acf(i)=dvar*(dfsum/dble(ndata-i))*(one-dble(i-1)/dble(ndata-2))
!    3 continue
!      endif
!c
!c     approximate acf here
!c
!      if(iopt.gt.1) then
!      rfact=two*sqrt(dble(klarge))/dble(nacfa-1)
!      wsum1=zero
!      wsum2=zero
!      warg=zero
!      do 5 i=1,nacfa-1
!      acf(i)=abs(acf(i))
!      if(acf(i).gt.zero) then
!      wterm=exp(weight*dble(nacfa-i)/dble(nacfa-1))
!      wsum1=wsum1+wterm
!      wsum2=wsum2+wterm*log(acf(i))/(rfact*dble(i))
!      endif
!    5 continue
!      if(abs(wsum1).gt.small) warg=wsum2/wsum1
!      do 55 i=1,nacfa-1
!      acf(i)=acf(i)-exp(warg*dble(i)*rfact)
!   55 continue
!      nfunc=iapprx
!      if(nfunc.le.0) nfunc=15
!      do 6 k=1,nfunc
!      do 7 l=1,nfunc
!      cdap(k,l)=zero
!      do 8 m=1,nacfa-1
!      xappm=dble(m)*rfact
!      cdawgt=exp(weight*dble(nacfa-m)/dble(nacfa-1))
!      cdap(k,l)=cdap(k,l)+cdawgt*exp((-xappm)*xappm)*
!     + approx(k)*approx(l)*(xappm**(kapprx(k)+kapprx(l)))
!    8 continue
!      cdaq(k,l)=cdap(k,l)
!    7 continue
!    6 continue
!      do 9 i=1,nfunc
!      cdad(i)=zero
!      do 10 j=1,nacfa-1
!      xappj=dble(j)*rfact
!      cdawgt=exp(weight*dble(nacfa-j)/dble(nacfa-1))
!      cdad(i)=cdad(i)+cdawgt*exp((-half)*xappj*xappj)*
!     + acf(j)*approx(i)*(xappj**kapprx(i))
!   10 continue
!    9 continue
!      do 11 k=1,nfunc
!      do 12 i=k,nfunc
!      cdaq(i,k)=cdap(i,k)
!      do 13 l=1,k-1
!      cdaq(i,k)=cdaq(i,k)-cdaq(i,l)*cdaq(l,k)
!   13 continue
!   12 continue
!      do 14 i=k+1,nfunc
!      cdaq(k,i)=cdap(k,i)
!      do 15 l=1,k-1
!      cdaq(k,i)=cdaq(k,i)-cdaq(k,l)*cdaq(l,i)
!   15 continue
!      cdaq(k,i)=cdaq(k,i)/cdaq(k,k)
!   14 continue
!   11 continue
!      do 16 j=1,nfunc
!      cdac(j)=cdad(j)
!      do 17 i=1,j-1
!      cdac(j)=cdac(j)-cdaq(j,i)*cdac(i)
!   17 continue
!      cdac(j)=cdac(j)/cdaq(j,j)
!   16 continue
!      do 18 k=1,nfunc
!      j=nfunc+1-k
!      do 19 i=j+1,nfunc
!      cdac(j)=cdac(j)-cdaq(j,i)*cdac(i)
!   19 continue
!   18 continue
!      do 20 i=1,nacfa-1
!      xappi=dble(i)*rfact
!      acf(i)=exp(warg*xappi)
!      do 21 j=1,nfunc
!      acf(i)=acf(i)+exp((-half)*xappi*xappi)*
!     + approx(j)*cdac(j)*(xappi**kapprx(j))
!   21 continue
!   20 continue
!      endif
!c
!c     sampling ratio
!c
!      ratio=one
!      if(iopt.gt.0) then
!      ratio=half
!      do 22 i=1,nacfa-1
!      ratio=ratio+(one-(dble(i)/dble(ndata))**2)*abs(acf(i))
!   22 continue
!      ratio=sqrt(two*abs(ratio))
!      if(ratio.lt.one) ratio=one
!      endif
!c
!c
!      corerr=ratio*stderr
!c
      return
      end
#if defined(CRAY_T3D) || defined(CRAY_T3E)
      integer function util_nint(x)
      real (kind=8) :: x
      util_nint=nint(x)
      return
      end
#endif
      logical function argos_zmat(x1,x2,x3,x4,d,a,t)

      implicit none

      external :: tps_rotate

      real (kind=8) :: x1(3),x2(3),x3(3),x4(3),d,a,t

      real (kind=8) :: p(3),q(3),s(3),v(3),w(3)
      real (kind=8) :: r
      integer :: i

      do 1 i=1,3
      p(i)=x3(i)-x2(i)
      q(i)=x4(i)-x2(i)
    1 continue
      v(1)=p(2)*q(3)-p(3)*q(2)+x2(1)
      v(2)=p(3)*q(1)-p(1)*q(3)+x2(2)
      v(3)=p(1)*q(2)-p(2)*q(1)+x2(3)

      r=sqrt(p(1)*p(1)+p(2)*p(2)+p(3)*p(3))
      do 2 i=1,3
      s(i)=d*p(i)/r+x2(i)
    2 continue

      call tps_rotate(x2,v,a,s,w)
      call tps_rotate(x3,x2,t,w,s)

      do 3 i=1,3
      x1(i)=s(i)
    3 continue

      argos_zmat=.true.
      return
      end

      subroutine tps_abort(message,id)
      implicit none
      character(*) :: message
      integer :: id
      integer (kind=4) :: ierror,ierr

      write(*,999) id,message
 999  format('Error code ',i2,': ',a)

      ierr=0
      ierror=0
!      call MPI_Abort(MPI_COMM_WORLD,ierror,ierr)
      stop "Aborting from tps_abort"
      return
      end
