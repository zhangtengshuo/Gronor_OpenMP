


      subroutine gfunct (l,m,a,b,p,t,g,n)

      implicit none

      integer :: l,m,n
      real (kind=8) :: a,b,p,t,g(7,3)

      integer :: ll,mm
      real (kind=8) :: temp,aa,bb,pp,ab,p2,t2,a2,f0,f1,f2,f3
      real (kind=8) :: a3,a1,f4,f5,b2

      ll=l+1
      mm=m+1
      goto (100,101,102,103),ll
  100 goto (110,111,112,113),mm
  101 goto (120,121,122,123),mm
  102 goto (130,131,132,133),mm
  103 goto (140,141,142,143),mm
!      00
  110 g(1,n)=1.
      goto 300
!      01
  111 g(1,n)=b
      g(2,n)=-p
      goto 300
!      02
  112 g(1,n)=b*b+0.5*t
      g(2,n)=-2.*b*p-0.5*t
      g(3,n)=p*p
      goto 300
!      03
  113 temp=a
      a=b
      b=temp
      goto 140
!      10
  120 g(1,n)=a
      g(2,n)=-p
      goto 300
!      11
  121 g(1,n)=a*b+0.5*t
      g(2,n)=-p*(a+b)-0.5*t
      g(3,n)=p*p
      goto 300
!      12
  122 g(1,n)=b*b*a+0.5*t*(a+2.*b)
      g(2,n)=-p*b*(2.*a+b)-0.5*t*((a+2.*b)+3.*p)
      g(3,n)=p*(p*(a+2.*b)+1.5*t)
      g(4,n)=-p*p*p
      goto 300
!      13
  123 temp=a
      a=b
      b=temp
      goto 141
!      20
  130 g(1,n)=a*a+0.5*t
      g(2,n)=-2.*a*p-0.5*t
      g(3,n)=p*p
      goto 300
!      21
  131 g(1,n)=a*a*b+0.5*t*(2.*a+b)
      g(2,n)=-p*a*(a+2.*b)-0.5*t*((2.*a+b)+3.*p)
      g(3,n)=p*(p*(2.*a+b)+1.5*t)
      g(4,n)=-p*p*p
      goto 300
!      22
  132 aa=a*a
      bb=b*b
      pp=p*p
      ab=4.*a*b
      g(1,n)=aa*bb+t*(0.5*(aa+ab+bb)+0.75*t)
      g(2,n)=-(2.*p*(aa*b+a*bb)+t*(0.5*(aa+ab+bb)+3.*(a+b)*p  +1.5*t))
      g(3,n)=pp*((aa+ab+bb)+3.*t)+t*(3.*p*(a+b)+0.75*t)
      g(4,n)=-(pp*(2.*p*(a+b)+3.*t))
      g(5,n)=pp*pp
      goto 300
!      23
  133 temp=a
      a=b
      b=temp
      goto 142
!      30
  140 g(1,n)=a*(a*a+1.5*t)
      g(2,n)=-3.*(a*(a*p+0.5*t)+0.5*p*t)
      g(3,n)=3.*p*(a*p+0.5*t)
      g(4,n)=-p*p*p
      goto 300
!      31
  141 p2=p*p
      t2=t*t
      a2=a*a
      ab=a*b
      f0=a2*ab
      f1=a2*(a+3.*b)
      f2=3.*a*(a+b)
      f3=3.*a+b
      g(1,n)=f0+.5*t*f2+.75*t2
      g(2,n)=-(p*f1+.5*t*f2+1.5*p*t*f3+1.5*t2)
      g(3,n)=p2*f2+1.5*p*t*f3+3.*t*p2+.75*t2
      g(4,n)=-(p*p2*f3+3.*t*p2)
      g(5,n)=p2*p2
      goto 300
!      32
  142 p2=p*p
      a2=a*a
      b2=b*b
      ab=a*b
      t2=t*t
      a3=3.*a2+6.*ab+b2
      a1=a2+6.*ab+3.*b2
      g(1,n)=a*a2*b2+0.5*t*a1*a+0.75*t2*(3.*a+2.*b)
      g(2,n)=-(ab*p*(2.*a2+3.*ab)+0.5*t*a*a1+1.5*t*a3*p+1.5*t2*(3.*a+2.*
     1b)+3.75*p*t2)
      g(3,n)=p2*a*a1+1.5*p*t*a3+3.*(3.*a+2.*b)*(t*p2+.25*t2)+7.5*p*t2
      g(4,n)=-(p2*(p*a3+t*( 9.*a+ 6.*b)+5.*p*t)+3.75*p*t2)
      g(5,n)=p*p2*(p*(3.*a+2.*b)+5.*t)
      g(6,n)=-p*p2*p2
      goto 300
!      33
  143 p2=p*p
      a2=a*a
      b2=b*b
      t2=t*t
      ab=a*b
      f0=a2*b2*ab
      f1=3.*a2*b2*(a+b)
      f2=3.*ab*(a2+3.*ab+b2)
      f3=a2*(a+9.*b)+b2*(b+9.*a)
      f4=3.*(a2+3.*ab+b2)
      f5=3.*(a+b)
      g(1,n)=f0+.5*t*f2+.75*t2*f4+1.875*t*t2
      g(2,n)=-(p*f1+.5*t*f2+1.5*t*p*f3+1.5*t2*f4+3.75*p*t2*f5+5.625*t*
     1t2)
      g(3,n)=p2*f2+1.5*p*t*f3+3.*f4*(t*p2+.25*t2)+7.5*p*t2*f5+11.25*t2*
     1(.5*t+p2)
      g(4,n)=-(p2*(p*f3+3.*t*f4)+5.*f5*p*t*(p2+.75*t)+t2*(22.50*p2+1.875
     1*t))
      g(5,n)=p2*(f4*p2+5.*p*t*f5+7.5*p2*t+11.25*t2)
      g(6,n)=-(p2*p2*(p*f5+7.5*t))
      g(7,n)=p2*p2*p2
  300 continue
      return
      end

      real (kind=8) function ovlap (l,m,a,b,t)
!
      implicit none
      integer :: l,m
      real (kind=8) :: a,b,t
!
      integer :: ll,mm
      real (kind=8) :: aa,tt,bb,ab,abab,t2,b2,a2,u,u2,a4

      ll=l+1
      mm=m+1
      goto (100,101,102,103,104,105),ll
  100 goto (110,111,112,113,200,200),mm
  101 goto (120,121,122,123,200,200),mm
  102 goto (130,131,132,133,200,200),mm
  103 goto (140,141,142,143,200,200),mm
  104 goto (150,151,152,153,200,200),mm
  105 goto (160,161,162,163,200,200),mm
  200 write(*,201) l,m,a,b,t
  201 format (//2i5,3f17.7,15h error in ovlap  )
      stop
!      00
  110 ovlap=1.
      goto 300
!      01
  111 ovlap=b
      goto 300
!      02
  112 ovlap=b*b+0.5*t
      goto 300
!      03
  113 ovlap=b*(b*b+1.5*t)
      goto 300
!      10
  120 ovlap=a
      goto 300
!      11
  121 ovlap=a*b+0.5*t
      goto 300
!      12
  122 ovlap=a*b*b+t*(b+0.5*a)
      goto 300
!      13
  123 ovlap=b*(b*(b*a+1.5*t)+1.5*a*t)+0.75*t*t
      goto 300
!      20
  130 ovlap=a*a+0.5*t
      goto 300
!      21
  131 ovlap=a*a*b+t*(a+0.5*b)
      goto 300
!      22
  132 ovlap=a*a*b*b+t*(0.5*(a*a+b*b)+2.*a*b+0.75*t)
      goto 300
!      23
  133 aa=a*a
      tt=t*t
      ovlap=b*(b*(b*(aa+0.5*t)+3.*a*t)+1.5*aa*t+2.25*tt)+1.5*a*tt
      goto 300
!      30
  140 ovlap=a*(a*a+1.5*t)
      goto 300
!      31
  141 ovlap=a*(a*(a*b+1.5*t)+1.5*b*t)+0.75*t*t
      goto 300
!      32
  142 bb=b*b
      tt=t*t
      ovlap=a*(a*(a*(bb+0.5*t)+3.*b*t)+1.5*bb*t+2.25*tt)+1.5*b*tt
      goto 300
!      33
  143 ab=a*b
      abab=a*a+3.*ab+b*b
      ovlap=ab*(ab*ab+1.5*t*abab)+t*t*(2.25*abab+1.875*t)
      goto 300
!      40
  150 ovlap=a*a*(a*a+3.*t)+0.75*t*t
      goto 300
!      41
  151 ovlap=a*(a*(a*(a*b+2.*t)+3.*b*t)+3.*t*t)+0.75*b*t*t
      goto 300
!      42
  152 bb=b*b
      tt=t*t
      ovlap=a*(a*(a*(a*(bb+0.5*t)+4.*b*t)+3.*bb*t+4.5*tt)+6.*b*tt)+
     1tt*(0.75*bb+1.875*t)
      goto 300
!      43
  153 t2=t*t
      b2=b*b
      ovlap=a*(a*(a*(a*(b*(b2+1.5*t))+3.*t*(2.*b2+t))+3.*b*t*(b2+4.5*t))
     1+t2*(9.*b2+7.5*t))+b*t2*(0.75*b2+5.625*t)
      goto 300
!      50
  160 a2=a*a
      ovlap=a*(5.*t*(0.75*t+a2)+a2*a2)
      goto 300
!      51
  161 t2=t*t
      ovlap=a*(a*(a*(a*(a*b+2.5*t)+5.*b*t)+7.5*t2)+3.75*b*t2)+1.875*t*t2
      goto 300
!      52
  162 t2=t*t
      b2=b*b
      ovlap=a*(a*(a*(a*(a*(0.5*t+b2)+5.*b*t)+5.*b2*t+7.5*t2)+
     115.*b*t2)+t2*(3.75*b2+9.375*t))+3.75*b*t*t2
      goto 300
!      53
  163 u=t/2.
      u2=u*u
      a2=a*a
      a4=a2*a2
      b2=b*b
      ab=a*b
      ovlap=ab*a4*b2+u*(3.*ab*a4+15.*a4*b2+10.*a2*b2*ab)+3.*u2*(
     15.*a4+30.*a2*ab+30.*a2*b2+5.*ab*b2)+15.*u*u2*(10.*a2+15.*ab+
     23.*b2)+105.*u2*u2
  300 continue
      return
      end

      real (kind=8) function fmch( m,x,y)

      use const

      implicit none

      real (kind=8) :: x,y
      integer :: m

!        this subroutine evaluates the integral from 0 to 1 of
!                       (u**(2*m))*expf(-x*(u**2))

      integer :: i,notrms
      real (kind=8) a,term,ptlsum,b,xd,approx,fimult,fiprop

      if(x-10.) 10,10,20
   10  a=m
       a=a+0.5
       term=1.0/a
       ptlsum=term
       do 11 i=2,50
       a=a+1.0
       term=term*x/a
       ptlsum=ptlsum+term
       if(term/ptlsum-0.00000001) 12, 11, 11
   11  continue
   12 fmch=0.5*ptlsum*y
      goto 150
   20  a=m
      b=a+0.5
       a=a-0.5
      xd=1./x
      approx=piroot*(sqrt(xd)*xd**m)
       if(m) 21, 23, 21
   21  do 22 i=1,m
       b=b-1.0
   22  approx=approx*b
   23 fimult=0.5*y*xd
      ptlsum=0.
      if(fimult) 421,25,421
  421 continue
       fiprop=fimult/approx
       term=1.0
       ptlsum=term
       notrms=x
       notrms=notrms+m
       do 24 i=2,notrms
       term=term*a*xd
       ptlsum=ptlsum+term
      if(abs (term*fiprop/ptlsum)-0.00000001)  25, 25, 24
   24  a=a-1.0
   25  fmch=approx-fimult*ptlsum
  150 return
  999  format (24h no convergence for fmch, i6, e16.9)
      end
