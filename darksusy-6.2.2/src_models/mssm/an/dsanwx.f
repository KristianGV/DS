*******************************************************************************
*** Function dsanwx provides the  WIMP self-annihilation invariant rate.    ***
***                                                                         ***
***  type : interface                                                       ***
***                                                                         ***
***  Input:                                                                 ***
***    p - initial cm momentum (real) for DM annihilations                  ***
***  Output:                                                                ***
***  BeginTex
***    \begin{displaymath}
***    W_{\rm{eff}} = \sum_{ij}\frac{p_{ij}}{p_{11}}
***    \frac{g_ig_j}{g_1^2} W_{ij} = 
***    \sum_{ij} \sqrt{\frac{[s-(m_{i}-m_{j})^2][s-(m_{i}+m_{j})^2]}
***    {s(s-4m_1^2)}} \frac{g_ig_j}{g_1^2} W_{ij}.
***    \end{displaymath}
***    where the $p$'s are the momenta, the $g$'s are the internal
***    degrees of freedom, the $m$'s are the masses and $W_{ij}$ is
***    the invariant annihilation rate for the included subprocess.
***  EndTex
***  uses dsabsq.
***  passed to dsrdens by dsrdomega.                                      *** 
***                                                                       ***
***  author: paolo gondolo (gondolo@lpthe.jussieu.fr) 1994
***  modified: joakim edsjo (edsjo@fysik.su.se) 97-09-09
***  mod TB 2016-02-08: added OMP compatibility
***  mod TB 2018-04-27: caught NAN output
*******************************************************************************
      real*8 function dsanwx(p)

      implicit none
      include 'dsmssm.h'
      include 'dsandwcom.h'
      include 'dsidtag.h'
      include 'dsio.h'
      real*8 p,dsanwxint,dsandwdcosopt,dsandwdcos
      real*8 epsabs,epsrel,result,abserr,alist(5000),blist(5000),
     &  elist(5000),rlist(5000)
      real*8 al,be,h1width
      integer ier,iord(5000),last,limit,neval

      logical dsisnan

      real*8 dsandwdcoss,eps,sum,y1,y2,dsandwdcosy
      external dsandwdcoss,dsandwdcosopt,dsandwdcos,dsandwdcosy

      real*8 dsandwdcosd
      external dsandwdcosd

      integer dsidnumber
      integer idold
      data idold/-123456789/
      save idold

c-----------------------------------------------------------------------
      real*8 pd
      common /gadint/ pd

      real*8 alph,bet
      common /yint/ alph,bet
      save /yint/
*$OMP THREADPRIVATE (/gadint/,/yint/)
c-----------------------------------------------------------------------


c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('MSSM','dsanwx')


c...temporarily change width of h1 higgs boson
      h1width=width(kh1)
      width(kh1)=max(width(kh1),0.1d0) ! je, dec 11, 1998

      if (idold.ne.dsidnumber()) then
        call dsanalbe(al,be)
        alph=al
        bet=be
        idold=dsidnumber()
      endif

      goto 90

      epsabs=1.d-2     !numerical accuracy
      epsrel=1.d-2
      limit=5000
      pd = p
      call dqagse(dsandwdcosd,-1.0d0,1.0d0,epsabs,epsrel,limit,result,
     &   abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      dsanwx=result/1.0d15
      goto 900


 90   eps=0.001d0
      sum=0.0d0
      pd = p

c      pd=16.064d0
c      call gewspecs(dsandwdcoss,-1.0d0,1.0d0,100,
c     & 'p1cth.dat                      ',33)
c      y1=(mco(1)**2+bet*pd**2)**(alph)
c      y2=(mco(1)**2)**(alph)
c      call gewspecs(dsandwdcosy,y1,y2,100,
c     & 'p1y.dat                        ',33)

c      pd=325.0d0
c      call gewspecs(dsandwdcoss,-1.0d0,1.0d0,100,
c     & 'p2cth.dat                      ',33)
c      y1=(mco(1)**2+bet*pd**2)**(alph)
c      y2=(mco(1)**2)**(alph)
c      call gewspecs(dsandwdcosy,y1,y2,100,
c     & 'p2y.dat                         ',33)

c      pd=1000.0d0
c      call gewspecs(dsandwdcoss,-1.0d0,1.0d0,100,
c     & 'p3cth.dat                      ',33)
c      y1=(mco(1)**2+bet*pd**2)**(alph)
c      y2=(mco(1)**2)**(alph)
c      call gewspecs(dsandwdcosy,y1,y2,100,
c     & 'p3y.dat                        ',33)
c      stop

      if (p.ge.2.5d0*mco(1)) then
c        call dsanalbe2(p,alph,be)
c        bet=be
        y1=(mco(1)**2+bet*pd**2)**(alph)
        y2=(mco(1)**2)**(alph)
c...changed to $OMP by TB
        call dgadapOMP(y1,y2,dsandwdcosy,eps,sum)
        dsanwx=sum/1.0d15
      else
c...changed to $OMP by TB
        call dgadapOMP(-1.0d0,1.0d0,dsandwdcoss,eps,sum)
        dsanwx=sum/1.0d15
      endif
c      write(*,*) rdtag,nr,p,dsanwx

c      call gadap(-0.99,0.99,dsandwdcoss,eps,sum)
c      dsanwx=dsanwx+dble(sum)/1.0d15
c      call gadap(0.99,1.0,dsandwdcoss,eps,sum)
c      dsanwx=dsanwx+dble(sum)/1.0d15
      goto 900

      if (p.lt.0.5d0*mco(1)) then
        dsanwx=dsanwxint(p,-1.0d0,1.0d0)
      else if ((p.ge.0.5d0*mco(1)).and.(p.lt.1.0d0*mco(1))) then
        dsanwx=dsanwxint(p,-1.0d0,-0.9d0)
        dsanwx=dsanwx+dsanwxint(p,-0.9d0,0.9d0)
        dsanwx=dsanwx+dsanwxint(p,0.9d0,1.0d0)
      else
        dsanwx=dsanwxint(p,-1.0d0,-0.999d0)
        dsanwx=dsanwx+dsanwxint(p,-0.999d0,-0.99d0)
        dsanwx=dsanwx+dsanwxint(p,-0.99d0,-0.97d0)
        dsanwx=dsanwx+dsanwxint(p,-0.97d0,-0.9d0)
        dsanwx=dsanwx+dsanwxint(p,-0.9d0,-0.5d0)
        dsanwx=dsanwx+dsanwxint(p,-0.5d0,0.5d0)
        dsanwx=dsanwx+dsanwxint(p,0.5d0,0.9d0)
        dsanwx=dsanwx+dsanwxint(p,0.9d0,0.97d0)
        dsanwx=dsanwx+dsanwxint(p,0.97d0,0.99d0)
        dsanwx=dsanwx+dsanwxint(p,0.99d0,0.999d0)
        dsanwx=dsanwx+dsanwxint(p,0.999d0,1.0d0)
      endif

 900  continue

      width(kh1)=h1width        ! change back h1 width

      if (dsisnan(dsanwx)) then
        if (prtlevel.ge.1) write(*,*) 'WARNING: dsanwx returned NaN !!!'
        if (prtlevel.ge.2) write(*,*) 'p = ',p
        dsanwx=1.0d-50
      endif

      return

      end
