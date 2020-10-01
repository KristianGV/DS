      subroutine dssuconst_yukawa_running
c_______________________________________________________________________
c  Calculate Yukawas with running masses
c  common:
c    'dsmssm.h' - file with susy common blocks
c  author: paolo gondolo 1994-1999
c  modified: 031105 neutrino's yukawa corrected (pg)
c  modified: 2015-10-29 switched of yukawa running for quarks if full NLO 
c                       result for cross section is used (torsten bringmann) 
c
c This routine replaces mass(kq) in the tree-level expressions with
c dsrmq(2m*mass(lsp),kq). This is always the correct scaling with 
c energy, but the overall *normalization* / renormalization condition
c will in general depend on the process.
c For neutralino-neutralino annihilation, e.g., the result must be
c multiplied by mass(kq)/dsrmq(2m*mass(kq),kq), see 1510.02473. 

c=======================================================================
      implicit none
      include 'dsmssm.h'
      include 'dsmpconst.h'
      real*8 aux,mscale
      real*8 dsrmq,dsralph3

      call dssuconst_yukawa ! non-running default yukawas

      mscale=2.d0*mass(lsp)
      alph3=dsralph3(mscale)
      g3stro=sqrt(4.d0*pi*alph3)
c      write(*,*) 'lsp = ',lsp
c      write(*,*) 'mass(lsp) = ',mass(lsp)
c      write(*,*) 'mscale=',mscale
c      write(*,*) 'alpha3=',alph3,' alph3mz=',alph3mz

      aux = g2weak/dsqrt(2.d0)/mass(kw)
      yukawa(ktau)= aux*dsrmq(mscale,ktau)/cosbe
      if (NLOoption.eq.'off') then
        yukawa(kqu(2))= aux*dsrmq(mscale,kc)/sinbe
        yukawa(kqu(3))= aux*dsrmq(mscale,kt)/sinbe
        yukawa(kqd(3))= aux*dsrmq(mscale,kb)/cosbe
      else ! default
        yukawa(kqu(2))= aux*mass(kc)/sinbe
        yukawa(kqu(3))= aux*mass(kt)/sinbe
        yukawa(kqd(3))= aux*mass(kb)/cosbe      
      endif

      end
