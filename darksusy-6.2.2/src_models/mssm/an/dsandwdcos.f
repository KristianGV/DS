      function dsandwdcos(p,costheta)
c_______________________________________________________________________
c  annihilation differential invariant rate.
c  input:
c    p - initial cm momentum (real) for lsp annihilations
c    costheta - cosine of c.m. annihilation angle
c  Output:
c  BeginTex
c   \begin{displaymath}
c     \frac{dW_{ij}}{d\cos\theta}
c   \end{displaymath}
c   where
c   \begin{displaymath}
c      W_{ij} = 4 p_{ij} \sqrt{s} \sigma_{ij} 
c      = 4 \sigma_{ij} \sqrt{(p_i \cdot p_j)^2 - m_i^2 m_j^2} 
c      = 4 E_{i} E_{j} \sigma_{ij} v_{ij}
c   \end{displaymath}
c  EndTex
c  The returned dW/dcos(theta) is dimensionless
c  uses dsandwdcosnn, dsandwdcoscn and dsandwdcoscc and
c  routines in src/as
c  called by dsanwx.
c  author: joakim edsjo (edsjo@fysik.su.se)
c  date: 96-02-21
c  modified: 97-05-12 Joakim Edsjo (edsjo@fysik.su.se)
c  modified: 01-01-30 paolo gondolo (paolo@mamma-mia.phys.cwru.edu)
c  modified: 02-03-09 Joakim Edsjo (edsjo@fysik.su.se)
c  modified: 06-02-22 Paolo Gondolo (paolo@physics.utah.edu)
c=======================================================================
      implicit none
      include 'dsandwcom.h'
      real*8 dsandwdcos
      real*8 mx
      real*8 w,p,costheta,pnew,s,mp1,mp2,f,tmp
      integer i,j

      real*8 dsandwdcosij,dw
      integer kkp1,kkp2,kkx

c      write(*,*) 'dsandwdcos called with p = ',p
c      write(*,*) '            costheta = ',costheta
c...initialize dsankinvar variables
      call dsankinvar(0.0d0,0.0d0,0,0,0,0,0)


c      do i=1,nco
c         write (*,*) "PG-DEBUG: i,kco,mco,mass,mdof,kdof",
c     &        i,kco(i),mco(i),mass(kco(i)),mdof(i),kdof(kco(i))
c      enddo

      kkx = kco(1)
      mx = mco(1)
      s = 4.0d0*(mx**2+p**2)
      w = 0.d0
      do i=1,nco
         do j=i,nco
            dw=0.d0
            kkp1=kco(i)
            kkp2=kco(j)
            if (kkp1.eq.kkx.and.kkp2.eq.kkx) then
               dw = dsandwdcosij(p,costheta,kkx,kkx)
            else
               if (p.gt.0.d0) then
                  mp1 = mco(i)
                  mp2 = mco(j)
                  tmp=(s-(mp1-mp2)**2)*(s-(mp1+mp2)**2)/(4.0d0*s)
                  if (tmp.gt.0.0d0) then
                     pnew = sqrt(tmp)
                     f=mdof(i)*mdof(j)/4.d0
                     ! f=sum(gi*gj/g1^2) with same mass and wij
                     if (kkp1.ne.kkp2) f=2.0d0*f
                     dw=f*pnew/p
     &                    *dsandwdcosij(pnew,costheta,kkp1,kkp2)
                  endif
               endif
            endif
            w=w+dw
c            if (newfirstprint.and.prtlevel.eq.-222) then
c               write (*,*) 'dsandwdcos: coann ',
c     &              kkp1,kkp2,p,costheta,dw
c            endif
         enddo
      enddo
c      newfirstprint=.false.

      dsandwdcos=w
c      write(33,*) costheta,w    ! je test

      end
