*******************************************************************
*** subroutine dsib2ffH3sVIB add to the common block array amp(i,j), i = 0, 
*** the helicity amplitudes of the s channel VIB diagrams. 
*** The H3 is emitted from the mediators (H3, Z), becoming H_{i}, i= 1, 2. H_{i}.
*** The sum is performed only over ih.
*** 
*** Author: Francesca Calore, 2014-02-20
*******************************************************************
      subroutine dsib2ffA0sVIB(f, ih) 
      implicit none
      include 'dsmssm.h'
      include 'dsib2com.h'

      integer f, ih

      integer kh(2),i 
      real*8 Mkz, Mh3, Mkw, Mhk(2)

      complex*16 deltaH3, deltaZ
           
      kh(1)=kh1
      kh(2)=kh2

c...  masses/mx mediators
      Mkz = mass(kz)/mx
      Mh3 = mass(kh3)/mx
      do i= 1,2
         Mhk(i) = mass(kh(i))/mx
      end do
      Mkw = mass(kw)/mx

c... denominators
c... H3 
      deltaH3 = (dcmplx(4*E1J**2 - Mhk(ih)**2,Mhk(ih)*width(kh(ih))/mx)*
     -     dcmplx(4*E1J*(E1J + EvJ),Mh3*width(kh3)/mx))
c... Z
      deltaZ = (dcmplx(4*E1J**2 - Mhk(ih)**2,Mhk(ih)*width(kh(ih))/mx)*
     -     dcmplx(4*E1J*(E1J + EvJ) - Mkz**2 + MB**2,Mkz*width(kz)/mx))

c... helicity amplitudes
      amp(0, 2) = amp(0, 2) + (
     -     (0,-4)*kJ*Mkw*gl(kh(ih),kh3,kh3)*
     -     (Conjg(gr(kh3,kn(1),kn(1))) - gr(kh3,kn(1),kn(1)))*
     -     gr(kh(ih),f,f)/deltaH3
     -     )                    ! H3 contribution
     -     + (                  ! Z contribution
     -     ((0,2)*kJ*(-4*E1J**2 + MB**2)*(4*E1J*(E1J + EvJ) - Mkz**2 + MB**2)*
     -     gl(kz,kh3,kh(ih))*(Conjg(gl(kz,kn(1),kn(1))) + gl(kz,kn(1),kn(1)))*
     -     gr(kh(ih),f,f))/Mkz**2/deltaZ
     -     )

      return
      end


