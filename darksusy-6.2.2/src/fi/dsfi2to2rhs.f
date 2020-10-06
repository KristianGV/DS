c_______________________________________________________________________
c  Function dsfi2to2rhs gives the right hand side of the Boltzmann equation
c   s*dY/dt=RHS. Done by integrating dsfi2to2int over s, and multiplying by
c   temperature dependent factors in the end.
c   
c   input:
c     T - temperature of heath bath
c
c
c
c  author: Kristian Gjestad Vangsnes (kristgva@uio.no)   2020-09-07
c=======================================================================
      real*8 function dsfi2to2rhs(lnT)
      implicit none
      include 'dsmpconst.h'
      include 'dsficom.h'
      include 'dsidtag.h'



      real*8 lnT,T,a,b,eps,geff,sqrtgstar,heff,s_ent,H,HPrime,
     &dsf_int2,abserr,resabs,resasc,c,d,sum1,sum2,sum3,dsrdthav,
     &dsbessek2,x,debug1,debug2,debug3,sum,sum_tmp,tmp_upper,
     &tmp_upper1
      real*8, external :: dsfi2to2int_simp,dsanwx
      integer narray,ncoann,nrs,nthr,i,j
      parameter (narray=1000)
      real*8 mcoann(narray),dof(narray),tm(narray),
     &rm(narray),rw(narray)

      T=exp(lnT)
      call dsrdset('dof','default')
      call dskdgeff(T,geff)
      call dsrddof(T,sqrtgstar,heff)

      s_ent=heff*2*pi*pi/45*T*T*T
      H=sqrt(4*pi**3*geff/45)*T*T/mpl
      HPrime=H*heff/sqrt(geff)/sqrtgstar 

      T_22=T;x1_22=m1_22/T;x2_22=m1_22/T
      x=mdm/T
c     Integration limits
      ! a=1/(100*mdm**2);eps=1.d-5
      ! b=1/(4*mdm**2)
      ! b=1
      if(stat.eq.0) then
c     TRIED TO USE <sv>            
!             if(x.eq.0) then
!                   dsfi2to2rhs=0
!             else
!                   dsfi2to2rhs=T**2*dsrdthav(x,dsanwx)*
!      &(dsbessek2(x)/exp(x))**2/HPrime/s_ent
!             end if



            call dsrdparticles(0,narray,
     &  selfcon,ncoann,mcoann,dof,nrs,rm,rw,nthr,tm)
            a=(4*mdm**2)+1.d-30
            b=(4*mdm**2)+1.d80;eps=1.d-4
            sum=0
            if(nthr.eq.0) then
                  do i=0,100
                        tmp_upper=a+1.d4
                        tmp_upper1=tmp_upper
                        call dgadap(a,tmp_upper,dsfi2to2int_simp,
     &eps,sum_tmp)
                        sum=sum+sum_tmp
                        a=tmp_upper1
                  end do
            else
            do i =0,nthr
                  if(i.eq.nthr) then
                        tmp_upper=tm(nthr)**2
                  else
                        tmp_upper=tm(i+1)**2
                  end if
                  tmp_upper1=tmp_upper
                  if(nrs.eq.0) then
                        if(i.eq.0) then
                              tmp_upper=tm(1)**2
                              tmp_upper1=tmp_upper
                              
                        else if(i.ne.0.and.i.ne.nthr) then
                              tmp_upper=tm(i+1)**2
                              tmp_upper1=tmp_upper
                        else
                              tmp_upper=b
                              tmp_upper1=tmp_upper
                        end if
                        call dgadap(a,a+1.d2,
     &dsfi2to2int_simp,eps,sum_tmp)
                        sum=sum+sum_tmp
                        call dgadap(a+1.d2,tmp_upper,
     &dsfi2to2int_simp,eps,sum_tmp)
                        sum=sum+sum_tmp
                        a=tmp_upper1                                           
                  else   
                  do j =1,nrs
                    if((i.eq.0).and.(rm(j).ne.real(0))) then
                      if((rm(j).lt.tm(i+1))) then
                        tmp_upper=rm(j)**2-10
                        tmp_upper1=tmp_upper
                        call dgadap(a,tmp_upper,
     &dsfi2to2int_simp,eps,sum_tmp)
                        sum=sum+sum_tmp
                        a=tmp_upper1

                        tmp_upper=rm(j)**2+10
                        tmp_upper1=tmp_upper
                        call dgadap(a,tmp_upper,
     &dsfi2to2int_simp,eps,sum_tmp)
                        sum=sum+sum_tmp
                        a=tmp_upper1
                      else if(j.eq.nrs) then
                        tmp_upper=tm(1)**2
                        tmp_upper1=tmp_upper
                        call dgadap(a,tmp_upper,
     &dsfi2to2int_simp,eps,sum_tmp)
                        sum=sum+sum_tmp
                        a=tmp_upper1
                      end if

                    else if((i.gt.0).and.(i.lt.nthr+1)) then
                      if((rm(j).lt.tm(i+1)).and.
     &(tm(i).lt.rm(j))) then
                        tmp_upper=rm(j)**2-1000
                        tmp_upper1=tmp_upper
                        call dgadap(a,tmp_upper,
     &dsfi2to2int_simp,eps,sum_tmp)
                        sum=sum+sum_tmp
                        a=tmp_upper1

                        tmp_upper=rm(j)**2+1000
                        tmp_upper1=tmp_upper
                        call dgadap(a,tmp_upper,
     &dsfi2to2int_simp,eps,sum_tmp)
                        sum=sum+sum_tmp
                        a=tmp_upper1
                      else if(j.eq.nrs) then
                        if(i.eq.nthr) then
                              tmp_upper=tm(nthr)**2
                        else
                              tmp_upper=tm(i+1)**2
                        end if
                        tmp_upper1=tmp_upper
                        call dgadap(a,tmp_upper,
     &dsfi2to2int_simp,eps,sum_tmp)
                        sum=sum+sum_tmp
                        a=tmp_upper1
                      end if
                    else if(i.eq.nthr+1) then
                      if((rm(j).lt.b).and.
     &(tm(nthr).lt.rm(j))) then
                        tmp_upper=rm(j)**2-1000
                        tmp_upper1=tmp_upper
                        call dgadap(a,tmp_upper,
     &dsfi2to2int_simp,eps,sum_tmp)
                        sum=sum+sum_tmp
                        a=tmp_upper1

                        tmp_upper=rm(j)**2+1000
                        tmp_upper1=tmp_upper
                        call dgadap(a,tmp_upper,
     &dsfi2to2int_simp,eps,sum_tmp)
                        sum=sum+sum_tmp
                        
                        a=tmp_upper1
                      else if(j.eq.nrs) then
                        call dgadap(a,b,dsfi2to2int_simp,
     &eps,sum_tmp)
                        sum=sum+sum_tmp
                        a=tmp_upper1
                      end if
                    end if
                  end do
                  end if
            end do
      end if
c -------------------------------------------------------
            dsfi2to2rhs=sum/HPrime/s_ent*T
      else 
            a=(4*mdm**2)
            b=(4*mdm**2)+1.d40;eps=1.d-4

            call dgadap(a,b,dsfi2to2int_simp,eps,sum1)
            dsfi2to2rhs=(sum1)/HPrime/s_ent*T
      end if

      return
      end
