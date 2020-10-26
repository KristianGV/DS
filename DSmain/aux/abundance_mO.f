      program abundance_mO
c----------------------------------------------------------------------
        implicit none
        include 'dsmpconst.h'
        include 'dsparticles.h'
        include 'dssm.h'
        include 'dsidtag.h'
        include 'dsficom.h'

        real*8 oh2,dsfi2to2oh2,TR,Tmin,inputlambda,dsrdomega,xf,omega

        integer ierr,iwar,nfc

        mdm=130.d0; TR=1.d6; Tmin=1.d-5; inputlambda=1.d-11

        stat=0; eta1_22=1;eta2_22=1;etaX_22=0;g1_22=1;g2_22=1
        c12_22=1.d0
        call dsinit

        call dsgivemodel_silveira_zee(inputlambda,mdm)
        call dsmodelsetup()
        ! omega=dsrdomega(0,0,xf,ierr,iwar,nfc)

        oh2=dsfi2to2oh2(Tmin,TR,m1_22,m2_22
     &,eta1_22,eta2_22,etaX_22,g1_22,g2_22,c12_22)

        write(*,*) "Oh2=", oh2
      end program
