      subroutine dsddset_mssm(cs)
c...set parameters for mssm scattering cross section
c...  cs - character string specifying choice to be made
c...author: paolo gondolo 2000-07-07
c...modified by Gintaras Duda 2007-06-27 for new FF options
c...modified by Paolo Gondolo 2008-02-18
c...modified by Paolo Gondolo 2013-11-29
c...modified by torsten bringmann 2018-06-13: changed default, and allowed to set  
c...                                          poles and DN independently
      implicit none
      include 'dsddmssmcom.h'
      character*(*) cs

      if (cs.eq.'help') then
         write (*,*) 
     &        'dsddset_mssm: use "call dsddset_mssm(value)" where'
         write (*,*) ' value=''help'' this help'
         write (*,*) ' value=''tree_pole'' for full tree-level ',
     &        'scattering amplitude'
         write (*,*) ' value=''tree_nopole'' to skip squark pole in ',
     &        'scattering amplitude as in Griest 1988 and Gondolo'
         write (*,*) ' value=''dn_pole'' to treet QCD corrections to ',
     &        'scattering amplitude as in Drees & Nojiri 1993'
         write (*,*) ' value=''dn_nopole'' to treet QCD corrections to ',
     &        'scattering amplitude as in Drees & Nojiri 1993, but explicitly ',
     &        'remove the poles [default]'
      else if (cs.eq.'dn_nopole'.or.cs.eq.'default') then
         dddn = 1
         ddpole = 0
      else if (cs.eq.'tree_pole') then
         dddn = 0
         ddpole = 1
      else if (cs.eq.'tree_nopole') then
         dddn = 0
         ddpole = 0
      else if (cs.eq.'dn_pole') then
         dddn = 1
         ddpole = 1
      else
         write (*,*) 'dsddset_mssm: unrecognized options ',cs
      endif
      
      return
      end
