c...Subroutine dssenu_set
c...set parameters for neutrino telescope routines
c...  c - character string specifying choice to be made
c...author: joakim edsjo, 2000-08-16
c...modified: joakim Edsjo, 2015      
c...Short description of options set here:
c...  secalcmet = 1 -> JKG approximation of capture rates
c...              2 -> JKG + similar approx for Earth, but w/ Gould
c...              4 -> Full numerical integration over velocity distribution
c...                   and sum over solar and terrestrial elements
c...              5 -> Full numerical integration over velocity distribution
c...                   and form factors (as taken from DD routines) and
c...                   sum over solar and terrestrial elements
c...              (3 -> Deprecated: Damour-Krauss population)
c...  setab = 0 -> no tables are used, instead the capture rates are
c...               calculated (integrated) directly
c...          1 -> the capture rates are read from tables (created from 
c...               numerical integrations) to speed things up. This option
c...               is not available for all scenarios (e.g. full Jupiter
c...               cut-off, sejup=3, where numerical integrations need to
c...               be performed model for model).
c...  sejup = 0 -> no suppression due to Jupiter
c...          1 -> simple approximate of effects of Jupiter on solar capture
c...               if the WIMPs reach out to 5.2 AU after first scatter they
c...               are considered lost as Jupiter will most likely throw the
c...               WIMPs out of the solar system before being able to scatter
c...               in the Sun again.
c...          2 -> more complete treatment of Jupiter effects on solar capture.
c...               The maximal distance from the Sun where WIMPs are considered
c...               captured after first scatter is here calculated as a function
c...               of the scattering cross sections, based on results by
c...               Annika Peter. For very high cross sections, Jupiter effects
c...               are smaller than option 1, whereas the effects are larger
c...               for small cross sections. This option requires numerical
c...               integrations model by model.
c...          sejup options are available for secalcmet 4 and 5.
c...  sesunacc = 1 -> high accuracy for solar integration, i.e. all elements
c...                  and isotopes are included for both SI and SD scattering.
c...             2 -> medium accuracy for solar integration, i.e. 
c...                  elements with very low abundances are ignored and some
c...                  heavier isotopes are clumped together element by element
c...             3 -> low accuracy for solar integration, i.e. elements with
c...                  low abundances are ignored and some
c...                  heavier isotopes are clumped together element by element
c...    In most cases, the low option is accurate enough (to within 1%).
c...  Note that the solar model used is set with a call to the routine
c...  dssem_sunset.
      subroutine dssenu_set(c)
      implicit none
      include 'dssecom.h'
      character*(*) c

c...Initialize tables
c      call dssem_sunread   ! this is now done in the routines that need it

c...Defaults
      sejup=0     ! no Jupiter effect

c...approximate formulae given in jungman, kamionkowski, griest.
      if (c.eq.'jkg') then
         secalcmet=1           ! calculate method

c...approximate as above for the sun, but the full formulae from
c...gould, apj 321 (1987) 571 for the earth
      else if (c.eq.'gould') then
         secalcmet=2

c...Use full expressions in Gould, ApJ 321 (1987) 571 for the Earth, but
c...instead of using analytic expressions for a Gaussian and approximate
c...positions within the Earth for the elements, the velocity profile
c...specified in dshmset.f (i.e. by the parameter veldfearth) is
c...integrated numerically and a full integration over the Earth's
c...interior is performed. For the Sun, the same full expressions are used
c...and the numerical profile as specified in dshmset.f (i.e. by the 
c...parameter veldf) is integrated numerically over a realistic Sun
c...profile.
      else if (c.eq.'num') then
         secalcmet=4
         setab=0
         sesunacc = 2 ! medium solar model accuracy

      else if (c.eq.'numcut') then
         secalcmet=4
         setab=0
         sesunacc = 2 ! medium solar model accuracy
         sejup=1
c...This option includes the effects of Jupiter by assuming that WIMPs that
c...reach out to Jupiter after first scatter will not be captured.

      else if (c.eq.'numcutfull') then
         secalcmet=4
         setab=0
         sesunacc = 2 ! medium solar accuracy
         sejup=2
c...This is the same as numcut above, i.e. that effects of Jupiter is
c...taken into account. This option is more accurate though in that it
c...determines the effects of Jupiter as a function of WIMP mass and
c...scattering cross section. The higher the cross section, the smaller
c...the effects of Jupiter are. The effects are typically larger than
c...for the simplified Jupiter approach above in numcut. As this option
c...requires integration over the velocity distribution for every single 
c...combination of WIMP mass and cross section, there is no corresponding
c...fast option available. veout is determined for each model in
c...dssenu_capsunnum/dssenu_capsunnumff


c...Same as 'num' or 'numcut' above, but use tabuled versions of the
c...numerical integrations.
c...If the tables do not exist, they are recreated with the numerical
c...routines as in option 'num' above.

      else if (c.eq.'numhi') then
         secalcmet=5
         setab=0
         sesunacc = 1 ! high solar model accuracy

      else if (c.eq.'numcuthi') then
         secalcmet=5
         setab=0
         sesunacc = 1 ! high solar model accuracy
         sejup=1

      else if (c.eq.'numlo') then 
         secalcmet=5
         setab=0
         sesunacc = 3 ! low solar model accuracy

      else if (c.eq.'numcutlo') then
         secalcmet=5
         setab=0
         sesunacc = 3 ! low solar model accuracy
         sejup=1

      else if (c.eq.'nummed') then
         secalcmet=5
         setab=0
         sesunacc = 2 ! low solar model accuracy

      else if (c.eq.'numcutmed') then
         secalcmet=5
         setab=0
         sesunacc = 2 ! low solar model accuracy
         sejup=1

      else if (c.eq.'numcutfullhi') then ! more accurate Jupiter capture
         secalcmet=5
         setab=0
         sesunacc = 1 ! high solar model accuracy
         sejup=2

      else if (c.eq.'numcutfullmed') then ! more accurate Jupiter capture
         secalcmet=5
         setab=0
         sesunacc = 2 ! low solar model accuracy
         sejup=2

      else if (c.eq.'numcutfulllo') then ! more accurate Jupiter capture
         secalcmet=5
         setab=0
         sesunacc = 3 ! low solar model accuracy
         sejup=2

      else if (c.eq.'numold') then ! old routines with default FR
         secalcmet=4
         sesunacc=2
         setab=0

      else if (c.eq.'tabhi') then
         secalcmet=5
         setab=1
         sesunacc = 1 ! high solar model accuracy

      else if (c.eq.'tabcuthi') then
         secalcmet=5
         setab=1
         sesunacc = 1 ! high solar model accuracy
         sejup=1

      else if (c.eq.'tablo') then 
         secalcmet=5
         setab=1
         sesunacc = 3 ! low solar model accuracy

      else if (c.eq.'tabcutlo') then
         secalcmet=5
         setab=1
         sesunacc = 3 ! low solar model accuracy
         sejup=1

      else if (c.eq.'tabmed'.or.c.eq.'default') then 
         secalcmet=5
         setab=1
         sesunacc = 2 ! low solar model accuracy

      else if (c.eq.'tabcutmed') then
         secalcmet=5
         setab=1
         sesunacc = 2 ! low solar model accuracy
         sejup=1

c...Same as above, but use tabuled versions of the numerical integrations.
c...If the tables do not exist, they are recreated with the numerical
c...routines as in option 'num' above.

      elseif (c.eq.'tab') then

        secalcmet=4
        setab=1

      else if (c.eq.'tabcut') then
         secalcmet=4
         setab=1
         sejup=1

c...invalid choice
      else
         write (*,*) 'dssenu_set: unrecognized option ',c
         stop
      endif

      return
      end




