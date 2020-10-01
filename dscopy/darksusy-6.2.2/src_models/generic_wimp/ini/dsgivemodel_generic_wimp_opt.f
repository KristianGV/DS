*******************************************************************************
***  subroutine dsgivemodel_generic_wimp_opt sets optional inputs not       ***
***  covered by the default set of input parameters given in                ***
***  dsgivemodel_generic_wimp (and overrides those specified there).        ***
***                                                                         ***
***  input:                                                                 ***
***                                                                         ***
***    key   - character string specifying input to be provided             ***
***    value - input value assigned (real)                                  ***
***                                                                         ***
***   This function should be called *after* dsgivemodel_generic_wimp (and) ***
***   before dsmodelsetup). For details about which values of "key" can be  ***
***   provided, see below.                                                  ***
***                                                                         ***
***                                                                         ***
*** author: Torsten Bringmann, torsten.bringmann@fys.uio.no                 ***
*** date: 2018-10-24                                                        ***
*******************************************************************************
      subroutine dsgivemodel_generic_wimp_opt(key,value)
      implicit none
      character*(*) key
      real*8 value
      include 'dsgeneric_wimp.h'
      include 'dsio.h'

      if (key.eq.'spin') then ! set DM spin; whether this is an acceptable 
                              ! assignment will only be tested after a new
                              ! call to dsmodelsetup
        spin(kdm)=value 
        kdof(kdm)=2*spin(kdm) + 1


      elseif (key.eq.'sigsip') then ! set spin-independent scattering  
                                    ! cross section on protons (pb)     
         sigsip = value*1.d-36 


      elseif (key.eq.'sigsin') then ! set spin-independent scattering  
                                    ! cross section on neutrons (pb)     
         sigsin = value*1.d-36 


      elseif (key.eq.'sigsdp') then ! set spin-dependent scattering  
                                    ! cross section on protons (pb)     
         sigsdp = value*1.d-36 


      elseif (key.eq.'sigsdn') then ! set spin-dependent scattering  
                                    ! cross section on neutrons (pb)     
         sigsdn = value*1.d-36 


      else 
        if (prtlevel.gt.0) write (*,*) 'dsgivemodel_generic_wimp_opt: '//
     &                    'unrecognized option ''',key,''' ''',value,''''
      endif
            
      return
      end 
