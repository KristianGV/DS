      subroutine dsdvmatch(itypein,nrein,rein,itypeout,nreout,reout)
c_______________________________________________________________________
c subroutine to match the input rein(nrein) dependent variables of type
c itypein onto the output reout(nreout) dependent variables of type
c itypeout.
c the possible matches are: 
c   itypeout = itypein -> you just copy rein(nrein) into reout(nreout)
c   itypeout spherical from itypein axisymmetric or triaxial
c   itypeout axisymmetric from itypein triaxial
c in all other cases a error message is returned and the program stops      
c_______________________________________________________________________
      implicit none
      include 'dsdvcom.h'
      integer itypein,nrein,itypeout,nreout 
      real*8 rein(nrein),reout(nreout)
      integer ii
ccc
      if(itypeout.eq.itypein) then
c
c just copy the input to the output        
c    
        do ii=1,nrein
          reout(ii)=rein(ii)
        enddo
        return
      endif
ccc      
      if(itypeout.eq.itydvsph) then
c
c if you are calling for a spherically symmetric function with
c axisymmetric or triaxial coordinates:       
c
        if(itypein.eq.itydvaxi) then
          reout(idvrsph)=dsqrt((rein(idvraxi))**2+(rein(idvzaxi))**2)
          return
        elseif(itypein.eq.itydvtri) then
          reout(idvrsph)=dsqrt((rein(idvx))**2+(rein(idvy))**2
     &     +(rein(idvz))**2)
          return
        endif
      elseif(itypeout.eq.itydvaxi) then
c
c if you are calling for a axisymmetric symmetric function with
c triaxial coordinates:       
c
        if(itypein.eq.itydvtri) then
          reout(idvraxi)=dsqrt((rein(idvx))**2+(rein(idvy))**2)
          reout(idvzaxi)=rein(idvz)
          return
        endif
      endif
c
c otherwise you called this subroutine without a possible match
c      
      write(*,*) 'DS: in dsdvmatch, wrong itypein - itypeout match:'
      write(*,*) 'DS: itypein = ',itypein
      write(*,*) 'DS: itypeout = ',itypeout
      stop
      end
