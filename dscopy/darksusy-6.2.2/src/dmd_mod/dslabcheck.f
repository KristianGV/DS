      subroutine dslabcheck(nlabin,labin,nlabsuff,labsuff,match)
c_______________________________________________________________________
c subroutine to search fot the character string 'labsuff' made of
c nlabsuff characters within the label 'labin' made of nlabin characters  
c if the string is found match=.true. is returned, otherwise
c match=.false. is returned
c_______________________________________________________________________
      implicit none
      character(*) labin,labsuff
      integer nlabin,nlabsuff
      logical  match
      integer i
      match=.false.
      i=0
 10   if(labin(1+i:nlabsuff+i).eq.labsuff) then
        match=.true.
        return
      endif  
      i=i+1
      if(nlabsuff+i.le.nlabin) goto 10
      return
      end
