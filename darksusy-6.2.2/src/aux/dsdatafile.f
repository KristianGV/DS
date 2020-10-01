**********************************************************************
*** subroutine dsdatafile generates a filename for files stored
*** in the data/ dirctory in the DarkSUSY root.
*** This routine is useful to get absolute filenames to be used
*** on a given system.      
*** Output: filename, full absolute path to data file
*** Input: filename_in (name of file in data/)
**********************************************************************      
      subroutine dsdatafile(filename_out,filename_in)
      character*(*) filename_out,filename_in
      integer dsi_trim
      include 'dsdir.h'
      filename_out=dsdatapath(:dsi_trim(dsdatapath))//
     &     filename_in
      return
      end
