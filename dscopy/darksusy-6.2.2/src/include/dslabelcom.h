*         -*- mode: fortran -*-

      integer npar     ! number of parameters
      real*8 par(100)    ! set of parameters (save integers as real*8)
      logical addlabel    ! option for authomatic generation of labels
      character*100 labelfile  ! file storing labels and associated parameters
c      character*10 label    ! label
      character*4 labelsuffix  ! suffix to authomatically generate a label
      common/dslabel2com/par,npar,addlabel,labelfile, !label,
     &  labelsuffix

      integer labin   ! option to pass parameters from labelinoutfun
      integer labout   ! option to pass parameters to labelinoutfun
      integer labset   ! option to pass parameters to labelinoutfun
      parameter(labin=1,labout=2,labset=3)

      integer labelunit
      parameter(labelunit=14)  !????      
      
      save /dslabel2com/
