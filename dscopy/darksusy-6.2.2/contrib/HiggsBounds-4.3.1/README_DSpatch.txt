CHANGES implemented to comply with f90 standards of gfortran >=6.3
======================================================

usefulbits.f90, S95table.f90, interpolate.f90, ...
--------------
* added additional space after stop:
  stop'error XYZ' -> stop 'error XYZ'

string_manip.f90
----------------
* make sure that integer that appears as length
  in array declaration is declared before array