In order to work smoothly together with DarkSUSY, we had in
some cases to implement minor changes to the contributed codes
contained in this directory.
[In cases, where more non-trivial changes were necessary, like
for Isajet, the respective subdirectory name contains the string 
"-for-darksusy"]

All these minor changes are marked with the string "DS mod", and thus
easily searchable.

The following provides a summary of which codes are affected, and what
has been changed:


HiggsSignals-1.4.0
==================
* suppressed excessive writing to standard output 
  (introduced a 'debug' flag which by default is set to false)

