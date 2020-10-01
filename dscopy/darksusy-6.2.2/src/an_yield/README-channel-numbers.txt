In the new DS setup, we only keep simulated (SM) final states in src/. The
scalars are model dependent and that part is in src_models.

To facilitate this, we now use a new channel coding scheme where we instead use
exactly the same channel numbers in haloann also in src.

On top of this, the haloann channel numbers have now also been updated
to also include the lighter quarks. Hence, the following is a translation
the different channel numbers

Old chi is what is implemented in DS up until March 2014
New chi is going to be implemented in April 2014


*** Ch No  Particles                 Old chi   New chi
*** -----  ---------                 -------   ----------
***  1     S1 S1                     -         
***  2     S1 S2                     -
***  3     S2 S2                     -
***  4     S3 S3                     -
***  5     S1 S3                     -
***  6     S2 S3                     -
***  7     S- S+                     -
***  8     S1 Z                      -
***  9     S2 Z                      -
*** 10     S3 Z	                     -
*** 11     W- S+ and W+ S-           -
*** 12     Z0 Z0 	             6          9
*** 13     W+ W-                     5          8
*** 14     nu_e nu_e-bar             -
*** 15     e+ e-                     -
*** 16     nu_mu nu_mu-bar           -
*** 17     mu+ mu-                   7         10
*** 18     nu_tau nu_tau-bar	     -
*** 19     tau+ tau-	             4         11
*** 20     u u-bar                   -          2
*** 21     d d-bar                   -          1
*** 22     c c-bar                   1          4
*** 23     s s-bar                   -          3
*** 24     t t-bar                   3          6
*** 25     b b-bar                   2          5
*** 26     gluon gluon               8          7
*** 27     q q gluon (not implemented yet, put to zero)
*** 28     gamma gamma (1-loop)
*** 29     Z0 gamma (1-loop)
