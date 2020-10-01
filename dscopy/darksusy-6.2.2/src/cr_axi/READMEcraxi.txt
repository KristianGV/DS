
[TB&JE] Update after long discussion: get rid of all common block logical 
parameters and rely instead fully on the concept of replaceable functions.
Move all functions from /templates to /user_replaceables to illustrate concept



-----------------------------------------


This subroutine contains the propagation of charged cosmic ray in an
axisymmetric two-zone propagation model (spatially constant diffusion
coefficient and energy loss term, ect. ect. - details of the model are 
in the accompaining tex file).

Steps in running these routines:

1) Define the axisymmetric profile of the source function: 
this is done by setting the function "dsnpaxi", and setting a 
truncation radius defining a spherical region inside which the source 
function is assumed to be constant. You can implement such truncation 
in view of the fact that it is unlikely that the approximation of 
spatially constant diffusion coefficient or positron energy loss term 
holds in the innest part of the galaxy. Two different truncation 
radius are allowed for pb/db and ep, respectively, "diffrcpb" and 
"diffrcep", stored in the common block "dscraxicom.h"; in practice
this matters only for very singular source functions. The source 
function is passed to the pb/db and ep routines through, respectively, 
"dspbnpairsaxi" and "dsepnpairsaxi". 

NOTE: dsnpaxi can be used to implement the most generic source 
function, including the case of dark matter pair annihilations and 
dark matter decay; the user should keep track only of factorizations
and the appropriate units. The recommended factorizations are:

- for wimp pair annihilations: dsnpaxi is the number density of wimp 
pairs, up to the factor 0.5d0/hamwimp**2, in units of GeV^2 cm^-6

- for dark matter decays: dsnpaxi is the number density of DM 
particles, up to the factor  1/hamwimp, in units of GeV cm^-3

2) Enter the propagation parameters: 

- All routines in this directory assumes a spatially constant diffusion 
coefficient, given by the function "dskdiff". The antiproton routines 
allow for a different diffusion coefficient in the gas disc; otherwise 
one can implement a generic dependence on rigidity and beta (=v/c). 
Setting the "dscraxicom.h" logical variable "kdiffdef" to .true., the 
fuction "dskdiff" links to a set of default hardcoded versions of the 
diffusion coefficient given in the function "dskdiff_def", which are
defined in terms of a set of parameters defined in the common block
"dscraxicom.h". If "kdiffdef" is set to .false. the fuction "dskdiff" 
links to the user provided function "dskdiff_user". NOTE: when linking
to dskdiff_user, if the user wants to have tabulations for 

- The other parameters defining the propagation volume and convection 
are all include in the common block "dscraxicom.h" and are: 
a) the half thickness of the diffusive halo "diffhh";
b) the radial size of the diffusive halo "diffRh";
c) the half thickness of the of the gas disc "diffhg";
d) the density of hydrogen in the disc "diffng";
e) the density of hydrogen in the halo "diffnh";
f) the galactic wind velocity, assumed to be constant and outgoing in 
   the whole halo.
Of these parameters only "diffhh" enters in the computation of the
ep diffusion, since in the ep propagation model the radial boundary is
neglected, there is no gas disc, the galactic wind cannot be included.

- The energy losses for ep/em are defined through the function 
"dseppdotm", which assumes a spatially constant energy loss rate. Two
option can be implemented: 
g) the energy loss is assumed to be in the approximate form: 
blossmean * pp^2, with pp the momentum of the ep/em, and blossmean the
mean energy loss rate at 1 GeV;
h) the energy loss links to a generic function of momentum with all 
energy losses terms - NOTE: at the moment this is not actually fully
implemented and only IC and synchrotron losses are included so far; 
these two terms are specified by the mean starlight density 
"starlightmean", and the mean magnetic field "bfieldmean" - ON MY
LIST OF THINGS TO IMPROVE.

3) Define preliminary steps and options to propagate ep/em and pb/db.

- The computation of the ep/em equilibrium number density is 
implemented thorugh the spatial integration of a green function, which
depends only on the ep/em diffusion length $\lambda_e$. This function is
tabulated and its tabulations needs to be reloaded whenever one changes:
a) the source function (step 1 above); 
b) the half thickness of the diffusive halo "diffhh" (one of the 
parameters in step 2 above) and the inner radial cutoff "diffrcep"; 
c) the radial/vertical coordinate where ep/em are computed.
Three labels are assigned, respectively, to a), b), and c): "nptag",
"epcraxitag", and "rzaxitag". For the latest two, labels and 
corresponding parameters are stored in repository files - respectively, 
"axirzlabel.dat" and "axieplabel.dat"; a new label can be added 
to the repositories with the dsaddlabel subroutines. In the subroutine 
"dsepgreenaxitab", for a given set of parameters, the corresponding 
labels are extracted from repository files and used for defining the
file from which to read (or on which to write) the corrresponding 
tabulation. If the currently used parameters do not match any of the
combinationations in the repository files, there are two possibilities,
depending on whether the user allows or not (setting the corresponding 
logical variable, respectively, to .true. or .false.) to modify 
interactively repository files. If you do not allow this possibility,  
an error message is printed and the program stopped, if you do allow 
it, a label is automatically created glueing in sequence a "a" (for
authomatic) a suffix (up to 4 characters) which is, eventually, defined 
by the user, and the number corresponding to the line number at which
the label will appear in the repository. If you make sure that the 
labels the user is adding to repositories do not start with a "a", you
can use the subroutine "dscleanlabel" to clean the a repository from 
automatically generated label and move them to another file (the same
routine can be used to take out any given label). WARNING: the label 
'nptag' and its number of characters 'nnptag' needs to be fixed by the 
user (in the common blocks of dsdmdcom.h) and there is no internal
consistency check on whether the label really corresponds to the 
currently used source function.

- The computation of the pb/db confinement time is performed by the 
functions "dspbtdaxi" and "dsdbtdaxi". They involve a double numerical
integral and a sum over zeros of bessel functions; if the input 
parameter "isin" is set to 1 the computation is performed assuming that 
the gas disc can be approximated with a delta function, if it is set 
to 2 the propagation model is the two-zone model with finite thickness 
of the disc - the difference in results between the two is in most cases
marginal, the second option is slightly more CPU consuming. The 
structure of the computation is through a nested set of functions:
  i) "dspbint_intgz" ("dspbinte_intgz") integrand in the z direction 
     for isin=1 case (isin=2 case);
  ii) "dspbint_intz" ("dspbinte_intz") function performing the integral
     in the z direction for isin=1 case (isin=2 case);
  iii) "dspbint_intgRis1" ("dspbint_intgRis2") integrand in the R 
     direction for isin=1 case (isin=2 case);
  iv) "dspbint_intgR2" function linking to either "dspbint_intgRis1" or
     "dspbint_intgRis2";
  v) "dspbintrs" function performing the integral in the R direction;
  vi) "dspbsums" function performing the sum over zeros of the bessel 
     function J0.
This set of fuction is the same for pb and db, the only difference is 
the link to the appropriate cross section for particle losses via
collisions and the appropriate definition of rigidity with the input
being the kinetic energy per nucleon. There is a little complication when
dealing with the computation, which is related to the fact that for 
source functions that are large in the central part of the system, or 
for the case when the position at which the confinement time is computed 
get close to the origin in the coordinate system, the sum of the over 
the zeros of the bessel function J0 oscillates around a mean, converging 
very slowly.

NOTE: "The convergence is checked by summing over one full 
oscillation and comparing the estimated mean to the one in the previous 
cycle; if these differ of less than the input parameter "prec", the last 
estimated mean value is returned as the output. If convergence is not found 
within a number of zeros correponding to the input parameter "nprec", a 
cylinder at the center of the system is cut out of the integral in the R 
and z directions, and replaced by a spatial integral of the point source 
green function (computed in the limit of no radial boundary conditions),
appropriately weighted over the source function in the same region; to 
do the latter, a series of tabulations for the green function
is needed, a procedure which is anyway CPU consuming. Also, once the
method of cutting out a central part is activated, it is mantained for
any following call (since, e.g. subsequent calls with the same source
function but different kinetic energies, are likely to have the same
convergence problem); in case the user is making subsequent calls with
very different source functions, it is advisable to reinitialize the
parameters defining the inner cylinder, and this can be done by hand
with a call to the subroutine "dspbcutini". Also, for singular source 
function, a fair guess is that at least one iteration of the inner 
cylinder setting will be needed and this can be done by the user with 
a call to the subroutine "dspbcylcut"."

IS NOW REPLACED BY: "The convergence is Improved by defining recursively
averages of cycles of oscillations, with a slight variant of the method
proposed in the PhD Thesis of Torsten Bringmann, pag. 72, Eq. 5.33, and 
checking the stability of such averages."
 

- The function "dspbtdaxi" and "dsdbtdaxi" are tabulated, respectively, 
by the functions "dspbtdaxitab" and "dsdbtdaxitab" (the tabulation goes 
from 0.1 GeV to 10 TeV in log scale, starting with 100 points, and 
adding points if values between two nearby points differ more than 20%
- both the 100 and the 20% at this stage are hardcoded in the code). These 
functions needs to be computed or reloaded from disk whenever one changes:
a) the source function (step 1 above) - the user provided label 'nptag'
needs to be properly set, as for the ep case above;
b) the radial/vertical coordinates corresponding to the point where 
pb/db are computed - the label "rzaxitag" takes care of keeping track 
of this as for the ep case above;
c) any of the parameters setting propagation in the axisymmetric model -
the label pbcraxitag keeps track of this;
d) the option isin is changed - the label islab keeps track of this.


4) compute the equilibrium number density and the flux.

- For ep/em the equilibrium number density is obtained by linking to the
function "dsepdndpaxi" up to a normalization factor (different for 
annihilating or decaying dark matter). At this level, you need to have 
specified the ep/em injection spectrum; this is defined in the function 
"dsepdndpi", where, if the option "epspecdef" is set to .true. there is
a link to the default function "dsepdndpi_def", which is just the wimp 
pair halo annihilation function "dshaloyield"; if the option "epspecdef"
is set to .false., "dsepdndpi" links to the user provided function 
"dsepdndpi_user" (in the standard DS release there is just a dummy 
version of this function. Note: to use these routines for the case of 
dark matter particles decays in the halo, the user needs to provide an
appropriate injection spectrum). In the function "dsepdndpaxi" you need 
also to choose the prescription to link momentum to the diffusion length
scale. For "ivopt=1" this is done through the generic (but spatially
homogeneous) ep/em energy loss rate function "dseppdotmmean" and any
given diffusion coefficient, using on the first call a numerical 
integrations and then linking to a tabulation; this tabulation is 
stored in the file "vofpnum.dat" and reloaded in "dsepdndpaxi" if the
energy loss rate or the diffusion coefficient have changed at 1 GeV or
50 GeV, and beta=0.5d0 - note: if you need to reload the tabulation by
hand you have to set both the variable "vofpset" in the common block
"vofenumstore" to zero and remove the file "vofpnum.dat". For "ivopt=2",
the diffusion length is computed using analytical integrals and 
inversions which are possible only for a specific form of the momentum 
loss rate, i.e. function "dseppdotmana" (see point g in step 2 above)
and a subset of the diffusion coefficients given in the default 
hardcoded version of the diffusion coefficient as given in "dskdiff_def"
- namely those corresponding to "nkdiff" = 1 or 4 and "betalabel" =
.false. (these are parameters in "dscraxicom.h"); a sample check that 
this is the case is inserted in the sample setup file "dscraxiset", 
however the program may complete the run even for an inconsistently-
defined "ivopt=2" setup, it is responsibility of the user to do a cross 
check if "ivopt=2" is used.

- The ep/em differential flux at the local galactocentric radius 
R=objdistance,z=0 is computed by the function: 
  "dsepdphidphaaxi" if the computation refers to wimp pair annihilations
                    in the halo; 
  "dsepdphidpdeaxi" if the computation refers to dark matter particle
                    decays in the halo.
Note: these functions assumes the recommended factorizations for 
dsnpaxi; if different factorizations are used, the user should 
implement his own link to "dsepnpairsaxi".

- The pb/db differential flux at the local galactocentric radius 
R=objdistance,z=0 is computed by the function: 
  "dspbdphidthaaxi/dsdbdphidthaaxi" if the computation refers to wimp 
                    pair annihilations in the halo; 
  "dspbdphidthdaxi/dsdbdphidthdaxi" if the computation refers to dark 
                    matter particle decays in the halo; the user needs to
                    set dspbdndti/dsdbdndti appropriately.
Note: as abpve these functions assumes the recommended factorizations for 
dsnpaxi; if different factorizations are used, the user should 
implement his own link to "dspbnpairsaxi" (this is assumed to be valid 
for both pb and db).
