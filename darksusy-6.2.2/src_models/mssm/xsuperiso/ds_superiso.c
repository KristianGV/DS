/* Small set of header functions that make calling superiso functions */
/* from Fortran easy                                                  */
/* Author: Joakim Edsjo, edsjo@fysik.su.se                            */
/* Date: 2013-05-01                                                   */

/*void myfunc(int *a,int *b);   */
/*extern "C" void myfunc_(int *a,int *b) { return myfunc(a,b);} */

int test_slha(char *name[]);
extern int test_slha_(char *name[]) { return test_slha(name);}

double bsgamma_calculator(char *name[]);
extern double bsgamma_calculator_(char *name[]) { return bsgamma_calculator(name);}

double delta0_calculator(char *name[]);
extern double delta0_calculator_(char *name[]) { return delta0_calculator(name);}

double Bsmumu_calculator(char *name[]);
extern double bsmumu_calculator_(char *name[]) { return Bsmumu_calculator(name);}

double Bsmumu_untag_calculator(char *name[]);
extern double bsmumu_untag_calculator_(char *name[]) { return Bsmumu_untag_calculator(name);}

double Btaunu_calculator(char *name[]);
extern double btaunu_calculator_(char *name[]) { return Btaunu_calculator(name);}

double BDtaunu_calculator(char *name[]);
extern double bdtaunu_calculator_(char *name[]) { return BDtaunu_calculator(name);}

double Rmu23_calculator(char *name[]);
extern double rmu23_calculator_(char *name[]) { return Rmu23_calculator(name);}

double muon_gm2_calculator(char *name[]);
extern double muon_gm2_calculator_(char *name[]) { return muon_gm2_calculator(name);}

