#include "include.h"

double chi2(int likelihood, char* names0[], double values0[], int nbobs0, struct parameters* param, int reset)
/* computes the chi2 of the observables whose names are contained in names0[], with predicted values in values0[]. The name of the observables can be found in experimental.h.
   nbobs0 is the number of input observables defining the length of the tables.
   likelihood gives the switch to use either B->K*0 mu mu results found using the likelihood method (=1) or the method of moments (=0).
   reset has to be set to 1 for the first run. In the subsequent runs, reset can be set to 0 to avoid recomputing the inverse covariance matrix and gain time, provided the observable list has not changed. */
{
	double prediction[NBOBS];
	static char* names[NBOBS];
	static double central[NBOBS],errors[NBOBS],covariance[NBOBS][NBOBS],invcov[NBOBS][NBOBS];
	static int onoff[NBOBS];
	
	int ie,je,ke,le;
		
	if(reset==1)
	{	
#include "experimental.h" /* contains the experimental measurements as well as the names of the available observables */

		for(ie=0;ie<NBOBS;ie++) onoff[ie]=0;
		for(ie=0;ie<NBOBS;ie++) prediction[ie]=0.;
					
		for(ie=0;ie<nbobs0;ie++)
		{
			je=0;
			while(strcmp(names0[ie],names[je])&&je<NBOBS-1) je++;
			if(!strcmp(names0[ie],names[je]))
			{
				prediction[je]=values0[ie];
				onoff[je]=1;
			}
			if((je==NBOBS-1)&&strcmp(names0[ie],names[je])) printf("Covariance matrix elements not available for %s                       Removed from chi2\n",names0[ie]);
		}
						
		for(ie=0;ie<NBOBS;ie++) 
		{
			if(errors[ie]==0.&&onoff[ie]==1)
			{
				onoff[ie]=0; 
				printf("No experimental data available for %s                                 Removed from chi2!\n",names[ie]);
			}
			
			if(likelihood==1&&onoff[ie]==1) if(ie>=116&&ie<=245)
			{
				onoff[ie]=0;
				printf("Using LHCb likelihood method, measurement not available for %s        Removed from chi2!\n",names[ie]);
			}
			
			if(likelihood==0&&onoff[ie]==1) if(ie>=26&&ie<=105)
			{
				onoff[ie]=0;
				printf("Using LHCb method of moments, measurement not available for %s        Removed from chi2!\n",names[ie]);
			}
			
			if(onoff[ie]==1) if((ie>=66&&ie<=75)||(ie>=186&&ie<=205)||(ie>=282&&ie<=283)||(ie>=311&&ie<=315)||(ie>=370&&ie<=375)||(ie>=441&&ie<=452))
			{
				onoff[ie]=0;
				printf("Theoretical prediction not reliable (due to resonances) for %s        Removed from chi2!\n",names[ie]);
			}
			
		}
		
		int nbobs=0; /* number of activated observables */
		for(ie=0;ie<NBOBS;ie++) nbobs+=onoff[ie]; 
		
		printf("n_obs=%d\n",nbobs);

		double covariance_reduced[nbobs][nbobs],inverted_covariance_reduced[nbobs][nbobs]; /* reduced covariance matrix containing only the activated observables */
		int matching[nbobs]; 

#include "covariance.h" /* contains the total covariance matrices */

		if(param->fullFF==0&&likelihood==0) for(ie=0;ie<NBOBS;ie++) for(je=0;je<NBOBS;je++) covariance[ie][je]=covariance_mom_soft[ie][je];
		if(param->fullFF==0&&likelihood==1) for(ie=0;ie<NBOBS;ie++) for(je=0;je<NBOBS;je++) covariance[ie][je]=covariance_like_soft[ie][je];
		if(param->fullFF==1&&likelihood==0) for(ie=0;ie<NBOBS;ie++) for(je=0;je<NBOBS;je++) covariance[ie][je]=covariance_mom_full[ie][je];
		if(param->fullFF==1&&likelihood==1) for(ie=0;ie<NBOBS;ie++) for(je=0;je<NBOBS;je++) covariance[ie][je]=covariance_like_full[ie][je];

		ke=0;
		for(ie=0;ie<NBOBS;ie++)
		{			
			if(onoff[ie]==1)
			{
				le=0;		
				for(je=0;je<NBOBS;je++)
				{
					if(onoff[je]==1)
					{
						covariance_reduced[ke][le]=covariance[ie][je];
						le++;
					}
				}	
				matching[ke]=ie;
				ke++;
			}
		}
		
		if(invert_matrix(nbobs,covariance_reduced,inverted_covariance_reduced)==0) return -1.; 
		
		for(ie=0;ie<NBOBS;ie++) for(je=0;je<NBOBS;je++) invcov[ie][je]=0.;
		for(ie=0;ie<nbobs;ie++) for(je=0;je<nbobs;je++) invcov[matching[ie]][matching[je]]=inverted_covariance_reduced[ie][je]; 
	}
	else
	{
		for(ie=0;ie<NBOBS;ie++) prediction[ie]=0.;
					
		for(ie=0;ie<nbobs0;ie++)
		{
			je=0;
			while(strcmp(names0[ie],names[je])&&je<NBOBS-1) je++;
			if(!strcmp(names0[ie],names[je]))
			{
				prediction[je]=values0[ie];
				onoff[je]=1;
			}
		}
	}
	
	double chi2=0.;
	for(ie=0;ie<NBOBS;ie++) for(je=0;je<NBOBS;je++) if(onoff[ie]*onoff[je]==1) chi2+=(prediction[ie]-central[ie])*invcov[ie][je]*(prediction[je]-central[je]); 

	return chi2;
}
