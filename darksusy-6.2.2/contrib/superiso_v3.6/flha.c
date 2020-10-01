#include "src/include.h"
/*--------------------------------------------------------*/
/* Calculation of the observables using a given FLHA file */
/*--------------------------------------------------------*/

int main(int argc,char** argv)
{ 
	struct parameters param;
	struct flhaparam paramflha;
	char name[50];
	int test;
		
  	if(argc<2) 
  	{ 
    		printf(" This program needs 1 parameter:\n"
           	"   name    name of the FLHA file\n");
      		exit(1); 
  	} 
	else 
  	{
  		sscanf(argv[1],"%s",name);
  	}
	
	printf("\n");
	
	printf("SuperIso v3.6 - F. Mahmoudi\n\n");
	printf("FLHA input file\n\n");

	Init_param(&param);
	
	Les_Houches_Reader(name,&param);
	slha_adjust(&param);
	param.SM=1;

	FLHA_Reader(name,&paramflha);
	
	printf("\n");
	
	printf("Observable\t\t\tValue\n\n");

	double C0b[11],C0spec[11],C1b[11],C1spec[11],C0w[11],C1w[11],C2w[11],C2b[11],Cpb[11];
	double complex CQpb[3],CQ0b[3],CQ1b[3];
	CQ0b[1]=CQ0b[2]=CQ1b[1]=CQ1b[2]=CQpb[1]=CQpb[2]=0.;
	double obs[Nobs_BKsll+1];

	double mu_W=2.*param.mass_W;
	if(paramflha.Q>50.) mu_W=paramflha.Q;
		
	CW_calculator(2,C0w,C1w,C2w,mu_W,&param);
	
	C0w[7]+=paramflha.dC7[0];
	C1w[7]+=paramflha.dC7[1];
	C2w[7]+=paramflha.dC7[2];
	C0w[8]+=paramflha.dC8[0];
	C1w[8]+=paramflha.dC8[1];
	C2w[8]+=paramflha.dC8[2];
	C0w[9]+=paramflha.dC9[0];
	C1w[9]+=paramflha.dC9[1];
	C2w[9]+=paramflha.dC9[2];
	C0w[10]+=paramflha.dC10[0];
	C1w[10]+=paramflha.dC10[1];
	C2w[10]+=paramflha.dC10[2];

	
	double mu_b=2.;
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	Cprime_calculator(2,Cpb,CQpb,mu_W,mu_b,&param);
	Cpb[7]+=paramflha.dC7p[0];
	Cpb[8]+=paramflha.dC8p[0];
	Cpb[9]+=paramflha.dC9p[0];
	Cpb[10]+=paramflha.dC10p[0];
	printf("BR(b->s gamma)\t\t\t%.3e\n",bsgamma(C0b,C1b,C2b,Cpb,mu_b,mu_W,&param));
	
	double lambda_h=0.5;
	double mu_spec=sqrt(lambda_h*param.mass_b);		
	C_calculator_base2(C0w,C1w,mu_W,C0b,C1b,mu_b,&param);
	C_calculator_base2(C0w,C1w,mu_W,C0spec,C1spec,mu_spec,&param);
	printf("delta0(B->K* gamma)\t\t%.3e\n\n",delta0(C0b,C0spec,C1b,C1spec,Cpb,&param,mu_b,mu_spec,lambda_h));
	
	mu_b=param.mass_b_pole;
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	Cprime_calculator(2,Cpb,CQpb,mu_W,mu_b,&param);
	Cpb[7]+=paramflha.dC7p[0];
	Cpb[8]+=paramflha.dC8p[0];
	Cpb[9]+=paramflha.dC9p[0];
	Cpb[10]+=paramflha.dC10p[0];
	
	param.SM=0;
	CQ_calculator(2,CQ0b,CQ1b,mu_W,mu_b,&param);
	param.SM=1;
	
	printf("BR(Bs->mu mu)_untag\t\t%.3e\n",Bsmumu_untag(C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b));
	printf("BR(Bd->mu mu)\t\t\t%.3e\n",Bdmumu(C0b,C1b,C2b,CQ0b,CQ1b,&param,mu_b));

	param.SM=0;
 	printf("BR(B->tau nu)\t\t\t%.3e\n",Btaunu(&param));
    printf("R(B->tau nu)\t\t\t%.3e\n\n",RBtaunu(&param));
    printf("BR(B->D tau nu)\t\t\t%.3e\n",BRBDlnu_full(3,1,obs,&param));
    printf("AFB(B->D tau nu)\t\t%.3e\n",obs[1]);
    printf("Ptau(B->D tau nu)\t\t%.3e\n",obs[2]);
	printf("BR(B->D tau nu)/BR(B->D e nu)\t%.3e\n\n",BDtaunu_BDenu(&param));
    printf("BR(B->D* tau nu)\t\t%.3e\n",BRBDstarlnu_full(3,1,obs,&param));
    printf("AFB(B->D* tau nu)\t\t%.3e\n",obs[1]);
    printf("Ptau(B->D* tau nu)\t\t%.3e\n",obs[2]);
    printf("PD*(B->D* tau nu)\t\t%.3e\n",obs[3]);
	printf("BR(B->D* tau nu)/BR(B->D* e nu)\t%.3e\n\n",BDstartaunu_BDstarenu(&param));
	printf("BR(Ds->tau nu)\t\t\t%.3e\n",Dstaunu(&param));
	printf("BR(Ds->mu nu)\t\t\t%.3e\n",Dsmunu(&param));
	printf("BR(D->mu nu)\t\t\t%.3e\n",Dmunu(&param));
	printf("BR(K->mu nu)/BR(pi->mu nu)\t%.3e\n",Kmunu_pimunu(&param));
	printf("Rmu23(K->mu nu)\t\t\t%.3e\n\n",Rmu23(&param));
	return 1;
}
