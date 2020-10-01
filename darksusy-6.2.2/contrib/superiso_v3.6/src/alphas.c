#include "include.h"


double alphas_running(double Q, double mtop, double mbot, struct parameters* param)
/* computes the QCD coupling constant alphas at the energy Q, fixing at mtop and mbot the matching scales between the flavours */
{
	double beta0,beta0_sq,beta1,beta1_sq,beta2,alphas_running,Lambda4,Lambda5,Lambda6,Lambda_min,Lambda_max,Lambda_moy,alphas_min,alphas_max,alphas_moy;
	int nf;

	double MZ=param->mass_Z;
	double alphas_MZ=param->alphas_MZ;

	nf=5;
	beta0 = 11.-2./3.*nf;
    beta0_sq = beta0*beta0;
	beta1=51.-19./3.*nf;
    beta1_sq = beta1*beta1;
	beta2=2857.-5033.*nf/9.+325./27.*nf*nf;

	if(param->Lambda5==-1.) return -1.;

	if(param->Lambda5==0.)
	{
		Lambda_min=1.e-3;
		Lambda_max=1.;
		alphas_min=0.;

		while((fabs(1.-alphas_min/alphas_MZ)>=1.e-4)&&(fabs(1.-Lambda_min/Lambda_max)>1.e-5))
		{
            double log_MZ_Lambda_min = log(pow(MZ/Lambda_min,2.));
            double log_MZ_Lambda_max = log(pow(MZ/Lambda_max,2.));
                        
			alphas_min=4.*pi/beta0/log_MZ_Lambda_min*(1.-2.*beta1/beta0_sq*log(log_MZ_Lambda_min)/log_MZ_Lambda_min+4.*beta1_sq/pow(beta0_sq*log_MZ_Lambda_min,2.)*(pow(log(log_MZ_Lambda_min)-1./2.,2.)+beta2*beta0/8./beta1_sq-5./4.));

			alphas_max=4.*pi/beta0/log_MZ_Lambda_max*(1.-2.*beta1/beta0_sq*log(log_MZ_Lambda_max)/log_MZ_Lambda_max+4.*beta1_sq/pow(beta0_sq*log_MZ_Lambda_max,2.)*(pow(log(log_MZ_Lambda_max)-1./2.,2.)+beta2*beta0/8./beta1_sq-5./4.));

			Lambda_moy=(Lambda_min+Lambda_max)/2.;
                        double log_MZ_Lambda_moy = log(pow(MZ/Lambda_moy,2.));

			alphas_moy=4.*pi/beta0/log_MZ_Lambda_moy*(1.-2.*beta1/beta0_sq*log(log_MZ_Lambda_moy)/log_MZ_Lambda_moy+4.*beta1_sq/pow(beta0_sq*log_MZ_Lambda_moy,2.)*(pow(log(log_MZ_Lambda_moy)-1./2.,2.)+beta2*beta0/8./beta1_sq-5./4.));

			if((alphas_MZ>=alphas_min)&&(alphas_MZ<=alphas_moy))
				Lambda_max=Lambda_moy;
			else Lambda_min=Lambda_moy;
		}

		Lambda5=Lambda_min;
		param->Lambda5=Lambda5;
		
		if(fabs(1.-Lambda_min/Lambda_max)<=1.e-5)
		{
			param->Lambda5=-1.;
			return -1.;
		}
	}
	else Lambda5=param->Lambda5;

	if((Q<=mtop)&&(Q>=mbot))
	/* 5 active flavours */
	{
		alphas_running=4.*pi/beta0/log(pow(Q/Lambda5,2.))*(1.-2.*beta1/beta0_sq*log(log(pow(Q/Lambda5,2.)))/log(pow(Q/Lambda5,2.))+4.*beta1_sq/pow(beta0_sq*log(pow(Q/Lambda5,2.)),2.)*(pow(log(log(pow(Q/Lambda5,2.)))-1./2.,2.)+beta2*beta0/8./beta1_sq-5./4.));
		return alphas_running;
	}
	else if((Q>mtop))
	/* 6 active flavours */
	{
		alphas_running=4.*pi/beta0/log(pow(mtop/Lambda5,2.))*(1.-2.*beta1/beta0_sq*log(log(pow(mtop/Lambda5,2.)))/log(pow(mtop/Lambda5,2.))+4.*beta1_sq/pow(beta0_sq*log(pow(mtop/Lambda5,2.)),2.)*(pow(log(log(pow(mtop/Lambda5,2.)))-1./2.,2.)+beta2*beta0/8./beta1_sq-5./4.));

		nf=6;
		beta0 = 11.-2./3.*nf;
                beta0_sq = beta0*beta0;
		beta1=51.-19./3.*nf;
                beta1_sq = beta1*beta1;
		beta2=2857.-5033.*nf/9.+325./27.*nf*nf;

		Lambda_min=1.e-3;
		Lambda_max=1.;
		alphas_min=0.;

		while((fabs(1.-alphas_min/alphas_running)>=1.e-4)&&(fabs(1.-Lambda_min/Lambda_max)>1.e-5))
		{
                        double log_mtop_Lambda_min = log(pow(mtop/Lambda_min,2.));
                        double log_mtop_Lambda_max = log(pow(mtop/Lambda_max,2.));
                        
			alphas_min=4.*pi/beta0/log_mtop_Lambda_min*(1.-2.*beta1/beta0_sq*log(log_mtop_Lambda_min)/log_mtop_Lambda_min+4.*beta1_sq/pow(beta0_sq*log_mtop_Lambda_min,2.)*(pow(log(log_mtop_Lambda_min)-1./2.,2.)+beta2*beta0/8./beta1_sq-5./4.));

			alphas_max=4.*pi/beta0/log_mtop_Lambda_max*(1.-2.*beta1/beta0_sq*log(log_mtop_Lambda_max)/log_mtop_Lambda_max+4.*beta1_sq/pow(beta0_sq*log_mtop_Lambda_max,2.)*(pow(log(log_mtop_Lambda_max)-1./2.,2.)+beta2*beta0/8./beta1_sq-5./4.));

			Lambda_moy=(Lambda_min+Lambda_max)/2.;
                        double log_mtop_Lambda_moy = log(pow(mtop/Lambda_moy,2.));
                        
			alphas_moy=4.*pi/beta0/log_mtop_Lambda_moy*(1.-2.*beta1/beta0_sq*log(log_mtop_Lambda_moy)/log_mtop_Lambda_moy+4.*beta1_sq/pow(beta0_sq*log_mtop_Lambda_moy,2.)*(pow(log(log_mtop_Lambda_moy)-1./2.,2.)+beta2*beta0/8./beta1_sq-5./4.));

			if((alphas_running>=alphas_min)&&(alphas_running<=alphas_moy))
				Lambda_max=Lambda_moy;
			else Lambda_min=Lambda_moy;
		}

		Lambda6=Lambda_min;
		
		if(fabs(1.-Lambda_min/Lambda_max)<=1.e-5)
		{
			param->Lambda5=-1.;
			return -1.;
		}
	alphas_running=4.*pi/beta0/log(pow(Q/Lambda6,2.))*(1.-2.*beta1/beta0_sq*log(log(pow(Q/Lambda6,2.)))/log(pow(Q/Lambda6,2.))+4.*beta1_sq/pow(beta0_sq*log(pow(Q/Lambda6,2.)),2.)*(pow(log(log(pow(Q/Lambda6,2.)))-1./2.,2.)+beta2*beta0/8./beta1_sq-5./4.));

		return alphas_running;
	}
	else
	/* 4 active flavours */
	{
		alphas_running=4.*pi/beta0/log(pow(mbot/Lambda5,2.))*(1.-2.*beta1/beta0_sq*log(log(pow(mbot/Lambda5,2.)))/log(pow(mbot/Lambda5,2.))+4.*beta1_sq/pow(beta0_sq*log(pow(mbot/Lambda5,2.)),2.)*(pow(log(log(pow(mbot/Lambda5,2.)))-1./2.,2.)+beta2*beta0/8./beta1_sq-5./4.));

		nf=4;
		beta0 = 11.-2./3.*nf;
                beta0_sq = beta0*beta0;
		beta1=51.-19./3.*nf;
                beta1_sq = beta1*beta1;
		beta2=2857.-5033.*nf/9.+325./27.*nf*nf;

		Lambda_min=1.e-3;
		Lambda_max=1.;
		alphas_min=0.;

		while((fabs(1.-alphas_min/alphas_running)>=1.e-4)&&(fabs(1.-Lambda_min/Lambda_max)>1.e-5))
		{
                        double log_mbot_Lambda_min = log(pow(mbot/Lambda_min,2.));
                        double log_mbot_Lambda_max = log(pow(mbot/Lambda_max,2.));
                        
			alphas_min=4.*pi/beta0/log_mbot_Lambda_min*(1.-2.*beta1/beta0_sq*log(log_mbot_Lambda_min)/log_mbot_Lambda_min+4.*beta1_sq/pow(beta0_sq*log_mbot_Lambda_min,2.)*(pow(log(log_mbot_Lambda_min)-1./2.,2.)+beta2*beta0/8./beta1_sq-5./4.));

			alphas_max=4.*pi/beta0/log_mbot_Lambda_max*(1.-2.*beta1/beta0_sq*log(log_mbot_Lambda_max)/log_mbot_Lambda_max+4.*beta1_sq/pow(beta0_sq*log_mbot_Lambda_max,2.)*(pow(log(log_mbot_Lambda_max)-1./2.,2.)+beta2*beta0/8./beta1_sq-5./4.));

			Lambda_moy=(Lambda_min+Lambda_max)/2.;
                        double log_mbot_Lambda_moy = log(pow(mbot/Lambda_moy,2.));
                        
			alphas_moy=4.*pi/beta0/log_mbot_Lambda_moy*(1.-2.*beta1/beta0_sq*log(log_mbot_Lambda_moy)/log_mbot_Lambda_moy+4.*beta1_sq/pow(beta0_sq*log_mbot_Lambda_moy,2.)*(pow(log(log_mbot_Lambda_moy)-1./2.,2.)+beta2*beta0/8./beta1_sq-5./4.));

			if((alphas_running>=alphas_min)&&(alphas_running<=alphas_moy))
				Lambda_max=Lambda_moy;
			else Lambda_min=Lambda_moy;
		}
		Lambda4=Lambda_min;
		
		if(fabs(1.-Lambda_min/Lambda_max)<=1.e-5)
		{
			param->Lambda5=-1.;
			return -1.;
		}
		
		if(Q<Lambda4+0.3) Q=Lambda4+0.3;
		
        double log_Q_lambda4 = log(pow(Q/Lambda4,2.));
	alphas_running=4.*pi/beta0/log_Q_lambda4*(1.-2.*beta1/beta0_sq*log(log_Q_lambda4)/log_Q_lambda4+4.*beta1_sq/pow(beta0_sq*log_Q_lambda4,2.)*(pow(log(log_Q_lambda4)-1./2.,2.)+beta2*beta0/8./beta1_sq-5./4.));

		return alphas_running;
	}
}
