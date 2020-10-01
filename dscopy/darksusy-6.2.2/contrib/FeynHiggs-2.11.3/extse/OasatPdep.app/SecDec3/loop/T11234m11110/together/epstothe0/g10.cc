#include "intfile.hh"

double Pr10(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
double y[45];
double FOUT;
y[1]=1./bi;
y[2]=em[0];
y[3]=x1*x1;
y[4]=esx[0];
y[5]=-x0;
y[6]=1.+y[5];
y[7]=y[1]*y[2];
y[8]=2.*x1*y[1]*y[2];
y[9]=y[1]*y[2]*y[3];
y[10]=2.*x2*y[1]*y[2];
y[11]=2.*x0*x2*y[1]*y[2];
y[12]=2.*x1*x2*y[1]*y[2];
y[13]=4.*x0*x1*x2*y[1]*y[2];
y[14]=2.*x0*x2*y[1]*y[2]*y[3];
y[15]=-(x1*y[1]*y[4]);
y[16]=-(x2*y[1]*y[4]);
y[17]=-2.*x0*x1*x2*y[1]*y[4];
y[18]=y[7]+y[8]+y[9]+y[10]+y[11]+y[12]+y[13]+y[14]+y[15]+y[16]+y[17];
y[19]=lrs[0];
y[20]=x0*x0;
y[21]=-x1;
y[22]=1.+y[21];
y[23]=2.*x0*y[1]*y[2];
y[24]=2.*x0*x1*y[1]*y[2];
y[25]=2.*x2*y[1]*y[2]*y[20];
y[26]=2.*x1*x2*y[1]*y[2]*y[20];
y[27]=-(x0*y[1]*y[4]);
y[28]=-(x2*y[1]*y[4]*y[20]);
y[29]=y[7]+y[11]+y[23]+y[24]+y[25]+y[26]+y[27]+y[28];
y[30]=lrs[1];
y[31]=-x2;
y[32]=1.+y[31];
y[33]=y[1]*y[2]*y[20];
y[34]=2.*x1*y[1]*y[2]*y[20];
y[35]=y[1]*y[2]*y[3]*y[20];
y[36]=-(x1*y[1]*y[4]*y[20]);
y[37]=y[7]+y[23]+y[24]+y[27]+y[33]+y[34]+y[35]+y[36];
y[38]=lrs[2];
y[39]=pow(y[6],2);
y[40]=pow(y[18],2);
y[41]=pow(y[19],2);
y[42]=pow(y[22],2);
y[43]=pow(y[29],2);
y[44]=pow(y[30],2);
FOUT=(2.*x0*x1*x2*y[1]*y[2]*y[6]*y[18]*y[19]*y[22]*y[29]*y[30]*y[32]*y[37]*y\
[38]+4.*x1*x2*y[1]*y[2]*y[6]*y[18]*y[19]*y[20]*y[22]*y[29]*y[30]*y[32]*y[37\
]*y[38]+4.*x2*y[1]*y[2]*y[3]*y[6]*y[18]*y[19]*y[20]*y[22]*y[29]*y[30]*y[32]\
*y[37]*y[38]-2.*x1*x2*y[1]*y[4]*y[6]*y[18]*y[19]*y[20]*y[22]*y[29]*y[30]*y[\
32]*y[37]*y[38]+2.*x1*x2*y[1]*y[2]*y[20]*y[22]*y[29]*y[30]*y[39]*y[40]*y[41\
]+2.*x2*y[1]*y[2]*y[3]*y[20]*y[22]*y[29]*y[30]*y[39]*y[40]*y[41]-x1*x2*y[1]\
*y[4]*y[20]*y[22]*y[29]*y[30]*y[39]*y[40]*y[41]+x2*y[1]*y[2]*y[20]*y[32]*y[\
37]*y[38]*y[39]*y[40]*y[41]+2.*x1*x2*y[1]*y[2]*y[20]*y[32]*y[37]*y[38]*y[39\
]*y[40]*y[41]+x2*y[1]*y[2]*y[3]*y[20]*y[32]*y[37]*y[38]*y[39]*y[40]*y[41]-x\
1*x2*y[1]*y[4]*y[20]*y[32]*y[37]*y[38]*y[39]*y[40]*y[41]+x0*y[1]*y[2]*y[3]*\
y[6]*y[18]*y[19]*y[42]*y[43]*y[44]+2.*x2*y[1]*y[2]*y[3]*y[6]*y[18]*y[19]*y[\
20]*y[42]*y[43]*y[44]+x2*y[1]*y[2]*y[3]*y[20]*y[32]*y[37]*y[38]*y[42]*y[43]\
*y[44])/(-(x0*y[1]*y[2]*y[6]*y[18]*y[19])-2.*x0*x1*y[1]*y[2]*y[6]*y[18]*y[1\
9]-2.*x0*x2*y[1]*y[2]*y[6]*y[18]*y[19]-2.*x0*x1*x2*y[1]*y[2]*y[6]*y[18]*y[1\
9]-x0*y[1]*y[2]*y[3]*y[6]*y[18]*y[19]+x0*x1*y[1]*y[4]*y[6]*y[18]*y[19]+x0*x\
2*y[1]*y[4]*y[6]*y[18]*y[19]-2.*x2*y[1]*y[2]*y[6]*y[18]*y[19]*y[20]-4.*x1*x\
2*y[1]*y[2]*y[6]*y[18]*y[19]*y[20]-2.*x2*y[1]*y[2]*y[3]*y[6]*y[18]*y[19]*y[\
20]+2.*x1*x2*y[1]*y[4]*y[6]*y[18]*y[19]*y[20]-x1*y[1]*y[2]*y[22]*y[29]*y[30\
]-2.*x0*x1*y[1]*y[2]*y[22]*y[29]*y[30]-2.*x0*x1*x2*y[1]*y[2]*y[22]*y[29]*y[\
30]-2.*x0*y[1]*y[2]*y[3]*y[22]*y[29]*y[30]+x0*x1*y[1]*y[4]*y[22]*y[29]*y[30\
]-2.*x1*x2*y[1]*y[2]*y[20]*y[22]*y[29]*y[30]-2.*x2*y[1]*y[2]*y[3]*y[20]*y[2\
2]*y[29]*y[30]+x1*x2*y[1]*y[4]*y[20]*y[22]*y[29]*y[30]-x2*y[1]*y[2]*y[32]*y\
[37]*y[38]-2.*x0*x2*y[1]*y[2]*y[32]*y[37]*y[38]-2.*x0*x1*x2*y[1]*y[2]*y[32]\
*y[37]*y[38]+x0*x2*y[1]*y[4]*y[32]*y[37]*y[38]-x2*y[1]*y[2]*y[20]*y[32]*y[3\
7]*y[38]-2.*x1*x2*y[1]*y[2]*y[20]*y[32]*y[37]*y[38]-x2*y[1]*y[2]*y[3]*y[20]\
*y[32]*y[37]*y[38]+x1*x2*y[1]*y[4]*y[20]*y[32]*y[37]*y[38]);
return (FOUT);
}
double Pm10(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
double y[45];
double FOUT;
y[1]=1./bi;
y[2]=em[0];
y[3]=x0*x0;
y[4]=x1*x1;
y[5]=esx[0];
y[6]=y[1]*y[2];
y[7]=2.*x0*y[1]*y[2];
y[8]=2.*x0*x1*y[1]*y[2];
y[9]=-(x0*y[1]*y[5]);
y[10]=2.*x0*x2*y[1]*y[2];
y[11]=-x0;
y[12]=1.+y[11];
y[13]=2.*x1*y[1]*y[2];
y[14]=y[1]*y[2]*y[4];
y[15]=2.*x2*y[1]*y[2];
y[16]=2.*x1*x2*y[1]*y[2];
y[17]=4.*x0*x1*x2*y[1]*y[2];
y[18]=2.*x0*x2*y[1]*y[2]*y[4];
y[19]=-(x1*y[1]*y[5]);
y[20]=-(x2*y[1]*y[5]);
y[21]=-2.*x0*x1*x2*y[1]*y[5];
y[22]=y[6]+y[10]+y[13]+y[14]+y[15]+y[16]+y[17]+y[18]+y[19]+y[20]+y[21];
y[23]=lrs[0];
y[24]=-x1;
y[25]=1.+y[24];
y[26]=2.*x2*y[1]*y[2]*y[3];
y[27]=2.*x1*x2*y[1]*y[2]*y[3];
y[28]=-(x2*y[1]*y[3]*y[5]);
y[29]=y[6]+y[7]+y[8]+y[9]+y[10]+y[26]+y[27]+y[28];
y[30]=lrs[1];
y[31]=-x2;
y[32]=1.+y[31];
y[33]=y[1]*y[2]*y[3];
y[34]=2.*x1*y[1]*y[2]*y[3];
y[35]=y[1]*y[2]*y[3]*y[4];
y[36]=-(x1*y[1]*y[3]*y[5]);
y[37]=y[6]+y[7]+y[8]+y[9]+y[33]+y[34]+y[35]+y[36];
y[38]=lrs[2];
y[39]=pow(y[12],2);
y[40]=pow(y[22],2);
y[41]=pow(y[23],2);
y[42]=pow(y[25],2);
y[43]=pow(y[29],2);
y[44]=pow(y[30],2);
FOUT=pow(lambda*(-(x0*y[1]*y[2]*y[12]*y[22]*y[23])-2.*x0*x1*y[1]*y[2]*y[12]*\
y[22]*y[23]-2.*x0*x2*y[1]*y[2]*y[12]*y[22]*y[23]-2.*x0*x1*x2*y[1]*y[2]*y[12\
]*y[22]*y[23]-2.*x2*y[1]*y[2]*y[3]*y[12]*y[22]*y[23]-4.*x1*x2*y[1]*y[2]*y[3\
]*y[12]*y[22]*y[23]-x0*y[1]*y[2]*y[4]*y[12]*y[22]*y[23]-2.*x2*y[1]*y[2]*y[3\
]*y[4]*y[12]*y[22]*y[23]+x0*x1*y[1]*y[5]*y[12]*y[22]*y[23]+x0*x2*y[1]*y[5]*\
y[12]*y[22]*y[23]+2.*x1*x2*y[1]*y[3]*y[5]*y[12]*y[22]*y[23]-x1*y[1]*y[2]*y[\
25]*y[29]*y[30]-2.*x0*x1*y[1]*y[2]*y[25]*y[29]*y[30]-2.*x0*x1*x2*y[1]*y[2]*\
y[25]*y[29]*y[30]-2.*x1*x2*y[1]*y[2]*y[3]*y[25]*y[29]*y[30]-2.*x0*y[1]*y[2]\
*y[4]*y[25]*y[29]*y[30]-2.*x2*y[1]*y[2]*y[3]*y[4]*y[25]*y[29]*y[30]+x0*x1*y\
[1]*y[5]*y[25]*y[29]*y[30]+x1*x2*y[1]*y[3]*y[5]*y[25]*y[29]*y[30]-x2*y[1]*y\
[2]*y[32]*y[37]*y[38]-2.*x0*x2*y[1]*y[2]*y[32]*y[37]*y[38]-2.*x0*x1*x2*y[1]\
*y[2]*y[32]*y[37]*y[38]-x2*y[1]*y[2]*y[3]*y[32]*y[37]*y[38]-2.*x1*x2*y[1]*y\
[2]*y[3]*y[32]*y[37]*y[38]-x2*y[1]*y[2]*y[3]*y[4]*y[32]*y[37]*y[38]+x0*x2*y\
[1]*y[5]*y[32]*y[37]*y[38]+x1*x2*y[1]*y[3]*y[5]*y[32]*y[37]*y[38])-x2*pow(l\
ambda,5)*y[1]*y[2]*y[3]*y[4]*y[32]*y[37]*y[38]*y[39]*y[40]*y[41]*y[42]*y[43\
]*y[44]+pow(lambda,3)*(2.*x0*x1*x2*y[1]*y[2]*y[12]*y[22]*y[23]*y[25]*y[29]*\
y[30]*y[32]*y[37]*y[38]+4.*x1*x2*y[1]*y[2]*y[3]*y[12]*y[22]*y[23]*y[25]*y[2\
9]*y[30]*y[32]*y[37]*y[38]+4.*x2*y[1]*y[2]*y[3]*y[4]*y[12]*y[22]*y[23]*y[25\
]*y[29]*y[30]*y[32]*y[37]*y[38]-2.*x1*x2*y[1]*y[3]*y[5]*y[12]*y[22]*y[23]*y\
[25]*y[29]*y[30]*y[32]*y[37]*y[38]+2.*x1*x2*y[1]*y[2]*y[3]*y[25]*y[29]*y[30\
]*y[39]*y[40]*y[41]+2.*x2*y[1]*y[2]*y[3]*y[4]*y[25]*y[29]*y[30]*y[39]*y[40]\
*y[41]-x1*x2*y[1]*y[3]*y[5]*y[25]*y[29]*y[30]*y[39]*y[40]*y[41]+x2*y[1]*y[2\
]*y[3]*y[32]*y[37]*y[38]*y[39]*y[40]*y[41]+2.*x1*x2*y[1]*y[2]*y[3]*y[32]*y[\
37]*y[38]*y[39]*y[40]*y[41]+x2*y[1]*y[2]*y[3]*y[4]*y[32]*y[37]*y[38]*y[39]*\
y[40]*y[41]-x1*x2*y[1]*y[3]*y[5]*y[32]*y[37]*y[38]*y[39]*y[40]*y[41]+x0*y[1\
]*y[2]*y[4]*y[12]*y[22]*y[23]*y[42]*y[43]*y[44]+2.*x2*y[1]*y[2]*y[3]*y[4]*y\
[12]*y[22]*y[23]*y[42]*y[43]*y[44]+x2*y[1]*y[2]*y[3]*y[4]*y[32]*y[37]*y[38]\
*y[42]*y[43]*y[44]),2)+pow(x0*y[1]*y[2]+x1*y[1]*y[2]+x2*y[1]*y[2]+2.*x0*x1*\
x2*y[1]*y[2]+x2*y[1]*y[2]*y[3]+x0*y[1]*y[2]*y[4]+x2*y[1]*y[2]*y[3]*y[4]-x0*\
x1*y[1]*y[5]-x0*x2*y[1]*y[5]-x1*x2*y[1]*y[3]*y[5]+y[6]+y[8]+y[10]+y[27]+lam\
bda*lambda*(-2.*x0*x1*y[1]*y[2]*y[12]*y[22]*y[23]*y[25]*y[29]*y[30]-2.*x0*x\
1*x2*y[1]*y[2]*y[12]*y[22]*y[23]*y[25]*y[29]*y[30]-4.*x1*x2*y[1]*y[2]*y[3]*\
y[12]*y[22]*y[23]*y[25]*y[29]*y[30]-2.*x0*y[1]*y[2]*y[4]*y[12]*y[22]*y[23]*\
y[25]*y[29]*y[30]-4.*x2*y[1]*y[2]*y[3]*y[4]*y[12]*y[22]*y[23]*y[25]*y[29]*y\
[30]+x0*x1*y[1]*y[5]*y[12]*y[22]*y[23]*y[25]*y[29]*y[30]+2.*x1*x2*y[1]*y[3]\
*y[5]*y[12]*y[22]*y[23]*y[25]*y[29]*y[30]-2.*x0*x2*y[1]*y[2]*y[12]*y[22]*y[\
23]*y[32]*y[37]*y[38]-2.*x0*x1*x2*y[1]*y[2]*y[12]*y[22]*y[23]*y[32]*y[37]*y\
[38]-2.*x2*y[1]*y[2]*y[3]*y[12]*y[22]*y[23]*y[32]*y[37]*y[38]-4.*x1*x2*y[1]\
*y[2]*y[3]*y[12]*y[22]*y[23]*y[32]*y[37]*y[38]-2.*x2*y[1]*y[2]*y[3]*y[4]*y[\
12]*y[22]*y[23]*y[32]*y[37]*y[38]+x0*x2*y[1]*y[5]*y[12]*y[22]*y[23]*y[32]*y\
[37]*y[38]+2.*x1*x2*y[1]*y[3]*y[5]*y[12]*y[22]*y[23]*y[32]*y[37]*y[38]-2.*x\
0*x1*x2*y[1]*y[2]*y[25]*y[29]*y[30]*y[32]*y[37]*y[38]-2.*x1*x2*y[1]*y[2]*y[\
3]*y[25]*y[29]*y[30]*y[32]*y[37]*y[38]-2.*x2*y[1]*y[2]*y[3]*y[4]*y[25]*y[29\
]*y[30]*y[32]*y[37]*y[38]+x1*x2*y[1]*y[3]*y[5]*y[25]*y[29]*y[30]*y[32]*y[37\
]*y[38]-x2*y[1]*y[2]*y[3]*y[39]*y[40]*y[41]-2.*x1*x2*y[1]*y[2]*y[3]*y[39]*y\
[40]*y[41]-x2*y[1]*y[2]*y[3]*y[4]*y[39]*y[40]*y[41]+x1*x2*y[1]*y[3]*y[5]*y[\
39]*y[40]*y[41]-x0*y[1]*y[2]*y[4]*y[42]*y[43]*y[44]-x2*y[1]*y[2]*y[3]*y[4]*\
y[42]*y[43]*y[44])+pow(lambda,4)*(2.*x1*x2*y[1]*y[2]*y[3]*y[25]*y[29]*y[30]\
*y[32]*y[37]*y[38]*y[39]*y[40]*y[41]+2.*x2*y[1]*y[2]*y[3]*y[4]*y[25]*y[29]*\
y[30]*y[32]*y[37]*y[38]*y[39]*y[40]*y[41]-x1*x2*y[1]*y[3]*y[5]*y[25]*y[29]*\
y[30]*y[32]*y[37]*y[38]*y[39]*y[40]*y[41]+2.*x2*y[1]*y[2]*y[3]*y[4]*y[12]*y\
[22]*y[23]*y[32]*y[37]*y[38]*y[42]*y[43]*y[44]+x2*y[1]*y[2]*y[3]*y[4]*y[39]\
*y[40]*y[41]*y[42]*y[43]*y[44]),2);
return (FOUT);
}
double Ps10(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
double y[45];
double FOUT;
y[1]=1./bi;
y[2]=em[0];
y[3]=x0*x0;
y[4]=x1*x1;
y[5]=esx[0];
y[6]=y[1]*y[2];
y[7]=2.*x0*y[1]*y[2];
y[8]=2.*x0*x1*y[1]*y[2];
y[9]=-(x0*y[1]*y[5]);
y[10]=2.*x0*x2*y[1]*y[2];
y[11]=-x0;
y[12]=1.+y[11];
y[13]=2.*x1*y[1]*y[2];
y[14]=y[1]*y[2]*y[4];
y[15]=2.*x2*y[1]*y[2];
y[16]=2.*x1*x2*y[1]*y[2];
y[17]=4.*x0*x1*x2*y[1]*y[2];
y[18]=2.*x0*x2*y[1]*y[2]*y[4];
y[19]=-(x1*y[1]*y[5]);
y[20]=-(x2*y[1]*y[5]);
y[21]=-2.*x0*x1*x2*y[1]*y[5];
y[22]=y[6]+y[10]+y[13]+y[14]+y[15]+y[16]+y[17]+y[18]+y[19]+y[20]+y[21];
y[23]=lrs[0];
y[24]=-x1;
y[25]=1.+y[24];
y[26]=2.*x2*y[1]*y[2]*y[3];
y[27]=2.*x1*x2*y[1]*y[2]*y[3];
y[28]=-(x2*y[1]*y[3]*y[5]);
y[29]=y[6]+y[7]+y[8]+y[9]+y[10]+y[26]+y[27]+y[28];
y[30]=lrs[1];
y[31]=-x2;
y[32]=1.+y[31];
y[33]=y[1]*y[2]*y[3];
y[34]=2.*x1*y[1]*y[2]*y[3];
y[35]=y[1]*y[2]*y[3]*y[4];
y[36]=-(x1*y[1]*y[3]*y[5]);
y[37]=y[6]+y[7]+y[8]+y[9]+y[33]+y[34]+y[35]+y[36];
y[38]=lrs[2];
y[39]=pow(y[12],2);
y[40]=pow(y[22],2);
y[41]=pow(y[23],2);
y[42]=pow(y[25],2);
y[43]=pow(y[29],2);
y[44]=pow(y[30],2);
FOUT=lambda*(-(x0*y[1]*y[2]*y[12]*y[22]*y[23])-2.*x0*x1*y[1]*y[2]*y[12]*y[22\
]*y[23]-2.*x0*x2*y[1]*y[2]*y[12]*y[22]*y[23]-2.*x0*x1*x2*y[1]*y[2]*y[12]*y[\
22]*y[23]-2.*x2*y[1]*y[2]*y[3]*y[12]*y[22]*y[23]-4.*x1*x2*y[1]*y[2]*y[3]*y[\
12]*y[22]*y[23]-x0*y[1]*y[2]*y[4]*y[12]*y[22]*y[23]-2.*x2*y[1]*y[2]*y[3]*y[\
4]*y[12]*y[22]*y[23]+x0*x1*y[1]*y[5]*y[12]*y[22]*y[23]+x0*x2*y[1]*y[5]*y[12\
]*y[22]*y[23]+2.*x1*x2*y[1]*y[3]*y[5]*y[12]*y[22]*y[23]-x1*y[1]*y[2]*y[25]*\
y[29]*y[30]-2.*x0*x1*y[1]*y[2]*y[25]*y[29]*y[30]-2.*x0*x1*x2*y[1]*y[2]*y[25\
]*y[29]*y[30]-2.*x1*x2*y[1]*y[2]*y[3]*y[25]*y[29]*y[30]-2.*x0*y[1]*y[2]*y[4\
]*y[25]*y[29]*y[30]-2.*x2*y[1]*y[2]*y[3]*y[4]*y[25]*y[29]*y[30]+x0*x1*y[1]*\
y[5]*y[25]*y[29]*y[30]+x1*x2*y[1]*y[3]*y[5]*y[25]*y[29]*y[30]-x2*y[1]*y[2]*\
y[32]*y[37]*y[38]-2.*x0*x2*y[1]*y[2]*y[32]*y[37]*y[38]-2.*x0*x1*x2*y[1]*y[2\
]*y[32]*y[37]*y[38]-x2*y[1]*y[2]*y[3]*y[32]*y[37]*y[38]-2.*x1*x2*y[1]*y[2]*\
y[3]*y[32]*y[37]*y[38]-x2*y[1]*y[2]*y[3]*y[4]*y[32]*y[37]*y[38]+x0*x2*y[1]*\
y[5]*y[32]*y[37]*y[38]+x1*x2*y[1]*y[3]*y[5]*y[32]*y[37]*y[38])-x2*pow(lambd\
a,5)*y[1]*y[2]*y[3]*y[4]*y[32]*y[37]*y[38]*y[39]*y[40]*y[41]*y[42]*y[43]*y[\
44]+pow(lambda,3)*(2.*x0*x1*x2*y[1]*y[2]*y[12]*y[22]*y[23]*y[25]*y[29]*y[30\
]*y[32]*y[37]*y[38]+4.*x1*x2*y[1]*y[2]*y[3]*y[12]*y[22]*y[23]*y[25]*y[29]*y\
[30]*y[32]*y[37]*y[38]+4.*x2*y[1]*y[2]*y[3]*y[4]*y[12]*y[22]*y[23]*y[25]*y[\
29]*y[30]*y[32]*y[37]*y[38]-2.*x1*x2*y[1]*y[3]*y[5]*y[12]*y[22]*y[23]*y[25]\
*y[29]*y[30]*y[32]*y[37]*y[38]+2.*x1*x2*y[1]*y[2]*y[3]*y[25]*y[29]*y[30]*y[\
39]*y[40]*y[41]+2.*x2*y[1]*y[2]*y[3]*y[4]*y[25]*y[29]*y[30]*y[39]*y[40]*y[4\
1]-x1*x2*y[1]*y[3]*y[5]*y[25]*y[29]*y[30]*y[39]*y[40]*y[41]+x2*y[1]*y[2]*y[\
3]*y[32]*y[37]*y[38]*y[39]*y[40]*y[41]+2.*x1*x2*y[1]*y[2]*y[3]*y[32]*y[37]*\
y[38]*y[39]*y[40]*y[41]+x2*y[1]*y[2]*y[3]*y[4]*y[32]*y[37]*y[38]*y[39]*y[40\
]*y[41]-x1*x2*y[1]*y[3]*y[5]*y[32]*y[37]*y[38]*y[39]*y[40]*y[41]+x0*y[1]*y[\
2]*y[4]*y[12]*y[22]*y[23]*y[42]*y[43]*y[44]+2.*x2*y[1]*y[2]*y[3]*y[4]*y[12]\
*y[22]*y[23]*y[42]*y[43]*y[44]+x2*y[1]*y[2]*y[3]*y[4]*y[32]*y[37]*y[38]*y[4\
2]*y[43]*y[44]);
return (FOUT);
}
double Pa10(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
double y[45];
double FOUT;
y[1]=1./bi;
y[2]=em[0];
y[3]=x0*x0;
y[4]=x1*x1;
y[5]=esx[0];
y[6]=y[1]*y[2];
y[7]=2.*x0*y[1]*y[2];
y[8]=2.*x0*x1*y[1]*y[2];
y[9]=-(x0*y[1]*y[5]);
y[10]=2.*x0*x2*y[1]*y[2];
y[11]=-x0;
y[12]=1.+y[11];
y[13]=2.*x1*y[1]*y[2];
y[14]=y[1]*y[2]*y[4];
y[15]=2.*x2*y[1]*y[2];
y[16]=2.*x1*x2*y[1]*y[2];
y[17]=4.*x0*x1*x2*y[1]*y[2];
y[18]=2.*x0*x2*y[1]*y[2]*y[4];
y[19]=-(x1*y[1]*y[5]);
y[20]=-(x2*y[1]*y[5]);
y[21]=-2.*x0*x1*x2*y[1]*y[5];
y[22]=y[6]+y[10]+y[13]+y[14]+y[15]+y[16]+y[17]+y[18]+y[19]+y[20]+y[21];
y[23]=lrs[0];
y[24]=-x1;
y[25]=1.+y[24];
y[26]=2.*x2*y[1]*y[2]*y[3];
y[27]=2.*x1*x2*y[1]*y[2]*y[3];
y[28]=-(x2*y[1]*y[3]*y[5]);
y[29]=y[6]+y[7]+y[8]+y[9]+y[10]+y[26]+y[27]+y[28];
y[30]=lrs[1];
y[31]=-x2;
y[32]=1.+y[31];
y[33]=y[1]*y[2]*y[3];
y[34]=2.*x1*y[1]*y[2]*y[3];
y[35]=y[1]*y[2]*y[3]*y[4];
y[36]=-(x1*y[1]*y[3]*y[5]);
y[37]=y[6]+y[7]+y[8]+y[9]+y[33]+y[34]+y[35]+y[36];
y[38]=lrs[2];
y[39]=pow(y[12],2);
y[40]=pow(y[22],2);
y[41]=pow(y[23],2);
y[42]=pow(y[25],2);
y[43]=pow(y[29],2);
y[44]=pow(y[30],2);
FOUT=(lambda*(-(x0*y[1]*y[2]*y[12]*y[22]*y[23])-2.*x0*x1*y[1]*y[2]*y[12]*y[2\
2]*y[23]-2.*x0*x2*y[1]*y[2]*y[12]*y[22]*y[23]-2.*x0*x1*x2*y[1]*y[2]*y[12]*y\
[22]*y[23]-2.*x2*y[1]*y[2]*y[3]*y[12]*y[22]*y[23]-4.*x1*x2*y[1]*y[2]*y[3]*y\
[12]*y[22]*y[23]-x0*y[1]*y[2]*y[4]*y[12]*y[22]*y[23]-2.*x2*y[1]*y[2]*y[3]*y\
[4]*y[12]*y[22]*y[23]+x0*x1*y[1]*y[5]*y[12]*y[22]*y[23]+x0*x2*y[1]*y[5]*y[1\
2]*y[22]*y[23]+2.*x1*x2*y[1]*y[3]*y[5]*y[12]*y[22]*y[23]-x1*y[1]*y[2]*y[25]\
*y[29]*y[30]-2.*x0*x1*y[1]*y[2]*y[25]*y[29]*y[30]-2.*x0*x1*x2*y[1]*y[2]*y[2\
5]*y[29]*y[30]-2.*x1*x2*y[1]*y[2]*y[3]*y[25]*y[29]*y[30]-2.*x0*y[1]*y[2]*y[\
4]*y[25]*y[29]*y[30]-2.*x2*y[1]*y[2]*y[3]*y[4]*y[25]*y[29]*y[30]+x0*x1*y[1]\
*y[5]*y[25]*y[29]*y[30]+x1*x2*y[1]*y[3]*y[5]*y[25]*y[29]*y[30]-x2*y[1]*y[2]\
*y[32]*y[37]*y[38]-2.*x0*x2*y[1]*y[2]*y[32]*y[37]*y[38]-2.*x0*x1*x2*y[1]*y[\
2]*y[32]*y[37]*y[38]-x2*y[1]*y[2]*y[3]*y[32]*y[37]*y[38]-2.*x1*x2*y[1]*y[2]\
*y[3]*y[32]*y[37]*y[38]-x2*y[1]*y[2]*y[3]*y[4]*y[32]*y[37]*y[38]+x0*x2*y[1]\
*y[5]*y[32]*y[37]*y[38]+x1*x2*y[1]*y[3]*y[5]*y[32]*y[37]*y[38])-x2*pow(lamb\
da,5)*y[1]*y[2]*y[3]*y[4]*y[32]*y[37]*y[38]*y[39]*y[40]*y[41]*y[42]*y[43]*y\
[44]+pow(lambda,3)*(2.*x0*x1*x2*y[1]*y[2]*y[12]*y[22]*y[23]*y[25]*y[29]*y[3\
0]*y[32]*y[37]*y[38]+4.*x1*x2*y[1]*y[2]*y[3]*y[12]*y[22]*y[23]*y[25]*y[29]*\
y[30]*y[32]*y[37]*y[38]+4.*x2*y[1]*y[2]*y[3]*y[4]*y[12]*y[22]*y[23]*y[25]*y\
[29]*y[30]*y[32]*y[37]*y[38]-2.*x1*x2*y[1]*y[3]*y[5]*y[12]*y[22]*y[23]*y[25\
]*y[29]*y[30]*y[32]*y[37]*y[38]+2.*x1*x2*y[1]*y[2]*y[3]*y[25]*y[29]*y[30]*y\
[39]*y[40]*y[41]+2.*x2*y[1]*y[2]*y[3]*y[4]*y[25]*y[29]*y[30]*y[39]*y[40]*y[\
41]-x1*x2*y[1]*y[3]*y[5]*y[25]*y[29]*y[30]*y[39]*y[40]*y[41]+x2*y[1]*y[2]*y\
[3]*y[32]*y[37]*y[38]*y[39]*y[40]*y[41]+2.*x1*x2*y[1]*y[2]*y[3]*y[32]*y[37]\
*y[38]*y[39]*y[40]*y[41]+x2*y[1]*y[2]*y[3]*y[4]*y[32]*y[37]*y[38]*y[39]*y[4\
0]*y[41]-x1*x2*y[1]*y[3]*y[5]*y[32]*y[37]*y[38]*y[39]*y[40]*y[41]+x0*y[1]*y\
[2]*y[4]*y[12]*y[22]*y[23]*y[42]*y[43]*y[44]+2.*x2*y[1]*y[2]*y[3]*y[4]*y[12\
]*y[22]*y[23]*y[42]*y[43]*y[44]+x2*y[1]*y[2]*y[3]*y[4]*y[32]*y[37]*y[38]*y[\
42]*y[43]*y[44]))/(lambda*(x0*y[1]*y[2]+x1*y[1]*y[2]+x2*y[1]*y[2]+2.*x0*x1*\
x2*y[1]*y[2]+x2*y[1]*y[2]*y[3]+x0*y[1]*y[2]*y[4]+x2*y[1]*y[2]*y[3]*y[4]-x0*\
x1*y[1]*y[5]-x0*x2*y[1]*y[5]-x1*x2*y[1]*y[3]*y[5]+y[6]+y[8]+y[10]+y[27]+lam\
bda*lambda*(-2.*x0*x1*y[1]*y[2]*y[12]*y[22]*y[23]*y[25]*y[29]*y[30]-2.*x0*x\
1*x2*y[1]*y[2]*y[12]*y[22]*y[23]*y[25]*y[29]*y[30]-4.*x1*x2*y[1]*y[2]*y[3]*\
y[12]*y[22]*y[23]*y[25]*y[29]*y[30]-2.*x0*y[1]*y[2]*y[4]*y[12]*y[22]*y[23]*\
y[25]*y[29]*y[30]-4.*x2*y[1]*y[2]*y[3]*y[4]*y[12]*y[22]*y[23]*y[25]*y[29]*y\
[30]+x0*x1*y[1]*y[5]*y[12]*y[22]*y[23]*y[25]*y[29]*y[30]+2.*x1*x2*y[1]*y[3]\
*y[5]*y[12]*y[22]*y[23]*y[25]*y[29]*y[30]-2.*x0*x2*y[1]*y[2]*y[12]*y[22]*y[\
23]*y[32]*y[37]*y[38]-2.*x0*x1*x2*y[1]*y[2]*y[12]*y[22]*y[23]*y[32]*y[37]*y\
[38]-2.*x2*y[1]*y[2]*y[3]*y[12]*y[22]*y[23]*y[32]*y[37]*y[38]-4.*x1*x2*y[1]\
*y[2]*y[3]*y[12]*y[22]*y[23]*y[32]*y[37]*y[38]-2.*x2*y[1]*y[2]*y[3]*y[4]*y[\
12]*y[22]*y[23]*y[32]*y[37]*y[38]+x0*x2*y[1]*y[5]*y[12]*y[22]*y[23]*y[32]*y\
[37]*y[38]+2.*x1*x2*y[1]*y[3]*y[5]*y[12]*y[22]*y[23]*y[32]*y[37]*y[38]-2.*x\
0*x1*x2*y[1]*y[2]*y[25]*y[29]*y[30]*y[32]*y[37]*y[38]-2.*x1*x2*y[1]*y[2]*y[\
3]*y[25]*y[29]*y[30]*y[32]*y[37]*y[38]-2.*x2*y[1]*y[2]*y[3]*y[4]*y[25]*y[29\
]*y[30]*y[32]*y[37]*y[38]+x1*x2*y[1]*y[3]*y[5]*y[25]*y[29]*y[30]*y[32]*y[37\
]*y[38]-x2*y[1]*y[2]*y[3]*y[39]*y[40]*y[41]-2.*x1*x2*y[1]*y[2]*y[3]*y[39]*y\
[40]*y[41]-x2*y[1]*y[2]*y[3]*y[4]*y[39]*y[40]*y[41]+x1*x2*y[1]*y[3]*y[5]*y[\
39]*y[40]*y[41]-x0*y[1]*y[2]*y[4]*y[42]*y[43]*y[44]-x2*y[1]*y[2]*y[3]*y[4]*\
y[42]*y[43]*y[44])+pow(lambda,4)*(2.*x1*x2*y[1]*y[2]*y[3]*y[25]*y[29]*y[30]\
*y[32]*y[37]*y[38]*y[39]*y[40]*y[41]+2.*x2*y[1]*y[2]*y[3]*y[4]*y[25]*y[29]*\
y[30]*y[32]*y[37]*y[38]*y[39]*y[40]*y[41]-x1*x2*y[1]*y[3]*y[5]*y[25]*y[29]*\
y[30]*y[32]*y[37]*y[38]*y[39]*y[40]*y[41]+2.*x2*y[1]*y[2]*y[3]*y[4]*y[12]*y\
[22]*y[23]*y[32]*y[37]*y[38]*y[42]*y[43]*y[44]+x2*y[1]*y[2]*y[3]*y[4]*y[39]\
*y[40]*y[41]*y[42]*y[43]*y[44])));
return (FOUT);
}
double Pt10t1(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
double y[5];
double FOUT;
y[1]=1./bi;
y[2]=em[0];
y[3]=x1*x1;
y[4]=esx[0];
FOUT=(1.-x0)*x0*(y[1]*y[2]+2.*x1*y[1]*y[2]+2.*x2*y[1]*y[2]+2.*x0*x2*y[1]*y[2\
]+2.*x1*x2*y[1]*y[2]+4.*x0*x1*x2*y[1]*y[2]+y[1]*y[2]*y[3]+2.*x0*x2*y[1]*y[2\
]*y[3]-x1*y[1]*y[4]-x2*y[1]*y[4]-2.*x0*x1*x2*y[1]*y[4]);
return (FOUT);
}
double Pt10t2(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
double y[5];
double FOUT;
y[1]=1./bi;
y[2]=em[0];
y[3]=x0*x0;
y[4]=esx[0];
FOUT=(1.-x1)*x1*(y[1]*y[2]+2.*x0*y[1]*y[2]+2.*x0*x1*y[1]*y[2]+2.*x0*x2*y[1]*\
y[2]+2.*x2*y[1]*y[2]*y[3]+2.*x1*x2*y[1]*y[2]*y[3]-x0*y[1]*y[4]-x2*y[1]*y[3]\
*y[4]);
return (FOUT);
}
double Pt10t3(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
double y[5];
double FOUT;
y[1]=1./bi;
y[2]=em[0];
y[3]=x0*x0;
y[4]=esx[0];
FOUT=(1.-x2)*x2*(y[1]*y[2]+2.*x0*y[1]*y[2]+2.*x0*x1*y[1]*y[2]+y[1]*y[2]*y[3]\
+2.*x1*y[1]*y[2]*y[3]+x1*x1*y[1]*y[2]*y[3]-x0*y[1]*y[4]-x1*y[1]*y[3]*y[4]);
return (FOUT);
}