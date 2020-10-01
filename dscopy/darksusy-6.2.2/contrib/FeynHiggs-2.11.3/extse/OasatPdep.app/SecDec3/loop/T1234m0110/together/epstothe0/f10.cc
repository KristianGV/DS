#include "intfile.hh"

dcmplx Pf10(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
dcmplx y[72];
dcmplx FOUT;
dcmplx MYI(0.,1.);
y[1]=1./bi;
y[2]=em[0];
y[3]=x0*x0;
y[4]=esx[0];
y[5]=y[1]*y[2];
y[6]=2.*x2*y[1]*y[2];
y[7]=-x1;
y[8]=1.+y[7];
y[9]=2.*x0*y[1]*y[2];
y[10]=y[1]*y[2]*y[3];
y[11]=2.*x2*y[1]*y[2]*y[3];
y[12]=-(x0*y[1]*y[4]);
y[13]=-(y[1]*y[3]*y[4]);
y[14]=y[9]+y[10]+y[11]+y[12]+y[13];
y[15]=-x0;
y[16]=1.+y[15];
y[17]=x2*x2;
y[18]=lrs[0];
y[19]=2.*x1*x2*y[1]*y[2];
y[20]=-(x2*y[1]*y[4]);
y[21]=x1*y[1]*y[2];
y[22]=x2*y[1]*y[2];
y[23]=2.*x0*x1*x2*y[1]*y[2];
y[24]=y[1]*y[2]*y[17];
y[25]=2.*x0*x1*y[1]*y[2]*y[17];
y[26]=-(x1*x2*y[1]*y[4]);
y[27]=-2.*x0*x1*x2*y[1]*y[4];
y[28]=y[19]+y[20]+y[21]+y[22]+y[23]+y[24]+y[25]+y[26]+y[27];
y[29]=lrs[1];
y[30]=-x2;
y[31]=1.+y[30];
y[32]=2.*x1*y[1]*y[2];
y[33]=2.*x0*x1*y[1]*y[2];
y[34]=4.*x0*x1*x2*y[1]*y[2];
y[35]=-(y[1]*y[4]);
y[36]=-(x1*y[1]*y[4]);
y[37]=-2.*x0*x1*y[1]*y[4];
y[38]=y[5]+y[6]+y[32]+y[33]+y[34]+y[35]+y[36]+y[37];
y[39]=lambda*lambda;
y[40]=2.*x0*x2*y[1]*y[2];
y[41]=2.*x0*y[1]*y[2]*y[17];
y[42]=-2.*x0*x2*y[1]*y[4];
y[43]=y[5]+y[6]+y[20]+y[40]+y[41]+y[42];
y[44]=x0*y[1]*y[2];
y[45]=x2*y[1]*y[2]*y[3];
y[46]=y[1]*y[2]*y[3]*y[17];
y[47]=-(x0*x2*y[1]*y[4]);
y[48]=-(x2*y[1]*y[3]*y[4]);
y[49]=y[5]+y[40]+y[44]+y[45]+y[46]+y[47]+y[48];
y[50]=lrs[2];
y[51]=2.*x1*y[1]*y[2]*y[17];
y[52]=-2.*x1*x2*y[1]*y[4];
y[53]=y[19]+y[51]+y[52];
y[54]=-(lambda*MYI*x0*y[16]*y[18]*y[53]);
y[55]=-(lambda*MYI*y[16]*y[18]*y[28]);
y[56]=lambda*MYI*x0*y[18]*y[28];
y[57]=1.+y[54]+y[55]+y[56];
y[58]=-(lambda*MYI*y[8]*y[29]*y[49]);
y[59]=lambda*MYI*x1*y[29]*y[49];
y[60]=1.+y[58]+y[59];
y[61]=x1*y[1]*y[2]*y[3];
y[62]=2.*x1*x2*y[1]*y[2]*y[3];
y[63]=-(x0*x1*y[1]*y[4]);
y[64]=-(x1*y[1]*y[3]*y[4]);
y[65]=y[5]+y[12]+y[33]+y[40]+y[44]+y[61]+y[62]+y[63]+y[64];
y[66]=-(lambda*MYI*x1*y[8]*y[29]*y[49]);
y[67]=x1+y[66];
y[68]=-(lambda*MYI*x0*y[16]*y[18]*y[28]);
y[69]=x0+y[68];
y[70]=-(lambda*MYI*x2*y[31]*y[50]*y[65]);
y[71]=x2+y[70];
FOUT=pow(bi,-2)*pow(y[1]+y[1]*y[67]+y[1]*y[67]*y[69]+y[1]*y[71]+y[1]*y[67]*y\
[69]*y[71],-2)*(lambda*MYI*x2*y[14]*y[31]*y[50]*(x0*x1*y[8]*y[16]*y[18]*y[2\
9]*y[38]*y[39]*y[43]-lambda*MYI*x1*y[8]*y[14]*y[29]*y[57])-lambda*MYI*x2*y[\
31]*y[38]*y[50]*(-(x0*x1*y[8]*y[14]*y[16]*y[18]*y[29]*y[39]*y[43])+lambda*M\
YI*x0*y[16]*y[18]*y[38]*y[60])+(x0*x1*pow(y[43],2)*y[8]*y[16]*y[18]*y[29]*y\
[39]+y[57]*y[60])*(1.-lambda*MYI*x2*(2.*x1*y[1]*y[2]*y[3]+y[9])*y[31]*y[50]\
+lambda*MYI*x2*y[50]*y[65]-lambda*MYI*y[31]*y[50]*y[65]));
return (FOUT);
}