#include "intfile.hh"

dcmplx Pf8(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
dcmplx y[74];
dcmplx FOUT;
dcmplx MYI(0.,1.);
y[1]=1./bi;
y[2]=em[0];
y[3]=x0*x0;
y[4]=esx[0];
y[5]=2.*y[1]*y[2];
y[6]=2.*x1*y[1]*y[2];
y[7]=4.*x0*x1*y[1]*y[2];
y[8]=2.*x2*y[1]*y[2];
y[9]=4.*x0*x2*y[1]*y[2];
y[10]=-x1;
y[11]=1.+y[10];
y[12]=2.*x0*y[1]*y[2];
y[13]=2.*y[1]*y[2]*y[3];
y[14]=-(x0*y[1]*y[4]);
y[15]=-(y[1]*y[3]*y[4]);
y[16]=y[12]+y[13]+y[14]+y[15];
y[17]=-x0;
y[18]=1.+y[17];
y[19]=lrs[0];
y[20]=x1*x1;
y[21]=x2*x2;
y[22]=-(x1*y[1]*y[4]);
y[23]=y[1]*y[2]*y[20];
y[24]=2.*x0*y[1]*y[2]*y[20];
y[25]=2.*x1*x2*y[1]*y[2];
y[26]=4.*x0*x1*x2*y[1]*y[2];
y[27]=y[1]*y[2]*y[21];
y[28]=2.*x0*y[1]*y[2]*y[21];
y[29]=-(x1*x2*y[1]*y[4]);
y[30]=-2.*x0*x1*x2*y[1]*y[4];
y[31]=y[6]+y[8]+y[22]+y[23]+y[24]+y[25]+y[26]+y[27]+y[28]+y[29]+y[30];
y[32]=lrs[1];
y[33]=-x2;
y[34]=1.+y[33];
y[35]=-2.*x0*x1*y[1]*y[4];
y[36]=y[5]+y[6]+y[7]+y[8]+y[9]+y[22]+y[35];
y[37]=lambda*lambda;
y[38]=-(y[1]*y[4]);
y[39]=-(x2*y[1]*y[4]);
y[40]=-2.*x0*x2*y[1]*y[4];
y[41]=y[5]+y[6]+y[7]+y[8]+y[9]+y[38]+y[39]+y[40];
y[42]=y[1]*y[2];
y[43]=2.*x0*x1*y[1]*y[2];
y[44]=2.*x1*y[1]*y[2]*y[3];
y[45]=2.*x0*x2*y[1]*y[2];
y[46]=2.*x2*y[1]*y[2]*y[3];
y[47]=-(x0*x2*y[1]*y[4]);
y[48]=-(x2*y[1]*y[3]*y[4]);
y[49]=y[12]+y[14]+y[42]+y[43]+y[44]+y[45]+y[46]+y[47]+y[48];
y[50]=lrs[2];
y[51]=2.*y[1]*y[2]*y[20];
y[52]=4.*x1*x2*y[1]*y[2];
y[53]=2.*y[1]*y[2]*y[21];
y[54]=-2.*x1*x2*y[1]*y[4];
y[55]=y[51]+y[52]+y[53]+y[54];
y[56]=-(lambda*MYI*x0*y[18]*y[19]*y[55]);
y[57]=-(lambda*MYI*y[18]*y[19]*y[31]);
y[58]=lambda*MYI*x0*y[19]*y[31];
y[59]=1.+y[56]+y[57]+y[58];
y[60]=y[12]+y[13];
y[61]=-(lambda*MYI*x1*y[11]*y[32]*y[60]);
y[62]=-(lambda*MYI*y[11]*y[32]*y[49]);
y[63]=lambda*MYI*x1*y[32]*y[49];
y[64]=1.+y[61]+y[62]+y[63];
y[65]=-(x0*x1*y[1]*y[4]);
y[66]=-(x1*y[1]*y[3]*y[4]);
y[67]=y[12]+y[42]+y[43]+y[44]+y[45]+y[46]+y[65]+y[66];
y[68]=-(lambda*MYI*x1*y[11]*y[32]*y[49]);
y[69]=x1+y[68];
y[70]=-(lambda*MYI*x0*y[18]*y[19]*y[31]);
y[71]=x0+y[70];
y[72]=-(lambda*MYI*x2*y[34]*y[50]*y[67]);
y[73]=x2+y[72];
FOUT=pow(bi,-2)*pow(y[1]+y[1]*y[69]+y[1]*y[69]*y[71]+y[1]*y[73]+y[1]*y[71]*y\
[73],-2)*(lambda*MYI*x2*y[16]*y[34]*y[50]*(x0*x1*y[11]*y[18]*y[19]*y[32]*y[\
36]*y[37]*y[41]-lambda*MYI*x1*y[11]*y[16]*y[32]*y[59])-lambda*MYI*x2*y[34]*\
y[36]*y[50]*(-(x0*x1*y[11]*y[16]*y[18]*y[19]*y[32]*y[37]*y[41])+lambda*MYI*\
x0*y[18]*y[19]*y[36]*y[64])+(x0*x1*pow(y[41],2)*y[11]*y[18]*y[19]*y[32]*y[3\
7]+y[59]*y[64])*(1.-lambda*MYI*x2*y[34]*y[50]*y[60]+lambda*MYI*x2*y[50]*y[6\
7]-lambda*MYI*y[34]*y[50]*y[67]));
return (FOUT);
}
