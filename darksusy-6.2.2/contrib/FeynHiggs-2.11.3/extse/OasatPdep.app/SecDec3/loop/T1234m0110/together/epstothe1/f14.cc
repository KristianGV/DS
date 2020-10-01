#include "intfile.hh"

dcmplx Pf14(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
dcmplx y[97];
dcmplx FOUT;
dcmplx MYI(0.,1.);
y[1]=1./bi;
y[2]=em[0];
y[3]=x1*x1;
y[4]=esx[0];
y[5]=x0*x0;
y[6]=-x0;
y[7]=1.+y[6];
y[8]=2.*x0*x1*y[1]*y[2];
y[9]=y[1]*y[2]*y[3];
y[10]=-(x1*y[1]*y[4]);
y[11]=y[1]*y[2];
y[12]=2.*x1*y[1]*y[2];
y[13]=2.*x1*x2*y[1]*y[2];
y[14]=lrs[0];
y[15]=-x1;
y[16]=1.+y[15];
y[17]=2.*x0*y[1]*y[2];
y[18]=y[1]*y[2]*y[5];
y[19]=2.*x1*y[1]*y[2]*y[5];
y[20]=-(x0*y[1]*y[4]);
y[21]=y[8]+y[17]+y[18]+y[19]+y[20];
y[22]=x1*y[1]*y[2];
y[23]=x2*y[1]*y[2];
y[24]=2.*x0*x1*x2*y[1]*y[2];
y[25]=x2*y[1]*y[2]*y[3];
y[26]=2.*x0*x2*y[1]*y[2]*y[3];
y[27]=-(x1*x2*y[1]*y[4]);
y[28]=y[9]+y[10]+y[13]+y[22]+y[23]+y[24]+y[25]+y[26]+y[27];
y[29]=-(lambda*MYI*y[7]*y[14]*y[28]);
y[30]=lrs[1];
y[31]=-x2;
y[32]=1.+y[31];
y[33]=2.*x0*y[1]*y[2]*y[3];
y[34]=y[8]+y[9]+y[10]+y[11]+y[12]+y[33];
y[35]=lambda*lambda;
y[36]=2.*x2*y[1]*y[2];
y[37]=2.*x0*x2*y[1]*y[2];
y[38]=4.*x0*x1*x2*y[1]*y[2];
y[39]=-(y[1]*y[4]);
y[40]=-(x2*y[1]*y[4]);
y[41]=y[11]+y[12]+y[13]+y[36]+y[37]+y[38]+y[39]+y[40];
y[42]=x0*y[1]*y[2];
y[43]=x2*y[1]*y[2]*y[5];
y[44]=2.*x1*x2*y[1]*y[2]*y[5];
y[45]=-(x0*x2*y[1]*y[4]);
y[46]=y[8]+y[11]+y[20]+y[24]+y[37]+y[42]+y[43]+y[44]+y[45];
y[47]=lrs[2];
y[48]=2.*x2*y[1]*y[2]*y[3];
y[49]=y[13]+y[48];
y[50]=-(lambda*MYI*x0*y[7]*y[14]*y[49]);
y[51]=lambda*MYI*x0*y[14]*y[28];
y[52]=1.+y[29]+y[50]+y[51];
y[53]=2.*x2*y[1]*y[2]*y[5];
y[54]=y[17]+y[37]+y[53];
y[55]=-(lambda*MYI*x1*y[16]*y[30]*y[54]);
y[56]=-(lambda*MYI*y[16]*y[30]*y[46]);
y[57]=lambda*MYI*x1*y[30]*y[46];
y[58]=1.+y[55]+y[56]+y[57];
y[59]=x1*y[1]*y[2]*y[5];
y[60]=x0*y[1]*y[2]*y[3];
y[61]=y[1]*y[2]*y[3]*y[5];
y[62]=-(x0*x1*y[1]*y[4]);
y[63]=y[8]+y[42]+y[59]+y[60]+y[61]+y[62];
y[64]=-(lambda*MYI*x0*y[7]*y[14]*y[28]);
y[65]=x0+y[64];
y[66]=-(lambda*MYI*x1*y[16]*y[30]*y[46]);
y[67]=x1+y[66];
y[68]=-(lambda*MYI*x2*y[32]*y[47]*y[63]);
y[69]=x2+y[68];
y[70]=pow(bi,-2);
y[71]=1.+y[29];
y[72]=x0*x1*y[7]*y[14]*y[16]*y[30]*y[34]*y[35]*y[41];
y[73]=-(lambda*MYI*x1*y[16]*y[21]*y[30]*y[52]);
y[74]=y[72]+y[73];
y[75]=lambda*MYI*x2*y[21]*y[32]*y[47]*y[74];
y[76]=-(x0*x1*y[7]*y[14]*y[16]*y[21]*y[30]*y[35]*y[41]);
y[77]=lambda*MYI*x0*y[7]*y[14]*y[34]*y[58];
y[78]=y[76]+y[77];
y[79]=-(lambda*MYI*x2*y[32]*y[34]*y[47]*y[78]);
y[80]=pow(y[41],2);
y[81]=x0*x1*y[7]*y[14]*y[16]*y[30]*y[35]*y[80];
y[82]=y[52]*y[58];
y[83]=y[81]+y[82];
y[84]=-(lambda*MYI*y[32]*y[47]*y[63]);
y[85]=lambda*MYI*x2*y[47]*y[63];
y[86]=1.+y[84]+y[85];
y[87]=y[83]*y[86];
y[88]=y[75]+y[79]+y[87];
y[89]=y[1]*y[65]*y[67];
y[90]=y[1]*y[65]*y[69];
y[91]=y[1]*y[65]*y[67]*y[69];
y[92]=pow(y[65],2);
y[93]=y[1]*y[67]*y[69]*y[92];
y[94]=y[1]+y[89]+y[90]+y[91]+y[93];
y[95]=pow(y[94],-2);
y[96]=pow(y[67],2);
FOUT=(-2.*x0*myLog(x0)-x0*myLog(x2))*y[70]*y[71]*y[88]*y[95]+x0*(myLog(bi)*y\
[70]*y[71]*y[88]*y[95]-2.*myLog(y[71])*y[70]*y[71]*y[88]*y[95]-myLog(1.+y[8\
4])*y[70]*y[71]*y[88]*y[95]+3.*myLog(y[94])*y[70]*y[71]*y[88]*y[95]-2.*myLo\
g(y[11]+y[39]+y[1]*y[2]*y[67]+y[1]*y[2]*y[65]*y[67]-y[1]*y[4]*y[65]*y[67]+y\
[1]*y[2]*y[65]*y[69]+2.*y[1]*y[2]*y[65]*y[67]*y[69]-y[1]*y[4]*y[65]*y[67]*y\
[69]+y[1]*y[2]*y[67]*y[69]*y[92]+y[1]*y[2]*y[65]*y[96]+y[1]*y[2]*y[65]*y[69\
]*y[96]+y[1]*y[2]*y[69]*y[92]*y[96])*y[70]*y[71]*y[88]*y[95]);
return (FOUT);
}
