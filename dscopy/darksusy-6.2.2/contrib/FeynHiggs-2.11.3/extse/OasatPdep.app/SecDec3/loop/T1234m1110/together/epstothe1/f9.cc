#include "intfile.hh"

dcmplx Pf9(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
dcmplx y[97];
dcmplx FOUT;
dcmplx MYI(0.,1.);
y[1]=1./bi;
y[2]=em[0];
y[3]=x0*x0;
y[4]=esx[0];
y[5]=2.*y[1]*y[2];
y[6]=2.*x0*y[1]*y[2];
y[7]=2.*x2*y[1]*y[2];
y[8]=-(y[1]*y[4]);
y[9]=-x1;
y[10]=1.+y[9];
y[11]=2.*y[1]*y[2]*y[3];
y[12]=2.*x2*y[1]*y[2]*y[3];
y[13]=-(y[1]*y[3]*y[4]);
y[14]=y[6]+y[11]+y[12]+y[13];
y[15]=-x0;
y[16]=1.+y[15];
y[17]=2.*x1*y[1]*y[2];
y[18]=x2*x2;
y[19]=lrs[0];
y[20]=4.*x0*x1*x2*y[1]*y[2];
y[21]=y[1]*y[2];
y[22]=2.*x0*x1*y[1]*y[2];
y[23]=2.*x1*x2*y[1]*y[2];
y[24]=y[1]*y[2]*y[18];
y[25]=2.*x0*x1*y[1]*y[2]*y[18];
y[26]=-(x1*y[1]*y[4]);
y[27]=-(x2*y[1]*y[4]);
y[28]=-2.*x0*x1*x2*y[1]*y[4];
y[29]=y[7]+y[17]+y[20]+y[21]+y[22]+y[23]+y[24]+y[25]+y[26]+y[27]+y[28];
y[30]=lrs[1];
y[31]=-x2;
y[32]=1.+y[31];
y[33]=4.*x0*x1*y[1]*y[2];
y[34]=-2.*x0*x1*y[1]*y[4];
y[35]=y[5]+y[7]+y[8]+y[17]+y[20]+y[33]+y[34];
y[36]=lambda*lambda;
y[37]=4.*x0*x2*y[1]*y[2];
y[38]=2.*x0*y[1]*y[2]*y[18];
y[39]=-2.*x0*x2*y[1]*y[4];
y[40]=y[5]+y[6]+y[7]+y[8]+y[37]+y[38]+y[39];
y[41]=y[1]*y[2]*y[3];
y[42]=2.*x0*x2*y[1]*y[2];
y[43]=y[1]*y[2]*y[3]*y[18];
y[44]=-(x0*y[1]*y[4]);
y[45]=-(x2*y[1]*y[3]*y[4]);
y[46]=y[6]+y[12]+y[21]+y[41]+y[42]+y[43]+y[44]+y[45];
y[47]=lrs[2];
y[48]=4.*x1*x2*y[1]*y[2];
y[49]=2.*x1*y[1]*y[2]*y[18];
y[50]=-2.*x1*x2*y[1]*y[4];
y[51]=y[17]+y[48]+y[49]+y[50];
y[52]=-(lambda*MYI*x0*y[16]*y[19]*y[51]);
y[53]=-(lambda*MYI*y[16]*y[19]*y[29]);
y[54]=lambda*MYI*x0*y[19]*y[29];
y[55]=1.+y[52]+y[53]+y[54];
y[56]=-(lambda*MYI*y[10]*y[30]*y[46]);
y[57]=lambda*MYI*x1*y[30]*y[46];
y[58]=1.+y[56]+y[57];
y[59]=2.*x1*y[1]*y[2]*y[3];
y[60]=2.*x1*x2*y[1]*y[2]*y[3];
y[61]=-(x1*y[1]*y[3]*y[4]);
y[62]=y[6]+y[21]+y[22]+y[42]+y[44]+y[59]+y[60]+y[61];
y[63]=-(lambda*MYI*x1*y[10]*y[30]*y[46]);
y[64]=x1+y[63];
y[65]=-(lambda*MYI*x0*y[16]*y[19]*y[29]);
y[66]=x0+y[65];
y[67]=-(lambda*MYI*x2*y[32]*y[47]*y[62]);
y[68]=x2+y[67];
y[69]=pow(bi,-2);
y[70]=x0*x1*y[10]*y[16]*y[19]*y[30]*y[35]*y[36]*y[40];
y[71]=-(lambda*MYI*x1*y[10]*y[14]*y[30]*y[55]);
y[72]=y[70]+y[71];
y[73]=lambda*MYI*x2*y[14]*y[32]*y[47]*y[72];
y[74]=-(x0*x1*y[10]*y[14]*y[16]*y[19]*y[30]*y[36]*y[40]);
y[75]=lambda*MYI*x0*y[16]*y[19]*y[35]*y[58];
y[76]=y[74]+y[75];
y[77]=-(lambda*MYI*x2*y[32]*y[35]*y[47]*y[76]);
y[78]=pow(y[40],2);
y[79]=x0*x1*y[10]*y[16]*y[19]*y[30]*y[36]*y[78];
y[80]=y[55]*y[58];
y[81]=y[79]+y[80];
y[82]=y[6]+y[59];
y[83]=-(lambda*MYI*x2*y[32]*y[47]*y[82]);
y[84]=-(lambda*MYI*y[32]*y[47]*y[62]);
y[85]=lambda*MYI*x2*y[47]*y[62];
y[86]=1.+y[83]+y[84]+y[85];
y[87]=y[81]*y[86];
y[88]=y[73]+y[77]+y[87];
y[89]=y[1]*y[64];
y[90]=y[1]*y[64]*y[66];
y[91]=y[1]*y[68];
y[92]=y[1]*y[64]*y[66]*y[68];
y[93]=y[1]+y[89]+y[90]+y[91]+y[92];
y[94]=pow(y[93],-2);
y[95]=pow(y[66],2);
y[96]=pow(y[68],2);
FOUT=myLog(bi)*y[69]*y[88]*y[94]+myLog(x0)*y[69]*y[88]*y[94]+myLog(1.+y[53])\
*y[69]*y[88]*y[94]+3.*myLog(y[93])*y[69]*y[88]*y[94]-2.*myLog(y[21]+y[1]*y[\
2]*y[64]+y[1]*y[2]*y[66]+2.*y[1]*y[2]*y[64]*y[66]-y[1]*y[4]*y[64]*y[66]+y[1\
]*y[2]*y[68]+2.*y[1]*y[2]*y[66]*y[68]-y[1]*y[4]*y[66]*y[68]+2.*y[1]*y[2]*y[\
64]*y[66]*y[68]+y[1]*y[2]*y[64]*y[95]+2.*y[1]*y[2]*y[64]*y[68]*y[95]-y[1]*y\
[4]*y[64]*y[68]*y[95]+y[1]*y[2]*y[66]*y[96]+y[1]*y[2]*y[64]*y[95]*y[96])*y[\
69]*y[88]*y[94];
return (FOUT);
}