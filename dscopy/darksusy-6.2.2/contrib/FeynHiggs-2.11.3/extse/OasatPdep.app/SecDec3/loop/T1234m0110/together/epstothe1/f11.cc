#include "intfile.hh"

dcmplx Pf11(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
dcmplx y[98];
dcmplx FOUT;
dcmplx MYI(0.,1.);
y[1]=1./bi;
y[2]=em[0];
y[3]=x1*x1;
y[4]=esx[0];
y[5]=x0*x0;
y[6]=-x0;
y[7]=1.+y[6];
y[8]=y[1]*y[2];
y[9]=-(x1*y[1]*y[4]);
y[10]=2.*x1*y[1]*y[2];
y[11]=x2*y[1]*y[2];
y[12]=lrs[0];
y[13]=-x1;
y[14]=1.+y[13];
y[15]=x0*y[1]*y[2];
y[16]=y[1]*y[2]*y[5];
y[17]=2.*x1*y[1]*y[2]*y[5];
y[18]=-(x0*y[1]*y[4]);
y[19]=y[15]+y[16]+y[17]+y[18];
y[20]=y[1]*y[2]*y[3];
y[21]=x1*x2*y[1]*y[2];
y[22]=2.*x0*x1*x2*y[1]*y[2];
y[23]=2.*x0*x2*y[1]*y[2]*y[3];
y[24]=-(x1*x2*y[1]*y[4]);
y[25]=y[8]+y[9]+y[10]+y[11]+y[20]+y[21]+y[22]+y[23]+y[24];
y[26]=-(lambda*MYI*y[7]*y[12]*y[25]);
y[27]=lrs[1];
y[28]=-x2;
y[29]=1.+y[28];
y[30]=x1*y[1]*y[2];
y[31]=2.*x0*x1*y[1]*y[2];
y[32]=2.*x0*y[1]*y[2]*y[3];
y[33]=y[8]+y[9]+y[30]+y[31]+y[32];
y[34]=lambda*lambda;
y[35]=2.*y[1]*y[2];
y[36]=2.*x0*x2*y[1]*y[2];
y[37]=4.*x0*x1*x2*y[1]*y[2];
y[38]=-(y[1]*y[4]);
y[39]=-(x2*y[1]*y[4]);
y[40]=y[10]+y[11]+y[35]+y[36]+y[37]+y[38]+y[39];
y[41]=2.*x0*y[1]*y[2];
y[42]=x0*x2*y[1]*y[2];
y[43]=x2*y[1]*y[2]*y[5];
y[44]=2.*x1*x2*y[1]*y[2]*y[5];
y[45]=-(x0*x2*y[1]*y[4]);
y[46]=y[8]+y[18]+y[31]+y[38]+y[41]+y[42]+y[43]+y[44]+y[45];
y[47]=lrs[2];
y[48]=2.*x1*x2*y[1]*y[2];
y[49]=2.*x2*y[1]*y[2]*y[3];
y[50]=y[48]+y[49];
y[51]=-(lambda*MYI*x0*y[7]*y[12]*y[50]);
y[52]=lambda*MYI*x0*y[12]*y[25];
y[53]=1.+y[26]+y[51]+y[52];
y[54]=2.*x2*y[1]*y[2]*y[5];
y[55]=y[41]+y[54];
y[56]=-(lambda*MYI*x1*y[14]*y[27]*y[55]);
y[57]=-(lambda*MYI*y[14]*y[27]*y[46]);
y[58]=lambda*MYI*x1*y[27]*y[46];
y[59]=1.+y[56]+y[57]+y[58];
y[60]=x0*x1*y[1]*y[2];
y[61]=x1*y[1]*y[2]*y[5];
y[62]=y[1]*y[2]*y[3]*y[5];
y[63]=-(x0*x1*y[1]*y[4]);
y[64]=y[15]+y[60]+y[61]+y[62]+y[63];
y[65]=-(lambda*MYI*x0*y[7]*y[12]*y[25]);
y[66]=x0+y[65];
y[67]=-(lambda*MYI*x1*y[14]*y[27]*y[46]);
y[68]=x1+y[67];
y[69]=-(lambda*MYI*x2*y[29]*y[47]*y[64]);
y[70]=x2+y[69];
y[71]=pow(bi,-2);
y[72]=1.+y[26];
y[73]=x0*x1*y[7]*y[12]*y[14]*y[27]*y[33]*y[34]*y[40];
y[74]=-(lambda*MYI*x1*y[14]*y[19]*y[27]*y[53]);
y[75]=y[73]+y[74];
y[76]=lambda*MYI*x2*y[19]*y[29]*y[47]*y[75];
y[77]=-(x0*x1*y[7]*y[12]*y[14]*y[19]*y[27]*y[34]*y[40]);
y[78]=lambda*MYI*x0*y[7]*y[12]*y[33]*y[59];
y[79]=y[77]+y[78];
y[80]=-(lambda*MYI*x2*y[29]*y[33]*y[47]*y[79]);
y[81]=pow(y[40],2);
y[82]=x0*x1*y[7]*y[12]*y[14]*y[27]*y[34]*y[81];
y[83]=y[53]*y[59];
y[84]=y[82]+y[83];
y[85]=-(lambda*MYI*y[29]*y[47]*y[64]);
y[86]=lambda*MYI*x2*y[47]*y[64];
y[87]=1.+y[85]+y[86];
y[88]=y[84]*y[87];
y[89]=y[76]+y[80]+y[88];
y[90]=y[1]*y[66];
y[91]=y[1]*y[66]*y[68];
y[92]=y[1]*y[66]*y[70];
y[93]=pow(y[66],2);
y[94]=y[1]*y[68]*y[70]*y[93];
y[95]=y[1]+y[90]+y[91]+y[92]+y[94];
y[96]=pow(y[95],-2);
y[97]=pow(y[68],2);
FOUT=(-2.*x0*myLog(x0)-x0*myLog(x2))*y[71]*y[72]*y[89]*y[96]+x0*(myLog(bi)*y\
[71]*y[72]*y[89]*y[96]-2.*myLog(y[72])*y[71]*y[72]*y[89]*y[96]-myLog(1.+y[8\
5])*y[71]*y[72]*y[89]*y[96]+3.*myLog(y[95])*y[71]*y[72]*y[89]*y[96]-2.*myLo\
g(y[8]+y[1]*y[2]*y[66]+y[1]*y[2]*y[68]-y[1]*y[4]*y[68]+2.*y[1]*y[2]*y[66]*y\
[68]-y[1]*y[4]*y[66]*y[68]+y[1]*y[2]*y[66]*y[70]+y[1]*y[2]*y[66]*y[68]*y[70\
]-y[1]*y[4]*y[66]*y[68]*y[70]+y[1]*y[2]*y[68]*y[70]*y[93]+y[1]*y[2]*y[66]*y\
[97]+y[1]*y[2]*y[70]*y[93]*y[97])*y[71]*y[72]*y[89]*y[96]);
return (FOUT);
}