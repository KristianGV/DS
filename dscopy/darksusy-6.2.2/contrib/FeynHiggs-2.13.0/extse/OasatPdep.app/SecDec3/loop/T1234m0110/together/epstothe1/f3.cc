#include "intfile.hh"

dcmplx Pf3(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
dcmplx y[113];
dcmplx FOUT;
dcmplx MYI(0.,1.);
y[1]=1./bi;
y[2]=-x0;
y[3]=1.+y[2];
y[4]=em[0];
y[5]=y[1]*y[4];
y[6]=esx[0];
y[7]=lrs[0];
y[8]=x1*y[1]*y[4];
y[9]=-(x1*y[1]*y[6]);
y[10]=y[5]+y[8]+y[9];
y[11]=-x1;
y[12]=1.+y[11];
y[13]=-(y[1]*y[6]);
y[14]=lrs[1];
y[15]=x0*y[1]*y[4];
y[16]=-(x0*y[1]*y[6]);
y[17]=y[5]+y[13]+y[15]+y[16];
y[18]=-(lambda*MYI*x0*y[3]*y[7]*y[10]);
y[19]=x0+y[18];
y[20]=-(lambda*MYI*x1*y[12]*y[14]*y[17]);
y[21]=x1+y[20];
y[22]=1./x2;
y[23]=pow(bi,-2);
y[24]=lambda*lambda;
y[25]=y[5]+y[13];
y[26]=pow(y[25],2);
y[27]=x0*x1*y[3]*y[7]*y[12]*y[14]*y[24]*y[26];
y[28]=-(lambda*MYI*y[3]*y[7]*y[10]);
y[29]=lambda*MYI*x0*y[7]*y[10];
y[30]=1.+y[28]+y[29];
y[31]=-(lambda*MYI*y[12]*y[14]*y[17]);
y[32]=lambda*MYI*x1*y[14]*y[17];
y[33]=1.+y[31]+y[32];
y[34]=y[30]*y[33];
y[35]=y[27]+y[34];
y[36]=y[1]*y[19];
y[37]=y[1]*y[21];
y[38]=y[1]*y[19]*y[21];
y[39]=y[1]+y[36]+y[37]+y[38];
y[40]=pow(y[39],-2);
y[41]=myLog(bi);
y[42]=2.*x0*y[1]*y[4];
y[43]=x0*x1*y[1]*y[4];
y[44]=y[5]+y[8]+y[16]+y[42]+y[43];
y[45]=lrs[2];
y[46]=-(lambda*MYI*y[44]*y[45]);
y[47]=1.+y[46];
y[48]=myLog(y[47]);
y[49]=myLog(y[39]);
y[50]=y[1]*y[4]*y[19];
y[51]=y[1]*y[4]*y[21];
y[52]=-(y[1]*y[6]*y[21]);
y[53]=y[1]*y[4]*y[19]*y[21];
y[54]=-(y[1]*y[6]*y[19]*y[21]);
y[55]=y[5]+y[50]+y[51]+y[52]+y[53]+y[54];
y[56]=myLog(y[55]);
y[57]=myLog(x2);
y[58]=-x2;
y[59]=1.+y[58];
y[60]=y[5]+y[15];
y[61]=2.*x2*y[1]*y[4];
y[62]=x1*x2*y[1]*y[4];
y[63]=x2*x2;
y[64]=y[1]*y[4]*y[63];
y[65]=-(x2*y[1]*y[6]);
y[66]=y[5]+y[8]+y[9]+y[61]+y[62]+y[64]+y[65];
y[67]=2.*y[1]*y[4];
y[68]=y[8]+y[13]+y[61]+y[67];
y[69]=x2*y[1]*y[4];
y[70]=y[5]+y[13]+y[69];
y[71]=x0*x2*y[1]*y[4];
y[72]=y[5]+y[13]+y[15]+y[16]+y[69]+y[71];
y[73]=-(lambda*MYI*y[3]*y[7]*y[66]);
y[74]=lambda*MYI*x0*y[7]*y[66];
y[75]=1.+y[73]+y[74];
y[76]=-(lambda*MYI*y[12]*y[14]*y[72]);
y[77]=lambda*MYI*x1*y[14]*y[72];
y[78]=1.+y[76]+y[77];
y[79]=2.*x0*x2*y[1]*y[4];
y[80]=y[5]+y[8]+y[16]+y[42]+y[43]+y[79];
y[81]=-(lambda*MYI*y[45]*y[59]*y[80]);
y[82]=-(lambda*MYI*x0*y[3]*y[7]*y[66]);
y[83]=x0+y[82];
y[84]=-(lambda*MYI*x1*y[12]*y[14]*y[72]);
y[85]=x1+y[84];
y[86]=1.+y[81];
y[87]=1./y[86];
y[88]=x0*x1*y[3]*y[7]*y[12]*y[14]*y[24]*y[68]*y[70];
y[89]=-(lambda*MYI*x1*y[12]*y[14]*y[60]*y[75]);
y[90]=y[88]+y[89];
y[91]=lambda*MYI*x2*y[45]*y[59]*y[60]*y[90];
y[92]=-(x0*x1*y[3]*y[7]*y[12]*y[14]*y[24]*y[60]*y[70]);
y[93]=lambda*MYI*x0*y[3]*y[7]*y[68]*y[78];
y[94]=y[92]+y[93];
y[95]=-(lambda*MYI*x2*y[45]*y[59]*y[68]*y[94]);
y[96]=pow(y[70],2);
y[97]=x0*x1*y[3]*y[7]*y[12]*y[14]*y[24]*y[96];
y[98]=y[75]*y[78];
y[99]=y[97]+y[98];
y[100]=-2.*lambda*MYI*x0*x2*y[1]*y[4]*y[45]*y[59];
y[101]=lambda*MYI*x2*y[45]*y[80];
y[102]=1.+y[81]+y[100]+y[101];
y[103]=y[99]*y[102];
y[104]=y[91]+y[95]+y[103];
y[105]=y[1]*y[83];
y[106]=y[1]*y[85];
y[107]=y[1]*y[83]*y[85];
y[108]=-(lambda*MYI*x2*y[45]*y[59]*y[80]);
y[109]=x2+y[108];
y[110]=y[1]*y[83]*y[109];
y[111]=y[1]+y[105]+y[106]+y[107]+y[110];
y[112]=pow(y[111],-2);
FOUT=-(y[22]*(y[23]*y[35]*y[40]*y[41]+y[23]*y[35]*y[40]*y[48]+3.*y[23]*y[35]\
*y[40]*y[49]-2.*y[23]*y[35]*y[40]*y[56]))+0.5*(y[35]*y[40]*(pow(y[41],2)*y[\
23]+pow(y[48],2)*y[23]+2.*y[23]*y[41]*y[48])+2.*(y[23]*y[35]*y[41]+y[23]*y[\
35]*y[48])*(3.*y[40]*y[49]-2.*y[40]*y[56])+y[23]*y[35]*(9.*pow(y[49],2)*y[4\
0]+4.*pow(y[56],2)*y[40]-12.*y[40]*y[49]*y[56]))-y[22]*y[23]*y[35]*y[40]*y[\
57]+y[22]*y[23]*y[57]*y[87]*y[104]*y[112]+y[22]*(myLog(y[86])*y[23]*y[87]*y\
[104]*y[112]-2.*myLog(y[5]+y[1]*y[4]*y[83]+pow(y[109],2)*y[1]*y[4]*y[83]+y[\
1]*y[4]*y[85]-y[1]*y[6]*y[85]+y[1]*y[4]*y[83]*y[85]-y[1]*y[6]*y[83]*y[85]+y\
[1]*y[4]*y[109]+2.*y[1]*y[4]*y[83]*y[109]-y[1]*y[6]*y[83]*y[109]+y[1]*y[4]*\
y[85]*y[109]+y[1]*y[4]*y[83]*y[85]*y[109])*y[23]*y[87]*y[104]*y[112]+3.*myL\
og(y[111])*y[23]*y[87]*y[104]*y[112]+y[23]*y[41]*y[87]*y[104]*y[112]);
return (FOUT);
}