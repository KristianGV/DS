#include "intfile.hh"

dcmplx Pf3(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
dcmplx y[89];
dcmplx FOUT;
dcmplx MYI(0.,1.);
y[1]=1./bi;
y[2]=em[0];
y[3]=esx[0];
y[4]=-x0;
y[5]=1.+y[4];
y[6]=2.*y[1]*y[2];
y[7]=2.*x0*y[1]*y[2];
y[8]=-(y[1]*y[3]);
y[9]=lrs[0];
y[10]=2.*x1*y[1]*y[2];
y[11]=2.*x0*x1*y[1]*y[2];
y[12]=-(x1*y[1]*y[3]);
y[13]=y[6]+y[7]+y[8]+y[10]+y[11]+y[12];
y[14]=-(lambda*MYI*y[5]*y[9]*y[13]);
y[15]=-x1;
y[16]=1.+y[15];
y[17]=lrs[1];
y[18]=y[1]*y[2];
y[19]=x0*x0;
y[20]=y[1]*y[2]*y[19];
y[21]=-(x0*y[1]*y[3]);
y[22]=y[7]+y[18]+y[20]+y[21];
y[23]=-(lambda*MYI*x0*y[5]*y[9]*y[13]);
y[24]=x0+y[23];
y[25]=-(lambda*MYI*x1*y[16]*y[17]*y[22]);
y[26]=x1+y[25];
y[27]=pow(y[24],2);
y[28]=pow(bi,-2);
y[29]=1.+y[14];
y[30]=lambda*lambda;
y[31]=y[6]+y[7]+y[8];
y[32]=pow(y[31],2);
y[33]=x0*x1*y[5]*y[9]*y[16]*y[17]*y[30]*y[32];
y[34]=y[6]+y[10];
y[35]=-(lambda*MYI*x0*y[5]*y[9]*y[34]);
y[36]=lambda*MYI*x0*y[9]*y[13];
y[37]=1.+y[14]+y[35]+y[36];
y[38]=-(lambda*MYI*y[16]*y[17]*y[22]);
y[39]=lambda*MYI*x1*y[17]*y[22];
y[40]=1.+y[38]+y[39];
y[41]=y[37]*y[40];
y[42]=y[33]+y[41];
y[43]=y[1]*y[24];
y[44]=y[1]*y[26];
y[45]=y[1]*y[24]*y[26];
y[46]=y[1]+y[43]+y[44]+y[45];
y[47]=1./y[46];
y[48]=2.*y[1]*y[2]*y[24];
y[49]=-(y[1]*y[3]*y[24]);
y[50]=y[1]*y[2]*y[27];
y[51]=y[1]*y[2]*y[26];
y[52]=2.*y[1]*y[2]*y[24]*y[26];
y[53]=-(y[1]*y[3]*y[24]*y[26]);
y[54]=y[1]*y[2]*y[26]*y[27];
y[55]=y[18]+y[48]+y[49]+y[50]+y[51]+y[52]+y[53]+y[54];
y[56]=1./y[55];
y[57]=1./x2;
y[58]=x0*y[1]*y[2];
y[59]=2.*x1*x2*y[1]*y[2];
y[60]=lrs[2];
y[61]=-x2;
y[62]=1.+y[61];
y[63]=2.*x2*y[1]*y[2];
y[64]=y[6]+y[7]+y[8]+y[63];
y[65]=x2*y[1]*y[2];
y[66]=y[6]+y[7]+y[8]+y[10]+y[11]+y[12]+y[59]+y[65];
y[67]=-(lambda*MYI*y[5]*y[9]*y[66]);
y[68]=y[10]+y[18];
y[69]=2.*x0*x2*y[1]*y[2];
y[70]=x2*x2;
y[71]=y[1]*y[2]*y[70];
y[72]=-(x2*y[1]*y[3]);
y[73]=y[7]+y[18]+y[20]+y[21]+y[63]+y[69]+y[71]+y[72];
y[74]=pow(y[64],2);
y[75]=lambda*MYI*x0*y[9]*y[66];
y[76]=1.+y[35]+y[67]+y[75];
y[77]=-(lambda*MYI*y[16]*y[17]*y[73]);
y[78]=lambda*MYI*x1*y[17]*y[73];
y[79]=1.+y[77]+y[78];
y[80]=y[10]+y[11]+y[12]+y[18]+y[58]+y[59];
y[81]=-(lambda*MYI*y[60]*y[62]*y[80]);
y[82]=-(lambda*MYI*x0*y[5]*y[9]*y[66]);
y[83]=x0+y[82];
y[84]=-(lambda*MYI*x1*y[16]*y[17]*y[73]);
y[85]=x1+y[84];
y[86]=pow(y[83],2);
y[87]=-(lambda*MYI*x2*y[60]*y[62]*y[80]);
y[88]=x2+y[87];
FOUT=x0*(myLog(bi)*y[28]*y[29]*y[42]*y[47]*y[56]+3.*myLog(y[46])*y[28]*y[29]\
*y[42]*y[47]*y[56]-2.*myLog(y[55])*y[28]*y[29]*y[42]*y[47]*y[56]+myLog(1.-l\
ambda*MYI*(y[10]+y[11]+y[12]+y[18]+y[58])*y[60])*y[28]*y[29]*y[42]*y[47]*y[\
56])-x0*y[28]*y[29]*y[42]*y[47]*y[56]*y[57]+(x0*y[28]*y[57]*(1.+y[67])*(lam\
bda*MYI*x2*y[60]*y[62]*y[64]*(x0*x1*y[5]*y[9]*y[16]*y[17]*y[30]*y[64]*y[68]\
-lambda*MYI*x1*y[16]*y[17]*y[64]*y[76])-lambda*MYI*x2*y[60]*y[62]*y[68]*(-(\
x0*x1*y[5]*y[9]*y[16]*y[17]*y[30]*y[74])+lambda*MYI*x0*y[5]*y[9]*y[68]*y[79\
])+(x0*x1*y[5]*y[9]*y[16]*y[17]*y[30]*y[74]+y[76]*y[79])*(1.-2.*lambda*MYI*\
x1*x2*y[1]*y[2]*y[60]*y[62]+lambda*MYI*x2*y[60]*y[80]+y[81])))/((1.+y[81])*\
(y[1]+y[1]*y[83]+y[1]*y[85]+y[1]*y[83]*y[85]+y[1]*y[85]*y[88])*(y[18]+2.*y[\
1]*y[2]*y[83]-y[1]*y[3]*y[83]+y[1]*y[2]*y[85]+pow(y[88],2)*y[1]*y[2]*y[85]+\
2.*y[1]*y[2]*y[83]*y[85]-y[1]*y[3]*y[83]*y[85]+y[1]*y[2]*y[86]+y[1]*y[2]*y[\
85]*y[86]+y[1]*y[2]*y[88]+y[1]*y[2]*y[83]*y[88]+2.*y[1]*y[2]*y[85]*y[88]-y[\
1]*y[3]*y[85]*y[88]+2.*y[1]*y[2]*y[83]*y[85]*y[88]));
return (FOUT);
}
