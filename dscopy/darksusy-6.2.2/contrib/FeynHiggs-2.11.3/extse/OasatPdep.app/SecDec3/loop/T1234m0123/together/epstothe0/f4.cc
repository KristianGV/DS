#include "intfile.hh"

dcmplx Pf4(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
dcmplx y[99];
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
y[22]=pow(bi,-2);
y[23]=lambda*lambda;
y[24]=y[5]+y[13];
y[25]=pow(y[24],2);
y[26]=x0*x1*y[3]*y[7]*y[12]*y[14]*y[23]*y[25];
y[27]=-(lambda*MYI*y[3]*y[7]*y[10]);
y[28]=lambda*MYI*x0*y[7]*y[10];
y[29]=1.+y[27]+y[28];
y[30]=-(lambda*MYI*y[12]*y[14]*y[17]);
y[31]=lambda*MYI*x1*y[14]*y[17];
y[32]=1.+y[30]+y[31];
y[33]=y[29]*y[32];
y[34]=y[26]+y[33];
y[35]=y[1]*y[19];
y[36]=y[1]*y[21];
y[37]=y[1]*y[19]*y[21];
y[38]=y[1]+y[35]+y[36]+y[37];
y[39]=pow(y[38],-2);
y[40]=em[1];
y[41]=em[2];
y[42]=x0*x0;
y[43]=1./x2;
y[44]=y[1]*y[40];
y[45]=x0*y[1]*y[40];
y[46]=x1*y[1]*y[40];
y[47]=x0*x1*y[1]*y[40];
y[48]=x0*y[1]*y[41];
y[49]=y[1]*y[41]*y[42];
y[50]=x0*x1*y[1]*y[41];
y[51]=x1*y[1]*y[41]*y[42];
y[52]=lrs[2];
y[53]=-x2;
y[54]=1.+y[53];
y[55]=y[44]+y[45]+y[48]+y[49];
y[56]=x2*y[1]*y[40];
y[57]=x2*x2;
y[58]=x2*y[1]*y[41];
y[59]=2.*x0*x2*y[1]*y[41];
y[60]=x2*y[1]*y[4];
y[61]=x1*x2*y[1]*y[40];
y[62]=y[1]*y[40]*y[57];
y[63]=x1*x2*y[1]*y[41];
y[64]=2.*x0*x1*x2*y[1]*y[41];
y[65]=2.*x0*y[1]*y[41]*y[57];
y[66]=-(x2*y[1]*y[6]);
y[67]=y[5]+y[8]+y[9]+y[56]+y[58]+y[59]+y[60]+y[61]+y[62]+y[63]+y[64]+y[65]+y\
[66];
y[68]=2.*x2*y[1]*y[40];
y[69]=y[1]*y[41];
y[70]=2.*x0*y[1]*y[41];
y[71]=x1*y[1]*y[41];
y[72]=2.*x0*x1*y[1]*y[41];
y[73]=4.*x0*x2*y[1]*y[41];
y[74]=y[5]+y[13]+y[44]+y[46]+y[68]+y[69]+y[70]+y[71]+y[72]+y[73];
y[75]=y[5]+y[13]+y[56]+y[58]+y[59];
y[76]=x0*x2*y[1]*y[40];
y[77]=x0*x2*y[1]*y[41];
y[78]=x2*y[1]*y[41]*y[42];
y[79]=y[5]+y[13]+y[15]+y[16]+y[56]+y[76]+y[77]+y[78];
y[80]=2.*x2*y[1]*y[41];
y[81]=2.*x1*x2*y[1]*y[41];
y[82]=2.*y[1]*y[41]*y[57];
y[83]=y[80]+y[81]+y[82];
y[84]=-(lambda*MYI*x0*y[3]*y[7]*y[83]);
y[85]=-(lambda*MYI*y[3]*y[7]*y[67]);
y[86]=lambda*MYI*x0*y[7]*y[67];
y[87]=1.+y[84]+y[85]+y[86];
y[88]=-(lambda*MYI*y[12]*y[14]*y[79]);
y[89]=lambda*MYI*x1*y[14]*y[79];
y[90]=1.+y[88]+y[89];
y[91]=2.*x0*x2*y[1]*y[40];
y[92]=2.*x2*y[1]*y[41]*y[42];
y[93]=y[15]+y[16]+y[44]+y[45]+y[46]+y[47]+y[48]+y[49]+y[50]+y[51]+y[91]+y[92\
];
y[94]=-(lambda*MYI*y[52]*y[54]*y[93]);
y[95]=-(lambda*MYI*x0*y[3]*y[7]*y[67]);
y[96]=x0+y[95];
y[97]=-(lambda*MYI*x1*y[12]*y[14]*y[79]);
y[98]=x1+y[97];
FOUT=myLog(bi)*y[22]*y[34]*y[39]-2.*myLog(y[5]+y[1]*y[4]*y[19]+y[1]*y[4]*y[2\
1]-y[1]*y[6]*y[21]+y[1]*y[4]*y[19]*y[21]-y[1]*y[6]*y[19]*y[21])*y[22]*y[34]\
*y[39]+3.*myLog(y[38])*y[22]*y[34]*y[39]+myLog(1.-lambda*MYI*(y[15]+y[16]+y\
[44]+y[45]+y[46]+y[47]+y[48]+y[49]+y[50]+y[51])*y[52])*y[22]*y[34]*y[39]-y[\
22]*y[34]*y[39]*y[43]+(pow(y[1]+y[1]*y[96]+y[1]*(x2-lambda*MYI*x2*y[52]*y[5\
4]*y[93])*y[96]+y[1]*y[98]+y[1]*y[96]*y[98],-2)*y[22]*y[43]*(lambda*MYI*x2*\
y[52]*y[54]*y[55]*(x0*x1*y[3]*y[7]*y[12]*y[14]*y[23]*y[74]*y[75]-lambda*MYI\
*x1*y[12]*y[14]*y[55]*y[87])-lambda*MYI*x2*y[52]*y[54]*y[74]*(-(x0*x1*y[3]*\
y[7]*y[12]*y[14]*y[23]*y[55]*y[75])+lambda*MYI*x0*y[3]*y[7]*y[74]*y[90])+(x\
0*x1*pow(y[75],2)*y[3]*y[7]*y[12]*y[14]*y[23]+y[87]*y[90])*(1.-lambda*MYI*x\
2*(2.*x0*y[1]*y[40]+2.*y[1]*y[41]*y[42])*y[52]*y[54]+lambda*MYI*x2*y[52]*y[\
93]+y[94])))/(1.+y[94]);
return (FOUT);
}
