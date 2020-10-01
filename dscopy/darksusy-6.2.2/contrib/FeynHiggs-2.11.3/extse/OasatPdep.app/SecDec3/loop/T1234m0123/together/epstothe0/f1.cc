#include "intfile.hh"

dcmplx Pf1(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
dcmplx y[110];
dcmplx FOUT;
dcmplx MYI(0.,1.);
y[1]=pow(bi,-2);
y[2]=-x0;
y[3]=1.+y[2];
y[4]=em[0];
y[5]=1./bi;
y[6]=lrs[0];
y[7]=y[4]*y[5];
y[8]=x1*y[4]*y[5];
y[9]=esx[0];
y[10]=-(y[5]*y[9]);
y[11]=y[7]+y[8]+y[10];
y[12]=-x1;
y[13]=1.+y[12];
y[14]=lrs[1];
y[15]=x0*y[4]*y[5];
y[16]=y[7]+y[15];
y[17]=-(lambda*MYI*x0*y[3]*y[6]*y[11]);
y[18]=x0+y[17];
y[19]=-(lambda*MYI*x1*y[13]*y[14]*y[16]);
y[20]=x1+y[19];
y[21]=lambda*lambda;
y[22]=pow(y[4],2);
y[23]=x0*x1*y[1]*y[3]*y[6]*y[13]*y[14]*y[21]*y[22];
y[24]=-(lambda*MYI*y[3]*y[6]*y[11]);
y[25]=lambda*MYI*x0*y[6]*y[11];
y[26]=1.+y[24]+y[25];
y[27]=-(lambda*MYI*y[13]*y[14]*y[16]);
y[28]=lambda*MYI*x1*y[14]*y[16];
y[29]=1.+y[27]+y[28];
y[30]=y[26]*y[29];
y[31]=y[23]+y[30];
y[32]=y[5]*y[18];
y[33]=y[5]*y[20];
y[34]=y[5]*y[18]*y[20];
y[35]=y[5]+y[32]+y[33]+y[34];
y[36]=pow(y[35],-2);
y[37]=em[1];
y[38]=x0*x0;
y[39]=em[2];
y[40]=1./x2;
y[41]=x0*x1*y[4]*y[5];
y[42]=x0*y[5]*y[37];
y[43]=y[5]*y[37]*y[38];
y[44]=x0*x1*y[5]*y[37];
y[45]=x1*y[5]*y[37]*y[38];
y[46]=y[5]*y[39];
y[47]=x0*y[5]*y[39];
y[48]=x1*y[5]*y[39];
y[49]=x0*x1*y[5]*y[39];
y[50]=-(x0*x1*y[5]*y[9]);
y[51]=lrs[2];
y[52]=-x2;
y[53]=1.+y[52];
y[54]=x2*x2;
y[55]=2.*x2*y[5]*y[37]*y[38];
y[56]=2.*x0*x2*y[5]*y[39];
y[57]=-(x0*y[5]*y[9]);
y[58]=y[15]+y[42]+y[43]+y[46]+y[47]+y[55]+y[56]+y[57];
y[59]=x2*y[5]*y[37];
y[60]=2.*x0*x2*y[5]*y[37];
y[61]=x2*y[5]*y[39];
y[62]=x1*x2*y[4]*y[5];
y[63]=x1*x2*y[5]*y[37];
y[64]=2.*x0*x1*x2*y[5]*y[37];
y[65]=2.*x0*x1*y[5]*y[37]*y[54];
y[66]=x1*x2*y[5]*y[39];
y[67]=x1*y[5]*y[39]*y[54];
y[68]=-(x1*x2*y[5]*y[9]);
y[69]=y[7]+y[8]+y[10]+y[59]+y[60]+y[61]+y[62]+y[63]+y[64]+y[65]+y[66]+y[67]+\
y[68];
y[70]=y[5]*y[37];
y[71]=2.*x0*y[5]*y[37];
y[72]=x1*y[5]*y[37];
y[73]=2.*x0*x1*y[5]*y[37];
y[74]=4.*x0*x1*x2*y[5]*y[37];
y[75]=2.*x1*x2*y[5]*y[39];
y[76]=-(x1*y[5]*y[9]);
y[77]=y[8]+y[46]+y[48]+y[70]+y[71]+y[72]+y[73]+y[74]+y[75]+y[76];
y[78]=x2*y[4]*y[5];
y[79]=2.*x0*y[5]*y[37]*y[54];
y[80]=y[5]*y[39]*y[54];
y[81]=-(x2*y[5]*y[9]);
y[82]=y[7]+y[59]+y[60]+y[61]+y[78]+y[79]+y[80]+y[81];
y[83]=x0*x2*y[4]*y[5];
y[84]=x0*x2*y[5]*y[37];
y[85]=x2*y[5]*y[37]*y[38];
y[86]=y[5]*y[37]*y[38]*y[54];
y[87]=x0*x2*y[5]*y[39];
y[88]=x0*y[5]*y[39]*y[54];
y[89]=-(x0*x2*y[5]*y[9]);
y[90]=y[7]+y[15]+y[61]+y[83]+y[84]+y[85]+y[86]+y[87]+y[88]+y[89];
y[91]=2.*x2*y[5]*y[37];
y[92]=2.*x1*x2*y[5]*y[37];
y[93]=2.*x1*y[5]*y[37]*y[54];
y[94]=y[91]+y[92]+y[93];
y[95]=-(lambda*MYI*x0*y[3]*y[6]*y[94]);
y[96]=-(lambda*MYI*y[3]*y[6]*y[69]);
y[97]=lambda*MYI*x0*y[6]*y[69];
y[98]=1.+y[95]+y[96]+y[97];
y[99]=-(lambda*MYI*y[13]*y[14]*y[90]);
y[100]=lambda*MYI*x1*y[14]*y[90];
y[101]=1.+y[99]+y[100];
y[102]=2.*x1*x2*y[5]*y[37]*y[38];
y[103]=2.*x0*x1*x2*y[5]*y[39];
y[104]=y[41]+y[42]+y[43]+y[44]+y[45]+y[46]+y[47]+y[48]+y[49]+y[50]+y[102]+y[\
103];
y[105]=-(lambda*MYI*y[51]*y[53]*y[104]);
y[106]=-(lambda*MYI*x0*y[3]*y[6]*y[69]);
y[107]=x0+y[106];
y[108]=-(lambda*MYI*x1*y[13]*y[14]*y[90]);
y[109]=x1+y[108];
FOUT=myLog(bi)*y[1]*y[31]*y[36]-myLog(x1)*y[1]*y[31]*y[36]-2.*myLog(y[7]+y[1\
0]+y[4]*y[5]*y[18]-y[5]*y[9]*y[18]+y[4]*y[5]*y[20]+y[4]*y[5]*y[18]*y[20])*y\
[1]*y[31]*y[36]-myLog(1.+y[27])*y[1]*y[31]*y[36]+3.*myLog(y[35])*y[1]*y[31]\
*y[36]+myLog(1.-lambda*MYI*(y[41]+y[42]+y[43]+y[44]+y[45]+y[46]+y[47]+y[48]\
+y[49]+y[50])*y[51])*y[1]*y[31]*y[36]-y[1]*y[31]*y[36]*y[40]+(pow(y[5]+y[5]\
*y[107]+y[5]*y[109]+y[5]*y[107]*y[109]+y[5]*(x2-lambda*MYI*x2*y[51]*y[53]*y\
[104])*y[107]*y[109],-2)*y[1]*y[40]*(lambda*MYI*x2*y[51]*y[53]*y[58]*(x0*x1\
*y[3]*y[6]*y[13]*y[14]*y[21]*y[77]*y[82]-lambda*MYI*x1*y[13]*y[14]*y[58]*y[\
98])-lambda*MYI*x2*y[51]*y[53]*y[77]*(-(x0*x1*y[3]*y[6]*y[13]*y[14]*y[21]*y\
[58]*y[82])+lambda*MYI*x0*y[3]*y[6]*y[77]*y[101])+(x0*x1*pow(y[82],2)*y[3]*\
y[6]*y[13]*y[14]*y[21]+y[98]*y[101])*(1.-lambda*MYI*x2*(2.*x1*y[5]*y[37]*y[\
38]+2.*x0*x1*y[5]*y[39])*y[51]*y[53]+lambda*MYI*x2*y[51]*y[104]+y[105])))/(\
1.+y[105]);
return (FOUT);
}
