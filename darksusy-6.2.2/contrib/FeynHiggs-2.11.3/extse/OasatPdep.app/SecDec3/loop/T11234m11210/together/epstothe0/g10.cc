#include "intfile.hh"

double Pr10(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
double y[54];
double FOUT;
y[1]=1./bi;
y[2]=em[0];
y[3]=x1*x1;
y[4]=em[1];
y[5]=esx[0];
y[6]=-x0;
y[7]=1.+y[6];
y[8]=x1*y[1]*y[2];
y[9]=y[1]*y[2]*y[3];
y[10]=x2*y[1]*y[2];
y[11]=2.*x1*x2*y[1]*y[2];
y[12]=2.*x0*x1*x2*y[1]*y[2];
y[13]=2.*x0*x2*y[1]*y[2]*y[3];
y[14]=y[1]*y[4];
y[15]=x1*y[1]*y[4];
y[16]=x2*y[1]*y[4];
y[17]=2.*x0*x2*y[1]*y[4];
y[18]=2.*x0*x1*x2*y[1]*y[4];
y[19]=-(x1*y[1]*y[5]);
y[20]=-(x2*y[1]*y[5]);
y[21]=-2.*x0*x1*x2*y[1]*y[5];
y[22]=y[8]+y[9]+y[10]+y[11]+y[12]+y[13]+y[14]+y[15]+y[16]+y[17]+y[18]+y[19]+\
y[20]+y[21];
y[23]=lrs[0];
y[24]=x0*x0;
y[25]=-x1;
y[26]=1.+y[25];
y[27]=y[1]*y[2];
y[28]=x0*y[1]*y[2];
y[29]=2.*x0*x1*y[1]*y[2];
y[30]=2.*x0*x2*y[1]*y[2];
y[31]=x2*y[1]*y[2]*y[24];
y[32]=2.*x1*x2*y[1]*y[2]*y[24];
y[33]=x0*y[1]*y[4];
y[34]=x2*y[1]*y[4]*y[24];
y[35]=-(x0*y[1]*y[5]);
y[36]=-(x2*y[1]*y[5]*y[24]);
y[37]=y[27]+y[28]+y[29]+y[30]+y[31]+y[32]+y[33]+y[34]+y[35]+y[36];
y[38]=lrs[1];
y[39]=-x2;
y[40]=1.+y[39];
y[41]=x1*y[1]*y[2]*y[24];
y[42]=y[1]*y[2]*y[3]*y[24];
y[43]=y[1]*y[4]*y[24];
y[44]=x1*y[1]*y[4]*y[24];
y[45]=-(x1*y[1]*y[5]*y[24]);
y[46]=y[27]+y[28]+y[29]+y[33]+y[35]+y[41]+y[42]+y[43]+y[44]+y[45];
y[47]=lrs[2];
y[48]=pow(y[7],2);
y[49]=pow(y[22],2);
y[50]=pow(y[23],2);
y[51]=pow(y[26],2);
y[52]=pow(y[37],2);
y[53]=pow(y[38],2);
FOUT=(2.*x0*x1*x2*y[1]*y[2]*y[7]*y[22]*y[23]*y[26]*y[37]*y[38]*y[40]*y[46]*y\
[47]+2.*x1*x2*y[1]*y[2]*y[7]*y[22]*y[23]*y[24]*y[26]*y[37]*y[38]*y[40]*y[46\
]*y[47]+4.*x2*y[1]*y[2]*y[3]*y[7]*y[22]*y[23]*y[24]*y[26]*y[37]*y[38]*y[40]\
*y[46]*y[47]+2.*x1*x2*y[1]*y[4]*y[7]*y[22]*y[23]*y[24]*y[26]*y[37]*y[38]*y[\
40]*y[46]*y[47]-2.*x1*x2*y[1]*y[5]*y[7]*y[22]*y[23]*y[24]*y[26]*y[37]*y[38]\
*y[40]*y[46]*y[47]+x1*x2*y[1]*y[2]*y[24]*y[26]*y[37]*y[38]*y[48]*y[49]*y[50\
]+2.*x2*y[1]*y[2]*y[3]*y[24]*y[26]*y[37]*y[38]*y[48]*y[49]*y[50]+x1*x2*y[1]\
*y[4]*y[24]*y[26]*y[37]*y[38]*y[48]*y[49]*y[50]-x1*x2*y[1]*y[5]*y[24]*y[26]\
*y[37]*y[38]*y[48]*y[49]*y[50]+x1*x2*y[1]*y[2]*y[24]*y[40]*y[46]*y[47]*y[48\
]*y[49]*y[50]+x2*y[1]*y[2]*y[3]*y[24]*y[40]*y[46]*y[47]*y[48]*y[49]*y[50]+x\
2*y[1]*y[4]*y[24]*y[40]*y[46]*y[47]*y[48]*y[49]*y[50]+x1*x2*y[1]*y[4]*y[24]\
*y[40]*y[46]*y[47]*y[48]*y[49]*y[50]-x1*x2*y[1]*y[5]*y[24]*y[40]*y[46]*y[47\
]*y[48]*y[49]*y[50]+x0*y[1]*y[2]*y[3]*y[7]*y[22]*y[23]*y[51]*y[52]*y[53]+2.\
*x2*y[1]*y[2]*y[3]*y[7]*y[22]*y[23]*y[24]*y[51]*y[52]*y[53]+x2*y[1]*y[2]*y[\
3]*y[24]*y[40]*y[46]*y[47]*y[51]*y[52]*y[53])/(-(x0*x1*y[1]*y[2]*y[7]*y[22]\
*y[23])-x0*x2*y[1]*y[2]*y[7]*y[22]*y[23]-2.*x0*x1*x2*y[1]*y[2]*y[7]*y[22]*y\
[23]-x0*y[1]*y[2]*y[3]*y[7]*y[22]*y[23]-x0*y[1]*y[4]*y[7]*y[22]*y[23]-x0*x1\
*y[1]*y[4]*y[7]*y[22]*y[23]-x0*x2*y[1]*y[4]*y[7]*y[22]*y[23]+x0*x1*y[1]*y[5\
]*y[7]*y[22]*y[23]+x0*x2*y[1]*y[5]*y[7]*y[22]*y[23]-2.*x1*x2*y[1]*y[2]*y[7]\
*y[22]*y[23]*y[24]-2.*x2*y[1]*y[2]*y[3]*y[7]*y[22]*y[23]*y[24]-2.*x2*y[1]*y\
[4]*y[7]*y[22]*y[23]*y[24]-2.*x1*x2*y[1]*y[4]*y[7]*y[22]*y[23]*y[24]+2.*x1*\
x2*y[1]*y[5]*y[7]*y[22]*y[23]*y[24]-x1*y[1]*y[2]*y[26]*y[37]*y[38]-x0*x1*y[\
1]*y[2]*y[26]*y[37]*y[38]-2.*x0*x1*x2*y[1]*y[2]*y[26]*y[37]*y[38]-2.*x0*y[1\
]*y[2]*y[3]*y[26]*y[37]*y[38]-x0*x1*y[1]*y[4]*y[26]*y[37]*y[38]+x0*x1*y[1]*\
y[5]*y[26]*y[37]*y[38]-x1*x2*y[1]*y[2]*y[24]*y[26]*y[37]*y[38]-2.*x2*y[1]*y\
[2]*y[3]*y[24]*y[26]*y[37]*y[38]-x1*x2*y[1]*y[4]*y[24]*y[26]*y[37]*y[38]+x1\
*x2*y[1]*y[5]*y[24]*y[26]*y[37]*y[38]-x2*y[1]*y[2]*y[40]*y[46]*y[47]-x0*x2*\
y[1]*y[2]*y[40]*y[46]*y[47]-2.*x0*x1*x2*y[1]*y[2]*y[40]*y[46]*y[47]-x0*x2*y\
[1]*y[4]*y[40]*y[46]*y[47]+x0*x2*y[1]*y[5]*y[40]*y[46]*y[47]-x1*x2*y[1]*y[2\
]*y[24]*y[40]*y[46]*y[47]-x2*y[1]*y[2]*y[3]*y[24]*y[40]*y[46]*y[47]-x2*y[1]\
*y[4]*y[24]*y[40]*y[46]*y[47]-x1*x2*y[1]*y[4]*y[24]*y[40]*y[46]*y[47]+x1*x2\
*y[1]*y[5]*y[24]*y[40]*y[46]*y[47]);
return (FOUT);
}
double Pm10(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
double y[54];
double FOUT;
y[1]=1./bi;
y[2]=em[0];
y[3]=x0*x0;
y[4]=x1*x1;
y[5]=em[1];
y[6]=esx[0];
y[7]=y[1]*y[2];
y[8]=x0*y[1]*y[2];
y[9]=2.*x0*x1*y[1]*y[2];
y[10]=x0*y[1]*y[5];
y[11]=-(x0*y[1]*y[6]);
y[12]=-x0;
y[13]=1.+y[12];
y[14]=x1*y[1]*y[2];
y[15]=y[1]*y[2]*y[4];
y[16]=x2*y[1]*y[2];
y[17]=2.*x1*x2*y[1]*y[2];
y[18]=2.*x0*x1*x2*y[1]*y[2];
y[19]=2.*x0*x2*y[1]*y[2]*y[4];
y[20]=y[1]*y[5];
y[21]=x1*y[1]*y[5];
y[22]=x2*y[1]*y[5];
y[23]=2.*x0*x2*y[1]*y[5];
y[24]=2.*x0*x1*x2*y[1]*y[5];
y[25]=-(x1*y[1]*y[6]);
y[26]=-(x2*y[1]*y[6]);
y[27]=-2.*x0*x1*x2*y[1]*y[6];
y[28]=y[14]+y[15]+y[16]+y[17]+y[18]+y[19]+y[20]+y[21]+y[22]+y[23]+y[24]+y[25\
]+y[26]+y[27];
y[29]=lrs[0];
y[30]=-x1;
y[31]=1.+y[30];
y[32]=2.*x0*x2*y[1]*y[2];
y[33]=x2*y[1]*y[2]*y[3];
y[34]=2.*x1*x2*y[1]*y[2]*y[3];
y[35]=x2*y[1]*y[3]*y[5];
y[36]=-(x2*y[1]*y[3]*y[6]);
y[37]=y[7]+y[8]+y[9]+y[10]+y[11]+y[32]+y[33]+y[34]+y[35]+y[36];
y[38]=lrs[1];
y[39]=-x2;
y[40]=1.+y[39];
y[41]=x1*y[1]*y[2]*y[3];
y[42]=y[1]*y[2]*y[3]*y[4];
y[43]=y[1]*y[3]*y[5];
y[44]=x1*y[1]*y[3]*y[5];
y[45]=-(x1*y[1]*y[3]*y[6]);
y[46]=y[7]+y[8]+y[9]+y[10]+y[11]+y[41]+y[42]+y[43]+y[44]+y[45];
y[47]=lrs[2];
y[48]=pow(y[13],2);
y[49]=pow(y[28],2);
y[50]=pow(y[29],2);
y[51]=pow(y[31],2);
y[52]=pow(y[37],2);
y[53]=pow(y[38],2);
FOUT=pow(lambda*(-(x0*x1*y[1]*y[2]*y[13]*y[28]*y[29])-x0*x2*y[1]*y[2]*y[13]*\
y[28]*y[29]-2.*x0*x1*x2*y[1]*y[2]*y[13]*y[28]*y[29]-2.*x1*x2*y[1]*y[2]*y[3]\
*y[13]*y[28]*y[29]-x0*y[1]*y[2]*y[4]*y[13]*y[28]*y[29]-2.*x2*y[1]*y[2]*y[3]\
*y[4]*y[13]*y[28]*y[29]-x0*y[1]*y[5]*y[13]*y[28]*y[29]-x0*x1*y[1]*y[5]*y[13\
]*y[28]*y[29]-x0*x2*y[1]*y[5]*y[13]*y[28]*y[29]-2.*x2*y[1]*y[3]*y[5]*y[13]*\
y[28]*y[29]-2.*x1*x2*y[1]*y[3]*y[5]*y[13]*y[28]*y[29]+x0*x1*y[1]*y[6]*y[13]\
*y[28]*y[29]+x0*x2*y[1]*y[6]*y[13]*y[28]*y[29]+2.*x1*x2*y[1]*y[3]*y[6]*y[13\
]*y[28]*y[29]-x1*y[1]*y[2]*y[31]*y[37]*y[38]-x0*x1*y[1]*y[2]*y[31]*y[37]*y[\
38]-2.*x0*x1*x2*y[1]*y[2]*y[31]*y[37]*y[38]-x1*x2*y[1]*y[2]*y[3]*y[31]*y[37\
]*y[38]-2.*x0*y[1]*y[2]*y[4]*y[31]*y[37]*y[38]-2.*x2*y[1]*y[2]*y[3]*y[4]*y[\
31]*y[37]*y[38]-x0*x1*y[1]*y[5]*y[31]*y[37]*y[38]-x1*x2*y[1]*y[3]*y[5]*y[31\
]*y[37]*y[38]+x0*x1*y[1]*y[6]*y[31]*y[37]*y[38]+x1*x2*y[1]*y[3]*y[6]*y[31]*\
y[37]*y[38]-x2*y[1]*y[2]*y[40]*y[46]*y[47]-x0*x2*y[1]*y[2]*y[40]*y[46]*y[47\
]-2.*x0*x1*x2*y[1]*y[2]*y[40]*y[46]*y[47]-x1*x2*y[1]*y[2]*y[3]*y[40]*y[46]*\
y[47]-x2*y[1]*y[2]*y[3]*y[4]*y[40]*y[46]*y[47]-x0*x2*y[1]*y[5]*y[40]*y[46]*\
y[47]-x2*y[1]*y[3]*y[5]*y[40]*y[46]*y[47]-x1*x2*y[1]*y[3]*y[5]*y[40]*y[46]*\
y[47]+x0*x2*y[1]*y[6]*y[40]*y[46]*y[47]+x1*x2*y[1]*y[3]*y[6]*y[40]*y[46]*y[\
47])-x2*pow(lambda,5)*y[1]*y[2]*y[3]*y[4]*y[40]*y[46]*y[47]*y[48]*y[49]*y[5\
0]*y[51]*y[52]*y[53]+pow(lambda,3)*(2.*x0*x1*x2*y[1]*y[2]*y[13]*y[28]*y[29]\
*y[31]*y[37]*y[38]*y[40]*y[46]*y[47]+2.*x1*x2*y[1]*y[2]*y[3]*y[13]*y[28]*y[\
29]*y[31]*y[37]*y[38]*y[40]*y[46]*y[47]+4.*x2*y[1]*y[2]*y[3]*y[4]*y[13]*y[2\
8]*y[29]*y[31]*y[37]*y[38]*y[40]*y[46]*y[47]+2.*x1*x2*y[1]*y[3]*y[5]*y[13]*\
y[28]*y[29]*y[31]*y[37]*y[38]*y[40]*y[46]*y[47]-2.*x1*x2*y[1]*y[3]*y[6]*y[1\
3]*y[28]*y[29]*y[31]*y[37]*y[38]*y[40]*y[46]*y[47]+x1*x2*y[1]*y[2]*y[3]*y[3\
1]*y[37]*y[38]*y[48]*y[49]*y[50]+2.*x2*y[1]*y[2]*y[3]*y[4]*y[31]*y[37]*y[38\
]*y[48]*y[49]*y[50]+x1*x2*y[1]*y[3]*y[5]*y[31]*y[37]*y[38]*y[48]*y[49]*y[50\
]-x1*x2*y[1]*y[3]*y[6]*y[31]*y[37]*y[38]*y[48]*y[49]*y[50]+x1*x2*y[1]*y[2]*\
y[3]*y[40]*y[46]*y[47]*y[48]*y[49]*y[50]+x2*y[1]*y[2]*y[3]*y[4]*y[40]*y[46]\
*y[47]*y[48]*y[49]*y[50]+x2*y[1]*y[3]*y[5]*y[40]*y[46]*y[47]*y[48]*y[49]*y[\
50]+x1*x2*y[1]*y[3]*y[5]*y[40]*y[46]*y[47]*y[48]*y[49]*y[50]-x1*x2*y[1]*y[3\
]*y[6]*y[40]*y[46]*y[47]*y[48]*y[49]*y[50]+x0*y[1]*y[2]*y[4]*y[13]*y[28]*y[\
29]*y[51]*y[52]*y[53]+2.*x2*y[1]*y[2]*y[3]*y[4]*y[13]*y[28]*y[29]*y[51]*y[5\
2]*y[53]+x2*y[1]*y[2]*y[3]*y[4]*y[40]*y[46]*y[47]*y[51]*y[52]*y[53]),2)+pow\
(x0*x1*y[1]*y[2]+x0*x2*y[1]*y[2]+x1*x2*y[1]*y[2]*y[3]+x0*y[1]*y[2]*y[4]+x2*\
y[1]*y[2]*y[3]*y[4]+x0*x1*y[1]*y[5]+x0*x2*y[1]*y[5]+x1*x2*y[1]*y[3]*y[5]-x0\
*x1*y[1]*y[6]-x0*x2*y[1]*y[6]-x1*x2*y[1]*y[3]*y[6]+y[7]+y[10]+y[14]+y[16]+y\
[18]+y[35]+lambda*lambda*(-(x0*x1*y[1]*y[2]*y[13]*y[28]*y[29]*y[31]*y[37]*y\
[38])-2.*x0*x1*x2*y[1]*y[2]*y[13]*y[28]*y[29]*y[31]*y[37]*y[38]-2.*x1*x2*y[\
1]*y[2]*y[3]*y[13]*y[28]*y[29]*y[31]*y[37]*y[38]-2.*x0*y[1]*y[2]*y[4]*y[13]\
*y[28]*y[29]*y[31]*y[37]*y[38]-4.*x2*y[1]*y[2]*y[3]*y[4]*y[13]*y[28]*y[29]*\
y[31]*y[37]*y[38]-x0*x1*y[1]*y[5]*y[13]*y[28]*y[29]*y[31]*y[37]*y[38]-2.*x1\
*x2*y[1]*y[3]*y[5]*y[13]*y[28]*y[29]*y[31]*y[37]*y[38]+x0*x1*y[1]*y[6]*y[13\
]*y[28]*y[29]*y[31]*y[37]*y[38]+2.*x1*x2*y[1]*y[3]*y[6]*y[13]*y[28]*y[29]*y\
[31]*y[37]*y[38]-x0*x2*y[1]*y[2]*y[13]*y[28]*y[29]*y[40]*y[46]*y[47]-2.*x0*\
x1*x2*y[1]*y[2]*y[13]*y[28]*y[29]*y[40]*y[46]*y[47]-2.*x1*x2*y[1]*y[2]*y[3]\
*y[13]*y[28]*y[29]*y[40]*y[46]*y[47]-2.*x2*y[1]*y[2]*y[3]*y[4]*y[13]*y[28]*\
y[29]*y[40]*y[46]*y[47]-x0*x2*y[1]*y[5]*y[13]*y[28]*y[29]*y[40]*y[46]*y[47]\
-2.*x2*y[1]*y[3]*y[5]*y[13]*y[28]*y[29]*y[40]*y[46]*y[47]-2.*x1*x2*y[1]*y[3\
]*y[5]*y[13]*y[28]*y[29]*y[40]*y[46]*y[47]+x0*x2*y[1]*y[6]*y[13]*y[28]*y[29\
]*y[40]*y[46]*y[47]+2.*x1*x2*y[1]*y[3]*y[6]*y[13]*y[28]*y[29]*y[40]*y[46]*y\
[47]-2.*x0*x1*x2*y[1]*y[2]*y[31]*y[37]*y[38]*y[40]*y[46]*y[47]-x1*x2*y[1]*y\
[2]*y[3]*y[31]*y[37]*y[38]*y[40]*y[46]*y[47]-2.*x2*y[1]*y[2]*y[3]*y[4]*y[31\
]*y[37]*y[38]*y[40]*y[46]*y[47]-x1*x2*y[1]*y[3]*y[5]*y[31]*y[37]*y[38]*y[40\
]*y[46]*y[47]+x1*x2*y[1]*y[3]*y[6]*y[31]*y[37]*y[38]*y[40]*y[46]*y[47]-x1*x\
2*y[1]*y[2]*y[3]*y[48]*y[49]*y[50]-x2*y[1]*y[2]*y[3]*y[4]*y[48]*y[49]*y[50]\
-x2*y[1]*y[3]*y[5]*y[48]*y[49]*y[50]-x1*x2*y[1]*y[3]*y[5]*y[48]*y[49]*y[50]\
+x1*x2*y[1]*y[3]*y[6]*y[48]*y[49]*y[50]-x0*y[1]*y[2]*y[4]*y[51]*y[52]*y[53]\
-x2*y[1]*y[2]*y[3]*y[4]*y[51]*y[52]*y[53])+pow(lambda,4)*(x1*x2*y[1]*y[2]*y\
[3]*y[31]*y[37]*y[38]*y[40]*y[46]*y[47]*y[48]*y[49]*y[50]+2.*x2*y[1]*y[2]*y\
[3]*y[4]*y[31]*y[37]*y[38]*y[40]*y[46]*y[47]*y[48]*y[49]*y[50]+x1*x2*y[1]*y\
[3]*y[5]*y[31]*y[37]*y[38]*y[40]*y[46]*y[47]*y[48]*y[49]*y[50]-x1*x2*y[1]*y\
[3]*y[6]*y[31]*y[37]*y[38]*y[40]*y[46]*y[47]*y[48]*y[49]*y[50]+2.*x2*y[1]*y\
[2]*y[3]*y[4]*y[13]*y[28]*y[29]*y[40]*y[46]*y[47]*y[51]*y[52]*y[53]+x2*y[1]\
*y[2]*y[3]*y[4]*y[48]*y[49]*y[50]*y[51]*y[52]*y[53]),2);
return (FOUT);
}
double Ps10(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
double y[54];
double FOUT;
y[1]=1./bi;
y[2]=em[0];
y[3]=x0*x0;
y[4]=x1*x1;
y[5]=em[1];
y[6]=esx[0];
y[7]=y[1]*y[2];
y[8]=x0*y[1]*y[2];
y[9]=2.*x0*x1*y[1]*y[2];
y[10]=x0*y[1]*y[5];
y[11]=-(x0*y[1]*y[6]);
y[12]=-x0;
y[13]=1.+y[12];
y[14]=x1*y[1]*y[2];
y[15]=y[1]*y[2]*y[4];
y[16]=x2*y[1]*y[2];
y[17]=2.*x1*x2*y[1]*y[2];
y[18]=2.*x0*x1*x2*y[1]*y[2];
y[19]=2.*x0*x2*y[1]*y[2]*y[4];
y[20]=y[1]*y[5];
y[21]=x1*y[1]*y[5];
y[22]=x2*y[1]*y[5];
y[23]=2.*x0*x2*y[1]*y[5];
y[24]=2.*x0*x1*x2*y[1]*y[5];
y[25]=-(x1*y[1]*y[6]);
y[26]=-(x2*y[1]*y[6]);
y[27]=-2.*x0*x1*x2*y[1]*y[6];
y[28]=y[14]+y[15]+y[16]+y[17]+y[18]+y[19]+y[20]+y[21]+y[22]+y[23]+y[24]+y[25\
]+y[26]+y[27];
y[29]=lrs[0];
y[30]=-x1;
y[31]=1.+y[30];
y[32]=2.*x0*x2*y[1]*y[2];
y[33]=x2*y[1]*y[2]*y[3];
y[34]=2.*x1*x2*y[1]*y[2]*y[3];
y[35]=x2*y[1]*y[3]*y[5];
y[36]=-(x2*y[1]*y[3]*y[6]);
y[37]=y[7]+y[8]+y[9]+y[10]+y[11]+y[32]+y[33]+y[34]+y[35]+y[36];
y[38]=lrs[1];
y[39]=-x2;
y[40]=1.+y[39];
y[41]=x1*y[1]*y[2]*y[3];
y[42]=y[1]*y[2]*y[3]*y[4];
y[43]=y[1]*y[3]*y[5];
y[44]=x1*y[1]*y[3]*y[5];
y[45]=-(x1*y[1]*y[3]*y[6]);
y[46]=y[7]+y[8]+y[9]+y[10]+y[11]+y[41]+y[42]+y[43]+y[44]+y[45];
y[47]=lrs[2];
y[48]=pow(y[13],2);
y[49]=pow(y[28],2);
y[50]=pow(y[29],2);
y[51]=pow(y[31],2);
y[52]=pow(y[37],2);
y[53]=pow(y[38],2);
FOUT=lambda*(-(x0*x1*y[1]*y[2]*y[13]*y[28]*y[29])-x0*x2*y[1]*y[2]*y[13]*y[28\
]*y[29]-2.*x0*x1*x2*y[1]*y[2]*y[13]*y[28]*y[29]-2.*x1*x2*y[1]*y[2]*y[3]*y[1\
3]*y[28]*y[29]-x0*y[1]*y[2]*y[4]*y[13]*y[28]*y[29]-2.*x2*y[1]*y[2]*y[3]*y[4\
]*y[13]*y[28]*y[29]-x0*y[1]*y[5]*y[13]*y[28]*y[29]-x0*x1*y[1]*y[5]*y[13]*y[\
28]*y[29]-x0*x2*y[1]*y[5]*y[13]*y[28]*y[29]-2.*x2*y[1]*y[3]*y[5]*y[13]*y[28\
]*y[29]-2.*x1*x2*y[1]*y[3]*y[5]*y[13]*y[28]*y[29]+x0*x1*y[1]*y[6]*y[13]*y[2\
8]*y[29]+x0*x2*y[1]*y[6]*y[13]*y[28]*y[29]+2.*x1*x2*y[1]*y[3]*y[6]*y[13]*y[\
28]*y[29]-x1*y[1]*y[2]*y[31]*y[37]*y[38]-x0*x1*y[1]*y[2]*y[31]*y[37]*y[38]-\
2.*x0*x1*x2*y[1]*y[2]*y[31]*y[37]*y[38]-x1*x2*y[1]*y[2]*y[3]*y[31]*y[37]*y[\
38]-2.*x0*y[1]*y[2]*y[4]*y[31]*y[37]*y[38]-2.*x2*y[1]*y[2]*y[3]*y[4]*y[31]*\
y[37]*y[38]-x0*x1*y[1]*y[5]*y[31]*y[37]*y[38]-x1*x2*y[1]*y[3]*y[5]*y[31]*y[\
37]*y[38]+x0*x1*y[1]*y[6]*y[31]*y[37]*y[38]+x1*x2*y[1]*y[3]*y[6]*y[31]*y[37\
]*y[38]-x2*y[1]*y[2]*y[40]*y[46]*y[47]-x0*x2*y[1]*y[2]*y[40]*y[46]*y[47]-2.\
*x0*x1*x2*y[1]*y[2]*y[40]*y[46]*y[47]-x1*x2*y[1]*y[2]*y[3]*y[40]*y[46]*y[47\
]-x2*y[1]*y[2]*y[3]*y[4]*y[40]*y[46]*y[47]-x0*x2*y[1]*y[5]*y[40]*y[46]*y[47\
]-x2*y[1]*y[3]*y[5]*y[40]*y[46]*y[47]-x1*x2*y[1]*y[3]*y[5]*y[40]*y[46]*y[47\
]+x0*x2*y[1]*y[6]*y[40]*y[46]*y[47]+x1*x2*y[1]*y[3]*y[6]*y[40]*y[46]*y[47])\
-x2*pow(lambda,5)*y[1]*y[2]*y[3]*y[4]*y[40]*y[46]*y[47]*y[48]*y[49]*y[50]*y\
[51]*y[52]*y[53]+pow(lambda,3)*(2.*x0*x1*x2*y[1]*y[2]*y[13]*y[28]*y[29]*y[3\
1]*y[37]*y[38]*y[40]*y[46]*y[47]+2.*x1*x2*y[1]*y[2]*y[3]*y[13]*y[28]*y[29]*\
y[31]*y[37]*y[38]*y[40]*y[46]*y[47]+4.*x2*y[1]*y[2]*y[3]*y[4]*y[13]*y[28]*y\
[29]*y[31]*y[37]*y[38]*y[40]*y[46]*y[47]+2.*x1*x2*y[1]*y[3]*y[5]*y[13]*y[28\
]*y[29]*y[31]*y[37]*y[38]*y[40]*y[46]*y[47]-2.*x1*x2*y[1]*y[3]*y[6]*y[13]*y\
[28]*y[29]*y[31]*y[37]*y[38]*y[40]*y[46]*y[47]+x1*x2*y[1]*y[2]*y[3]*y[31]*y\
[37]*y[38]*y[48]*y[49]*y[50]+2.*x2*y[1]*y[2]*y[3]*y[4]*y[31]*y[37]*y[38]*y[\
48]*y[49]*y[50]+x1*x2*y[1]*y[3]*y[5]*y[31]*y[37]*y[38]*y[48]*y[49]*y[50]-x1\
*x2*y[1]*y[3]*y[6]*y[31]*y[37]*y[38]*y[48]*y[49]*y[50]+x1*x2*y[1]*y[2]*y[3]\
*y[40]*y[46]*y[47]*y[48]*y[49]*y[50]+x2*y[1]*y[2]*y[3]*y[4]*y[40]*y[46]*y[4\
7]*y[48]*y[49]*y[50]+x2*y[1]*y[3]*y[5]*y[40]*y[46]*y[47]*y[48]*y[49]*y[50]+\
x1*x2*y[1]*y[3]*y[5]*y[40]*y[46]*y[47]*y[48]*y[49]*y[50]-x1*x2*y[1]*y[3]*y[\
6]*y[40]*y[46]*y[47]*y[48]*y[49]*y[50]+x0*y[1]*y[2]*y[4]*y[13]*y[28]*y[29]*\
y[51]*y[52]*y[53]+2.*x2*y[1]*y[2]*y[3]*y[4]*y[13]*y[28]*y[29]*y[51]*y[52]*y\
[53]+x2*y[1]*y[2]*y[3]*y[4]*y[40]*y[46]*y[47]*y[51]*y[52]*y[53]);
return (FOUT);
}
double Pa10(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
double y[54];
double FOUT;
y[1]=1./bi;
y[2]=em[0];
y[3]=x0*x0;
y[4]=x1*x1;
y[5]=em[1];
y[6]=esx[0];
y[7]=y[1]*y[2];
y[8]=x0*y[1]*y[2];
y[9]=2.*x0*x1*y[1]*y[2];
y[10]=x0*y[1]*y[5];
y[11]=-(x0*y[1]*y[6]);
y[12]=-x0;
y[13]=1.+y[12];
y[14]=x1*y[1]*y[2];
y[15]=y[1]*y[2]*y[4];
y[16]=x2*y[1]*y[2];
y[17]=2.*x1*x2*y[1]*y[2];
y[18]=2.*x0*x1*x2*y[1]*y[2];
y[19]=2.*x0*x2*y[1]*y[2]*y[4];
y[20]=y[1]*y[5];
y[21]=x1*y[1]*y[5];
y[22]=x2*y[1]*y[5];
y[23]=2.*x0*x2*y[1]*y[5];
y[24]=2.*x0*x1*x2*y[1]*y[5];
y[25]=-(x1*y[1]*y[6]);
y[26]=-(x2*y[1]*y[6]);
y[27]=-2.*x0*x1*x2*y[1]*y[6];
y[28]=y[14]+y[15]+y[16]+y[17]+y[18]+y[19]+y[20]+y[21]+y[22]+y[23]+y[24]+y[25\
]+y[26]+y[27];
y[29]=lrs[0];
y[30]=-x1;
y[31]=1.+y[30];
y[32]=2.*x0*x2*y[1]*y[2];
y[33]=x2*y[1]*y[2]*y[3];
y[34]=2.*x1*x2*y[1]*y[2]*y[3];
y[35]=x2*y[1]*y[3]*y[5];
y[36]=-(x2*y[1]*y[3]*y[6]);
y[37]=y[7]+y[8]+y[9]+y[10]+y[11]+y[32]+y[33]+y[34]+y[35]+y[36];
y[38]=lrs[1];
y[39]=-x2;
y[40]=1.+y[39];
y[41]=x1*y[1]*y[2]*y[3];
y[42]=y[1]*y[2]*y[3]*y[4];
y[43]=y[1]*y[3]*y[5];
y[44]=x1*y[1]*y[3]*y[5];
y[45]=-(x1*y[1]*y[3]*y[6]);
y[46]=y[7]+y[8]+y[9]+y[10]+y[11]+y[41]+y[42]+y[43]+y[44]+y[45];
y[47]=lrs[2];
y[48]=pow(y[13],2);
y[49]=pow(y[28],2);
y[50]=pow(y[29],2);
y[51]=pow(y[31],2);
y[52]=pow(y[37],2);
y[53]=pow(y[38],2);
FOUT=(lambda*(-(x0*x1*y[1]*y[2]*y[13]*y[28]*y[29])-x0*x2*y[1]*y[2]*y[13]*y[2\
8]*y[29]-2.*x0*x1*x2*y[1]*y[2]*y[13]*y[28]*y[29]-2.*x1*x2*y[1]*y[2]*y[3]*y[\
13]*y[28]*y[29]-x0*y[1]*y[2]*y[4]*y[13]*y[28]*y[29]-2.*x2*y[1]*y[2]*y[3]*y[\
4]*y[13]*y[28]*y[29]-x0*y[1]*y[5]*y[13]*y[28]*y[29]-x0*x1*y[1]*y[5]*y[13]*y\
[28]*y[29]-x0*x2*y[1]*y[5]*y[13]*y[28]*y[29]-2.*x2*y[1]*y[3]*y[5]*y[13]*y[2\
8]*y[29]-2.*x1*x2*y[1]*y[3]*y[5]*y[13]*y[28]*y[29]+x0*x1*y[1]*y[6]*y[13]*y[\
28]*y[29]+x0*x2*y[1]*y[6]*y[13]*y[28]*y[29]+2.*x1*x2*y[1]*y[3]*y[6]*y[13]*y\
[28]*y[29]-x1*y[1]*y[2]*y[31]*y[37]*y[38]-x0*x1*y[1]*y[2]*y[31]*y[37]*y[38]\
-2.*x0*x1*x2*y[1]*y[2]*y[31]*y[37]*y[38]-x1*x2*y[1]*y[2]*y[3]*y[31]*y[37]*y\
[38]-2.*x0*y[1]*y[2]*y[4]*y[31]*y[37]*y[38]-2.*x2*y[1]*y[2]*y[3]*y[4]*y[31]\
*y[37]*y[38]-x0*x1*y[1]*y[5]*y[31]*y[37]*y[38]-x1*x2*y[1]*y[3]*y[5]*y[31]*y\
[37]*y[38]+x0*x1*y[1]*y[6]*y[31]*y[37]*y[38]+x1*x2*y[1]*y[3]*y[6]*y[31]*y[3\
7]*y[38]-x2*y[1]*y[2]*y[40]*y[46]*y[47]-x0*x2*y[1]*y[2]*y[40]*y[46]*y[47]-2\
.*x0*x1*x2*y[1]*y[2]*y[40]*y[46]*y[47]-x1*x2*y[1]*y[2]*y[3]*y[40]*y[46]*y[4\
7]-x2*y[1]*y[2]*y[3]*y[4]*y[40]*y[46]*y[47]-x0*x2*y[1]*y[5]*y[40]*y[46]*y[4\
7]-x2*y[1]*y[3]*y[5]*y[40]*y[46]*y[47]-x1*x2*y[1]*y[3]*y[5]*y[40]*y[46]*y[4\
7]+x0*x2*y[1]*y[6]*y[40]*y[46]*y[47]+x1*x2*y[1]*y[3]*y[6]*y[40]*y[46]*y[47]\
)-x2*pow(lambda,5)*y[1]*y[2]*y[3]*y[4]*y[40]*y[46]*y[47]*y[48]*y[49]*y[50]*\
y[51]*y[52]*y[53]+pow(lambda,3)*(2.*x0*x1*x2*y[1]*y[2]*y[13]*y[28]*y[29]*y[\
31]*y[37]*y[38]*y[40]*y[46]*y[47]+2.*x1*x2*y[1]*y[2]*y[3]*y[13]*y[28]*y[29]\
*y[31]*y[37]*y[38]*y[40]*y[46]*y[47]+4.*x2*y[1]*y[2]*y[3]*y[4]*y[13]*y[28]*\
y[29]*y[31]*y[37]*y[38]*y[40]*y[46]*y[47]+2.*x1*x2*y[1]*y[3]*y[5]*y[13]*y[2\
8]*y[29]*y[31]*y[37]*y[38]*y[40]*y[46]*y[47]-2.*x1*x2*y[1]*y[3]*y[6]*y[13]*\
y[28]*y[29]*y[31]*y[37]*y[38]*y[40]*y[46]*y[47]+x1*x2*y[1]*y[2]*y[3]*y[31]*\
y[37]*y[38]*y[48]*y[49]*y[50]+2.*x2*y[1]*y[2]*y[3]*y[4]*y[31]*y[37]*y[38]*y\
[48]*y[49]*y[50]+x1*x2*y[1]*y[3]*y[5]*y[31]*y[37]*y[38]*y[48]*y[49]*y[50]-x\
1*x2*y[1]*y[3]*y[6]*y[31]*y[37]*y[38]*y[48]*y[49]*y[50]+x1*x2*y[1]*y[2]*y[3\
]*y[40]*y[46]*y[47]*y[48]*y[49]*y[50]+x2*y[1]*y[2]*y[3]*y[4]*y[40]*y[46]*y[\
47]*y[48]*y[49]*y[50]+x2*y[1]*y[3]*y[5]*y[40]*y[46]*y[47]*y[48]*y[49]*y[50]\
+x1*x2*y[1]*y[3]*y[5]*y[40]*y[46]*y[47]*y[48]*y[49]*y[50]-x1*x2*y[1]*y[3]*y\
[6]*y[40]*y[46]*y[47]*y[48]*y[49]*y[50]+x0*y[1]*y[2]*y[4]*y[13]*y[28]*y[29]\
*y[51]*y[52]*y[53]+2.*x2*y[1]*y[2]*y[3]*y[4]*y[13]*y[28]*y[29]*y[51]*y[52]*\
y[53]+x2*y[1]*y[2]*y[3]*y[4]*y[40]*y[46]*y[47]*y[51]*y[52]*y[53]))/(lambda*\
(x0*x1*y[1]*y[2]+x0*x2*y[1]*y[2]+x1*x2*y[1]*y[2]*y[3]+x0*y[1]*y[2]*y[4]+x2*\
y[1]*y[2]*y[3]*y[4]+x0*x1*y[1]*y[5]+x0*x2*y[1]*y[5]+x1*x2*y[1]*y[3]*y[5]-x0\
*x1*y[1]*y[6]-x0*x2*y[1]*y[6]-x1*x2*y[1]*y[3]*y[6]+y[7]+y[10]+y[14]+y[16]+y\
[18]+y[35]+lambda*lambda*(-(x0*x1*y[1]*y[2]*y[13]*y[28]*y[29]*y[31]*y[37]*y\
[38])-2.*x0*x1*x2*y[1]*y[2]*y[13]*y[28]*y[29]*y[31]*y[37]*y[38]-2.*x1*x2*y[\
1]*y[2]*y[3]*y[13]*y[28]*y[29]*y[31]*y[37]*y[38]-2.*x0*y[1]*y[2]*y[4]*y[13]\
*y[28]*y[29]*y[31]*y[37]*y[38]-4.*x2*y[1]*y[2]*y[3]*y[4]*y[13]*y[28]*y[29]*\
y[31]*y[37]*y[38]-x0*x1*y[1]*y[5]*y[13]*y[28]*y[29]*y[31]*y[37]*y[38]-2.*x1\
*x2*y[1]*y[3]*y[5]*y[13]*y[28]*y[29]*y[31]*y[37]*y[38]+x0*x1*y[1]*y[6]*y[13\
]*y[28]*y[29]*y[31]*y[37]*y[38]+2.*x1*x2*y[1]*y[3]*y[6]*y[13]*y[28]*y[29]*y\
[31]*y[37]*y[38]-x0*x2*y[1]*y[2]*y[13]*y[28]*y[29]*y[40]*y[46]*y[47]-2.*x0*\
x1*x2*y[1]*y[2]*y[13]*y[28]*y[29]*y[40]*y[46]*y[47]-2.*x1*x2*y[1]*y[2]*y[3]\
*y[13]*y[28]*y[29]*y[40]*y[46]*y[47]-2.*x2*y[1]*y[2]*y[3]*y[4]*y[13]*y[28]*\
y[29]*y[40]*y[46]*y[47]-x0*x2*y[1]*y[5]*y[13]*y[28]*y[29]*y[40]*y[46]*y[47]\
-2.*x2*y[1]*y[3]*y[5]*y[13]*y[28]*y[29]*y[40]*y[46]*y[47]-2.*x1*x2*y[1]*y[3\
]*y[5]*y[13]*y[28]*y[29]*y[40]*y[46]*y[47]+x0*x2*y[1]*y[6]*y[13]*y[28]*y[29\
]*y[40]*y[46]*y[47]+2.*x1*x2*y[1]*y[3]*y[6]*y[13]*y[28]*y[29]*y[40]*y[46]*y\
[47]-2.*x0*x1*x2*y[1]*y[2]*y[31]*y[37]*y[38]*y[40]*y[46]*y[47]-x1*x2*y[1]*y\
[2]*y[3]*y[31]*y[37]*y[38]*y[40]*y[46]*y[47]-2.*x2*y[1]*y[2]*y[3]*y[4]*y[31\
]*y[37]*y[38]*y[40]*y[46]*y[47]-x1*x2*y[1]*y[3]*y[5]*y[31]*y[37]*y[38]*y[40\
]*y[46]*y[47]+x1*x2*y[1]*y[3]*y[6]*y[31]*y[37]*y[38]*y[40]*y[46]*y[47]-x1*x\
2*y[1]*y[2]*y[3]*y[48]*y[49]*y[50]-x2*y[1]*y[2]*y[3]*y[4]*y[48]*y[49]*y[50]\
-x2*y[1]*y[3]*y[5]*y[48]*y[49]*y[50]-x1*x2*y[1]*y[3]*y[5]*y[48]*y[49]*y[50]\
+x1*x2*y[1]*y[3]*y[6]*y[48]*y[49]*y[50]-x0*y[1]*y[2]*y[4]*y[51]*y[52]*y[53]\
-x2*y[1]*y[2]*y[3]*y[4]*y[51]*y[52]*y[53])+pow(lambda,4)*(x1*x2*y[1]*y[2]*y\
[3]*y[31]*y[37]*y[38]*y[40]*y[46]*y[47]*y[48]*y[49]*y[50]+2.*x2*y[1]*y[2]*y\
[3]*y[4]*y[31]*y[37]*y[38]*y[40]*y[46]*y[47]*y[48]*y[49]*y[50]+x1*x2*y[1]*y\
[3]*y[5]*y[31]*y[37]*y[38]*y[40]*y[46]*y[47]*y[48]*y[49]*y[50]-x1*x2*y[1]*y\
[3]*y[6]*y[31]*y[37]*y[38]*y[40]*y[46]*y[47]*y[48]*y[49]*y[50]+2.*x2*y[1]*y\
[2]*y[3]*y[4]*y[13]*y[28]*y[29]*y[40]*y[46]*y[47]*y[51]*y[52]*y[53]+x2*y[1]\
*y[2]*y[3]*y[4]*y[48]*y[49]*y[50]*y[51]*y[52]*y[53])));
return (FOUT);
}
double Pt10t1(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
double y[6];
double FOUT;
y[1]=1./bi;
y[2]=em[0];
y[3]=x1*x1;
y[4]=em[1];
y[5]=esx[0];
FOUT=(1.-x0)*x0*(x1*y[1]*y[2]+x2*y[1]*y[2]+2.*x1*x2*y[1]*y[2]+2.*x0*x1*x2*y[\
1]*y[2]+y[1]*y[2]*y[3]+2.*x0*x2*y[1]*y[2]*y[3]+y[1]*y[4]+x1*y[1]*y[4]+x2*y[\
1]*y[4]+2.*x0*x2*y[1]*y[4]+2.*x0*x1*x2*y[1]*y[4]-x1*y[1]*y[5]-x2*y[1]*y[5]-\
2.*x0*x1*x2*y[1]*y[5]);
return (FOUT);
}
double Pt10t2(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
double y[6];
double FOUT;
y[1]=1./bi;
y[2]=em[0];
y[3]=x0*x0;
y[4]=em[1];
y[5]=esx[0];
FOUT=(1.-x1)*x1*(y[1]*y[2]+x0*y[1]*y[2]+2.*x0*x1*y[1]*y[2]+2.*x0*x2*y[1]*y[2\
]+x2*y[1]*y[2]*y[3]+2.*x1*x2*y[1]*y[2]*y[3]+x0*y[1]*y[4]+x2*y[1]*y[3]*y[4]-\
x0*y[1]*y[5]-x2*y[1]*y[3]*y[5]);
return (FOUT);
}
double Pt10t3(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
double y[6];
double FOUT;
y[1]=1./bi;
y[2]=em[0];
y[3]=x0*x0;
y[4]=em[1];
y[5]=esx[0];
FOUT=(1.-x2)*x2*(y[1]*y[2]+x0*y[1]*y[2]+2.*x0*x1*y[1]*y[2]+x1*y[1]*y[2]*y[3]\
+x1*x1*y[1]*y[2]*y[3]+x0*y[1]*y[4]+y[1]*y[3]*y[4]+x1*y[1]*y[3]*y[4]-x0*y[1]\
*y[5]-x1*y[1]*y[3]*y[5]);
return (FOUT);
}
