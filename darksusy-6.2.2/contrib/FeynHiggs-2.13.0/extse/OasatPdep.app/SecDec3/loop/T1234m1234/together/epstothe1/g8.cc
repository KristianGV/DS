#include "intfile.hh"

double Pr8(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
double y[77];
double FOUT;
y[1]=1./bi;
y[2]=em[0];
y[3]=x2*x2;
y[4]=em[1];
y[5]=em[2];
y[6]=em[3];
y[7]=x1*x1;
y[8]=esx[0];
y[9]=-x0;
y[10]=1.+y[9];
y[11]=x2*y[1]*y[2];
y[12]=x1*x2*y[1]*y[2];
y[13]=2.*x0*x1*x2*y[1]*y[2];
y[14]=y[1]*y[2]*y[3];
y[15]=2.*x0*x1*y[1]*y[2]*y[3];
y[16]=y[1]*y[4];
y[17]=x1*y[1]*y[4];
y[18]=2.*x0*x1*y[1]*y[4];
y[19]=x2*y[1]*y[4];
y[20]=2.*x0*x1*x2*y[1]*y[4];
y[21]=x1*y[1]*y[5];
y[22]=x1*x2*y[1]*y[5];
y[23]=x1*y[1]*y[6];
y[24]=y[1]*y[6]*y[7];
y[25]=2.*x0*y[1]*y[6]*y[7];
y[26]=x1*x2*y[1]*y[6];
y[27]=2.*x0*x2*y[1]*y[6]*y[7];
y[28]=-(x1*y[1]*y[8]);
y[29]=-(x2*y[1]*y[8]);
y[30]=-2.*x0*x1*x2*y[1]*y[8];
y[31]=y[11]+y[12]+y[13]+y[14]+y[15]+y[16]+y[17]+y[18]+y[19]+y[20]+y[21]+y[22\
]+y[23]+y[24]+y[25]+y[26]+y[27]+y[28]+y[29]+y[30];
y[32]=lrs[0];
y[33]=x0*x0;
y[34]=-x1;
y[35]=1.+y[34];
y[36]=x0*x2*y[1]*y[2];
y[37]=x2*y[1]*y[2]*y[33];
y[38]=y[1]*y[2]*y[3]*y[33];
y[39]=x0*y[1]*y[4];
y[40]=y[1]*y[4]*y[33];
y[41]=x2*y[1]*y[4]*y[33];
y[42]=y[1]*y[5];
y[43]=x0*y[1]*y[5];
y[44]=x0*x2*y[1]*y[5];
y[45]=x0*y[1]*y[6];
y[46]=2.*x0*x1*y[1]*y[6];
y[47]=2.*x1*y[1]*y[6]*y[33];
y[48]=x0*x2*y[1]*y[6];
y[49]=2.*x1*x2*y[1]*y[6]*y[33];
y[50]=-(x0*y[1]*y[8]);
y[51]=-(x2*y[1]*y[8]*y[33]);
y[52]=y[36]+y[37]+y[38]+y[39]+y[40]+y[41]+y[42]+y[43]+y[44]+y[45]+y[46]+y[47\
]+y[48]+y[49]+y[50]+y[51];
y[53]=lrs[1];
y[54]=-x2;
y[55]=1.+y[54];
y[56]=x0*y[1]*y[2];
y[57]=x0*x1*y[1]*y[2];
y[58]=x1*y[1]*y[2]*y[33];
y[59]=2.*x0*x2*y[1]*y[2];
y[60]=2.*x1*x2*y[1]*y[2]*y[33];
y[61]=x1*y[1]*y[4]*y[33];
y[62]=x0*x1*y[1]*y[5];
y[63]=x0*x1*y[1]*y[6];
y[64]=y[1]*y[6]*y[7]*y[33];
y[65]=-(x1*y[1]*y[8]*y[33]);
y[66]=y[39]+y[42]+y[50]+y[56]+y[57]+y[58]+y[59]+y[60]+y[61]+y[62]+y[63]+y[64\
]+y[65];
y[67]=lrs[2];
y[68]=pow(y[10],2);
y[69]=pow(y[31],2);
y[70]=pow(y[32],2);
y[71]=pow(y[35],2);
y[72]=pow(y[52],2);
y[73]=pow(y[53],2);
y[74]=pow(y[55],2);
y[75]=pow(y[66],2);
y[76]=pow(y[67],2);
FOUT=(x0*x1*x2*y[1]*y[2]*y[10]*y[31]*y[32]*y[35]*y[52]*y[53]*y[55]*y[66]*y[6\
7]+x0*x1*x2*y[1]*y[5]*y[10]*y[31]*y[32]*y[35]*y[52]*y[53]*y[55]*y[66]*y[67]\
+x0*x1*x2*y[1]*y[6]*y[10]*y[31]*y[32]*y[35]*y[52]*y[53]*y[55]*y[66]*y[67]+2\
.*x1*x2*y[1]*y[2]*y[10]*y[31]*y[32]*y[33]*y[35]*y[52]*y[53]*y[55]*y[66]*y[6\
7]+4.*x1*y[1]*y[2]*y[3]*y[10]*y[31]*y[32]*y[33]*y[35]*y[52]*y[53]*y[55]*y[6\
6]*y[67]+2.*x1*x2*y[1]*y[4]*y[10]*y[31]*y[32]*y[33]*y[35]*y[52]*y[53]*y[55]\
*y[66]*y[67]+4.*x2*y[1]*y[6]*y[7]*y[10]*y[31]*y[32]*y[33]*y[35]*y[52]*y[53]\
*y[55]*y[66]*y[67]-2.*x1*x2*y[1]*y[8]*y[10]*y[31]*y[32]*y[33]*y[35]*y[52]*y\
[53]*y[55]*y[66]*y[67]+x1*x2*y[1]*y[2]*y[33]*y[35]*y[52]*y[53]*y[68]*y[69]*\
y[70]+x1*y[1]*y[2]*y[3]*y[33]*y[35]*y[52]*y[53]*y[68]*y[69]*y[70]+x1*y[1]*y\
[4]*y[33]*y[35]*y[52]*y[53]*y[68]*y[69]*y[70]+x1*x2*y[1]*y[4]*y[33]*y[35]*y\
[52]*y[53]*y[68]*y[69]*y[70]+2.*y[1]*y[6]*y[7]*y[33]*y[35]*y[52]*y[53]*y[68\
]*y[69]*y[70]+2.*x2*y[1]*y[6]*y[7]*y[33]*y[35]*y[52]*y[53]*y[68]*y[69]*y[70\
]-x1*x2*y[1]*y[8]*y[33]*y[35]*y[52]*y[53]*y[68]*y[69]*y[70]+x1*x2*y[1]*y[2]\
*y[33]*y[55]*y[66]*y[67]*y[68]*y[69]*y[70]+2.*x1*y[1]*y[2]*y[3]*y[33]*y[55]\
*y[66]*y[67]*y[68]*y[69]*y[70]+x1*x2*y[1]*y[4]*y[33]*y[55]*y[66]*y[67]*y[68\
]*y[69]*y[70]+x2*y[1]*y[6]*y[7]*y[33]*y[55]*y[66]*y[67]*y[68]*y[69]*y[70]-x\
1*x2*y[1]*y[8]*y[33]*y[55]*y[66]*y[67]*y[68]*y[69]*y[70]+x0*y[1]*y[6]*y[7]*\
y[10]*y[31]*y[32]*y[71]*y[72]*y[73]+2.*y[1]*y[6]*y[7]*y[10]*y[31]*y[32]*y[3\
3]*y[71]*y[72]*y[73]+2.*x2*y[1]*y[6]*y[7]*y[10]*y[31]*y[32]*y[33]*y[71]*y[7\
2]*y[73]+x2*y[1]*y[6]*y[7]*y[33]*y[55]*y[66]*y[67]*y[71]*y[72]*y[73]+x0*y[1\
]*y[2]*y[3]*y[10]*y[31]*y[32]*y[74]*y[75]*y[76]+2.*x1*y[1]*y[2]*y[3]*y[10]*\
y[31]*y[32]*y[33]*y[74]*y[75]*y[76]+x1*y[1]*y[2]*y[3]*y[33]*y[35]*y[52]*y[5\
3]*y[74]*y[75]*y[76])/(-(x0*x2*y[1]*y[2]*y[10]*y[31]*y[32])-x0*x1*x2*y[1]*y\
[2]*y[10]*y[31]*y[32]-x0*y[1]*y[2]*y[3]*y[10]*y[31]*y[32]-x0*y[1]*y[4]*y[10\
]*y[31]*y[32]-x0*x1*y[1]*y[4]*y[10]*y[31]*y[32]-x0*x2*y[1]*y[4]*y[10]*y[31]\
*y[32]-x0*x1*y[1]*y[5]*y[10]*y[31]*y[32]-x0*x1*x2*y[1]*y[5]*y[10]*y[31]*y[3\
2]-x0*x1*y[1]*y[6]*y[10]*y[31]*y[32]-x0*x1*x2*y[1]*y[6]*y[10]*y[31]*y[32]-x\
0*y[1]*y[6]*y[7]*y[10]*y[31]*y[32]+x0*x1*y[1]*y[8]*y[10]*y[31]*y[32]+x0*x2*\
y[1]*y[8]*y[10]*y[31]*y[32]-2.*x1*x2*y[1]*y[2]*y[10]*y[31]*y[32]*y[33]-2.*x\
1*y[1]*y[2]*y[3]*y[10]*y[31]*y[32]*y[33]-2.*x1*y[1]*y[4]*y[10]*y[31]*y[32]*\
y[33]-2.*x1*x2*y[1]*y[4]*y[10]*y[31]*y[32]*y[33]-2.*y[1]*y[6]*y[7]*y[10]*y[\
31]*y[32]*y[33]-2.*x2*y[1]*y[6]*y[7]*y[10]*y[31]*y[32]*y[33]+2.*x1*x2*y[1]*\
y[8]*y[10]*y[31]*y[32]*y[33]-x0*x1*x2*y[1]*y[2]*y[35]*y[52]*y[53]-x0*x1*y[1\
]*y[4]*y[35]*y[52]*y[53]-x1*y[1]*y[5]*y[35]*y[52]*y[53]-x0*x1*y[1]*y[5]*y[3\
5]*y[52]*y[53]-x0*x1*x2*y[1]*y[5]*y[35]*y[52]*y[53]-x0*x1*y[1]*y[6]*y[35]*y\
[52]*y[53]-x0*x1*x2*y[1]*y[6]*y[35]*y[52]*y[53]-2.*x0*y[1]*y[6]*y[7]*y[35]*\
y[52]*y[53]+x0*x1*y[1]*y[8]*y[35]*y[52]*y[53]-x1*x2*y[1]*y[2]*y[33]*y[35]*y\
[52]*y[53]-x1*y[1]*y[2]*y[3]*y[33]*y[35]*y[52]*y[53]-x1*y[1]*y[4]*y[33]*y[3\
5]*y[52]*y[53]-x1*x2*y[1]*y[4]*y[33]*y[35]*y[52]*y[53]-2.*y[1]*y[6]*y[7]*y[\
33]*y[35]*y[52]*y[53]-2.*x2*y[1]*y[6]*y[7]*y[33]*y[35]*y[52]*y[53]+x1*x2*y[\
1]*y[8]*y[33]*y[35]*y[52]*y[53]-x0*x2*y[1]*y[2]*y[55]*y[66]*y[67]-x0*x1*x2*\
y[1]*y[2]*y[55]*y[66]*y[67]-2.*x0*y[1]*y[2]*y[3]*y[55]*y[66]*y[67]-x0*x2*y[\
1]*y[4]*y[55]*y[66]*y[67]-x2*y[1]*y[5]*y[55]*y[66]*y[67]-x0*x1*x2*y[1]*y[5]\
*y[55]*y[66]*y[67]-x0*x1*x2*y[1]*y[6]*y[55]*y[66]*y[67]+x0*x2*y[1]*y[8]*y[5\
5]*y[66]*y[67]-x1*x2*y[1]*y[2]*y[33]*y[55]*y[66]*y[67]-2.*x1*y[1]*y[2]*y[3]\
*y[33]*y[55]*y[66]*y[67]-x1*x2*y[1]*y[4]*y[33]*y[55]*y[66]*y[67]-x2*y[1]*y[\
6]*y[7]*y[33]*y[55]*y[66]*y[67]+x1*x2*y[1]*y[8]*y[33]*y[55]*y[66]*y[67]);
return (FOUT);
}
double Pm8(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
double y[77];
double FOUT;
y[1]=1./bi;
y[2]=em[0];
y[3]=x0*x0;
y[4]=x2*x2;
y[5]=em[1];
y[6]=em[2];
y[7]=em[3];
y[8]=x1*x1;
y[9]=esx[0];
y[10]=x1*y[1]*y[6];
y[11]=-x0;
y[12]=1.+y[11];
y[13]=pow(y[12],2);
y[14]=x2*y[1]*y[2];
y[15]=x1*x2*y[1]*y[2];
y[16]=2.*x0*x1*x2*y[1]*y[2];
y[17]=y[1]*y[2]*y[4];
y[18]=2.*x0*x1*y[1]*y[2]*y[4];
y[19]=y[1]*y[5];
y[20]=x1*y[1]*y[5];
y[21]=2.*x0*x1*y[1]*y[5];
y[22]=x2*y[1]*y[5];
y[23]=2.*x0*x1*x2*y[1]*y[5];
y[24]=x1*x2*y[1]*y[6];
y[25]=x1*y[1]*y[7];
y[26]=y[1]*y[7]*y[8];
y[27]=2.*x0*y[1]*y[7]*y[8];
y[28]=x1*x2*y[1]*y[7];
y[29]=2.*x0*x2*y[1]*y[7]*y[8];
y[30]=-(x1*y[1]*y[9]);
y[31]=-(x2*y[1]*y[9]);
y[32]=-2.*x0*x1*x2*y[1]*y[9];
y[33]=y[10]+y[14]+y[15]+y[16]+y[17]+y[18]+y[19]+y[20]+y[21]+y[22]+y[23]+y[24\
]+y[25]+y[26]+y[27]+y[28]+y[29]+y[30]+y[31]+y[32];
y[34]=pow(y[33],2);
y[35]=lrs[0];
y[36]=pow(y[35],2);
y[37]=x0*x2*y[1]*y[2];
y[38]=x0*y[1]*y[5];
y[39]=y[1]*y[6];
y[40]=-x1;
y[41]=1.+y[40];
y[42]=x2*y[1]*y[2]*y[3];
y[43]=y[1]*y[2]*y[3]*y[4];
y[44]=y[1]*y[3]*y[5];
y[45]=x2*y[1]*y[3]*y[5];
y[46]=x0*y[1]*y[6];
y[47]=x0*x2*y[1]*y[6];
y[48]=x0*y[1]*y[7];
y[49]=2.*x0*x1*y[1]*y[7];
y[50]=2.*x1*y[1]*y[3]*y[7];
y[51]=x0*x2*y[1]*y[7];
y[52]=2.*x1*x2*y[1]*y[3]*y[7];
y[53]=-(x0*y[1]*y[9]);
y[54]=-(x2*y[1]*y[3]*y[9]);
y[55]=y[37]+y[38]+y[39]+y[42]+y[43]+y[44]+y[45]+y[46]+y[47]+y[48]+y[49]+y[50\
]+y[51]+y[52]+y[53]+y[54];
y[56]=lrs[1];
y[57]=pow(y[41],2);
y[58]=pow(y[55],2);
y[59]=pow(y[56],2);
y[60]=x1*y[1]*y[3]*y[5];
y[61]=x0*x1*y[1]*y[6];
y[62]=x0*x1*y[1]*y[7];
y[63]=y[1]*y[3]*y[7]*y[8];
y[64]=-x2;
y[65]=1.+y[64];
y[66]=x0*y[1]*y[2];
y[67]=x0*x1*y[1]*y[2];
y[68]=x1*y[1]*y[2]*y[3];
y[69]=2.*x0*x2*y[1]*y[2];
y[70]=2.*x1*x2*y[1]*y[2]*y[3];
y[71]=-(x1*y[1]*y[3]*y[9]);
y[72]=y[38]+y[39]+y[53]+y[60]+y[61]+y[62]+y[63]+y[66]+y[67]+y[68]+y[69]+y[70\
]+y[71];
y[73]=lrs[2];
y[74]=pow(y[65],2);
y[75]=pow(y[72],2);
y[76]=pow(y[73],2);
FOUT=pow(x0*x1*x2*y[1]*y[2]+x1*x2*y[1]*y[2]*y[3]+x0*y[1]*y[2]*y[4]+x1*y[1]*y\
[2]*y[3]*y[4]+x0*x1*y[1]*y[5]+x0*x2*y[1]*y[5]+x1*x2*y[1]*y[3]*y[5]+x2*y[1]*\
y[6]+x0*x1*x2*y[1]*y[6]+x0*x1*x2*y[1]*y[7]+x0*y[1]*y[7]*y[8]+x2*y[1]*y[3]*y\
[7]*y[8]-x0*x1*y[1]*y[9]-x0*x2*y[1]*y[9]-x1*x2*y[1]*y[3]*y[9]+y[10]+y[37]+y\
[38]+y[39]+y[60]+y[61]+y[62]+y[63]+lambda*lambda*(-(x1*x2*y[1]*y[2]*y[3]*y[\
13]*y[34]*y[36])-x1*y[1]*y[2]*y[3]*y[4]*y[13]*y[34]*y[36]-x1*y[1]*y[3]*y[5]\
*y[13]*y[34]*y[36]-x1*x2*y[1]*y[3]*y[5]*y[13]*y[34]*y[36]-y[1]*y[3]*y[7]*y[\
8]*y[13]*y[34]*y[36]-x2*y[1]*y[3]*y[7]*y[8]*y[13]*y[34]*y[36]+x1*x2*y[1]*y[\
3]*y[9]*y[13]*y[34]*y[36]-x0*x1*x2*y[1]*y[2]*y[12]*y[33]*y[35]*y[41]*y[55]*\
y[56]-2.*x1*x2*y[1]*y[2]*y[3]*y[12]*y[33]*y[35]*y[41]*y[55]*y[56]-2.*x1*y[1\
]*y[2]*y[3]*y[4]*y[12]*y[33]*y[35]*y[41]*y[55]*y[56]-x0*x1*y[1]*y[5]*y[12]*\
y[33]*y[35]*y[41]*y[55]*y[56]-2.*x1*y[1]*y[3]*y[5]*y[12]*y[33]*y[35]*y[41]*\
y[55]*y[56]-2.*x1*x2*y[1]*y[3]*y[5]*y[12]*y[33]*y[35]*y[41]*y[55]*y[56]-x0*\
x1*y[1]*y[6]*y[12]*y[33]*y[35]*y[41]*y[55]*y[56]-x0*x1*x2*y[1]*y[6]*y[12]*y\
[33]*y[35]*y[41]*y[55]*y[56]-x0*x1*y[1]*y[7]*y[12]*y[33]*y[35]*y[41]*y[55]*\
y[56]-x0*x1*x2*y[1]*y[7]*y[12]*y[33]*y[35]*y[41]*y[55]*y[56]-2.*x0*y[1]*y[7\
]*y[8]*y[12]*y[33]*y[35]*y[41]*y[55]*y[56]-4.*y[1]*y[3]*y[7]*y[8]*y[12]*y[3\
3]*y[35]*y[41]*y[55]*y[56]-4.*x2*y[1]*y[3]*y[7]*y[8]*y[12]*y[33]*y[35]*y[41\
]*y[55]*y[56]+x0*x1*y[1]*y[9]*y[12]*y[33]*y[35]*y[41]*y[55]*y[56]+2.*x1*x2*\
y[1]*y[3]*y[9]*y[12]*y[33]*y[35]*y[41]*y[55]*y[56]-x0*y[1]*y[7]*y[8]*y[57]*\
y[58]*y[59]-y[1]*y[3]*y[7]*y[8]*y[57]*y[58]*y[59]-x2*y[1]*y[3]*y[7]*y[8]*y[\
57]*y[58]*y[59]-x0*x2*y[1]*y[2]*y[12]*y[33]*y[35]*y[65]*y[72]*y[73]-x0*x1*x\
2*y[1]*y[2]*y[12]*y[33]*y[35]*y[65]*y[72]*y[73]-2.*x1*x2*y[1]*y[2]*y[3]*y[1\
2]*y[33]*y[35]*y[65]*y[72]*y[73]-2.*x0*y[1]*y[2]*y[4]*y[12]*y[33]*y[35]*y[6\
5]*y[72]*y[73]-4.*x1*y[1]*y[2]*y[3]*y[4]*y[12]*y[33]*y[35]*y[65]*y[72]*y[73\
]-x0*x2*y[1]*y[5]*y[12]*y[33]*y[35]*y[65]*y[72]*y[73]-2.*x1*x2*y[1]*y[3]*y[\
5]*y[12]*y[33]*y[35]*y[65]*y[72]*y[73]-x0*x1*x2*y[1]*y[6]*y[12]*y[33]*y[35]\
*y[65]*y[72]*y[73]-x0*x1*x2*y[1]*y[7]*y[12]*y[33]*y[35]*y[65]*y[72]*y[73]-2\
.*x2*y[1]*y[3]*y[7]*y[8]*y[12]*y[33]*y[35]*y[65]*y[72]*y[73]+x0*x2*y[1]*y[9\
]*y[12]*y[33]*y[35]*y[65]*y[72]*y[73]+2.*x1*x2*y[1]*y[3]*y[9]*y[12]*y[33]*y\
[35]*y[65]*y[72]*y[73]-x0*x1*x2*y[1]*y[2]*y[41]*y[55]*y[56]*y[65]*y[72]*y[7\
3]-x1*x2*y[1]*y[2]*y[3]*y[41]*y[55]*y[56]*y[65]*y[72]*y[73]-2.*x1*y[1]*y[2]\
*y[3]*y[4]*y[41]*y[55]*y[56]*y[65]*y[72]*y[73]-x1*x2*y[1]*y[3]*y[5]*y[41]*y\
[55]*y[56]*y[65]*y[72]*y[73]-x0*x1*x2*y[1]*y[6]*y[41]*y[55]*y[56]*y[65]*y[7\
2]*y[73]-x0*x1*x2*y[1]*y[7]*y[41]*y[55]*y[56]*y[65]*y[72]*y[73]-2.*x2*y[1]*\
y[3]*y[7]*y[8]*y[41]*y[55]*y[56]*y[65]*y[72]*y[73]+x1*x2*y[1]*y[3]*y[9]*y[4\
1]*y[55]*y[56]*y[65]*y[72]*y[73]-x0*y[1]*y[2]*y[4]*y[74]*y[75]*y[76]-x1*y[1\
]*y[2]*y[3]*y[4]*y[74]*y[75]*y[76])+pow(lambda,4)*(y[1]*y[3]*y[7]*y[8]*y[13\
]*y[34]*y[36]*y[57]*y[58]*y[59]+x2*y[1]*y[3]*y[7]*y[8]*y[13]*y[34]*y[36]*y[\
57]*y[58]*y[59]+x1*x2*y[1]*y[2]*y[3]*y[13]*y[34]*y[36]*y[41]*y[55]*y[56]*y[\
65]*y[72]*y[73]+2.*x1*y[1]*y[2]*y[3]*y[4]*y[13]*y[34]*y[36]*y[41]*y[55]*y[5\
6]*y[65]*y[72]*y[73]+x1*x2*y[1]*y[3]*y[5]*y[13]*y[34]*y[36]*y[41]*y[55]*y[5\
6]*y[65]*y[72]*y[73]+2.*x2*y[1]*y[3]*y[7]*y[8]*y[13]*y[34]*y[36]*y[41]*y[55\
]*y[56]*y[65]*y[72]*y[73]-x1*x2*y[1]*y[3]*y[9]*y[13]*y[34]*y[36]*y[41]*y[55\
]*y[56]*y[65]*y[72]*y[73]+2.*x2*y[1]*y[3]*y[7]*y[8]*y[12]*y[33]*y[35]*y[57]\
*y[58]*y[59]*y[65]*y[72]*y[73]+x1*y[1]*y[2]*y[3]*y[4]*y[13]*y[34]*y[36]*y[7\
4]*y[75]*y[76]+2.*x1*y[1]*y[2]*y[3]*y[4]*y[12]*y[33]*y[35]*y[41]*y[55]*y[56\
]*y[74]*y[75]*y[76]),2)+pow(lambda*(-(x0*x2*y[1]*y[2]*y[12]*y[33]*y[35])-x0\
*x1*x2*y[1]*y[2]*y[12]*y[33]*y[35]-2.*x1*x2*y[1]*y[2]*y[3]*y[12]*y[33]*y[35\
]-x0*y[1]*y[2]*y[4]*y[12]*y[33]*y[35]-2.*x1*y[1]*y[2]*y[3]*y[4]*y[12]*y[33]\
*y[35]-x0*y[1]*y[5]*y[12]*y[33]*y[35]-x0*x1*y[1]*y[5]*y[12]*y[33]*y[35]-x0*\
x2*y[1]*y[5]*y[12]*y[33]*y[35]-2.*x1*y[1]*y[3]*y[5]*y[12]*y[33]*y[35]-2.*x1\
*x2*y[1]*y[3]*y[5]*y[12]*y[33]*y[35]-x0*x1*y[1]*y[6]*y[12]*y[33]*y[35]-x0*x\
1*x2*y[1]*y[6]*y[12]*y[33]*y[35]-x0*x1*y[1]*y[7]*y[12]*y[33]*y[35]-x0*x1*x2\
*y[1]*y[7]*y[12]*y[33]*y[35]-x0*y[1]*y[7]*y[8]*y[12]*y[33]*y[35]-2.*y[1]*y[\
3]*y[7]*y[8]*y[12]*y[33]*y[35]-2.*x2*y[1]*y[3]*y[7]*y[8]*y[12]*y[33]*y[35]+\
x0*x1*y[1]*y[9]*y[12]*y[33]*y[35]+x0*x2*y[1]*y[9]*y[12]*y[33]*y[35]+2.*x1*x\
2*y[1]*y[3]*y[9]*y[12]*y[33]*y[35]-x0*x1*x2*y[1]*y[2]*y[41]*y[55]*y[56]-x1*\
x2*y[1]*y[2]*y[3]*y[41]*y[55]*y[56]-x1*y[1]*y[2]*y[3]*y[4]*y[41]*y[55]*y[56\
]-x0*x1*y[1]*y[5]*y[41]*y[55]*y[56]-x1*y[1]*y[3]*y[5]*y[41]*y[55]*y[56]-x1*\
x2*y[1]*y[3]*y[5]*y[41]*y[55]*y[56]-x1*y[1]*y[6]*y[41]*y[55]*y[56]-x0*x1*y[\
1]*y[6]*y[41]*y[55]*y[56]-x0*x1*x2*y[1]*y[6]*y[41]*y[55]*y[56]-x0*x1*y[1]*y\
[7]*y[41]*y[55]*y[56]-x0*x1*x2*y[1]*y[7]*y[41]*y[55]*y[56]-2.*x0*y[1]*y[7]*\
y[8]*y[41]*y[55]*y[56]-2.*y[1]*y[3]*y[7]*y[8]*y[41]*y[55]*y[56]-2.*x2*y[1]*\
y[3]*y[7]*y[8]*y[41]*y[55]*y[56]+x0*x1*y[1]*y[9]*y[41]*y[55]*y[56]+x1*x2*y[\
1]*y[3]*y[9]*y[41]*y[55]*y[56]-x0*x2*y[1]*y[2]*y[65]*y[72]*y[73]-x0*x1*x2*y\
[1]*y[2]*y[65]*y[72]*y[73]-x1*x2*y[1]*y[2]*y[3]*y[65]*y[72]*y[73]-2.*x0*y[1\
]*y[2]*y[4]*y[65]*y[72]*y[73]-2.*x1*y[1]*y[2]*y[3]*y[4]*y[65]*y[72]*y[73]-x\
0*x2*y[1]*y[5]*y[65]*y[72]*y[73]-x1*x2*y[1]*y[3]*y[5]*y[65]*y[72]*y[73]-x2*\
y[1]*y[6]*y[65]*y[72]*y[73]-x0*x1*x2*y[1]*y[6]*y[65]*y[72]*y[73]-x0*x1*x2*y\
[1]*y[7]*y[65]*y[72]*y[73]-x2*y[1]*y[3]*y[7]*y[8]*y[65]*y[72]*y[73]+x0*x2*y\
[1]*y[9]*y[65]*y[72]*y[73]+x1*x2*y[1]*y[3]*y[9]*y[65]*y[72]*y[73])+pow(lamb\
da,3)*(x1*x2*y[1]*y[2]*y[3]*y[13]*y[34]*y[36]*y[41]*y[55]*y[56]+x1*y[1]*y[2\
]*y[3]*y[4]*y[13]*y[34]*y[36]*y[41]*y[55]*y[56]+x1*y[1]*y[3]*y[5]*y[13]*y[3\
4]*y[36]*y[41]*y[55]*y[56]+x1*x2*y[1]*y[3]*y[5]*y[13]*y[34]*y[36]*y[41]*y[5\
5]*y[56]+2.*y[1]*y[3]*y[7]*y[8]*y[13]*y[34]*y[36]*y[41]*y[55]*y[56]+2.*x2*y\
[1]*y[3]*y[7]*y[8]*y[13]*y[34]*y[36]*y[41]*y[55]*y[56]-x1*x2*y[1]*y[3]*y[9]\
*y[13]*y[34]*y[36]*y[41]*y[55]*y[56]+x0*y[1]*y[7]*y[8]*y[12]*y[33]*y[35]*y[\
57]*y[58]*y[59]+2.*y[1]*y[3]*y[7]*y[8]*y[12]*y[33]*y[35]*y[57]*y[58]*y[59]+\
2.*x2*y[1]*y[3]*y[7]*y[8]*y[12]*y[33]*y[35]*y[57]*y[58]*y[59]+x1*x2*y[1]*y[\
2]*y[3]*y[13]*y[34]*y[36]*y[65]*y[72]*y[73]+2.*x1*y[1]*y[2]*y[3]*y[4]*y[13]\
*y[34]*y[36]*y[65]*y[72]*y[73]+x1*x2*y[1]*y[3]*y[5]*y[13]*y[34]*y[36]*y[65]\
*y[72]*y[73]+x2*y[1]*y[3]*y[7]*y[8]*y[13]*y[34]*y[36]*y[65]*y[72]*y[73]-x1*\
x2*y[1]*y[3]*y[9]*y[13]*y[34]*y[36]*y[65]*y[72]*y[73]+x0*x1*x2*y[1]*y[2]*y[\
12]*y[33]*y[35]*y[41]*y[55]*y[56]*y[65]*y[72]*y[73]+2.*x1*x2*y[1]*y[2]*y[3]\
*y[12]*y[33]*y[35]*y[41]*y[55]*y[56]*y[65]*y[72]*y[73]+4.*x1*y[1]*y[2]*y[3]\
*y[4]*y[12]*y[33]*y[35]*y[41]*y[55]*y[56]*y[65]*y[72]*y[73]+2.*x1*x2*y[1]*y\
[3]*y[5]*y[12]*y[33]*y[35]*y[41]*y[55]*y[56]*y[65]*y[72]*y[73]+x0*x1*x2*y[1\
]*y[6]*y[12]*y[33]*y[35]*y[41]*y[55]*y[56]*y[65]*y[72]*y[73]+x0*x1*x2*y[1]*\
y[7]*y[12]*y[33]*y[35]*y[41]*y[55]*y[56]*y[65]*y[72]*y[73]+4.*x2*y[1]*y[3]*\
y[7]*y[8]*y[12]*y[33]*y[35]*y[41]*y[55]*y[56]*y[65]*y[72]*y[73]-2.*x1*x2*y[\
1]*y[3]*y[9]*y[12]*y[33]*y[35]*y[41]*y[55]*y[56]*y[65]*y[72]*y[73]+x2*y[1]*\
y[3]*y[7]*y[8]*y[57]*y[58]*y[59]*y[65]*y[72]*y[73]+x0*y[1]*y[2]*y[4]*y[12]*\
y[33]*y[35]*y[74]*y[75]*y[76]+2.*x1*y[1]*y[2]*y[3]*y[4]*y[12]*y[33]*y[35]*y\
[74]*y[75]*y[76]+x1*y[1]*y[2]*y[3]*y[4]*y[41]*y[55]*y[56]*y[74]*y[75]*y[76]\
)+pow(lambda,5)*(-(x2*y[1]*y[3]*y[7]*y[8]*y[13]*y[34]*y[36]*y[57]*y[58]*y[5\
9]*y[65]*y[72]*y[73])-x1*y[1]*y[2]*y[3]*y[4]*y[13]*y[34]*y[36]*y[41]*y[55]*\
y[56]*y[74]*y[75]*y[76]),2);
return (FOUT);
}
double Ps8(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
double y[77];
double FOUT;
y[1]=1./bi;
y[2]=em[0];
y[3]=x2*x2;
y[4]=em[1];
y[5]=em[2];
y[6]=em[3];
y[7]=x1*x1;
y[8]=esx[0];
y[9]=-x0;
y[10]=1.+y[9];
y[11]=x2*y[1]*y[2];
y[12]=x1*x2*y[1]*y[2];
y[13]=2.*x0*x1*x2*y[1]*y[2];
y[14]=y[1]*y[2]*y[3];
y[15]=2.*x0*x1*y[1]*y[2]*y[3];
y[16]=y[1]*y[4];
y[17]=x1*y[1]*y[4];
y[18]=2.*x0*x1*y[1]*y[4];
y[19]=x2*y[1]*y[4];
y[20]=2.*x0*x1*x2*y[1]*y[4];
y[21]=x1*y[1]*y[5];
y[22]=x1*x2*y[1]*y[5];
y[23]=x1*y[1]*y[6];
y[24]=y[1]*y[6]*y[7];
y[25]=2.*x0*y[1]*y[6]*y[7];
y[26]=x1*x2*y[1]*y[6];
y[27]=2.*x0*x2*y[1]*y[6]*y[7];
y[28]=-(x1*y[1]*y[8]);
y[29]=-(x2*y[1]*y[8]);
y[30]=-2.*x0*x1*x2*y[1]*y[8];
y[31]=y[11]+y[12]+y[13]+y[14]+y[15]+y[16]+y[17]+y[18]+y[19]+y[20]+y[21]+y[22\
]+y[23]+y[24]+y[25]+y[26]+y[27]+y[28]+y[29]+y[30];
y[32]=lrs[0];
y[33]=x0*x0;
y[34]=-x1;
y[35]=1.+y[34];
y[36]=x0*x2*y[1]*y[2];
y[37]=x2*y[1]*y[2]*y[33];
y[38]=y[1]*y[2]*y[3]*y[33];
y[39]=x0*y[1]*y[4];
y[40]=y[1]*y[4]*y[33];
y[41]=x2*y[1]*y[4]*y[33];
y[42]=y[1]*y[5];
y[43]=x0*y[1]*y[5];
y[44]=x0*x2*y[1]*y[5];
y[45]=x0*y[1]*y[6];
y[46]=2.*x0*x1*y[1]*y[6];
y[47]=2.*x1*y[1]*y[6]*y[33];
y[48]=x0*x2*y[1]*y[6];
y[49]=2.*x1*x2*y[1]*y[6]*y[33];
y[50]=-(x0*y[1]*y[8]);
y[51]=-(x2*y[1]*y[8]*y[33]);
y[52]=y[36]+y[37]+y[38]+y[39]+y[40]+y[41]+y[42]+y[43]+y[44]+y[45]+y[46]+y[47\
]+y[48]+y[49]+y[50]+y[51];
y[53]=lrs[1];
y[54]=-x2;
y[55]=1.+y[54];
y[56]=x0*y[1]*y[2];
y[57]=x0*x1*y[1]*y[2];
y[58]=x1*y[1]*y[2]*y[33];
y[59]=2.*x0*x2*y[1]*y[2];
y[60]=2.*x1*x2*y[1]*y[2]*y[33];
y[61]=x1*y[1]*y[4]*y[33];
y[62]=x0*x1*y[1]*y[5];
y[63]=x0*x1*y[1]*y[6];
y[64]=y[1]*y[6]*y[7]*y[33];
y[65]=-(x1*y[1]*y[8]*y[33]);
y[66]=y[39]+y[42]+y[50]+y[56]+y[57]+y[58]+y[59]+y[60]+y[61]+y[62]+y[63]+y[64\
]+y[65];
y[67]=lrs[2];
y[68]=pow(y[10],2);
y[69]=pow(y[31],2);
y[70]=pow(y[32],2);
y[71]=pow(y[35],2);
y[72]=pow(y[52],2);
y[73]=pow(y[53],2);
y[74]=pow(y[55],2);
y[75]=pow(y[66],2);
y[76]=pow(y[67],2);
FOUT=lambda*(-(x0*x2*y[1]*y[2]*y[10]*y[31]*y[32])-x0*x1*x2*y[1]*y[2]*y[10]*y\
[31]*y[32]-x0*y[1]*y[2]*y[3]*y[10]*y[31]*y[32]-x0*y[1]*y[4]*y[10]*y[31]*y[3\
2]-x0*x1*y[1]*y[4]*y[10]*y[31]*y[32]-x0*x2*y[1]*y[4]*y[10]*y[31]*y[32]-x0*x\
1*y[1]*y[5]*y[10]*y[31]*y[32]-x0*x1*x2*y[1]*y[5]*y[10]*y[31]*y[32]-x0*x1*y[\
1]*y[6]*y[10]*y[31]*y[32]-x0*x1*x2*y[1]*y[6]*y[10]*y[31]*y[32]-x0*y[1]*y[6]\
*y[7]*y[10]*y[31]*y[32]+x0*x1*y[1]*y[8]*y[10]*y[31]*y[32]+x0*x2*y[1]*y[8]*y\
[10]*y[31]*y[32]-2.*x1*x2*y[1]*y[2]*y[10]*y[31]*y[32]*y[33]-2.*x1*y[1]*y[2]\
*y[3]*y[10]*y[31]*y[32]*y[33]-2.*x1*y[1]*y[4]*y[10]*y[31]*y[32]*y[33]-2.*x1\
*x2*y[1]*y[4]*y[10]*y[31]*y[32]*y[33]-2.*y[1]*y[6]*y[7]*y[10]*y[31]*y[32]*y\
[33]-2.*x2*y[1]*y[6]*y[7]*y[10]*y[31]*y[32]*y[33]+2.*x1*x2*y[1]*y[8]*y[10]*\
y[31]*y[32]*y[33]-x0*x1*x2*y[1]*y[2]*y[35]*y[52]*y[53]-x0*x1*y[1]*y[4]*y[35\
]*y[52]*y[53]-x1*y[1]*y[5]*y[35]*y[52]*y[53]-x0*x1*y[1]*y[5]*y[35]*y[52]*y[\
53]-x0*x1*x2*y[1]*y[5]*y[35]*y[52]*y[53]-x0*x1*y[1]*y[6]*y[35]*y[52]*y[53]-\
x0*x1*x2*y[1]*y[6]*y[35]*y[52]*y[53]-2.*x0*y[1]*y[6]*y[7]*y[35]*y[52]*y[53]\
+x0*x1*y[1]*y[8]*y[35]*y[52]*y[53]-x1*x2*y[1]*y[2]*y[33]*y[35]*y[52]*y[53]-\
x1*y[1]*y[2]*y[3]*y[33]*y[35]*y[52]*y[53]-x1*y[1]*y[4]*y[33]*y[35]*y[52]*y[\
53]-x1*x2*y[1]*y[4]*y[33]*y[35]*y[52]*y[53]-2.*y[1]*y[6]*y[7]*y[33]*y[35]*y\
[52]*y[53]-2.*x2*y[1]*y[6]*y[7]*y[33]*y[35]*y[52]*y[53]+x1*x2*y[1]*y[8]*y[3\
3]*y[35]*y[52]*y[53]-x0*x2*y[1]*y[2]*y[55]*y[66]*y[67]-x0*x1*x2*y[1]*y[2]*y\
[55]*y[66]*y[67]-2.*x0*y[1]*y[2]*y[3]*y[55]*y[66]*y[67]-x0*x2*y[1]*y[4]*y[5\
5]*y[66]*y[67]-x2*y[1]*y[5]*y[55]*y[66]*y[67]-x0*x1*x2*y[1]*y[5]*y[55]*y[66\
]*y[67]-x0*x1*x2*y[1]*y[6]*y[55]*y[66]*y[67]+x0*x2*y[1]*y[8]*y[55]*y[66]*y[\
67]-x1*x2*y[1]*y[2]*y[33]*y[55]*y[66]*y[67]-2.*x1*y[1]*y[2]*y[3]*y[33]*y[55\
]*y[66]*y[67]-x1*x2*y[1]*y[4]*y[33]*y[55]*y[66]*y[67]-x2*y[1]*y[6]*y[7]*y[3\
3]*y[55]*y[66]*y[67]+x1*x2*y[1]*y[8]*y[33]*y[55]*y[66]*y[67])+pow(lambda,3)\
*(x0*x1*x2*y[1]*y[2]*y[10]*y[31]*y[32]*y[35]*y[52]*y[53]*y[55]*y[66]*y[67]+\
x0*x1*x2*y[1]*y[5]*y[10]*y[31]*y[32]*y[35]*y[52]*y[53]*y[55]*y[66]*y[67]+x0\
*x1*x2*y[1]*y[6]*y[10]*y[31]*y[32]*y[35]*y[52]*y[53]*y[55]*y[66]*y[67]+2.*x\
1*x2*y[1]*y[2]*y[10]*y[31]*y[32]*y[33]*y[35]*y[52]*y[53]*y[55]*y[66]*y[67]+\
4.*x1*y[1]*y[2]*y[3]*y[10]*y[31]*y[32]*y[33]*y[35]*y[52]*y[53]*y[55]*y[66]*\
y[67]+2.*x1*x2*y[1]*y[4]*y[10]*y[31]*y[32]*y[33]*y[35]*y[52]*y[53]*y[55]*y[\
66]*y[67]+4.*x2*y[1]*y[6]*y[7]*y[10]*y[31]*y[32]*y[33]*y[35]*y[52]*y[53]*y[\
55]*y[66]*y[67]-2.*x1*x2*y[1]*y[8]*y[10]*y[31]*y[32]*y[33]*y[35]*y[52]*y[53\
]*y[55]*y[66]*y[67]+x1*x2*y[1]*y[2]*y[33]*y[35]*y[52]*y[53]*y[68]*y[69]*y[7\
0]+x1*y[1]*y[2]*y[3]*y[33]*y[35]*y[52]*y[53]*y[68]*y[69]*y[70]+x1*y[1]*y[4]\
*y[33]*y[35]*y[52]*y[53]*y[68]*y[69]*y[70]+x1*x2*y[1]*y[4]*y[33]*y[35]*y[52\
]*y[53]*y[68]*y[69]*y[70]+2.*y[1]*y[6]*y[7]*y[33]*y[35]*y[52]*y[53]*y[68]*y\
[69]*y[70]+2.*x2*y[1]*y[6]*y[7]*y[33]*y[35]*y[52]*y[53]*y[68]*y[69]*y[70]-x\
1*x2*y[1]*y[8]*y[33]*y[35]*y[52]*y[53]*y[68]*y[69]*y[70]+x1*x2*y[1]*y[2]*y[\
33]*y[55]*y[66]*y[67]*y[68]*y[69]*y[70]+2.*x1*y[1]*y[2]*y[3]*y[33]*y[55]*y[\
66]*y[67]*y[68]*y[69]*y[70]+x1*x2*y[1]*y[4]*y[33]*y[55]*y[66]*y[67]*y[68]*y\
[69]*y[70]+x2*y[1]*y[6]*y[7]*y[33]*y[55]*y[66]*y[67]*y[68]*y[69]*y[70]-x1*x\
2*y[1]*y[8]*y[33]*y[55]*y[66]*y[67]*y[68]*y[69]*y[70]+x0*y[1]*y[6]*y[7]*y[1\
0]*y[31]*y[32]*y[71]*y[72]*y[73]+2.*y[1]*y[6]*y[7]*y[10]*y[31]*y[32]*y[33]*\
y[71]*y[72]*y[73]+2.*x2*y[1]*y[6]*y[7]*y[10]*y[31]*y[32]*y[33]*y[71]*y[72]*\
y[73]+x2*y[1]*y[6]*y[7]*y[33]*y[55]*y[66]*y[67]*y[71]*y[72]*y[73]+x0*y[1]*y\
[2]*y[3]*y[10]*y[31]*y[32]*y[74]*y[75]*y[76]+2.*x1*y[1]*y[2]*y[3]*y[10]*y[3\
1]*y[32]*y[33]*y[74]*y[75]*y[76]+x1*y[1]*y[2]*y[3]*y[33]*y[35]*y[52]*y[53]*\
y[74]*y[75]*y[76])+pow(lambda,5)*(-(x2*y[1]*y[6]*y[7]*y[33]*y[55]*y[66]*y[6\
7]*y[68]*y[69]*y[70]*y[71]*y[72]*y[73])-x1*y[1]*y[2]*y[3]*y[33]*y[35]*y[52]\
*y[53]*y[68]*y[69]*y[70]*y[74]*y[75]*y[76]);
return (FOUT);
}
double Pa8(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
double y[77];
double FOUT;
y[1]=1./bi;
y[2]=em[0];
y[3]=x0*x0;
y[4]=x2*x2;
y[5]=em[1];
y[6]=em[2];
y[7]=em[3];
y[8]=x1*x1;
y[9]=esx[0];
y[10]=x1*y[1]*y[6];
y[11]=-x0;
y[12]=1.+y[11];
y[13]=pow(y[12],2);
y[14]=x2*y[1]*y[2];
y[15]=x1*x2*y[1]*y[2];
y[16]=2.*x0*x1*x2*y[1]*y[2];
y[17]=y[1]*y[2]*y[4];
y[18]=2.*x0*x1*y[1]*y[2]*y[4];
y[19]=y[1]*y[5];
y[20]=x1*y[1]*y[5];
y[21]=2.*x0*x1*y[1]*y[5];
y[22]=x2*y[1]*y[5];
y[23]=2.*x0*x1*x2*y[1]*y[5];
y[24]=x1*x2*y[1]*y[6];
y[25]=x1*y[1]*y[7];
y[26]=y[1]*y[7]*y[8];
y[27]=2.*x0*y[1]*y[7]*y[8];
y[28]=x1*x2*y[1]*y[7];
y[29]=2.*x0*x2*y[1]*y[7]*y[8];
y[30]=-(x1*y[1]*y[9]);
y[31]=-(x2*y[1]*y[9]);
y[32]=-2.*x0*x1*x2*y[1]*y[9];
y[33]=y[10]+y[14]+y[15]+y[16]+y[17]+y[18]+y[19]+y[20]+y[21]+y[22]+y[23]+y[24\
]+y[25]+y[26]+y[27]+y[28]+y[29]+y[30]+y[31]+y[32];
y[34]=pow(y[33],2);
y[35]=lrs[0];
y[36]=pow(y[35],2);
y[37]=x0*x2*y[1]*y[2];
y[38]=x0*y[1]*y[5];
y[39]=y[1]*y[6];
y[40]=-x1;
y[41]=1.+y[40];
y[42]=x2*y[1]*y[2]*y[3];
y[43]=y[1]*y[2]*y[3]*y[4];
y[44]=y[1]*y[3]*y[5];
y[45]=x2*y[1]*y[3]*y[5];
y[46]=x0*y[1]*y[6];
y[47]=x0*x2*y[1]*y[6];
y[48]=x0*y[1]*y[7];
y[49]=2.*x0*x1*y[1]*y[7];
y[50]=2.*x1*y[1]*y[3]*y[7];
y[51]=x0*x2*y[1]*y[7];
y[52]=2.*x1*x2*y[1]*y[3]*y[7];
y[53]=-(x0*y[1]*y[9]);
y[54]=-(x2*y[1]*y[3]*y[9]);
y[55]=y[37]+y[38]+y[39]+y[42]+y[43]+y[44]+y[45]+y[46]+y[47]+y[48]+y[49]+y[50\
]+y[51]+y[52]+y[53]+y[54];
y[56]=lrs[1];
y[57]=pow(y[41],2);
y[58]=pow(y[55],2);
y[59]=pow(y[56],2);
y[60]=x1*y[1]*y[3]*y[5];
y[61]=x0*x1*y[1]*y[6];
y[62]=x0*x1*y[1]*y[7];
y[63]=y[1]*y[3]*y[7]*y[8];
y[64]=-x2;
y[65]=1.+y[64];
y[66]=x0*y[1]*y[2];
y[67]=x0*x1*y[1]*y[2];
y[68]=x1*y[1]*y[2]*y[3];
y[69]=2.*x0*x2*y[1]*y[2];
y[70]=2.*x1*x2*y[1]*y[2]*y[3];
y[71]=-(x1*y[1]*y[3]*y[9]);
y[72]=y[38]+y[39]+y[53]+y[60]+y[61]+y[62]+y[63]+y[66]+y[67]+y[68]+y[69]+y[70\
]+y[71];
y[73]=lrs[2];
y[74]=pow(y[65],2);
y[75]=pow(y[72],2);
y[76]=pow(y[73],2);
FOUT=(lambda*(-(x0*x2*y[1]*y[2]*y[12]*y[33]*y[35])-x0*x1*x2*y[1]*y[2]*y[12]*\
y[33]*y[35]-2.*x1*x2*y[1]*y[2]*y[3]*y[12]*y[33]*y[35]-x0*y[1]*y[2]*y[4]*y[1\
2]*y[33]*y[35]-2.*x1*y[1]*y[2]*y[3]*y[4]*y[12]*y[33]*y[35]-x0*y[1]*y[5]*y[1\
2]*y[33]*y[35]-x0*x1*y[1]*y[5]*y[12]*y[33]*y[35]-x0*x2*y[1]*y[5]*y[12]*y[33\
]*y[35]-2.*x1*y[1]*y[3]*y[5]*y[12]*y[33]*y[35]-2.*x1*x2*y[1]*y[3]*y[5]*y[12\
]*y[33]*y[35]-x0*x1*y[1]*y[6]*y[12]*y[33]*y[35]-x0*x1*x2*y[1]*y[6]*y[12]*y[\
33]*y[35]-x0*x1*y[1]*y[7]*y[12]*y[33]*y[35]-x0*x1*x2*y[1]*y[7]*y[12]*y[33]*\
y[35]-x0*y[1]*y[7]*y[8]*y[12]*y[33]*y[35]-2.*y[1]*y[3]*y[7]*y[8]*y[12]*y[33\
]*y[35]-2.*x2*y[1]*y[3]*y[7]*y[8]*y[12]*y[33]*y[35]+x0*x1*y[1]*y[9]*y[12]*y\
[33]*y[35]+x0*x2*y[1]*y[9]*y[12]*y[33]*y[35]+2.*x1*x2*y[1]*y[3]*y[9]*y[12]*\
y[33]*y[35]-x0*x1*x2*y[1]*y[2]*y[41]*y[55]*y[56]-x1*x2*y[1]*y[2]*y[3]*y[41]\
*y[55]*y[56]-x1*y[1]*y[2]*y[3]*y[4]*y[41]*y[55]*y[56]-x0*x1*y[1]*y[5]*y[41]\
*y[55]*y[56]-x1*y[1]*y[3]*y[5]*y[41]*y[55]*y[56]-x1*x2*y[1]*y[3]*y[5]*y[41]\
*y[55]*y[56]-x1*y[1]*y[6]*y[41]*y[55]*y[56]-x0*x1*y[1]*y[6]*y[41]*y[55]*y[5\
6]-x0*x1*x2*y[1]*y[6]*y[41]*y[55]*y[56]-x0*x1*y[1]*y[7]*y[41]*y[55]*y[56]-x\
0*x1*x2*y[1]*y[7]*y[41]*y[55]*y[56]-2.*x0*y[1]*y[7]*y[8]*y[41]*y[55]*y[56]-\
2.*y[1]*y[3]*y[7]*y[8]*y[41]*y[55]*y[56]-2.*x2*y[1]*y[3]*y[7]*y[8]*y[41]*y[\
55]*y[56]+x0*x1*y[1]*y[9]*y[41]*y[55]*y[56]+x1*x2*y[1]*y[3]*y[9]*y[41]*y[55\
]*y[56]-x0*x2*y[1]*y[2]*y[65]*y[72]*y[73]-x0*x1*x2*y[1]*y[2]*y[65]*y[72]*y[\
73]-x1*x2*y[1]*y[2]*y[3]*y[65]*y[72]*y[73]-2.*x0*y[1]*y[2]*y[4]*y[65]*y[72]\
*y[73]-2.*x1*y[1]*y[2]*y[3]*y[4]*y[65]*y[72]*y[73]-x0*x2*y[1]*y[5]*y[65]*y[\
72]*y[73]-x1*x2*y[1]*y[3]*y[5]*y[65]*y[72]*y[73]-x2*y[1]*y[6]*y[65]*y[72]*y\
[73]-x0*x1*x2*y[1]*y[6]*y[65]*y[72]*y[73]-x0*x1*x2*y[1]*y[7]*y[65]*y[72]*y[\
73]-x2*y[1]*y[3]*y[7]*y[8]*y[65]*y[72]*y[73]+x0*x2*y[1]*y[9]*y[65]*y[72]*y[\
73]+x1*x2*y[1]*y[3]*y[9]*y[65]*y[72]*y[73])+pow(lambda,3)*(x1*x2*y[1]*y[2]*\
y[3]*y[13]*y[34]*y[36]*y[41]*y[55]*y[56]+x1*y[1]*y[2]*y[3]*y[4]*y[13]*y[34]\
*y[36]*y[41]*y[55]*y[56]+x1*y[1]*y[3]*y[5]*y[13]*y[34]*y[36]*y[41]*y[55]*y[\
56]+x1*x2*y[1]*y[3]*y[5]*y[13]*y[34]*y[36]*y[41]*y[55]*y[56]+2.*y[1]*y[3]*y\
[7]*y[8]*y[13]*y[34]*y[36]*y[41]*y[55]*y[56]+2.*x2*y[1]*y[3]*y[7]*y[8]*y[13\
]*y[34]*y[36]*y[41]*y[55]*y[56]-x1*x2*y[1]*y[3]*y[9]*y[13]*y[34]*y[36]*y[41\
]*y[55]*y[56]+x0*y[1]*y[7]*y[8]*y[12]*y[33]*y[35]*y[57]*y[58]*y[59]+2.*y[1]\
*y[3]*y[7]*y[8]*y[12]*y[33]*y[35]*y[57]*y[58]*y[59]+2.*x2*y[1]*y[3]*y[7]*y[\
8]*y[12]*y[33]*y[35]*y[57]*y[58]*y[59]+x1*x2*y[1]*y[2]*y[3]*y[13]*y[34]*y[3\
6]*y[65]*y[72]*y[73]+2.*x1*y[1]*y[2]*y[3]*y[4]*y[13]*y[34]*y[36]*y[65]*y[72\
]*y[73]+x1*x2*y[1]*y[3]*y[5]*y[13]*y[34]*y[36]*y[65]*y[72]*y[73]+x2*y[1]*y[\
3]*y[7]*y[8]*y[13]*y[34]*y[36]*y[65]*y[72]*y[73]-x1*x2*y[1]*y[3]*y[9]*y[13]\
*y[34]*y[36]*y[65]*y[72]*y[73]+x0*x1*x2*y[1]*y[2]*y[12]*y[33]*y[35]*y[41]*y\
[55]*y[56]*y[65]*y[72]*y[73]+2.*x1*x2*y[1]*y[2]*y[3]*y[12]*y[33]*y[35]*y[41\
]*y[55]*y[56]*y[65]*y[72]*y[73]+4.*x1*y[1]*y[2]*y[3]*y[4]*y[12]*y[33]*y[35]\
*y[41]*y[55]*y[56]*y[65]*y[72]*y[73]+2.*x1*x2*y[1]*y[3]*y[5]*y[12]*y[33]*y[\
35]*y[41]*y[55]*y[56]*y[65]*y[72]*y[73]+x0*x1*x2*y[1]*y[6]*y[12]*y[33]*y[35\
]*y[41]*y[55]*y[56]*y[65]*y[72]*y[73]+x0*x1*x2*y[1]*y[7]*y[12]*y[33]*y[35]*\
y[41]*y[55]*y[56]*y[65]*y[72]*y[73]+4.*x2*y[1]*y[3]*y[7]*y[8]*y[12]*y[33]*y\
[35]*y[41]*y[55]*y[56]*y[65]*y[72]*y[73]-2.*x1*x2*y[1]*y[3]*y[9]*y[12]*y[33\
]*y[35]*y[41]*y[55]*y[56]*y[65]*y[72]*y[73]+x2*y[1]*y[3]*y[7]*y[8]*y[57]*y[\
58]*y[59]*y[65]*y[72]*y[73]+x0*y[1]*y[2]*y[4]*y[12]*y[33]*y[35]*y[74]*y[75]\
*y[76]+2.*x1*y[1]*y[2]*y[3]*y[4]*y[12]*y[33]*y[35]*y[74]*y[75]*y[76]+x1*y[1\
]*y[2]*y[3]*y[4]*y[41]*y[55]*y[56]*y[74]*y[75]*y[76])+pow(lambda,5)*(-(x2*y\
[1]*y[3]*y[7]*y[8]*y[13]*y[34]*y[36]*y[57]*y[58]*y[59]*y[65]*y[72]*y[73])-x\
1*y[1]*y[2]*y[3]*y[4]*y[13]*y[34]*y[36]*y[41]*y[55]*y[56]*y[74]*y[75]*y[76]\
))/(lambda*(x0*x1*x2*y[1]*y[2]+x1*x2*y[1]*y[2]*y[3]+x0*y[1]*y[2]*y[4]+x1*y[\
1]*y[2]*y[3]*y[4]+x0*x1*y[1]*y[5]+x0*x2*y[1]*y[5]+x1*x2*y[1]*y[3]*y[5]+x2*y\
[1]*y[6]+x0*x1*x2*y[1]*y[6]+x0*x1*x2*y[1]*y[7]+x0*y[1]*y[7]*y[8]+x2*y[1]*y[\
3]*y[7]*y[8]-x0*x1*y[1]*y[9]-x0*x2*y[1]*y[9]-x1*x2*y[1]*y[3]*y[9]+y[10]+y[3\
7]+y[38]+y[39]+y[60]+y[61]+y[62]+y[63]+lambda*lambda*(-(x1*x2*y[1]*y[2]*y[3\
]*y[13]*y[34]*y[36])-x1*y[1]*y[2]*y[3]*y[4]*y[13]*y[34]*y[36]-x1*y[1]*y[3]*\
y[5]*y[13]*y[34]*y[36]-x1*x2*y[1]*y[3]*y[5]*y[13]*y[34]*y[36]-y[1]*y[3]*y[7\
]*y[8]*y[13]*y[34]*y[36]-x2*y[1]*y[3]*y[7]*y[8]*y[13]*y[34]*y[36]+x1*x2*y[1\
]*y[3]*y[9]*y[13]*y[34]*y[36]-x0*x1*x2*y[1]*y[2]*y[12]*y[33]*y[35]*y[41]*y[\
55]*y[56]-2.*x1*x2*y[1]*y[2]*y[3]*y[12]*y[33]*y[35]*y[41]*y[55]*y[56]-2.*x1\
*y[1]*y[2]*y[3]*y[4]*y[12]*y[33]*y[35]*y[41]*y[55]*y[56]-x0*x1*y[1]*y[5]*y[\
12]*y[33]*y[35]*y[41]*y[55]*y[56]-2.*x1*y[1]*y[3]*y[5]*y[12]*y[33]*y[35]*y[\
41]*y[55]*y[56]-2.*x1*x2*y[1]*y[3]*y[5]*y[12]*y[33]*y[35]*y[41]*y[55]*y[56]\
-x0*x1*y[1]*y[6]*y[12]*y[33]*y[35]*y[41]*y[55]*y[56]-x0*x1*x2*y[1]*y[6]*y[1\
2]*y[33]*y[35]*y[41]*y[55]*y[56]-x0*x1*y[1]*y[7]*y[12]*y[33]*y[35]*y[41]*y[\
55]*y[56]-x0*x1*x2*y[1]*y[7]*y[12]*y[33]*y[35]*y[41]*y[55]*y[56]-2.*x0*y[1]\
*y[7]*y[8]*y[12]*y[33]*y[35]*y[41]*y[55]*y[56]-4.*y[1]*y[3]*y[7]*y[8]*y[12]\
*y[33]*y[35]*y[41]*y[55]*y[56]-4.*x2*y[1]*y[3]*y[7]*y[8]*y[12]*y[33]*y[35]*\
y[41]*y[55]*y[56]+x0*x1*y[1]*y[9]*y[12]*y[33]*y[35]*y[41]*y[55]*y[56]+2.*x1\
*x2*y[1]*y[3]*y[9]*y[12]*y[33]*y[35]*y[41]*y[55]*y[56]-x0*y[1]*y[7]*y[8]*y[\
57]*y[58]*y[59]-y[1]*y[3]*y[7]*y[8]*y[57]*y[58]*y[59]-x2*y[1]*y[3]*y[7]*y[8\
]*y[57]*y[58]*y[59]-x0*x2*y[1]*y[2]*y[12]*y[33]*y[35]*y[65]*y[72]*y[73]-x0*\
x1*x2*y[1]*y[2]*y[12]*y[33]*y[35]*y[65]*y[72]*y[73]-2.*x1*x2*y[1]*y[2]*y[3]\
*y[12]*y[33]*y[35]*y[65]*y[72]*y[73]-2.*x0*y[1]*y[2]*y[4]*y[12]*y[33]*y[35]\
*y[65]*y[72]*y[73]-4.*x1*y[1]*y[2]*y[3]*y[4]*y[12]*y[33]*y[35]*y[65]*y[72]*\
y[73]-x0*x2*y[1]*y[5]*y[12]*y[33]*y[35]*y[65]*y[72]*y[73]-2.*x1*x2*y[1]*y[3\
]*y[5]*y[12]*y[33]*y[35]*y[65]*y[72]*y[73]-x0*x1*x2*y[1]*y[6]*y[12]*y[33]*y\
[35]*y[65]*y[72]*y[73]-x0*x1*x2*y[1]*y[7]*y[12]*y[33]*y[35]*y[65]*y[72]*y[7\
3]-2.*x2*y[1]*y[3]*y[7]*y[8]*y[12]*y[33]*y[35]*y[65]*y[72]*y[73]+x0*x2*y[1]\
*y[9]*y[12]*y[33]*y[35]*y[65]*y[72]*y[73]+2.*x1*x2*y[1]*y[3]*y[9]*y[12]*y[3\
3]*y[35]*y[65]*y[72]*y[73]-x0*x1*x2*y[1]*y[2]*y[41]*y[55]*y[56]*y[65]*y[72]\
*y[73]-x1*x2*y[1]*y[2]*y[3]*y[41]*y[55]*y[56]*y[65]*y[72]*y[73]-2.*x1*y[1]*\
y[2]*y[3]*y[4]*y[41]*y[55]*y[56]*y[65]*y[72]*y[73]-x1*x2*y[1]*y[3]*y[5]*y[4\
1]*y[55]*y[56]*y[65]*y[72]*y[73]-x0*x1*x2*y[1]*y[6]*y[41]*y[55]*y[56]*y[65]\
*y[72]*y[73]-x0*x1*x2*y[1]*y[7]*y[41]*y[55]*y[56]*y[65]*y[72]*y[73]-2.*x2*y\
[1]*y[3]*y[7]*y[8]*y[41]*y[55]*y[56]*y[65]*y[72]*y[73]+x1*x2*y[1]*y[3]*y[9]\
*y[41]*y[55]*y[56]*y[65]*y[72]*y[73]-x0*y[1]*y[2]*y[4]*y[74]*y[75]*y[76]-x1\
*y[1]*y[2]*y[3]*y[4]*y[74]*y[75]*y[76])+pow(lambda,4)*(y[1]*y[3]*y[7]*y[8]*\
y[13]*y[34]*y[36]*y[57]*y[58]*y[59]+x2*y[1]*y[3]*y[7]*y[8]*y[13]*y[34]*y[36\
]*y[57]*y[58]*y[59]+x1*x2*y[1]*y[2]*y[3]*y[13]*y[34]*y[36]*y[41]*y[55]*y[56\
]*y[65]*y[72]*y[73]+2.*x1*y[1]*y[2]*y[3]*y[4]*y[13]*y[34]*y[36]*y[41]*y[55]\
*y[56]*y[65]*y[72]*y[73]+x1*x2*y[1]*y[3]*y[5]*y[13]*y[34]*y[36]*y[41]*y[55]\
*y[56]*y[65]*y[72]*y[73]+2.*x2*y[1]*y[3]*y[7]*y[8]*y[13]*y[34]*y[36]*y[41]*\
y[55]*y[56]*y[65]*y[72]*y[73]-x1*x2*y[1]*y[3]*y[9]*y[13]*y[34]*y[36]*y[41]*\
y[55]*y[56]*y[65]*y[72]*y[73]+2.*x2*y[1]*y[3]*y[7]*y[8]*y[12]*y[33]*y[35]*y\
[57]*y[58]*y[59]*y[65]*y[72]*y[73]+x1*y[1]*y[2]*y[3]*y[4]*y[13]*y[34]*y[36]\
*y[74]*y[75]*y[76]+2.*x1*y[1]*y[2]*y[3]*y[4]*y[12]*y[33]*y[35]*y[41]*y[55]*\
y[56]*y[74]*y[75]*y[76])));
return (FOUT);
}
double Pt8t1(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
double y[9];
double FOUT;
y[1]=1./bi;
y[2]=em[0];
y[3]=x2*x2;
y[4]=em[1];
y[5]=em[2];
y[6]=em[3];
y[7]=x1*x1;
y[8]=esx[0];
FOUT=(1.-x0)*x0*(x2*y[1]*y[2]+x1*x2*y[1]*y[2]+2.*x0*x1*x2*y[1]*y[2]+y[1]*y[2\
]*y[3]+2.*x0*x1*y[1]*y[2]*y[3]+y[1]*y[4]+x1*y[1]*y[4]+2.*x0*x1*y[1]*y[4]+x2\
*y[1]*y[4]+2.*x0*x1*x2*y[1]*y[4]+x1*y[1]*y[5]+x1*x2*y[1]*y[5]+x1*y[1]*y[6]+\
x1*x2*y[1]*y[6]+y[1]*y[6]*y[7]+2.*x0*y[1]*y[6]*y[7]+2.*x0*x2*y[1]*y[6]*y[7]\
-x1*y[1]*y[8]-x2*y[1]*y[8]-2.*x0*x1*x2*y[1]*y[8]);
return (FOUT);
}
double Pt8t2(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
double y[8];
double FOUT;
y[1]=1./bi;
y[2]=em[0];
y[3]=x0*x0;
y[4]=em[1];
y[5]=em[2];
y[6]=em[3];
y[7]=esx[0];
FOUT=(1.-x1)*x1*(x0*x2*y[1]*y[2]+x2*y[1]*y[2]*y[3]+x2*x2*y[1]*y[2]*y[3]+x0*y\
[1]*y[4]+y[1]*y[3]*y[4]+x2*y[1]*y[3]*y[4]+y[1]*y[5]+x0*y[1]*y[5]+x0*x2*y[1]\
*y[5]+x0*y[1]*y[6]+2.*x0*x1*y[1]*y[6]+x0*x2*y[1]*y[6]+2.*x1*y[1]*y[3]*y[6]+\
2.*x1*x2*y[1]*y[3]*y[6]-x0*y[1]*y[7]-x2*y[1]*y[3]*y[7]);
return (FOUT);
}
double Pt8t3(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
double y[8];
double FOUT;
y[1]=1./bi;
y[2]=em[0];
y[3]=x0*x0;
y[4]=em[1];
y[5]=em[2];
y[6]=em[3];
y[7]=esx[0];
FOUT=(1.-x2)*x2*(x0*y[1]*y[2]+x0*x1*y[1]*y[2]+2.*x0*x2*y[1]*y[2]+x1*y[1]*y[2\
]*y[3]+2.*x1*x2*y[1]*y[2]*y[3]+x0*y[1]*y[4]+x1*y[1]*y[3]*y[4]+y[1]*y[5]+x0*\
x1*y[1]*y[5]+x0*x1*y[1]*y[6]+x1*x1*y[1]*y[3]*y[6]-x0*y[1]*y[7]-x1*y[1]*y[3]\
*y[7]);
return (FOUT);
}
