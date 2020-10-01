#include "intfile.hh"

dcmplx Pf3(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
dcmplx y[97];
dcmplx FOUT;
dcmplx MYI(0.,1.);
y[1]=1./bi;
y[2]=em[0];
y[3]=em[1];
y[4]=esx[0];
y[5]=-x0;
y[6]=1.+y[5];
y[7]=y[1]*y[2];
y[8]=2.*x0*y[1]*y[2];
y[9]=y[1]*y[3];
y[10]=-(y[1]*y[4]);
y[11]=lrs[0];
y[12]=x1*y[1]*y[2];
y[13]=2.*x0*x1*y[1]*y[2];
y[14]=x1*y[1]*y[3];
y[15]=-(x1*y[1]*y[4]);
y[16]=y[7]+y[8]+y[9]+y[10]+y[12]+y[13]+y[14]+y[15];
y[17]=-(lambda*MYI*y[6]*y[11]*y[16]);
y[18]=-x1;
y[19]=1.+y[18];
y[20]=lrs[1];
y[21]=x0*y[1]*y[2];
y[22]=x0*x0;
y[23]=y[1]*y[2]*y[22];
y[24]=x0*y[1]*y[3];
y[25]=-(x0*y[1]*y[4]);
y[26]=y[9]+y[21]+y[23]+y[24]+y[25];
y[27]=-(lambda*MYI*x0*y[6]*y[11]*y[16]);
y[28]=x0+y[27];
y[29]=-(lambda*MYI*x1*y[19]*y[20]*y[26]);
y[30]=x1+y[29];
y[31]=pow(y[28],2);
y[32]=pow(bi,-2);
y[33]=1.+y[17];
y[34]=lambda*lambda;
y[35]=y[7]+y[8]+y[9]+y[10];
y[36]=pow(y[35],2);
y[37]=x0*x1*y[6]*y[11]*y[19]*y[20]*y[34]*y[36];
y[38]=2.*y[1]*y[2];
y[39]=2.*x1*y[1]*y[2];
y[40]=y[38]+y[39];
y[41]=-(lambda*MYI*x0*y[6]*y[11]*y[40]);
y[42]=lambda*MYI*x0*y[11]*y[16];
y[43]=1.+y[17]+y[41]+y[42];
y[44]=-(lambda*MYI*y[19]*y[20]*y[26]);
y[45]=lambda*MYI*x1*y[20]*y[26];
y[46]=1.+y[44]+y[45];
y[47]=y[43]*y[46];
y[48]=y[37]+y[47];
y[49]=y[1]*y[28];
y[50]=y[1]*y[30];
y[51]=y[1]*y[28]*y[30];
y[52]=y[1]+y[49]+y[50]+y[51];
y[53]=1./y[52];
y[54]=y[1]*y[2]*y[28];
y[55]=y[1]*y[3]*y[28];
y[56]=-(y[1]*y[4]*y[28]);
y[57]=y[1]*y[2]*y[31];
y[58]=y[1]*y[3]*y[30];
y[59]=y[1]*y[2]*y[28]*y[30];
y[60]=y[1]*y[3]*y[28]*y[30];
y[61]=-(y[1]*y[4]*y[28]*y[30]);
y[62]=y[1]*y[2]*y[30]*y[31];
y[63]=y[9]+y[54]+y[55]+y[56]+y[57]+y[58]+y[59]+y[60]+y[61]+y[62];
y[64]=1./y[63];
y[65]=1./x2;
y[66]=2.*x1*x2*y[1]*y[2];
y[67]=lrs[2];
y[68]=-x2;
y[69]=1.+y[68];
y[70]=2.*x2*y[1]*y[2];
y[71]=y[7]+y[8]+y[9]+y[10]+y[70];
y[72]=x2*y[1]*y[2];
y[73]=y[7]+y[8]+y[9]+y[10]+y[12]+y[13]+y[14]+y[15]+y[66]+y[72];
y[74]=-(lambda*MYI*y[6]*y[11]*y[73]);
y[75]=y[7]+y[39];
y[76]=2.*x0*x2*y[1]*y[2];
y[77]=x2*x2;
y[78]=y[1]*y[2]*y[77];
y[79]=x2*y[1]*y[3];
y[80]=-(x2*y[1]*y[4]);
y[81]=y[9]+y[21]+y[23]+y[24]+y[25]+y[72]+y[76]+y[78]+y[79]+y[80];
y[82]=pow(y[71],2);
y[83]=lambda*MYI*x0*y[11]*y[73];
y[84]=1.+y[41]+y[74]+y[83];
y[85]=-(lambda*MYI*y[19]*y[20]*y[81]);
y[86]=lambda*MYI*x1*y[20]*y[81];
y[87]=1.+y[85]+y[86];
y[88]=y[7]+y[12]+y[13]+y[14]+y[15]+y[21]+y[66];
y[89]=-(lambda*MYI*y[67]*y[69]*y[88]);
y[90]=-(lambda*MYI*x0*y[6]*y[11]*y[73]);
y[91]=x0+y[90];
y[92]=-(lambda*MYI*x1*y[19]*y[20]*y[81]);
y[93]=x1+y[92];
y[94]=pow(y[91],2);
y[95]=-(lambda*MYI*x2*y[67]*y[69]*y[88]);
y[96]=x2+y[95];
FOUT=x0*(myLog(bi)*y[32]*y[33]*y[48]*y[53]*y[64]+3.*myLog(y[52])*y[32]*y[33]\
*y[48]*y[53]*y[64]-2.*myLog(y[63])*y[32]*y[33]*y[48]*y[53]*y[64]+myLog(1.-l\
ambda*MYI*(y[7]+y[12]+y[13]+y[14]+y[15]+y[21])*y[67])*y[32]*y[33]*y[48]*y[5\
3]*y[64])-x0*y[32]*y[33]*y[48]*y[53]*y[64]*y[65]+(x0*y[32]*y[65]*(1.+y[74])\
*(lambda*MYI*x2*y[67]*y[69]*y[71]*(x0*x1*y[6]*y[11]*y[19]*y[20]*y[34]*y[71]\
*y[75]-lambda*MYI*x1*y[19]*y[20]*y[71]*y[84])-lambda*MYI*x2*y[67]*y[69]*y[7\
5]*(-(x0*x1*y[6]*y[11]*y[19]*y[20]*y[34]*y[82])+lambda*MYI*x0*y[6]*y[11]*y[\
75]*y[87])+(x0*x1*y[6]*y[11]*y[19]*y[20]*y[34]*y[82]+y[84]*y[87])*(1.-2.*la\
mbda*MYI*x1*x2*y[1]*y[2]*y[67]*y[69]+lambda*MYI*x2*y[67]*y[88]+y[89])))/((1\
.+y[89])*(y[1]+y[1]*y[91]+y[1]*y[93]+y[1]*y[91]*y[93]+y[1]*y[93]*y[96])*(y[\
9]+y[1]*y[2]*y[91]+y[1]*y[3]*y[91]-y[1]*y[4]*y[91]+pow(y[96],2)*y[1]*y[2]*y\
[93]+y[1]*y[3]*y[93]+y[1]*y[2]*y[91]*y[93]+y[1]*y[3]*y[91]*y[93]-y[1]*y[4]*\
y[91]*y[93]+y[1]*y[2]*y[94]+y[1]*y[2]*y[93]*y[94]+y[1]*y[2]*y[96]+y[1]*y[2]\
*y[91]*y[96]+y[1]*y[2]*y[93]*y[96]+y[1]*y[3]*y[93]*y[96]-y[1]*y[4]*y[93]*y[\
96]+2.*y[1]*y[2]*y[91]*y[93]*y[96]));
return (FOUT);
}