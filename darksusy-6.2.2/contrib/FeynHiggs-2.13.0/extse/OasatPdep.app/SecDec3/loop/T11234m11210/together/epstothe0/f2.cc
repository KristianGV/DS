#include "intfile.hh"

dcmplx Pf2(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
dcmplx y[106];
dcmplx FOUT;
dcmplx MYI(0.,1.);
y[1]=1./bi;
y[2]=em[1];
y[3]=-x0;
y[4]=1.+y[3];
y[5]=lrs[0];
y[6]=em[0];
y[7]=y[1]*y[6];
y[8]=y[1]*y[2];
y[9]=2.*x0*y[1]*y[2];
y[10]=esx[0];
y[11]=-(y[1]*y[10]);
y[12]=x1*y[1]*y[6];
y[13]=x1*y[1]*y[2];
y[14]=2.*x0*x1*y[1]*y[2];
y[15]=-(x1*y[1]*y[10]);
y[16]=y[7]+y[8]+y[9]+y[11]+y[12]+y[13]+y[14]+y[15];
y[17]=-x1;
y[18]=1.+y[17];
y[19]=lrs[1];
y[20]=x0*y[1]*y[6];
y[21]=x0*y[1]*y[2];
y[22]=x0*x0;
y[23]=y[1]*y[2]*y[22];
y[24]=-(x0*y[1]*y[10]);
y[25]=y[7]+y[20]+y[21]+y[23]+y[24];
y[26]=-(lambda*MYI*x0*y[4]*y[5]*y[16]);
y[27]=x0+y[26];
y[28]=-(lambda*MYI*x1*y[18]*y[19]*y[25]);
y[29]=x1+y[28];
y[30]=pow(y[27],2);
y[31]=pow(bi,-2);
y[32]=lambda*lambda;
y[33]=y[7]+y[8]+y[9]+y[11];
y[34]=pow(y[33],2);
y[35]=x0*x1*y[4]*y[5]*y[18]*y[19]*y[32]*y[34];
y[36]=2.*y[1]*y[2];
y[37]=2.*x1*y[1]*y[2];
y[38]=y[36]+y[37];
y[39]=-(lambda*MYI*x0*y[4]*y[5]*y[38]);
y[40]=-(lambda*MYI*y[4]*y[5]*y[16]);
y[41]=lambda*MYI*x0*y[5]*y[16];
y[42]=1.+y[39]+y[40]+y[41];
y[43]=-(lambda*MYI*y[18]*y[19]*y[25]);
y[44]=lambda*MYI*x1*y[19]*y[25];
y[45]=1.+y[43]+y[44];
y[46]=y[42]*y[45];
y[47]=y[35]+y[46];
y[48]=y[1]*y[27];
y[49]=y[1]*y[29];
y[50]=y[1]*y[27]*y[29];
y[51]=y[1]+y[48]+y[49]+y[50];
y[52]=1./y[51];
y[53]=y[1]*y[6]*y[27];
y[54]=y[1]*y[2]*y[27];
y[55]=-(y[1]*y[10]*y[27]);
y[56]=y[1]*y[2]*y[30];
y[57]=y[1]*y[6]*y[29];
y[58]=y[1]*y[6]*y[27]*y[29];
y[59]=y[1]*y[2]*y[27]*y[29];
y[60]=-(y[1]*y[10]*y[27]*y[29]);
y[61]=y[1]*y[2]*y[29]*y[30];
y[62]=y[7]+y[53]+y[54]+y[55]+y[56]+y[57]+y[58]+y[59]+y[60]+y[61];
y[63]=1./y[62];
y[64]=1./x2;
y[65]=2.*x1*y[1]*y[6];
y[66]=x0*x1*y[1]*y[6];
y[67]=x0*x1*y[1]*y[2];
y[68]=-(x0*x1*y[1]*y[10]);
y[69]=lrs[2];
y[70]=-x2;
y[71]=1.+y[70];
y[72]=2.*y[1]*y[6];
y[73]=2.*x2*y[1]*y[6];
y[74]=y[20]+y[21]+y[24]+y[72]+y[73];
y[75]=x2*y[1]*y[6];
y[76]=x1*x2*y[1]*y[6];
y[77]=x1*x2*y[1]*y[2];
y[78]=-(x1*x2*y[1]*y[10]);
y[79]=y[7]+y[8]+y[9]+y[11]+y[12]+y[13]+y[14]+y[15]+y[75]+y[76]+y[77]+y[78];
y[80]=y[7]+y[12]+y[13]+y[15];
y[81]=x2*y[1]*y[2];
y[82]=-(x2*y[1]*y[10]);
y[83]=y[7]+y[8]+y[9]+y[11]+y[75]+y[81]+y[82];
y[84]=x0*x2*y[1]*y[6];
y[85]=x2*x2;
y[86]=y[1]*y[6]*y[85];
y[87]=x0*x2*y[1]*y[2];
y[88]=-(x0*x2*y[1]*y[10]);
y[89]=y[7]+y[20]+y[21]+y[23]+y[24]+y[73]+y[84]+y[86]+y[87]+y[88];
y[90]=-(lambda*MYI*y[4]*y[5]*y[79]);
y[91]=lambda*MYI*x0*y[5]*y[79];
y[92]=1.+y[39]+y[90]+y[91];
y[93]=-(lambda*MYI*y[18]*y[19]*y[89]);
y[94]=lambda*MYI*x1*y[19]*y[89];
y[95]=1.+y[93]+y[94];
y[96]=2.*x1*x2*y[1]*y[6];
y[97]=y[7]+y[20]+y[65]+y[66]+y[67]+y[68]+y[96];
y[98]=-(lambda*MYI*y[69]*y[71]*y[97]);
y[99]=-(lambda*MYI*x0*y[4]*y[5]*y[79]);
y[100]=x0+y[99];
y[101]=-(lambda*MYI*x1*y[18]*y[19]*y[89]);
y[102]=x1+y[101];
y[103]=pow(y[100],2);
y[104]=-(lambda*MYI*x2*y[69]*y[71]*y[97]);
y[105]=x2+y[104];
FOUT=myLog(bi)*y[31]*y[47]*y[52]*y[63]+3.*myLog(y[51])*y[31]*y[47]*y[52]*y[6\
3]-2.*myLog(y[62])*y[31]*y[47]*y[52]*y[63]+myLog(1.-lambda*MYI*(y[7]+y[20]+\
y[65]+y[66]+y[67]+y[68])*y[69])*y[31]*y[47]*y[52]*y[63]-y[31]*y[47]*y[52]*y\
[63]*y[64]+(y[31]*y[64]*(lambda*MYI*x2*y[69]*y[71]*y[74]*(x0*x1*y[4]*y[5]*y\
[18]*y[19]*y[32]*y[80]*y[83]-lambda*MYI*x1*y[18]*y[19]*y[74]*y[92])-lambda*\
MYI*x2*y[69]*y[71]*y[80]*(-(x0*x1*y[4]*y[5]*y[18]*y[19]*y[32]*y[74]*y[83])+\
lambda*MYI*x0*y[4]*y[5]*y[80]*y[95])+(x0*x1*pow(y[83],2)*y[4]*y[5]*y[18]*y[\
19]*y[32]+y[92]*y[95])*(1.-2.*lambda*MYI*x1*x2*y[1]*y[6]*y[69]*y[71]+lambda\
*MYI*x2*y[69]*y[97]+y[98])))/((1.+y[98])*(y[1]+y[1]*y[100]+y[1]*y[102]+y[1]\
*y[100]*y[102]+y[1]*y[102]*y[105])*(y[7]+y[1]*y[2]*y[100]+y[1]*y[6]*y[100]-\
y[1]*y[10]*y[100]+y[1]*y[6]*y[102]+pow(y[105],2)*y[1]*y[6]*y[102]+y[1]*y[2]\
*y[100]*y[102]+y[1]*y[6]*y[100]*y[102]-y[1]*y[10]*y[100]*y[102]+y[1]*y[2]*y\
[103]+y[1]*y[2]*y[102]*y[103]+y[1]*y[6]*y[105]+y[1]*y[6]*y[100]*y[105]+2.*y\
[1]*y[6]*y[102]*y[105]+y[1]*y[2]*y[100]*y[102]*y[105]+y[1]*y[6]*y[100]*y[10\
2]*y[105]-y[1]*y[10]*y[100]*y[102]*y[105]));
return (FOUT);
}