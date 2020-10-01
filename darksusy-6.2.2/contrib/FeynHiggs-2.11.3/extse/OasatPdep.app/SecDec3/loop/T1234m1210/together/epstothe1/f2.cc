#include "intfile.hh"

dcmplx Pf2(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
dcmplx y[139];
dcmplx FOUT;
dcmplx MYI(0.,1.);
y[1]=1./bi;
y[2]=em[0];
y[3]=-x0;
y[4]=1.+y[3];
y[5]=em[1];
y[6]=y[1]*y[5];
y[7]=esx[0];
y[8]=lrs[0];
y[9]=x1*y[1]*y[2];
y[10]=x1*x1;
y[11]=y[1]*y[2]*y[10];
y[12]=x1*y[1]*y[5];
y[13]=-(x1*y[1]*y[7]);
y[14]=y[6]+y[9]+y[11]+y[12]+y[13];
y[15]=-x1;
y[16]=1.+y[15];
y[17]=lrs[1];
y[18]=y[1]*y[2];
y[19]=2.*x1*y[1]*y[2];
y[20]=-(y[1]*y[7]);
y[21]=x0*y[1]*y[2];
y[22]=2.*x0*x1*y[1]*y[2];
y[23]=x0*y[1]*y[5];
y[24]=-(x0*y[1]*y[7]);
y[25]=y[6]+y[18]+y[19]+y[20]+y[21]+y[22]+y[23]+y[24];
y[26]=-(lambda*MYI*x0*y[4]*y[8]*y[14]);
y[27]=x0+y[26];
y[28]=-(lambda*MYI*x1*y[16]*y[17]*y[25]);
y[29]=x1+y[28];
y[30]=1./x2;
y[31]=pow(bi,-2);
y[32]=lambda*lambda;
y[33]=y[6]+y[18]+y[19]+y[20];
y[34]=pow(y[33],2);
y[35]=x0*x1*y[4]*y[8]*y[16]*y[17]*y[32]*y[34];
y[36]=-(lambda*MYI*y[4]*y[8]*y[14]);
y[37]=lambda*MYI*x0*y[8]*y[14];
y[38]=1.+y[36]+y[37];
y[39]=2.*y[1]*y[2];
y[40]=2.*x0*y[1]*y[2];
y[41]=y[39]+y[40];
y[42]=-(lambda*MYI*x1*y[16]*y[17]*y[41]);
y[43]=-(lambda*MYI*y[16]*y[17]*y[25]);
y[44]=lambda*MYI*x1*y[17]*y[25];
y[45]=1.+y[42]+y[43]+y[44];
y[46]=y[38]*y[45];
y[47]=y[35]+y[46];
y[48]=y[1]*y[27];
y[49]=y[1]*y[29];
y[50]=y[1]*y[27]*y[29];
y[51]=y[1]+y[48]+y[49]+y[50];
y[52]=pow(y[51],-2);
y[53]=pow(y[29],2);
y[54]=x0*x0;
y[55]=myLog(bi);
y[56]=y[1]*y[2]*y[54];
y[57]=x1*y[1]*y[2]*y[54];
y[58]=y[21]+y[22]+y[23]+y[24]+y[56]+y[57];
y[59]=lrs[2];
y[60]=-(lambda*MYI*y[58]*y[59]);
y[61]=1.+y[60];
y[62]=myLog(y[61]);
y[63]=myLog(y[51]);
y[64]=y[1]*y[5]*y[27];
y[65]=y[1]*y[2]*y[29];
y[66]=y[1]*y[5]*y[29];
y[67]=-(y[1]*y[7]*y[29]);
y[68]=y[1]*y[2]*y[27]*y[29];
y[69]=y[1]*y[5]*y[27]*y[29];
y[70]=-(y[1]*y[7]*y[27]*y[29]);
y[71]=y[1]*y[2]*y[53];
y[72]=y[1]*y[2]*y[27]*y[53];
y[73]=y[6]+y[64]+y[65]+y[66]+y[67]+y[68]+y[69]+y[70]+y[71]+y[72];
y[74]=myLog(y[73]);
y[75]=myLog(x2);
y[76]=-x2;
y[77]=1.+y[76];
y[78]=y[40]+y[56];
y[79]=2.*x2*y[1]*y[2];
y[80]=2.*x0*x2*y[1]*y[2];
y[81]=2.*x1*x2*y[1]*y[2];
y[82]=x2*x2;
y[83]=x2*y[1]*y[2];
y[84]=2.*x0*x1*x2*y[1]*y[2];
y[85]=2.*x0*y[1]*y[2]*y[82];
y[86]=x2*y[1]*y[5];
y[87]=-(x2*y[1]*y[7]);
y[88]=y[6]+y[9]+y[11]+y[12]+y[13]+y[80]+y[81]+y[83]+y[84]+y[85]+y[86]+y[87];
y[89]=4.*x0*x2*y[1]*y[2];
y[90]=y[6]+y[18]+y[19]+y[20]+y[22]+y[40]+y[89];
y[91]=y[6]+y[18]+y[19]+y[20]+y[79]+y[80];
y[92]=x2*y[1]*y[2]*y[54];
y[93]=y[6]+y[18]+y[19]+y[20]+y[21]+y[22]+y[23]+y[24]+y[80]+y[92];
y[94]=2.*y[1]*y[2]*y[82];
y[95]=y[79]+y[81]+y[94];
y[96]=-(lambda*MYI*x0*y[4]*y[8]*y[95]);
y[97]=-(lambda*MYI*y[4]*y[8]*y[88]);
y[98]=lambda*MYI*x0*y[8]*y[88];
y[99]=1.+y[96]+y[97]+y[98];
y[100]=-(lambda*MYI*y[16]*y[17]*y[93]);
y[101]=lambda*MYI*x1*y[17]*y[93];
y[102]=1.+y[42]+y[100]+y[101];
y[103]=2.*x2*y[1]*y[2]*y[54];
y[104]=y[21]+y[22]+y[23]+y[24]+y[56]+y[57]+y[103];
y[105]=-(lambda*MYI*y[59]*y[77]*y[104]);
y[106]=-(lambda*MYI*x0*y[4]*y[8]*y[88]);
y[107]=x0+y[106];
y[108]=-(lambda*MYI*x1*y[16]*y[17]*y[93]);
y[109]=x1+y[108];
y[110]=1.+y[105];
y[111]=1./y[110];
y[112]=x0*x1*y[4]*y[8]*y[16]*y[17]*y[32]*y[90]*y[91];
y[113]=-(lambda*MYI*x1*y[16]*y[17]*y[78]*y[99]);
y[114]=y[112]+y[113];
y[115]=lambda*MYI*x2*y[59]*y[77]*y[78]*y[114];
y[116]=-(x0*x1*y[4]*y[8]*y[16]*y[17]*y[32]*y[78]*y[91]);
y[117]=lambda*MYI*x0*y[4]*y[8]*y[90]*y[102];
y[118]=y[116]+y[117];
y[119]=-(lambda*MYI*x2*y[59]*y[77]*y[90]*y[118]);
y[120]=pow(y[91],2);
y[121]=x0*x1*y[4]*y[8]*y[16]*y[17]*y[32]*y[120];
y[122]=y[99]*y[102];
y[123]=y[121]+y[122];
y[124]=-2.*lambda*MYI*x2*y[1]*y[2]*y[54]*y[59]*y[77];
y[125]=lambda*MYI*x2*y[59]*y[104];
y[126]=1.+y[105]+y[124]+y[125];
y[127]=y[123]*y[126];
y[128]=y[115]+y[119]+y[127];
y[129]=y[1]*y[107];
y[130]=y[1]*y[109];
y[131]=y[1]*y[107]*y[109];
y[132]=-(lambda*MYI*x2*y[59]*y[77]*y[104]);
y[133]=x2+y[132];
y[134]=y[1]*y[107]*y[133];
y[135]=y[1]+y[129]+y[130]+y[131]+y[134];
y[136]=pow(y[135],-2);
y[137]=pow(y[109],2);
y[138]=pow(y[107],2);
FOUT=-(y[30]*(y[31]*y[47]*y[52]*y[55]+y[31]*y[47]*y[52]*y[62]+3.*y[31]*y[47]\
*y[52]*y[63]-2.*y[31]*y[47]*y[52]*y[74]))+0.5*(y[47]*y[52]*(pow(y[55],2)*y[\
31]+pow(y[62],2)*y[31]+2.*y[31]*y[55]*y[62])+2.*(y[31]*y[47]*y[55]+y[31]*y[\
47]*y[62])*(3.*y[52]*y[63]-2.*y[52]*y[74])+y[31]*y[47]*(9.*pow(y[63],2)*y[5\
2]+4.*pow(y[74],2)*y[52]-12.*y[52]*y[63]*y[74]))-y[30]*y[31]*y[47]*y[52]*y[\
75]+y[30]*y[31]*y[75]*y[111]*y[128]*y[136]+y[30]*(myLog(y[110])*y[31]*y[111\
]*y[128]*y[136]+3.*myLog(y[135])*y[31]*y[111]*y[128]*y[136]-2.*myLog(y[6]+y\
[1]*y[5]*y[107]+y[1]*y[2]*y[109]+y[1]*y[5]*y[109]-y[1]*y[7]*y[109]+y[1]*y[2\
]*y[107]*y[109]+y[1]*y[5]*y[107]*y[109]-y[1]*y[7]*y[107]*y[109]+y[1]*y[2]*y\
[107]*y[133]+y[1]*y[5]*y[107]*y[133]-y[1]*y[7]*y[107]*y[133]+2.*y[1]*y[2]*y\
[107]*y[109]*y[133]+y[1]*y[2]*y[137]+y[1]*y[2]*y[107]*y[137]+pow(y[133],2)*\
y[1]*y[2]*y[138]+y[1]*y[2]*y[133]*y[138]+y[1]*y[2]*y[109]*y[133]*y[138])*y[\
31]*y[111]*y[128]*y[136]+y[31]*y[55]*y[111]*y[128]*y[136]);
return (FOUT);
}
