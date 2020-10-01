#include "intfile.hh"

dcmplx Pf1(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
dcmplx y[167];
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
y[54]=em[2];
y[55]=x0*x0;
y[56]=em[3];
y[57]=myLog(bi);
y[58]=x0*x1*y[1]*y[2];
y[59]=x0*y[1]*y[54];
y[60]=y[1]*y[54]*y[55];
y[61]=x0*x1*y[1]*y[54];
y[62]=x1*y[1]*y[54]*y[55];
y[63]=y[1]*y[56];
y[64]=x0*y[1]*y[56];
y[65]=x1*y[1]*y[56];
y[66]=x0*x1*y[1]*y[56];
y[67]=y[23]+y[24]+y[58]+y[59]+y[60]+y[61]+y[62]+y[63]+y[64]+y[65]+y[66];
y[68]=lrs[2];
y[69]=-(lambda*MYI*y[67]*y[68]);
y[70]=1.+y[69];
y[71]=myLog(y[70]);
y[72]=myLog(y[51]);
y[73]=y[1]*y[5]*y[27];
y[74]=y[1]*y[2]*y[29];
y[75]=y[1]*y[5]*y[29];
y[76]=-(y[1]*y[7]*y[29]);
y[77]=y[1]*y[2]*y[27]*y[29];
y[78]=y[1]*y[5]*y[27]*y[29];
y[79]=-(y[1]*y[7]*y[27]*y[29]);
y[80]=y[1]*y[2]*y[53];
y[81]=y[1]*y[2]*y[27]*y[53];
y[82]=y[6]+y[73]+y[74]+y[75]+y[76]+y[77]+y[78]+y[79]+y[80]+y[81];
y[83]=myLog(y[82]);
y[84]=myLog(x2);
y[85]=-x2;
y[86]=1.+y[85];
y[87]=y[21]+y[59]+y[60]+y[63]+y[64];
y[88]=x2*y[1]*y[54];
y[89]=2.*x0*x2*y[1]*y[54];
y[90]=x2*x2;
y[91]=x2*y[1]*y[56];
y[92]=x1*x2*y[1]*y[2];
y[93]=x2*y[1]*y[5];
y[94]=x1*x2*y[1]*y[54];
y[95]=2.*x0*x1*x2*y[1]*y[54];
y[96]=2.*x0*y[1]*y[54]*y[90];
y[97]=x1*x2*y[1]*y[56];
y[98]=y[1]*y[56]*y[90];
y[99]=-(x2*y[1]*y[7]);
y[100]=y[6]+y[9]+y[11]+y[12]+y[13]+y[88]+y[89]+y[91]+y[92]+y[93]+y[94]+y[95]\
+y[96]+y[97]+y[98]+y[99];
y[101]=y[1]*y[54];
y[102]=2.*x0*y[1]*y[54];
y[103]=x1*y[1]*y[54];
y[104]=2.*x0*x1*y[1]*y[54];
y[105]=4.*x0*x2*y[1]*y[54];
y[106]=2.*x2*y[1]*y[56];
y[107]=y[6]+y[9]+y[20]+y[63]+y[65]+y[101]+y[102]+y[103]+y[104]+y[105]+y[106]\
;
y[108]=x2*y[1]*y[2];
y[109]=y[6]+y[18]+y[19]+y[20]+y[88]+y[89]+y[91]+y[108];
y[110]=x0*x2*y[1]*y[2];
y[111]=x0*x2*y[1]*y[54];
y[112]=x2*y[1]*y[54]*y[55];
y[113]=x0*x2*y[1]*y[56];
y[114]=y[6]+y[18]+y[19]+y[20]+y[21]+y[22]+y[23]+y[24]+y[91]+y[110]+y[111]+y[\
112]+y[113];
y[115]=2.*x2*y[1]*y[54];
y[116]=2.*x1*x2*y[1]*y[54];
y[117]=2.*y[1]*y[54]*y[90];
y[118]=y[115]+y[116]+y[117];
y[119]=-(lambda*MYI*x0*y[4]*y[8]*y[118]);
y[120]=-(lambda*MYI*y[4]*y[8]*y[100]);
y[121]=lambda*MYI*x0*y[8]*y[100];
y[122]=1.+y[119]+y[120]+y[121];
y[123]=-(lambda*MYI*y[16]*y[17]*y[114]);
y[124]=lambda*MYI*x1*y[17]*y[114];
y[125]=1.+y[42]+y[123]+y[124];
y[126]=2.*x2*y[1]*y[54]*y[55];
y[127]=2.*x0*x2*y[1]*y[56];
y[128]=y[23]+y[24]+y[58]+y[59]+y[60]+y[61]+y[62]+y[63]+y[64]+y[65]+y[66]+y[1\
26]+y[127];
y[129]=-(lambda*MYI*y[68]*y[86]*y[128]);
y[130]=-(lambda*MYI*x0*y[4]*y[8]*y[100]);
y[131]=x0+y[130];
y[132]=-(lambda*MYI*x1*y[16]*y[17]*y[114]);
y[133]=x1+y[132];
y[134]=1.+y[129];
y[135]=1./y[134];
y[136]=x0*x1*y[4]*y[8]*y[16]*y[17]*y[32]*y[107]*y[109];
y[137]=-(lambda*MYI*x1*y[16]*y[17]*y[87]*y[122]);
y[138]=y[136]+y[137];
y[139]=lambda*MYI*x2*y[68]*y[86]*y[87]*y[138];
y[140]=-(x0*x1*y[4]*y[8]*y[16]*y[17]*y[32]*y[87]*y[109]);
y[141]=lambda*MYI*x0*y[4]*y[8]*y[107]*y[125];
y[142]=y[140]+y[141];
y[143]=-(lambda*MYI*x2*y[68]*y[86]*y[107]*y[142]);
y[144]=pow(y[109],2);
y[145]=x0*x1*y[4]*y[8]*y[16]*y[17]*y[32]*y[144];
y[146]=y[122]*y[125];
y[147]=y[145]+y[146];
y[148]=2.*y[1]*y[54]*y[55];
y[149]=2.*x0*y[1]*y[56];
y[150]=y[148]+y[149];
y[151]=-(lambda*MYI*x2*y[68]*y[86]*y[150]);
y[152]=lambda*MYI*x2*y[68]*y[128];
y[153]=1.+y[129]+y[151]+y[152];
y[154]=y[147]*y[153];
y[155]=y[139]+y[143]+y[154];
y[156]=y[1]*y[131];
y[157]=y[1]*y[133];
y[158]=y[1]*y[131]*y[133];
y[159]=-(lambda*MYI*x2*y[68]*y[86]*y[128]);
y[160]=x2+y[159];
y[161]=y[1]*y[131]*y[160];
y[162]=y[1]+y[156]+y[157]+y[158]+y[161];
y[163]=pow(y[162],-2);
y[164]=pow(y[133],2);
y[165]=pow(y[131],2);
y[166]=pow(y[160],2);
FOUT=-(y[30]*(y[31]*y[47]*y[52]*y[57]+y[31]*y[47]*y[52]*y[71]+3.*y[31]*y[47]\
*y[52]*y[72]-2.*y[31]*y[47]*y[52]*y[83]))+0.5*(y[47]*y[52]*(pow(y[57],2)*y[\
31]+pow(y[71],2)*y[31]+2.*y[31]*y[57]*y[71])+2.*(y[31]*y[47]*y[57]+y[31]*y[\
47]*y[71])*(3.*y[52]*y[72]-2.*y[52]*y[83])+y[31]*y[47]*(9.*pow(y[72],2)*y[5\
2]+4.*pow(y[83],2)*y[52]-12.*y[52]*y[72]*y[83]))-y[30]*y[31]*y[47]*y[52]*y[\
84]+y[30]*y[31]*y[84]*y[135]*y[155]*y[163]+y[30]*(myLog(y[134])*y[31]*y[135\
]*y[155]*y[163]+3.*myLog(y[162])*y[31]*y[135]*y[155]*y[163]-2.*myLog(y[6]+y\
[1]*y[5]*y[131]+y[1]*y[2]*y[133]+y[1]*y[5]*y[133]-y[1]*y[7]*y[133]+y[1]*y[2\
]*y[131]*y[133]+y[1]*y[5]*y[131]*y[133]-y[1]*y[7]*y[131]*y[133]+y[1]*y[56]*\
y[160]+y[1]*y[5]*y[131]*y[160]-y[1]*y[7]*y[131]*y[160]+y[1]*y[54]*y[131]*y[\
160]+y[1]*y[56]*y[131]*y[160]+y[1]*y[56]*y[133]*y[160]+y[1]*y[2]*y[131]*y[1\
33]*y[160]+y[1]*y[54]*y[131]*y[133]*y[160]+y[1]*y[56]*y[131]*y[133]*y[160]+\
y[1]*y[2]*y[164]+y[1]*y[2]*y[131]*y[164]+y[1]*y[54]*y[160]*y[165]+y[1]*y[54\
]*y[133]*y[160]*y[165]+y[1]*y[56]*y[131]*y[166]+y[1]*y[54]*y[165]*y[166])*y\
[31]*y[135]*y[155]*y[163]+y[31]*y[57]*y[135]*y[155]*y[163]);
return (FOUT);
}
