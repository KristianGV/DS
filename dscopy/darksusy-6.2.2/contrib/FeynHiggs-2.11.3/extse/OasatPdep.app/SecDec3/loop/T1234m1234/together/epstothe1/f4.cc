#include "intfile.hh"

dcmplx Pf4(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
dcmplx y[170];
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
y[30]=1./x2;
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
y[52]=pow(y[51],-2);
y[53]=pow(y[27],2);
y[54]=em[2];
y[55]=em[3];
y[56]=x1*x1;
y[57]=myLog(bi);
y[58]=x0*x1*y[1]*y[2];
y[59]=y[1]*y[54];
y[60]=x0*y[1]*y[54];
y[61]=x1*y[1]*y[54];
y[62]=x0*x1*y[1]*y[54];
y[63]=x1*y[1]*y[55];
y[64]=x0*x1*y[1]*y[55];
y[65]=y[1]*y[55]*y[56];
y[66]=x0*y[1]*y[55]*y[56];
y[67]=-(x0*x1*y[1]*y[10]);
y[68]=y[12]+y[58]+y[59]+y[60]+y[61]+y[62]+y[63]+y[64]+y[65]+y[66]+y[67];
y[69]=lrs[2];
y[70]=-(lambda*MYI*y[68]*y[69]);
y[71]=1.+y[70];
y[72]=myLog(y[71]);
y[73]=myLog(y[51]);
y[74]=y[1]*y[6]*y[27];
y[75]=y[1]*y[2]*y[27];
y[76]=-(y[1]*y[10]*y[27]);
y[77]=y[1]*y[2]*y[53];
y[78]=y[1]*y[6]*y[29];
y[79]=y[1]*y[6]*y[27]*y[29];
y[80]=y[1]*y[2]*y[27]*y[29];
y[81]=-(y[1]*y[10]*y[27]*y[29]);
y[82]=y[1]*y[2]*y[29]*y[53];
y[83]=y[7]+y[74]+y[75]+y[76]+y[77]+y[78]+y[79]+y[80]+y[81]+y[82];
y[84]=myLog(y[83]);
y[85]=myLog(x2);
y[86]=-x2;
y[87]=1.+y[86];
y[88]=2.*x2*y[1]*y[54];
y[89]=y[1]*y[55];
y[90]=x0*y[1]*y[55];
y[91]=2.*x1*y[1]*y[55];
y[92]=2.*x0*x1*y[1]*y[55];
y[93]=4.*x1*x2*y[1]*y[55];
y[94]=y[7]+y[21]+y[24]+y[59]+y[60]+y[88]+y[89]+y[90]+y[91]+y[92]+y[93];
y[95]=x2*y[1]*y[54];
y[96]=x1*x2*y[1]*y[2];
y[97]=x1*x2*y[1]*y[54];
y[98]=x1*x2*y[1]*y[55];
y[99]=x2*y[1]*y[55]*y[56];
y[100]=-(x1*x2*y[1]*y[10]);
y[101]=y[7]+y[8]+y[9]+y[11]+y[12]+y[13]+y[14]+y[15]+y[95]+y[96]+y[97]+y[98]+\
y[99]+y[100];
y[102]=y[13]+y[15]+y[59]+y[61]+y[63]+y[65];
y[103]=x2*y[1]*y[2];
y[104]=x2*y[1]*y[55];
y[105]=2.*x1*x2*y[1]*y[55];
y[106]=-(x2*y[1]*y[10]);
y[107]=y[7]+y[8]+y[9]+y[11]+y[95]+y[103]+y[104]+y[105]+y[106];
y[108]=x2*x2;
y[109]=x2*y[1]*y[6];
y[110]=x0*x2*y[1]*y[2];
y[111]=x0*x2*y[1]*y[54];
y[112]=y[1]*y[54]*y[108];
y[113]=x0*x2*y[1]*y[55];
y[114]=2.*x0*x1*x2*y[1]*y[55];
y[115]=2.*x1*y[1]*y[55]*y[108];
y[116]=-(x0*x2*y[1]*y[10]);
y[117]=y[7]+y[20]+y[21]+y[23]+y[24]+y[95]+y[104]+y[105]+y[109]+y[110]+y[111]\
+y[112]+y[113]+y[114]+y[115]+y[116];
y[118]=-(lambda*MYI*y[4]*y[5]*y[101]);
y[119]=lambda*MYI*x0*y[5]*y[101];
y[120]=1.+y[39]+y[118]+y[119];
y[121]=2.*x2*y[1]*y[55];
y[122]=2.*x0*x2*y[1]*y[55];
y[123]=2.*y[1]*y[55]*y[108];
y[124]=y[121]+y[122]+y[123];
y[125]=-(lambda*MYI*x1*y[18]*y[19]*y[124]);
y[126]=-(lambda*MYI*y[18]*y[19]*y[117]);
y[127]=lambda*MYI*x1*y[19]*y[117];
y[128]=1.+y[125]+y[126]+y[127];
y[129]=2.*x1*x2*y[1]*y[54];
y[130]=2.*x2*y[1]*y[55]*y[56];
y[131]=y[12]+y[58]+y[59]+y[60]+y[61]+y[62]+y[63]+y[64]+y[65]+y[66]+y[67]+y[1\
29]+y[130];
y[132]=-(lambda*MYI*y[69]*y[87]*y[131]);
y[133]=-(lambda*MYI*x0*y[4]*y[5]*y[101]);
y[134]=x0+y[133];
y[135]=-(lambda*MYI*x1*y[18]*y[19]*y[117]);
y[136]=x1+y[135];
y[137]=1.+y[132];
y[138]=1./y[137];
y[139]=x0*x1*y[4]*y[5]*y[18]*y[19]*y[32]*y[102]*y[107];
y[140]=-(lambda*MYI*x1*y[18]*y[19]*y[94]*y[120]);
y[141]=y[139]+y[140];
y[142]=lambda*MYI*x2*y[69]*y[87]*y[94]*y[141];
y[143]=-(x0*x1*y[4]*y[5]*y[18]*y[19]*y[32]*y[94]*y[107]);
y[144]=lambda*MYI*x0*y[4]*y[5]*y[102]*y[128];
y[145]=y[143]+y[144];
y[146]=-(lambda*MYI*x2*y[69]*y[87]*y[102]*y[145]);
y[147]=pow(y[107],2);
y[148]=x0*x1*y[4]*y[5]*y[18]*y[19]*y[32]*y[147];
y[149]=y[120]*y[128];
y[150]=y[148]+y[149];
y[151]=2.*x1*y[1]*y[54];
y[152]=2.*y[1]*y[55]*y[56];
y[153]=y[151]+y[152];
y[154]=-(lambda*MYI*x2*y[69]*y[87]*y[153]);
y[155]=lambda*MYI*x2*y[69]*y[131];
y[156]=1.+y[132]+y[154]+y[155];
y[157]=y[150]*y[156];
y[158]=y[142]+y[146]+y[157];
y[159]=y[1]*y[134];
y[160]=y[1]*y[136];
y[161]=y[1]*y[134]*y[136];
y[162]=-(lambda*MYI*x2*y[69]*y[87]*y[131]);
y[163]=x2+y[162];
y[164]=y[1]*y[136]*y[163];
y[165]=y[1]+y[159]+y[160]+y[161]+y[164];
y[166]=pow(y[165],-2);
y[167]=pow(y[134],2);
y[168]=pow(y[136],2);
y[169]=pow(y[163],2);
FOUT=-(y[30]*(y[31]*y[47]*y[52]*y[57]+y[31]*y[47]*y[52]*y[72]+3.*y[31]*y[47]\
*y[52]*y[73]-2.*y[31]*y[47]*y[52]*y[84]))+0.5*(y[47]*y[52]*(pow(y[57],2)*y[\
31]+pow(y[72],2)*y[31]+2.*y[31]*y[57]*y[72])+2.*(y[31]*y[47]*y[57]+y[31]*y[\
47]*y[72])*(3.*y[52]*y[73]-2.*y[52]*y[84])+y[31]*y[47]*(9.*pow(y[73],2)*y[5\
2]+4.*pow(y[84],2)*y[52]-12.*y[52]*y[73]*y[84]))-y[30]*y[31]*y[47]*y[52]*y[\
85]+y[30]*y[31]*y[85]*y[138]*y[158]*y[166]+y[30]*(myLog(y[137])*y[31]*y[138\
]*y[158]*y[166]+3.*myLog(y[165])*y[31]*y[138]*y[158]*y[166]-2.*myLog(y[7]+y\
[1]*y[2]*y[134]+y[1]*y[6]*y[134]-y[1]*y[10]*y[134]+y[1]*y[6]*y[136]+y[1]*y[\
2]*y[134]*y[136]+y[1]*y[6]*y[134]*y[136]-y[1]*y[10]*y[134]*y[136]+y[1]*y[54\
]*y[163]+y[1]*y[54]*y[134]*y[163]+y[1]*y[6]*y[136]*y[163]+y[1]*y[54]*y[136]\
*y[163]+y[1]*y[55]*y[136]*y[163]+y[1]*y[2]*y[134]*y[136]*y[163]-y[1]*y[10]*\
y[134]*y[136]*y[163]+y[1]*y[54]*y[134]*y[136]*y[163]+y[1]*y[55]*y[134]*y[13\
6]*y[163]+y[1]*y[2]*y[167]+y[1]*y[2]*y[136]*y[167]+y[1]*y[55]*y[163]*y[168]\
+y[1]*y[55]*y[134]*y[163]*y[168]+y[1]*y[54]*y[136]*y[169]+y[1]*y[55]*y[168]\
*y[169])*y[31]*y[138]*y[158]*y[166]+y[31]*y[57]*y[138]*y[158]*y[166]);
return (FOUT);
}
