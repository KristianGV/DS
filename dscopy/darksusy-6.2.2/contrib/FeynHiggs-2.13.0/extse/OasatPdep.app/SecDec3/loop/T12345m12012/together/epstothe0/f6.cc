#include "intfile.hh"

dcmplx Pf6(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
double x3=x[3];
dcmplx y[174];
dcmplx FOUT;
dcmplx MYI(0.,1.);
y[1]=1./bi;
y[2]=em[0];
y[3]=x0*x0;
y[4]=em[1];
y[5]=esx[0];
y[6]=2.*x0*x2*y[1]*y[4];
y[7]=-x1;
y[8]=1.+y[7];
y[9]=y[1]*y[2];
y[10]=x0*y[1]*y[2];
y[11]=y[1]*y[4];
y[12]=x0*y[1]*y[4];
y[13]=2.*x1*y[1]*y[4];
y[14]=-(y[1]*y[5]);
y[15]=-(x0*y[1]*y[5]);
y[16]=-x0;
y[17]=1.+y[16];
y[18]=x3*x3;
y[19]=x2*x2;
y[20]=lrs[0];
y[21]=x2*y[1]*y[2];
y[22]=2.*x3*y[1]*y[2];
y[23]=x2*x3*y[1]*y[2];
y[24]=y[1]*y[2]*y[18];
y[25]=x2*y[1]*y[4];
y[26]=y[1]*y[4]*y[19];
y[27]=x2*x3*y[1]*y[4];
y[28]=-(x2*y[1]*y[5]);
y[29]=-(x2*x3*y[1]*y[5]);
y[30]=x1*x3*y[1]*y[2];
y[31]=2.*x0*x2*x3*y[1]*y[2];
y[32]=x1*x2*x3*y[1]*y[2];
y[33]=2.*x0*y[1]*y[2]*y[18];
y[34]=x1*y[1]*y[2]*y[18];
y[35]=2.*x1*x2*y[1]*y[4];
y[36]=2.*x0*y[1]*y[4]*y[19];
y[37]=x1*y[1]*y[4]*y[19];
y[38]=x1*x3*y[1]*y[4];
y[39]=2.*x0*x2*x3*y[1]*y[4];
y[40]=x1*x2*x3*y[1]*y[4];
y[41]=-(x1*x3*y[1]*y[5]);
y[42]=-2.*x0*x2*x3*y[1]*y[5];
y[43]=-(x1*x2*x3*y[1]*y[5]);
y[44]=y[21]+y[22]+y[23]+y[24]+y[25]+y[26]+y[27]+y[28]+y[29]+y[30]+y[31]+y[32\
]+y[33]+y[34]+y[35]+y[36]+y[37]+y[38]+y[39]+y[40]+y[41]+y[42]+y[43];
y[45]=lrs[1];
y[46]=-x2;
y[47]=1.+y[46];
y[48]=x3*y[1]*y[2];
y[49]=2.*x0*x3*y[1]*y[2];
y[50]=2.*x2*y[1]*y[4];
y[51]=x3*y[1]*y[4];
y[52]=-(x3*y[1]*y[5]);
y[53]=lambda*lambda;
y[54]=x0*x2*y[1]*y[2];
y[55]=x0*x2*y[1]*y[4];
y[56]=-(x0*x2*y[1]*y[5]);
y[57]=y[9]+y[10]+y[11]+y[12]+y[13]+y[14]+y[15]+y[49]+y[54]+y[55]+y[56];
y[58]=y[23]+y[24]+y[26]+y[27]+y[29]+y[48]+y[50]+y[51]+y[52];
y[59]=2.*y[1]*y[2];
y[60]=x1*y[1]*y[2];
y[61]=2.*x0*x2*y[1]*y[2];
y[62]=x1*x2*y[1]*y[2];
y[63]=4.*x0*x3*y[1]*y[2];
y[64]=2.*x1*x3*y[1]*y[2];
y[65]=x1*y[1]*y[4];
y[66]=x1*x2*y[1]*y[4];
y[67]=-(x1*y[1]*y[5]);
y[68]=-2.*x0*x2*y[1]*y[5];
y[69]=-(x1*x2*y[1]*y[5]);
y[70]=y[6]+y[21]+y[22]+y[25]+y[28]+y[59]+y[60]+y[61]+y[62]+y[63]+y[64]+y[65]\
+y[66]+y[67]+y[68]+y[69];
y[71]=x0*x3*y[1]*y[2];
y[72]=x0*x3*y[1]*y[4];
y[73]=-(x0*x3*y[1]*y[5]);
y[74]=x0*x2*x3*y[1]*y[2];
y[75]=x0*y[1]*y[2]*y[18];
y[76]=x0*y[1]*y[4]*y[19];
y[77]=2.*x1*x3*y[1]*y[4];
y[78]=x0*x2*x3*y[1]*y[4];
y[79]=-(x0*x2*x3*y[1]*y[5]);
y[80]=y[6]+y[9]+y[11]+y[13]+y[14]+y[21]+y[25]+y[28]+y[35]+y[48]+y[51]+y[52]+\
y[71]+y[72]+y[73]+y[74]+y[75]+y[76]+y[77]+y[78]+y[79];
y[81]=lrs[2];
y[82]=y[1]*y[2]*y[3];
y[83]=x0*x1*y[1]*y[2];
y[84]=y[1]*y[3]*y[4];
y[85]=x0*x1*y[1]*y[4];
y[86]=-(y[1]*y[3]*y[5]);
y[87]=-(x0*x1*y[1]*y[5]);
y[88]=y[10]+y[12]+y[15]+y[82]+y[83]+y[84]+y[85]+y[86]+y[87];
y[89]=2.*x2*x3*y[1]*y[2];
y[90]=2.*y[1]*y[2]*y[18];
y[91]=2.*y[1]*y[4]*y[19];
y[92]=2.*x2*x3*y[1]*y[4];
y[93]=-2.*x2*x3*y[1]*y[5];
y[94]=y[89]+y[90]+y[91]+y[92]+y[93];
y[95]=-(lambda*MYI*x0*y[17]*y[20]*y[94]);
y[96]=-(lambda*MYI*y[17]*y[20]*y[44]);
y[97]=lambda*MYI*x0*y[20]*y[44];
y[98]=1.+y[95]+y[96]+y[97];
y[99]=2.*y[1]*y[4];
y[100]=2.*x3*y[1]*y[4];
y[101]=y[50]+y[99]+y[100];
y[102]=-(lambda*MYI*x1*y[8]*y[45]*y[101]);
y[103]=-(lambda*MYI*y[8]*y[45]*y[80]);
y[104]=lambda*MYI*x1*y[45]*y[80];
y[105]=1.+y[102]+y[103]+y[104];
y[106]=-x3;
y[107]=1.+y[106];
y[108]=4.*x0*x2*y[1]*y[4];
y[109]=2.*x0*x3*y[1]*y[4];
y[110]=-2.*x0*x3*y[1]*y[5];
y[111]=y[9]+y[11]+y[13]+y[14]+y[30]+y[35]+y[38]+y[41]+y[48]+y[49]+y[50]+y[51\
]+y[52]+y[108]+y[109]+y[110];
y[112]=2.*x0*y[1]*y[4];
y[113]=y[6]+y[9]+y[11]+y[13]+y[14]+y[71]+y[72]+y[73]+y[112];
y[114]=x0*x1*y[8]*y[17]*y[20]*y[45]*y[53]*y[58]*y[70];
y[115]=-(lambda*MYI*x1*y[8]*y[45]*y[57]*y[98]);
y[116]=y[114]+y[115];
y[117]=2.*x0*x1*y[1]*y[4];
y[118]=x3*y[1]*y[2]*y[3];
y[119]=x0*x1*x3*y[1]*y[2];
y[120]=x1*x1;
y[121]=y[1]*y[4]*y[120];
y[122]=2.*x2*y[1]*y[3]*y[4];
y[123]=2.*x0*x1*x2*y[1]*y[4];
y[124]=x3*y[1]*y[3]*y[4];
y[125]=x0*x1*x3*y[1]*y[4];
y[126]=-(x3*y[1]*y[3]*y[5]);
y[127]=-(x0*x1*x3*y[1]*y[5]);
y[128]=y[6]+y[9]+y[10]+y[12]+y[15]+y[60]+y[65]+y[67]+y[71]+y[72]+y[73]+y[117\
]+y[118]+y[119]+y[121]+y[122]+y[123]+y[124]+y[125]+y[126]+y[127];
y[129]=lrs[3];
y[130]=x0*x1*y[8]*y[17]*y[20]*y[45]*y[53]*y[70]*y[113];
y[131]=-(x0*x1*y[8]*y[17]*y[20]*y[45]*y[53]*y[57]*y[111]);
y[132]=y[130]+y[131];
y[133]=-(x0*x1*y[8]*y[17]*y[20]*y[45]*y[53]*y[57]*y[58]);
y[134]=lambda*MYI*x0*y[17]*y[20]*y[70]*y[105];
y[135]=y[133]+y[134];
y[136]=2.*y[1]*y[3]*y[4];
y[137]=y[112]+y[117]+y[136];
y[138]=-(lambda*MYI*x2*y[47]*y[81]*y[137]);
y[139]=-(lambda*MYI*y[47]*y[81]*y[128]);
y[140]=lambda*MYI*x2*y[81]*y[128];
y[141]=1.+y[138]+y[139]+y[140];
y[142]=x0*x1*y[8]*y[17]*y[20]*y[45]*y[53]*y[58]*y[111];
y[143]=-(lambda*MYI*x1*y[8]*y[45]*y[98]*y[113]);
y[144]=y[142]+y[143];
y[145]=-(x0*x1*y[8]*y[17]*y[20]*y[45]*y[53]*y[58]*y[113]);
y[146]=lambda*MYI*x0*y[17]*y[20]*y[105]*y[111];
y[147]=y[145]+y[146];
y[148]=pow(y[58],2);
y[149]=x0*x1*y[8]*y[17]*y[20]*y[45]*y[53]*y[148];
y[150]=y[98]*y[105];
y[151]=y[149]+y[150];
y[152]=2.*x0*y[1]*y[2];
y[153]=x2*y[1]*y[2]*y[3];
y[154]=x0*x1*x2*y[1]*y[2];
y[155]=2.*x3*y[1]*y[2]*y[3];
y[156]=2.*x0*x1*x3*y[1]*y[2];
y[157]=x2*y[1]*y[3]*y[4];
y[158]=x0*x1*x2*y[1]*y[4];
y[159]=-(x2*y[1]*y[3]*y[5]);
y[160]=-(x0*x1*x2*y[1]*y[5]);
y[161]=y[9]+y[49]+y[54]+y[55]+y[56]+y[60]+y[65]+y[67]+y[83]+y[85]+y[87]+y[12\
1]+y[152]+y[153]+y[154]+y[155]+y[156]+y[157]+y[158]+y[159]+y[160];
y[162]=-(lambda*MYI*x2*y[47]*y[81]*y[128]);
y[163]=x2+y[162];
y[164]=-(lambda*MYI*x1*y[8]*y[45]*y[80]);
y[165]=x1+y[164];
y[166]=-(lambda*MYI*x0*y[17]*y[20]*y[44]);
y[167]=x0+y[166];
y[168]=-(lambda*MYI*x3*y[107]*y[129]*y[161]);
y[169]=x3+y[168];
y[170]=pow(y[165],2);
y[171]=pow(y[163],2);
y[172]=pow(y[167],2);
y[173]=pow(y[169],2);
FOUT=(pow(bi,-2)*(-(lambda*MYI*x3*y[57]*y[107]*y[129]*(-(lambda*MYI*x2*y[47]\
*y[81]*y[111]*y[132])-y[116]*y[141]-lambda*MYI*x2*y[47]*y[81]*y[88]*y[144])\
)+lambda*MYI*x3*y[70]*y[107]*y[129]*(-(lambda*MYI*x2*y[47]*y[81]*y[113]*y[1\
32])-y[135]*y[141]-lambda*MYI*x2*y[47]*y[81]*y[88]*y[147])+lambda*MYI*x3*y[\
88]*y[107]*y[129]*(lambda*MYI*x2*y[47]*y[81]*y[113]*y[116]-lambda*MYI*x2*y[\
47]*y[81]*y[111]*y[135]-lambda*MYI*x2*y[47]*y[81]*y[88]*y[151])+(lambda*MYI\
*x2*y[47]*y[81]*y[113]*y[144]-lambda*MYI*x2*y[47]*y[81]*y[111]*y[147]+y[141\
]*y[151])*(1.-lambda*MYI*x3*y[107]*y[129]*(2.*x0*x1*y[1]*y[2]+2.*y[1]*y[2]*\
y[3]+y[152])+lambda*MYI*x3*y[129]*y[161]-lambda*MYI*y[107]*y[129]*y[161])))\
/((y[1]+y[1]*y[163]+y[1]*y[165]+y[1]*y[163]*y[165]+y[1]*y[163]*y[167]+y[1]*\
y[169]+y[1]*y[165]*y[169]+y[1]*y[167]*y[169])*(y[9]+y[1]*y[2]*y[163]+y[1]*y\
[2]*y[165]+y[1]*y[4]*y[165]-y[1]*y[5]*y[165]+y[1]*y[2]*y[163]*y[165]+y[1]*y\
[4]*y[163]*y[165]-y[1]*y[5]*y[163]*y[165]+y[1]*y[2]*y[163]*y[167]+y[1]*y[4]\
*y[163]*y[167]-y[1]*y[5]*y[163]*y[167]+2.*y[1]*y[4]*y[163]*y[165]*y[167]+y[\
1]*y[2]*y[169]+y[1]*y[2]*y[165]*y[169]+y[1]*y[4]*y[165]*y[169]-y[1]*y[5]*y[\
165]*y[169]+2.*y[1]*y[2]*y[167]*y[169]+y[1]*y[2]*y[163]*y[167]*y[169]+y[1]*\
y[4]*y[163]*y[167]*y[169]-y[1]*y[5]*y[163]*y[167]*y[169]+y[1]*y[2]*y[165]*y\
[167]*y[169]+y[1]*y[4]*y[165]*y[167]*y[169]-y[1]*y[5]*y[165]*y[167]*y[169]+\
y[1]*y[2]*y[163]*y[165]*y[167]*y[169]+y[1]*y[4]*y[163]*y[165]*y[167]*y[169]\
-y[1]*y[5]*y[163]*y[165]*y[167]*y[169]+y[1]*y[4]*y[170]+y[1]*y[4]*y[163]*y[\
170]+y[1]*y[4]*y[169]*y[170]+y[1]*y[4]*y[167]*y[171]+y[1]*y[4]*y[165]*y[167\
]*y[171]+y[1]*y[2]*y[163]*y[169]*y[172]+y[1]*y[4]*y[163]*y[169]*y[172]-y[1]\
*y[5]*y[163]*y[169]*y[172]+y[1]*y[4]*y[171]*y[172]+y[1]*y[2]*y[167]*y[173]+\
y[1]*y[2]*y[165]*y[167]*y[173]+y[1]*y[2]*y[172]*y[173]));
return (FOUT);
}
