#include "intfile.hh"

dcmplx Pf3(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
dcmplx y[196];
dcmplx FOUT;
dcmplx MYI(0.,1.);
y[1]=1./bi;
y[2]=em[2];
y[3]=-x0;
y[4]=1.+y[3];
y[5]=lrs[0];
y[6]=-(lambda*MYI*x0*y[1]*y[2]*y[4]*y[5]);
y[7]=x0+y[6];
y[8]=em[0];
y[9]=em[1];
y[10]=lrs[1];
y[11]=x0*x0;
y[12]=esx[0];
y[13]=pow(bi,-2);
y[14]=lambda*lambda;
y[15]=pow(y[2],2);
y[16]=1./x1;
y[17]=myLog(x1);
y[18]=y[1]*y[8];
y[19]=x0*y[1]*y[8];
y[20]=x0*y[1]*y[9];
y[21]=y[1]*y[9]*y[11];
y[22]=x0*y[1]*y[2];
y[23]=-(x0*y[1]*y[12]);
y[24]=x1*x1;
y[25]=-(lambda*MYI*y[1]*y[2]*y[5]);
y[26]=2.*lambda*MYI*x0*y[1]*y[2]*y[5];
y[27]=-(lambda*MYI*y[1]*y[8]*y[10]);
y[28]=-(lambda*MYI*x0*y[1]*y[8]*y[10]);
y[29]=-(lambda*MYI*x0*y[1]*y[9]*y[10]);
y[30]=-(lambda*MYI*y[1]*y[9]*y[10]*y[11]);
y[31]=-(lambda*MYI*x0*y[1]*y[2]*y[10]);
y[32]=lambda*MYI*x0*y[1]*y[10]*y[12];
y[33]=pow(y[8],2);
y[34]=pow(x1,3);
y[35]=pow(x1,4);
y[36]=pow(x0,3);
y[37]=pow(y[9],2);
y[38]=pow(x0,4);
y[39]=-(y[2]*y[5]*y[8]*y[10]*y[13]*y[14]);
y[40]=x0*y[2]*y[5]*y[8]*y[10]*y[13]*y[14];
y[41]=2.*y[2]*y[5]*y[8]*y[10]*y[11]*y[13]*y[14];
y[42]=-(x0*y[2]*y[5]*y[9]*y[10]*y[13]*y[14]);
y[43]=y[2]*y[5]*y[9]*y[10]*y[11]*y[13]*y[14];
y[44]=2.*y[2]*y[5]*y[9]*y[10]*y[13]*y[14]*y[36];
y[45]=-(x0*y[5]*y[10]*y[13]*y[14]*y[15]);
y[46]=2.*y[5]*y[10]*y[11]*y[13]*y[14]*y[15];
y[47]=x0*y[2]*y[5]*y[10]*y[12]*y[13]*y[14];
y[48]=-2.*y[2]*y[5]*y[10]*y[11]*y[12]*y[13]*y[14];
y[49]=pow(y[12],2);
y[50]=y[1]*y[2];
y[51]=x1*y[1]*y[8];
y[52]=y[1]*y[8]*y[24];
y[53]=x1*y[1]*y[9];
y[54]=2.*x0*x1*y[1]*y[9];
y[55]=2.*x0*y[1]*y[9]*y[24];
y[56]=x1*y[1]*y[2];
y[57]=-(x1*y[1]*y[12]);
y[58]=y[50]+y[51]+y[52]+y[53]+y[54]+y[55]+y[56]+y[57];
y[59]=-(lambda*MYI*x0*y[4]*y[5]*y[58]);
y[60]=x0+y[59];
y[61]=-x1;
y[62]=1.+y[61];
y[63]=2.*x0*x1*y[1]*y[8];
y[64]=2.*x1*y[1]*y[9]*y[11];
y[65]=y[18]+y[19]+y[20]+y[21]+y[22]+y[23]+y[63]+y[64];
y[66]=-(lambda*MYI*x1*y[10]*y[62]*y[65]);
y[67]=x1+y[66];
y[68]=pow(y[60],2);
y[69]=pow(y[67],2);
y[70]=y[1]*y[7];
y[71]=y[1]+y[70];
y[72]=pow(y[71],-3);
y[73]=y[1]*y[2]*y[7];
y[74]=y[50]+y[73];
y[75]=y[18]+y[19]+y[20]+y[21]+y[22]+y[23];
y[76]=-(lambda*MYI*y[10]*y[75]);
y[77]=1.+y[76];
y[78]=1./y[77];
y[79]=1.+y[25]+y[26]+y[27]+y[28]+y[29]+y[30]+y[31]+y[32]+y[39]+y[40]+y[41]+y\
[42]+y[43]+y[44]+y[45]+y[46]+y[47]+y[48];
y[80]=myLog(bi);
y[81]=-(lambda*MYI*y[10]*y[62]*y[65]);
y[82]=1.+y[81];
y[83]=1./y[82];
y[84]=-(lambda*MYI*x1*y[1]*y[5]*y[8]);
y[85]=2.*lambda*MYI*x0*x1*y[1]*y[5]*y[8];
y[86]=-(lambda*MYI*y[1]*y[5]*y[8]*y[24]);
y[87]=2.*lambda*MYI*x0*y[1]*y[5]*y[8]*y[24];
y[88]=-(lambda*MYI*x1*y[1]*y[5]*y[9]);
y[89]=-2.*lambda*MYI*x0*x1*y[1]*y[5]*y[9];
y[90]=6.*lambda*MYI*x1*y[1]*y[5]*y[9]*y[11];
y[91]=-4.*lambda*MYI*x0*y[1]*y[5]*y[9]*y[24];
y[92]=6.*lambda*MYI*y[1]*y[5]*y[9]*y[11]*y[24];
y[93]=-(lambda*MYI*x1*y[1]*y[2]*y[5]);
y[94]=2.*lambda*MYI*x0*x1*y[1]*y[2]*y[5];
y[95]=lambda*MYI*x1*y[1]*y[5]*y[12];
y[96]=-2.*lambda*MYI*x0*x1*y[1]*y[5]*y[12];
y[97]=2.*lambda*MYI*x1*y[1]*y[8]*y[10];
y[98]=-2.*lambda*MYI*x0*x1*y[1]*y[8]*y[10];
y[99]=6.*lambda*MYI*x0*y[1]*y[8]*y[10]*y[24];
y[100]=2.*lambda*MYI*x0*x1*y[1]*y[9]*y[10];
y[101]=-2.*lambda*MYI*x1*y[1]*y[9]*y[10]*y[11];
y[102]=6.*lambda*MYI*y[1]*y[9]*y[10]*y[11]*y[24];
y[103]=2.*lambda*MYI*x0*x1*y[1]*y[2]*y[10];
y[104]=-2.*lambda*MYI*x0*x1*y[1]*y[10]*y[12];
y[105]=-(x1*y[5]*y[10]*y[13]*y[14]*y[33]);
y[106]=2.*x0*x1*y[5]*y[10]*y[13]*y[14]*y[33];
y[107]=x1*y[5]*y[10]*y[11]*y[13]*y[14]*y[33];
y[108]=y[5]*y[10]*y[13]*y[14]*y[24]*y[33];
y[109]=-2.*x0*y[5]*y[10]*y[13]*y[14]*y[24]*y[33];
y[110]=3.*y[5]*y[10]*y[11]*y[13]*y[14]*y[24]*y[33];
y[111]=2.*y[5]*y[10]*y[13]*y[14]*y[33]*y[34];
y[112]=-8.*y[5]*y[10]*y[11]*y[13]*y[14]*y[33]*y[34];
y[113]=2.*x0*y[5]*y[10]*y[13]*y[14]*y[33]*y[35];
y[114]=-8.*y[5]*y[10]*y[11]*y[13]*y[14]*y[33]*y[35];
y[115]=-(x1*y[5]*y[8]*y[9]*y[10]*y[13]*y[14]);
y[116]=-2.*x0*x1*y[5]*y[8]*y[9]*y[10]*y[13]*y[14];
y[117]=7.*x1*y[5]*y[8]*y[9]*y[10]*y[11]*y[13]*y[14];
y[118]=4.*x1*y[5]*y[8]*y[9]*y[10]*y[13]*y[14]*y[36];
y[119]=2.*y[5]*y[8]*y[9]*y[10]*y[13]*y[14]*y[24];
y[120]=x0*y[5]*y[8]*y[9]*y[10]*y[13]*y[14]*y[24];
y[121]=-9.*y[5]*y[8]*y[9]*y[10]*y[11]*y[13]*y[14]*y[24];
y[122]=12.*y[5]*y[8]*y[9]*y[10]*y[13]*y[14]*y[24]*y[36];
y[123]=12.*x0*y[5]*y[8]*y[9]*y[10]*y[13]*y[14]*y[34];
y[124]=-4.*y[5]*y[8]*y[9]*y[10]*y[11]*y[13]*y[14]*y[34];
y[125]=-32.*y[5]*y[8]*y[9]*y[10]*y[13]*y[14]*y[34]*y[36];
y[126]=14.*y[5]*y[8]*y[9]*y[10]*y[11]*y[13]*y[14]*y[35];
y[127]=-32.*y[5]*y[8]*y[9]*y[10]*y[13]*y[14]*y[35]*y[36];
y[128]=4.*x1*y[5]*y[10]*y[13]*y[14]*y[36]*y[37];
y[129]=2.*x1*y[5]*y[10]*y[13]*y[14]*y[37]*y[38];
y[130]=x0*y[5]*y[10]*y[13]*y[14]*y[24]*y[37];
y[131]=3.*y[5]*y[10]*y[11]*y[13]*y[14]*y[24]*y[37];
y[132]=-6.*y[5]*y[10]*y[13]*y[14]*y[24]*y[36]*y[37];
y[133]=6.*y[5]*y[10]*y[13]*y[14]*y[24]*y[37]*y[38];
y[134]=6.*y[5]*y[10]*y[11]*y[13]*y[14]*y[34]*y[37];
y[135]=-24.*y[5]*y[10]*y[13]*y[14]*y[34]*y[37]*y[38];
y[136]=8.*y[5]*y[10]*y[13]*y[14]*y[35]*y[36]*y[37];
y[137]=-20.*y[5]*y[10]*y[13]*y[14]*y[35]*y[37]*y[38];
y[138]=x1*y[2]*y[5]*y[8]*y[10]*y[13]*y[14];
y[139]=-4.*x0*x1*y[2]*y[5]*y[8]*y[10]*y[13]*y[14];
y[140]=6.*x1*y[2]*y[5]*y[8]*y[10]*y[11]*y[13]*y[14];
y[141]=2.*y[2]*y[5]*y[8]*y[10]*y[13]*y[14]*y[24];
y[142]=3.*x0*y[2]*y[5]*y[8]*y[10]*y[13]*y[14]*y[24];
y[143]=-12.*y[2]*y[5]*y[8]*y[10]*y[11]*y[13]*y[14]*y[24];
y[144]=4.*x0*y[2]*y[5]*y[8]*y[10]*y[13]*y[14]*y[34];
y[145]=-12.*y[2]*y[5]*y[8]*y[10]*y[11]*y[13]*y[14]*y[34];
y[146]=2.*x0*x1*y[2]*y[5]*y[9]*y[10]*y[13]*y[14];
y[147]=-5.*x1*y[2]*y[5]*y[9]*y[10]*y[11]*y[13]*y[14];
y[148]=8.*x1*y[2]*y[5]*y[9]*y[10]*y[13]*y[14]*y[36];
y[149]=2.*x0*y[2]*y[5]*y[9]*y[10]*y[13]*y[14]*y[24];
y[150]=6.*y[2]*y[5]*y[9]*y[10]*y[11]*y[13]*y[14]*y[24];
y[151]=-18.*y[2]*y[5]*y[9]*y[10]*y[13]*y[14]*y[24]*y[36];
y[152]=6.*y[2]*y[5]*y[9]*y[10]*y[11]*y[13]*y[14]*y[34];
y[153]=-16.*y[2]*y[5]*y[9]*y[10]*y[13]*y[14]*y[34]*y[36];
y[154]=2.*x0*x1*y[5]*y[10]*y[13]*y[14]*y[15];
y[155]=-3.*x1*y[5]*y[10]*y[11]*y[13]*y[14]*y[15];
y[156]=x0*y[5]*y[10]*y[13]*y[14]*y[15]*y[24];
y[157]=-3.*y[5]*y[10]*y[11]*y[13]*y[14]*y[15]*y[24];
y[158]=x1*y[5]*y[8]*y[10]*y[12]*y[13]*y[14];
y[159]=-2.*x0*x1*y[5]*y[8]*y[10]*y[12]*y[13]*y[14];
y[160]=-2.*x1*y[5]*y[8]*y[10]*y[11]*y[12]*y[13]*y[14];
y[161]=-2.*y[5]*y[8]*y[10]*y[12]*y[13]*y[14]*y[24];
y[162]=3.*x0*y[5]*y[8]*y[10]*y[12]*y[13]*y[14]*y[24];
y[163]=-4.*x0*y[5]*y[8]*y[10]*y[12]*y[13]*y[14]*y[34];
y[164]=12.*y[5]*y[8]*y[10]*y[11]*y[12]*y[13]*y[14]*y[34];
y[165]=-(x1*y[5]*y[9]*y[10]*y[11]*y[12]*y[13]*y[14]);
y[166]=-4.*x1*y[5]*y[9]*y[10]*y[12]*y[13]*y[14]*y[36];
y[167]=-2.*x0*y[5]*y[9]*y[10]*y[12]*y[13]*y[14]*y[24];
y[168]=6.*y[5]*y[9]*y[10]*y[12]*y[13]*y[14]*y[24]*y[36];
y[169]=-6.*y[5]*y[9]*y[10]*y[11]*y[12]*y[13]*y[14]*y[34];
y[170]=16.*y[5]*y[9]*y[10]*y[12]*y[13]*y[14]*y[34]*y[36];
y[171]=-2.*x0*x1*y[2]*y[5]*y[10]*y[12]*y[13]*y[14];
y[172]=2.*x1*y[2]*y[5]*y[10]*y[11]*y[12]*y[13]*y[14];
y[173]=-2.*x0*y[2]*y[5]*y[10]*y[12]*y[13]*y[14]*y[24];
y[174]=6.*y[2]*y[5]*y[10]*y[11]*y[12]*y[13]*y[14]*y[24];
y[175]=x1*y[5]*y[10]*y[11]*y[13]*y[14]*y[49];
y[176]=x0*y[5]*y[10]*y[13]*y[14]*y[24]*y[49];
y[177]=-3.*y[5]*y[10]*y[11]*y[13]*y[14]*y[24]*y[49];
y[178]=1.+y[25]+y[26]+y[27]+y[28]+y[29]+y[30]+y[31]+y[32]+y[39]+y[40]+y[41]+\
y[42]+y[43]+y[44]+y[45]+y[46]+y[47]+y[48]+y[84]+y[85]+y[86]+y[87]+y[88]+y[8\
9]+y[90]+y[91]+y[92]+y[93]+y[94]+y[95]+y[96]+y[97]+y[98]+y[99]+y[100]+y[101\
]+y[102]+y[103]+y[104]+y[105]+y[106]+y[107]+y[108]+y[109]+y[110]+y[111]+y[1\
12]+y[113]+y[114]+y[115]+y[116]+y[117]+y[118]+y[119]+y[120]+y[121]+y[122]+y\
[123]+y[124]+y[125]+y[126]+y[127]+y[128]+y[129]+y[130]+y[131]+y[132]+y[133]\
+y[134]+y[135]+y[136]+y[137]+y[138]+y[139]+y[140]+y[141]+y[142]+y[143]+y[14\
4]+y[145]+y[146]+y[147]+y[148]+y[149]+y[150]+y[151]+y[152]+y[153]+y[154]+y[\
155]+y[156]+y[157]+y[158]+y[159]+y[160]+y[161]+y[162]+y[163]+y[164]+y[165]+\
y[166]+y[167]+y[168]+y[169]+y[170]+y[171]+y[172]+y[173]+y[174]+y[175]+y[176\
]+y[177];
y[179]=y[1]*y[60];
y[180]=y[1]*y[60]*y[67];
y[181]=y[1]+y[179]+y[180];
y[182]=pow(y[181],-3);
y[183]=y[1]*y[2]*y[60];
y[184]=y[1]*y[8]*y[67];
y[185]=y[1]*y[8]*y[60]*y[67];
y[186]=y[1]*y[9]*y[60]*y[67];
y[187]=y[1]*y[2]*y[60]*y[67];
y[188]=-(y[1]*y[12]*y[60]*y[67]);
y[189]=y[1]*y[9]*y[67]*y[68];
y[190]=y[1]*y[8]*y[60]*y[69];
y[191]=y[1]*y[9]*y[68]*y[69];
y[192]=y[50]+y[183]+y[184]+y[185]+y[186]+y[187]+y[188]+y[189]+y[190]+y[191];
y[193]=myLog(y[71]);
y[194]=myLog(y[74]);
y[195]=myLog(y[77]);
FOUT=-(y[13]*y[16]*y[17]*y[72]*y[74]*y[78]*y[79])+y[13]*y[16]*y[17]*y[83]*y[\
178]*y[182]*y[192]+y[16]*(myLog(y[82])*y[13]*y[83]*y[178]*y[182]*y[192]+3.*\
myLog(y[181])*y[13]*y[83]*y[178]*y[182]*y[192]-2.*myLog(y[192])*y[13]*y[83]\
*y[178]*y[182]*y[192]+y[13]*y[80]*y[83]*y[178]*y[182]*y[192])-y[16]*(y[13]*\
y[72]*y[74]*y[78]*y[79]*y[80]+3.*y[13]*y[72]*y[74]*y[78]*y[79]*y[193]-2.*y[\
13]*y[72]*y[74]*y[78]*y[79]*y[194]+y[13]*y[72]*y[74]*y[78]*y[79]*y[195])+0.\
5*(y[13]*y[78]*y[79]*(9.*pow(y[193],2)*y[72]*y[74]+4.*pow(y[194],2)*y[72]*y\
[74]-12.*y[72]*y[74]*y[193]*y[194])+2.*(3.*y[72]*y[74]*y[193]-2.*y[72]*y[74\
]*y[194])*(y[13]*y[78]*y[79]*y[80]+y[13]*y[78]*y[79]*y[195])+y[72]*y[74]*y[\
79]*(pow(y[80],2)*y[13]*y[78]+pow(y[195],2)*y[13]*y[78]+2.*y[13]*y[78]*y[80\
]*y[195]));
return (FOUT);
}
