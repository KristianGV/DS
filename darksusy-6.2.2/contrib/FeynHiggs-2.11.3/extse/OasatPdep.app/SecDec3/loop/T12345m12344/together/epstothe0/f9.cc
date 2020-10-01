#include "intfile.hh"

dcmplx Pf9(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
double x3=x[3];
dcmplx y[209];
dcmplx FOUT;
dcmplx MYI(0.,1.);
y[1]=1./bi;
y[2]=em[0];
y[3]=x0*x0;
y[4]=em[2];
y[5]=em[3];
y[6]=esx[0];
y[7]=em[1];
y[8]=x0*y[1]*y[4];
y[9]=x0*y[1]*y[5];
y[10]=-(x0*y[1]*y[6]);
y[11]=-(y[1]*y[6]);
y[12]=y[1]*y[7];
y[13]=x2*y[1]*y[4];
y[14]=x2*x2;
y[15]=x2*y[1]*y[5];
y[16]=-(x2*y[1]*y[6]);
y[17]=-x1;
y[18]=1.+y[17];
y[19]=x0*y[1]*y[2];
y[20]=x0*y[1]*y[7];
y[21]=2.*y[1]*y[5];
y[22]=2.*x1*y[1]*y[5];
y[23]=-x0;
y[24]=1.+y[23];
y[25]=x3*x3;
y[26]=lrs[0];
y[27]=x3*y[1]*y[2];
y[28]=x2*x3*y[1]*y[2];
y[29]=y[1]*y[2]*y[25];
y[30]=x1*y[1]*y[7];
y[31]=x2*y[1]*y[7];
y[32]=2.*x0*x2*y[1]*y[7];
y[33]=x3*y[1]*y[7];
y[34]=x1*x2*y[1]*y[4];
y[35]=y[1]*y[4]*y[14];
y[36]=2.*x0*y[1]*y[4]*y[14];
y[37]=x2*x3*y[1]*y[4];
y[38]=x1*x2*y[1]*y[5];
y[39]=x2*x3*y[1]*y[5];
y[40]=-(x3*y[1]*y[6]);
y[41]=x1*x3*y[1]*y[2];
y[42]=2.*x0*x2*x3*y[1]*y[2];
y[43]=x1*x2*x3*y[1]*y[2];
y[44]=x1*y[1]*y[2]*y[25];
y[45]=2.*x0*x2*y[1]*y[2]*y[25];
y[46]=x1*x2*y[1]*y[7];
y[47]=x1*x3*y[1]*y[7];
y[48]=2.*x0*x2*x3*y[1]*y[7];
y[49]=x1*y[1]*y[4]*y[14];
y[50]=x1*x2*x3*y[1]*y[4];
y[51]=2.*x0*x3*y[1]*y[4]*y[14];
y[52]=x1*x2*x3*y[1]*y[5];
y[53]=-(x1*x2*y[1]*y[6]);
y[54]=-(x1*x3*y[1]*y[6]);
y[55]=-(x2*x3*y[1]*y[6]);
y[56]=-2.*x0*x2*x3*y[1]*y[6];
y[57]=y[12]+y[13]+y[15]+y[27]+y[28]+y[29]+y[30]+y[31]+y[32]+y[33]+y[34]+y[35\
]+y[36]+y[37]+y[38]+y[39]+y[40]+y[41]+y[42]+y[43]+y[44]+y[45]+y[46]+y[47]+y\
[48]+y[49]+y[50]+y[51]+y[52]+y[53]+y[54]+y[55]+y[56];
y[58]=lrs[1];
y[59]=-x2;
y[60]=1.+y[59];
y[61]=2.*x0*x3*y[1]*y[2];
y[62]=-(x1*y[1]*y[6]);
y[63]=lambda*lambda;
y[64]=x0*x2*y[1]*y[2];
y[65]=x0*x2*y[1]*y[4];
y[66]=x0*x2*y[1]*y[5];
y[67]=y[10]+y[11]+y[19]+y[20]+y[21]+y[22]+y[61]+y[64]+y[65]+y[66];
y[68]=y[12]+y[13]+y[15]+y[16]+y[27]+y[28]+y[29]+y[31]+y[33]+y[35]+y[37]+y[39\
]+y[40];
y[69]=y[1]*y[2];
y[70]=x1*y[1]*y[2];
y[71]=x2*y[1]*y[2];
y[72]=2.*x0*x2*y[1]*y[2];
y[73]=x1*x2*y[1]*y[2];
y[74]=2.*x3*y[1]*y[2];
y[75]=2.*x1*x3*y[1]*y[2];
y[76]=4.*x0*x2*x3*y[1]*y[2];
y[77]=-2.*x0*x2*y[1]*y[6];
y[78]=y[11]+y[12]+y[13]+y[15]+y[16]+y[30]+y[32]+y[34]+y[36]+y[38]+y[62]+y[69\
]+y[70]+y[71]+y[72]+y[73]+y[74]+y[75]+y[76]+y[77];
y[79]=x0*x3*y[1]*y[2];
y[80]=2.*x2*y[1]*y[5];
y[81]=2.*x3*y[1]*y[5];
y[82]=x0*x2*x3*y[1]*y[2];
y[83]=x0*y[1]*y[2]*y[25];
y[84]=x0*x2*y[1]*y[7];
y[85]=x0*x3*y[1]*y[7];
y[86]=x0*y[1]*y[4]*y[14];
y[87]=x0*x2*x3*y[1]*y[4];
y[88]=2.*x1*x2*y[1]*y[5];
y[89]=2.*x1*x3*y[1]*y[5];
y[90]=x0*x2*x3*y[1]*y[5];
y[91]=-(x0*x2*y[1]*y[6]);
y[92]=-(x0*x3*y[1]*y[6]);
y[93]=y[11]+y[16]+y[20]+y[21]+y[22]+y[40]+y[65]+y[66]+y[79]+y[80]+y[81]+y[82\
]+y[83]+y[84]+y[85]+y[86]+y[87]+y[88]+y[89]+y[90]+y[91]+y[92];
y[94]=lrs[2];
y[95]=y[1]*y[2]*y[3];
y[96]=x0*x1*y[1]*y[2];
y[97]=2.*x3*y[1]*y[2]*y[3];
y[98]=y[1]*y[3]*y[7];
y[99]=x0*x1*y[1]*y[4];
y[100]=2.*x2*y[1]*y[3]*y[4];
y[101]=x0*x1*y[1]*y[5];
y[102]=-(y[1]*y[3]*y[6]);
y[103]=y[8]+y[9]+y[10]+y[19]+y[95]+y[96]+y[97]+y[98]+y[99]+y[100]+y[101]+y[1\
02];
y[104]=2.*x2*x3*y[1]*y[2];
y[105]=2.*x2*y[1]*y[2]*y[25];
y[106]=2.*x2*y[1]*y[7];
y[107]=2.*x2*x3*y[1]*y[7];
y[108]=2.*y[1]*y[4]*y[14];
y[109]=2.*x3*y[1]*y[4]*y[14];
y[110]=-2.*x2*x3*y[1]*y[6];
y[111]=y[104]+y[105]+y[106]+y[107]+y[108]+y[109]+y[110];
y[112]=-(lambda*MYI*x0*y[24]*y[26]*y[111]);
y[113]=-(lambda*MYI*y[24]*y[26]*y[57]);
y[114]=lambda*MYI*x0*y[26]*y[57];
y[115]=1.+y[112]+y[113]+y[114];
y[116]=y[21]+y[80]+y[81];
y[117]=-(lambda*MYI*x1*y[18]*y[58]*y[116]);
y[118]=-(lambda*MYI*y[18]*y[58]*y[93]);
y[119]=lambda*MYI*x1*y[58]*y[93];
y[120]=1.+y[117]+y[118]+y[119];
y[121]=-x3;
y[122]=1.+y[121];
y[123]=2.*x0*y[1]*y[2]*y[25];
y[124]=2.*x0*y[1]*y[7];
y[125]=2.*x0*x3*y[1]*y[7];
y[126]=y[1]*y[4];
y[127]=x1*y[1]*y[4];
y[128]=2.*x2*y[1]*y[4];
y[129]=4.*x0*x2*y[1]*y[4];
y[130]=2.*x1*x2*y[1]*y[4];
y[131]=x3*y[1]*y[4];
y[132]=x1*x3*y[1]*y[4];
y[133]=4.*x0*x2*x3*y[1]*y[4];
y[134]=y[1]*y[5];
y[135]=x1*y[1]*y[5];
y[136]=x3*y[1]*y[5];
y[137]=x1*x3*y[1]*y[5];
y[138]=-2.*x0*x3*y[1]*y[6];
y[139]=y[12]+y[27]+y[30]+y[40]+y[41]+y[61]+y[62]+y[123]+y[124]+y[125]+y[126]\
+y[127]+y[128]+y[129]+y[130]+y[131]+y[132]+y[133]+y[134]+y[135]+y[136]+y[13\
7]+y[138];
y[140]=2.*x0*x2*y[1]*y[4];
y[141]=x0*x3*y[1]*y[4];
y[142]=x0*x3*y[1]*y[5];
y[143]=y[8]+y[9]+y[10]+y[11]+y[20]+y[21]+y[22]+y[79]+y[140]+y[141]+y[142];
y[144]=x0*x1*y[18]*y[24]*y[26]*y[58]*y[63]*y[68]*y[78];
y[145]=-(lambda*MYI*x1*y[18]*y[58]*y[67]*y[115]);
y[146]=y[144]+y[145];
y[147]=x3*y[1]*y[2]*y[3];
y[148]=x0*x1*x3*y[1]*y[2];
y[149]=y[1]*y[2]*y[3]*y[25];
y[150]=x0*x1*y[1]*y[7];
y[151]=x3*y[1]*y[3]*y[7];
y[152]=2.*x0*x1*x2*y[1]*y[4];
y[153]=x0*x1*x3*y[1]*y[4];
y[154]=2.*x2*x3*y[1]*y[3]*y[4];
y[155]=x1*x1;
y[156]=y[1]*y[5]*y[155];
y[157]=x0*x1*x3*y[1]*y[5];
y[158]=-(x0*x1*y[1]*y[6]);
y[159]=-(x3*y[1]*y[3]*y[6]);
y[160]=y[8]+y[9]+y[20]+y[22]+y[62]+y[79]+y[92]+y[98]+y[99]+y[100]+y[101]+y[1\
34]+y[140]+y[141]+y[142]+y[147]+y[148]+y[149]+y[150]+y[151]+y[152]+y[153]+y\
[154]+y[156]+y[157]+y[158]+y[159];
y[161]=lrs[3];
y[162]=x0*x1*y[18]*y[24]*y[26]*y[58]*y[63]*y[78]*y[143];
y[163]=-(x0*x1*y[18]*y[24]*y[26]*y[58]*y[63]*y[67]*y[139]);
y[164]=y[162]+y[163];
y[165]=-(x0*x1*y[18]*y[24]*y[26]*y[58]*y[63]*y[67]*y[68]);
y[166]=lambda*MYI*x0*y[24]*y[26]*y[78]*y[120];
y[167]=y[165]+y[166];
y[168]=2.*x0*y[1]*y[4];
y[169]=2.*y[1]*y[3]*y[4];
y[170]=2.*x0*x1*y[1]*y[4];
y[171]=2.*x3*y[1]*y[3]*y[4];
y[172]=y[168]+y[169]+y[170]+y[171];
y[173]=-(lambda*MYI*x2*y[60]*y[94]*y[172]);
y[174]=-(lambda*MYI*y[60]*y[94]*y[160]);
y[175]=lambda*MYI*x2*y[94]*y[160];
y[176]=1.+y[173]+y[174]+y[175];
y[177]=x0*x1*y[18]*y[24]*y[26]*y[58]*y[63]*y[68]*y[139];
y[178]=-(lambda*MYI*x1*y[18]*y[58]*y[115]*y[143]);
y[179]=y[177]+y[178];
y[180]=-(x0*x1*y[18]*y[24]*y[26]*y[58]*y[63]*y[68]*y[143]);
y[181]=lambda*MYI*x0*y[24]*y[26]*y[120]*y[139];
y[182]=y[180]+y[181];
y[183]=pow(y[68],2);
y[184]=x0*x1*y[18]*y[24]*y[26]*y[58]*y[63]*y[183];
y[185]=y[115]*y[120];
y[186]=y[184]+y[185];
y[187]=x2*y[1]*y[2]*y[3];
y[188]=x0*x1*x2*y[1]*y[2];
y[189]=2.*x0*x1*x3*y[1]*y[2];
y[190]=2.*x2*x3*y[1]*y[2]*y[3];
y[191]=x2*y[1]*y[3]*y[7];
y[192]=x0*x1*x2*y[1]*y[4];
y[193]=y[1]*y[3]*y[4]*y[14];
y[194]=x0*x1*x2*y[1]*y[5];
y[195]=-(x2*y[1]*y[3]*y[6]);
y[196]=y[10]+y[19]+y[20]+y[22]+y[61]+y[62]+y[64]+y[65]+y[66]+y[91]+y[96]+y[1\
34]+y[150]+y[156]+y[158]+y[187]+y[188]+y[189]+y[190]+y[191]+y[192]+y[193]+y\
[194]+y[195];
y[197]=-(lambda*MYI*x2*y[60]*y[94]*y[160]);
y[198]=x2+y[197];
y[199]=-(lambda*MYI*x1*y[18]*y[58]*y[93]);
y[200]=x1+y[199];
y[201]=-(lambda*MYI*x3*y[122]*y[161]*y[196]);
y[202]=x3+y[201];
y[203]=-(lambda*MYI*x0*y[24]*y[26]*y[57]);
y[204]=x0+y[203];
y[205]=pow(y[200],2);
y[206]=pow(y[204],2);
y[207]=pow(y[198],2);
y[208]=pow(y[202],2);
FOUT=(pow(bi,-2)*(-(lambda*MYI*x3*y[67]*y[122]*y[161]*(-(lambda*MYI*x2*y[60]\
*y[94]*y[139]*y[164])-y[146]*y[176]-lambda*MYI*x2*y[60]*y[94]*y[103]*y[179]\
))+lambda*MYI*x3*y[78]*y[122]*y[161]*(-(lambda*MYI*x2*y[60]*y[94]*y[143]*y[\
164])-y[167]*y[176]-lambda*MYI*x2*y[60]*y[94]*y[103]*y[182])+lambda*MYI*x3*\
y[103]*y[122]*y[161]*(lambda*MYI*x2*y[60]*y[94]*y[143]*y[146]-lambda*MYI*x2\
*y[60]*y[94]*y[139]*y[167]-lambda*MYI*x2*y[60]*y[94]*y[103]*y[186])+(lambda\
*MYI*x2*y[60]*y[94]*y[143]*y[179]-lambda*MYI*x2*y[60]*y[94]*y[139]*y[182]+y\
[176]*y[186])*(1.-lambda*MYI*x3*(2.*x0*y[1]*y[2]+2.*x0*x1*y[1]*y[2]+2.*x2*y\
[1]*y[2]*y[3])*y[122]*y[161]+lambda*MYI*x3*y[161]*y[196]-lambda*MYI*y[122]*\
y[161]*y[196])))/((y[1]+y[1]*y[198]+y[1]*y[200]+y[1]*y[198]*y[200]+y[1]*y[2\
02]+y[1]*y[200]*y[202]+y[1]*y[198]*y[204]+y[1]*y[198]*y[202]*y[204])*(y[134\
]+y[1]*y[5]*y[198]+2.*y[1]*y[5]*y[200]-y[1]*y[6]*y[200]+2.*y[1]*y[5]*y[198]\
*y[200]-y[1]*y[6]*y[198]*y[200]+y[1]*y[5]*y[202]+2.*y[1]*y[5]*y[200]*y[202]\
-y[1]*y[6]*y[200]*y[202]+y[1]*y[7]*y[204]+y[1]*y[4]*y[198]*y[204]+y[1]*y[5]\
*y[198]*y[204]+y[1]*y[7]*y[198]*y[204]+y[1]*y[7]*y[200]*y[204]+y[1]*y[4]*y[\
198]*y[200]*y[204]+y[1]*y[5]*y[198]*y[200]*y[204]-y[1]*y[6]*y[198]*y[200]*y\
[204]+y[1]*y[7]*y[198]*y[200]*y[204]+y[1]*y[2]*y[202]*y[204]-y[1]*y[6]*y[20\
2]*y[204]+y[1]*y[7]*y[202]*y[204]+y[1]*y[2]*y[198]*y[202]*y[204]+y[1]*y[4]*\
y[198]*y[202]*y[204]+y[1]*y[5]*y[198]*y[202]*y[204]-y[1]*y[6]*y[198]*y[202]\
*y[204]+y[1]*y[2]*y[200]*y[202]*y[204]-y[1]*y[6]*y[200]*y[202]*y[204]+y[1]*\
y[7]*y[200]*y[202]*y[204]+y[1]*y[2]*y[198]*y[200]*y[202]*y[204]+y[1]*y[4]*y\
[198]*y[200]*y[202]*y[204]+y[1]*y[5]*y[198]*y[200]*y[202]*y[204]+y[1]*y[5]*\
y[205]+y[1]*y[5]*y[198]*y[205]+y[1]*y[5]*y[202]*y[205]+y[1]*y[7]*y[198]*y[2\
06]+y[1]*y[2]*y[198]*y[202]*y[206]-y[1]*y[6]*y[198]*y[202]*y[206]+y[1]*y[7]\
*y[198]*y[202]*y[206]+y[1]*y[4]*y[204]*y[207]+y[1]*y[4]*y[200]*y[204]*y[207\
]+y[1]*y[4]*y[206]*y[207]+y[1]*y[4]*y[202]*y[206]*y[207]+y[1]*y[2]*y[204]*y\
[208]+y[1]*y[2]*y[200]*y[204]*y[208]+y[1]*y[2]*y[198]*y[206]*y[208]));
return (FOUT);
}
