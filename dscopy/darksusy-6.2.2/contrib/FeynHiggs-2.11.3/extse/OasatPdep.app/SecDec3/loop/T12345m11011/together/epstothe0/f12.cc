#include "intfile.hh"

dcmplx Pf12(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
double x3=x[3];
dcmplx y[126];
dcmplx FOUT;
dcmplx MYI(0.,1.);
y[1]=1./bi;
y[2]=em[0];
y[3]=esx[0];
y[4]=2.*y[1]*y[2];
y[5]=2.*x0*y[1]*y[2];
y[6]=2.*x2*y[1]*y[2];
y[7]=2.*x0*x3*y[1]*y[2];
y[8]=-(y[1]*y[3]);
y[9]=-(x0*y[1]*y[3]);
y[10]=2.*x3*y[1]*y[2];
y[11]=2.*x2*x3*y[1]*y[2];
y[12]=-(x2*y[1]*y[3]);
y[13]=-x1;
y[14]=1.+y[13];
y[15]=x0*x0;
y[16]=-x0;
y[17]=1.+y[16];
y[18]=2.*x1*y[1]*y[2];
y[19]=x3*x3;
y[20]=lrs[0];
y[21]=2.*x0*x1*y[1]*y[2];
y[22]=2.*x1*x2*y[1]*y[2];
y[23]=4.*x0*x1*x3*y[1]*y[2];
y[24]=-(x3*y[1]*y[3]);
y[25]=y[1]*y[2];
y[26]=x2*y[1]*y[2];
y[27]=2.*x1*x3*y[1]*y[2];
y[28]=2.*x1*x2*x3*y[1]*y[2];
y[29]=y[1]*y[2]*y[19];
y[30]=2.*x0*x1*y[1]*y[2]*y[19];
y[31]=x2*y[1]*y[2]*y[19];
y[32]=-(x1*x2*y[1]*y[3]);
y[33]=-(x1*x3*y[1]*y[3]);
y[34]=-2.*x0*x1*x3*y[1]*y[3];
y[35]=-(x2*x3*y[1]*y[3]);
y[36]=y[10]+y[11]+y[18]+y[21]+y[22]+y[23]+y[24]+y[25]+y[26]+y[27]+y[28]+y[29\
]+y[30]+y[31]+y[32]+y[33]+y[34]+y[35];
y[37]=lrs[1];
y[38]=-x2;
y[39]=1.+y[38];
y[40]=-(x1*y[1]*y[3]);
y[41]=lambda*lambda;
y[42]=2.*y[1]*y[2]*y[15];
y[43]=2.*x0*x2*y[1]*y[2];
y[44]=2.*x3*y[1]*y[2]*y[15];
y[45]=-(y[1]*y[3]*y[15]);
y[46]=y[5]+y[9]+y[42]+y[43]+y[44]+y[45];
y[47]=4.*x0*x3*y[1]*y[2];
y[48]=2.*x0*y[1]*y[2]*y[19];
y[49]=-2.*x0*x3*y[1]*y[3];
y[50]=y[4]+y[5]+y[6]+y[10]+y[11]+y[12]+y[24]+y[47]+y[48]+y[49];
y[51]=4.*x0*x1*y[1]*y[2];
y[52]=-2.*x0*x1*y[1]*y[3];
y[53]=y[4]+y[6]+y[8]+y[10]+y[11]+y[12]+y[18]+y[22]+y[23]+y[40]+y[51]+y[52];
y[54]=y[1]*y[2]*y[15];
y[55]=x2*x2;
y[56]=y[1]*y[2]*y[55];
y[57]=2.*x0*x2*x3*y[1]*y[2];
y[58]=y[1]*y[2]*y[15]*y[19];
y[59]=-(x0*x2*y[1]*y[3]);
y[60]=-(x0*x3*y[1]*y[3]);
y[61]=-(x3*y[1]*y[3]*y[15]);
y[62]=y[5]+y[6]+y[7]+y[12]+y[25]+y[43]+y[44]+y[54]+y[56]+y[57]+y[58]+y[59]+y\
[60]+y[61];
y[63]=lrs[2];
y[64]=y[4]+y[5]+y[6]+y[7]+y[8]+y[9]+y[21];
y[65]=4.*x1*x3*y[1]*y[2];
y[66]=2.*x1*y[1]*y[2]*y[19];
y[67]=-2.*x1*x3*y[1]*y[3];
y[68]=y[18]+y[65]+y[66]+y[67];
y[69]=-(lambda*MYI*x0*y[17]*y[20]*y[68]);
y[70]=-(lambda*MYI*y[17]*y[20]*y[36]);
y[71]=lambda*MYI*x0*y[20]*y[36];
y[72]=1.+y[69]+y[70]+y[71];
y[73]=-(lambda*MYI*y[14]*y[37]*y[62]);
y[74]=lambda*MYI*x1*y[37]*y[62];
y[75]=1.+y[73]+y[74];
y[76]=-x3;
y[77]=1.+y[76];
y[78]=y[10]+y[18]+y[24]+y[25]+y[27]+y[29]+y[40];
y[79]=y[4]+y[5]+y[6]+y[7]+y[8]+y[9];
y[80]=x0*x1*y[14]*y[17]*y[20]*y[37]*y[41]*y[50]*y[53];
y[81]=-(lambda*MYI*x1*y[14]*y[37]*y[46]*y[72]);
y[82]=y[80]+y[81];
y[83]=x0*y[1]*y[2];
y[84]=2.*x0*x1*x3*y[1]*y[2];
y[85]=x0*y[1]*y[2]*y[19];
y[86]=-(x0*x1*y[1]*y[3]);
y[87]=y[4]+y[6]+y[7]+y[8]+y[10]+y[11]+y[18]+y[21]+y[22]+y[24]+y[40]+y[60]+y[\
83]+y[84]+y[85]+y[86];
y[88]=lrs[3];
y[89]=x0*x1*y[14]*y[17]*y[20]*y[37]*y[41]*y[53]*y[79];
y[90]=-(x0*x1*y[14]*y[17]*y[20]*y[37]*y[41]*y[46]*y[78]);
y[91]=y[89]+y[90];
y[92]=-(x0*x1*y[14]*y[17]*y[20]*y[37]*y[41]*y[46]*y[50]);
y[93]=lambda*MYI*x0*y[17]*y[20]*y[53]*y[75];
y[94]=y[92]+y[93];
y[95]=y[4]+y[10]+y[18];
y[96]=-(lambda*MYI*x2*y[39]*y[63]*y[95]);
y[97]=-(lambda*MYI*y[39]*y[63]*y[87]);
y[98]=lambda*MYI*x2*y[63]*y[87];
y[99]=1.+y[96]+y[97]+y[98];
y[100]=x0*x1*y[14]*y[17]*y[20]*y[37]*y[41]*y[50]*y[78];
y[101]=-(lambda*MYI*x1*y[14]*y[37]*y[72]*y[79]);
y[102]=y[100]+y[101];
y[103]=-(x0*x1*y[14]*y[17]*y[20]*y[37]*y[41]*y[50]*y[79]);
y[104]=lambda*MYI*x0*y[17]*y[20]*y[75]*y[78];
y[105]=y[103]+y[104];
y[106]=pow(y[50],2);
y[107]=x0*x1*y[14]*y[17]*y[20]*y[37]*y[41]*y[106];
y[108]=y[72]*y[75];
y[109]=y[107]+y[108];
y[110]=2.*x1*y[1]*y[2]*y[15];
y[111]=2.*x0*x1*x2*y[1]*y[2];
y[112]=2.*x1*x3*y[1]*y[2]*y[15];
y[113]=-(x1*y[1]*y[3]*y[15]);
y[114]=y[5]+y[6]+y[7]+y[9]+y[12]+y[21]+y[25]+y[43]+y[56]+y[57]+y[59]+y[86]+y\
[110]+y[111]+y[112]+y[113];
y[115]=-(lambda*MYI*x1*y[14]*y[37]*y[62]);
y[116]=x1+y[115];
y[117]=-(lambda*MYI*x2*y[39]*y[63]*y[87]);
y[118]=x2+y[117];
y[119]=-(lambda*MYI*x0*y[17]*y[20]*y[36]);
y[120]=x0+y[119];
y[121]=-(lambda*MYI*x3*y[77]*y[88]*y[114]);
y[122]=x3+y[121];
y[123]=pow(y[118],2);
y[124]=pow(y[120],2);
y[125]=pow(y[122],2);
FOUT=(pow(bi,-2)*(-(lambda*MYI*x3*y[46]*y[77]*y[88]*(-(lambda*MYI*x2*y[39]*y\
[63]*y[78]*y[91])-y[82]*y[99]-lambda*MYI*x2*y[39]*y[63]*y[64]*y[102]))+lamb\
da*MYI*x3*y[53]*y[77]*y[88]*(-(lambda*MYI*x2*y[39]*y[63]*y[79]*y[91])-y[94]\
*y[99]-lambda*MYI*x2*y[39]*y[63]*y[64]*y[105])+lambda*MYI*x3*y[64]*y[77]*y[\
88]*(lambda*MYI*x2*y[39]*y[63]*y[79]*y[82]-lambda*MYI*x2*y[39]*y[63]*y[78]*\
y[94]-lambda*MYI*x2*y[39]*y[63]*y[64]*y[109])+(lambda*MYI*x2*y[39]*y[63]*y[\
79]*y[102]-lambda*MYI*x2*y[39]*y[63]*y[78]*y[105]+y[99]*y[109])*(1.-lambda*\
MYI*x3*y[77]*y[88]*(y[5]+y[43]+y[110])+lambda*MYI*x3*y[88]*y[114]-lambda*MY\
I*y[77]*y[88]*y[114])))/((y[1]+y[1]*y[116]+y[1]*y[118]+y[1]*y[116]*y[118]+y\
[1]*y[116]*y[120]+y[1]*y[122]+y[1]*y[118]*y[122]+y[1]*y[116]*y[120]*y[122])\
*(y[25]+y[1]*y[2]*y[116]+2.*y[1]*y[2]*y[118]-y[1]*y[3]*y[118]+2.*y[1]*y[2]*\
y[116]*y[118]-y[1]*y[3]*y[116]*y[118]+y[1]*y[2]*y[120]+2.*y[1]*y[2]*y[116]*\
y[120]+y[1]*y[2]*y[118]*y[120]+2.*y[1]*y[2]*y[116]*y[118]*y[120]-y[1]*y[3]*\
y[116]*y[118]*y[120]+y[1]*y[2]*y[122]+2.*y[1]*y[2]*y[118]*y[122]-y[1]*y[3]*\
y[118]*y[122]+2.*y[1]*y[2]*y[120]*y[122]-y[1]*y[3]*y[120]*y[122]+2.*y[1]*y[\
2]*y[116]*y[120]*y[122]-y[1]*y[3]*y[116]*y[120]*y[122]+2.*y[1]*y[2]*y[118]*\
y[120]*y[122]-y[1]*y[3]*y[118]*y[120]*y[122]+2.*y[1]*y[2]*y[116]*y[118]*y[1\
20]*y[122]+y[1]*y[2]*y[123]+y[1]*y[2]*y[116]*y[123]+y[1]*y[2]*y[122]*y[123]\
+y[1]*y[2]*y[116]*y[124]+2.*y[1]*y[2]*y[116]*y[122]*y[124]-y[1]*y[3]*y[116]\
*y[122]*y[124]+y[1]*y[2]*y[120]*y[125]+y[1]*y[2]*y[118]*y[120]*y[125]+y[1]*\
y[2]*y[116]*y[124]*y[125]));
return (FOUT);
}
