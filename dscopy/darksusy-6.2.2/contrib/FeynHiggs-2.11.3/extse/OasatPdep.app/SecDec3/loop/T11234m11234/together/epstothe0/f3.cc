#include "intfile.hh"

dcmplx Pf3(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
dcmplx y[133];
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
y[65]=em[2];
y[66]=x1*x1;
y[67]=em[3];
y[68]=1./x2;
y[69]=x0*x1*y[1]*y[2];
y[70]=x1*y[1]*y[65];
y[71]=x0*x1*y[1]*y[65];
y[72]=y[1]*y[65]*y[66];
y[73]=x0*y[1]*y[65]*y[66];
y[74]=y[1]*y[67];
y[75]=x0*y[1]*y[67];
y[76]=x1*y[1]*y[67];
y[77]=x0*x1*y[1]*y[67];
y[78]=lrs[2];
y[79]=-x2;
y[80]=1.+y[79];
y[81]=x2*y[1]*y[67];
y[82]=y[1]*y[65];
y[83]=x0*y[1]*y[65];
y[84]=2.*x1*y[1]*y[65];
y[85]=2.*x0*x1*y[1]*y[65];
y[86]=4.*x1*x2*y[1]*y[65];
y[87]=2.*x2*y[1]*y[67];
y[88]=y[9]+y[10]+y[21]+y[74]+y[75]+y[82]+y[83]+y[84]+y[85]+y[86]+y[87];
y[89]=x1*x2*y[1]*y[2];
y[90]=x1*x2*y[1]*y[65];
y[91]=x2*y[1]*y[65]*y[66];
y[92]=x1*x2*y[1]*y[67];
y[93]=y[7]+y[8]+y[9]+y[10]+y[12]+y[13]+y[14]+y[15]+y[81]+y[89]+y[90]+y[91]+y\
[92];
y[94]=-(lambda*MYI*y[6]*y[11]*y[93]);
y[95]=y[12]+y[70]+y[72]+y[74]+y[76];
y[96]=x2*y[1]*y[2];
y[97]=x2*y[1]*y[65];
y[98]=2.*x1*x2*y[1]*y[65];
y[99]=y[7]+y[8]+y[9]+y[10]+y[81]+y[96]+y[97]+y[98];
y[100]=x2*x2;
y[101]=x0*x2*y[1]*y[2];
y[102]=x2*y[1]*y[3];
y[103]=x0*x2*y[1]*y[65];
y[104]=2.*x0*x1*x2*y[1]*y[65];
y[105]=2.*x1*y[1]*y[65]*y[100];
y[106]=x0*x2*y[1]*y[67];
y[107]=y[1]*y[67]*y[100];
y[108]=-(x2*y[1]*y[4]);
y[109]=y[9]+y[21]+y[23]+y[24]+y[25]+y[81]+y[97]+y[98]+y[101]+y[102]+y[103]+y\
[104]+y[105]+y[106]+y[107]+y[108];
y[110]=lambda*MYI*x0*y[11]*y[93];
y[111]=1.+y[41]+y[94]+y[110];
y[112]=2.*x2*y[1]*y[65];
y[113]=2.*x0*x2*y[1]*y[65];
y[114]=2.*y[1]*y[65]*y[100];
y[115]=y[112]+y[113]+y[114];
y[116]=-(lambda*MYI*x1*y[19]*y[20]*y[115]);
y[117]=-(lambda*MYI*y[19]*y[20]*y[109]);
y[118]=lambda*MYI*x1*y[20]*y[109];
y[119]=1.+y[116]+y[117]+y[118];
y[120]=2.*x2*y[1]*y[65]*y[66];
y[121]=2.*x1*x2*y[1]*y[67];
y[122]=y[14]+y[15]+y[69]+y[70]+y[71]+y[72]+y[73]+y[74]+y[75]+y[76]+y[77]+y[1\
20]+y[121];
y[123]=-(lambda*MYI*y[78]*y[80]*y[122]);
y[124]=-(lambda*MYI*x0*y[6]*y[11]*y[93]);
y[125]=x0+y[124];
y[126]=-(lambda*MYI*x1*y[19]*y[20]*y[109]);
y[127]=x1+y[126];
y[128]=pow(y[125],2);
y[129]=-(lambda*MYI*x2*y[78]*y[80]*y[122]);
y[130]=x2+y[129];
y[131]=pow(y[127],2);
y[132]=pow(y[130],2);
FOUT=x0*(myLog(bi)*y[32]*y[33]*y[48]*y[53]*y[64]+3.*myLog(y[52])*y[32]*y[33]\
*y[48]*y[53]*y[64]-2.*myLog(y[63])*y[32]*y[33]*y[48]*y[53]*y[64]+myLog(1.-l\
ambda*MYI*(y[14]+y[15]+y[69]+y[70]+y[71]+y[72]+y[73]+y[74]+y[75]+y[76]+y[77\
])*y[78])*y[32]*y[33]*y[48]*y[53]*y[64])-x0*y[32]*y[33]*y[48]*y[53]*y[64]*y\
[68]+(x0*y[32]*y[68]*(1.+y[94])*(lambda*MYI*x2*y[78]*y[80]*y[88]*(x0*x1*y[6\
]*y[11]*y[19]*y[20]*y[34]*y[95]*y[99]-lambda*MYI*x1*y[19]*y[20]*y[88]*y[111\
])-lambda*MYI*x2*y[78]*y[80]*y[95]*(-(x0*x1*y[6]*y[11]*y[19]*y[20]*y[34]*y[\
88]*y[99])+lambda*MYI*x0*y[6]*y[11]*y[95]*y[119])+(x0*x1*pow(y[99],2)*y[6]*\
y[11]*y[19]*y[20]*y[34]+y[111]*y[119])*(1.-lambda*MYI*x2*(2.*y[1]*y[65]*y[6\
6]+2.*x1*y[1]*y[67])*y[78]*y[80]+lambda*MYI*x2*y[78]*y[122]+y[123])))/((1.+\
y[123])*(y[1]+y[1]*y[125]+y[1]*y[127]+y[1]*y[125]*y[127]+y[1]*y[127]*y[130]\
)*(y[9]+y[1]*y[2]*y[125]+y[1]*y[3]*y[125]-y[1]*y[4]*y[125]+y[1]*y[3]*y[127]\
+y[1]*y[2]*y[125]*y[127]+y[1]*y[3]*y[125]*y[127]-y[1]*y[4]*y[125]*y[127]+y[\
1]*y[2]*y[128]+y[1]*y[2]*y[127]*y[128]+y[1]*y[67]*y[130]+y[1]*y[67]*y[125]*\
y[130]+y[1]*y[3]*y[127]*y[130]-y[1]*y[4]*y[127]*y[130]+y[1]*y[65]*y[127]*y[\
130]+y[1]*y[67]*y[127]*y[130]+y[1]*y[2]*y[125]*y[127]*y[130]+y[1]*y[65]*y[1\
25]*y[127]*y[130]+y[1]*y[67]*y[125]*y[127]*y[130]+y[1]*y[65]*y[130]*y[131]+\
y[1]*y[65]*y[125]*y[130]*y[131]+y[1]*y[67]*y[127]*y[132]+y[1]*y[65]*y[131]*\
y[132]));
return (FOUT);
}
