#include "intfile.hh"

dcmplx Pf10(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
dcmplx y[152];
dcmplx FOUT;
dcmplx MYI(0.,1.);
y[1]=1./bi;
y[2]=x0*x0;
y[3]=em[1];
y[4]=em[3];
y[5]=esx[0];
y[6]=em[0];
y[7]=em[2];
y[8]=y[1]*y[6];
y[9]=-x1;
y[10]=1.+y[9];
y[11]=y[1]*y[2]*y[6];
y[12]=x0*y[1]*y[3];
y[13]=y[1]*y[2]*y[3];
y[14]=2.*x2*y[1]*y[2]*y[3];
y[15]=x0*y[1]*y[7];
y[16]=x0*y[1]*y[4];
y[17]=2.*x1*y[1]*y[2]*y[4];
y[18]=-(x0*y[1]*y[5]);
y[19]=-(y[1]*y[2]*y[5]);
y[20]=y[11]+y[12]+y[13]+y[14]+y[15]+y[16]+y[17]+y[18]+y[19];
y[21]=-x0;
y[22]=1.+y[21];
y[23]=x2*x2;
y[24]=x1*x1;
y[25]=lrs[0];
y[26]=2.*x0*x1*y[1]*y[6];
y[27]=x2*y[1]*y[3];
y[28]=x1*y[1]*y[7];
y[29]=x1*y[1]*y[4];
y[30]=2.*x0*y[1]*y[4]*y[24];
y[31]=-(x2*y[1]*y[5]);
y[32]=x1*y[1]*y[6];
y[33]=x2*y[1]*y[6];
y[34]=2.*x0*x1*x2*y[1]*y[6];
y[35]=x1*x2*y[1]*y[3];
y[36]=2.*x0*x1*x2*y[1]*y[3];
y[37]=y[1]*y[3]*y[23];
y[38]=2.*x0*x1*y[1]*y[3]*y[23];
y[39]=x1*x2*y[1]*y[7];
y[40]=y[1]*y[4]*y[24];
y[41]=x1*x2*y[1]*y[4];
y[42]=2.*x0*x2*y[1]*y[4]*y[24];
y[43]=-(x1*x2*y[1]*y[5]);
y[44]=-2.*x0*x1*x2*y[1]*y[5];
y[45]=y[8]+y[26]+y[27]+y[28]+y[29]+y[30]+y[31]+y[32]+y[33]+y[34]+y[35]+y[36]\
+y[37]+y[38]+y[39]+y[40]+y[41]+y[42]+y[43]+y[44];
y[46]=lrs[1];
y[47]=-x2;
y[48]=1.+y[47];
y[49]=y[1]*y[3];
y[50]=x1*y[1]*y[3];
y[51]=2.*x0*x1*y[1]*y[3];
y[52]=2.*x2*y[1]*y[3];
y[53]=4.*x0*x1*x2*y[1]*y[3];
y[54]=-(y[1]*y[5]);
y[55]=-(x1*y[1]*y[5]);
y[56]=-2.*x0*x1*y[1]*y[5];
y[57]=y[8]+y[26]+y[28]+y[29]+y[30]+y[49]+y[50]+y[51]+y[52]+y[53]+y[54]+y[55]\
+y[56];
y[58]=lambda*lambda;
y[59]=2.*x0*y[1]*y[6];
y[60]=2.*x0*x2*y[1]*y[6];
y[61]=2.*x0*x2*y[1]*y[3];
y[62]=2.*x0*y[1]*y[3]*y[23];
y[63]=y[1]*y[7];
y[64]=x2*y[1]*y[7];
y[65]=y[1]*y[4];
y[66]=2.*x1*y[1]*y[4];
y[67]=4.*x0*x1*y[1]*y[4];
y[68]=x2*y[1]*y[4];
y[69]=4.*x0*x1*x2*y[1]*y[4];
y[70]=-2.*x0*x2*y[1]*y[5];
y[71]=y[8]+y[27]+y[31]+y[59]+y[60]+y[61]+y[62]+y[63]+y[64]+y[65]+y[66]+y[67]\
+y[68]+y[69]+y[70];
y[72]=x0*y[1]*y[6];
y[73]=x2*y[1]*y[2]*y[6];
y[74]=x0*x2*y[1]*y[3];
y[75]=x2*y[1]*y[2]*y[3];
y[76]=y[1]*y[2]*y[3]*y[23];
y[77]=x0*x2*y[1]*y[7];
y[78]=2.*x0*x1*y[1]*y[4];
y[79]=x0*x2*y[1]*y[4];
y[80]=2.*x1*x2*y[1]*y[2]*y[4];
y[81]=-(x0*x2*y[1]*y[5]);
y[82]=-(x2*y[1]*y[2]*y[5]);
y[83]=y[11]+y[15]+y[16]+y[17]+y[63]+y[72]+y[73]+y[74]+y[75]+y[76]+y[77]+y[78\
]+y[79]+y[80]+y[81]+y[82];
y[84]=lrs[2];
y[85]=2.*x1*y[1]*y[6];
y[86]=2.*x1*x2*y[1]*y[6];
y[87]=2.*x1*x2*y[1]*y[3];
y[88]=2.*x1*y[1]*y[3]*y[23];
y[89]=2.*y[1]*y[4]*y[24];
y[90]=2.*x2*y[1]*y[4]*y[24];
y[91]=-2.*x1*x2*y[1]*y[5];
y[92]=y[85]+y[86]+y[87]+y[88]+y[89]+y[90]+y[91];
y[93]=-(lambda*MYI*x0*y[22]*y[25]*y[92]);
y[94]=-(lambda*MYI*y[22]*y[25]*y[45]);
y[95]=lambda*MYI*x0*y[25]*y[45];
y[96]=1.+y[93]+y[94]+y[95];
y[97]=2.*x0*y[1]*y[4];
y[98]=2.*y[1]*y[2]*y[4];
y[99]=2.*x2*y[1]*y[2]*y[4];
y[100]=y[97]+y[98]+y[99];
y[101]=-(lambda*MYI*x1*y[10]*y[46]*y[100]);
y[102]=-(lambda*MYI*y[10]*y[46]*y[83]);
y[103]=lambda*MYI*x1*y[46]*y[83];
y[104]=1.+y[101]+y[102]+y[103];
y[105]=x1*y[1]*y[2]*y[6];
y[106]=x0*x1*y[1]*y[3];
y[107]=x1*y[1]*y[2]*y[3];
y[108]=2.*x1*x2*y[1]*y[2]*y[3];
y[109]=x0*x1*y[1]*y[7];
y[110]=x0*x1*y[1]*y[4];
y[111]=y[1]*y[2]*y[4]*y[24];
y[112]=-(x0*x1*y[1]*y[5]);
y[113]=-(x1*y[1]*y[2]*y[5]);
y[114]=y[12]+y[18]+y[61]+y[63]+y[72]+y[105]+y[106]+y[107]+y[108]+y[109]+y[11\
0]+y[111]+y[112]+y[113];
y[115]=-(lambda*MYI*x1*y[10]*y[46]*y[83]);
y[116]=x1+y[115];
y[117]=-(lambda*MYI*x0*y[22]*y[25]*y[45]);
y[118]=x0+y[117];
y[119]=-(lambda*MYI*x2*y[48]*y[84]*y[114]);
y[120]=x2+y[119];
y[121]=pow(bi,-2);
y[122]=x0*x1*y[10]*y[22]*y[25]*y[46]*y[57]*y[58]*y[71];
y[123]=-(lambda*MYI*x1*y[10]*y[20]*y[46]*y[96]);
y[124]=y[122]+y[123];
y[125]=lambda*MYI*x2*y[20]*y[48]*y[84]*y[124];
y[126]=-(x0*x1*y[10]*y[20]*y[22]*y[25]*y[46]*y[58]*y[71]);
y[127]=lambda*MYI*x0*y[22]*y[25]*y[57]*y[104];
y[128]=y[126]+y[127];
y[129]=-(lambda*MYI*x2*y[48]*y[57]*y[84]*y[128]);
y[130]=pow(y[71],2);
y[131]=x0*x1*y[10]*y[22]*y[25]*y[46]*y[58]*y[130];
y[132]=y[96]*y[104];
y[133]=y[131]+y[132];
y[134]=2.*x0*y[1]*y[3];
y[135]=2.*x1*y[1]*y[2]*y[3];
y[136]=y[134]+y[135];
y[137]=-(lambda*MYI*x2*y[48]*y[84]*y[136]);
y[138]=-(lambda*MYI*y[48]*y[84]*y[114]);
y[139]=lambda*MYI*x2*y[84]*y[114];
y[140]=1.+y[137]+y[138]+y[139];
y[141]=y[133]*y[140];
y[142]=y[125]+y[129]+y[141];
y[143]=y[1]*y[116];
y[144]=y[1]*y[116]*y[118];
y[145]=y[1]*y[120];
y[146]=y[1]*y[116]*y[118]*y[120];
y[147]=y[1]+y[143]+y[144]+y[145]+y[146];
y[148]=pow(y[147],-2);
y[149]=pow(y[118],2);
y[150]=pow(y[116],2);
y[151]=pow(y[120],2);
FOUT=myLog(bi)*y[121]*y[142]*y[148]+myLog(x0)*y[121]*y[142]*y[148]+myLog(1.+\
y[94])*y[121]*y[142]*y[148]+3.*myLog(y[147])*y[121]*y[142]*y[148]-2.*myLog(\
y[63]+y[1]*y[7]*y[116]+y[1]*y[6]*y[118]+y[1]*y[4]*y[116]*y[118]+y[1]*y[6]*y\
[116]*y[118]+y[1]*y[7]*y[116]*y[118]+y[1]*y[7]*y[120]+y[1]*y[3]*y[118]*y[12\
0]-y[1]*y[5]*y[118]*y[120]+y[1]*y[6]*y[118]*y[120]+y[1]*y[3]*y[116]*y[118]*\
y[120]+y[1]*y[4]*y[116]*y[118]*y[120]-y[1]*y[5]*y[116]*y[118]*y[120]+y[1]*y\
[7]*y[116]*y[118]*y[120]+y[1]*y[6]*y[116]*y[149]+y[1]*y[3]*y[116]*y[120]*y[\
149]-y[1]*y[5]*y[116]*y[120]*y[149]+y[1]*y[6]*y[116]*y[120]*y[149]+y[1]*y[4\
]*y[118]*y[150]+y[1]*y[4]*y[149]*y[150]+y[1]*y[4]*y[120]*y[149]*y[150]+y[1]\
*y[3]*y[118]*y[151]+y[1]*y[3]*y[116]*y[149]*y[151])*y[121]*y[142]*y[148];
return (FOUT);
}
