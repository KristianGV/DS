#include "intfile.hh"

dcmplx Pf1(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
dcmplx y[127];
dcmplx FOUT;
dcmplx MYI(0.,1.);
y[1]=pow(bi,-2);
y[2]=-x0;
y[3]=1.+y[2];
y[4]=em[0];
y[5]=1./bi;
y[6]=lrs[0];
y[7]=y[4]*y[5];
y[8]=x1*y[4]*y[5];
y[9]=esx[0];
y[10]=-(y[5]*y[9]);
y[11]=y[7]+y[8]+y[10];
y[12]=-x1;
y[13]=1.+y[12];
y[14]=lrs[1];
y[15]=x0*y[4]*y[5];
y[16]=y[7]+y[15];
y[17]=-(lambda*MYI*x0*y[3]*y[6]*y[11]);
y[18]=x0+y[17];
y[19]=-(lambda*MYI*x1*y[13]*y[14]*y[16]);
y[20]=x1+y[19];
y[21]=myLog(x1);
y[22]=1./x2;
y[23]=lambda*lambda;
y[24]=pow(y[4],2);
y[25]=x0*x1*y[1]*y[3]*y[6]*y[13]*y[14]*y[23]*y[24];
y[26]=-(lambda*MYI*y[3]*y[6]*y[11]);
y[27]=lambda*MYI*x0*y[6]*y[11];
y[28]=1.+y[26]+y[27];
y[29]=-(lambda*MYI*y[13]*y[14]*y[16]);
y[30]=lambda*MYI*x1*y[14]*y[16];
y[31]=1.+y[29]+y[30];
y[32]=y[28]*y[31];
y[33]=y[25]+y[32];
y[34]=y[5]*y[18];
y[35]=y[5]*y[20];
y[36]=y[5]*y[18]*y[20];
y[37]=y[5]+y[34]+y[35]+y[36];
y[38]=pow(y[37],-2);
y[39]=myLog(bi);
y[40]=y[1]*y[33]*y[38]*y[39];
y[41]=1.+y[29];
y[42]=myLog(y[41]);
y[43]=-(y[1]*y[33]*y[38]*y[42]);
y[44]=myLog(y[37]);
y[45]=3.*y[1]*y[33]*y[38]*y[44];
y[46]=y[4]*y[5]*y[18];
y[47]=-(y[5]*y[9]*y[18]);
y[48]=y[4]*y[5]*y[20];
y[49]=y[4]*y[5]*y[18]*y[20];
y[50]=y[7]+y[10]+y[46]+y[47]+y[48]+y[49];
y[51]=myLog(y[50]);
y[52]=-2.*y[1]*y[33]*y[38]*y[51];
y[53]=2.*x0*x1*y[4]*y[5];
y[54]=-(x0*x1*y[5]*y[9]);
y[55]=y[7]+y[8]+y[15]+y[53]+y[54];
y[56]=lrs[2];
y[57]=-(lambda*MYI*y[55]*y[56]);
y[58]=1.+y[57];
y[59]=myLog(y[58]);
y[60]=y[1]*y[33]*y[38]*y[59];
y[61]=y[40]+y[43]+y[45]+y[52]+y[60];
y[62]=y[1]*y[39];
y[63]=-(y[1]*y[42]);
y[64]=myLog(x2);
y[65]=-x2;
y[66]=1.+y[65];
y[67]=2.*x0*y[4]*y[5];
y[68]=2.*x0*x2*y[4]*y[5];
y[69]=-(x0*y[5]*y[9]);
y[70]=y[7]+y[67]+y[68]+y[69];
y[71]=2.*x1*x2*y[4]*y[5];
y[72]=x2*x2;
y[73]=x2*y[4]*y[5];
y[74]=x1*y[4]*y[5]*y[72];
y[75]=-(x1*x2*y[5]*y[9]);
y[76]=y[7]+y[8]+y[10]+y[71]+y[73]+y[74]+y[75];
y[77]=2.*x1*y[4]*y[5];
y[78]=-(x1*y[5]*y[9]);
y[79]=y[7]+y[71]+y[77]+y[78];
y[80]=2.*x2*y[4]*y[5];
y[81]=y[4]*y[5]*y[72];
y[82]=-(x2*y[5]*y[9]);
y[83]=y[7]+y[80]+y[81]+y[82];
y[84]=x0*y[4]*y[5]*y[72];
y[85]=-(x0*x2*y[5]*y[9]);
y[86]=y[7]+y[15]+y[68]+y[73]+y[84]+y[85];
y[87]=-(lambda*MYI*y[3]*y[6]*y[76]);
y[88]=lambda*MYI*x0*y[6]*y[76];
y[89]=1.+y[87]+y[88];
y[90]=-(lambda*MYI*y[13]*y[14]*y[86]);
y[91]=lambda*MYI*x1*y[14]*y[86];
y[92]=1.+y[90]+y[91];
y[93]=2.*x0*x1*x2*y[4]*y[5];
y[94]=y[7]+y[8]+y[15]+y[53]+y[54]+y[93];
y[95]=-(lambda*MYI*y[56]*y[66]*y[94]);
y[96]=-(lambda*MYI*x0*y[3]*y[6]*y[76]);
y[97]=x0+y[96];
y[98]=-(lambda*MYI*x1*y[13]*y[14]*y[86]);
y[99]=x1+y[98];
y[100]=1.+y[95];
y[101]=1./y[100];
y[102]=x0*x1*y[3]*y[6]*y[13]*y[14]*y[23]*y[79]*y[83];
y[103]=-(lambda*MYI*x1*y[13]*y[14]*y[70]*y[89]);
y[104]=y[102]+y[103];
y[105]=lambda*MYI*x2*y[56]*y[66]*y[70]*y[104];
y[106]=-(x0*x1*y[3]*y[6]*y[13]*y[14]*y[23]*y[70]*y[83]);
y[107]=lambda*MYI*x0*y[3]*y[6]*y[79]*y[92];
y[108]=y[106]+y[107];
y[109]=-(lambda*MYI*x2*y[56]*y[66]*y[79]*y[108]);
y[110]=pow(y[83],2);
y[111]=x0*x1*y[3]*y[6]*y[13]*y[14]*y[23]*y[110];
y[112]=y[89]*y[92];
y[113]=y[111]+y[112];
y[114]=-2.*lambda*MYI*x0*x1*x2*y[4]*y[5]*y[56]*y[66];
y[115]=lambda*MYI*x2*y[56]*y[94];
y[116]=1.+y[95]+y[114]+y[115];
y[117]=y[113]*y[116];
y[118]=y[105]+y[109]+y[117];
y[119]=y[5]*y[97];
y[120]=y[5]*y[99];
y[121]=y[5]*y[97]*y[99];
y[122]=-(lambda*MYI*x2*y[56]*y[66]*y[94]);
y[123]=x2+y[122];
y[124]=y[5]*y[97]*y[99]*y[123];
y[125]=y[5]+y[119]+y[120]+y[121]+y[124];
y[126]=pow(y[125],-2);
FOUT=0.5*pow(y[21],2)*y[1]*y[33]*y[38]-y[21]*y[61]-y[22]*y[61]+0.5*(y[1]*(9.\
*pow(y[44],2)*y[33]*y[38]+4.*pow(y[51],2)*y[33]*y[38]-12.*y[33]*y[38]*y[44]\
*y[51])+2.*(3.*y[33]*y[38]*y[44]-2.*y[33]*y[38]*y[51])*(y[1]*y[59]+y[62]+y[\
63])+y[33]*y[38]*(pow(y[39],2)*y[1]+pow(y[42],2)*y[1]+pow(y[59],2)*y[1]-2.*\
y[1]*y[39]*y[42]+2.*y[59]*(y[62]+y[63])))+y[1]*y[33]*y[38]*(y[21]*y[22]-y[2\
2]*y[64])+y[1]*(-(y[21]*y[22])+y[22]*y[64])*y[101]*y[118]*y[126]+y[22]*(-(m\
yLog(1.+y[90])*y[1]*y[101]*y[118]*y[126])+myLog(y[100])*y[1]*y[101]*y[118]*\
y[126]-2.*myLog(y[7]+y[10]+y[4]*y[5]*y[97]-y[5]*y[9]*y[97]+y[4]*y[5]*y[99]+\
y[4]*y[5]*y[97]*y[99]+pow(y[123],2)*y[4]*y[5]*y[97]*y[99]+y[4]*y[5]*y[123]+\
y[4]*y[5]*y[97]*y[123]+y[4]*y[5]*y[99]*y[123]+2.*y[4]*y[5]*y[97]*y[99]*y[12\
3]-y[5]*y[9]*y[97]*y[99]*y[123])*y[1]*y[101]*y[118]*y[126]+3.*myLog(y[125])\
*y[1]*y[101]*y[118]*y[126]+y[1]*y[39]*y[101]*y[118]*y[126]);
return (FOUT);
}