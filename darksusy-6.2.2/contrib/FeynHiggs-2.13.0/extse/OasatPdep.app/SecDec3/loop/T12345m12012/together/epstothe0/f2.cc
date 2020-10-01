#include "intfile.hh"

dcmplx Pf2(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
double x3=x[3];
dcmplx y[149];
dcmplx FOUT;
dcmplx MYI(0.,1.);
y[1]=1./bi;
y[2]=em[0];
y[3]=em[1];
y[4]=esx[0];
y[5]=y[1]*y[2];
y[6]=x1*y[1]*y[2];
y[7]=x1*y[1]*y[3];
y[8]=-(x1*y[1]*y[4]);
y[9]=x3*y[1]*y[2];
y[10]=x3*y[1]*y[3];
y[11]=-(x3*y[1]*y[4]);
y[12]=-x1;
y[13]=1.+y[12];
y[14]=x0*y[1]*y[2];
y[15]=x2*y[1]*y[2];
y[16]=x0*y[1]*y[3];
y[17]=x0*x0;
y[18]=y[1]*y[3]*y[17];
y[19]=2.*x0*x1*y[1]*y[3];
y[20]=x2*y[1]*y[3];
y[21]=-(x0*y[1]*y[4]);
y[22]=-(x2*y[1]*y[4]);
y[23]=-x0;
y[24]=1.+y[23];
y[25]=2.*y[1]*y[3];
y[26]=2.*x1*x3*y[1]*y[3];
y[27]=lrs[0];
y[28]=x2*x3*y[1]*y[2];
y[29]=x2*x2;
y[30]=y[1]*y[3];
y[31]=2.*x1*y[1]*y[3];
y[32]=x1*x1;
y[33]=x2*x3*y[1]*y[3];
y[34]=-(y[1]*y[4]);
y[35]=-(x2*x3*y[1]*y[4]);
y[36]=x1*x3*y[1]*y[2];
y[37]=x1*x2*x3*y[1]*y[2];
y[38]=x3*y[1]*y[2]*y[29];
y[39]=2.*x0*y[1]*y[3];
y[40]=x1*x3*y[1]*y[3];
y[41]=2.*x0*x1*x3*y[1]*y[3];
y[42]=x3*y[1]*y[3]*y[32];
y[43]=2.*x0*x2*x3*y[1]*y[3];
y[44]=x1*x2*x3*y[1]*y[3];
y[45]=-(x1*x3*y[1]*y[4]);
y[46]=-(x1*x2*x3*y[1]*y[4]);
y[47]=y[5]+y[15]+y[20]+y[22]+y[28]+y[30]+y[31]+y[33]+y[34]+y[35]+y[36]+y[37]\
+y[38]+y[39]+y[40]+y[41]+y[42]+y[43]+y[44]+y[45]+y[46];
y[48]=lrs[1];
y[49]=-x2;
y[50]=1.+y[49];
y[51]=2.*x0*x3*y[1]*y[3];
y[52]=lambda*lambda;
y[53]=x0*x2*y[1]*y[2];
y[54]=x0*x2*y[1]*y[3];
y[55]=-(x0*x2*y[1]*y[4]);
y[56]=y[5]+y[14]+y[15]+y[16]+y[18]+y[19]+y[20]+y[21]+y[22]+y[31]+y[53]+y[54]\
+y[55];
y[57]=y[9]+y[10]+y[11]+y[25]+y[26]+y[28]+y[33]+y[35]+y[51];
y[58]=x1*x2*y[1]*y[2];
y[59]=y[1]*y[2]*y[29];
y[60]=y[1]*y[3]*y[32];
y[61]=2.*x0*x2*y[1]*y[3];
y[62]=x1*x2*y[1]*y[3];
y[63]=-(x1*x2*y[1]*y[4]);
y[64]=y[6]+y[7]+y[8]+y[15]+y[19]+y[20]+y[22]+y[58]+y[59]+y[60]+y[61]+y[62]+y\
[63];
y[65]=x0*x3*y[1]*y[2];
y[66]=x0*x3*y[1]*y[3];
y[67]=-(x0*x3*y[1]*y[4]);
y[68]=x0*x2*x3*y[1]*y[2];
y[69]=x3*y[1]*y[3]*y[17];
y[70]=x0*x2*x3*y[1]*y[3];
y[71]=-(x0*x2*x3*y[1]*y[4]);
y[72]=y[5]+y[9]+y[15]+y[20]+y[22]+y[26]+y[28]+y[30]+y[31]+y[33]+y[34]+y[35]+\
y[39]+y[41]+y[65]+y[66]+y[67]+y[68]+y[69]+y[70]+y[71];
y[73]=lrs[2];
y[74]=x0*x1*y[1]*y[2];
y[75]=2.*x2*y[1]*y[2];
y[76]=2.*x0*x2*y[1]*y[2];
y[77]=x0*x1*y[1]*y[3];
y[78]=-(x0*x1*y[1]*y[4]);
y[79]=y[5]+y[6]+y[7]+y[8]+y[14]+y[16]+y[18]+y[21]+y[74]+y[75]+y[76]+y[77]+y[\
78];
y[80]=2.*x2*x3*y[1]*y[3];
y[81]=y[25]+y[26]+y[80];
y[82]=-(lambda*MYI*x0*y[24]*y[27]*y[81]);
y[83]=-(lambda*MYI*y[24]*y[27]*y[47]);
y[84]=lambda*MYI*x0*y[27]*y[47];
y[85]=1.+y[82]+y[83]+y[84];
y[86]=2.*x3*y[1]*y[3];
y[87]=y[25]+y[51]+y[86];
y[88]=-(lambda*MYI*x1*y[13]*y[48]*y[87]);
y[89]=-(lambda*MYI*y[13]*y[48]*y[72]);
y[90]=lambda*MYI*x1*y[48]*y[72];
y[91]=1.+y[88]+y[89]+y[90];
y[92]=-x3;
y[93]=1.+y[92];
y[94]=2.*x2*x3*y[1]*y[2];
y[95]=y[5]+y[9]+y[10]+y[11]+y[30]+y[34]+y[36]+y[40]+y[45]+y[51]+y[94];
y[96]=y[5]+y[9]+y[10]+y[11]+y[30]+y[34]+y[65]+y[66]+y[67];
y[97]=x0*x1*y[13]*y[24]*y[27]*y[48]*y[52]*y[57]*y[64];
y[98]=-(lambda*MYI*x1*y[13]*y[48]*y[56]*y[85]);
y[99]=y[97]+y[98];
y[100]=2.*y[1]*y[2];
y[101]=x0*x1*x3*y[1]*y[2];
y[102]=2.*x0*x2*x3*y[1]*y[2];
y[103]=x0*x1*x3*y[1]*y[3];
y[104]=-(x0*x1*x3*y[1]*y[4]);
y[105]=y[6]+y[7]+y[8]+y[9]+y[14]+y[16]+y[21]+y[36]+y[40]+y[45]+y[65]+y[66]+y\
[67]+y[69]+y[75]+y[94]+y[100]+y[101]+y[102]+y[103]+y[104];
y[106]=lrs[3];
y[107]=x0*x1*y[13]*y[24]*y[27]*y[48]*y[52]*y[64]*y[96];
y[108]=-(x0*x1*y[13]*y[24]*y[27]*y[48]*y[52]*y[56]*y[95]);
y[109]=y[107]+y[108];
y[110]=-(x0*x1*y[13]*y[24]*y[27]*y[48]*y[52]*y[56]*y[57]);
y[111]=lambda*MYI*x0*y[24]*y[27]*y[64]*y[91];
y[112]=y[110]+y[111];
y[113]=2.*x3*y[1]*y[2];
y[114]=2.*x0*x3*y[1]*y[2];
y[115]=y[100]+y[113]+y[114];
y[116]=-(lambda*MYI*x2*y[50]*y[73]*y[115]);
y[117]=-(lambda*MYI*y[50]*y[73]*y[105]);
y[118]=lambda*MYI*x2*y[73]*y[105];
y[119]=1.+y[116]+y[117]+y[118];
y[120]=x0*x1*y[13]*y[24]*y[27]*y[48]*y[52]*y[57]*y[95];
y[121]=-(lambda*MYI*x1*y[13]*y[48]*y[85]*y[96]);
y[122]=y[120]+y[121];
y[123]=-(x0*x1*y[13]*y[24]*y[27]*y[48]*y[52]*y[57]*y[96]);
y[124]=lambda*MYI*x0*y[24]*y[27]*y[91]*y[95];
y[125]=y[123]+y[124];
y[126]=pow(y[57],2);
y[127]=x0*x1*y[13]*y[24]*y[27]*y[48]*y[52]*y[126];
y[128]=y[85]*y[91];
y[129]=y[127]+y[128];
y[130]=x0*x1*x2*y[1]*y[2];
y[131]=x0*y[1]*y[2]*y[29];
y[132]=x1*y[1]*y[3]*y[17];
y[133]=x0*y[1]*y[3]*y[32];
y[134]=x2*y[1]*y[3]*y[17];
y[135]=x0*x1*x2*y[1]*y[3];
y[136]=-(x0*x1*x2*y[1]*y[4]);
y[137]=y[6]+y[15]+y[53]+y[54]+y[55]+y[58]+y[59]+y[60]+y[62]+y[63]+y[74]+y[77\
]+y[78]+y[130]+y[131]+y[132]+y[133]+y[134]+y[135]+y[136];
y[138]=-(lambda*MYI*x1*y[13]*y[48]*y[72]);
y[139]=x1+y[138];
y[140]=-(lambda*MYI*x0*y[24]*y[27]*y[47]);
y[141]=x0+y[140];
y[142]=-(lambda*MYI*x3*y[93]*y[106]*y[137]);
y[143]=x3+y[142];
y[144]=-(lambda*MYI*x2*y[50]*y[73]*y[105]);
y[145]=x2+y[144];
y[146]=pow(y[141],2);
y[147]=pow(y[139],2);
y[148]=pow(y[145],2);
FOUT=(pow(bi,-2)*(-(lambda*MYI*x3*y[56]*y[93]*y[106]*(-(lambda*MYI*x2*y[50]*\
y[73]*y[95]*y[109])-y[99]*y[119]-lambda*MYI*x2*y[50]*y[73]*y[79]*y[122]))+l\
ambda*MYI*x3*y[64]*y[93]*y[106]*(-(lambda*MYI*x2*y[50]*y[73]*y[96]*y[109])-\
y[112]*y[119]-lambda*MYI*x2*y[50]*y[73]*y[79]*y[125])+lambda*MYI*x3*y[79]*y\
[93]*y[106]*(lambda*MYI*x2*y[50]*y[73]*y[96]*y[99]-lambda*MYI*x2*y[50]*y[73\
]*y[95]*y[112]-lambda*MYI*x2*y[50]*y[73]*y[79]*y[129])+(lambda*MYI*x2*y[50]\
*y[73]*y[96]*y[122]-lambda*MYI*x2*y[50]*y[73]*y[95]*y[125]+y[119]*y[129])*(\
1.+lambda*MYI*x3*y[106]*y[137]-lambda*MYI*y[93]*y[106]*y[137])))/((y[1]+y[1\
]*y[139]+y[1]*y[141]+y[1]*y[139]*y[143]+y[1]*y[139]*y[141]*y[143]+y[1]*y[14\
5]+y[1]*y[143]*y[145]+y[1]*y[141]*y[143]*y[145])*(y[5]+y[1]*y[2]*y[139]+y[1\
]*y[3]*y[139]-y[1]*y[4]*y[139]+y[1]*y[2]*y[141]+y[1]*y[3]*y[141]-y[1]*y[4]*\
y[141]+2.*y[1]*y[3]*y[139]*y[141]+y[1]*y[2]*y[139]*y[143]+y[1]*y[2]*y[139]*\
y[141]*y[143]+y[1]*y[3]*y[139]*y[141]*y[143]-y[1]*y[4]*y[139]*y[141]*y[143]\
+2.*y[1]*y[2]*y[145]+y[1]*y[2]*y[139]*y[145]+y[1]*y[3]*y[139]*y[145]-y[1]*y\
[4]*y[139]*y[145]+y[1]*y[2]*y[141]*y[145]+y[1]*y[3]*y[141]*y[145]-y[1]*y[4]\
*y[141]*y[145]+y[1]*y[2]*y[143]*y[145]+y[1]*y[2]*y[139]*y[143]*y[145]+y[1]*\
y[3]*y[139]*y[143]*y[145]-y[1]*y[4]*y[139]*y[143]*y[145]+y[1]*y[2]*y[141]*y\
[143]*y[145]+y[1]*y[3]*y[141]*y[143]*y[145]-y[1]*y[4]*y[141]*y[143]*y[145]+\
y[1]*y[2]*y[139]*y[141]*y[143]*y[145]+y[1]*y[3]*y[139]*y[141]*y[143]*y[145]\
-y[1]*y[4]*y[139]*y[141]*y[143]*y[145]+y[1]*y[3]*y[146]+y[1]*y[3]*y[139]*y[\
143]*y[146]+y[1]*y[3]*y[143]*y[145]*y[146]+y[1]*y[3]*y[147]+y[1]*y[3]*y[143\
]*y[147]+y[1]*y[3]*y[141]*y[143]*y[147]+y[1]*y[2]*y[148]+y[1]*y[2]*y[143]*y\
[148]+y[1]*y[2]*y[141]*y[143]*y[148]));
return (FOUT);
}