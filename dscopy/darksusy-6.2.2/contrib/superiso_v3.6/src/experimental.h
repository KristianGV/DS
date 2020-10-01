central[0]=5.2e-2;
errors[0]=2.6e-2;
names[0]="delta0                      ";

central[1]=3.32e-4;
errors[1]=0.15e-4;
names[1]="bsgamma                     ";

central[2]=3.0e-9;
errors[2]=sqrt(0.6e-9*0.6e-9+0.3e-9*0.3e-9);
names[2]="Bsmumu_untag                ";

central[3]=-0.;
errors[3]=-0.;
names[3]="Bsee_untag                  ";

central[4]=1.5e-10;
errors[4]=sqrt(1.2e-10*1.2e-10+0.2e-10*0.2e-10);
names[4]="Bdmumu                      ";

central[5]=0.66e-6;
errors[5]=0.88e-6;
names[5]="BRBXsmumu_lowq2             ";

central[6]=0.6e-6;
errors[6]=0.31e-6;
names[6]="BRBXsmumu_highq2            ";

central[7]=-0.;
errors[7]=-0.;
names[7]="BRBXsmumu_full              ";

central[8]=1.93e-6;
errors[8]=0.55e-6;
names[8]="BRBXsee_lowq2               ";

central[9]=0.56e-6;
errors[9]=0.19e-6;
names[9]="BRBXsee_highq2              ";

central[10]=-0.;
errors[10]=-0.;
names[10]="BRBXsee_full                ";

central[11]=0.6e-7;
errors[11]=0.1e-7;
names[11]="dG_B0Kstar0mumu_0.1-2       ";

central[12]=0.3e-7;
errors[12]=0.05e-7;
names[12]="dG_B0Kstar0mumu_2-4.3       ";

central[13]=0.49e-7;
errors[13]=0.08e-7;
names[13]="dG_B0Kstar0mumu_4.3-8.68    ";

central[14]=0.56e-7;
errors[14]=0.1e-7;
names[14]="dG_B0Kstar0mumu_14.18-16    ";

central[15]=0.41e-7;
errors[15]=0.07e-7;
names[15]="dG_B0Kstar0mumu_16-19       ";

central[16]=-0.;
errors[16]=-0.;
names[16]="dG_B0Kstar0mumu_0.1-0.98    ";

if(likelihood==1)
{
central[17]=0.263;
errors[17]=sqrt(0.046*0.046+0.017*0.017);
names[17]="FL_B0Kstar0mumu_0.1-0.98    ";

central[18]=-0.003;
errors[18]=sqrt(0.059*0.059+0.008*0.008);
names[18]="AFB_B0Kstar0mumu_0.1-0.98   ";

central[19]=-0.036;
errors[19]=sqrt(0.063*0.063+0.005*0.005);
names[19]="S3_B0Kstar0mumu_0.1-0.98    ";

central[20]=-0.082;
errors[20]=sqrt(0.07*0.07+0.009*0.009);
names[20]="S4_B0Kstar0mumu_0.1-0.98    ";

central[21]=0.17;
errors[21]=sqrt(0.06*0.06+0.018*0.018);
names[21]="S5_B0Kstar0mumu_0.1-0.98    ";

central[22]=0.015;
errors[22]=sqrt(0.059*0.059+0.006*0.006);
names[22]="S7_B0Kstar0mumu_0.1-0.98    ";

central[23]=0.079;
errors[23]=sqrt(0.078*0.078+0.007*0.007);
names[23]="S8_B0Kstar0mumu_0.1-0.98    ";

central[24]=-0.083;
errors[24]=sqrt(0.06*0.06+0.004*0.004);
names[24]="S9_B0Kstar0mumu_0.1-0.98    ";

central[25]=0.387;
errors[25]=sqrt(0.141*0.141+0.052*0.052);
names[25]="P5p_B0Kstar0mumu_0.1-0.98   ";
}
else
{
central[17]=0.242;
errors[17]=sqrt(0.058*0.058+0.026*0.026);
names[17]="FL_B0Kstar0mumu_0.1-0.98    ";

central[18]=-0.138;
errors[18]=sqrt(0.095*0.095+0.072*0.072);
names[18]="AFB_B0Kstar0mumu_0.1-0.98   ";

central[19]=-0.014;
errors[19]=sqrt(0.060*0.060+0.008*0.008);
names[19]="S3_B0Kstar0mumu_0.1-0.98    ";

central[20]=-0.039;
errors[20]=sqrt(0.091*0.091+0.015*0.015);
names[20]="S4_B0Kstar0mumu_0.1-0.98    ";

central[21]=0.129;
errors[21]=sqrt(0.068*0.068+0.011*0.011);
names[21]="S5_B0Kstar0mumu_0.1-0.98    ";

central[22]=0.038;
errors[22]=sqrt(0.063*0.063+0.009*0.009);
names[22]="S7_B0Kstar0mumu_0.1-0.98    ";

central[23]=0.063;
errors[23]=sqrt(0.080*0.080+0.009*0.009);
names[23]="S8_B0Kstar0mumu_0.1-0.98    ";

central[24]=-0.113;
errors[24]=sqrt(0.063*0.063+0.004*0.004);
names[24]="S9_B0Kstar0mumu_0.1-0.98    ";

central[25]=0.300;
errors[25]=sqrt(0.171*0.171+0.023*0.023);
names[25]="P5p_B0Kstar0mumu_0.1-0.98   ";
}

central[26]=-0.;
errors[26]=-0.;
names[26]="dG_B0Kstar0mumu_1.1-2.5     ";

central[27]=0.66;
errors[27]=sqrt(0.088*0.088+0.022*0.022);
names[27]="FL_B0Kstar0mumu_1.1-2.5     ";

central[28]=-0.191;
errors[28]=sqrt(0.078*0.078+0.012*0.012);
names[28]="AFB_B0Kstar0mumu_1.1-2.5    ";

central[29]=-0.077;
errors[29]=sqrt(0.104*0.104+0.005*0.005);
names[29]="S3_B0Kstar0mumu_1.1-2.5     ";

central[30]=0.077;
errors[30]=sqrt(0.112*0.112+0.005*0.005);
names[30]="S4_B0Kstar0mumu_1.1-2.5     ";

central[31]=0.137;
errors[31]=sqrt(0.098*0.098+0.009*0.009);
names[31]="S5_B0Kstar0mumu_1.1-2.5     ";

central[32]=-0.219;
errors[32]=sqrt(0.105*0.105+0.003*0.003);
names[32]="S7_B0Kstar0mumu_1.1-2.5     ";

central[33]=-0.098;
errors[33]=sqrt(0.122*0.122+0.005*0.005);
names[33]="S8_B0Kstar0mumu_1.1-2.5     ";

central[34]=-0.119;
errors[34]=sqrt(0.101+0.101+0.005*0.005);
names[34]="S9_B0Kstar0mumu_1.1-2.5     ";

central[35]=0.289;
errors[35]=sqrt(0.216*0.216+0.023*0.023);
names[35]="P5p_B0Kstar0mumu_1.1-2.5    ";

central[36]=-0.;
errors[36]=-0.;
names[36]="dG_B0Kstar0mumu_2.5-4       ";

central[37]=0.877;
errors[37]=sqrt(0.096*0.096+0.017*0.017);
names[37]="FL_B0Kstar0mumu_2.5-4       ";

central[38]=-0.118;
errors[38]=sqrt(0.088*0.088+0.007*0.007);
names[38]="AFB_B0Kstar0mumu_2.5-4      ";

central[39]=0.035;
errors[39]=sqrt(0.101*0.101+0.006*0.006);
names[39]="S3_B0Kstar0mumu_2.5-4       ";

central[40]=0.234;
errors[40]=sqrt(0.144*0.144+0.006*0.006);
names[40]="S4_B0Kstar0mumu_2.5-4       ";

central[41]=-0.022;
errors[41]=sqrt(0.11*0.11+0.008*0.008);
names[41]="S5_B0Kstar0mumu_2.5-4       ";

central[42]=0.068;
errors[42]=sqrt(0.119*0.119+0.005*0.005);
names[42]="S7_B0Kstar0mumu_2.5-4       ";

central[43]=0.03;
errors[43]=sqrt(0.127*0.127+0.006*0.006);
names[43]="S8_B0Kstar0mumu_2.5-4       ";

central[44]=-0.092;
errors[44]=sqrt(0.125*0.125+0.007*0.007);
names[44]="S9_B0Kstar0mumu_2.5-4       ";

central[45]=-0.066;
errors[45]=sqrt(0.36*0.36+0.023*0.023);
names[45]="P5p_B0Kstar0mumu_2.5-4      ";

central[46]=-0.;
errors[46]=-0.;
names[46]="dG_B0Kstar0mumu_4-6         ";

central[47]=0.611;
errors[47]=sqrt(0.054*0.054+0.017*0.017);
names[47]="FL_B0Kstar0mumu_4-6         ";

central[48]=0.025;
errors[48]=sqrt(0.05*0.05+0.004*0.004);
names[48]="AFB_B0Kstar0mumu_4-6        ";

central[49]=0.035;
errors[49]=sqrt(0.069*0.069+0.006*0.006);
names[49]="S3_B0Kstar0mumu_4-6         ";

central[50]=0.219;
errors[50]=sqrt(0.086*0.086+0.008*0.008);
names[50]="S4_B0Kstar0mumu_4-6         ";

central[51]=-0.146;
errors[51]=sqrt(0.078*0.078+0.011*0.011);
names[51]="S5_B0Kstar0mumu_4-6         ";

central[52]=-0.016;
errors[52]=sqrt(0.083*0.083+0.004*0.004);
names[52]="S7_B0Kstar0mumu_4-6         ";

central[53]=0.167;
errors[53]=sqrt(0.097*0.097+0.004*0.004);
names[53]="S8_B0Kstar0mumu_4-6         ";

central[54]=-0.032;
errors[54]=sqrt(0.073*0.073+0.004*0.004);
names[54]="S9_B0Kstar0mumu_4-6         ";

central[55]=-0.3;
errors[55]=sqrt(0.156*0.156+0.023*0.023);
names[55]="P5p_B0Kstar0mumu_4-6        ";

central[56]=-0.;
errors[56]=-0.;
names[56]="dG_B0Kstar0mumu_6-8         ";

central[57]=0.579;
errors[57]=sqrt(0.047*0.047+0.015*0.015);
names[57]="FL_B0Kstar0mumu_6-8         ";

central[58]=0.152;
errors[58]=sqrt(0.04*0.04+0.008*0.008);
names[58]="AFB_B0Kstar0mumu_6-8        ";

central[59]=-0.042;
errors[59]=sqrt(0.058*0.058+0.011*0.011);
names[59]="S3_B0Kstar0mumu_6-8         ";

central[60]=0.296;
errors[60]=sqrt(0.065*0.065+0.011*0.011);
names[60]="S4_B0Kstar0mumu_6-8         ";

central[61]=-0.249;
errors[61]=sqrt(0.062*0.062+0.012*0.012);
names[61]="S5_B0Kstar0mumu_6-8         ";

central[62]=-0.047;
errors[62]=sqrt(0.066*0.066+0.003*0.003);
names[62]="S7_B0Kstar0mumu_6-8         ";

central[63]=-0.085;
errors[63]=sqrt(0.073*0.073+0.006*0.006);
names[63]="S8_B0Kstar0mumu_6-8         ";

central[64]=-0.024;
errors[64]=sqrt(0.062*0.062+0.005*0.005);
names[64]="S9_B0Kstar0mumu_6-8         ";

central[65]=-0.505;
errors[65]=sqrt(0.118*0.118+0.024*0.024);
names[65]="P5p_B0Kstar0mumu_6-8        ";

central[66]=-0.;
errors[66]=-0.;
names[66]="dG_B0Kstar0mumu_11-12.5     ";

central[67]=0.493;
errors[67]=sqrt(0.051*0.051+0.013*0.013);
names[67]="FL_B0Kstar0mumu_11-12.5     ";

central[68]=0.318;
errors[68]=sqrt(0.04*0.04+0.009*0.009);
names[68]="AFB_B0Kstar0mumu_11-12.5    ";

central[69]=-0.189;
errors[69]=sqrt(0.053*0.053+0.005*0.005);
names[69]="S3_B0Kstar0mumu_11-12.5     ";

central[70]=0.283;
errors[70]=sqrt(0.092*0.092+0.009*0.009);
names[70]="S4_B0Kstar0mumu_11-12.5     ";

central[71]=-0.327;
errors[71]=sqrt(0.08*0.08+0.009*0.009);
names[71]="S5_B0Kstar0mumu_11-12.5     ";

central[72]=-0.141;
errors[72]=sqrt(0.075*0.075+0.005*0.005);
names[72]="S7_B0Kstar0mumu_11-12.5     ";

central[73]=-0.007;
errors[73]=sqrt(0.073*0.073+0.005*0.005);
names[73]="S8_B0Kstar0mumu_11-12.5     ";

central[74]=-0.004;
errors[74]=sqrt(0.073*0.073+0.006*0.006);
names[74]="S9_B0Kstar0mumu_11-12.5     ";

central[75]=-0.654;
errors[75]=sqrt(0.153*0.153+0.015*0.015);
names[75]="P5p_B0Kstar0mumu_11-12.5    ";

central[76]=-0.;
errors[76]=-0.;
names[76]="dG_B0Kstar0mumu_15-17       ";

central[77]=0.349;
errors[77]=sqrt(0.04*0.04+0.009*0.009);
names[77]="FL_B0Kstar0mumu_15-17       ";

central[78]=0.411;
errors[78]=sqrt(0.04*0.04+0.008*0.008);
names[78]="AFB_B0Kstar0mumu_15-17      ";

central[79]=-0.142;
errors[79]=sqrt(0.047*0.047+0.007*0.007);
names[79]="S3_B0Kstar0mumu_15-17       ";

central[80]=0.321;
errors[80]=sqrt(0.078*0.078+0.007*0.007);
names[80]="S4_B0Kstar0mumu_15-17       ";

central[81]=-0.316;
errors[81]=sqrt(0.058*0.058+0.009*0.009);
names[81]="S5_B0Kstar0mumu_15-17       ";

central[82]=0.061;
errors[82]=sqrt(0.06*0.06+0.005*0.005);
names[82]="S7_B0Kstar0mumu_15-17       ";

central[83]=0.003;
errors[83]=sqrt(0.06*0.06+0.003*0.003);
names[83]="S8_B0Kstar0mumu_15-17       ";

central[84]=-0.019;
errors[84]=sqrt(0.057*0.057+0.004*0.004);
names[84]="S9_B0Kstar0mumu_15-17       ";

central[85]=-0.662;
errors[85]=sqrt(0.126*0.126+0.017*0.017);
names[85]="P5p_B0Kstar0mumu_15-17      ";

central[86]=-0.;
errors[86]=-0.;
names[86]="dG_B0Kstar0mumu_17-19       ";

central[87]=0.354;
errors[87]=sqrt(0.048*0.048+0.025*0.025);
names[87]="FL_B0Kstar0mumu_17-19       ";

central[88]=0.305;
errors[88]=sqrt(0.048*0.048+0.013*0.013);
names[88]="AFB_B0Kstar0mumu_17-19      ";

central[89]=-0.188;
errors[89]=sqrt(0.086*0.086+0.017*0.017);
names[89]="S3_B0Kstar0mumu_17-19       ";

central[90]=0.266;
errors[90]=sqrt(0.071*0.071+0.01*0.01);
names[90]="S4_B0Kstar0mumu_17-19       ";

central[91]=-0.323;
errors[91]=sqrt(0.069*0.069+0.009*0.009);
names[91]="S5_B0Kstar0mumu_17-19       ";

central[92]=0.044;
errors[92]=sqrt(0.073*0.073+0.013*0.013);
names[92]="S7_B0Kstar0mumu_17-19       ";

central[93]=0.013;
errors[93]=sqrt(0.071*0.071+0.005*0.005);
names[93]="S8_B0Kstar0mumu_17-19       ";

central[94]=-0.094;
errors[94]=sqrt(0.069*0.069+0.004*0.004);
names[94]="S9_B0Kstar0mumu_17-19       ";

central[95]=-0.675;
errors[95]=sqrt(0.152*0.152+0.017*0.017);
names[95]="P5p_B0Kstar0mumu_17-19      ";

central[96]=0.;
errors[96]=0.;
names[96]="dG_B0Kstar0mumu_1.1-6       ";

central[97]=0.690;
errors[97]=sqrt(0.037*0.037+0.017*0.017);
names[97]="FL_B0Kstar0mumu_1.1-6       ";

central[98]=-0.075;
errors[98]=sqrt(0.035*0.035+0.007*0.007);
names[98]="AFB_B0Kstar0mumu_1.1-6      ";

central[99]=0.012;
errors[99]=sqrt(0.038*0.038+0.004*0.004);
names[99]="S3_B0Kstar0mumu_1.1-6       ";

central[100]=0.155;
errors[100]=sqrt(0.058*0.058+0.004*0.004);
names[100]="S4_B0Kstar0mumu_1.1-6       ";

central[101]=-0.023;
errors[101]=sqrt(0.052*0.052+0.005*0.005);
names[101]="S5_B0Kstar0mumu_1.1-6       ";

central[102]=-0.077;
errors[102]=sqrt(0.052*0.052+0.006*0.006);
names[102]="S7_B0Kstar0mumu_1.1-6       ";

central[103]=0.028;
errors[103]=sqrt(0.058*0.058+0.008*0.008);
names[103]="S8_B0Kstar0mumu_1.1-6       ";

central[104]=-0.064;
errors[104]=sqrt(0.042*0.042+0.004*0.004);
names[104]="S9_B0Kstar0mumu_1.1-6       ";

central[105]=-0.049;
errors[105]=sqrt(0.109*0.109+0.014*0.014);
names[105]="P5p_B0Kstar0mumu_1.1-6      ";

central[106]=0.;
errors[106]=0.;
names[106]="dG_B0Kstar0mumu_15-19       ";

if(likelihood==1)
{
central[107]=0.344;
errors[107]=sqrt(0.030*0.030+0.008*0.008);
names[107]="FL_B0Kstar0mumu_15-19       ";

central[108]=0.355;
errors[108]=sqrt(0.028*0.028+0.009*0.009);
names[108]="AFB_B0Kstar0mumu_15-19      ";

central[109]=-0.163;
errors[109]=sqrt(0.035*0.035+0.009*0.009);
names[109]="S3_B0Kstar0mumu_15-19       ";

central[110]=0.284;
errors[110]=sqrt(0.040*0.040+0.007*0.007);
names[110]="S4_B0Kstar0mumu_15-19       ";

central[111]=-0.325;
errors[111]=sqrt(0.038*0.038+0.009*0.009);
names[111]="S5_B0Kstar0mumu_15-19       ";

central[112]=0.048;
errors[112]=sqrt(0.044*0.044+0.006*0.006);
names[112]="S7_B0Kstar0mumu_15-19       ";

central[113]=0.028;
errors[113]=sqrt(0.045*0.045+0.003*0.003);
names[113]="S8_B0Kstar0mumu_15-19       ";

central[114]=-0.053;
errors[114]=sqrt(0.040*0.040+0.002*0.002);
names[114]="S9_B0Kstar0mumu_15-19       ";

central[115]=-0.684;
errors[115]=sqrt(0.083*0.083+0.020*0.020);
names[115]="P5p_B0Kstar0mumu_15-19      ";
}
else
{
central[107]=0.357;
errors[107]=sqrt(0.035*0.035+0.011*0.011);
names[107]="FL_B0Kstar0mumu_15-19       ";

central[108]=0.367;
errors[108]=sqrt(0.037*0.037+0.007*0.007);
names[108]="AFB_B0Kstar0mumu_15-19      ";

central[109]=-0.135;
errors[109]=sqrt(0.050*0.050+0.011*0.011);
names[109]="S3_B0Kstar0mumu_15-19       ";

central[110]=0.314;
errors[110]=sqrt(0.054*0.054+0.027*0.027);
names[110]="S4_B0Kstar0mumu_15-19       ";

central[111]=-0.335;
errors[111]=sqrt(0.047*0.047+0.007*0.007);
names[111]="S5_B0Kstar0mumu_15-19       ";

central[112]=0.066;
errors[112]=sqrt(0.049*0.049+0.014*0.014);
names[112]="S7_B0Kstar0mumu_15-19       ";

central[113]=0.024;
errors[113]=sqrt(0.048*0.048+0.009*0.009);
names[113]="S8_B0Kstar0mumu_15-19       ";

central[114]=-0.056;
errors[114]=sqrt(0.047*0.047+0.014*0.014);
names[114]="S9_B0Kstar0mumu_15-19       ";

central[115]=-0.709;
errors[115]=sqrt(0.093*0.093+0.016*0.016);
names[115]="P5p_B0Kstar0mumu_15-19      ";
}

central[116]=0.;
errors[116]=0.;
names[116]="dG_B0Kstar0mumu_1.1-2       ";

central[117]=0.768;
errors[117]=sqrt(0.141*0.141+0.025*0.025);
names[117]="FL_B0Kstar0mumu_1.1-2       ";

central[118]=-0.333;
errors[118]=sqrt(0.130*0.130+0.012*0.012);
names[118]="AFB_B0Kstar0mumu_1.1-2      ";

central[119]=0.065;
errors[119]=sqrt(0.137*0.137+0.007*0.007);
names[119]="S3_B0Kstar0mumu_1.1-2       ";

central[120]=-0.127;
errors[120]=sqrt(0.190*0.190+0.027*0.027);
names[120]="S4_B0Kstar0mumu_1.1-2       ";

central[121]=0.286;
errors[121]=sqrt(0.172*0.172+0.009*0.009);
names[121]="S5_B0Kstar0mumu_1.1-2       ";

central[122]=-0.293;
errors[122]=sqrt(0.180*0.180+0.005*0.005);
names[122]="S7_B0Kstar0mumu_1.1-2       ";

central[123]=-0.114;
errors[123]=sqrt(0.196*0.196+0.006*0.006);
names[123]="S8_B0Kstar0mumu_1.1-2       ";

central[124]=-0.110;
errors[124]=sqrt(0.140*0.140+0.001*0.001);
names[124]="S9_B0Kstar0mumu_1.1-2       ";

central[125]=0.606;
errors[125]=sqrt(0.769*0.769+0.017*0.017);
names[125]="P5p_B0Kstar0mumu_1.1-2      ";

central[126]=0.;
errors[126]=0.;
names[126]="dG_B0Kstar0mumu_2-3         ";

central[127]=0.690;
errors[127]=sqrt(0.113*0.113+0.023*0.023);
names[127]="FL_B0Kstar0mumu_2-3         ";

central[128]=-0.158;
errors[128]=sqrt(0.090*0.090+0.008*0.008);
names[128]="AFB_B0Kstar0mumu_2-3        ";

central[129]=0.006;
errors[129]=sqrt(0.100*0.100+0.007*0.007);
names[129]="S3_B0Kstar0mumu_2-3         ";

central[130]=0.339;
errors[130]=sqrt(0.140*0.140+0.041*0.041);
names[130]="S4_B0Kstar0mumu_2-3         ";

central[131]=0.206;
errors[131]=sqrt(0.131*0.131+0.009*0.009);
names[131]="S5_B0Kstar0mumu_2-3         ";

central[132]=-0.252;
errors[132]=sqrt(0.151*0.151+0.002*0.002);
names[132]="S7_B0Kstar0mumu_2-3         ";

central[133]=-0.176;
errors[133]=sqrt(0.165*0.165+0.006*0.006);
names[133]="S8_B0Kstar0mumu_2-3         ";

central[134]=0.000;
errors[134]=sqrt(0.102*0.102+0.003*0.003);
names[134]="S9_B0Kstar0mumu_2-3         ";

central[135]=0.461;
errors[135]=sqrt(0.313*0.313+0.019*0.019);
names[135]="P5p_B0Kstar0mumu_2-3        ";

central[136]=0.;
errors[136]=0.;
names[136]="dG_B0Kstar0mumu_3-4         ";

central[137]=0.873;
errors[137]=sqrt(0.154*0.154+0.023*0.023);
names[137]="FL_B0Kstar0mumu_3-4         ";

central[138]=-0.041;
errors[138]=sqrt(0.091*0.091+0.002*0.002);
names[138]="AFB_B0Kstar0mumu_3-4        ";

central[139]=0.078;
errors[139]=sqrt(0.131*0.131+0.008*0.008);
names[139]="S3_B0Kstar0mumu_3-4         ";

central[140]=0.046;
errors[140]=sqrt(0.196*0.196+0.046*0.046);
names[140]="S4_B0Kstar0mumu_3-4         ";

central[141]=-0.110;
errors[141]=sqrt(0.169*0.169+0.004*0.004);
names[141]="S5_B0Kstar0mumu_3-4         ";

central[142]=0.171;
errors[142]=sqrt(0.175*0.175+0.002*0.002);
names[142]="S7_B0Kstar0mumu_3-4         ";

central[143]=0.097;
errors[143]=sqrt(0.189*0.189+0.002*0.002);
names[143]="S8_B0Kstar0mumu_3-4         ";

central[144]=-0.203;
errors[144]=sqrt(0.132*0.132+0.002*0.002);
names[144]="S9_B0Kstar0mumu_3-4         ";

central[145]=-0.295;
errors[145]=sqrt(7.112*7.112+0.023*0.023);
names[145]="P5p_B0Kstar0mumu_3-4        ";

central[146]=0.;
errors[146]=0.;
names[146]="dG_B0Kstar0mumu_4-5         ";

central[147]=0.899;
errors[147]=sqrt(0.106*0.106+0.023*0.023);
names[147]="FL_B0Kstar0mumu_4-5         ";

central[148]=0.052;
errors[148]=sqrt(0.080*0.080+0.004*0.004);
names[148]="AFB_B0Kstar0mumu_4-5        ";

central[149]=0.200;
errors[149]=sqrt(0.101*0.101+0.007*0.007);
names[149]="S3_B0Kstar0mumu_4-5         ";

central[150]=0.148;
errors[150]=sqrt(0.154*0.154+0.047*0.047);
names[150]="S4_B0Kstar0mumu_4-5         ";

central[151]=-0.306;
errors[151]=sqrt(0.141*0.141+0.004*0.004);
names[151]="S5_B0Kstar0mumu_4-5         ";

central[152]=-0.082;
errors[152]=sqrt(0.129*0.129+0.001*0.001);
names[152]="S7_B0Kstar0mumu_4-5         ";

central[153]=0.107;
errors[153]=sqrt(0.146*0.146+0.003*0.003);
names[153]="S8_B0Kstar0mumu_4-5         ";

central[154]=0.181;
errors[154]=sqrt(0.105*0.105+0.001*0.001);
names[154]="S9_B0Kstar0mumu_4-5         ";

central[155]=-0.799;
errors[155]=sqrt(18.19*18.19+0.022*0.022);
names[155]="P5p_B0Kstar0mumu_4-5        ";

central[156]=0.;
errors[156]=0.;
names[156]="dG_B0Kstar0mumu_5-6         ";

central[157]=0.644;
errors[157]=sqrt(0.130*0.130+0.025*0.025);
names[157]="FL_B0Kstar0mumu_5-6         ";

central[158]=0.057;
errors[158]=sqrt(0.094*0.094+0.006*0.006);
names[158]="AFB_B0Kstar0mumu_5-6        ";

central[159]=-0.122;
errors[159]=sqrt(0.126*0.126+0.009*0.009);
names[159]="S3_B0Kstar0mumu_5-6         ";

central[160]=0.273;
errors[160]=sqrt(0.184*0.184+0.048*0.048);
names[160]="S4_B0Kstar0mumu_5-6         ";

central[161]=-0.095;
errors[161]=sqrt(0.142*0.142+0.004*0.004);
names[161]="S5_B0Kstar0mumu_5-6         ";

central[162]=0.038;
errors[162]=sqrt(0.135*0.135+0.002*0.002);
names[162]="S7_B0Kstar0mumu_5-6         ";

central[163]=-0.037;
errors[163]=sqrt(0.160*0.160+0.003*0.003);
names[163]="S8_B0Kstar0mumu_5-6         ";

central[164]=-0.080;
errors[164]=sqrt(0.120*0.120+0.001*0.001);
names[164]="S9_B0Kstar0mumu_5-6         ";

central[165]=-0.197;
errors[165]=sqrt(0.334*0.334+0.018*0.018);
names[165]="P5p_B0Kstar0mumu_5-6        ";

central[166]=0.;
errors[166]=0.;
names[166]="dG_B0Kstar0mumu_6-7         ";

central[167]=0.644;
errors[167]=sqrt(0.089*0.089+0.025*0.025);
names[167]="FL_B0Kstar0mumu_6-7         ";

central[168]=0.058;
errors[168]=sqrt(0.064*0.064+0.009*0.009);
names[168]="AFB_B0Kstar0mumu_6-7        ";

central[169]=-0.069;
errors[169]=sqrt(0.091*0.091+0.004*0.004);
names[169]="S3_B0Kstar0mumu_6-7         ";

central[170]=0.311;
errors[170]=sqrt(0.118*0.118+0.052*0.052);
names[170]="S4_B0Kstar0mumu_6-7         ";

central[171]=-0.339;
errors[171]=sqrt(0.114*0.114+0.008*0.008);
names[171]="S5_B0Kstar0mumu_6-7         ";

central[172]=0.009;
errors[172]=sqrt(0.124*0.124+0.004*0.004);
names[172]="S7_B0Kstar0mumu_6-7         ";

central[173]=0.080;
errors[173]=sqrt(0.131*0.131+0.002*0.002);
names[173]="S8_B0Kstar0mumu_6-7         ";

central[174]=0.061;
errors[174]=sqrt(0.091*0.091+0.001*0.001);
names[174]="S9_B0Kstar0mumu_6-7         ";

central[175]=-0.713;
errors[175]=sqrt(0.268*0.268+0.015*0.015);
names[175]="P5p_B0Kstar0mumu_6-7        ";

central[176]=0.;
errors[176]=0.;
names[176]="dG_B0Kstar0mumu_7-8         ";

central[177]=0.609;
errors[177]=sqrt(0.103*0.103+0.025*0.025);
names[177]="FL_B0Kstar0mumu_7-8         ";

central[178]=0.241;
errors[178]=sqrt(0.080*0.080+0.0012*0.012);
names[178]="AFB_B0Kstar0mumu_7-8        ";

central[179]=-0.054;
errors[179]=sqrt(0.099*0.099+0.005*0.005);
names[179]="S3_B0Kstar0mumu_7-8         ";

central[180]=0.236;
errors[180]=sqrt(0.136*0.136+0.058*0.058);
names[180]="S4_B0Kstar0mumu_7-8         ";

central[181]=-0.386;
errors[181]=sqrt(0.135*0.135+0.007*0.007);
names[181]="S5_B0Kstar0mumu_7-8         ";

central[182]=-0.094;
errors[182]=sqrt(0.130*0.130+0.003*0.003);
names[182]="S7_B0Kstar0mumu_7-8         ";

central[183]=-0.295;
errors[183]=sqrt(0.139*0.139+0.002*0.002);
names[183]="S8_B0Kstar0mumu_7-8         ";

central[184]=0.030;
errors[184]=sqrt(0.100*0.100+0.001*0.001);
names[184]="S9_B0Kstar0mumu_7-8         ";

central[185]=-0.808;
errors[185]=sqrt(0.303*0.303+0.010*0.010);
names[185]="P5p_B0Kstar0mumu_7-8        ";

central[186]=0.;
errors[186]=0.;
names[186]="dG_B0Kstar0mumu_11-11.75    ";

central[187]=0.502;
errors[187]=sqrt(0.090*0.090+0.022*0.022);
names[187]="FL_B0Kstar0mumu_11-11.75    ";

central[188]=0.370;
errors[188]=sqrt(0.076*0.076+0.015*0.015);
names[188]="AFB_B0Kstar0mumu_11-11.75   ";

central[189]=-0.217;
errors[189]=sqrt(0.090*0.090+0.008*0.008);
names[189]="S3_B0Kstar0mumu_11-11.75    ";

central[190]=0.252;
errors[190]=sqrt(0.113*0.113+0.063*0.063);
names[190]="S4_B0Kstar0mumu_11-11.75    ";

central[191]=-0.235;
errors[191]=sqrt(0.115*0.115+0.013*0.013);
names[191]="S5_B0Kstar0mumu_11-11.75    ";

central[192]=-0.110;
errors[192]=sqrt(0.114*0.114+0.002*0.002);
names[192]="S7_B0Kstar0mumu_11-11.75    ";

central[193]=-0.079;
errors[193]=sqrt(0.122*0.122+0.003*0.003);
names[193]="S8_B0Kstar0mumu_11-11.75    ";

central[194]=-0.084;
errors[194]=sqrt(0.102*0.102+0.003*0.003);
names[194]="S9_B0Kstar0mumu_11-11.75    ";

central[195]=-0.485;
errors[195]=sqrt(0.224*0.224+0.028*0.028);
names[195]="P5p_B0Kstar0mumu_11-11.75   ";

central[196]=0.;
errors[196]=0.;
names[196]="dG_B0Kstar0mumu_11.75-12.5  ";

central[197]=0.734;
errors[197]=sqrt(0.107*0.107+0.018*0.018);
names[197]="FL_B0Kstar0mumu_11.75-12.5  ";

central[198]=0.293;
errors[198]=sqrt(0.064*0.064+0.014*0.014);
names[198]="AFB_B0Kstar0mumu_11.75-12.5 ";

central[199]=-0.157;
errors[199]=sqrt(0.098*0.098+0.008*0.008);
names[199]="S3_B0Kstar0mumu_11.75-12.5  ";

central[200]=0.309;
errors[200]=sqrt(0.111*0.111+0.056*0.056);
names[200]="S4_B0Kstar0mumu_11.75-12.5  ";

central[201]=-0.366;
errors[201]=sqrt(0.112*0.112+0.012*0.012);
names[201]="S5_B0Kstar0mumu_11.75-12.5  ";

central[202]=-0.212;
errors[202]=sqrt(0.118*0.118+0.002*0.002);
names[202]="S7_B0Kstar0mumu_11.75-12.5  ";

central[203]=-0.090;
errors[203]=sqrt(0.111*0.111+0.003*0.003);
names[203]="S8_B0Kstar0mumu_11.75-12.5  ";

central[204]=0.030;
errors[204]=sqrt(0.093*0.093+0.002*0.002);
names[204]="S9_B0Kstar0mumu_11.75-12.5  ";

central[205]=-0.827;
errors[205]=sqrt(0.357*0.357+0.026*0.026);
names[205]="P5p_B0Kstar0mumu_11.75-12.5 ";

central[206]=0.;
errors[206]=0.;
names[206]="dG_B0Kstar0mumu_15-16       ";

central[207]=0.385;
errors[207]=sqrt(0.067*0.067+0.013*0.013);
names[207]="FL_B0Kstar0mumu_15-16       ";

central[208]=0.396;
errors[208]=sqrt(0.068*0.068+0.009*0.009);
names[208]="AFB_B0Kstar0mumu_15-16      ";

central[209]=-0.060;
errors[209]=sqrt(0.085*0.085+0.006*0.006);
names[209]="S3_B0Kstar0mumu_15-16       ";

central[210]=0.321;
errors[210]=sqrt(0.099*0.099+0.007*0.007);
names[210]="S4_B0Kstar0mumu_15-16       ";

central[211]=-0.360;
errors[211]=sqrt(0.092*0.092+0.006*0.006);
names[211]="S5_B0Kstar0mumu_15-16       ";

central[212]=0.040;
errors[212]=sqrt(0.092*0.092+0.002*0.002);
names[212]="S7_B0Kstar0mumu_15-16       ";

central[213]=-0.057;
errors[213]=sqrt(0.095*0.095+0.005*0.005);
names[213]="S8_B0Kstar0mumu_15-16       ";

central[214]=-0.054;
errors[214]=sqrt(0.087*0.087+0.005*0.005);
names[214]="S9_B0Kstar0mumu_15-16       ";

central[215]=-0.758;
errors[215]=sqrt(0.179*0.179+0.013*0.013);
names[215]="P5p_B0Kstar0mumu_15-16      ";

central[216]=0.;
errors[216]=0.;
names[216]="dG_B0Kstar0mumu_16-17       ";

central[217]=0.295;
errors[217]=sqrt(0.062*0.062+0.013*0.013);
names[217]="FL_B0Kstar0mumu_16-17       ";

central[218]=0.451;
errors[218]=sqrt(0.071*0.071+0.007*0.007);
names[218]="AFB_B0Kstar0mumu_16-17      ";

central[219]=-0.250;
errors[219]=sqrt(0.092*0.092+0.007*0.007);
names[219]="S3_B0Kstar0mumu_16-17       ";

central[220]=0.246;
errors[220]=sqrt(0.096*0.096+0.029*0.029);
names[220]="S4_B0Kstar0mumu_16-17       ";

central[221]=-0.254;
errors[221]=sqrt(0.081*0.081+0.010*0.010);
names[221]="S5_B0Kstar0mumu_16-17       ";

central[222]=0.144;
errors[222]=sqrt(0.091*0.091+0.005*0.005);
names[222]="S7_B0Kstar0mumu_16-17       ";

central[223]=0.055;
errors[223]=sqrt(0.090*0.090+0.005*0.005);
names[223]="S8_B0Kstar0mumu_16-17       ";

central[224]=-0.014;
errors[224]=sqrt(0.086*0.086+0.004*0.004);
names[224]="S9_B0Kstar0mumu_16-17       ";

central[225]=-0.567;
errors[225]=sqrt(0.186*0.186+0.014*0.014);
names[225]="P5p_B0Kstar0mumu_16-17      ";

central[226]=0.;
errors[226]=0.;
names[226]="dG_B0Kstar0mumu_17-18       ";

central[227]=0.363;
errors[227]=sqrt(0.073*0.073+0.017*0.017);
names[227]="FL_B0Kstar0mumu_17-18       ";

central[228]=0.274;
errors[228]=sqrt(0.069*0.069+0.008*0.008);
names[228]="AFB_B0Kstar0mumu_17-18      ";

central[229]=-0.099;
errors[229]=sqrt(0.092*0.092+0.011*0.011);
names[229]="S3_B0Kstar0mumu_17-18       ";

central[230]=0.229;
errors[230]=sqrt(0.096*0.096+0.045*0.045);
names[230]="S4_B0Kstar0mumu_17-18       ";

central[231]=-0.305;
errors[231]=sqrt(0.088*0.088+0.015*0.015);
names[231]="S5_B0Kstar0mumu_17-18       ";

central[232]=0.022;
errors[232]=sqrt(0.094*0.094+0.011*0.011);
names[232]="S7_B0Kstar0mumu_17-18       ";

central[233]=-0.007;
errors[233]=sqrt(0.098*0.098+0.001*0.001);
names[233]="S8_B0Kstar0mumu_17-18       ";

central[234]=-0.090;
errors[234]=sqrt(0.095*0.095+0.002*0.002);
names[234]="S9_B0Kstar0mumu_17-18       ";

central[235]=-0.646;
errors[235]=sqrt(0.176*0.176+0.027*0.027);
names[235]="P5p_B0Kstar0mumu_17-18      ";

central[236]=0.;
errors[236]=0.;
names[236]="dG_B0Kstar0mumu_18-19       ";

central[237]=0.421;
errors[237]=sqrt(0.100*0.100+0.013*0.013);
names[237]="FL_B0Kstar0mumu_18-19       ";

central[238]=0.354;
errors[238]=sqrt(0.111*0.111+0.012*0.012);
names[238]="AFB_B0Kstar0mumu_18-19      ";

central[239]=-0.131;
errors[239]=sqrt(0.130*0.130+0.012*0.012);
names[239]="S3_B0Kstar0mumu_18-19       ";

central[240]=0.607;
errors[240]=sqrt(0.170*0.170+0.059*0.059);
names[240]="S4_B0Kstar0mumu_18-19       ";

central[241]=-0.534;
errors[241]=sqrt(0.150*0.150+0.015*0.015);
names[241]="S5_B0Kstar0mumu_18-19       ";

central[242]=0.058;
errors[242]=sqrt(0.124*0.124+0.006*0.006);
names[242]="S7_B0Kstar0mumu_18-19       ";

central[243]=0.149;
errors[243]=sqrt(0.139*0.139+0.010*0.010);
names[243]="S8_B0Kstar0mumu_18-19       ";

central[244]=-0.079;
errors[244]=sqrt(0.122*0.122+0.007*0.007);
names[244]="S9_B0Kstar0mumu_18-19       ";

central[245]=-1.070;
errors[245]=sqrt(0.349*0.349+0.029*0.029);
names[245]="P5p_B0Kstar0mumu_18-19      ";

central[246]=59.2e-9;
errors[246]=sqrt(14.4*14.4+4.0*4.0)*1.e-9;
names[246]="dG_B+Kstar+mumu_0.1-2       ";

central[247]=55.9e-9;
errors[247]=sqrt(15.9*15.9+3.8*3.8)*1.e-9;
names[247]="dG_B+Kstar+mumu_2-4         ";

central[248]=24.9e-9;
errors[248]=sqrt(11.0*11.0+1.7*1.7)*1.e-9;
names[248]="dG_B+Kstar+mumu_4-6         ";

central[249]=33.0e-9;
errors[249]=sqrt(11.3*11.3+2.3*2.3)*1.e-9;
names[249]="dG_B+Kstar+mumu_6-8         ";

central[250]=82.8e-9;
errors[250]=sqrt(15.8*15.8+5.6*5.6)*1.e-9;
names[250]="dG_B+Kstar+mumu_11-12.5     ";

central[251]=64.4e-9;
errors[251]=sqrt(12.9*12.9+4.4*4.4)*1.e-9;
names[251]="dG_B+Kstar+mumu_15-17       ";

central[252]=11.6e-9;
errors[252]=sqrt(9.1*9.1+0.8*0.8)*1.e-9;
names[252]="dG_B+Kstar+mumu_17-19       ";

central[253]=36.6e-9;
errors[253]=8.7e-9;
names[253]="dG_B+Kstar+mumu_1.1-6       ";

central[254]=39.5e-9;
errors[254]=8.5e-9;
names[254]="dG_B+Kstar+mumu_15-19       ";

central[255]=1.03e-6;
errors[255]=0.19e-6;
names[255]="BR_B0Kstar0ee_full          ";

central[256]=3.1e-7;
errors[256]=0.9e-7;
names[256]="dG_B0Kstar0ee_0.0009-1      ";

central[257]=0.;
errors[257]=0.;
names[257]="dG_B0Kstar0ee_0.002-1.12    ";

central[258]=0.16;
errors[258]=sqrt(0.06*0.06+0.03*0.03);
names[258]="FL_B0Kstar0ee_0.002_1.12    ";

central[259]=-0.23;
errors[259]=sqrt(0.23*0.23+0.05*0.05);
names[259]="AT2_B0Kstar0ee_0.002_1.12   ";

central[260]=0.10;
errors[260]=sqrt(0.18*0.18+0.05*0.05);
names[260]="ATRe_B0Kstar0ee_0.002_1.12  ";

central[261]=0.;
errors[261]=0.;
names[261]="dG_B0Kstar0ee_1.1-6         ";

central[262]=0.;
errors[262]=0.;
names[262]="FL_B0Kstar0ee_1.1_6         ";

central[263]=0.;
errors[263]=0.;
names[263]="P1_B0Kstar0ee_1.1_6         ";

central[264]=0.;
errors[264]=0.;
names[264]="ATRe_B0Kstar0ee_1.1_6       ";

central[265]=12.2e-9;
errors[265]=sqrt(5.9*5.9+0.6*0.6)*1.e-9;
names[265]="dG_B0K0mumu_0.1-2           ";

central[266]=18.7e-9;
errors[266]=sqrt(5.5*5.5+0.9*0.9)*1.e-9;
names[266]="dG_B0K0mumu_2-4             ";

central[267]=17.3e-9;
errors[267]=sqrt(5.3*5.3+0.9*0.9)*1.e-9;
names[267]="dG_B0K0mumu_4-6             ";

central[268]=27.0e-9;
errors[268]=sqrt(5.8*5.8+1.4*1.4)*1.e-9;
names[268]="dG_B0K0mumu_6-8             ";

central[269]=12.7e-9;
errors[269]=sqrt(4.5*4.5+0.6*0.6)*1.e-9;
names[269]="dG_B0K0mumu_11-12.5         ";

central[270]=14.3e-9;
errors[270]=sqrt(3.5*3.5+0.7*0.7)*1.e-9;
names[270]="dG_B0K0mumu_15-17           ";

central[271]=7.8e-9;
errors[271]=sqrt(1.7*1.7+0.4*0.4)*1.e-9;
names[271]="dG_B0K0mumu_17-22           ";

central[272]=18.7e-9;
errors[272]=3.6e-9;
names[272]="dG_B0K0mumu_1.1-6           ";

central[273]=9.5e-9;
errors[273]=1.7e-9;
names[273]="dG_B0K0mumu_15-22           ";

central[274]=33.2e-9;
errors[274]=sqrt(1.8*1.8+1.7*1.7)*1.e-9;
names[274]="dG_B+K+mumu_0.1-0.98         ";

central[275]=23.3e-9;
errors[275]=sqrt(1.5*1.5+1.2*1.2)*1.e-9;
names[275]="dG_B+K+mumu_1.1-2           ";

central[276]=28.2e-9;
errors[276]=sqrt(1.6*1.6+1.4*1.4)*1.e-9;
names[276]="dG_B+K+mumu_2-3             ";

central[277]=25.4e-9;
errors[277]=sqrt(1.5*1.5+1.3*1.3)*1.e-9;
names[277]="dG_B+K+mumu_3-4             ";

central[278]=22.1e-9;
errors[278]=sqrt(1.4*1.4+1.1*1.1)*1.e-9;
names[278]="dG_B+K+mumu_4-5             ";

central[279]=23.1e-9;
errors[279]=sqrt(1.4*1.4+1.2*1.2)*1.e-9;
names[279]="dG_B+K+mumu_5-6             ";

central[280]=24.5e-9;
errors[280]=sqrt(1.4*1.4+1.2*1.2)*1.e-9;
names[280]="dG_B+K+mumu_6-7             ";

central[281]=23.1e-9;
errors[281]=sqrt(1.4*1.4+1.2*1.2)*1.e-9;
names[281]="dG_B+K+mumu_7-8             ";

central[282]=17.7e-9;
errors[282]=sqrt(1.3*1.3+0.9*0.9)*1.e-9;
names[282]="dG_B+K+mumu_11-11.8         ";

central[283]=19.3e-9;
errors[283]=sqrt(1.2*1.2+1.0*1.0)*1.e-9;
names[283]="dG_B+K+mumu_11.8-12.5       ";

central[284]=16.1e-9;
errors[284]=sqrt(1.0*1.0+0.8*0.8)*1.e-9;
names[284]="dG_B+K+mumu_15-16           ";

central[285]=16.4e-9;
errors[285]=sqrt(1.0*1.0+0.8*0.8)*1.e-9;
names[285]="dG_B+K+mumu_16-17           ";

central[286]=20.6e-9;
errors[286]=sqrt(1.1*1.1+1.0*1.0)*1.e-9;
names[286]="dG_B+K+mumu_17-18           ";

central[287]=13.7e-9;
errors[287]=sqrt(1.0*1.0+0.7*0.7)*1.e-9;
names[287]="dG_B+K+mumu_18-19           ";

central[288]=7.4e-9;
errors[288]=sqrt(0.8*0.8+0.4*0.4)*1.e-9;
names[288]="dG_B+K+mumu_19-20           ";

central[289]=5.9e-9;
errors[289]=sqrt(0.7*0.7+0.3*0.3)*1.e-9;
names[289]="dG_B+K+mumu_20-21           ";

central[290]=4.3e-9;
errors[290]=sqrt(0.7*0.7+0.2*0.2)*1.e-9;
names[290]="dG_B+K+mumu_21-22           ";

central[291]=24.2e-9;
errors[291]=1.4e-9;
names[291]="dG_B+K+mumu_1.1-6           ";

central[292]=12.1e-9;
errors[292]=0.7e-9;
names[292]="dG_B+K+mumu_15-22           ";

central[293]=-0.;
errors[293]=-0.;
names[293]="dG_B+K+mumu_1-6             ";

central[294]=1.56e-7;
errors[294]=sqrt(0.19*0.19+0.06*0.06);
names[294]="dG_B+K+ee_1-6               ";

central[295]=0.745-1.;
errors[295]=sqrt(0.09*0.09+0.036*0.036);
names[295]="RK-1                        ";

central[296]=5.85e-8;
errors[296]=sqrt(0.73*0.73+0.14*0.14+0.44*0.44)*1.e-8;
names[296]="dG_Bsphimumu_0.1-2          ";

central[297]=0.20;
errors[297]=sqrt(0.09*0.09+0.02*0.02);
names[297]="FL_Bsphimumu_0.1-2          ";

central[298]=-0.05;
errors[298]=sqrt(0.13*0.13+0.01*0.01);
names[298]="S3_Bsphimumu_0.1-2          ";

central[299]=-0.27;
errors[299]=sqrt(0.28*0.28+0.01*0.01);
names[299]="S4_Bsphimumu_0.1-2          ";

central[300]=0.04;
errors[300]=0.12;
names[300]="S7_Bsphimumu_0.1-2          ";

central[301]=2.56e-8;
errors[301]=sqrt(0.42*0.42+0.06*0.06+0.19*0.19)*1.e-8;
names[301]="dG_Bsphimumu_2-5            ";

central[302]=0.68;
errors[302]=sqrt(0.16*0.16+0.03*0.03);
names[302]="FL_Bsphimumu_2-5            ";

central[303]=-0.06;
errors[303]=sqrt(0.23*0.23+0.01*0.01);
names[303]="S3_Bsphimumu_2-5            ";

central[304]=0.47;
errors[304]=sqrt(0.44*0.44+0.01*0.01);
names[304]="S4_Bsphimumu_2-5            ";

central[305]=-0.03;
errors[305]=sqrt(0.23*0.23+0.01*0.01);
names[305]="S7_Bsphimumu_2-5            ";

central[306]=3.21e-8;
errors[306]=sqrt(0.44*0.44+0.08*0.08+0.24*0.24)*1.e-8;
names[306]="dG_Bsphimumu_5-8            ";

central[307]=0.54;
errors[307]=sqrt(0.10*0.10+0.02*0.02);
names[307]="FL_Bsphimumu_5-8            ";

central[308]=-0.10;
errors[308]=sqrt(0.29*0.29+0.01*0.01);
names[308]="S3_Bsphimumu_5-8            ";

central[309]=0.10;
errors[309]=sqrt(0.18*0.18+0.01*0.01);
names[309]="S4_Bsphimumu_5-8            ";

central[310]=0.04;
errors[310]=sqrt(0.20*0.20+0.01*0.01);
names[310]="S7_Bsphimumu_5-8            ";

central[311]=4.71e-8;
errors[311]=sqrt(0.69*0.69+0.15*0.15+0.36*0.36)*1.e-8;
names[311]="dG_Bsphimumu_11-12.5        ";

central[312]=0.29;
errors[312]=sqrt(0.11*0.11+0.04*0.04);
names[312]="FL_Bsphimumu_11-12.5        ";

central[313]=-0.19;
errors[313]=sqrt(0.23*0.23+0.01*0.01);
names[313]="S3_Bsphimumu_11-12.5        ";

central[314]=0.47;
errors[314]=sqrt(0.29*0.29+0.01*0.01);
names[314]="S4_Bsphimumu_11-12.5        ";

central[315]=0.;
errors[315]=sqrt(0.17*0.17+0.01*0.01);
names[315]="S7_Bsphimumu_11-12.5        ";

central[316]=4.52e-8;
errors[316]=sqrt(0.57*0.57+0.12*0.12+0.34*0.34)*1.e-8;
names[316]="dG_Bsphimumu_15-17          ";

central[317]=0.23;
errors[317]=sqrt(0.09*0.09+0.02*0.02);
names[317]="FL_Bsphimumu_15-17          ";

central[318]=-0.06;
errors[318]=sqrt(0.19*0.19+0.01*0.01);
names[318]="S3_Bsphimumu_15-17          ";

central[319]=0.03;
errors[319]=sqrt(0.15*0.15+0.01*0.01);
names[319]="S4_Bsphimumu_15-17          ";

central[320]=0.12;
errors[320]=sqrt(0.16*0.16+0.01*0.01);
names[320]="S7_Bsphimumu_15-17          ";

central[321]=3.96e-8;
errors[321]=sqrt(0.57*0.57+0.14*0.14+0.30*0.30)*1.e-8;
names[321]="dG_Bsphimumu_17-19          ";

central[322]=0.40;
errors[322]=sqrt(0.15*0.15+0.02*0.02);
names[322]="FL_Bsphimumu_17-19          ";

central[323]=-0.07;
errors[323]=sqrt(0.27*0.27+0.02*0.02);
names[323]="S3_Bsphimumu_17-19          ";

central[324]=0.39;
errors[324]=sqrt(0.34*0.34+0.02*0.02);
names[324]="S4_Bsphimumu_17-19          ";

central[325]=0.20;
errors[325]=sqrt(0.29*0.29+0.01*0.01);
names[325]="S7_Bsphimumu_17-19          ";

central[326]=2.58e-8;
errors[326]=sqrt(0.33*0.33+0.08*0.08+0.19*0.19)*1.e-8;
names[326]="dG_Bsphimumu_1-6            ";

central[327]=0.63;
errors[327]=sqrt(0.09*0.09+0.03*0.03);
names[327]="FL_Bsphimumu_1-6            ";

central[328]=-0.02;
errors[328]=sqrt(0.13*0.13+0.01*0.01);
names[328]="S3_Bsphimumu_1-6            ";

central[329]=0.19;
errors[329]=sqrt(0.14*0.14+0.01*0.01);
names[329]="S4_Bsphimumu_1-6            ";

central[330]=-0.03;
errors[330]=0.14;
names[330]="S7_Bsphimumu_1-6            ";

central[331]=4.04e-8;
errors[331]=sqrt(0.39*0.39+0.13*0.13+0.30*0.30)*1.e-8;
names[331]="dG_Bsphimumu_15-19          ";

central[332]=0.29;
errors[332]=sqrt(0.07*0.07+0.02*0.02);
names[332]="FL_Bsphimumu_15-19          ";

central[333]=-0.09;
errors[333]=sqrt(0.12*0.12+0.01*0.01);
names[333]="S3_Bsphimumu_15-19          ";

central[334]=0.14;
errors[334]=sqrt(0.11*0.11+0.01*0.01);
names[334]="S4_Bsphimumu_15-19          ";

central[335]=0.13;
errors[335]=sqrt(0.11*0.11+0.01*0.01);
names[335]="S7_Bsphimumu_15-19          ";

central[336]=0.080;
errors[336]=sqrt(0.248*0.248+0.044*0.044);
names[336]="P1_B0Kstar0mumu_1.1_6       ";

central[337]=-0.323;
errors[337]=0.159;
names[337]="ATRe_B0Kstar0mumu_1.1_6     ";

if(likelihood==1)
{
central[338]=-0.099;
errors[338]=sqrt(0.168*0.168+0.014*0.014);
names[338]="P1_B0Kstar0mumu_0.1_0.98   ";
}
else
{
central[338]=-0.038;
errors[338]=sqrt(0.157*0.157+0.020*0.020);
names[338]="P1_B0Kstar0mumu_0.1_0.98   ";
}

central[339]=0.;
errors[339]=0.;
names[339]="ATRe_B0Kstar0mumu_0.1_0.98  ";

central[340]=0.;
errors[340]=0.;
names[340]="dG_B0Kstar0ee_1-6           ";

if(likelihood==1)
{
central[341]=-0.003;
errors[341]=sqrt(0.052*0.052+0.007*0.007);
names[341]="P2_B0Kstar0mumu_0.1-0.98    ";

central[342]=0.113;
errors[342]=sqrt(0.079*0.079+0.006*0.006);
names[342]="P3_B0Kstar0mumu_0.1-0.98    ";

central[343]=0.185;
errors[343]=sqrt(0.158*0.158+0.023*0.023);
names[343]="P4p_B0Kstar0mumu_0.1-0.98   ";

central[344]=0.034;
errors[344]=sqrt(0.135*0.135+0.015*0.015);
names[344]="P6p_B0Kstar0mumu_0.1-0.98   ";

central[345]=0.180;
errors[345]=sqrt(0.174*0.174+0.007*0.007);
names[345]="P8p_B0Kstar0mumu_0.1-0.98   ";
}
else
{
central[341]=-0.119;
errors[341]=sqrt(0.081*0.081+0.063*0.063);
names[341]="P2_B0Kstar0mumu_0.1-0.98    ";

central[342]=0.147;
errors[342]=sqrt(0.086*0.086+0.005*0.005);
names[342]="P3_B0Kstar0mumu_0.1-0.98    ";

central[343]=0.086;
errors[343]=sqrt(0.221*0.221+0.026*0.026);
names[343]="P4p_B0Kstar0mumu_0.1-0.98   ";

central[344]=0.086;
errors[344]=sqrt(0.152*0.152+0.024*0.024);
names[344]="P6p_B0Kstar0mumu_0.1-0.98   ";

central[345]=0.143;
errors[345]=sqrt(0.195*0.195+0.022*0.022);
names[345]="P8p_B0Kstar0mumu_0.1-0.98   ";
}

central[346]=-0.451;
errors[346]=sqrt(0.636*0.636+0.038*0.038);
names[346]="P1_B0Kstar0mumu_1.1-2.5     ";

central[347]=-0.373;
errors[347]=sqrt(0.199*0.199+0.027*0.027);
names[347]="P2_B0Kstar0mumu_1.1-2.5     ";

central[348]=0.350;
errors[348]=sqrt(0.330*0.330+0.015*0.015);
names[348]="P3_B0Kstar0mumu_1.1-2.5     ";

central[349]=-0.163;
errors[349]=sqrt(0.240*0.240+0.021*0.021);
names[349]="P4p_B0Kstar0mumu_1.1-2.5    ";

central[350]=-0.463;
errors[350]=sqrt(0.221*0.221+0.012*0.012);
names[350]="P6p_B0Kstar0mumu_1.1-2.5    ";

central[351]=-0.208;
errors[351]=sqrt(0.270*0.270+0.024*0.024);
names[351]="P8p_B0Kstar0mumu_1.1-2.5    ";

central[352]=0.571;
errors[352]=sqrt(2.404*2.404+0.045*0.045);
names[352]="P1_B0Kstar0mumu_2.5-4       ";

central[353]=-0.636;
errors[353]=sqrt(1.735*1.735+0.015*0.015);
names[353]="P2_B0Kstar0mumu_2.5-4       ";

central[354]=0.745;
errors[354]=sqrt(2.587*2.587+0.030*0.030);
names[354]="P3_B0Kstar0mumu_2.5-4       ";

central[355]=-0.713;
errors[355]=sqrt(1.305*1.305+0.024*0.024);
names[355]="P4p_B0Kstar0mumu_2.5-4      ";

central[356]=0.205;
errors[356]=sqrt(0.964*0.964+0.013*0.013);
names[356]="P6p_B0Kstar0mumu_2.5-4      ";

central[357]=0.091;
errors[357]=sqrt(0.650*0.650+0.025*0.025);
names[357]="P8p_B0Kstar0mumu_2.5-4      ";

central[358]=0.180;
errors[358]=sqrt(0.364*0.364+0.027*0.027);
names[358]="P1_B0Kstar0mumu_4-6         ";

central[359]=0.042;
errors[359]=sqrt(0.088*0.088+0.011*0.011);
names[359]="P2_B0Kstar0mumu_4-6         ";

central[360]=0.083;
errors[360]=sqrt(0.187*0.187+0.023*0.023);
names[360]="P3_B0Kstar0mumu_4-6         ";

central[361]=-0.448;
errors[361]=sqrt(0.172*0.172+0.020*0.020);
names[361]="P4p_B0Kstar0mumu_4-6        ";

central[362]=-0.032;
errors[362]=sqrt(0.167*0.167+0.007*0.007);
names[362]="P6p_B0Kstar0mumu_4-6        ";

central[363]=0.342;
errors[363]=sqrt(0.188*0.188+0.009*0.009);
names[363]="P8p_B0Kstar0mumu_4-6        ";

central[364]=-0.199;
errors[364]=sqrt(0.281*0.281+0.025*0.025);
names[364]="P1_B0Kstar0mumu_6-8         ";

central[365]=0.241;
errors[365]=sqrt(0.062*0.062+0.013*0.013);
names[365]="P2_B0Kstar0mumu_6-8         ";

central[366]=0.057;
errors[366]=sqrt(0.148*0.148+0.013*0.013);
names[366]="P3_B0Kstar0mumu_6-8         ";

central[367]=-0.599;
errors[367]=sqrt(0.135*0.135+0.010*0.010);
names[367]="P4p_B0Kstar0mumu_6-8        ";

central[368]=-0.095;
errors[368]=sqrt(0.135*0.135+0.011*0.011);
names[368]="P6p_B0Kstar0mumu_6-8        ";

central[369]=-0.171;
errors[369]=sqrt(0.143*0.143+0.006*0.006);
names[369]="P8p_B0Kstar0mumu_6-8        ";

central[370]=-0.745;
errors[370]=sqrt(0.230*0.230+0.015*0.015);
names[370]="P1_B0Kstar0mumu_11-12.5     ";

central[371]=0.418;
errors[371]=sqrt(0.053*0.053+0.005*0.005);
names[371]="P2_B0Kstar0mumu_11-12.5     ";

central[372]=0.007;
errors[372]=sqrt(0.141*0.141+0.010*0.010);
names[372]="P3_B0Kstar0mumu_11-12.5     ";

central[373]=-0.567;
errors[373]=sqrt(0.187*0.187+0.014*0.014);
names[373]="P4p_B0Kstar0mumu_11-12.5    ";

central[374]=-0.282;
errors[374]=sqrt(0.151*0.151+0.007*0.007);
names[374]="P6p_B0Kstar0mumu_11-12.5    ";

central[375]=-0.015;
errors[375]=sqrt(0.145*0.145+0.005*0.005);
names[375]="P8p_B0Kstar0mumu_11-12.5    ";

central[376]=-0.436;
errors[376]=sqrt(0.147*0.147+0.018*0.018);
names[376]="P1_B0Kstar0mumu_15-17       ";

central[377]=0.421;
errors[377]=sqrt(0.042*0.042+0.005*0.005);
names[377]="P2_B0Kstar0mumu_15-17       ";

central[378]=0.029;
errors[378]=sqrt(0.084*0.084+0.006*0.006);
names[378]="P3_B0Kstar0mumu_15-17       ";

central[379]=-0.672;
errors[379]=sqrt(0.151*0.151+0.016*0.016);
names[379]="P4p_B0Kstar0mumu_15-17      ";

central[380]=0.127;
errors[380]=sqrt(0.122*0.122+0.006*0.006);
names[380]="P6p_B0Kstar0mumu_15-17      ";

central[381]=0.007;
errors[381]=sqrt(0.129*0.129+0.005*0.005);
names[381]="P8p_B0Kstar0mumu_15-17      ";

central[382]=-0.581;
errors[382]=sqrt(0.263*0.263+0.037*0.037);
names[382]="P1_B0Kstar0mumu_17-19       ";

central[383]=0.314;
errors[383]=sqrt(0.048*0.048+0.007*0.007);
names[383]="P2_B0Kstar0mumu_17-19       ";

central[384]=0.145;
errors[384]=sqrt(0.107*0.107+0.008*0.008);
names[384]="P3_B0Kstar0mumu_17-19       ";

central[385]=-0.556;
errors[385]=sqrt(0.156*0.156+0.016*0.016);
names[385]="P4p_B0Kstar0mumu_17-19      ";

central[386]=0.092;
errors[386]=sqrt(0.152*0.152+0.025*0.025);
names[386]="P6p_B0Kstar0mumu_17-19      ";

central[387]=0.027;
errors[387]=sqrt(0.147*0.147+0.009*0.009);
names[387]="P8p_B0Kstar0mumu_17-19      ";

central[388]=-0.162;
errors[388]=sqrt(0.073*0.073+0.010*0.010);
names[388]="P2_B0Kstar0mumu_1.1_6       ";

central[389]=0.205;
errors[389]=sqrt(0.135*0.135+0.017*0.017);
names[389]="P3_B0Kstar0mumu_1.1_6       ";

central[390]=-0.336;
errors[390]=sqrt(0.124*0.124+0.012*0.012);
names[390]="P4p_B0Kstar0mumu_1.1_6      ";

central[391]=-0.166;
errors[391]=sqrt(0.108*0.108+0.021*0.021);
names[391]="P6p_B0Kstar0mumu_1.1_6      ";

central[392]=0.060;
errors[392]=sqrt(0.124*0.124+0.009*0.009);
names[392]="P8p_B0Kstar0mumu_1.1_6      ";

if(likelihood==1)
{
central[393]=-0.497;
errors[393]=sqrt(0.102*0.102+0.027*0.027);
names[393]="P1_B0Kstar0mumu_15-19       ";

central[394]=0.361;
errors[394]=sqrt(0.026*0.026+0.010*0.010);
names[394]="P2_B0Kstar0mumu_15-19       ";

central[395]=0.081;
errors[395]=sqrt(0.060*0.060+0.005*0.005);
names[395]="P3_B0Kstar0mumu_15-19       ";

central[396]=-0.597;
errors[396]=sqrt(0.085*0.085+0.015*0.015);
names[396]="P4p_B0Kstar0mumu_15-19      ";

central[397]=0.101;
errors[397]=sqrt(0.092*0.092+0.011*0.011);
names[397]="P6p_B0Kstar0mumu_15-19      ";

central[398]=0.059;
errors[398]=sqrt(0.094*0.094+0.008*0.008);
names[398]="P8p_B0Kstar0mumu_15-19      ";
}
else
{
central[393]=-0.424;
errors[393]=sqrt(0.150*0.150+0.028*0.028);
names[393]="P1_B0Kstar0mumu_15-19       ";

central[394]=0.385;
errors[394]=sqrt(0.036*0.036+0.010*0.010);
names[394]="P2_B0Kstar0mumu_15-19       ";

central[395]=0.089;
errors[395]=sqrt(0.072*0.072+0.019*0.019);
names[395]="P3_B0Kstar0mumu_15-19       ";

central[396]=-0.663;
errors[396]=sqrt(0.105*0.105+0.055*0.055);
names[396]="P4p_B0Kstar0mumu_15-19      ";

central[397]=0.140;
errors[397]=sqrt(0.101*0.101+0.032*0.032);
names[397]="P6p_B0Kstar0mumu_15-19      ";

central[398]=0.049;
errors[398]=sqrt(0.106*0.106+0.021*0.021);
names[398]="P8p_B0Kstar0mumu_15-19      ";
}

central[399]=0.439;
errors[399]=sqrt(1.916*1.916+0.012*0.012);
names[399]="P1_B0Kstar0mumu_1.1-2       ";

central[400]=-0.667;
errors[400]=sqrt(1.939*1.939+0.017*0.017);
names[400]="P2_B0Kstar0mumu_1.1-2       ";

central[401]=0.363;
errors[401]=sqrt(1.088*1.088+0.001*0.001);
names[401]="P3_B0Kstar0mumu_1.1-2       ";

central[402]=0.266;
errors[402]=sqrt(0.648*0.648+0.057*0.057);
names[402]="P4p_B0Kstar0mumu_1.1-2      ";

central[403]=-0.632;
errors[403]=sqrt(0.753*0.753+0.009*0.009);
names[403]="P6p_B0Kstar0mumu_1.1-2      ";

central[404]=-0.244;
errors[404]=sqrt(0.645*0.645+0.012*0.012);
names[404]="P8p_B0Kstar0mumu_1.1-2      ";

central[405]=0.055;
errors[405]=sqrt(0.756*0.756+0.007*0.007);
names[405]="P1_B0Kstar0mumu_2-3         ";

central[406]=-0.323;
errors[406]=sqrt(0.316*0.316+0.033*0.033);
names[406]="P2_B0Kstar0mumu_2-3         ";

central[407]=0.005;
errors[407]=sqrt(0.364*0.364+0.012*0.012);
names[407]="P3_B0Kstar0mumu_2-3         ";

central[408]=-0.765;
errors[408]=sqrt(0.359*0.359+0.099*0.099);
names[408]="P4p_B0Kstar0mumu_2-3        ";

central[409]=-0.549;
errors[409]=sqrt(0.393*0.393+0.005*0.005);
names[409]="P6p_B0Kstar0mumu_2-3        ";

central[410]=-0.393;
errors[410]=sqrt(0.388*0.388+0.002*0.002);
names[410]="P8p_B0Kstar0mumu_2-3        ";

central[411]=0.421;
errors[411]=sqrt(18.35*18.25+0.018*0.018);
names[411]="P1_B0Kstar0mumu_3-4         ";

central[412]=-0.117;
errors[412]=sqrt(4.435*4.435+0.015*0.015);
names[412]="P2_B0Kstar0mumu_3-4         ";

central[413]=0.905;
errors[413]=sqrt(17.51*17.51+0.009*0.009);
names[413]="P3_B0Kstar0mumu_3-4         ";

central[414]=-0.134;
errors[414]=sqrt(1.343*1.343+0.108*0.108);
names[414]="P4p_B0Kstar0mumu_3-4        ";

central[415]=0.449;
errors[415]=sqrt(19.04*19.04+0.007*0.007);
names[415]="P6p_B0Kstar0mumu_3-4        ";

central[416]=0.303;
errors[416]=sqrt(1.394*1.394+0.006*0.006);
names[416]="P8p_B0Kstar0mumu_3-4        ";

central[417]=2.296;
errors[417]=sqrt(17.71*17.71+0.024*0.024);
names[417]="P1_B0Kstar0mumu_4-5         ";

central[418]=0.174;
errors[418]=sqrt(3.034*3.034+0.010*0.010);
names[418]="P2_B0Kstar0mumu_4-5         ";

central[419]=-0.801;
errors[419]=sqrt(17.42*17.42+0.007*0.007);
names[419]="P3_B0Kstar0mumu_4-5         ";

central[420]=-0.415;
errors[420]=sqrt(1.911*1.911+0.104*0.104);
names[420]="P4p_B0Kstar0mumu_4-5        ";

central[421]=-0.215;
errors[421]=sqrt(1.243*1.243+0.006*0.006);
names[421]="P6p_B0Kstar0mumu_4-5        ";

central[422]=0.293;
errors[422]=sqrt(1.522*1.522+0.006*0.006);
names[422]="P8p_B0Kstar0mumu_4-5        ";

central[423]=-0.540;
errors[423]=sqrt(1.100*1.100+0.025*0.025);
names[423]="P1_B0Kstar0mumu_5-6         ";

central[424]=0.089;
errors[424]=sqrt(0.227*0.227+0.012*0.012);
names[424]="P2_B0Kstar0mumu_5-6         ";

central[425]=0.178;
errors[425]=sqrt(0.465*0.465+0.007*0.007);
names[425]="P3_B0Kstar0mumu_5-6         ";

central[426]=-0.561;
errors[426]=sqrt(0.465*0.465+0.101*0.101);
names[426]="P4p_B0Kstar0mumu_5-6        ";

central[427]=0.074;
errors[427]=sqrt(0.309*0.309+0.005*0.005);
names[427]="P6p_B0Kstar0mumu_5-6        ";

central[428]=-0.068;
errors[428]=sqrt(0.372*0.372+0.006*0.006);
names[428]="P8p_B0Kstar0mumu_5-6        ";

central[429]=-0.353;
errors[429]=sqrt(0.602*0.602+0.026*0.026);
names[429]="P1_B0Kstar0mumu_6-7         ";

central[430]=0.104;
errors[430]=sqrt(0.136*0.136+0.013*0.013);
names[430]="P2_B0Kstar0mumu_6-7         ";

central[431]=-0.161;
errors[431]=sqrt(0.291*0.291+0.001*0.001);
names[431]="P3_B0Kstar0mumu_6-7         ";

central[432]=-0.641;
errors[432]=sqrt(0.294*0.294+0.106*0.106);
names[432]="P4p_B0Kstar0mumu_6-7        ";

central[433]=0.017;
errors[433]=sqrt(0.267*0.267+0.007*0.007);
names[433]="P6p_B0Kstar0mumu_6-7        ";

central[434]=0.162;
errors[434]=sqrt(0.289*0.289+0.005*0.005);
names[434]="P8p_B0Kstar0mumu_6-7        ";

central[435]=-0.284;
errors[435]=sqrt(0.548*0.548+0.025*0.025);
names[435]="P1_B0Kstar0mumu_7-8         ";

central[436]=0.393;
errors[436]=sqrt(0.231*0.231+0.013*0.013);
names[436]="P2_B0Kstar0mumu_7-8         ";

central[437]=-0.063;
errors[437]=sqrt(0.298*0.298+0.002*0.002);
names[437]="P3_B0Kstar0mumu_7-8         ";

central[438]=-0.503;
errors[438]=sqrt(0.288*0.288+0.118*0.118);
names[438]="P4p_B0Kstar0mumu_7-8        ";

central[439]=-0.201;
errors[439]=sqrt(0.274*0.274+0.007*0.007);
names[439]="P6p_B0Kstar0mumu_7-8        ";

central[440]=-0.623;
errors[440]=sqrt(0.295*0.295+0.005*0.005);
names[440]="P8p_B0Kstar0mumu_7-8        ";

central[441]=-0.869;
errors[441]=sqrt(0.408*0.408+0.030*0.030);
names[441]="P1_B0Kstar0mumu_11-11.75    ";

central[442]=0.494;
errors[442]=sqrt(0.134*0.134+0.013*0.013);
names[442]="P2_B0Kstar0mumu_11-11.75    ";

central[443]=0.166;
errors[443]=sqrt(0.221*0.221+0.005*0.005);
names[443]="P3_B0Kstar0mumu_11-11.75    ";

central[444]=-0.522;
errors[444]=sqrt(0.222*0.222+0.128*0.128);
names[444]="P4p_B0Kstar0mumu_11-11.75   ";

central[445]=-0.233;
errors[445]=sqrt(0.227*0.227+0.004*0.004);
names[445]="P6p_B0Kstar0mumu_11-11.75   ";

central[446]=-0.159;
errors[446]=sqrt(0.250*0.250+0.007*0.007);
names[446]="P8p_B0Kstar0mumu_11-11.75   ";

central[447]=-1.002;
errors[447]=sqrt(1.360*1.360+0.030*0.030);
names[447]="P1_B0Kstar0mumu_11.75-12.5  ";

central[448]=0.637;
errors[448]=sqrt(0.599*0.599+0.008*0.008);
names[448]="P2_B0Kstar0mumu_11.75-12.5  ";

central[449]=-0.105;
errors[449]=sqrt(0.42*0.42+0.004*0.004);
names[449]="P3_B0Kstar0mumu_11.75-12.5  ";

central[450]=-0.701;
errors[450]=sqrt(0.342*0.342+0.114*0.114);
names[450]="P4p_B0Kstar0mumu_11.75-12.5 ";

central[451]=-0.473;
errors[451]=sqrt(0.344*0.344+0.004*0.004);
names[451]="P6p_B0Kstar0mumu_11.75-12.5 ";

central[452]=-0.211;
errors[452]=sqrt(0.274*0.274+0.007*0.007);
names[452]="P8p_B0Kstar0mumu_11.75-12.5 ";

central[453]=-0.199;
errors[453]=sqrt(0.285*0.285+0.014*0.014);
names[453]="P1_B0Kstar0mumu_15-16       ";

central[454]=0.433;
errors[454]=sqrt(0.074*0.074+0.005*0.005);
names[454]="P2_B0Kstar0mumu_15-16       ";

central[455]=0.087;
errors[455]=sqrt(0.144*0.144+0.007*0.007);
names[455]="P3_B0Kstar0mumu_15-16       ";

central[456]=-0.673;
errors[456]=sqrt(0.199*0.199+0.013*0.013);
names[456]="P4p_B0Kstar0mumu_15-16      ";

central[457]=0.083;
errors[457]=sqrt(0.189*0.189+0.004*0.004);
names[457]="P6p_B0Kstar0mumu_15-16      ";

central[458]=-0.120;
errors[458]=sqrt(0.198*0.198+0.010*0.010);
names[458]="P8p_B0Kstar0mumu_15-16      ";

central[459]=-0.726;
errors[459]=sqrt(0.241*0.241+0.014*0.014);
names[459]="P1_B0Kstar0mumu_16-17       ";

central[460]=0.430;
errors[460]=sqrt(0.063*0.063+0.007*0.007);
names[460]="P2_B0Kstar0mumu_16-17       ";

central[461]=0.019;
errors[461]=sqrt(0.122*0.122+0.006*0.006);
names[461]="P3_B0Kstar0mumu_16-17       ";

central[462]=-0.552;
errors[462]=sqrt(0.213*0.213+0.014*0.014);
names[462]="P4p_B0Kstar0mumu_16-17      ";

central[463]=0.328;
errors[463]=sqrt(0.195*0.195+0.012*0.012);
names[463]="P6p_B0Kstar0mumu_16-17      ";

central[464]=0.122;
errors[464]=sqrt(0.199*0.199+0.010*0.010);
names[464]="P8p_B0Kstar0mumu_16-17      ";

central[465]=-0.313;
errors[465]=sqrt(0.293*0.293+0.019*0.019);
names[465]="P1_B0Kstar0mumu_17-18       ";

central[466]=0.288;
errors[466]=sqrt(0.075*0.075+0.006*0.006);
names[466]="P2_B0Kstar0mumu_17-18       ";

central[467]=0.144;
errors[467]=sqrt(0.149*0.149+0.002*0.002);
names[467]="P3_B0Kstar0mumu_17-18       ";

central[468]=-0.486;
errors[468]=sqrt(0.200*0.200+0.092*0.092);
names[468]="P4p_B0Kstar0mumu_17-18      ";

central[469]=0.047;
errors[469]=sqrt(0.198*0.198+0.023*0.023);
names[469]="P6p_B0Kstar0mumu_17-18      ";

central[470]=-0.006;
errors[470]=sqrt(0.215*0.215+0.001*0.001);
names[470]="P8p_B0Kstar0mumu_17-18      ";

central[471]=-0.450;
errors[471]=sqrt(0.447*0.447+0.022*0.022);
names[471]="P1_B0Kstar0mumu_18-19       ";

central[472]=0.393;
errors[472]=sqrt(0.159*0.159+0.011*0.011);
names[472]="P2_B0Kstar0mumu_18-19       ";

central[473]=0.134;
errors[473]=sqrt(0.219*0.219+0.010*0.010);
names[473]="P3_B0Kstar0mumu_18-19       ";

central[474]=-1.221;
errors[474]=sqrt(0.388*0.388+0.119*0.119);
names[474]="P4p_B0Kstar0mumu_18-19      ";

central[475]=0.128;
errors[475]=sqrt(0.265*0.265+0.012*0.012);
names[475]="P6p_B0Kstar0mumu_18-19      ";

central[476]=0.300;
errors[476]=sqrt(0.297*0.297+0.022*0.022);
names[476]="P8p_B0Kstar0mumu_18-19      ";


