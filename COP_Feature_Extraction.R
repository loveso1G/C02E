#
rm(list=ls())
gc()
# install.packages("pracma")#"#" must be removed for the first installation
# install.packages("EMD")#"#" must be removed for the first installation
# signal: COPx, COPy
t1=Sys.time();print(t1)
rawdata=read.csv("BDS00001.csv",header=T);dim(rawdata);#Just change BDS00001 and make changes according to the input data
tt=rawdata[,1];copx=rawdata[,8];
copy=rawdata[,9];
# EMD
library(EMD)
E1=emd(copx,tt,boundary="wave");
E1.no=E1$nimf;
E1.imf=E1$imf
E1.residue=E1$residue

E2=emd(copy,tt,boundary="wave");
E2.no=E2$nimf;
E2.imf=E2$imf
E2.residue=E2$residue

# time domain features
Result_all=matrix(0,1,80);
	# Feature-1: RMSD
m0E1=sqrt(mean((copx-mean(copx))^2))            ;Result_all[1]=m0E1 #x.0.1
m1E1=sqrt(mean((E1.imf[,1]-mean(E1.imf[,1]))^2));Result_all[2]=m1E1 #x.1.1
m2E1=sqrt(mean((E1.imf[,2]-mean(E1.imf[,2]))^2));Result_all[3]=m2E1 #x.2.1
m3E1=sqrt(mean((E1.imf[,3]-mean(E1.imf[,3]))^2));Result_all[4]=m3E1 #x.3.1
m4E1=sqrt(mean((E1.imf[,4]-mean(E1.imf[,4]))^2));Result_all[5]=m4E1 #x.4.1
m5E1=sqrt(mean((E1.imf[,5]-mean(E1.imf[,5]))^2));Result_all[6]=m5E1 #x.5.1
m6E1=sqrt(mean((E1.imf[,6]-mean(E1.imf[,6]))^2));Result_all[7]=m6E1 #x.6.1

m0E2=sqrt(mean((copy-mean(copy))^2))            ;Result_all[9]=m0E2 #y.0.1
m1E2=sqrt(mean((E2.imf[,1]-mean(E2.imf[,1]))^2));Result_all[10]=m1E2#y.1.1
m2E2=sqrt(mean((E2.imf[,2]-mean(E2.imf[,2]))^2));Result_all[11]=m2E2#y.2.1
m3E2=sqrt(mean((E2.imf[,3]-mean(E2.imf[,3]))^2));Result_all[12]=m3E2#y.3.1
m4E2=sqrt(mean((E2.imf[,4]-mean(E2.imf[,4]))^2));Result_all[13]=m4E2#y.4.1
m5E2=sqrt(mean((E2.imf[,5]-mean(E2.imf[,5]))^2));Result_all[14]=m5E2#y.5.1
m6E2=sqrt(mean((E2.imf[,6]-mean(E2.imf[,6]))^2));Result_all[15]=m6E2#y.6.1

SD0E1=sd(copx);	  	
SD1E1=sd(E1.imf[,1]);
SD2E1=sd(E1.imf[,2]);
SD3E1=sd(E1.imf[,3]);
SD4E1=sd(E1.imf[,4]);
SD5E1=sd(E1.imf[,5]);
SD6E1=sd(E1.imf[,6]);

SD0E2=sd(copy)
SD1E2=sd(E2.imf[,1]);
SD2E2=sd(E2.imf[,2]);
SD3E2=sd(E2.imf[,3]);
SD4E2=sd(E2.imf[,4]);
SD5E2=sd(E2.imf[,5]);
SD6E2=sd(E2.imf[,6]);

# http://127.0.0.1:17734/library/pracma/html/entropy.html
	#Feature-2:  Approximate Entropy 
 library(pracma)
print("Approximate Entropy E1");Sys.time()
AP0E1=approx_entropy(copx,edim=2,r=0.2*SD0E1,elag = 1);	 Result_all[17]=AP0E1; #x.0.4
AP1E1=approx_entropy(E1.imf[,1],edim=2,r=0.2*SD1E1,elag = 1);Result_all[18]=AP1E1; #x.1.4
AP2E1=approx_entropy(E1.imf[,2],edim=2,r=0.2*SD2E1,elag = 1);Result_all[19]=AP2E1; #x.2.4
AP3E1=approx_entropy(E1.imf[,3],edim=2,r=0.2*SD3E1,elag = 1);Result_all[20]=AP3E1; #x.3.4
AP4E1=approx_entropy(E1.imf[,4],edim=2,r=0.2*SD4E1,elag = 1);Result_all[21]=AP4E1; #x.4.4
AP5E1=approx_entropy(E1.imf[,5],edim=2,r=0.2*SD5E1,elag = 1);Result_all[22]=AP5E1; #x.5.4
AP6E1=approx_entropy(E1.imf[,6],edim=2,r=0.2*SD6E1,elag = 1);Result_all[23]=AP6E1; #x.6.4


print("Approximate Entropy E2");Sys.time()
AP0E2=approx_entropy(copy,edim=2,r=0.2*SD0E2,elag = 1);	 Result_all[25]=AP0E2; #y.0.4
AP1E2=approx_entropy(E2.imf[,1],edim=2,r=0.2*SD1E2,elag = 1);Result_all[26]=AP1E2; #y.1.4
AP2E2=approx_entropy(E2.imf[,2],edim=2,r=0.2*SD2E2,elag = 1);Result_all[27]=AP2E2; #y.2.4
AP3E2=approx_entropy(E2.imf[,3],edim=2,r=0.2*SD3E2,elag = 1);Result_all[28]=AP3E2; #y.3.4
AP4E2=approx_entropy(E2.imf[,4],edim=2,r=0.2*SD4E2,elag = 1);Result_all[29]=AP4E2; #y.4.4
AP5E2=approx_entropy(E2.imf[,5],edim=2,r=0.2*SD5E2,elag = 1);Result_all[30]=AP5E2; #y.5.4
AP6E2=approx_entropy(E2.imf[,6],edim=2,r=0.2*SD6E2,elag = 1);Result_all[31]=AP6E2; #y.6.4


	#Feature-3:  Sample entropy

print("Sample Entropy E1");Sys.time()
SA0E1=sample_entropy(copx,edim=2,r=0.2*SD1E1,tau = 1);	Result_all[33]=SA0E1; #x.0.5
SA1E1=sample_entropy(E1.imf[,1],edim=2,r=0.2*SD1E1,tau = 1);Result_all[34]=SA1E1; #x.1.5
SA2E1=sample_entropy(E1.imf[,2],edim=2,r=0.2*SD2E1,tau = 1);Result_all[35]=SA2E1; #x.2.5
SA3E1=sample_entropy(E1.imf[,3],edim=2,r=0.2*SD3E1,tau = 1);Result_all[36]=SA3E1; #x.3.5
SA4E1=sample_entropy(E1.imf[,4],edim=2,r=0.2*SD4E1,tau = 1);Result_all[37]=SA4E1; #x.4.5
SA5E1=sample_entropy(E1.imf[,5],edim=2,r=0.2*SD5E1,tau = 1);Result_all[38]=SA5E1; #x.5.5
SA6E1=sample_entropy(E1.imf[,6],edim=2,r=0.2*SD6E1,tau = 1);Result_all[39]=SA6E1; #x.6.5

print("Sample Entropy E2");Sys.time()
SA0E2=sample_entropy(copy,edim=2,r=0.2*SD1E1,tau = 1);	Result_all[41]=SA0E2; #y.0.5
SA1E2=sample_entropy(E2.imf[,1],edim=2,r=0.2*SD1E2,tau = 1);Result_all[42]=SA1E2; #y.1.5
SA2E2=sample_entropy(E2.imf[,2],edim=2,r=0.2*SD2E2,tau = 1);Result_all[43]=SA2E2; #y.2.5
SA3E2=sample_entropy(E2.imf[,3],edim=2,r=0.2*SD3E2,tau = 1);Result_all[44]=SA3E2; #y.3.5
SA4E2=sample_entropy(E2.imf[,4],edim=2,r=0.2*SD4E2,tau = 1);Result_all[45]=SA4E2; #y.4.5
SA5E2=sample_entropy(E2.imf[,5],edim=2,r=0.2*SD5E2,tau = 1);Result_all[46]=SA5E2; #y.5.5
SA6E2=sample_entropy(E2.imf[,6],edim=2,r=0.2*SD6E2,tau = 1);Result_all[47]=SA6E2; #y.6.5

	# Feature-4 Frequency power 
# install.packages("psd")
print("Frequency domain features");Sys.time()
library(psd)
f0E1=pspectrum(copx);
f1E1=pspectrum(E1.imf[,1]);
f2E1=pspectrum(E1.imf[,2]);
f3E1=pspectrum(E1.imf[,3]);
f4E1=pspectrum(E1.imf[,4]);
f5E1=pspectrum(E1.imf[,5]);
f6E1=pspectrum(E1.imf[,6]);

P0E1=sum(f0E1$spec); Result_all[49]=P0E1; #x.0.3
P1E1=sum(f1E1$spec); Result_all[50]=P1E1; #x.1.3
P2E1=sum(f2E1$spec); Result_all[51]=P2E1; #x.2.3
P3E1=sum(f3E1$spec); Result_all[52]=P3E1; #x.3.3
P4E1=sum(f4E1$spec); Result_all[53]=P4E1; #x.4.3
P5E1=sum(f5E1$spec); Result_all[54]=P5E1; #x.5.3
P6E1=sum(f6E1$spec); Result_all[55]=P6E1; #x.6.3

f0E2=pspectrum(copy);
f1E2=pspectrum(E2.imf[,1]);
f2E2=pspectrum(E2.imf[,2]);
f3E2=pspectrum(E2.imf[,3]);
f4E2=pspectrum(E2.imf[,4]);
f5E2=pspectrum(E2.imf[,5]);
f6E2=pspectrum(E2.imf[,6]);

P0E2=sum(f0E2$spec); Result_all[57]=P0E2; #y.0.3
P1E2=sum(f1E2$spec); Result_all[58]=P1E2; #y.1.3
P2E2=sum(f2E2$spec); Result_all[59]=P2E2; #y.2.3
P3E2=sum(f3E2$spec); Result_all[60]=P3E2; #y.3.3
P4E2=sum(f4E2$spec); Result_all[61]=P4E2; #y.4.3
P5E2=sum(f5E2$spec); Result_all[62]=P5E2; #y.5.3
P6E2=sum(f6E2$spec); Result_all[63]=P6E2; #y.6.3

	# Feature-5 median frequency, unit=Hz
c0E1=diff(sign(cumsum(f0E1$spec)-(P0E1/2)));Result_all[65]=which.max(c0E1)/60; #x.0.2
c1E1=diff(sign(cumsum(f1E1$spec)-(P1E1/2)));Result_all[66]=which.max(c1E1)/60; #x.1.2
c2E1=diff(sign(cumsum(f2E1$spec)-(P2E1/2)));Result_all[67]=which.max(c2E1)/60; #x.2.2
c3E1=diff(sign(cumsum(f3E1$spec)-(P3E1/2)));Result_all[68]=which.max(c3E1)/60; #x.3.2
c4E1=diff(sign(cumsum(f4E1$spec)-(P4E1/2)));Result_all[69]=which.max(c4E1)/60; #x.4.2
c5E1=diff(sign(cumsum(f5E1$spec)-(P5E1/2)));Result_all[70]=which.max(c5E1)/60; #x.5.2
c6E1=diff(sign(cumsum(f6E1$spec)-(P6E1/2)));Result_all[71]=which.max(c6E1)/60; #x.6.2

c0E2=diff(sign(cumsum(f0E2$spec)-(P0E2/2)));Result_all[73]=which.max(c0E2)/60; #y.0.2
c1E2=diff(sign(cumsum(f1E2$spec)-(P1E2/2)));Result_all[74]=which.max(c1E2)/60; #y.1.2
c2E2=diff(sign(cumsum(f2E2$spec)-(P2E2/2)));Result_all[75]=which.max(c2E2)/60; #y.2.2
c3E2=diff(sign(cumsum(f3E2$spec)-(P3E2/2)));Result_all[76]=which.max(c3E2)/60; #y.3.2
c4E2=diff(sign(cumsum(f4E2$spec)-(P4E2/2)));Result_all[77]=which.max(c4E2)/60; #y.4.2
c5E2=diff(sign(cumsum(f5E2$spec)-(P5E2/2)));Result_all[78]=which.max(c5E2)/60; #y.5.2
c6E2=diff(sign(cumsum(f6E2$spec)-(P6E2/2)));Result_all[79]=which.max(c6E2)/60; #y.6.2

print(Result_all);
print("Sample Entropy E2");Sys.time()