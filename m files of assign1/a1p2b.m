clc ;
clear all;
close all;
Dit=1e12;
ep=1.05e-12; q=1.6e-19;
Na=5e17;
Phim=4.1;
ni=1.5*(10^10);
Eg=1.1;
Vt=0.026;
chi=4.05;
Phis=(chi+(Eg/2)+(Vt*log(Na/ni)));
Phi_ms=Phim-Phis;
Qf=2*10^(-7);
tox1 =2e-7 ; 
eps_siO2 = 4*8.854*1e-14;
Cox=(eps_siO2/tox1);
Vfb=Phi_ms-(Qf/Cox);
Cox=0.345e-6;
phi_f=0.39;
phi_t=0.026;
a1=sqrt(2*ep*Na*q);
psi_s=-0.82:0.00001:1.9;
n=length(psi_s);
Cc1_prime=[];
Cch_prime=[];
y1=[];
y2=[];
Vgb=[];
Cgb_prime=[];
Cgbh_prime=[];
Qc1_prime=[];
for i=1:n
if Vgb==Vfb
Cgb_prime=sqrt(q*ep*Na/phi_t);
else
Qc_prime(1,i)=sign(psi_s(1,i))*a1*sqrt((phi_t*exp(-psi_s(1,i)/phi_t))+psi_s(1,i)-phi_t+(exp(-2*phi_f/phi_t))*(phi_t*exp(psi_s(1,i)/phi_t)-psi_s(1,i)-phi_t));
y1(1,i)=2*sqrt((phi_t*exp(-psi_s(1,i)/phi_t))+psi_s(1,i)-phi_t+(exp(-2*phi_f/phi_t))*(phi_t*exp(psi_s(1,i)/phi_t)-psi_s(1,i)-phi_t));
y2(1,i)=2*sqrt((phi_t*exp(-psi_s(1,i)/phi_t))+psi_s(1,i)-phi_t);
Cc1_prime(1,i)=sign(psi_s(1,i))*a1*((1-exp(-psi_s(1,i)/phi_t)+exp(-2*phi_f/phi_t)*(exp(psi_s(1,i)/phi_t)-1))/y1(1,i));
Cc1h_prime(1,i)=sign(psi_s(1,i))*a1*((1-exp(-psi_s(1,i)/phi_t))/y2(1,i));
Vgb(1,i)=Vfb+psi_s(1,i)+(-sign(psi_s(1,i))*a1*sqrt((phi_t*exp(-psi_s(1,i)/phi_t))+psi_s(1,i)-phi_t-(exp((-2*phi_f)/phi_t))*(phi_t*exp(psi_s(1,i)/phi_t)-psi_s(1,i)-phi_t)));
Cgb_prime(1,i)=Cox*Cc1_prime(1,i)/(Cox+Cc1_prime(1,i));
Cgbh_prime(1,i)=Cox*Cc1h_prime(1,i)/(Cox+Cc1h_prime(1,i));
Qit_prime(1,i)=q*Dit*(phi_f-psi_s(1,i));
Cit=q*Dit;
Cgb1_prime(1,i)=Cox*(Cc1_prime(1,i)+Cit)/(Cit+Cox+Cc1_prime(1,i));
Vgb2(1,i)=Vfb+psi_s(1,i)-(Qit_prime(1,i)/Cox)+(-sign(psi_s(1,i))*a1*sqrt((phi_t*exp(-psi_s(1,i)/phi_t))+psi_s(1,i)-phi_t-(exp((-2*phi_f)/phi_t))*(phi_t*exp(psi_s(1,i)/phi_t)-psi_s(1,i)-phi_t)));
end
end
% Cgb' and Vgb'
plot(psi_s,Cgb_prime,psi_s,Cgbh_prime,psi_s,Cgb1_prime,psi_s,Cgbh_prime,'LineWidth',3);
title(' Ctotal Vs Vgb');
xlabel('Vgb in volts ');
ylabel('Ctotal in F/cm^2');
% Cgb' and Vgb' With Dit HF and LF
plot(psi_s,Cgb_prime,psi_s,Cgbh_prime,'LineWidth',3);
title(' Ctotal Vs Vgb');
xlabel('Vgb in volts');
ylabel('Ctotal in F/cm^2');
legend('low freq without Dit','High freq without Dit','low freq withDit','High freq with Dit');