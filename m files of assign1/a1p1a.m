clc
clear
na=5*10^17;
tox=2*10^-7;
q=1.6*10^-19; 
ep_si=106.2480*10^-14; 
ep_ox=35.4160*10^-14;
c_ox=ep_ox/tox;
eg=1.1;
ki=4.05;
ni=1.5*10^10;
fi_m=4.04;
pi_f=0.026*log(na/ni);
fi_sc=ki+(eg/2)+pi_f;
fi_ms=fi_m-fi_sc;
Q_fix=2*10^-7;
%vt=(fi_ms-(Q_fix/c_ox)+(2*pi_f))+(((2*ep_si*q*na*2*pi_f)^0.5)/c_ox);
vt=0.026;
x1 =-0.925 :0.0001:0.001; % surface potential variation from -5*vt to 2 volts
y1=1*((2*ep_si*q*na)^0.5)*(((vt*exp(-x1/vt)+x1-vt)+(exp((-2*pi_f)/vt)*((vt*exp(x1/vt))-x1-vt))));
x2 =0.001 :0.0001:2; % surface potential variation from -5*vt to 2 volts
y2=-1*((2*ep_si*q*na)^0.5)*(((vt*exp(-x2/vt)+x2-vt)+(exp((-2*pi_f)/vt)*((vt*exp(x2/vt))-x2-vt))));
x=[x1,x2];
y=[y1,y2];
t = tiledlayout(1,2);
nexttile
plot(x2,y2,"Marker",".","Color",'b')
xlabel('surface potential phi_sin volts')
ylabel('semiconductor charge Qs in columb/cm^2')
title('semiconductor charge vs surface potential')
grid on
nexttile
plot(x,log(y),"Marker",".","Color",'b')
xlabel('surface potential  of semiconductor in volts')
ylabel('semiconductor charge Qs in log scale columb/cm^2')
title('semiconductor charge Qs vs surface potential')
grid on