%for vcb=1

clc;
epsilon_s = 11.6*8.854187817*(10^(-14)); %F/cm
epsilon_ox = 3.9*8.854187817*(10^(-14)); %F/cm
q = 1.6*(10^(-19));     %C
Na = 5*(10^17);     %cm-3
ni = 1.5*(10^10); %cm-3
tox = 2*(10^(-7));  %cm
Qf = 2*(10^(-7)); %C/cm2
Cox = epsilon_ox/tox;  %F/cm2
Vfb = Qf/Cox;   %Volts
Vt = .026;     %threshold voltage(V)
k = (2*q*epsilon_s*Na)^(1/2); %constant

V_surface1 =-5*Vt:.00001:2 ;%Volts
s = sign(V_surface1);    %sign of the surface potential
k1 = (Na/ni)^(-2);
Qs1 = -s.*k.*(  (V_surface1-Vt+(Vt*exp(-V_surface1/Vt)))+  k1.*(-V_surface1-Vt*(exp(-1/Vt))+(Vt*exp((V_surface1-1)/Vt)))).^(1/2);
subplot(2,1,1)
plot(V_surface1,Qs1)
xlabel("Surface Potential(Vsurface)")
ylabel("Surface Charge(Qs)")
title("Plot of Qs vs Vsurface")

subplot(2,1,2)
plot(V_surface1,log(Qs1))
xlabel("Surface Potential(Vsurface)")
ylabel("Surface Chargelog scale log(Qs)")
title("Plot of Qs vs Vsurface")


