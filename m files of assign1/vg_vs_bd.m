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
k1 = (Na/ni)^(-2);
Vcb = 0; 
Vgb = linspace(-2, 2, 100);
for i=1:100
Vs0 = .1;
    func = @(Vs) (Vfb+ ((sign(Vs)).*(k.*((Vs-Vt+(Vt*exp(-Vs/Vt)))+  k1.*(-Vs-(Vt*(exp(-Vcb/Vt)))+(Vt*exp((Vs-Vcb)/Vt)))).^(1/2)))/Cox + Vs -Vgb(i));
    z(1,i) = fzero(func,Vs0);
end
plot(Vgb,z)
xlabel("Gate to Body Voltage Vgb(V)");
ylabel("Surface potential Vs(V)");
title("Vs vs Vgb");
