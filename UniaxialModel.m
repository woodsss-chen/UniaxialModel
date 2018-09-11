clear;
%% Input constants
% h = 16e-9; %m
global Kfc F Ef m Qc R T A sr Ys;

Kfc = 0.76e6; %Pa m^1/2
F = 2;
Ef = 52e9;
v = 0.21; %plane strain!
Ef = Ef/(1-v);
sr = 3e-4; %strain rate

m = 6.7;
Qc = 36e3; %J/mol
R = 8.31; %J/mol/K
T = 298; %K
A = 3.5e5;

Ys = sr^(1/m)*A*exp(Qc/m/R/T);
Ys = 2.28e6;

h = 60e-10:1e-10:(400*1e-10);
e_c = zeros(size(h));
d_su = e_c;
d_sd = e_c;
for i=1:size(h,2)
    e_c(i) = CriticalStrain(h(i));
    d_su(i) = 2*sqrt(3)*e_c(i)*Ef*h(i)/Ys;
    d_sd(i) = sqrt(3)*e_c(i)*Ef*h(i)/Ys;
end
figure;
plot(h,1./d_su); hold on;
plot(h,1./d_sd); % hold off;

d_fit = 1.416e-3./h;
plot(h,d_fit,'--');

plot(80e-10,0.3e6,'*');
plot(120e-10,0.15e6,'*');
plot(160e-10,0.16e6,'o');
plot(160e-10,0.088e6,'*');
plot(240e-10,0.06e6,'*');
plot(320e-10,0.044e6,'*');
xlabel('Film Thickness [m]');
ylabel('Saturation Crack Density [m^{-1}]');

function e_c = CriticalStrain(h)
global Kfc F Ef m Qc R T A sr Ys;
%% Solve for sigma_c
a = zeros(1,4);
a(1) = 1;
a(2) = sqrt(3)*pi*Ys*F;
a(4) = -sqrt(3)*Ys*Kfc^2/h;
sigma_c = roots(a);
sigma_c = sigma_c(imag(sigma_c)==0);
sigma_c = sigma_c(sigma_c>=0);
e_c = sigma_c/Ef;
end