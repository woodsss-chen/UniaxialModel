clear;
%% Input constants
% h = 16e-9; %m
global Kfc F Ef m Qc R T A sr Ys;

Kfc = 0.76e6; %Pa m^1/2
F = 2;
Ef = 52e9;
time = 1e3;
s = 0:0.2e-2:60e-2;
sr = s./time; %strain rate

m = 6.7;
Qc = 36e3; %J/mol
R = 8.31; %J/mol/K
T = 298; %K
A = 3.5e5;

Ys = sr.^(1/m)*A*exp(Qc/m/R/T);
Ys = Ys';

h = 16e-9;

%% Solve for sigma_c
a = zeros(size(Ys,1),4);
a(:,1) = 1;
a(:,2) = sqrt(3)*pi.*Ys*F;
a(:,4) = -sqrt(3).*Ys*Kfc^2/h;
sigma_c = zeros(size(a,1),3);

for i = 1:size(a,1)
    sigma_c(i,:) = roots(a(i,:));
end
sigma_c = sigma_c(imag(sigma_c)==0);
sigma_c = sigma_c(sigma_c>=0);
e_c = sigma_c/Ef;

d_su = zeros(size(Ys));
d_sd = zeros(size(Ys));
for i=1:size(Ys,1)
    d_su(i) = 2*sqrt(3)*e_c(i)*Ef*h/Ys(i);
    d_sd(i) = sqrt(3).*e_c(i)*Ef*h./Ys(i);
end


figure;
plot(s,1./d_su); hold on;
plot(s,1./d_sd); % hold off;
