clc
clear all

%% ec = EE - EA = EE - 0.25c
c = 1.4;               %[m]
rho = 1.2;             %[kg/m3]
ec = 0.35*c-c/4;
a = 2*(ec/c - 1/4);
b = c/2;

%k = f*2*pi*(c/2)/U
k = 0.05;
Ck = theodorsen2(k)
G = imag(Ck);
F = real(Ck);

% componente da matriz de rigidez
Lz = 2*pi*(-k^2/2 - G*k);
bLteta = b*2*pi*((k^2)*(a/2) + F - G*k*(0.5-a));
bMz = b*2*pi*(-k^2*a/2 - k*(a+0.5)*G);
b2Mteta = (b^2)*2*pi*((k^2/2)*(1/8+a^2)+F*(a+0.5)-k*G*(a+0.5)*(0.5-a));
C = [Lz bLteta;bMz b2Mteta]
%rhoU2C = rho*(U^2)*C

% componente matriz de amortecimento
bLz_dot = b*2*pi*F;
b2Lteta_dot = b^2*2*pi*(0.5+F*(0.5-a)+G/k);
b2Mz_dot = b^2*2*pi*(a+0.5)*F;
b3Mteta_dot = b^3*2*pi*(-0.5*k*(0.5-a)+k*F*(a+0.5)*(0.5-a)+(G/k)*(a+0.5));
B = [bLz_dot b2Lteta_dot;b2Mz_dot b3Mteta_dot]
%rhoUB = rho*U*B

%% Funcao de Theodorsen
function C = theodorsen2(k)
  C = besselk(1,1i*k)/(besselk(0,1i*k)+besselk(1,1i*k));
end