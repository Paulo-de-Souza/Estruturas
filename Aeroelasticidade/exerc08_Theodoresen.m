clc
clear all
close all

k = 0.01:0.01:1.99;
for n=1:length(k)
    Ck(n) = theodorsen2(k(n));
end

subplot(1,2,1)
plot(k,real(Ck),'k','LineWidth',1.5)
title("Parte real da funcao de Theodorsen")
xlabel("Frequencia Adimensional k")
ylabel("F(k) = Re[C(k)]")
subplot(1,2,2)
plot(k,imag(Ck),'k','LineWidth',1.5)
title("Parte imaginaria da funcao de Theodorsen")
xlabel("Frequencia Adimensional k")
ylabel("G(k) = Im[C(k)]")

function C = theodorsen2(k)
  C = besselk(1,1i*k)/(besselk(0,1i*k)+besselk(1,1i*k));
end