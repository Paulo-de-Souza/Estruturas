clc
clear all
close all

X = 0.0065;    %K/m
rho_0 = 1.225; %Kg/m3
T_0 = 288.15;  %K

h = linspace(0,10000); %variacao altitude em [m]
T = T_0 - X.*h;
rho = rho_0*((T/T_0).^4.2524);
rho_aprox = rho_0*((20 - h./1000)./(20+h./1000));
erro = (1-(rho_aprox./rho))*100;
subplot(1,2,1)
plot(h,rho,'b','LineWidth',1.5,'DisplayName','Geral')
hold on
plot(h,rho,'r--','LineWidth',1.2,'DisplayName','Aproximado')
xlabel("Altura h em [m]")
ylabel("Massa especifica [kg/m^3]")
legend()
hold off
subplot(1,2,2)
plot(h,erro,'k','DisplayName','Erro')
xlabel("Altura h em [m]")
ylabel('Diferenca [%]')
legend()
hold off