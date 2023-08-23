clc
clear all
close all

%% EXERCICIO 04 - IMPLEMENTACAO DO MODELO DINAMICO
% Utilizacao do modelo desenvolvido em sala

m = 1000;               %[kg] massa do modelo
c = 1;                  %[m] corda
b = c/2;                %[m] semi-corda
a = -0.5;               %[-]
d = 0.1;                %[-]
Iee = 100;              %[m^4] Momento de Inercia
kf = 20;                %[] rigidez flexional
kt = 10;                 %[] rigidez torsional

A = [m -(a-d)*b*m;-(a-d)*b*m Iee] % matriz de inercia
E = [kf 0;0 kt]                   % matriz de rigidez

%% EXERCICIO 05 - FREQUENCIAS NATURAIS E COEFS. DE AMORTECIMENTO
% Definindo a matriz de amortecimento
lamb1 = 0.05;   %parametro de ajuste 1
lamb2 = 0.01;   %parametro de ajuste 2
C = lamb1*A + lamb2*E

% Definindo as componentes e a matriz para o problema de autovalores
phi_11 = zeros(2);
phi_12 = eye(2);
phi_21 = -inv(A)*E;
phi_22 = -inv(A)*C;
PHI = [phi_11 phi_12;phi_21 phi_22]

% Autovetores (q) e autovalores (lamb) para a matriz PHI
[q,lamb] = eig(PHI)

% coletando os coeficientes de amortecimento 
for np=1:length(lamb)
    zeta_j(np) = -real(lamb(np,np))/(norm(lamb(np,np)));
    omega_j(np) = real(lamb(np,np))/(-zeta_j(np));
end

% Mudanca da posicao e velocidade para algum modo

tempo = 0:0.1:150;

for modo = 1:4
% -zeta_j(1)*omega_j(1)+1i*omega_j(1)*sqrt(1-zeta_j(1)^2)
    for nt=1:length(tempo)
    z(nt) = q(1,modo)*exp(lamb(modo,modo)*tempo(nt));
    theta(nt) = q(2,modo)*exp(lamb(modo,modo)*tempo(nt));
    z_dot(nt) = q(3,modo)*exp(lamb(modo,modo)*tempo(nt));
    theta_dot(nt) = q(4,modo)*exp(lamb(modo,modo)*tempo(nt));
    end
    
    subplot(2,2,1)
    plot(tempo,z)
    ylabel('Posicao Z')
    hold on
    subplot(2,2,2)
    plot(tempo,theta)
    ylabel('Posicao Theta')
    hold on
    subplot(2,2,3)
    plot(tempo,z_dot)
    ylabel('Velocidade Z')
    hold on
    subplot(2,2,4)
    plot(tempo,theta_dot)
    ylabel('Velocidade Theta')
    hold on
    sgtitle("Resolucao Analitica")
end
%% EXERCICIO 6 - METODO DE EULER E PREVISOR/CORRETOR
% Condicoes iniciais para ambos os metodos
z0 = -0.7;        % posicao z inicial
zdot0 = -0.2;     % velocidade z inicial
theta0 = -0.5;     % posicao theta inicial
thetadot0 = 0.1;  % velocidade theta inicial
deltaT = 0.001;   % passo no tempo
Np = 150/deltaT;  % numero de iteracoes para 150 segundos

% Metodo de Euler
for ktt = 1:Np
    if ktt==1
        XDOT(:,ktt) = [zdot0;thetadot0];
        X(:,ktt) = [z0;theta0];
    else
        X2DOT(:,ktt) = inv(A)*(-C*XDOT(:,ktt-1)-E*X(:,ktt-1));
        XDOT(:,ktt) = XDOT(:,ktt-1) + X2DOT(:,ktt)*deltaT;
        X(:,ktt) = X(:,ktt-1) + XDOT(:,ktt-1)*deltaT;
    end
end
tempo_euler = deltaT:deltaT:150;
figure(2)
subplot(2,2,1)
plot(tempo_euler,X(1,:))
ylabel('Posicao Z')

subplot(2,2,2)
plot(tempo_euler,X(2,:))
ylabel('Posicao Theta')

subplot(2,2,3)
plot(tempo_euler,XDOT(1,:))
ylabel('Velocidade Z')

subplot(2,2,4)
plot(tempo_euler,XDOT(2,:))
ylabel('Velocidade Theta')
sgtitle("Metodo Euler")

%%
% Metodo Previsor/Corretor

for ktt = 1:Np
    if ktt==1
        XDOTc(:,ktt) = [zdot0;thetadot0];
        Xc(:,ktt) = [z0;theta0];
        % aceleracao inicial
        X2DOTc(:,ktt) = inv(A)*(-C*XDOTc(:,ktt)-E*Xc(:,ktt));
    else
        % previsao aceleracao, velocidade e posicao
        XPDOTc(:,ktt) = XDOTc(:,ktt-1) + X2DOTc(:,ktt-1)*deltaT;
        XPc(:,ktt) = Xc(:,ktt-1) + (XDOTc(:,ktt-1)+XPDOTc(:,ktt))*0.5*deltaT;
        X2PDOTc(:,ktt) = inv(A)*(-C*XPDOTc(:,ktt)-E*XPc(:,ktt));
        % correcao aceleracao, velocidade e posicao
        XDOTc(:,ktt) = XDOTc(:,ktt-1) + (X2DOTc(:,ktt-1)+X2PDOTc(:,ktt))*0.5*deltaT;
        Xc(:,ktt) = Xc(:,ktt-1) + (XDOTc(:,ktt-1)+XPDOTc(:,ktt))*0.5*deltaT;
        X2DOTc(:,ktt) = inv(A)*(-C*XDOTc(:,ktt)-E*Xc(:,ktt));
    end
    
end
figure(3)
subplot(2,2,1)
plot(tempo_euler,Xc(1,:))
ylabel('Posicao Z')

subplot(2,2,2)
plot(tempo_euler,Xc(2,:))
ylabel('Posicao Theta')

subplot(2,2,3)
plot(tempo_euler,XDOTc(1,:))
ylabel('Velocidade Z')

subplot(2,2,4)
plot(tempo_euler,XDOTc(2,:))
ylabel('Velocidade Theta')
sgtitle("Metodo Previsor/Corretor")

fprintf("Amortecimentos analiticos: ")
disp(zeta_j)
fprintf("Frequencias analiticas: ")
disp(omega_j)
%% Exercicio 7 - Obtencao 
% figure(4)
% 
% % USANDO METODO NUMERICO
% Fs = 1/max(tempo_euler);
% L = length(tempo_euler);
% Xomega = abs(fft(Xc(2,:)));
% 
% % ANALITICO
% % Fs = 1/max(tempo);
% % L = length(tempo);
% % Xomega = abs(fft(z(1,:)));
% 
% P2 = abs(Xomega/L);
% P1 = P2(:,1:L/2+1);
% P1(:,2:end-1) = 2*P1(:,2:end-1);
% omega = 2*pi*(0:(Fs/L):(Fs/2-Fs/L));
% plot(omega,P1(1,1:L/2),'LineWidth',1.5)
% ylabel("|X(w)|")
% 
% max(Xomega/(1*L))



