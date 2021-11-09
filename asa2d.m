clc 
clear all
close all
u = linspace(0,1,60);
x = u.^2;
y = 5*0.12*(0.2969*u - 0.1260*u.^2 - 0.3516*u.^4 + 0.2843*u.^6 - 0.1015*u.^8);

hold on
plot(x,y,'m')
plot(x,-y,'m');
xlabel('x/c')
axis equal;
ax = gca;
grid on
set(ax,'Color','k')
ax.GridColor = [0.52, 0.52, 0.52];