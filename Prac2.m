% Ramirez Rojas Gabriel Alejandro %
% 220477725     INRO   %

% ACT: Integración Sistemas EDO con Matlab %
% Archivo: dinamica_prac2.m
function dx = din_prac2(t, x)

    % Parámetros físicos
    Ip = 0.0079;      
    Mc = 0.7031;      
    lp = 0.3302;      
    Mp = 0.23;        
    Fc = 0;           
    Beq= 4.3;         
    g  = 9.81;        
    Bp = 0.0024;      

    xc   = x(1); 
    vc   = x(2);             
    th   = x(3);             
    w    = x(4);              
% Denominador común
    D = (Mc+Mp)*Ip + Mc*Mp*lp^2 + Mp^2*lp^2*sin(th)^2;

    acc_carro = ((Ip+Mp*lp^2)*Fc ...
                + Mp^2*lp^2*g*cos(th)*sin(th) ...
                - (Ip+Mp*lp^2)*Beq*vc ...
                - (Ip*Mp*lp - Mp^2*lp^3)*w^2*sin(th) ...
                - Mp*lp*w*cos(th)*Bp) / D;
    acc_pend = ((Mc+Mp)*Mp*g*lp*sin(th) ...
               - (Mc+Mp)*Bp*w ...
               + Fc*Mp*lp*cos(th) ...
               - Mp^2*lp^2*w^2*sin(th)*cos(th) ...
               - Beq*Mp*lp*vc*cos(th)) / D;

    % Vector de derivadas
    dx = [vc; acc_carro; w; acc_pend];
end
% ------------------------------------

CODIGO HECHO Y EJECUTADO EN UNA VENTANA APARTE EN MATLAB ODE 45      
% Archivo principal: main_prac2.m
clc; clear; close all;


[t, z] = ode45(@din_prac2, [0, 10], [0, 0, 0.1, 0]);

figure(1);
subplot(2, 1, 1);
plot(t, z(:,1), 'b', 'LineWidth', 1.5); 
hold on;
plot(t, z(:,2), 'r', 'LineWidth', 1.5);
hold off;
xlabel('Tiempo (s)');
ylabel('Carrito');
title('Variables del Carrito');
legend('Posición (x_c)', 'Velocidad (dx_c/dt)');
grid on;

subplot(2, 1, 2);
plot(t, z(:,3), 'b', 'LineWidth', 1.5); 
hold on;
plot(t, z(:,4), 'r', 'LineWidth', 1.5);
hold off;
xlabel('Tiempo (s)');
ylabel('Péndulo');
title('Variables del Péndulo');
legend('\alpha(t)', 'd\alpha/dt');
grid on;
