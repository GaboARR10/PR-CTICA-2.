% Ramirez Rojas Gabriel Alejandro %
% 220477725     INRO   %

% ACT: Integración Sistemas EDO con Matlab %
% Archivo: dinamica_prac2.m
function dx = din_prac2(t, x)

    % Parámetros físicos
    Ip  = 0.0079;      
    Mc  = 0.7031;      
    lp  = 0.3302;      
    Mp  = 0.23;        
    Fc  = 0;           
    Beq = 4.3;         
    g   = 9.81;        
    Bp  = 0.0024;

    dx = zeros(4, 1); 

    % Estados
    X_1      = x(1);
    X_2      = x(2);
    alpha    = x(3);
    alpha2_2 = x(4);

    % Pre-cálculos
    sen_alpha = sin(alpha);
    cos_alpha = cos(alpha);
    alpha_pot = alpha2_2^2;

    % Ecuaciones de estado
    dx(1) = X_2;

    dx(2) = (1 / ((Mc+Mp)*Ip + Mc*Mp*lp^2 + Mp^2*lp^2*sen_alpha^2)) * ...
            ( (Ip+Mp*lp^2)*Fc + Mp^2*lp^2*g*cos_alpha*sen_alpha ...
            - (Ip+Mp*lp^2)*Beq*X_2 ...
            - (Ip*Mp*lp - Mp^2*lp^3)*alpha_pot*sen_alpha ...
            - Mp*lp*alpha2_2*cos_alpha*Bp );

    dx(3) = alpha2_2;

    dx(4) = (1 / ((Mc+Mp)*Ip + Mc*Mp*lp^2 + Mp^2*lp^2*sen_alpha^2)) * ...
            ( (Mc+Mp)*Mp*g*lp*sen_alpha ...
            - (Mc+Mp)*Bp*alpha2_2 ...
            + Fc*Mp*lp*cos_alpha ...
            - Mp^2*lp^2*alpha_pot*sen_alpha*cos_alpha ...
            - Beq*Mp*lp*X_2*cos_alpha );

end
% ------------------------------------

CODIGO HECHO Y EJECUTADO EN UNA VENTANA APARTE EN MATLAB ODE 45      
% Condiciones iniciales
x0 = [0; 0; pi/180; 0];   % [posición; velocidad; ángulo(rad); vel. angular]
tspan = [0 10];           % intervalo de simulación

% Resolver con ODE45
[t, x] = ode45(@din_prac2, tspan, x0);

% Graficar resultados
figure;

subplot(2,1,1);
plot(t, x(:,1), 'y', 'LineWidth', 2);
xlabel('Tiempo [s]');
ylabel('Posición [m]');
title('Movimiento del carro');
grid on;

subplot(2,1,2);
plot(t, x(:,3)* 50*pi/180, 'y', 'LineWidth', 2);
xlabel('Tiempo [s]');
ylabel('Ángulo [°]');
title('Movimiento del péndulo');
grid on;
