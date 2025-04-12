clear all
close all
clc

% Static_Char()

F = [1.5 2 2.5 3 3.5 4 4.5 5];
F = F(1);

n = 100000;
h = 0.1;
    
T0 = 15;
tau = 19.08*F^(-0.4293) - 4.042;
tau0 = 11.93*F^(-0.78) + 2.73;

x = zeros(1,3);
x(1) = T0;
x(3) = T0;
y = zeros(n/10,1) + T0;
y_prev = 0;
u_vec = zeros(n/10,1);

T_Setpoint = 30;
u = (T_Setpoint - T0) * (5*4200*F)/(3*12000);

k = 0;
for i=1:n
    if k == 0
        u = (T_Setpoint - 3 - T0) * (5*4200*F)/(3*12000);
    elseif k > 0 && k < 10
        u = (T_Setpoint + 3 - T0) * (5*4200*F)/(3*12000);
    else
        break
    end
    u_vec(i) = u;

    G = -0.0002347*u + 1.012;
    x = RungeKutta(x, T0, F, u, G, tau, h);

    err = x(1) - y_prev;
    if err < 0.000017
        k = k + 1;
    end
    y_prev = x(1);
    y(i + round(tau0/h)) = x(1);
end
y = y(1:length(u_vec));

data = iddata(y, u_vec, 0.1);
model_fopdt = procest(data, 'P1D');

plot(y);

function [x] = RungeKutta(x, x0, F, u, G, tau, h)
    k11 = h*StateEquationObject(x(3), x0, F, 1.6, u, 12000, 4200, 1);
    k21 = h*StateEquationObject(x(3) + 0.5*k11, x0, F, 1.6, u, 12000, 4200, 1);
    k31 = h*StateEquationObject(x(3) + 0.5*k21, x0, F, 1.6, u, 12000, 4200, 1);
    k41 = h*StateEquationObject(x(3) + k31, x0, F, 1.6, u, 12000, 4200, 1);

    x(3) = x(3) + (k11 + 2*k21 + 2*k31 + k41)/6;

    k12 = h*StateEquationSensor(x(1:2), x(3), G, tau);
    k22 = h*StateEquationSensor(x(1:2) + 0.5*k12, x(3), G, tau);
    k32 = h*StateEquationSensor(x(1:2) + 0.5*k22, x(3), G, tau);
    k42 = h*StateEquationSensor(x(1:2) + k32, x(3), G, tau);

    x(1:2) = x(1:2) + (k12 + 2*k22 + 2*k32 + k42)/6;
end

function [dxdt] = StateEquationObject(x, x0, F, V, Ph, Pnom, cw, ro)
    dxdt = -F/(60*V)*x + F/(60*V)*x0 + (Ph*Pnom)/(100*cw*ro*V);
end

function [dxdt] = StateEquationSensor(x, u, G, tau)
    dxdt(1) = x(2);
    dxdt(2) = -2/tau*x(2) - 1/tau^2*x(1) + G*u/tau^2;
end

function Static_Char()
    % Tworzenie siatki dla P i F
    [P, F] = meshgrid(2:0.1:10, 0:0.1:100);
    
    % Parametry
    Tin = 15; % Parametr Tin
    Pnom = 12000; % Parametr Pnom
    cw = 4200; % Parametr cw
    ro = 1; % Parametr ro
    
    % Funkcja nieliniowa T0
    T0 = @(P, F) -0.0002347 * P * Tin + 1.012 * Tin + (3 * Pnom) ./ (5 * cw * ro) .* (-0.0002347 * P.^2 ./ F + 1.012 * P ./ F);
    
    % Tworzenie siatki dla pełnego zakresu P i F
    [P_grid_full, F_grid_full] = meshgrid(linspace(0, 100, 30), linspace(2, 10, 30));
    
    % Obliczenie wartości funkcji T0 na siatce
    T0_full = T0(P_grid_full, F_grid_full);
    
    % Wykres 3D
    figure;
    hold on;
    
    % Powierzchnia funkcji nieliniowej
    surf(P_grid_full, F_grid_full, T0_full, 'FaceAlpha', 0.5); % Przezroczystość 50%
    
    xlabel('P');
    ylabel('F');
    zlabel('T0');
    legend({'Funkcja nieliniowa'}, 'Location', 'best');
    colorbar;
    shading interp;
    grid on;
    hold off;
end