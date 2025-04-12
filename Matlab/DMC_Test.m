clear all
close all
clc

n = 10000;

T_in = 15;
F_max = 12;
Ph = 80;
V = 1.6;
Pnom = 12000;
cw = 4200;
ro = 1;

x = zeros(1,3);
x(1) = T_in;
x(3) = T_in;
y = zeros(2*n,1) + T_in;

alpha = F_max/(100*60*V);
beta = (Ph*Pnom)/(100*cw*ro*V);
G = -0.0002347*Ph + 1.012;

T_OP = 40;
OR = 1/(alpha/(beta*G)*T_OP - alpha/beta*T_in);

for i=1:2*n
    if i < n
        OR = 1/(alpha/(beta*G)*(T_OP-5) - alpha/beta*T_in);
        u_norm(1) = OR;
    else
        OR = 1/(alpha/(beta*G)*(T_OP+5) - alpha/beta*T_in);
        u_norm(2) = OR;
    end
    F = OR*F_max/100;
    tau = 19.08*F^-0.4293 - 4.042;
    tau0 = 11.93*F^-0.78 + 2.37;
    x = RungeKutta(x, T_in, OR, Ph, tau);
    y(i + round(tau0*10):i+round(tau0*10)+n/10) = x(1);
end

Tp = 1;
N1 = round(tau0)/Tp;
s_raw = (y(n:end)-y(n))/(u_norm(2)-u_norm(1));
s = s_raw(1:Tp*10:end);

% time = 0:0.1:(length(y)-1)/10;
% plot(time, y)

%%
D = [1 40 50 75 100 126 200 300 500 1000];
N = [2 5 10 30 56 75 100 150 200 300];
Nu = [3 1 2 5 10 20 50];
l = [5 0.1 1 2.78 5 10 30 100 200 500 1000];

params = [126, 59, 2, N1, 2.78];

p = D;
params_end = repmat(params, length(p)-1, 1);
params_end(:, p(1)) = p(2:end)';

figure

for k=1:length(p)-1

    D = params_end(k,1);
    N = params_end(k,2);
    Nu = params_end(k,3);
    N1 = params_end(k,4);
    l = params_end(k,5);

    MP = zeros(N-N1+1, D-1);
    
    for i=1:N-N1+1
        for j=1:D-1
            if i+j+N1-1 < D
                MP(i,j) = s(i+j+N1-1) - s(j);
            else
                MP(i,j) = s(D) - s(j);
            end
        end
    end
    
    M = zeros(N-N1+1, Nu);
    
    for i=1:Nu
        for j=i:N-N1+1
            M(j,i) = s(j+N1-i);
        end
    end
    
    K = ((M'*M+l*eye(Nu))^-1)*M';
    
    Ku = K(1,:)*MP;
    Ke = sum(K(1,:));
    
    n = 20000;
    
    Ph = 50;
    F_max = 12;
    T0 = 15;
    
    x = zeros(1,3);
    x(1) = T0;
    x(3) = T0;
    y = zeros(n,1) + T0;
    u_vec = zeros(n,1);
    
    % T_Setpoint = 50;
    % u = (T_Setpoint - T0) * 21000/12000;
    % for i=1:n/2
    %     % Calculate Gain and process Runge-Kutta method for output approx
    %     G = -0.0002347*u + 1.012;
    %     x = RungeKutta(x, T0, u, G, tau);
    %     y(i + round(tau0*10)) = x(1);
    %     u_vec(i) = u;
    % end
    
    % Control law
    T_OP = 35;
    u_subtruct = 0;
    OR = 50;
    dOR = 0;
    dOR_prev = zeros(D-1,1);
    k = 0;
    for i=1:n
        if k == 10*Tp
            for j=1:D-1
                u_subtruct = u_subtruct + Ku(j)*dOR_prev(j);
            end
            dOR = Ke*(T_OP - y(i-1)) - u_subtruct;
            u_subtruct = 0;
    
            if dOR > 10
                dOR = 10;
            elseif dOR < -10
                dOR = -10;
            end
    
            OR = OR + dOR;
    
            if OR > 100
                OR = 100;
            elseif OR < 10
                OR = 20;
            end
    
            for j=D-1:-1:2
                dOR_prev(j) = dOR_prev(j-1);
            end
            dOR_prev(1) = dOR;
            k = 0;
        end
    
        if i == 15000
            T_OP = 40;
        end
        
        F = OR*F_max/100;
        tau = 19.08*F^-0.4293 - 4.042;
        tau0 = 11.93*F^-0.78 + 2.37;

        x = RungeKutta(x, T0, OR, Ph, tau);
        y(i + round(tau0*10):i+round(tau0*10)+30) = x(1);
        u_vec(i) = OR;
        k = k + 1;
    end
    
    % disp(max(y))
    
    [size_y, ~] = size(y);
    time = (0:0.1:(size_y-1)/10)';
    time_u = (0:0.1:(n-1)/10);
    
    T = 37.05;
    T0 = 21.55;
    gain = 0.5588;
    T_IN = 17;
    
    T_SP = 45;
    
    Kp = 0.6*T/(gain*T0);
    Ti = 0.8*T0+0.5*T;

    F_vec = u_vec*F_max/100;

    subplot(2,1,1)
    plot(time, y, 'LineWidth', 1)
    xlim([n*3/40 n/10])
    hold on
    grid on
    legend(D)
    subplot(2,1,2)
    stairs(time_u, F_vec, 'LineWidth', 1)
    xlim([n*3/40 n/10])
    hold on
    grid on
end
hold off

% sim('Model/Model_2022b.slx', 600);
% 
% figure
% subplot(2,1,1)
% plot(time, y, 'LineWidth', 1.5, 'Color', 'b')
% xlim([0 n/10])
% hold on
% plot(output.Time, output.Data, 'LineWidth', 1.5, 'Color', 'r')
% grid on
% hold off
% subplot(2,1,2)
% stairs(time_u, u_vec, 'LineWidth', 1.5, 'Color', 'b')
% hold on
% stairs(cv.Time, cv.Data, 'LineWidth', 1.5, 'Color', 'r')
% hold off
% grid on

function [x] = RungeKutta(x, x0, OR, Ph, tau)
    h = 0.1;

    k11 = h*StateEquationObject(x(3), x0, OR, 12, 1.6, Ph, 12000, 4200, 1);
    k21 = h*StateEquationObject(x(3) + 0.5*k11, x0, OR, 12, 1.6, Ph, 12000, 4200, 1);
    k31 = h*StateEquationObject(x(3) + 0.5*k21, x0, OR, 12, 1.6, Ph, 12000, 4200, 1);
    k41 = h*StateEquationObject(x(3) + k31, x0, OR, 12, 1.6, Ph, 12000, 4200, 1);

    x(3) = x(3) + (k11 + 2*k21 + 2*k31 + k41)/6;

    G = -0.0002347*Ph + 1.012;

    k12 = h*StateEquationSensor(x(1:2), x(3), G, tau);
    k22 = h*StateEquationSensor(x(1:2) + 0.5*k12, x(3), G, tau);
    k32 = h*StateEquationSensor(x(1:2) + 0.5*k22, x(3), G, tau);
    k42 = h*StateEquationSensor(x(1:2) + k32, x(3), G, tau);

    x(1:2) = x(1:2) + (k12 + 2*k22 + 2*k32 + k42)/6;
end

function [dxdt] = StateEquationObject(x, x0, OR, F_max, V, Ph, Pnom, cw, ro)
    dxdt = -OR*F_max/(100*60*V)*x + OR*F_max/(100*60*V)*x0 + (Ph*Pnom)/(100*cw*ro*V);
end

function [dxdt] = StateEquationSensor(x, u, G, tau)
    dxdt(1) = x(2);
    dxdt(2) = -2/tau*x(2) - 1/tau^2*x(1) + G*u/tau^2;
end

