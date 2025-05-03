%% Acquire step response for all operating points
clear all
close all
clc

n = 10000;

T_in = 15;
F_max = 12;
V = 1.6;
Pnom = 12000;
cw = 4200;
ro = 1;

x = zeros(1,3);
x(1) = T_in;
x(3) = T_in;
y = zeros(2*n,1) + T_in;

% Ph_array = [30 40 50 60 70 80 90 100];
% T_OP_array = [30 40 50 60 70 80];
Ph_array = 80;
T_OP_array = 50;

u = zeros(length(T_OP_array), length(Ph_array));

for k=1:length(Ph_array)
    for j=1:length(T_OP_array)
        Ph = Ph_array(k);
        G = -0.0002347*Ph + 1.012;
        alpha = F_max/(100*60*V);
        beta = (Ph*Pnom)/(100*cw*ro*V);
        y = [];
    
        T_OP = T_OP_array(j);
        OR = 1/(alpha/(beta*G)*T_OP - alpha/beta*T_in);
        for i=1:2*n
            if i < n
                OR = 1/(alpha/(beta*G)*(T_OP-10) - alpha/beta*T_in);
                u_norm(1) = OR;
            else
                OR = 1/(alpha/(beta*G)*(T_OP+10) - alpha/beta*T_in);
                u_norm(2) = OR;
            end
            F = OR*F_max/100;
            tau = 19.08*F^-0.4293 - 4.042;
            tau0 = 11.93*F^-0.78 + 2.37;
            x = RungeKutta(x, T_in, OR, Ph, tau);
            y(i + round(tau0*10):i+round(tau0*10)+n/10) = x(1);
        end
    
        Tp = 1.4;
        s_raw = (y(n:end)-y(n))/(u_norm(2)-u_norm(1));
        
        step_response(j,k,:) = s_raw(1:1.1*n);

        if j > 1 || k > 1
            s_temp = s_raw(1:Tp*10:end);
            if length(s_temp) < length(s(1,1,:))
                s_temp(end:length(s(1,1,:))) = s_temp(end);
            end
            s_temp = s_temp(1:length(s(1,1,:)));
            s(j,k,:) = s_temp;
        else
            s(j,k,:) = s_raw(1:Tp*10:end);
        end
    end
end

%% Calculating FOPDT model parameters for all operating points [k, T, T0]

FOPDT_params = zeros(length(T_OP_array), length(Ph_array), 3);
FOPDT_params(:,:,1) = s(:,:,end);

for k=1:length(Ph_array)
    for j=1:length(T_OP_array)
        t1 = find(abs(step_response(j,k,:)) >= abs(0.283*step_response(j,k,end)), 1)/10;
        t2 = find(abs(step_response(j,k,:)) >= abs(0.632*step_response(j,k,end)), 1)/10;
        FOPDT_params(j,k,2) = 1.5*(t2-t1);
        FOPDT_params(j,k,3) = t2 - FOPDT_params(j,k,2);
    end
end

%% Calculating sampling time Tp

T_temp = FOPDT_params(:,:,2);
Tp = 0.1*floor(min(T_temp(:)));

%% Calculating DMC parameters for all operating points [D, N, Nu, N1, lambda]

DMC_params = zeros(length(T_OP_array), length(Ph_array), 5);

DMC_params(:,:,1) = max(max(round(3*FOPDT_params(:,:,2)/Tp + FOPDT_params(:,:,3)/Tp)));
DMC_params(:,:,2) = max(max(round(FOPDT_params(:,:,2)/Tp + FOPDT_params(:,:,3)/Tp)));
DMC_params(:,:,3) = 2;
DMC_params(:,:,4) = round(FOPDT_params(:,:,3)/Tp);
DMC_params(:,:,5) = 0.005*FOPDT_params(:,:,1).^2 .* DMC_params(:,:,2);

%% Calculating DMC parameters for all operating points

Ke_params = zeros(length(T_OP_array), length(Ph_array));
Ku_params = zeros(length(T_OP_array), length(Ph_array), DMC_params(1,1,1));

for k=1:length(Ph_array)
    for j=1:length(T_OP_array)

        D = DMC_params(j,k,1);
        N = DMC_params(j,k,2);
        Nu = DMC_params(j,k,3);
        N1 = DMC_params(j,k,4);
        l = DMC_params(j,k,5);
        
        MP = zeros(N-N1+1, D);
        
        for i=1:N-N1+1
            for m=1:D
                if i+m+N1-1 < D
                    MP(i,m) = s(j,k,i+m+N1-1) - s(j,k,m);
                else
                    MP(i,m) = s(j,k,D) - s(j,k,m);
                end
            end
        end
        
        M = zeros(N-N1+1, Nu);
        
        for i=1:Nu
            for m=i:N-N1+1
                M(m,i) = s(j,k,m+N1-i);
            end
        end
        
        K = ((M'*M+l*eye(Nu))^-1)*M';
        
        Ku_params(j,k,:) = K(1,:)*MP;
        Ke_params(j,k) = sum(K(1,:));
    end
end

%% Adaptive DMC test

n = 20000;

F_max = 12;
T0 = 15;

x = zeros(1,3);
x(1) = T0;
x(3) = T0;
y = zeros(n,1) + T0;
u_vec = zeros(n,1);

Ph = 80;
T_OP = 45;
T_OP_prev = T_OP;
T_OP_next = 50;

Ph_idx = find(Ph_array >= Ph, 1);
T_OP_idx = find(T_OP_array >= T_OP_next, 1);

% Ke = (Ke_params(T_OP_idx,Ph_idx)-Ke_params(T_OP_idx,Ph_idx-1))/(Ph_array(Ph_idx)-Ph_array(Ph_idx-1))*Ph + ...
%     (Ke_params(T_OP_idx,Ph_idx-1)*Ph_array(Ph_idx)-Ke_params(T_OP_idx,Ph_idx)*Ph_array(Ph_idx-1))/(Ph_array(Ph_idx)-Ph_array(Ph_idx-1));
% 
% Ku = zeros(length(Ku_params(1,1,:)), 1);
% for k=1:length(Ku_params(1,1,:))
%     Ku(k) = (Ku_params(T_OP_idx,Ph_idx,k)-Ku_params(T_OP_idx,Ph_idx-1,k))/(Ph_array(Ph_idx)-Ph_array(Ph_idx-1))*Ph + ...
%     (Ku_params(T_OP_idx,Ph_idx-1,k)*Ph_array(Ph_idx)-Ku_params(T_OP_idx,Ph_idx,k)*Ph_array(Ph_idx-1))/(Ph_array(Ph_idx)-Ph_array(Ph_idx-1));
% end

x1 = Ph_array(Ph_idx-1);
x2 = Ph_array(Ph_idx);
y1 = T_OP_array(T_OP_idx-1);
y2 = T_OP_array(T_OP_idx);
Q11 = Ke_params(T_OP_idx-1,Ph_idx-1);
Q12 = Ke_params(T_OP_idx-1,Ph_idx);
Q21 = Ke_params(T_OP_idx,Ph_idx-1);
Q22 = Ke_params(T_OP_idx,Ph_idx);
Ke = 1/((x2-x1)*(y2-y1)) * (Q22*abs(x1-Ph)*abs(y1-T_OP_next) + Q21*abs(x2-Ph)*abs(y1-T_OP_next) + Q12*abs(x1-Ph)*abs(y2-T_OP_next) + Q11*abs(x2-Ph)*abs(y2-T_OP_next));

Ku = zeros(length(Ku_params(1,1,:)), 1);
for k=1:length(Ku_params(1,1,:))
    Q11 = Ku_params(T_OP_idx-1,Ph_idx-1,k);
    Q12 = Ku_params(T_OP_idx-1,Ph_idx,k);
    Q21 = Ku_params(T_OP_idx,Ph_idx-1,k);
    Q22 = Ku_params(T_OP_idx,Ph_idx,k);
    Ku(k) = 1/((x2-x1)*(y2-y1)) * (Q22*abs(x1-Ph)*abs(y1-T_OP_next) + Q21*abs(x2-Ph)*abs(y1-T_OP_next) + Q12*abs(x1-Ph)*abs(y2-T_OP_next) + Q11*abs(x2-Ph)*abs(y2-T_OP_next));
end

u_subtruct = 0;
OR = 50;
dOR = 0;
dOR_prev = zeros(D-1,1);
k = 0;
for i=1:n
    if k == round(10*Tp)
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
            OR = 10;
        end

        for j=D-1:-1:2
            dOR_prev(j) = dOR_prev(j-1);
        end
        dOR_prev(1) = dOR;
        k = 0;
    end

    if i == 15000
        T_OP = T_OP_next;

        dOR_prev(:) = 0;

        % T_OP_idx = find(T_OP_array >= T_OP, 1);
        % Ke = (Ke_params(T_OP_idx,Ph_idx)-Ke_params(T_OP_idx-1,Ph_idx))/(T_OP_array(T_OP_idx)-T_OP_array(T_OP_idx-1))*T_OP + ...
        %     (Ke_params(T_OP_idx-1,Ph_idx)*T_OP_array(T_OP_idx)-Ke_params(T_OP_idx,Ph_idx)*T_OP_array(T_OP_idx-1))/(T_OP_array(T_OP_idx)-T_OP_array(T_OP_idx-1));
        % 
        % Ku = zeros(length(Ku_params(1,1,:)), 1);
        % for j=1:length(Ku_params(1,1,:))
        %     Ku(j) = (Ku_params(T_OP_idx,Ph_idx,j)-Ku_params(T_OP_idx-1,Ph_idx,j))/(T_OP_array(T_OP_idx)-T_OP_array(T_OP_idx-1))*T_OP + ...
        %     (Ku_params(T_OP_idx-1,Ph_idx,j)*T_OP_array(T_OP_idx)-Ku_params(T_OP_idx,Ph_idx,j)*T_OP_array(T_OP_idx-1))/(T_OP_array(T_OP_idx)-T_OP_array(T_OP_idx-1));
        % end
    end
    
    F = OR*F_max/100;
    tau = 19.08*F^-0.4293 - 4.042;
    tau0 = 11.93*F^-0.78 + 2.73;

    x = RungeKutta(x, T0, OR, Ph, tau);
    y(i + round(tau0*10):i+round(tau0*10)+30) = x(1);
    u_vec(i) = OR;
    k = k + 1;
end

[size_y, ~] = size(y);
time = (0:0.1:(size_y-15001+50)/10)';
time_u = (0:0.1:(n-15001+50)/10);

T_OP_vec = zeros(length(time),1);
T_OP_vec(1:50) = T_OP_prev;
T_OP_vec(51:end) = T_OP_next;
Ph_vec = zeros(length(time_u),1) + Ph;
F_vec = u_vec(15000-49:end)*F_max/100;

f = figure;
f.Position = [100 100 800 600];

subplot(2,1,1)
stairs(time, T_OP_vec, '--', 'LineWidth', 1.3, 'Color', '#FF7F00');
yl = sprintf("T, %cC", char(176));
ylabel(yl);
xlabel("t, s")
hold on
plot(time, y(15000-49:end), 'LineWidth', 1.5, 'Color', 'b')
xlim([0 400]);
grid on
subplot(2,1,2)
yyaxis left
stairs(time_u, F_vec, 'LineWidth', 1.5, 'Color', '#7E2F8E');
xlabel("t, s")
ylabel("F, l/min")
hold on
yyaxis right
stairs(time_u, Ph_vec, 'LineWidth', 1, 'Color', '#34831B');
ylim([40 60])
xa = gca;
xa.YAxis(1).Color = [0 0 0];
xa.YAxis(2).Color = '#34831B';
ylabel("P_h, %")
xlim([0 400]);
grid on

%% Whole experiment plot
[size_y, ~] = size(y);
time = (0:0.1:(size_y-1)/10)';
time_u = (0:0.1:(n-1)/10);

f = figure;
f.Position = [100 100 800 600];

subplot(2,1,1)
stairs(time, T_OP_vec, '--', 'LineWidth', 1.3, 'Color', '#FF7F00');
xlim([3*n/40-5 1850])
yl = sprintf("T, %cC", char(176));
ylabel(yl);
xlabel("t, s")
hold on
plot(time, y, 'LineWidth', 1.5, 'Color', 'b')
grid on
subplot(2,1,2)
yyaxis left
stairs(time_u, F_vec, 'LineWidth', 1.5, 'Color', '#7E2F8E');
xlim([3*n/40-5 1850])
xlabel("t, s")
ylabel("F, l/min")
hold on
yyaxis right
stairs(time_u, Ph_vec, 'LineWidth', 1, 'Color', '#34831B');
ylim([40 60])
xa = gca;
xa.YAxis(1).Color = [0 0 0];
xa.YAxis(2).Color = '#34831B';
ylabel("P_h, %")
grid on

%% Utilities

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

