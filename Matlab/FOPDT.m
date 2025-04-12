clear all
n = 8000;

h = 0.1;

F_max = 12;
V = 1.6;
Pnom = 12000;
Ph = 50;
cw = 4200;
ro = 1;

T_in = 17;

G = -0.0002347*Ph + 1.012;

x = zeros(1,3);
x(1) = T_in;
x(3) = T_in;
y = zeros(n,1) + T_in;
u_vector = zeros(n,1) + T_in;

T_Setpoint_Begin = 35;
T_Setpoint_End = 45;

alpha = F_max/(100*60*V);
beta = (Ph*Pnom)/(100*cw*ro*V);

OR = 1/(alpha/(beta*G)*T_Setpoint_Begin - alpha/beta*T_in);
u_prev(1) = OR;

for i=1:1:n
    if i == 4000
        OR = 1/(alpha/(beta*G)*T_Setpoint_End - alpha/beta*T_in);
        u_prev(2) = OR;
    end
    F = OR*F_max/100;
    tau = 19.08*F^(-0.4293) - 4.042;
    tau0 = 11.93*F^(-0.78) + 2.37;
    tau0 = round(tau0/h);
    x = RungeKutta(h, x, T_in, G, tau, F, V, Ph, Pnom, cw, ro);
    y(i + tau0:i + tau0 + 100) = x(1);
    u_vector(i) = F;
end

sys = tf([-0.7051], [35.1 1], 'InputDelay', 21);
t = 399.9:0.1:799.9;
stepresponse = step(sys, t, RespConfig(Amplitude = -14.1819));
stepresponse = stepresponse + 35;
t = (0:0.1:792.5)';
stepresponse = [zeros(3999,1) + 35; stepresponse];

[size_u, ~] = size(u_vector);
y = y(1:size_u);
time = (0:h:(size_u-1)*h)';

y = y(3999:end);
time = (0:0.1:(length(y)-1)/10)';
stepresponse = stepresponse(3999:end);

plot(time, y, 'LineWidth', 1.5, 'color', 'r');
xlim([0 max(time)]);
% ylim([0.9*T0 1.05*max(y)])
hold on
% plot(time, u_vector)
plot(time, stepresponse, 'LineWidth', 1.5, 'color', 'b');
yl = sprintf("T, %cC", char(176));
ylabel(yl);
xlabel("t, s");
legend('Postać równań różniczkowych', 'Postać przybliżenia FOPDT');
grid on

%% Utilities
function [x] = RungeKutta(h, x, x0, G, tau, F, V, Ph, Pnom, cw, ro)
    k11 = h*StateEquationObject(x(3), x0, F, V, Ph, Pnom, cw, ro);
    k21 = h*StateEquationObject(x(3) + 0.5*k11, x0, F, V, Ph, Pnom, cw, ro);
    k31 = h*StateEquationObject(x(3) + 0.5*k21, x0, F, V, Ph, Pnom, cw, ro);
    k41 = h*StateEquationObject(x(3) + k31, x0, F, V, Ph, Pnom, cw, ro);

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