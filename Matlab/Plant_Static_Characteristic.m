%% 2D - Full Flow Range and specific Heat
clear all
close all
clc

% Zakres F i stałe P = 50
F = linspace(2, 12, 1000);
P = 50 * ones(size(F));

% Parametry
Tin = 15;
Pnom = 12000;
cw = 4200;
ro = 1;

% Funkcja nieliniowa T0
T0_func = @(P, Tin, F, Pnom, cw, ro) -0.0002347*P.*Tin + 1.012*Tin + ...
    (3*Pnom)./(5*cw*ro).*(-0.0002347*P.^2./F + 1.012*P./F);

% Obliczanie T0 dla P=50 i różnych F
T0_values = T0_func(P, Tin, F, Pnom, cw, ro);

% Znalezienie punktu, w którym T0 = 50
F_symbolic = sym('F');
P_fixed = 50;
T0_sym = -0.0002347*P_fixed*Tin + 1.012*Tin + ...
    (3*Pnom)/(5*cw*ro)*(-0.0002347*P_fixed^2/F_symbolic + 1.012*P_fixed/F_symbolic);

% Rozwiązanie T0 = 50
F0 = double(vpasolve(T0_sym == 40, F_symbolic, [2 12]));

% Obliczenie T0 i pochodnej w tym punkcie
T0_at_F0 = double(subs(T0_sym, F_symbolic, F0));
dT0_dF = diff(T0_sym, F_symbolic);
dT0_dF_at_F0 = double(subs(dT0_dF, F_symbolic, F0));

% Równanie stycznej
T0_tangent = T0_at_F0 + dT0_dF_at_F0 * (F(1:500) - F0);

% Wykres
figure;
plot(F, T0_values, 'b-', 'LineWidth', 2); hold on;
plot(F(1:500), T0_tangent, 'r--', 'LineWidth', 2);
plot(F0, T0_at_F0, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
xlabel('F_0, l/min');
yl = sprintf("T_0, %cC", char(176));
ylabel(yl);
title('Charakterystyka statyczna obiektu dla zadanej mocy podgrzewacza P_h');
legend('Temperatura wyjściowa', 'Temperatura wyjściowa po linearyzacji w punkcie pracy');
xlim([2 max(F)])
grid on;

%% 3D - Full Power and Flow Range
clear all
close all
clc

[F, P] = meshgrid(0.5:0.1:10, 0:0.1:100);

% Parametry
Tin = 15; % Parametr Tin
Pnom = 12000; % Parametr Pnom
cw = 4200; % Parametr cw
ro = 1; % Parametr ro
T0_target = 50; % Punkt pracy (T0 = 50)
F_target = 3; % Punkt pracy (F = 3)

% Funkcja nieliniowa T0
T0 = @(P, Tin, F, Pnom, cw, ro) -0.0002347*P*Tin + 1.012*Tin + (3*Pnom)/(5*cw*ro)*(-0.0002347*P.^2./F + 1.012*P./F);

% Funkcja do obliczenia P w punkcie pracy
f_P = @(P) T0(P, Tin, F_target, Pnom, cw, ro) - T0_target; % Przyjmujemy Tin = 15

% Używamy fsolve do rozwiązania równania T0 = 50 dla P
P_target = fsolve(f_P, 10); % Rozwiązywanie równania dla P, przybliżenie początkowe P = 10

% Obliczanie wartości funkcji nieliniowej w punkcie pracy
T0_at_point = T0(P_target, Tin, F_target, Pnom, cw, ro); % T0 w punkcie pracy

% Obliczanie pochodnych cząstkowych funkcji T0 względem P i F
syms P F
T0_symbolic = -0.0002347*P*Tin + 1.012*Tin + (3*Pnom)/(5*cw*ro)*(-0.0002347*P.^2./F + 1.012*P./F);

% Pochodna cząstkowa względem P
dT0_dP = diff(T0_symbolic, P); 

% Pochodna cząstkowa względem F
dT0_dF = diff(T0_symbolic, F); 

% Obliczanie pochodnych w punkcie pracy (P_target, F_target)
dT0_dP_at_point = double(subs(dT0_dP, {P, F}, {P_target, F_target})); % Pochodna względem P w punkcie pracy
dT0_dF_at_point = double(subs(dT0_dF, {P, F}, {P_target, F_target})); % Pochodna względem F w punkcie pracy

% Funkcja zlinearyzowana wokół punktu (P_target, F_target)
T0_lin = T0_at_point + dT0_dP_at_point * (P - P_target) + dT0_dF_at_point * (F - F_target);

% Tworzenie siatki dla funkcji nieliniowej (pełny zakres dla P i F)
[P_grid_full, F_grid_full] = meshgrid(linspace(0, 100, 100), linspace(0.5, 12, 100));
T0_full = double(subs(T0_symbolic, {P, F}, {P_grid_full, F_grid_full}));

% Tworzenie siatki dla funkcji zlinearyzowanej (w zakresie linearyzacji)
P_range_lin = linspace(P_target - 10, P_target + 10, 30);
F_range_lin = linspace(F_target - 1, F_target + 1, 30);
[P_grid_lin, F_grid_lin] = meshgrid(P_range_lin, F_range_lin);
T0_lin_grid = double(subs(T0_lin, {P, F}, {P_grid_lin, F_grid_lin}));

% Wykres 3D
figure;
hold on;

% Funkcja nieliniowa na pełnym zakresie
surf(P_grid_full, F_grid_full, T0_full, 'FaceAlpha', 0.5); % Powierzchnia funkcji nieliniowej (przezroczysta)

% % Funkcja zlinearyzowana w wąskim zakresie
% surf(P_grid_lin, F_grid_lin, T0_lin_grid, 'FaceColor', 'r'); % Powierzchnia funkcji zlinearyzowanej

xlabel('P_0, %');
ylabel('F_0, l/min');
zl = sprintf("T_0, %cC", char(176));
zlabel(zl);
title('Charakterystyka statyczna obiektu cieplnego');
% legend({'Funkcja nieliniowa', 'Funkcja zlinearyzowana'}, 'Location', 'best');
shading interp;
grid on;
hold off;
