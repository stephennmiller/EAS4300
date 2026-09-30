%% EAS 4300 HW #7 Question 2
% Stephen Miller
%
% Created 3/28/22.
%
% Adding an afterburner to question 1. The maximum stagnation temperature
% downstream of the afterburner is 2000 K. The afterburner combustion
% efficiency eta_ab is 0.95 and the total pressure ratio rab is 0.97.
% All other efficiencies remain the same. Use a gamma value of 1.4 up to
% the primary burner, and a value of 1.3 for the rest of the engine.
% Assume R is 0.287 kJ/kg*K throughout the engine. Plot the specific
% thrust, TSFC, eta_th, eta_p, and eta_o as a function of rc, the total
% pressure ratio across the compressor. Consider a range of rc from 2 to
% 60. Assume the exhaust is ideally expanded. Also plot the nozzle area
% ratio as a function of rc. Comment on the changes that occur due to the
% addition of the afterburner (i.e., compare the plots to that of problem
% 1). Note that the overall fuel-to-air ratio cannot exceed the
% stoichiometric value (fb + fab <= fst).

clc; clear; close all;

%% Given Parameters
M             = 1.8;
T04_max       = 1500;       % [K], max turbine inlet temperature
T06_ab_max    = 2000;       % [K], max afterburner stagnation temperature
delta_h_c     = 43124;      % [kJ/kg], fuel lower heating value
eta_d         = 0.9;        % diffuser efficiency
eta_c         = 0.9;        % compressor efficiency
eta_b         = 0.98;       % main burner efficiency
eta_t         = 0.92;       % turbine efficiency
eta_n         = 0.98;       % nozzle efficiency
eta_ab        = 0.95;       % afterburner combustion efficiency
rb            = 0.97;       % burner pressure ratio (P04 / P03)
rab           = 0.97;       % afterburner pressure ratio (P06 / P05)
f_st          = 0.06;       % stoichiometric fuel-to-air ratio
k             = [1.4, 1.3]; % gamma values: 1.4 (up to burner), 1.3 (rest)
R             = 0.287;      % [kJ/kg·K], gas constant
Ta            = 216.65;     % [K], ambient temperature at 15 km
Pa            = 1.2112e4;   % [Pa], ambient pressure at 15 km
delta         = 100;        % number of data points

%% Preallocate
rc       = linspace(2, 60, delta)';  % compressor pressure ratio range
n        = length(rc);

% Non-afterburner arrays
f_b      = zeros(n,1);
T04      = zeros(n,1);
T05      = zeros(n,1);
P05      = zeros(n,1);
I        = zeros(n,1);
TSFC     = zeros(n,1);
eta_th   = zeros(n,1);
eta_p    = zeros(n,1);
eta_0    = zeros(n,1);
Ratio    = zeros(n,1);

% Afterburner arrays
f_ab     = zeros(n,1);
f        = zeros(n,1);  % total fuel (main + afterburner)
T06_ab   = zeros(n,1);
I_ab     = zeros(n,1);
TSFC_ab  = zeros(n,1);
eta_p_ab = zeros(n,1);
eta_th_ab= zeros(n,1);
eta_0_ab = zeros(n,1);
Ratio_ab = zeros(n,1);

%% Inlet / Freestream
% Stagnation temperature and flight velocity:
%
% $$T_{0a} = T_a\left(1 + \frac{\gamma_1 - 1}{2}M^2\right), \qquad u = M\sqrt{\gamma_1 R T_a}$$
%
% The diffuser is adiabatic, so $T_{02} = T_{0a}$. Its efficiency sets the
% pressure recovery:
%
% $$T_{02s} = T_a + \eta_d\,(T_{0a} - T_a), \qquad P_{02} = P_a\left(\frac{T_{02s}}{T_a}\right)^{\gamma_1/(\gamma_1 - 1)}$$
%
% Specific heats: $c_{p1} = \frac{\gamma_1}{\gamma_1 - 1}R$ for air ($\gamma_1 = 1.4$)
% and $c_{p2} = \frac{\gamma_2}{\gamma_2 - 1}R$ for burned gas ($\gamma_2 = 1.3$).
T0a = Ta * (1 + ((k(1) - 1)/2) * M^2);
u   = M * sqrt(R * 1000 * k(1) * Ta);

% Ideal diffuser (for T02s), then real P02
T02s = eta_d * (T0a - Ta) + Ta;
P02  = Pa * (T02s / Ta)^(k(1)/(k(1) - 1));

% Specific heats [kJ/kg·K]: air (k = 1.4) and burned gas (k = 1.3)
Cp = (k ./ (k - 1)) * R;

%% Turbojet Cycle With and Without Afterburner
% *Compressor* (efficiency $\eta_c$):
%
% $$T_{03s} = T_{02}\,r_c^{(\gamma_1-1)/\gamma_1}, \qquad T_{03} = T_{02} + \frac{T_{03s} - T_{02}}{\eta_c}, \qquad w_c = c_{p1}(T_{03} - T_{02})$$
%
% *Main burner* energy balance, with the stoichiometric limit $f_b \le f_{st}$:
%
% $$c_{p1}T_{03} + f_b\,\eta_b\,\Delta h_c = (1 + f_b)\,c_{p2}T_{04} \quad\Rightarrow\quad f_b = \frac{c_{p2}T_{04} - c_{p1}T_{03}}{\eta_b\,\Delta h_c - c_{p2}T_{04}}, \qquad P_{04} = r_b\,r_c\,P_{02}$$
%
% If $f_b$ would exceed $f_{st}$, the burner runs at $f_{st}$ and $T_{04}$ is
% solved from the same balance.
%
% *Turbine* work drives the compressor (efficiency $\eta_t$):
%
% $$(1 + f_b)\,c_{p2}(T_{04} - T_{05}) = w_c, \qquad T_{05s} = T_{04} - \frac{T_{04} - T_{05}}{\eta_t}$$
%
% $$P_{05} = P_{04}\left(\frac{T_{05s}}{T_{04}}\right)^{\gamma_2/(\gamma_2-1)}$$
%
% *Afterburner* energy balance, reheating the turbine exhaust to $T_{06}$:
%
% $$(1 + f_b)\,c_{p2}T_{05} + f_{ab}\,\eta_{ab}\,\Delta h_c = (1 + f_b + f_{ab})\,c_{p2}T_{06}$$
%
% $$f_{ab} = \frac{(1 + f_b)(T_{06} - T_{05})}{\eta_{ab}\,\Delta h_c / c_{p2} - T_{06}}, \qquad P_{06} = r_{ab}\,P_{05}$$
%
% If $f_b + f_{ab}$ would exceed $f_{st}$, the afterburner burns only
% $f_{ab} = f_{st} - f_b$ and $T_{06}$ is solved from the same balance.
%
% *Nozzle* (function |nozzle_exit|) expands ideally to $P_a$ from
% $(T_0, P_0) = (T_{05}, P_{05})$ without the afterburner, or
% $(T_{06}, P_{06})$ with it:
%
% $$T_{es} = T_0\left(\frac{P_a}{P_0}\right)^{(\gamma_2-1)/\gamma_2}, \qquad T_e = T_0 - \eta_n(T_0 - T_{es})$$
%
% $$M_e = \sqrt{\frac{2}{\gamma_2 - 1}\left(\frac{T_0}{T_e} - 1\right)}, \qquad u_e = M_e\sqrt{\gamma_2 R T_e}$$
%
% $$\frac{A_e}{A^*} = \frac{1}{M_e}\left[\frac{2}{\gamma_2+1}\left(1 + \frac{\gamma_2-1}{2}M_e^2\right)\right]^{\frac{\gamma_2+1}{2(\gamma_2-1)}}$$
%
% *Performance*, with $f = f_b$ without the afterburner and
% $f = f_b + f_{ab}$ with it:
%
% $$I = (1 + f)\,u_e - u, \qquad \mathrm{TSFC} = \frac{f}{I}$$
%
% $$\eta_{th} = \frac{(1 + f)\,u_e^2 - u^2}{2 f\,\Delta h_c}, \qquad \eta_p = \frac{I\,u}{\frac{1}{2}\left[(1 + f)\,u_e^2 - u^2\right]}, \qquad \eta_0 = \eta_{th}\,\eta_p$$
for j = 1:n
    % Compressor
    T03s = T0a * rc(j)^((k(1) - 1)/k(1));
    T03  = ((T03s - T0a) / eta_c) + T0a;
    W    = Cp(1) * (T03 - T0a);

    % Main burner energy balance: Cp(1)*T03 + f*eta_b*dhc = (1 + f)*Cp(2)*T04
    f_b(j) = (Cp(2)*T04_max - Cp(1)*T03) / (eta_b*delta_h_c - Cp(2)*T04_max);
    T04(j) = T04_max;
    if f_b(j) > f_st
        % Stoichiometric limit: burn at f_st and find the lower T04
        f_b(j) = f_st;
        T04(j) = (Cp(1)*T03 + f_st*eta_b*delta_h_c) / ((1 + f_st)*Cp(2));
    end
    P04 = rb * rc(j) * P02;

    % Turbine exit (T05 = T06 = T07 without the afterburner)
    T05(j) = T04(j) - (W/(Cp(2)*(1 + f_b(j))));
    T05s   = ((T05(j) - T04(j))/eta_t) + T04(j);
    P05(j) = P04 * (T05s/T04(j))^(k(2)/(k(2) - 1));

    % ---- Without afterburner: nozzle expansion to ambient ----
    [Ue, Me]  = nozzle_exit(T05(j), P05(j), Pa, eta_n, k(2), R);
    I(j)      = (1 + f_b(j))*Ue - u;
    TSFC(j)   = f_b(j) / I(j);
    eta_p(j)  = (I(j)*u) / ((1 + f_b(j))*(Ue^2)/2 - (u^2)/2);
    eta_th(j) = ((1 + f_b(j))*Ue^2 - u^2)/(2*f_b(j)*delta_h_c*1000);
    eta_0(j)  = eta_th(j) * eta_p(j);
    Ratio(j)  = area_ratio(Me, k(2));

    % ---- With afterburner ----
    % Energy balance: (1 + f_b)*Cp(2)*T05 + f_ab*eta_ab*dhc = (1 + f_b + f_ab)*Cp(2)*T06
    f_ab(j)   = ((1 + f_b(j))*(T06_ab_max - T05(j))) / ...
                (((eta_ab*delta_h_c)/Cp(2)) - T06_ab_max);
    T06_ab(j) = T06_ab_max;
    if f_b(j) + f_ab(j) > f_st
        % Overall fuel is capped at f_st: burn the rest and find the lower T06
        f_ab(j)   = f_st - f_b(j);
        T06_ab(j) = ((1 + f_b(j))*T05(j) + f_ab(j)*eta_ab*delta_h_c/Cp(2)) / (1 + f_st);
    end
    f(j) = f_b(j) + f_ab(j);

    % Afterburner exit and nozzle expansion to ambient
    P06_ab       = rab * P05(j);
    [Ue_ab, Me_ab] = nozzle_exit(T06_ab(j), P06_ab, Pa, eta_n, k(2), R);
    I_ab(j)      = (1 + f(j))*Ue_ab - u;
    TSFC_ab(j)   = f(j) / I_ab(j);
    eta_th_ab(j) = ((1 + f(j))*Ue_ab^2 - u^2) / (2*f(j)*delta_h_c*1000);
    eta_p_ab(j)  = (I_ab(j)*u) / ((1 + f(j))*(Ue_ab^2)/2 - (u^2)/2);
    eta_0_ab(j)  = eta_th_ab(j)*eta_p_ab(j);
    Ratio_ab(j)  = area_ratio(Me_ab, k(2));
end

%% Optimum rc
[I_max, iI]         = max(I);
[TSFC_min, iT]      = min(TSFC);
[I_ab_max, iIab]    = max(I_ab);
[TSFC_ab_min, iTab] = min(TSFC_ab);
fprintf('Without afterburner: max I = %.1f m/s at rc = %.2f; min TSFC = %.4e s/m at rc = %.2f\n', ...
    I_max, rc(iI), TSFC_min, rc(iT));
fprintf('With afterburner:    max I = %.1f m/s at rc = %.2f; min TSFC = %.4e s/m at rc = %.2f\n', ...
    I_ab_max, rc(iIab), TSFC_ab_min, rc(iTab));
fprintf('Largest total fuel-to-air ratio: %.4f (f_st = %.2f)\n', max(f), f_st);

%% Figure 1: Specific Thrust vs. rc
figure('Position', [100, 100, 600, 400]);
plot(rc, I, 'DisplayName', 'Non-afterburner', 'LineWidth', 1.5);
hold on
plot(rc, I_ab, '--', 'Color', 'b', 'DisplayName', 'Afterburner', 'LineWidth', 1.5);
hold off
ylim([300 1200]);
legend('Location','southwest');
xlabel('r_c [dim]', 'FontWeight','bold');
ylabel('I [m/s]', 'FontWeight','bold');
title('Specific Thrust vs. r_c');
grid on;

%% Figure 2: TSFC vs. rc
figure('Position', [120, 120, 600, 400]);
plot(rc, TSFC, 'DisplayName','Non-afterburner','LineWidth',1.5);
hold on
plot(rc, TSFC_ab, '--', 'Color','b','DisplayName','Afterburner','LineWidth',1.5);
hold off
legend('Location','northeast');
ax = gca; ax.YAxis.Exponent = 0;  % no scientific notation
xlabel('r_c [dim]', 'FontWeight','bold');
ylabel('TSFC [s/m]', 'FontWeight','bold');
title('TSFC vs. r_c');
grid on;

%% Figure 3: Efficiencies vs. rc
figure('Position', [140, 140, 600, 400]);
plot(rc, eta_th,  'Color','b','LineWidth',1.5, 'DisplayName','\eta_{th} Non-AB');
hold on
plot(rc, eta_0,   'Color','k','LineWidth',1.5, 'DisplayName','\eta_{o} Non-AB');
plot(rc, eta_p,   'Color','r','LineWidth',1.5, 'DisplayName','\eta_{p} Non-AB');
plot(rc, eta_th_ab, '--','Color','b','LineWidth',1.5,'DisplayName','\eta_{th} AB');
plot(rc, eta_0_ab,  '--','Color','k','LineWidth',1.5,'DisplayName','\eta_{o} AB');
plot(rc, eta_p_ab,  '--','Color','r','LineWidth',1.5,'DisplayName','\eta_{p} AB');
hold off
legend('Location','northwest','FontSize',9);
xlabel('r_c [dim]', 'FontWeight','bold');
ylabel('Efficiency', 'FontWeight','bold');
title('Efficiencies vs. r_c');
grid on;

%% Figure 4: Nozzle Area Ratio vs. rc
figure('Position', [160, 160, 600, 400]);
plot(rc, Ratio, 'DisplayName','Non-afterburner','LineWidth',1.5);
hold on
plot(rc, Ratio_ab, '--','Color','b','DisplayName','Afterburner','LineWidth',1.5);
hold off
legend('Location','northeast');
xlabel('r_c [dim]', 'FontWeight','bold');
ylabel('Area Ratio', 'FontWeight','bold');
title('Nozzle Area Ratio vs. r_c');
grid on;

%% Save Data to CSV
data = table( ...
    rc, ...
    f_b,       f_ab,      f,        T06_ab, ...
    I,         I_ab, ...
    TSFC,      TSFC_ab, ...
    eta_th,    eta_th_ab, ...
    eta_p,     eta_p_ab, ...
    eta_0,     eta_0_ab, ...
    Ratio,     Ratio_ab, ...
    'VariableNames', { ...
    'rc', ...
    'f_b_nonAB','f_ab_AB','f_total','T06_AB', ...
    'I_nonAB',  'I_AB', ...
    'TSFC_nonAB','TSFC_AB', ...
    'eta_th_nonAB','eta_th_AB', ...
    'eta_p_nonAB','eta_p_AB', ...
    'eta_0_nonAB','eta_0_AB', ...
    'AR_nonAB','AR_AB'});

writetable(data, 'HW7_Q2_Data.csv');
disp('Results saved to HW7_Q2_Data.csv');

%% Helper Functions
function [Ue, Me] = nozzle_exit(T0, P0, Pa, eta_n, k, R)
    % Expand from stagnation state (T0, P0) to ambient pressure Pa
    Ts = T0 / (P0/Pa)^((k - 1)/k);      % Isentropic exit temperature
    Te = eta_n * (Ts - T0) + T0;        % Actual exit temperature
    Me = sqrt(((T0/Te) - 1)/((k - 1)/2));
    Ue = Me * sqrt(k*R*1000*Te);        % Exit velocity [m/s]
end

function AR = area_ratio(Me, k)
    % Nozzle exit-to-throat area ratio for exit Mach number Me
    AR = (1/Me) * ((2/(k+1)) * (1 + ((k-1)/2)*Me^2))^((k+1)/(2*(k-1)));
end
