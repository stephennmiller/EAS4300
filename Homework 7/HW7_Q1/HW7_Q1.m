%% EAS 4300 HW #7 Question 1
% Stephen Miller
%
% Created 3/28/22.
%
% A non-afterburning turbojet is being designed for operation at an
% altitude of 15 km and a Mach number of 1.8. The maximum stagnation
% temperature at the inlet of the turbine is 1500 K. The fuel is a jet fuel
% having a LHV of 43124 kJ/kg and fst is 0.06. The following efficiencies
% apply at this Mach number: eta_d = 0.9, eta_c = 0.9, eta_b = 0.98,
% rb = 0.97, eta_t = 0.92, eta_n = 0.98. Use a gamma value of 1.4 up to the
% burner, and a value of 1.3 for the rest of the engine. Assume R is
% 0.287 kJ/kgK throughout the engine. Plot the specific thrust, TSFC,
% eta_th, eta_p, and eta_o as a function of rc, the total pressure ratio
% across the compressor. Is there an optimum rc that minimizes TSFC? Is
% there an optimum rc that maximizes specific thrust? Consider a range of
% rc from 2 to 60. Assume the exhaust is ideally expanded. Also plot the
% nozzle area ratio as a function of rc.

clc; clear; close all;

%% Given Parameters
M        = 1.8;                         % Mach number
T04_max  = 1500;                        % Maximum turbine inlet temperature [K]
delta_h_c= 43124;                       % Fuel lower heating value [kJ/kg]
f_st     = 0.06;                        % Stoichiometric fuel-to-air ratio
eta_d    = 0.9;                         % Diffuser efficiency
eta_c    = 0.9;                         % Compressor efficiency
eta_b    = 0.98;                        % Burner efficiency
eta_t    = 0.92;                        % Turbine efficiency
eta_n    = 0.98;                        % Nozzle efficiency
rb       = 0.97;                        % Burner pressure ratio (P04/P03)
k        = [1.4 1.3];                   % Specific heat ratios (before/after burner)
R        = 0.287;                       % Gas constant [kJ/kg·K]
Ta       = 216.65;                      % Ambient temperature at 15 km [K]
Pa       = 1.2112e4;                    % Ambient pressure at 15 km [Pa]
delta    = 58;                          % Number of data points for rc
rc       = linspace(2, 60, delta)';     % Compressor pressure ratio range

%% Preallocate Variables
n        = length(rc);
T03s     = zeros(n, 1);
T03      = zeros(n, 1);
W        = zeros(n, 1);
f_b      = zeros(n, 1);
T04      = zeros(n, 1);
P04      = zeros(n, 1);
T05      = zeros(n, 1);
T05s     = zeros(n, 1);
P05      = zeros(n, 1);
T7s      = zeros(n, 1);
T7       = zeros(n, 1);
Me       = zeros(n, 1);
Ue       = zeros(n, 1);
I        = zeros(n, 1);
TSFC     = zeros(n, 1);
eta_th   = zeros(n, 1);
eta_p    = zeros(n, 1);
eta_0    = zeros(n, 1);
Ratio    = zeros(n, 1);

%% Inlet and Diffuser
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
T0a      = Ta * (1 + ((k(1) - 1)/2) * M^2);                    % T0a = T02
u        = M * sqrt(R * 1000 * k(1) * Ta);                     % Flight velocity [m/s]
T02s     = eta_d * (T0a - Ta) + Ta;                            % Isentropic-equivalent T02
P02      = Pa * (T02s / Ta)^(k(1) / (k(1) - 1));                % Diffuser exit pressure

% Specific heats [kJ/kg·K]: air (k = 1.4) and burned gas (k = 1.3)
Cp = (k ./ (k - 1)) * R;

%% Turbojet Cycle for Each rc
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
% *Nozzle* expands ideally to $P_a$ (efficiency $\eta_n$), with
% $(T_0, P_0) = (T_{05}, P_{05})$ and $T_e = T_7$:
%
% $$T_{es} = T_0\left(\frac{P_a}{P_0}\right)^{(\gamma_2-1)/\gamma_2}, \qquad T_e = T_0 - \eta_n(T_0 - T_{es})$$
%
% $$M_e = \sqrt{\frac{2}{\gamma_2 - 1}\left(\frac{T_0}{T_e} - 1\right)}, \qquad u_e = M_e\sqrt{\gamma_2 R T_e}$$
%
% $$\frac{A_e}{A^*} = \frac{1}{M_e}\left[\frac{2}{\gamma_2+1}\left(1 + \frac{\gamma_2-1}{2}M_e^2\right)\right]^{\frac{\gamma_2+1}{2(\gamma_2-1)}}$$
%
% *Performance*, with $f = f_b$:
%
% $$I = (1 + f)\,u_e - u, \qquad \mathrm{TSFC} = \frac{f}{I}$$
%
% $$\eta_{th} = \frac{(1 + f)\,u_e^2 - u^2}{2 f\,\Delta h_c}, \qquad \eta_p = \frac{I\,u}{\frac{1}{2}\left[(1 + f)\,u_e^2 - u^2\right]}, \qquad \eta_0 = \eta_{th}\,\eta_p$$
for j = 1:n
    % Compressor
    T03s(j) = T0a * rc(j)^((k(1) - 1) / k(1));          % Ideal compressor exit temperature
    T03(j)  = ((T03s(j) - T0a) / eta_c) + T0a;          % Actual compressor exit temperature
    W(j)    = Cp(1) * (T03(j) - T0a);                   % Compressor work [kJ/kg]

    % Burner energy balance: Cp(1)*T03 + f*eta_b*dhc = (1 + f)*Cp(2)*T04
    f_b(j)  = (Cp(2) * T04_max - Cp(1) * T03(j)) / (eta_b * delta_h_c - Cp(2) * T04_max);
    T04(j)  = T04_max;
    if f_b(j) > f_st
        % Stoichiometric limit: burn at f_st and find the lower T04
        f_b(j) = f_st;
        T04(j) = (Cp(1) * T03(j) + f_st * eta_b * delta_h_c) / ((1 + f_st) * Cp(2));
    end
    P04(j)  = rb * rc(j) * P02;

    % Turbine drives the compressor (T05 = T06 = T07)
    T05(j)  = T04(j) - (W(j) / (Cp(2) * (1 + f_b(j))));
    T05s(j) = ((T05(j) - T04(j)) / eta_t) + T04(j);
    P05(j)  = P04(j) * (T05s(j) / T04(j))^(k(2) / (k(2) - 1)); % P05 = P06

    % Nozzle: ideal expansion to Pa
    T7s(j)  = T05(j) / (P05(j) / Pa)^((k(2) - 1) / k(2));
    T7(j)   = eta_n * (T7s(j) - T05(j)) + T05(j);
    Me(j)   = sqrt(((T05(j) / T7(j)) - 1) / ((k(2) - 1) / 2));
    Ue(j)   = Me(j) * sqrt(k(2) * R * 1000 * T7(j));    % Exhaust velocity [m/s]

    % Performance
    I(j)      = (1 + f_b(j)) * Ue(j) - u;               % Specific thrust [m/s]
    TSFC(j)   = f_b(j) / I(j);                          % Thrust-specific fuel consumption [s/m]
    eta_p(j)  = (I(j) * u) / ((1 + f_b(j)) * (Ue(j)^2) / 2 - (u^2) / 2);
    eta_th(j) = ((1 + f_b(j)) * Ue(j)^2 - u^2) / (2 * f_b(j) * delta_h_c * 1000);
    eta_0(j)  = eta_th(j) * eta_p(j);                   % Overall efficiency
    Ratio(j)  = (1 / Me(j)) * ((2 / (k(2) + 1)) * (1 + ((k(2) - 1) / 2) * Me(j)^2))^((k(2) + 1) / (2 * (k(2) - 1)));
end

%% Optimum rc
[I_max, iI]      = max(I);
[TSFC_min, iT]   = min(TSFC);
fprintf('Maximum specific thrust: I = %.1f m/s at rc = %.2f\n', I_max, rc(iI));
fprintf('Minimum TSFC: %.4e s/m at rc = %.2f\n', TSFC_min, rc(iT));

%% Figure 1: Specific Thrust vs. rc
figure('Position', [100, 100, 600, 400]);
plot(rc, I, 'LineWidth', 1.5);
hold on;
plot(rc(iI), I_max, 'o', 'MarkerSize', 6);
text(rc(iI) + 1, I_max, sprintf('max: (%.2f, %.1f)', rc(iI), I_max), 'FontSize', 8, 'VerticalAlignment', 'bottom');
hold off;
xlabel('r_c [dim]', 'FontWeight', 'bold', 'FontSize', 10);
ylabel('I [m/s]', 'FontWeight', 'bold', 'FontSize', 10);
title('Specific Thrust vs. r_c', 'FontSize', 12);
grid on;

%% Figure 2: TSFC vs. rc
figure('Position', [150, 150, 600, 400]);
plot(rc, TSFC, 'LineWidth', 1.5);
hold on;
plot(rc(iT), TSFC_min, 'o', 'MarkerSize', 6);
text(rc(iT), TSFC_min, sprintf('min: (%.2f, %.3e)', rc(iT), TSFC_min), 'FontSize', 8, ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'center');
hold off;
ax = gca;
ax.YAxis.Exponent = 0;
xlabel('r_c [dim]', 'FontWeight', 'bold', 'FontSize', 10);
ylabel('TSFC [s/m]', 'FontWeight', 'bold', 'FontSize', 10);
title('TSFC vs. r_c', 'FontSize', 12);
grid on;

%% Figure 3: Efficiencies vs. rc
figure('Position', [200, 200, 600, 400]);
plot(rc, eta_th, 'LineWidth', 1.5, 'DisplayName', '\eta_{th}');
hold on;
plot(rc, eta_0, 'LineWidth', 1.5, 'DisplayName', '\eta_{o}');
plot(rc, eta_p, 'LineWidth', 1.5, 'Color', 'k', 'DisplayName', '\eta_p');
hold off;
legend('Location', 'southeast', 'FontSize', 9);
xlabel('r_c [dim]', 'FontWeight', 'bold', 'FontSize', 10);
ylabel('Efficiency', 'FontWeight', 'bold', 'FontSize', 10);
title('Efficiencies vs. r_c', 'FontSize', 12);
grid on;

%% Figure 4: Nozzle Area Ratio vs. rc
figure('Position', [250, 250, 600, 400]);
plot(rc, Ratio, 'LineWidth', 1.5);
ylim([0.9 * min(Ratio), 1.1 * max(Ratio)]);
xlabel('r_c [dim]', 'FontWeight', 'bold', 'FontSize', 10);
ylabel('Area Ratio', 'FontWeight', 'bold', 'FontSize', 10);
title('Nozzle Area Ratio vs. r_c', 'FontSize', 12);
grid on;

%% Save Data to CSV
data = table(rc, f_b, T03, T04, T05, P04, P05, T7, Me, Ue, W, I, TSFC, eta_th, eta_p, eta_0, Ratio, ...
    'VariableNames', {'rc', 'f_b', 'T03', 'T04', 'T05', 'P04', 'P05', 'T7', 'Me', 'Ue', 'W', ...
    'I', 'TSFC', 'eta_th', 'eta_p', 'eta_0', 'AreaRatio'});
writetable(data, 'HW7_Q1_Data.csv');
disp('Results saved to HW7_Q1_Data.csv');
