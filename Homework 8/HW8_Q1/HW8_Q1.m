%% EAS 4300 HW #8 Question 1
% Stephen Miller
%
% Created 4/13/22.
%
% You are to specify the bypass ratio of a turbofan engine for a comercial
% aircraft that is to cruise at a Mach number of 0.8. The exhaust from two
% streams are not internally mixed and both streams are ideally expanded.
% The design altitude is 10,000 m. The fuel has a heat of combustion of
% 43,000 kJ/kg, and a stoichiometric fuel-to-air ratio of 0.06. Assume the
% gamma is constant at 1.4 throughout the fan, 1.4 in the core up to the
% entrance of the burner, and 1.35 for the rest of the core (burner,
% turbine, and core nozzle). Assume the gas constant is fixed at
% R = 287 J/kgK. The maximum temperature at the inlet of the turbine is
% 1500 K. The compressor pressure ratio is 24, and the fan pressure ratio
% is 2.0. Plot I, TSFC, eta_o, eta_p, eta_th vs. bypass ratio. Is there an
% optimum bypass ratio that maximizes the overall efficiency? Use the
% following efficiencies: eta_d = 0.94, eta_c = 0.87, eta_f = 0.92,
% eta_b = 0.98, rb = 0.97, eta_t = 0.85, eta_cn = 0.97, eta_fn = 0.98.

clc; clear; close all;

%% Given Parameters
M_f = 0.8;                          % Flight Mach number
Ta = 223.252;                       % Ambient temperature at 10 km [K]
Pa = 26500;                         % Ambient pressure at 10 km [Pa]
R = 0.287;                          % Gas constant [kJ/(kg*K)]
T04_max = 1500;                     % Maximum turbine inlet temperature [K]
dhc = 43000;                        % Fuel heat of combustion [kJ/kg]
f_stoic = 0.06;                     % Stoichiometric fuel-to-air ratio

% Pressure ratios
rc = 24;       % Compressor pressure ratio
rb = 0.97;     % Burner pressure loss factor
rf = 2;        % Fan pressure ratio

% Efficiencies
eta_d = 0.94;  % Diffuser efficiency
eta_c = 0.87;  % Compressor efficiency
eta_b = 0.98;  % Burner efficiency
eta_t = 0.85;  % Turbine efficiency
eta_f = 0.92;  % Fan efficiency
eta_cn = 0.97; % Core nozzle efficiency
eta_fn = 0.98; % Fan nozzle efficiency

% Specific heat ratios
% k(1): inlet, diffuser, compressor, fan
% k(2): burner, turbine, core nozzle
k = [1.4 1.35];
Cp = (k ./ (k - 1)) * R;            % Specific heats [kJ/(kg*K)]

% Bypass ratio range. Above beta = 7.26 the turbine can no longer drive
% the fan and still leave enough pressure to expand the core flow to Pa.
beta = linspace(2, 7.26, 100)';
n = length(beta);

%% Inlet, Compressor, Fan, and Burner
% These do not depend on the bypass ratio.
%
% *Inlet and diffuser* (adiabatic, so $T_{02} = T_{0a}$):
%
% $$T_{0a} = T_a\left(1 + \frac{\gamma_1 - 1}{2}M^2\right), \qquad u = M\sqrt{\gamma_1 R T_a}$$
%
% $$T_{02s} = T_a + \eta_d\,(T_{0a} - T_a), \qquad P_{02} = P_a\left(\frac{T_{02s}}{T_a}\right)^{\gamma_1/(\gamma_1 - 1)}$$
%
% *Compressor* and *fan*:
%
% $$T_{03} = T_{02} + \frac{T_{02}\,r_c^{(\gamma_1-1)/\gamma_1} - T_{02}}{\eta_c}, \qquad w_c = c_{p1}(T_{03} - T_{02}), \qquad P_{03} = r_c P_{02}$$
%
% $$T_{07} = T_{02} + \frac{T_{02}\,r_f^{(\gamma_1-1)/\gamma_1} - T_{02}}{\eta_f}, \qquad P_{07} = r_f P_{02}$$
%
% *Burner* energy balance, with $c_{p1}$ for the incoming air and
% $c_{p2}$ for the products:
%
% $$c_{p1}T_{03} + f\,\eta_b\,\Delta h_c = (1 + f)\,c_{p2}T_{04} \quad\Rightarrow\quad f = \frac{c_{p2}T_{04} - c_{p1}T_{03}}{\eta_b\,\Delta h_c - c_{p2}T_{04}}, \qquad P_{04} = r_b P_{03}$$
%
% If $f$ would exceed $f_{st}$, the burner runs at $f_{st}$ and $T_{04}$ is
% solved from the same balance.

% Inlet (diffuser) conditions
T0a = Ta * (1 + ((k(1)-1)/2)*M_f^2);
T02 = T0a;
T02s = eta_d * (T02 - Ta) + Ta;
P02 = Pa * (T02s/Ta)^(k(1)/(k(1)-1));
U = M_f * sqrt(k(1)*R*1000*Ta);     % Freestream velocity [m/s]

% Compressor
T03s = T02 * rc^((k(1)-1)/k(1));
T03 = (T03s - T02)/eta_c + T02;
Wc = Cp(1) * (T03 - T02);           % Compressor work per kg of core air [kJ/kg]
P03 = rc * P02;

% Fan
T07s = T02 * rf^((k(1)-1)/k(1));
T07 = (T07s - T02)/eta_f + T02;
P07 = rf * P02;

% Burner
Fb = (Cp(2)*T04_max - Cp(1)*T03) / (eta_b*dhc - Cp(2)*T04_max);
T04 = T04_max;
if Fb > f_stoic
    Fb = f_stoic;
    T04 = (Cp(1)*T03 + f_stoic*eta_b*dhc) / ((1 + f_stoic)*Cp(2));
end
P04 = rb * P03;

%% Turbine, Nozzles, and Performance for Each Bypass Ratio
% *Turbine* drives the compressor and the fan ($\beta$ kg of bypass air
% per kg of core air):
%
% $$(1 + f)\,c_{p2}(T_{04} - T_{05}) = w_c + \beta\,c_{p1}(T_{07} - T_{02}), \qquad T_{05s} = T_{04} - \frac{T_{04} - T_{05}}{\eta_t}, \qquad P_{05} = P_{04}\left(\frac{T_{05s}}{T_{04}}\right)^{\gamma_2/(\gamma_2-1)}$$
%
% *Core nozzle* (hot) and *fan nozzle* (cold), both ideally expanded to
% $P_a$. For each, with inlet state $(T_0, P_0)$, efficiency $\eta_n$ and
% $\gamma$:
%
% $$T_{es} = T_0\left(\frac{P_a}{P_0}\right)^{(\gamma-1)/\gamma}, \qquad T_e = T_0 - \eta_n(T_0 - T_{es}), \qquad M_e = \sqrt{\frac{2}{\gamma - 1}\left(\frac{T_0}{T_e} - 1\right)}, \qquad u_e = M_e\sqrt{\gamma R T_e}$$
%
% using $(T_{05}, P_{05}, \eta_{cn}, \gamma_2)$ for the core and
% $(T_{07}, P_{07}, \eta_{fn}, \gamma_1)$ for the fan.
%
% *Performance* (per kg of core air):
%
% $$I = (1 + f)\,u_{e,hot} - u + \beta\,(u_{e,cold} - u), \qquad \mathrm{TSFC} = \frac{f}{I}$$
%
% $$\eta_{th} = \frac{(1 + f)\,u_{e,hot}^2 + \beta\,u_{e,cold}^2 - (1 + \beta)\,u^2}{2 f\,\Delta h_c}, \qquad \eta_p = \frac{2\,I\,u}{(1 + f)\,u_{e,hot}^2 + \beta\,u_{e,cold}^2 - (1 + \beta)\,u^2}, \qquad \eta_0 = \eta_{th}\,\eta_p$$

Wf = zeros(n,1);
T05 = zeros(n,1);
P05 = zeros(n,1);
Me = zeros(n,1);
M8 = zeros(n,1);
Ue_hot = zeros(n,1);
Ue_cold = zeros(n,1);
I = zeros(n,1);
TSFC = zeros(n,1);
eta_th = zeros(n,1);
eta_p = zeros(n,1);
eta_0 = zeros(n,1);

% Fan nozzle (same for every beta)
T8s = T07 / (P07/Pa)^((k(1)-1)/k(1));
T8 = T07 - eta_fn * (T07 - T8s);
M8(:) = sqrt((2*(T07/T8 - 1))/(k(1)-1));
Ue_cold(:) = M8(1) * sqrt(k(1)*R*1000*T8);

for i = 1:n
    % Turbine: work extraction drives the compressor and fan
    Wf(i) = beta(i) * Cp(1) * (T07 - T02);
    T05(i) = T04 - (Wc + Wf(i)) / ((1 + Fb) * Cp(2));
    T05s = T04 - (T04 - T05(i)) / eta_t;
    P05(i) = P04 * (T05s/T04)^(k(2)/(k(2)-1));
    if P05(i) <= Pa
        error('At beta = %.2f the turbine exit pressure is below Pa; lower the beta range.', beta(i));
    end

    % Core nozzle: expansion to ambient
    T6s = T05(i) / (P05(i)/Pa)^((k(2)-1)/k(2));
    T6 = T05(i) - eta_cn * (T05(i) - T6s);
    Me(i) = sqrt((2*(T05(i)/T6 - 1))/(k(2)-1));
    Ue_hot(i) = Me(i) * sqrt(k(2)*R*1000*T6);

    % Specific thrust and TSFC
    I(i) = (1 + Fb) * Ue_hot(i) - U + beta(i)*(Ue_cold(i) - U);
    TSFC(i) = Fb / I(i);

    % Efficiencies
    KE = (1 + Fb)*Ue_hot(i)^2 + beta(i)*Ue_cold(i)^2 - (beta(i)+1)*U^2;
    eta_th(i) = KE / (2*Fb*dhc*1000);
    eta_p(i) = 2*I(i)*U / KE;
    eta_0(i) = eta_th(i) * eta_p(i);
end

%% Optimum Bypass Ratio
[eta_0_max, i0] = max(eta_0);
[I_max, iI] = max(I);
[TSFC_min, iT] = min(TSFC);
fprintf('Fuel-to-air ratio: f = %.4f\n', Fb);
fprintf('Maximum overall efficiency: eta_0 = %.4f at beta = %.2f\n', eta_0_max, beta(i0));
fprintf('Maximum specific thrust: I = %.1f m/s at beta = %.2f\n', I_max, beta(iI));
fprintf('Minimum TSFC: %.4e s/m at beta = %.2f\n', TSFC_min, beta(iT));

%% Figure 1: Specific Thrust vs. Bypass Ratio
figure;
plot(beta, I, 'LineWidth', 1.5)
hold on
plot(beta(iI), I_max, 'o', 'MarkerSize', 6)
hold off
xlabel('\beta', 'FontWeight', 'bold')
ylabel('Specific Thrust (m/s)', 'FontWeight', 'bold')
title('Specific Thrust vs. Bypass Ratio')
grid on

%% Figure 2: TSFC vs. Bypass Ratio
figure;
plot(beta, TSFC, 'LineWidth', 1.5)
hold on
plot(beta(iT), TSFC_min, 'o', 'MarkerSize', 6)
hold off
ax = gca;
ax.YAxis.Exponent = 0;
xlabel('\beta', 'FontWeight', 'bold')
ylabel('TSFC (s/m)', 'FontWeight', 'bold')
title('TSFC vs. Bypass Ratio')
grid on

%% Figure 3: Efficiencies vs. Bypass Ratio
figure;
plot(beta, eta_th, 'LineWidth', 1.5)
hold on
plot(beta, eta_p, 'LineWidth', 1.5)
plot(beta, eta_0, 'LineWidth', 1.5)
plot(beta(i0), eta_0_max, 'ko', 'MarkerSize', 6)
text(beta(i0), eta_0_max, sprintf('  max \\eta_0 = %.3f at \\beta = %.2f', eta_0_max, beta(i0)), ...
    'VerticalAlignment', 'bottom', 'FontSize', 8)
hold off
xlabel('\beta', 'FontWeight', 'bold')
ylabel('Efficiency', 'FontWeight', 'bold')
title('Efficiencies vs. Bypass Ratio')
legend('\eta_{th}', '\eta_{p}', '\eta_{0}', 'Location', 'best')
grid on

%% Figure 4: Exit Mach Numbers vs. Bypass Ratio
figure;
plot(beta, Me, 'LineWidth', 1.5)
hold on
plot(beta, M8, 'LineWidth', 1.5)
hold off
xlabel('\beta', 'FontWeight', 'bold')
ylabel('Mach Number', 'FontWeight', 'bold')
title('Exit Mach Numbers vs. Bypass Ratio')
legend('Core', 'Bypass', 'Location', 'best')
grid on

%% Save Data to CSV
dataTable = table(beta, I, TSFC, eta_th, eta_p, eta_0, Me, M8, repmat(Fb, n, 1), ...
    'VariableNames', {'BypassRatio','SpecificThrust','TSFC','ThermalEff','PropEff', ...
    'OverallEff','CoreMach','BypassMach','FuelAirRatio'});
writetable(dataTable, 'HW8_Q1_data.csv');
disp('Data saved to HW8_Q1_data.csv');
