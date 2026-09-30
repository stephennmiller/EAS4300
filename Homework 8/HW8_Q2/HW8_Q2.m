%% EAS 4300 HW #8 Question 2
% Stephen Miller
%
% Use the model developed in problem 1 to further optimize rc and rf at
% the optimum bypass ratio determined in problem 1. Do this using a single
% table where these values are varied across the following ranges:
% 1.5 < rf < 2.2, 20 < rc < 28. Make a table with 900 rows. Vary rf over
% the range using the "repeat pattern" across 30 rows. Vary rc over the
% range using the "apply pattern" across 30 rows. From this table
% solution, you can make an x-y-z plot, with eta_o(rc, rf). For the
% conditions and efficiencies listed in problem 1, determine if a
% different combination of rc and rf (other than 24 and 2) result in a
% further increase in the overall efficiency.

clear; clc; close all;

%% Given Parameters
M_f     = 0.8;      % Flight Mach number
Ta      = 223.252;  % Ambient temperature at 10 km [K]
Pa      = 26500;    % Ambient pressure at 10 km [Pa]
R       = 0.287;    % Gas constant [kJ/(kg*K)]
T04_max = 1500;     % Turbine inlet temperature [K]
dhc     = 43000;    % Fuel heat of combustion [kJ/kg]
f_stoic = 0.06;     % Stoichiometric fuel-to-air ratio
beta    = 5.88;     % Optimum bypass ratio from Question 1
rb      = 0.97;     % Burner pressure loss factor

% Efficiencies
eta_d  = 0.94;  % Diffuser efficiency
eta_c  = 0.87;  % Compressor efficiency
eta_b  = 0.98;  % Burner efficiency
eta_t  = 0.85;  % Turbine efficiency
eta_f  = 0.92;  % Fan efficiency
eta_cn = 0.97;  % Core nozzle efficiency
eta_fn = 0.98;  % Fan nozzle efficiency

% Specific heat ratios
% k(1) for inlet, diffuser, compressor, fan; k(2) for burner, turbine, nozzle
k = [1.4, 1.35];
Cp = (k ./ (k - 1)) * R;   % Specific heats [kJ/(kg*K)]

%% 900-Row Table of rf and rc
% rf takes 30 values that repeat every 30 rows ("repeat pattern"); rc takes
% 30 values, each held for 30 rows ("apply pattern").
N     = 30;
rfVec = linspace(1.5, 2.2, N)';
rcVec = linspace(20, 28, N)';
rf    = repmat(rfVec, N, 1);       % 1.5 ... 2.2, 1.5 ... 2.2, ...
rc    = repelem(rcVec, N);         % 20 (x30), 20.28 (x30), ...
nRows = numel(rf);

%% Inlet Conditions
% Same inlet and diffuser as Question 1 ($T_{02} = T_{0a}$, adiabatic):
%
% $$T_{0a} = T_a\left(1 + \frac{\gamma_1 - 1}{2}M^2\right), \qquad T_{02s} = T_a + \eta_d\,(T_{0a} - T_a), \qquad P_{02} = P_a\left(\frac{T_{02s}}{T_a}\right)^{\gamma_1/(\gamma_1 - 1)}$$

T0a  = Ta * (1 + (k(1)-1)/2 * M_f^2);
T02  = T0a;
T02s = eta_d * (T0a - Ta) + Ta;
P02  = Pa * (T02s/Ta)^(k(1)/(k(1)-1));
U    = M_f * sqrt(k(1) * R * 1000 * Ta);   % Freestream velocity [m/s]

%% Turbofan Cycle for Each Row
% The cycle is the Question 1 model, evaluated at $\beta = 5.88$ for each
% $(r_c, r_f)$ pair:
%
% $$T_{03} = T_{02} + \frac{T_{02}\,r_c^{(\gamma_1-1)/\gamma_1} - T_{02}}{\eta_c}, \qquad T_{07} = T_{02} + \frac{T_{02}\,r_f^{(\gamma_1-1)/\gamma_1} - T_{02}}{\eta_f}$$
%
% $$f = \frac{c_{p2}T_{04} - c_{p1}T_{03}}{\eta_b\,\Delta h_c - c_{p2}T_{04}}, \qquad (1 + f)\,c_{p2}(T_{04} - T_{05}) = c_{p1}(T_{03} - T_{02}) + \beta\,c_{p1}(T_{07} - T_{02})$$
%
% $$P_{05} = r_b\,r_c\,P_{02}\left(\frac{T_{05s}}{T_{04}}\right)^{\gamma_2/(\gamma_2-1)}, \qquad T_{05s} = T_{04} - \frac{T_{04} - T_{05}}{\eta_t}$$
%
% Both nozzles expand ideally to $P_a$ (see Question 1), and
%
% $$I = (1 + f)\,u_{e,hot} - u + \beta\,(u_{e,cold} - u), \qquad \eta_0 = \eta_{th}\,\eta_p = \frac{I\,u}{f\,\Delta h_c}$$

eta_0 = zeros(nRows, 1);
I     = zeros(nRows, 1);
TSFC  = zeros(nRows, 1);

for j = 1:nRows
    % Compressor
    T03 = (T02 * rc(j)^((k(1)-1)/k(1)) - T02)/eta_c + T02;
    Wc  = Cp(1) * (T03 - T02);
    P04 = rb * rc(j) * P02;

    % Burner: Cp(1)*T03 + f*eta_b*dhc = (1 + f)*Cp(2)*T04
    Fb  = (Cp(2)*T04_max - Cp(1)*T03) / (eta_b*dhc - Cp(2)*T04_max);
    T04 = T04_max;
    if Fb > f_stoic
        Fb  = f_stoic;
        T04 = (Cp(1)*T03 + f_stoic*eta_b*dhc) / ((1 + f_stoic)*Cp(2));
    end

    % Fan (bypass)
    T07 = (T02 * rf(j)^((k(1)-1)/k(1)) - T02)/eta_f + T02;
    Wf  = beta * Cp(1) * (T07 - T02);
    P07 = rf(j) * P02;

    % Turbine drives compressor and fan
    T05  = T04 - (Wc + Wf) / ((1 + Fb)*Cp(2));
    T05s = T04 - (T04 - T05)/eta_t;
    P05  = P04 * (T05s/T04)^(k(2)/(k(2)-1));

    % Core and bypass nozzles, expanded to ambient
    T6  = T05 - eta_cn * (T05 - T05 / (P05/Pa)^((k(2)-1)/k(2)));
    Me  = sqrt((2*(T05/T6 - 1))/(k(2)-1));
    T8  = T07 - eta_fn * (T07 - T07 / (P07/Pa)^((k(1)-1)/k(1)));
    M8  = sqrt((2*(T07/T8 - 1))/(k(1)-1));
    Ue_hot  = Me * sqrt(k(2)*R*1000*T6);
    Ue_cold = M8 * sqrt(k(1)*R*1000*T8);

    % Performance
    I(j)     = (1 + Fb)*Ue_hot - U + beta*(Ue_cold - U);
    TSFC(j)  = Fb / I(j);
    KE       = (1 + Fb)*Ue_hot^2 + beta*Ue_cold^2 - (beta + 1)*U^2;
    eta_th   = KE / (2*Fb*dhc*1000);
    eta_p    = 2*I(j)*U / KE;
    eta_0(j) = eta_th * eta_p;
end

%% Best Combination of rc and rf
[eta_0_max, iBest] = max(eta_0);
[~, iBase] = min(abs(rc - 24) + abs(rf - 2));   % grid point closest to rc = 24, rf = 2
fprintf('Maximum overall efficiency: eta_0 = %.4f at rc = %.2f, rf = %.3f\n', ...
    eta_0_max, rc(iBest), rf(iBest));
fprintf('Near the Question 1 design (rc = %.2f, rf = %.3f): eta_0 = %.4f\n', ...
    rc(iBase), rf(iBase), eta_0(iBase));

% Reshape the table into a 30 x 30 grid: rows follow rf, columns follow rc
ETA0 = reshape(eta_0, N, N);
[RC, RF] = meshgrid(rcVec, rfVec);

%% Figure 1: Overall Efficiency Surface
figure;
surf(RC, RF, ETA0);
xlabel('Compressor Pressure Ratio (r_c)');
ylabel('Fan Pressure Ratio (r_f)');
zlabel('\eta_0');
title('Overall Efficiency \eta_0(r_c, r_f) at \beta = 5.88');
shading interp;
colorbar;
grid on;

%% Figure 2: Overall Efficiency Contours
figure;
[C, h] = contour(RC, RF, ETA0, 12);
clabel(C, h, 'FontSize', 8);
hold on
plot(rc(iBest), rf(iBest), 'r*', 'MarkerSize', 10);
plot(24, 2, 'ko', 'MarkerSize', 8);
hold off
legend('\eta_0 contours', 'Maximum', 'r_c = 24, r_f = 2', 'Location', 'southwest');
xlabel('Compressor Pressure Ratio (r_c)');
ylabel('Fan Pressure Ratio (r_f)');
title('Overall Efficiency Contours at \beta = 5.88');
grid on;

%% Save the 900-Row Table to CSV
dataTable = table(rf, rc, eta_0, I, TSFC, ...
    'VariableNames', {'FanPressureRatio', 'CompPressureRatio', 'OverallEfficiency', ...
    'SpecificThrust', 'TSFC'});
writetable(dataTable, 'HW8_Q2_data.csv');
disp('Data saved to HW8_Q2_data.csv');
