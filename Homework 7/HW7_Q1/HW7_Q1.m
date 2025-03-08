% HW7_Q1.m
% A non-afterburning turbojet analysis for operation at 15 km altitude and Mach 1.8.
% Plots specific thrust, TSFC, efficiencies, and nozzle area ratio vs. compressor pressure ratio (rc).
% Created by Stephen Miller on 3/28/22.
% Modified to produce multiple figures, each appended to a single PDF, with smaller figure sizes.

clc; clear; clear vars;

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
rb       = 0.97;                        % Burner pressure ratio (P03/P02)
k        = [1.4 1.3];                  % Specific heat ratios (before/after burner)
R        = 0.287;                       % Gas constant [kJ/kg·K]
Ta       = 216.65;                      % Ambient temperature at 15 km [K]
Pa       = 1.2112e4;                    % Ambient pressure at 15 km [Pa]
delta    = 58;                          % Number of data points for rc
rc       = linspace(2, 60, delta);      % Compressor pressure ratio range

%% Preallocate Variables
n        = length(rc);
num      = length(k);
Cp       = zeros(num, 1);
T03s     = zeros(n, 1);
T03      = zeros(n, 1);
W        = zeros(n, 1);
f_b      = zeros(n, 1);
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
T04      = zeros(n, 1);

%% Initial Conditions
% Stagnation temperature and pressure at inlet
T0a      = Ta * (1 + ((k(1) - 1)/2) * M^2);                    % T0a = T02
P0a      = Pa * (1 + ((k(1) - 1)/2) * M^2)^(k(1)/(k(1) - 1));  % P0a = P02
u        = M * sqrt(R * 1000 * k(1) * Ta);                     % Flight velocity [m/s]
T02s     = eta_d * (T0a - Ta) + Ta;                            % Ideal T02
P02      = Pa * (T02s / Ta)^(k(1) / (k(1) - 1));                % Inlet pressure

% Specific heat calculations
for i = 1:num
    Cp(i) = (k(i) / (k(i) - 1)) * R;
end

%% Main Loop: Calculate Parameters for Each rc
for j = 1:n
    % Compressor
    T03s(j) = T0a * rc(j)^((k(1) - 1) / k(1));          % Ideal compressor exit temperature
    T03(j)  = ((T03s(j) - T0a) / eta_c) + T0a;          % Actual compressor exit temperature
    W(j)    = Cp(1) * (T03(j) - T0a);                   % Compressor work [kJ/kg]

    % Burner: Fuel-air ratio to reach T04_max
    f_b(j)  = ((T04_max / T03(j)) - 1) / ...
              (eta_b * delta_h_c / (Cp(2) * T03(j)) - (T04_max / T03(j)));

    if f_b(j) > f_st
        % If fuel-air ratio exceeds stoichiometric limit
        T04(j) = (1 / (1 + f_st)) * ((eta_b * delta_h_c) / (Cp(2)) + T03);
        P04(j) = rb * rc(j) * P02;
        T05(j) = T04(j) - (W(j) / (Cp(2) * (1 + f_b(j)))); % T05 = T06 = T07
        T05s(j)= ((T05(j) - T04(j)) / eta_t) + T04(j);
        P05(j) = P04(j) * (T05s(j) / T04(j))^(k(2) / (k(2) - 1)); % P05 = P06
        T7s(j) = T05(j) / (P05(j) / Pa)^((k(2) - 1) / k(2));
        T7(j)  = eta_n * (T7s(j) - T05(j)) + T05(j);
        Me(j)  = sqrt(((T05(j) / T7(j)) - 1) / ((k(2) - 1) / 2));
        Ue(j)  = Me(j) * sqrt(k(2) * R * 1000 * T7(j)); % Exhaust velocity [m/s]
        I(j)   = (1 + f_b(j)) * Ue(j) - u;              % Specific thrust [m/s]
        TSFC(j)= f_b(j) / I(j);                         % Thrust-specific fuel consumption [s/m]
        eta_p(j)= (I(j) * u) / ((1 + f_b(j)) * (Ue(j)^2) / 2 - (u^2) / 2); 
        eta_th(j)= ((1 + f_b(j)) * Ue(j)^2 - u^2) / (2 * f_b(j) * delta_h_c * 1000); 
        eta_0(j)= eta_th(j) * eta_p(j);                 % Overall efficiency
        Ratio(j)= (1 / Me(j)) * ((2 / (k(2) + 1)) * (1 + ((k(2) - 1) / 2) * Me(j)^2))^((k(2) + 1) / (2 * (k(2) - 1)));
    else
        % Normal operation within stoichiometric limit
        P04(j) = rb * rc(j) * P02;
        T05(j) = T04_max - (W(j) / (Cp(2) * (1 + f_b(j)))); % T05 = T06 = T07
        T05s(j)= ((T05(j) - T04_max) / eta_t) + T04_max;
        P05(j) = P04(j) * (T05s(j) / T04_max)^(k(2) / (k(2) - 1)); % P05 = P06
        T7s(j) = T05(j) / (P05(j) / Pa)^((k(2) - 1) / k(2));
        T7(j)  = eta_n * (T7s(j) - T05(j)) + T05(j);
        Me(j)  = sqrt(((T05(j) / T7(j)) - 1) / ((k(2) - 1) / 2));
        Ue(j)  = Me(j) * sqrt(k(2) * R * 1000 * T7(j)); % Exhaust velocity [m/s]
        I(j)   = (1 + f_b(j)) * Ue(j) - u;              % Specific thrust [m/s]
        TSFC(j)= f_b(j) / I(j);                         % Thrust-specific fuel consumption [s/m]
        eta_p(j)= (I(j) * u) / ((1 + f_b(j)) * (Ue(j)^2) / 2 - (u^2) / 2);
        eta_th(j)= ((1 + f_b(j)) * Ue(j)^2 - u^2) / (2 * f_b(j) * delta_h_c * 1000);
        eta_0(j)= eta_th(j) * eta_p(j);
        Ratio(j)= (1 / Me(j)) * ((2 / (k(2) + 1)) * (1 + ((k(2) - 1) / 2) * Me(j)^2))^((k(2) + 1) / (2 * (k(2) - 1)));
    end
end

%% Plotting (Four Separate Figures, Appended to One PDF)
% (Requires MATLAB R2020a or newer for exportgraphics with 'Append' option)

% 1) Specific Thrust vs. r_c
figure('Position', [100, 100, 600, 400]);
plot(rc, I, 'LineWidth', 1.5);
hold on;
[Peak, PeakIdx] = findpeaks(I);
rc_Peak = rc(PeakIdx);
if ~isempty(rc_Peak)
    text(rc_Peak, 715, sprintf('(%g, %4.4f)', rc_Peak, Peak), 'FontSize', 8);
    plot(rc_Peak, Peak, '-o', 'MarkerSize', 6);
end
hold off;
xlabel('r_c [dim]', 'FontWeight', 'bold', 'FontSize', 10);
ylabel('I [m/s]', 'FontWeight', 'bold', 'FontSize', 10);
title('Specific Thrust vs. r_c', 'FontSize', 12);
grid on;
exportgraphics(gcf, 'HW7_Q1_Plots.pdf', 'Append', false);

% 2) TSFC vs. r_c
figure('Position', [150, 150, 600, 400]);
plot(rc, TSFC, 'LineWidth', 1.5);
hold on;
plot(60, min(TSFC), '-o', 'MarkerSize', 6);
text(54, 0.000022, sprintf('(60, %f)', min(TSFC)), 'FontSize', 8);
ax = gca;
ax.YAxis.Exponent = 0;
xlim([0 65]);
hold off;
xlabel('r_c [dim]', 'FontWeight', 'bold', 'FontSize', 10);
ylabel('TSFC [s/m]', 'FontWeight', 'bold', 'FontSize', 10);
title('TSFC vs. r_c', 'FontSize', 12);
grid on;
exportgraphics(gcf, 'HW7_Q1_Plots.pdf', 'Append', true);

% 3) Efficiencies vs. r_c
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
exportgraphics(gcf, 'HW7_Q1_Plots.pdf', 'Append', true);

% 4) Nozzle Area Ratio vs. r_c
figure('Position', [250, 250, 600, 400]);
plot(rc, Ratio, 'LineWidth', 1.5);

% Dynamically set y-limits so the curve is fully visible
ratioMin = min(Ratio);
ratioMax = max(Ratio);
ylim([0.9*ratioMin, 1.1*ratioMax]);

xlabel('r_c [dim]', 'FontWeight', 'bold', 'FontSize', 10);
ylabel('Area Ratio', 'FontWeight', 'bold', 'FontSize', 10);
title('Nozzle Area Ratio vs. r_c', 'FontSize', 12);
grid on;

exportgraphics(gcf, 'HW7_Q1_Plots.pdf', 'Append', true);


%% Save Data to CSV
data = table(f_b, I, TSFC, Me, P04, P05, T03, T03s, T05, T05s, T7, T7, T7s, TSFC, Ue, W, ...
    'VariableNames', {'f_b', 'I', 'TSFC', 'Me', 'P04', 'P05', 'T03', 'T03s', 'T05', 'T05s', 'T7', 'T7_1', 'T7s', 'TSFC_1', 'Ue', 'W'});
writetable(data, 'HW7_Q1_Data.csv');

%% Display Data (Optional)
disp(data);
