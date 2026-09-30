%% EAS 4300 HW #6 Question 4
% Stephen Miller
%
% Ramjet performance sensitivity study: normalized sensitivity of specific
% thrust (I) and TSFC to combustion efficiency (eta_b) and nozzle total
% pressure ratio (r_n), for flight Mach numbers 1 to 6.

% Clear workspace and close all figures
clear all; close all;

%% Constants
gamma = 1.4; % Specific heat ratio
R = 287; % Gas constant (J/kg·K)
cp = 1000; % Specific heat at constant pressure (J/kg·K)
p0 = 8500; % Ambient pressure (Pa)
T0 = 220; % Ambient temperature (K)
hf = 43000000; % Fuel heat of combustion (J/kg)
f_st = 0.06; % Stoichiometric fuel-to-air ratio
Tt4 = 2540; % Turbine inlet temperature (K)

% Flight Mach number range
M_flight = 1:0.1:6;
nM = length(M_flight);

% Initialize arrays for baseline and perturbed cases
I_base = zeros(1, nM); % Thrust
TSFC_base = zeros(1, nM); % Thrust-specific fuel consumption
I_deta = zeros(1, nM); % Thrust for perturbed eta_b
TSFC_deta = zeros(1, nM); % TSFC for perturbed eta_b
I_dr = zeros(1, nM); % Thrust for perturbed r_n
TSFC_dr = zeros(1, nM); % TSFC for perturbed r_n
f = zeros(1, nM); % Fuel-to-air ratio at each Mach number

% Baseline eta_b and r_n
eta_b = 1;
r_n = 1;
delta = 0.01; % Perturbation

%% Baseline and Perturbed Performance, M_flight = 1 to 6
for i = 1:nM
    M0 = M_flight(i);
    
    % Ambient conditions
    a0 = sqrt(gamma * R * T0);
    u0 = M0 * a0;
    
    % Stagnation conditions
    Tt0 = T0 * (1 + (gamma - 1)/2 * M0^2);
    pt0 = p0 * (1 + (gamma - 1)/2 * M0^2)^(gamma/(gamma - 1));

    % Fuel-to-air ratio that reaches Tt4 at the baseline eta_b (f << 1);
    % held fixed for the perturbed cases, so a change in eta_b changes Tt4
    f(i) = cp * (Tt4 - Tt0) / (eta_b * hf);
    if f(i) > f_st
        warning('f = %g exceeds f_st at M = %g', f(i), M0);
    end

    % Baseline case
    [I_base(i), TSFC_base(i)] = compute_performance(M0, eta_b, r_n, Tt0, pt0, p0, T0, f(i), hf, gamma, R, cp);

    % Perturbed eta_b
    [I_deta(i), TSFC_deta(i)] = compute_performance(M0, eta_b + delta, r_n, Tt0, pt0, p0, T0, f(i), hf, gamma, R, cp);

    % Perturbed r_n
    [I_dr(i), TSFC_dr(i)] = compute_performance(M0, eta_b, r_n + delta, Tt0, pt0, p0, T0, f(i), hf, gamma, R, cp);
end

%% Sensitivity Derivatives and Normalized Coefficients
% Derivatives approximated with finite differences (delta = 0.01)
dI_deta = (I_deta - I_base) / delta;
dI_dr = (I_dr - I_base) / delta;
dTSFC_deta = (TSFC_deta - TSFC_base) / delta;
dTSFC_dr = (TSFC_dr - TSFC_base) / delta;

% Normalize by the baseline values
norm_dI_deta = (eta_b ./ I_base) .* dI_deta;
norm_dI_dr = (r_n ./ I_base) .* dI_dr;
norm_dTSFC_deta = (eta_b ./ TSFC_base) .* dTSFC_deta;
norm_dTSFC_dr = (r_n ./ TSFC_base) .* dTSFC_dr;

%% Figure 1: Sensitivity Derivatives vs. Flight Mach Number
figure;
subplot(2, 2, 1);
plot(M_flight, dI_deta, 'b-', 'LineWidth', 2);
xlabel('Flight Mach Number');
ylabel('dI/d\eta_b (m/s)');
title('Thrust Sensitivity to \eta_b');
xlim([1 6]);
grid on;

subplot(2, 2, 2);
plot(M_flight, dI_dr, 'r-', 'LineWidth', 2);
xlabel('Flight Mach Number');
ylabel('dI/dr_n (m/s)');
title('Thrust Sensitivity to r_n');
xlim([1 6]);
grid on;

subplot(2, 2, 3);
plot(M_flight, dTSFC_deta, 'b-', 'LineWidth', 2);
xlabel('Flight Mach Number');
ylabel('dTSFC/d\eta_b (kg/N·hr)');
title('TSFC Sensitivity to \eta_b');
xlim([1 6]);
grid on;

subplot(2, 2, 4);
plot(M_flight, dTSFC_dr, 'r-', 'LineWidth', 2);
xlabel('Flight Mach Number');
ylabel('dTSFC/dr_n (kg/N·hr)');
title('TSFC Sensitivity to r_n');
xlim([1 6]);
grid on;

%% Figure 2: Normalized Sensitivities vs. Flight Mach Number
figure;
subplot(2, 2, 1);
plot(M_flight, norm_dI_deta, 'b-', 'LineWidth', 2);
xlabel('Flight Mach Number');
ylabel('(1/I) dI/d\eta_b');
title('Normalized Sensitivity of Thrust to \eta_b');
xlim([1 6]);
grid on;

subplot(2, 2, 2);
plot(M_flight, norm_dI_dr, 'r-', 'LineWidth', 2);
xlabel('Flight Mach Number');
ylabel('(1/I) dI/dr_n');
title('Normalized Sensitivity of Thrust to r_n');
xlim([1 6]);
ylim([0 1.2]);
grid on;

subplot(2, 2, 3);
plot(M_flight, norm_dTSFC_deta, 'b-', 'LineWidth', 2);
xlabel('Flight Mach Number');
ylabel('(1/TSFC) dTSFC/d\eta_b');
title('Normalized Sensitivity of TSFC to \eta_b');
xlim([1 6]);
grid on;

subplot(2, 2, 4);
plot(M_flight, norm_dTSFC_dr, 'r-', 'LineWidth', 2);
xlabel('Flight Mach Number');
ylabel('(1/TSFC) dTSFC/dr_n');
title('Normalized Sensitivity of TSFC to r_n');
xlim([1 6]);
ylim([-1.2 0]);
grid on;

%% Save Results
data = [M_flight' f' I_base' TSFC_base' dI_deta' dI_dr' dTSFC_deta' dTSFC_dr' norm_dI_deta' norm_dI_dr' norm_dTSFC_deta' norm_dTSFC_dr'];
headers = {'Mach_Number', 'f', 'Thrust_I', 'TSFC', 'dI_deta', 'dI_dr', 'dTSFC_deta', 'dTSFC_dr', 'Norm_dI_deta', 'Norm_dI_dr', 'Norm_dTSFC_deta', 'Norm_dTSFC_dr'};
T = array2table(data, 'VariableNames', headers);
writetable(T, 'HW6_Q4_results.csv');

disp('Results saved to HW6_Q4_results.csv');

%% Helper Function: Ideal Ramjet Performance
function [I, TSFC] = compute_performance(M0, eta_b, r_n, Tt0, pt0, p0, T0, f, hf, gamma, R, cp)
    % Flight velocity
    a0 = sqrt(gamma * R * T0);
    u0 = M0 * a0;
    
    % Post-combustion temperature (adjusted for eta_b)
    Tt4 = Tt0 + eta_b * f * hf / cp;
    
    % Nozzle total pressure (adjusted for r_n)
    pt9 = pt0 * r_n;
    
    % Nozzle exit conditions
    T9 = Tt4 * (p0 / pt9)^((gamma - 1)/gamma);
    a9 = sqrt(gamma * R * T9);
    M9 = sqrt((2/(gamma - 1)) * ((pt9/p0)^((gamma - 1)/gamma) - 1));
    u9 = M9 * a9;
    
    % Specific thrust
    I = u9 - u0;
    
    % TSFC
    TSFC = f / I * 3600; % Convert to kg/N·hr
end