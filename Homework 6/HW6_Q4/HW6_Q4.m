% HW6_Q4 - Ramjet Performance Sensitivity Study
clear all; close all;

% Constants
gamma = 1.4; % Specific heat ratio
R = 287; % Gas constant (J/kg·K)
cp = 1000; % Specific heat at constant pressure (J/kg·K)
p0 = 8500; % Ambient pressure (Pa)
T0 = 220; % Ambient temperature (K)
hf = 43000000; % Fuel heat of combustion (J/kg)
f = 0.06; % Fuel-to-air ratio
Tt4 = 2540; % Turbine inlet temperature (K)

% Flight Mach number range
M_flight = 1:0.1:6;
nM = length(M_flight);

% Initialize arrays for baseline and perturbed cases
I_base = zeros(1, nM); % Thrust
TSFC_base = zeros(1, nM); % Thrust-specific fuel consumption
I_deta = zeros(1, nM); % Thrust for perturbed eta_b
TSFC_deta = zeros(1, nM); % TSFC for perturbed eta_b
I_dpi = zeros(1, nM); % Thrust for perturbed pi_n
TSFC_dpi = zeros(1, nM); % TSFC for perturbed pi_n

% Baseline eta_b and pi_n
eta_b = 1;
pi_n = 1;
delta = 0.01; % Perturbation

% Loop over Mach numbers
for i = 1:nM
    M0 = M_flight(i);
    
    % Ambient conditions
    a0 = sqrt(gamma * R * T0);
    u0 = M0 * a0;
    
    % Stagnation conditions
    Tt0 = T0 * (1 + (gamma - 1)/2 * M0^2);
    pt0 = p0 * (1 + (gamma - 1)/2 * M0^2)^(gamma/(gamma - 1));
    
    % Baseline case
    [I_base(i), TSFC_base(i)] = compute_performance(M0, eta_b, pi_n, Tt0, pt0, p0, T0, f, hf, gamma, R, cp);
    
    % Perturbed eta_b
    [I_deta(i), TSFC_deta(i)] = compute_performance(M0, eta_b + delta, pi_n, Tt0, pt0, p0, T0, f, hf, gamma, R, cp);
    
    % Perturbed pi_n
    [I_dpi(i), TSFC_dpi(i)] = compute_performance(M0, eta_b, pi_n + delta, Tt0, pt0, p0, T0, f, hf, gamma, R, cp);
end

% Compute derivatives using finite differences
dI_deta = (I_deta - I_base) / delta;
dI_dpi = (I_dpi - I_base) / delta;
dTSFC_deta = (TSFC_deta - TSFC_base) / delta;
dTSFC_dpi = (TSFC_dpi - TSFC_base) / delta;

% Normalized sensitivity coefficients
norm_dI_deta = (eta_b ./ I_base) .* dI_deta;
norm_dI_dpi = (pi_n ./ I_base) .* dI_dpi;
norm_dTSFC_deta = (eta_b ./ TSFC_base) .* dTSFC_deta;
norm_dTSFC_dpi = (pi_n ./ TSFC_base) .* dTSFC_dpi;

% Plotting
figure;
subplot(2, 2, 1);
plot(M_flight, norm_dI_deta, 'b-', 'LineWidth', 2);
xlabel('Flight Mach Number');
ylabel('(1/I) dI/d\eta_b');
title('Normalized Sensitivity of Thrust to \eta_b');
grid on;

subplot(2, 2, 2);
plot(M_flight, norm_dI_dpi, 'r-', 'LineWidth', 2);
xlabel('Flight Mach Number');
ylabel('(1/I) dI/d\pi_n');
title('Normalized Sensitivity of Thrust to \pi_n');
grid on;

subplot(2, 2, 3);
plot(M_flight, norm_dTSFC_deta, 'b-', 'LineWidth', 2);
xlabel('Flight Mach Number');
ylabel('(1/TSFC) dTSFC/d\eta_b');
title('Normalized Sensitivity of TSFC to \eta_b');
grid on;

subplot(2, 2, 4);
plot(M_flight, norm_dTSFC_dpi, 'r-', 'LineWidth', 2);
xlabel('Flight Mach Number');
ylabel('(1/TSFC) dTSFC/d\pi_n');
title('Normalized Sensitivity of TSFC to \pi_n');
grid on;

% Save Figure 1 as PDF
set(gcf, 'PaperPositionMode', 'auto'); % Use screen size for PDF
print('-dpdf', 'HW6_Q4_Figure1.pdf', '-bestfit'); % Save the entire figure as PDF

% Save results to CSV file
data = [M_flight' I_base' TSFC_base' norm_dI_deta' norm_dI_dpi' norm_dTSFC_deta' norm_dTSFC_dpi'];
headers = {'Mach_Number', 'Thrust_I', 'TSFC', 'Norm_dI_deta', 'Norm_dI_dpi', 'Norm_dTSFC_deta', 'Norm_dTSFC_dpi'};
T = array2table(data, 'VariableNames', headers);
writetable(T, 'HW6_Q4_results.csv');

disp('Results saved to HW6_Q4_results.csv and HW6_Q4_Figure1.pdf');

% Helper function to compute performance
function [I, TSFC] = compute_performance(M0, eta_b, pi_n, Tt0, pt0, p0, T0, f, hf, gamma, R, cp)
    % Flight velocity
    a0 = sqrt(gamma * R * T0);
    u0 = M0 * a0;
    
    % Post-combustion temperature (adjusted for eta_b)
    Tt4 = Tt0 + eta_b * f * hf / cp;
    
    % Nozzle total pressure (adjusted for pi_n)
    pt9 = pt0 * pi_n;
    
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