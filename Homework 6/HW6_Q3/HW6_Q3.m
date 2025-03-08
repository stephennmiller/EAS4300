%% EAS 4300 HW #6 Question 3
% Stephen Miller
% 
% 3-22-21 (Updated: 03-07-2025)
% 
% All rights reserved.

% Clear workspace, close all figures, and clear command window
clear all; close all; clc;

% Constants
f_st = 0.06;                              % Stoichiometric fuel-air ratio
T0_4_max = 2600;                          % Maximum combustor temperature (K)
delta_H_c = 43000;                        % Heat of combustion (kJ/kg)
k = [1.4 1.33];                           % Specific heat ratios (inlet: 1.4, combustor: 1.33)
R = 0.287;                                % Gas constant (kJ/kg·K)
delta = 1000;                             % Number of points for Mach range

% Ambient conditions at 10,000 meters
T_a = 223.252;                            % Ambient temperature (K)
P_a = 2.65e4;                             % Ambient pressure (Pa)

% Symbolic values
syms T0_4                           
syms M_e

% Specific heat at constant pressure for the combustor
Cp = (k(2)/(k(2)-1)) * R;                 % kJ/kg·K

% Preallocation of variables
M_flight = linspace(1, 6, delta);
n = length(M_flight);
P0_a = zeros(n,1);
Me = zeros(n,1);
Me_t = zeros(n,1);
A_exit_A_throat = zeros(n,1);             % Renamed for clarity
T0_a = zeros(n,1);
f = zeros(n,1);
T04 = zeros(n,1);
T6 = zeros(n,1);
u_e = zeros(n,1);
u = zeros(n,1);
I = zeros(n,1);                           % Specific thrust (m/s)
TSFC = zeros(n,1);                        % Thrust specific fuel consumption (kg/N·s)
eta_th = zeros(n,1);
eta_p = zeros(n,1);
eta_0 = zeros(n,1);

% Loop through Flight Mach numbers from 1 -> 6
for i = 1:n
    % Stagnation pressure at inlet
    P0_a(i) = P_a * (1 + ((k(1) - 1)/2) * M_flight(i)^2)^(k(1) / (k(1) - 1));
    
    % Flight velocity
    u(i) = M_flight(i) * sqrt(k(1) * R * 1000 * T_a); % m/s
    
    % Solve for exit Mach number
    Me(i) = vpasolve(P_a * (1 + ((k(2)-1)/2) * (M_e)^2)^(k(2)/(k(2)-1)) == P0_a(i), M_e, [0 Inf]);
    Me_t(i) = double(Me(i)); % Convert to numeric
    
    % Area ratio (A_exit / A_throat)
    A_exit_A_throat(i) = (1/Me_t(i)) * ((2/(k(2)+1)) * (1 + ((k(2)-1)/2) * Me_t(i)^2))^((k(2)+1)/(2*(k(2)-1)));
    
    % Stagnation temperature at inlet
    T0_a(i) = T_a * (1 + ((k(1)-1)/2) * M_flight(i)^2);
    
    % Initial fuel-to-air ratio calculation
    f(i) = ((T0_4_max / T0_a(i)) - 1) / (delta_H_c / (Cp * T0_a(i)) - (T0_4_max / T0_a(i)));
  
    % Special condition if f exceeds stoichiometric value
    if (f(i) > f_st)
        num = (T0_4 / T0_a(i)) - 1;
        den = (delta_H_c / (Cp * T0_a(i)) - (T0_4 / T0_a(i)));
        T04_sol = vpasolve(num / den == f_st, T0_4); % Solve symbolically
        T04(i) = double(T04_sol); % Convert to numeric
    else
        T04(i) = T0_4_max; % Assign numeric value
    end
    
    % Exit temperature
    T6(i) = T04(i) / (1 + ((k(2)-1)/2) * Me_t(i)^2);
    
    % Exit velocity
    u_e(i) = Me_t(i) * sqrt(k(2) * R * 1000 * T6(i)); % m/s
    
    % Specific thrust
    I(i) = (1 + f(i)) * u_e(i) - u(i); % m/s
    
    % Thrust specific fuel consumption
    TSFC(i) = f(i) / I(i); % kg/(N·s)
    
    % Efficiencies
    eta_p(i) = (I(i) * u(i)) / (((1 + f(i)) * (u_e(i)^2) / 2) - (u(i)^2) / 2);
    eta_th(i) = ((1 + f(i)) * u_e(i)^2 - u(i)^2) / (f(i) * delta_H_c * 1000); % Corrected thermal efficiency
    eta_0(i) = eta_th(i) * eta_p(i);
end

% Ensure all vectors are column vectors and numeric
M_flight = double(M_flight(:)); % Ensure column vector
Me_t = double(Me_t(:));
u = double(u(:));
u_e = double(u_e(:));
f = double(f(:));
I = double(I(:));
TSFC = double(TSFC(:));
T04 = double(T04(:));
A_exit_A_throat = double(A_exit_A_throat(:));
eta_th = double(eta_th(:));
eta_p = double(eta_p(:));
eta_0 = double(eta_0(:));

% Debugging: Check for NaN or Inf values
disp('Checking for NaN or Inf values:');
disp(['M_flight: ' num2str(any(isnan(M_flight) | isinf(M_flight)))]);
disp(['Me_t: ' num2str(any(isnan(Me_t) | isinf(Me_t)))]);
disp(['u: ' num2str(any(isnan(u) | isinf(u)))]);
disp(['u_e: ' num2str(any(isnan(u_e) | isinf(u_e)))]);
disp(['f: ' num2str(any(isnan(f) | isinf(f)))]);
disp(['I: ' num2str(any(isnan(I) | isinf(I)))]);
disp(['TSFC: ' num2str(any(isnan(TSFC) | isinf(TSFC)))]);
disp(['T04: ' num2str(any(isnan(T04) | isinf(T04)))]);
disp(['A_exit_A_throat: ' num2str(any(isnan(A_exit_A_throat) | isinf(A_exit_A_throat)))]);
disp(['eta_th: ' num2str(any(isnan(eta_th) | isinf(eta_th)))]);
disp(['eta_p: ' num2str(any(isnan(eta_p) | isinf(eta_p)))]);
disp(['eta_0: ' num2str(any(isnan(eta_0) | isinf(eta_0)))]);

% Plotting
% Figure 1: Specific Thrust vs. M_flight
figure(1)
plot(M_flight, I, 'LineWidth', 1.5)
xticks(1:6)
ylim([400 1200])
xlabel('Flight Mach Number (M_{flight})')
ylabel('Specific Thrust (m/s)')
title('Specific Thrust vs. Flight Mach Number')

% Figure 2: TSFC vs. M_flight
figure(2)
plot(M_flight, TSFC, 'LineWidth', 1.5)
xticks(1:6)
ax = gca;
ax.YAxis.Exponent = 0;
ylim([0.00004 0.0001])
xlabel('Flight Mach Number (M_{flight})')
ylabel('TSFC (kg/N·s)')
title('Thrust Specific Fuel Consumption vs. Flight Mach Number')

% Figure 3: T04 vs. M_flight
figure(3)
plot(M_flight, T04, 'LineWidth', 1.5)
xticks(1:6)
ylim([2350 2700])
xlabel('Flight Mach Number (M_{flight})')
ylabel('Combustor Exit Temperature (T_{0_4}, K)')
title('Combustor Exit Temperature vs. Flight Mach Number')

% Figure 4: A_exit/A_throat vs. M_flight
figure(4)
plot(M_flight, A_exit_A_throat, 'LineWidth', 1.5)
xticks(1:6)
xlabel('Flight Mach Number (M_{flight})')
ylabel('Area Ratio (A_{exit}/A_{throat})')
title('Area Ratio vs. Flight Mach Number')

% Figure 5: Efficiencies vs. M_flight
figure(5)
plot(M_flight, eta_0, 'b-', 'LineWidth', 1.5, 'DisplayName', '\eta_{0}')
hold on
plot(M_flight, eta_th, 'r-', 'LineWidth', 1.5, 'DisplayName', '\eta_{th}')
plot(M_flight, eta_p, 'g-', 'LineWidth', 1.5, 'DisplayName', '\eta_{p}')
hold off
xlabel('Flight Mach Number (M_{flight})')
ylabel('Efficiency')
title('Efficiencies vs. Flight Mach Number')
legend('Location', 'best')

% Figure 6: Fuel-to-air ratio vs. M_flight
figure(6)
plot(M_flight, f, 'LineWidth', 1.5)
xticks(1:6)
xlabel('Flight Mach Number (M_{flight})')
ylabel('Fuel-to-Air Ratio (f)')
yline(f_st, '--', 'LineWidth', 1.5)
text(5, 0.062, 'f_{st} = 0.06', 'FontWeight', 'bold', 'FontSize', 12)
title('Fuel-to-Air Ratio vs. Flight Mach Number')

% Data table
data = table(M_flight, Me_t, u, u_e, f, I, TSFC, T04, A_exit_A_throat, eta_th, eta_p, eta_0, ...
    'VariableNames', {'M_flight', 'Me_t', 'u', 'u_e', 'f', 'Specific_Thrust', 'TSFC', 'T04', 'A_exit_A_throat', 'eta_th', 'eta_p', 'eta_0'});

% Save data to a file
writetable(data, 'hw6_q3_results.csv');

disp('Calculation and plotting complete. Results saved to hw6_q3_results.csv.');