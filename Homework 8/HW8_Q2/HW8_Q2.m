%% EAS 4300 HW8 Q2
% Stephen Miller
% 2D Parametric Study of Turbofan Engine Overall Efficiency
%
% This script performs a 2D parametric sweep over the fan pressure ratio (rf)
% and compressor pressure ratio (rc) to compute the overall efficiency (η₀)
% of a turbofan engine. It exports the computed data to a CSV file and
% saves the 3D surface plot as a PDF file.

clear; clc; close all;

%% Given Parameters
M_f     = 0.8;      % Flight Mach number
Ta      = 223.252;  % Ambient temperature [K]
Pa      = 26500;    % Ambient pressure [Pa]
R       = 0.287;    % Gas constant [kJ/(kg*K)]
T04_max = 1500;     % Turbine inlet temperature [K]
dhc     = 43000;    % Fuel heat of combustion [kJ/kg]
beta    = 5.88;     % Bypass ratio (fixed for this study)
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

% Convert R from kJ/(kg*K) to J/(kg*K) for velocity calculations
RkJ = R;          % in kJ/(kg*K)
R = R * 1000;     % in J/(kg*K)

%% Define Parameter Ranges
rfVec = linspace(1.5, 2.2, 20);  % Fan pressure ratio: 20 points from 1.5 to 2.2
rcVec = linspace(20, 28, 20);    % Compressor pressure ratio: 20 points from 20 to 28

% Create a meshgrid for the parametric study
[RF, RC] = meshgrid(rfVec, rcVec);

% Preallocate a matrix for overall efficiency
eta_0_mat = zeros(size(RF));

%% Compute Constant Inlet Conditions
% Inlet total temperature and pressure (using flight Mach number)
T0a = Ta * (1 + (k(1)-1)/2 * M_f^2);
P0a = Pa * (1 + (k(1)-1)/2 * M_f^2)^(k(1)/(k(1)-1));

% Adjusted diffuser conditions
T02s = eta_d * (T0a - Ta) + Ta;
P02  = Pa * (T02s/Ta)^(k(1)/(k(1)-1));

% Specific heats [kJ/(kg*K)]
Cp_inlet = (k(1)/(k(1)-1)) * RkJ;  % For inlet, compressor, fan
Cp_core  = (k(2)/(k(2)-1)) * RkJ;  % For burner, turbine, nozzle

% Freestream velocity [m/s]
U = M_f * sqrt(k(1) * R * Ta);

%% Parametric Study: Nested Loops over rf and rc
for i = 1:length(rfVec)
    for j = 1:length(rcVec)
        rf_local = rfVec(i);
        rc_local = rcVec(j);
        
        % Compressor stage
        T03s = T02s * rc_local^((k(1)-1)/k(1));
        T03  = (T03s - T02s)/eta_c + T02s;
        Wc   = Cp_inlet * (T03 - T02s);
        P03  = rc_local * P02;
        
        % Burner: Fuel-to-air ratio (energy balance)
        Fb = (T04_max - T03) / ((eta_b*dhc)/Cp_core - T04_max);
        
        % Burner pressure drop
        P04 = rb * P03;
        
        % Fan stage (bypass)
        T07s = T02s * rf_local^((k(1)-1)/k(1));
        T07  = (T07s - T02s)/eta_f + T02s;
        Wf   = beta * Cp_inlet * (T07 - T02s);
        
        % Turbine: Drives compressor and fan
        T05  = T04_max - (Wc + Wf) / ((1+Fb)*Cp_core);
        T05s = T04_max - (T04_max - T05)/eta_t;
        P05  = P04 * (T05s/T04_max)^(k(2)/(k(2)-1));
        
        % Core nozzle expansion to ambient
        T6s  = T05 / (P05/Pa)^((k(2)-1)/k(2));
        T6   = T05 - eta_cn * (T05 - T6s);
        Me   = sqrt((2*(T05/T6 - 1))/(k(2)-1));
        
        % Bypass nozzle expansion to ambient
        P07  = rf_local * P02;
        T8s  = T07 / (P07/Pa)^((k(1)-1)/k(1));
        T8   = T07 - eta_fn * (T07 - T8s);
        M8   = sqrt((2*(T07/T8 - 1))/(k(1)-1));
        
        % Exit velocities [m/s]
        Ue_hot  = Me * sqrt(k(2)*R*T6);
        Ue_cold = M8 * sqrt(k(1)*R*T8);
        
        % Efficiency Calculations
        % Thermal efficiency
        numerator   = (1+Fb)*Ue_hot^2 + beta*Ue_cold^2 - (beta+1)*U^2;
        denominator = 2*Fb*dhc*1000;
        eta_th_val  = numerator / denominator;
        
        % Propulsive efficiency
        thrust_core   = (1+Fb)*Ue_hot - U;
        thrust_bypass = beta*(Ue_cold - U);
        numerator_p   = 2 * (thrust_core + thrust_bypass) * U;
        denominator_p = (1+Fb)*Ue_hot^2 + beta*Ue_cold^2 - (beta+1)*U^2;
        eta_p_val     = numerator_p / denominator_p;
        
        % Overall efficiency
        eta_0_val = eta_th_val * eta_p_val;
        
        % Store result in matrix
        eta_0_mat(j,i) = eta_0_val;
    end
end

%% Save Data to CSV File
% Convert meshgrid data to column vectors for table creation
dataTable = table(RF(:), RC(:), eta_0_mat(:), ...
    'VariableNames', {'FanPressureRatio', 'CompPressureRatio', 'OverallEfficiency'});
writetable(dataTable, 'HW8_Q2_data.csv');
disp('Data saved to HW8_Q2_data.csv');

%% Create 3D Surface Plot and Save as PDF
figure;
surf(RF, RC, eta_0_mat);
xlabel('Fan Pressure Ratio (r_f)');
ylabel('Compressor Pressure Ratio (r_c)');
zlabel('\eta_0');
title('Overall Efficiency as a Function of r_f and r_c');
shading interp;  % Smooth color transitions
colorbar;
grid on;

% Save the figure as a PDF file
exportgraphics(gcf, 'HW8_Q2_surface_plot.pdf', 'ContentType','vector');
disp('Surface plot saved to HW8_Q2_surface_plot.pdf');