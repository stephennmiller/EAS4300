clc; clear; close all;

%% Given Parameters
M_f = 0.8;
Ta = 223.252;                       % Ambient temperature [K]
Pa = 26500;                         % Ambient pressure [Pa]
R = 0.287;                          % Gas constant [kJ/(kg*K)]
T04_max = 1500;                     % Maximum turbine inlet temperature [K]
dhc = 43000;                        % Fuel heat of combustion [kJ/kg]
f_stoic = 0.06;                     % Stoichiometric fuel-to-air ratio

%% Pressure Ratios
rc = 24;       % Compressor pressure ratio
rb = 0.97;     % Burner pressure loss factor
rf = 2;        % Fan pressure ratio

%% Efficiencies
eta_d = 0.94;  % Diffuser efficiency
eta_c = 0.87;  % Compressor efficiency
eta_b = 0.98;  % Burner efficiency
eta_t = 0.85;  % Turbine efficiency
eta_f = 0.92;  % Fan efficiency
eta_cn = 0.97; % Core nozzle efficiency
eta_fn = 0.98; % Fan nozzle efficiency

%% Specific Heat Ratios for Engine Sections
% k(1): Inlet, diffuser, compressor, fan
% k(2): Burner, turbine, core nozzle
k = [1.4 1.35];

%% Bypass Ratio Range and Preallocation
beta = linspace(2,7.26);
n = length(beta);
num = length(k);
Cp = zeros(num,1);
Wf = zeros(n,1);
T05 = zeros(n,1);
T05s = zeros(n,1);
P05 = zeros(n,1);
T6s = zeros(n,1);
T6 = zeros(n,1);
T06 = zeros(n,1);
Me = zeros(n,1);
T08 = zeros(n,1);
T8s = zeros(n,1);
T8 = zeros(n,1);
M8 = zeros(n,1);
Ue_hot = zeros(n,1);
Ue_cold = zeros(n,1);
I = zeros(n,1);
TSFC = zeros(n,1);
eta_th = zeros(n,1);
eta_p = zeros(n,1);
eta_0 = zeros(n,1);

%% Calculate Specific Heats [kJ/(kg*K)]
for j = 1:num
    Cp(j) = (k(j)/(k(j)-1)) * R;  
end

%% Main Engine Cycle Calculations Over Bypass Ratio
for i = 1:n
    % Inlet (diffuser) conditions: stagnation temperature and pressure
    T0a = Ta * (1 + ((k(1)-1)/2)*M_f^2);
    P0a = Pa * (1 + ((k(1)-1)/2)*M_f^2)^(k(1)/(k(1)-1));
    T02 = T0a;
    T02s = eta_d * (T02 - Ta) + Ta;
    P02 = Pa * (T02s/Ta)^(k(1)/(k(1)-1));
    
    % Compressor Stage: Ideal and actual outlet temperature
    T03s = T02 * rc^((k(1)-1)/k(1));
    T03 = (T03s - T02)/eta_c + T02;
    Wc = Cp(1) * (T03 - T02);
    P03 = rc * P02;
    
    % Burner: Compute fuel-to-air ratio based on energy balance
    Fb_calc = (T04_max - T03) / ((eta_b * dhc / Cp(2)) - T04_max);
    Fb = min(Fb_calc, f_stoic);
    
    % Burner pressure loss
    P04 = rb * P03;
    
    % Fan Stage: Ideal and actual outlet temperature for bypass stream
    T07s = T02 * rf^((k(1)-1)/k(1));
    T07 = (T07s - T02)/eta_f + T02;
    Wf(i) = beta(i) * Cp(1) * (T07 - T02);
    
    % Turbine: Energy balance for work extraction (driving compressor and fan)
    T05(i) = T04_max - (Wc + Wf(i)) / ((1 + Fb) * Cp(2));
    T05s(i) = T04_max - (T04_max - T05(i)) / eta_t;
    P05(i) = P04 * (T05s(i)/T04_max)^(k(2)/(k(2)-1));
    
    % Core Nozzle: Expansion to ambient conditions
    P6 = Pa;
    T6s(i) = T05(i) / (P05(i)/P6)^((k(2)-1)/k(2));
    T6(i) = T05(i) - eta_cn * (T05(i) - T6s(i));
    T06(i) = T05(i);
    Me(i) = sqrt((2*(T06(i)/T6(i) - 1))/(k(2)-1));
    
    % Bypass Nozzle: Expansion to ambient conditions
    T08 = T07;
    P07 = rf * P02;
    P8 = Pa;
    T8s(i) = T07 / (P07/P8)^((k(1)-1)/k(1));
    T8(i) = T07 - eta_fn * (T07 - T8s(i));
    M8(i) = sqrt((2*(T08/T8(i) - 1))/(k(1)-1));
    
    % Exit velocities [m/s]
    Ue_hot(i) = Me(i) * sqrt(k(2)*R*1000*T6(i));
    Ue_cold(i) = M8(i) * sqrt(k(1)*R*1000*T8(i));
    
    % Freestream velocity [m/s]
    U = M_f * sqrt(k(1)*R*1000*Ta);
    
    % Specific Thrust and TSFC
    I(i) = (1 + Fb) * Ue_hot(i) - U + beta(i)*(Ue_cold(i) - U);
    TSFC(i) = Fb / I(i);
    
    % Efficiency Calculations
    eta_th(i) = ((1 + Fb)*Ue_hot(i)^2 + beta(i)*Ue_cold(i)^2 - (beta(i)+1)*U^2) / (2*Fb*dhc*1000);
    eta_p(i) = 2*(beta(i)*(Ue_cold(i) - U) + ((1 + Fb)*Ue_hot(i) - U)) / ...
               ((1 + Fb)*Ue_hot(i)^2 + beta(i)*Ue_cold(i)^2 - (beta(i)+1)*U^2)*U;
    eta_0(i) = eta_th(i) * eta_p(i);
end

%% Save Data to CSV
dataTable = table(beta(:), I(:), TSFC(:), eta_th(:), eta_p(:), eta_0(:), Me(:), M8(:), ...
    'VariableNames',{'BypassRatio','SpecificThrust','TSFC','ThermalEff','PropEff','OverallEff','CoreMach','BypassMach'});
writetable(dataTable, 'HW8_Q1_data.csv');
disp('Data saved to HW8_Q1_data.csv');

%% Set up figure properties for small plots (6x4 inches)
paperWidth  = 6;
paperHeight = 4;
figPos = [1 1 paperWidth paperHeight];

% Delete old PDF if it exists
if exist('HW8_Q1_plots.pdf','file')
    delete('HW8_Q1_plots.pdf')
end

%% Create and save Figure 1: Specific Thrust vs. Bypass Ratio
fig1 = figure;
set(fig1, 'Units', 'inches', 'Position', figPos);
plot(beta, I, 'LineWidth',1.5)
xlabel('\beta','FontWeight','bold')
ylabel('Specific Thrust (m/s)','FontWeight','bold')
title('Specific Thrust vs. Bypass Ratio')
grid on
exportgraphics(gcf, 'HW8_Q1_plots.pdf', 'ContentType','vector'); % first page

%% Create and append Figure 2: TSFC vs. Bypass Ratio
fig2 = figure;
set(fig2, 'Units', 'inches', 'Position', figPos);
plot(beta, TSFC, 'LineWidth',1.5)
xlabel('\beta','FontWeight','bold')
ylabel('TSFC','FontWeight','bold')
title('TSFC vs. Bypass Ratio')
grid on
exportgraphics(gcf, 'HW8_Q1_plots.pdf', 'Append',true, 'ContentType','vector');

%% Create and append Figure 3: Efficiencies vs. Bypass Ratio
fig3 = figure;
set(fig3, 'Units', 'inches', 'Position', figPos);
plot(beta, eta_th, 'LineWidth',1.5)
hold on
plot(beta, eta_p, 'LineWidth',1.5)
plot(beta, eta_0, 'LineWidth',1.5)
hold off
xlabel('\beta','FontWeight','bold')
ylabel('Efficiency','FontWeight','bold')
title('Efficiencies vs. Bypass Ratio')
legend('\eta_{th}','\eta_{p}','\eta_{0}','Location','best')
grid on
exportgraphics(gcf, 'HW8_Q1_plots.pdf', 'Append',true, 'ContentType','vector');

%% Create and append Figure 4: Mach Numbers vs. Bypass Ratio (Core vs. Bypass)
fig4 = figure;
set(fig4, 'Units', 'inches', 'Position', figPos);
plot(beta, Me, 'LineWidth',1.5)
hold on 
plot(beta, M8, 'LineWidth',1.5)
hold off
xlabel('\beta','FontWeight','bold')
ylabel('Mach Number','FontWeight','bold')
title('Mach Numbers vs. Bypass Ratio')
legend('Core','Bypass','Location','best')
grid on
exportgraphics(gcf, 'HW8_Q1_plots.pdf', 'Append',true, 'ContentType','vector');

disp('Four figures saved to HW8_Q1_plots.pdf (each on a separate page).');
