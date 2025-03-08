%% EAS 4300 HW#7 Question 2
% Created by Stephen Miller on 3/28/22.
% Modified to produce multiple figures in a single PDF and fix minor logic issue.
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

clc; clear; clear vars;

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
rb            = 0.97;       % burner pressure ratio (P03 / P02)
rab           = 0.97;       % afterburner pressure ratio (P06 / P05)
f_st          = 0.06;       % stoichiometric fuel-to-air ratio
k             = [1.4, 1.3]; % gamma values: 1.4 (up to burner), 1.3 (rest)
R             = 0.287;      % [kJ/kg·K], gas constant
Ta            = 216.65;     % [K], ambient temperature at 15 km
Pa            = 1.2112e4;   % [Pa], ambient pressure at 15 km
delta         = 100;        % number of data points
afterbuner    = 1;          % flag to run afterburner logic (kept spelling)

%% Preallocate
rc       = linspace(2, 60, delta);  % compressor pressure ratio range
n        = length(rc);
num      = length(k);
Cp       = zeros(num,1);

% Non-afterburner arrays
T03s     = zeros(n,1);
T03      = zeros(n,1);
W        = zeros(n,1);
f_b      = zeros(n,1);
P04      = zeros(n,1);
T05      = zeros(n,1);
T05s     = zeros(n,1);
P05      = zeros(n,1);
T7s      = zeros(n,1);
T7       = zeros(n,1);
Me       = zeros(n,1);
Ue       = zeros(n,1);
I        = zeros(n,1);
TSFC     = zeros(n,1);
eta_th   = zeros(n,1);
eta_p    = zeros(n,1);
eta_0    = zeros(n,1);
Ratio    = zeros(n,1);

% Afterburner arrays
f_ab     = zeros(n,1);
f        = zeros(n,1);  % total fuel (main + afterburner)
P06_ab   = zeros(n,1);
T7_as    = zeros(n,1);
T7_ab    = zeros(n,1);
M_e_ab   = zeros(n,1);
u_e_ab   = zeros(n,1);
I_ab     = zeros(n,1);
TSFC_ab  = zeros(n,1);
eta_p_ab = zeros(n,1);
eta_th_ab= zeros(n,1);
eta_0_ab = zeros(n,1);
Ratio_ab = zeros(n,1);

%% Inlet / Freestream
T0a = Ta * (1 + ((k(1) - 1)/2) * M^2);                     
P0a = Pa * (1 + ((k(1) - 1)/2) * M^2)^(k(1)/(k(1) - 1));   
u   = M * sqrt(R * 1000 * k(1) * Ta);

% Ideal diffuser (for T02s), then real P02
T02s = eta_d * (T0a - Ta) + Ta;
P02  = Pa * (T02s / Ta)^(k(1)/(k(1) - 1));

% Cp calculations
for i = 1:num
    Cp(i) = (k(i)/(k(i) - 1)) * R;
end

%% ---- Non-Afterburner Performance ----
for j = 1:n
    % Compressor
    T03s(j) = T0a * rc(j)^((k(1) - 1)/k(1));
    T03(j)  = ((T03s(j) - T0a) / eta_c) + T0a;
    W(j)    = Cp(1) * (T03(j) - T0a);

    % Main burner fuel requirement
    f_b(j) = ((T04_max/T03(j)) - 1) / ...
             (eta_b * delta_h_c/(Cp(2)*T03(j)) - (T04_max/T03(j)));

    % Pressures
    P04(j)  = rb * rc(j) * P02;

    % Turbine exit (T05 = T06 = T07 in a simple turbojet)
    T05(j)  = T04_max - (W(j)/(Cp(2)*(1 + f_b(j))));
    T05s(j) = ((T05(j) - T04_max)/eta_t) + T04_max;
    P05(j)  = P04(j) * (T05s(j)/T04_max)^(k(2)/(k(2) - 1));

    % Nozzle expansion to ambient
    T7s(j) = T05(j) / (P05(j)/Pa)^((k(2) - 1)/k(2));
    T7(j)  = eta_n * (T7s(j) - T05(j)) + T05(j);
    Me(j)  = sqrt(((T05(j)/T7(j)) - 1)/((k(2) - 1)/2));
    Ue(j)  = Me(j)*sqrt(k(2)*R*1000*T7(j));

    % Thrust, TSFC, efficiencies
    I(j)      = (1 + f_b(j))*Ue(j) - u;
    TSFC(j)   = f_b(j) / I(j);
    eta_p(j)  = (I(j)*u) / ((1 + f_b(j))*(Ue(j)^2)/2 - (u^2)/2);
    eta_th(j) = ((1 + f_b(j))*Ue(j)^2 - u^2)/(2*f_b(j)*delta_h_c*1000);
    eta_0(j)  = eta_th(j) * eta_p(j);

    % Nozzle area ratio
    Ratio(j) = (1/Me(j)) * ...
        ((2/(k(2)+1)) * (1 + ((k(2) - 1)/2)*Me(j)^2))^((k(2)+1)/(2*(k(2)-1)));
end

%% ---- Afterburner Performance ----
if afterbuner == 1
    for j = 1:n
        % Recompute T03, W, etc., for the same rc
        T03s(j) = T0a * rc(j)^((k(1) - 1)/k(1));
        T03(j)  = ((T03s(j) - T0a) / eta_c) + T0a;
        W(j)    = Cp(1) * (T03(j) - T0a);

        % Main burner with efficiency (CHANGED HERE to match above logic)
        f_b(j) = ((T04_max/T03(j)) - 1) / ...
                 (eta_b * delta_h_c/(Cp(2)*T03(j)) - (T04_max/T03(j)));

        % Turbine exit
        P04(j)  = rb * rc(j) * P02;
        T05(j)  = T04_max - (W(j)/(Cp(2)*(1 + f_b(j))));
        T05s(j) = ((T05(j) - T04_max)/eta_t) + T04_max;
        P05(j)  = P04(j) * (T05s(j)/T04_max)^(k(2)/(k(2) - 1));

        % Afterburner fuel
        % f_ab formula references (1 + f_b) mass flow and T06_ab_max
        f_ab(j) = ((1 + f_b(j))*(T06_ab_max - T05(j))) / ...
                  (((eta_ab*delta_h_c)/Cp(2)) - T06_ab_max);

        % Total fuel
        f(j) = f_b(j) + f_ab(j); 
        % (If needed, you could clamp f(j) <= f_st here.)

        % Afterburner exit
        P06_ab(j) = rab * P05(j);
        T7_as(j)  = T06_ab_max / (P06_ab(j)/Pa)^((k(2) - 1)/k(2));
        T7_ab(j)  = eta_n * (T7_as(j) - T06_ab_max) + T06_ab_max;
        M_e_ab(j) = sqrt(((T06_ab_max/T7_ab(j)) - 1)/((k(2) - 1)/2));
        u_e_ab(j) = M_e_ab(j)*sqrt(k(2)*R*1000*T7_ab(j));

        % Afterburning thrust & TSFC
        I_ab(j)      = (1 + f(j))*u_e_ab(j) - u;
        TSFC_ab(j)   = f(j) / I_ab(j);
        eta_th_ab(j) = ((1 + f(j))*u_e_ab(j)^2 - u^2) / ...
                       (2*f(j)*delta_h_c*1000);
        eta_p_ab(j)  = (I_ab(j)*u) / ((1 + f(j))*(u_e_ab(j)^2)/2 - (u^2)/2);
        eta_0_ab(j)  = eta_th_ab(j)*eta_p_ab(j);

        % Nozzle area ratio (with afterburner)
        Ratio_ab(j) = (1/M_e_ab(j)) * ...
            ((2/(k(2)+1)) * (1 + ((k(2)-1)/2)*M_e_ab(j)^2))^((k(2)+1)/(2*(k(2)-1)));
    end
end

%% ---- Plotting (Four Figures -> Single PDF) ----
% Requires MATLAB R2020a+ for 'Append' option

% Figure 1: Specific Thrust
figure('Position', [100, 100, 600, 400]);
plot(rc, I, 'DisplayName', 'Non-afterburner', 'LineWidth', 1.5);
hold on
plot(rc, I_ab, '--', 'Color', 'b', 'DisplayName', 'Afterburner', 'LineWidth', 1.5);
hold off
ylim([300 1200]);  % Adjust as desired
legend('Location','southwest');
xlabel('r_c [dim]', 'FontWeight','bold');
ylabel('I [m/s]', 'FontWeight','bold');
title('Specific Thrust vs. r_c');
grid on;
exportgraphics(gcf, 'HW7_Q2_Plots.pdf', 'Append', false);

% Figure 2: TSFC
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
exportgraphics(gcf, 'HW7_Q2_Plots.pdf', 'Append', true);

% Figure 3: Efficiencies
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
exportgraphics(gcf, 'HW7_Q2_Plots.pdf', 'Append', true);

% Figure 4: Nozzle Area Ratio
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
exportgraphics(gcf, 'HW7_Q2_Plots.pdf', 'Append', true);

%% ---- Save Key Data to CSV ----
% Collect main and afterburner data into a single table
% (Feel free to add/remove columns as needed.)

data = table( ...
    rc', ...
    f_b,       f_ab,      f, ...
    I,         I_ab, ...
    TSFC,      TSFC_ab, ...
    eta_th,    eta_th_ab, ...
    eta_p,     eta_p_ab, ...
    eta_0,     eta_0_ab, ...
    Ratio,     Ratio_ab, ...
    'VariableNames', { ...
    'rc', ...
    'f_b_nonAB','f_ab_AB','f_total', ...
    'I_nonAB',  'I_AB', ...
    'TSFC_nonAB','TSFC_AB', ...
    'eta_th_nonAB','eta_th_AB', ...
    'eta_p_nonAB','eta_p_AB', ...
    'eta_0_nonAB','eta_0_AB', ...
    'AR_nonAB','AR_AB'});

writetable(data, 'HW7_Q2_Data.csv');

%% (Optional) Display a snippet of the table
disp(data(1:5,:));
