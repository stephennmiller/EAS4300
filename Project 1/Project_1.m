clc;
clear;
close all;

% Constants
k = 1.4;  % Specific heat ratio (gamma)
R = 287;  % Gas constant (J/kg·K)
T_stp = 293.15;  % Standard temperature (K)
Pb = 101325;  % Back pressure (Pa)
P0_tank = 4500 * 6894.76;  % Tank pressure (Pa)
frictioncoeff = 0.02;  % Friction coefficient
time = 120;  % Time in seconds (2 minutes)
tankVolume = 49 / 1000;  % Tank volume (m^3)

% Question 1: Determine the A_nozzle exit / A* for Mach 5 flow
M1 = 5;  % Design Mach number at nozzle exit
A_e = (45^2) * (1 / 1000)^2;  % Exit area (m^2)
Ratio = (1 / M1) * (((2 / (k + 1)) * (1 + ((k - 1) / 2) * M1^2)))^((k + 1) / (2 * (k - 1)));  % Corrected area ratio
A_t = A_e / Ratio;  % Throat area (m^2)

fprintf("Question 1:\n");
fprintf("The desired Mach number at the exit is : M = %d\n", M1);
fprintf("The exit area is : A_e = %g m^2\n", A_e);
fprintf("Area ratio = %g\n\n", Ratio);

% Question 2: Stagnation pressure needed for normal shock at exit plane
ductLength = 225 / 1000;  % Channel length (m)
% Hydraulic diameter for square duct (45 mm x 45 mm)
D_h = (4 * A_e) / (4 * (45 / 1000));  % Correct hydraulic diameter
% The handout gives f = 0.02 without saying whether the Fanno equation is
% written with f*L/D or 4f*L/D (Fanning). Set the multiplier to match the course.
fannoMultiplier = 1;  % 1 -> f*L/D, 4 -> 4f*L/D
fannoParameter = fannoMultiplier * (frictioncoeff * ductLength) / D_h;

% Fanno function: f*L*/D from Mach M to the choking point (M = 1)
fannoLstar = @(M) (1 - M^2) / (k * M^2) + ((k + 1) / (2 * k)) * log((k + 1) * M^2 / (2 + (k - 1) * M^2));

% Supersonic flow slows toward M = 1 along the duct: fL*/D(M2) = fL*/D(M1) - fL/D
fannoTarget = fannoLstar(M1) - fannoParameter;
if fannoTarget <= 0
    error('Duct is longer than L* for M1 = %g; the flow chokes inside the channel.', M1);
end
M2 = fzero(@(M) fannoLstar(M) - fannoTarget, [1 + 1e-6, M1]);  % Supersonic branch

% Normal shock relations
M3 = sqrt((1 + ((k - 1) / 2) * M2^2) / (k * M2^2 - ((k - 1) / 2)));

% Stagnation pressure ratio across shock
P02_P01 = (M1 / M2) * ((1 + ((k - 1) / 2) * M2^2) / (1 + ((k - 1) / 2) * M1^2))^((k + 1) / (2 * (k - 1)));  % Fanno loss
P03_P02 = (((k + 1) * M2^2 / 2) / (1 + ((k - 1) / 2) * M2^2))^(k / (k - 1)) * ((2 * k * M2^2 - (k - 1)) / (k + 1))^(-1 / (k - 1));  % Shock loss

% Pressures
P2 = Pb / (1 + (2 * k / (k + 1)) * (M2^2 - 1));  % Pressure before shock
P02 = P2 * (1 + ((k - 1) / 2) * M2^2)^(k / (k - 1));  % Stagnation pressure before shock
P03 = Pb * (1 + ((k - 1) / 2) * M3^2)^(k / (k - 1));  % Stagnation pressure after shock
P01 = P03 / (P03_P02 * P02_P01);  % Total stagnation pressure at nozzle inlet

fprintf("Question 2:\n");
fprintf("Diameter of the channel is : %g m\n", D_h);
fprintf("The length of the channel is : %g m\n", ductLength);
fprintf("Fanno Parameter is : %g\n", fannoParameter);
fprintf("The Mach number at the end of the channel is : %g\n", M2);
fprintf("The pressure before shock is : %g Pa\n", P2);
fprintf("P02 = %g Pa\n", P02);
fprintf("P03 = %g Pa\n", P03);
fprintf("M3 is %g\n", M3);
fprintf("P03 / P02 = %g\n\n", P03 / P02);
fprintf("P01 = %g Pa\n\n", P01);

% Question 3: Number of tanks required for 2 minutes
% Mass flow rate based on corrected P01
m_dot = A_t * (P01 / sqrt(T_stp)) * sqrt(k / R) * ((2 / (k + 1))^((k + 1) / (2 * (k - 1))));
expelledMass = m_dot * time;

% Each tank feeds the regulator through a choked 1/8" throat, and the tank
% gas is treated as isothermal at T_stp. A tank stops being useful once
% either (a) its throat no longer chokes against the regulated P01, or
% (b) the n tanks together can no longer pass m_dot through their throats.
tankThroatDiameter = 0.125 * 0.0254;  % m
A_tankThroat = pi * tankThroatDiameter^2 / 4;  % m^2
chokedFlux = sqrt(k / R) * ((2 / (k + 1))^((k + 1) / (2 * (k - 1)))) / sqrt(T_stp);  % kg/s per (Pa*m^2)
P_chokeMin = P01 * ((k + 1) / 2)^(k / (k - 1));  % Lowest tank pressure that still chokes the throat
W_f = P0_tank * tankVolume / (R * T_stp);  % Mass of full tank
if P_chokeMin >= P0_tank
    error('Tanks at %g Pa cannot choke their throats against P01 = %g Pa; no usable air.', P0_tank, P01);
end

n = 0;
usableMass = 0;
while usableMass < expelledMass
    n = n + 1;
    P_flowMin = m_dot / (n * A_tankThroat * chokedFlux);  % Tank pressure where n throats pass m_dot
    P_end = max(P_chokeMin, P_flowMin);  % Tank pressure at end of run
    usableMass = n * (P0_tank - P_end) * tankVolume / (R * T_stp);
end

fprintf("Question 3:\n");
fprintf("P0_tank = %g Pa\n", P0_tank);
fprintf("A_t = %g m^2\n", A_t);
fprintf("Tank volume = %g m^3\n", tankVolume);
fprintf("Mass flow rate is %g kg/s\n", m_dot);
fprintf("Expelled mass = %g kg\n", expelledMass);
fprintf("W_f = %g kg\n", W_f);
fprintf("Tank throat area = %g m^2\n", A_tankThroat);
fprintf("Min tank pressure to choke the throat = %g Pa\n", P_chokeMin);
fprintf("Min tank pressure for %d throats to pass m_dot = %g Pa\n", n, P_flowMin);
fprintf("Tank pressure at end of run = %g Pa\n", P_end);
fprintf("Usable mass with %d tanks = %g kg\n", n, usableMass);
fprintf("Number of tanks = %g\n", n);

% Cross-check (2022 method): flow through one choked 1/8" tank throat at
% full tank pressure, run for 2 minutes, divided by one tank's usable mass
% down to 1 atm.
m_dot_tankThroat = A_tankThroat * P0_tank * chokedFlux;
W_t_atm = (P0_tank - Pb) * tankVolume / (R * T_stp);
n_throatMethod = ceil(m_dot_tankThroat * time / W_t_atm);
fprintf("Cross-check: tank-throat flow at 4500 psi = %g kg/s\n", m_dot_tankThroat);
fprintf("Cross-check: number of tanks = %g\n\n", n_throatMethod);

% Question 4: Heat release for stoichiometric combustion of hydrogen
dhc = 241820000;  % Heat of combustion (J/kmol), converted from kJ/kmol
molar_mass_H2 = 2.016;  % kg/kmol
dhc_per_kg = dhc / molar_mass_H2;  % J/kg
% H2 + 0.5 (O2 + 3.76 N2) -> H2O + 1.88 N2
air_fuel_ratio = 0.5 * (32 + 3.76 * 28.013) / molar_mass_H2;  % Stoichiometric air-fuel ratio (kg air / kg H2)
m_dot_H2 = m_dot / air_fuel_ratio;  % Mass flow rate of H2
heat_release = m_dot_H2 * dhc_per_kg;  % Heat release rate (W)

% H2 + 2.38095 (0.79 N2 + 0.21 O2) -> H2O + 1.88095 N2; only H2O has a
% nonzero heat of formation, so dhc = -h_f(H2O) = 241,820 kJ/kmol H2.
fprintf("Question 4:\n");
fprintf("Heat release = %g kJ/kmol H2\n", dhc / 1000);
fprintf("Heat of combustion = %g J/kg\n", dhc_per_kg);
fprintf("Mass flow rate of H2 = %g kg/s\n", m_dot_H2);
fprintf("Heat release rate at the Q2 mass flow = %g MW\n\n", heat_release / 1e6);

% Question 5: Compare stagnation pressures
% Ideally expanded: the flow leaves the channel at M2 with P2 = Pb (no
% shock), with the same Fanno friction loss through the channel as in Q2.
P02_ideal = Pb * (1 + ((k - 1) / 2) * M2^2)^(k / (k - 1));
P01_ideal = P02_ideal / P02_P01;

fprintf("Question 5:\n");
fprintf("Ideally expanded stagnation pressure : P0 = %g Pa\n", P01_ideal);
fprintf("Stagnation pressure to maintain shock : P0 = %g Pa\n", P01);
fprintf("Ideally expanded / shock case = %g\n", P01_ideal / P01);