clc;
clear vars;
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
fannoParameter = (frictioncoeff * ductLength) / D_h;

% Fanno flow to solve for M2 (approximate due to symbolic solution limitation)
% Using iterative approach or table-based approximation (simplified here)
M2_initial = 5;  % Initial guess
fannoFunc = @(M) ((k + 1) / (2 * k)) * log((1 + ((k - 1) / 2) * M^2) / ((k + 1) * M^2 / 2)) - (1 / (k * M^2)) - fannoParameter;
M2 = fzero(fannoFunc, M2_initial);  % Numerical solution for M2

% Normal shock relations
M3 = sqrt((1 + ((k - 1) / 2) * M2^2) / (k * M2^2 - ((k - 1) / 2)));

% Stagnation pressure ratio across shock
P02_P01 = (M1 / M2) * ((1 + ((k - 1) / 2) * M2^2) / (1 + ((k - 1) / 2) * M1^2))^((k + 1) / (2 * (k - 1)));  % Fanno loss
P03_P02 = ((k + 1) * M2^2 / 2) / (1 + ((k - 1) / 2) * M2^2) * (k / (k - 1)) * ((2 * k * M2^2 - (k - 1)) / (k + 1))^(-1 / (k - 1));  % Shock loss

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
rho_f = P0_tank / (R * T_stp);  % Density of full tank
rho_e = Pb / (R * T_stp);  % Density at 1 atm
W_f = rho_f * tankVolume;  % Mass of full tank
W_e = rho_e * tankVolume;  % Mass at empty (1 atm)
W_t = W_f - W_e;  % Usable mass per tank
expelledMass = m_dot * time;
n = ceil(expelledMass / W_t);  % Round up to nearest integer

fprintf("Question 3:\n");
fprintf("P0_tank = %g Pa\n", P0_tank);
fprintf("A_t = %g m^2\n", A_t);
fprintf("Tank volume = %g m^3\n", tankVolume);
fprintf("Mass flow rate is %g kg/s\n", m_dot);
fprintf("Rho_f = %g kg/m^3\n", rho_f);
fprintf("Rho_e = %g kg/m^3\n", rho_e);
fprintf("W_f = %g kg\n", W_f);
fprintf("W_e = %g kg\n", W_e);
fprintf("W_t = %g kg\n", W_t);
fprintf("Expelled mass = %g kg\n", expelledMass);
fprintf("Number of tanks = %g\n\n", n);

% Question 4: Heat release for stoichiometric combustion of hydrogen
dhc = 241820000;  % Heat of combustion (J/kmol), converted from kJ/kmol
molar_mass_H2 = 2.016;  % kg/kmol
dhc_per_kg = dhc / molar_mass_H2;  % J/kg
air_fuel_ratio = 2.38095 * 0.21 * (32 / 2.016);  % Stoichiometric air-fuel ratio
m_dot_H2 = m_dot / air_fuel_ratio;  % Mass flow rate of H2
heat_release = m_dot_H2 * dhc_per_kg;  % Heat release rate (W)

fprintf("Question 4:\n");
fprintf("Heat of combustion = %g J/kg\n", dhc_per_kg);
fprintf("Mass flow rate of H2 = %g kg/s\n", m_dot_H2);
fprintf("Heat release = %g W\n", heat_release);
fprintf("Heat release = %g MW\n\n", heat_release / 1e6);

% Question 5: Compare stagnation pressures
P01_ideal = Pb * (1 + ((k - 1) / 2) * M1^2)^(k / (k - 1));  % Ideally expanded at M = 5

fprintf("Question 5:\n");
fprintf("Ideally expanded stagnation pressure : P0 = %g Pa\n", P01_ideal);
fprintf("Stagnation pressure to maintain shock : P0 = %g Pa\n", P01);