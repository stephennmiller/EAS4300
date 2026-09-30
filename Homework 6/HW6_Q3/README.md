# EAS 4300 Homework 6 - Question 3

## Overview
This script (`HW6_Q3.m`) solves HW6 Question 3, analyzing the performance of a ramjet engine at an altitude of 10,000 meters over a flight Mach number range from 1 to 6. The analysis includes key performance parameters such as specific thrust, thrust specific fuel consumption (TSFC), combustor exit temperature, area ratio, efficiencies, and fuel-to-air ratio. The results are computed for 1000 evenly spaced Mach numbers (`delta = 1000`) to provide smooth plots and detailed data.

### Problem Description
The ramjet operates at 10,000 meters with the following conditions:
- Ambient temperature: 223.252 K
- Ambient pressure: 26,500 Pa
- Stoichiometric fuel-air ratio: 0.06
- Maximum combustor temperature: 2600 K
- Heat of combustion: 43,000 kJ/kg
- Specific heat ratios: 1.4 (inlet), 1.33 (combustor)
- Gas constant: 0.287 kJ/kg·K

The script calculates and plots the following:
- Specific Thrust vs. Flight Mach Number
- Thrust Specific Fuel Consumption (TSFC) vs. Flight Mach Number
- Combustor Exit Temperature (T04) vs. Flight Mach Number
- Area Ratio (A_exit/A_throat) vs. Flight Mach Number
- Efficiencies (overall, thermal, propulsive) vs. Flight Mach Number
- Fuel-to-Air Ratio (f) vs. Flight Mach Number

Results are saved in `hw6_q3_results.csv`. The published report `HW6_Q3.pdf` (MATLAB `publish`) contains the code, output, and all six plots.

## Methodology
1. **Setup**: Define constants and ambient conditions at 10,000 meters. Use `linspace` to create an array of 1000 Mach numbers from 1 to 6.
2. **Calculations**:
   - Compute stagnation pressure and temperature at the inlet.
   - Compute the exit Mach number for ideal expansion back to ambient pressure (closed-form isentropic relation).
   - Calculate the area ratio (A_exit/A_throat) based on the exit Mach number.
   - Determine the fuel-to-air ratio from the burner energy balance, `Cp_air*T0a + f*Δh_c = (1+f)*Cp*T04`, using the inlet Cp (k = 1.4) for the incoming air and the combustor Cp (k = 1.33) for the products. Where reaching 2600 K would need more than the stoichiometric limit, f is held at f_st = 0.06 and T04 is solved from the same balance.
   - Compute exit temperature, exit velocity, specific thrust, TSFC, and efficiencies (thermal, propulsive, overall).
3. **Output**:
   - Save results to `hw6_q3_results.csv`.
   - Generate six plots (Figures 1–6) with grid lines for readability; each appears under its own section in the published PDF.

## Key Findings
- **Specific Thrust**: Rises from 641 m/s at M=1 to a peak of about 1103 m/s at M≈2.93, where the engine first reaches the 2600 K limit, then falls to 515 m/s at M=6.
- **TSFC**: Falls from 9.36e-5 kg/N·s at M=1 to a minimum of 5.28e-5 kg/N·s near M=3.94, then rises slightly to 5.67e-5 kg/N·s at M=6 as thrust drops faster than fuel flow.
- **Combustor Exit Temperature (T04)**: Rises from 2324 K at M=1 (fuel-limited at f_st) to the 2600 K maximum at M≈2.93, and stays there for higher Mach numbers.
- **Area Ratio (A_exit/A_throat)**: Increases from 1.0003 at M=1 to 65.68 at M=6, as the nozzle expands a much higher stagnation pressure back to ambient.
- **Efficiencies**:
  - Overall efficiency (η_0) increases from 0.074 at M=1 to 0.737 at M=6.
  - Thermal efficiency (η_th) increases from 0.144 to 0.782, staying below the ideal Brayton limit (1 − T_a/T_0a = 0.878 at M=6).
  - Propulsive efficiency (η_p) increases from 0.516 to 0.942 as the exhaust velocity approaches the flight velocity.
- **Fuel-to-Air Ratio (f)**: Held at the stoichiometric limit f_st = 0.06 from M=1 to M≈2.93, then decreases to 0.0292 at M=6 as ram heating does more of the work of reaching 2600 K.

## Plots
All six plots are in the published report, [HW6_Q3.pdf](HW6_Q3.pdf), each under its own heading:
1. Specific Thrust vs. Flight Mach Number
2. TSFC vs. Flight Mach Number
3. Combustor Exit Temperature vs. Flight Mach Number
4. Nozzle Area Ratio vs. Flight Mach Number
5. Efficiencies vs. Flight Mach Number
6. Fuel-to-Air Ratio vs. Flight Mach Number

## Files
- `HW6_Q3.m`: MATLAB script that performs the calculations and generates plots.
- `hw6_q3_results.csv`: Output data table with 1000 rows, containing M_flight, Me_t, u, u_e, f, Specific_Thrust, TSFC, T04, A_exit_A_throat, eta_th, eta_p, and eta_0.
- `HW6_Q3.pdf`: Published report (code, output, and Figures 1–6).

## Notes
- The script includes clearing commands (`clear all; close all; clc;`) to ensure a fresh start for each run.
- Grid lines are added to all plots for better readability.
- The `delta = 1000` setting provides smooth plots and detailed data; with the closed-form relations the script runs in well under a second.
- To regenerate the report, run `publish('HW6_Q3.m', 'format', 'pdf', 'outputDir', pwd);` from this folder. This also rewrites the CSV.