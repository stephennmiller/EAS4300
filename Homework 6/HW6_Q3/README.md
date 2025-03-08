# EAS 4300 Homework 6 - Question 3

## Overview
This script (`q3.m`) solves HW6 Question 3, analyzing the performance of a ramjet engine at an altitude of 10,000 meters over a flight Mach number range from 1 to 6. The analysis includes key performance parameters such as specific thrust, thrust specific fuel consumption (TSFC), combustor exit temperature, area ratio, efficiencies, and fuel-to-air ratio. The results are computed for 1000 evenly spaced Mach numbers (`delta = 1000`) to provide smooth plots and detailed data.

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

Results are saved in `hw6_q3_results.csv`, and plots are saved as PDF files (`HW6_Q3_Figure1.pdf` to `HW6_Q3_Figure6.pdf`).

## Methodology
1. **Setup**: Define constants and ambient conditions at 10,000 meters. Use `linspace` to create an array of 1000 Mach numbers from 1 to 6.
2. **Calculations**:
   - Compute stagnation pressure and temperature at the inlet.
   - Solve for the exit Mach number using symbolic math (`vpasolve`).
   - Calculate the area ratio (A_exit/A_throat) based on the exit Mach number.
   - Determine the fuel-to-air ratio, constrained by the stoichiometric limit (f_st = 0.06) and maximum combustor temperature (2600 K).
   - Compute exit temperature, exit velocity, specific thrust, TSFC, and efficiencies (thermal, propulsive, overall).
3. **Output**:
   - Save results to `hw6_q3_results.csv`.
   - Generate six plots (Figures 1–6) with grid lines for readability and save them as PDF files.

## Key Findings
- **Specific Thrust**: Peaks around M=2.67 at approximately 1104 m/s and decreases to 499 m/s at M=6, reflecting the ramjet’s optimal operating range.
- **TSFC**: Decreases from 0.000103 kg/N·s at M=1 to 0.0000445 kg/N·s at M=6, indicating improved fuel efficiency at higher Mach numbers.
- **Combustor Exit Temperature (T04)**: Increases from 2357 K at M=1 to 2600 K (the maximum limit) around M=2.62 and remains constant thereafter due to the stoichiometric constraint.
- **Area Ratio (A_exit/A_throat)**: Increases exponentially from 1.0003 at M=1 to 65.68 at M=6, as expected for a ramjet to accommodate expanding exhaust gases at higher speeds.
- **Efficiencies**:
  - Overall efficiency (\eta_0) increases from 0.135 at M=1 to 0.469 at M=6.
  - Thermal efficiency (\eta_th) increases from 0.263 to 0.505.
  - Propulsive efficiency (\eta_p) increases from 0.514 to 0.930, showing better energy conversion at higher speeds.
- **Fuel-to-Air Ratio (f)**: Decreases from 0.0675 at M=1 to 0.0223 at M=6, constrained by the stoichiometric limit (f_st = 0.06) at lower Mach numbers.

## Plots
- [Specific Thrust vs. Flight Mach Number](HW6_Q3_Figure1.pdf)
- [TSFC vs. Flight Mach Number](HW6_Q3_Figure2.pdf)
- [Combustor Exit Temperature vs. Flight Mach Number](HW6_Q3_Figure3.pdf)
- [Area Ratio vs. Flight Mach Number](HW6_Q3_Figure4.pdf)
- [Efficiencies vs. Flight Mach Number](HW6_Q3_Figure5.pdf)
- [Fuel-to-Air Ratio vs. Flight Mach Number](HW6_Q3_Figure6.pdf)

These are downloadable PDF files for high-quality viewing and printing.

## Files
- `q3.m`: MATLAB script that performs the calculations and generates plots.
- `hw6_q3_results.csv`: Output data table with 1000 rows, containing M_flight, Me_t, u, u_e, f, Specific_Thrust, TSFC, T04, A_exit_A_throat, eta_th, eta_p, and eta_0.
- `HW6_Q3_Figure1.pdf`: Specific Thrust vs. Flight Mach Number.
- `HW6_Q3_Figure2.pdf`: TSFC vs. Flight Mach Number.
- `HW6_Q3_Figure3.pdf`: Combustor Exit Temperature vs. Flight Mach Number.
- `HW6_Q3_Figure4.pdf`: Area Ratio vs. Flight Mach Number.
- `HW6_Q3_Figure5.pdf`: Efficiencies vs. Flight Mach Number.
- `HW6_Q3_Figure6.pdf`: Fuel-to-Air Ratio vs. Flight Mach Number.

## Notes
- The script includes clearing commands (`clear all; close all; clc;`) to ensure a fresh start for each run.
- Grid lines are added to all plots for better readability.
- The `delta = 1000` setting provides smooth plots and detailed data, but increases computation time slightly compared to `delta = 10`.
- Plots are saved as PDFs for high-quality output, downloadable from the links above.