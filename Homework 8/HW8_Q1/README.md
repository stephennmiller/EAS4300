# HW8_Q1 - Turbofan Engine Performance Analysis

MATLAB analysis for EAS 4300 Homework #8, Question 1: choosing the bypass ratio of an unmixed-exhaust turbofan for a commercial aircraft cruising at Mach 0.8 at 10,000 m. Specific thrust, TSFC, thermal, propulsive, and overall efficiency, and both exit Mach numbers are computed vs. bypass ratio $\beta$.

## Problem Statement

Cruise at $M = 0.8$, 10,000 m. Fuel $\Delta h_c = 43{,}000$ kJ/kg, $f_{st} = 0.06$. $\gamma = 1.4$ in the fan and in the core up to the burner, $\gamma = 1.35$ in the burner, turbine, and core nozzle; $R = 287$ J/kg·K. $T_{04,max} = 1500$ K, $r_c = 24$, $r_f = 2.0$. Efficiencies: $\eta_d = 0.94$, $\eta_c = 0.87$, $\eta_f = 0.92$, $\eta_b = 0.98$, $r_b = 0.97$, $\eta_t = 0.85$, $\eta_{cn} = 0.97$, $\eta_{fn} = 0.98$. Both streams are ideally expanded.

**Is there an optimum bypass ratio that maximizes the overall efficiency?**

## Results

- **Yes:** $\eta_0$ peaks at **0.281 at $\beta \approx 5.88$**. Specific thrust peaks at the same bypass ratio (1218 m/s), and TSFC is lowest there (1.98e-5 s/m).
- Below the optimum, adding bypass air lowers the jet velocity and raises propulsive efficiency. Above it, the turbine has to extract so much work for the fan that the core jet collapses: the core exit Mach number falls from 1.74 at $\beta = 2$ to 0.12 at $\beta = 7.26$.
- $\beta = 7.26$ is the practical upper limit of the plots: beyond about 7.29 the turbine cannot drive the fan and still leave the core flow enough pressure to expand to ambient.
- The bypass exit Mach number is 1.33 for every $\beta$, since the fan stream does not depend on $\beta$.
- Fuel-to-air ratio: $f = 0.0241$, independent of $\beta$.

## Method Notes

- Burner energy balance: $c_{p1}T_{03} + f\,\eta_b\,\Delta h_c = (1+f)\,c_{p2}T_{04}$, with $c_{p1}$ ($\gamma = 1.4$) for the incoming air and $c_{p2}$ ($\gamma = 1.35$) for the products, the same treatment used in Homework 6 and 7. The original 2022 version used $c_{p2}$ for both, which gave $\eta_0 = 0.302$ at the same optimum $\beta$.
- The turbine drives both the compressor and the fan: $(1+f)\,c_{p2}(T_{04} - T_{05}) = c_{p1}(T_{03} - T_{02}) + \beta\,c_{p1}(T_{07} - T_{02})$.
- The published report shows every equation above the code that uses it.

## Files

- **HW8_Q1.m** – MATLAB script, organized into sections for MATLAB `publish`.
- **HW8_Q1_data.csv** – one row per $\beta$: `BypassRatio, SpecificThrust, TSFC, ThermalEff, PropEff, OverallEff, CoreMach, BypassMach, FuelAirRatio`.
- **HW8_Q1.pdf** – published report with the code, equations, optimum output, and four figures (specific thrust, TSFC, efficiencies, exit Mach numbers vs. $\beta$).

## How to Run

From this folder in MATLAB:

```matlab
publish('HW8_Q1.m', 'format', 'pdf', 'outputDir', pwd);
```

This runs the analysis, writes `HW8_Q1_data.csv`, and creates `HW8_Q1.pdf`.

## Acknowledgments

- Original code was created by Stephen Miller on 4/13/22.
- This project was developed for educational purposes in the EAS 4300 course.
