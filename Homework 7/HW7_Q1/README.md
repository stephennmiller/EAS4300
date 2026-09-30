# EAS4300 – HW7 Q1: Turbojet Performance Analysis

This folder contains the code and outputs for a **non-afterburning turbojet** operating at 15 km altitude and Mach $M = 1.8$. The MATLAB script computes:

- **Specific Thrust** $I$  
- **Thrust-Specific Fuel Consumption (TSFC)**  
- **Thermal**, **Propulsive**, and **Overall** Efficiencies: $\eta_{th}$, $\eta_{p}$, and $\eta_{0}$  
- **Nozzle Area Ratio**  

and plots them against the **compressor pressure ratio** $r_c$ from 2 to 60.

---

## Problem Statement

A non-afterburning turbojet is designed for 15 km altitude and Mach 1.8, with a maximum turbine inlet stagnation temperature of 1500 K. The fuel has an LHV of 43,124 kJ/kg and $f_{st} = 0.06$. Efficiencies: $\eta_d = 0.9$, $\eta_c = 0.9$, $\eta_b = 0.98$, $r_b = 0.97$, $\eta_t = 0.92$, $\eta_n = 0.98$. Use $\gamma = 1.4$ up to the burner and $\gamma = 1.3$ for the rest of the engine, with $R = 0.287$ kJ/kg·K. The exhaust is ideally expanded.

- Is there an optimum $r_c$ that minimizes TSFC?
- Is there an optimum $r_c$ that maximizes specific thrust?

---

## Results

- **Maximum specific thrust:** yes, about **708 m/s at $r_c \approx 6.1$**. Above that, the extra compressor work leaves less energy for the exhaust, and thrust falls to 421 m/s at $r_c = 60$.
- **Minimum TSFC:** yes, about **3.59e-5 s/m at $r_c \approx 43$**. TSFC falls from 5.35e-5 s/m at $r_c = 2$, reaches its minimum near $r_c = 43$, then rises slightly to 3.66e-5 s/m at $r_c = 60$. With non-ideal components, the gain in cycle efficiency at very high $r_c$ no longer outweighs the drop in thrust.
- **Efficiencies:** $\eta_{th}$ peaks at 0.51, $\eta_p$ ranges from 0.62 to 0.73, and $\eta_0$ peaks at 0.34.
- **Nozzle area ratio:** 1.74 to 3.02.
- The main burner needs $f_b$ = 0.035 at $r_c = 2$, falling to 0.015 at $r_c = 60$, well below $f_{st}$.

---

## Method Notes

- **Burner energy balance:** $C_{p,air} T_{03} + f\,\eta_b\,\Delta h_c = (1+f)\,C_{p,gas}\,T_{04}$, with $C_{p,air}$ from $\gamma = 1.4$ for the incoming air and $C_{p,gas}$ from $\gamma = 1.3$ for the products. If reaching 1500 K ever needed more than $f_{st}$, the script would burn at $f_{st}$ and solve for the lower $T_{04}$ (this never happens for $r_c$ = 2–60).
- The optimum $r_c$ values are found with `max`/`min` on the computed curves, so no toolbox is required.

---

## Files

1. [**HW7_Q1.m**](HW7_Q1.m)  
   Main MATLAB script, organized into sections for MATLAB `publish`.

2. [**HW7_Q1_Data.csv**](HW7_Q1_Data.csv)  
   One row per $r_c$: `rc, f_b, T03, T04, T05, P04, P05, T7, Me, Ue, W, I, TSFC, eta_th, eta_p, eta_0, AreaRatio`.

3. [**HW7_Q1.pdf**](HW7_Q1.pdf)  
   Published report with the code, the optimum-$r_c$ output, and all four figures:
   - Figure 1: Specific Thrust vs. $r_c$ (maximum marked)
   - Figure 2: TSFC vs. $r_c$ (minimum marked)
   - Figure 3: Efficiencies $\eta_{th}$, $\eta_{p}$, $\eta_{0}$ vs. $r_c$
   - Figure 4: Nozzle Area Ratio vs. $r_c$

---

## Usage

From this folder in MATLAB, run:

```matlab
publish('HW7_Q1.m', 'format', 'pdf', 'outputDir', pwd);
```

This runs the analysis, writes `HW7_Q1_Data.csv`, and creates `HW7_Q1.pdf`.
