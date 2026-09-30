# EAS 4300 HW8 Q2 - Optimizing the Compressor and Fan Pressure Ratios

This question uses the turbofan model from [HW8_Q1](../HW8_Q1/HW8_Q1.m) at its optimum bypass ratio, $\beta = 5.88$, and varies the fan pressure ratio $r_f$ and the compressor pressure ratio $r_c$ to see whether a combination other than $r_c = 24$, $r_f = 2$ gives a higher overall efficiency $\eta_0$.

## Problem Statement

Using the model from problem 1 at the optimum bypass ratio, vary $1.5 < r_f < 2.2$ and $20 < r_c < 28$ in a single 900-row table: $r_f$ uses the "repeat pattern" across 30 rows and $r_c$ uses the "apply pattern" across 30 rows. Plot $\eta_0(r_c, r_f)$ and determine whether a different combination of $r_c$ and $r_f$ increases the overall efficiency.

## Results

- **Yes.** The highest overall efficiency in the table is **$\eta_0 = 0.287$ at $r_c = 28$, $r_f \approx 1.89$**, compared with 0.281 at $r_c = 24$, $r_f = 2$.
- $\eta_0$ keeps rising with $r_c$ all the way to the edge of the range, so the best $r_c$ is limited by the range studied (28), not by the cycle.
- For every $r_c$, $\eta_0$ has a clear optimum fan pressure ratio between about $r_f = 1.89$ and $1.93$. A higher $r_f$ takes too much turbine work from the core, and a lower one gives too little bypass thrust.

## Table Layout

`HW8_Q2_data.csv` has 900 rows. $r_f$ steps through its 30 values and then repeats; $r_c$ holds each of its 30 values for 30 consecutive rows. Columns: `FanPressureRatio, CompPressureRatio, OverallEfficiency, SpecificThrust, TSFC`.

## Files

- **HW8_Q2.m** – MATLAB script, organized into sections for MATLAB `publish`.
- **HW8_Q2_data.csv** – the 900-row table described above.
- **HW8_Q2.pdf** – published report with the code, equations, the best combination, and two figures: a surface plot and a contour plot of $\eta_0(r_c, r_f)$ with the optimum and the $r_c = 24$, $r_f = 2$ design marked.

## How to Run

From this folder in MATLAB:

```matlab
publish('HW8_Q2.m', 'format', 'pdf', 'outputDir', pwd);
```

This runs the analysis, writes `HW8_Q2_data.csv`, and creates `HW8_Q2.pdf`.

## Method Notes

- Same cycle and burner energy balance as Question 1, with $T_{02} = T_{0a}$ as the compressor and fan inlet temperature. An earlier revision of this script used the isentropic-equivalent $T_{02s}$ there, which is only used to find $P_{02}$.

## Acknowledgments

- Original code was created by Stephen Miller.
- This project was developed for educational purposes in the EAS 4300 course.
