# HW6_Q4 - Ramjet Performance Sensitivity Study

## Purpose
This MATLAB script (`HW6_Q4.m`) is designed to solve a ramjet propulsion problem by conducting a sensitivity study of the ramjet engine's performance. The study evaluates the impact of combustor efficiency ($\eta_b$) and exhaust nozzle total pressure ratio ($r_n$) on thrust ($I$) and thrust-specific fuel consumption (TSFC) across a range of flight Mach numbers from 1 to 6. The code computes baseline performance metrics, perturbs the key parameters, calculates normalized sensitivity coefficients, generates plots, and saves the results to a `.csv` file and a `.pdf` figure for further analysis.

The problem is based on a scenario where the ramjet operates under specific ambient conditions and uses jet fuel with given properties. The script leverages finite difference approximations to derive sensitivity coefficients and visualizes the results in two figures (the four derivatives, then the four normalized coefficients), both included in the published report `HW6_Q4.pdf`.

## Problem Statements

### Problem 2
A ramjet is to propel an aircraft at Mach 3 at high altitude where the ambient pressure is 8.5 kPa and the ambient temperature $T_0$ is 220 K. The turbine inlet temperature $T_t$ is 2540 K.
- a. The thermal efficiency,
- b. The propulsion efficiency,
- c. The overall efficiency.
Assume all components of the engine are ideal—that is, frictionless—determine the above efficiencies. The specific heat ratio by $\gamma = 1.4$ and make the approximations appropriate to $f \ll 1$.

### Problem 4
(15 points) Conduct a sensitivity study of the performance of the ramjet engine with respect to combustion efficiency ($\eta_b$) and the exhaust nozzle total pressure ratio ($r_n$). Using the conditions and properties listed for problem 2 as a starting point and jet fuel used has a heat of combustion of 43,000 kJ/kg, and a stoichiometric fuel-to-air ratio of 0.06, construct plots of

$$ \frac{d(I)}{d\eta_b}, \quad \frac{d(I)}{dr_n}, \quad \frac{d(\text{TSFC})}{d\eta_b}, \quad \frac{d(\text{TSFC})}{dr_n}, $$

as a function of flight Mach number from 1 to 6. The derivatives can be approximated using small but discrete changes, such as

$$ \frac{d(I)}{d\eta_b} \approx \frac{\Delta I}{\Delta \eta_b}. $$

Use baseline $\eta_b$ and $r_n$ values of one, and use delta values of 0.01. Make plots of the normalized sensitivity coefficient (e.g., derivative) as a function of flight Mach number:

1. $\frac{1}{I} \frac{dI}{d\eta_b}$ vs. $M_{\text{flight}}$
2. $\frac{1}{I} \frac{dI}{dr_n}$ vs. $M_{\text{flight}}$
3. $\frac{1}{\text{TSFC}} \frac{d(\text{TSFC})}{d\eta_b}$ vs. $M_{\text{flight}}$
4. $\frac{1}{\text{TSFC}} \frac{d(\text{TSFC})}{dr_n}$ vs. $M_{\text{flight}}$

(15 points) Conduct a sensitivity study with respect to the performance of the ramjet engine with respect to combustion efficiency ($\eta_b$) and the exhaust nozzle total pressure ratio ($r_n$).

## Files
- `HW6_Q4.m`: The main MATLAB script that performs the calculations, generates plots, and saves the results.
- `HW6_Q4_results.csv`: Contains the computed data including Mach number, fuel-to-air ratio ($f$), thrust ($I$), TSFC (kg/N·hr), the four sensitivity derivatives, and the four normalized sensitivity coefficients.
- `HW6_Q4.pdf`: Published report (MATLAB `publish`) with the code, output, and two 2×2 figures vs. flight Mach number: Figure 1 shows the derivatives $dI/d\eta_b$, $dI/dr_n$, $d(\text{TSFC})/d\eta_b$, $d(\text{TSFC})/dr_n$; Figure 2 shows the normalized sensitivity coefficients.

## How to Use
1. **Publish the Script**: From this folder, run `publish('HW6_Q4.m', 'format', 'pdf', 'outputDir', pwd);` in MATLAB. This runs the calculations, writes `HW6_Q4_results.csv`, and creates `HW6_Q4.pdf`.
2. **Access Data**: `HW6_Q4_results.csv` can be opened in any spreadsheet software, and `HW6_Q4.pdf` in any PDF reader.
3. **Analysis**: Use the `.csv` file and `.pdf` figure for reporting or additional analysis (e.g., with an AI tool like Grok 3 for trend identification).

## Notes
- The script assumes ideal components ($\eta_b = 1$, $r_n = 1$) as a baseline, with perturbations of $\Delta = 0.01$.
- At each Mach number, the fuel-to-air ratio is chosen so the baseline ($\eta_b = 1$) reaches $T_{t4} = 2540$ K from problem 2: $f = c_p (T_{t4} - T_{t0}) / (\eta_b h_f)$. This gives $f = 0.053$ at M=1 down to $0.017$ at M=6, always below the stoichiometric limit of 0.06.
- $f$ is then held fixed for the perturbed cases, so a less efficient combustor reaches a lower $T_{t4}$ with the same fuel flow. Because TSFC $= f/I$ with $f$ fixed, each TSFC sensitivity is (to first order) the negative of the matching thrust sensitivity.
- Following problem 2, the $f \ll 1$ approximation is used: thrust is $I = u_9 - u_0$ without the $(1+f)$ factor.
- Trends: thrust sensitivity to $\eta_b$ grows from 0.66 at M=1 to 0.92 at M=6, while sensitivity to $r_n$ drops from 1.04 at M=1 to about 0.10–0.13 above M=4. Combustor efficiency matters most at high Mach numbers; nozzle losses matter most at low Mach numbers.

## Acknowledgments
This work was developed as part of a homework assignment (HW6_Q4) to explore ramjet performance sensitivity. The analysis was enhanced with assistance from Grok 3, an AI tool by xAI, for data interpretation and documentation.