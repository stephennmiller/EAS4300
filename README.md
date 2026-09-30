# EAS 4300 – Propulsion (UCF)

MATLAB coursework for EAS 4300: gas dynamics of a Mach 5 test facility (Project 1) and cycle analysis of ramjets (HW6), turbojets with and without afterburning (HW7), and turbofans (HW8).

Each homework question has its own folder with a MATLAB script, a CSV of results, a published PDF report (code, equations, output, and figures), and a README with the problem statement and discussion.

## Contents

| Folder | Topic | Report |
|---|---|---|
| [Project 1](Project%201) | Mach 5 nozzle, Fanno channel, normal shock, tank blowdown, H₂ combustion | [Project_1.m](Project%201/Project_1.m) |
| [Homework 6 / HW6_Q3](Homework%206/HW6_Q3) | Ramjet performance vs. Mach number (M = 1–6) | [HW6_Q3.pdf](Homework%206/HW6_Q3/HW6_Q3.pdf) |
| [Homework 6 / HW6_Q4](Homework%206/HW6_Q4) | Ramjet sensitivity to $\eta_b$ and $r_n$ | [HW6_Q4.pdf](Homework%206/HW6_Q4/HW6_Q4.pdf) |
| [Homework 7 / HW7_Q1](Homework%207/HW7_Q1) | Turbojet performance vs. compressor pressure ratio | [HW7_Q1.pdf](Homework%207/HW7_Q1/HW7_Q1.pdf) |
| [Homework 7 / HW7_Q2](Homework%207/HW7_Q2) | Turbojet with afterburner | [HW7_Q2.pdf](Homework%207/HW7_Q2/HW7_Q2.pdf) |
| [Homework 8 / HW8_Q1](Homework%208/HW8_Q1) | Turbofan: optimum bypass ratio | [HW8_Q1.pdf](Homework%208/HW8_Q1/HW8_Q1.pdf) |
| [Homework 8 / HW8_Q2](Homework%208/HW8_Q2) | Turbofan: optimum $r_c$ and $r_f$ | [HW8_Q2.pdf](Homework%208/HW8_Q2/HW8_Q2.pdf) |

## Key Results

### Project 1 – Mach 5 Test Facility

- **Q1:** Nozzle area ratio $A_e/A^* = 25$ for $M = 5$, giving a throat area of $8.1\times10^{-5}$ m².
- **Q2:** Fanno friction in the 0.225 m channel ($fL/D = 0.1$) slows the flow from M = 5 to M = 3.57. Holding a normal shock at the exit then needs $P_{01} \approx 2.00$ MPa, with $P_{03}/P_{02} = 0.20$ across the shock.
- **Q3:** Mass flow is 0.383 kg/s, so a 2-minute run expels 46 kg of air. **4 tanks** are needed: the tank throats must stay choked and still pass the nozzle's mass flow, which ends the run at about 5.1 MPa tank pressure. A throat-flow cross-check also gives 4 tanks.
- **Q4:** Stoichiometric H₂–air combustion releases 241,820 kJ/kmol H₂ (120 MJ/kg). That is 1.35 MW at the Q2 mass flow.
- **Q5:** Ideally expanding the flow at the channel exit would need $P_0 \approx 29.4$ MPa, **14.7×** the 2.00 MPa needed to hold the shock.

### Homework 6 – Ramjet

- **Q3 (M = 1–6 at 10 km, $T_{04,max}$ = 2600 K):**
  - Specific thrust peaks at **1103 m/s at M ≈ 2.93**, where the engine first reaches 2600 K.
  - TSFC is lowest at **5.28e-5 kg/N·s near M ≈ 3.94**.
  - The engine runs at $f_{st}$ = 0.06 below M ≈ 2.93.
  - $\eta_0$ rises to 0.737 at M = 6, and $\eta_{th}$ to 0.782, below the ideal Brayton limit of 0.878.
- **Q4 (sensitivity):**
  - Thrust sensitivity to combustor efficiency, $(1/I)\,dI/d\eta_b$, grows from 0.66 at M = 1 to 0.92 at M = 6.
  - Sensitivity to nozzle pressure ratio, $(1/I)\,dI/dr_n$, falls from 1.04 to about 0.1 above M = 4.
  - Nozzle losses matter most at low Mach numbers, and combustor efficiency at high ones.

### Homework 7 – Turbojet (15 km, M = 1.8, $T_{04}$ = 1500 K)

| | Without afterburner | With afterburner ($T_{06}$ = 2000 K) |
|---|---|---|
| Max specific thrust | **708 m/s at $r_c \approx 6.1$** | **1106 m/s at $r_c \approx 21$** |
| Min TSFC | **3.59e-5 s/m at $r_c \approx 43$** | **4.93e-5 s/m at $r_c \approx 21$** |
| Max $\eta_0$ | 0.34 | 0.25 |
| Total fuel-to-air ratio | 0.015–0.035 | 0.054–0.055 (below $f_{st}$ = 0.06) |

The afterburner raises thrust by 1.4–2.4×, but it also raises TSFC by 14–46%. The reason is that it burns its fuel at low pressure, and the faster exhaust lowers propulsive efficiency.

### Homework 8 – Turbofan (10 km, M = 0.8, $T_{04}$ = 1500 K)

- **Q1:** Overall efficiency peaks at **$\eta_0$ = 0.281 at bypass ratio β ≈ 5.88**. At the same point:
  - specific thrust is at its maximum, 1218 m/s;
  - TSFC is at its minimum, 1.98e-5 s/m;
  - $f$ = 0.0241.

  Above β ≈ 7.3, the turbine can no longer drive the fan and still let the core expand to ambient pressure.
- **Q2:** At β = 5.88, a different combination does better: **$\eta_0$ = 0.287 at $r_c$ = 28, $r_f \approx 1.89$**, vs. 0.281 at $r_c$ = 24, $r_f$ = 2. $\eta_0$ keeps rising with $r_c$ up to the edge of the range studied. The best fan pressure ratio stays between 1.89 and 1.93 for every $r_c$.

## Modeling Notes

- **Burner energy balance:** all cycle analyses use $c_{p,air}T_{in} + f\,\eta_b\,\Delta h_c = (1+f)\,c_{p,gas}T_{out}$. $c_{p,air}$ (γ = 1.4) applies to the incoming air and $c_{p,gas}$ to the combustion products.
  - Using one $c_p$ for both overstates the fuel energy, and in HW6 and HW7 it gave thermal efficiencies above the ideal Brayton limit.
- **Stoichiometric limit:** if reaching the target temperature would need more fuel than $f_{st}$, the script burns $f_{st}$ and solves for the lower temperature it reaches.
  - For the afterburner, the limit applies to the total, $f_b + f_{ab} \le f_{st}$.
- Optima are found with `max`/`min` on the computed curves. No toolboxes are required.

## Running the Scripts

Open a question folder in MATLAB and publish its script, for example:

```matlab
publish('HW7_Q1.m', 'format', 'pdf', 'outputDir', pwd);
```

This runs the analysis, writes the CSV, and regenerates the PDF report. Project 1 prints its results to the Command Window, so run it with `run('Project_1.m')`.

## Author

Stephen Miller, University of Central Florida. All analysis is my own work; the problem statements were provided by the course instructor.
