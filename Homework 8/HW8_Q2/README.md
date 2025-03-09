# EAS 4300 HW8 Q2 - 2D Parametric Study of Turbofan Engine Performance

This repository contains MATLAB code that performs a two-dimensional parametric study of a turbofan engine’s performance. The code varies the fan pressure ratio $r_f\$ and the compressor pressure ratio $r_c$ to compute the overall efficiency $eta_0$ of the engine. A surface plot is generated to visualize the relationship between these parameters and the efficiency.

## Files

- **HW8_Q2.m**  
  The main MATLAB script that:
  - Defines engine and flight parameters.
  - Computes performance metrics (including overall efficiency) over a grid of $r_f$ and $r_c$ values.
  - Generates a surface plot of overall efficiency.
  
- **HW8_Q2_data.csv**  
Data file from script

- **README.md**  
  This file.

## Requirements

- MATLAB (preferably R2020a or later for full functionality)
- Basic familiarity with MATLAB and thermodynamic cycle analysis

## How to Run

1. **Clone or Download the Repository**  
   Ensure all files are in your working directory.

2. **Run the MATLAB Script**  
   Open MATLAB and navigate to the repository directory. Execute the script by typing:
   ```matlab
   HW8_Q2

The script will:

- Perform a full parametric sweep over $r_f $ and $r_c$.
- Compute the overall efficiency $eta_0$.
- Save the computed data into `HW8_Q2_data.csv`.
- Generate and export a 3D surface plot as `HW8_Q2_surface_plot.pdf`.

**Review the Output:**

Open `HW8_Q2_data.csv` to inspect the numerical results and view `HW8_Q2_surface_plot.pdf` for the graphical representation of the results.

## Graph Explanation

The generated graph is a 3D surface plot where:
- **X-axis:** Represents the fan pressure ratio $r_f$, varied between 1.5 and 2.2.  
- **Y-axis:** Represents the compressor pressure ratio $r_c$, varied between 20 and 28.
- **Z-axis:** Represents the overall efficiency $eta_0$ of the turbofan engine.

The surface plot allows you to visualize how changes in both the fan and compressor pressure ratios affect the engine’s overall efficiency. A color bar is provided to indicate the range of efficiency values. This visualization helps in identifying the combinations of $r_f$ and $r_c$ that yield optimal performance.

## Key Equations

The analysis is based on several key thermodynamic equations. For example:

**Inlet Total Temperature:**
$$T_{0a} = T_a \left(1 + \frac{k-1}{2} M_f^2\right)$$

**Fuel-to-Air Ratio (Energy Balance):**
$$F_b = \frac{T_{04_{\text{max}}} - T_3}{\left(\frac{\eta_b \, dh_c}{C_{p2}}\right) - T_{04_{\text{max}}}}$$

**Overall Efficiency:**
$$\eta_0 = \eta_{th} \times \eta_p$$

## Customization

You can adjust several parameters in the script:
- **Engine and Flight Conditions:**  
  - Modify ambient conditions:  
    - Temperature: $T_a$  
    - Pressure: $P_a$  
    - Mach number: $M_f$  
    - Turbine inlet temperature: $T_{04_{\text{max}}}$  

- **Pressure Ratio Ranges:**  
  - Change the range and resolution of:  
    - Fan pressure ratio: $r_f$  
    - Compressor pressure ratio: $r_c$  
  - Adjust the vectors `rfVec` and `rcVec`.

- **Efficiency Parameters:**  
  Update the efficiency values (e.g., diffuser, compressor, turbine) as needed.

## Acknowledgments

- Original code was created by Stephen Miller.
- This project was developed for educational purposes in the EAS 4300 course.
