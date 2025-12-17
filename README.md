# JournelDataDrivenControlWithCBF

This repository contains the simulation code for our paper **"Data-Driven Control Based on Control Barrier Functions with Recursive Feasibility Guarantee"**.

## **Related Tools and Software**

### **Simulation Platform**
- MATLAB/Simulink **2024b** (recommended)
  
### **Required MATLAB Toolboxes**
The following MATLAB toolboxes are required to run the simulations:
- Algorithms require:
  **Optimization Toolbox**
- Simulations require:
  **Control System Toolbox**, **Parallel Computing Toolbox**, **Robust Control Toolbox**, **Automated Driving Toolbox**

### **Additional Required Toolboxes**
Please manually install the following toolboxes:
- Programming algorithms require:
  **YALMIP**: [Download Here](https://yalmip.github.io/download/), **MOSEK**: [Download Here](https://www.mosek.com/downloads/)
- Convex Polyhedron require:
  **MPT3**: [Installation Guide](https://www.mpt3.org/pmwiki.php/Main/Installation)

---

## **Setup Instructions**
Follow the steps below to set up the simulation environment:

1. Install **MATLAB/Simulink** (**2024b** recommended) along with the required toolboxes listed above.
2. Install the additional toolboxes (YALMIP, MOSEK, and MPT3) by following the provided links.
3. Download the folders `LK` and `ACC`.
4. Navigate to the `LK` folder and make sure it is the working directory. Then run `LK_Run.mlx` (The plotting code is in a separate mlx file, and `LK_Run` will call it).
5. Navigate to the `ACC` folder and make sure it is the working directory. Then run `ACC_Run.mlx` (The plotting code is in a separate mlx file, and `ACC_Run` will call it).
   
  (Note: Both the `LK` and `ACC` folders contain a `File.mat` file storing the simulation data from the paper. You can directly load these files and use the provided plotting scripts to generate figures, which are located in the `Figures` subfolder of each respective directory.)

---
## **Features**
### **Simulation Results For Lane-Keeping Problem**

1. Lateral displacement and front wheel steering angle (Case 1 with Control Invariant Set, Case 2 without). 
<table align="center">
    <tr>
        <td align="center" style="background-color: white;">
            <img src="https://raw.githubusercontent.com/aicpslab/DDControlWithCBF/main/LK/Figures/LK3.jpg" width="200"><br>
            <b>Case 1: Lateral Displacement</b>
        </td>
        <td align="center" style="background-color: white;">
            <img src="https://raw.githubusercontent.com/aicpslab/DDControlWithCBF/main/LK/Figures/LK5.jpg" width="200"><br>
            <b>Case 1: Steering Angle</b>
        </td>
    </tr>
    <tr>
        <td align="center" style="background-color: white;">
            <img src="https://raw.githubusercontent.com/aicpslab/DDControlWithCBF/main/LK/Figures/LK4.jpg" width="200"><br>
            <b>Case 2: Lateral Displacement</b>
        </td>
        <td align="center" style="background-color: white;">
            <img src="https://raw.githubusercontent.com/aicpslab/DDControlWithCBF/main/LK/Figures/LK6.jpg" width="200"><br>
            <b>Case 2: Steering Angle</b>
        </td>
    </tr>
</table>


2. State variations of other states, \( v, \psi, r \). 

<table>
    <tr>
        <td align="center">
            <img src="https://raw.githubusercontent.com/aicpslab/DDControlWithCBF/main/LK/Figures/LK1.jpg" width="500"><br>
            <b>With Control Invariant Set</b>
        </td>
        <td align="center">
            <img src="https://raw.githubusercontent.com/aicpslab/DDControlWithCBF/main/LK/Figures/LK2.jpg" width="500"><br>
            <b>Without Control Invariant Set</b>
        </td>
    </tr>
</table>

### **Simulation Results For Adaptive Cruise Problem**

1.Iteration Process and Robust Control Invariant Set. 
<p align="center">
    <img src="https://raw.githubusercontent.com/aicpslab/DDControlWithCBF/main/ACC/Figures/ACC1.jpg" width="500"><br>
    <b>Iteration Process and Robust Control Invariant Set.</b>
</p>

2.Velocity, Traction Force, and Safety Function Curves for Autonomous Vehicles. 
<p align="center">
    <img src="https://raw.githubusercontent.com/aicpslab/DDControlWithCBF/main/ACC/Figures/ACC2.jpg" width="800"><br>
    <b>Velocity, Traction Force, and Safety Function Curves for Autonomous Vehicles.</b>
</p>

---

## **Additional Information**
If you encounter any issues, please check the respective toolbox documentation or contact the authors for assistance.
