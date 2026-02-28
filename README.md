# Multiphysics Optimization and FEA Validation of AFPM Motors

This repository contains a Python-based design and optimization engine for Axial Flux Permanent Magnet (AFPM) synchronous machines. The project focuses on bridging the gap between analytical pre-design and numerical validation (FEA) by integrating seven distinct physical domains.

## FEA Validation Heatmap
Below is the generated magnetic flux density distribution of the optimized design, automatically validated through the automated Python-FEMM bridge. It clearly illustrates the fringing effects, non-linear saturation levels, and full magnetic circuit completion in the M-19 steel stator teeth and rotor yoke.

![Optimized AFPM Motor Flux Density Heatmap](femm_screen.png)

## Technical Overview
The core algorithm utilizes a Sequential Least Squares Programming (SLSQP) optimizer to navigate a non-linear design space. Unlike simplified magnetic circuit models, this engine accounts for critical electromagnetic and mechanical second-order effects that are often neglected in early-stage prototyping.

### Implemented Analytical Modules
The optimization is governed by the following physics-driven modules:
* **Fringing & Leakage Control:** Implementation of Carter’s coefficient ($k_c$) for precise air-gap permeance modeling.
* **High-Frequency Loss Evaluation:** Utilization of Dowell's equations to predict skin and proximity effects in stator windings.
* **Non-linear Magnetic Modeling:** 9th-order power-law fit for silicon steel B-H characteristics to monitor $K_{sat}$.
* **Aerodynamic Analysis:** Reynolds-based windage loss estimation for high-RPM disc rotors.
* **Thermal Magnet Safety:** Calculation of eddy current losses in NdFeB magnets based on slot-passing harmonics.
* **Structural Stability:** Maxwell Stress Tensor analysis coupled with plate theory to calculate axial pull and rotor deflection.
* **NVH Optimization:** Minimization of peak-to-peak cogging torque through LCM (Least Common Multiple) pole-slot selection.

## Project Structure
* `/src`: Core Python scripts and physics modules.
* `/docs`: IEEE-format technical report and analytical derivations.
* `/results`: FEA heatmaps and flux density plots generated via the FEMM bridge.

## Performance Metrics
The finalized designs achieve an efficiency of >92% while maintaining a saturation factor below 1.5 and a thermal ceiling of 120°C. Structural integrity is ensured by keeping rotor deflection within 10% of the nominal air-gap.

## Author
**Mehmet Ege Göz** Electronics and Communication Engineering  
Contact: [mehmetegegoz@gmail.com](mailto:mehmetegegoz@gmail.com) | [GitHub](https://github.com/mege97)
