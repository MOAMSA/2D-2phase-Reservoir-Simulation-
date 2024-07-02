![Solving 2phase 2dimention Reservoir Simulation with Aquifer ](readme.png)

Reservoir studies are conducted to predict future reservoir performance based on past data and to enhance the final hydrocarbon recovery. Reservoir simulators are commonly used for this purpose.

A reservoir simulator is a computer program that solves a system of differential equations representing multiphase flow in a porous medium. Simulators are categorized based on:
* Number of phases and components in the reservoir.
* Type of reservoir processing.
* Flow direction.
* Equations used to solve flow equations.

Based on the number of phases, simulators can be single-phase, two-phase, or three-phase. They can handle from one to multiple components. Regarding reservoir processing, simulators can be classified as Black Oil, Compositional, or Thermal Oil Recovery. Flow direction can be one-dimensional, two-dimensional, or three-dimensional. For solving equations, methods like IMPES, Fully Implicit, or Adaptive Implicit Model can be used. Depending on the reservoir model, the simulator can be Single Porosity, Dual Porosity, or Dual Porosity-Dual Permeability.
The written program is a reservoir simulator (Black Oil model) operating in MATLAB, designed for Single Porosity in two phases (water and oil) and two directions (x and y). It uses the IMPES method, assuming Darcy's flow in two dimensions under isothermal conditions with fluid properties dependent on pressure. The program supports various production mechanisms like water drive, rock and fluid expansion, and gravitational drainage.

The IMPES method employs implicit pressure and explicit water saturation calculations. This method iteratively solves the equations for pressure distribution at each time step, followed by fluid saturation calculation. It's faster and requires less memory compared to Fully Implicit methods but may become unstable with rapid saturation changes, which can be managed by using smaller time steps.

Program capabilities include:

* Handling unlimited wells with different conditions (production, injection, fixed rate, fixed bottom-hole pressure, etc.).
* Modifying various input parameters, including matrices for depth, permeability, dimensions, initial oil pressure, and saturation. It can interpolate values from different tables.
* Automatically performing Material Balance Checks to validate results.
* Allowing input of different time intervals, densities, dimensions, and adjusting boundaries (Aquifer, No Flow Boundary).
* Providing outputs such as total oil produced, recovery factor, final pressure, saturation levels, water cut, and total water injected.

