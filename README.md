# Repository Use
Welcome to my code repository for NE 571 Nuclear Reactor Theory at The University of Tennessee Knoxville. This repository exists to help me keep track of previous assignments as well as turn in my code. Any code found in this repository is mine alone (i.e., written by referring to original documentation and without the use of LLMs, community forums, or other people) unless state otherwise.
Feel free to use any code here for your own project.

# Scripts
## Project 1
### Assignment Description
The following equation describes a 2-D, 1-speed, homogenous, finite slab with dimensions $x$ and $y$:
```math
-D\frac{d^2\phi}{d^2x}-d\frac{d^2\phi}{d^2y}+\Sigma_a\phi=\frac{1}{k}\nu\Sigma_f\phi
```
### Deliverables
Write a short report that solves the equation numerically and analytically for any $N_x$ by $N_y$ arrangement of nodes. The flux should be normalized to power for the sake of comparison.
### Code Use Instructions
The code for this project is in `project1.m`. The size and node count of the slab can be altered by varying the `node_x`, `node_y`, `width_x`, and `width_y` variables.
Running the script will result in 2 figures, one for the numerical solution and the other for the analytical solution.

## Project 2
### Assignment Description
*This assignment is a continuation of assignmet one*
### Deliverables
Expand your code from Project 1 to allow for a flux with 2 energy groups. Your code must work independently of the number of nodes, and you must provide the time needed fro your code to find a solution. Calculate the dimentions for a critical reactor with no reflector, then add a reflector with increasing thicknesses to find the maximum reflector savings. Provide plots of the flux of each group and the criticality of the core to the nearest pcm. The flux must me normalized to power of $3000$ $\text{MW}_\text{th}$.
### Code Use Instructions
The code for this project is in `project2.m` which relies on `CreateFissMat.m`, `CreateLossMat.m`, `CreateSctrMat.m`, `SpatialFlux.m`, and `GetFissNormFactor.m` to run. New materials can be added using structs and the the core's layout can be directly assigned or can be quickly generated as a square slab in a larger square slab. The code automatically detects the number of energy groups based on the material array-struct used. The "Quick Square-in-Square" allows the internal core material's width to be assigned instead to the total width. This allows the reflector savings to be found without needed to alter the total width, as this is calculated for you.
