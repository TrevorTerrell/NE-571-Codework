# Repository Use
Welcome to my code repository for NE 571. This repository exists to help me keep track of previous assignments as well as turn in my code.
Any code found in this repository is mine alone (i.e., written by referring to original documentation and without the use of LLMs, community forums, or other people) unless state otherwise.
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

