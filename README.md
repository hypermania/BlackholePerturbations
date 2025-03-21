# BlackholePerturbations
My numerical tools for studying blackhole perturbations. Currently implements two equations:

1. The Regge-Wheeler equation with a customizable source.
2. The Regge-Wheeler equations for all real harmonics modes `l <= l_max`, coupled via a cubic self-interaction.

These equations were used to produce the numerical results in paper [arXiv:XXXX.XXXXX](https://arxiv.org).

## Plotting
Some data are already included in the codebase under `/output`. Run the `plot.nb` Mathematica notebook to plot these data. New data can be generated from the code.

## Running the code
The `src/main.cpp` file has two functions:

1. `run_sourced_eqn()` sets up parameters for the sourced RW equation, and run them in parallel. See section 2.3 of [arXiv:XXXX.XXXXX](https://arxiv.org) for details.
2. `run_coupled_eqn()` sets up parameters for the coupled RW equations corresponding to a scalar field with cubic self-interaction. See section 3.3 of [arXiv:XXXX.XXXXX](https://arxiv.org) for details.

## Compilation
Compiler requirement: a C++ compiler supporting C++20. (I used [g++ 12.2.0](https://gcc.gnu.org/).)

Compilation should be as easy as running `make` at the project directory. Note that the compilation could take a while.

Makes use of `boost` and `Eigen` library, which are included in `/external`.
