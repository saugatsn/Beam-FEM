# Beam Analysis using Finite Element Method

A Python script for analyzing beams using the Finite Element Method (FEM). This tool allows users to:

- Define beam properties (length, number of elements)
- Specify support conditions (fixed, roller, hinge)
- Apply various load types (point loads, uniformly distributed loads, moments)
- Calculate and visualize:
  - Displacement
  - Rotation
  - Shear force
  - Bending moment

Features include ability to solve indeterminate beams, interactive input, support for multiple load cases, and matplotlib visualizations of results.

## Recommended: 100 number of elements or greater

## Key Components:
- Stiffness matrix generation
- Global assembly
- Boundary condition application
- Load vector calculation
- Equation solving using SymPy
- Result visualization with Matplotlib
