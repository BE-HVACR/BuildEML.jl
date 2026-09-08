# BuildEML

BuildEML is a Julia package for equation-based modeling and simulation of building and HVAC systems. It is built on ModelingToolkit.jl, so models can be used together with the Julia SciML ecosystem, including its differential-equation solvers, automatic differentiation, and optimization tooling.


## Status

This is an early public release under active development. The current release
provides a building energy system modeling prototype that includes single-zone
building models, hydronic and air-side HVAC components, and weather-file
handling.

## Installation

You can install it from GitHub.


## Package Structure

The package currently includes the following main modules:

- `Utils`
  Smoothing functions, unit conversions, time utilities, and interpolation-based helper components.
- `Media`
  Air and water connectors together with supporting thermophysical property functions.
- `FundamentalComponents`
  Basic boundary, and mixing-volume components for air and water systems.
- `HVACComponents`
  Movers, actuators, pressure-drop elements, heat exchangers, and thermal source models.
- `BuildingComponents`
  Assembled single-zone SimpleHouse building models.
- `Disturbances`
  Weather-file handling, TMYx files, and schedule (occupant/interal heat gain) generation.

## Verification examples

Example workflows are available in the [`examples/`](examples) directory,
including:

- `examples/simplehouse0to6_verification`
  Reference model: `Buildings.Examples.Tutorial.SimpleHouse` (MBL v11.0.0)
- `examples/simplehouse_standalone_verification`
  Reference model: `Buildings.Examples.SimpleHouse` (MBL v11.0.0)

## Citation

If you use BuildEML in your research, please cite the following accepted
conference paper, which describes the modeling prototype and its verification:

```
R. Song, M. Liu, Z. Yang, and Z. O'Neill. "Toward Equation-Based Building
Energy System Modeling in Julia: A Modular Prototype and Verification Study."
Proceedings of the American Modelica & FMI Conference 2026, Atlanta, GA, USA,
October 12–14, 2026. To appear.
```

## Julia Version

The package currently targets Julia `1.12`.
The package currently depends on ModelingToolkit.jl `v10.21.0`.

## License

This project is released under the MIT License.
