# BuildEML

BuildEML is a Julia package for equation-based building/HVAC system energy modeling and simulation.


## Status

This repository is an early public release of the package and remains under
active development.

At its current stage, BuildEML is a prototype-oriented package. Some of its
current functionality is inspired by the Modelica Buildings Library (MBL)
SimpleHouse examples, which are also used here as verification references,
though the implementation here is not a one-to-one replication.

## Installation

You can install it from GitHub.


## Package Structure

The package currently includes the following main modules:

- `Utils`
  Smoothing functions, unit conversions, time utilities, and interpolation-based helper components.
- `Media`
  Air and water connectors together with supporting thermophysical property functions.
- `FundamentalComponents`
  Basic boundary, sensing, and mixing-volume components for air and water systems.
- `HVACComponents`
  Movers, actuators, pressure-drop elements, heat exchangers, and thermal source models.
- `BuildingComponents`
  Assembled single-zone SimpleHouse building models.
- `Disturbances`
  Weather-file handling, weather buses, and schedule generation/read-in utilities.

## Examples

Example workflows are available in the [`examples/`](examples) directory,
including:

- `examples/simplehouse0to6_verification`
  Reference model: `Buildings.Examples.Tutorial.SimpleHouse` (MBL v11.0.0)
- `examples/simplehouse_standalone_verification`
  Reference model: `Buildings.Examples.SimpleHouse` (MBL v11.0.0)

## Citation

If BuildEML is helpful in your research, please consider citing this
repository. A related paper citation will be added here once the associated
manuscript is published.

## Julia Version

The package currently targets Julia `1.12`.
The package currently depends on ModelingToolkit.jl `v10.21.0`.

## License

This project is released under the MIT License.
