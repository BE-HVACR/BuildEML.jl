# BuildEML

BuildEML is a Julia package for equation-based building/HVAC system energy modeling and simulation.


## Status

This repository is an early public release of the package and is currently under development.

At its current stage, BuildEML is a prototype-oriented package. Some of its
current functionality is inspired by the Modelica Buildings Library (MBL)
SimpleHouse examples, which are also used here as verification references,
though the implementation here is not a one-to-one replication.

## Installation

You can install it from GitHub.


## Package Structure

The package currently includes the following main modules:

- `Utils`
- `Media`
- `FundamentalComponents`
- `HVACComponents`
- `Disturbances`
- `BuildingComponents`

## Examples

Example workflows are available in the [`examples/`](examples) directory,
including:

- `examples/simplehouse0to6_verification`
- `examples/simplehouse_standalone_verification`

## Citation

If BuildEML is helpful in your research, please consider citing this
repository. A related paper citation will be added here once the associated
manuscript is published.

## Julia Version

The package currently targets Julia `1.12`.

## License

This project is released under the MIT License.
