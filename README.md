# Welcome to Héricendre

## Installation

Héricendre uses [vcpkg](https://vcpkg.io/en/) for dependency management.
Simply install vcpkg and identify the path to the `vcpkg.cmake` file.
It should be located in `$VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake`.

To only build the `hericendre` binary, go to the Héricendre root directory and do:

```console
$ cmake -Bbuild/ -DCMAKE_TOOLCHAIN_FILE=/path/to/vcpkg.cmake .
$ cd build
$ make
```

To build the `hericendre` python package:

```console
$ uv build -Ccmake.define.CMAKE_TOOLCHAIN_FILE=/path/to/vcpkg.cmake
```

## Use

The compute the decay of a single nuclide over 10 years:

```toml
name = "Example"

[Settings]
chain = "data/chain_casl_sfr.xml"
results = "results/results_cram.h5"
solver = "CRAM48"

[Material]
concentrations = { Pu239 = 1 }

[Time]
unit = "y"
timestamps = [0, "linspace 1 10 9", 10]
```

To compute the decay of a material containing 1 single unit of every nuclide in the chain:

```toml
name = "Example"

[Settings]
chain = "data/chain_casl_sfr.xml"
results = "results/results_cram.h5"
solver = "CRAM48"

[Material]
uniform = 1.0

[Time]
unit = "y"
timestamps = [0, "linspace 1 10 9", 10]
```
