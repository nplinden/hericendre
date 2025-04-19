# Welcome to Héricendre

A work in progress depletion code, for now supporting radioactive decay with CRAM and analytical solutions.

## Installation

You can install Héricendre with `pip`:

```console
pip install hericendre
```

This would install the Héricendre binary and library in your python environment directory.

Or add it to your `uv` project:

```console
uv add hericendre
```

This makes the Héricendre library available for your project.
To run the binary simply call

```console
uv run hericendre
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
