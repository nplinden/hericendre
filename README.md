# Welcome to Héricendre

Héricendre is a standalone Bateman equation solver. Its purpose is to be easy to use through human readable `yaml` input
file as well as (in the future) a python API. Its secondary purpose is to teach myself `c++`.

It uses depletion chain in the OpenMC `xml` format.

## Roadmap

Depletion:

- Add HDF5 result files
- Support neutron reactions (no fission)
  - Find a suitable file format to store multigroup cross-sections
- Support secondary particle for neutron reactions
- Support fission

Utilities:

- Support chain reduction

An example:

```toml
name = "An example input file"

[Settings]
chain = "data/chain_casl_sfr.xml"
results = "results/results_cram.csv"
solver = "CRAM48"

[Material]
uniform = 1.0

[Time]
unit = "y"
timestamps = [0, "linspace 1 10 9", 10]
```
