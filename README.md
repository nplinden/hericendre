# Welcome to Héricendre
Héricendre is a standalone Bateman equation solver. Its purpose is to be easy to use through easily readable `yaml` input file as well as (in the future) as python API. Its secondary purpose is to teach myself `c++`.

It uses depletion chain in the OpenMC `xml` format.

## Roadmap

Depletion:
- Support secondary particle for decay reactions
- Add HDF5 result files
- Support neutron reactions (no fission)
    - Find a suitable file format to store multigroup cross-sections
- Support secondary particle for neutron reactions
- Support fission

Utilities:
- Support chain reduction

