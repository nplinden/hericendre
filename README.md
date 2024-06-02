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

```yaml
Solver: Decay
Chain: /home/nlinden/workspace/hericendre/data/chain-jeff33-long.xml
Results: results_decay.csv
TimeMode: Timestamps
Time: >
  logspace 1e-10 1e17 271 ;
  1e17
ConcentrationMode: Uniform
Concentrations: 1.
```

Another:

```yaml
Solver: CRAM48
Chain: /home/nlinden/workspace/hericendre/data/chain-jeff33-long.xml
Results: results_decay.csv
TimeMode: Timestamps
Time: >
  logspace 1e-10 1e17 271 ;
  1e17
ConcentrationMode: Explicit
Concentrations: >
  Pu239 10.
  U235 10.
```
