# Héricendre Input File Reference

This document describes the TOML input file format for Héricendre simulations.

## File Structure

The input file uses [TOML](https://toml.io/) format and is organized into four main sections:
- Top-level settings
- `[Settings]` - Simulation configuration
- `[Time]` - Time grid definition
- `[Material]` - Initial nuclide concentrations

---

## Top-Level Settings

### `name`
**Type:** String (optional)  
**Description:** Name of the simulation run.

**Example:**
```toml
name = "U235 Decay Simulation"
```

---

## `[Settings]` Section

### `secondaries`
**Type:** Boolean (optional)  
**Description:** Include secondary nuclides in the chain. Defaults to `true` if not specified.

**Example:**
```toml
[Settings]
secondaries = true
```

### `chain`
**Type:** String (required)  
**Description:** Path to the depletion chain XML file containing nuclide and decay data.

**Example:**
```toml
[Settings]
chain = "/path/to/chain.xml"
```

### `results`
**Type:** String (required)  
**Description:** Output file path for simulation results. File extension determines format:
- `.csv` - CSV format
- `.h5` or `.hdf5` - HDF5 format

**Example:**
```toml
[Settings]
results = "output.h5"
```

### `save_matrix`
**Type:** String (optional)  
**Description:** Path to save the decay matrix in CSV format.

**Example:**
```toml
[Settings]
save_matrix = "decay_matrix.csv"
```

### `solver`
**Type:** String (optional)  
**Description:** Solver algorithm to use. If not specified, defaults to `CRAM48`.

**Allowed values:**
- `"CRAM48"` - Chebyshev Rational Approximation Method (48th order)
- `"Decay"` - Analytical decay solver

**Example:**
```toml
[Settings]
solver = "CRAM48"
```

---

## `[Time]` Section

### `timestamps`
**Type:** Array (required)  
**Description:** Time points for the simulation. Can contain:
- **Numbers** (integers or floats) - Explicit time points
- **Strings** - Time functions that generate multiple points

**Time Functions:**

#### `linspace`
Generates linearly spaced points.

**Format:** `"linspace START STOP NPOINTS"`

**Example:**
```toml
[Time]
timestamps = ["linspace 0 100 11"]  # Generates: 0, 10, 20, ..., 100
```

#### `logspace`
Generates logarithmically spaced points (base 10).

**Format:** `"logspace LOG_START LOG_STOP NPOINTS"`

**Example:**
```toml
[Time]
timestamps = ["logspace 0 6 7"]  # Generates: 10^0, 10^1, ..., 10^6
```

**Mixed Example:**
```toml
[Time]
timestamps = [0, 1e3, 1e4, "linspace 1e5 1e6 10", "logspace 6 9 4"]
```

### `unit`
**Type:** String or Table (optional)  
**Description:** Time unit for the timestamps. If not specified, defaults to seconds.

#### Option 1: Predefined Units (String)

**Allowed values:**
- `"s"` or `"second"` - seconds (default)
- `"h"` or `"hour"` - hours (3600 seconds)
- `"d"` or `"day"` - days (86400 seconds)
- `"y"` or `"year"` - years (365.25 days)

**Example:**
```toml
[Time]
timestamps = [0, 1, 2, 5, 10]
unit = "year"
```

#### Option 2: Custom Unit (Table)

Define a custom unit with a name and magnitude (conversion factor to seconds).

**Example:**
```toml
[Time]
timestamps = [0, 100, 200]
unit = { name = "minute", magnitude = 60.0 }
```

**Complete Time Example:**
```toml
[Time]
timestamps = ["linspace 0 10 11", "logspace 1 5 5"]
unit = "year"
```

---

## `[Material]` Section

Defines initial nuclide concentrations and optional microscopic cross sections.

### `microxs`
**Type:** String (optional)  
**Description:** Path to microscopic cross section data file. Required when using `CRAM48` solver with neutron-induced reactions. Cannot be used with `Decay` solver.

**Example:**
```toml
[Material]
microxs = "microxs_data.h5"
```

### Initial Concentrations

Use **either** `uniform` **or** `concentrations`, not both.

#### Option 1: `uniform`
**Type:** Number  
**Description:** Sets the same concentration for all nuclides in the chain.

**Example:**
```toml
[Material]
uniform = 1.0
```

#### Option 2: `concentrations`
**Type:** Table  
**Description:** Specifies individual concentrations for each nuclide by name.

**Example:**
```toml
[Material]
concentrations = { U235 = 1.0, U238 = 0.5, Pu239 = 0.1 }
```

**Multi-line format:**
```toml
[Material.concentrations]
U235 = 1.0
U238 = 0.5
Pu239 = 0.1
```

---

## Complete Example

```toml
name = "Uranium Decay Analysis"

[Settings]
chain = "chains/endfb8.xml"
results = "results/uranium_decay.h5"
solver = "CRAM48"
secondaries = true

[Time]
timestamps = [0, "linspace 1e3 1e6 100", "logspace 6 9 10"]
unit = "year"

[Material]
microxs = "data/microxs.h5"

[Material.concentrations]
U235 = 1.0
U238 = 0.007
```

---

## Validation Rules

1. **Timestamps must be sorted** - Time values must be in ascending order
2. **Unique section keys** - Cannot define both `uniform` and `concentrations` in `[Material]`
3. **Required fields** - `chain`, `results`, and `timestamps` are mandatory
4. **File existence** - The chain file must exist and be a valid depletion chain XML
5. **Nuclide names** - Concentrations must reference nuclides that exist in the chain
6. **Solver compatibility** - Cannot use `Decay` solver with `microxs` data

---

## Error Messages

Common errors and their meanings:

- `"No depletion chain file was specified"` - Missing `chain` in `[Settings]`
- `"Depletion chain file {} does not exist"` - Invalid path to chain file
- `"No result file was specified"` - Missing `results` in `[Settings]`
- `"No timestamps provided"` - Missing `timestamps` in `[Time]`
- `"Time vector is not sorted"` - Timestamps are not in ascending order
- `"Invalid solver type"` - Solver must be "CRAM48" or "Decay"
- `"uniform and concentration can't be defined at the same time"` - Use only one in `[Material]`
- `"Can't use \"Decay\" solver with microxs data"` - Remove microxs when using Decay solver
- `"Invalid time function in timestamps definition"` - Time function format is incorrect (needs at least 4 parts)
- `"Invalid unit name"` - Time unit not recognized
- `"Invalid type in timestamps definition"` - Timestamp must be number or string
- `"Concentration must be a number"` - Concentration values must be numeric
- `"Invalid time unit name"` - Custom unit must have a name field
- `"Invalid time unit magnitude"` - Custom unit must have a numeric magnitude field
- `"Time unit must be either a string or a table"` - Unit format is invalid
