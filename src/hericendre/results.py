import polars as pl
import h5py
import numpy as np


class Results:
    def __init__(self, path):
        self.path = path
        with h5py.File(path) as f:
            self.nuclides = np.array([n.decode("utf-8") for n in f["NUCLIDES"][...]])
            self.concentrations = f["CONCENTRATIONS"][...]
            self.times = f["TIMES"][...]

        self.df = pl.DataFrame(
            np.concatenate([self.times[None, :], self.concentrations.T]),
            schema=[("Time", pl.Float64)] + [(n, pl.Float64) for n in self.nuclides],
        )

    def __str__(self):
        return self.df.__str__()

    def __repr__(self):
        return self.df.__repr__()
