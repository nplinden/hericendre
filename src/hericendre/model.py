from pathlib import Path
import pint
from .units import ur
import subprocess
import shutil
import os


class Model:
    def __init__(self, name=""):
        self.name = name
        self.settings = Settings()
        self.time = Time()
        self.material = Material()

        if self.name:
            self.settings.results = f"{name}.h5"
        else:
            self.settings.results = "hericendre_results.h5"

    def __str__(self):
        concentrations = ".".join(
            [f"{k} = {v}" for k, v in self.material.concentrations.items()]
        )
        s = f"""name = "{self.name}"

[Settings]
chain = "{self.settings.chain.absolute()}"
results = "{self.settings.results.absolute()}"
solver = "{self.settings.solver}"

[Time]
unit_name = {{"name" = "{self.time.unit.u}", "magnitude" = {self.time.unit.to("s").m}}}
timestamps = {[float(t) for t in self.time.timestamps]}

[Material]
concentrations = {{{concentrations}}}
"""
        return s

    def __repr__(self):
        return self.__str__()

    def export_to_toml(self, path="model.toml"):
        with open(path, "w") as f:
            print(self, file=f)

    def run(self):
        self.export_to_toml()
        if "HERICENDRE" in os.environ:
            subprocess.run([os.environ["HERICENDRE"], "model.toml"])
        else:
            print("No HERICENDRE environment variable is set.")


class Settings:
    def __init__(self):
        self.solver = "CRAM48"

    @property
    def solver(self):
        return self._solver

    @solver.setter
    def solver(self, rhs):
        self._solver = rhs

    @property
    def results(self):
        return self._results

    @results.setter
    def results(self, rhs):
        self._results = Path(rhs)

    @property
    def chain(self):
        return self._chain

    @chain.setter
    def chain(self, rhs):
        self._chain = Path(rhs)


class Time:
    def __init__(self):
        self.unit = ur("s")

    @property
    def unit(self):
        return self._unit

    @unit.setter
    def unit(self, rhs):
        if isinstance(rhs, pint.Quantity):
            self._unit = rhs
        elif isinstance(rhs, str):
            self._unit = ur(rhs)
        else:
            raise ValueError("Can't parse unit")

    @property
    def unit_magnitude(self):
        return self.unit.to("s").m

    @property
    def unit_name(self):
        return str(self.unit.u)

    @property
    def timestamps(self):
        return self._timestamps

    @timestamps.setter
    def timestamps(self, rhs):
        self._timestamps = rhs


class Material:
    def __init__(self):
        pass

    @property
    def concentrations(self):
        return self._concentrations

    @concentrations.setter
    def concentrations(self, rhs):
        self._concentrations = rhs
