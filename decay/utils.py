"""Function or variable used for setup."""
import importlib

ROOT_PATH = importlib.util.find_spec("decay").submodule_search_locations[0]