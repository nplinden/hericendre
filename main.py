from hericendre import Model
import numpy as np

model = Model(name="Test")
model.settings.chain = "playground/data/chain_casl_sfr.xml"

model.time.unit = "a"
model.time.timestamps = np.linspace(0, 10, 50)

model.material.concentrations = {"Pu239": 1}

print(model)
model.run()
