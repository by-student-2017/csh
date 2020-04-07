import pybel
import RASPA2
import pandas as pd
import matplotlib.pyplot as plt

# Set up
gas = "H2"
pressures = [1e4 * 10**(0.1 * i) for i in range(21)]

#pressures = [350e6, 700e6]
# Use pybel to parse, fill, and charge cif structure
my_structure = pybel.readfile("cif", "IRMOF-1.cif").next()
#my_structure.unitcell.FillUnitCell(mol.OBMol)
#my_structure.calccharges("eqeq")

# Run
results = [RASPA2.run(my_structure, gas, temperature=298, pressure=pressure) for pressure in pressures]

# Parse
uptakes = [r["Number of molecules"][gas]["Average loading absolute [cm^3 (STP)/cm^3 framework]"][0] for r in results]

# print
z = [pressures,uptakes]
df = pd.DataFrame(z).T
df.to_csv("Pa_vs_NumAbsGas.csv")

# Plot
plt.plot(pressures, uptakes)
plt.show()
