import pybel

# Use pybel to parse, fill, and charge cif structure
mol = pybel.readfile("cif", "ztc_CCDC1580143.cif").next()
#mol.unitcell.FillUnitCell(mol.OBMol)

output = pybel.Outputfile("cif", "streamed.cif")
output.write(mol)
output.close()