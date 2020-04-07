import RASPA2
input_script = RASPA2.create_script("H2",temperature=298, pressure=1000000, helium_void_fraction=1.0, input_file_type="cif")
print input_script