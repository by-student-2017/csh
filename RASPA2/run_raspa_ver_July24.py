# Connect to database and get 3000 sample MOFs.
# Let's assume that these contain an id, a charged RASPA molfile, and a helium void fraction.
import copy
import pymongo
db = pymongo.MongoClient().sample
mofs = db.mof.find().limit(200)

# On each core, we will get the v/v methane uptake at 65 bar.
def f(mof):
    import RASPA2
    output = RASPA2.run(mof["cif"], "H2", temperature=298, pressure=65e5, helium_void_fraction=mof["helium void fraction"], input_file_type="cif")
    return output["Number of molecules"]["H2"]["Average loading absolute [milligram/gram framework]"][0]

# Distribute the jobs on the amazon cloud cores. Wait for the jobs to finish.
#import cloud
#jids = cloud.map(f, mofs)
#uptakes = cloud.result(mofs)
mofs_data = copy.deepcopy(mofs)
uptakes = map(f, mofs)

# Let's print the results for postprocessing, and additionally save them back into the database.
# By saving into a central database, you can share your simulation results with everyone else. (Useful for data mining!
for mof, uptake in zip(mofs_data, uptakes):
    print mof["_id"], uptake
    db.mof.update({"_id": mof["_id"]}, {"$set": {"H2": uptake}})
