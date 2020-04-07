# Connect to database and get 3000 sample MOFs.
# Let's assume that these contain an id, a charged RASPA molfile, and a helium void fraction.
import pymongo
db = pymongo.MongoClient().sample
mofs = db.mof.find().limit(3000)

print "start"
for doc in mofs:
    print "mofs", doc
print "end"

# On each core, we will get the v/v methane uptake at 65 bar.
def f(mof):
    import RASPA2
    output = RASPA2.run(mof["charged mol"], "CH4", temperature=298, pressure=65e5, helium_void_fraction=mof["helium void fraction"], input_file_type="mol")
    print 'mfs', mof, output
    return output["Number of molecules"]["CH4"]["Average loading absolute [cm^3 (STP)/cm^3 framework]"][0]

# Distribute the jobs on the amazon cloud cores. Wait for the jobs to finish.
#import cloud
#jids = cloud.map(f, mofs)
#uptakes = cloud.result(mofs)
#map(f, mofs)
#uptakes = result(mofs)

# Let's print the results for postprocessing, and additionally save them back into the database.
# By saving into a central database, you can share your simulation results with everyone else. (Useful for data mining!)
#for mof, uptake in zip(mofs, uptakes):
#    print mof["_id"], uptake
#    db.mof.update({"_id": mof["_id"]}, {"$set": {"CH4": uptake}})
