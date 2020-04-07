import glob
import pybel
import RASPA2
from pymongo import MongoClient
from datetime import datetime

# get cif
cif_list = glob.glob('*.cif')

for cif_file in cif_list:
    print cif_file   # python2
    # print(cif_file)    # python3
    #
    # Use pybel to parse, fill, and charge cif structure
    mol = pybel.readfile("cif", cif_file).next()
    mol.unitcell.FillUnitCell(mol.OBMol)
    print mol;
    #mol.calccharges("eqeq")
    #
    # Mongo setting
    #client_1 = MongoClient()
    #client_2 = MongoClient('localhost', 27017)
    #client_3 = MongoClient('mongodb://localhost:27017/') 
    #
    class Mongo_cif_DB(object):
        
        def __init__(self):
            self.clint = MongoClient()
            self.db = self.clint['mof']
        
        def add_one(self):
            """insert data"""
            post = {
                    'title': cif_file,
                    'mol': mol,
                    'created_at': datetime.now()
                    }
            return self.db.mof.insert_one(post)
    
    def main():
        obj = Mongo_cif_DB()
        rest = obj.add_one()
        print(rest)
    
    if __name__ == '__main__':
        main()

