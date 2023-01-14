### Get coordinates of the binding site center in Chain A (mean of the PDB_4xv2 ligand Dabrafenib)

from Bio.PDB.PDBParser import PDBParser
import statistics

parser = PDBParser(PERMISSIVE=1)
structure = parser.get_structure("Dabrafenib", "./files/Ligand_4xv2_ChA.pdb")
model = structure[0]
ligand = model.get_list()[0].get_list()[0]

def get_coords(PDB_ligand):
    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure("Dabrafenib", PDB_ligand)
    model = structure[0]
    ligand = model.get_list()[0].get_list()[0]
    
    coord_x=[]
    coord_y=[]
    coord_z=[]

    for atom in ligand:
        coord = (atom.get_coord())
        coord_x.append(coord[0])
        coord_y.append(coord[1])
        coord_z.append(coord[2])
                            
    mean_x = statistics.mean(coord_x)
    mean_y = statistics.mean(coord_y)
    mean_z = statistics.mean(coord_z)

    return(mean_x, mean_y, mean_z)

print("Binding center Chain A:", get_coords("./files/Ligand_4xv2_ChA.pdb"))
print("Binding center Chain B:", get_coords("./files/Ligand_4xv2_ChB.pdb"))