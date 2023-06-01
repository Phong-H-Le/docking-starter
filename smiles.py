import os
import chimera
from chimera import runCommand as rc
from chimera import replyobj
from chimera import openModels, Molecule
from functions import *

os.chdir("/Users/phong/Desktop/Lynch Lab/PBP2a Docking")

# ligands_dir = "./ligands"
# file_names = [fn for fn in os.listdir(ligands_dir) if fn.endswith(".sdf")]
receptors_dir = "./"
receptor_names = [fn for fn in os.listdir(receptors_dir) if fn.endswith(".pdb")]
vina_location = "'C:/Program Files (x86)/The Scripps Research Institute/Vina/vina.exe'"
output_location = "'C:/Users/phong/Desktop/Lynch Lab/PBP2a Docking/output2/test"
xyz_cords = "30,28,86"
xyz_size = "20,28,20"
steps = 2000

def get_atom(model_number, sub_model_number, residue_number, atom_names):
    model_number = int(model_number)
    sub_model_number = int(sub_model_number)
    models = chimera.openModels.list(modelTypes=[Molecule])
    model = next((m for m in models if m.id == model_number and m.subid == sub_model_number), None)
    if model is None:
        raise ValueError("Model with number {model_number}.{sub_model_number} not found")

    residue = None
    for res in model.residues:
        if res.id.position == residue_number:
            residue = res
            break

    if residue is None:
        raise ValueError("Residue {residue_number} not found in model {model_number}.{sub_model_number}")

    for atom_name in atom_names:
        atom = residue.findAtom(atom_name)
        if atom is not None:
            break

    if atom is None:
        raise ValueError(
            "None of the atoms {atom_names} found in residue {residue_number} of model {model_number}.{sub_model_number}"
        )

    return atom, atom_name

def build_molecule_from_smiles(smiles_str):
    rc("open smiles:" + smiles_str)
    rc("addh")

import csv

csv_file = 'naringenin_analogues.txt'
smiles_list = []

with open(csv_file, 'r') as file:
    reader = csv.reader(file)
    entries = next(reader)  # Read the first row containing the entries

    for i, entry in enumerate(entries):
        if i == 0 or (i) % 11 == 0:  # Gather entries at indices 0, 11, 23, etc.
            index = str(i)  # Convert the integer index to string format
            smiles_list.append(entry.strip().strip("'"))
print(smiles_list)

for index, smiles_str in enumerate(smiles_list, start=1591):
    name = str(index)
    build_molecule_from_smiles(smiles_str)
    minimize(steps)
    open_receptor(receptors_dir, receptor_names[0])
    perform_docking(vina_location, xyz_cords, xyz_size, output_location, name)

#     rc("color green #2.1")
#     # rc("open " + receptors_dir + '/mev kinase with mva.pdb')
#     # rc("focus #1")

#     # receptor, receptor_atom_name = get_chain("1", "0", 314, "B", "CA")  # Change the 3rd argument for residue number of receptor atom
#     # ligand, ligand_atom_name = get_atom("2", "1", 1, ["C1", "C2", "C3"])
# # Measure distance
#     # point1 = receptor.xformCoord()
#     # point2 = ligand.xformCoord()
#     # distance = chimera.distance(point1, point2)
#     #select #0:314.A
#     # Save image
#     rc("focus #2.1")
#     rc("turn y 110")
#     rc("scale 0.35")

    # # Check if distance is greater than 10 Angstroms
    # if distance > 10:
    #     # Show distance measurement
    #     rc("distance #1:314.B #2.1@" + str(ligand_atom_name))
    #     rc("copy file " + "RERUN" + name + "_" + os.path.splitext(receptor_names[0])[0] + ".png")
    # else:
    #     rc("distance #1:314.B #2.1@" + str(ligand_atom_name))
    #     rc("copy file " + "test" + name + ".png")
        # rc("save " + name + "_" + os.path.splitext(receptor_names[0])[0])
    # rc("copy file " + "test" + name + ".png")
    # Close all models
    rc("close session")