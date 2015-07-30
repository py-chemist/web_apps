from rdkit import Chem
from rdkit.Chem import Descriptors

m = Chem.MolFromMolFile('molecule.mol')
smiles = Chem.MolToSmiles(m)
print smiles