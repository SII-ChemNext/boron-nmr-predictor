from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd

# Analyze the first 5 molecules containing Morgan bit 690
smiles_list = [
    "C[SiH](C)C1=C(B2C3CCCC2CCC3)CCC[Si]1(C)C",
    "C/C(=C1\C[C@@H]2C[C@@H]3CB1C[C@@H](C3)C2)[Si](C)(C)C",
    "C[Si](C)(C)C=C=C1C=C([Si](C)(C)C)C23C[C@@H]4C[C@H](CB12)C[C@@H](C4)C3",
    "C=CCB1CCC[Si](C)(C)C([SiH](C)C)=C1CC=C",
    "C[Si](C)(C)[C@H]1[C@H]([Si](C)(C)C)B(Cl)B(Cl)[C@@H]1[Si](C)(C)C"
]

print("Analyzing atoms that trigger Morgan bit 690 across molecules:\n")

for i, smiles in enumerate(smiles_list, 1):
    mol = Chem.MolFromSmiles(smiles)
    bit_info = {}
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024, bitInfo=bit_info)
    
    print(f"Molecule {i}: {smiles[:50]}...")
    
    if 690 in bit_info:
        for atom_idx, radius in bit_info[690]:
            atom = mol.GetAtomWithIdx(atom_idx)
            print(f"  Bit 690 triggered at: atom {atom_idx} ({atom.GetSymbol()}), radius={radius}")
            
            # Check the neighbors of this atom
            neighbors = [mol.GetAtomWithIdx(n.GetIdx()).GetSymbol() for n in atom.GetNeighbors()]
            print(f"    Neighbor atoms: {neighbors}")
    else:
        print(f"  ⚠ Bit 690 not found")
    print()
