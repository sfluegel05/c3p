"""
Classifies: CHEBI:35359 carboxamidine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_carboxamidine(smiles: str):
    """
    Determines if a molecule is a carboxamidine based on its SMILES string.
    A carboxamidine has the structure RC(=NR)NR2, where R can be H.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carboxamidine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Core carboxamidine pattern
    core_pattern = Chem.MolFromSmarts("C(=N)[N]")
    matches = mol.GetSubstructMatches(core_pattern)
    
    if not matches:
        return False, "Core carboxamidine structure C(=N)N not found."
    
    for match in matches:
        # Get the indices of the matched atoms
        c_idx, n1_idx, n2_idx = match

        # check n2 substitution
        n2_atom = mol.GetAtomWithIdx(n2_idx)
        n2_bonds = n2_atom.GetTotalValence()
        if n2_bonds > 3:
            return False, f"Nitrogen {n2_idx} is bonded to more than 3 atoms"
        
        # check c substitution
        c_atom = mol.GetAtomWithIdx(c_idx)
        c_bonds = c_atom.GetTotalValence()
        if c_bonds > 4:
           return False, f"Carbon {c_idx} is bonded to more than 4 atoms"

    return True, "Molecule matches the core structure of a carboxamidine"