"""
Classifies: CHEBI:37141 organobromine compound
"""
"""
Classifies: CHEBI:38241 organobromine compound
"""
from rdkit import Chem

def is_organobromine_compound(smiles: str):
    """
    Determines if a molecule is an organobromine compound based on its SMILES string.
    An organobromine compound is a compound containing at least one carbon-bromine bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organobromine compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for bromine atoms
    bromine_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 35]
    if not bromine_atoms:
        return False, "No bromine atoms present"

    # Check for carbon-bromine bonds (direct and in rings)
    has_c_br_bond = False
    for br_atom in bromine_atoms:
        for neighbor in br_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # Carbon atom
                has_c_br_bond = True
                break
        if has_c_br_bond:
            break

    if has_c_br_bond:
        return True, "Contains at least one carbon-bromine bond"
    else:
        return False, "No carbon-bromine bonds found"