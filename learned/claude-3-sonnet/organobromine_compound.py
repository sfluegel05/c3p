"""
Classifies: CHEBI:37141 organobromine compound
"""
"""
Classifies: CHEBI:33636 Organobromine compounds
An organobromine compound is a compound containing at least one carbon-bromine bond.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_organobromine_compound(smiles: str):
    """
    Determines if a molecule is an organobromine compound based on its SMILES string.

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
    
    # Check for presence of carbon-bromine bonds
    has_c_br_bond = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 35:  # Bromine
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6:  # Carbon
                    has_c_br_bond = True
                    break
            if has_c_br_bond:
                break
    
    if has_c_br_bond:
        return True, "Contains at least one carbon-bromine bond"
    else:
        return False, "No carbon-bromine bonds found"