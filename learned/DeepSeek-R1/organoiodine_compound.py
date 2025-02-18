"""
Classifies: CHEBI:37142 organoiodine compound
"""
"""
Classifies organoiodine compounds (CHEBI:37139) - compounds containing at least one carbon-iodine bond.
"""
from rdkit import Chem

def is_organoiodine_compound(smiles: str):
    """
    Determines if a molecule is an organoiodine compound based on the presence of at least one carbon-iodine bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains C-I bond, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check all iodine atoms to see if any are bonded to carbon
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 53:  # Iodine
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6:  # Carbon
                    return True, "Contains carbon-iodine bond"
    
    # If no C-I bonds found, return False
    return False, "No carbon-iodine bonds found"