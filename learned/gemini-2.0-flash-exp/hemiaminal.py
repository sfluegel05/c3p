"""
Classifies: CHEBI:73080 hemiaminal
"""
from rdkit import Chem

def is_hemiaminal(smiles: str):
    """
    Determines if a molecule is a hemiaminal based on its SMILES string.
    A hemiaminal has a carbon bonded to both a hydroxyl (-OH) and an amino (-NH2 or substituted amine) group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hemiaminal, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate over all carbon atoms
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Check if it is a carbon atom
            
            has_oh = False
            has_nhx = False

            # Check the neighbors for -OH and -NHx
            for neighbor in atom.GetNeighbors():
                neighbor_atomic_num = neighbor.GetAtomicNum()
                if neighbor_atomic_num == 8: #check if it is oxygen
                    has_oh = True
                elif neighbor_atomic_num == 7: #check if it is nitrogen
                    has_nhx = True

            if has_oh and has_nhx:
                return True, "Found a carbon bonded to both a hydroxyl and an amine group"

    return False, "No hemiaminal structure found"