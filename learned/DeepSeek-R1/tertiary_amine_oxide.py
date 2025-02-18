"""
Classifies: CHEBI:134363 tertiary amine oxide
"""
"""
Classifies: CHEBI:134687 tertiary amine oxide
"""
from rdkit import Chem

def is_tertiary_amine_oxide(smiles: str):
    """
    Determines if a molecule is a tertiary amine oxide based on its SMILES string.
    A tertiary amine oxide has a nitrogen atom with three organic groups and an oxygen atom with a negative charge.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tertiary amine oxide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find all nitrogen atoms with +1 charge
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7 and atom.GetFormalCharge() == 1:
            # Check if bonded to an oxygen with -1 charge
            o_found = False
            organic_substituents = 0
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8 and neighbor.GetFormalCharge() == -1:
                    o_found = True
                else:
                    # Check if the neighbor is part of an organic group (carbon-based)
                    if neighbor.GetAtomicNum() == 6:
                        organic_substituents += 1
            if o_found and organic_substituents == 3:
                return True, "Nitrogen with +1 charge bonded to O- and three organic groups"
    
    return False, "No tertiary amine oxide group detected"