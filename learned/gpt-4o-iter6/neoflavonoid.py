"""
Classifies: CHEBI:71971 neoflavonoid
"""
from rdkit import Chem

def is_neoflavonoid(smiles: str):
    """
    Determines if a molecule is a neoflavonoid based on its SMILES string.
    A neoflavonoid is any 1-benzopyran with an aryl substituent at position 4.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a neoflavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define 1-benzopyran pattern 
    benzopyran_pattern = Chem.MolFromSmarts("c1cc2ccocc2c1")
    if not mol.HasSubstructMatch(benzopyran_pattern):
        return False, "No 1-benzopyran structure found"

    # Identify atoms that match the 1-benzopyran pattern
    benzopyran_matches = mol.GetSubstructMatches(benzopyran_pattern)

    # For each match, check if there is an aryl ring (aromatic six-membered ring) substituent connected at the 4-position
    for match in benzopyran_matches:
        # The atom at index 4 in the SMARTS string corresponds to the position where aryl group should be attached
        aryl_position_idx = match[4]
        aryl_position_atom = mol.GetAtomWithIdx(aryl_position_idx)

        # Scan the neighbors for an aryl group
        for neighbor in aryl_position_atom.GetNeighbors():
            # Define an aryl group as an aromatic six-membered ring
            if Chem.MolFromSmiles(Chem.MolToSmiles(neighbor.GetOwningMol())).HasSubstructMatch(Chem.MolFromSmarts('c1ccccc1')):
                return True, "Aryl substituent found at position 4"

    return False, "No aryl substitution at position 4 in the 1-benzopyran"