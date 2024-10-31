from rdkit import Chem
from rdkit.Chem import AllChem

def is_hexol(smiles: str):
    """
    Determines if a molecule is a hexol (contains exactly 6 hydroxy groups).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hexol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Count number of OH groups
    oh_pattern = Chem.MolFromSmarts('[OH]')
    oh_matches = mol.GetSubstructMatches(oh_pattern)
    
    num_oh = len(oh_matches)
    
    if num_oh != 6:
        return False, f"Contains {num_oh} hydroxy groups, not 6"
    
    # Get the indices of all O atoms that are part of OH groups
    oh_o_indices = [match[0] for match in oh_matches]
    
    # Check if the OH groups are bonded to carbons
    for o_idx in oh_o_indices:
        o_atom = mol.GetAtomWithIdx(o_idx)
        neighbors = [n for n in o_atom.GetNeighbors() if n.GetAtomicNum() != 1] # exclude H
        if not neighbors or neighbors[0].GetAtomicNum() != 6:  # C has atomic number 6
            return False, "At least one OH group is not bonded to carbon"
            
    return True, "Contains exactly 6 hydroxy groups bonded to carbon atoms"
# Pr=1.0
# Recall=1.0