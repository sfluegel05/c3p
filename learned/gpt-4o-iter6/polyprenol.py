"""
Classifies: CHEBI:26199 polyprenol
"""
from rdkit import Chem

def is_polyprenol(smiles: str):
    """
    Determines if a molecule is a polyprenol based on its SMILES string.
    A polyprenol has more than one isoprene unit with a terminal hydroxyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyprenol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern for detecting segments of isoprene units: C(=C/C)C, consider trans configurations
    isoprene_pattern = Chem.MolFromSmarts("C(=C)C")  
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    
    if len(isoprene_matches) <= 1:
        return False, f"Found {len(isoprene_matches)} isoprene units, need more than 1"
    
    # SMARTS pattern for terminal hydroxyl group: Simplified for terminal carbon with OH
    oh_pattern = Chem.MolFromSmarts("[CH2]O")
    terminal_matches = mol.GetSubstructMatches(oh_pattern)
    
    if not terminal_matches:
        return False, "No appropriate terminal hydroxyl group"
    
    # Additional validation: check if the OH is connected to the end of a chain of isoprene units
    for oh_match in terminal_matches:
        oh_atom_idx = oh_match[1]
        is_terminal = True
        
        # Check if the hydroxyl group connects directly to a continuous chain of isoprene units
        oh_atom = mol.GetAtomWithIdx(oh_atom_idx)
        for neighbor in oh_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # Carbon check
                connected_isoprene = [match for match in isoprene_matches if neighbor.GetIdx() in match]
                if not connected_isoprene:
                    is_terminal = False
        
        if is_terminal:
            return True, "Verified: Contains more than one isoprene unit with a terminal hydroxyl group"
    
    # If no valid terminal hydroxyl connection found
    return False, "Hydroxyl is not positioned terminally relative to isoprene chain"