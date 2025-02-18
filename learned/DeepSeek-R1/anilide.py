"""
Classifies: CHEBI:13248 anilide
"""
"""
Classifies: CHEBI:22314 anilide
An aromatic amide obtained by acylation of aniline.
"""
from rdkit import Chem

def is_anilide(smiles: str):
    """
    Determines if a molecule is an anilide based on its SMILES string.
    An anilide has an amide group (-CONH-) where the nitrogen is directly attached to an aromatic ring.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an anilide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Find all amide groups (N connected to C=O)
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if not amide_matches:
        return False, "No amide group found"
    
    # Check each amide's nitrogen for direct aromatic attachment
    for match in amide_matches:
        nitrogen_idx = match[0]  # First atom in match is the N
        nitrogen = mol.GetAtomWithIdx(nitrogen_idx)
        
        # Check if nitrogen is bonded to an aromatic carbon
        for neighbor in nitrogen.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIsAromatic():
                return True, "Amide nitrogen attached to aromatic ring"
    
    return False, "Amide group(s) present but not attached to aromatic ring"