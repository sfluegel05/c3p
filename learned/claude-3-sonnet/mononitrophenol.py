"""
Classifies: CHEBI:39362 mononitrophenol
"""
"""
Classifies: CHEBI:51088 mononitrophenol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_mononitrophenol(smiles: str):
    """
    Determines if a molecule is a mononitrophenol based on its SMILES string.
    A mononitrophenol is a phenol carrying a single nitro substituent at unspecified position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mononitrophenol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # First find all aromatic rings with OH substituents (phenols)
    # This pattern matches any aromatic ring with an OH group
    phenol_pattern = Chem.MolFromSmarts("c1([OH])ccccc1")
    phenol_matches = mol.GetSubstructMatches(phenol_pattern)
    
    if not phenol_matches:
        return False, "No phenol group found"

    # Check for nitro groups - both resonance forms
    nitro_patterns = [
        Chem.MolFromSmarts("[N+](=[O-])=O"),  # Charged form
        Chem.MolFromSmarts("[N+]([O-])=O"),    # Alternative charged form
        Chem.MolFromSmarts("N(=O)=O")          # Uncharged form
    ]
    
    total_nitro = 0
    nitro_matches = []
    for pattern in nitro_patterns:
        if pattern is not None:
            matches = mol.GetSubstructMatches(pattern)
            total_nitro += len(matches)
            nitro_matches.extend(matches)

    if total_nitro == 0:
        return False, "No nitro group found"
    elif total_nitro > 1:
        return False, f"Found {total_nitro} nitro groups, must have exactly one"

    # For each phenol ring found, check if it has exactly one nitro group attached
    for phenol_match in phenol_matches:
        phenol_atoms = set(phenol_match)
        
        # Check if the nitro group is attached to this specific phenol ring
        for nitro_match in nitro_matches:
            nitro_N = nitro_match[0]  # Get the nitrogen atom index
            nitro_neighbors = mol.GetAtomWithIdx(nitro_N).GetNeighbors()
            
            for neighbor in nitro_neighbors:
                if neighbor.GetIdx() in phenol_atoms:
                    # Found a nitro group attached to this phenol ring
                    return True, "Contains phenol ring with single nitro substituent"

    return False, "Nitro group not properly attached to phenol ring"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:51088',
        'name': 'mononitrophenol',
        'definition': 'A nitrophenol that is phenol carrying a single nitro substituent at unspecified position.',
        'parents': ['CHEBI:33622']
    }
}