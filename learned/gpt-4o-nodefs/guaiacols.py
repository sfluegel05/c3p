"""
Classifies: CHEBI:134251 guaiacols
"""
"""
Classifies: Guaiacols
"""
from rdkit import Chem

def is_guaiacols(smiles: str):
    """
    Determines if a molecule is a guaiacol based on its SMILES string.
    Guaiacols typically have a methoxy group and a hydroxyl group attached to a benzene ring,
    possibly in ortho, meta, or para positions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a guaiacol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Expanded guaiacol pattern allowing for ortho, meta, and para hydroxy-methoxy pairings
    patterns = [
        "c1cc(OC)c(O)cc1",     # ortho (positions 1-2 on benzene)
        "c1c(OC)cc(O)cc1",     # meta (positions 1-3 on benzene)
        "c1cc(O)cc(OC)c1"      # para (positions 1-4 on benzene)
    ]
    
    # Check for multiple patterns to allow flexibility in structural identification
    for pattern in patterns:
        mol_pattern = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(mol_pattern):
            return True, "Contains guaiacol core structure"
    
    return False, "No guaiacol core structure found"