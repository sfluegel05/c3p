"""
Classifies: CHEBI:26244 prenols
"""
from rdkit import Chem

def is_prenols(smiles: str):
    """
    Determines if a molecule is a prenol based on its SMILES string.
    A prenol is an alcohol possessing one or more isoprene units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prenol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # More flexible isoprene pattern: [CH2-C(Me)=CH-CH2], dealing with traditional isoprene and conjugated units
    isoprene_patterns = [
        Chem.MolFromSmarts("C(C)=C(C)"),      # General isoprene core
        Chem.MolFromSmarts("C=C(C)C"),        # Alternative representation
        Chem.MolFromSmarts("[CH2]C(=C)C"),    # Left-projected isoprene
        Chem.MolFromSmarts("C(=C)C[CH2]")     # Right-projected isoprene
    ]
    isoprene_matches = any(mol.HasSubstructMatch(pattern) for pattern in isoprene_patterns)
    
    # Check for isoprene units
    if not isoprene_matches:
        return False, "No isoprene units found"

    # Verify the presence of a terminal alcohol group (-OH)
    if not mol.HasSubstructMatch(Chem.MolFromSmarts("[CX4][OH]")):
        return False, "Missing terminal or appropriate OH group"

    # Checking linearity (no cyclic structure in the backbone)
    if not mol.GetRingInfo().IsStraightChain():
        return False, "Should not contain cyclic structures"

    return True, "Contains isoprene units with terminal alcohol group"