"""
Classifies: CHEBI:18379 nitrile
"""
"""
Classifies: Nitrile
"""
from rdkit import Chem

def is_nitrile(smiles: str):
    """
    Determines if a molecule is a nitrile based on its SMILES string.
    A nitrile contains a carbon atom triple-bonded to a nitrogen atom (C#N),
    which is not part of a more complex group like amides, cyanates, or isocyanates.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nitrile, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define basic nitrile pattern: carbon triple-bonded to nitrogen
    nitrile_pattern = Chem.MolFromSmarts("[C]#[N]")
    
    # Check for known exclusion patterns: amide, cyanate, isocyanate, etc.
    exclusion_patterns = [
        Chem.MolFromSmarts("[$([NX3]=[CX3]-[CX2]#[N])]-[OX2]"),  # Carbamates and imides
        Chem.MolFromSmarts("[CX3]=[NX2]=[CX2]#[N]"),              # Cyanates
        Chem.MolFromSmarts("[NX3]=[CX2]#[N]"),                   # General exclusion for isocyanates
    ]
    
    # Find nitrile substructures
    matches = mol.GetSubstructMatches(nitrile_pattern)
    if not matches:
        return False, "No nitrile (C#N) group found"
    
    for excl_pattern in exclusion_patterns:
        if mol.HasSubstructMatch(excl_pattern):
            return False, "Nitrile (C#N) exists but as part of an excluded group (amide, cyanate, etc.)"

    return True, "Contains a nitrile (C#N) group not part of an excluded functional group"