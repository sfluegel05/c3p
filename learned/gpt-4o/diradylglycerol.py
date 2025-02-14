"""
Classifies: CHEBI:76578 diradylglycerol
"""
"""
Classifies: CHEBI:17977 diradylglycerol
"""
from rdkit import Chem

def is_diradylglycerol(smiles: str):
    """
    Determines if a molecule is a diradylglycerol based on its SMILES string.
    A diradylglycerol is defined as a glycerol backbone (three-carbon chain) with two acyl, alkyl, or alk-1-enyl linkages.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diradylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a pattern for glycerol backbone where backbone can have various linkages
    glycerol_pattern = Chem.MolFromSmarts("C(O)(O)C")
    
    # Check for glycerol backbone
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Look for acyl, alkyl, and alk-1-enyl linkages
    acyl_pattern = Chem.MolFromSmarts("[C](=O)O")  # ester linkage (acyl chain)
    alkyl_pattern = Chem.MolFromSmarts("[O][CX4]")  # ether linkage (alkyl chain)
    alkenyl_pattern = Chem.MolFromSmarts("[CH]=[CH]-[C]") # alk-1-enyl ether

    # Count matches
    acyl_matches = len(mol.GetSubstructMatches(acyl_pattern))
    alkyl_matches = len(mol.GetSubstructMatches(alkyl_pattern))
    alkenyl_matches = len(mol.GetSubstructMatches(alkenyl_pattern))
    
    linkage_count = acyl_matches + alkyl_matches + alkenyl_matches
    if linkage_count < 2:
        return False, f"Found only {linkage_count} applicable linkages, need at least 2"

    return True, "Contains glycerol backbone with two acyl, alkyl, or alk-1-enyl linkages"