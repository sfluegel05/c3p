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
    A diradylglycerol is characterized by a glycerol backbone with two acyl, alkyl, or alk-1-enyl linkages.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a diradylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a broader pattern for the glycerol backbone with potential stereochemical configurations
    glycerol_pattern = Chem.MolFromSmarts("[C@H](O)[C@H](O)CO | [C@@H](O)[C@@H](O)CO | O[C@H](CO)[C@H]O")
    
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Consider more comprehensive patterns for acyl, alkyl, and alk-1-enyl linkages
    acyl_pattern = Chem.MolFromSmarts("C(=O)[O;R0]")  # more complex ester types
    alkyl_pattern = Chem.MolFromSmarts("O[CX4,CX3]")  # broader ether patterns
    alkenyl_pattern = Chem.MolFromSmarts("O=C[C]=[C]")  # conjugated systems across linkages

    # Count linkages
    acyl_count = len(mol.GetSubstructMatches(acyl_pattern))
    alkyl_count = len(mol.GetSubstructMatches(alkyl_pattern))
    alkenyl_count = len(mol.GetSubstructMatches(alkenyl_pattern))
    
    linkage_count = acyl_count + alkyl_count + alkenyl_count

    if linkage_count < 2:
        return False, f"Found only {linkage_count} applicable linkages, need at least 2."

    return True, "Contains glycerol backbone with two acyl, alkyl, or alk-1-enyl linkages"