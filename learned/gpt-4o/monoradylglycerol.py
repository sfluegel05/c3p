"""
Classifies: CHEBI:76575 monoradylglycerol
"""
from rdkit import Chem

def is_monoradylglycerol(smiles: str):
    """
    Determines if a molecule is a monoradylglycerol based on its SMILES string.
    A monoradylglycerol has a glycerol backbone with a single acyl, alkyl, or alk-1-enyl 
    substituent at an unspecified position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoradylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string 
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Identify glycerol backbone pattern: C-C-C with 2 or 3 hydroxyl groups allowed
    glycerol_pattern = Chem.MolFromSmarts("C(CO)CO")
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    
    if not glycerol_matches:
        return False, "No glycerol backbone found"
    
    # Verify presence of exactly one substituent
    # Patterns for substituents
    acyl_pattern = Chem.MolFromSmarts("C(=O)O")
    alkyl_pattern = Chem.MolFromSmarts("CC")
    alk1_enyl_pattern = Chem.MolFromSmarts("C=CC")
    
    # Check for the number of matching substituents
    acyl_count = len(mol.GetSubstructMatches(acyl_pattern))
    alkyl_count = len(mol.GetSubstructMatches(alkyl_pattern))
    alk1_enyl_count = len(mol.GetSubstructMatches(alk1_enyl_pattern))
    
    # Ensure only one type of substituent is detected
    if sum([acyl_count, alkyl_count, alk1_enyl_count]) != 1:
        return False, f"Expected one substituent, found {sum([acyl_count, alkyl_count, alk1_enyl_count])}"
    
    # Ensure only one total substituent
    if acyl_count > 1 or alkyl_count > 1 or alk1_enyl_count > 1:
        return False, f"Found multiple substituents of the same type"
    
    return True, "Contains a glycerol backbone with exactly one lipid substituent"