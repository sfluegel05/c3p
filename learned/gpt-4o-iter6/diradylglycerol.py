"""
Classifies: CHEBI:76578 diradylglycerol
"""
from rdkit import Chem

def is_diradylglycerol(smiles: str):
    """
    Determines if a molecule is a diradylglycerol based on its SMILES string.
    A diradylglycerol has a glycerol backbone with exactly two substituent groups (acyl, alkyl, or alk-1-enyl).

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
    
    # Define improved glycerol backbone pattern to account for stereochemistry
    glycerol_pattern = Chem.MolFromSmarts("[CH2][CH](O)[CH2]O")
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    if not glycerol_matches:
        return False, "Glycerol backbone pattern not matched"

    # Identify acyl groups (R-C(=O)-O-)
    acyl_pattern = Chem.MolFromSmarts("C(=O)O[CH2]")
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)

    # Identify alkyl groups (R-O-)
    alkyl_pattern = Chem.MolFromSmarts("[CH2]O")
    alkyl_matches = mol.GetSubstructMatches(alkyl_pattern)

    # Identify alk-1-enyl groups
    alk1enyl_pattern = Chem.MolFromSmarts("C=C[CH2]O")
    alk1enyl_matches = mol.GetSubstructMatches(alk1enyl_pattern)

    # Combine the indices from matches ensuring unique attachment checking
    substituents_count = len(set(
        [match[2] for match in acyl_matches] + 
        [match[0] for match in alkyl_matches] +
        [match[0] for match in alk1enyl_matches]
    ))
    
    if substituents_count != 2:
        return False, f"Expected exactly 2 substituent groups, found {substituents_count}"

    return True, "Valid diradylglycerol with two substituents connected via acyl, alkyl, or alk-1-enyl bonds"