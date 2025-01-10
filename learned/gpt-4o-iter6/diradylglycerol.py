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

    # Look for a generalized glycerol backbone pattern
    glycerol_pattern = Chem.MolFromSmarts("C(CO)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "Glycerol backbone pattern not matched"

    # Identify acyl groups (R-C(=O)-O-)
    acyl_pattern = Chem.MolFromSmarts("C(=O)O")
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)

    # Identify alkyl groups (R-O-) that aren't part of acyl groups
    alkyl_pattern = Chem.MolFromSmarts("CO")
    alkyl_matches = mol.GetSubstructMatches(alkyl_pattern)

    # Ensure that ether link is not double-counted by excluding acyl's matched terminal carbon
    unique_alkyl_positions = set(m[0] for m in alkyl_matches) - set(m[0] for m in acyl_matches)

    # Identify alk-1-enyl groups
    alk1enyl_pattern = Chem.MolFromSmarts("C=C-O")
    alk1enyl_matches = mol.GetSubstructMatches(alk1enyl_pattern)

    # Count substituents before applying stereochemistry considerations
    total_unique_substituents = (
        len(set([match[0] for match in acyl_matches])) +
        len(unique_alkyl_positions) +
        len(set([match[0] for match in alk1enyl_matches]))
    )

    if total_unique_substituents != 2:
        return False, f"Expected exactly 2 substituent groups, found {total_unique_substituents}"

    return True, "Valid diradylglycerol with two substituents connected via acyl, alkyl, or alk-1-enyl bonds"