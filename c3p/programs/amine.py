"""
Classifies: CHEBI:32952 amine
"""
from rdkit import Chem

def is_amine(smiles: str):
    """
    Determines if a molecule is an amine based on its SMILES string.
    An amine is a compound formally derived from ammonia by replacing one,
    two, or three hydrogen atoms by hydrocarbyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS patterns for primary, secondary, tertiary amines:
    primary_amine_pattern = Chem.MolFromSmarts("[NX3;H2][#6]")
    secondary_amine_pattern = Chem.MolFromSmarts("[NX3;H1]([#6])[#6]")
    tertiary_amine_pattern = Chem.MolFromSmarts("[NX3]([#6])([#6])[#6]")
    
    # Exclusion patterns
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=O)")
    imine_pattern = Chem.MolFromSmarts("[NX2]=[CX3]")
    nitrile_pattern = Chem.MolFromSmarts("[NX1]#[CX2]")
    
    # If a nitrogen is amide, imine, or nitrile, we exclude it
    if mol.HasSubstructMatch(amide_pattern) or mol.HasSubstructMatch(imine_pattern) or mol.HasSubstructMatch(nitrile_pattern):
        return False, "Structure matches amide, imine, or nitrile, excluding it as an amine"

    # Functions to check if there is any amine pattern match
    def has_amine_structure(pattern):
        matches = mol.GetSubstructMatches(pattern)
        for match in matches:
            if not mol.GetRingInfo().IsAtomInRingOfSize(match[0], 5):
                # Check if the matched nitrogen is part of a ring
                if not mol.GetRingInfo().NumAtomRings(match[0]):
                    return True
        return False
    
    # Check if any primary, secondary, or tertiary amine structure matches
    if has_amine_structure(primary_amine_pattern) or has_amine_structure(secondary_amine_pattern) or has_amine_structure(tertiary_amine_pattern):
        return True, "Contains amine group not part of excluded structures"
    
    return False, "No amine group found, or nitrogen not bonded as required"