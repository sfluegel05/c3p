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
    
    # Refined SMARTS patterns for primary, secondary, and tertiary amines:
    primary_amine_pattern = Chem.MolFromSmarts("[NX3;H2][#6]")
    secondary_amine_pattern = Chem.MolFromSmarts("[NX3;H1]([#6])[#6]")
    tertiary_amine_pattern = Chem.MolFromSmarts("[NX3]([#6])([#6])[#6]")
    
    # Check the molecule for amine patterns excluding amides, imines, nitriles, etc.
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=O)")
    imine_pattern = Chem.MolFromSmarts("[NX2]=[CX3]")
    nitrile_pattern = Chem.MolFromSmarts("[NX1]#[CX2]")
    
    # Checking for amine but not amides, imines, nitriles and not part of a ring
    if (not mol.HasSubstructMatch(amide_pattern) and
        not mol.HasSubstructMatch(imine_pattern) and
        not mol.HasSubstructMatch(nitrile_pattern) and
        not mol.GetRingInfo().IsAtomInRingOfSize(mol.GetSubstructMatches(primary_amine_pattern)[0][0], 5) and
        (mol.HasSubstructMatch(primary_amine_pattern) or
         mol.HasSubstructMatch(secondary_amine_pattern) or
         mol.HasSubstructMatch(tertiary_amine_pattern))):
        return True, "Contains amine group"
    
    return False, "No amine group found, or nitrogen not bonded as required for amines"