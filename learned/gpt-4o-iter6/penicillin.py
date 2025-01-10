"""
Classifies: CHEBI:17334 penicillin
"""
from rdkit import Chem

def is_penicillin(smiles: str):
    """
    Determines if a molecule is a penicillin based on its SMILES string.
    A penicillin refers to the presence of a 4-thia-1-azabicyclo[3.2.0]heptane
    core with specific substituents: two methyl groups at position 2,
    a carboxylate at position 3, and a carboxamido group at position 6.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a penicillin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for the core penicillin scaffold: 4-thia-1-azabicyclo[3.2.0]heptane
    penicillin_core_pattern = Chem.MolFromSmarts("C1([C@@H]2SC3N2C(=O)[C@@H](C(=O)O)C3)C")
    if penicillin_core_pattern is None:
        return None, None

    # Check for the core penicillin structure and required stereochemistry
    if not mol.HasSubstructMatch(penicillin_core_pattern):
        return False, "Core penicillin structure not found"

    # Verify two methyl groups at position 2 of the penicillin core
    methyl_pattern = Chem.MolFromSmarts("[C@@]([C@](C)(C)C1([C@@H]2SC3N2)[N](C(=O))[CH](C3)C(=O)O)=O")
    if not mol.HasSubstructMatch(methyl_pattern):
        return False, "Methyl groups at position 2 are required but not found"

    # Carboxylate group should be present: C(=O)O moiety at the outward-position
    carboxylate_pattern = Chem.MolFromSmarts("C(=O)[O-]")
    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "Required carboxylate group not found"

    # Carboxamido group at position 6 of the penicillin backbone
    carboxamido_pattern = Chem.MolFromSmarts("NC(=O)")
    if not mol.HasSubstructMatch(carboxamido_pattern):
        return False, "Required carboxamido group not found"

    return True, "Matches penicillin structure requirements"