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

    # Define SMARTS pattern for the 4-thia-1-azabicyclo[3.2.0]heptane core in penicillins
    penicillin_core_pattern = Chem.MolFromSmarts("N1C(=O)[C@@H]2[C@@H](S1)[C@]2(C(=O)O)C")
    if not mol.HasSubstructMatch(penicillin_core_pattern):
        return False, "Core penicillin structure not found"

    # Check for two methyl substituents at position 2 of the core
    methyl_pattern = Chem.MolFromSmarts("C(C)(C)C1CNC(=O)[C@@H](S1)")
    if not mol.HasSubstructMatch(methyl_pattern):
        return False, "Methyl groups at position 2 are missing"

    # Check for a carboxylate group at position 3
    carboxylate_pattern = Chem.MolFromSmarts("C(=O)[O-]")
    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "Required carboxylate group not found at position 3"

    # Check for a carboxamido group at position 6
    carboxamido_pattern = Chem.MolFromSmarts("NC(=O)")
    if not mol.HasSubstructMatch(carboxamido_pattern):
        return False, "Required carboxamido group not found at position 6"

    return True, "Matches penicillin structure requirements"