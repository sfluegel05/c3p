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

    # Define SMARTS pattern for more precise stereochemistry in the core structure
    penicillin_core_pattern = Chem.MolFromSmarts("[C@@H]1([S][C@@H]2C[N@]1C(=O)[C@H]2C(=O)O)C")
    if not mol.HasSubstructMatch(penicillin_core_pattern):
        return False, "Core penicillin structure not found"

    # Check for two methyl substituents attached to the penam core (position 2)
    methyl_pattern = Chem.MolFromSmarts("[C@@](C)(C)C1([S][C@@H]2C[N@]1C(=O)[C@H]2C(=O)O)")
    if not mol.HasSubstructMatch(methyl_pattern):
        return False, "Methyl groups at position 2 are missing"

    # Check for carboxylate group presence - flexible matching to accommodate ionization states
    carboxylate_pattern = Chem.MolFromSmarts("C(=O)[O]")
    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "Required carboxylate group not found at position 3"

    # Check for a carboxamido group in a reasonable proximity - flexible approach
    carboxamido_pattern = Chem.MolFromSmarts("NC(=O)")
    if not mol.HasSubstructMatch(carboxamido_pattern):
        return False, "Required carboxamido group not found at position 6"

    return True, "Matches penicillin structure requirements"