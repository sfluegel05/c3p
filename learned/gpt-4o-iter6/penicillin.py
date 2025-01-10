"""
Classifies: CHEBI:17334 penicillin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_penicillin(smiles: str):
    """
    Determines if a molecule is a penicillin based on its SMILES string.
    A penicillin is defined as a substituted penam with specific structural features.

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

    # The basic penicillin scaffold: 4-thia-1-azabicyclo[3.2.0]heptane.
    # With two methyl groups at position 2, a carboxylate at position 3, and a carboxamido at position 6.
    penicillin_scaffold_pattern = Chem.MolFromSmarts("[C@@H]1([C@@H]2N([C@H]1C(=O)O)C(=O)S2(C)C)")
    if not mol.HasSubstructMatch(penicillin_scaffold_pattern):
        return False, "Penicillin scaffold not found"
    
    # Ensure methylation at position 2
    methylation_pattern = Chem.MolFromSmarts("S1[C@@]2([C@H]1C)C")
    methylation_matches = mol.GetSubstructMatches(methylation_pattern)
    if len(methylation_matches) == 0:
        return False, "Required methyl groups at position 2 not found"
    
    # Carboxylate group must be present on position 3
    carboxylate_pattern = Chem.MolFromSmarts("C(=O)O")
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_pattern)
    if len(carboxylate_matches) == 0:
        return False, "Carboxylate group not found at required position"
    
    # Carboxamido group must be present on position 6
    carboxamido_pattern = Chem.MolFromSmarts("NC(=O)")
    carboxamido_matches = mol.GetSubstructMatches(carboxamido_pattern)
    if len(carboxamido_matches) == 0:
        return False, "Carboxamido group not found at required position"

    return True, "SMILES string represents a penicillin"