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

    # Improved penicillin scaffold: 4-thia-1-azabicyclo[3.2.0]heptane core
    penicillin_scaffold_pattern = Chem.MolFromSmarts("C1([C@H]2SC3N2C(=O)C([C@@H]3C1(=O)O)(C)C)C(=O)OX")  # X here represents possible esters or amides
    if not mol.HasSubstructMatch(penicillin_scaffold_pattern):
        return False, "Penicillin scaffold not found"

    # Verify the presence of the two methyl groups at position 2
    methylation_pattern = Chem.MolFromSmarts("SC(C)(C)C")
    if not mol.HasSubstructMatch(methylation_pattern):
        return False, "Required methyl groups at position 2 not found"

    # Carboxylate should be present, assuming variable position due to diversity
    carboxylate_pattern = Chem.MolFromSmarts("C(=O)[O-]")
    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "Carboxylate group not found"

    # Carboxamido group as part of the peptidic bond
    carboxamido_pattern = Chem.MolFromSmarts("NC(=O)")
    if not mol.HasSubstructMatch(carboxamido_pattern):
        return False, "Carboxamido group not found"

    return True, "SMILES string represents a penicillin"