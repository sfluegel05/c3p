"""
Classifies: CHEBI:26493 quinic acid
"""
from rdkit import Chem

def is_quinic_acid(smiles: str):
    """
    Determines if a molecule is quinic acid or one of its derivatives based on its SMILES string.
    Quinic acid is characterized as a cyclitol carboxylic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as quinic acid or a derivative, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the pattern for a hydroxylated cyclohexane with stereochemistry considerations
    hydroxy_cyclohexane_pattern = Chem.MolFromSmarts("[C@H]1([C@H](O)[C@H](O)[C@H](O)[C@H]1)C(=O)O")
    
    if not mol.HasSubstructMatch(hydroxy_cyclohexane_pattern):
        return False, "No cyclohexane with correct stereochemistry and hydroxyl groups found"

    # Look for presence of ester or ether linkages connected to aromatic systems
    ester_linkage_pattern = Chem.MolFromSmarts("[C](=O)O[C@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1")
    if mol.HasSubstructMatch(ester_linkage_pattern):
        return True, "Quinic acid derivative with ester linkage detected"

    # Assert cyclohexanecarboxylic acid main backbone
    cyclohexanecarboxylic_acid_pattern = Chem.MolFromSmarts("[C@H]1([C@H](O)[C@H](O)[C@H](O)[C@H]1)C(=O)[O]C")
    if mol.HasSubstructMatch(cyclohexanecarboxylic_acid_pattern):
        return True, "Cyclohexanecarboxylic acid backbone identified"

    return False, "Does not meet quinic acid structural criteria"