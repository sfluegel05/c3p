"""
Classifies: CHEBI:48039 dihydroflavonols
"""
from rdkit import Chem

def is_dihydroflavonols(smiles: str):
    """
    Determines if a molecule is a dihydroflavonol based on its SMILES string.
    A dihydroflavonol is a hydroxyflavanone with a hydroxyl group at position 3 of the heterocyclic ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dihydroflavonol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define refined SMARTS pattern for dihydroflavonol
    # This pattern captures the chromanone-like core with hydroxyl group at C3
    dihydroflavonol_pattern = Chem.MolFromSmarts("O[C@H]1[C@@H](O)C(=O)Cc2c1cccc2")

    # Check if the molecule has a match
    if mol.HasSubstructMatch(dihydroflavonol_pattern):
        return True, "Contains dihydroflavonol structure"

    return False, "Does not fit dihydroflavonol pattern adequately"