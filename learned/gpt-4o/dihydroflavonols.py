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

    # Define the dihydroflavonol pattern: bicyclic core with specific hydroxylation
    # Representative pattern for dihydroflavonol structural backbone with hydroxyl on position 3
    core_pattern = "[C@H]1([O])C(=O)CC2C(C=CC=C2O)O1"  # Simplified pattern focusing on key features

    # Compile SMARTS pattern to match dihydroflavonol core
    dihydroflavonol_pattern = Chem.MolFromSmarts(core_pattern)
    
    # Check if the molecule has a match
    if mol.HasSubstructMatch(dihydroflavonol_pattern):
        return True, "Contains dihydroflavonol structure"

    return False, "Does not fit dihydroflavonol pattern adequately"

# Note: This is an example pattern for educational purposes. Real-world applications
# would refine the pattern based on precise stereochemistry and functional group placement.