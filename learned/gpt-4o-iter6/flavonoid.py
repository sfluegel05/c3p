"""
Classifies: CHEBI:47916 flavonoid
"""
from rdkit import Chem

def is_flavonoid(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    A flavonoid is characterized by a 1-benzopyran structure with an aryl group at position 2.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Adjusted pattern to include variations for 1-benzopyran (flavonoid core)
    benzopyran_pattern = Chem.MolFromSmarts("c1c2ccccc2oc1")
    if not mol.HasSubstructMatch(benzopyran_pattern):
        return False, "No flavonoid core (1-benzopyran structure) found"

    # Alternative approach to find the 1-benzopyran core without specific aryl substitution matching
    # Allow for any aryl substitution as long as 1-benzopyran core is intact
    # Here, we are assuming any group attached to position 2 of benzopyran is acceptable as a substitution
    aryl_substituted_pattern = Chem.MolFromSmarts("C1=COc2ccccc2C1-*")  # Using * to represent any possible substitution
    if not mol.HasSubstructMatch(aryl_substituted_pattern):
        return False, "Aryl substitution not found at position 2"

    # Hydroxyl groups are commonly found but not mandatory
    # hydroxyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() == 1)
    # if hydroxyl_count < 1:
    #     return False, "Insufficient hydroxyl groups for a typical flavonoid structure"

    return True, "Identified flavonoid structure with 1-benzopyran core and aryl substitution"