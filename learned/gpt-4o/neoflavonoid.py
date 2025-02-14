"""
Classifies: CHEBI:71971 neoflavonoid
"""
from rdkit import Chem

def is_neoflavonoid(smiles: str):
    """
    Determines if a molecule is a neoflavonoid based on its SMILES string.
    A neoflavonoid is defined as any 1-benzopyran with an aryl substituent at position 4.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a neoflavonoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1-benzopyran core pattern (basic chromene structure), flexible substitution allowed
    benzopyran_pattern = Chem.MolFromSmarts('c1ccc2c(c1)occ2')  # Core without specificity on the substitution
    if not mol.HasSubstructMatch(benzopyran_pattern):
        return False, "1-benzopyran core not found"

    # Locate the benzopyran core if found, then check for aryl at position 4
    # Use a defined substructure to ensure position matching
    # 4-position in benzopyran is typically part of the C-ring in flavonoids
    # e.g., Aryl directly connected to chromene/coumarin system
    aryl_position_pattern = Chem.MolFromSmarts('c1ccc(-c2ccc3c(c2)occ3)o1')  # Ensure direct connection for aryl
    if not mol.HasSubstructMatch(aryl_position_pattern):
        return False, "Aryl group not found at position 4 of benzopyran"

    return True, "Contains a 1-benzopyran core with an aryl group at position 4"