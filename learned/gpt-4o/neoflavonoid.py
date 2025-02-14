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
        bool: True if molecule is a neoflavonoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Refined 1-benzopyran core pattern
    benzopyran_pattern = Chem.MolFromSmarts('c1cc2oc(C)c(C=CC2=c1)')  # Open valences allow detection of potential site for aryl group
    if not mol.HasSubstructMatch(benzopyran_pattern):
        return False, "1-benzopyran core not found"
    
    # Refined pattern for aryl group at position 4
    aryl_pattern = Chem.MolFromSmarts('c2ccc(-c1coc3ccccc3c1)c2')  # General pattern for an aryl group connected to the core
    if not mol.HasSubstructMatch(aryl_pattern):
        return False, "Aryl group not found at position 4"
    
    return True, "Contains a 1-benzopyran core with an aryl group at position 4"