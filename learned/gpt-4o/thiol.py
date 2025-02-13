"""
Classifies: CHEBI:29256 thiol
"""
from rdkit import Chem

def is_thiol(smiles: str):
    """
    Determines if a molecule is a thiol based on its SMILES string.
    A thiol is an organosulfur compound with a thiol group (-SH) attached to a carbon atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a thiol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for thiol group pattern, -SH
    thiol_pattern = Chem.MolFromSmarts("[CX4,SX3][SH]")
    if mol.HasSubstructMatch(thiol_pattern):
        return True, "Contains a thiol group (-SH) attached to a carbon atom"
    
    return False, "No thiol group found"

# Examples
print(is_thiol("C(CC(CCO)(S)C)C"))  # 3-methyl-3-sulfanylhexan-1-ol
print(is_thiol("SCCCCCC"))          # 1-Hexanethiol