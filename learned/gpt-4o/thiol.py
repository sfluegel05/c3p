"""
Classifies: CHEBI:29256 thiol
"""
from rdkit import Chem

def is_thiol(smiles: str):
    """
    Determines if a molecule is a thiol based on its SMILES string.
    A thiol is an organosulfur compound in which a thiol group (-SH) is attached to a carbon atom.

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
    
    # Define thiol group pattern: Thiol group (-SH) attached to carbon [C]-[S]
    thiol_pattern = Chem.MolFromSmarts("[C]-[SH]")
    
    if mol.HasSubstructMatch(thiol_pattern):
        return True, "Contains a thiol group (-SH) attached to a carbon atom"
    
    return False, "No thiol group attached to carbon found"

# Example usage:
# result, reason = is_thiol("CC(C)(COP(O)(O)=O)[C@@H](O)C(=O)NCCC(=O)NCCS")
# print(result, reason)