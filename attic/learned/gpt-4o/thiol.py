"""
Classifies: CHEBI:29256 thiol
"""
from rdkit import Chem

def is_thiol(smiles: str):
    """
    Determines if a molecule is a thiol based on its SMILES string.
    A thiol is an organosulfur compound in which a thiol group (-SH) is attached to a carbon atom
    of an aliphatic or aromatic moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a thiol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
    
    # Correct the pattern to detect thiol groups in both aliphatic and aromatic moieties
    # Thiol attached to non-aromatic carbon
    aliphatic_thiol_pattern = Chem.MolFromSmarts("[CX4,CX3,CX2]-[SH]")
    
    # Thiol which could be part of an aromatic ring
    aromatic_thiol_pattern = Chem.MolFromSmarts("[c]-[SH]")
    
    # Check if either pattern is a substructure within the molecule
    if mol.HasSubstructMatch(aliphatic_thiol_pattern) or mol.HasSubstructMatch(aromatic_thiol_pattern):
        return True, "Contains a thiol group (-SH) attached to a carbon atom"
    
    return False, "No thiol group attached to carbon found"

# Example usage:
# result, reason = is_thiol("CC(C)(COP(O)(O)=O)[C@@H](O)C(=O)NCCC(=O)NCCS")
# print(result, reason)