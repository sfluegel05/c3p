"""
Classifies: CHEBI:47908 alkanethiol
"""
from rdkit import Chem

def is_alkanethiol(smiles: str):
    """
    Determines if a molecule is an alkanethiol based on its SMILES string.
    An alkanethiol has a sulfanyl group (-SH) attached directly to an alkyl group (carbon).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an alkanethiol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Thiol group pattern - sulfur with hydrogen and a bonding carbon
    # This searches for sulfur with one attached hydrogen and a connection to a carbon atom
    thiol_pattern = Chem.MolFromSmarts("[SX2][H]C")
    
    if not mol.HasSubstructMatch(thiol_pattern):
        return False, "No thiol (-SH) group directly attached to carbon found"
    
    return True, "Contains alkanethiol pattern with -SH group attached to carbon"

# Example usage
print(is_alkanethiol("SCC(CC)C"))  # Example SMILES for 2-Methyl-1-butanethiol