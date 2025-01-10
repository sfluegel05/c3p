"""
Classifies: CHEBI:47908 alkanethiol
"""
from rdkit import Chem

def is_alkanethiol(smiles: str):
    """
    Determines if a molecule is an alkanethiol based on its SMILES string.
    An alkanethiol is defined as a compound where a sulfanyl group (-SH) is attached to an alkyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: Whether the molecule is an alkanethiol
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern correctly for alkanethiol (simple alkyl group + thiol group)
    # A simple chain or branch ending with a thiol group
    thiol_pattern = Chem.MolFromSmarts("[CX4][SX2H1]")  # Aliphatic carbon directly attached to sulfur with hydrogen

    # We assume we're looking for an aliphatic (sp3) carbon bonded to the thiol group
    if not mol.HasSubstructMatch(thiol_pattern):
        return False, "No alkanethiol with a thiol group (-SH) attached directly to an alkyl group found"
    
    # Further checks can be added here to ensure no unexpected heteroatoms or rings in simple alkane context
    # Such constraints help narrow down false positives from complex molecules
    
    return True, "Contains alkanethiol pattern with -SH group specifically bound to an aliphatic carbon"

# Example usage
print(is_alkanethiol("SCC(CC)C"))  # Example SMILES for 2-Methyl-1-butanethiol