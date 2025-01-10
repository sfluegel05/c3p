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

    # Define the SMARTS pattern for a thiol group: sulfur connected to carbon and with one hydrogen
    thiol_pattern = Chem.MolFromSmarts("[SX2H1][#6]")  # sulfur with single hydrogen bound to carbon
    
    if not mol.HasSubstructMatch(thiol_pattern):
        return False, "No alkanethiol with thiol group (-SH) attached to a carbon group found"
    
    return True, "Contains alkanethiol pattern with -SH group specifically bound to a carbon group"

# Example usage
print(is_alkanethiol("SCC(CC)C"))  # Example SMILES for 2-Methyl-1-butanethiol