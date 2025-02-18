"""
Classifies: CHEBI:47908 alkanethiol
"""
"""
Classifies: CHEBI:46643 alkanethiol
"""
from rdkit import Chem

def is_alkanethiol(smiles: str):
    """
    Determines if a molecule is an alkanethiol based on its SMILES string.
    An alkanethiol has at least one sulfanyl group (-SH) attached to an alkyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkanethiol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for at least one -SH group attached to a non-aromatic carbon (alkyl)
    pattern = Chem.MolFromSmarts('[C;!a]-[SH]')
    if mol.HasSubstructMatch(pattern):
        return True, "Contains a sulfanyl group (-SH) attached to an alkyl group"
    else:
        return False, "No sulfanyl group attached to an alkyl group"