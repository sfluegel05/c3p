"""
Classifies: CHEBI:47908 alkanethiol
"""
"""
Classifies: CHEBI:48477 alkanethiol
An alkanethiol is a compound in which a sulfanyl group, -SH, is attached to an alkyl group.
"""
from rdkit import Chem

def is_alkanethiol(smiles: str):
    """
    Determines if a molecule is an alkanethiol based on its SMILES string.
    An alkanethiol is a compound in which a sulfanyl group, -SH, is attached to an alkyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkanethiol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for at least one sulfanyl group (-SH)
    sulfanyl_pattern = Chem.MolFromSmarts("[SH]")
    sulfanyl_matches = mol.GetSubstructMatches(sulfanyl_pattern)
    if len(sulfanyl_matches) == 0:
        return False, "No sulfanyl group (-SH) found"

    # Check if the sulfanyl group is attached to an alkyl group (carbon chain)
    # We look for a carbon atom connected to the sulfur atom
    for match in sulfanyl_matches:
        sulfur_idx = match[0]
        sulfur_atom = mol.GetAtomWithIdx(sulfur_idx)
        neighbors = sulfur_atom.GetNeighbors()
        for neighbor in neighbors:
            if neighbor.GetAtomicNum() == 6:  # Carbon atom
                # Check if the carbon is part of a simple alkyl chain
                # We allow the carbon to be part of a chain or ring, and we don't restrict the rest of the molecule
                return True, "Sulfanyl group (-SH) attached to an alkyl group"

    return False, "Sulfanyl group (-SH) not attached to an alkyl group"