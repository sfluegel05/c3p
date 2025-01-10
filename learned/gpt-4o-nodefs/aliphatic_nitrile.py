"""
Classifies: CHEBI:80291 aliphatic nitrile
"""
from rdkit import Chem

def is_aliphatic_nitrile(smiles: str):
    """
    Determines if a molecule is an aliphatic nitrile based on its SMILES string.
    An aliphatic nitrile contains a nitrile group (C#N) connected to an aliphatic (non-aromatic)
    carbon chain or environment.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an aliphatic nitrile, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for nitrile group pattern (-C#N)
    nitrile_pattern = Chem.MolFromSmarts("[CX2]#N")
    nitrile_matches = mol.GetSubstructMatches(nitrile_pattern)
    
    if not nitrile_matches:
        return False, "No nitrile group found"

    # Ensuring nitrile group is connected to an aliphatic chain
    for match in nitrile_matches:
        nitrile_carbon = mol.GetAtomWithIdx(match[0])

        # Consider the carbon attached to the nitrile group
        for neighbor in nitrile_carbon.GetNeighbors():
            # Check if neighbor is a non-aromatic carbon
            if neighbor.GetAtomicNum() == 6 and not neighbor.GetIsAromatic():
                # Further check if this neighbor is part of a non-aromatic carbon chain
                if all(not n.GetIsAromatic() and n.GetAtomicNum() == 6 for n in neighbor.GetNeighbors() if n.GetIdx() != nitrile_carbon.GetIdx()):
                    return True, "Contains an aliphatic nitrile group"

    return False, "Nitrile group is not connected to an aliphatic chain"