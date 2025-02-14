"""
Classifies: CHEBI:80291 aliphatic nitrile
"""
"""
Classifies: CHEBI:33569 aliphatic nitrile
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_aliphatic_nitrile(smiles: str):
    """
    Determines if a molecule is an aliphatic nitrile based on its SMILES string.
    An aliphatic nitrile is a compound containing a nitrile (-C≡N) group attached
    to an aliphatic (non-aromatic) carbon chain.

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

    # Look for nitrile group pattern (C-C#N)
    nitrile_pattern = Chem.MolFromSmarts("[C-;!#1]C#N")
    nitrile_matches = mol.GetSubstructMatches(nitrile_pattern)

    if not nitrile_matches:
        return False, "No nitrile group found"

    # Check if the carbon attached to the nitrile is aliphatic
    for match in nitrile_matches:
        c_idx = match[0]  # Index of the carbon atom attached to the nitrile
        if mol.GetAtomWithIdx(c_idx).GetIsAromatic():
            return False, "Nitrile group is attached to an aromatic carbon"

    # Check if the molecule has any aromatic rings
    if mol.GetAromaticRings():
        return False, "Molecule contains aromatic rings"

    return True, "Molecule contains a nitrile group attached to an aliphatic carbon chain"