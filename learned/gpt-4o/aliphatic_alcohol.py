"""
Classifies: CHEBI:2571 aliphatic alcohol
"""
from rdkit import Chem

def is_aliphatic_alcohol(smiles: str):
    """
    Determines if a molecule is an aliphatic alcohol based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an aliphatic alcohol, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for hydroxyl group attached to an aliphatic carbon
    # [CX4;H] for sp3 hybridized carbon (saturated), or [CX3;H] for sp2 hybridized (unsaturated, but not aromatic)
    aliphatic_alcohol_pattern = Chem.MolFromSmarts("[CX4;!r][OH] | [CX3;!r][OH]")
    
    if mol.HasSubstructMatch(aliphatic_alcohol_pattern):
        return True, "Molecule contains a hydroxyl group attached to an aliphatic carbon"

    return False, "No aliphatic alcohol pattern found"