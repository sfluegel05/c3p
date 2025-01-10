"""
Classifies: CHEBI:80291 aliphatic nitrile
"""
"""
Classifies: aliphatic nitrile
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_aliphatic_nitrile(smiles: str):
    """
    Determines if a molecule is an aliphatic nitrile based on its SMILES string.
    An aliphatic nitrile is any nitrile derived from an aliphatic compound.

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
    
    # Define SMARTS pattern for nitrile group attached to aliphatic carbon
    aliphatic_nitrile_pattern = Chem.MolFromSmarts("[C;!a]#N")
    matches = mol.GetSubstructMatches(aliphatic_nitrile_pattern)
    
    if not matches:
        return False, "No aliphatic nitrile group found"
    
    # Confirm that nitrile group is not attached to any aromatic system
    # Although the pattern ensures the carbon is aliphatic, we check for entire molecule
    # being derived from an aliphatic compound (no aromatic rings)
    is_aliphatic = True
    for bond in mol.GetBonds():
        if bond.IsAromatic():
            is_aliphatic = False
            break

    if not is_aliphatic:
        return True, "Contains aliphatic nitrile group, but molecule has aromatic rings"
    else:
        return True, "Contains aliphatic nitrile group and molecule is aliphatic"

__metadata__ = {   'chemical_class': {   'name': 'aliphatic nitrile',
                              'definition': 'Any nitrile derived from an aliphatic compound.'},
        'config': {   'llm_model_name': 'lbl/claude-sonnet',
                      'f1_threshold': 0.8}}