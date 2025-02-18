"""
Classifies: CHEBI:86315 methyl sulfide
"""
"""
Classifies: CHEBI:39194 methyl sulfide
"""
from rdkit import Chem
from rdkit.Chem import MolFromSmiles, MolFromSmarts

def is_methyl_sulfide(smiles: str):
    """
    Determines if a molecule is a methyl sulfide based on its SMILES string.
    A methyl sulfide is an aliphatic sulfide where at least one organyl group attached to sulfur is methyl.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a methyl sulfide, False otherwise
        str: Reason for classification
    """
    mol = MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Aliphatic sulfide (non-aromatic sulfur) with at least one methyl group attached
    # SMARTS: [SX2;!a] connected to [CH3;!a] and any carbon (including aromatic)
    methyl_sulfide_pattern = MolFromSmarts('[SX2;!a]([CH3;!a])[#6]')
    if mol.HasSubstructMatch(methyl_sulfide_pattern):
        return True, "Aliphatic sulfide with methyl group attached to sulfur"

    # Check for dimethyl sulfide (both groups are methyl)
    dimethyl_pattern = MolFromSmarts('[SX2;!a]([CH3;!a])([CH3;!a])')
    if mol.HasSubstructMatch(dimethyl_pattern):
        return True, "Dimethyl sulfide (both groups are methyl)"

    return False, "No aliphatic sulfide with methyl group found"