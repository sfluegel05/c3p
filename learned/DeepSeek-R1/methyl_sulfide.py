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

    # Check for aliphatic sulfide (thioether) with at least one methyl group attached to sulfur
    # Sulfur must have exactly two single bonds (SX2) and not be aromatic
    # At least one substituent is a methyl group (CH3)
    sulfide_pattern = MolFromSmarts('[SX2;!a]([CH3])[#6]')
    if mol.HasSubstructMatch(sulfide_pattern):
        return True, "Aliphatic sulfide with methyl group attached to sulfur"

    # Check for dimethyl sulfide case (both substituents are methyl)
    dimethyl_pattern = MolFromSmarts('[SX2;!a]([CH3])([CH3])')
    if mol.HasSubstructMatch(dimethyl_pattern):
        return True, "Dimethyl sulfide (both groups are methyl)"

    return False, "No aliphatic sulfide with methyl group found"