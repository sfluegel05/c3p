"""
Classifies: CHEBI:59238 cyclic fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cyclic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a cyclic fatty acid based on its SMILES string.
    A cyclic fatty acid contains a ring of atoms in its structure and retains characteristics of fatty acids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyclic fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify ring structures
    has_ring = mol.GetRingInfo().NumRings() > 0
    if not has_ring:
        return False, "No ring structure found"

    # Look for carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for long hydrocarbon chains (consider 6+ contiguous carbons as 'long')
    chain_pattern = Chem.MolFromSmarts("[CH2]~[CH2]~[CH2]~[CH2]~[CH2]")
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No long hydrocarbon chains found"

    return True, "Contains ring structure with fatty acid traits (carboxylic acid group and long hydrocarbon chain)"