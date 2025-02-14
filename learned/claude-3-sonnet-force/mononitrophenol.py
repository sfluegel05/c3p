"""
Classifies: CHEBI:39362 mononitrophenol
"""
"""
Classifies: CHEBI:16581 mononitrophenol

A mononitrophenol is a phenol carrying a single nitro substituent at an unspecified position.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_mononitrophenol(smiles: str):
    """
    Determines if a molecule is a mononitrophenol based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mononitrophenol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for phenol ring
    phenol_pattern = Chem.MolFromSmarts("c1ccc(O)cc1")
    if not mol.HasSubstructMatch(phenol_pattern):
        return False, "No phenol ring found"

    # Check for single nitro group
    nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
    nitro_matches = mol.GetSubstructMatches(nitro_pattern)
    if len(nitro_matches) != 1:
        return False, f"Found {len(nitro_matches)} nitro groups, expected 1"

    # Check if nitro group is attached to phenol ring
    nitro_atom_idx = nitro_matches[0][0]
    phenol_ring_atoms = mol.GetAtomWithIdx(nitro_atom_idx).GetNeighbors()
    if not any(mol.GetAtomWithIdx(idx).IsInRingSize(6) and mol.GetAtomWithIdx(idx).GetAtomicNum() == 8 for idx in [x.GetIdx() for x in phenol_ring_atoms]):
        return False, "Nitro group not attached to phenol ring"

    return True, "Contains a phenol ring with a single nitro substituent"