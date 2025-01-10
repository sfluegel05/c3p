"""
Classifies: CHEBI:15693 aldose
"""
from rdkit import Chem

def is_aldose(smiles: str):
    """
    Determines if a molecule is an aldose based on its SMILES string.
    An aldose is a sugar containing an aldehyde group and several hydroxyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldose, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for aldehyde presence: R-CHO (R-single atom group in sugar contexts)
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H]=O")
    if not mol.HasSubstructMatch(aldehyde_pattern):
        return False, "No aldehyde group found"

    # Check for multiple hydroxyl groups: -OH
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 2:
        return False, "Not enough hydroxyl groups for an aldose"

    # Check for cyclic hemiacetal formation: ring systems like furanose or pyranose
    # The example SMILES can help inform this needed pattern

    # Hemiacetal generically is [C@]1O[C@@H](O) other hydroxylated carbon chains ...
    ring_found = any(len(ring) == 5 or len(ring) == 6 for ring in mol.GetRingInfo().AtomRings())
    if not ring_found:
        return False, "No suitable cyclic hemiacetal structure found"

    return True, "Contains aldehyde and multiple hydroxyl groups, with possible cyclic hemiacetal structure"