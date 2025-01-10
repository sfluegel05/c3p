"""
Classifies: CHEBI:15693 aldose
"""
from rdkit import Chem

def is_aldose(smiles: str):
    """
    Determines if a molecule is an aldose based on its SMILES string.
    An aldose is a sugar containing an aldehydic parent group and its potential cyclic forms.

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

    # Check if it's an open chain form with an aldehyde group.
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H](O)=O")
    if mol.HasSubstructMatch(aldehyde_pattern):
        found_aldehyde = True
        reason = "Contains an aldehyde group."
    else:
        found_aldehyde = False
        reason = "No free aldehyde group found."

    # Check for cyclic forms: furanoses (5-membered rings) and pyranoses (6-membered rings) with oxygen
    ring_info = mol.GetRingInfo()
    furanose_ring = any(len(ring) == 5 and any(mol.GetAtomWithIdx(idx).GetSymbol() == 'O' for idx in ring) for ring in ring_info.AtomRings())
    pyranose_ring = any(len(ring) == 6 and any(mol.GetAtomWithIdx(idx).GetSymbol() == 'O' for idx in ring) for ring in ring_info.AtomRings())

    if furanose_ring:
        found_cyclic = True
        reason += " Contains a furanose (5-membered) ring."
    elif pyranose_ring:
        found_cyclic = True
        reason += " Contains a pyranose (6-membered) ring."
    else:
        found_cyclic = False
        reason += " No furanose or pyranose ring detected."

    # Check for multiple hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    min_hydroxyl_groups = 2 if found_aldehyde else 3

    if len(hydroxyl_matches) < min_hydroxyl_groups:
        return False, f"Insufficient hydroxyl groups: found {len(hydroxyl_matches)}, but requires at least {min_hydroxyl_groups}."

    if found_aldehyde or found_cyclic:
        return True, reason

    return False, "Did not match criteria for aldose structure."