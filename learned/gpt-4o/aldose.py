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

    # Check for aldehyde group or potential cyclic form
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H]=O")
    if mol.HasSubstructMatch(aldehyde_pattern):
        found_aldehyde = True
        reason = "Found free aldehyde group"
    else:
        found_aldehyde = False
        reason = "No free aldehyde group found,"

    # Check cyclic hemiacetal forms: look for furanose or pyranose structures
    # These are represented by 5- and 6- atom rings respectively, often involving oxygen
    ring_info = mol.GetRingInfo()
    ring_size_furanose = 5
    ring_size_pyranose = 6
    
    furanose_ring = any(len(ring) == ring_size_furanose and Chem.MolFromSmiles('O') in [mol.GetAtomWithIdx(idx) for idx in ring] for ring in ring_info.AtomRings())
    pyranose_ring = any(len(ring) == ring_size_pyranose and Chem.MolFromSmiles('O') in [mol.GetAtomWithIdx(idx) for idx in ring] for ring in ring_info.AtomRings())
    
    if furanose_ring or pyranose_ring:
        found_cyclic = True
        reason += " and potential furanose or pyranose ring structure found."
    else:
        found_cyclic = False
        reason += " no furanose or pyranose ring structure detected."

    # Check for multiple hydroxyl groups appropriately
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    min_hydroxyl_groups = 2 if found_aldehyde else 3
    if len(hydroxyl_matches) < min_hydroxyl_groups:
        return False, "Insufficient hydroxyl groups: requires at least " + str(min_hydroxyl_groups)

    if found_aldehyde or found_cyclic:
        return True, reason

    return False, "Did not match criteria for aldose structure."