"""
Classifies: CHEBI:15693 aldose
"""
from rdkit import Chem

def is_aldose(smiles: str):
    """
    Determines if a molecule is an aldose based on its SMILES string.
    An aldose is a sugar containing an aldehydic group and its cyclic forms known as furanoses and pyranoses.

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

    # Check for an open chain form with an aldehyde group: [CX3H]=O
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H](=O)")
    if mol.HasSubstructMatch(aldehyde_pattern):
        found_aldehyde = True
        reason = "Contains an aldehyde group."
    else:
        found_aldehyde = False
        reason = "No free aldehyde group found."

    # Check for cyclic forms: furanoses (5-membered rings) and pyranoses (6-membered rings) with oxygen
    ring_info = mol.GetRingInfo()
    furanose_ring = False
    pyranose_ring = False
    for ring in ring_info.AtomRings():
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        if len(ring) == 5 and any(atom.GetSymbol() == 'O' for atom in ring_atoms):
            furanose_ring = True
        elif len(ring) == 6 and any(atom.GetSymbol() == 'O' for atom in ring_atoms):
            pyranose_ring = True

    if furanose_ring:
        found_cyclic = True
        reason += " Contains a furanose (5-membered) ring."
    elif pyranose_ring:
        found_cyclic = True
        reason += " Contains a pyranose (6-membered) ring."
    else:
        found_cyclic = False
        reason += " No furanose or pyranose ring detected."

    # Check for sufficient hydroxyl groups, root structure is important
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    min_hydroxyl_groups = 3  # Have multiple hydroxyls in cyclic and acyclic forms

    if len(hydroxyl_matches) < min_hydroxyl_groups:
        return False, f"Insufficient hydroxyl groups: found {len(hydroxyl_matches)}, need at least {min_hydroxyl_groups}."

    if found_aldehyde or found_cyclic:
        return True, reason

    return False, "Did not match criteria for aldose structure."