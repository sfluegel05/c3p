"""
Classifies: CHEBI:33447 phospho sugar
"""
"""
Classifies: phospho sugar
"""
from rdkit import Chem

def is_phospho_sugar(smiles: str):
    """
    Determines if a molecule is a phospho sugar based on its SMILES string.
    A phospho sugar is any monosaccharide containing an alcoholic hydroxy group esterified with phosphoric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phospho sugar, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for monosaccharide ring (5 or 6-membered ring with one oxygen)
    monosaccharide_ring_smarts = "[O;R][C;R][C;R][C;R][C;R]"  # 5-membered ring with one oxygen
    monosaccharide_ring = Chem.MolFromSmarts(monosaccharide_ring_smarts)
    ring_matches_5 = mol.GetSubstructMatches(monosaccharide_ring)

    monosaccharide_ring_smarts_6 = "[O;R][C;R][C;R][C;R][C;R][C;R]"  # 6-membered ring with one oxygen
    monosaccharide_ring_6 = Chem.MolFromSmarts(monosaccharide_ring_smarts_6)
    ring_matches_6 = mol.GetSubstructMatches(monosaccharide_ring_6)

    # Combine ring matches
    ring_matches = list(ring_matches_5) + list(ring_matches_6)

    # Define SMARTS pattern for acyclic monosaccharide (chain of carbons with hydroxyl groups)
    acyclic_sugar_smarts = "[CH2][CH](O)[CH](O)[CH](O)"  # Simplified pattern for linear tetrose
    acyclic_sugar = Chem.MolFromSmarts(acyclic_sugar_smarts)
    acyclic_matches = mol.GetSubstructMatches(acyclic_sugar)

    if len(ring_matches) == 0 and len(acyclic_matches) == 0:
        return False, "No monosaccharide unit found"

    # Define SMARTS pattern for phosphate ester group attached via oxygen
    phosphate_ester_smarts = "[O]-P(=O)([O])[O]"  # Ester linkage to phosphate group
    phosphate_ester = Chem.MolFromSmarts(phosphate_ester_smarts)
    phosphate_matches = mol.GetSubstructMatches(phosphate_ester)

    if len(phosphate_matches) == 0:
        return False, "No phosphate ester group found"

    # Collect atoms in monosaccharide units
    sugar_atoms = set()
    for match in ring_matches + acyclic_matches:
        sugar_atoms.update(match)

    # Check if any phosphate ester is attached to the monosaccharide
    found_phosphate_on_sugar = False
    for match in phosphate_matches:
        # The first atom in the match is the oxygen connected to phosphate
        oxygen_idx = match[0]
        # Check if this oxygen is connected to a sugar atom
        for neighbor in mol.GetAtomWithIdx(oxygen_idx).GetNeighbors():
            if neighbor.GetIdx() in sugar_atoms:
                found_phosphate_on_sugar = True
                break
        if found_phosphate_on_sugar:
            break

    if not found_phosphate_on_sugar:
        return False, "Phosphate ester group not attached to monosaccharide unit"

    return True, "Contains monosaccharide unit with phosphate ester group attached to hydroxyl group"

__metadata__ = {
    'chemical_class': {
        'name': 'phospho sugar',
        'definition': 'Any monosaccharide containing an alcoholic hydroxy group esterified with phosphoric acid.'
    }
}