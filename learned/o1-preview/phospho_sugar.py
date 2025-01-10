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
    A phospho sugar is a monosaccharide containing an alcoholic hydroxy group esterified with phosphoric acid.

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

    # Define SMARTS patterns for pyranose and furanose rings
    pyranose_smarts = "[O;R][C;R][C;R][C;R][C;R][C;R]"  # 6-membered ring with one oxygen and five carbons
    furanose_smarts = "[O;R][C;R][C;R][C;R][C;R]"       # 5-membered ring with one oxygen and four carbons

    pyranose = Chem.MolFromSmarts(pyranose_smarts)
    furanose = Chem.MolFromSmarts(furanose_smarts)

    is_pyranose = mol.HasSubstructMatch(pyranose)
    is_furanose = mol.HasSubstructMatch(furanose)

    if not (is_pyranose or is_furanose):
        return False, "No monosaccharide ring (pyranose or furanose) found"

    # Define SMARTS pattern for phosphate ester group attached via oxygen
    phosphate_ester_smarts = "[O]-[P](=O)([O])[O]"  # Ester linkage to phosphate group
    phosphate_ester = Chem.MolFromSmarts(phosphate_ester_smarts)
    phosphate_matches = mol.GetSubstructMatches(phosphate_ester)

    if len(phosphate_matches) == 0:
        return False, "No phosphate ester group found"

    # Verify that phosphate is attached to the sugar ring
    # For each phosphate match, check if the oxygen is connected to the sugar ring
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    found_phosphate_on_sugar = False
    for match in phosphate_matches:
        # The first atom in the match is the oxygen connected to phosphate
        oxygen_idx = match[0]
        # Check if this oxygen is connected to a ring atom
        for ring in atom_rings:
            if oxygen_idx in ring:
                found_phosphate_on_sugar = True
                break
        if found_phosphate_on_sugar:
            break

    if not found_phosphate_on_sugar:
        return False, "Phosphate ester group not attached to sugar ring"

    return True, "Contains monosaccharide ring with phosphate ester group attached to hydroxyl group"

__metadata__ = {
    'chemical_class': {
        'name': 'phospho sugar',
        'definition': 'Any monosaccharide containing an alcoholic hydroxy group esterified with phosphoric acid.'
    }
}