"""
Classifies: CHEBI:1624 3-oxo-5beta-steroid
"""
"""
Classifies: CHEBI:XXXXX 3-oxo-5beta-steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_oxo_5beta_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5beta-steroid based on its SMILES string.
    A 3-oxo-5beta-steroid is a steroid with a ketone at position 3 and beta-configuration at position 5.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-5beta-steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Assign stereochemistry
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True, flagPossibleStereoCenters=True)

    # Define the steroid nucleus SMARTS pattern
    steroid_core_smarts = """
    [#6]1CC[C@H]2C[C@H]3C[C@@H](C)C[C@H]4CC(=O)CC[C@]34C[C@H]2C1
    """
    steroid_core = Chem.MolFromSmarts(steroid_core_smarts)
    if steroid_core is None:
        return False, "Invalid steroid core SMARTS pattern"

    # Check for steroid core match
    matches = mol.GetSubstructMatches(steroid_core)
    if not matches:
        return False, "Steroid core not found"

    # Define SMARTS pattern for ketone at position 3
    ketone_smarts = "[C;R1]=O"
    ketone_pattern = Chem.MolFromSmarts(ketone_smarts)
    if ketone_pattern is None:
        return False, "Invalid ketone SMARTS pattern"

    # Check for ketone at position 3
    ketone_matches = []
    for match in mol.GetSubstructMatches(ketone_pattern):
        atom = mol.GetAtomWithIdx(match[0])
        if atom.IsInRing() and atom.GetIdx() in matches[0]:
            # Assuming atom indices correspond to position 3
            ketone_matches.append(match)
            break
    if not ketone_matches:
        return False, "Ketone group at position 3 not found"

    # Identify chiral center at position 5
    # Position 5 corresponds to a specific atom in the steroid nucleus
    # We need to map the SMARTS pattern to actual atom indices
    steroid_match = matches[0]
    # The index of position 5 in the SMARTS pattern (adjusted for zero-based indexing)
    position5_idx_in_smarts = 4  # Adjust this index based on the SMARTS pattern
    position5_atom_idx = steroid_match[position5_idx_in_smarts]

    position5_atom = mol.GetAtomWithIdx(position5_atom_idx)
    if position5_atom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
        return False, "Chiral center at position 5 not found"

    # Check stereochemistry at position 5
    # Beta-configuration corresponds to 'S' absolute configuration for position 5
    cip_code = position5_atom.GetProp('_CIPCode') if position5_atom.HasProp('_CIPCode') else None
    if cip_code != 'S':
        return False, f"Position 5 is not in beta-configuration (found {cip_code})"

    return True, "Molecule is a 3-oxo-5beta-steroid"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:XXXXX',
        'name': '3-oxo-5beta-steroid',
        'definition': "Any 3-oxo steroid that has beta- configuration at position 5.",
        'parents': []
    },
    'config': {},
    'message': None,
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
}