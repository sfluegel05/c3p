"""
Classifies: CHEBI:23044 carotenoid
"""
"""
Classifies: CHEBI:23044 carotenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_carotenoid(smiles: str):
    """
    Determines if a molecule is a carotenoid based on its SMILES string.
    A carotenoid is a tetraterpenoid (C40), derived from psi,psi-carotene by various modifications
    such as hydrogenation, dehydrogenation, cyclization, oxidation, or combinations thereof.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carotenoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 35 or c_count > 45:
        return False, f"Molecule has {c_count} carbon atoms, expected around 40 for carotenoids"

    # Check for conjugated double bonds (polyene chain)
    # Count the number of conjugated double bonds
    num_conj_double_bonds = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE and bond.GetIsConjugated():
            num_conj_double_bonds += 1

    if num_conj_double_bonds < 7:
        return False, f"Only {num_conj_double_bonds} conjugated double bonds found, expected at least 7"

    # Check for isoprene units
    # Isoprene unit SMARTS pattern: C[C]=C[C]C
    isoprene_pattern = Chem.MolFromSmarts("C(C)=C(C)C")
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    if len(isoprene_matches) < 6:
        return False, f"Found {len(isoprene_matches)} isoprene units, expected at least 6"

    # Check for rings (cyclization)
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings > 6:
        return False, f"Molecule has {num_rings} rings, which is more than typically found in carotenoids"

    # Optional: Check for common functional groups (hydroxyl, keto, epoxide)
    # But carotenoids can lack these groups, so we don't exclude molecules without them

    return True, "Molecule matches criteria for a carotenoid"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:23044',
                          'name': 'carotenoid',
                          'definition': 'One of a class of tetraterpenoids (C40), formally derived from the acyclic parent, psi,psi-carotene by hydrogenation, dehydrogenation, cyclization, oxidation, or combination of these processes. This class includes carotenes, xanthophylls and certain compounds that arise from rearrangement of the skeleton of psi,psi-carotene or by loss of part of this structure. Retinoids are excluded.',
                          'parents': ['CHEBI:26859']},
    'config': {   'llm_model_name': 'lbl/claude-sonnet',
                  'f1_threshold': 0.8,
                  'max_attempts': 5,
                  'max_positive_instances': None,
                  'max_positive_to_test': None,
                  'max_negative_to_test': None,
                  'max_positive_in_prompt': 50,
                  'max_negative_in_prompt': 20,
                  'max_instances_in_prompt': 100,
                  'test_proportion': 0.1},
    'message': None,
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': None,
    'num_false_positives': None,
    'num_true_negatives': None,
    'num_false_negatives': None,
    'num_negatives': None,
    'precision': None,
    'recall': None,
    'f1': None,
    'accuracy': None }