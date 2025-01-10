"""
Classifies: CHEBI:61703 nonclassic icosanoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nonclassic_icosanoid(smiles: str):
    """
    Determines if a molecule is a nonclassic icosanoid based on its SMILES string.
    A nonclassic icosanoid is any biologically active signalling molecule made by oxygenation of C20 fatty acids
    other than the classic icosanoids (the leukotrienes and the prostanoids).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nonclassic icosanoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count the number of carbon atoms
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < 18 or num_carbons > 24:
        return False, f"Molecule has {num_carbons} carbons, expected between 18 and 24 carbons"

    # Check for terminal carboxylic acid group (-COOH) or derivative
    carboxylic_acid_pattern = Chem.MolFromSmarts('C(=O)[O;H1,H0]')
    ester_pattern = Chem.MolFromSmarts('C(=O)O[C]')
    if not mol.HasSubstructMatch(carboxylic_acid_pattern) and not mol.HasSubstructMatch(ester_pattern):
        return False, "No terminal carboxylic acid group or ester found"

    # Check for oxygenated functional groups (hydroxyl, epoxide, peroxide, ether, ketone)
    oxy_funcs = ['[OX2H]', '[OX2]', '[OX1]', 'C=O']
    has_oxygenation = False
    for func in oxy_funcs:
        pattern = Chem.MolFromSmarts(func)
        if mol.HasSubstructMatch(pattern):
            has_oxygenation = True
            break
    if not has_oxygenation:
        return False, "No oxygenated functional groups found"

    # Ensure there are additional oxygen atoms besides carboxylic acid or ester
    num_oxygen_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    num_oxygen_required = 3  # At least 3 oxygen atoms including carboxyl
    if num_oxygen_atoms < num_oxygen_required:
        return False, f"Molecule has {num_oxygen_atoms} oxygen atoms, expected at least {num_oxygen_required}"

    # Check that molecule is mostly linear (no large rings except possible epoxides)
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings > 0:
        # Allow only 3-membered (epoxide) rings
        for ring in ring_info.AtomRings():
            if len(ring) > 3:
                return False, "Molecule contains rings larger than epoxides"

    # Exclude prostanoids (contains cyclopentane ring fused into the chain)
    prostanoid_pattern = Chem.MolFromSmarts('C1CCCC1')
    if mol.HasSubstructMatch(prostanoid_pattern):
        return False, "Molecule matches prostanoid pattern (possible prostanoid)"

    # Exclude leukotrienes (specific conjugated triene system)
    leukotriene_pattern = Chem.MolFromSmarts('C=CC=CC=CC=CC(=O)')
    if mol.HasSubstructMatch(leukotriene_pattern):
        return False, "Molecule matches leukotriene pattern (possible leukotriene)"

    # Check for long hydrocarbon chain
    chain_pattern = Chem.MolFromSmarts('C' + ('C' * 10))
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "Molecule does not contain a long hydrocarbon chain"

    return True, "Molecule is a nonclassic icosanoid"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'nonclassic icosanoid',
        'definition': 'Any biologically active signalling molecule made by oxygenation of C20 fatty acids other than the classic icosanoids (the leukotrienes and the prostanoids).',
        'parents': []
    }
}