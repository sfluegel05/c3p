"""
Classifies: CHEBI:23044 carotenoid
"""
"""
Classifies: CHEBI:27687 Carotenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_carotenoid(smiles: str):
    """
    Determines if a molecule is a carotenoid based on its SMILES string.
    A carotenoid is a tetraterpenoid (C40) derived from the acyclic parent psi,psi-carotene
    by hydrogenation, dehydrogenation, cyclization, oxidation, or a combination of these processes.
    This class includes carotenes, xanthophylls, and certain compounds arising from
    rearrangement or loss of part of the psi,psi-carotene structure.

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

    # Count carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 18 or c_count > 40:
        return False, "Carbon atom count outside typical carotenoid range"

    # Calculate effective backbone length
    backbone_length = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            backbone_length += 2
        elif bond.GetBondType() == Chem.BondType.SINGLE:
            backbone_length += 1

    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) == 6:  # Cyclohexene or cyclohexadiene ring
            backbone_length += 6

    if backbone_length < 18 or backbone_length > 24:
        return False, "Backbone length outside typical carotenoid range"

    # Check for common functional groups (hydroxyl, keto, epoxy, glycosidic)
    functional_group_patterns = [
        Chem.MolFromSmarts("O"),  # Hydroxyl
        Chem.MolFromSmarts("O=C"),  # Keto
        Chem.MolFromSmarts("C1OC1"),  # Epoxy
        Chem.MolFromSmarts("O[C@H]1[C@H]([C@@H]([C@H](O[C@@H]2[C@@H]([C@H]([C@@H]([C@@H]2O)O)O)O)O1)O"),  # Glycosidic
    ]
    has_functional_groups = any(mol.HasSubstructMatch(pattern) for pattern in functional_group_patterns)

    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400 or mol_wt > 700:
        return False, "Molecular weight outside typical carotenoid range"

    # Classify as carotenoid if it has an appropriate backbone length and functional groups
    if has_functional_groups:
        return True, "Contains carotenoid-like backbone and functional groups"
    else:
        return False, "No carotenoid-like features found"