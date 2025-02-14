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
    if c_count < 35:
        return False, "Too few carbon atoms for a carotenoid"

    # Check for linear carbon backbone (at least 10 contiguous C=C bonds)
    backbone_pattern = Chem.MolFromSmarts("[C]=[C]=[C]=[C]=[C]=[C]=[C]=[C]=[C]=[C]")
    backbone_matches = mol.GetSubstructMatches(backbone_pattern)
    if not backbone_matches:
        return False, "No linear carbon backbone found"

    # Check for common functional groups (hydroxyl, keto, epoxy, glycosidic)
    functional_group_patterns = [
        Chem.MolFromSmarts("O"),  # Hydroxyl
        Chem.MolFromSmarts("O=C"),  # Keto
        Chem.MolFromSmarts("C1OC1"),  # Epoxy
        Chem.MolFromSmarts("O[C@H]1[C@H]([C@@H]([C@H](O[C@@H]2[C@@H]([C@H]([C@@H]([C@@H]2O)O)O)O)O1)O"),  # Glycosidic
    ]
    has_functional_groups = any(mol.HasSubstructMatch(pattern) for pattern in functional_group_patterns)

    # Check for common ring systems (cyclohexene, cyclohexadiene)
    ring_patterns = [
        Chem.MolFromSmarts("C1=CCCCC1"),  # Cyclohexene
        Chem.MolFromSmarts("C1=CC=CCC1"),  # Cyclohexadiene
    ]
    has_ring_systems = any(mol.HasSubstructMatch(pattern) for pattern in ring_patterns)

    # Classify as carotenoid if it has a linear backbone, functional groups, and/or ring systems
    if has_functional_groups or has_ring_systems:
        return True, "Contains linear carbon backbone and carotenoid-like features"
    else:
        return False, "No carotenoid-like features found"